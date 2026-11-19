"""
Draw point labels on a matplotlib plot without letting them overlap.

Labels are first drawn where they belong. Any label that overlaps one placed
before it is moved up or right (whichever is the shorter move) until it is
clear of the others, and an arrow is drawn back to its point.

    from jme.jupy_tools.smart_labels import annotate_points
    annotate_points(ax, [("label 1", (x1, y1)), ("label 2", (x2, y2))])

Each annotation is a ("text", (x, y)) tuple, optionally with a third item: a
dict of keyword arguments for ax.annotate() for that label alone.
"""
import logging

LOGGER = logging.getLogger(__name__)

DEFAULT_ARROWPROPS = {
    "width": 2,
    "headlength": 2,
    "headwidth": 2,
    "shrink": 0,
    "color": "grey",
    "alpha": 0.66,
}


class TooManyAdjustmentsError(RuntimeError):
    """ Raised when a label can't be placed after a reasonable number of moves """


def bb_overlaps(box0, box1, padding=0):
    """
    Return True if two bounding boxes overlap.

    Boxes are anything with x0, x1, y0, and y1 attributes (like
    matplotlib.transforms.Bbox). Boxes closer than padding count as overlapping.
    Boxes that only touch along an edge do not overlap.
    """
    overlaps_y = box0.y0 < (box1.y1 + padding) and box1.y0 < (box0.y1 + padding)
    overlaps_x = box0.x0 < (box1.x1 + padding) and box1.x0 < (box0.x1 + padding)
    return overlaps_y and overlaps_x


def find_overlaps(bounds, other_bounds, padding=0):
    """
    Return the first box in other_bounds that overlaps bounds (or None).

    See bb_overlaps() for what counts as an overlap.
    """
    for other_box in other_bounds:
        if bb_overlaps(bounds, other_box, padding=padding):
            return other_box
    return None


def move_bounds(box, ovl_box, padding=0, axis=None):
    """
    Return the (delta_x, delta_y) that moves box clear of ovl_box.

    Only one axis is moved and only in the positive direction (up or right).
    By default, the axis with the shorter move is used (up wins ties). Set
    axis to 'x' or 'y' to force the direction. The move includes the padding.
    """
    if axis is not None:
        axis = axis.strip().lower()
        if axis == "":
            # empty axis is as good as None
            axis = None
        elif axis[0] not in ("x", "y"):
            raise ValueError(f"axis must be 'x', 'y', or None, not {axis!r}")

    delta_up = ovl_box.y1 - box.y0
    delta_right = ovl_box.x1 - box.x0

    if axis is None:
        move_up = delta_up <= delta_right
    else:
        move_up = axis.startswith("y")

    if move_up:
        return 0, delta_up + padding
    return delta_right + padding, 0


def _annotate_point(ax, annot_tuple, xytext=None, **kwargs):
    """
    Draw a single annotation tuple with ax.annotate() and return the Annotation.

    annot_tuple is (text, (x, y)) or (text, (x, y), {extra kwargs}). The extra
    kwargs override the shared kwargs for this label only.
    """
    text, position = annot_tuple[:2]
    if len(annot_tuple) > 2:
        annot_kwargs = {**kwargs, **annot_tuple[2]}
    else:
        annot_kwargs = kwargs
    return ax.annotate(text, position, xytext=xytext, **annot_kwargs)


def annotate_points(
    ax,
    annotations,
    padding=0,
    axis=None,
    arrowprops=None,
    max_moves=None,
    **kwargs,
):
    """
    Draw point labels, moving any that overlap previously drawn labels.

    Parameters:
        ax: the matplotlib Axes to draw on
        annotations: iterable of ("text", (x, y)) tuples. A third item, a dict
            of ax.annotate() keyword arguments, applies to that label alone.
        padding: minimum gap (in data units) to leave between labels
        axis: restrict movement to 'x' or 'y'. None (default): use the shorter
            move for each overlap.
        arrowprops: style of the arrows pointing from moved labels back to
            their points (see matplotlib's Axes.annotate). Defaults to
            DEFAULT_ARROWPROPS.
        max_moves: give up (raising TooManyAdjustmentsError) if one label has
            to move more than this many times. Defaults to twice the number
            of labels.
        **kwargs: passed on to ax.annotate for every label

    Returns:
        list of the matplotlib Annotation objects (in the input order)

    Labels are placed in order, so earlier labels never move. Each label is
    drawn in memory to measure it, and only re-drawn if it had to move.
    Use logging at DEBUG level for details of each move.
    """
    if arrowprops is None:
        arrowprops = DEFAULT_ARROWPROPS
    annotations = list(annotations)
    if max_moves is None:
        max_moves = 2 * len(annotations)

    # Place every label naively first. There are no arrows yet, because
    # everything is still at its original position.
    drawn_labels = [
        _annotate_point(ax, annot_tuple, **kwargs) for annot_tuple in annotations
    ]

    # draw the figure (in memory) so the labels have sizes
    figure = ax.get_figure()
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    display_to_data = ax.transData.inverted()

    # moved labels get arrows pointing back to their points
    moved_kwargs = dict(kwargs, arrowprops=arrowprops)
    placed_bounds = []
    final_labels = []
    for annot_tuple, label in zip(annotations, drawn_labels):
        text_bounds = label.get_window_extent(renderer=renderer).transformed(
            display_to_data
        )

        ovl_box = find_overlaps(text_bounds, placed_bounds, padding=padding)
        move_count = 0
        total_dx, total_dy = 0, 0
        while ovl_box is not None:
            move_count += 1
            if move_count > max_moves:
                raise TooManyAdjustmentsError(
                    f"Could not place label {label.get_text()!r} at "
                    f"{label.get_position()} after {max_moves} moves. "
                    f"Last overlap: {text_bounds} with {ovl_box}"
                )
            delta_x, delta_y = move_bounds(text_bounds, ovl_box, padding, axis)
            LOGGER.debug("moving %r by (%s, %s) to avoid %s",
                         label.get_text(), delta_x, delta_y, ovl_box)
            text_bounds = text_bounds.translated(delta_x, delta_y)
            total_dx += delta_x
            total_dy += delta_y
            ovl_box = find_overlaps(text_bounds, placed_bounds, padding=padding)

        if move_count > 0:
            # Replace the original label with an arrowed one, shifted by the
            # total move (so any text alignment settings are respected)
            text_x, text_y = label.get_position()
            label.remove()
            label = _annotate_point(
                ax,
                annot_tuple,
                xytext=(text_x + total_dx, text_y + total_dy),
                **moved_kwargs,
            )

        placed_bounds.append(text_bounds)
        final_labels.append(label)

    return final_labels
