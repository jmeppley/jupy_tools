def find_overlaps(bounds, other_bounds, padding=0):
    """
    Look for an overlap between the given rectangle "bounds" and
    all the rectangles in "other_bounds".
    
    Return first overlapping box
    """
    for b in other_bounds:
        if bb_overlaps(bounds, b, padding=padding):
            return b

    return None


def bb_overlaps(box0, box1, padding=0):
    """
    Look for an overlap between the given rectangles
    
    return True or False
    """
    ovl_y = box0.y0 < (box1.y1 + padding) and box1.y0 < (box0.y1 + padding)
    ovl_x = box0.x0 < (box1.x1 + padding) and box1.x0 < (box0.x1 + padding)
    return ovl_y and ovl_x


def _annotate_point(ax, annot_tuple, xytext=None, **kwargs):
    """
    Use ax to annotate a given tuple.
    
    annot_tuple should be (text, (x,y)) 
    but can have additional postional or keyword args
    EG: (text, (x,y), {'rotation': 45})
    
    Extra args will be merged with kwargs. This allows for global settings to be 
    overridden by a sigle label.
    """
    text, pos = annot_tuple[:2]

    if len(annot_tuple) > 2:
        # shallow copy keyword arguments
        annot_kwargs = dict(kwargs)
        annot_kwargs.update(annot_tuple[2])
    else:
        annot_kwargs = kwargs

    return ax.annotate(text, pos, xytext=xytext, **annot_kwargs)


def annotate_points(
    ax,
    annotations,
    padding=0,
    debug=False,
    axis=None,
    arrowprops={
        "width": 2,
        "headlength": 2,
        "headwidth": 2,
        "shrink": 0,
        "color": "grey",
        "alpha": 0.66,
    },
    **kwargs,
):
    """
    Draws non-overlapping point labels

    Parameters:
        ax: the matplotlib axes object for drawing
        annotations: iterable of ("text", (x, y)) nested tuples
        debug: True -> print lots of messages
        axis: restrict movement to 'x' or 'y'. None (defaul): choose shortest
        arrowprops: what the lines pointing back to the points look like
            see: https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.annotate.html?highlight=annotate#matplotlib.axes.Axes.annotate
        **kwargs: all other params passed on to ax.annotate
        
        
    Outline:
     First, text is drawn naively (in memory only) to get the full size of the label
     Then, each annotation is deleted and re-added somewhere else if it overlaps with something
    """

    # simple logging mechanism
    if debug:

        def printd(msg):
            print(msg)

    else:

        def printd(msg):
            pass

    # First, place all labels naively
    # (we don't use the arrowprops here, because everthing is still at the original position)
    annots = {
        _annotate_point(ax, a_tuple, **kwargs): a_tuple for a_tuple in annotations
    }

    # draw plot (in memory) to get text positions
    fig = ax.get_figure()
    fig.canvas.draw()

    # we'll need this later to decode positions
    transf = ax.transData.inverted()

    ## loop over all annotations and adjust any overlaps
    # adjusted positions willneed arrows
    kwargs["arrowprops"] = arrowprops
    # keep track of each label so we can see if others overlap
    saved_positions = []

    # loop
    for annot in annots:
        printd(f"adjusting {annot.get_text()}")
        # get the drawn position in x,y coordinates
        text_bounds = annot.get_window_extent(renderer=fig.canvas.renderer).transformed(
            transf
        )

        # find an overlap
        ovl_box = find_overlaps(text_bounds, saved_positions, padding=padding)
        if ovl_box is not None:
            # there is an overlap, we'll have to move
            # track original position
            x, y = annot.get_position()

            move_count = 0
            while ovl_box is not None:
                printd(f"{text_bounds} overlaps {ovl_box}")

                # if we have to move a label too much, there is something wrong with the code
                move_count += 1
                if move_count > 2 * len(annots):
                    print(annot.get_text())
                    print(annot.get_position())
                    print(text_bounds)
                    print(ovl_box)
                    raise Exception("Too many adjustments!")

                # find shortest move to avoid this overlap
                delta_x, delta_y = move_bounds(text_bounds, ovl_box, padding, axis)
                printd(f"moving {delta_x}, {delta_y}")
                text_bounds.y0 += delta_y
                text_bounds.y1 += delta_y
                text_bounds.x0 += delta_x
                text_bounds.x1 += delta_x

                # check for ovelaps again
                ovl_box = find_overlaps(text_bounds, saved_positions, padding=padding)

            ## Reposition the actual annotation
            # get the original annotation tuple
            annot_tuple = annots[annot]
            # delete old annotation
            annot.remove()
            del annot
            # redraw (this time with the arrow propoerties)
            text_pos = (text_bounds.x0, text_bounds.y0)
            annot = _annotate_point(ax, annot_tuple, xytext=text_pos, **kwargs)
        # save position
        saved_positions.append(text_bounds)


def move_bounds(box, ovl_box, padding=0, axis=None):
    """ 
    given two boxes, return the x and y deltas needed to move the first box to avoid the second.
    
    for simplicity, we chose only one axis, whichever is shorter. Also, only allow positive moves.
    
    can choose axis with axis='y' or axis='x'
    """
    if axis is not None and len(axis.strip()) == 0:
        # empty axis is as good as None
        axis = None

    # calculate moves in each axis
    delta_up = ovl_box.y1 - box.y0
    delta_right = ovl_box.x1 - box.x0

    # find shortest
    min_delta = min(delta_up, delta_right)

    # return shortest
    if axis.lower().startswith("y") or (axis is None and delta_up == min_delta):
        return 0, delta_up + padding
    else:
        return delta_right + padding, 0
