from collections import namedtuple
from itertools import combinations

import pytest
from matplotlib.text import Text

from jme.jupy_tools import smart_labels

Box = namedtuple("Box", ["x0", "y0", "x1", "y1"])


@pytest.mark.parametrize("box0, box1, padding, expected", [
    (Box(0, 0, 2, 2), Box(1, 1, 3, 3), 0, True),
    (Box(0, 0, 2, 2), Box(0.5, 0.5, 1, 1), 0, True),  # contained
    (Box(0, 0, 1, 1), Box(2, 2, 3, 3), 0, False),
    (Box(0, 0, 1, 1), Box(0, 2, 1, 3), 0, False),  # same x, different y
    (Box(0, 0, 1, 1), Box(2, 0, 3, 1), 0, False),  # same y, different x
    (Box(0, 0, 1, 1), Box(1, 0, 2, 1), 0, False),  # touching edges
    (Box(0, 0, 1, 1), Box(1.5, 0, 2, 1), 1, True),  # within padding
])
def test_bb_overlaps(box0, box1, padding, expected):
    assert smart_labels.bb_overlaps(box0, box1, padding) is expected
    assert smart_labels.bb_overlaps(box1, box0, padding) is expected


def test_find_overlaps():
    others = [Box(5, 5, 6, 6), Box(0, 0, 1, 1), Box(0.5, 0.5, 2, 2)]
    assert smart_labels.find_overlaps(Box(0.8, 0.8, 1.2, 1.2), others) == others[1]
    assert smart_labels.find_overlaps(Box(3, 3, 4, 4), others) is None
    assert smart_labels.find_overlaps(Box(3, 3, 4, 4), []) is None


def test_move_bounds():
    box = Box(0, 0, 2, 2)
    # shorter to move up
    assert smart_labels.move_bounds(box, Box(1, -1, 5, 1)) == (0, 1)
    # shorter to move right
    assert smart_labels.move_bounds(box, Box(-1, 1, 1, 5)) == (1, 0)
    # padding is added
    assert smart_labels.move_bounds(box, Box(1, -1, 5, 1), padding=0.5) == (0, 1.5)


@pytest.mark.parametrize("axis, expected", [
    (None, (0, 1)),
    ("", (0, 1)),
    ("  ", (0, 1)),
    ("y", (0, 1)),
    ("Y", (0, 1)),
    ("x", (4, 0)),
    ("X", (4, 0)),
])
def test_move_bounds_axis(axis, expected):
    assert smart_labels.move_bounds(Box(0, 0, 2, 2), Box(1, -1, 4, 1), axis=axis) == expected


def test_move_bounds_bad_axis():
    with pytest.raises(ValueError):
        smart_labels.move_bounds(Box(0, 0, 2, 2), Box(1, -1, 4, 1), axis="z")


def _label_bounds(ax, labels):
    """
    get the bounding boxes of the label text (not the arrows) in data
    coordinates
    """
    figure = ax.get_figure()
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    to_data = ax.transData.inverted()
    return [
        Text.get_window_extent(label, renderer=renderer).transformed(to_data)
        for label in labels
    ]


def _setup_axes(ax):
    """ fix the axis limits so labels have a predictable size """
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)


def test_annotate_points_no_overlaps(ax):
    _setup_axes(ax)
    points = [(5, 5), (5.1, 5), (5, 5.1), (5.2, 5.2), (1, 1)]
    annotations = [(f"label {i}", point) for i, point in enumerate(points)]

    labels = smart_labels.annotate_points(ax, annotations)

    assert [label.get_text() for label in labels] == \
           [text for text, _ in annotations]
    # only the current labels are on the axes
    assert len(ax.texts) == len(annotations)
    for box0, box1 in combinations(_label_bounds(ax, labels), 2):
        assert not smart_labels.bb_overlaps(box0, box1)

    # the first and lone labels didn't move and have no arrows
    for index in (0, 4):
        assert labels[index].get_position() == points[index]
        assert labels[index].arrow_patch is None
    # the moved labels have arrows and still point at their points
    for index in (1, 2, 3):
        assert labels[index].xy == points[index]
        assert labels[index].get_position() != points[index]
        assert labels[index].arrow_patch is not None


@pytest.mark.parametrize("axis", ["x", "y"])
def test_annotate_points_axis(ax, axis):
    _setup_axes(ax)
    annotations = [("first", (5, 5)), ("second", (5, 5))]
    labels = smart_labels.annotate_points(ax, annotations, axis=axis)
    x_pos, y_pos = labels[1].get_position()
    if axis == "x":
        assert y_pos == 5
        assert x_pos > 5
    else:
        assert x_pos == 5
        assert y_pos > 5


def test_annotate_points_kwargs(ax):
    _setup_axes(ax)
    annotations = [
        ("first", (5, 5)),
        ("second", (5, 5), {"color": "red"}),
    ]
    labels = smart_labels.annotate_points(
        ax, annotations, color="blue", fontsize=6,
        arrowprops={"arrowstyle": "-"},
    )
    assert labels[0].get_color() == "blue"
    assert labels[1].get_color() == "red"
    assert all(label.get_fontsize() == 6 for label in labels)
    assert labels[0].arrow_patch is None
    assert labels[1].arrow_patch is not None


def test_annotate_points_padding(ax):
    _setup_axes(ax)
    annotations = [("first", (5, 5)), ("second", (5, 5))]
    labels = smart_labels.annotate_points(ax, annotations, padding=1, axis="y")
    box0, box1 = _label_bounds(ax, labels)
    assert box1.y0 >= box0.y1 + 1 - 1e-9


def test_annotate_points_too_many_moves(ax):
    _setup_axes(ax)
    annotations = [(f"label {i}", (5, 5)) for i in range(3)]
    with pytest.raises(smart_labels.TooManyAdjustmentsError):
        smart_labels.annotate_points(ax, annotations, max_moves=1)


def test_annotate_points_does_not_modify_defaults(ax):
    _setup_axes(ax)
    defaults = dict(smart_labels.DEFAULT_ARROWPROPS)
    smart_labels.annotate_points(ax, [("a", (5, 5)), ("b", (5, 5))])
    assert smart_labels.DEFAULT_ARROWPROPS == defaults


def test_annotate_points_empty(ax):
    assert smart_labels.annotate_points(ax, []) == []
