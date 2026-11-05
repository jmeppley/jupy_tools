"""
hit_map_plots needs edl.blastm8 (from jmeppley::py_metagenomics, conda only)
at import time. If it's not installed, a minimal stand-in is used so the rest
of the module can still be tested.
"""
import importlib
import sys
import types
from collections import namedtuple

import numpy
import pandas
import pytest


def _install_fake_edl():
    try:
        import edl.blastm8  # noqa: F401
        return
    except ImportError:
        pass
    edl = types.ModuleType("edl")
    blastm8 = types.ModuleType("edl.blastm8")
    blastm8.SAM = "sam"

    def generate_hits(*args, **kwargs):
        raise NotImplementedError("stand-in for edl.blastm8.generate_hits")
    blastm8.generate_hits = generate_hits
    edl.blastm8 = blastm8
    sys.modules["edl"] = edl
    sys.modules["edl.blastm8"] = blastm8


_install_fake_edl()
from jme.jupy_tools import hit_map_plots  # noqa: E402


Hit = namedtuple('Hit', ['read', 'hit', 'qstart', 'qend', 'hstart', 'hend',
                         'line'])


def make_hit(read, hstart, hend, qstart=1, qend=None, contig='c1', line=None):
    if qend is None:
        qend = qstart + abs(hend - hstart)
    return Hit(read, contig, qstart, qend, hstart, hend, line)


###
# missing dependency


def test_missing_edl_message(monkeypatch):
    monkeypatch.setitem(sys.modules, 'edl', None)
    monkeypatch.setitem(sys.modules, 'edl.blastm8', None)
    monkeypatch.delitem(sys.modules, 'jme.jupy_tools.hit_map_plots')
    with pytest.raises(ModuleNotFoundError, match="py_metagenomics") as error:
        importlib.import_module('jme.jupy_tools.hit_map_plots')
    assert isinstance(error.value.__cause__, ImportError)
    assert error.value.name == 'edl'


###
# range helpers


@pytest.mark.parametrize("r1, r2, expected", [
    ((1, 10), (5, 15), True),
    ((5, 15), (1, 10), True),
    ((1, 10), (10, 20), True),
    ((1, 10), (11, 20), False),
    ((1, 100), (40, 50), True),
])
def test_ranges_intersect(r1, r2, expected):
    assert hit_map_plots.ranges_intersect(r1, r2) == expected


def test_hits_intersect():
    assert hit_map_plots.hits_intersect(make_hit('a', 1, 10),
                                        make_hit('b', 10, 20))
    assert not hit_map_plots.hits_intersect(make_hit('a', 1, 10),
                                            make_hit('b', 11, 20))


def test_hit_list_intersects_range():
    hits = [make_hit('a', 1, 10), make_hit('a', 50, 60)]
    assert hit_map_plots.hit_list_intersects_range(hits, (55, 100))
    assert not hit_map_plots.hit_list_intersects_range(hits, (20, 40))
    assert not hit_map_plots.hit_list_intersects_range([], (20, 40))


def test_get_hit_range():
    hits = [make_hit('a', 50, 60), make_hit('a', 30, 5), make_hit('a', 20, 40)]
    assert hit_map_plots.get_hit_range(hits) == (5, 60)
    assert hit_map_plots.get_hit_range([]) == (None, None)


def test_get_y_pos():
    rows = []
    assert hit_map_plots.get_y_pos((1, 10), rows) == 0
    # doesn't overlap: same row
    assert hit_map_plots.get_y_pos((20, 30), rows) == 0
    # overlaps: new row
    assert hit_map_plots.get_y_pos((5, 25), rows) == 1
    # fits in the first row
    assert hit_map_plots.get_y_pos((12, 18), rows) == 0
    # overlaps both rows
    assert hit_map_plots.get_y_pos((8, 22), rows) == 2
    assert rows == [[(1, 10), (20, 30), (12, 18)], [(5, 25)], [(8, 22)]]


###
# removing overlaps


def test_remove_overlapping_hits():
    hits = [make_hit('a', 1, 100), make_hit('b', 90, 200),
            make_hit('c', 101, 200), make_hit('d', 150, 160)]
    kept = [h.read for h in hit_map_plots.remove_overlapping_hits(hits)]
    assert kept == ['a', 'c']

    # buffer allows small overlaps (b), and hits no longer than the buffer
    #  are always kept (d)
    kept = [h.read for h in
            hit_map_plots.remove_overlapping_hits(hits, buffer=10)]
    assert kept == ['a', 'b', 'd']


def test_remove_overlapping_hits_on_query():
    hits = [make_hit('a', 1, 100, qstart=1, qend=100),
            make_hit('b', 1, 100, qstart=200, qend=101),
            make_hit('c', 1000, 1100, qstart=50, qend=150)]
    kept = [h.read for h in
            hit_map_plots.remove_overlapping_hits(hits, on_hit=False)]
    assert kept == ['a', 'b']


def test_remove_overlapping_hits_short_hits_kept():
    hits = [make_hit('a', 1, 100), make_hit('b', 50, 55)]
    kept = [h.read for h in
            hit_map_plots.remove_overlapping_hits(hits, buffer=10)]
    assert kept == ['a', 'b']


###
# collecting hits by contig


def test_parse_contig_hits(monkeypatch):
    calls = []

    def fake_generate_hits(map_file, format=None, **kwargs):
        calls.append((map_file, format, kwargs))
        yield 'r1', [make_hit('r1', 1, 10, contig='c1'),
                     make_hit('r1', 20, 30, contig='c2'),
                     make_hit('r1', 40, 50, contig='c1')]
        yield 'r2', [make_hit('r2', 5, 15, contig='c2')]

    monkeypatch.setattr(hit_map_plots.blastm8, 'generate_hits',
                        fake_generate_hits)
    contig_hits = hit_map_plots.parse_contig_hits('reads.sam', filter_top_pct=0)
    assert set(contig_hits) == {'c1', 'c2'}
    assert [h.hstart for h in contig_hits['c1']['r1']] == [1, 40]
    assert set(contig_hits['c2']) == {'r1', 'r2'}
    assert calls == [('reads.sam', hit_map_plots.blastm8.SAM,
                      {'filter_top_pct': 0})]


def test_contig_hits_from_table():
    hits = pandas.DataFrame({
        'query': ['r1', 'r1', 'r2', 'r3'],
        'hit': ['c1', 'c1', 'c1', 'c2'],
        'hstart': [1, 50, 10, 1],
        'hend': [20, 70, 30, 10],
    })
    hit_dict = hit_map_plots.contig_hits_from_table(hits)
    assert set(hit_dict) == {'c1', 'c2'}
    assert set(hit_dict['c1']) == {'r1', 'r2'}
    assert [h.hstart for h in hit_dict['c1']['r1']] == [1, 50]
    assert hit_dict['c2']['r3'][0].hend == 10

    # custom column names
    renamed = hits.rename(columns={'query': 'read', 'hit': 'contig'})
    hit_dict = hit_map_plots.contig_hits_from_table(renamed,
                                                    contig_col='contig',
                                                    read_col='read')
    assert set(hit_dict['c1']) == {'r1', 'r2'}


###
# alignment details


REFERENCE = "ACGTACGTAC"


def test_get_diffs():
    # 2 soft clipped, 3 match, 1 insertion, 2 match, 1 deletion, 2 match
    # (with one mismatch in the last block)
    query = "NNACGTTAGA"
    diffs = list(hit_map_plots.get_diffs("1H2S3M1I2M1D2M", query, 1, 8,
                                         REFERENCE))
    assert diffs == [(5, 'T', 3, ''), (8, '', 5, 'C'), (9, 'A', 7, 'T')]


def test_get_diffs_perfect():
    assert list(hit_map_plots.get_diffs("4M", "GTAC", 3, 6, REFERENCE)) == []


def test_get_diffs_bad_cigar():
    with pytest.raises(Exception, match="Unknown cigar character: X"):
        list(hit_map_plots.get_diffs("2M1X", "ACG", 1, 3, REFERENCE))


def test_get_hit_pctids_no_diffs():
    pctids = list(hit_map_plots.get_hit_pctids(1, 100, [], 10))
    assert pctids == [100] * 10


def test_get_hit_pctids_counts_all_diffs():
    diffs = [(0, 'A', 5, 'C'), (0, 'A', 15, 'C'), (0, 'A', 16, 'C')]
    pctids = list(hit_map_plots.get_hit_pctids(1, 100, diffs, 10))
    assert len(pctids) == 10
    assert min(pctids) < 100


def test_get_hit_pctids_positions():
    # 1 diff in region 0 (positions 0-9), 2 in region 1 (10-19), 1 in the last
    diffs = [(0, 'A', 5, 'C'), (0, 'A', 15, 'C'), (0, 'A', 16, 'C'),
             (0, 'A', 95, 'C')]
    pctids = list(hit_map_plots.get_hit_pctids(1, 100, diffs, 10))
    assert pctids == [90, 80] + [100] * 7 + [90]


def test_get_hit_pctids_hit_offset():
    """ regions start at the hit start (diff positions are 0-based) """
    # hit covers 0-based positions 500-599, one diff at the start of each half
    diffs = [(0, 'A', 500, 'C'), (0, 'A', 550, 'C')]
    pctids = list(hit_map_plots.get_hit_pctids(501, 600, diffs, 2))
    assert pctids == [98, 98]
    # a diff on the last base of the first half
    pctids = list(hit_map_plots.get_hit_pctids(501, 600, [(0, 'A', 549, 'C')], 2))
    assert pctids == [98, 100]


def test_get_hit_pctids_fractional_regions():
    """ with a fractional N, the last region is shorter (and scored as such) """
    # 100bp in 2.5 regions: 0-40, 40-80, 80-100
    diffs = [(0, 'A', 10, 'C'), (0, 'A', 90, 'C')]
    pctids = list(hit_map_plots.get_hit_pctids(1, 100, diffs, 2.5))
    assert pctids == [pytest.approx(97.5), 100, pytest.approx(95)]


@pytest.mark.parametrize("N", [3, 3.0, 300 / 100, 0.1 * 3 * 10, 7])
def test_get_hit_pctids_region_count(N):
    """ whole numbers give exactly N regions, even with float rounding """
    pctids = list(hit_map_plots.get_hit_pctids(1, 100, [], N))
    assert len(pctids) == round(N)
    assert pctids == [pytest.approx(100)] * round(N)


def test_get_hit_pctids_less_than_one_region():
    pctids = list(hit_map_plots.get_hit_pctids(1, 100, [(0, 'A', 4, 'C')], 0.5))
    assert pctids == [99]


def test_get_hit_pctids_from_cigar():
    """ diffs from get_diffs line up with the regions """
    reference = "A" * 40
    # 20bp hit at reference positions 11-30, mismatch at the 5th base
    query = "A" * 4 + "C" + "A" * 15
    diffs = hit_map_plots.get_diffs("20M", query, 11, 30, reference)
    assert list(hit_map_plots.get_hit_pctids(11, 30, diffs, 2)) == [90, 100]
    # an insertion (in the second half) counts as a diff, too
    query = "A" * 15 + "G" + "A" * 5
    diffs = hit_map_plots.get_diffs("15M1I5M", query, 11, 30, reference)
    assert list(hit_map_plots.get_hit_pctids(11, 30, diffs, 2)) == [100, 90]


###
# plotting


def test_get_dot_count(ax):
    fig = ax.get_figure()
    pos = ax.get_position()
    x_dots = fig.dpi * fig.get_size_inches()[0] * (pos.x1 - pos.x0)
    y_dots = fig.dpi * fig.get_size_inches()[1] * (pos.y1 - pos.y0)

    assert hit_map_plots.get_dot_count((0, 50), (0, 100), ax) == \
        pytest.approx(x_dots / 2)
    assert hit_map_plots.get_dot_count_x((10, 20), (0, 100), ax) == \
        pytest.approx(x_dots / 10)
    assert hit_map_plots.get_dot_count_y((0, 100), (0, 100), ax) == \
        pytest.approx(y_dots)
    assert hit_map_plots.get_dot_count((0, 25), (0, 100), ax, axis='y') == \
        pytest.approx(y_dots / 4)


@pytest.fixture
def hits_dict():
    """ reads mapped to a 1000bp contig """
    return {
        'r1': [make_hit('r1', 1, 100, qstart=1, qend=100),
               make_hit('r1', 150, 300, qstart=101, qend=250)],
        'r2': [make_hit('r2', 200, 400)],
        'r3': [make_hit('r3', 500, 700)],
        'r4': [make_hit('r4', 900, 1000, qstart=200, qend=100)],
    }


def test_plot_hits(ax, hits_dict):
    positions = hit_map_plots.plot_hits(ax, hits_dict, reference_length=1000)
    # r1 and r2 overlap, so r2 gets bumped up a row
    assert positions == {'r1': 1, 'r2': 2, 'r3': 1, 'r4': 1}
    # reference + one read span per read + one bar per hit
    assert len(ax.patches) == 1 + 4 + 5
    assert list(ax.get_yticks()) == [1, 2]


def test_plot_hits_range(ax, hits_dict):
    positions = hit_map_plots.plot_hits(ax, hits_dict, plot_range=(450, 950))
    assert positions == {'r3': 1, 'r4': 1}
    # reference bar covers the range
    ref_bar = ax.patches[0]
    assert ref_bar.get_x() == 450
    assert ref_bar.get_width() == 501


def test_plot_hits_range_keep_pairs(ax, hits_dict):
    # r1 has one hit in range, but both are drawn with hide_missing_pairs=False
    positions = hit_map_plots.plot_hits(ax, hits_dict, plot_range=(1, 120),
                                        hide_missing_pairs=False)
    assert positions == {'r1': 1}
    assert len(ax.patches) == 1 + 1 + 2


def test_plot_hits_colors(ax, hits_dict):
    from matplotlib.colors import to_rgba
    hit_map_plots.plot_hits(
        ax, hits_dict, reference_length=1000,
        bar_color=lambda hit: 'green' if hit.read == 'r3' else 'orange',
        sort_key=lambda k: -int(k[1:]))
    colors = [p.get_facecolor() for p in ax.patches]
    assert colors.count(to_rgba('green')) == 1
    assert colors.count(to_rgba('orange')) == 4


def test_plot_hits_needs_length(ax, hits_dict):
    with pytest.raises(Exception, match="reference length or plot range"):
        hit_map_plots.plot_hits(ax, hits_dict)


def test_plot_hits_colored_py_pos(ax, hits_dict):
    positions = hit_map_plots.plot_hits_colored_py_pos(ax, hits_dict,
                                                       reference_length=1000)
    assert positions == {'r1': 1, 'r2': 2, 'r3': 1, 'r4': 1}
    # one image per hit
    assert len(ax.images) == 5
    # all on the same color scale (longest read)
    assert {im.get_clim() for im in ax.images} == {(0, 250)}
    assert ax.get_xlim() == (0, 1000)
    assert ax.get_ylim() == (-0.5, 2.5)


def test_plot_hits_colored_py_pos_options(ax, hits_dict):
    positions = hit_map_plots.plot_hits_colored_py_pos(
        ax, hits_dict, plot_range=(150, 750), norm_read_color=True)
    assert set(positions) == {'r1', 'r2', 'r3'}
    assert ax.get_xlim() == (149, 750)
    # each read scaled to its own max
    clims = sorted(im.get_clim()[1] for im in ax.images)
    assert clims == [201, 201, 250]


def test_plot_hits_colored_py_pos_read_lengths(ax, hits_dict):
    hit_map_plots.plot_hits_colored_py_pos(
        ax, hits_dict, reference_length=1000,
        norm_read_color={'r1': 1000, 'r2': 2000, 'r3': 3000, 'r4': 4000})
    clims = sorted(im.get_clim()[1] for im in ax.images)
    assert clims == [1000, 1000, 2000, 3000, 4000]


def test_plot_hits_colored_py_pos_needs_length(ax, hits_dict):
    with pytest.raises(Exception, match="reference length or plot range"):
        hit_map_plots.plot_hits_colored_py_pos(ax, hits_dict)


def _sam_line(read, pos, cigar, seq):
    return "\t".join([read, "0", "c1", str(pos), "60", cigar, "*", "0", "0",
                      seq, "*"])


def test_plot_hits_colored_pctid(ax):
    reference = "ACGT" * 25
    hits_dict = {
        'r1': [Hit('r1', 'c1', 1, 8, 1, 8, _sam_line('r1', 1, '8M', 'ACGTACGT'))],
        # secondary hit without a sequence reuses the primary's
        'r2': [Hit('r2', 'c1', 1, 8, 41, 48,
                   _sam_line('r2', 41, '8M', 'ACGTTCGT')),
               Hit('r2', 'c1', 1, 8, 61, 68, _sam_line('r2', 61, '8M', '*'))],
    }
    positions = hit_map_plots.plot_hits_colored_pctid(ax, hits_dict, reference)
    assert positions == {'r1': 1, 'r2': 1}
    assert len(ax.images) == 3
    assert ax.get_xlim() == (0, 100)

    ax.clear()
    positions = hit_map_plots.plot_hits_colored_pctid(
        ax, hits_dict, reference, plot_range=(30, 50), pctid_range=[50, 100])
    assert positions == {'r2': 1}
    assert len(ax.images) == 1
    assert ax.images[0].get_clim() == (50, 100)
