import numpy
import pandas
import pytest

from jme.jupy_tools import mags


@pytest.fixture
def mag_hits():
    """
    Hits between the contigs of two MAGs:
        hit MAG: h1 (1000bp), h2 (500bp)
        query MAG: qA (400bp), qB (300bp)
    """
    return pandas.DataFrame([
        # query, qlen, qstart, qend, hit, hlen, hstart, hend, mlen, pctid
        ['qA', 400, 1, 50, 'h1', 1000, 901, 950, 50, 80.0],
        ['qA', 400, 1, 400, 'h2', 500, 51, 450, 400, 99.0],
        ['qB', 300, 1, 300, 'h1', 1000, 101, 400, 300, 95.0],
    ], columns=['query', 'qlen', 'qstart', 'qend', 'hit', 'hlen', 'hstart',
                'hend', 'mlen', 'pctid'])


def test_get_contig_starts():
    lens = {'a': 100, 'b': 50, 'c': 10}
    assert mags.get_contig_starts(['a', 'b', 'c'], lens) == \
        {'a': 0, 'b': 600, 'c': 1150}
    assert mags.get_contig_starts(['c', 'a'], lens, buffer=0) == \
        {'c': 0, 'a': 10}
    assert mags.get_contig_starts([], lens) == {}


def test_random_color():
    color = mags.random_color()
    assert color.shape == (3,)
    assert all(0 <= c <= 1 for c in color)


###
# pseudo_mauve


def test_pseudo_mauve(ax, mag_hits):
    min_pctid = mags.pseudo_mauve(mag_hits, ax)
    assert min_pctid == 80.0
    # one line per contig
    assert len(ax.get_lines()) == 4
    # one polygon per hit
    assert len(ax.patches) == 3
    # hits are on the bottom, queries on top
    assert ax.get_ylim() == (-1, 21)
    # x axis covers both genomes:
    #   hits (sorted by length): h1 + 500 + h2 + 500
    assert ax.get_xlim() == (0, 1000 + 500 + 500 + 500)


def test_pseudo_mauve_min_pctid(ax, mag_hits):
    assert mags.pseudo_mauve(mag_hits, ax, min_pctid=50) == 50
    # negative values are relative to the observed minimum
    assert mags.pseudo_mauve(mag_hits, ax, min_pctid=-10) == 70


def test_pseudo_mauve_layout(ax, mag_hits):
    mags.pseudo_mauve(mag_hits, ax, hit_y=5, query_y=10, buffer_x=0,
                      connection_colormap='viridis',
                      contig_colors=['red'])
    assert ax.get_ylim() == (4, 11)
    assert ax.get_xlim() == (0, 1500)
    assert all(l.get_color() == 'red' for l in ax.get_lines())

    # polygons are drawn with the highest pctid last (on top) and colored
    #  on a scale from the min pctid (80) to 100
    first_poly, last_poly = ax.patches[0], ax.patches[-1]
    from matplotlib import pyplot as plt
    assert tuple(last_poly.get_facecolor()[:3]) == \
        pytest.approx(plt.get_cmap('viridis')((99 - 80) / 20)[:3])
    assert tuple(first_poly.get_facecolor()[:3]) == \
        pytest.approx(plt.get_cmap('viridis')(0.0)[:3])


###
# sorted hits and dot plots


def test_get_sorted_hits(mag_hits):
    sorted_hits = mags.get_sorted_hits(mag_hits, other_cols=['pctid'])
    assert list(sorted_hits.columns) == ['hit', 'hstart', 'hend', 'query',
                                         'qstart', 'qend', 'pctid']
    # h1 is longer, so it's first. h2 starts after h1 + 500bp buffer
    h2_hit = sorted_hits[sorted_hits.hit == 'h2'].iloc[0]
    assert (h2_hit.hstart, h2_hit.hend) == (1551, 1950)
    h1_hit = sorted_hits[sorted_hits['query'] == 'qB'].iloc[0]
    assert (h1_hit.hstart, h1_hit.hend) == (101, 400)


def test_get_sorted_hits_query_order(mag_hits):
    sorted_hits = mags.get_sorted_hits(mag_hits, other_cols=['pctid'])
    # qB's best hit is on h1 (first hit contig), qA's best is on h2,
    #  so qB should come first in the query genome
    qB = sorted_hits[sorted_hits['query'] == 'qB'].iloc[0]
    assert qB.qstart == 1
    # qA starts after qB (300bp) + 500bp buffer
    qA = sorted_hits[sorted_hits['query'] == 'qA']
    assert set(qA.qstart) == {801}


@pytest.mark.parametrize("mlen, pctid, best", [
    # longer alignment wins at similar identity
    ((400, 50), (95.0, 99.0), 'h2'),
    ((50, 400), (99.0, 95.0), 'h1'),
    # aligned length * pctid: 100 * 99 < 120 * 85
    ((100, 120), (99.0, 85.0), 'h1'),
])
def test_get_sorted_hits_best_hit(mlen, pctid, best):
    """ each query is placed by its hit with the highest mlen * pctid """
    hits = pandas.DataFrame([
        # qX hits both contigs, qY only hits h1, qZ only hits h2
        ['qX', 500, 1, 400, 'h2', 1000, 1, 400, mlen[0], pctid[0]],
        ['qX', 500, 1, 400, 'h1', 2000, 1, 400, mlen[1], pctid[1]],
        ['qY', 500, 1, 400, 'h1', 2000, 1, 400, 400, 99.0],
        ['qZ', 500, 1, 400, 'h2', 1000, 1, 400, 400, 99.0],
    ], columns=['query', 'qlen', 'qstart', 'qend', 'hit', 'hlen', 'hstart',
                'hend', 'mlen', 'pctid'])
    sorted_hits = mags.get_sorted_hits(hits)
    query_starts = sorted_hits.groupby('query').qstart.min()
    # h1 is longer, so it's first. Queries placed on h1 come first
    if best == 'h1':
        assert query_starts['qX'] < query_starts['qZ']
    else:
        assert query_starts['qX'] > query_starts['qY']


def test_get_sorted_hits_score_column(mag_hits):
    """ an existing score column (eg from BLAST) is not used or changed """
    hits = mag_hits.assign(score=[1000.0, 1.0, 1.0])
    sorted_hits = mags.get_sorted_hits(hits, other_cols=['score'])
    assert list(sorted_hits.score) == [1000.0, 1.0, 1.0]
    qB = sorted_hits[sorted_hits['query'] == 'qB'].iloc[0]
    assert qB.qstart == 1


def test_get_sorted_hits_default_cols(mag_hits):
    sorted_hits = mags.get_sorted_hits(mag_hits)
    assert list(sorted_hits.columns) == ['hit', 'hstart', 'hend', 'query',
                                         'qstart', 'qend']
    assert len(sorted_hits) == 3


def test_get_sorted_hits_does_not_modify_input(mag_hits):
    original = mag_hits.copy()
    mags.get_sorted_hits(mag_hits, other_cols=['pctid'])
    pandas.testing.assert_frame_equal(mag_hits, original)


def test_get_sorted_hits_slice():
    """ a filtered slice of a bigger table works without a pandas warning """
    import warnings
    hits = pandas.DataFrame([
        ['qA', 400, 1, 400, 'h1', 1000, 51, 450, 400, 99.0],
        ['qB', 300, 1, 300, 'h1', 1000, 501, 800, 300, 95.0],
        ['qC', 300, 1, 300, 'h1', 1000, 1, 300, 300, 50.0],
    ], columns=['query', 'qlen', 'qstart', 'qend', 'hit', 'hlen', 'hstart',
                'hend', 'mlen', 'pctid'])
    good_hits = hits[hits.pctid > 90]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sorted_hits = mags.get_sorted_hits(good_hits, other_cols=['pctid'])
    assert list(sorted_hits['query']) == ['qA', 'qB']
    assert list(good_hits.columns) == list(hits.columns)


def test_sorted_dot_plot(ax, mag_hits):
    mags.sorted_dot_plot(mag_hits.copy(), ax=ax, min_pctid=50)
    assert len(ax.get_lines()) == 3


def test_sorted_dot_plot_no_ax(close_figures, mag_hits):
    from matplotlib import pyplot as plt
    mags.sorted_dot_plot(mag_hits.copy(), min_pctid=50, colormap='viridis')
    assert len(plt.gca().get_lines()) == 3


@pytest.mark.xfail(strict=True, raises=NameError, reason=(
    "sorted_dot_plot refers to an undefined 'hits' variable when min_pctid "
    "is None or negative"))
def test_sorted_dot_plot_default_min_pctid(ax, mag_hits):
    mags.sorted_dot_plot(mag_hits.copy(), ax=ax)
    assert len(ax.get_lines()) == 3
