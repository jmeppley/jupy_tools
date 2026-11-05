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
    sorted_hits = mags.get_sorted_hits(mag_hits.copy(), other_cols=['pctid'])
    assert list(sorted_hits.columns) == ['hit', 'hstart', 'hend', 'query',
                                         'qstart', 'qend', 'pctid']
    # h1 is longer, so it's first. h2 starts after h1 + 500bp buffer
    h2_hit = sorted_hits[sorted_hits.hit == 'h2'].iloc[0]
    assert (h2_hit.hstart, h2_hit.hend) == (1551, 1950)
    h1_hit = sorted_hits[sorted_hits['query'] == 'qB'].iloc[0]
    assert (h1_hit.hstart, h1_hit.hend) == (101, 400)


@pytest.mark.xfail(strict=True, reason=(
    "get_sorted_hits says it orders queries by their best hit, but it sorts "
    "scores ascending and uses the first (worst) hit"))
def test_get_sorted_hits_query_order(mag_hits):
    sorted_hits = mags.get_sorted_hits(mag_hits.copy(), other_cols=['pctid'])
    # qB's best hit is on h1 (first hit contig), qA's best is on h2,
    #  so qB should come first in the query genome
    qB = sorted_hits[sorted_hits['query'] == 'qB'].iloc[0]
    assert qB.qstart == 1


@pytest.mark.xfail(strict=True, raises=TypeError, reason=(
    "get_sorted_hits defaults other_cols to None, then adds it to a list"))
def test_get_sorted_hits_default_cols(mag_hits):
    sorted_hits = mags.get_sorted_hits(mag_hits.copy())
    assert len(sorted_hits.columns) == 6


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
