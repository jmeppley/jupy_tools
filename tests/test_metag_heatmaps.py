"""
Tests for jme.jupy_tools.experimental.metag_heatmaps, using small synthetic
tables (the real data lives on the lab servers).
"""
import os

import numpy
import pandas
import pytest

polars = pytest.importorskip("polars")
pytest.importorskip("scipy")

from jme.jupy_tools.experimental import metag_heatmaps  # noqa: E402

V, C = metag_heatmaps.ALOHA_VIRAL, metag_heatmaps.ALOHA_CELLULAR


@pytest.fixture
def md():
    """ 2 depths x 2 filters x 2 replicates """
    return polars.DataFrame(
        [(f"{f}_{d}_{r}", d, r + 1, f)
         for f in (V, C) for d in (25, 500) for r in range(2)],
        schema=['sample', 'depth', 'month', 'filter'], orient='row')


@pytest.fixture
def hits():
    rows = [
        # g1: only in one of the two viral 25m samples
        (f"{V}_25_0", 'g1', 4.0),
        # g2: in both viral 25m samples and one cellular 500m sample
        (f"{V}_25_0", 'g2', 1.0),
        (f"{V}_25_1", 'g2', 3.0),
        (f"{C}_500_0", 'g2', 10.0),
        # g3: cellular only
        (f"{C}_25_0", 'g3', 2.0),
        (f"{C}_25_1", 'g3', 2.0),
        (f"{C}_500_1", 'g3', 6.0),
        # a sample that's not in the metadata
        ("other", 'g1', 100.0),
    ]
    return polars.DataFrame(rows, schema=['sample', 'genome',
                                          'cov_per_bn_reads'], orient='row')


def test_normalize_coverage():
    hits = polars.DataFrame({'sample': ['a', 'a', 'b', 'c'],
                             'genome': ['g1', 'g2', 'g1', 'g1'],
                             'cov': [1.0, 2.0, 3.0, 4.0]})
    out = metag_heatmaps.normalize_coverage(hits, {'a': 1e9, 'b': 2e9})
    # c has no read count
    assert out.height == 3
    assert out['cov_per_bn_reads'].to_list() == [1.0, 2.0, 1.5]


def test_long_from_wide():
    wide = pandas.DataFrame({'g1': [1.0, 0.0], 'g2': [2.0, 3.0]},
                            index=['s1', 's2'])
    long = metag_heatmaps.long_from_wide(wide, ['g1', 'g2'])
    assert long.columns == ['sample', 'genome', 'cov_per_bn_reads']
    assert long.height == 4


def test_aggregate_counts_absent_as_zero(hits, md):
    agg = metag_heatmaps.aggregate(hits, md, 'depth',
                                   genomes=['g1', 'g2', 'g3'])
    assert list(agg.index) == [25, 500]
    assert list(agg.columns) == ['g1', 'g2', 'g3']
    # 4 samples at 25m (2 viral, 2 cellular), g1 in only one of them
    assert agg.loc[25, 'g1'] == pytest.approx(1.0)
    assert agg.loc[25, 'g2'] == pytest.approx(1.0)
    assert agg.loc[500, 'g2'] == pytest.approx(2.5)
    # the sample not in the metadata is ignored
    assert agg.loc[500, 'g1'] == 0


def test_aggregate_signed_split(hits, md):
    agg = metag_heatmaps.aggregate(hits, md, 'depth',
                                   genomes=['g1', 'g2', 'g3'],
                                   split_col='filter', split_negative=C)
    assert isinstance(agg.columns, pandas.MultiIndex)
    # viral first then cellular within each genome
    assert list(agg.columns) == [('g1', V), ('g1', C), ('g2', V), ('g2', C),
                                 ('g3', V), ('g3', C)]
    # mean over both viral 25m samples, g1 is absent in one
    assert agg.loc[25, ('g1', V)] == pytest.approx(2.0)
    assert agg.loc[25, ('g2', V)] == pytest.approx(2.0)
    # cellular values are negative
    assert agg.loc[500, ('g2', C)] == pytest.approx(-5.0)
    assert agg.loc[25, ('g3', C)] == pytest.approx(-2.0)
    assert (agg.xs(C, axis=1, level=1) <= 0).all().all()
    assert (agg.xs(V, axis=1, level=1) >= 0).all().all()


def test_aggregate_order(hits, md):
    agg = metag_heatmaps.aggregate(hits, md, 'month', genomes=['g2'],
                                   order=[2, 1, 3])
    assert list(agg.index) == [2, 1, 3]
    # no samples in month 3
    assert numpy.isnan(agg.loc[3, 'g2'])


def test_aggregate_pandas_input(hits, md):
    """ pandas tables work too """
    hits_pd = pandas.DataFrame(hits.to_dict(as_series=False))
    md_pd = pandas.DataFrame(md.to_dict(as_series=False))
    a = metag_heatmaps.aggregate(hits_pd, md_pd, 'depth', genomes=['g1', 'g2'])
    b = metag_heatmaps.aggregate(hits, md, 'depth', genomes=['g1', 'g2'])
    pandas.testing.assert_frame_equal(a, b)


def test_paired_frame():
    pos = pandas.DataFrame({'a': [1.0, 2.0]}, index=[1, 2])
    neg = pandas.DataFrame({'a': [3.0], 'b': [4.0]}, index=[2])
    paired = metag_heatmaps.paired_frame(pos, neg)
    assert list(paired.columns) == [('a', 'vf'), ('a', 'cf'),
                                    ('b', 'vf'), ('b', 'cf')]
    assert paired.loc[2, ('a', 'cf')] == -3.0
    assert paired.loc[1, ('b', 'cf')] == 0


def test_top_genomes_union_keeps_order(hits, md):
    # ('other' sample with g1=100 is excluded by md)
    top = metag_heatmaps.top_genomes(hits, {V: 1, C: 2}, group_col='filter',
                                     samples_md=md)
    # top viral: g1 and g2 both sum to 4, ties broken by name -> g1
    # top 2 cellular: g3 (10), g2 (10) -> tie -> g2, g3
    assert top == ['g1', 'g2', 'g3']

    top = metag_heatmaps.top_genomes(hits, {C: 1, V: 2}, group_col='filter',
                                     samples_md=md)
    assert top == ['g2', 'g1']

    # without metadata, 'other' sample counts
    assert metag_heatmaps.top_genomes(hits, 1) == ['g1']
    # with an int n per group
    assert set(metag_heatmaps.top_genomes(hits, 1, group_col='filter',
                                          samples_md=md)) == {'g1', 'g2'}


def test_top_genomes_group_col_in_hits():
    hits = polars.DataFrame({
        'sample': ['a', 'b', 'c', 'd'],
        'genome': ['x', 'y', 'z', 'x'],
        'cov_per_bn_reads': [1.0, 2.0, 5.0, 1.0],
        'dataset': ['ALOHA', 'ALOHA', 'TARA', 'TARA']})
    assert metag_heatmaps.top_genomes(hits, {'ALOHA': 1, 'TARA': 2},
                                      group_col='dataset') == ['y', 'z', 'x']


def test_signed_log():
    a = pandas.DataFrame({'x': [0.0, 1.0, -4.0], 'y': [2.0, 0.0, -0.5]})
    b = pandas.DataFrame({'z': [8.0, 0.0]})
    (la, lb), zero_val = metag_heatmaps.signed_log([a, b])
    assert zero_val == pytest.approx(0.25)
    # signs kept and zero stays zero
    assert (numpy.sign(la) == numpy.sign(a)).all().all()
    assert (numpy.sign(lb) == numpy.sign(b)).all().all()
    assert la.loc[0, 'x'] == 0
    # monotonic in |x|
    assert lb.loc[0, 'z'] > abs(la.loc[2, 'x']) > la.loc[1, 'x']
    # invertible
    back = metag_heatmaps.signed_exp(la.values, zero_val)
    numpy.testing.assert_allclose(back, a.values)

    # single frame in, single frame out
    single, zv = metag_heatmaps.signed_log(a, zero_frac=.1)
    assert isinstance(single, pandas.DataFrame)
    assert zv == pytest.approx(.05)


def test_parse_tara_samples():
    md = metag_heatmaps.parse_tara_samples(
        ['Station85_MES_0d2-3', 'Station7_SUR_0d22-1d6',
         'not_tara', 'Station100_DCM_5-20'])
    assert md.columns == ['sample', 'station', 'depth', 'type']
    assert md['station'].to_list() == [85, 7, 100]
    assert md['depth'].to_list() == ['MES', 'SUR', 'DCM']
    assert md['type'].to_list() == ['0d2-3', '0d22-1d6', '5-20']


def test_color_fn_from_dict():
    hosts = {'g1': 'Pelagibacter', 'g2': 'Oceanicaulis', 'g3': 'Weird'}
    fn = metag_heatmaps.color_fn_from_dict(
        {'Pelagibacter': 'pink', 'Alpha': 'purple'}, key_fn=hosts.get,
        translate={'Oceanicaulis': 'Alpha'}, missing='white', other='grey')
    assert [fn(g) for g in ['g1', 'g2', 'g3', 'g4']] == \
        ['pink', 'purple', 'grey', 'white']


def _synthetic_panels(n_genomes=8, seed=0):
    rng = numpy.random.default_rng(seed)
    genomes = [f"g{i}" for i in range(n_genomes)]
    md = polars.DataFrame(
        [(f"{f}_{d}_{m}", d, m, f) for f in (V, C)
         for d in (25, 75, 150, 500) for m in range(1, 13)],
        schema=['sample', 'depth', 'month', 'filter'], orient='row')
    rows = [(s, g, float(rng.lognormal(0, 1)))
            for s in md['sample'] for g in genomes if rng.uniform() < .6]
    hits = polars.DataFrame(rows, schema=['sample', 'genome',
                                          'cov_per_bn_reads'], orient='row')
    tara_md = metag_heatmaps.parse_tara_samples(
        [f"Station{s}_{d}_x"
         for s in range(5) for d in metag_heatmaps.TARA_DEPTHS])
    trows = [(s, g, float(rng.lognormal(0, 1)))
             for s in tara_md['sample'] for g in genomes[2:]]
    tara_hits = polars.DataFrame(trows, schema=['sample', 'genome',
                                                'cov_per_bn_reads'],
                                 orient='row')
    split = dict(split_col='filter', split_negative=C, split_order=[V, C])
    by_depth = metag_heatmaps.aggregate(hits, md, 'depth', genomes=genomes,
                                        **split)
    by_month = metag_heatmaps.aggregate(hits, md, 'month', genomes=genomes,
                                        order=range(1, 13), **split)
    by_tara = metag_heatmaps.aggregate(
        tara_hits, tara_md, 'depth', genomes=genomes,
        order=list(metag_heatmaps.TARA_DEPTHS))
    return genomes, by_depth, by_month, by_tara


def _ticklabels(ax):
    return [t.get_text() for t in ax.get_xticklabels()]


def test_abundance_heatmap(close_figures, tmp_path):
    genomes, by_depth, by_month, by_tara = _synthetic_panels()
    panels, zero_val = metag_heatmaps.signed_log(
        [by_depth, by_month, by_tara])
    hm = metag_heatmaps.abundance_heatmap(
        [metag_heatmaps.Panel(panels[0], 'Depth (m)\nALOHA'),
         metag_heatmaps.Panel(panels[1], 'Month\nALOHA',
                              row_labels=metag_heatmaps.MONTHS),
         metag_heatmaps.Panel(
             panels[2], 'Depth\nTARA',
             row_labels=list(metag_heatmaps.TARA_DEPTHS.values()))],
        strips=[metag_heatmaps.MetadataStrip('Type', lambda g: 'blue'),
                metag_heatmaps.MetadataStrip('Host', {'g1': 'pink'})],
        label_color_fn=lambda g: 'red' if g == 'g3' else 'black',
        label_rename_fn=str.upper,
        log_zero_val=zero_val,
    )
    assert sorted(hm.genome_order) == sorted(genomes)
    assert hm.linkage is not None

    # the column order matches the dendrogram for every panel
    from scipy.cluster.hierarchy import dendrogram
    leaves = dendrogram(hm.linkage, no_plot=True)['leaves']
    assert hm.genome_order == [genomes[i] for i in leaves]
    for panel in hm.panels:
        ordered = panel.ordered(hm.genome_order)
        top_level = list(dict.fromkeys(
            ordered.columns.get_level_values(0))) if panel.paired \
            else list(ordered.columns)
        assert top_level == hm.genome_order

    # each panel's drawn values match the data in that order
    for i, panel in enumerate(hm.panels):
        mesh = hm.axes[f'panel_{i}'].collections[0]
        drawn = numpy.ma.filled(mesh.get_array(), numpy.nan).reshape(
            panel.data.shape)
        numpy.testing.assert_allclose(
            drawn, panel.ordered(hm.genome_order).values)

    # paired columns: viral (positive) then cellular (negative)
    depth = hm.panels[0].ordered(hm.genome_order)
    assert (depth.values[:, ::2] >= 0).all()
    assert (depth.values[:, 1::2] <= 0).all()

    # all panels and strips share an x axis with one unit per genome
    for name in ['panel_0', 'panel_1', 'panel_2', 'strips']:
        assert hm.axes[name].get_xlim() == (0, len(genomes))

    # genome labels under the strips, renamed and colored
    labels = _ticklabels(hm.axes['strips'])
    assert labels == [g.upper() for g in hm.genome_order]
    colors = {t.get_text(): t.get_color()
              for t in hm.axes['strips'].get_xticklabels()}
    assert colors['G3'] == 'red'

    assert [t.get_text() for t in hm.axes['panel_1'].get_yticklabels()] \
        == metag_heatmaps.MONTHS
    assert 'colorbar' in hm.axes
    assert [t.get_text() for t in hm.axes['colorbar'].get_xticklabels()] \
        == list(metag_heatmaps.FRACTION_LABELS)
    # first color bar tick is 0
    assert hm.axes['colorbar'].get_yticklabels()[0].get_text() == '0'

    # panel letters
    texts = [t.get_text() for t in hm.fig.texts]
    assert texts == ['(A)', '(B)', '(C)', '(D)']

    pdf = tmp_path / 'fig.pdf'
    metag_heatmaps.save_pdf(hm, str(pdf))
    assert os.path.getsize(pdf) > 0


def test_abundance_heatmap_fixed_order_single_panel(close_figures):
    genomes, by_depth, by_month, by_tara = _synthetic_panels()
    order = genomes[::-1]
    hm = metag_heatmaps.abundance_heatmap([by_tara], genome_order=order,
                                          panel_letters=False, fontsize=8)
    assert hm.genome_order == order
    assert hm.linkage is None
    assert 'dendrogram' not in hm.axes
    assert _ticklabels(hm.axes['panel_0']) == order
    # unpaired -> standard color bar
    assert hm.axes['colorbar'].get_xticklabels() == [] or \
        not any(t.get_text() for t in hm.axes['colorbar'].get_xticklabels())
    assert hm.fig.texts == []


def test_order_genomes_handles_constant_columns():
    data = pandas.DataFrame({'a': [1.0, 2.0, 3.0], 'b': [0.0, 0.0, 0.0],
                             'c': [3.0, 2.0, 1.0], 'd': [1.1, 2.1, 2.9]})
    order, link = metag_heatmaps.order_genomes([metag_heatmaps.Panel(data)])
    assert sorted(order) == ['a', 'b', 'c', 'd']
    # a and d are perfectly correlated, so they should be neighbors
    assert abs(order.index('a') - order.index('d')) == 1


def test_double_color_bar_linear(ax):
    metag_heatmaps.double_color_bar(ax, 'PuOr', 10,
                                    fraction_labels=('pos', 'neg'))
    ticks = [float(t.get_text()) for t in ax.get_yticklabels()]
    assert ticks[0] == 0 and ticks[-1] <= 10
    assert _ticklabels(ax) == ['pos', 'neg']


def _colorbar_tick_values(hm):
    cax = hm.axes['colorbar']
    return [float(t.get_text().replace('−', '-'))
            for t in cax.get_yticklabels() if t.get_text()]


def test_abundance_heatmap_one_sided(close_figures):
    """ unpaired, all positive data: only the positive half of the scale """
    genomes, by_depth, by_month, by_tara = _synthetic_panels()
    (tara_log,), zero_val = metag_heatmaps.signed_log([by_tara])
    hm = metag_heatmaps.abundance_heatmap([tara_log], log_zero_val=zero_val,
                                          panel_letters=False)
    ticks = _colorbar_tick_values(hm)
    assert min(ticks) == 0 and max(ticks) > 0
    mesh = hm.axes['panel_0'].collections[0]
    assert mesh.norm.vmin == 0 and mesh.norm.vmax == hm.vmax
    # zero is drawn white-ish (the middle of PuOr), the max is dark purple
    assert numpy.allclose(mesh.cmap(0.0), plt_cmap('PuOr')(.5), atol=.02)
    assert numpy.allclose(mesh.cmap(1.0), plt_cmap('PuOr')(1.0), atol=.02)

    # negated data uses the other half
    hm = metag_heatmaps.abundance_heatmap([-tara_log], log_zero_val=zero_val,
                                          panel_letters=False)
    ticks = _colorbar_tick_values(hm)
    assert max(ticks) == 0 and min(ticks) < 0

    # sign='both' keeps the full diverging scale
    hm = metag_heatmaps.abundance_heatmap([tara_log], log_zero_val=zero_val,
                                          panel_letters=False, sign='both')
    ticks = _colorbar_tick_values(hm)
    assert min(ticks) < 0 < max(ticks)


def plt_cmap(name):
    from matplotlib import pyplot as plt
    return plt.get_cmap(name)


def test_log_ticks_not_crowded_at_zero():
    # tiny zero_val (like real coverages) used to give a 1e-4 tick on top of 0
    zero_val = 5e-5
    vmax = float(numpy.log(1e4 + zero_val) - numpy.log(zero_val))
    ticks = metag_heatmaps._log_ticks(vmax, zero_val)
    assert ticks[0] == 0
    positions = numpy.log(numpy.array(ticks[1:]) + zero_val) - \
        numpy.log(zero_val)
    assert (positions >= .08 * vmax).all()
