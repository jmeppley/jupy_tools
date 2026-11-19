"""
Tests for the plotting modules moved out of experimental: plot_metadata,
gene_maps, tree_plots and circular_trees. All data is synthetic and figures
are drawn with the Agg backend (see conftest.py).
"""
import json
import warnings
from io import StringIO

import numpy
import pandas
import pytest

from Bio import Phylo
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba

from jme.jupy_tools import circular_trees, gene_maps, plot_metadata, \
    tree_plots, utils

NEWICK = ("((A:0.1,B:0.2)90:0.3,((C:0.1,D:0.1)75:0.2,"
          "(E:0.3,F:0.1)40:0.1)99:0.2);")
GENOMES = list("ABCDEF")
ANNOTS = ['VOG1', 'VOG2', 'VOG3', 'VOG4', 'VOG5']
DESCS = {'VOG1': 'terminase large subunit', 'VOG2': 'portal protein',
         'VOG3': 'major capsid protein', 'VOG4': 'hypothetical protein',
         'VOG5': 'tail fiber'}


@pytest.fixture
def tree():
    return Phylo.read(StringIO(NEWICK), 'newick')


@pytest.fixture
def md_df():
    """ genome metadata: pandas, indexed by genome """
    return pandas.DataFrame({
        'genome': GENOMES,
        'Length': [10000, 12000, 30000, 31000, 29000, 8000],
        'GC': [.31, .33, .45, .44, .50, .38],
        'seq_type': ['Ref', 'Ref', 'Contig', 'Contig', 'Contig', 'Ref'],
        'Host': ['Cyano', 'Cyano', None, 'SAR11', None, 'SAR11'],
        'is_novel': [False, False, True, True, True, False],
    }).set_index('genome')


@pytest.fixture
def gene_df(md_df):
    """ 5 genes per genome, one per annotation, D on the - strand """
    rows = []
    for genome, length in md_df.Length.items():
        step = length // 6
        strand = -1 if genome == 'D' else 1
        for i, annot in enumerate(ANNOTS):
            rows.append((f"{genome}_{i + 1}", genome, i + 1, i * step + 10,
                         i * step + step - 10, strand, annot, DESCS[annot]))
    return pandas.DataFrame(
        rows, columns=['gene', 'genome', 'gene_no', 'start', 'end', 'strand',
                       'annot', 'desc']).set_index('gene')


@pytest.fixture
def config_dict():
    """ a plot_metadata dict as it comes out of json.load """
    return {
        "genomes": {
            "metadata": [
                {"column": "seq_type", "label": "Type",
                 "color_dict": {"Ref": "black", "Contig": "blue",
                                "null": [1.0, 1.0, 1.0, 0.0]}},
                {"column": "is_novel", "label": "New",
                 "color_dict": {"true": "red", "false": "white"}},
                None,
                {"column": "Length", "log": True,
                 "color_map": "nipy_spectral_r"},
                {"column": "GC", "color_map": "brg"},
                {"column": "Host", "use_tree_seqs": True,
                 "color_dict": {"Cyano": "cyan", "SAR11": "red"}},
            ],
            "table_locs": {"length_col": "Length", "hostname_col": "Host",
                           "seq_type_col": "seq_type", "ref_seq_type": "Ref",
                           "contig_seq_type": "Contig", "clusters_col": None},
        },
        "genes": {"n_genes": 4, "base_colors": ["red", "blue", "green"],
                  "type_patterns": [["terminase", "terminase"],
                                    ["portal", "portal"]]},
    }


def to_polars(df):
    pytest.importorskip("polars")
    return utils.to_polars(df.reset_index())


##########
# plot_metadata

def test_load_plot_metadata_json_keys(tmp_path, config_dict):
    path = tmp_path / 'plot_metadata.json'
    path.write_text(json.dumps(config_dict))
    config = plot_metadata.load_plot_metadata(str(path))
    assert len(config.genome_metadata) == 6
    assert config.genome_metadata[2] is None
    new = config.get_spec('New')
    # JSON can't have bool keys: they're aliased
    assert new.color_dict[True] == 'red'
    assert new.color_dict[False] == 'white'
    assert config.get_spec('Type').color_dict[None] == [1.0, 1.0, 1.0, 0.0]
    assert config.genes['n_genes'] == 4
    # passing a config or dict through gives the same thing
    assert plot_metadata.load_plot_metadata(config) is config
    assert plot_metadata.load_plot_metadata(None).genome_metadata == []


def test_plot_config_round_trip(tmp_path, config_dict):
    config = plot_metadata.PlotConfig.from_dict(config_dict)
    path = tmp_path / 'out.json'
    config.to_json(str(path))
    reloaded = json.loads(path.read_text())
    # bool aliases collapse back to one string key
    assert reloaded['genomes']['metadata'][1]['color_dict'] == \
        {"true": "red", "false": "white"}
    assert reloaded['genomes']['table_locs'] == \
        config_dict['genomes']['table_locs']
    assert reloaded['genes'] == config_dict['genes']
    assert plot_metadata.PlotConfig.from_dict(reloaded).to_dict() \
        == config.to_dict()


def test_spec_build_categorical(md_df):
    spec = plot_metadata.MetadataSpec(column='seq_type', label='Type',
                                      color_dict={'Ref': 'k', 'Contig': 'b'})
    md = spec.build(md_df)
    assert isinstance(md, plot_metadata.TreeMetadata)
    assert md.label == 'Type'
    assert md.get_item_color('A') == 'k'
    assert md.get_item_color('C') == 'b'
    assert dict(md.get_category_colors(['A', 'B'])) == {'Ref': 'k'}


def test_spec_callable_colors_and_polars(md_df):
    spec = plot_metadata.MetadataSpec.from_dict(
        {'column': 'Host',
         'color_dict': lambda host: 'cyan' if host == 'Cyano' else 'grey'})
    fn = spec.color_fn(to_polars(md_df))
    assert fn('A') == 'cyan'
    assert fn('D') == 'grey'
    # missing values get the (transparent) null color
    assert to_rgba(fn('C'))[3] == 0


def test_spec_continuous(md_df):
    spec = plot_metadata.MetadataSpec(column='Length', log=True)
    md = spec.build(md_df)
    assert md.label == 'Length'
    c_small, c_big = to_rgba(md.get_item_color('F')), \
        to_rgba(md.get_item_color('D'))
    assert c_small != c_big
    # bar plot when color_map is a single color
    bars = plot_metadata.MetadataSpec(column='GC', color_map='red').build(
        md_df)
    assert bars.size_norm is not None
    assert bars.get_item_size('E') == pytest.approx(1.0)


def test_tree_metadata_collapsed_clade(tree, md_df):
    md = plot_metadata.MetadataSpec(column='seq_type').build(md_df)
    clade_ab = tree.common_ancestor('A', 'B')
    # single category: that color, fully opaque
    assert md.get_item_color(clade_ab) == to_rgba(md.get_color('Ref'))
    size_md = plot_metadata.MetadataSpec(column='GC',
                                         color_map='red').build(md_df)
    # (this used an undefined variable before the move)
    assert 0 <= size_md.get_item_size(clade_ab) <= 1


def test_color_fn_from_dict():
    fn = plot_metadata.color_fn_from_dict({'x': 'red'}, key_fn=str.lower,
                                          other='grey')
    assert fn('X') == 'red'
    assert fn('Y') == 'grey'


def test_get_min_max():
    assert plot_metadata.get_min_max([3, 1, 2]) == (1, 3)


##########
# gene_maps

def test_import_genes(tmp_path):
    faa = tmp_path / 'genes.faa'
    faa.write_text(
        ">contigA_1 # 1 # 300 # 1 # ID=1_1;partial=00\nMKV\n"
        ">contigA_2 # 400 # 900 # -1 # ID=1_2;partial=00\nMKV\n"
        ">contigB_1 # 5 # 200 # 1 # ID=2_1;partial=00\nMKV\n")
    genes = gene_maps.import_genes(str(faa), {'contigA_2': 'VOG9'})
    assert list(genes.index) == ['contigA_1', 'contigA_2', 'contigB_1']
    assert list(genes.genome) == ['contigA', 'contigA', 'contigB']
    assert genes.loc['contigA_2', 'strand'] == -1
    assert genes.loc['contigA_2', 'annot'] == 'VOG9'


def test_get_cluster_genes_flips_strand(gene_df, md_df):
    genes = gene_maps.get_cluster_genes(gene_df, ['C', 'D'],
                                        md_df.Length.to_dict())
    assert set(genes.genome) == {'C', 'D'}
    assert (genes.strand == 1).all()
    # and accepts polars
    genes_pl = gene_maps.get_cluster_genes(to_polars(gene_df), ['C', 'D'],
                                           md_df.Length.to_dict())
    pandas.testing.assert_frame_equal(genes, genes_pl)


def test_draw_genes_and_synteny_plot(gene_df, md_df, close_figures):
    genes = gene_maps.get_cluster_genes(gene_df, GENOMES,
                                        md_df.Length.to_dict())
    colors = {'VOG1': 'red', 'VOG2': 'blue'}
    fig, ax = plt.subplots()
    gene_maps.draw_genes(genes, colors, ax=ax)
    assert len(ax.patches) + len(ax.collections) > 0

    fig, (ax1, ax2) = plt.subplots(2, 1)
    out = gene_maps.cluster_synteny_plot(
        {g: i for i, g in enumerate(GENOMES)}, ax1, ax2,
        md_df.Length.to_dict(), gene_df, max_colored_genes=3)
    assert out is not None


def test_get_subcluster_colors():
    colors = gene_maps.get_subcluster_colors({'a': 1, 'b': 1, 'c': 2})
    assert colors['a'] == colors['b'] != colors['c']


##########
# tree_plots

def test_get_tree_seqs_and_read_tree(tree, tmp_path):
    assert tree_plots.get_tree_seqs(tree) == set(GENOMES)
    path = tmp_path / 't.nwk'
    path.write_text(NEWICK)
    assert tree_plots.get_tree_seqs(tree_plots.read_tree(str(path))) == \
        set(GENOMES)


def test_synteny_tree_plotter_init(md_df, gene_df, config_dict):
    plotter = tree_plots.SyntenyTreePlotter(md_df, gene_df, config_dict,
                                            annot_col='annot')
    # caller's dict isn't modified, config is parsed
    assert 'leaf_color' not in config_dict['genomes']
    assert isinstance(plotter.plot_config, plot_metadata.PlotConfig)
    assert plotter.type_patterns == config_dict['genes']['type_patterns']
    assert plotter.get_type('terminase large subunit') == 'terminase'
    # the old name still works
    assert tree_plots.TreePlotterVLP is tree_plots.SyntenyTreePlotter


@pytest.mark.parametrize('gene_method', [tree_plots.GM_TOP,
                                         tree_plots.GM_TYPES])
def test_draw_tree(tree, md_df, gene_df, config_dict, close_figures,
                   gene_method):
    plotter = tree_plots.SyntenyTreePlotter(md_df, gene_df, config_dict,
                                            annot_col='annot')
    fig = plotter.draw_tree('test tree', tree, collapse=False,
                            gene_method=gene_method)
    assert fig is not None
    assert len(fig.axes) >= 3


def test_draw_tree_polars_and_collapse(tree, md_df, gene_df, config_dict,
                                       close_figures):
    plotter = tree_plots.SyntenyTreePlotter(to_polars(md_df),
                                            to_polars(gene_df), config_dict,
                                            annot_col='annot')

    def check_clade(clade):
        seqs = tree_plots.get_tree_seqs(clade)
        if seqs == {'A', 'B'}:
            return "2 refs", "Ref"
        return None

    fig = plotter.draw_tree('collapsed', tree,
                            collapse=dict(check_clade_fn=check_clade),
                            length_hist=True)
    assert fig is not None


def test_draw_tree_no_genes(tree, md_df, gene_df, config_dict,
                            close_figures):
    plotter = tree_plots.SyntenyTreePlotter(md_df, gene_df, config_dict,
                                            annot_col='annot')
    fig, md_list = plotter.draw_tree('no genes', tree, collapse=False,
                                     draw_genes=False, draw_legends=False)
    assert [md.label if md else None for md in md_list][:2] == \
        ['Type', 'New']


def test_draw_plain_phylogram(tree, ax):
    tree_plots.draw(tree, axes=ax, do_show=False)
    assert len(ax.collections) + len(ax.lines) > 0


##########
# circular_trees

def test_ring_order_and_clade_styles(tree):
    order = circular_trees.get_ring_order(tree.root, {})
    assert [c.name for c in order] == GENOMES

    clade_ab = tree.common_ancestor('A', 'B')
    collapsed = {clade_ab: 3}
    order = circular_trees.get_ring_order(tree.root, collapsed)
    assert order[:3] == [clade_ab] * 3

    types = {'A': 'x', 'B': 'x', 'C': 'y', 'D': 'y', 'E': 'y', 'F': 'z'}
    styles = {}
    circular_trees.get_clade_style_by_type(
        tree.root, styles, types.get,
        {'x': {'fc': 'red'}, 'y': {'fc': 'blue'}, 'z': {'fc': 'green'}})
    assert styles[clade_ab] == {'fc': 'red'}
    assert styles[tree.common_ancestor('C', 'D')] == {'fc': 'blue'}


def test_draw_circular_tree(tree, md_df, gene_df, config_dict,
                            close_figures):
    plotter = tree_plots.SyntenyTreePlotter(md_df, gene_df, config_dict,
                                            annot_col='annot')
    clade_ab = tree.common_ancestor('A', 'B')
    fig, md_list = plotter.draw_circular_tree(
        'circle', tree, collapsed_dict={clade_ab: dict(width=2, fc='grey')},
        bg_styles_dict={tree.common_ancestor('C', 'D'): dict(fc='yellow')},
        support_md=circular_trees.SupportMetadata(),
    )
    assert len(md_list) == 6
    ax = fig.axes[0]
    assert len(ax.patches) + len(ax.collections) > 0
    xmin, xmax = ax.get_xlim()
    assert xmin == -xmax


def test_circular_tree_without_leaf_arc(tree, ax):
    """ the leaf_arc=None branch referenced undefined names before """
    radius = circular_trees.draw_circular_tree_recursive(
        ax, tree.root, collapsed_dict={'nothing': 2})
    assert radius > 0


##########
# old experimental imports

def test_experimental_shims_warn_and_reexport():
    import importlib
    from jme.jupy_tools.experimental import cluster_trees, gene_synteny, \
        metag_heatmaps, vlp_plots_2
    from jme.jupy_tools import abundance_heatmaps
    for module in (cluster_trees, gene_synteny, metag_heatmaps, vlp_plots_2):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            importlib.reload(module)
        assert any(issubclass(w.category, DeprecationWarning)
                   for w in caught), module.__name__

    assert vlp_plots_2.TreePlotterVLP is tree_plots.SyntenyTreePlotter
    assert vlp_plots_2.GM_TYPES == tree_plots.GM_TYPES
    for name in ['get_tree_seqs', 'get_clade_style_by_type', 'get_min_max',
                 'add_cat_md_legend', 'get_tree_metadata', 'draw_md_ring',
                 'get_ring_order', 'SupportMetadata']:
        assert hasattr(vlp_plots_2, name), name
    assert gene_synteny.import_genes is gene_maps.import_genes
    assert gene_synteny.get_N_colors is not None
    assert cluster_trees.TreeMetadata is plot_metadata.TreeMetadata
    assert cluster_trees.TreePlotter is tree_plots.TreePlotter
    assert cluster_trees.cluster_synteny_plot is gene_maps.cluster_synteny_plot
    # tree building stays in experimental
    assert cluster_trees.TreeBuilder.__module__ == \
        'jme.jupy_tools.experimental.cluster_trees'
    assert metag_heatmaps.abundance_heatmap is \
        abundance_heatmaps.abundance_heatmap
    assert metag_heatmaps._log_ticks is abundance_heatmaps._log_ticks


def test_gene_type_colors(md_df, gene_df, config_dict):
    """ GM_TYPES colours genes by the config's type_patterns """
    plotter = tree_plots.SyntenyTreePlotter(md_df, gene_df, config_dict,
                                            annot_col='annot')
    colors, descs = plotter.get_gene_colors_and_descs(
        set(GENOMES), tree_plots.GM_TYPES)
    assert colors['VOG1'] == 'red'    # terminase -> first base color
    assert colors['VOG2'] == 'blue'   # portal -> second
    assert len(colors) == config_dict['genes']['n_genes']
    assert 'terminase' in descs['VOG1']
