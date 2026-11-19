"""
Shared genome metadata and colour model for the tree, gene map and heatmap
plots (tree_plots, circular_trees, gene_maps, abundance_heatmaps).

Data
====

Genome metadata table
    One row per genome/sequence. A pandas table is indexed by genome name; a
    polars table needs a ``genome`` column (see to_pandas_indexed). Column
    names are free form. The ones the tree plots need (length, seq type,
    host, cluster) are named in PlotConfig.table_locs.

Metadata spec (MetadataSpec)
    How to turn one metadata column into colours. Its fields are the keys
    already used in the plot_metadata JSON files::

        column, label, md_type ('cat' or 'cont'), color_dict, color_map,
        log, null_value, null_color, rename, use_tree_seqs, bar_plot

    color_dict may be a dict (value -> color) or a callable (value -> color).
    color_map is a matplotlib colormap name/object (or, for continuous data,
    a single color to draw bars). A ``None`` in a list of specs is a spacer.

    spec.build(md_table) -> TreeMetadata (what the tree plots draw)
    spec.color_fn(md_table) -> callable genome -> color (eg for heatmap strips)

Plot config (PlotConfig)
    The whole plot_metadata JSON::

        {"genomes": {"metadata": [spec, None, spec, ...],
                     "table_locs": {"length_col": ..., ...},
                     "leaf_color": spec},
         "genes": {"n_genes": 20, "base_colors": [...],
                   "type_patterns": [[name, regex], ...]},
         "draw_args": {...}}

    load_plot_metadata(path_or_dict) reads it. Because JSON can't have
    boolean or null keys, "true"/"false"/"null" keys in color dicts are
    also made available as True/False/None.

Moved from experimental.cluster_trees (TreeMetadata), experimental.vlp_plots_2
(get_tree_metadata and legend helpers) and experimental.metag_heatmaps
(color_fn_from_dict).
"""
import json
from collections import Counter
from dataclasses import dataclass, field, fields
from typing import Any, Dict, List, Optional

import attr
import numpy
import pandas
from Bio import Phylo
from pandas.api.types import is_numeric_dtype, is_bool_dtype

import matplotlib.collections as mcollections
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba

from jme.jupy_tools.utils import get_N_colors, to_pandas_indexed


@attr.s
class TreeMetadata():
    """
    init params:
        label: :str: (required)
        data_dict: :dict: (required)
            a dictionary mapping names to metadata values
        null_value: (default: False)
            value to always map to white
        color_dict: :dict:
            direct map from values to colors (overrides other color options)
        color_map: str or object
            matplotlib colormap name or object (default: "binary")
        color_method: str or callable
            either "category" or "cat" (default) or "fixed" 
            (chooses evenly spaced colors from color_map for each unique value in data_dict)
            or "continuous" or "con"
            (data_dict values msut be in range [0.-1.] or [-255])
            or a callable transformation to convert values to a cmap-suitable range (implyies continuous)
        size_norm: a function normalizing values to the 0. to 1. interval
            Reterun a scaling value for bar plots metadata.
            Defaults to None: (all the same size)
        null_color: :str: or :tuple:
            color for missing values (default: transparent)            
        rename: :dict: or callable or None (default)
            for 'category' metadata, translates category names in
            get_category_colors() (if not None)
    """

    label = attr.ib()
    data_dict = attr.ib()
    null_value = attr.ib(default=False)
    color_dict = attr.ib(default=None)
    color_map = attr.ib(default="binary")
    color_method = attr.ib(default="category")
    size_norm = attr.ib(default=None)
    null_color = attr.ib(default=(1., 1., 1., 0.))
    relative_size = attr.ib(default=1)
    collapsed_clade_labels = attr.ib(default=None)
    rename = attr.ib(default=None)

    def __attrs_post_init__(self):
        """ using the color attributes, set up the get_color method """
        if self.color_dict:
            self._use_color_dict()
            # make a copy
            self.color_dict = dict(self.color_dict)
        else:
            if isinstance(self.color_method, str):
                c_method = self.color_method.lower()[:3]
            elif callable(self.color_method):
                c_method = "con"
            if c_method == 'cat':
                # build a categorical color dict
                categories = set(self.data_dict.values())
                self.color_dict = dict(zip(categories, 
                                           get_N_colors(len(categories),
                                                                     self.color_map)))
                self._use_color_dict()
            elif c_method == 'con':
                # was a transform provided?
                if callable(self.color_method):
                    self._transform = self.color_method
                else:
                    self._transform = lambda x: x
                
                # make sure colormap is callable:
                if not callable(self.color_map):
                    self.color_map = plt.get_cmap(self.color_map)
                    
                self.get_color = self._get_color_from_map
            elif c_method == 'fix':
                # a fixed color, usually used with sized bars

                # either get the "1.0" value from the color map or interpret
                # the color_map argumetn as a color
                # make sure colormap is callable:
                color_thing, thing_type = TreeMetadata._check_color_map(self.color_map)
                if thing_type == 'MAP':
                    fixed_color = color_thing(1.0)
                else:
                    fixed_color = color_thing

                self.get_color = lambda x: (
                    self.null_color 
                    if x == self.null_value
                    else fixed_color
                )

    def _check_color_map(color_map):
        try:
            color = to_rgba(color_map)
            return color, 'COLOR'
        except ValueError:
            # it was not a color, so it's a color map
            if not callable(color_map):
                try:
                    color_map = plt.get_cmap(color_map)
                except ValueError:
                    raise ValueError("For a fixed color Metadata, pass"
                                     " either a color or a color_map.")
            return color_map, 'MAP'

    def _use_color_dict(self):
        self.get_color = self._get_color_from_dict
        # add null value to color dict
        self.color_dict[self.null_value] = self.null_color

    def get_category_colors(self, seqs=None):
        """ yield uples of value, color pairs for all the possible value/color
        pairs represented by this metadata object.
        
        If a collection of seqs is given, only use values/colors returnable by
        itmes in that collection
        """
        # If a rename function or dict given, rename the values returned
        if callable(self.rename):
            rename_fn = self.rename
        elif self.rename is None:
            def rename_fn(value):
                return value
        else:
            def rename_fn(value):
                return self.rename.get(value, value)
        
        # get values (either from seqs or all values)
        if seqs:
            unique_values = set(self.data_dict.get(s, self.null_value)
                                for s in seqs)
        else:
            if self.color_dict is not None:
                unique_values = self.color_dict.keys()
            else:
                unique_values = set(self.data_dict.keys())

        # generate tuples of values paired with the corrsponding color
        for value in unique_values:
            yield (rename_fn(value), self.color_dict[value])

    def get_color(self):
        pass
    
    def _get_color_from_dict(self, value):
        return self.color_dict[value]

    def _get_color_from_map(self, value):
        if value == self.null_value:
            return self.null_color
        return self.color_map(self._transform(value))

    def get_item_value(self, item):
        if isinstance(item, Phylo.Newick.Clade):
            if item.is_terminal():
                return self.get_item_value(item.name)
            else:
                # handle the case where we get a collapsed clade
                if self.color_dict:
                    # catergorical md, retun all the values unless:
                    # is there a special entry for this clade?
                    if item in self.data_dict:
                        return self.data_dict[item]

                    # return all values for this clade
                    return Counter(self.data_dict.get(t.name, self.null_value) 
                                           for t in item.get_terminals())
                else:
                    # for continuous color maps, average the values
                    try:
                        mean_val = numpy.mean(
                            [self.data_dict.get(t.name, self.null_value)
                             for t in item.get_terminals()]
                        )
                    except TypeError:
                        mean_val = self.null_value
                    return mean_val
                    
        # if it's not a clade object, assume it's in the data dict if it has a value
        return self.data_dict.get(item, self.null_value)

    def get_item_color(self, item):
        item_value = self.get_item_value(item)
        if (
            isinstance(item_value, Counter)
            and
            isinstance(item, Phylo.Newick.Clade)
       ):
            # this clade was not in the data dict, check to see if it only has
            # one value
            categories = item_value
            #categories = Counter(self.data_dict.get(t.name, self.null_value) 
            #                       for t in item.get_terminals())
            non_null = []
            total, total_non_null = 0, 0
            for category, count in categories.items():
                total += count
                if category != self.null_value:
                    non_null.append(category)
                    total_non_null += count
            if len(non_null) == 1:
                # get the RGB of the mapped color
                color = to_rgba(self.get_color(non_null[0]))[:3]
                # make transparent in proportion to fraction with non-null value
                alpha = total_non_null / total                        
                return tuple(list(color) + [alpha,])
            
            # too many categories or all were null, use null color
            return self.null_color

        else:
            # we should have a single value to apply to the color dict or color
            # function
            return self.get_color(item_value)

    def get_item_size(self, item):
        item_value = self.get_item_value(item)
        if (
            isinstance(item_value, Counter)
            and
            isinstance(item, Phylo.Newick.Clade)
       ):
            # this clade was not in the data dict, do a weighted mean
            total, n_non_null = 0, 0
            for category, count in item_value.items():
                if category == self.null_value or pandas.isna(category):
                    continue
                total += count * category
                n_non_null += count
            item_value = total / n_non_null
        return self.size_norm(item_value)


def add_cat_md_legend(label, values_colors, fig, **legend_kwargs):
    patches = []
    names = []
    for value, color in values_colors:
        patches.append(mpatches.Patch(color=color, ec='k', lw=.25))
        names.append(str(value))

    for k,v in dict(bbox_to_anchor=[.74, .81],
                    loc='upper left',
                    facecolor='white',
                    title=label).items():
        legend_kwargs.setdefault(k, v)
    return fig.legend(patches, names, **legend_kwargs)

def get_tree_metadata(
    column,        # A named pandas.Series from the md_df DataFrame
    seqs=None,
    md_type=None,
    **kwargs,
):
    """
    Creates a TreeMetadata object from a metadata column.

    Can be configured with options:

        md_type: one of categorical(or 'cat') or continuous (or 'cont')
                 default: auto from datatype
        label: Any string. Default: column.name
        other args passed to get_continuous_md() or get_cat_md()
    """
    data_series = column.dropna()

    """ limit seqs to only those requested that have non-null values """
    if seqs is None:
        seqs = set(data_series.index)
    else:
        seqs = set(seqs).intersection(data_series.index)

    if md_type is None:
        md_type = 'cont' \
            if (is_numeric_dtype(data_series) 
                and 
                not is_bool_dtype(data_series)
               )\
            else 'cat'
    if md_type.lower().startswith('cont'):
        return get_continuous_md(data_series, seqs, **kwargs)
    else:
        return get_cat_md(data_series, seqs, **kwargs)

def get_cat_md(data_series,
               seqs,
               label=None,
               null_value=None,
               color_dict=None,
               color_map=None,
               **kwargs,
              ):
    if label is None:
        label = data_series.name
    data_dict = get_md_dict_from_series(data_series, seqs)
    if callable(color_dict):
        # a function from value to color: evaluate it on the values we have
        color_fn = color_dict
        color_dict = {v: color_fn(v) for v in set(data_dict.values())}
    if color_dict is None:
        uniq = list(set(data_dict.values()))
        N = len(uniq)
        if color_map is None:
            if N <= 10:
                color_map = 'tab10'
            elif N <= 20:
                color_map = 'tab20'
            else:
                color_map = 'nipy_spectral'
        color_dict = dict(zip(
            uniq,
            get_N_colors(N, color_map)
        ))
        
        # if there is no value, it should be transparent
        kwargs.setdefault('null_color', [1., 1., 1., 0.])
    return TreeMetadata(
        label,
        data_dict=data_dict,
        null_value=null_value,
        color_dict=color_dict,
        color_method='cat',
        **kwargs,
    )

def get_md_dict_from_series(data, seqs):
    return {k:v
            for k,v in data.items()
            if pandas.notna(v) and k in seqs
           }

def get_continuous_md(
    data_series,
    seqs,
    label=None,
    color_map='binary',
    log=False,
    log_shift=1e-12,
    bar_plot=None,
    **kwargs,
):
    """
    Generate a continuos TreeMetedata object from a data series

    Calculates a Normalization function to map all values into [0., 1.0]

    If bar_plot is True:
        * normalization is passed as the size_norm to make a bar plot
    If bar_plot is None and color_map is a color (not a map):
        * normalization is passed as the size_norm and the color method is set
        to "fixed"
    If bar_plot is False and the color_map is a color:
        * the color method is dset to 'fixed' and you'll get a useless md
        object
    """
    name = data_series.name if label is None else label
    " construct Metadata object from a column of numeric data "

    color_thing, thing_type = \
        TreeMetadata._check_color_map(color_map)
    kwargs.setdefault('color_map', color_thing)
    kwargs.setdefault('null_color', [1.,1.,1.,0.,])

    # if no specifics are supploied, do some magic
    if ('color_method' not in kwargs) and ('size_norm' not in kwargs):
        if log:
            values = {s:data_series[s] + log_shift
                      for s in seqs}
            vmin, vmax = min(values.values()), max(values.values())
            v_norm = mcolors.LogNorm(vmin=vmin, vmax=vmax, clip=True)
            null_val = log_shift
        else:
            values = {s:data_series[s] for s in seqs}
            vmin, vmax = min(values.values()), max(values.values())
            v_norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
            null_val = None

        if thing_type == "COLOR":
            color_method = 'fixed'
            if bar_plot is None:
                bar_plot = True
        else:
            color_method = v_norm

        if bar_plot is True:
            size_norm = v_norm
        else:
            size_norm = None

        for k,v in dict(
            null_value=null_val,
            color_method=color_method,
            size_norm=size_norm,
        ).items():
            kwargs.setdefault(k, v)
    else:
        values = {s:data_series[s] for s in seqs}
        kwargs.setdefault('size_norm', None)

    return TreeMetadata(
        name,
        values,
        **kwargs,
    )


def add_legend(fig, colors_by_label, size=15, **legend_kwargs):
    patches = []
    names = []
    for label, color in colors_by_label:
        patches.append(mcollections.CircleCollection([size,], facecolors=[color,]))
        names.append(f"{label}")

    for k,v in dict(loc='upper left',
                    facecolor='white',
                    title='Colors').items():
        legend_kwargs.setdefault(k, v)
    return fig.legend(patches, names,
                      **legend_kwargs)


def get_min_max(values):
    mn, mx = None, None
    for value in values:
        if mn is None and mx is None:
            mn, mx = value, value
        else:
            mn = min(mn, value)
            mx = max(mx, value)
    return mn, mx


def color_fn_from_dict(colors, key_fn=None, translate=None,
                       missing=(1., 1., 1., 1.), other=None):
    """
    Build a genome -> color function from a dict.

    colors: {category: color}
    key_fn: genome -> category (or None if unknown). Default: identity.
            Eg known_hosts.get or lambda g: g[:3]
    translate: optional {category: category} to apply first
               (eg {'Synechococcaceae': 'Cyanobacteria'})
    missing: color when key_fn returns None
    other: color for categories not in colors (default: same as missing)
    """
    translate = translate or {}
    other = missing if other is None else other

    def color_fn(genome):
        key = genome if key_fn is None else key_fn(genome)
        if key is None:
            return missing
        key = translate.get(key, key)
        return colors.get(key, other)
    return color_fn


##########
# Specs and config (the shared data model)

# JSON can't have boolean or null keys
_JSON_KEY_ALIASES = {'true': True, 'false': False, 'null': None}


def _alias_json_keys(color_dict):
    """ add True/False/None keys for "true"/"false"/"null" (keeps both) """
    if not isinstance(color_dict, dict):
        return color_dict
    color_dict = dict(color_dict)
    for str_key, key in _JSON_KEY_ALIASES.items():
        if str_key in color_dict and key not in color_dict:
            color_dict[key] = color_dict[str_key]
    return color_dict


def _jsonable_color_dict(color_dict):
    """ inverse of _alias_json_keys: one string key per entry """
    out = {}
    for key, color in color_dict.items():
        if key is True or key is False or key is None:
            key = json.dumps(key)
        else:
            key = str(key)
        if isinstance(color, tuple):
            color = list(color)
        out.setdefault(key, color)
    return out


@dataclass
class MetadataSpec:
    """
    How to draw one genome metadata column (see module docstring).

    Any keys not listed as fields are kept in ``extra`` and passed on to
    TreeMetadata (eg relative_size, color_method, size_norm).
    """
    column: Optional[str] = None
    label: Optional[str] = None
    md_type: Optional[str] = None
    color_dict: Any = None
    color_map: Any = None
    log: bool = False
    null_value: Any = None
    null_color: Any = None
    rename: Any = None
    use_tree_seqs: bool = False
    bar_plot: Any = None
    extra: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, config):
        """ None (a spacer) stays None, specs pass through """
        if config is None or isinstance(config, cls):
            return config
        names = {f.name for f in fields(cls)} - {'extra'}
        kwargs = {k: v for k, v in config.items() if k in names}
        extra = {k: v for k, v in config.items() if k not in names}
        if 'color_dict' in kwargs:
            kwargs['color_dict'] = _alias_json_keys(kwargs['color_dict'])
        return cls(**kwargs, extra=extra)

    def to_dict(self):
        """
        The kwargs for get_tree_metadata (and the JSON form): only fields
        that differ from the defaults.
        """
        out = {}
        for f in fields(self):
            if f.name == 'extra':
                continue
            value = getattr(self, f.name)
            if value != f.default:
                out[f.name] = value
        out.update(self.extra)
        return out

    def to_json_dict(self):
        out = self.to_dict()
        if isinstance(out.get('color_dict'), dict):
            out['color_dict'] = _jsonable_color_dict(out['color_dict'])
        elif callable(out.get('color_dict')):
            raise ValueError(f"Can't save a callable color_dict for "
                             f"{self.label or self.column} to JSON")
        return out

    def build(self, md_table, seqs=None):
        """
        Make a TreeMetadata object from this spec and a genome metadata table
        (pandas indexed by genome, or polars with a 'genome' column).
        If seqs is given, only those genomes are included.
        """
        kwargs = self.to_dict()
        column = kwargs.pop('column')
        kwargs.pop('use_tree_seqs', None)
        if isinstance(column, str):
            md_table = to_pandas_indexed(md_table, 'genome')
            column = md_table[column]
        return get_tree_metadata(column, seqs=seqs, **kwargs)

    def color_fn(self, md_table, seqs=None):
        """ a callable from genome name to color """
        return self.build(md_table, seqs=seqs).get_item_color


@dataclass
class PlotConfig:
    """
    A parsed plot_metadata JSON. See module docstring.

    Unknown keys are kept (in genomes_extra and extra) so that to_dict() and
    to_json() round trip.
    """
    genome_metadata: List[Optional[MetadataSpec]] = field(default_factory=list)
    table_locs: Dict[str, Any] = field(default_factory=dict)
    leaf_color: Optional[MetadataSpec] = None
    genes: Dict[str, Any] = field(default_factory=dict)
    draw_args: Dict[str, Any] = field(default_factory=dict)
    genomes_extra: Dict[str, Any] = field(default_factory=dict)
    extra: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, config):
        config = dict(config)
        genomes = dict(config.pop('genomes', {}))
        metadata = [MetadataSpec.from_dict(m)
                    for m in genomes.pop('metadata', [])]
        leaf_color = MetadataSpec.from_dict(genomes.pop('leaf_color', None))
        table_locs = genomes.pop('table_locs', {})
        genes = config.pop('genes', {})
        draw_args = config.pop('draw_args', {})
        return cls(genome_metadata=metadata, table_locs=dict(table_locs),
                   leaf_color=leaf_color, genes=dict(genes),
                   draw_args=dict(draw_args), genomes_extra=genomes,
                   extra=config)

    def to_dict(self, _json=False):
        """ the plot_metadata dict that SyntenyTreePlotter reads """
        spec_dict = (lambda s: s.to_json_dict()) if _json \
            else (lambda s: s.to_dict())
        genomes = dict(self.genomes_extra)
        genomes['metadata'] = [None if s is None else spec_dict(s)
                               for s in self.genome_metadata]
        genomes['table_locs'] = dict(self.table_locs)
        if self.leaf_color is not None:
            genomes['leaf_color'] = spec_dict(self.leaf_color)
        out = dict(self.extra)
        out['genomes'] = genomes
        if self.genes:
            out['genes'] = dict(self.genes)
        if self.draw_args:
            out['draw_args'] = dict(self.draw_args)
        return out

    def to_json(self, path=None, **json_kwargs):
        """ write to path (if given) and return the JSON string """
        json_kwargs.setdefault('indent', 1)
        text = json.dumps(self.to_dict(_json=True), **json_kwargs)
        if path is not None:
            with open(path, 'w') as JSON:
                JSON.write(text)
        return text

    def get_spec(self, label_or_column):
        """ first spec with this label (or column) """
        for spec in self.genome_metadata:
            if spec is not None and label_or_column in (spec.label,
                                                        spec.column):
                return spec
        raise KeyError(label_or_column)

    def build_metadata(self, md_table, seqs=None):
        """ TreeMetadata objects for each spec (None for spacers) """
        return [None if s is None else s.build(md_table, seqs=seqs)
                for s in self.genome_metadata]


def load_plot_metadata(plot_metadata):
    """
    Read a plot_metadata JSON file (or an already loaded dict, or a
    PlotConfig) into a PlotConfig. None gives an empty config.
    """
    if plot_metadata is None:
        return PlotConfig()
    if isinstance(plot_metadata, PlotConfig):
        return plot_metadata
    if isinstance(plot_metadata, str):
        with open(plot_metadata) as JSON:
            plot_metadata = json.load(JSON)
    return PlotConfig.from_dict(plot_metadata)
