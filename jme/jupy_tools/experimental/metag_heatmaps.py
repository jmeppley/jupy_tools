"""
Clustered heatmaps of genome abundances in metagenomes, aggregated by sample
metadata (eg depth, month).

Generalized from the nanopore-2020_biller notebooks (phage mgs_002-004 and
VEIME Aloha_abundance_012-016). The typical figure is:

 * a genome dendrogram on top (one column per genome)
 * one or more heatmap panels below it, all in the dendrogram's column order
   (eg ALOHA by depth, ALOHA by month, TARA by depth)
 * panels can be "paired": each genome column is split by a sample category
   (eg size fraction) and one category is made negative so a single diverging
   colormap shows both (viral fraction positive, cellular fraction negative)
 * panels that aren't split use a sequential colormap (magma_r) by default
   (see color_scale in abundance_heatmap)
 * optional colored metadata strips (type, module, host, ...) under the panels
 * a "double" color bar with one side per size fraction
 * an optional signed log transform to bring out the low values

Data prep and plotting are kept separate. Data prep works on a long table of
coverages (polars or pandas) with columns: sample, genome, value and a sample
metadata table. Plotting takes pandas tables (rows = metadata levels,
columns = genomes, or a (genome, fraction) MultiIndex for paired panels).

All coloring is passed in as callables or dicts. Nothing here knows about
VEIMEs, hosts or modules.

Example::

    import polars as pl
    from jme.jupy_tools.experimental import metag_heatmaps
    from jme.jupy_tools.utils import save_fig_to_pdf

    # long table of hits: sample, genome, cov
    hits = pl.read_csv('phage.covs.tsv', separator='\\t')
    hits = metag_heatmaps.normalize_coverage(hits, read_counts)
    # -> adds cov_per_bn_reads

    # sample metadata: sample, depth, month, filter, ...
    md = pl.read_csv('aloha.sample.md.tsv', separator='\\t')
    md = md.filter(pl.col('cruise') < 298)  # eg public cruises

    # top 15 cellular + top 35 viral genomes (union, in rank order)
    genomes = metag_heatmaps.top_genomes(
        hits, {'00-20': 15, '20-02': 35}, group_col='filter', samples_md=md)

    # mean over ALL samples in each depth (absent = 0), split by size fraction
    #  with the cellular fraction negated
    split = dict(split_col='filter', split_negative='00-20',
                 split_order=['20-02', '00-20'])
    by_depth = metag_heatmaps.aggregate(hits, md, 'depth', genomes=genomes,
                                        **split)
    by_month = metag_heatmaps.aggregate(hits, md, 'month', genomes=genomes,
                                        order=range(1, 13), **split)

    # TARA: parse sample names to get depth classes
    tara_md = metag_heatmaps.parse_tara_samples(tara_hits['sample'].unique())
    by_tara = metag_heatmaps.aggregate(
        tara_hits, tara_md, 'depth', genomes=genomes,
        order=list(metag_heatmaps.TARA_DEPTHS))

    # optional log scale (keeps sign, 0 stays 0)
    (by_depth, by_month, by_tara), zero_val = metag_heatmaps.signed_log(
        [by_depth, by_month, by_tara])

    host_color_fn = metag_heatmaps.color_fn_from_dict(
        host_colors, key_fn=known_hosts.get, other='#808080')
    hm = metag_heatmaps.abundance_heatmap(
        [metag_heatmaps.Panel(by_depth, 'Depth (m)\\nALOHA'),
         metag_heatmaps.Panel(by_month, 'Month\\nALOHA',
                              row_labels=metag_heatmaps.MONTHS),
         metag_heatmaps.Panel(
             by_tara, 'Depth\\nTARA',
             row_labels=list(metag_heatmaps.TARA_DEPTHS.values()))],
        strips=[metag_heatmaps.MetadataStrip('Sequence Type', type_colors_fn),
                metag_heatmaps.MetadataStrip('Mapped Host', host_color_fn)],
        log_zero_val=zero_val,
        fraction_labels=('particles < 0.2 \\u03bcm',
                         'particles > 0.2 \\u03bcm'),
    )
    save_fig_to_pdf(hm.fig, 'Figure.4.pdf')  # (also takes a list of figs)

    hm.genome_order    # genomes in dendrogram order
    hm.axes['panel_0'] # matplotlib axes for further tweaks
"""
from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

import numpy
import pandas
import polars

MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
          'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# TARA depth class codes (in plotting order) and display names
TARA_DEPTHS = {'SUR': 'Surface', 'DCM': 'DCM', 'MES': 'Mesopelagic'}
TARA_SAMPLE_RE = r'^Station(\d+)_(MES|SUR|DCM)_(.+)$'

# ALOHA size fractions: viral is positive, cellular negative by convention
ALOHA_VIRAL = '20-02'
ALOHA_CELLULAR = '00-20'
FRACTION_LABELS = ('particles < 0.2 μm', 'particles > 0.2 μm')

DEFAULT_VALUE_COL = 'cov_per_bn_reads'

# default colormaps for abundance_heatmap(color_scale=...)
DEFAULT_DIVERGING_CMAP = 'PuOr'
# light at 0, since most cells are 0 or close to it
DEFAULT_SEQUENTIAL_CMAP = 'magma_r'


##########
# Data prep

def _to_polars(df):
    """ accept polars or pandas tables """
    if isinstance(df, polars.DataFrame):
        return df
    if isinstance(df, pandas.DataFrame):
        # (avoid polars.from_pandas, it needs pyarrow for string columns)
        return polars.DataFrame({c: df[c].tolist() for c in df.columns})
    raise TypeError(f"Expected a polars or pandas DataFrame, got {type(df)}")


def _to_pandas(df):
    """ polars -> pandas without needing pyarrow """
    return pandas.DataFrame(df.to_dict(as_series=False))


def normalize_coverage(hits, read_counts, value_col='cov',
                       sample_col='sample', reads_col='reads', per=1e9,
                       out_col=DEFAULT_VALUE_COL):
    """
    Normalize coverage by the number of reads in each sample.

    hits: long table with sample_col and value_col
    read_counts: dict or pandas.Series (sample -> reads) or a table with
                 sample_col and reads_col
    per: scale factor (1e9 gives coverage per billion reads)

    Returns a polars table with a new out_col. Samples without a read count
    are dropped.
    """
    hits = _to_polars(hits)
    if isinstance(read_counts, dict):
        read_counts = pandas.Series(read_counts)
    if isinstance(read_counts, pandas.Series):
        read_counts = polars.DataFrame({
            sample_col: [str(s) for s in read_counts.index],
            reads_col: read_counts.values,
        })
    read_counts = _to_polars(read_counts).select(sample_col, reads_col)
    read_counts = read_counts.with_columns(
        polars.col(sample_col).cast(hits.schema[sample_col]))
    return (
        hits.join(read_counts, on=sample_col, how='inner')
            .with_columns(
                (polars.col(value_col) / polars.col(reads_col) * per)
                .alias(out_col))
            .drop(reads_col)
    )


def long_from_wide(df, genome_cols, sample_col='sample',
                   value_col=DEFAULT_VALUE_COL, genome_col='genome'):
    """
    Convert a wide table (one row per sample, one column per genome; like the
    older *.bam_norm_covs.tsv files) to the long format used here.

    If df is a pandas table with no sample_col, the index is used.
    """
    if isinstance(df, pandas.DataFrame) and sample_col not in df.columns:
        df = df.rename_axis(sample_col).reset_index()
    df = _to_polars(df)
    return df.unpivot(index=sample_col, on=list(genome_cols),
                      variable_name=genome_col, value_name=value_col)


def _sample_groups(samples_md, group_col, sample_col):
    """ return {group: set(samples)} """
    md = _to_polars(samples_md).select(sample_col, group_col) \
        .drop_nulls()
    groups = {}
    for sample, group in md.iter_rows():
        groups.setdefault(group, set()).add(sample)
    return groups


def top_genomes(hits, n, group_col=None, samples_md=None, score='sum',
                sample_col='sample', genome_col='genome',
                value_col=DEFAULT_VALUE_COL):
    """
    Pick the most abundant genomes.

    n: an int, or a dict {group: n} to take the top n from each group of
       samples and return the union (eg {'00-20': 15, '20-02': 35})
    group_col: column to group samples by (in samples_md if given, else in
               hits). If n is an int and group_col is given, the top n of
               every group are taken.
    samples_md: sample metadata table. If given, only samples in it are used.
    score: how to rank genomes within a group: 'sum', 'mean' (absent counted
           as zero, so it ranks the same as sum), 'max' or 'median'
           (over samples with hits).

    Returns a list of genomes: the first group's list in rank order, then any
    new genomes from the next group, and so on. Groups are taken in the order
    of the n dict, or sorted if n is an int.
    """
    hits = _to_polars(hits)
    if samples_md is not None:
        md = _to_polars(samples_md)
        keep = md[sample_col].unique().to_list()
        hits = hits.filter(polars.col(sample_col).is_in(keep))
    else:
        md = None

    agg = {
        'sum': polars.col(value_col).sum(),
        'mean': polars.col(value_col).sum(),
        'max': polars.col(value_col).max(),
        'median': polars.col(value_col).median(),
    }[score]

    def _top(table, count):
        return (table.group_by(genome_col).agg(agg.alias('_score'))
                     .sort(['_score', genome_col], descending=[True, False])
                     .head(count)[genome_col].to_list())

    if group_col is None:
        if isinstance(n, dict):
            raise ValueError("n can only be a dict if group_col is given")
        return _top(hits, n)

    if md is not None and group_col in md.columns:
        groups = _sample_groups(md, group_col, sample_col)
        group_hits = {g: hits.filter(polars.col(sample_col).is_in(list(s)))
                      for g, s in groups.items()}
    else:
        group_hits = {g: hits.filter(polars.col(group_col) == g)
                      for g in hits[group_col].drop_nulls().unique()}

    if isinstance(n, dict):
        n_by_group = n
    else:
        n_by_group = {g: n for g in sorted(group_hits)}

    genomes = []
    seen = set()
    for group, count in n_by_group.items():
        if group not in group_hits:
            continue
        for genome in _top(group_hits[group], count):
            if genome not in seen:
                seen.add(genome)
                genomes.append(genome)
    return genomes


def aggregate(hits, samples_md, by, genomes=None, split_col=None,
              split_negative=None, split_order=None, order=None,
              sample_col='sample', genome_col='genome',
              value_col=DEFAULT_VALUE_COL, agg='mean'):
    """
    Aggregate per-sample values by a sample metadata column.

    The mean is taken over ALL samples in samples_md with that value of `by`,
    counting samples with no hits as 0 (the table is pivoted to
    sample x genome and filled with 0 before averaging).

    hits: long table (sample, genome, value)
    samples_md: one row per sample. Only these samples are used.
    by: metadata column to aggregate on (eg 'depth', 'month')
    genomes: genomes (columns) to keep, in order. Default: all genomes in
             hits, ranked by total value.
    split_col: optional metadata column (eg 'filter') to split columns on.
               The result then has (genome, split_value) MultiIndex columns.
    split_negative: split value (or list of values) to multiply by -1
    split_order: order of split values within each genome
                 (default: non-negated values first, then sorted)
    order: row order (and the set of rows). Levels with no samples are NaN.
           Default: sorted levels.
    agg: 'mean' (default) or 'sum', or 'median'

    Returns a pandas.DataFrame with levels of `by` as the index.
    """
    hits = _to_polars(hits)
    md = _to_polars(samples_md)
    md_cols = [sample_col, by] + ([split_col] if split_col else [])
    md = md.select(md_cols).drop_nulls().unique(subset=[sample_col])
    md_pd = _to_pandas(md).set_index(sample_col)

    # sample x genome table, zero filled for all samples in md
    summed = _to_pandas(
        hits.filter(polars.col(sample_col).is_in(md[sample_col].to_list()))
            .group_by([sample_col, genome_col])
            .agg(polars.col(value_col).sum()))
    wide = summed.pivot(index=sample_col, columns=genome_col,
                        values=value_col)
    if genomes is None:
        genomes = list(wide.sum().sort_values(ascending=False).index)
    wide = wide.reindex(index=md_pd.index, columns=list(genomes)) \
        .fillna(0.0)
    wide.columns.name = genome_col

    if split_col is None:
        result = wide.groupby(md_pd[by]).agg(agg)
    else:
        if split_negative is None:
            negatives = set()
        elif isinstance(split_negative, (list, tuple, set)):
            negatives = set(split_negative)
        else:
            negatives = {split_negative}
        if split_order is None:
            levels = sorted(md_pd[split_col].unique(), key=str)
            split_order = [s for s in levels if s not in negatives] + \
                          [s for s in levels if s in negatives]
        parts = {}
        for split_value in split_order:
            in_split = md_pd[split_col] == split_value
            part = wide[in_split].groupby(md_pd.loc[in_split, by]).agg(agg)
            if split_value in negatives:
                part = -1 * part
            parts[split_value] = part
        result = pandas.concat(parts, axis=1, names=['_split', genome_col])
        # we want (genome, split) columns ordered by genome then split
        result = result.swaplevel(axis=1)
        result.columns.names = [genome_col, split_col]
        result = result.reindex(
            columns=pandas.MultiIndex.from_product(
                [list(genomes), list(split_order)],
                names=[genome_col, split_col]))
        # a level missing from a split (eg no cellular samples at 1000m)
        result = result.fillna(0.0)

    if order is not None:
        result = result.reindex(list(order))
    else:
        result = result.sort_index()
    result.index.name = by
    return result


def paired_frame(positive, negative, names=('vf', 'cf'), negate=True):
    """
    Combine two already aggregated tables (rows = levels, cols = genomes)
    into one paired table with (genome, fraction) MultiIndex columns.

    The negative table is multiplied by -1 if negate is True.
    Missing genomes or levels are filled with 0.
    """
    genomes = list(positive.columns) + \
        [g for g in negative.columns if g not in set(positive.columns)]
    index = positive.index.union(negative.index, sort=False)
    pos = positive.reindex(index=index, columns=genomes).fillna(0.0)
    neg = negative.reindex(index=index, columns=genomes).fillna(0.0)
    if negate:
        neg = -1 * neg
    result = pandas.concat({names[0]: pos, names[1]: neg}, axis=1) \
        .swaplevel(axis=1)
    return result.reindex(
        columns=pandas.MultiIndex.from_product([genomes, list(names)]))


def signed_log(dfs, zero_frac=0.5, zero_val=None):
    """
    Log transform that keeps the sign of the values and maps 0 to 0:

        y = sign(x) * (log(|x| + zero_val) - log(zero_val))

    where zero_val is zero_frac times the smallest non-zero |x| in all of dfs
    (unless zero_val is given). Use the same zero_val for all panels of a
    figure, so they share one color scale.

    dfs: a DataFrame or a list of DataFrames

    Returns (transformed, zero_val), where transformed is a DataFrame or list
    to match the input. Use signed_exp() to go back.
    """
    single = isinstance(dfs, pandas.DataFrame)
    if single:
        dfs = [dfs]
    if zero_val is None:
        nonzero = [numpy.abs(df.values[numpy.isfinite(df.values)])
                   for df in dfs]
        nonzero = numpy.concatenate([v[v > 0] for v in nonzero])
        if len(nonzero) == 0:
            raise ValueError("No non-zero values to transform")
        zero_val = nonzero.min() * zero_frac
    out = [numpy.sign(df) * (numpy.log(df.abs() + zero_val) -
                             numpy.log(zero_val))
           for df in dfs]
    return (out[0] if single else out), zero_val


def signed_exp(values, zero_val):
    """ the inverse of signed_log """
    values = numpy.asarray(values, dtype=float)
    return numpy.sign(values) * \
        (numpy.exp(numpy.abs(values) + numpy.log(zero_val)) - zero_val)


def parse_tara_samples(samples, regex=TARA_SAMPLE_RE,
                       columns=('station', 'depth', 'type'),
                       sample_col='sample'):
    """
    Parse TARA sample names (eg Station85_MES_0d2-3) into a metadata table
    with columns: sample, station (int), depth, type.

    Samples that don't match the regex are dropped.
    """
    if isinstance(samples, (polars.Series, pandas.Series)):
        samples = list(samples)
    import re
    rexp = re.compile(regex)
    rows = []
    for sample in samples:
        match = rexp.search(str(sample))
        if match is not None:
            rows.append((sample,) + match.groups())
    md = polars.DataFrame(rows, schema=[sample_col] + list(columns),
                          orient='row')
    if 'station' in columns:
        md = md.with_columns(polars.col('station').cast(polars.Int64))
    return md


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
# Plotting

@dataclass
class Panel:
    """
    One heatmap panel.

    data: rows = levels (eg depths), columns = genomes, or a
          (genome, fraction) MultiIndex for a paired panel (as returned by
          aggregate(..., split_col=...))
    label: y axis label (eg 'Depth (m)\\nALOHA')
    row_labels: tick labels for the rows (default: the index)
    height: relative height (default: number of rows)
    cluster: use this panel to build the genome dendrogram
    """
    data: pandas.DataFrame
    label: str = ''
    row_labels: Optional[Sequence] = None
    height: Optional[float] = None
    cluster: bool = True

    @property
    def paired(self):
        return isinstance(self.data.columns, pandas.MultiIndex)

    @property
    def fractions(self):
        if not self.paired:
            return []
        return list(dict.fromkeys(self.data.columns.get_level_values(1)))

    @property
    def genomes(self):
        if self.paired:
            return list(dict.fromkeys(self.data.columns.get_level_values(0)))
        return list(self.data.columns)

    def cluster_table(self):
        """ rows = levels (x fractions), cols = genomes """
        if not self.paired:
            return self.data
        return pandas.concat({f: self.data.xs(f, axis=1, level=1)
                              for f in self.fractions}, axis=0)

    def ordered(self, genome_order):
        """
        data with columns in genome_order. Genomes missing from this panel
        are 0. (Rows that are all NaN, ie levels with no samples, stay NaN.)
        """
        if self.paired:
            columns = pandas.MultiIndex.from_product(
                [list(genome_order), self.fractions],
                names=self.data.columns.names)
        else:
            columns = pandas.Index(list(genome_order))
        data = self.data.reindex(columns=columns).astype(float)
        missing = ~columns.isin(self.data.columns)
        if missing.any():
            data.loc[:, missing] = 0.0
        return data


@dataclass
class MetadataStrip:
    """
    A row of colored boxes under the heatmaps, one per genome.

    color_fn: genome -> color, or a dict {genome: color}
    default: color for genomes missing from a dict
    """
    label: str
    color_fn: Union[Callable[[str], Any], Dict[str, Any]]
    default: Any = (1., 1., 1., 1.)

    def color(self, genome):
        if callable(self.color_fn):
            return self.color_fn(genome)
        return self.color_fn.get(genome, self.default)


@dataclass
class HeatmapFigure:
    """ what abundance_heatmap returns """
    fig: Any
    axes: Dict[str, Any]
    genome_order: List[str]
    linkage: Optional[numpy.ndarray] = None
    vmax: float = None
    panels: List[Panel] = field(default_factory=list)


def order_genomes(panels, metric='correlation', method='average',
                  genomes=None):
    """
    Cluster genomes on the panels' values (all panels with cluster=True,
    paired panels contribute one row per level and fraction).

    Returns (genome_order, linkage) (linkage is None if < 3 genomes).

    Genomes with constant values (eg all zero) have undefined correlation;
    those distances are set to the maximum distance in the matrix.
    """
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import pdist

    if genomes is None:
        genomes = list(dict.fromkeys(g for p in panels for g in p.genomes))
    tables = [p.cluster_table().reindex(columns=genomes)
              for p in panels if p.cluster]
    if not tables:
        tables = [p.cluster_table().reindex(columns=genomes)
                  for p in panels]
    table = pandas.concat(tables, axis=0, ignore_index=True).fillna(0.0)
    if len(genomes) < 3:
        return list(genomes), None
    with numpy.errstate(divide='ignore', invalid='ignore'):
        dists = pdist(table.T.values, metric=metric)
    bad = ~numpy.isfinite(dists)
    if bad.any():
        good = dists[~bad]
        dists[bad] = good.max() if len(good) else 1.0
    dists = numpy.clip(dists, 0, None)
    link = linkage(dists, method=method)
    return [genomes[i] for i in leaves_list(link)], link


def _draw_dendrogram(ax, link):
    """ draw dendrogram so leaf i is centered at x = i + .5 """
    from scipy.cluster.hierarchy import dendrogram
    dendro = dendrogram(link, no_plot=True)
    for xs, ys in zip(dendro['icoord'], dendro['dcoord']):
        ax.plot(numpy.array(xs) / 10, ys, color='k', lw=.75)
    ax.set_xlim(0, len(dendro['leaves']))
    ax.set_ylim(0, max(max(y) for y in dendro['dcoord']) * 1.02)
    ax.axis('off')
    return dendro


def _draw_panel(ax, panel, genome_order, cmap, norm, fontsize):
    data = panel.ordered(genome_order)
    n_genomes = len(genome_order)
    n_rows = data.shape[0]
    n_sub = max(1, len(panel.fractions))
    x_edges = numpy.linspace(0, n_genomes, n_genomes * n_sub + 1)
    y_edges = numpy.arange(n_rows + 1)
    values = numpy.ma.masked_invalid(data.values.astype(float))
    mesh = ax.pcolormesh(x_edges, y_edges, values, cmap=cmap, norm=norm)
    ax.set_xlim(0, n_genomes)
    ax.set_ylim(n_rows, 0)
    ax.set_yticks(numpy.arange(n_rows) + .5)
    row_labels = panel.row_labels if panel.row_labels is not None \
        else list(data.index)
    ax.set_yticklabels(row_labels, fontsize=fontsize, rotation=0)
    ax.set_xticks(numpy.arange(n_genomes) + .5)
    ax.tick_params(axis='x', labelbottom=False, length=0)
    ax.set_ylabel(panel.label, fontsize=fontsize)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return mesh


def _draw_strips(ax, strips, genome_order, fontsize):
    n = len(genome_order)
    X = numpy.arange(n) + .5
    for i, strip in enumerate(strips):
        ax.bar(X, numpy.ones(n), width=1, bottom=i,
               color=[strip.color(g) for g in genome_order])
    ax.set_xlim(0, n)
    ax.set_ylim(len(strips), 0)
    ax.set_yticks(numpy.arange(len(strips)) + .5)
    ax.set_yticklabels([s.label for s in strips], fontsize=fontsize)
    ax.set_xticks(X)
    ax.tick_params(axis='x', labelbottom=False, length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _log_ticks(vmax, zero_val, max_ticks=6, min_gap=.08):
    """
    nice real values (1, 2, 5 x 10^k) for a signed_log scale

    Ticks closer to 0 than min_gap (as a fraction of the bar) are dropped so
    their labels don't overlap the 0 label.
    """
    real_max = float(signed_exp(vmax, zero_val))
    if real_max <= 0:
        return [0.0]
    lo = numpy.floor(numpy.log10(zero_val))
    hi = numpy.ceil(numpy.log10(real_max))
    for mults in ([1], [1, 3], [1, 2, 5]):
        ticks = [m * 10 ** k for k in numpy.arange(lo, hi + 1)
                 for m in mults]
        ticks = [t for t in ticks
                 if zero_val < t <= real_max and
                 numpy.log(t + zero_val) - numpy.log(zero_val)
                 >= min_gap * vmax]
        if len(ticks) >= 3:
            break
    if len(ticks) > max_ticks - 1:
        step = int(numpy.ceil(len(ticks) / (max_ticks - 1)))
        ticks = ticks[::-1][::step][::-1]
    return [0.0] + ticks


def _tick_text(value):
    if value == 0:
        return '0'
    if abs(value) >= 1000 or abs(value) < .01:
        return f"{value:.0e}".replace('e+0', 'e').replace('e-0', 'e-')
    return f"{value:g}"


def double_color_bar(ax, cmap, vmax, n_steps=100, zero_val=None,
                     fraction_labels=FRACTION_LABELS, label=None,
                     fontsize=None, ticks=None):
    """
    Draw a two column color bar: the left column goes from 0 to vmax and the
    right one from 0 to -vmax, both drawn upward from 0, so each size fraction
    gets its own scale. Tick labels are absolute values.

    zero_val: if the data went through signed_log(), pass its zero_val so
              the ticks show untransformed values
    ticks: tick values (in real, untransformed units)
    """
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize
    cmap = plt.get_cmap(cmap)
    levels = numpy.linspace(0, vmax, n_steps + 1)
    mids = (levels[:-1] + levels[1:]) / 2
    grid = numpy.column_stack([mids, -mids])
    ax.pcolormesh([0, 1, 2], levels, grid, cmap=cmap,
                  norm=Normalize(-vmax, vmax))
    ax.set_xlim(0, 2)
    ax.set_ylim(0, vmax)

    if ticks is None:
        if zero_val is not None:
            ticks = _log_ticks(vmax, zero_val)
        else:
            from matplotlib.ticker import MaxNLocator
            ticks = [t for t in MaxNLocator(5).tick_values(0, vmax)
                     if 0 <= t <= vmax]
    ticks = list(ticks)
    if zero_val is not None:
        positions = signed_log(pandas.DataFrame({'t': ticks}),
                               zero_val=zero_val)[0]['t'].values
    else:
        positions = ticks
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')
    ax.set_yticks(positions)
    ax.set_yticklabels([_tick_text(t) for t in ticks], fontsize=fontsize)
    ax.set_xticks([.5, 1.5])
    ax.set_xticklabels(list(fraction_labels), rotation=90, fontsize=fontsize)
    if label:
        ax.set_ylabel(label, fontsize=fontsize)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return ax


def _data_sign(panels):
    """ 'positive', 'negative' or 'both' """
    if any(p.paired for p in panels):
        return 'both'
    values = numpy.concatenate([p.data.values.astype(float).ravel()
                                for p in panels])
    values = values[numpy.isfinite(values)]
    if (values >= 0).all():
        return 'positive'
    if (values <= 0).all():
        return 'negative'
    return 'both'


def _half_cmap(cmap, sign, n=256):
    """ the upper (positive) or lower (negative) half of a diverging cmap """
    from matplotlib import pyplot as plt
    from matplotlib.colors import ListedColormap
    cmap = plt.get_cmap(cmap)
    lo, hi = (.5, 1) if sign == 'positive' else (0, .5)
    return ListedColormap(cmap(numpy.linspace(lo, hi, n)),
                          name=f"{cmap.name}_{sign}")


def abundance_heatmap(panels, strips=(), cmap=None, vmax=None,
                      genome_order=None, metric='correlation',
                      method='average', dendrogram=True,
                      label_color_fn=None, label_rename_fn=None,
                      genome_labels=True, colorbar='auto',
                      color_scale='auto',
                      colorbar_label='Mean coverage\nper billion reads',
                      log_zero_val=None, fraction_labels=FRACTION_LABELS,
                      panel_letters=True, separators=True,
                      figsize=None, row_height=.25, dendrogram_height=4,
                      colorbar_width=None, colorbar_height=.6,
                      fontsize=14, title=None):
    """
    Draw stacked, column-aligned heatmap panels under one genome dendrogram.

    panels: list of Panel (or DataFrames). All panels share one x axis:
            genome i spans x = [i, i+1]; paired panels split that space
            into one sub-column per fraction.
    strips: list of MetadataStrip drawn under the panels
    cmap: colormap. Default: DEFAULT_DIVERGING_CMAP ('PuOr') for a
          diverging color_scale, DEFAULT_SEQUENTIAL_CMAP ('magma_r') for a
          sequential one.
    vmax: default: max |value| over all panels
    genome_order: use this order instead of clustering
    metric, method: for scipy clustering of genomes
    dendrogram: draw the dendrogram (if clustering)
    label_color_fn: genome -> color for the genome labels
    label_rename_fn: genome -> text for the genome labels
    genome_labels: label genomes under the bottom panel
    color_scale: 'sequential', 'diverging' or 'auto' (the default).
          'sequential': the whole cmap from 0 (light end of the default) to
              vmax. Data must not have mixed signs; all-negative data (eg a
              negated cellular fraction) is plotted as absolute values.
          'diverging': 0 is the middle of cmap, the scale is -vmax..vmax.
              If no panel is paired and all values have one sign, only that
              half of cmap is used (and shown in the color bar), so a
              single size fraction looks like its side of a paired figure.
          'auto': 'diverging' if any panel is paired (ie made with
              aggregate(..., split_col=...)) or the data has both signs,
              else 'sequential'.
    colorbar: 'double', 'single', None or 'auto' (double if any panel is
              paired, else single)
    log_zero_val: zero_val from signed_log, to label the color bar in real
                  units
    fraction_labels: (positive, negative) labels for the double color bar
    panel_letters: True for (A), (B), ... or a list of strings, or False
    separators: draw thin lines between genomes
    figsize: default from the number of genomes and rows
    row_height: inches per table row (when figsize is not given)
    dendrogram_height: dendrogram height, in table rows
    colorbar_width: inches for the color bar column (incl. its labels)
    colorbar_height: color bar height as a fraction of the panels' height

    Returns a HeatmapFigure with .fig, .axes (dict), .genome_order
    """
    from matplotlib import pyplot as plt
    from matplotlib.colors import Normalize

    panels = [p if isinstance(p, Panel) else Panel(p) for p in panels]
    strips = list(strips)

    link = None
    if genome_order is None:
        genome_order, link = order_genomes(panels, metric=metric,
                                           method=method)
    genome_order = list(genome_order)
    n_genomes = len(genome_order)
    show_dendro = dendrogram and link is not None

    if vmax is None:
        vmax = max(float(numpy.nanmax(numpy.abs(p.data.values)))
                   for p in panels)
        if not vmax > 0:
            vmax = 1.0
    # sign: which part of a signed scale the data needs: 'both',
    #  'positive' or 'negative'
    sign = _data_sign(panels)
    if color_scale == 'auto':
        color_scale = 'sequential' if sign != 'both' else 'diverging'
    if color_scale == 'sequential':
        if sign == 'both':
            raise ValueError("A sequential color_scale needs data with one "
                             "sign (or paired panels need 'diverging')")
        if sign == 'negative':
            panels = [replace(p, data=p.data.abs()) for p in panels]
            sign = 'positive'
        cmap = DEFAULT_SEQUENTIAL_CMAP if cmap is None else cmap
        norm = Normalize(0, vmax)
    elif color_scale == 'diverging':
        cmap = DEFAULT_DIVERGING_CMAP if cmap is None else cmap
        if sign == 'both':
            norm = Normalize(-vmax, vmax)
        else:
            cmap = _half_cmap(cmap, sign)
            norm = Normalize(0, vmax) if sign == 'positive' else \
                Normalize(-vmax, 0)
    else:
        raise ValueError(f"color_scale must be 'sequential', 'diverging' or "
                         f"'auto', not {color_scale!r}")

    if colorbar == 'auto':
        colorbar = 'double' if any(p.paired for p in panels) else 'single'

    # layout: one gridspec row per element, heights in "table rows"
    heights = []
    names = []
    if show_dendro:
        heights.append(dendrogram_height)
        names.append('dendrogram')
    for i, panel in enumerate(panels):
        heights.append(panel.height if panel.height is not None
                       else panel.data.shape[0])
        names.append(f'panel_{i}')
    if strips:
        heights.append(len(strips))
        names.append('strips')

    if figsize is None:
        width = max(6, .3 * n_genomes + 3)
        height = max(3, sum(heights) * row_height + 1.5)
        figsize = (width, height)
    if colorbar_width is None:
        # inches for the bar plus its tick labels and axis label
        colorbar_width = (1.0 if colorbar == 'double' else .6) + \
            .35 * fontsize / 10
    fig = plt.figure(figsize=figsize)
    cbar_cols = 1 if colorbar else 0
    gs = fig.add_gridspec(
        len(heights), 1 + cbar_cols, height_ratios=heights,
        width_ratios=[1] + ([colorbar_width / figsize[0]]
                            if colorbar else []),
        hspace=.08, wspace=.02 if colorbar else 0)

    axes = {}
    main_axes = []
    first = None
    for row, name in enumerate(names):
        ax = fig.add_subplot(gs[row, 0], sharex=first)
        first = first or ax
        axes[name] = ax
        if name != 'dendrogram':
            main_axes.append(ax)

    if show_dendro:
        _draw_dendrogram(axes['dendrogram'], link)

    for i, panel in enumerate(panels):
        ax = axes[f'panel_{i}']
        _draw_panel(ax, panel, genome_order, cmap, norm, fontsize)
        if separators:
            for x in range(1, n_genomes):
                ax.axvline(x, color='k', lw=.4)

    if strips:
        _draw_strips(axes['strips'], strips, genome_order, fontsize)
        if separators:
            for x in range(1, n_genomes):
                axes['strips'].axvline(x, color='k', lw=.4)

    # genome labels under the bottom axes
    bottom = main_axes[-1]
    if genome_labels:
        texts = [label_rename_fn(g) if label_rename_fn else str(g)
                 for g in genome_order]
        bottom.tick_params(axis='x', labelbottom=True, length=0)
        bottom.set_xticks(numpy.arange(n_genomes) + .5)
        bottom.set_xticklabels(texts, rotation=90, fontsize=fontsize)
        if label_color_fn is not None:
            for genome, label in zip(genome_order,
                                     bottom.get_xticklabels()):
                label.set_color(label_color_fn(genome))

    # line up the y labels of all panels
    fig.align_ylabels(main_axes)

    if colorbar:
        panel_rows = [r for r, n in enumerate(names) if n.startswith('panel')]
        holder = fig.add_subplot(gs[panel_rows[0]:panel_rows[-1] + 1, 1])
        holder.axis('off')
        # bar is the left part of its column (labels go to the right) and
        #  colorbar_height of the panels' height, centered
        bar_frac = (.4 if colorbar == 'double' else .25)
        cax = holder.inset_axes([0, (1 - colorbar_height) / 2,
                                 bar_frac, colorbar_height])
        axes['colorbar'] = cax
        if colorbar == 'double':
            double_color_bar(cax, cmap, vmax, zero_val=log_zero_val,
                             fraction_labels=fraction_labels,
                             label=colorbar_label, fontsize=fontsize)
        else:
            from matplotlib.cm import ScalarMappable
            cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
            if colorbar_label:
                cb.set_label(colorbar_label, fontsize=fontsize)
            if log_zero_val is not None:
                ticks = _log_ticks(vmax, log_zero_val)
                if sign == 'both':
                    ticks = [-t for t in ticks[:0:-1]] + ticks
                elif sign == 'negative':
                    ticks = [-t for t in ticks[::-1]]
                cb.set_ticks(signed_log(pandas.DataFrame({'t': ticks}),
                                        zero_val=log_zero_val)[0]['t'].values)
                cb.set_ticklabels([_tick_text(t) for t in ticks])
            cax.tick_params(labelsize=fontsize)

    if panel_letters:
        letters = panel_letters if not isinstance(panel_letters, bool) else \
            [f"({chr(65 + i)})" for i in range(len(main_axes))]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        # put letters left of the left-most y label
        lefts = [ax.get_tightbbox(renderer).x0 for ax in main_axes]
        x = fig.transFigure.inverted().transform((min(lefts), 0))[0] - .01
        for ax, letter in zip(main_axes, letters):
            y = ax.get_position().y1
            fig.text(x, y, letter, fontsize=fontsize, ha='right', va='top')

    if title:
        fig.suptitle(title, fontsize=fontsize)

    return HeatmapFigure(fig=fig, axes=axes, genome_order=genome_order,
                         linkage=link, vmax=vmax, panels=panels)

