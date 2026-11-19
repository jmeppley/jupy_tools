"""
Phylogenetic trees drawn next to genome metadata columns and gene synteny
maps.

Moved from jme.jupy_tools.experimental.vlp_plots_2 (TreePlotterVLP, now
SyntenyTreePlotter) and experimental.cluster_trees (TreePlotter, draw,
clade collapsing).

Inputs (see plot_metadata for the shared data model):

 * genome metadata table: one row per genome (pandas indexed by genome, or
   polars with a 'genome' column)
 * gene table: one row per gene with genome, gene_no, start, end, strand,
   annot[, desc] (see gene_maps.import_genes; polars needs a 'gene' column)
 * plot config: a plot_metadata JSON path, dict or PlotConfig listing the
   metadata columns to draw and their colours, plus gene colouring options
 * a Bio.Phylo tree

Example::

    from jme.jupy_tools import tree_plots

    plotter = tree_plots.SyntenyTreePlotter(md_df, genes_df,
                                            'plot_metadata.json',
                                            desc_col='desc')
    tree = tree_plots.read_tree('concat.treefile')
    # rectangular tree + metadata + gene maps
    fig = plotter.draw_tree('MCL 12', tree,
                            collapse=dict(check_clade_fn=my_label_fn),
                            gene_method=tree_plots.GM_TYPES)
    # circular tree with metadata rings
    fig, md_list = plotter.draw_circular_tree('MCL 12', tree)

Lower level:

 * TreePlotter(tree, seq_lens, gene_table, taxonomy).master_plot(...)
 * draw(tree, axes=ax): a plain rectangular phylogram (from Bio.Phylo.draw)
 * get_tree_seqs(tree_or_clade): set of leaf names
 * get_collapsed_nodes(tree, check_clade_fn, leaf_color_md)
"""
import json
import re
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import chain, combinations

import numpy
import pandas
import matplotlib.patches as mpatches
from matplotlib import patches
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import to_rgba
from scipy.cluster.hierarchy import dendrogram, linkage

from Bio import Entrez, MissingPythonDependencyError, Phylo, SeqIO
from Bio.Phylo import BaseTree

from jme.jupy_tools import gene_maps
from jme.jupy_tools.circular_trees import (
    _draw_circular_tree_recursive, draw_md_ring, get_ring_order,
)
from jme.jupy_tools.gene_maps import cluster_synteny_plot
from jme.jupy_tools.plot_metadata import (
    TreeMetadata, add_cat_md_legend, get_min_max, get_tree_metadata,
    load_plot_metadata,
)
from jme.jupy_tools.utils import first, read_tsv, to_pandas, to_pandas_indexed

Entrez.email = 'jmeppley@hawaii.edu'

GM_TOP = 'top'
GM_TYPES = 'types'
GM_MERGE = 'merge'

CN_REFS = 'refs'
CN_HOSTS = 'host'
CN_CUSTOM = 'cust'


class SyntenyTreePlotter():
    def __init__(
        self,
        md_tsv,
        gene_tsv,
        metadata_json=None,
        hosts_json=None,
        clusters_col='VC_cl',
        annot_col='annot_e',
        desc_col=None,
        length_col='Length',
        hostname_col='Host',
        seq_type_col='seq_type',
        ref_seq_type='Ref',
        contig_seq_type='Contig',
        override_gene_descs=None,
        parents=None,
        type_patterns=None,
        debug=False,
    ):
        self.debug = debug

        # load TSVs
        #  (if filenames are not passed, assume its a DataFrame)
        #  (polars tables are converted to pandas: genomes indexed by
        #   'genome' and genes by 'gene')
        self.md_df = read_tsv(md_tsv) if isinstance(md_tsv, str) \
                    else to_pandas_indexed(md_tsv, 'genome')
        self.gene_df = read_tsv(gene_tsv) if isinstance(gene_tsv, str) \
                    else _gene_table(gene_tsv)

        # we will look in here for metadata colors
        #  (a JSON file, a loaded dict, or a PlotConfig. This makes a copy, so
        #   the caller's dict isn't modified, and fixes "true"/"false"/"null"
        #   color_dict keys)
        self.plot_config = load_plot_metadata(metadata_json)
        self.plot_metadata = self.plot_config.to_dict()

        # apply custom column choices
        column_configs = self.plot_metadata.setdefault('genomes',
                                                       {}).setdefault('table_locs',
                                                                      {})
        hostname_col = column_configs.setdefault('hostname_col', hostname_col)
        length_col = column_configs.setdefault('length_col', length_col)
        seq_type_col = column_configs.setdefault('seq_type_col', seq_type_col)
        ref_seq_type = column_configs.setdefault('ref_seq_type', ref_seq_type)
        contig_seq_type = column_configs.setdefault('contig_seq_type', contig_seq_type)
        clusters_col = column_configs.setdefault('clusters_col', clusters_col)
        parents = column_configs.setdefault('parents', parents)

        # get seq lengths
        self.seq_lens = self.md_df[length_col].to_dict()

        # get reference geome name map witout any host assignments for novel
        # seqs
        self.ref_names = \
            self.md_df[
                self.md_df[seq_type_col] == ref_seq_type
            ][hostname_col].dropna()
        self.seq_type_col = seq_type_col
        self.contig_seq_type = contig_seq_type
        self.ref_seq_type = ref_seq_type
        self.parents = parents
        # gene type patterns: [(type_name, regex), ...]. Default to the ones
        #  in the plot config (the experimental version ignored those, so
        #  GM_TYPES never matched any types unless they were passed here)
        if not type_patterns:
            type_patterns = self.plot_metadata.get('genes', {}) \
                                .get('type_patterns', [])
        self.type_patterns = type_patterns if type_patterns else []
        if annot_col != "annot":
            if "annot" in self.gene_df.columns:
                self.gene_df = self.gene_df.drop('annot', axis=1)
            self.gene_df = self.gene_df.rename({annot_col: 'annot'}, axis=1)
        if desc_col is not None and desc_col != "desc":
            if "desc" in self.gene_df.columns:
                self.gene_df = self.gene_df.drop('desc', axis=1)
            self.gene_df = self.gene_df.rename({desc_col: 'desc'}, axis=1)

        self.clusters_col = clusters_col
        self.clusters = None

        # if there is no metadata requested, put some defaults
        self.genome_metadata_list = self.plot_metadata['genomes'].setdefault('metadata', [])
        if len(self.genome_metadata_list) == 0:
            self.genome_metadata_list.extend([
                            dict(column=length_col, log=True),
                            dict(column=seq_type_col),
            ])
            if hostname_col is not None:
                self.genome_metadata_list.append(dict(column=hostname_col))

        # make metadata objects
        self.metadata_metadata_all = [
            self.get_tree_metadata(**md_props)
            if md_props is not None
            else None
            for md_props in self.genome_metadata_list
        ]

        ## Color phylogram labels
        if 'leaf_color' in self.plot_metadata['genomes']:
            # User can speicy a TreeMetadatra config in
            # plot_metadaat['genomes']['leaf_color']
            md_params = self.plot_metadata['genomes']['leaf_color']
            if debug:
                print("Making leaf color md from:\n" + str(md_params))
            self.leaf_color_md = \
                self.get_tree_metadata(**md_params)
        else:
            # Otherwise we look for a seq_type metadata in the genome metadata
            leaf_color_col = seq_type_col
            # take the first MD that matches
            for i, md_props in enumerate(self.genome_metadata_list):
                if (
                    md_props is not None
                    and
                    md_props['column'] == leaf_color_col
                ):
                    leaf_md_index = i
                    break
            else:
                leaf_md_index = None
            if leaf_md_index is not None:
                # make a copy of the metadata object that defaults to black
                leaf_md_props = dict(md_props)
                leaf_md_props['null_color'] = 'black'
                if debug:
                    print("Making leaf color md from seq type params:\n" +
                          str(leaf_md_props))
                self.leaf_color_md = self.get_tree_metadata(**leaf_md_props)
            else:
                # if we find nothing, let everything be written in black
                if debug:
                    print("No leaf meatadarta found...making everything black")
                self.leaf_color_md = None
        leaf_color_fn = self.leaf_color_md.get_item_color \
                        if self.leaf_color_md is not None \
                        else None

        # merge metadata with actual data
        self.plot_data = \
            dict(
                seq_lens=self.seq_lens,
                gene_table=self.gene_df,
                taxonomy=None,
                ref_names=self.ref_names.to_dict(),
                leaf_color_fn=leaf_color_fn,
            )

        # desc overrides can come from JSON or paameters
        self.override_gene_descs = {} \
                if override_gene_descs is None \
                else override_gene_descs


    def debug_on(self):
        self.debug = True

    def debug_off(self):
        self.debug = False

    def get_plot_data(self):
        return {k: self.plot_data[k]
                for k in ['seq_lens',
                          'gene_table',
                          'taxonomy',
                          'ref_names',
                          'leaf_color_fn']}


    def get_clusters(self):
        if self.clusters is None:
            self.clusters = self.md_df \
                .groupby(self.clusters_col) \
                .apply(lambda df: set(df.index)) \
                .to_dict()
        return self.clusters

    def get_top_genes(self, tree, N=20):
        try:
            seqs = get_tree_seqs(tree)
        except AttributeError:
            # lets assume the user passed a collection of seq names
            seqs = tree
        genes = gene_maps.get_cluster_genes(self.gene_df,
                                               seqs,
                                               self.seq_lens)
        top_genes, _ = \
            gene_maps.plot_annot_positions(genes, ax=False,
                                              max_plotted_genes=N,
                                              max_colored_genes=N,
                                             )
        return top_genes

    def get_collapsed_nodes(
        self, tree,
        **kwargs
    ):
        """
        There are 2 ways to collapse nodes:
        
         * By default, collapse ref clades, breaking up by cluster if a clusters_col is set
         * Custom: if a check_clade_fn is passed, collapse only if that function returns a label
         
         """
        if "check_clade_fn" in kwargs:
            check_clade_fn = kwargs.pop("check_clade_fn")
            collapsed_nodes = get_collapsed_nodes(tree,
                                                  check_clade_fn,
                                                  self.leaf_color_md,
                                                  debug=self.debug,
                                                  **kwargs)
        else:
            collapsed_nodes = self.get_collapsed_nodes_refs(tree, **kwargs)
        
        if collapsed_nodes is not None:
            # set_clade_heights
            biggest_clade = max(len(c.get_terminals()) for c in collapsed_nodes)
            collapsed_node_heights = {
                c:numpy.interp(len(c.get_terminals()), [2, biggest_clade], [2,5])
                for c in collapsed_nodes
            }
            return collapsed_nodes, collapsed_node_heights
        
        # if there were no collapsed nodes
        return None, None
        
    def get_collapsed_nodes_refs(
        self, tree, 
        max_d_cutoff=1.5, 
        non_ref_label='Novel',
        min_ref_size=3,
        min_novel_size=8,
    ):
        # lets collapse:
        #  - any clade over 10 seqs with no ref in it
        #  - any clade with more than 3 refs all the same host
        collapsed_nodes = {}

        # constrain to clusters if there's a cluster col set
        cl_id_map = self.md_df[self.clusters_col] \
                    if self.clusters_col is not None \
                    else None

        if self.debug:
            counts = Counter()

        # loop over clades (exclude single terminals)
        for clade in tree.get_nonterminals():
            if self.debug:
                counts['clades_all'] += 1

            # no nested collapsed clades
            for c in collapsed_nodes:
                if c.is_parent_of(clade):
                    skip = True
                    break
            else:
                skip = False
            if skip:
                continue

            if self.debug:
                counts['clades_checked'] += 1

            # get the seq names from this clade
            terminals = get_tree_seqs(clade)
            # get the cluster of the first seq, if we have a cluster column set
            cl_id = cl_id_map.get(first(terminals), -1) \
                    if cl_id_map is not None \
                    else None

            # if there is a reference sequence in the clade
            if any(self.md_df[self.seq_type_col][t] == self.ref_seq_type for t in terminals):
                # ccollapse if:
                #  - the whole clade is ref seqs
                #  - (all in the same cluster)
                #  - at least 3 seqs
                if (
                    all(
                        (t in self.ref_names
                         and (
                                cl_id is None or
                                cl_id_map.get(t, -1) == cl_id
                         )
                        )
                        for t in terminals
                    )
                    and (len(terminals) >= min_ref_size)
                ):
                    if self.debug:
                        counts['collapsed_w_ref'] += 1
                    host_counts = get_host_counts(terminals,
                                                  self.ref_names,
                                                  self.parents,
                                                 )
                    count_str = "\n".join(f"{c} {h}" for h,c in sorted(host_counts.items()))
                    label = f"Published phage:\n{count_str}"
                    collapsed_nodes[clade] = label

                # don't collapse any clades with a mix of refs and novel
                continue

            # don't collapse novel clades w fewer than 8 seqs
            if len(terminals) < min_novel_size:
                continue

            depths = clade.depths()
            if not max(depths.values()):
                depths = clade.depths(unit_branch_lengths=True)
            max_d = max(depths[t] for t in clade.get_terminals())
            # collapse novel clades if:
            #  - they don't differ too much (deepest branch less than max_d_cutoff)
            #  - (they are all in the same cluster)
            if (
                max_d < max_d_cutoff
                and
                (
                    cl_id_map is None
                    or
                    all(cl_id_map.get(t, -1) == cl_id for t in terminals)
                )
            ):
                if self.debug:
                    counts['collapsed_no_ref'] += 1
                seq_type_counts = Counter(self.md_df[self.seq_type_col][t]
                                          for t in terminals)
                new_label = "\n".join(f"{count} novel {seq_type}s"
                                      for seq_type, count
                                      in sorted(seq_type_counts.items()))
                collapsed_nodes[clade] = new_label
                if len(seq_type_counts) > 1 and self.leaf_color_md is not None:
                    self.leaf_color_md.data_dict[clade] = non_ref_label

        if self.debug:
            print(counts)

        # return Nones if set of clades is empty or a single clade that is the whole tree
        if (
            (len(collapsed_nodes) == 0)
            or
            (
                (len(collapsed_nodes) == 1)
                and
                (len(next(iter(collapsed_nodes)).get_terminals()) == len(tree.get_terminals()))
            )
        ):
            return None

        return collapsed_nodes

    def get_type(self, gene_annot):
        return get_type(gene_annot, self.type_patterns)

    def standardize_descs(self, gene_df, annot_col='annot', desc_col='desc'):
        """
        for each gene/desc pair in the table:
            * replace the desc if it's in self.override_gene_descs
            * replace the desc if it matches one of the patterns in
                self.plot_metadata['type_patterns']
        """
        def get_new_desc(annot, desc):
            desc = self.override_gene_descs.get(annot, desc)
            gene_type = get_type(desc,
                                 self.plot_metadata.get('genes',
                                                        {}).get('type_patterns'))
            if gene_type is None:
                return desc
            return gene_type

        gene_df[desc_col] = [
            get_new_desc(annot, desc)
            for annot, desc
            in gene_df[[annot_col, desc_col]].values
        ]

        return gene_df

    def get_gene_colors_and_descs(self, seqs, gene_method=GM_TOP,
                                  gene_footnotes=None):
        """
        gene_method:
            GM_TOP (default):
                simply color by the top genes in this tree
            GM_MERGE:
                simplify gene descs with type patterns and use simpified names
                as annots, merging them together before laying out the plot.
            GM_TYPES:
                keep exising annots and descriptions, but color by
                type_patterns, trying to keep the same colors across plots
        """
        gene_colors = {}

        ## go through all genes in cluster and pull out any with the primary types
        # first ge tthe genes for just these seqs
        cluster_genes = gene_maps.get_cluster_genes(self.gene_df,
                                                       seqs,
                                                       self.seq_lens)

        # if merging by descs, standardize descs with type pattern names
        if gene_method == GM_MERGE:
            # use type_patterns and override_gene_descs!
            cluster_genes = self.standardize_descs(cluster_genes,)
            gene_col = 'desc'
        else:
            gene_col = 'annot'
            annot_descs = dict(tuple(v) for v in cluster_genes[[gene_col, 'desc']].values)
            if self.override_gene_descs is not None:
                annot_descs.update(self.override_gene_descs)

        # count the annotations
        non_na_genome_counts = Counter(
            a
            for g,a in set(tuple(p) for p in cluster_genes[['genome', gene_col]].values)
            if pandas.notna(a)
        )

        def generate_colors(starting_colors, remainder='black'):
            for color in starting_colors:
                yield color
            while True:
                yield remainder

        n_genes = self.plot_metadata['genes'].get('n_genes', 20)
        base_colors = self.plot_metadata['genes']['base_colors']
        n_colored_genes = min(n_genes, len(base_colors))
        n_base_colors = len(base_colors)
        type_colors = dict(zip(
            (p[0] for p in self.plot_metadata['genes']['type_patterns']),
            generate_colors(base_colors),
        ))

        if gene_method == GM_TYPES:
            # first, find the most common annotations with any of the primary types
            for gene, count in non_na_genome_counts.most_common():
                gene_annot = annot_descs[gene]
                gene_type = self.get_type(gene_annot)
                if gene_type is not None:
                    gene_colors[gene] = type_colors[gene_type]

                if len(gene_colors) >= n_base_colors:
                    break
            # use left over base and type colors before using black for
            # everything else
            color_iter = generate_colors([c for c in reversed(base_colors)
                                          if c not in gene_colors.values()])
        if gene_method == GM_MERGE:
            for (type_name, color) in type_colors.items():
                if type_name in non_na_genome_counts:
                    gene_colors[type_name] = color
            # use left over base and type colors before using black for
            # everything else
            color_iter = generate_colors([c for c in reversed(base_colors)
                                          if c not in gene_colors.values()])
        else:
            color_iter = generate_colors(base_colors)

        # Add from top_genes until we reach n_genes
        tree_top_genes = self.get_top_genes(seqs)
        top_gene_iter = iter(dict(non_na_genome_counts \
                                    .most_common()).keys())
        while len(gene_colors) < n_genes:
            try:
                top_gene = next(top_gene_iter)
            except StopIteration:
                break
            if top_gene not in gene_colors:
                gene_colors[top_gene] = next(color_iter)

        if gene_method == GM_MERGE:
            # return gene table with desc as annot
            #  (gene_col is 'desc' and make_gene_desc will skip including the
            #  annot if the gene name is found in the description)
            gene_labels = []
            gene_label_colors = {}
            for gene in cluster_genes[gene_col]:
                if gene in gene_colors:
                    #gene_label = make_gene_desc(gene, gene, non_na_genome_counts[gene])
                    gene_label = gene
                    #print(gene, gene_label, gene_colors[geene])
                    gene_label_colors[gene_label] = gene_colors[gene]
                else:
                    gene_label = gene
                gene_labels.append(gene_label)
            cluster_genes['annot'] = gene_labels
            # return gene table instead of descriptions dict
            return gene_label_colors, cluster_genes

        # otherwise return dict of descrtiptions for labelling genes
        if gene_footnotes is None:
            gene_footnotes = {}
            max_footnotes = 0
        else:
            max_footnotes = max(len(v) for v in gene_footnotes.values())
        gene_descs = {
            gene: make_gene_desc(gene, annot_descs[gene],
                                 non_na_genome_counts[gene],
                                 gene_footnotes.get(gene, []),
                                 max_footnotes,
                                )
            for gene in gene_colors
        }
        return gene_colors, gene_descs

    def get_lens_by_type(self, seqs):
        lens_by_type = defaultdict(list)
        for seq in seqs:
            seq_type = self.md_df[self.seq_type_col][seq]
            lens_by_type[seq_type].append(self.md_df.Length[seq])
        return lens_by_type

    def filter_tree_seqs(self, tree_seqs, cutoff=10, multiplier=3, buffer=.05):
        """
        Return the set of sequences such that:
            If there are at least 10 trusted sequences (AFVG or Ref):
                All non truseted seqs are:
                    within 3sigma of the trusted length distribtion
                    or within {buffer} of the trusted length max and min
            Otherwse, return all seqs
        """
        lens_by_source = self.get_lens_by_type(tree_seqs)
        good_lens = list(chain(*(lens_by_source[t]
                                 for t in lens_by_source
                                 if t != self.contig_seq_type)))
        if len(good_lens) < cutoff:
            return tree_seqs

        three_sigma_range = get_sigma_range(good_lens, multiplier)
        length_range = [min(three_sigma_range[0], min(good_lens) * (1-buffer)),
                        max(three_sigma_range[1], max(good_lens) * (1+buffer))]

        return [seq
                for seq in tree_seqs
                if (self.md_df[self.seq_type_col][seq] != self.contig_seq_type
                    or
                    length_range[0] <= self.md_df.Length[seq] <= length_range[1]
                   )
               ]

    def prune_outlier_contigs(self, tree):
        tree_seqs = get_tree_seqs(tree)
        filter_seqs = self.filter_tree_seqs(tree_seqs)
        for dropped in set(tree_seqs).difference(filter_seqs):
            tree.prune(dropped)

    def prune_unknowns(self, tree):
        """ remove refs with no known hosts """
        seqs = get_tree_seqs(tree)
        known_hosts = self.ref_names
        for seq in seqs:
            if seq in known_hosts and known_hosts[seq] == "Unknown":
                tree.prune(seq)

    def prune_duplicates(self, tree, cutoff=1e-06):
        known_hosts = self.ref_names
        for clade in get_redundant_clades(tree, known_hosts, cutoff):
            keep = None
            for leaf in sorted(clade.get_terminals(),
                               key=lambda l: l.name):
                if known_hosts.get(leaf.name, "Unknown") != "Unknown":
                    keep = leaf
                    break
            else:
                # move on if we can't find a good keeper
                continue

            for leaf in clade.get_terminals():
                if leaf != keep:
                    tree.prune(leaf)

    def prune_all(self, tree):
        self.prune_outlier_contigs(tree)
        self.prune_unknowns(tree)
        self.prune_duplicates(tree)

    def get_metadata_metadata(self, seqs, cl_id=None, drop_all_nas=True):
        """
        Prepends a cluster metadaata column if a specific cluster is flagged

        Pull metadata columns from metadata_metadata_all replacing any flagged
        columns with metadata rebuilt on just the tree sequences
        """
        metadata_list = []

        if self.clusters_col:
            # if a specific cluster is specified
            if cl_id is not None and cl_id != -1:
                # highlight indicated cluster with black
                cl_str = f'Cluster {cl_id}'
                self.in_cluster_dict = {
                    s:(cl_str if c == cl_id
                       else (None if c < 0
                             else "Other")
                      )
                    for s in seqs
                    for c in [self.md_df[self.clusters_col][s],]
                }
                # if any are NOT in the original cluster, add MD column with cluster status
                if any(v != cl_str for v in self.in_cluster_dict.values()):
                    metadata_list.append(TreeMetadata(
                        cl_str,
                        data_dict=self.in_cluster_dict,
                        null_value=None,
                        color_dict={cl_str: 'black', "Other": 'lightgrey', None: 'white'},
                    ))


        # if any metadata columns are flagged with "use_tree_seqs"
        # rebuild the metadata with only the clusters draw to limit the
        # number of colors needed
        #
        # also skip any MD that are all null for this set of seqs
        for i, md_props in enumerate(self.genome_metadata_list):
            md = self.metadata_metadata_all[i]

            # check for non null values
            try:
                if drop_all_nas and (md is not None) and (
                    all(md.get_item_value(s) in [None, md.null_value]
                        for s in seqs)
                ):
                    # if all the values for seqs are null, then skip this MD
                    continue
            except TypeError as e:
                print(md.label, len(seqs), md.null_value)
                print(dict(Counter(type(s) for s in seqs).most_common()))
                raise e

            # filter seqs if asked
            if (
                md_props is not None
                and
                md_props.get('use_tree_seqs', False)
            ):
                # by passing seqs, we limit the number of unique clusters
                # drawn
                metadata_list.append(
                    self.get_tree_metadata(seqs=seqs, **md_props)
                )
            else:
                metadata_list.append(md)

        return metadata_list

    def get_tree_metadata(self, column, use_tree_seqs=False, **kwargs):
        """
        A wrapper for get_tree_metadata that
        pulls out the column string from kwargs and replaces with
        the actual column from md_df

        Note, the param use_tree_seqs is ignored here, but its pulled out so it
        doesn't pass on to the TreeMetadata constructor. It should be used
        upstream  of this method.
        """
        if isinstance(column, str):
            column = self.md_df[column]
        return get_tree_metadata(column, **kwargs)

    def draw_tree(self, tree_name, tree,
                  collapse=True, cl_id=None, draw_genes=None,
                  gene_method=GM_TOP,
                  length_hist=False,
                  gene_footnotes=None,
                  tree_kwargs={},
                  draw_legends=True,
                  debug=False,
                  **kwargs):
        """
        Given a phylogenetic tree and genome metadata, make a synteny tree plot (using cluster_trees.master_plot)

        draw_genes: is on by default with a sorting function to prioritize drawing genomes closest
                    to the mean length. Set to True to just draw the first few genomes in each clade
                    or a different sorting function. Set to false to just draw the tree and metadata
                    with no gene locations.

        collapse: if True (default) try to collapse clades of all novel or all ref sequences
                  can also be a dict of one of two types:
                      * map from node to labels
                      * set of parameters for the get_collapsed_nodes() method
                  or just pass a cutoff value (max_d_cutoff) for the method

        gene_method: if "top" (default), draw the most common genes
                     if "types", prioritize genes that match the configured
                        type patterns
                     if "merge", merge genes with the same types and
                        descriptoins
        """

        # root at midpoint (ignore error if we got a clade that can't be re-rooted)
        try:
            tree.root_at_midpoint()
        except AttributeError:
            pass

        tree_seqs = get_tree_seqs(tree)
        if draw_genes is None:
            # default is sort by distance from mean length
            draw_genes = get_distance_to_mean_length_function(tree_seqs, self.seq_lens)

        if collapse:
            # is it a collapsed node label dict?
            if (
                isinstance(collapse, dict)
                and
                isinstance(next(iter(collapse)), BaseTree.Clade)
            ):
                # then we just add the heights
                collapsed_nodes = collapse
                # set clade heights and check that nodes do not overlap
                used_nodes = set()
                biggest_clade = max(len(c.get_terminals()) for c in collapsed_nodes)
                if (
                    "collapsed_clade_heights" in kwargs or
                    "collapsed_clade_heights" in tree_kwargs
                ):
                    if (
                        "collapsed_clade_heights" in kwargs and
                        "collapsed_clade_heights" in tree_kwargs
                    ):
                        print("WARNING: using clade heights form kwargs, not"
                              " tree_kwargs!!")
                    if "collapsed_clade_heights" in kwargs:
                        collapsed_node_heights = kwargs.pop("collapsed_clade_heights")
                        tree_kwargs['collapsed_clade_heights'] = collapsed_node_heights
                    else:
                        collapsed_node_heights = tree_kwargs.pop("collapsed_clade_heights")
                else:
                    collapsed_node_heights = {}
                    for c in collapsed_nodes:
                        seqs = get_tree_seqs(c)
                        if any(s in used_nodes for s in seqs):
                            raise Exception("Collapsed Clades cannnot overlap with each other!")
                        used_nodes.update(seqs)
                        collapsed_node_heights[c] = numpy.interp(len(seqs), [2, biggest_clade], [2,5])
            else:
                # we run get_collapsed_nodes with requested parameters
                if isinstance(collapse, dict):
                    # they passed a dict of params to use
                    collapse_kwargs = collapse
                elif collapse is True:
                    # use the default params
                    collapse_kwargs = {}
                else:
                    # use the collapse value as the cutoff
                    collapse_kwargs = {'max_d_cutoff': collapse}
 
                collapsed_nodes, collapsed_node_heights = \
                    self.get_collapsed_nodes(tree, **collapse_kwargs)
        else:
            collapsed_nodes, collapsed_node_heights = None, None

        # get list of nodes in this tree
        tree_nodes = tree_seqs
        if collapsed_nodes:
            # remove hidden nodes and replace with clade if collapsed
            for clade in collapsed_nodes:
                tree_nodes = tree_nodes.difference({t.name for t in
                                                    clade.get_terminals()})
                tree_nodes.update(clade)
        md_md = self.get_metadata_metadata(tree_nodes, cl_id)

        # pull master_plot args from plot_metadata, but allow kwarg overrides
        for k,v in self.plot_metadata.get('draw_args', {}).items():
            kwargs.setdefault(k, v)

        plot_data = self.get_plot_data()

        if draw_genes:
            # gene colors (TODO: This needs to be updated!)
            if (
                'desc' in self.gene_df
                and
                'base_colors' in self.plot_metadata.get('genes', {})
            ):
                # if gene_method is merge, gene descs will be a modified gene_table
                # otherwise, it's a dict from gene to description
                gene_colors, gene_descs = self.get_gene_colors_and_descs(tree_seqs,
                                                                         gene_method,
                                                                         gene_footnotes,
                                                                        )
                if gene_method == GM_MERGE:
                    plot_data['gene_table'] = gene_descs
                else:
                    kwargs.setdefault('add_gene_descs', gene_descs)
                kwargs.setdefault('gene_colors', gene_colors)
            elif 'color_map' in self.plot_metadata.get('genes', {}):
                kwargs.setdefault('gene_colors', self.plot_metadata['genes']['color_map'])

            n_genes = self.plot_metadata.get('genes', {}).get('n_genes', 20)
            kwargs.setdefault('max_colored_genes', n_genes)
            kwargs.setdefault('max_plotted_genes', n_genes)

        tree_plotter = TreePlotter(tree, **plot_data)
        fig = tree_plotter.master_plot(tree_name,
                                       md_md,
                                       collapsed_clade_labels=collapsed_nodes,
                                       draw_genes=draw_genes,
                                       tree_kwargs=tree_kwargs,
                                       **kwargs,
                                      )

        # make some legends (this will fail for draw_genes=False)
        if draw_genes == False or draw_genes is None:
            if draw_legends:
                print("WARNING: I don't know how to make legnds when draw_genes is False!!")
            return fig, md_md
        
        if not draw_legends:
            return fig, md_md
        else:
            # figure out where the plots are
            try:
                _ = fig.axes[0].bbox._bbox.get_points()
                def get_bounds(ax):
                    x1, x2 = ax.bbox._bbox.get_points()
                    return (*x1, x2[0] - x1[0], x2[1] - x1[1])
            except AttributeError:
                get_bounds = lambda ax: ax.figbox.bounds

            x, y, w, h = get_bounds(fig.axes[0])
            x1, y1, w1, h1 = get_bounds(fig.axes[1])
            x2, y2, w2, h2 = (
                get_bounds(fig.axes[2])
                if 'clade_ratio_dicts' not in kwargs
                else
                get_bounds(fig.axes[4])
            )
            top = y2 + h2

            # separate categorical and continuous MDs
            cat_mds = []
            cont_mds = []

            for md in md_md:
                if md is None:
                    continue
                label = md.label
                try:
                    # If draw_legends is a container, on draw things in it
                    if label not in draw_legends:
                        continue
                except TypeError:
                    # draw_legends was just a boolean
                    pass
                if md.color_method == 'cat':
                    vcs = list(md.get_category_colors(tree_nodes))
                    cat_mds.append((label, vcs))
                else:
                    cont_mds.append((label, md))

            N_rows = len(cont_mds) + 2
            cont_ax_step = 1 / N_rows
            cont_ax_h = cont_ax_step / 5

            lc_count = 0
            for label, md in cont_mds:
                cax = fig.add_axes([x/2, top, w * .4, h2 * cont_ax_h], )
                cax.patch.set_alpha(0)
                #_ = cax.set_title(label)
                _ = cax.yaxis.set_label_position('right')
                _ = cax.set_ylabel(label, rotation=0, va='top', ha='left')
                cb = fig.colorbar(ScalarMappable(norm=md.color_method, cmap=md.color_map,),
                                 cax=cax,
                                 orientation='horizontal')

                # set cbar limit to just values in this plot
                minv, maxv = get_min_max(md.data_dict[s]
                                         for s in tree_seqs
                                         if s in md.data_dict)
                _ = cax.set_xlim(minv, maxv)
                # erase the outline of the full (not clipped) colorbar
                cb.outline.set_visible(False)

                top -= h2 * cont_ax_step
                cax.set_clip_on(True)

            lax = None
            lax_x = x/2
            lax_w = w * .4 / 3
            top -= h2 * cont_ax_step
            for label, values_colors in sorted(cat_mds, key=lambda lvc:
                                               len(lvc[1]), reverse=True):
                if len(values_colors) == 1:
                    value, color = next(iter(values_colors))
                    if value in {None, "None", md.null_value}:
                        continue
                lax = fig.add_axes([lax_x, top, lax_w, h2 * cont_ax_step], )
                lax.patch.set_alpha(0)
                for spine in lax.spines.values():
                    spine.set_visible(False)
                lax.set_yticks([])
                lax.set_xticks([])
                add_cat_md_legend(
                           label,
                           values_colors,
                           lax,
                           bbox_to_anchor=None,
                           loc='upper center')
                lax_x += lax_w

            if length_hist:
                # make an axis in the upper right
                lax = fig.add_axes([x1, y2 + (h2/2), w1, h2/2])
                _ = lax.yaxis.tick_right()
                _ = lax.set_title('Genome Lengths')
                if isinstance(length_hist, dict):
                    # interperet as a dict of length sets to plot multiple
                    # histograms from
                    hist_range = get_min_max(chain(*length_hist.values()))
                    for label, seq_lens in length_hist.items():
                        h = lax.hist(seq_lens, bins=50, log=True, label=label,
                                     alpha=.66, range=hist_range,
                                     histtype='step')
                    _ =  lax.legend()
                else:
                    # other wise just make a hist of tree seq lengths
                    tree_seq_lens = [self.seq_lens[s] for s in tree_seqs]
                    h = lax.hist(tree_seq_lens, bins=25, log=True)

        return fig

    def draw_circular_tree(
        self,
        tree_name,
        tree, 
        figsize=[12,12], 
        bg_styles_dict={},
        line_styles_dict={},
        collapsed_dict={},
        scale_lines_step=1,
        ring_frac=.025,
        min_collapsed_size=2,
        support_md=None
    ):

        # get metadata
        tree_seqs = get_tree_seqs(tree)
        md_md = self.get_metadata_metadata(tree_seqs)

        fig, ax = plt.subplots(1,1,figsize=figsize)

        ring_order = get_ring_order(tree.root, collapsed_dict,
                                    min_collapsed_size)
        leaf_arc = 355 / len(ring_order)
        r = _draw_circular_tree_recursive(ax, tree.root, leaf_arc=leaf_arc, 
                           collapsed_dict=collapsed_dict,
                           support_md=support_md,
                           min_collapsed_size=min_collapsed_size,
                           bg_styles_dict=bg_styles_dict,
                           line_styles_dict=line_styles_dict)


        # we drew a tree out to a radius of r
        # let's add some scale lines in the form of circles of unit radii
        for i in range(0, int(numpy.ceil(r)), scale_lines_step):
            ax.add_patch(patches.Circle((0,0),
                                        radius=i, 
                                        lw=.5, ec='grey', alpha=.7,
                                        fill=None, zorder=-2))



        # ring width should be ring_frac of tree radius
        ring_width = r * ring_frac
        
        for mdmd in md_md:
            if mdmd is not None:
                actual_ring_width = draw_md_ring(
                    ax, ring_order, mdmd, 
                    radius=r, 
                    width=ring_width, 
                    leaf_arc=leaf_arc
                )
                r += actual_ring_width
            else:
                r += ring_width / 2

        _ = ax.set_xlim(-r,r)
        _ = ax.set_ylim(-r,r)
        _ = ax.set_xticks([])
        _ = ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        _ = ax.set_title(tree_name, pad=30)
        
        return fig, md_md


def get_distance_to_mean_length_function(seqs, seq_lens):
    """ return function that calculates:
            negative (-) of the abs distance from mean length for given sequence name or tree leaf

        use returned funtion as key for sorting by disrtance from mean length
    """
    seqs_mean_length = numpy.mean([seq_lens[s] for s in seqs])

    # dynamically create function with mean length and seq_len dict already embedded
    def distance_to_mean_length(seq):
        try:
            seq = seq.name
        except:
            pass
        return -1 * numpy.absolute(seq_lens[seq] - seqs_mean_length)

    return distance_to_mean_length

def generate_annot_desc_pairs(gene_table, override_gene_descs={}):
    for annot, desc in gene_table[['annot', 'desc']].values:
        desc = override_gene_descs.get(annot, desc)
        yield annot, desc

def get_type_mask(type_filter, cluster_genes, override_gene_descs={}):
    has_type = {}
    if not callable(type_filter):
        if isinstance(type_filter, str):
            type_filter = re.compile(type_filter)
        # assume it's a compiled pattern if not collable or string
        def is_type(annot, desc):
            if pandas.isna(desc):
                return False
            return type_filter.search(desc) is not None
    else:
        # assume callable filter takes two arguments
        def is_type(annot, desc):
            try:
                return type_filter(annot, desc)
            except:
                # if the data i null, assume that's the cause of the error
                if pandas.isna(annot) or pandas.isna(desc):
                    return False
                else:
                    # otherwise pass the error through
                    raise

    for genome, genes in cluster_genes.groupby('genome'):
        if any(is_type(a,d) for a,d
               in generate_annot_desc_pairs(genes,
                                            override_gene_descs)
              ):
            has_type[genome] = True
    return has_type

def get_host_counts(refs, known_hosts, parents=None):
    host_counts = Counter(re.sub('^Other ', '', known_hosts[r]) for r in refs)

    # merge any child taxon counts with parent counts if there are other members of the parent clade
    if parents is not None and len(parents) > 0:
        for host, count in list(host_counts.items()):
            # if a host has a parent and there are no other counted hosts with the same parent
            if (
                host in parents
                and
                any(parents[host] == parents.get(h,h)
                    for h in host_counts
                    if h != host)
            ):
                # collect counts under the parent
                host_counts[parents[host]] += count
                del (host_counts[host])

    if 'Unknown' in host_counts:
        del host_counts['Unknown']

    return host_counts

def make_gene_desc(gene, desc, count, footnotes, footnote_size):
    gene_id_str = f' ({gene})' if re.search(gene, desc) is None else ''
    desc = f'{desc}{gene_id_str} (in {count} genomes)'
    for n in range(footnote_size, -1, -1):
        char = " " if n >= len(footnotes) else footnotes[n]
        desc += char
    return desc

def get_tree_seqs(tree):
    return {t.name for t in tree.get_terminals()}

def get_type(desc, type_patterns):
    """ classify gene based on it's description """
    if pandas.isna(desc):
        return None
    for name, patt in type_patterns:
        if re.search(patt, desc, re.I) is not None:
            return name
    else:
        return None

def get_sigma_range(lens, multiplier=3):
    lstd = numpy.std(lens)
    lmean = numpy.mean(lens)
    return (lmean - (multiplier*lstd),
            lmean + (multiplier*lstd))

def get_redundant_clades(tree, known_hosts=None, cutoff=1e-06):
    small_clades = []
    for clade in tree.get_nonterminals():
        if all(clade.distance(t) < cutoff
               for t in clade.get_terminals()):
            if known_hosts and all(t.name in known_hosts
                                   for t in clade.get_terminals()):
                small_clades.append(clade)
    return small_clades


def intersection(row, feature):
    return min(row.end - feature.location.start, feature.location.end - row.start) + 1

def get_desc_overrides(seqs, known_hosts, gene_df):
    try:
        seqs = get_tree_seqs(seqs)
    except:
        pass

    from urllib.request import HTTPError

    matches = defaultdict(Counter)
    failures, attempts = 0, 0
    for seq in seqs:
        if seq in known_hosts:
            attempts += 1
            try:
                r = Entrez.efetch(db='nucleotide', id=seq, rettype='gb', retmode='text')
                record = first(SeqIO.parse(r, 'gb'))
            except HTTPError:
                failures += 1
                continue

            for gene, row in gene_df.query(f'genome == "{seq}"').iterrows():
                if pandas.isna(row.annot):
                    continue
                if re.search(r'^(Hypothetical|EggNOG cluster|[Uu]ncharacterized)',
                             row.desc) is None:
                    continue
                for feature in record.features:
                    if feature.type == 'CDS':
                        ovl = intersection(row, feature)
                        if ovl > 50:
                            for d in feature.qualifiers['product']:
                                if re.search('hypothetical protein', d) is None:
                                    if ovl * 2 > numpy.abs(row.end - row.start):
                                        matches[row.annot][d] += 1

    if failures:
        print(f"WARNING: {failures} of {attempts} searches failed!")
    return {a:m.most_common()[0][0] for a,m in matches.items() if len(m) == 1}


def get_collapsed_nodes(
    tree, 
    check_clade_fn,
    leaf_color_md,
    debug=False,
):
    """
    Collapse clades.

    The params are:

     * check_clade_fn: a callable that returns:
         None if the clade should not be collapsed
         If collapsed, then either:
             label: just the name for the calde
             label, color_label: the name and the value to look up in leaf_color_md
    """
        
    collapsed_nodes = {}

    if debug:
        counts = Counter()

    # loop over clades (exclude single terminals)

    for clade in tree.get_nonterminals():
        if debug:
            counts['clades_all'] += 1

        # no nested collapsed clades
        for c in collapsed_nodes:
            if c.is_parent_of(clade):
                skip = True
                break
        else:
            skip = False
        if skip:
            continue

        # check the max depth
        # (commented out for now: max_d isn't used yet, and the depths are
        #  slow to compute for every clade of a big tree)
        # depths = clade.depths()
        # if not max(depths.values()):
        #     depths = clade.depths(unit_branch_lengths=True)
        # max_d = max(depths[t] for t in clade.get_terminals())

        if debug:
            counts['clades_checked'] += 1

        # get the seq names from this clade
        label = check_clade_fn(clade)
        if label is None:
            continue
        
        if debug:
            counts['clades_collapsed'] += 1
            counts['leaves_collapsed'] += clade.count_terminals()
    
        if isinstance(label, tuple) and len(label) == 2:
            label, leaf_color_value = label
            
            # add clade to leaf_color_md object
            if leaf_color_md is not None:
                leaf_color_md.data_dict[clade] = leaf_color_value
            
        collapsed_nodes[clade] = label
                    
    if debug:
        print(counts)

    # return None if set of clades is empty or a single clade that is the whole tree
    if (
        (len(collapsed_nodes) == 0)
        or
        (
            (len(collapsed_nodes) == 1)
            and
            (len(next(iter(collapsed_nodes)).get_terminals()) == len(tree.get_terminals()))
        )
    ):
        return None

    return collapsed_nodes


def draw_tree_metadata(ax, y_posns, metadata_metadata, md_font_size=7, x_labels_on_top=False, thickness=1, collapsed_clade_heights=defaultdict(lambda: 2)):
    """
    Draw metadata columns
    
    params:
        ax: matplotlib axes object to draw in
        y_posns: dict
            map from Phylo tree nodes or clades to y locations
        metadata_metadata: list of TreeMetadata objects
            
    """

    # plot metadata as heatmap (using bars())
    bottoms, heights, tops, nodes = [], [], [], []
    for node in y_posns:
        y = y_posns[node]
        top = y + thickness/2
        if isinstance(node, str):
            height = thickness
        elif node.is_terminal():
            height = thickness
            node = node.name
        else:
            height = thickness * collapsed_clade_heights[node]
        nodes.append(node)
        heights.append(height)
        bottom = top - height
        bottoms.append(bottom)
        tops.append(top)

    x_pos = -1
    x_labels = {}
    for metadata in metadata_metadata:
        if metadata is None:
            x_pos += .5
            continue
        
        x_pos += 1
        x_labels[x_pos] = metadata.label

        colors = [metadata.get_item_color(g) for g in nodes]
        _ = ax.bar(x_pos, 
                   height=heights, 
                   width=1, 
                   bottom=bottoms, 
                   ec='grey',
                   linewidth=.1,
                   color=colors)
    
    # x axis labels
    _ = ax.set_xlim([-1,max(x_labels) + .5])
    _ = ax.set_xticks(list(x_labels.keys()))
    if x_labels_on_top:
        _ = ax.set_xticklabels(x_labels.values(), rotation=-90, fontsize=md_font_size)
        _ = ax.xaxis.tick_top()
    else:
        _ = ax.set_xticklabels(x_labels.values(), rotation=90, fontsize=md_font_size)
        
    # clear y axis labels    
    ticks = []
    _ = ax.set_yticks(ticks)

    
def add_label_refs(ax, y_positions, y_axis_names=None, **kwargs):
    """ put selected names on the right size of a plot 
    
        ax: matplotlib axes
        y_positions: map from node objects to y positions
        y_axis_names: optinal map from object to name or (sub)set of names to plot
            terminal nodes are expected to be indexed by name. Non-terminal nodes, by object.
        kwargs: arguments (eg fontsize) for set_yticklabels.
    """
    if y_axis_names:
        if isinstance(y_axis_names, dict):
            def get_label(genome):
                return y_axis_names[genome]
        else:
            def get_label(genome):
                return genome
        def use_genome(genome):
            return genome in y_axis_names
    else:
        def get_label(genome):
            return genome
        def use_genome(genome):
            return True
        

    ticks, labels = [], []
    for node in y_positions:
        if isinstance(node, str):
            height = 1
            genome = node
        elif node.is_terminal():
            height = 1
            genome = node.name
        else:
            height = 2
            genome = node
        if use_genome(genome):
            y = y_positions[node]
            top = y + .5
            bottom = top - height
            
            ticks.append((bottom + top) / 2)
            labels.append(get_label(genome))

    _ = ax.yaxis.tick_right()
    _ = ax.set_yticks(ticks)
    kwargs.setdefault('fontsize', 7)
    _ = ax.set_yticklabels(labels, **kwargs)

# draw a tree with the option of collapsed nodes
# (modified from Bio.Phylo.draw())
def draw(
    tree,
    label_func=str,
    do_show=True,
    show_confidence=True,
    # For power users
    axes=None,
    branch_labels=None,
    branch_label_x_delta=0.025,
    branch_label_y_delta=0.25,
    label_colors=None,
    collapsed_clade_labels=None,
    collapsed_clade_heights=None,
    *args,
    **kwargs
):
    """Plot the given tree using matplotlib (or pylab).

    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.

    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).

    Example using the pyplot options 'axhspan' and 'axvline'::

        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})

    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).

    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            A dict should map node labels to colors.
            A function should expect clade objects.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.
        collapsed_clade_labels: dict
            A map from clade objects to the label for the collapsed clade. Labels 
            default to first node plus a count if only height given for a clade.
        collapsed_clade_heights: dict
            A map from clade to the height of the collapsed clade marker. Defaults
            to 1 if only a label given for a clade.

    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MissingPythonDependencyError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    import matplotlib.collections as mpcollections

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []
    collapsed_polycollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)

    if not branch_labels:
        if show_confidence:

            def format_branch_label(clade):
                try:
                    confidences = clade.confidences
                    # phyloXML supports multiple confidences
                except AttributeError:
                    pass
                else:
                    return "/".join(conf2str(cnf.value) for cnf in confidences)
                if clade.confidence is not None:
                    return conf2str(clade.confidence)
                return None

        else:

            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(clade):
            return branch_labels.get(clade)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels

    # options for displaying label colors.
    if label_colors:
        if callable(label_colors):

            def get_label_color(clade):
                return label_colors(clade)

        else:
            # label_colors is presumed to be a dict
            def get_label_color(clade):
                return label_colors.get(label_func(clade), "black")

    else:

        def get_label_color(clade):
            # if label_colors is not specified, use black
            return "black"

    # make sure collapsed node dicts agree
    if collapsed_clade_heights or collapsed_clade_labels:
        # flag any sub clades to be skipped
        all_collapsed_clades = set()
        if collapsed_clade_heights:
            all_collapsed_clades.update(collapsed_clade_heights)
            def get_cc_height(clade):
                try:
                    return collapsed_clade_heights[clade]
                except:
                    return 1
        else:
            def get_cc_height(clade):
                return 1
        if collapsed_clade_labels:
            all_collapsed_clades.update(collapsed_clade_labels)
            def get_cc_label(clade):
                try:
                    return collapsed_clade_labels[clade]
                except:
                    return f"{clade.count_terminals()} tips"
        else:
            def get_cc_label(clade):
                return f"{clade.count_terminals()} tips"

        redundant_collapsed_clades = set()
        for n1 in all_collapsed_clades:
            for n2 in all_collapsed_clades:
                if n1 != n2:
                    if n1.is_parent_of(n2):
                        redundant_collapsed_clades.add(n2)
                    elif n2.is_parent_of(n1):
                        redundant_collapsed_clades.add(n1)
        good_collapsed_clades = all_collapsed_clades.difference(redundant_collapsed_clades)

        collapsed_clade_labels = {c:get_cc_label(c) for c in good_collapsed_clades}
        collapsed_clade_heights = {c:get_cc_height(c) for c in good_collapsed_clades}
    else:
        collapsed_clade_labels = {}
        collapsed_clade_heights = {}
        

    # Layout

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths


        """
        ## TODO: this is never reached!
        # map terminals to any collapsed parent clades
        # and count how many rows the collapsed clades will take up

        # did we get a mapping or collection?
        try:
            # if it has .items(), assume values() are row heights
            collapsed_clade_heights = {c:h for c,h in collapsed_clades.items() 
                                       if c not in redundant_collapsed_clades}
        except AttributeError:
            # otherwise all collapsed clades take up one row
            collapsed_clade_heights = {c:1 for c in collapsed_clades 
                                       if c not in redundant_collapsed_clades}
        """

    def get_y_positions(tree, collapsed_clade_heights):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        if collapsed_clade_heights:
            # flag all terminals to be hidden
            collapsed_terminals = {}
            cumulative_collapsed_height = 0
            for clade, height in collapsed_clade_heights.items():
                for tip in clade.get_terminals():
                    collapsed_terminals[tip] = clade
                cumulative_collapsed_height += height

            # total height is # all leaves - hidden leaves + collapsed clades
            maxheight = (tree.count_terminals()
                         + cumulative_collapsed_height
                         - len(collapsed_terminals))

            # Rows are defined by the tips
            heights = {}
            row = 0
            seen_clades = set()
            for tip in tree.get_terminals():
                if tip in collapsed_terminals:
                    clade = collapsed_terminals[tip]
                    if clade not in seen_clades:
                        seen_clades.add(clade)
                        row += 1 + collapsed_clade_heights[clade]
                        heights[clade] = row
                else:
                    row += 1
                    heights[tip] = row
        else:
            # with no collapsed clades, total height is number of leaves
            maxheight = tree.count_terminals()
            # Rows are defined by the tips
            heights = {
                tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
            }

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    x_posns = get_x_positions(tree)
    y_posns  = get_y_positions(tree, collapsed_clade_heights)
    terminals = []
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for axes: %s" % axes)

    def draw_clade_triangle(
        x_left=0,
        x_right=1,
        y_bot=0,
        y_top=1,
        color="black",
        lw=.1,
    ):
        """
        Create a triangle representing a collapsed clade
        """
        verts = [
            [x_right, y_bot],
            [x_left, y_top],
            [x_right, y_top],
        ]
        collapsed_polycollections.append(
            mpcollections.PolyCollection(
                verts=[verts,],
                closed=True,
                edgecolors=[color,],
                linewidths=[lw,],
                facecolors=[(0.,0.,0.,0.),],
            )
        )

    def draw_clade_lines(
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw
                )
            )
        elif orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw
                )
            )


    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        if clade not in collapsed_clade_heights:
            # Add label above the branch (optional)
            conf_label = format_branch_label(clade)
            if conf_label:
                axes.text(
                    #0.5 * (x_start + x_here),
                    x_here - branch_label_x_delta,
                    y_here - branch_label_y_delta,
                    conf_label,
                    fontsize="small",
                    #horizontalalignment="center",
                    horizontalalignment="right",
                    verticalalignment='bottom',
                )
            # Add node/taxon labels
            label = label_func(clade)
            if label not in (None, clade.__class__.__name__):
                axes.text(
                    x_here,
                    y_here,
                    " %s" % label,
                    verticalalignment="center",
                    color=get_label_color(clade),
                )
            if clade.clades:
                # Draw a vertical line connecting all children
                y_top = y_posns[clade.clades[0]]
                y_bot = y_posns[clade.clades[-1]]
                # Only apply widths to horizontal lines, like Archaeopteryx
                draw_clade_lines(
                    orientation="vertical",
                    x_here=x_here,
                    y_bot=y_bot,
                    y_top=y_top,
                    color=color,
                    lw=lw,
                )
                # Draw descendents
                for child in clade:
                    draw_clade(child, x_here, color, lw)
            else:
                terminals.append(clade)
        else:
            # draw triangle for collapsed clade
            terminals.append(clade)
            max_x = max(x_posns[t] for t in clade.get_terminals())
            d_y = collapsed_clade_heights[clade]
            draw_clade_triangle(
                x_left=x_here,
                x_right=max_x,
                y_bot=y_here - d_y,
                y_top=y_here,
                color=color,
                lw=lw,
            )
            # Add clade labels
            label = collapsed_clade_labels[clade]
            if label not in (None, ""):
                axes.text(
                    max_x,
                    y_here - (d_y / 2),
                    " %s" % label,
                    verticalalignment="center",
                    color=get_label_color(clade),
                )

    draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)
    for i in collapsed_polycollections:
        axes.add_collection(i)

    # Aesthetics

    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            axes.set_title(name)
    axes.set_xlabel("branch length")
    axes.set_ylabel("taxa")
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if do_show:
        plt.show()
        
    return {t:y_posns[t] for t in terminals}

# turn taxonmy into linkage for matplotlib
def _get_tax_linkage_recursive(node, linkage_list, node_index_dict):
    i = max(node_index_dict.get(node, -1), node_index_dict.get(node.name,-1))
    if i >= 0:
        return i
    
    if len(node.children) > 0:
        # if this is not a leaf
        # process the children if we get this far
        links = []
        for child in node.children:
            if child == node:
                continue
            i = _get_tax_linkage_recursive(child, linkage_list, node_index_dict)
            if i >= 0:
                links.append(i)
                
                             

        if len(links) == 0:
            return -1
    
        if len(links) == 1:
            return links[0]

        # there are at least two nodes to connect
        N=len(node_index_dict)
        d0 = numpy.double(1)
        while len(links) > 1:
            i1 = links.pop()
            i2 = links.pop()
            clust_size = 0
            d = 0
            for i in [i1, i2]:
                if i < N:
                    clust_size += 1
                else:
                    clust_size += linkage_list[i - N][3]
                    d += linkage_list[i - N][2]
            j = N + len(linkage_list)
            linkage_list.append([i1, i2, d + d0, clust_size])
            links.append(j)
            
            d0 = 0
        return j
    return -1

def get_tax_linkage(root_node, leaves):
    """ given a taxonom startig at root_node, build a linkage to the leaves (given as node names)
        returns a tuple: (linkage, leaf_index_dict), where:
            linkage: is a linkage table suiatble for scipy.cluster.hierarchy.dendrogram
            leaf_index_dict: is a map from leaf names to index in the linkage
    """
    taxlinkage = []
    leaf_index_dict = {leaf:i for i, leaf in enumerate(leaves)}
    _get_tax_linkage_recursive(root_node, taxlinkage, leaf_index_dict)
    return numpy.array(taxlinkage), leaf_index_dict

def generate_collapsed_nodes(tree, cluster_dict, 
                             max_node_size=0,
                             exceptions={}):
    """
    Pick clades to collapse in tree. Finds maximal clades where all terminals are in the same cluster.
    
    params:
        tree: a Bio.Phylo tree
        cluster_dict: map from leaf node names to cluster ids
        exceptions: set of node names to never collapse (eg reference seqs).
    
    yields: 2-tuple
        (clade object to collapse, cluster id)
    """
    collapsed_node_labels = {}
    collapsed_leaves = set()
    for node in tree.get_nonterminals():
        leaves = node.get_terminals()
        if len(leaves) == tree.count_terminals():
            # don't collapse the whole tree
            continue
        if max_node_size and len(leaves) > max_node_size:
            # skip nodes that are too big
            continue
        # skip clades that contain already used leaves or excluded leaves
        if any(((t in collapsed_leaves) or (t.name in exceptions))
               for t in leaves):
            continue

        # get the cluster id of all leaves in clade
        cluster_ids = set(cluster_dict.get(t.name, -1) for t in leaves)
        if len(cluster_ids) == 1:
            # single cluster clade, use it!
            yield (node, next(iter(cluster_ids)))
            
            # don't use its leaves in any other returned clades
            collapsed_leaves.update(leaves)


def get_tree_plotter_data(genomes_tsv, 
                          genes_tsv,
                          taxdump_dir,
                          refs_tsv,
                          leaf_color_fn=lambda l: 'black',
                         ):
    from edl.taxon import readTaxonomy
    genome_df = read_tsv(genomes_tsv)
    gene_table = read_tsv(genes_tsv)
    taxonomy = readTaxonomy(taxdump_dir)
    ref_taxa = read_tsv(refs_tsv)

    return dict(seq_lens=genome_df.length.to_dict(),
                gene_table=gene_table,
                ref_taxa=ref_taxa,
                genome_df=genome_df,
                taxonomy=taxonomy,
                ref_names=ref_taxa.name.to_dict(),
                leaf_color_fn=leaf_color_fn,
               )

class TreePlotter():
    def __init__(self, tree, 
                 seq_lens,
                 gene_table,
                 taxonomy,
                 ref_names={},
                 leaf_color_fn=lambda l: 'black',
                 **kwargs
                ):
        '''If tree is a dictionary, it is assumed to be a distance matrix 
           and a hierarical clustering tree will be generated
           
           (tree dict should look like {(genome1, genome2): dist12, (genome1, genome3): dist13, ...})
           
           '''
        self.tree = tree
        self.is_hier = isinstance(tree, dict)
        if self.is_hier:
            self.leaves = set(chain(*tree.keys()))
        else:
            self.leaves = set(self.tree.get_terminals())
        self.seq_lens = seq_lens
        self.gene_table = gene_table
        self.taxonomy = taxonomy
        self.gray_map = plt.get_cmap('binary')
        self.highlight_map = plt.get_cmap('Greens')    
        self.ref_names = ref_names
        self.leaf_color_fn = leaf_color_fn

    @lru_cache(maxsize=None)
    def get_taxon(self, name):
        for taxon in self.taxonomy.idMap.values():
            if taxon.name == name:
                return taxon
        return None

    def master_plot(self, title, 
                    metadata_metadata=None,
                    clade_ratio_dicts=None,
                    draw_genes=False,
                    max_plotted_genes=None,
                    max_colored_genes=10,
                    collapsed_clade_labels=None,
                    fig_width=None, 
                    tree_font_size=7,
                    axis_font_size=12,
                    md_font_size=9, 
                    gene_colors=None,
                    add_gene_descs='desc',
                    synteny_anchor_annot=None,
                    synteny_anchor_side=None,
                    synteny_anchor_strand=None,
                    tax_tree_height=1/3,
                    md_width_factor=0.02,
                    items_per_inch=6,
                    gene_kwargs={},
                    tree_kwargs={},
             ):

        """
        Do the work!

        params:
            draw_genes:
                False (default): only draw tree and metadata
                True: Draw synteny plot over aligned genome maps
                callabale: pick representatives of clades to draw using draw_genes() as the sorting key. This let's you chose the best representatives for collapsed clades.
        """

        # set some defualts
        if max_plotted_genes is None:
            if isinstance(gene_colors, dict):
                max_plotted_genes = len(gene_colors)
            else:
                max_plotted_genes = max_colored_genes
        n_genes = max_plotted_genes if isinstance(max_plotted_genes, int) \
                  else len(max_plotted_genes)

        # calculate the figure and subplot sizes

        # default figure width is bigger if we're drawing genes
        fig_width = fig_width if fig_width else (20 if draw_genes else 12)
        
        # how much vertical space for the primary tree
        if collapsed_clade_labels:
            # are there custom clade heights?
            if "collapsed_clade_heights" in tree_kwargs:
                N = len(self.leaves)
                for c in collapsed_clade_labels:
                    n = c.count_terminals()
                    h = tree_kwargs['collapsed_clade_heights'].get(c, 2)
                    N = N - n + h
            else:
                N = len(self.leaves) \
                     - sum(c.count_terminals() 
                          for c in collapsed_clade_labels) \
                     + 2 * len(collapsed_clade_labels)
        else:
            N = len(self.leaves)
        tree_height = N/items_per_inch
        
        # how much vertical space for genes/tax tree?

        # zero out tree height if no med plotted
        if clade_ratio_dicts is None:
            tax_tree_height = 0
        
        gene_height = n_genes/items_per_inch if draw_genes else 0         
        fig_height = tree_height + max(gene_height, tax_tree_height)   
        hr_tax_tree = tax_tree_height / tree_height
        hr_gene = 0 \
            if (tax_tree_height >= gene_height) \
            else (gene_height - tax_tree_height) / tree_height
            
        # genes take up zero width if we're not darwing them
        wr_genes = 1 if draw_genes else 0
        
        # relative width of metadata bars
        wr_a = md_width_factor * (len(metadata_metadata) + 1) if metadata_metadata else 0
        wr_b = md_width_factor * (len(clade_ratio_dicts) + 1) if clade_ratio_dicts else 0

        # set up grid
        fig = plt.figure(figsize=[fig_width, fig_height], constrained_layout=False)
        if self.is_hier: 
            # hier tree goes on the right
            gs = fig.add_gridspec(3, 4, 
                                  width_ratios=[wr_genes, wr_a, wr_b, .15],
                                  height_ratios=[hr_gene, hr_tax_tree, 1])

            # create only the subplots we're going to use 
            ax_tree = fig.add_subplot(gs[2,3])
            ax_md1 = fig.add_subplot(gs[2,1]) if metadata_metadata else None
            ax_md2 = fig.add_subplot(gs[2,2]) if clade_ratio_dicts else None
            ax_tax = fig.add_subplot(gs[1,2]) if clade_ratio_dicts else None
            ax_gene = fig.add_subplot(gs[0:2,0]) if draw_genes else None
            ax_synt = fig.add_subplot(gs[2,0]) if draw_genes else None 
        else:
            # tax tree goes on the left (it contains labels)
            gs = fig.add_gridspec(3, 4, 
                                  width_ratios=[1, wr_genes, wr_a, wr_b],
                                  height_ratios=[hr_gene, hr_tax_tree, 1])

            # create only the subplots we're going to use 
            ax_tree = fig.add_subplot(gs[2,0])
            ax_md1 = fig.add_subplot(gs[2,2]) if metadata_metadata else None
            ax_md2 = fig.add_subplot(gs[2,3]) if clade_ratio_dicts else None
            ax_tax = fig.add_subplot(gs[1,3]) if clade_ratio_dicts else None
            ax_gene = fig.add_subplot(gs[0:2,1]) if draw_genes else None
            ax_synt = fig.add_subplot(gs[2,1]) if draw_genes else None
        
        # display the title (differnt subplot depending on layout)
        title_ax = ax_gene if draw_genes else ax_tree        
        _ = title_ax.set_title(title, fontsize=axis_font_size)

        # this is a tight layout
        plt.subplots_adjust(wspace=0, hspace=0)

        y_posns, thickness = self.tree_plot(
            ax_tree, ax_tax,
            ax_md1, ax_md2, 
            metadata_metadata,
            clade_ratio_dicts,
            collapsed_clade_labels=collapsed_clade_labels,
            md_font_size=md_font_size,
            tree_font_size=tree_font_size,
            axis_font_size=axis_font_size,
            **tree_kwargs,
        )

        if draw_genes:
            if not callable(draw_genes):
                draw_genes = lambda genome: True
            
            if self.is_hier:
                genome_ys = y_posns
            else:
                genome_ys = {}
                for clade, y in y_posns.items():
                    if clade.is_terminal():
                        # it's a terminal, eg a single genomes
                        if draw_genes(clade) or draw_genes(clade.name):
                            genome_ys[clade.name] = y
                    else:
                        clade_height = tree_kwargs.get('collapsed_clade_heights',
                                                       defaultdict(lambda: 2))[clade]
                        leaves_to_draw = sorted(
                            [leaf.name
                             for leaf in clade.get_terminals()],
                            key=draw_genes,
                            reverse=True,
                        )[:int(numpy.floor(clade_height))]

                        pad = (clade_height - len(leaves_to_draw)) / 2
                        next_y = y - pad - .5
                        for l in leaves_to_draw:
                            genome_ys[l] = next_y
                            next_y -= 1
        
            cluster_synteny_plot(genome_ys,
                                 ax_gene, ax_synt, 
                                 genome_lens=self.seq_lens,
                                 gene_table=self.gene_table,
                                 max_plotted_genes=max_plotted_genes,
                                 max_colored_genes=max_colored_genes,
                                 max_genomes=len(genome_ys),
                                 add_gene_descs=add_gene_descs,
                                 set_titles=False,
                                 label_genomes=self.is_hier,
                                 gene_colors=gene_colors,
                                 anchor_annot=synteny_anchor_annot,
                                 anchor_side=synteny_anchor_side,
                                 anchor_strand=synteny_anchor_strand,
                                 thickness=thickness / 2,
                                 **gene_kwargs,
                                )
            _ = ax_synt.set_ylim(*ax_tree.get_ylim())
            
        return fig

    def _find_shadow_clades(self, clade_names):
        # collect all the ancestors of the clades we have counts for
        ancestors = defaultdict(list)
        clades = set()
        for clade_name in clade_names:
            clade = self.get_taxon(clade_name)
            clades.add(clade)
            for ancestor in clade.getLineage()[:-1]:
                ancestors[ancestor].append(clade)

        clade_shadows = {}
        shadows_to_clade = {}
        # are any clades ancestors of other clades?
        for clade_name in clade_names:
            clade = self.get_taxon(clade_name)
            # if this clade is an ancestor of another clade
            if clade in ancestors:
                # then: find a child that's not a targeted clade or an ancestor of one
                for child in clade.children:
                    if child not in ancestors and child not in clades:
                        # use as shadow clade for collecting other counts
                        clade_shadows[clade_name] = child.name
                        shadows_to_clade[child.name] = clade_name
                        break
                else:
                    dummy_name = f".{clade_name}.shadow"
                    from edl.taxon import TaxNode
                    dummy = TaxNode(
                        taxid=max(self.taxonomy.idMap.keys()) + 1,
                        parentid=clade.id,
                        rank=None
                    )
                    dummy.name = dummy_name
                    dummy.setParent(clade)
                    self.taxonomy.idMap[dummy.id] = dummy
                    print(f'WARNING: there is no unused child taxon for'
                          f'{clade_name} to use as a shadow. We are '
                          f'creating a dummy child taxon ({dummy_name}) '
                          'for making the tree.')
                    # use as shadow clade for collecting other counts
                    clade_shadows[clade_name] = dummy_name
                    shadows_to_clade[dummy_name] = clade_name
                    
        return clade_shadows, shadows_to_clade
            
    def tree_plot(self, ax_tree, ax_tax, ax_md1, ax_md2,
                  metadata_metadata, clade_ratio_dicts,
                  collapsed_clade_labels=None,
                  tree_font_size=7,
                  axis_font_size=10,
                  md_font_size=9,
                  internal_ref_labels=None,
                  **kwargs,
                 ):

        if internal_ref_labels is None:
            if self.is_hier:
                # default for hierachical trees is external ref labels
                internal_ref_labels = False
            else:
                # default for phylo trees is internal ref labels
                internal_ref_labels = True

        if internal_ref_labels:
            def label_func(node):
                if node.name:
                    if node.name in self.ref_names:
                        return f"{self.ref_names[node.name]} ({node.name})"                    
                    else:
                        return node.name
                return None
        else:
            label_func = str

        # draw tree
        if self.is_hier: 
            # hierarchical clustering
            
            genomes = self.leaves
            distances = self.tree
            # restructure for numpy
            genome_dists = [
                distances.get((g1,g2), distances.get((g2,g1), 1))
                for g1, g2
                in combinations(genomes, 2)
            ]
            # draw
            genome_linkage = linkage(genome_dists, method='ward')
            dn = dendrogram(genome_linkage, orientation='right', ax=ax_tree)
            x0,x1 = ax_tree.get_ylim()
            step = (x1-x0)/len(genomes)
            positions = numpy.arange(step/2, x1, step)
            y_posns = {g:positions[dn['leaves'].index(i)] 
                       for i, g in enumerate(genomes)}
        else:
            # Bio.Phylo tree
            plt.rcParams["font.size"] = tree_font_size
            y_posns = draw(self.tree,
                           axes=ax_tree,
                           do_show=False,
                           label_colors=self.leaf_color_fn,
                           collapsed_clade_labels=collapsed_clade_labels,
                           label_func=label_func, **kwargs)
            step = 1
            _ = ax_tree.set_ylabel("Genome", fontsize=axis_font_size)
            _ = ax_tree.set_xlabel("Branch Length", fontsize=axis_font_size)
            #  This is already done in the draw() function and it's less buggy
            #max_depth = max(self.tree.depths().values())
            #_ = ax_tree.set_xlim(-.5, max_depth + 1.5)
            
        _ = ax_tree.set_yticks([])

        # this may be overridden next if we have metadata
        right_most_ax = ax_tree
        
        # draw simple metadata
        if metadata_metadata:
            draw_tree_metadata(ax_md1, y_posns,
                               metadata_metadata,
                               md_font_size=md_font_size,
                               x_labels_on_top=not bool(clade_ratio_dicts),
                               thickness=step,
                               collapsed_clade_heights=kwargs.get(
                                   "collapsed_clade_heights",
                                   defaultdict(lambda: 2)
                               ),
                              )
            _ = ax_md1.set_ylim(*ax_tree.get_ylim())
            if not self.is_hier: 
                right_most_ax = ax_md1

        # draw taxonomic metadata
        if clade_ratio_dicts:
            clade_to_shadow, shadow_to_clade = self._find_shadow_clades(clade_ratio_dicts.keys())
            taxlinkage, leaf_index_dict = get_tax_linkage(self.taxonomy.root, 
                                                       [clade_to_shadow.get(c,c) 
                                                        for c in clade_ratio_dicts.keys()])
            index_leaf_dict = {i:l for l,i in leaf_index_dict.items()}
            d = dendrogram(taxlinkage, ax=ax_tax)
            order_order = [index_leaf_dict[i]
                           for i in d['leaves']]

            x0 = min(ax_tax.get_xticks())
            x1 = max(ax_tax.get_xticks())
            dx = x1 - x0
            dmd = len(order_order)
            buffer = dx * (1/dmd)
            _ = ax_tax.set_xlim(x0 - buffer, x1 + buffer)
            _ = ax_tax.set_yticks([])
            _ = ax_tax.set_xticks([])


            max_ratio = max(max(clade_ratio_dicts[o].values()) for o in clade_ratio_dicts.keys())
            metadata_metadata_2 = []
            for o in order_order:
                if o in {}:
                    cmap = self.highlight_map
                else:
                    cmap = self.gray_map
                if o in shadow_to_clade:
                    clade_dict = clade_ratio_dicts[shadow_to_clade[o]]
                    o = "Other " + shadow_to_clade[o]
                else:
                    clade_dict = clade_ratio_dicts[o]
                metadata_metadata_2.append(
                    TreeMetadata(o, clade_dict, 
                                 null_value=0, 
                                 color_method=lambda v: v/max_ratio)
                )

            draw_tree_metadata(ax_md2, y_posns,
                                           metadata_metadata_2,
                                           md_font_size=md_font_size,
                                   thickness=step,
                                   collapsed_clade_heights=kwargs.get(
                                       "collapsed_clade_heights",
                                       defaultdict(lambda: 2)
                                   ),
                                          )



            _ = ax_md2.set_ylim(*ax_tree.get_ylim())
            if not self.is_hier: 
                right_most_ax = ax_md2
            
        # label left side of left-most subplot with ref names
        if self.ref_names and not internal_ref_labels:
            add_label_refs(right_most_ax, y_posns, self.ref_names)
        

        return y_posns, step
    
def get_clade_lens(clade, seq_lens, filter_condition=lambda x: True):
    return numpy.array([
        seq_lens[s] 
        for s in (t.name for t in clade.get_terminals())
        if filter_condition(s)
    ])
    
def get_clade_ratio_dicts_dict(gene_hits, seqs=None, columns=None):
    columns = columns if columns else ['clade', 'clade_v', 'clade_3']
    
    seq_genes = \
        gene_hits[[q in seqs for q in gene_hits['query']]] \
        if seqs \
        else gene_hits

    clade_dicts_dict = {}
    for clade_col in columns:
        clade_counts = \
            seq_genes \
                [['query', 'domain', clade_col]] \
                .query('domain != "Viruses"') \
                .groupby(['query', clade_col]) \
                .agg(len)

        clade_counts_pivot = clade_counts \
            .reset_index() \
            .pivot('query', clade_col, 'domain') \
            .fillna(0)

        o_sum = clade_counts_pivot.sum(axis=1)
        clade_ratios = clade_counts_pivot.divide(o_sum, axis='rows')
        clade_dicts_dict[clade_col] = {o:clade_ratios[o].to_dict() for o in clade_ratios.columns}    
    
    return clade_dicts_dict


# draw a tree in polar space
# (modified from Bio.Phylo.draw())
def draw_polar(
    tree,
    label_func=str,
    do_show=True,
    show_confidence=True,
    # For power users
    axes=None,
    branch_labels=None,
    branch_label_x_delta=0,
    branch_label_y_delta=-.25,
    label_colors=None,
    *args,
    **kwargs
):
    """Plot the given tree using matplotlib (or pylab).

    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.

    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).

    Example using the pyplot options 'axhspan' and 'axvline'::

        from Bio import Phylo, AlignIO
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        constructor = DistanceTreeConstructor()
        aln = AlignIO.read(open('TreeConstruction/msa.phy'), 'phylip')
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        tree = constructor.upgma(dm)
        Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
        ... axvline={'x':0, 'ymin':0, 'ymax':1})

    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).

    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
        label_colors : dict or callable
            A function or a dictionary specifying the color of the tip label.
            If the tip label can't be found in the dict or label_colors is
            None, the label will be shown in black.

    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        try:
            import pylab as plt
        except ImportError:
            raise MissingPythonDependencyError(
                "Install matplotlib or pylab if you want to use draw."
            ) from None

    import matplotlib.collections as mpcollections

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)

    if not branch_labels:
        if show_confidence:

            def format_branch_label(clade):
                try:
                    confidences = clade.confidences
                    # phyloXML supports multiple confidences
                except AttributeError:
                    pass
                else:
                    return "/".join(conf2str(cnf.value) for cnf in confidences)
                if clade.confidence is not None:
                    return conf2str(clade.confidence)
                return None

        else:

            def format_branch_label(clade):
                return None

    elif isinstance(branch_labels, dict):

        def format_branch_label(clade):
            return branch_labels.get(clade)

    else:
        if not callable(branch_labels):
            raise TypeError(
                "branch_labels must be either a dict or a callable (function)"
            )
        format_branch_label = branch_labels

    # options for displaying label colors.
    if label_colors:
        if callable(label_colors):

            def get_label_color(clade):
                return label_colors(clade)

        else:
            # label_colors is presumed to be a dict
            def get_label_color(clade):
                return label_colors.get(label_func(clade), "black")

    else:

        def get_label_color(clade):
            # if label_colors is not specified, use black
            return "black"

    # Layout

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths


    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        # total height is number of leaves
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
        }

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights

    x_posns = get_x_positions(tree)
    y_posns  = get_y_positions(tree)
    terminals = []
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1, subplot_kw={'projection': 'polar'})
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for axes: %s" % axes)

    def draw_clade_lines(
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw=".1",
    ):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if orientation == "horizontal":
            horizontal_linecollections.append(
                mpcollections.LineCollection(
                    [[(y_here, x_start), (y_here, x_here)]], color=color, lw=lw
                )
            )
        elif orientation == "vertical":
            vertical_linecollections.append(
                mpcollections.LineCollection(
                    [[(y_bot, x_here), (y_top, x_here)]], color=color, lw=lw
                )
            )

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, "color") and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, "width") and clade.width is not None:
            lw = clade.width * plt.rcParams["lines.linewidth"]
        # Draw a horizontal line from start to here
        draw_clade_lines(
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            ## TODO: polar!
            axes.text(
                #0.5 * (x_start + x_here),
                x_here - branch_label_x_delta,
                y_here + branch_label_y_delta,
                conf_label,
                fontsize="small",
                #horizontalalignment="center",
                horizontalalignment="right",
                verticalalignment='bottom',
            )
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            ## TODO: polar!
            axes.text(
                x_here,
                y_here,
                " %s" % label,
                verticalalignment="center",
                color=get_label_color(clade),
            )
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]]
            y_bot = y_posns[clade.clades[-1]]
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(
                orientation="vertical",
                x_here=x_here,
                y_bot=y_bot,
                y_top=y_top,
                color=color,
                lw=lw,
            )
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)
        else:
            terminals.append(clade)


    draw_clade(tree.root, 0, "k", plt.rcParams["lines.linewidth"])

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)

    # Aesthetics

    try:
        name = tree.name
    except AttributeError:
        pass
    else:
        if name:
            axes.set_title(name)
    axes.set_xlabel("branch length")
    axes.set_ylabel("taxa")
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_ylim(-0.05 * xmax, 1.25 * xmax)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    axes.set_xlim(0.2, max(y_posns.values()) + 0.8)

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            list(value)
        except TypeError:
            raise ValueError(
                'Keyword argument "%s=%s" is not in the format '
                "pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),"
                " or pyplot_option_name=(dict) " % (key, value)
            ) from None
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif isinstance(value[0], tuple):
            getattr(plt, str(key))(*value[0], **dict(value[1]))

    if do_show:
        plt.show()
        
    return {t:y_posns[t] for t in terminals}


##########
# helpers added in the move out of experimental

def _gene_table(gene_table):
    """ polars gene tables get indexed by 'gene' (if present) """
    if gene_table is None or isinstance(gene_table, pandas.DataFrame):
        return gene_table
    index_col = 'gene' if 'gene' in gene_table.columns else None
    if index_col is None:
        return to_pandas(gene_table)
    return to_pandas_indexed(gene_table, index_col)


def read_tree(tree_file, fmt='newick', midpoint=True):
    """ read the first tree from a file and (by default) midpoint root it """
    tree = first(Phylo.parse(tree_file, fmt))
    if midpoint:
        tree.root_at_midpoint()
    return tree


# old name
TreePlotterVLP = SyntenyTreePlotter
