"""
This is a second version of vlp_plots, for use with the phage depth profile
analyses. Hopefully it will prove to be a bit more future proof than the
original.

== Synteny plots ==

class TreePlotterVLP()

    hides all the mess

    Initialize it with:
         * table of genome metadata
         * table of gene annotations
         * json dict of metadata colors

functions:

def draw_tree(tree_name, tree, collapse=True, cl_id=None, **kwargs):

Given a phylogenetic tree and genome metadata, make a synteny tree plot (using cluster_trees.master_plot)

=== Changes from prev version ===

All metadata should be in the metadata table. Only colors are in the
plot_metadata json:

The main arguments can either be paths to a file or the loaded data object.

The various ways to group hosts should be different columns in the metadat
table. There are now methods provided to translate host names into broader
categories if you really want to do it on the fly.

This is also true for the type_patters code. We no wexpect to apply groupings
before calling the plot drawing code.

== Single Gene Trees ==
Code for drawing circularized trees of genes

"""

import json
import numpy
import pandas

from pandas.api.types import is_numeric_dtype

from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable

from Bio import Entrez
from Bio.Phylo import BaseTree
Entrez.email = 'jmeppley@hawaii.edu'

from jme.jupy_tools.faves import *
from jme.jupy_tools.utils import get_N_colors
from jme.jupy_tools.experimental import cluster_trees, gene_synteny

GM_TOP = 'top'
GM_TYPES = 'types'
GM_MERGE = 'merge'

class TreePlotterVLP():
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
    ):
        self.debug = False

        # load TSVs
        #  (if filenames are not passed, assume its a DataFrame)
        self.md_df = read_tsv(md_tsv) if isinstance(md_tsv, str) else md_tsv
        self.gene_df = read_tsv(gene_tsv) if isinstance(gene_tsv, str) \
                    else gene_tsv

        # we will look in here for metadata colors
        if isinstance(metadata_json, str):
            with open(metadata_json) as JSON:
                self.plot_metadata = json.load(JSON)
        else:
            if metadata_json is not None:
                self.plot_metadata = metadata_json
            else:
                # put a blank dict there
                self.plot_metadata = {}

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
        self.type_patterns = type_patterns if type_patterns else {}
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
                            dict(column='length_col', log=True),
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
            self.leaf_color_md = \
                self.get_tree_metadata(**self.plot_metadata['genomes']['leaf_color'])
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
                self.leaf_color_md = self.metadata_metadata_all[leaf_md_index]
            else:
                # if we find nothing, let everything be written in black
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
        self.debyg = False

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
        genes = gene_synteny.get_cluster_genes(self.gene_df,
                                               seqs,
                                               self.seq_lens)
        top_genes, _ = \
            gene_synteny.plot_annot_positions(genes, ax=False,
                                              max_plotted_genes=N,
                                              max_colored_genes=N,
                                             )
        return top_genes

    def get_collapsed_nodes(
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
            return None, None

        # set_clade_heights
        biggest_clade = max(len(c.get_terminals()) for c in collapsed_nodes)
        collapsed_node_heights = {
            c:numpy.interp(len(c.get_terminals()), [2, biggest_clade], [2,5])
            for c in collapsed_nodes
        }

        return collapsed_nodes, collapsed_node_heights

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
        cluster_genes = gene_synteny.get_cluster_genes(self.gene_df,
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
                    metadata_list.append(cluster_trees.TreeMetadata(
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
                if "collapsed_clade_heights" in kwargs:
                    collapsed_node_heights = kwargs.pop("collapsed_clade_heights")
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

        tree_plotter = cluster_trees.TreePlotter(tree, **plot_data)
        fig = tree_plotter.master_plot(tree_name,
                                       md_md,
                                       collapsed_clade_labels=collapsed_nodes,
                                       draw_genes=draw_genes,
                                      tree_kwargs=dict(
                                          collapsed_clade_heights=collapsed_node_heights),
                                       **kwargs,
                                      )

        # make some legends (this will fail for draw_genes=False)
        if draw_genes == False:
            print("WARNING: I don't know how to make legnds when draw_genes is False!!")
            return fig
        
        # figure out where the plots are
        x, y, w, h = fig.axes[0].figbox.bounds
        x1, y1, w1, h1 = fig.axes[1].figbox.bounds
        x2, y2, w2, h2 = (
            fig.axes[2].figbox.bounds
            if 'clade_ratio_dicts' not in kwargs
            else
            fig.axes[4].figbox.bounds
        )
        top = y2 + h2

        # separate categorical and continuous MDs
        cat_mds = []
        cont_mds = []

        for md in md_md:
            if md is None:
                continue
            label = md.label
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
        ring_width=.1,
    ):

        # get metadata
        tree_seqs = get_tree_seqs(tree)
        md_md = self.get_metadata_metadata(tree_seqs)

        fig, ax = plt.subplots(1,1,figsize=figsize)

        ring_order = get_ring_order(tree.root)
        leaf_arc = 355 / len(ring_order)
        r = recursive_draw(ax, tree.root, leaf_arc=leaf_arc, 
                           bg_styles_dict=bg_styles_dict)

        
        for mdmd in md_md:
            if mdmd is not None:
                draw_md_ring(ax, ring_order, mdmd, radius=r, width=ring_width, leaf_arc=leaf_arc)
                r += ring_width
            else:
                r += ring_width / 2

        _ = ax.set_xlim(-r,r)
        _ = ax.set_ylim(-r,r)
        _ = ax.set_xticks([])
        _ = ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        _ = ax.set_title(tree_name, pad=30)
        
        return fig

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
        md_type = 'cont' if is_numeric_dtype(data_series) else 'cat'
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
    return cluster_trees.TreeMetadata(
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
):
    name = data_series.name if label is None else label
    " construct Metadata object from a column of numeric data "
    if log:
        values = {s:data_series[s] + log_shift
                  for s in seqs}
        vmin, vmax = min(values.values()), max(values.values())
        v_norm = colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)
        null_val = log_shift
    else:
        values = {s:data_series[s] for s in seqs}
        vmin, vmax = min(values.values()), max(values.values())
        v_norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        null_val = None

    v_cmap = plt.get_cmap(color_map) \
        if isinstance(color_map, str) \
        else color_map

    return cluster_trees.TreeMetadata(
        name,
        values,
        null_value=null_val,
        null_color=[1.,1.,1.,0.,],
        color_map=v_cmap,
        color_method=v_norm,
    )

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

## cicular gene trees
from matplotlib import patches, collections

def plot_polar_in_cart(ax, r, th, offset=[0,0], from_degrees=False, **kwargs):
    """ convert r,theta to x,y and plot on cartesian axis """
    x, y = polar_to_cart(r, th, from_degrees)
    ax.plot(x + offset[0],
            y + offset[1],
            **kwargs)

def polar_to_cart(r, th, from_degrees=False):
    if from_degrees:
        th = numpy.radians(th)
        return  (r * numpy.cos(th),
                 r * numpy.sin(th))
    return  (r * numpy.cos(th),
             r * numpy.sin(th))

def get_clade_style_by_type(clade, style_dict, get_leaf_type, type_styles):
    """
    Depth-first crawl of a tree to find clades of all the same type
    PArams:
        node: Bio.Phylo clade
        style_dict: empty dict to save styles to
        get_leaf_type: callable on a leaf node, returns a node type
        type_styles: dict-like mapping types to matplotlib patch kwargs dict
    """
    node_types = set()

    # jest get the type of a leaf node
    if clade.is_terminal():
        return [get_leaf_type(clade.name),]

    # get types of all children
    children_with_one_type = {}
    for child in clade.clades:
        child_types = get_clade_style_by_type(child, style_dict, get_leaf_type, type_styles)
        node_types.update(child_types)
        if len(child_types) == 1:
            children_with_one_type[child] = next(iter(child_types))

    # if this node has multiple types, set the the style for any children with one type
    if len(node_types) > 1:
        for child, child_type in children_with_one_type.items():
            style_dict[child] = type_styles[child_type]

    return node_types

def get_ring_order(clade):
    if clade.is_terminal():
        return [clade,]
    return list(chain(*(get_ring_order(child)
                        for child in clade.clades)))

def add_legend(fig, colors_by_label, size=15, **legend_kwargs):
    patches = []
    names = []
    for label, color in colors_by_label:
        patches.append(collections.CircleCollection([size,], facecolors=[color,]))
        names.append(f"{label}")

    for k,v in dict(loc='upper left',
                    facecolor='white',
                    title='Colors').items():
        legend_kwargs.setdefault(k, v)
    return fig.legend(patches, names,
                      **legend_kwargs)

from matplotlib.colors import to_rgb, to_rgba

def get_collapsed_clades_source_host(
    tree,
    get_leaf_type,
    known_hosts,
    source_colors,
):
    """
    returns collapsed clades colored by source, but split by host
    """
    collapsed_leaves = set()
    collapsed_clades = {}
    for clade in tree.get_nonterminals():
        consensus = None
        for leaf in clade.get_terminals():
            if leaf in collapsed_leaves:
                # skip if parent clade is collapsed
                break

            source = get_leaf_type(leaf.name)
            host = known_hosts.get(leaf.name, None)

            # only collapse if same source and host (or host is None)
            if consensus is None:
                consensus = (source, host)
            elif consensus[1] is None:
                if source == consensus[0]:
                    consensus = (source, host)
                else:
                    break
            elif host is None:
                if consensus[0] != source:
                    break
            elif consensus != (source, host):
                break
        else:
            # if we get here, ccollapse the clade
            color = source_colors[consensus[0]]
            collapsed_clades[clade] = dict(
                width=clade.count_terminals(),
                ec=color,
                fc=to_rgba(color, .5),
                lw=.5,
            )
            collapsed_leaves.update(clade.get_terminals())
    return collapsed_clades

default_line_style = dict(color='k', lw=.5)
def recursive_draw(ax, clade,
                   distance=0,
                   ax_patches=None,
                   last_distance=None,
                   center=180,
                   total_arc_span=340,
                   min_collapsed_size=2,
                   line_style=None,
                   style_dicts=None,
                   line_styles_dict=None,
                   bg_styles_dict=None,
                   collapsed_dict=None,
                   offset=(0,0),
                   max_distance=None,
                   unit_branch_lengths=None,
                   leaf_arc=None):

    # some things get set automatically based on the first (root) clade
    if max_distance is None:
        max_distance = max(clade.depths().values())
        if max_distance == 0:
            max_distance = max(clade.depths(unit_branch_lengths=True).values())
            unit_branch_lengths = True
        else:
            unit_branch_lengths = False
    if style_dicts is None:
        style_dicts = defaultdict(dict)
        if line_styles_dict is not None:
            style_dicts['lines'] = line_styles_dict
        if bg_styles_dict is not None:
            style_dicts['bg'] = bg_styles_dict
        if collapsed_dict is not None:
            style_dicts['collapsed'] = collapsed_dict
    if leaf_arc is None:
        # how many arc segments total
        total_segs = clade.count_terminals()

        # for collapsed nodes, replace leaf count with collapsed size
        skip_nodes = set()
        for cl, cdata in style_dicts['collapsed'].items():
            try:
                c_size = float(cdata)
            except (ValueError, TypeError):
                c_size = cdata.get('width', min_collapsed_size)
            # don't consider clades inside other collapsed clades
            if not any(t in skip_nodes for t in cl.get_terminals()):
                c_size = max(int(c_size), min_collapsed_size)
                total_segs += c_size - cl.count_terminals()
                skip_nodes.update(cl.get_terminals())

        # get fraction of requested span
        leaf_arc = total_arc_span / total_segs

    if unit_branch_lengths:
        def get_clade_distance(node):
            return 1
    else:
        get_clade_distance = clade.distance

    if ax_patches is None:
        add_patches_on_exit = True
        ax_patches = defaultdict(list)
    else:
        add_patches_on_exit = False

    # get the line style. Use clade specific if available, fall back to parent style, then to default
    if clade in style_dicts['lines']:
        # keep default settings for anything not overriden by the clade style
        line_style = dict(default_line_style)
        line_style.update(style_dicts['lines'][clade])
    elif line_style is None:
        line_style = default_line_style

    # draw radial line from parent arc
    if last_distance is not None:
        tot_distance = last_distance + distance
        plot_polar_in_cart(ax,
                           (last_distance, tot_distance),
                           (center, center),
                           offset=offset,
                           from_degrees=True,
                           **line_style)
    else:
        tot_distance = distance

    # figure out where this clade is
    N = clade.count_terminals()
    start = center - leaf_arc * (N/2)
    end = center + leaf_arc * (N/2)

    # if clade has a background set
    bg_style = style_dicts['bg'].get(clade, None)
    if bg_style is not None:
        ax_patches['wedges'].append(patches.Wedge(offset,
                                   r=tot_distance,
                                   theta1=start,
                                   theta2=end,
                                   width=tot_distance - max_distance,
                                   **bg_style,
                                  )
                    )

    # is this clade collapsed?
    collapsed = style_dicts['collapsed'].get(clade, None)

    # Yes if there is a nonnull or False value. Empty dict is OK (but would resolve to False)
    if collapsed is not None and collapsed is not False:
        # get default style
        # use bg settings with line style overrides
        triangle_style = dict(line_style)
        if bg_style is not None:
            triangle_style.update(bg_style)

        # collapsed value could be one of: True, an integer, or a style
        if isinstance(collapsed, bool):
            triangle_width = min_collapsed_size
        else:
            try:
                triangle_width = max(float(collapsed), min_collapsed_size)
            except (ValueError, TypeError):
                # assume it is a dict
                collapsed_style = dict(collapsed)
                # pull width from dict (if it's there)
                try:
                    triangle_width = max(collapsed_style.pop('width'), min_collapsed_size)
                except KeyError:
                    triangle_width = min_collapsed_size
                # use rest of dict to update style
                triangle_style.update(collapsed_style)

        # get triangle shape
        distances = numpy.array([get_clade_distance(t) for t in clade.get_terminals()])
        # set up triangle in polar corrdinates
        triangle_verts_polar_unzipped = [
            (tot_distance, tot_distance + distances.max(), tot_distance + distances.min()),
            (center, start, end)
        ]
        # convert to cartesian pairs
        triangle_verts = list(zip(*polar_to_cart(*triangle_verts_polar_unzipped,
                                                 from_degrees=True)))


        ax.add_patch(patches.Polygon(triangle_verts, **triangle_style))

        # Alternative approach to use PolyCollection to speed things up
        #ax_patches['collapsed'].append(dict(
        #    verts=triangle_verts,
        #    **triangle_style
        #))
        depths = distances + tot_distance
    else:
        # figure out where children go and draw them
        last_end = start
        centers = []
        depths = [tot_distance,]
        for child in clade.clades:
            d = get_clade_distance(child)    # radial distance
            n = child.count_terminals()  # how many leaves
            w2 = n * leaf_arc / 2        # 1/2 angle covered
            c = last_end + w2            # center (half from previous end)
            centers.append(c)
            d = recursive_draw(ax, child,
                               ax_patches=ax_patches,
                               distance=d,
                               last_distance=tot_distance,
                               center=c,
                               offset=offset,
                               line_style=line_style,
                               style_dicts=style_dicts,
                               leaf_arc=leaf_arc,
                               max_distance=max_distance,
                               unit_branch_lengths=unit_branch_lengths,
                              )
            depths.append(d)
            last_end = c + w2

        # connect children with arc
        if centers and tot_distance > 0:
            diameter = tot_distance * 2
            arc = \
                patches.Arc(
                    offset,
                    width=diameter, height=diameter,
                    theta1=min(centers),
                    theta2=max(centers),
                    **line_style
                )
            #ax_patches['arcs'].append(arc)
            ax.add_patch(arc)


    if add_patches_on_exit:
        ax.add_collection(collections.PatchCollection(ax_patches['wedges'], match_original=True))
        #ax.add_collection(collections.PatchCollection(ax_patches['arcs']))

        """
        # collapsed triangles take a little massaging
        if 'collapsed' in ax_patches:
            triangle_data = defaultdict(list)
            all_keys = set(chain(*(t.keys() for t in ax_patches['collapsed'])))
            for triangle in ax_patches['collapsed']:
                for key in all_keys:
                    if key not in triangle:
                        raise Exception(f'All collapsed clade styles must have the same keys. Key {key} is missing.')
                    triangle_data[TRIANGLE_KEY_MAP[key]].append(triangle[key])
            triangle_data['closed'] = True
            ax.add_collection(collections.PolyCollection(**triangle_data))
        """
    return max(depths)

def get_min_max(values):
    mn, mx = None, None
    for value in values:
        if mn is None and mx is None:
            mn, mx = value, value
        else:
            mn = min(mn, value)
            mx = max(mx, value)
    return mn, mx

# Set up collapsed clades

def draw_md_ring(ax, ring_order, mdmd, radius=3, width=.1, center=180, leaf_arc=None, offset=[0,0]):
    N = len(ring_order)
    if leaf_arc is None:
        leaf_arc = 340 / N

    ax_patches = []
    start = center - leaf_arc * (N/2)
    for leaf in ring_order:
        color = mdmd.get_item_color(leaf.name)
        end = start + leaf_arc
        ax_patches.append(patches.Wedge(offset,
                                   r=radius,
                                   theta1=start,
                                   theta2=end,
                                   width=-width,
                                   fc=color,
                                   ec=None,
                                  )
            )
        start = end

    ax.add_collection(collections.PatchCollection(ax_patches, match_original=True))

