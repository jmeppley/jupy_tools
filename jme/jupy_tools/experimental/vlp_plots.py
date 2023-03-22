"""
This is a very specifc set of tools for the HOT319 tailed phage analysis plots. It probably won't hsve uch general usefulness without massive refactoring, but I'm using it too much not to put it into a module

== Synteny plots ==

class TreePlotterVLP()

    hides all the mess

    Initialize it with:
         * table of genome metadata
         * table of gene annotations
         * json dict maping reference genome ids to host taxon
         * json dict of assorted plot metadata (this is the REALLY messy one)

functions:

def draw_tree(tree_name, tree, collapse=True, cl_id=None, **kwargs):

Given a phylogenetic tree and genome metadata, make a synteny tree plot (using cluster_trees.master_plot)

== Single Gene Trees ==
Code for drawing circularized trees of genes

"""

import json
import numpy
import pandas

from collections import defaultdict

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable

from Bio import Entrez
Entrez.email = 'jmeppley@hawaii.edu'

from jme.jupy_tools.faves import *
from jme.jupy_tools.experimental import cluster_trees, gene_synteny

GM_TOP = 'top'
GM_TYPES = 'types'
GM_MERGE = 'merge'

class TreePlotterVLP():
    def __init__(
        self,
        md_tsv,
        gene_tsv,
        metadata_json,
        hosts_json=None,
        clusters_col='mcl_5_I2',
        annot_col='annot_e',
        desc_col='desc',
        length_col='Length',
        override_gene_descs=None,
    ):

        # load TSVs
        self.md_df = read_tsv(md_tsv)
        self.gene_df = read_tsv(gene_tsv)

        # apply custom column choices
        self.seq_lens = self.md_df[length_col].to_dict()
        if annot_col != "annot":
            if "annot" in self.gene_df.columns:
                self.gene_df = self.gene_df.drop('annot', axis=1)
            self.gene_df = self.gene_df.rename({annot_col: 'annot'}, axis=1)
        if desc_col != "desc":
            if "desc" in self.gene_df.columns:
                self.gene_df = self.gene_df.drop('desc', axis=1)
            self.gene_df = self.gene_df.rename({desc_col: 'desc'}, axis=1)

        self.clusters_col = clusters_col
        self.clusters = None

        # this should have keys: base_colors, type_patterns, parents, host_colors, source_colors
        # it can have known_hosts or (known_hosts_text and knwon_hosts_color)
        # or neither
        with open(metadata_json) as JSON:
            self.plot_metadata = json.load(JSON)
        missing_keys = {k for k
                        in ('base_colors',
                            'parents', 'host_colors', 'source_colors')
                        if k not in self.plot_metadata}

        # apply default types
        if 'type_patterns' not in self.plot_metadata:
            self.plot_metadata['type_patterns'] = [
                ['Terminase', 'terminase'],
                 ['Tail protein', 'tail'],
                 ['Helicase', 'helicase'],
                 ['Endonuclease', 'endonuclease'],
                 ['Primase', 'primase'],
                 ['DNA Polymerase', 'dna.polymerase'],
                 ['Polymerase', 'polymerase'],
                 ['Exonuclease', 'exonuclease'],
                 ['Integrase', 'integra(se|tion)'],
                 ['Pilus', 'pilus']
            ]

        # check the possible combinations of known_hosts keys
        if 'known_hosts_text' in self.plot_metadata.keys():
            if 'known_hosts_color' in self.plot_metadata.keys():
                # all good
                pass
            else:
                if 'known_hosts' in self.plot_metadata.keys():
                    self.plot_metadata['known_hosts_color'] = \
                            self.plot_metadata['known_hosts']
                else:
                    self.plot_metadata['known_hosts_color'] = \
                            self.plot_metadata['known_hosts_text']
        else:
            if 'known_hosts_color' in self.plot_metadata.keys():
                if 'known_hosts' in self.plot_metadata.keys():
                    self.plot_metadata['known_hosts_color'] = \
                            self.plot_metadata['known_hosts']
                else:
                    self.plot_metadata['known_hosts_color'] = \
                            self.plot_metadata['known_hosts_text']
            else:
                # neither is here, do we have known hosts?
                if 'known_hosts' in self.plot_metadata.keys():
                    self.plot_metadata['known_hosts_text'] = \
                            self.plot_metadata['known_hosts']
                    self.plot_metadata['known_hosts_color'] = \
                            self.plot_metadata['known_hosts']
                else:
                    if hosts_json is not None:
                        # if neither known hosts in metadata, load from separate JSON
                        with open(hosts_json) as JSON:
                            known_hosts = json.load(JSON)
                        self.plot_metadata['known_hosts_text'] = known_hosts
                        self.plot_metadata['known_hosts_color'] = known_hosts
                    else:
                        missing_keys.add('known_hosts')

        if len(missing_keys) > 0:
            raise Exception('ERROR: plot metdata is missing the following keys: ' +
                            str(missing_keys))

        # merge metadata with actual data
        self.plot_metadata.update(
            dict(
                seq_lens=self.seq_lens,
                gene_table=self.gene_df,
                taxonomy=None,
                ref_names=self.plot_metadata['known_hosts_text'],
                leaf_color_fn=self.get_leaf_color,
            )
        )


        # desc overrides can come from JSON or paameters
        # start with JSON
        self.override_gene_descs = \
            self.plot_metadata.get('override_gene_descs', {})
        # anything passed in params takes precedent
        if override_gene_descs:
            self.override_gene_descs.update(override_gene_descs)

        self.mds_cat = [
            ('Seq. Type', 
             ('Source_pretty' 
              if 'Source_pretty' in self.md_df.columns
              else 'Source'),
              None,
             'source_colors'
            ),
        ]
        #self.check_source_colors()
        self.mds_gene = []
        self.add_metadata('gene', 'Integrase',)
        
        self.mds_cont = [
            ('Abundance', self.md_df.LongReadHits),
            ('Element length', self.md_df[length_col], {'log':True}),
            ('GC', self.md_df.GC),
        ]

    def get_leaf_color(self, name):
        try:
            return self.plot_metadata['source_colors'][self.md_df.Source[name]]
        except:
            if re.search('novel', name):
                return self.plot_metadata['source_colors'].get('Novel',
                                                               'rebeccapurple')
            return 'black'


    def get_plot_data(self):
        return {k: self.plot_metadata[k]
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

    def get_collapsed_nodes(self, tree, max_d_cutoff=1.5):
        # lets collapse:
        #  - any clade over 10 seqs with no ref in it
        #  - any clade with more than 3 refs all the same host
        collapsed_nodes = {}
        cl_id_map = self.md_df[self.clusters_col]
        for clade in tree.get_nonterminals():

            # no nested collapsed clades
            for c in collapsed_nodes:
                if c.is_parent_of(clade):
                    skip = True
                    break
            else:
                skip = False
            if skip:
                continue

            terminals = {t.name for t in clade.get_terminals()}
            cl_id = cl_id_map.get(first(terminals), -1)

            if any(self.md_df.Source[t] == 'Refs' for t in terminals):
                if (
                    all(
                        (self.md_df.Source[t] == 'Refs'
                         and cl_id_map.get(t, -1) == cl_id
                        )
                        for t in terminals
                    )
                    and (len(terminals) > 2)
                ):
                    host_counts = get_host_counts(terminals,
                                                  self.plot_metadata['known_hosts_color'],
                                                  self.plot_metadata['parents'],
                                                 )
                    count_str = "\n".join(f"{c} {h}" for h,c in sorted(host_counts.items()))
                    label = f"Published phage:\n{count_str}"
                    collapsed_nodes[clade] = label
                continue

            if len(terminals) < 8:
                continue

            depths = clade.depths()
            max_d = max(depths[t] for t in clade.get_terminals())
            if max_d < max_d_cutoff and all(cl_id_map.get(t, -1) == cl_id for t in terminals):
                collapsed_nodes[clade] = f"{len(clade.get_terminals())} novel viral genomes"

        # don't collapse if set of clades is empty or a single clade that is the whole tree
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
        return get_type(gene_annot, self.plot_metadata['type_patterns'])

    def standardize_descs(self, gene_df, annot_col='annot', desc_col='desc'):
        """
        for each gene/desc pair in the table:
            * replace the desc if it's in self.override_gene_descs
            * replace the desc if it matches one of the patterns in 
                self.plot_metadata['type_patterns']
        """
        def get_new_desc(annot, desc):
            desc = self.override_gene_descs.get(annot, desc)
            gene_type = get_type(desc, self.plot_metadata['type_patterns'])
            if gene_type is None:
                return desc
            return gene_type

        gene_df[desc_col] = [
            get_new_desc(annot, desc)
            for annot, desc
            in gene_df[[annot_col, desc_col]].values
        ]

        return gene_df

    def get_gene_colors_and_descs(self, seqs, gene_method=GM_TOP):
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
            annot_descs.update(self.override_gene_descs)

        # count the annotations
        non_na_genome_counts = Counter(
            a
            for g,a in set(tuple(p) for p in cluster_genes[['genome', gene_col]].values)
            if pandas.notna(a)
        )

        base_colors = self.plot_metadata['base_colors']
        type_colors = dict(zip(
            (p[0] for p in self.plot_metadata['type_patterns']),
            base_colors,
        ))
        n_types = len(type_colors)

        if gene_method == GM_TYPES:
            # first, find the most common annotations with any of the primary types
            for gene, count in non_na_genome_counts.most_common():
                gene_annot = annot_descs[gene]
                gene_type = self.get_type(gene_annot)
                if gene_type is not None:
                    gene_colors[gene] = type_colors[gene_type]

                if len(gene_colors) >= 18:
                    break
            # prioritize non-type colors for the remaining genes
            color_iter = iter(base_colors[n_types:] + base_colors[:n_types])
        if gene_method == GM_MERGE:
            for (type_name, color) in type_colors.items():
                if type_name in non_na_genome_counts:
                    gene_colors[type_name] = color
            # prioritize non-type colors for the remaining genes
            color_iter = iter(base_colors[n_types:] + base_colors[:n_types])
        else:
            color_iter = iter(base_colors)

        # Add from top_genes until we reach 18
        tree_top_genes = self.get_top_genes(seqs)
        top_gene_iter = iter(dict(non_na_genome_counts \
                                    .most_common()).keys())
        while len(gene_colors) < 18:
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
                    gene_label = make_gene_desc(gene, gene, non_na_genome_counts[gene])
                    #print(gene, gene_label, gene_colors[geene])
                    gene_label_colors[gene_label] = gene_colors[gene]
                else:
                    gene_label = gene
                gene_labels.append(gene_label)
            cluster_genes['annot'] = gene_labels
            # return gene table instead of descriptions dict
            return gene_label_colors, cluster_genes

        # otherwise return dict of descrtiptions for labelling genes
        gene_descs = {
            gene: make_gene_desc(gene, annot_descs[gene], non_na_genome_counts[gene])
            for gene in gene_colors
        }
        return gene_colors, gene_descs

    def get_lens_by_source(self, seqs):
        lens_by_source = defaultdict(list)
        for seq in seqs:
            source = self.md_df.Source[seq]
            source = 'AFVG' if source.startswith('AFVG') else source
            lens_by_source[source].append(self.seq_lens[seq])
        return lens_by_source

    def filter_tree_seqs(self, tree_seqs, cutoff=10, multiplier=3, buffer=.05):
        """
        Return the set of sequences such that:
            If there are at least 10 trusted sequences (AFVG or Ref):
                All non truseted seqs are:
                    within 3sigma of the trusted length distribtion
                    or within {buffer} of the trusted length max and min
            Otherwse, return all seqs
        """
        lens_by_source = self.get_lens_by_source(tree_seqs)
        good_lens = lens_by_source['AFVG'] + lens_by_source['Refs']
        if len(good_lens) < cutoff:
            return tree_seqs

        three_sigma_range = get_sigma_range(good_lens, multiplier)
        length_range = [min(three_sigma_range[0], min(good_lens) * (1-buffer)),
                        max(three_sigma_range[1], max(good_lens) * (1+buffer))]

        return [seq
                for seq in tree_seqs
                if (self.md_df.Source[seq] != 'Opera'
                    or
                    length_range[0] <= self.seq_lens[seq] <= length_range[1]
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
        known_hosts = self.plot_metadata['known_hosts_text']
        for seq in seqs:
            if seq in known_hosts and known_hosts[seq] == "Unknown":
                tree.prune(seq)

    def prune_duplicates(self, tree, cutoff=1e-06):
        known_hosts = self.plot_metadata['known_hosts_text']
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

    def add_metadata(self, md_type, params):
        """
        Add custom metadata to the plots. There are three types:

        Categorical (md_type="cat"):
            pass a 4-tuple with:
                (label: string name for metadata
                 data_dict: dictionary mapping genome ids to cateories
                            or column name in md_df
                 null_value: what to pass to color dict if id not in val_dict
                 color_dict: dict from values to colors 
                            or key in plot_metadata
                )

        Gene Presence (md_type="gene"):
            pass 2- or 3- tuple with:
                (Label,
                 gene_filter: regex pattern (string o compiled)
                        or callable on (anot, desc) that returns True/False
                        (make sure callable can handle null genes and descs)
                 colors: dict mapping True/False to colors  (defaults to black/white)
                )
            or pass a string that matches a type_patterns name

        Continnuous (md_type="cont"):
            pass a 2- or 3-tuple with:
                (label: a string name for the metadata
                 data_dict: a dicttionaey from ids to values
                 kwargs: dict of extended args for get_continuous_md():
                     color_map: a matplotlib colormap name (optional)
                     log: use log scale if True
                )
        """
        if md_type == 'cat':
            self.mds_cat.append(params)
        elif md_type =='cont':
            self.mds_cont.append(params)
        elif md_type == 'gene':
            if len(params) == 1:
                params = params[0]
            if isinstance(params, str):
                for tname, tpattern in self.plot_metadata['type_patterns']:
                    if tname == params:
                        params = ('Has' + params, tpattern)
                        break
                else:
                    raise Exception("Cannot find a type with name: " +
                                    repr(params))
            self.mds_gene.append(params)

    def remove_metadata(self, label):
        for md_list in (self.mds_cat, self.mds_cont, self.mds_gene):
            for i, data in list(enumerate(reversed(md_list))):
                if data[0] == label:
                    md_list.pop(i)

    def md_dict_from_col(self, column):
        """ make a dict from a md column, filtering out the nans """
        return {
            k:v
            for k,v in self.md_df[column].items()
            if pandas.notna(v)
        }

    def get_metadata_metadata(self, seqs, cl_id=None):
        cluster_genes = gene_synteny.get_cluster_genes(self.gene_df,
                                                       seqs,
                                                       self.seq_lens)

        mdmd_data = []

        cl_str = f'Cluster {cl_id}'
        self.in_cluster_dict = {
            s:(cl_str if c == cl_id
               else (None if c < 0
                     else "Other")
              )
            for s in seqs
            for c in [self.md_df[self.clusters_col][s],]
        }
        if cl_id is not None and cl_id != -1:
            # if any are NOT in the original cluster, add MD column with cluster status
            if any(v != cl_str for v in self.in_cluster_dict.values()):
                mdmd_data.append(
                        (f'Cluster',
                         self.in_cluster_dict,
                         None,
                         {cl_str: 'black', "Other": 'lightgrey', None: 'white'},
                        )
                )

        # if there are any refs, add host column
        _hosts_dict = (
            self.md_df['Host']
            if 'Host' in self.md_df.columns
            else self.plot_metadata['known_hosts_color']
        )
        host_dict = {
            s:_hosts_dict[s]
            for s in seqs
            if pandas.notna(_hosts_dict.get(s, None))
        }
        if len(host_dict) > 0:
            mdmd_data.append(
                ('Host',
                 host_dict,
                 None,
                 'host_colors',
                )
            )

        # insert any user added MDs
        mdmd_data.extend(self.mds_cat)

        # if there are any interesting gene types
        # (Just integrase by deafult, but theuser can add others)
        default_gene_color_dict = {True: 'k'}
        for name, *params in self.mds_gene:
            type_filter = params[0]
            color_dict = default_gene_color_dict \
                            if len(params) == 1 \
                            else params[1]
            has_type = get_type_mask(type_filter, cluster_genes,
                                     self.override_gene_descs)
            if any(has_type.values()):
                mdmd_data.append(
                    (name, has_type, False, {True: 'k'}, )
                )

        # build metadata
        metadata_metadata = {
            label:cluster_trees.TreeMetadata(label,
                                             data_dict=(
                                                 self.md_dict_from_col(data)
                                                 if isinstance(data, str)
                                                 else data),
                                             null_value=null_val,
                                             color_dict=(
                                                 self.plot_metadata[c_dict]
                                                 if isinstance(c_dict, str)
                                                 else c_dict),
                                             color_method='cat')
            for label, data, null_val, c_dict in mdmd_data
        }

        # add continuous metadata
        for params in self.mds_cont:
            name, value_dict = params[:2]
            if len(params) > 2:
                # 3rd param should be kwargs dict
                gcmd_args = params[2]
                # (bakwards compat: if it's a string, it's the color map)
                if isinstance(gcmd_args, str):
                    gcmd_args = dict(color_map=gcmd_args)
            else:
                gcmd_args = {}
            metadata_metadata[name] = \
                get_continuous_md(name, value_dict, seqs, **gcmd_args)

        return metadata_metadata

    def draw_tree(self, tree_name, tree,
                  collapse=True, cl_id=None, draw_genes=None,
                  gene_method=GM_TOP,
                  **kwargs):
        """
        Given a phylogenetic tree and genome metadata, make a synteny tree plot (using cluster_trees.master_plot)

        draw_genes: is on by default with a sorting function to prioritize drawing genomes closest
                    to the mean length. Set to True to just draw the first few genomes in each clade
                    or a different sorting function. Set to false to just draw the tree and metadata
                    with no gene locations.

        collapse: if True (default) try to collapse clades of all novel or all ref sequences

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

        # if gene_method is merger, gene descs will be a modified gene_table
        # otherwise, it's a dict from gene to description
        gene_colors, gene_descs = self.get_gene_colors_and_descs(tree_seqs,
                                                                 gene_method)
        metadata_metadata = self.get_metadata_metadata(tree_seqs, cl_id)
        md_md = metadata_metadata.values()

        if collapse:
            if collapse is True:
                collapsed_nodes, collapsed_node_heights = self.get_collapsed_nodes(tree)
            else:
                collapsed_nodes, collapsed_node_heights = \
                    self.get_collapsed_nodes(tree, collapse)
        else:
            collapsed_nodes, collapsed_node_heights = None, None

        # provide defaults if not given through kwargs
        for k,v in dict(
            fig_width=20,
            items_per_inch=6,
            max_colored_genes=18,
            max_plotted_genes=18,
            add_gene_descs=(
                None if gene_method == GM_MERGE
                else gene_descs),
        ).items():
            kwargs.setdefault(k, v)

        plot_data = self.get_plot_data()
        if gene_method == GM_MERGE:
            plot_data['gene_table'] = gene_descs
        tree_plotter = cluster_trees.TreePlotter(tree, **plot_data)
        fig = tree_plotter.master_plot(tree_name,
                                       md_md,
                                       collapsed_clade_labels=collapsed_nodes,
                                       draw_genes=draw_genes,
                                       gene_colors=gene_colors,
                                      tree_kwargs=dict(
                                          collapsed_clade_heights=collapsed_node_heights),
                                       **kwargs,
                                      )

        # make some legends (this will fail for draw_genes=False)
        # figure out where the plots are
        x, y, w, h = fig.axes[0].figbox.bounds
        x2, y2, w2, h2 = fig.axes[2].figbox.bounds
        top = y2 + h2

        for md_params in self.mds_cont:
            md_key = md_params[0]
            md = metadata_metadata[md_key]
            cax = fig.add_axes([x/2, top, w * .4, h2 * .05], )
            _ = cax.set_title(md_key)
            _ = fig.colorbar(ScalarMappable(norm=md.color_method, cmap=md.color_map,),
                             cax=cax,
                             orientation='horizontal')
            top -= h2 * .25

        user_md_names = {m[0] for m in self.mds_cat}
        user_md_names.update([m[0] for m in self.mds_gene])
        lax = None
        lax_x = x/2
        lax_w = w * .4 / 3
        for name in sorted(metadata_metadata):
            # "Has xxx" or "Host"
            if (name in user_md_names 
                or name == 'Host' 
                or name.startswith('Cluster')
            ):
                lax = fig.add_axes([lax_x, top - h2 * .05, lax_w, h2 * .2], )
                for spine in lax.spines.values():
                    spine.set_visible(False)
                lax.set_yticks([])
                lax.set_xticks([])
                add_md_legend(tree_seqs,
                           metadata_metadata[name],
                           lax,
                           bbox_to_anchor=None,
                           loc='upper center')
                lax_x += lax_w

        return fig

def add_md_legend(seqs, metadata, fig, label=None, **legend_kwargs):
    if label is None:
        label = metadata.label

    patches = []
    names = []
    for value in set(metadata.data_dict.get(s, metadata.null_value)
                     for s in seqs):
        patches.append(mpatches.Patch(color=metadata.color_dict[value], ec='k', lw=.25))
        names.append(str(value))

    for k,v in dict(bbox_to_anchor=[.74, .81],
                    loc='upper left',
                    facecolor='white',
                    title=label).items():
        legend_kwargs.setdefault(k, v)
    return fig.legend(patches, names, **legend_kwargs)

def get_continuous_md(
    name,
    value_dict,
    seqs,
    color_map='binary',
    log=False,
    log_shift=1e-12,
):
    " construct Metadata object from a column of numeric data "
    if log:
        values = {s:value_dict[s] + log_shift
                  for s in seqs}
        vmin, vmax = min(values.values()), max(values.values())
        v_norm = colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)
        null_val = log_shift
    else:
        values = {s:value_dict[s] for s in seqs}
        vmin, vmax = min(values.values()), max(values.values())
        v_norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        null_val = None

    v_cmap = plt.get_cmap(color_map) if isinstance(color_map, str) else color_map

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

def get_host_counts(refs, known_hosts, parents):
    host_counts = Counter(re.sub('^Other ', '', known_hosts[r]) for r in refs)

    for host, count in list(host_counts.items()):
        if host in parents and parents[host] in host_counts:
            host_counts[parents[host]] += count
            del (host_counts[host])

    return host_counts

def make_gene_desc(gene, desc, count):
    gene_id_str = f' ({gene})' if re.search(gene, desc) is None else ''
    return f'{desc}{gene_id_str} (in {count} genomes)'

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
                   leaf_arc=None):
    
    # some things get set automatically based on the first (root) clade
    if max_distance is None:
        max_distance = max(clade.distance(t) for t in clade.get_terminals())
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
        distances = numpy.array([clade.distance(t) for t in clade.get_terminals()])
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
            d = clade.distance(child)    # radial distance
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

