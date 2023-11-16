"""

Collection of tools for making my multi-panel tree and gene-synteny plots.

Includes methods for making trees, both phylogenitic and distance based dendrograms.

def cat_tree_dir_alignment(aln_files_by_vog, **kwargs):
    concatenates multiple single-gene protein alignemnts into one big alignment


def filter_cat_alignment(cat_alignment, genome_residue_frac=.05, residue_genome_frac=.5):
    Filters genomes and residues from an MSA:
        First, drops residues with < 50% of genomes
        Srcond, drops genomes with < 5% of residues

def get_genes_from_genomes(gene_table, genomes, column='genome'):
    Filters given pandas DataFrame to only include listed items in a colomn

class TreeBuilder()
    class for picking shared genes for tree building from a set of genomes
    def get_tree_genes()
        pick genes for a set of genomes
    def expand_tree_genomes()
        look for other genomes with these genes to include in a tree
    def pull_vogs()
        get the protein sequences from the selected genes and genomes into faa files

class TreeMetadata()
    holds metadata for right-hand-side metadata in TreePlotter plots

def draw()
    modified version of Bio.Phylo tree drawing. only used internally

def generate_collapsed_nodes()
    Automatically select nodes to collapse (to fit a big tree into one plot)

def cluster_synteny_plot()
    draws the two-panel synenty plots at the core of my bigger plots
        top: locations for each gene
        bottom: gene map, one geneome per line

def get_tree_plotter_data()
    convenience function for reading much of the data for the tree plots

class TreePlotter()
    makes the big plots with gene synteny tied to a genome tree and metadata bars
    def master_plot()

def draw_polar()
    a modified Bio.Phylo.draw() that can make a circular tree with collapsed branches

"""

import re, numpy, pandas, os
import attr
from Bio import Align, Seq, SeqRecord, AlignIO, SeqIO, Phylo
from collections import defaultdict, Counter
from functools import lru_cache
from itertools import combinations, islice, chain
from jme.jupy_tools.utils import read_tsv
from jme.jupy_tools.experimental import gene_synteny
from edl.taxon import readTaxonomy
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.colors import to_rgba
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree

def cat_tree_dir_alignment(aln_files_by_vog, format='stockholm'):
    """
    Given: dict from gene name to stockholm aln file,
    
    Returns a Bio.MultipleSequenceAlignment built from multiple MSA files
    """
    alignments = {}
    gene_dict = defaultdict(dict)
    total_length = 0

    # parse alignments
    for vog, align_file in aln_files_by_vog.items():
        aln_format = format if isinstance(format, str) else format[vog]
        vog_alignment = AlignIO.read(align_file, format=aln_format)

        # keep a record of where genomes aligned
        for i, gene in enumerate(vog_alignment):
            gene_id = gene.id
            genome = re.sub(r"_\d+$", "", gene_id)
            gene_dict[genome].setdefault(vog, []).append((i, gene_id))

        alignments[vog] = vog_alignment
        total_length += len(vog_alignment)


    # concat the alignments
    records = []
    vog_spacers = {v:"-"*a.get_alignment_length() for v,a in alignments.items()}
    for genome in gene_dict.keys():
        genome_aligned_sequence = ""
        for vog, vog_alignment in alignments.items():
            # get alignment(s) for the current genome from this vog alignment
            genome_vog_alignments = gene_dict[genome].get(vog, [])

            if len(genome_vog_alignments) == 0:
                # no alignment, fill with all gaps
                genome_aligned_sequence += vog_spacers[vog]
            elif len(genome_vog_alignments) == 1:
                # get the  only alignment
                a_index, gene_id = genome_vog_alignments.pop()
                genome_aligned_sequence += str(vog_alignment[a_index].seq)
            else:
                # merge the alignments 
                def pick_aa(res_chars):
                    " Take the first non-gap in the set of chars "
                    for c in res_chars:
                        if c not in {'.','-'}:
                            return c
                    return res_chars[0]
                res_char_generator = (pick_aa(res_chars)
                                      for res_chars 
                                      in zip(*[str(vog_alignment[ai].seq) 
                                               for ai, gi in genome_vog_alignments])
                                     )
                genome_aligned_sequence += "".join(res_char_generator)

        # make into a record and save
        records.append(SeqRecord.SeqRecord(Seq.Seq(genome_aligned_sequence), id=genome))
    # stack records into alignment object
    final_alignment = Align.MultipleSeqAlignment(records, )

    return final_alignment

def filter_cat_alignment(cat_alignment, genome_residue_frac=.05, residue_genome_frac=.5):
    """
    Filters genomes and residues from an MSA:
    
    First, drops residues with < 50% of genomes
    
    Srcond, drops genomes with < 5% of residues
    
    Returns a new Bio.MultipleSequenceAlignment object
    """
    if residue_genome_frac > 0:
        # count gapped genomes at each residue
        gapped_genome_counts = defaultdict(int)
        for j, alignment in enumerate(cat_alignment):
            for i, c in enumerate(alignment.seq):
                if c in {'.', '-'}:
                    gapped_genome_counts[i] += 1  # add 1 for this residue

        # flag residues with <50% of genomes
        cutoff = len(cat_alignment) * residue_genome_frac
        residues_to_drop = {i
                            for i, gap_count in gapped_genome_counts.items()
                            if (len(cat_alignment) - gap_count) < cutoff}
    else:
        residues_to_drop = {}

    # build new alignment dropping flagged residues and skipping genomes with too many gaps
    alignments = []

    for j, alignment in enumerate(cat_alignment):
        new_seq = ""
        align_gap_count = 0
        for i,c in enumerate(str(alignment.seq)):
            if i not in residues_to_drop:
                if c in {'.', '-'}:
                    align_gap_count += 1
                new_seq += c
        if len(new_seq) == 0 or ((1 - genome_residue_frac) < (align_gap_count / len(new_seq))):
            # skip if genomes has < 5% of residues
            continue

        alignments.append(SeqRecord.SeqRecord(Seq.Seq(new_seq), 
                                              id=alignment.id))
    return Align.MultipleSeqAlignment(alignments,)

def get_genes_from_genomes(gene_table, genomes, column='genome'):
    """ subset table off the genome column """
    genomes = set(genomes)
    return gene_table[[g in genomes 
                       for g in gene_table.genome.values]]

class TreeBuilder():
    " loads and holds dat for a clustering "
    def __init__(self, gene_table):
        self.gene_table = gene_table
    
    @lru_cache(maxsize=None)
    def get_cl_genes(self, seqs):
        return get_genes_from_genomes(self.gene_table, 
                                                   seqs)
    
    @lru_cache(maxsize=None)
    def get_cl_50_genes(self, seqs, min_hmm_frac=.5, min_hmm_len=100):
        # return genes in these seqs that are over 50% of the HMM and over 100AA
        return self.get_cl_genes(seqs) \
                   .query(f'hmm_frac >= {min_hmm_frac} and hmm_start + {min_hmm_len} <= hmm_end')

    @lru_cache(maxsize=None)
    def get_seq_vog_counts(self, seqs, cutoff=.1, max_copy_num=1.2, **cl_genes_kwargs):
        min_count = len(seqs) * cutoff
        cl_genes_50 = self.get_cl_50_genes(seqs, **cl_genes_kwargs)
        copy_num = cl_genes_50.groupby(['genome', 'vog']).agg({'gene_no': get_copy_no})
        mean_copy_num = copy_num.reset_index().groupby('vog').agg({'gene_no': numpy.mean})
        sc_genes = set(mean_copy_num.query(f'gene_no <= {max_copy_num}').index)

        return copy_num \
            .reset_index()[[vog in sc_genes for genome, vog in copy_num.index]][['genome','vog']] \
            .groupby('vog').agg(len) \
            .query(f'genome >= {min_count}') \
            .sort_values('genome', ascending=False)

    def get_tree_genes(self,
                       seqs,
                       n_genes=None,
                       gene_seq_frac=.1,
                       seq_gene_frac=.05,
                       max_search_combos=5000,
                       required_genomes=None,
                       **kwargs
                      ):
        """
        Pick a set of genes to build a tree from with the following constraints:
        
        seqs: collection of genome names
        n_genes: choose this many genes (default: None -> as many as possible)
        gene_seq_frac: chosen genes must be present in this fraction of the genomes
        seq_gene_frac: genomes must have at least this fraction of the genes
        required_genomes: one of:
            :None: no constraints
            :Collection/Iterable: given genomes must be included
            :Callable: chose genes to maximize the return value of this function
        
        other params:
            max_search_combos:
                when chosing between multiple possible gene sets, try this many combinations
            kwargs: params for get_cl_50_genes
        """
        # make sure seq set is hashable and consistent for caching
        seqs = hashable(seqs)
        n_seqs = len(seqs)
        cl_genes_50 = self.get_cl_50_genes(seqs, **kwargs)
        genome_vog_counts = self.get_seq_vog_counts(seqs, cutoff=gene_seq_frac)

        if n_genes is None or len(genome_vog_counts) <= n_genes:
            n_genes = len(genome_vog_counts)
            

        # get the biggest set of seqs using these genes
        best_score = 0
        all_scores = []
        best_gene_sets = [(n_genes, 0, [], [])]
        
        # if no max given, do them all
        if max_search_combos is None or max_search_combos <= 0:
            max_search_combos = len(n_genes) * len(n_genes)
            
        # set up genome set check
        if required_genomes:
            if not callable(required_genomes):
                required_genome_set = set(required_genomes)
                n_required_genomes = len(required_genome_set)
                # return # genomes or 0 if not all of the reuired genomes are present
                score_genomes = lambda genomes: len(genomes) if len(required_genome_set.intersection(genomes)) == n_required_genomes else 0
            else:
                # it's already a callable, assume (or hope) the user knows what they're doing
                score_genomes = required_genomes
        else:
            # choose the biggest set of genomes
            score_genomes = lambda genomes: len(genomes)
                

        for gene_set in islice(combinations(genome_vog_counts.index, 
                                            n_genes), 
                               max_search_combos):
            gene_set = set(gene_set)
            # number of genes per genome
            genomes_hit = Counter(cl_genes_50[[(v in gene_set) 
                                               for v in cl_genes_50.vog]].genome)
            # total genes in all genomes
            score = score_genomes([g 
                                   for g,c in genomes_hit.items() 
                                   if c >= n_genes * seq_gene_frac])
            all_scores.append(score)
            if score == 0 or score < best_score:
                continue
                
            key = tuple(sorted(gene_set))
            data_tuple = (n_genes,
                          score, 
                          gene_set, 
                          genomes_hit)
            if score > best_score:
                best_gene_sets = [data_tuple,]
            else:
                best_gene_sets.append(data_tuple)
                
        tree_data = dict(zip(['n_genes','score', 'genes', 'genomes'],
                             next(iter(reversed(sorted(best_gene_sets))))))
        tree_data['all_scores'] = all_scores
        return tree_data
    
    @lru_cache(maxsize=None)
    def get_vog_hmm_lens(self):
        return dict(tuple(a) for a in self.gene_table[['vog','hmm_len']].values)
        
        
    @lru_cache(maxsize=None)
    def get_genome_vog_mlens(self):
        """ we can calculate this once for all genome/vog pairs and re use them in combination each time """
        return self.gene_table.eval('mlen = hmm_end + 1 - hmm_start').groupby(['genome','vog']).agg({'mlen':sum}).reset_index()
    
    def get_genomes_from_genes(self, genes, cutoff=0):
        """
        Get the toal hmm alignment length to the given genes for all genomes (with a nonzero value)
        """
        genome_vog_mlens = self.get_genome_vog_mlens()
        gene_set = set(genes)
        genome_mlen_totals = genome_vog_mlens[[(v in gene_set) for v in genome_vog_mlens.vog]].groupby('genome').agg({'mlen':sum})
        if cutoff <= 0:
            return genome_mlen_totals
        
        vog_lens = self.get_vog_hmm_lens()
        genes_len_total = sum(vog_lens[v] for v in gene_set)
        return genome_mlen_totals.query(f"mlen >= {cutoff * genes_len_total}")


    def expand_tree_genomes(self, tbd, min_gene_frac=1/3, max_new_genome_frac=1/4, new_tbd=True, **kwargs):
        """
        add new genomes to a given tree
        
         * tbd: tree building data retuurned from get_tree_genes
         * min_gene_frac: new genomes need at least this frac of the gene set in tbd
         * max_new_genome_frac: don't add more than this frac of the genomes already in tbd
         * new_tbd: if True return new tbd data structure, otherwise just new seq list
         
        returns tbd dict or new seq list (if new_tbd==False)
        """
        other_genomes_gene_counts = \
            self.gene_table[[(v in tbd['genes'] 
                              and g not in tbd['genomes']) 
                             for g, v in self.gene_table[['genome','vog']].values]] \
                .groupby(['genome','vog']).agg(len) \
                .reset_index()[['genome','vog']] \
                .groupby('genome').agg(len)

        cutoff = min_gene_frac * len(tbd['genes'])
        max_ext_genomes = int(numpy.floor(len(tbd['genomes']) * max_new_genome_frac))
        new_seqs = list(other_genomes_gene_counts \
                            .query(f'vog >= {cutoff}') \
                            .vog 
                            .sort_values() \
                            .tail(max_ext_genomes) \
                            .index)
        if new_tbd:
            seqs = tuple(sorted(list(tbd['genomes']) + new_seqs))
            return self.get_tree_genes(seqs, **kwargs)
        else:
            return new_seqs

    def pull_vogs(self, tree_dir, tree_genomes, gene_set, faa_file='genomes.prodigal.faa'):
        # get list of genes to pull out (as dict from gene_id to vog)
        gene_set = set(gene_set)
        tree_genes = \
            self.gene_table[[((vog in gene_set) and (genome in tree_genomes))
                            for genome, vog 
                            in self.gene_table[['genome', 'vog']].values]] \
                .vog.to_dict()

        # load genes
        vog_genes = defaultdict(list)
        for gene in SeqIO.parse(faa_file, 'fasta'):
            if gene.id not in tree_genes:
                continue
            vog = tree_genes[gene.id]

            vog_genes[vog].append(gene)

        for vog, genes in vog_genes.items():
            with open(f'{tree_dir}/{vog}.faa', 'wt') as vog_faa:
                SeqIO.write(genes, vog_faa, 'fasta')

def get_copy_no(gene_nos):
    copy = 1
    for i, num in enumerate(gene_nos[1:], start=1):
        if gene_nos[i - 1] + 1 < num:
            copy += 1
    return copy

def hashable(collection):
    return tuple(sorted(set(collection)))

## Plotting functions
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
            either "category" or "cat" (default) 
            (chooses evenly spaced colors from color_map for each unique value in data_dict)
            or "continuous" or "con"
            (data_dict values msut be in range [0.-1.] or [-255])
            or a callable transformation to convert values to a cmap-suitable range (implyies continuous)
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
    null_color = attr.ib(default=(1., 1., 1., 0.))
    collapsed_clade_labels = attr.ib(default=None)
    rename = attr.ib(default=None)
    
    def __attrs_post_init__(self):
        """ using the color attributes, set up the get_color method """
        if self.color_dict:
            self._use_color_dict()
        else:
            if isinstance(self.color_method, str):
                c_method = self.color_method.lower()[:3]
            elif callable(self.color_method):
                c_method = "con"
            if c_method == 'cat':
                # build a categorical color dict
                categories = set(self.data_dict.values())
                self.color_dict = dict(zip(categories, 
                                           gene_synteny.get_N_colors(len(categories),
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
            unique_values = self.color_dict.keys()

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
            
def cluster_synteny_plot(genome_y_vals,
                         ax1, ax2, 
                         genome_lens, 
                         gene_table,
                         genome_2_subcluster=None,
                         subcluster_colors=None,
                         gene_colors=None,
                         max_colored_genes=8,
                         max_plotted_genes=8,
                         max_genomes=40,
                         add_gene_descs=None,
                         set_titles=True,
                         label_genomes=True,
                         anchor_annot=None,
                         anchor_side=None,
                         anchor_strand=None,
                         gene_plot_range=None,
                         **kwargs,
                        ):
    
    my_genes = gene_synteny.get_cluster_genes(gene_table, set(genome_y_vals.keys()), genome_lens)
    shifted_genes = gene_synteny.align_gene_locs(my_genes,
                                                 anchor_annot,
                                                 anchor_side,
                                                 anchor_strand,
                                                 genome_lens,
                                                )

    """
    # sync X axes
    ax1.get_shared_x_axes().join(ax1, ax2)
    ax1.set_xticklabels([])
    # ax2.autoscale() ## call autoscale if needed
    """
 
    if isinstance(gene_colors, dict):
        if not isinstance(max_plotted_genes, int):
            # instead of a number of genes, it can be a list/set
            # to subset the color dict
            gene_synteny_colors = {
                g:c
                for g,c in gene_colors.items()
                if g in max_plotted_genes
            }
            max_plotted_genes = len(gene_synteny_colors)
        else:
            gene_synteny_colors = gene_colors
        syn_kwargs = dict(gene_color_dict=gene_synteny_colors)
    else:
        syn_kwargs = dict(gene_cmap=gene_colors)
        
    top_annots, calculated_gene_colors = \
        gene_synteny.plot_annot_positions(
            shifted_genes, 
            ax=ax1, 
            max_colored_genes=max_colored_genes,
            max_plotted_genes=max_plotted_genes,
            **syn_kwargs        
    )
    if not isinstance(gene_colors, dict):
        gene_colors = calculated_gene_colors
    n_genes = len(ax1.get_yticks())

    if set_titles:
        ax1.set_title(f"top {n_genes} annots of {len(set(shifted_genes.annot.values)) - 1}")

    if add_gene_descs is True:
        desc_columns = [c for c in my_genes.columns 
                       if re.search(r'^desc(ription)?$', c, flags=re.I)]
        if len(desc_columns) >= 1:
            add_gene_descs = desc_columns[-1]
    if isinstance(add_gene_descs, str) or isinstance(add_gene_descs, dict) or callable(add_gene_descs):
        gene_synteny.decorate_gene_labels(ax1, my_genes, add_gene_descs)
    
    plotted_reads = gene_synteny.draw_genes(shifted_genes, gene_colors, 
                                            ax=ax2, 
                                            genome_order=genome_y_vals,
                                            max_genomes=max_genomes, **kwargs)
    if set_titles:
        _ = ax2.set_title(f"top {len(plotted_reads)} genomes of {len(set(my_genes.genome))}")

    if label_genomes:
        if genome_2_subcluster or subcluster_colors:
            subcluster_colors = subcluster_colors if subcluster_colors is not None \
                                else get_subcluster_colors(genome_2_subcluster)
            for label in ax2.get_yticklabels():
                g = label.get_text()
                label.set_color(subcluster_colors.get(g,'black'))
    else:
        _ = ax2.set_yticks([])
 
    if gene_plot_range is None or gene_plot_range == 'genomes':
        # default: use genome extents
        gene_plot_range = ax2.get_xlim()
    elif gene_plot_range == 'genes':
        gene_plot_range = ax1.get_xlim()
    else:
        # just assume it's manually supplied interval
        pass

    ax1.set_xlim(gene_plot_range)
    ax2.set_xlim(gene_plot_range)

    # make a second copy of the genome coords at the top
    ax3 = ax2.twiny()
    ax3.set_xticks(ax2.get_xticks())
    ax3.set_xlim(gene_plot_range)

    return shifted_genes, gene_colors

def get_tree_plotter_data(genomes_tsv, 
                          genes_tsv,
                          taxdump_dir,
                          refs_tsv,
                          leaf_color_fn=lambda l: 'black',
                         ):
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

