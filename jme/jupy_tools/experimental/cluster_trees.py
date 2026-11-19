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

Tree building (TreeBuilder, cat_tree_dir_alignment, filter_cat_alignment,
...) lives here still. The plotting code has moved:

 * TreeMetadata -> jme.jupy_tools.plot_metadata
 * TreePlotter, draw, generate_collapsed_nodes, ... -> jme.jupy_tools.tree_plots
 * cluster_synteny_plot -> jme.jupy_tools.gene_maps

and is re-exported here so old notebooks keep working.
"""
import re, numpy, pandas, os
import attr
from Bio import Align, Seq, SeqRecord, AlignIO, SeqIO, Phylo
from collections import defaultdict, Counter
from functools import lru_cache
from itertools import combinations, islice, chain
from jme.jupy_tools.utils import read_tsv
from jme.jupy_tools import gene_maps as gene_synteny
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


## Plotting functions (moved)
import warnings as _warnings
from jme.jupy_tools import plot_metadata as _pm, gene_maps as _gm, \
    tree_plots as _tp

_warnings.warn(
    "The plotting code in jme.jupy_tools.experimental.cluster_trees has "
    "moved to jme.jupy_tools.tree_plots, plot_metadata and gene_maps",
    DeprecationWarning, stacklevel=2)

_TREE_BUILDING = {k for k in list(globals()) if not k.startswith('__')}
for _module in [_pm, _gm, _tp]:
    globals().update({_k: _v for _k, _v in vars(_module).items()
                      if not (_k.startswith('__') and _k.endswith('__'))
                      and _k not in _TREE_BUILDING})
del _module
