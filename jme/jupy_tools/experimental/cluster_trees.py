
import re, numpy, pandas, os
from Bio import Align, Seq, SeqRecord, AlignIO, SeqIO
from collections import defaultdict, Counter
from functools import lru_cache
from itertools import combinations, islice
from jme.jupy_tools.experimental import gene_synteny

def cat_tree_dir_alignment(aln_files_by_vog):
    """
    Given: dict from gene name to stockholm aln file,
    
    Returns a Bio.MultipleSequenceAlignment built from multiple MSA files
    """
    alignments = {}
    gene_dict = defaultdict(dict)
    total_length = 0

    # parse alignments
    for vog, align_file in aln_files_by_vog.items():
        vog_alignment = AlignIO.read(align_file, format='stockholm')

        # keep a record of where genomes aligned
        for i, gene in enumerate(vog_alignment):
            gene_id = gene.id
            genome = re.sub(r"_\d+$", "", gene_id)
            gene_dict[genome][vog] = (i, gene_id)

        alignments[vog] = vog_alignment
        total_length += len(vog_alignment)


    if len(alignments) == 1:
        # no need to cat, just write a fasta version
        final_alignment = next(iter(alignments.values()))

    else:
        # concat the alignments
        records = []
        vog_spacers = {v:"-"*a.get_alignment_length() for v,a in alignments.items()}
        for genome in gene_dict.keys():
            genome_aligned_sequence = ""
            for vog, vog_alignment in alignments.items():
                a_index, gene_id = gene_dict[genome].get(vog, (-1,None))
                if gene_id is None:
                    genome_aligned_sequence += vog_spacers[vog]
                else:
                    genome_aligned_sequence += str(vog_alignment[a_index].seq)

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
    # count gapped genomes at each residue
    gapped_genome_counts = defaultdict(int)
    for j, alignment in enumerate(cat_alignment):
        for i, c in enumerate(alignment.seq):
            if c in {'.', '-'}:
                gapped_genome_counts[i] += 1  # add 1 for this residue

    # flag residues with <50% of genomes
    cutoff = len(cat_alignment) * residue_genome_frac
    residues_to_drop = {i
                        for i, count in gapped_genome_counts.items()
                        if count < cutoff}

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
        if .95 < align_gap_count / len(new_seq):
            # skip if genomes has < 5% of residues
            continue

        alignments.append(SeqRecord.SeqRecord(Seq.Seq(new_seq), 
                                              id=alignment.id))
    return Align.MultipleSeqAlignment(alignments,)

class TreeBuilder():
    " loads and holds dat for a clustering "
    def __init__(self, gene_table, genome_lengths):
        self.gene_table = gene_table
        self.seq_lens = genome_lengths
    
    @lru_cache(maxsize=None)
    def get_cl_50_genes(self, seqs, min_hmm_frac=.5, min_hmm_len=100):
        # return genes in these seqs that are over 50% of the HMM and over 100AA
        return gene_synteny.get_cluster_genes(self.gene_table, 
                                              seqs, 
                                              self.seq_lens) \
                    .query(f'hmm_frac >= {min_hmm_frac} and hmm_start + {min_hmm_len} <= hmm_end')

    @lru_cache(maxsize=None)
    def get_seq_vog_counts(self, seqs, cutoff=.1):
        min_count = len(seqs) * cutoff
        cl_genes_50 = self.get_cl_50_genes(seqs)
        copy_num = cl_genes_50.groupby(['genome', 'vog']).agg({'gene_no': get_copy_no})
        mean_copy_num = copy_num.reset_index().groupby('vog').agg({'gene_no': numpy.mean})
        sc_genes = set(mean_copy_num.query('gene_no <= 1.2').index)

        return copy_num \
            .reset_index()[[vog in sc_genes for genome, vog in copy_num.index]][['genome','vog']] \
            .groupby('vog').agg(len) \
            .query(f'genome >= {min_count}') \
            .sort_values('genome', ascending=False)

    @lru_cache(maxsize=None)
    def get_tree_genes(self,
                       seqs,
                       n_genes=None,
                       gene_seq_frac=.1,
                       seq_gene_frac=.05,
                       max_search_combos=5000
                      ):
        # make sure seq set is hashable and consistent for caching
        n_seqs = len(seqs)
        cl_genes_50 = self.get_cl_50_genes(seqs)
        genome_vog_counts = self.get_seq_vog_counts(seqs, cutoff=gene_seq_frac)

        if n_genes is None or len(genome_vog_counts) <= n_genes:
            n_genes = len(genome_vog_counts)
        # get the biggest set of seqs using these genes
        best_n = 0
        all_ns = []
        best_gene_sets = []
        for gene_set in islice(combinations(genome_vog_counts.index, n_genes), max_search_combos):
            gene_set = set(gene_set)
            # number of genes per genome
            genomes_hit = Counter(cl_genes_50[[(v in gene_set) for v in cl_genes_50.vog]].genome)
            # total genes in all genomes
            n = sum(c for c in genomes_hit.values() if c >= n_genes * seq_gene_frac)
            all_ns.append(n)
            if n < best_n:
                continue
            key = tuple(sorted(gene_set))
            data_tuple = (n, 
                          gene_set, 
                          genomes_hit)
            if n > best_n:
                best_gene_sets = [data_tuple,]
            else:
                best_gene_sets.append(data_tuple)

        return dict(zip(['n_genes', 'gene_set', 'cl_genomes_filt'],
                        next(iter(reversed(sorted(best_gene_sets))))))

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
            self.gene_table[[(v in tbd['gene_set'] 
                              and g not in tbd['cl_genomes_filt']) 
                             for g, v in self.gene_table[['genome','vog']].values]] \
                .groupby(['genome','vog']).agg(len) \
                .reset_index()[['genome','vog']] \
                .groupby('genome').agg(len)

        cutoff = min_gene_frac * len(tbd['gene_set'])
        max_ext_genomes = int(numpy.floor(len(tbd['cl_genomes_filt']) * max_new_genome_frac))
        new_seqs = list(other_genomes_gene_counts \
                            .query(f'vog >= {cutoff}') \
                            .vog \
                            .sort_values() \
                            .tail(max_ext_genomes) \
                            .index)
        if new_tbd:
            seqs = tuple(sorted(list(tbd['cl_genomes_filt']) + new_seqs))
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
