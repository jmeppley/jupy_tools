"""
multiprocessing based tools for speeding up  distance calculations between thousands of items

includes genome-specific methods, but the core should be generic


"""
import numpy
from itertools import combinations, product
from functools import partial
from multiprocessing import Pool

from jme.jupy_tools.utils import first

def multithreaded_distances(
    data,
    dist_fn,
    chunk_root=5,
    threads=None,
    get_key_fn=first,
    mirror=True,
):
    """
    if mirror==True, save value for (key1, key2) and (key2, key1)
    """
    chunks = chunkify_data(data,
                           n=chunk_root)
    if threads is None or threads <= 0:
        threads = chunk_root * chunk_root

    _process_chunk = partial(process_chunk, process_element_pair_fn=dist_fn, get_key_fn=get_key_fn)

    results = {}
    with Pool(processes=threads) as pool:
        for chunk_results in pool.imap_unordered(_process_chunk, chunks):
            for g1, g2, score in chunk_results:
                results[(g1, g2)] = score
                if mirror:
                    results[(g2, g1)] = score

    return results

def chunkify_data(data, n=5):
    """
    Given an array of data (of any type), split into n^2 chunks for multithreaded all-v-all processing:

        * If N is the size of the array, calculate m as int(ceil(N/n))
        * break up the NxN comparisons into n^2 chunks, each m x m
        * return the chunks along the diagonal (will only process top half of each)
        * return only the chunks above the diagnol, but split each in half
    
    params:
        * data: array or list of arbitrary data
        * n: dimension of chunking

    yields:
        * n^2 chunks that should have roughly equal processing time
        * each chunk is two arrays of data to compare to each other
        * except, for the n chunks on the diagonal, the second array is None
    """
    end = 0
    size2 = int(numpy.ceil(len(data)/(n*2)))
    size1 = size2 * 2
                
    for i in range(n):
        start = end
        end = start + size1
        
        slice_a = data[start:end]
        
        # compare to self
        yield (slice_a, None)
        
        # compare to later slices
        end_b = end
        while end_b < len(data):
            start_b = end_b
            end_b = start_b + size2
            yield (slice_a, data[start_b:end_b])

def process_chunk(chunk, process_element_pair_fn, get_key_fn=first):
    """ given a pair of arrays:

         * run process_element_pair_fn on each combination of values
         * if second array is None, run all-v-all on first array

    returns:
        * list of 3-tuples: (key1, key2, values)
        * where value is the result of running process_element_pair_fn(element1, element2)
        * where key1 is the result of running get_key_fn on element1
    """
    slice_a, slice_b = chunk
    
    if slice_b is None:
        # all unique pairs in slice a
        pair_iter = combinations(slice_a, 2)
    else:
        # all combos of pairs of items from slices a and b
        pair_iter = product(slice_a, slice_b)
    
    results = []
    for element1, element2 in pair_iter:
        results.append(
            (get_key_fn(element1),
             get_key_fn(element2),
             process_element_pair_fn(element1, element2),
            )
        )
    return results



## specific methods for shared gene content in genomes
def calc_shared_gene_length_ratio(genome1_data, genome2_data):
    """
    ASsumes each element is a 3-tuple:
        * genome_id
        * dict of gene_type to gene length in genome
        * genome length

    returns the ratio of shared gene length (total in both genomes) to sum of genome lengths
    """
    (genome1, genes_lens_1, tot_len_1) = genome1_data
    (genome2, genes_lens_2, tot_len_2) = genome2_data

    # what genes are in both
    shared_genes = set(genes_lens_1).intersection(genes_lens_2)
    # what is the total length of these shared genes?
    tot_shared_gene_len = sum((genes_lens_1[g] + genes_lens_2[g] 
                               for g in shared_genes if g != "Unknown"))
    # compare to total of all genes for score
    tot_gene_len = tot_len_1 + tot_len_2

    return tot_shared_gene_len / tot_gene_len

def calculate_shared_gene_length_distances(
    gene_df,
    dist_fn=calc_shared_gene_length_ratio,
    data_kws={},
    **kwargs,
):
    """
    Cacluates the shared gene length ratio betwween every pair of genomes
    
    uses chunk_root * chunk_root threads unless specified

    Extra arguments are passed to gene_annots_to_data_array() via data_kws

    **kwargs passed to multithreaded_distances()
    """
    data = gene_annots_to_data_array(gene_df, **data_kws)
    return multithreaded_distances(data, dist_fn=dist_fn, **kwargs)

def gene_annots_to_data_array(
    gene_df,
    annot_col='annot',
    genome_col='genome',
    gene_start_col='start',
    gene_end_col='end',
    gene_len_col=None,
):
    """
    Preprocess gene annotations into data array for distance calculations

    params:
        * gene_df: pandas data frame with columns:
                * genome: genome id
                * annot: gene annotation
                * start: gtart position of gene in genome
                * end: end position of gene in genome

    if gene_len_col given, that is used instead of start/end

    returns:
        numpy array of 3-element arrays where each has:
            * genome id
            * dict of gene annotation to length in genome (may be summed over multiple copies/fragments)
            * total length of genes in genome
    """
    if gene_len_col is not None:
        data = gene_df[[gene_len_col,genome_col,annot_col]]
    else:
        gene_len_col = 'gene_len'
        data = gene_df[[gene_end_col,gene_start_col,genome_col,annot_col]]
        data[gene_len_col] = data[gene_end_col] + 1 - data[gene_start_col]

    return data \
        .groupby([genome_col, annot_col]) \
        .agg({gene_len_col:sum}) \
        .reset_index() \
        .set_index(annot_col) \
        [[genome_col, gene_len_col]] \
        .groupby(genome_col) \
        .agg({gene_len_col:(dict, sum)}) \
        .reset_index() \
        .values
