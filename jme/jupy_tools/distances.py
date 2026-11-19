"""
Calculate pairwise scores between thousands of items with multiprocessing.

The core (multithreaded_distances, chunkify_data, process_chunk) works with any
items and any pairwise function. There are also helpers for one specific
score: the fraction of two genomes' gene length that is in shared gene
annotations (calculate_shared_gene_length_distances).

    from jme.jupy_tools.distances import multithreaded_distances
    scores = multithreaded_distances(items, my_pair_function, threads=4)
    scores[(key1, key2)]

The pair function (and get_key_fn) must be picklable (defined at the top level
of a module, not a lambda) when using more than one thread.
"""
import math
import os
from functools import partial
from itertools import combinations, product
from multiprocessing import Pool

from jme.jupy_tools.utils import first

# Never use more than this many processes unless asked to
MAX_DEFAULT_THREADS = 8


def get_default_thread_count(max_threads=MAX_DEFAULT_THREADS):
    """
    Return a sensible number of worker processes: the number of available
    CPUs, but no more than max_threads.
    """
    try:
        available = len(os.sched_getaffinity(0))
    except AttributeError:
        # not available on all platforms (eg macOS)
        available = os.cpu_count() or 1
    return max(1, min(available, max_threads))


def multithreaded_distances(
    data,
    dist_fn,
    chunk_root=5,
    threads=None,
    get_key_fn=first,
    mirror=True,
):
    """
    Run dist_fn on every pair of items in data, using multiple processes.

    Parameters:
        data: a sequence (list, tuple, array) of items. Must support slicing.
        dist_fn: function taking two items and returning a score
        chunk_root: the work is split into about chunk_root^2 chunks of
            similar size
        threads: number of worker processes. None or <= 0: use
            get_default_thread_count() (all CPUs, capped at
            MAX_DEFAULT_THREADS). 1: run in this process (no multiprocessing).
        get_key_fn: function that gets an item's key (default: the first
            element of the item)
        mirror: if True, save each score as both (key1, key2) and (key2, key1)

    Returns:
        dict mapping (key1, key2) tuples to scores. Items are never compared
        to themselves.
    """
    if threads is None or threads <= 0:
        threads = get_default_thread_count()

    chunks = chunkify_data(data, n=chunk_root)
    _process_chunk = partial(
        process_chunk, process_element_pair_fn=dist_fn, get_key_fn=get_key_fn
    )

    if threads == 1:
        chunk_results_iter = map(_process_chunk, chunks)
        return _collect_results(chunk_results_iter, mirror)

    with Pool(processes=threads) as pool:
        return _collect_results(pool.imap_unordered(_process_chunk, chunks), mirror)


def _collect_results(chunk_results_iter, mirror):
    """ Merge lists of (key1, key2, score) tuples into one dict of scores """
    results = {}
    for chunk_results in chunk_results_iter:
        for key1, key2, score in chunk_results:
            results[(key1, key2)] = score
            if mirror:
                results[(key2, key1)] = score
    return results


def chunkify_data(data, n=5):
    """
    Split an all-v-all comparison of data into roughly n^2 equal-sized chunks.

    With N items, the data is cut into n "row" slices of about N/n items.
    Each row slice is yielded once paired with None (meaning: compare the
    slice to itself) and then paired with each half-size slice of all the
    data that comes after it. Every pair of items is covered exactly once.

    Parameters:
        data: sequence of arbitrary items (must support len() and slicing)
        n: the chunking dimension

    Yields:
        (slice_a, slice_b) tuples. slice_b is None for the self-comparison
        chunks.
    """
    data_len = len(data)
    if data_len == 0:
        return

    col_size = max(1, math.ceil(data_len / (n * 2)))
    row_size = col_size * 2

    for row_start in range(0, data_len, row_size):
        row_end = row_start + row_size
        slice_a = data[row_start:row_end]

        # compare to self
        yield (slice_a, None)

        # compare to everything after this row
        for col_start in range(row_end, data_len, col_size):
            yield (slice_a, data[col_start:col_start + col_size])


def process_chunk(chunk, process_element_pair_fn, get_key_fn=first):
    """
    Score every pair of items in a chunk from chunkify_data().

    If the second slice is None, every unique pair within the first slice is
    scored. Otherwise every item in the first slice is scored against every
    item in the second.

    Returns:
        list of (key1, key2, score) tuples, where the keys come from
        get_key_fn(item) and score is process_element_pair_fn(item1, item2)
    """
    slice_a, slice_b = chunk

    if slice_b is None:
        pair_iter = combinations(slice_a, 2)
    else:
        pair_iter = product(slice_a, slice_b)

    return [
        (
            get_key_fn(element1),
            get_key_fn(element2),
            process_element_pair_fn(element1, element2),
        )
        for element1, element2 in pair_iter
    ]


## specific methods for shared gene content in genomes
def calc_shared_gene_length_ratio(genome1_data, genome2_data, ignore_annots=("Unknown",)):
    """
    Return the fraction of two genomes' total gene length that is in genes
    with annotations found in both genomes.

    Each genome's data is a 3-tuple (as made by gene_annots_to_data_array):
        * genome_id
        * dict of gene annotation to total gene length in that genome
        * total gene length in the genome

    Annotations in ignore_annots are never counted as shared.

    Returns a value from 0 (nothing shared) to 1 (all genes shared). This is a
    similarity, so use (1 - value) if you need a distance.
    """
    (_genome1, gene_lens_1, tot_len_1) = genome1_data
    (_genome2, gene_lens_2, tot_len_2) = genome2_data

    tot_gene_len = tot_len_1 + tot_len_2
    if tot_gene_len == 0:
        return 0.0

    shared_annots = set(gene_lens_1).intersection(gene_lens_2).difference(ignore_annots)
    tot_shared_gene_len = sum(
        gene_lens_1[annot] + gene_lens_2[annot] for annot in shared_annots
    )
    return tot_shared_gene_len / tot_gene_len


def calculate_shared_gene_length_distances(
    gene_df,
    dist_fn=calc_shared_gene_length_ratio,
    data_kws=None,
    **kwargs,
):
    """
    Calculate the shared gene length ratio between every pair of genomes.

    Parameters:
        gene_df: pandas DataFrame of genes (see gene_annots_to_data_array)
        dist_fn: the pairwise function (default: calc_shared_gene_length_ratio)
        data_kws: dict of extra arguments for gene_annots_to_data_array()
        **kwargs: passed to multithreaded_distances() (eg threads, chunk_root)

    Returns:
        dict mapping (genome1, genome2) to the ratio
    """
    data = gene_annots_to_data_array(gene_df, **(data_kws or {}))
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
    Summarize gene annotations by genome for calculate_shared_gene_length_distances().

    Parameters:
        gene_df: pandas DataFrame with one row per gene and columns for:
            * the genome id (genome_col)
            * the gene annotation (annot_col)
            * the gene start and end (1-based, inclusive) positions, or
            * the gene length (gene_len_col, if given)

    Returns:
        list of (genome_id, {annotation: total_length}, total_gene_length)
        tuples, one per genome (sorted by genome id). Lengths from multiple
        genes with the same annotation are summed.
    """
    if gene_len_col is not None:
        gene_lens = gene_df[gene_len_col]
    else:
        gene_lens = gene_df[gene_end_col] + 1 - gene_df[gene_start_col]

    annot_lens = (
        gene_lens
        .groupby([gene_df[genome_col], gene_df[annot_col]])
        .sum()
    )

    return [
        (genome, genome_lens.droplevel(0).to_dict(), genome_lens.sum())
        for genome, genome_lens in annot_lens.groupby(level=0)
    ]
