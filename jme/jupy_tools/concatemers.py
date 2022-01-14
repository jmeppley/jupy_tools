import numpy, pandas

from collections import defaultdict
from itertools import chain
from functools import reduce

from matplotlib import pyplot as plt
from jme.jupy_tools.utils import get_best_fit
from jme.jupy_tools.hit_tables import parse_blast_m8, BLAST, BLAST_PLUS, PAF

FFT = 'fft'
CLUST = 'cluster'
JB = 'jb'

def dotplot(qhits, ax=None, subplots_kws=None, seq_length=None):
    if ax is None:
        if subplots_kws is None:
            subplots_kws = {}
        subplots_kws.setdefault("figsize", [8, 8])
        fig, ax = plt.subplots(1, 1, **subplots_kws)
    if seq_length is None:
        seq_length = qhits.mlen.max()

    for i, hit in qhits.iterrows():
        ax.plot(
            [hit.qstart / seq_length, hit.qend / seq_length],
            [hit.hstart / seq_length, hit.hend / seq_length],
            color="black",
            alpha=0.66,
        )


def get_qhits(last_out, mlen_cutoff=0, pctid_cutoff=0, table_format=BLAST, blast_plus=False, skiprows=0):
    cutoffs = {k:v
               for k,v in (('mlen', mlen_cutoff),
                           ('pctid', pctid_cutoff))
               if v > 0}
    
    # backwards compatability
    if blast_plus:
        print("DEPRECATED! use format=BLAST_PLUS instead of blast_plus=True")
        if table_format == BLAST:
            table_format = BLAST_PLUS
        elif table_format != BLAST_PLUS:
            print("WARNING: ignoring blast_plus and using format: " + table_format)

    # parse with cutoffs, and only keep parallel hits
    hit_table = parse_blast_m8(last_out, format=table_format, skiprows=skiprows, **cutoffs)
    if 'strand' in hit_table.columns:
        return hit_table.query('strand == "+"')
    else:
        return hit_table.query("((hend - hstart) > 0) == ((qend - qstart) > 0)")

def insert_self_hit(qhits, format=PAF):
    """ add a self it for tables that don't have it"""
    if format == PAF:
        query, qlength = qhits[['query','qlen']].values[0]
        new_row = pandas.Series({'query': query,
                                 'qlen': qlength,
                                 'qstart': 1,
                                 'qend': qlength,
                                 'strand': "+",
                                 'hit': query,
                                 'hlen': qlength,
                                 'hstart': 1,
                                 'hend': qlength,
                                 'matches': qlength,
                                 'mlen': qlength,
                                 'mapqv': 0,
                                 'pctid': 1,
                                }, name=-1)
        return qhits.append(new_row)
    else:
        raise Exception("I don't know how to add a self hit for format = " + format)

def evaluate_repeat(repeat_size, score, qlen, score_cutoff=0.5, frac_cutoff=0.6):
    return score > score_cutoff and repeat_size / qlen <= frac_cutoff


def get_ranks(qhits):
    # transform matches to group by position in the dotplot
    return qhits.eval("((hstart - qstart) + (hend - qend))" " / 2")


def get_rank_sizes(qhits):
    rank_sizes = defaultdict(int)
    for hit_index, rank in get_ranks(qhits).items():
        rank_sizes[rank] += qhits.loc[hit_index, "mlen"]
    return rank_sizes


def weighted_ranks(rank_sizes, factor=750):
    return numpy.array(
        list(chain(*([r,] * int(s / factor) for r, s in rank_sizes.items())))
    )


def find_cmer_stats_fft(
    qhits, resolution=10000, debug=False, kde_bw=0.05
):
    """ use KDE and FFT to get the repeat size """
    from scipy import stats, fftpack
    
    # just use positive dir hits
    qhits = qhits.query("((hend - hstart) > 0) == ((qend - qstart) > 0)")

    if qhits.query("qstart > hstart").shape[0] < 2:
        # there should be at least on long hit on each side
        return None

    # transform matches to group by position in the dotplot, add add up total hit length at each
    rank_sizes = get_rank_sizes(qhits)

    # make into a distribution using KDE
    xmin = min(rank_sizes)
    xmax = max(rank_sizes)
    size = xmax - xmin
    dx = size / resolution
    X = numpy.arange(xmin, xmax, dx)
    try:
        kernel = stats.gaussian_kde(weighted_ranks(rank_sizes), bw_method=kde_bw)
    except ValueError:
        # happens if no ranks over cutoff
        print("sequence has no self hits. Some tools, like minimap don't always return the obvious full length self hit. Try re-running using insert_self_hit() or add_full_hit=True.")
        raise
    except numpy.linalg.LinAlgError:
        # happens if only one rank over cutoff
        return None
    Z = numpy.reshape(kernel(X).T, X.shape)

    # use FFT to get dominat freq of KDE distrobution
    signal = Z - Z.mean()
    freq = numpy.linspace(0.0, 0.5 / dx, int(resolution / 2))
    signalFFT = fftpack.fft(signal) / resolution
    signalFFT = abs(signalFFT[: resolution // 2])
    peak_index = abs(signalFFT).argmax()

    repeat_size = 1 / freq[peak_index]
    if debug:
        return (
            repeat_size,
            rank_sizes,
            X,
            signal,
            abs(signalFFT[: resolution // 2]),
            freq,
        )
    else:
        return repeat_size


def find_cmer_stats_clust(
    qhits,
    seq_length=None,
    hit_cutoff=500,
    epsilon=1500,
    linkage="average",
    min_repeat_frac=0.5,
):
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.preprocessing import StandardScaler

    # just use positive dir hits
    qhits = qhits.query("((hend - hstart) > 0) == ((qend - qstart) > 0)")

    if qhits.query("qstart > hstart and mlen >= @hit_cutoff").shape[0] < 2:
        # there should be at least on long hit on each side
        return None

    # if we get here, caclulate repeat size

    # get seq length
    qlen = get_seq_length(qhits, seq_length)
    
    # transform matches to group by position in the dotplot
    ranks = get_ranks(qhits).values

    # squeeze data between 1 and -1 (roughly)
    rank_range = max(ranks) - min(ranks)
    scale = 2 / rank_range
    eps = scale * epsilon
    X = [[scale * r,] for r in ranks]

    # use Agglomerative clustering
    clustering = AgglomerativeClustering(
        linkage=linkage, n_clusters=None, distance_threshold=eps
    )
    clustering.fit(X)

    # group rank values by cluster
    clusters = {}
    for rank, label in zip(ranks, clustering.labels_):
        clusters.setdefault(label, []).append(rank)

    if len(clusters) <= 2:
        return None

    # average rank values in each cluster
    cluster_centers = sorted(
        numpy.mean([r for r in c_ranks]) for c_ranks in clusters.values()
    )

    # get diff between every cluster
    cc_diffs = [
        cluster_centers[i] - cluster_centers[i - 1]
        for i in range(1, len(cluster_centers))
    ]

    # throw out any that seem to have skipped a rank
    min_cc_diffs = min(cc_diffs)
    repeat_size = numpy.median([p for p in cc_diffs if p < 1.8 * min_cc_diffs])

    return repeat_size


def find_cmer_stats_jb(qhits, seq_length=None, hit_cutoff=3000):
    """ John B's quick method """
    qlen = get_seq_length(qhits, seq_length)
    
    # filter and sort
    qhits = qhits.query(f"mlen >= .9 * {hit_cutoff}") \
                 .sort_values('qstart')
    
    # need at least 3 hits
    if qhits.shape[0] > 2:
        
        # find median of difference of sorted hit starts
        starts = qhits.qstart.values
        diffs = [starts[i+1] - starts[i] for i in range(len(starts) - 1)]
        repeat_size = np.median(diffs)
        return repeat_size
    
def get_seq_length(qhits, seq_length):
    if seq_length is not None:
        return seq_length
    if 'qlen' in qhits.columns:
        return next(iter(qhits.qlen.values))
    else:
        return qhits.qend.values.max()

def make_presence_array(hit_range, size):
    start, end = sorted(hit_range)
    start = 1 if start <= 1 else int(start)
    if start > size:
        start = size
    end = int(end) if end < size else size

    out = numpy.empty(size)
    try:
        numpy.concatenate(
            (
                numpy.zeros(start - 1),
                numpy.ones(end + 1 - start),
                numpy.zeros(size - end),
            ),
            out=out,
        )
    except ValueError:
        print(start, end, size)
        raise
    return out


def get_score(qhits, repeat_size, qlen=None, width=0.1):
    if qlen is None:
        qlen = qhits.mlen.max()
    qhits = qhits.query("((hend - hstart) > 0) == ((qend - qstart) > 0)")
    size = int(qlen)
    upper = repeat_size * (1 + width / 2)
    lower = repeat_size * (1 - width / 2)
    coverage = reduce(
        lambda x, y: x + make_presence_array(y, size),
        qhits.eval("rank = ((hstart - qstart) + (hend - qend)) / 2")
        .query("rank <= @upper and rank >= @lower")
        .eval("mean_start = (hstart+qstart-rank)/2")
        .eval("mean_end = (hend+qend-rank)/2")[["mean_start", "mean_end"]]
        .values,
        numpy.zeros(size),
    )
    return sum(coverage > 0) / (qlen - repeat_size)


def generate_xys_from_clusters(clusters):
    for c, points in clusters.items():
        for point in points:
            yield c, point


def check_repeat_size(ranks, repeat_size, qlen, refine=True, epsilon=1500):
    clusters = {}
    for rank in ranks:
        if rank == 0:
            clusters.setdefault(0, []).append(rank)
        else:
            ranknum = int(numpy.round(rank / repeat_size))
            if abs(rank) % repeat_size < epsilon / 2:
                clusters.setdefault(ranknum, []).append(rank)

    found_clusters = len(clusters)
    expected_clusters = 1 + 2 * (int(numpy.round(repeat_size / qlen)) - 1)

    if refine:
        x, y = zip(*generate_xys_from_clusters(clusters))
        refined_rs = get_best_fit(x, y, force_intercept=True)
        return found_clusters, expected_clusters, refined_rs
    else:
        return found_clusters, expected_clusters

def build_all_v_all_cmer_table(hit_table, cmer_table_out=None,
                               self_filter=True,
                               method=FFT,
                               final_filter='notna',
                               add_full_hit=False,
                               qh_cutoff=750, **kwargs):
    """
    Loads an all v all blast hit table and makes a table of sequences that are tandem repeats

    params:
        cmer_table_out: output file. Return dataframe if None
        self_filter: if true filter incoming table for only self hits
        final_filter: what to include in the output cmer_table:
            'all': results for all sequences
            'notna': results with a non-null repeat size
            'state': results where state==True
            float(): results where score >= float value
            **kwargs are assed to get_qhits()
    """
    print(f"running {method} on {repr(hit_table)[:100]} and saving to {cmer_table_out}")
    
    if not isinstance(hit_table, pandas.DataFrame):
        # should be a dataframe or m8 file
        # load the full table
        hit_table = get_qhits(hit_table, mlen_cutoff=qh_cutoff, **kwargs)

    # filter for only self hits
    if self_filter:
        hit_table = hit_table.query('query == hit')
        
    # minimap2 doesn't return self hits
    fix_table = insert_self_hit if add_full_hit else lambda x: x

    # build the table
    cmer_table = pandas.DataFrame(
        (((q,) + check_cmer(fix_table(h), qh_cutoff=qh_cutoff, method=method))
         for q,h in hit_table.groupby('query')),
        columns=['read','repeat_size','copies','repeat_score','state','qlen']
    )

    # filter the table
    if final_filter == 'notna':
        cmer_table = cmer_table[cmer_table.repeat_size.notna()]
    elif final_filter == 'state':
        cmer_table = cmer_table.query('state')
    else:
        try:
            final_filter = float(final_filter)
            cmer_table = cmer_table[cmer_table['score'] >= final_filter]
        except ValueError:
            pass

    if cmer_table_out is not None:
        #save the table
        cmer_table.to_csv(cmer_table_out, index=None, sep='\t')
    return cmer_table


def check_cmer(qhits, method=FFT, qh_cutoff=750, **kwargs):
    if not isinstance(qhits, pandas.DataFrame):
        # should be a dataframe or m8 file
        qhits = get_qhits(qhits, mlen_cutoff=qh_cutoff, **kwargs)

    qlen = qhits.mlen.max()
    if method.lower() == FFT:
        repeat_size = find_cmer_stats_fft(qhits)
    elif method.lower() == CLUST:
        repeat_size = find_cmer_stats_clust(qhits)
    elif method.lower() == JB:
        repeat_size = find_cmer_stats_clust(qhits)
    else:
        raise Exception("Unkown cmer method: " + method)
        
    if repeat_size is not None:
        copies = qlen / repeat_size
        repeat_score = get_score(qhits, repeat_size, qlen)
        state = evaluate_repeat(repeat_size, repeat_score, qlen)
    else:
        repeat_score = -1
        state = False
        copies = None
    return repeat_size, copies, repeat_score, state, qlen

def get_cmer_check_line(input, output, candidate='', **kwargs):
    repeat_size, copies, score, state, length = check_cmer(input, **kwargs)
    with open(output, 'wt') as values_out:
        values_out.write(f"{candidate}\t"
                         f"{repeat_size}\t"
                         f"{copies}\t"
                         f"{score}\t"
                         f"{state}\t"
                         f"{length}\n")

# Scriptify for snakemake
def main():
    """ translate command line arguments into a function call: 
    the first arg is the function name
    subsequent args are passed as *args or **kwargs to the function as follows:
     * anything with a "=" is split into key/value pair for kwargs
     * all args and values are attempted to be translated into ints and floats    
    """
    
    # get args without the name of this module/script
    import sys, re
    cl_args = list(sys.argv[1:])
    print(f'processing arguments: {repr(cl_args)}')

    # get the function name and only allow word characters (\w in re)
    if len(cl_args) < 1:
        raise Exception("please supply a function name and arguments")

    function = cl_args.pop(0)
    if not re.match(r'^\w+$', function, flags=re.A):
        raise Exception("Function name can only contain word characters (a-zA-Z0-9_)")

    try:
        fn = eval(function)
    except NameError:
        print("Unknown function name: " + function)
        raise

    # parse the remaining arguments
    args = []
    kwargs = {}
    for arg in cl_args:
        split_arg = arg.split("=", 1)
        if len(split_arg) == 1:
            # single argument
            args.append(_try_to_convert_arg(arg))
        else:
            key, value = split_arg
            kwargs[key] = (_try_to_convert_arg(value))

    # run it
    print(f"running {function} with args {repr(args)} and {repr(kwargs)}")
    fn(*args, **kwargs)
    
def _try_to_convert_arg(argument):
    for cast in (int, float):
        try:
            return cast(argument)
        except ValueError:
            pass
    if argument == 'True':
        return True
    if argument == 'False':
        return False
    return argument
    

if __name__ == '__main__':
    main()
