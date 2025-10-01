import logging
import re
import numpy
import pandas
from itertools import islice
from functools import partial
from contextlib import contextmanager
from collections import defaultdict
from jme.jupy_tools.utils import iterable_to_stream, first

LOGGER = logging.getLogger(__name__)

# params for getting lastal to behave like blastx (from E Ottesen)
LASTP_PARAMS = "-b 1 -x 15 -y 7 -z 25"

BLAST_COLUMNS = [
        "query",
        "hit",
        "pctid",
        "mlen",
        "mismatches",
        "gaps",
        "qstart",
        "qend",
        "hstart",
        "hend",
        "evalue",
        "score",
]
BLAST_PLUS_COLUMNS = BLAST_COLUMNS + [
    'qlen', 'hlen', 'raw_score'
]
PAF_COLUMNS = ['query','qlen','qstart','qend','strand',
               'hit', 'hlen','hstart','hend','matches',
               'mlen','mapqv']

BLAST = 'BlastTab'.upper()
BLAST_PLUS = 'BlastTab+'.upper()
PAF = 'Paf'.upper()
PAF_ALL = 'Paf+'.upper()
LASTAL = 'lastal'.upper()
LASTAL_COLUMNS = ['score',
                  'hit', 'hstart', 'hmlen', 'hstrand', 'hlen',
                  'query', 'qstart', 'qmlen', 'qstrand', 'qlen',
                  'match_string', 'eg2', 'e']


HIT_TABLE_COLUMNS = {BLAST: BLAST_COLUMNS,
                     BLAST_PLUS: BLAST_PLUS_COLUMNS,
                     PAF: PAF_COLUMNS,
                     PAF_ALL: PAF_COLUMNS + ['tp','cm','dv','rl'],
                     LASTAL: LASTAL_COLUMNS,
                    }

paf_rexp = re.compile(r"""^(?P<query>\S+)\t
                          (?P<qlen>\d+)\t
                          (?P<qstart>\d+)\t
                          (?P<qend>\d+)\t
                          (?P<strand>[+-])\t
                          (?P<hit>\S+)\t
                          (?P<hlen>\d+)\t
                          (?P<hstart>\d+)\t
                          (?P<hend>\d+)\t
                          """,
                      re.VERBOSE)

blast_plus_rexp = re.compile(r"""^
                        (?P<query>\S+)\t
                        (?P<hit>\S+)\t
                        (?P<pctid>[0-9.]+)\t
                        (?P<mlen>\d+)\t
                        (?P<mismatches>\d+)\t
                        (?P<gaps>\d+)\t
                        (?P<qstart>\d+)\t
                        (?P<qend>\d+)\t
                        (?P<hstart>\d+)\t
                        (?P<hend>\d+)\t
                        (?P<evalue>\S+)\t
                        (?P<score>\S+)\t
                        (?P<qlen>\d+)\t
                        (?P<hlen>\d+)\t
                        (?P<raw_score>[0-9.])
                        """,
                             flags=re.VERBOSE)

blast_tab_rexp = re.compile(r"""^
                        (?P<query>\S+)\t
                        (?P<hit>\S+)\t
                        (?P<pctid>[0-9.]+)\t
                        (?P<mlen>\d+)\t
                        (?P<mismatches>\d+)\t
                        (?P<gaps>\d+)\t
                        (?P<qstart>\d+)\t
                        (?P<qend>\d+)\t
                        (?P<hstart>\d+)\t
                        (?P<hend>\d+)\t
                        (?P<evalue>\S+)\t
                        (?P<score>\S+)
                        """,
                             flags=re.VERBOSE)

HIT_TABLE_REXPS = {PAF: paf_rexp, BLAST: blast_tab_rexp, BLAST_PLUS: blast_plus_rexp}

def computeLastHitValues(blocks):
    """
    given the query length and 'blocks' string from a last hit
    return the:
        match length

    the blocks string looks something like this:
        "73,0:1,15,0:1,13,0:1,9"
        where integer elements indicate lenghts of matches and
        colon separated elements indicate lengths of gaps
    """
    matchLen = 0
    for segment in blocks.split(','):
        try:
            matchLen += int(segment)
        except ValueError:
            (hml, qml) = segment.split(":")
            mml = max(int(hml), int(qml))
            matchLen += mml

    return matchLen


def parse_blast_m8(
    hit_table,
    format=BLAST,
    skiprows=0,
    pair_overlap_buffer=-1,
    use_polars=False,
    **cutoffs
):
    """ utility for quickly loading a hit table into pandas 
    
        cutoff keys should match column names. Use negative cutoff value to accept all values.
        
        if pair_overlap_buffer is set to a nonzero value, hits that overlap previous hits 
        (in the hit_table order) between a given query/hit pair are dropped. This happens
        before any other curofs are applied.

        if use_polars is not False, use polars to load and filter. If
        use_polars is set to "lazy", return the LazyFrame (ie do not collect)
    """

    if pair_overlap_buffer >= 0:
        # call this function again after filtering for overlaps
        if format not in HIT_TABLE_REXPS:
            raise Exception("Overlap filtering is not supported for format: " + format)
        with remove_pair_overlaps(hit_table, buffer=pair_overlap_buffer, format=format) as new_table:
            return parse_blast_m8(new_table, format=format, skiprows=skiprows,
                                  use_polars=use_polars, **cutoffs)

    column_names = HIT_TABLE_COLUMNS[format]
    use_cols = list(range(len(column_names)))

    if use_polars:
        return _parse_blast_m8_polars(
            hit_table, format, skiprows, column_names,
            lazy=(isinstance(use_polars, str)
                  and (use_polars.upper() == "LAZY")),
            **cutoffs,
        )

    hits = \
        pandas.read_csv(
            hit_table,
            sep="\t",
            header=None,
            comment="#",
            usecols=use_cols,
            names=column_names,
            skiprows=skiprows,
        )

    if hits.shape[0] == 0:
        ## if there are no hits, return None
        return None

    # format specific tweaks
    if format == LASTAL:
        hits['qsmult'] = [1 if qs == '+' else -1 
                          for qs in hits.qstrand]
        hits = \
            hits.eval('hend = hstart + hmlen - 1') \
                .eval('qend = qstart + ((qmlen - 1) * qsmult)')
        hits['mlen'] = [computeLastHitValues(blocks)
                        for blocks in hits.match_string]
        hits['evalue'] = [float(re.sub('E=','',str(e)))
                          for e in hits['e']]
        
    if format in [PAF, PAF_ALL]:
        # calculate pctid
        hits = hits.eval('pctid = 100 * matches / mlen')
    
    query = " and ".join(f"{column} >= {value}" 
                         if value >=0 
                         else f"{column} <= {-1*value}"
                         for column, value in cutoffs.items()
                         if numpy.abs(value) > 0
                        )
    if len(query) > 0:
        return hits.query(query)
    return hits

def _parse_blast_m8_polars(
    hit_table,
    format,
    skiprows,
    column_names,
    lazy=False,
    scan_csv_args=None,
    **cutoffs,
):
    import polars as pl
    if scan_csv_args is None:
        scan_csv_args = {} 
    hits = \
        pl.scan_csv(
            hit_table,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            skip_rows=skiprows,
            new_columns=column_names,
            **scan_csv_args,
        ).select(
            *(
                pl.col(c) for c in column_names
            )
        )

    # format specific tweaks
    if format == LASTAL:
        hits = hits.with_columns(
            pl.col('qstrand').map_elements(
                lambda qs: 1 if qs == '+' else -1,
                return_dtype(pl.Int8)
            ),
            (
                pl.col('hstart') + pl.col('hmlen') + 1
            ).alias('hend'),
            (
                pl.col('qstart') + (
                    pl.col('qmlen') -1
                ) * pl.col('qsmult')
            ).alias('qend'),
            pl.col('match_string').map_elements(
                computeLastHitValues,
                return_dtype=pl.UInt32,
            ).alias('mlen'),
            pl.col('e').map_elements(
                lambda e: float(re.sub('E=','',str(e))),
                return_dtype=pl.Float64,
            ).alias('evalue')
        )

    if format in [PAF, PAF_ALL]:
        # calculate pctid
        hits = hits.with_columns(
            (
                pl.col('matches') * 100 / pl.col('mlen')
            ).alias('pctid')
        )
    
    filters = [
        (
            (pl.col(column) >= value)
            if value >= 0
            else
            (pl.col(column) <= -1 * value)

        )
        for column, value in cutoffs.items()
        if numpy.abs(value) > 0
    ]
    if len(filters) > 0:
        hits = hits.filter(*filters)

    if lazy:
        return hits

    # collect
    hits = hits.collect()

    if hits.shape[0] == 0:
        ## if there are no hits, return None
        return None

    return hits


@contextmanager
def remove_pair_overlaps(hit_table, **params):
    """ wrapper to  generate_nonoverlapping_lines that (1) returns a file-like
    object instead of an iterator and (2) serves as a context manager """
    yield iterable_to_stream(generate_nonoverlapping_lines(hit_table, **params))


def generate_nonoverlapping_lines(hit_file, format=BLAST_PLUS, buffer=0):
    """
    Loops over lines in a hit_table keeping track of hit and query positions. Currently
    supported formats are paf and blasttab+

    Only yields lines that don't overlap with previous hits between a query/hit pair
    """
    
    if isinstance(hit_file, str):
        # open file handle and call this function again
        with open(hit_file) as lines:
            yield from generate_nonoverlapping_lines(lines, format=format, buffer=buffer)
    else:
        hit_rexp = HIT_TABLE_REXPS[format]
        hit_regions_by_pair = defaultdict(list)
        query_regions_by_pair = defaultdict(list)

        for line in hit_file:
            if line.startswith("#"):
                # skip comment lines
                continue
            
            try:
                hit_data = hit_rexp.search(line).groupdict()
                query, qstart, qend, hit, hstart, hend = \
                    [hit_data[k] for k in ['query', 'qstart', 'qend', 'hit', 'hstart', 'hend']]
            except AttributeError:
                print("ERROR on line: \n" + line)
                raise
            qstart = int(qstart)
            qend = int(qend)
            hstart = int(hstart)
            hend = int(hend)
            
            pair_key = (hit, query)
            hit_regions = hit_regions_by_pair[pair_key]
            query_regions = query_regions_by_pair[pair_key]
            
            LOGGER.debug('processing line: \n%s', line)

            for regions, start, end in [(hit_regions, hstart, hend),
                                        (query_regions, qstart, qend)]:
                if start > end:
                    start, end = end, start
                    
                LOGGER.debug('checking for %d - %d for hit in \n%r', start, end, regions)

                # if hit is shorter than buffer, then keep it
                if end - start <= buffer:
                    continue

                # adjust for buffer
                buf_start, buf_end = start + buffer, end - buffer

                # check overlap for all ranges
                for already_occ_range in regions:
                    if (buf_start >= already_occ_range[1] 
                        or buf_end <= already_occ_range[0]):                    
                        # does not overlap this range (try next range) 
                        continue
                    else:
                        # overlaps a previous hit...move on
                        LOGGER.debug("overlaps %r", already_occ_range)
                        break
                else:
                    # there was no overlap, try next pair or exit loop
                    LOGGER.debug("no overlap here")
                    continue
                # there was some sort of collision
                break
            else:
                # the hit was accepted by both hstart,hend and qstart,qend
                LOGGER.debug("no overlaps")
                hit_regions.append(tuple(sorted((hstart,hend))))
                query_regions.append(tuple(sorted((qstart, qend))))
                yield line

def _check_columns(frame, col, is_lazy=False):
    if is_lazy:
        return col in frame.collect_schema()
    return col in frame.columns

def agg_hit_polars_df(non_ovl_hits, max_hit_fragments=0, keep_extent=False):
    import polars as pl
    
    if isinstance(non_ovl_hits, pl.LazyFrame):
        lazy = True
    elif isinstance(non_ovl_hits, pl.DataFrame):
        lazy = False
    else:
        raise Exception("I only know polars and panads dataframes. I can't handle type: " + 
                        str(type(non_ovl_hits)))

    columns = set(non_ovl_hits.collect_schema() if lazy else non_ovl_hits.columns)
    
    if 'matches' not in columns:
        non_ovl_hits = non_ovl_hits.with_columns(
            ((pl.col('pctid') * pl.col('mlen')) / 100).alias('matches')
        )

    # just summing the first N fragments is not yet supported for polars
    if max_hit_fragments > 0:
        raise Exception("Limiting hit fragments in aggregation with max_hit_fragments is not yet supported.")

    # add up the match counts and match lengths
    agg_col_exprs = [
        pl.col('matches').sum(),
        pl.col('mlen').sum()
    ]

    # keep one (the first) if either query or hit length is present
    for col in ['qlen', 'hlen']:
        if col in columns:
            agg_col_exprs.append(
                pl.col(col).first()
            )

    # get lengths and extents if we can and if asked
    updated_column_exprs = []
    for pref in ['q','h']:
        # do for both query and hit start/end pairs
        scol, ecol = [pref+col for col in ['start','end']]

        # only if the columns exist
        if scol in columns and ecol in columns:

            # force start/end pairs to be sorted (simplifies next 2 things)
            updated_column_exprs.extend([
                pl.when(pl.col(scol) > pl.col(ecol)).then(pl.col(ecol)).otherwise(pl.col(scol)).alias(scol),
                pl.when(pl.col(scol) > pl.col(ecol)).then(pl.col(scol)).otherwise(pl.col(ecol)).alias(ecol),
            ])

            # get max extent for all hits if requested
            if keep_extent:
                # now the start is the smallest start and the end is the biggest end
                agg_col_exprs.extend([
                    pl.col(scol).min(),
                    pl.col(ecol).max()
                ])

            # calculate query and hit specific alignment lengths
            agg_col_exprs.append(
                (pl.col(ecol) + 1 - pl.col(scol)).sum().alias(f"{pref}mlen")
            )

    post_agg_cols = [
        (pl.col('matches') * 100 / pl.col('mlen')).alias('pctid')
    ]
    if 'qlen' in columns and 'hlen' in columns:
        post_agg_cols.append(
            (pl.col('matches') * 200 / (pl.col('hlen') + pl.col('qlen'))).alias('mfrac')
        )
    
    agg_hits = non_ovl_hits \
        .with_columns(updated_column_exprs) \
        .group_by(['query','hit']) \
        .agg(agg_col_exprs) \
        .with_columns(post_agg_cols)

    return agg_hits
    
def agg_hit_df(non_ovl_hits, max_hit_fragments=0, keep_extent=False):
    if not isinstance(non_ovl_hits, pandas.DataFrame):
        # try polars
        return agg_hit_polars_df(non_ovl_hits, max_hit_fragments, keep_extent)
        
    # aggregate all hits by hit/query pair
    if 'matches' not in non_ovl_hits.columns:
        non_ovl_hits = non_ovl_hits.eval('matches = pctid * mlen / 100')
        
    # allow for the option to just sum the first N fragments
    if max_hit_fragments > 0:
        def sum_first_N(S):
            return sum(islice(S, max_hit_fragments))
        non_ovl_hits = non_ovl_hits.sort_values('matches', ascending=False)
    else:
        sum_first_N = sum

    # addup match ids and match lengths
    agg_dict = {'matches':sum_first_N,
                'mlen':sum_first_N}

    # also save qlen and hlen if present
    for col in ['qlen', 'hlen']:
        if col in non_ovl_hits.columns:
            agg_dict[col] = first

    # get lengths and extents if we can and if asked
    for pref in ['q','h']:
        # do for both query and hit start/end pairs
        scol, ecol = [pref+col for col in ['start','end']]

        # only if the columns exist
        if scol in non_ovl_hits.columns and ecol in non_ovl_hits.columns:

            # force start/end pairs to be sorted (simplifies next 2 things)
            non_ovl_hits[[scol, ecol]] = [
                sorted(v) for v in non_ovl_hits[[scol, ecol]].values
            ]

            # get max extent for all hits if requested
            if keep_extent:
                # now the start is the smallest start and the end is the biggest end
                agg_dict.update({scol: min, ecol: max})

            # calculate query and hit specific alignment lengths
            non_ovl_hits = non_ovl_hits.eval(
                f"{pref}mlen = 1 + {ecol} - {scol}"
            )
            agg_dict[pref + "mlen"] = sum
    
    agg_hits = non_ovl_hits \
        .groupby(['query','hit']) \
        .agg(agg_dict) \
        .eval('pctid=100 * matches / mlen') \
        .reset_index()
    
    # calculate mfrac if we can
    if 'qlen' in agg_hits.columns and 'hlen' in agg_hits.columns:
        agg_hits = agg_hits.eval('mfrac=200 * matches / (hlen + qlen)')
        
    return agg_hits    


def agg_hit_df_polars(non_ovl_hits, max_hit_fragments=0, keep_extent=False):
    # aggregate all hits by hit/query pair
    import polars as pl

    columns = set(
        non_ovl_hits.collect_schema().names()
        if isinstance(non_ovl_hits, pl.LazyFrame)
        else non_ovl_hits.schema.names()
    )
    if 'matches' not in columns:
        non_ovl_hits = non_ovl_hits.with_columns(
                    (pl.col('pctid') * pl.col('mlen') /
                     100).alias('matches')
                )

    # allow for the option to just sum the first N fragments
    if max_hit_fragments > 0:
        non_ovl_hits = non_ovl_hits.sort(
            pl.col('matches'), descending=True
        ).group_by(
            pl.col('query'), pl.col('hit')     # for each query/hit pair
        ).head(
            max_hit_fragments   # take the first N hits
        )

    # addup match ids and match lengths
    agg_args = [
        pl.col('matches').sum(),
        pl.col('mlen').sum(),
    ]

    # also save qlen and hlen if present
    for col in ['qlen', 'hlen']:
        if col in columns:
            agg_args.append(pl.col(col).first())
    
    
    # get lengths and extents if we can and if asked                                                
    for pref in ['q','h']:
        # do for both query and hit start/end pairs
        scol, ecol = [pref+col for col in ['start','end']]                                          
        
        # only if the columns exist
        if scol in columns and ecol in columns:                           
            agg_args.append(
                ((pl.col(ecol) - pl.col(scol)).abs() + 1).alias(f"{pref}mlen").sum()
            )

            # get max extent for all hits if requested                                              
            if keep_extent:
                # now the start is the smallest start and the end is the biggest end                
                agg_args.extend([
                    pl.min_horizontal(pl.col(scol), pl.col(ecol)).min(),
                    pl.max_horizontal(pl.col(scol), pl.col(ecol)).max(),
                ])                                                         
    
    agg_hits = non_ovl_hits.group_by(
        pl.col('query'), 
        pl.col('hit'),
    ).agg(
        *agg_args
    ).with_columns(
        (pl.col('matches') * 100 / pl.col('mlen')).alias('pctid')
    )
    

    # calculate mfrac if we can
    if ('qlen' in columns) and ('hlen' in columns):
        agg_hits = agg_hits.with_columns(
            (
                pl.col('matches') * 200 / (
                    pl.col('hlen') + pl.col('qlen')
                )
            ).alias('mfrac')            
        )

    return agg_hits


def agg_hit_table(hit_table, ovl_buffer=0, max_hit_fragments=0,
                  keep_extent=False, **parse_args):
    """
    for each hit/query pair return one line of data by merging multiple hit fragments
    
    params:
    
     * ovl_buffer: allow merged hit fragments to overlap by this much (0)
     * max_hit_fragments: only merge the N best fragments (unlimited)
    """
    non_ovl_hits = parse_blast_m8(hit_table, pair_overlap_buffer=ovl_buffer, **parse_args)
    
    if non_ovl_hits is None:
        return None

    if ('use_polars' in parse_args) and (parse_args['use_polars']):
        return agg_hit_df_polars(non_ovl_hits, max_hit_fragments, keep_extent)

    return agg_hit_df(non_ovl_hits, max_hit_fragments, keep_extent)

## SAM flags
bit_names = [
    'paired', 'mapped in pair', 'unmapped', 'mate unmapped', 
    'read rev strand', 'mate rev strand', 
    'first in pair', 'second in pair', 'not primary alig', 'fails checks', 
    'duplicate', 'supplementary',
]

def decode_sam_flags(flag):
    # hack using the binary representation of the flag integer
    # I couldn't get it to pad zeros properly, so I'm adding a high bit that I cut off at the end
    bits = list(reversed(f"{bin(int(flag) + 4096)}"))[:-3]
    assert len(bits) == len(bit_names)
    return {
        name: (bit == "1")
        for name, bit 
        in zip(bit_names, bits)
    }

def is_sam_flag_primary_hit(flag):
    flag_bits = decode_sam_flags(flag)
    return (
        (not flag_bits["unmapped"])
        and
        (not (
            flag_bits["supplementary"]
            or
            flag_bits["not primary alig"]
        ))
    )

def is_sam_flag_first_in_pair(flag):
    return decode_sam_flags(flag)['first in pair']

