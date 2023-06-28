"""
Collection of functions for processing viral and satellite genomes from nanopore sequencing runs

"""
import os, re
import pandas, numpy

from collections import defaultdict

from jme.jupy_tools import hit_tables

def get_genome_covs(hits, depth_file):
    """ Given the hit table retunred by get_best_genome_hits
    
    And a smatools depth file with base-by-base depth for each contig
    
    Return the "mean_cov" for each genome
    """
    # get coverage depths by position and index on contig+position
    depths_df = pandas.read_csv(
        depth_file, 
        sep='\t',
        index_col=None, 
        header=None, 
        names=['contig','pos','depth']
    ).set_index(['contig','pos'])

    for genome, genome_hits in hits.groupby('query'):
        tot_bases_of_cov = 0
        positions_covered = set()
        tot_cov_by_base = defaultdict(lambda: 0)
        for contig, start, end, qstart, qend in genome_hits[['hit', 'hstart', 'hend', 'qstart', 'qend']].values:
            # get the total cov of the hit
            for pos in numpy.arange(start, end + 1):
                depth = depths_df['depth'].get((contig, pos), 0)
                tot_bases_of_cov += depth
                
            # which genome positions are covered
            positions_covered.update(numpy.arange(qstart, qend + 1))
            
        # get mean depth value over genome positions
        mean_cov = tot_bases_of_cov / len(positions_covered)
        yield genome, mean_cov

def hits_to_bed(hit_table, out_file, contig_col='hit', start_col=None, end_col=None, name_col='query'):
    """
    Given a hit table with hits as contigs, write out a bed file with locations of the hits.
    
    By default contig names come from 'hit' and location from 'hstart', 'hend'. If contig_col
    is set to something else, start/end columns will us the first letter of the contig_col name
    if not also set manually.
    
    If name_col is not None, a fourth "name" column is added to the bed file pulling from that col.
    (DEfault is the query name)
    """
    if start_col is None:
        start_col = contig_col[0] + 'start'
    if end_col is None:
        end_col = contig_col[0] + 'end'
    if name_col is None:
        name_col = contig_col

    with open(out_file, 'wt') as bed_out:
        for contig, start, end, name in hit_table[[contig_col, start_col, end_col, name_col]].values:
            bed_out.write(f"{contig}\t{start}\t{end}\t{name}\n")

def get_best_any_hits(paf_file, min_pctid=80, min_hit_frac=.8):
    """
    Given a paf result mapping genomes to contigs: 
    
    return a condensed table of best non-opverlapping contig hits for each genome
    
    Steps:
     * aggregate hits by genome-contig pairs
     * filter hits (80% ID over at least 80% of the contig)
     * drop hits where contig also matches a different genome better
     
     NOTE: we expect the genomes to be complete and the contigs to be fragmented, so 
     we're assuming genomes are longer and multiple contigs can hit the same genome 
     
     """
    best_hits = []
    if os.path.getsize(paf_file) > 0:
        # aggregate all alignments for each query/hit pair
        hits = hit_tables.agg_hit_table(paf_file, format=hit_tables.PAF, keep_extent=True)

        # only keep hits that cover at least 80% of the query with at least 80% ID
        hits['mmfrac'] = hits.mlen / hits[['qlen', 'hlen']].min(axis=1)
        good_hits = hits.query(f'mmfrac >= {min_hit_frac} and pctid >= {min_pctid}')

        # contigs can only be used once each
        used_contigs = set()
        cols = ['query', 'hit', 'qstart', 'qend', 'hstart', 'hend', 'matches', 'mlen', 'pctid']
        for values in good_hits.sort_values('mmfrac', ascending=False)[cols].values:
            genome, contig, qstart, qend, hstart, hend, matches, mlen, pctid = values
            if contig in used_contigs:
                continue
            used_contigs.add(contig)

            best_hits.append(
                [genome, contig, matches, mlen, pctid, qstart, qend, hstart, hend]
            )
                        
    return pandas.DataFrame(
        best_hits,
        columns=['query', 'hit', 'matches', 'mlen', 'pctid', 'qstart', 'qend', 'hstart', 'hend'],
    )

def get_best_genome_hits(paf_file, min_pctid=80, min_hit_frac=.8):
    """
    Given a paf result mapping genomes to contigs: 
    
    return a condensed table of best non-opverlapping contig hits for each genome
    
    Steps:
     * aggregate hits by genome-contig pairs
     * filter hits (80% ID over at least 80% of the contig)
     * drop hits where contig also matches a different genome better
     
     NOTE: we expect the genomes to be complete and the contigs to be fragmented, so 
     we're assuming genomes are longer and multiple contigs can hit the same genome 
     
     """
    best_hits = []
    if os.path.getsize(paf_file) > 0:
        # aggregate all alignments for each query/hit pair
        hits = hit_tables.agg_hit_table(paf_file, format=hit_tables.PAF, keep_extent=True)

        # only keep hits that cover at least 80% of the query with at least 80% ID
        good_hits = hits \
            .eval(f'hmfrac = hmlen / hlen') \
            .query(f'hmfrac >= {min_hit_frac} and pctid >= {min_pctid}')

        # contigs can only be used once each
        used_contigs = set()
        cols = ['query', 'hit', 'qstart', 'qend', 'hstart', 'hend', 'matches', 'mlen', 'pctid']
        for values in good_hits.sort_values('hmfrac', ascending=False)[cols].values:
            genome, contig, qstart, qend, hstart, hend, matches, mlen, pctid = values
            if contig in used_contigs:
                continue
            used_contigs.add(contig)

            best_hits.append(
                [genome, contig, matches, mlen, pctid, qstart, qend, hstart, hend]
            )
                        
    return pandas.DataFrame(
        best_hits,
        columns=['query', 'hit', 'matches', 'mlen', 'pctid', 'qstart', 'qend', 'hstart', 'hend'],
    )

                        
def get_best_veime_hits(paf_file, min_pctid=80, min_hit_frac=.8):
    """
    Given a paf result mapping VEIMEs to contigs: 
    
    return a condensed table of best non-opverlapping contig hits for each VEIME
    
    Steps:
     * aggregate hits by VEIME-contig pairs
     * filter hits (80% ID over at least 80% of the VEIME)
     * drop hits to that overlap better hits to different contigs or VEIMEs
     
     NOTE: We expect VEIMEs and contigs to be the same size or the contig to be larger (for integrated VEIMEs)
     so multiple VEIMEs can hit the same contig as long as it's not the same location.
     
     We also allow a VEIME to match to multiple contigs if there are no better VEIME matches for either.
     """
    best_hits = []
    if os.path.getsize(paf_file) > 0:
        # aggregate all alignments for each query/hit pair
        hits = hit_tables.agg_hit_table(paf_file, format=hit_tables.PAF, keep_extent=True)

        # only keep hits that cover at least 80% of the query with at least 80% ID
        good_hits = hits \
            .eval(f'qmfrac = qmlen / qlen') \
            .query(f'qmfrac >= {min_hit_frac} and pctid >= {min_pctid}')

        # filter to best non-overlapping hits on each VEIMEs
        contig_hit_ranges = defaultdict(list)
        cols = ['query', 'hit', 'qstart', 'qend', 'hstart', 'hend', 'matches', 'mlen', 'pctid']
        for values in good_hits.sort_values('qmfrac', ascending=False)[cols].values:
            veime, contig, qstart, qend, hstart, hend, matches, mlen, pctid = values
            for hit_range in contig_hit_ranges[contig]:
                if intersects(hit_range, (hstart, hend)):
                    break
                else:
                    # this range oesn't overlap any others, keep it
                    contig_hit_ranges[contig].append((start, end))
                    best_hits.append(
                        [veime, contig, matches, mlen, pctid, qstart, qend, hstart, hend]
                    )
                        
    return pandas.DataFrame(
        best_hits,
        columns=['query', 'hit', 'matches', 'mlen', 'pctid', 'qstart', 'qend', 'hstart', 'hend'],
    )


def get_read_counts_from_report(report_file):
    read_counts = 0
    with open(report_file) as lines:
        line_iter = iter(lines)
        line = next(line_iter)
        while not line.strip().startswith('reads:'):
            line = next(line_iter)

        for line in line_iter:                                                                  
            if line.lstrip().startswith("###"):
                break
            m = re.search( 
                r'^\s*(\S+)\s+clean\s+reads:\s*(\d+)',                                          
                line,                                                                           
            )
            if m:
                sample, read_count = m.groups()                                                 
                read_counts += int(read_count)
    return read_counts
