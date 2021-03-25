"""
Utilities for inspecting and comparing MAGs
"""
import numpy
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from itertools import cycle

def get_contig_starts(contig_order, contig_lens, buffer=500):
    cumlen = 0
    cumlens = {}
    for contig in contig_order:
        cumlens[contig] = cumlen
        cumlen += contig_lens[contig] + buffer
    return cumlens

def random_color():
    return numpy.random.rand(3,)


def pseudo_mauve(hits, ax,
                 contig_line_props=dict(linewidth=5),
                 connection_alpha=.6,
                 connection_colormap='cool',
                 min_pctid=None,
                 contig_colors=['darkgrey', 'lightgrey'],
                 hit_y=0,
                 query_y=20,
                 buffer_x=500):

    # color connections by pctid
    if min_pctid is None:
        # get the true minimum %ID from data
        min_pctid = min(hits['pctid'])
    elif min_pctid < 0:
        # if negative, go under true min by that much
        min_pctid += min(hits['pctid'])
        
    pctid_range = 100 - min_pctid
    if isinstance(connection_colormap, str):
        connection_colormap = plt.get_cmap(connection_colormap)
    def get_conn_color(pctid):
        return connection_colormap((pctid - min_pctid) / pctid_range)

    # contig coloring
    query_color_iter = cycle(contig_colors)
    hit_color_iter = cycle(contig_colors)

    max_x = 0

    # sort hits by length
    hit_lens = {h:hl for h, hl in hits[['hit','hlen']].values}
    hit_order = sorted(hit_lens, key=lambda h: hit_lens[h], reverse=True)
    hit_starts = get_contig_starts(hit_order, hit_lens, buffer_x)

    # sort queries
    query_lens = {q:ql for q, ql in hits[['query','qlen']].values}
    
    # sort queries by mean hit start
    query_mean_hstarts = {query: numpy.mean([hstart + hit_starts[hit]
                                             for hit, hstart in data[['hit','hstart']].values])
                          for query, data in hits.groupby('query')}
    
    # sort queries by weighted hit starts
    """
    query_mean_hstarts = {}
    for query, data in hits.groupby('query'):
        weighted_sum, weigths_sum = 0, 0
        for hit, hstart, mlen in data[['hit','hstart','mlen']].values:
            weigths_sum += mlen
            weighted_sum += hstart + hit_starts[hit]
        query_mean_hstarts[query] = weighted_sum / weigths_sum
    """
    
    query_order = sorted(query_mean_hstarts, key=lambda q: query_mean_hstarts[q])
    
    query_starts = get_contig_starts(query_order, query_lens, buffer_x)
    
    # draw all queries
    for query, q_zero in query_starts.items():
        qlen = query_lens[query]
        ax.plot((q_zero, q_zero + qlen), (query_y, query_y), 
                color=next(query_color_iter),
                **contig_line_props)
        max_x = max(max_x, q_zero + qlen + buffer_x)
        
    # draw all hits
    for hit in hit_order:
        hlen = hit_lens[hit]
        h_zero = hit_starts[hit]
        ax.plot((h_zero, h_zero + hlen), (hit_y, hit_y), color=next(hit_color_iter), **contig_line_props)
        max_x = max(max_x, h_zero + hlen + buffer_x)

    # draw connections with highest pctid on top
    for query, qstart, qend, hit, hstart, hend, pctid in hits \
                [['query', 'qstart', 'qend', 'hit', 'hstart', 'hend', 'pctid']] \
                .sort_values('pctid') \
                .values:
        q_zero = query_starts[query]
        h_zero = hit_starts[hit]
        
        xy = numpy.array([(hstart + h_zero, hit_y),
                          (qstart + q_zero, query_y),
                          (qend + q_zero, query_y),
                          (hend + h_zero, hit_y)])
        ax.add_patch(Polygon(numpy.array(xy), 
                             color=get_conn_color(pctid), 
                             alpha=connection_alpha))

    xl = ax.set_xlim(0, max_x)
    yl = ax.set_ylim(hit_y - 1, query_y + 1)

    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    ax.tick_params(axis='both', which='both', length=0)
    return min_pctid

def get_sorted_hits(hit_table, other_cols=None):
    
    # get contig lengths from hit table (qlen and hlen columns)
    hit_lens = {h:hl for h, hl in hit_table[['hit','hlen']].values}
    # sort hits by length
    hit_order = sorted(hit_lens, key=lambda h: hit_lens[h], reverse=True)
    hit_starts = get_contig_starts(hit_order, hit_lens)

    # we sort on the best hit here (instead of mean) to get a better diagonal
    query_lens = {q:ql for q, ql in hit_table[['query','qlen']].values}
    query_hit_order = {q:hit_order.index(data.eval('score = mlen * pctid').sort_values('score')['hit'].values[0])
                       for q, data in hit_table.groupby('query')}
    query_order = sorted(query_lens, key=lambda q: query_hit_order[q])
    
    # using above orders, get start pos of each contig in genome
    query_starts = get_contig_starts(query_order, query_lens)

    # calculate genome start and end for each hit
    cols = []
    mcols = []
    for mag, starts in [('hit', hit_starts), ('query', query_starts)]:
        start_col = mag + "_start"
        hit_table[start_col] = [starts[m] for m in hit_table[mag]]
        mcols.append(mag)
        cols.append(mag)
        for col in ['start', 'end']:
            mag_col = mag[0] + col
            new_col = mag_col + "_mag"
            cols.append(new_col)
            mcols.append(mag_col)
            hit_table[new_col] = hit_table[mag_col] + hit_table[start_col]

            
    mag_hits = hit_table[cols + other_cols]
    mag_hits.columns = mcols + other_cols
    return mag_hits

def sorted_dot_plot(hit_table, ax=None, colormap='cool', min_pctid=None):
    mag_hits = get_sorted_hits(hit_table, other_cols=['pctid',])
    if isinstance(colormap, str):
        colormap = plt.get_cmap(colormap)
        
    # color connections by pctid
    if min_pctid is None:
        # get the true minimum %ID from data
        min_pctid = min(hits['pctid'])
    elif min_pctid < 0:
        # if negative, go under true min by that much
        min_pctid += min(hits['pctid'])

    def get_line_color(pctid, colormap=colormap, min_pctid=min_pctid):
        return colormap((pctid - min_pctid) / (100 - min_pctid))
    
    if ax is None:
        ax = plt
    
    for i, row in mag_hits.iterrows():
        ax.plot((row.hstart, row.hend), (row.qstart, row.qend), color=get_line_color(row.pctid), alpha=.9)
        
    
        
