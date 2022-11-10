"""
Plot read hit positions over reference sequences
"""
import re, numpy
from matplotlib.pyplot import cm
from edl import blastm8

def remove_overlapping_hits(hits, on_hit=True, buffer=0):
    regions = []
    for hit in hits:
        start, end = (hit.hstart, hit.hend) \
                     if on_hit else \
                     sorted((hit.qstart, hit.qend))
        
        # if hit is shorter than buffer, then keep it
        if end - start <= buffer:
            regions.append((start, end))
            yield hit
            continue
        
        # adjust for buffer
        buf_start, buf_end = start + buffer, end - buffer

        # check overlap for all ranges
        for i, already_occ_range in enumerate(regions):
            if (buf_start >= already_occ_range[1] 
                or buf_end <= already_occ_range[0]):                    
                # does not overlap this range (try next range) 
                continue
            else:
                # overlaps a previous hit...move on
                break
        else:
            # there was no overlap, save range and yield this hit
            regions.append((start, end))
            yield hit


def parse_contig_hits(map_file, format=blastm8.SAM, **kwargs):
    """
    Parse a hit table (default format is SAM) but group
    hits by contig then read in nested dict
    
    kwargs passed to blastm8.FilterParams
    """
    contig_hits = {}
    for read, hits in blastm8.generate_hits(map_file, format=format, **kwargs):
        read_contig_hits = {}
        for hit in hits:
            contig = hit.hit
            read_contig_hits.setdefault(contig, []).append(hit)
        for contig, chits in read_contig_hits.items():
            contig_hits.setdefault(contig, {})[read] = chits
    return contig_hits

def contig_hits_from_table(hits, contig_col='hit', read_col='query'):
    """
    A kludge to get tabular hits (from util.parse_blast_m8) into
    the proper format for the plot functions. Only tested with PAF
    """
    hit_dict = {} 
    for contig, c_hits in hits.groupby(contig_col):
        contig_hits = hit_dict.setdefault(contig, {})
        for query, q_hits in c_hits.groupby(read_col):
            hit_list = contig_hits.setdefault(query, [])
            for i, row in q_hits.iterrows():
                hit_list.append(row)
    return hit_dict

def hits_intersect(h1, h2):
    return h1.hstart <= h2.hend and h2.hstart <= h1.hend


def ranges_intersect(r1, r2):
    return r1[0] <= r2[1] and r2[0] <= r1[1]


def hit_list_intersects_range(hits, r):
    for hit in hits:
        if hit.hstart <= r[1] and r[0] <= hit.hend:
            return True
    else:
        return False


def get_y_pos(q_range, range_list):
    """ new version flips the script:
            ranges are now grouped by row
            we search rows until we find an empty one
        NB: this could be sped up by sorting ranges and dong a binary search!
    """
    for pos, ranges in enumerate(range_list):
        for r in ranges:
            if ranges_intersect(q_range, r):
                break
        else:
            # if we're still here, return row
            ranges.append(q_range)
            return pos

        # go to the next row if we collided
        continue

    # if we get out entirely, add a new row
    range_list.append(
        [q_range,]
    )
    return len(range_list) - 1


def get_hit_range(hits):
    start = None
    end = None
    for hit in hits:
        if start is None:
            start, end = sorted([hit.hstart, hit.hend])
        else:
            hs, he = sorted([hit.hstart, hit.hend])
            start = min(start, hs)
            end = max(end, he)
    return (start, end)


def plot_hits(
    ax,
    hits_dict,
    reference_length=None,
    plot_range=None,
    hide_missing_pairs=True,
    bar_height=0.8,
    bar_color="blue",
    sort_key=lambda i: i,
):
    """ Groups hits by query, and places above reference
        hits_dict: map from query_id to list of hit objects
        bar_color (def: "blue"): a color or a function on hit returning a color
    """
    if plot_range is None:
        if reference_length is None:
            raise Exception("please provide reference length or plot range!")

        ax.barh(0, reference_length, bar_height, 0, color="black")
        keys = sorted(list(hits_dict.keys()), key=sort_key)
    else:
        ax.barh(
            0,
            plot_range[1] + 1 - plot_range[0],
            bar_height,
            plot_range[0],
            color="black",
        )
        keys = sorted(
            list(
                [
                    k
                    for k, h in hits_dict.items()
                    if hit_list_intersects_range(h, plot_range)
                ]
            ),
            key=sort_key,
        )

    # set up bar color
    if callable(bar_color):
        get_color = bar_color
    else:

        def get_color(hit):
            return bar_color

    ranges = []
    key_positions = {}
    for i, key in enumerate(keys):
        if plot_range and hide_missing_pairs:
            hits = [
                h
                for h in hits_dict[key]
                if ranges_intersect((h.hstart, h.hend), plot_range)
            ]
        else:
            hits = hits_dict[key]
        q_range = get_hit_range(hits)
        pos = get_y_pos(q_range, ranges) + 1
        ax.barh(
            pos, q_range[1] + 1 - q_range[0], bar_height / 3, q_range[0], color="red"
        )
        for hit in hits:
            hit_color = get_color(hit)
            ax.barh(
                pos, hit.hend - hit.hstart + 1, bar_height, hit.hstart, color=hit_color,
            )

        key_positions[key] = pos

    ax.set_yticks(sorted(set(key_positions.values())))
    return key_positions


def plot_hits_colored_py_pos(
    ax,
    hits_dict,
    reference_length=None,
    bar_height=0.8,
    plot_range=None,
    hide_missing_pairs=True,
    bar_cmap=cm.gist_rainbow,
    norm_read_color=False,
    read_lengths=None,
):
    """ Groups hits by query, and places above reference
        colors hit bars by read position,
        hits_dict: map from query_id to list of hit objects
        norm_read_color: False (default): same color scale for all reads.
                         True: max is read end, not longest read.
                         dict: look up read length in dict as max per read
    """
    if plot_range is None:
        if reference_length is None:
            raise Exception("please provide reference length or plot range!")
        x0 = 0
        xw = reference_length
        keys = sorted(list(hits_dict.keys()))
    else:
        x0 = plot_range[0] - 1
        xw = plot_range[1] + 1 - plot_range[0]
        keys = sorted(
            list(
                [
                    k
                    for k, h in hits_dict.items()
                    if hit_list_intersects_range(h, plot_range)
                ]
            )
        )
    ax.barh(0, xw, bar_height, x0, color="black")

    # find longest read (largest of all qstarts and qends)
    if norm_read_color is False:
        read_max = 0
        for k in keys:
            for hit in hits_dict[k]:
                read_max = max(read_max, max(hit.qstart, hit.qend))
        # dummy for imshow
        X = [[0.6, 0.6], [0.7, 0.7]]

    ranges = []
    key_positions = {}
    for i, key in enumerate(keys):
        if plot_range and hide_missing_pairs:
            hits = [
                h
                for h in hits_dict[key]
                if ranges_intersect((h.hstart, h.hend), plot_range)
            ]
        else:
            hits = hits_dict[key]
        q_range = get_hit_range(hits)
        pos = get_y_pos(q_range, ranges) + 1
        ax.barh(
            pos - bar_height / 2,
            q_range[1] + 1 - q_range[0],
            bar_height / 4,
            q_range[0],
            color="grey",
        )
        if norm_read_color is not False:
            if isinstance(norm_read_color, dict):
                read_max = norm_read_color[next(iter(hits)).read]
            else:
                read_max = 0
                for hit in hits:
                    read_max = max(read_max, max(hit.qstart, hit.qend))
        for hit in hits:
            ax.imshow(
                2 * [[hit.qstart, hit.qend]],
                interpolation="bicubic",
                aspect="auto",
                cmap=bar_cmap,
                vmin=0,
                vmax=read_max,
                extent=[
                    hit.hstart,
                    hit.hend,
                    pos - bar_height / 2,
                    pos + bar_height / 2,
                ],
            )

        key_positions[key] = pos

    ax.set_xlim(x0, x0 + xw)
    ax.set_ylim(-0.5, max(key_positions.values()) + 0.5)
    ax.set_yticks(sorted(set(key_positions.values())))
    return key_positions

cigar_segment_rexp = re.compile(r'(\d+)([A-Z])')

def get_diffs(cigar_string, query_seq, hit_start, hit_end, hit_seq, debug=False):
    """ generates a tuple for each position in the alignment where the sequences differ 
            (query_position, query_base, hit_position, hit_base)
        For insertions/deletions, the missing base will be an empty string
        
        NB: diff postions are 0 indexed!
    """
    qpos = 0
    hpos = hit_start - 1
    for match in cigar_segment_rexp.finditer(cigar_string):
        base_count, segment_type = match.groups()
        base_count = int(base_count)
        if segment_type == "H":
            pass
        elif segment_type == "S":
            qpos += base_count
        elif segment_type == "I":
            for i in range(base_count):
                qbase = query_seq[qpos]
                yield (qpos, qbase, hpos, "")
                qpos += 1
        elif segment_type == "D":
            for i in range(base_count):
                hbase = hit_seq[hpos]
                yield (qpos, "", hpos, hbase)
                hpos += 1
        elif segment_type == "M":
            for i in range(base_count):
                qbase = query_seq[qpos]
                hbase = hit_seq[hpos]
                if qbase != hbase:
                    yield (qpos, qbase, hpos, hbase)
                qpos += 1
                hpos += 1
        else:
            raise Exception("Unknown cigar character: " + segment_type)


def get_hit_pctids(hit_start, hit_end, diffs, N):
    diff_iter = iter(diffs)
    try:
        next_diff = next(diff_iter)
    except StopIteration:
        next_diff = None
    mlen = hit_end + 1 - hit_start
    region_length = mlen / N  # NOT a integer!
    for hpos in numpy.arange(hit_start - 1, hit_end, region_length):
        diff_count = 0
        while next_diff is not None and next_diff[2] < hpos:
            diff_count += 1
            try:
                next_diff = next(diff_iter)
            except StopIteration:
                next_diff = None
        yield 100 * (region_length - diff_count) / region_length


def get_dot_count(span, axrange, ax, axis="x"):
    # dpi times the figure size
    fig_dots = (
        ax.get_figure().dpi * ax.get_figure().get_size_inches()[0 if axis == "x" else 1]
    )

    # ax postion in given in fraction of fig, so just multiply
    ax_pos = ax.get_position()
    if axis == "x":
        ax_dots = fig_dots * (ax_pos.x1 - ax_pos.x0)
    else:
        ax_dots = fig_dots * (ax_pos.y1 - ax_pos.y0)

    # ratio of span in question to full axis range
    return ax_dots * (span[1] - span[0]) / (axrange[1] - axrange[0])


def get_dot_count_y(span, xrange, ax):
    return get_dot_count(span, xrange, ax, axis="y")


def get_dot_count_x(span, xrange, ax):
    return get_dot_count(span, xrange, ax, axis="x")


def plot_hits_colored_pctid(
    ax,
    hits_dict,
    reference_sequence,
    bar_height=0.8,
    plot_range=None,
    hide_missing_pairs=True,
    hit_cmap=cm.gray,
    pctid_range=[0, 100],
):
    """ Groups hits by query, and places above reference
        colors hit bars via custom callable:
            takes hit as argumetns and returns
                hit: edl.blastm8.Hit object
                pos: y position
        hits_dict: map from query_id to list of hit objects
        norm_read_color: False (default): same color scale for all reads. True: max is read end, not longest read.
    """
    if plot_range is None:
        x0 = 0
        xw = len(reference_sequence)
        keys = sorted(list(hits_dict.keys()))
    else:
        x0 = plot_range[0] - 1
        xw = plot_range[1] + 1 - plot_range[0]
        keys = sorted(
            list(
                [
                    k
                    for k, h in hits_dict.items()
                    if hit_list_intersects_range(h, plot_range)
                ]
            )
        )
    x_range = (x0, x0 + xw)
    ax.barh(0, xw, bar_height, x0, color="black")

    ranges = []
    key_positions = {}
    for i, key in enumerate(keys):
        if plot_range and hide_missing_pairs:
            hits = [
                h
                for h in hits_dict[key]
                if ranges_intersect((h.hstart, h.hend), plot_range)
            ]
        else:
            hits = hits_dict[key]
        q_range = get_hit_range(hits)
        pos = get_y_pos(q_range, ranges) + 1
        ax.barh(
            pos - bar_height / 2,
            q_range[1] + 1 - q_range[0],
            bar_height / 4,
            q_range[0],
            color="grey",
        )

        # scan hits for hit with read seq
        for hit in hits:
            hit_seq = hit.line.split("\t")[9]
            if hit_seq != "*":
                read_seq = hit_seq

        read_seq = None
        for hit in hits:
            hit_seq = hit.line.split("\t")[9]
            if read_seq is None or len(hit_seq) > len(read_seq):
                read_seq = hit_seq
            elif hit_seq == "*":
                hit_seq = read_seq
            # get diffeernces from alignment
            cigar_string = hit.line.split("\t")[5]
            diffs = get_diffs(
                cigar_string, hit_seq, hit.hstart, hit.hend, reference_sequence
            )
            # how many dots will this hit span in the image
            dots = get_dot_count_x((hit.hstart, hit.hend), x_range, ax)
            # calc the pctid for each dot
            pctids = get_hit_pctids(hit.hstart, hit.hend, diffs, dots)
            # make into a rectangle
            ax.imshow(
                2 * [list(pctids),],
                interpolation="bicubic",
                aspect="auto",
                cmap=hit_cmap,
                vmin=pctid_range[0],
                vmax=pctid_range[1],
                extent=[
                    hit.hstart,
                    hit.hend,
                    pos - bar_height / 2,
                    pos + bar_height / 2,
                ],
            )

        key_positions[key] = pos

    ax.set_xlim(*x_range)
    ax.set_ylim(-0.5, max(key_positions.values()) + 0.5)
    ax.set_yticks(sorted(set(key_positions.values())))
    return key_positions
