"""
Functions modified from the np_read_clustering workflow that plot gene annotations

get_cluster_genes(gene_data, cluster, genome_lengths)

 * pulls out the annotations for the given cluster and flips annotations of any genomes that are predominantly in the reverse strand
 
plot_annot_positions(gene_data, genome_lengths, ax=None, **kwargs):

 * Generates a summary of where the most common annotations are in a set of genomes


"""
from collections import defaultdict, Counter
import re
import numpy
import pandas
from matplotlib import pyplot as plt
from jme.jupy_tools.utils import get_N_colors

BLACK = (0, 0, 0, 1)

def import_genes(faa_file, annot_dict, genome_name_map=None):
    """
    params:
        a faa file from prodigal (with start # end # strand in the header)
        a dict from gene id to annotation tag (PFAM, EggNOG, etc)
        an optional dict renaming contig or genomes
    returns:
        a dataframe of gene annotations suitable for the following functions
    """
    # we'll have to parse the gff or faa file
    gene_rows = []
    with open(faa_file) as faa_lines:
        for line in faa_lines:
            if line.startswith(">"):
                gene, start, end, strand, *_ = [s.strip() for s in line[1:].split("#")]
                genome, gene_no = re.search(r'^(.+)_(\d+)$', gene).groups()
                start, end, strand, gene_no = [int(x) for x in [start, end, strand, gene_no]]
                if genome_name_map is not None:
                    genome = genome_name_map[genome]
                gene_rows.append((gene, genome, gene_no, start, end, strand, annot_dict.get(gene, "")))
    return pandas.DataFrame(gene_rows, columns=['gene','genome', 'gene_no', 'start', 'end', 'strand', 'annot']).set_index('gene')

def get_cluster_genes(gene_data, cluster, genome_lens):
    """
    pulls out the annotations for the given cluster and flips annotations of any genomes that are predominantly in the reverse strand
    
    params:
        gene_data:
            DataFrame with columns: gene(index), genome, start, end, strand, annot
        cluster:
            collection of genome names
        genome_lens:
            dict of genome name to genome length
    """
    # select genes
    cluster_set = set(cluster)
    cluster_genes = gene_data[[g in cluster_set for g in gene_data.genome]].copy()
    
    # flip genome anot strand if needed
    reverse = {genome:(cluster_genes.query(f'genome == "{genome}"') \
                         .eval('glen = strand * (end - start)') \
                         .glen.sum() < 1)
               for genome in cluster}
    cluster_genes[['start','end','strand']] = \
        [(genome_lens[genome] - end, genome_lens[genome] - start, -strand) if reverse[genome] else (start, end, strand)
         for genome, start, end, strand in cluster_genes[['genome', 'start', 'end', 'strand']].values]

    return cluster_genes

def get_shared_genes(gene_data, N=10):
    # get the positions of the named genes
    annot_positions = defaultdict(list)
    for start, end, annot in gene_data[['start', 'end', 'annot']].values:
        if pandas.isna(annot) or len(annot.strip()) == 0:
            continue
        # add mean post to list for this annot
        annot_positions[annot].append((end + start) / 2)
        
    # chose which genes to color
    sorted_genes = sorted(annot_positions.keys(), key=lambda k: len(annot_positions[k]), reverse=True)
    top_N_annots = sorted_genes[:N]
    return top_N_annots

def plot_annot_positions(gene_data, ax=None, gene_color_dict=None, **kwargs):
    """
    Generates a summary of where the most common annotations are in a set of genomes
    
    params:
        gene_data:
            DataFrame of annotations for the selected genomes
            columns: gene(index), genome, start, end, strand, annot
            
        gene_color_dict:
            dict with genes to plot and their colors. Otherwise, use the most common
    """
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=kwargs.get('figsize', [8,6]))

    y_buff = kwargs.get('y_buffer', .5)    

    # get the positions of the named genes
    annot_positions = defaultdict(list)
    for start, end, annot in gene_data[['start', 'end', 'annot']].values:
        if pandas.isna(annot) or len(annot.strip()) == 0:
            continue
        # add mean post to list for this annot
        annot_positions[annot].append((end + start) / 2)

    # chose which genes to color
    if gene_color_dict is None:
        N = kwargs.get('max_colored_genes', 10)
        sorted_genes = sorted(annot_positions.keys(), key=lambda k: len(annot_positions[k]), reverse=True)
        top_N_annots = sorted_genes[:N]
        gene_color_dict = dict(zip(top_N_annots, get_N_colors(N, cmap_name=kwargs.get('gene_cmap', 'Dark2'))))
    else:
        N = len(gene_color_dict)
        n = N
        top_N_annots = sorted(gene_color_dict.keys(),
                              key=lambda k: len(annot_positions[k]),
                              reverse=True)
        sorted_genes = top_N_annots

    # sort ploted genes by median position
    n = kwargs.get('max_plotted_genes', 18)
    sorted_annot = sorted([p for p in sorted_genes if len(annot_positions[p]) > 1][:n], 
                       key=lambda p: numpy.median(list(annot_positions[p])))
    
    # scatter positions, one row per gene
    for i, p in enumerate(sorted_annot):
        x,y = zip(*((gp,i) for gp in annot_positions[p]))
        ax.scatter(x,y, 
                   c=len(y) * [gene_color_dict.get(p, BLACK)], 
                   ec=None, alpha=.5)
    _ = ax.set_yticks(range(len(sorted_annot)))
    ytl = ax.set_yticklabels(sorted_annot)
    for label in ytl:
        label.set_color(gene_color_dict.get(label.get_text(), BLACK))
    _ = ax.set_ylim(-y_buff, i + y_buff)
    
    return top_N_annots, gene_color_dict

def decorate_gene_labels(ax, gene_annots, desc_col):
    # add info to the gene labels
    genome_count = \
        gene_annots.groupby(['genome','annot']).agg({'gene_no':len}).reset_index().groupby('annot').agg({'genome':len}).genome.to_dict()
    if isinstance(desc_col, dict):
        vog_descs=desc_col
    else:
        vog_descs = {}
        for v,d in gene_annots[['annot',desc_col]].values:
            if pandas.notna(v):
                vog_descs[v] = f'{d} ({v}) (in {genome_count.get(v, "unknown")} genomes)'
    new_labels = []
    label_colors = []
    for label in ax.get_yticklabels():
        gene = label.get_text()
        new_labels.append(vog_descs[gene])
        label_colors.append(label.get_color())
        
    # This is a hack, but it works
    _ = ax.set_yticklabels(new_labels)
    for label, color in zip(ax.get_yticklabels(), label_colors):
        label.set_color(color)

def get_first_midpoint(gene_table):
    start, end = gene_table.sort_values('start')[['start','end']].values[0]
    return (start + end) / 2

def get_mean_midpoint(annot, gene_table):
    annot_midpoints = \
        gene_table.query(f"annot == '{annot}'") \
                  .groupby('genome') \
                  .apply(get_first_midpoint)
    if annot_midpoints.shape[0] == 0:
        return None
    return numpy.mean(annot_midpoints)

def align_gene_locs(my_genes, anchor_on_annot=None):
    """ Adjust annotation location in each genome to minimize difference 
    
    # TODO: some genomes need to be flipped, can we do this automatically?
    """

    cluster = set(my_genes.genome.values)
    
    # group by genome and annot to discount duplicate annots in one genome
    most_common_annots = \
        Counter(my_genes.query("annot != ''") \
                        .groupby(['genome', 'annot']) \
                        .agg({'strand': len}) \
                        .reset_index()\
                        .annot).most_common()
    if not anchor_on_annot:
        anchor_on_annot = most_common_annots[0][0]

    annot_midpoints = \
        my_genes.query(f"annot == '{anchor_on_annot}'") \
                .groupby('genome') \
                .apply(get_first_midpoint)

    shifted_genes = None
    shifted_genomes = set()
    for genome, annot_midpoint in annot_midpoints.items():
        shift = -1 * annot_midpoint
        shifted_genome_genes = \
            my_genes.query(f'genome == "{genome}"') \
                    .eval(f'start = start + {shift}') \
                    .eval(f'end = end + {shift}')
        if shifted_genes is not None:
            shifted_genes = shifted_genes.append(shifted_genome_genes)
        else:
            shifted_genes = shifted_genome_genes
        shifted_genomes.add(genome)

    # loop over unshifted genomes (didn't have the first annot)
    shifted_midpoints = {}
    for genome in cluster.difference(shifted_genomes):
        genome_genes = my_genes.query(f'genome == "{genome}"')
        genome_annots = set(genome_genes.annot)

        # get the average shift to allign genes with already shifted genomes
        shifts = []
        for annot, count in most_common_annots:
            if annot in genome_annots:
                if annot in shifted_midpoints:
                    shifted_midpoint = shifted_midpoints[annot]
                else:
                    shifted_midpoint = get_mean_midpoint(annot, shifted_genes)
                    if shifted_midpoint is None:
                        continue
                    shifted_midpoints[annot] = shifted_midpoint

                annot_midpoint = get_first_midpoint(genome_genes.query(f"annot == '{annot}'"))
                shifts.append(shifted_midpoint - annot_midpoint)
        if len(shifts) == 0:
            raise Exception(f"No common annots in genome: {genome}")
            
        shift = numpy.mean(numpy.array(shifts))

        shifted_genome_genes = \
            genome_genes \
                    .eval(f'start = start + {shift}') \
                    .eval(f'end = end + {shift}')
        shifted_genes = shifted_genes.append(shifted_genome_genes)
        #shifted_genomes.add(genome)
        
    return shifted_genes

def draw_genes(gene_data, annot_colors,
               genome_lens=None, 
               ax=None, 
               genome_order=None,
               **kwargs):
    """
    Draws annotations for a set of genomes. If more genomes are given than
    fit in the plot, the top few (by how many of the most common genes they have)
    are ploted
    
    params:
        gene_data:
            DataFrame of annotations for the selected genomes
            columns: gene(index), genome, start, end, strand, annot
        annot_colors:
            dict from annot tag to color. These are the most common genes
        genome_lens:
            dict from genome name to genome length. If given, draw spines 
            behind the annotations to show the genome extent
        max_genomes:
            cap on number of genomes plotted
        genome_order:
            If this is a dict, values are interpreted as Y positions
            If this is a list, use as model order.
            If set to "genes" or any string, the order in the genes table is used.
            Default is to use the number of shared genes.
        
    """
    if ax is None:
        fig, ax = plt.subplots(1,1, figsize=kwargs.get('figsize', [8,6]))

    # chose which genomes to draw
    M = kwargs.get('max_genomes', 20)
    top_M_genomes = gene_data.groupby('genome') \
        .apply(lambda S: len(set(S['annot']) \
        .intersection(annot_colors))) \
        .sort_values(ascending=False) \
        .head(M) \
        .index
    m = len(top_M_genomes)
    
    y_buff = kwargs.get('y_buffer', .5)        # y buffer for bottom plot
    thickness = kwargs.get('thickness', .5)    # arrow thickness
    
    # turn genome order into y positions
    if not isinstance(genome_order, dict):
        # sort top_M if asked (defaults to most shared genes)
        if genome_order is not None:
            # re-order the top_M
            if isinstance(genome_order, str):
                # use the order from the input table
                top_M_ordered = []
                used_genomes = set()
                for g in gene_data.genome.values:
                    if g not in used_genomes:
                        used_genomes.add(g)
                        top_M_ordered .append(g)
            else:
                # get order from list
                top_M_ordered = [g for g in genome_order if g in top_M_genomes]
        else:
            top_M_ordered = top_M_genomes
        genome_y_posns = list(range(m))
        max_y = m - 1
        min_y = 0
    else:
        top_genomes = set(top_M_genomes)
        top_M_ordered, genome_y_posns = \
            zip(*[(g,y) for g, y in genome_order.items() if g in top_genomes])
        min_y = min(genome_y_posns)
        max_y = max(genome_y_posns)

    # calculate the sizes necessary to draw genes using the matplotlib arrow function
    x,y,w,h = ax.bbox.bounds      # w,h are ax size in pixels
    
    y_ax_range = (max_y - min_y) + thickness + (2 * y_buff)
    if genome_lens is None:
        # get range from genes
        min_x = gene_data[['start', 'end']].min().min()
        max_x = gene_data[['start', 'end']].max().max()
        x_ax_range = max_x + 1 - min_x
    else:
        raise Exception("genome spines are not yet implemented")
        
    y_pix_size = y_ax_range / h   # value per pixel
    x_pix_size = x_ax_range / w 

    arrow_half_width = (thickness / 2) / y_pix_size  # y pixels in 1/2 an arrow
    head_length = arrow_half_width * x_pix_size      # value to get the same # of pixels in x axis

    prev_genome = None
    for name, y in zip(top_M_ordered, genome_y_posns):
        ## TODO: plot genome spines
        gene_table = gene_data.query(f'genome == "{name}"')

        # draw genes
        for start, end, strand, annot in gene_table[['start','end','strand','annot']].values:
            strand = int(strand)
            hl = min(head_length, end-start)
            al = max((end - start) - hl, .0001) * strand
            ast = start if al > 0 else end
            color = annot_colors.get(annot, 'k')
            ax.arrow(ast, y, al, 0, fc=color, ec=color, 
                      lw=0,
                      width=thickness, head_width=thickness, 
                      head_length=hl, 
                      head_starts_at_zero=(int(strand) > 0))

    y = ax.set_ylim(min_y-y_buff, max_y + thickness + y_buff)

    ax.set_yticks(genome_y_posns)
    ax.set_yticklabels(top_M_ordered)
    ax.set_xlabel('relative genome position')    
    return top_M_ordered