import re
from collections import namedtuple

CdhitClusters = namedtuple('CdhitClusters', 
                           ['clusters',
                            'reps',
                            'lookup'])

def parse_cdhit_clusters(cluster_file):
    """
    Parses cdhit output into three collections in a named tuple:
      clusters: list of lists of gene ids.
      reps: list of representative gene for each cluster
      lookup: dict mapping from gene names to cluster index

    In this setup, cluster ids are the position in either of the
     first two lists.
    """

    # re-call with file-like object if we are give an path
    if isinstance(cluster_file, str):
        with open(cluster_file) as cluster_handle:
            return parse_cdhit_clusters(cluster_handle)

    # initialize final containers
    clusters = []
    cluster_reps = []
    cluster_lookup = {}
    # expression for parsing cluster line (captures gene name and alignment)
    gene_expr = re.compile(r"\s>(\S+)\.\.\.\s\s*(.+)\s*$")

    # loop over lines
    for line in cluster_file:
        if line.startswith(">"):
            # create a new cluster
            cluster = []
            cluster_id = len(clusters)
            clusters.append(cluster)
            continue
        # parse gene name from line
        gene, alignment = gene_expr.search(line).groups()
        if alignment.strip() == "*":
            cluster_reps.append(gene)
        cluster_lookup[gene] = cluster_id
        cluster.append(gene)

    # done
    return CdhitClusters(clusters, cluster_reps, cluster_lookup)


def parse_multiple_clusters(cluster_files):
    """
    Merge multiple cluster files, must be given in reverse order. 
                                  IE 2D cluster file before individual cluster files.
    """
    cluster_data = None
    for cf in cluster_files:
        cf_data = parse_cdhit_clusters(cf)
        print("read {} clusters from {}".format(len(cf_data[0]), cf))
        if cluster_data is None:
            cluster_data = cf_data
        else:
            merge_clusters(cluster_data, cf_data)
    return cluster_data


def merge_clusters(clusters, new_data):
    for genes, rep in zip(*new_data[:2]):
        # is the rep already part of a cluster?
        try:
            cluster = clusters[2][rep]
        except KeyError as e:
            # no? add cluster to list
            cluster = len(clusters[0])
            clusters[0].append(genes)
            clusters[1].append(rep)
        else:
            # yes? then all genes go in that cluster
            clusters[0][cluster].extend(genes)
            rep = clusters[1][cluster]

        # either way, add entry to rep lookup
        for gene in genes:
            clusters[2][gene] = rep
