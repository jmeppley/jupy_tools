import io

import pytest

from jme.jupy_tools import cdhit


def test_parse_clust():
    temp_file = io.StringIO()
    temp_file.write(TEST_SNIPPET)
    temp_file.seek(0)

    clusters = cdhit.parse_cdhit_clusters(temp_file)

    assert len(clusters) == 3
    assert len(clusters.clusters) == len(clusters.reps)
    assert len(clusters.reps) == 7
    assert len(clusters.lookup) == 43
    assert clusters[0][clusters[2][clusters[1][1]]][3] == 'tp1k5d0.1n5ad455b48'


def test_parse_clust_details():
    clusters = cdhit.parse_cdhit_clusters(io.StringIO(TEST_SNIPPET))
    assert isinstance(clusters, cdhit.CdhitClusters)
    assert [len(c) for c in clusters.clusters] == [15, 4, 1, 1, 5, 1, 16]
    # the rep is the line marked with "*", even if it's not first
    assert clusters.reps[0] == 'tp1k6d0n15ecab2ccd'
    assert clusters.reps[1] == 'tp1k4d0.5n15ad455b48'
    assert clusters.clusters[1][0] == 'tp1k5d0n5ad455b48'
    # every gene maps to the cluster that contains it
    for i, genes in enumerate(clusters.clusters):
        for gene in genes:
            assert clusters.lookup[gene] == i
    for i, rep in enumerate(clusters.reps):
        assert rep in clusters.clusters[i]


def test_parse_clust_from_path(tmp_path):
    path = tmp_path / "test.clstr"
    path.write_text(TEST_SNIPPET)
    clusters = cdhit.parse_cdhit_clusters(str(path))
    assert len(clusters.reps) == 7
    assert len(clusters.lookup) == 43


def test_parse_clust_protein():
    """ protein clusters use "aa" instead of "nt" and report % identity """
    text = (
        ">Cluster 0\n"
        "0\t300aa, >geneA... *\n"
        "1\t280aa, >geneB... at 95.00%\n"
        ">Cluster 1\n"
        "0\t150aa, >geneC|with|pipes... *\n"
    )
    clusters = cdhit.parse_cdhit_clusters(io.StringIO(text))
    assert clusters.clusters == [['geneA', 'geneB'], ['geneC|with|pipes']]
    assert clusters.reps == ['geneA', 'geneC|with|pipes']


# For merging: a primary clustering, and a second file whose clusters either
# share a rep with the first (and get merged) or are new
PRIMARY_CLUSTERS = """>Cluster 0
0\t300aa, >a... *
1\t290aa, >b... at 90.00%
>Cluster 1
0\t200aa, >c... *
"""
SECONDARY_CLUSTERS = """>Cluster 0
0\t300aa, >a... *
1\t250aa, >d... at 80.00%
2\t240aa, >e... at 80.00%
>Cluster 1
0\t100aa, >f... *
1\t100aa, >g... at 99.00%
"""


@pytest.fixture
def cluster_files(tmp_path):
    paths = []
    for name, text in [('primary.clstr', PRIMARY_CLUSTERS),
                       ('secondary.clstr', SECONDARY_CLUSTERS)]:
        path = tmp_path / name
        path.write_text(text)
        paths.append(str(path))
    return paths


def test_merge_clusters():
    merged = cdhit.parse_cdhit_clusters(io.StringIO(PRIMARY_CLUSTERS))
    new = cdhit.parse_cdhit_clusters(io.StringIO(SECONDARY_CLUSTERS))
    cdhit.merge_clusters(merged, new)

    # clusters with a known rep are added to that cluster
    assert {'a', 'b', 'd', 'e'} == set(merged.clusters[0])
    assert merged.clusters[1] == ['c']
    # clusters with a new rep are appended
    assert merged.clusters[2] == ['f', 'g']
    assert merged.reps == ['a', 'c', 'f']
    for gene in ['a', 'b', 'c', 'd', 'e', 'f', 'g']:
        assert gene in merged.lookup


@pytest.mark.xfail(strict=True, reason=(
    "merge_clusters stores rep names in the lookup, but parse_cdhit_clusters "
    "stores cluster indices, so merged lookups are a mix of both"))
def test_merge_clusters_lookup_is_consistent():
    merged = cdhit.parse_cdhit_clusters(io.StringIO(PRIMARY_CLUSTERS))
    new = cdhit.parse_cdhit_clusters(io.StringIO(SECONDARY_CLUSTERS))
    cdhit.merge_clusters(merged, new)
    for i, genes in enumerate(merged.clusters):
        for gene in genes:
            assert merged.lookup[gene] == i


def test_parse_multiple_clusters(cluster_files, capsys):
    merged = cdhit.parse_multiple_clusters(cluster_files)
    assert merged.reps == ['a', 'c', 'f']
    assert len(merged.clusters) == 3
    out = capsys.readouterr().out
    assert "read 2 clusters from " + cluster_files[0] in out
    assert "read 2 clusters from " + cluster_files[1] in out


def test_parse_multiple_clusters_single(cluster_files, capsys):
    merged = cdhit.parse_multiple_clusters(cluster_files[:1])
    assert merged == cdhit.parse_cdhit_clusters(cluster_files[0])


def test_parse_multiple_clusters_empty():
    assert cdhit.parse_multiple_clusters([]) is None


TEST_SNIPPET = """>Cluster 0
0       95993nt, >tp1k6d0n15ecab2ccd... *
1       95993nt, >tp1k4d0n20ecab2ccd... at 1:95993:1:95993/+/100.00%
2       95993nt, >tp1k4d0n15ecab2ccd... at 1:95993:1:95993/+/100.00%
3       95993nt, >tp1k4d0.5n20ecab2ccd... at 1:95993:1:95993/+/100.00%
4       95993nt, >tp1k4d0.1n20ecab2ccd... at 1:95993:1:95993/+/100.00%
5       95993nt, >tp1k6d0.1n10ecab2ccd... at 1:95993:1:95993/+/100.00%
6       95993nt, >tp1k5d0.1n20ecab2ccd... at 1:95993:1:95993/+/100.00%
7       95993nt, >tp1k4d0.1n10ecab2ccd... at 1:95993:1:95993/+/100.00%
8       95993nt, >tp1k5d0n20ecab2ccd... at 1:95993:1:95993/+/100.00%
9       95993nt, >tp1k4d0n5ecab2ccd... at 1:95993:1:95993/+/100.00%
10      95993nt, >tp1k4d0n10ecab2ccd... at 1:95993:1:95993/+/100.00%
11      95993nt, >tp1k6d0.1n5ecab2ccd... at 1:95993:1:95993/+/100.00%
12      95993nt, >tp1k5d0.1n15ecab2ccd... at 1:95993:1:95993/+/100.00%
13      95993nt, >tp1k4d0.1n5ecab2ccd... at 1:95993:1:95993/+/100.00%
14      95993nt, >tp1k5d0n15ecab2ccd... at 1:95993:1:95993/+/100.00%
>Cluster 1
0       95983nt, >tp1k5d0n5ad455b48... at 1:95983:1:95991/+/99.35%
1       95992nt, >tp1k4d0.5n15ad455b48... *
2       95983nt, >tp1k5d0n10ad455b48... at 1:95983:1:95991/+/99.35%
3       95983nt, >tp1k5d0.1n5ad455b48... at 1:95983:1:95991/+/99.35%
>Cluster 2
0       95963nt, >tp1k4d0.5n10ad455b48... *
>Cluster 3
0       91018nt, >tp2k5d0.1n547bda91f... *
>Cluster 4
0       90980nt, >tp2k5d0n547bda91f... *
1       90980nt, >tp2k4d0.1n1047bda91f... at 1:90980:1:90980/+/100.00%
2       90980nt, >tp2k5d0.1n1547bda91f... at 1:90980:1:90980/+/100.00%
3       90980nt, >tp2k6d0n2047bda91f... at 1:90980:1:90980/+/100.00%
4       90980nt, >tp2k5d0.1n1047bda91f... at 1:90980:1:90980/+/100.00%
>Cluster 5
0       86298nt, >tp1k4d0.1n1578d616c9... *
>Cluster 6
0       86298nt, >tp1k6d0n1578d616c9... *
1       86298nt, >tp1k5d0n578d616c9... at 1:86298:1:86298/+/100.00%
2       86298nt, >tp1k4d0n2078d616c9... at 1:86298:1:86298/+/100.00%
3       86298nt, >tp1k4d0n1578d616c9... at 1:86298:1:86298/+/100.00%
4       86298nt, >tp1k4d0.1n2078d616c9... at 1:86298:1:86298/+/100.00%
5       86298nt, >tp1k5d0.1n1078d616c9... at 1:86298:1:86298/+/100.00%
6       86298nt, >tp1k6d0n2078d616c9... at 1:86298:1:86298/+/100.00%
7       86298nt, >tp1k6d0.5n2078d616c9... at 1:86298:1:86298/+/100.00%
8       86298nt, >tp1k6d0.1n1578d616c9... at 1:86298:1:86298/+/100.00%
9       86298nt, >tp1k6d0n578d616c9... at 1:86298:1:86298/+/100.00%
10      86298nt, >tp1k6d0.5n1078d616c9... at 1:86298:1:86298/+/100.00%
11      86298nt, >tp1k5d0.1n2078d616c9... at 1:86298:1:86298/+/100.00%
12      86298nt, >tp1k4d0.1n1078d616c9... at 1:86298:1:86298/+/100.00%
13      86298nt, >tp1k5d0n2078d616c9... at 1:86298:1:86298/+/100.00%
14      86298nt, >tp1k5d0.5n1578d616c9... at 1:86298:1:86298/+/100.00%
15      86298nt, >tp1k4d0n1078d616c9... at 1:86298:1:86298/+/100.00%"""
