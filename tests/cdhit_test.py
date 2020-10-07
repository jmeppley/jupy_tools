from jme.jupy_tools import cdhit

def test_parse_clust():
    import io
    temp_file = io.StringIO()
    temp_file.write(TEST_SNIPPET)
    temp_file.seek(0)

    clusters = cdhit.parse_cdhit_clusters(temp_file)
    
    assert len(clusters) == 3
    assert len(clusters.clusters) == len(clusters.reps)
    assert len(clusters.reps) == 7
    assert len(clusters.lookup) == 43
    assert clusters[0][clusters[2][clusters[1][1]]][3] == 'tp1k5d0.1n5ad455b48'

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
