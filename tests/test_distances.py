from itertools import combinations

import pandas
import pytest

from jme.jupy_tools import distances


def product_of_values(item1, item2):
    """ a picklable pair function for (key, value) items """
    return item1[1] * item2[1]


def identity(item):
    """ a picklable key function """
    return item


@pytest.mark.parametrize("data_len", [0, 1, 2, 3, 7, 10, 23, 100])
@pytest.mark.parametrize("n", [1, 2, 3, 5])
def test_chunkify_data_covers_all_pairs_once(data_len, n):
    data = list(range(data_len))
    pairs = []
    for slice_a, slice_b in distances.chunkify_data(data, n=n):
        assert len(slice_a) > 0
        assert slice_b is None or len(slice_b) > 0
        chunk_results = distances.process_chunk(
            (slice_a, slice_b), lambda a, b: None, get_key_fn=identity
        )
        pairs.extend(tuple(sorted((key1, key2))) for key1, key2, _ in chunk_results)

    assert len(pairs) == len(set(pairs))
    assert set(pairs) == set(combinations(data, 2))


def test_chunkify_data_chunk_count():
    chunks = list(distances.chunkify_data(list(range(100)), n=5))
    # 5 self chunks + (2 * 4 + 2 * 3 + 2 * 2 + 2 * 1) others
    assert len(chunks) == 25
    assert sum(1 for _, slice_b in chunks if slice_b is None) == 5


def test_process_chunk():
    data = [("a", 1), ("b", 2), ("c", 3)]
    results = distances.process_chunk((data, None), product_of_values)
    assert results == [("a", "b", 2), ("a", "c", 3), ("b", "c", 6)]

    results = distances.process_chunk((data[:1], data[1:]), product_of_values)
    assert results == [("a", "b", 2), ("a", "c", 3)]


def _expected_products(data, mirror=True):
    """ calculate the expected results serially """
    expected = {}
    for (key1, val1), (key2, val2) in combinations(data, 2):
        expected[(key1, key2)] = val1 * val2
        if mirror:
            expected[(key2, key1)] = val1 * val2
    return expected


@pytest.mark.parametrize("threads", [1, 2])
def test_multithreaded_distances(threads):
    data = [(f"item{i}", i) for i in range(20)]
    results = distances.multithreaded_distances(
        data, product_of_values, chunk_root=3, threads=threads
    )
    assert results == _expected_products(data)


def test_multithreaded_distances_no_mirror():
    data = [(f"item{i}", i) for i in range(10)]
    results = distances.multithreaded_distances(
        data, product_of_values, threads=1, mirror=False
    )
    assert len(results) == 45
    # every pair is in the results in one direction or the other
    expected = _expected_products(data)
    for (key1, key2), score in results.items():
        assert expected[(key1, key2)] == score


def test_multithreaded_distances_default_threads(monkeypatch):
    """ the default thread count should be capped """
    requested = {}

    class FakePool:
        """ records the number of processes requested and runs serially """

        def __init__(self, processes):
            """ save the process count """
            requested["processes"] = processes

        def __enter__(self):
            """ return self as the pool """
            return self

        def __exit__(self, *args):
            """ nothing to clean up """
            return False

        def imap_unordered(self, function, items):
            """ run in this process """
            return map(function, items)

    monkeypatch.setattr(distances, "Pool", FakePool)
    data = [(f"item{i}", i) for i in range(5)]
    results = distances.multithreaded_distances(data, product_of_values, threads=None)
    assert results == _expected_products(data)
    assert 1 <= requested["processes"] <= distances.MAX_DEFAULT_THREADS


def test_get_default_thread_count():
    assert 1 <= distances.get_default_thread_count() <= distances.MAX_DEFAULT_THREADS
    assert distances.get_default_thread_count(max_threads=1) == 1


GENES = pandas.DataFrame(
    [
        # genome, annot, start, end
        ("g1", "A", 1, 100),
        ("g1", "B", 101, 200),
        ("g1", "A", 201, 250),
        ("g2", "A", 1, 200),
        ("g2", "C", 201, 300),
        ("g3", "Unknown", 1, 100),
        ("g3", "D", 101, 200),
        ("g4", "Unknown", 1, 50),
        ("g4", "D", 51, 100),
    ],
    columns=["genome", "annot", "start", "end"],
)


def test_gene_annots_to_data_array():
    data = distances.gene_annots_to_data_array(GENES)
    assert data == [
        ("g1", {"A": 150, "B": 100}, 250),
        ("g2", {"A": 200, "C": 100}, 300),
        ("g3", {"Unknown": 100, "D": 100}, 200),
        ("g4", {"Unknown": 50, "D": 50}, 100),
    ]


def test_gene_annots_to_data_array_len_col():
    genes = GENES.assign(length=GENES.end + 1 - GENES.start) \
                 .drop(columns=["start", "end"]) \
                 .rename(columns={"genome": "gid", "annot": "og"})
    data = distances.gene_annots_to_data_array(
        genes, annot_col="og", genome_col="gid", gene_len_col="length"
    )
    assert data == distances.gene_annots_to_data_array(GENES)
    # input is not modified
    assert list(genes.columns) == ["gid", "og", "length"]


def test_calc_shared_gene_length_ratio():
    data = dict((d[0], d) for d in distances.gene_annots_to_data_array(GENES))
    ratio = distances.calc_shared_gene_length_ratio
    # g1 and g2 share A: (150 + 200) / (250 + 300)
    assert ratio(data["g1"], data["g2"]) == pytest.approx(350 / 550)
    # nothing shared
    assert ratio(data["g1"], data["g3"]) == 0
    # Unknown is not counted as shared by default
    assert ratio(data["g3"], data["g4"]) == pytest.approx(150 / 300)
    assert ratio(data["g3"], data["g4"], ignore_annots=()) == 1
    # empty genomes
    assert ratio(("x", {}, 0), ("y", {}, 0)) == 0


def test_calculate_shared_gene_length_distances():
    results = distances.calculate_shared_gene_length_distances(GENES, threads=1)
    assert len(results) == 12
    assert results[("g1", "g2")] == pytest.approx(350 / 550)
    assert results[("g2", "g1")] == pytest.approx(350 / 550)
    assert results[("g4", "g3")] == pytest.approx(0.5)

    renamed = GENES.rename(columns={"annot": "og"})
    assert results == distances.calculate_shared_gene_length_distances(
        renamed, data_kws={"annot_col": "og"}, threads=1
    )
