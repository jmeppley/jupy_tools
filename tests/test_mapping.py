import shutil

import numpy
import pandas
import pytest

from jme.jupy_tools import mapping
from jme.jupy_tools.mapping import SAMFlag

SAM_TEXT = "\n".join([
    "@HD\tVN:1.6\tSO:unsorted",
    "@SQ\tSN:contig1\tLN:20",
    "@SQ\tSN:contig2\tLN:10",
    # 5M covers 3-7
    "read1\t0\tcontig1\t3\t60\t5M\t*\t0\t0\tACGTA\tIIIII",
    # 2M1I2M2D1M: 2 + 2 + 2 + 1 = 7 ref bases, 5-11 (reverse strand is fine)
    "read2\t16\tcontig1\t5\t60\t2M1I2M2D1M\t*\t0\t0\tACGTAC\tIIIIII",
    # soft clips don't consume the reference: 3 ref bases, 1-3
    "read3\t0\tcontig2\t1\t60\t2S3M\t*\t0\t0\tACGTA\tIIIII",
    # secondary alignment: skipped by default
    "read3\t256\tcontig1\t15\t0\t3M\t*\t0\t0\t*\t*",
    # supplementary alignment: skipped by default
    "read1\t2048\tcontig2\t8\t0\t3M\t*\t0\t0\t*\t*",
    # unmapped
    "read4\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
]) + "\n"


@pytest.mark.parametrize("cigar, length", [
    ("10M", 10),
    ("5S10M5S", 10),
    ("3M2I3M", 6),
    ("3M2D3M", 8),
    ("3M100N3M", 106),
    ("4=1X4=", 9),
    ("5H3M", 3),
])
def test_get_cigar_ref_length(cigar, length):
    assert mapping.get_cigar_ref_length(cigar) == length


def test_get_cigar_ref_length_empty():
    with pytest.raises(ValueError):
        mapping.get_cigar_ref_length("*")


def test_parse_sam_alignments():
    alignments = list(mapping.parse_sam_alignments(SAM_TEXT.splitlines(True)))
    assert alignments == [
        ("read1", "contig1", 3, 7),
        ("read2", "contig1", 5, 11),
        ("read3", "contig2", 1, 3),
    ]


def test_parse_sam_alignments_flags():
    # keep everything that is mapped
    alignments = list(mapping.parse_sam_alignments(
        SAM_TEXT.splitlines(), exclude_flags=SAMFlag.UNMAPPED,
    ))
    assert len(alignments) == 5
    assert ("read3", "contig1", 15, 17) in alignments
    assert ("read1", "contig2", 8, 10) in alignments

    # drop reverse strand alignments too
    alignments = list(mapping.parse_sam_alignments(
        SAM_TEXT.splitlines(),
        exclude_flags=mapping.DEFAULT_EXCLUDE_FLAGS | SAMFlag.REV_COMP,
    ))
    assert [a[0] for a in alignments] == ["read1", "read3"]


def test_get_deltas():
    deltas = mapping.get_deltas([
        ("r1", "c1", 3, 5),
        ("r2", "c1", 8, 4),  # reversed coordinates are fine
        ("r3", "c2", 1, 1),
    ])
    assert deltas == {
        "c1": {3: 1, 6: -1, 4: 1, 9: -1},
        "c2": {1: 1, 2: -1},
    }
    # plain dicts, so missing references raise KeyError
    assert type(deltas["c1"]) is dict


def test_get_deltas_portion():
    deltas = mapping.get_deltas(
        [("r1", "c1", 1, 2), ("r1", "c2", 1, 2), ("r2", "c1", 1, 2)],
        portion=True,
    )
    assert deltas["c1"][1] == pytest.approx(1.5)
    assert deltas["c2"][1] == pytest.approx(0.5)


def test_get_deltas_empty():
    assert mapping.get_deltas([]) == {}


def test_deltas_to_cov():
    deltas = {3: 1, 6: -1, 4: 1, 9: -1}
    coverage = mapping.deltas_to_cov(deltas, nan_for_zero=False)
    numpy.testing.assert_array_equal(coverage, [0, 0, 1, 2, 2, 1, 1, 1])

    coverage = mapping.deltas_to_cov(deltas)
    numpy.testing.assert_array_equal(
        coverage, [numpy.nan, numpy.nan, 1, 2, 2, 1, 1, 1]
    )


def test_deltas_to_cov_x_max():
    deltas = {3: 1, 6: -1, 4: 1, 9: -1}
    # padded to the reference length
    coverage = mapping.deltas_to_cov(deltas, x_max=10, nan_for_zero=False)
    numpy.testing.assert_array_equal(coverage, [0, 0, 1, 2, 2, 1, 1, 1, 0, 0])

    # truncated to the reference length
    coverage = mapping.deltas_to_cov(deltas, x_max=4, nan_for_zero=False)
    numpy.testing.assert_array_equal(coverage, [0, 0, 1, 2])

    # alignment running right up to the end
    coverage = mapping.deltas_to_cov({1: 1, 6: -1}, x_max=5)
    numpy.testing.assert_array_equal(coverage, [1, 1, 1, 1, 1])


def test_deltas_to_cov_empty():
    assert len(mapping.deltas_to_cov({})) == 0
    numpy.testing.assert_array_equal(
        mapping.deltas_to_cov({}, x_max=3, nan_for_zero=False), [0, 0, 0]
    )


def test_deltas_to_cov_portion_rounding():
    """ fractional coverage that sums back to zero should be treated as zero """
    deltas = mapping.get_deltas(
        [("r1", "c1", 1, 3), ("r1", "c1", 1, 3), ("r1", "c1", 1, 3)],
        portion=True,
    )
    coverage = mapping.deltas_to_cov(deltas["c1"], x_max=5)
    numpy.testing.assert_allclose(coverage[:3], [1, 1, 1])
    assert numpy.isnan(coverage[3:]).all()


def test_parse_deltas_from_sam(tmp_path):
    sam_file = tmp_path / "test.sam"
    sam_file.write_text(SAM_TEXT)

    deltas = mapping.parse_deltas_from_sam(sam_file)
    assert set(deltas) == {"contig1", "contig2"}
    coverage = mapping.deltas_to_cov(deltas["contig1"], x_max=20, nan_for_zero=False)
    expected = numpy.zeros(20)
    expected[2:7] += 1   # read1: 3-7
    expected[4:11] += 1  # read2: 5-11
    numpy.testing.assert_array_equal(coverage, expected)


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")
def test_parse_deltas_from_sam_with_samtools(tmp_path):
    sam_file = tmp_path / "test.sam"
    sam_file.write_text(SAM_TEXT)
    with_samtools = mapping.parse_deltas_from_sam(sam_file, use_samtools=True)
    assert with_samtools == mapping.parse_deltas_from_sam(sam_file)


def test_samtools_failure(tmp_path, monkeypatch):
    """ a failing samtools command should raise an error, not return nothing """
    monkeypatch.setattr(mapping, "SAMTOOLS_COMMAND", "false")
    with pytest.raises(Exception):
        mapping.parse_deltas_from_sam(tmp_path / "missing.bam")


def test_get_deltas_from_hit_table():
    hits = pandas.DataFrame({
        "query": ["r1", "r2", "r2"],
        "hit": ["c1", "c1", "c2"],
        "hstart": [3, 8, 1],
        "hend": [5, 4, 2],
    })
    assert mapping.get_deltas_from_hit_table(hits) == {
        "c1": {3: 1, 6: -1, 4: 1, 9: -1},
        "c2": {1: 1, 3: -1},
    }

    deltas = mapping.get_deltas_from_hit_table(hits, portion=True)
    assert deltas["c2"] == {1: 0.5, 3: -0.5}


def test_get_deltas_from_hit_table_zero_indexed():
    """ PAF/BED style [start, end) coordinates """
    hits = pandas.DataFrame({
        "read": ["r1"], "contig": ["c1"], "start": [2], "end": [5],
    })
    deltas = mapping.get_deltas_from_hit_table(
        hits, read_col="read", ref_col="contig", start_col="start",
        end_col="end", zero_indexed=True,
    )
    # 0-based [2, 5) is 1-based 3-5
    assert deltas == {"c1": {3: 1, 6: -1}}


def test_samflag():
    flag = SAMFlag(16 + 64 + 1)
    assert SAMFlag.REV_COMP in flag
    assert SAMFlag.READ_1 in flag
    assert SAMFlag.PAIRED in flag
    assert SAMFlag.READ_2 not in flag
    assert int(mapping.DEFAULT_EXCLUDE_FLAGS) == 2308
