"""
Turn read mapping results (SAM/BAM files or hit tables) into per-base coverage.

Coverage is built in two steps. First, alignments are collected into
"coverage deltas": for each reference sequence, a dict of
{position: change in coverage at that position}. Then deltas_to_cov() turns
one reference's deltas into a numpy array of coverage by position.

    from jme.jupy_tools import mapping
    deltas_by_ref = mapping.parse_deltas_from_sam("reads.v.contigs.bam")
    coverage = mapping.deltas_to_cov(deltas_by_ref["contig_1"], x_max=contig_len)

Positions are 1-based throughout: coverage[0] is the coverage of position 1.

Functions:
    parse_deltas_from_sam: SAM/BAM file -> deltas by reference
    get_deltas_from_hit_table: hit table DataFrame (see hit_tables) -> deltas
    get_deltas: (read, reference, start, end) tuples -> deltas
    deltas_to_cov: deltas for one reference -> coverage array
    parse_sam_alignments: SAM lines -> (read, reference, start, end) tuples
    SAMFlag: decompose SAM flags into named bits
"""
import enum
import logging
import re
import subprocess
from collections import defaultdict

import numpy

LOGGER = logging.getLogger(__name__)


class SAMFlag(enum.IntFlag):
    """
    Bits of the SAM FLAG field. Combine and test them like any IntFlag:

        >>> SAMFlag.REV_COMP in SAMFlag(16 + 256)
        True

    From the SAM format specification:
    (int)   (meaning)
    1       template having multiple segments in sequencing (read is paired)
    2       each segment properly aligned according to the aligner (proper pair)
    4       segment unmapped
    8       next segment in the template unmapped (mate unmapped)
    16      SEQ being reverse complemented
    32      SEQ of the next segment in the template being reverse complemented
    64      the first segment in the template (read 1)
    128     the last segment in the template (read 2)
    256     secondary alignment (not primary)
    512     not passing filters, such as platform/vendor quality controls
    1024    PCR or optical duplicate
    2048    supplementary alignment (eg a portion of a split read)
    """

    PAIRED = 1
    PROPER_PAIR = 2
    UNMAPPED = 4
    NEXT_UNMAPPED = 8
    REV_COMP = 16
    NEXT_REV_COMP = 32
    READ_1 = 64
    READ_2 = 128
    NON_PRIMARY = 256
    LOW_Q = 512
    DUPLICATE = 1024
    SUPPLEMENTAL = 2048


# By default, skip alignments that aren't mapped, primary, or full alignments.
# (Same as "samtools view -F 2308")
DEFAULT_EXCLUDE_FLAGS = SAMFlag.UNMAPPED | SAMFlag.NON_PRIMARY | SAMFlag.SUPPLEMENTAL

SAMTOOLS_COMMAND = "samtools"

# CIGAR operations that consume reference bases
_REF_CONSUMING_OPS = set("MDN=X")
_CIGAR_REXP = re.compile(r"(\d+)([MIDNSHP=X])")


def get_cigar_ref_length(cigar):
    """
    Return the number of reference bases covered by a CIGAR string
    (the total of the M, D, N, =, and X operations).
    """
    if cigar == "*":
        raise ValueError("Cannot get alignment length from an empty CIGAR ('*')")
    return sum(
        int(length)
        for length, operation in _CIGAR_REXP.findall(cigar)
        if operation in _REF_CONSUMING_OPS
    )


def parse_sam_alignments(sam_lines, exclude_flags=DEFAULT_EXCLUDE_FLAGS):
    """
    Generate (read, reference, start, end) tuples from SAM formatted lines.

    Parameters:
        sam_lines: iterable of SAM lines (str). Header lines are skipped.
        exclude_flags: skip alignments with any of these flag bits set
            (default: unmapped, secondary, and supplementary alignments)

    Positions are 1-based and inclusive. The end comes from the CIGAR string.
    """
    exclude_flags = int(exclude_flags)
    for line in sam_lines:
        if line.startswith("@") or not line.strip():
            continue
        fields = line.split("\t", 6)
        read, flag, reference, position, _mapq, cigar = fields[:6]
        if int(flag) & exclude_flags:
            continue
        if reference == "*" or cigar == "*":
            # unmapped, even though the flag didn't say so
            continue
        start = int(position)
        end = start + get_cigar_ref_length(cigar) - 1
        yield read, reference, start, end


def _open_sam_lines(sam_file, use_samtools):
    """
    Return a context manager yielding the lines of a SAM or BAM file.

    BAM and CRAM files (and any file if use_samtools is True) are read
    through "samtools view -h". Other files are opened as text.
    """
    if use_samtools is None:
        use_samtools = str(sam_file).lower().endswith((".bam", ".cram"))

    if not use_samtools:
        return open(sam_file, "rt")

    return _SamtoolsLines(sam_file)


class _SamtoolsLines:
    """ Context manager for the text output of 'samtools view -h' """

    def __init__(self, sam_file):
        """ Save the file name for samtools to read """
        self.sam_file = str(sam_file)
        self.process = None

    def __enter__(self):
        """ Start samtools and return its stdout as a stream of str lines """
        self.process = subprocess.Popen(
            [SAMTOOLS_COMMAND, "view", "-h", self.sam_file],
            stdout=subprocess.PIPE,
            text=True,
        )
        return self.process.stdout

    def __exit__(self, exc_type, exc_value, traceback):
        """ Close the pipe, wait for samtools, and check its exit status """
        self.process.stdout.close()
        return_code = self.process.wait()
        if exc_type is None and return_code != 0:
            raise subprocess.CalledProcessError(
                return_code, [SAMTOOLS_COMMAND, "view", "-h", self.sam_file]
            )
        return False


def parse_deltas_from_sam(
    sam_file,
    exclude_flags=DEFAULT_EXCLUDE_FLAGS,
    use_samtools=None,
    portion=False,
):
    """
    Parse a SAM or BAM file into coverage deltas by reference sequence.

    Parameters:
        sam_file: path to a SAM or BAM file
        exclude_flags: skip alignments with any of these SAM flag bits set
            (default: unmapped, secondary, and supplementary alignments)
        use_samtools: if True, read the file with samtools. None (default):
            use samtools only for .bam and .cram files.
        portion: see get_deltas()

    Returns:
        see get_deltas()
    """
    with _open_sam_lines(sam_file, use_samtools) as sam_lines:
        return get_deltas(
            parse_sam_alignments(sam_lines, exclude_flags=exclude_flags),
            portion=portion,
        )


def get_deltas_from_hit_table(
    hit_table,
    read_col="query",
    ref_col="hit",
    start_col="hstart",
    end_col="hend",
    zero_indexed=False,
    portion=False,
):
    """
    Get coverage deltas from a pandas DataFrame of hits (eg from
    hit_tables.parse_blast_m8).

    Parameters:
        hit_table: pandas DataFrame with one row per alignment
        read_col, ref_col, start_col, end_col: names of the columns with the
            read name, reference name, and alignment start and end positions
            on the reference. Start and end can be in either order.
        zero_indexed: set to True if positions are 0-based half-open (BED/PAF
            style). By default positions are 1-based and inclusive.
        portion: see get_deltas()

    Returns:
        see get_deltas()
    """
    columns = [read_col, ref_col, start_col, end_col]
    alignments = hit_table[columns].itertuples(index=False, name=None)
    if zero_indexed:
        # convert [start, end) to [start + 1, end]
        alignments = (
            (read, ref, min(start, end) + 1, max(start, end))
            for read, ref, start, end in alignments
        )
    return get_deltas(alignments, portion=portion)


def get_deltas(alignments, portion=False):
    """
    Collect alignments into coverage deltas by reference sequence.

    Parameters:
        alignments: iterable of (read, reference, start, end) tuples, with
            1-based, inclusive positions (start and end can be in either
            order)
        portion: if True, a read with N alignments adds 1/N coverage for
            each of them (instead of 1)

    Returns:
        dict of {reference: {position: change_in_coverage}}. Coverage goes up
        at the start of each alignment and down just after its end.
    """
    if portion:
        # we need to know how many alignments each read has first
        alignments = list(alignments)
        hit_counts = defaultdict(int)
        for read, *_ in alignments:
            hit_counts[read] += 1

    deltas_by_ref = defaultdict(lambda: defaultdict(int))
    hit_count = 0
    for read, reference, start, end in alignments:
        hit_count += 1
        start, end = sorted((int(start), int(end)))
        factor = 1 / hit_counts[read] if portion else 1
        ref_deltas = deltas_by_ref[reference]
        ref_deltas[start] += factor
        ref_deltas[end + 1] -= factor

    LOGGER.debug(
        "parsed deltas for %d references from %d alignments",
        len(deltas_by_ref),
        hit_count,
    )
    return {reference: dict(deltas) for reference, deltas in deltas_by_ref.items()}


def deltas_to_cov(cov_deltas, x_max=None, nan_for_zero=True):
    """
    Convert one reference's coverage deltas into an array of coverage values.

    Parameters:
        cov_deltas: {position: change_in_coverage} (from get_deltas)
        x_max: length of the reference sequence. If None, the array ends at
            the last covered position.
        nan_for_zero: replace 0 coverage with NaN (so line plots have gaps
            where there is no coverage)

    Returns:
        numpy array of floats where element i is the coverage of position
        i + 1.
    """
    if x_max is None:
        x_max = max(cov_deltas) - 1 if cov_deltas else 0

    # numpy.add.at handles positions that appear more than once
    changes = numpy.zeros(x_max + 1)
    positions = numpy.fromiter(cov_deltas.keys(), dtype=int, count=len(cov_deltas))
    values = numpy.fromiter(cov_deltas.values(), dtype=float, count=len(cov_deltas))
    in_range = (positions >= 1) & (positions <= x_max)
    numpy.add.at(changes, positions[in_range] - 1, values[in_range])

    coverage = numpy.cumsum(changes)[:x_max]
    # fractional deltas (from portion=True) can leave rounding errors
    coverage[numpy.isclose(coverage, 0)] = numpy.nan if nan_for_zero else 0
    return coverage
