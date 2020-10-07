"""
methods for processing mapping results in SAM/BAM format

def parse_deltas(sam_file, ...):
    parse a sam/bam file into dicts of coverage changes by position
    
def deltas_to_cov(cov_deltas, x_max=None, nan_for_zero=True):
    convert coverage deltas into coverage array

class SAMFlag(IntFlag):
    class for decomposing SAM flags into bits for easier understanding
"""
import numpy
import subprocess
import enum
from collections import defaultdict
from edl import blastm8
import logging

SAMTOOLS_CMD_TEMPLATE = """samtools view -F 2308 {sam_file}"""


def parse_deltas(sam_file, samtools_cmd_template=SAMTOOLS_CMD_TEMPLATE, **kwargs):
    """
    Parse a SAM file into a collection of coverage deltas by reference sequences
    
    by default it uses samtools to strip all but the primary alignments first
    
    then hits are parsed with any kwargs sent to blastm8.FilterParams()
    """
    samtools_cmd = samtools_cmd_template.format(sam_file=sam_file)

    # dict of dicts of counts
    deltas_by_ref = defaultdict(lambda: defaultdict(int))
    read_count, hit_count = 0, 0
    with subprocess.Popen(samtools_cmd, shell=True, stdout=subprocess.PIPE,) as process:
        string_lines = (l.decode() for l in iter(process.stdout))
        for read, hits in blastm8.generate_hits(
            string_lines, format=blastm8.SAM, **kwargs
        ):
            read_count += 1
            for hit in hits:
                hit_count += 1
                start, end = sorted((hit.hstart, hit.hend))
                deltas = deltas_by_ref[hit.hit]
                deltas[start] += 1
                deltas[end + 1] -= 1

    logging.debug(
        "parsed deltas for %d contigs from %d reads and %d hits",
        len(deltas_by_ref),
        read_count,
        hit_count,
    )
    return deltas_by_ref


def deltas_to_cov(cov_deltas, x_max=None, nan_for_zero=True):
    """ converts dict of coverage deltas into array of coverage values
    
    cov_deltas: {pos: cov_change, ...}
    x_max: length of reference sequence (otherwise use last delta position)
    zero: replace 0 coverage with NAN so plots are discontinuous
    """
    zero = numpy.NAN if nan_for_zero else 0
    sorted_keys = sorted(cov_deltas)

    cov_segments = []
    last_pos = 1
    current_coverage = 0
    for pos in sorted_keys:
        delta = cov_deltas[pos]
        cov_segments.append(
            numpy.full(
                pos - last_pos, current_coverage if current_coverage != 0 else zero
            )
        )
        current_coverage += delta
        last_pos = pos
    cov_segments.append(
        numpy.full(
            1 if x_max is None else x_max + 1 - last_pos,
            current_coverage if current_coverage != 0 else zero,
        )
    )
    return numpy.concatenate(cov_segments)


class SAMFlag(enum.IntFlag):
    """
    From Wikipedia SAM Format entry:
    (int)	(binary)    	(meaning)
    1   	000000000001	template having multiple templates in sequencing (read is paired)
    2   	000000000010	each segment properly aligned according to the aligner (read mapped in proper pair)
    4   	000000000100	segment unmapped (read1 unmapped)
    8   	000000001000	next segment in the template unmapped (read2 unmapped)
    16  	000000010000	SEQ being reverse complemented (read1 reverse complemented)
    32  	000000100000	SEQ of the next segment in the template being reverse complemented (read2 reverse complemented)
    64  	000001000000	the first segment in the template (is read1)
    128 	000010000000	the last segment in the template (is read2)
    256 	000100000000	not primary alignment
    512 	001000000000	alignment fails quality checks
    1024	010000000000	PCR or optical duplicate
    2048	100000000000	supplementary alignment (e.g. aligner specific, could be a portion of a split read or a tied region)
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
