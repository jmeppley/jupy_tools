import importlib
import io
import sys

import pytest


@pytest.fixture
def biopython():
    pytest.importorskip("Bio")
    from jme.jupy_tools import biopython
    return biopython


FASTA = """>seq1 a description
ACGTACGTAC
GTACG
>seq2
AC
>seq3
"""

FASTQ = """@read1
ACGTA
+
IIIII
@read2
AC
+
II
"""


def test_get_seq_lens_fasta(biopython, tmp_path):
    path = tmp_path / "seqs.fasta"
    path.write_text(FASTA)
    assert biopython.get_seq_lens(str(path)) == \
        {'seq1': 15, 'seq2': 2, 'seq3': 0}


def test_get_seq_lens_handle(biopython):
    assert biopython.get_seq_lens(io.StringIO(FASTA))['seq1'] == 15


def test_get_seq_lens_fastq(biopython):
    assert biopython.get_seq_lens(io.StringIO(FASTQ), format='fastq') == \
        {'read1': 5, 'read2': 2}


@pytest.mark.xfail(strict=True, raises=NameError, reason=(
    "the import guard catches ModuleNotFoundException, which doesn't exist, "
    "so a missing biopython gives a NameError instead of an ImportError"))
def test_missing_biopython(monkeypatch):
    monkeypatch.setitem(sys.modules, 'Bio', None)
    monkeypatch.delitem(sys.modules, 'jme.jupy_tools.biopython',
                        raising=False)
    try:
        with pytest.raises(ImportError):
            importlib.import_module('jme.jupy_tools.biopython')
    finally:
        sys.modules.pop('jme.jupy_tools.biopython', None)
