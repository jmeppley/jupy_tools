try:
    from Bio import SeqIO
except ImportError as e:
    raise ModuleNotFoundError(
        "Error loading Bio.SeqIO. Install biopython to use this module.",
        name=e.name,
    ) from e

def get_seq_lens(sequence_file, format='fasta'):
    return {r.id:len(r) for r in SeqIO.parse(sequence_file, format=format)}

