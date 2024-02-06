try:
    from  Bio import SeqIO
except ModuleNotFoundException:
    print("Error loading Bio.SeqIO. Install biopython to use this module.")
    raise

def get_seq_lens(sequence_file, format='fasta'):
    return {r.id:len(r) for r in SeqIO.parse(sequence_file, format=format)}

