import sys


def call_primer3(in_file: str, out_file: str):
    """
    Calls the MAFFT algorithm to align the input DNA sequences of 'in_file' - a fasta format file. The aligned sequences
    are then stored in a fasta format file 'out_file'.

    :param in_file: fasta format file of sequences to be aligned
    :param out_file: fasta format file of aligned sequences
    """
    mafft_cline = MafftCommandline(cmd=sys.exec_prefix + "/bin/mafft", input=in_file)
    stdout, stderr = mafft_cline()
    with open(out_file, "w") as handle:
        handle.write(stdout)
