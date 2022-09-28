"""MAFFT and ProtDist functions"""
__author__ = 'GH'

# Unix version
import subprocess
import sys
from typing import List
import os
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO  # app to read fasta files (added by Udi)
import globals


########################################################################################################################
# the pipeline: make FASTA file, run mafft on it to get the output in a fasta file, convert this FASTA to PHYLIP,
# and run PHYLIP's tool for making distance matrix##
########################################################################################################################

def make_fasta_file(names_lst: List, seq_lst: List, outfile: str):
    """
    this function creates a fasta format file from given gene names and their sequences. gene names are written to
    individual lines starting with ">", and the gene sequence is written to the following line. e.g.:
    ">GG0"
    "TGTCGATGAACCCGGTGGCAATCCCCATAGACGAAGGACCTAGTGGCCACGA..."
    ">GG2"
    "TGGCCGACGACGACGAGATCGCTCT..."

    :param names_lst: list of gene names
    :param seq_lst: list of gene sequences
    :param outfile: the fasta format file to which the results are stored
    """
    f = open(outfile, 'w')
    for i in range(len(names_lst)):
        f.write('>' + names_lst[i] + '\n')
        f.write(seq_lst[i] + '\n')
    f.close()


def call_mafft(in_file: str, out_file: str):
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


def call_protdist(outpath: str):
    """
    calls the protDist algorithm to create a distance matrix from a phylip format file. the input file for protdist
    algorithm has to be a phylip format, in the 'outpath' and under the name 'infile' (these will be asserted in
    'create_protdist' and 'fasta_to_phylip' functions). the output distance matrix will be stored in 'outpath' under the
    name 'outfile'.

    :param outpath: the path of the input and output file for protdist
    """
    os.chdir(outpath)
    # when running the algorithm on Unix operating systems
    if globals.protdist_path is None:
        os.system(f'echo "Y\r\n" | "{sys.exec_prefix}/bin/protdist"')
    # when running the algorithm from Microsoft Windows using POSIX:
    else:
        subprocess.run([globals.protdist_path], input=b"Y")


def fasta_to_phylip(in_f: str, out_f: str):
    """
	this function changes the format of the alignments from fasta to phylip using biopython. written by Udi 25/01/22

	:param in_f: a fasta format file of aligned sequences
	:param out_f: a phylip format file of aligned sequences
	"""
    aligned_genes = list(SeqIO.parse(in_f, "fasta"))
    SeqIO.write(aligned_genes, out_f, "phylip")


def create_protdist(names_lst: List, seq_lst: List, out_path: str):
    """
    this function creates a protDist file for creating a distance matrix, using the MAFFT algorithm
    and PHYLIP's protDist algorithm. Given lists of gene names and their sequences this function aligns the sequences
    and returns a file with the distances between them according to the homology of the proteins the genes are encoding.
    :param names_lst: a list of names of genes to align
    :param seq_lst: a list of genes to align
    :param out_path: the output_path where the files created using this function will be stored
    """
    genes_fasta_for_mafft = out_path + "/genes_fasta_for_mafft.fa"
    mafft_output_aligned_fasta = out_path + "/mafft_output_aligned_fasta.fa"
    phylip_file = out_path + "/infile"
    make_fasta_file(names_lst, seq_lst, genes_fasta_for_mafft)
    call_mafft(genes_fasta_for_mafft, mafft_output_aligned_fasta)
    fasta_to_phylip(mafft_output_aligned_fasta, phylip_file)
    call_protdist(out_path)
