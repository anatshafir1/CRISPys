
from typing import List

import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Amplicon_construction.GetSequences import genes_fasta_to_list


def get_genomic_sites(out_path: str, fasta_file: str) -> List[str]:
    """

    :param out_path: the path to which the algorithm will store the results
    :param fasta_file: path to input FASTA format file of the genome
    :return:
    """
    # create BED format file

    bed_file = out_path + "/genomic_sites.bed"
    # run bedtools
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-s'],
                         stdout=subprocess.PIPE)
    sites_list = seq.stdout.decode().split()
    return sites_list


# #### TEST IF INDICES ARE CORRECT #### #
genome_fasta_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Banana_GAL.Phased_Scaffolds.fasta"
output_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/test"
# remember to reduce start indices by 1 for getfasta
sites = get_genomic_sites(output_path, genome_fasta_path)


sequences = []
for i in range(0, len(sites), 2):
    sequences += [SeqRecord(Seq(sites[i+1]), id=sites[i])]
SeqIO.write(sequences, output_path + "/seqs_to_align", "fasta")


def test_regions(aligned_seq_path, out_path):
    exon_region_aligned_lst = genes_fasta_to_list(out_path+aligned_seq_path)
    return {1: exon_region_aligned_lst}, {1: 1}
