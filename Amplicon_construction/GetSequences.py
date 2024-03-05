import subprocess
import sys
from typing import List

from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Crispys import globals
import pandas as pd
from pandas import DataFrame


# noinspection PyTypeChecker


def call_mafft(in_file: str, out_file: str):
    """
    Calls the MAFFT algorithm to align the input DNA sequences of 'in_file' - a fasta format file. The aligned sequences
    are then stored in a fasta format file 'out_file'.

    :param in_file: path to fasta format file of sequences to be aligned
    :param out_file: path to output of mafft - aligned sequences fasta file
    """
    mafft_cline = MafftCommandline(cmd=sys.exec_prefix + "/bin/mafft", input=in_file)
    stdout, stderr = mafft_cline()
    with open(out_file, "w") as handle:
        handle.write(stdout)


def add_id_parent_columns(filtered_by_feature):
    # add column with object ID
    start_char = 'ID='
    end_char = ';P'
    filtered_by_feature['start_position'] = filtered_by_feature['attribute'].str.find(start_char) + 3
    filtered_by_feature['end_position'] = filtered_by_feature['attribute'].str.find(end_char)
    filtered_by_feature['ID'] = filtered_by_feature.apply(
        lambda row: row['attribute'][row['start_position']:row['end_position']], axis=1)
    start_char = 'Parent='
    end_char = ';T'
    filtered_by_feature['start_position'] = filtered_by_feature['attribute'].str.find(start_char) + 7
    filtered_by_feature['end_position'] = filtered_by_feature['attribute'].str.find(end_char)
    filtered_by_feature['parent'] = filtered_by_feature.apply(
        lambda row: row['attribute'][row['start_position']:row['end_position']], axis=1)
    filtered_by_feature.drop(['start_position', 'end_position'], axis=1, inplace=True)


def annotations_to_lst_df(max_amplicon_len: int, primer_length: int, target_len: int,
                          annotations_file_path: str, distinct_alleles_num: int) -> List[DataFrame]:
    """
    Given a gene annotations file path, construct a list of dataframes, each representing the same exon from all the
    different alleles. The start and end columns of each exon will include a region around the exon from which an
    Amplicon can be later constructed.

    :param max_amplicon_len:
    :param primer_length:
    :param target_len:
    :param annotations_file_path:
    :param distinct_alleles_num:
    :return:
    """
    annotations_df = pd.read_csv(annotations_file_path, sep='\s{2,}')  # sep='\s{2,}' uses a regular expression to match two or more whitespace characters as the separator
    filtered_by_exon = annotations_df[annotations_df['feature'] == 'exon']
    add_id_parent_columns(filtered_by_exon)  # add columns of ID and Parent
    exon_surrounding_seq_len = max_amplicon_len - primer_length - target_len - globals.safety_padding_around_target * 2
    filtered_by_exon['new_start'] = filtered_by_exon.apply(lambda x: x['start'] - exon_surrounding_seq_len - 1, axis=1)
    filtered_by_exon['new_end'] = filtered_by_exon.apply(lambda x: x['end'] + exon_surrounding_seq_len, axis=1)
    filtered_by_exon.reset_index(drop=True, inplace=True)
    # Group exons by mRNA ID
    num_of_exons = int(len(filtered_by_exon)/distinct_alleles_num)
    exon_groups = filtered_by_exon[filtered_by_exon['feature'] == 'exon'].groupby('parent')
    exons_dfs_list = [pd.DataFrame(columns=filtered_by_exon.columns) for _ in range(num_of_exons)]
    # Iterate over the groups of exons
    for _, group in exon_groups:
        # Sort exons by start index
        group.reset_index(drop=True, inplace=True)
        group_sorted = group.sort_values(by='new_start') if group['strand'].iloc[0] == "+" else group.sort_values(by='new_start', ascending=False)
        group_sorted.reset_index(drop=True, inplace=True)
        for index, row in group_sorted.iterrows():
            exons_dfs_list[index].loc[len(exons_dfs_list[index].index)] = row
    return exons_dfs_list


def get_genomic_sites(out_path: str, fasta_file: str, filtered_gene_df) -> List[str]:
    """

    :param out_path: the path to which the algorithm will store the results
    :param fasta_file: path to input FASTA format file of the genome
    :param filtered_gene_df:
    :return:
    """
    # create BED format file
    filtered_gene_df.to_csv(out_path + '/genomic_sites.bed', sep='\t',
                            columns=['seqname', 'new_start', 'new_end', 'ID', 'score', 'strand'],
                            header=False, index=False)
    bed_file = out_path + "/genomic_sites.bed"
    # run bedtools
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-nameOnly', '-s'],
                         stdout=subprocess.PIPE)
    sites_list = seq.stdout.decode().split()
    return sites_list


def genomic_sites_list_to_fasta(genomic_sites_list: List, out_path: str):
    """

    :param genomic_sites_list:
    :param out_path:
    """
    sequences = []
    zipped_seqs = zip([genomic_sites_list[i] for i in range(0, len(genomic_sites_list), 2)],
                      [genomic_sites_list[i] for i in range(1, len(genomic_sites_list), 2)])
    for seq_id, seq in zipped_seqs:
        sequences += [SeqRecord(Seq(seq), id=seq_id)]
    SeqIO.write(sequences, out_path, "fasta")


def genes_fasta_to_list(aligned_genes: str) -> List[str]:
    """

    :param aligned_genes:
    :return:
    """
    gene_sequences = []

    # Iterate over the sequences in the FASTA file
    for record in SeqIO.parse(aligned_genes, "fasta"):
        gene_sequences.append(str(record.id))
        gene_sequences.append(str(record.seq).upper())
    return gene_sequences


def extract_exons_regions(max_amplicon_len: int, primer_length: int, target_len: int, annotations_file_path,
                          out_path: str, genome_fasta_file: str, distinct_alleles_num: int):
    """

    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param target_len: length of the sgRNA target sequence in the gene, defined by user, depending on the CAS used.
    :param annotations_file_path:
    :param out_path:
    :param genome_fasta_file:
    :param distinct_alleles_num: number of distinct alleles of the gene
    :return:
    """

    aligned_exons_regions_dict = {}
    exon_regions_path = out_path + "/exons_regions.fasta"
    aligned_exons_regions_path = out_path + "/aligned_exons_regions.fasta"
    exon_regions_lst_df = annotations_to_lst_df(max_amplicon_len, primer_length, target_len, annotations_file_path, distinct_alleles_num)
    for exon_num, exon_region in enumerate(exon_regions_lst_df):
        genomic_sites_list = get_genomic_sites(out_path, genome_fasta_file, exon_region)
        genomic_sites_list_to_fasta(genomic_sites_list, exon_regions_path)
        call_mafft(exon_regions_path, aligned_exons_regions_path)
        exon_region_aligned_lst = genes_fasta_to_list(aligned_exons_regions_path)
        aligned_exons_regions_dict[exon_num+1] = exon_region_aligned_lst
    return aligned_exons_regions_dict
