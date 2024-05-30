import math
import subprocess
import sys
from typing import List, Tuple, Dict

from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
from pandas import DataFrame
import warnings
warnings.filterwarnings("ignore")


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


def add_id_parent_columns(filtered_by_feature: DataFrame):
    """
    Given a DataFrame made from annotations GFF the function adds columns of object ID and object parent. This is
    done "in place", so the function returns None.

    :param filtered_by_feature: DataFrame made from annotations GFF
    """
    # add column with object ID
    start_char = 'ID='
    end_char = ';P'
    filtered_by_feature['start_position'] = filtered_by_feature['attribute'].str.find(start_char) + 3
    filtered_by_feature['end_position'] = filtered_by_feature['attribute'].str.find(end_char)
    filtered_by_feature['ID'] = filtered_by_feature.apply(
        lambda row: row['attribute'][row['start_position']:row['end_position']], axis=1)
    # add column with object parent
    start_char = 'Parent='
    end_char = ';T'
    filtered_by_feature['start_position'] = filtered_by_feature['attribute'].str.find(start_char) + 7
    filtered_by_feature['end_position'] = filtered_by_feature['attribute'].str.find(end_char)
    filtered_by_feature['parent'] = filtered_by_feature.apply(
        lambda row: row['attribute'][row['start_position']:row['end_position']], axis=1)
    filtered_by_feature.drop(['start_position', 'end_position'], axis=1, inplace=True)


def annotations_to_lst_df(max_amplicon_len: int, primer_length: int, cut_location: int,
                          annotations_file_path: str, distinct_alleles_num: int, target_surrounding_region: int) -> List[DataFrame]:
    """
    Given a gene annotations file path, construct a list of dataframes, with annotations of all the
    different alleles of a single exon. The start and end columns of each exon will include a region around the exon
    (part of the surrounding introns) from which an Amplicon can be later constructed.

    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run
    :param cut_location: length of the sgRNA target sequence in the gene, defined by user, depending on the CAS used.
    :param annotations_file_path: path to GFF file with annotations of the genome
    :param distinct_alleles_num: number of copies of each chromosome in the genome
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :return: a list of DataFrames, each with annotations of all the alleles of a single exon in the gene
    """
    annotations_df = pd.read_csv(annotations_file_path,
                                 sep='\s{2,}|\t')  # sep='\s{2,}' uses a regular expression to match two or more whitespace characters as the separator
    filtered_by_exon = annotations_df[annotations_df['feature'] == 'exon']
    add_id_parent_columns(filtered_by_exon)  # add columns of ID and Parent
    # calculate the maximum range upstream and downstream the exon, which are intron sites that are allowed be used to construct an amplicon
    exon_surrounding_seq_len = max_amplicon_len - cut_location - target_surrounding_region - primer_length
    # calculate the new start and end indexes of the sequence to extract from the genome fasta file
    filtered_by_exon['new_start'] = filtered_by_exon.apply(lambda x: x['start'] - exon_surrounding_seq_len - 1, axis=1)  # indices in genome fasta start from 1. getfasta calculates from 0. therefor "start" will be with -1
    filtered_by_exon['new_end'] = filtered_by_exon.apply(lambda x: x['end'] + exon_surrounding_seq_len, axis=1)
    filtered_by_exon.reset_index(drop=True, inplace=True)
    num_of_exons = int(len(filtered_by_exon) / distinct_alleles_num)  # calculate number of exons in the gene
    # Group exons by mRNA ID
    exon_groups = filtered_by_exon.groupby('parent')
    exons_dfs_list = [pd.DataFrame(columns=filtered_by_exon.columns) for _ in
                      range(num_of_exons)]  # a list of empty DataFrames with the columns for the annotations

    for _, group in exon_groups:
        # Sort exons by start index
        group.reset_index(drop=True, inplace=True)
        group_sorted = group.sort_values(by='new_start') if group['strand'].iloc[0] == "+" else group.sort_values(
            by='new_start', ascending=False)
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
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-s'],
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


def genes_fasta_to_list(aligned_genes: str) -> List[Tuple[str, str]]:
    """

    :param aligned_genes:
    :return:
    """
    gene_sequences = []

    # Iterate over the sequences in the FASTA file and save them to a list
    for record in SeqIO.parse(aligned_genes, "fasta"):
        gene_sequences.append((str(record.id), str(record.seq).upper()))
    return gene_sequences


def extract_exons_regions(max_amplicon_len: int, primer_length: int, target_surrounding_region: int, cut_location: int, annotations_file_path,
                          out_path: str, genome_fasta_file: str, distinct_alleles_num: int) -> Dict[int, List[Tuple[str, str]]]:
    """

    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param annotations_file_path: path to GFF file with annotations of the genome
    :param out_path: path to output directory where algorithm results will be saved
    :param genome_fasta_file: path to input FASTA format file of the genome
    :param distinct_alleles_num: number of copies of each chromosome in the genome
    :return: dictionary of exon number -> list of tuples of allele IDs and their sequences
    """

    aligned_exons_regions_dict = {}
    # create list of DataFrames of every exon in the gene with annotations of their alleles
    exon_regions_lst_df = annotations_to_lst_df(max_amplicon_len, primer_length, cut_location, annotations_file_path,
                                                distinct_alleles_num, target_surrounding_region)
    sliced_exon_regions_lst_df = exon_regions_lst_df[:math.floor(len(exon_regions_lst_df)*2/3)]  # ###TAKE ONLY FIRST 2/3 EXONS###
    for exon_num, exon_region in enumerate(sliced_exon_regions_lst_df):
        exon_regions_path = out_path + f"/exon_{exon_num + 1}_regions.fasta"
        aligned_exons_regions_path = out_path + f"/aligned_exon_{exon_num + 1}_regions.fasta"
        genomic_sites_list = get_genomic_sites(out_path, genome_fasta_file, exon_region)  # extract exon regions from genome FASTA
        genomic_sites_list_to_fasta(genomic_sites_list, exon_regions_path)  # save exon regions in FASTA file
        call_mafft(exon_regions_path, aligned_exons_regions_path)  # create an MSA of the alleles of the exon
        exon_region_aligned_lst = genes_fasta_to_list(aligned_exons_regions_path)  # save the aligned exon regions in a list
        aligned_exons_regions_dict[exon_num + 1] = exon_region_aligned_lst  # add the list of aligned exon regions to a dictionary
    return aligned_exons_regions_dict
