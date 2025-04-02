import subprocess
from typing import List, Dict, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pandas as pd
from pandas import DataFrame
import warnings

from Amplicon_construction.GetSequences import extract_exons_regions

warnings.filterwarnings("ignore")


# noinspection PyTypeChecker


def annotations_to_lst_df(annotations_file_path: str) -> List[DataFrame]:
    """
    Given a gene annotations file path, construct a list of dataframes, with annotations of all the
    different alleles of a single exon. The start and end columns of each exon will include a region around the exon
    (part of the surrounding introns) from which an Amplicon can be later constructed.

    :param annotations_file_path: path to GFF file with annotations of the genome
    :return: a list of DataFrames, each with annotations of all the alleles of a single exon in the gene
    """
    annotations_df = pd.read_csv(annotations_file_path,
                                 sep='\s{2,}|\t')  # sep='\s{2,}' uses a regular expression to match two or more whitespace characters as the separator
    filtered_by_exon = annotations_df[annotations_df['feature'] == 'exon']

    filtered_by_exon['new_start'] = filtered_by_exon.apply(lambda x: x['start'] - 1,
                                                           axis=1)  # indices in genome fasta start from 1. getfasta calculates from 0. therefor "start" will be with -1
    filtered_by_exon.reset_index(drop=True, inplace=True)

    # add column of sequence IDs
    filtered_by_exon["seq_id"] = filtered_by_exon["genename"] + "_" + filtered_by_exon["seqname"]
    # Group exons by scaffold
    allele_groups = filtered_by_exon.groupby('seq_id')
    exons_dfs_list = []

    for _, group in allele_groups:
        # Sort exons by start index
        group.reset_index(drop=True, inplace=True)
        group_sorted = group.sort_values(by='new_start') if group['strand'].iloc[0] == "+" else group.sort_values(
            by='new_start', ascending=False)
        group_sorted.reset_index(drop=True, inplace=True)
        exons_dfs_list.append(group_sorted)
    return exons_dfs_list


def get_genomic_sites(out_path: str, fasta_file: str, filtered_allele_df: DataFrame) -> List[str]:
    """
    use Bedtools Getfasta to extract the sequences from the genome fasta file.

    :param out_path: the path to which the algorithm will store the results
    :param fasta_file: path to input FASTA format file of the genome
    :param filtered_allele_df: DataFrame of exon sequences and their parameters
    :return: list of strings where even indices are scaffold names and odd indices are sequences
    """
    # create BED format file
    filtered_allele_df.to_csv(out_path + '/exon_sites.bed', sep='\t',
                              columns=['seqname', 'new_start', 'end', 'attribute', 'score', 'strand'],
                              header=False, index=False)
    bed_file = out_path + "/exon_sites.bed"
    # run bedtools
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-s'],
                         stdout=subprocess.PIPE)
    sites_str = seq.stdout.decode()
    sites_list = sites_str.split()
    return sites_list


def genomic_sites_dict_to_fasta(gene_seqs, out_path: str):
    """

    :param gene_seqs: Dict of allele ID -> concatenated exons sequences or List of exon region sequences
    :param out_path:
    """
    sequences = []
    if isinstance(gene_seqs, dict):
        for seq_id, seq in gene_seqs.items():
            sequences += [SeqRecord(Seq(seq), id=seq_id, description='')]

    elif isinstance(gene_seqs, list):
        zipped_seqs = zip([gene_seqs[i] for i in range(0, len(gene_seqs), 2)],
                          [gene_seqs[i] for i in range(1, len(gene_seqs), 2)])
        for seq_id, seq in zipped_seqs:
            sequences += [SeqRecord(Seq(seq), id=seq_id, description='')]
    else:
        return
    SeqIO.write(sequences, out_path, "fasta")


def create_crispys_input_fasta(annotations_file_path: str, out_path: str, genome_fasta_file: str) -> str:
    """

    :param annotations_file_path: path to GFF file with annotations of the genome
    :param out_path: path to output directory where algorithm results will be saved
    :param genome_fasta_file: path to input FASTA format file of the genome
    :return: dictionary of exon number -> list of tuples of allele IDs and their sequences
    """
    print("extracting exon regions for crispys".upper().center(40, "#"))
    # Create list of DataFrames each representing an allele and its exons
    alleles_df_lst = annotations_to_lst_df(annotations_file_path)
    gene_seqs_dict = dict()  # keys are allele IDs, values are strings of concatenated CDSs sequences
    exon_indices_dict = dict()  # keys are allele IDs, values are zips of exon start and end indices
    allele_strand_dict = dict()
    concat_CDSs_path = out_path + "/concat_CDSs_crispys.fasta"

    for allele_df in alleles_df_lst:
        exons_seqs_lst = get_genomic_sites(out_path, genome_fasta_file, allele_df)  # extract exons from genome FASTA
        allele_exons_indices_zip = zip(allele_df['start'].to_list(),
                                       allele_df['end'].to_list())  # save start and end indices of exons by their order
        allele_exons_indices_lst = [(int(start), int(end)) for start, end in allele_exons_indices_zip]
        allele_ID = allele_df['seq_id'].iloc[0]
        allele_str = "".join([exons_seqs_lst[i] for i in range(1, len(exons_seqs_lst),
                                                               2)])  # create a single string of all the exon sequences concatenated
        gene_seqs_dict[allele_ID] = allele_str
        exon_indices_dict[allele_ID] = allele_exons_indices_lst
        allele_strand_dict[allele_ID] = allele_df['strand'].iloc[0]
    genomic_sites_dict_to_fasta(gene_seqs_dict, concat_CDSs_path)  # save alleles in FASTA file
    return concat_CDSs_path


def get_sequences_dict(max_amplicon_len: int, primer_length: int, target_surrounding_region: int, cut_location: int,
                       annotations_file_path: str, out_path: str, genome_fasta_file: str) -> Dict[str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]]:
    sequences_dict = {}
    annotations_df = pd.read_csv(annotations_file_path, sep='\s{2,}|\t')  # sep='\s{2,}' uses a regular expression to match two or more whitespace characters as the separator
    genes_list = annotations_df['genename'].unique().tolist()
    for gene in genes_list:
        gene_annotations_df = annotations_df[annotations_df["genename"] == gene]
        gene_annotations_df_path = out_path + f"/{gene}_annotations.txt"
        columns = gene_annotations_df.columns.tolist().remove("genename")
        gene_annotations_df.to_csv(gene_annotations_df_path, sep='\t', columns=columns, index=False)
        aligned_exons_regions_dict, aligned_to_original_exon_num_dict = extract_exons_regions(max_amplicon_len, primer_length, target_surrounding_region, cut_location,
                                               annotations_file_path, out_path, genome_fasta_file)
        sequences_dict[gene] = aligned_exons_regions_dict, aligned_to_original_exon_num_dict
    return sequences_dict
