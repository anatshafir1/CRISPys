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
    mafft_cline = MafftCommandline(cmd=sys.exec_prefix + "/bin/mafft", input=in_file, genafpair=True)
    stdout, stderr = mafft_cline()
    with open(out_file, "w") as handle:
        handle.write(stdout)


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

    filtered_by_exon['new_start'] = filtered_by_exon.apply(lambda x: x['start'] - 1, axis=1)  # indices in genome fasta start from 1. getfasta calculates from 0. therefor "start" will be with -1
    filtered_by_exon.reset_index(drop=True, inplace=True)
    # Group exons by scaffold
    allele_groups = filtered_by_exon.groupby('seqname')
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
    sites_list = seq.stdout.decode().split()
    return sites_list


def genomic_sites_dict_to_fasta(genomic_sites_dict, out_path: str):
    """

    :param genomic_sites_dict:
    :param out_path:
    """
    sequences = []
    if isinstance(genomic_sites_dict, dict):
        for seq_id, seq in genomic_sites_dict.items():
            sequences += [SeqRecord(Seq(seq), id=seq_id)]

    elif isinstance(genomic_sites_dict, list):
        zipped_seqs = zip([genomic_sites_dict[i] for i in range(0, len(genomic_sites_dict), 2)],
                          [genomic_sites_dict[i] for i in range(1, len(genomic_sites_dict), 2)])
        for seq_id, seq in zipped_seqs:
            sequences += [SeqRecord(Seq(seq), id=seq_id)]
    else:
        return
    SeqIO.write(sequences, out_path, "fasta")


def genes_fasta_to_dict(aligned_genes: str) -> Tuple[Dict[str, str], int]:
    """

    :param aligned_genes:
    :return:
    """
    gene_sequences = dict()
    alignment_len = 0
    # Iterate over the sequences in the FASTA file and save them to a list
    for record in SeqIO.parse(aligned_genes, "fasta"):
        gene_sequences[str(record.id)] = str(record.seq).upper()
        alignment_len = len(str(record.seq))
    return gene_sequences, alignment_len


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


def get_overlapping_exons(exon_indices_dict: Dict[str, List[Tuple[int, int]]], alleles_exons_aligned_dict: Dict[str, str],
                          alignment_len: int, min_nucs_for_amplicon: int) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
    """
    Using aligned sequences of concatenated exons and exon start,end indices (given by annotations file) calculate
    which exons are properly aligned to enable potential amplicon construction.

    :param exon_indices_dict: dictionary of allele scaffold -> list of tuples of exons start,end indices
    :param alleles_exons_aligned_dict: dictionary of allele scaffold -> sequences of alignment of the alleles of concatenated exons
    :param alignment_len: length of the alignment of the alleles of concatenated exons
    :param min_nucs_for_amplicon: minimum necessary number of nucleotides on exon sequence to construct a potential amplicon
    :return:
    """
    aligned_overlapping_exons_dict = {allele: [] for allele in exon_indices_dict}
    prev_allele_idx = 0
    curr_allele_idx = 0
    alleles_to_loop_over_again = dict()
    nuc_count_dict = {allele: 0 for allele in exon_indices_dict}
    exon_num_dict = {allele: 0 for allele in exon_indices_dict}
    num_of_exons = 0

    while curr_allele_idx < alignment_len:

        for allele in alleles_exons_aligned_dict:
            curr_allele_seq = alleles_exons_aligned_dict[allele]
            curr_exon_num = exon_num_dict[allele]
            if curr_allele_seq[curr_allele_idx] != "-":  # make sure the aligned sequence doesn't start with a gap - 'real' start of exon sequence
                curr_exon_len = abs(exon_indices_dict[allele][curr_exon_num][1] - exon_indices_dict[allele][curr_exon_num][0]) + 1
                while nuc_count_dict[allele] < curr_exon_len:  # calculate the end index of the current exon in the alignment sequence
                    if curr_allele_seq[curr_allele_idx] != "-":
                        nuc_count_dict[allele] += 1
                    curr_allele_idx += 1
                for other_allele in alleles_exons_aligned_dict:  # add the other alleles to dictionary to be looped over
                    if other_allele != allele and other_allele not in alleles_to_loop_over_again:
                        alleles_to_loop_over_again[other_allele] = alleles_exons_aligned_dict[other_allele]
                if curr_allele_idx <= alignment_len:
                    exon_num_dict[allele] += 1
                    nuc_count_dict[allele] = 0
                break
            else:  # sequence starts with "-" (gap) and will be looped over again
                alleles_to_loop_over_again[allele] = curr_allele_seq

        legit_exon_dict = {allele: False for allele in alleles_to_loop_over_again}
        for allele in alleles_to_loop_over_again:  # loop over the other alleles and check if they are properly aligned to the first allele
            curr_allele_seq = alleles_exons_aligned_dict[allele]
            curr_exon_num = exon_num_dict[allele]
            curr_exon_len = abs(exon_indices_dict[allele][curr_exon_num][1] - exon_indices_dict[allele][curr_exon_num][0]) + 1
            indel_cnt = curr_allele_seq[prev_allele_idx:curr_allele_idx].count("-")
            nuc_cnt = curr_allele_idx - prev_allele_idx - indel_cnt
            if nuc_cnt > min_nucs_for_amplicon:  # enough nucleotides to use current exon allele to construct potential amplicons
                legit_exon_dict[allele] = True
            nuc_count_dict[allele] += nuc_cnt
            if nuc_count_dict[allele] >= curr_exon_len and curr_allele_idx <= alignment_len:  # update exons and nucleotides count for current allele
                exon_num_dict[allele] += 1
                nuc_count_dict[allele] = nuc_count_dict[allele] - curr_exon_len
        if all(legit_exon_dict[allele] is True for allele in legit_exon_dict):  # all exon alleles in the alignment are properly aligned
            for allele in aligned_overlapping_exons_dict:  # update result dictionary with the exon indices which are aligned properly
                aligned_overlapping_exons_dict[allele].append(exon_indices_dict[allele][exon_num_dict[allele] - 1])
            num_of_exons += 1
        prev_allele_idx = curr_allele_idx

    return aligned_overlapping_exons_dict, num_of_exons


def get_legit_exons_regions(annotations_file_path: str, out_path: str, genome_fasta_file: str, min_nucs_for_amplicon: int) -> \
        Tuple[Dict[str, List[Tuple[int, int]]], int, Dict[str, str]]:
    """

    :param annotations_file_path: path to GFF file with annotations of the genome
    :param out_path: path to output directory where algorithm results will be saved
    :param genome_fasta_file: path to input FASTA format file of the genome
    :param min_nucs_for_amplicon: minimum necessary number of nucleotides on exon sequence to construct a potential amplicon
    :return: dictionary of exon number -> list of tuples of allele IDs and their sequences
    """

    # Create list of DataFrames each representing an allele and its exons
    alleles_df_lst = annotations_to_lst_df(annotations_file_path)
    gene_seqs_dict = dict()  # keys are allele IDs, values are strings of concatenated CDSs sequences
    exon_indices_dict = dict()  # keys are allele IDs, values are zips of exon start and end indices
    allele_strand_dict = dict()
    concat_CDSs_path = out_path + "/concat_CDSs.fasta"
    aligned_concat_CDSs_path = out_path + "/aligned_concat_CDSs.fasta"

    for allele_df in alleles_df_lst:
        exons_seqs_lst = get_genomic_sites(out_path, genome_fasta_file, allele_df)  # extract exons from genome FASTA
        allele_exons_indices_zip = zip(allele_df['start'].to_list(),
                                       allele_df['end'].to_list())  # save start and end indices of exons by their order
        allele_exons_indices_lst = [(int(start), int(end)) for start, end in allele_exons_indices_zip]
        allele_ID = allele_df['seqname'].iloc[0]
        allele_str = "".join([exons_seqs_lst[i] for i in range(1, len(exons_seqs_lst), 2)])  # create a single string of all the exon sequences concatenated
        gene_seqs_dict[allele_ID] = allele_str
        exon_indices_dict[allele_ID] = allele_exons_indices_lst
        allele_strand_dict[allele_ID] = allele_df['strand'].iloc[0]
    genomic_sites_dict_to_fasta(gene_seqs_dict, concat_CDSs_path)  # save alleles in FASTA file
    call_mafft(concat_CDSs_path, aligned_concat_CDSs_path)  # create an MSA of the alleles of the gene
    alleles_exons_aligned_dict, alignment_len = genes_fasta_to_dict(aligned_concat_CDSs_path)  # save the aligned CDSs to a dict
    aligned_overlapping_exons_dict, num_of_exons = get_overlapping_exons(exon_indices_dict, alleles_exons_aligned_dict, alignment_len, min_nucs_for_amplicon)
    return aligned_overlapping_exons_dict, num_of_exons, allele_strand_dict


def get_exon_dict(aligned_overlapping_exons_dict: Dict[str, List[Tuple[int, int]]], exon_num: int,
                  exon_surrounding_seq_len, allele_strand_dict) -> Dict[str, List[str]]:
    """
    create a dictionary of overlapping exons with their attributes

    :param aligned_overlapping_exons_dict: dictionary of allele scaffold -> list of tuples exon start,end indices
    :param exon_num: current exon number
    :param exon_surrounding_seq_len: number of nucleotides from each side of the exon
    :param allele_strand_dict: dictionary of allele scaffold -> strand in genome fasta file
    :return: dictionary of overlapping exons with their attributes'
    """
    exon_dict = {'seqname': [], 'new_start': [], 'end': [], 'attribute': [], 'score': [], 'strand': []}
    for allele in aligned_overlapping_exons_dict:
        exon_dict['seqname'].append(allele)
        exon_dict['new_start'].append(
            aligned_overlapping_exons_dict[allele][exon_num][0] - exon_surrounding_seq_len - 1)
        exon_dict['end'].append(aligned_overlapping_exons_dict[allele][exon_num][1] + exon_surrounding_seq_len)
        exon_dict['attribute'].append(exon_num)
        exon_dict['score'].append(".")
        exon_dict['strand'].append(allele_strand_dict[allele])
    return exon_dict


def extract_exons_regions(max_amplicon_len: int, primer_length: int, target_surrounding_region: int, cut_location: int,
                          annotations_file_path,
                          out_path: str, genome_fasta_file: str) -> Dict[int, List[Tuple[str, str]]]:
    """

    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param annotations_file_path: path to GFF file with annotations of the genome
    :param out_path: path to output directory where algorithm results will be saved
    :param genome_fasta_file: path to input FASTA format file of the genome
    :return: dictionary of exon number -> list of tuples of allele IDs and their sequences
    """
    min_nucs_for_amplicon = primer_length + target_surrounding_region + cut_location
    aligned_overlapping_exons_dict, num_of_exons, allele_strand_dict = get_legit_exons_regions(annotations_file_path,
                                                                                               out_path,
                                                                                               genome_fasta_file, min_nucs_for_amplicon)
    exon_surrounding_seq_len = max_amplicon_len - cut_location - target_surrounding_region - primer_length
    aligned_exons_regions_dict = {}

    for exon_num in range(num_of_exons):
        exon_regions_path = out_path + f"/exon_{exon_num + 1}_regions.fasta"
        aligned_exons_regions_path = out_path + f"/aligned_exon_{exon_num + 1}_regions.fasta"
        exon_dict = get_exon_dict(aligned_overlapping_exons_dict, exon_num, exon_surrounding_seq_len, allele_strand_dict)
        exon_region = pd.DataFrame.from_dict(exon_dict)
        genomic_sites_list = get_genomic_sites(out_path, genome_fasta_file, exon_region)  # extract exon regions from genome FASTA
        genomic_sites_dict_to_fasta(genomic_sites_list, exon_regions_path)  # save exon regions in FASTA file
        call_mafft(exon_regions_path, aligned_exons_regions_path)  # create an MSA of the alleles of the exon
        exon_region_aligned_lst = genes_fasta_to_list(aligned_exons_regions_path)  # save the aligned exon regions in a list
        aligned_exons_regions_dict[exon_num + 1] = exon_region_aligned_lst  # add the list of aligned exon regions to a dictionary
    return aligned_exons_regions_dict
