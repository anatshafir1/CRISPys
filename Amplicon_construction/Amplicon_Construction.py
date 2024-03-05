"""Constructing Amplicons from targets and SNPs"""
from typing import List, Tuple

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.FindTargets import get_targets
from Amplicon_construction.GetSequences import extract_exons_regions
from Amplicon_construction.Get_SNPs import get_snps
from Amplicon_construction.SNP_Obj import SNP_Obj
from Crispys import globals


def possible_snp(max_dist: int, curr_snp: SNP_Obj, curr_snps_list: List[SNP_Obj], distinct_alleles_num) -> bool:
    """
    Check if there is another SNP in SNPs list such that it will be "close enough" (to construct an amplicon of the
    maximum amplicon size) and will cover the number of alleles needed (to union of the different_alleles_sets of the
    SNPs should be at least the number of alleles minus 1).

    :param max_dist:
    :param curr_snp:
    :param curr_snps_list:
    :param distinct_alleles_num:
    :return:
    """

    curr_alleles_set = curr_snp.different_alleles_set
    for snp in curr_snps_list:
        if curr_snp.position_in_sequence + max_dist < snp.position_in_sequence:
            return False
        else:
            if len(curr_alleles_set.union(snp.different_alleles_set)) >= distinct_alleles_num - 1:
                return True
            else:
                curr_alleles_set = curr_alleles_set.union(snp.different_alleles_set)
    return False


def get_relevant_targets(targets_list: List):
    pass


def get_minimum_snps(curr_snp: SNP_Obj, snps_list: List[SNP_Obj], max_dist, distinct_alleles_num):


    snps_lists = []
    curr_alleles_set = curr_snp.different_alleles_set
    i = 0
    while curr_snp.position_in_sequence + max_dist < snps_list[i].position_in_sequence:
        curr_min_snps_list = [curr_snp]
        if len(curr_alleles_set.union(snps_list[i].different_alleles_set)) >= distinct_alleles_num - 1:
            curr_min_snps_list += [snps_list[i]]
            snps_lists += curr_min_snps_list
            i += 1
            curr_min_snps_list = [curr_snp]
        elif curr_alleles_set.union(snps_list[i].different_alleles_set) == curr_alleles_set:
            i += 1

    return False


def construct_amplicons(gene_snps_dict, max_amplicon_len: int, primer_length: int, target_len: int, targets_list: List,
                       snps_list: List[SNP_Obj], distinct_alleles_num: int) -> List[Amplicon_Obj]:
    """
    Given a list of sgRNA targets and a list of SNP of a gene - for every target find SNPs so that all the distinct
    alleles of the gene will be represented by an SNP with the shortest possible size of the Amplicon constructed from
    the target and it's SNPs.

    :param gene_snps_dict:
    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param target_len: length of the sgRNA target sequence in the gene, defined by user, depending on the CAS used.
    :param targets_list: list of sgRNA target of a gene
    :param snps_list: list of all the SNPs of the target's gene
    :param distinct_alleles_num: number of distinct alleles of the gene
    :return: list of Amplicons - the shortest amplicon found for every target
    """
    amplicons_list = []
    max_dist_snp = max_amplicon_len - primer_length * 2 - 1
    for exon in gene_snps_dict:
        # loop over SNPs in SNPs list
        for i in range(len(gene_snps_dict[exon])):
            if possible_snp(max_dist_snp, snps_list[i], snps_list[i+1:], distinct_alleles_num):
                min_snp_list = get_minimum_snps(snps_list[i], snps_list[i+1:])
                relevant_targets_list = get_relevant_targets(targets_list)
            else:
                continue

    return amplicons_list


def run_all(max_amplicon_len: int, primer_length: int, target_len: int, annotations_file_path,
            out_path: str, genome_fasta_file: str, distinct_alleles_num: int, pams: Tuple):
    """

    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param target_len: length of the sgRNA target sequence in the gene, defined by user, depending on the CAS used.
    :param annotations_file_path:
    :param out_path:
    :param genome_fasta_file:
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param pams:
    :return:
    """

    gene_sequences_dict = extract_exons_regions(max_amplicon_len, primer_length, target_len, annotations_file_path,
                                                out_path, genome_fasta_file, distinct_alleles_num)
    gene_snps_dict = get_snps(gene_sequences_dict, distinct_alleles_num)
    gene_targets_dict = get_targets(gene_sequences_dict, pams, max_amplicon_len, primer_length, target_len)
    return gene_snps_dict, gene_targets_dict
