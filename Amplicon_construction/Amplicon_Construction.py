"""Constructing Amplicons from targets and SNPs"""
from typing import List

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.SNP_Obj import SNP_Obj
from Crispys.Candidate import Candidate
from Crispys import globals


def calc_amplicon_size(target: Candidate, first_snp: tuple[SNP_Obj, str], last_snp: tuple[SNP_Obj, str]) -> int:
    """
    Calculate the Initial Amplicon sequence length from the farthest SNP upstream the target to the farthest SNP downstream
    the target, or to the limits of the safe zone around the target, if there are only upstream or only downstream SNPs.

    :param target: sgRNA target of the Amplicon, as a Candidate object
    :param first_snp: first SNP upstream or downstream the target
    :param last_snp: last SNP upstream or downstream the target
    :return: Initial Amplicon sequence length
    """
    if first_snp[1] == "up":
        if last_snp[1] == "up":
            amp_size = target.position - first_snp[0].position_in_sequence + len(target.seq) + globals.safety_padding_around_target
        else:
            amp_size = last_snp[0].position_in_sequence - first_snp[0].position_in_sequence
    else:
        amp_size = last_snp[0].position_in_sequence - target.position + globals.safety_padding_around_target
    return amp_size


def create_relevant_snps_list(target: Candidate, snps_list: List[SNP_Obj]) -> List[(SNP_Obj, str)]:
    """
    Given a sgRNA target and a list of all the SNPs of the target's gene create a list of SNPs that are less than 921
    nucleotides upstream and downstream to the target. Also exclude SNPs which are in the safe zone around the target
    and the target sequence itself. The returned list has tuples with the SNP object and string "up" or "down" to
    indicate whether the SNP is upstream or downstream to the target.

    :param target: sgRNA target of the Amplicon, as a Candidate object
    :param snps_list: list of all the SNPs of the target's gene
    :return: list of SNPs with upstream/downstream notation as tuples
    """
    relevant_snp_list = list()
    # loop over SNPs list and create lists of SNPs for upstream and downstream seqs
    lower_limit = target.position - globals.safety_padding_around_target
    upper_limit = target.position + len(target.seq) + globals.safety_padding_around_target
    for snp in snps_list:
        if 1 <= (target.position - snp.position_in_sequence) <= 921 and not (lower_limit <= snp.position_in_sequence <= upper_limit):
            relevant_snp_list += [(snp, "up")]
        if -921 <= (target.position - snp.position_in_sequence) <= -1 and not (lower_limit <= snp.position_in_sequence <= upper_limit):
            relevant_snp_list += [(snp, "down")]
    return relevant_snp_list


def find_closest_snp(target: Candidate, snps_list: List[SNP_Obj]) -> bool:
    """
    Check if the current target has any SNP that is closer than 921 nucleotides upstream or downstream.

    :param target: sgRNA target of the Amplicon, as a Candidate object
    :param snps_list: list of all the SNPs of the target's gene
    :return: False if no relevant SNP found. True otherwise
    """
    max_distance = False
    for snp in snps_list:
        if abs(target.position - snp.position_in_sequence) <= 921:
            max_distance = True
            break
    return max_distance


def construct_amplicon(primer_length, candidates_list: List[Candidate], snps_list: List[SNP_Obj], distinct_alleles_num) -> List[Amplicon_Obj]:
    """
    Given a list of sgRNA targets and a list of SNP of a gene - for every target find SNPs so that all the distinct
    alleles of the gene will be represented by an SNP with the shortest possible size of the Amplicon constructed from
    the target and it's SNPs.

    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param candidates_list: list of sgRNA target of a gene, as Candidate objects
    :param snps_list: list of all the SNPs of the target's gene
    :param distinct_alleles_num: number of distinct alleles of the gene
    :return: list of Amplicons - the shortest amplicon found for every target
    """
    amplicons_list = []
    # loop over targets in candidates list
    for target in candidates_list:
        # check if the distance of the closest SNP to the target is 921 or less.
        if not find_closest_snp(target, snps_list):
            continue
        # create list of relevant SNPs for the target
        relevant_snp_list = create_relevant_snps_list(target, snps_list)
        # Initiate alleles set to be the set of the first SNP
        cur_alleles_set = set()
        # Initiate SNPs list with the first SNP
        cur_snps_list = list()
        min_amp_size = float("inf")
        shortest_amplicon = None
        for i in range(len(relevant_snp_list)):
            cur_alleles_set = cur_alleles_set.union(relevant_snp_list[i][0].different_alleles_set)
            cur_snps_list += [relevant_snp_list[i]]
            if len(cur_alleles_set) >= distinct_alleles_num-1:
                cur_amplicon_size = calc_amplicon_size(target, cur_snps_list[0], cur_snps_list[-1])
                if cur_amplicon_size < min_amp_size:
                    cur_amplicon = Amplicon_Obj(target.chromosome, target.gene, cur_amplicon_size, target.seq,
                                                primer_length, cur_snps_list)
                    shortest_amplicon = cur_amplicon
                    min_amp_size = cur_amplicon_size
                    cur_alleles_set = set()
                    cur_snps_list = list()
                    continue
            k = i + 1
            for j in range(k, len(relevant_snp_list)):
                if relevant_snp_list[j][0].different_alleles_set.issubset(cur_alleles_set):
                    continue
                cur_alleles_set = cur_alleles_set.union(relevant_snp_list[j][0].different_alleles_set)
                cur_snps_list += [relevant_snp_list[j]]
                if len(cur_alleles_set) >= distinct_alleles_num-1:
                    cur_amplicon_size = calc_amplicon_size(target, cur_snps_list[0], cur_snps_list[-1])
                    if cur_amplicon_size < min_amp_size:
                        cur_amplicon = Amplicon_Obj(target.chromosome, target.gene, cur_amplicon_size, target.seq,
                                                    primer_length, cur_snps_list)
                        shortest_amplicon = cur_amplicon
                        min_amp_size = cur_amplicon_size
                        cur_alleles_set = set()
                        cur_snps_list = list()
                        break

        amplicons_list += [shortest_amplicon]

    return amplicons_list
