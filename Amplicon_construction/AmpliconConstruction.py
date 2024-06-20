"""Constructing Amplicons from targets and SNPs"""
import statistics
from typing import List, Tuple, Dict

import pandas as pd

from Amplicon_Obj import Amplicon_Obj
from FindTargets import get_targets
from alleles_similarity_evaluation import extract_exons_regions
from Get_SNPs import get_snps
from primer3 import get_primers
from SNP_Obj import SNP_Obj
from Target_Obj import Target_Obj

import warnings
warnings.filterwarnings("ignore")


def get_relevant_targets(gene_targets_dict: Dict[int, List[Target_Obj]], gene_snps_dict: Dict[int, List[SNP_Obj]]) -> Dict[int, List[Target_Obj]]:
    """
    loop over every target in every exon of the gene, and check if any SNPs of that exon are on the target sequence.

    :param gene_targets_dict:
    :param gene_snps_dict:
    :return:
    """
    new_targets_dict = {}
    for exon in gene_snps_dict:
        new_exon_targets_lst = []
        for target in gene_targets_dict[exon]:
            if all(snp.position_in_sequence not in range(target.start_idx, target.end_idx) for snp in gene_snps_dict[exon]):
                new_exon_targets_lst.append(target)
        new_targets_dict[exon] = new_exon_targets_lst

    return new_targets_dict


def valid_amplicon(i: int, j: int, exon_snps_list: List[SNP_Obj], distinct_alleles_num: int) -> bool:
    """
    Check if SNPs i to j in exon_snps_list are enough to distinguish between the different alleles. The union of the sets
    of the alleles that each SNP covers should be at least the size of the number of different alleles minus 1 to create
    a set that will cover all the different alleles.

    :param i: index of current SNP
    :param j: index of SNP downstream of i
    :param exon_snps_list: list of the SNPs of the current exon
    :param distinct_alleles_num: number of distinct alleles of the gene
    :return:
    """
    allele_set = exon_snps_list[i].different_alleles_set
    for k in range(i + 1, j + 1):
        allele_set = allele_set.union(exon_snps_list[k].different_alleles_set)
        if len(allele_set) >= distinct_alleles_num - 1:
            return True
    return False


def valid_sgrna_target(i: int, j: int, exon_snps_list: List[SNP_Obj], distinct_alleles_num: int,
                       curr_target: Target_Obj, max_amplicon_len: int, primer_length: int, target_surrounding_region: int) -> List[SNP_Obj]:
    """
    Check if current target can be used to construct an amplicon.

    :param i: index of current SNP
    :param j: index of SNP downstream of i
    :param exon_snps_list: list of the SNPs of the current exon
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param curr_target:
    :param max_amplicon_len:
    :param primer_length:
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :return:
    """
    valid_snps_list = []
    # Make sure that target is not too far upstream or downstream from first and last SNPs
    max_range_from_snp = max_amplicon_len - primer_length * 2 - 1
    first_snp = exon_snps_list[i]
    last_snp = exon_snps_list[j]
    if not (curr_target.end_idx < first_snp.position_in_sequence + max_range_from_snp and curr_target.start_idx >
            last_snp.position_in_sequence - max_range_from_snp):
        return []
    # Create a list of SNPs that can be used to calculate the Amplicon statistics (SNPs that don't fall on the region around the target)
    for k in range(i, j + 1):
        if not curr_target.start_idx - target_surrounding_region <= exon_snps_list[k].position_in_sequence <= \
               curr_target.end_idx + target_surrounding_region:
            valid_snps_list.append(exon_snps_list[k])
    if len(valid_snps_list) == 0:
        return []
    # Make sure that the new SNPs list in enough for a valid amplicon
    if valid_amplicon(0, len(valid_snps_list) - 1, valid_snps_list, distinct_alleles_num):
        return valid_snps_list
    return []


def calculate_snps_statistics(valid_snps_for_target, distinct_alleles_num) -> Tuple[float, float]:
    """

    :param valid_snps_for_target:
    :param distinct_alleles_num:
    :return:
    """
    counts_lst = [0 for _ in range(distinct_alleles_num)]
    for snp in valid_snps_for_target:
        for allele_num in snp.different_alleles_set:
            counts_lst[allele_num] += 1
    snps_median = statistics.median(counts_lst)
    snps_mean = statistics.mean(counts_lst)

    return snps_median, snps_mean


# noinspection PyTypeChecker
def create_candidate_amplicon(snps_median: float, snps_mean: float, candidate_amplicon_snps_lst: List[SNP_Obj],
                              target: Target_Obj, max_amplicon_len: int, primer_length: int, exon_num: int,
                              current_exon_region_id: str, target_surrounding_region: int, min_amplicon_len: int) -> Amplicon_Obj:
    """

    :param snps_median:
    :param snps_mean:
    :param candidate_amplicon_snps_lst:
    :param target:
    :param max_amplicon_len:
    :param primer_length:
    :param exon_num:
    :param current_exon_region_id:
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param min_amplicon_len: minimum length of the amplicon, defined by user
    :return:
    """
    first_snp = candidate_amplicon_snps_lst[0]
    last_snp = candidate_amplicon_snps_lst[-1]
    max_range_from_snp = max_amplicon_len - primer_length - 1
    max_range_from_target = max_amplicon_len - target_surrounding_region - primer_length
    min_amp_start_index = target.end_idx - max_range_from_target if target.end_idx > last_snp.position_in_sequence else last_snp.position_in_sequence - max_range_from_snp
    max_amp_end_index = target.start_idx + max_range_from_target if target.end_idx < first_snp.position_in_sequence else first_snp.position_in_sequence + max_range_from_snp
    if max_amp_end_index - min_amp_start_index < min_amplicon_len:
        return Amplicon_Obj(exon_num, current_exon_region_id, "", "", 0, 0, 0.0, 0.0, target, candidate_amplicon_snps_lst, None)
    return Amplicon_Obj(exon_num, current_exon_region_id, "", "", min_amp_start_index, max_amp_end_index, snps_median, snps_mean, target, candidate_amplicon_snps_lst, None)


def get_candidate_amplicons(i: int, j: int, exon_snps_lst: List[SNP_Obj], distinct_alleles_num: int,
                            current_targets_lst: List[Target_Obj], max_amplicon_len: int,
                            primer_length: int, exon_num: int, current_exon_region_id: str,
                            target_surrounding_region: int, min_amplicon_len: int) -> List[Amplicon_Obj]:
    """

    :param i: index of current SNP
    :param j: index of SNP downstream of i
    :param exon_snps_lst: list of the SNPs of the current exon
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param current_targets_lst:
    :param max_amplicon_len:
    :param primer_length:
    :param exon_num:
    :param current_exon_region_id:
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param min_amplicon_len: minimum length of the amplicon, defined by user
    :return:
    """

    candidate_amplicons_list = []
    for target in current_targets_lst:
        valid_snps_for_target = valid_sgrna_target(i, j, exon_snps_lst, distinct_alleles_num, target,
                                                   max_amplicon_len, primer_length, target_surrounding_region)  # SNPs for Amplicon statistics
        if len(valid_snps_for_target) > 0:
            candidate_amplicon_snps_lst = exon_snps_lst[i:j + 1]  # all the SNPs of the Amplicon Candidate
            snps_median, snps_mean = calculate_snps_statistics(valid_snps_for_target, distinct_alleles_num)
            candidate_amplicon = create_candidate_amplicon(snps_median, snps_mean, candidate_amplicon_snps_lst, target,
                                                           max_amplicon_len, primer_length, exon_num, current_exon_region_id, target_surrounding_region, min_amplicon_len)
            if candidate_amplicon.snps_median != 0:
                candidate_amplicons_list.append(candidate_amplicon)

    return candidate_amplicons_list


def construct_amplicons(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                        gene_snps_dict: Dict[int, List[SNP_Obj]], gene_targets_dict: Dict[int, List[Target_Obj]],
                        max_amplicon_len: int, primer_length: int, distinct_alleles_num: int,
                        target_surrounding_region: int, min_amplicon_len: int) -> List[Amplicon_Obj]:
    """
    Given a dictionary of sgRNA targets, a dictionary of SNP of a gene and a dictionary of exon region sequences -
    for every SNP of every exon find possible amplicons that can be constructed using that SNP.

    :param gene_exon_regions_seqs_dict: dictionary of exon numbers -> a list of allele sequences and their description
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region
    :param gene_targets_dict: dictionary of exon numbers -> a list sgRNA targets of the exon
    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param min_amplicon_len: minimum length of the amplicon, defined by user
    :return: list of Amplicons
    """
    candidate_amplicons_list = []
    max_dist_snp = max_amplicon_len - primer_length * 2 - 1
    for exon in gene_snps_dict:
        exon_targets_lst = gene_targets_dict[exon]
        current_exon_region = gene_exon_regions_seqs_dict[exon][0][1]
        current_exon_region_id = gene_exon_regions_seqs_dict[exon][0][0]
        exon_snps_lst = gene_snps_dict[exon]
        max_range_from_snp = min(len(current_exon_region), max_dist_snp)
        for i in range(len(exon_snps_lst)):
            j = i + 1
            current_snp = exon_snps_lst[i]
            while j < len(exon_snps_lst) and exon_snps_lst[j].position_in_sequence < current_snp.position_in_sequence + max_range_from_snp:
                if valid_amplicon(i, j, exon_snps_lst, distinct_alleles_num):  # check if snps i to j are enough do distinct between different alleles
                    current_snp_candidate_amplicons = get_candidate_amplicons(i, j, exon_snps_lst,
                                                                              distinct_alleles_num, exon_targets_lst,
                                                                              max_amplicon_len, primer_length, exon, current_exon_region_id, target_surrounding_region, min_amplicon_len)
                    candidate_amplicons_list.extend(current_snp_candidate_amplicons)
                j += 1

    return candidate_amplicons_list


def save_results_to_csv(res_amplicons_lst, out_path):

    df = pd.DataFrame([amplicon.to_dict() for amplicon in res_amplicons_lst])
    df.to_csv(out_path + "/results.csv", index=False)


def get_amplicons(max_amplicon_len_category: int, primer_length: int, target_surrounding_region: int, cut_location: int, annotations_file_path: str,
                  out_path: str, genome_fasta_file: str, distinct_alleles_num: int, pams: Tuple, target_len: int,
                  primer3_core_path: str, n: int, genome_chroms_path: str, filter_off_targets: int):
    """


    :param max_amplicon_len_category: category of maximum length of the amplicon, defined by user
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param annotations_file_path: path to GFF file with annotations of the genome
    :param out_path: path to output directory where algorithm results will be saved
    :param genome_fasta_file: path to input FASTA format file of the genome
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param pams:
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :param primer3_core_path:
    :param n: desired maximum number of amplicons to return
    :param genome_chroms_path: path to the directory in which the genome separated by scaffold fasta files are.
    :param filter_off_targets: choose whether to filter amplicons with 'strong' off-targets for their gRNAs, or return
    them in the results.
    :return:
    """
    amplicon_ranges = [(200, 300), (300, 500), (500, 1000)]
    max_amplicon_len = max(amplicon_ranges[max_amplicon_len_category - 1])
    min_amplicon_len = min(amplicon_ranges[max_amplicon_len_category - 1])
    gene_exon_regions_seqs_dict = extract_exons_regions(max_amplicon_len, primer_length, target_surrounding_region, cut_location,
                                                        annotations_file_path, out_path, genome_fasta_file)
    gene_snps_dict = get_snps(gene_exon_regions_seqs_dict, distinct_alleles_num)
    gene_targets_dict = get_targets(gene_exon_regions_seqs_dict, pams, max_amplicon_len, primer_length, cut_location,
                                    target_surrounding_region, target_len)
    relevant_gene_targets_dict = get_relevant_targets(gene_targets_dict, gene_snps_dict)
    candidate_amplicons_list = construct_amplicons(gene_exon_regions_seqs_dict, gene_snps_dict,
                                                   relevant_gene_targets_dict, max_amplicon_len, primer_length,
                                                   distinct_alleles_num, target_surrounding_region, min_amplicon_len)
    sorted_candidate_amplicons = sorted(candidate_amplicons_list, key=lambda amplicon: (amplicon.snps_median, amplicon.snps_mean), reverse=True)
    amplicon_obj_with_primers = get_primers(gene_exon_regions_seqs_dict, sorted_candidate_amplicons, out_path,
                                            primer3_core_path, n, amplicon_ranges[max_amplicon_len_category - 1],
                                            distinct_alleles_num, target_surrounding_region, filter_off_targets, genome_chroms_path)
    if len(amplicon_obj_with_primers) > 0:
        save_results_to_csv(amplicon_obj_with_primers, out_path)
        return amplicon_obj_with_primers
    else:
        if max_amplicon_len_category < 3:
            get_amplicons(max_amplicon_len_category + 1, primer_length, target_surrounding_region, cut_location, annotations_file_path,
                          out_path, genome_fasta_file, distinct_alleles_num, pams, target_len,
                          primer3_core_path, n, genome_chroms_path, filter_off_targets)
        else:
            print("No amplicons found")
