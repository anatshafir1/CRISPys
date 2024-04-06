"""Constructing Amplicons from targets and SNPs"""
import statistics
from typing import List, Tuple, Dict

from Amplicon_construction.FindTargets import get_targets
from Amplicon_construction.GetSequences import extract_exons_regions
from Amplicon_construction.Get_SNPs import get_snps
from Amplicon_construction.SNP_Obj import SNP_Obj
from Crispys import globals


def get_relevant_targets(gene_targets_dict: Dict, gene_snps_dict: Dict) -> Dict:
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
            if all(snp.position_in_sequence not in range(target[1], target[2]) for snp in gene_snps_dict[exon]):
                new_exon_targets_lst.append(target)
        new_targets_dict[exon] = new_exon_targets_lst

    return new_targets_dict


def valid_amplicon(i: int, j: int, curr_snps_list: List[SNP_Obj], distinct_alleles_num: int) -> bool:
    """
    Check if SNPs i to j in curr_snps_list are enough to distinguish between the different alleles. The union of the sets
    of the alleles that each SNP covers should be at least the size of the number of different alleles minus 1 to create
    a set that will cover all the different alleles.

    :param i: index of current SNP
    :param j: index of SNP downstream of i
    :param curr_snps_list: list of the SNPs of the current exon
    :param distinct_alleles_num: number of distinct alleles of the gene
    :return:
    """
    allele_set = curr_snps_list[i].different_alleles_set
    for k in range(i + 1, j + 1):
        allele_set = allele_set.union(curr_snps_list[k].different_alleles_set)
        if len(allele_set) >= distinct_alleles_num - 1:
            return True
    return False


def valid_sgrna_target(i: int, j: int, curr_snps_list: List[SNP_Obj], distinct_alleles_num: int,
                       curr_target: Tuple[str, int, int, str], max_amplicon_len: int, primer_length: int) -> List:
    """
    Check if current target can be used to construct an amplicon.

    :param i: index of current SNP
    :param j: index of SNP downstream of i
    :param curr_snps_list: list of the SNPs of the current exon
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param curr_target:
    :param max_amplicon_len:
    :param primer_length:
    :return:
    """
    valid_snps_list = []
    max_range_from_snp = max_amplicon_len - primer_length*2 - 1
    if not curr_target[2] < curr_snps_list[0].position_in_sequence + max_range_from_snp:
        return []
    for k in range(i, j + 1):
        if not curr_target[1] - globals.safety_padding_around_target <= curr_snps_list[k].position_in_sequence <= \
               curr_target[2] + globals.safety_padding_around_target:
            valid_snps_list.append(curr_snps_list[k])
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


def create_candidate_amplicon(snps_median: float, snps_mean: float, candidate_amplicon_snps_lst: List[SNP_Obj],
                              target: Tuple[str, int, int, str], max_amplicon_len: int, primer_length: int) -> Tuple:
    """

    :param snps_median:
    :param snps_mean:
    :param candidate_amplicon_snps_lst:
    :param target:
    :param max_amplicon_len:
    :param primer_length:
    :return:
    """
    max_range_from_snp = max_amplicon_len - primer_length - 1
    max_range_from_target = max_amplicon_len - globals.safety_padding_around_target - primer_length
    min_amp_start_index = max(candidate_amplicon_snps_lst[-1].position_in_sequence - max_range_from_snp, target[2] - max_range_from_target)
    max_amp_end_index = min(candidate_amplicon_snps_lst[0].position_in_sequence + max_range_from_snp, target[1] + max_range_from_target)
    return min_amp_start_index, max_amp_end_index, snps_median, snps_mean, target, candidate_amplicon_snps_lst


def get_candidate_amplicons(i: int, j: int, current_snps_lst: List[SNP_Obj], distinct_alleles_num: int,
                            current_targets_lst: List[Tuple[str, int, int, str]], max_amplicon_len: int,
                        primer_length: int) -> List:
    """

    :param i: index of current SNP
    :param j: index of SNP downstream of i
    :param current_snps_lst: list of the SNPs of the current exon
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param current_targets_lst:
    :param max_amplicon_len:
    :param primer_length:
    :return:
    """
    candidate_amplicons_list = []
    for k in range(len(current_targets_lst)):
        valid_snps_for_target = valid_sgrna_target(i, j, current_snps_lst, distinct_alleles_num, current_targets_lst[k], max_amplicon_len, primer_length)
        if valid_snps_for_target:
            candidate_amplicon_snps_lst = current_snps_lst[i:j+1]
            snps_median, snps_mean = calculate_snps_statistics(valid_snps_for_target, distinct_alleles_num)
            candidate_amplicon = create_candidate_amplicon(snps_median, snps_mean, candidate_amplicon_snps_lst, current_targets_lst[k], max_amplicon_len, primer_length)
            candidate_amplicons_list.append(candidate_amplicon)
    return candidate_amplicons_list


def construct_amplicons(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                        gene_snps_dict: Dict[int, List[SNP_Obj]],
                        gene_targets_dict: Dict[int, List[Tuple[str, int, int, str]]], max_amplicon_len: int,
                        primer_length: int,
                        distinct_alleles_num: int) -> List:
    """
    Given a dictionary of sgRNA targets, a dictionary of SNP of a gene and a dictionary of exon region sequences -
    for every SNP of every exon find possible amplicons that can be constructed using that SNP.

    :param gene_exon_regions_seqs_dict: dictionary of exon numbers -> a list of allele sequences and their description
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region
    :param gene_targets_dict: dictionary of exon numbers -> a list sgRNA targets of the exon
    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: length of the primer sequence, defined by the user in the algorithm run
    :param distinct_alleles_num: number of distinct alleles of the gene
    :return: list of Amplicons
    """
    candidate_amplicons_list = []
    max_dist_snp = max_amplicon_len - primer_length * 2 - 1
    for exon in gene_snps_dict:
        current_targets_lst = gene_targets_dict[exon]
        current_exon_region = gene_exon_regions_seqs_dict[exon][1][1]
        current_snps_lst = gene_snps_dict[exon]
        max_range_from_snp = min(len(current_exon_region), max_dist_snp)
        for i in range(len(current_snps_lst)):
            j = i + 1
            current_snp = current_snps_lst[i]
            while j < len(current_snps_lst) and current_snps_lst[j].position_in_sequence < current_snp.position_in_sequence + max_range_from_snp:
                if valid_amplicon(i, j, current_snps_lst, distinct_alleles_num):  # check if snps i to j are enough do distinct between different alleles
                    current_snp_candidate_amplicons = get_candidate_amplicons(i, j, current_snps_lst, distinct_alleles_num, current_targets_lst, max_amplicon_len, primer_length)
                    candidate_amplicons_list.append(current_snp_candidate_amplicons)
                j += 1

    return candidate_amplicons_list


def run_all(max_amplicon_len: int, primer_length: int, cut_location: int, annotations_file_path: str,
            out_path: str, genome_fasta_file: str, distinct_alleles_num: int, pams: Tuple, target_len: int):
    """

    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param annotations_file_path: path to GFF file with annotations of the genome
    :param out_path: path to output directory where algorithm results will be saved
    :param genome_fasta_file: path to input FASTA format file of the genome
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param pams:
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :return:
    """

    gene_exon_regions_seqs_dict = extract_exons_regions(max_amplicon_len, primer_length, cut_location,
                                                        annotations_file_path,
                                                        out_path, genome_fasta_file, distinct_alleles_num)
    gene_snps_dict = get_snps(gene_exon_regions_seqs_dict, distinct_alleles_num)
    gene_targets_dict = get_targets(gene_exon_regions_seqs_dict, pams, max_amplicon_len, primer_length, cut_location,
                                    target_len)
    relevant_gene_targets_dict = get_relevant_targets(gene_targets_dict, gene_snps_dict)
    candidate_amplicons_list = construct_amplicons(gene_exon_regions_seqs_dict, gene_snps_dict, relevant_gene_targets_dict, max_amplicon_len, primer_length, distinct_alleles_num)
    return gene_snps_dict, gene_targets_dict, relevant_gene_targets_dict, candidate_amplicons_list
