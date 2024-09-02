
from typing import List, Dict, Tuple
from collections import Counter

from SNP_Obj import SNP_Obj


def create_snps_dict(exon_sequences_lst: List[Tuple[str, str]], distinct_alleles_num: int, primer_length: int) -> List[SNP_Obj]:
    """

    :param exon_sequences_lst: list of current exon's alleles sequences as tuples of (allele_ID, sequence)
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run
    :return:
    """
    snps_list = []
    number_to_scaffold_dict = {i: seq_tup[0].split(":")[0][1:] for i, seq_tup in enumerate(exon_sequences_lst)}
    all_alleles_set = set(seq_tup[0].split(":")[0][1:] for seq_tup in exon_sequences_lst)
    only_sequences_lst = [exon_sequences_lst[i][1] for i in range(distinct_alleles_num)]
    index = 0
    while index < len(only_sequences_lst[0]):
        curr_index_nuc_list = [only_sequences_lst[allele][index] for allele in range(distinct_alleles_num)]
        if all(curr_index_nuc_list[0] == allele for allele in curr_index_nuc_list):  # no SNP in current index
            index += 1
        elif all(allele != "-" for allele in curr_index_nuc_list):  # no indel in current index
            counts = Counter(curr_index_nuc_list)
            unique_alleles = {number_to_scaffold_dict[i] for i, val in enumerate(curr_index_nuc_list) if counts[val] == 1}
            snp = SNP_Obj(index, unique_alleles)
            if snp.position > primer_length:
                snps_list.append(snp)
            index += 1
        else:  # at least one indel in an allele at current index
            allele_set = set()
            index_run = index + 1
            gap_length = 1
            for j in range(distinct_alleles_num):  # find which allele starts with an indel
                if only_sequences_lst[j][index] == "-":
                    while index_run < len(only_sequences_lst[0]) and only_sequences_lst[j][index_run] == "-":  # calculate gap length
                        index_run += 1
                        gap_length += 1
                    if gap_length > 1:
                        break
            if gap_length > 1:
                for j in range(distinct_alleles_num):
                    if only_sequences_lst[j][index: index_run].count("-") == gap_length:
                        allele_set.add(number_to_scaffold_dict[j])
                if len(allele_set) > 1:
                    distinct_allele = all_alleles_set.difference(allele_set)
                    snp = SNP_Obj(index, distinct_allele, gap_length)
                    if snp.position > primer_length:
                        snps_list.append(snp)
                    index += gap_length
                else:
                    snp = SNP_Obj(index, allele_set, gap_length)
                    if snp.position > primer_length:
                        snps_list.append(snp)
                    index += gap_length
            else:
                counts = Counter(curr_index_nuc_list)
                unique_alleles = {number_to_scaffold_dict[i] for i, val in enumerate(curr_index_nuc_list) if
                                  counts[val] == 1}
                snp = SNP_Obj(index, unique_alleles)
                if snp.position > primer_length:
                    snps_list.append(snp)
                index += 1
    return snps_list


def get_snps(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], distinct_alleles_num: int,
             primer_length: int) -> Dict[int, List[SNP_Obj]]:
    """

    :param gene_sequences_dict: dictionary of exon number to list of exon's alleles sequences as tuples of (allele_ID, sequence)
    :param distinct_alleles_num: number of distinct alleles of the gene
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run
    :return:
    """
    snps_dict = {}
    for exon_region in gene_sequences_dict:
        curr_exon_region_snps_list = create_snps_dict(gene_sequences_dict[exon_region], distinct_alleles_num, primer_length)
        snps_dict[exon_region] = curr_exon_region_snps_list
    return snps_dict
