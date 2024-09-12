
from typing import List, Dict, Tuple
from collections import Counter

from SNP_Obj import SNP_Obj


def create_snps_dict(scaffold_to_seq_dict: Dict[str, str], distinct_alleles_num: int, primer_length: int) -> List[SNP_Obj]:
    """

    :param scaffold_to_seq_dict:
    :param distinct_alleles_num:
    :param primer_length:
    :return:
    """
    result_list = []
    all_alleles_ids_set = set(scaffold_to_seq_dict.keys())
    all_alleles_ids_list = list(scaffold_to_seq_dict.keys())
    seq_len = len(list(scaffold_to_seq_dict.values())[0])
    index = 0
    while index < seq_len:
        curr_index_nuc_list = [scaffold_to_seq_dict[allele][index] for allele in scaffold_to_seq_dict.keys()]
        counts = Counter(curr_index_nuc_list)
        if all(curr_index_nuc_list[0] == nuc for nuc in curr_index_nuc_list):  # no SNP in current index
            index += 1
        elif all(curr_index_nuc_list[i] != "-" for i in range(distinct_alleles_num)):  # no indel in current index
            unique_alleles = {all_alleles_ids_list[i] for i, val in enumerate(curr_index_nuc_list) if counts[val] == 1}
            snp = SNP_Obj(index, unique_alleles)
            if snp.position > primer_length:
                result_list.append(snp)
            index += 1
        elif all(counts[val] == 1 for val in counts):  # all nucleotides are different in current index. ex. ('A', '-', 'G')
            snp = SNP_Obj(index, all_alleles_ids_set)
            if snp.position > primer_length:
                result_list.append(snp)
            index += 1
        else:  # indel/s in current position with same nucleotide/s in other alleles
            gap_length = 1
            if distinct_alleles_num == 2:
                curr_distinct_allele = all_alleles_ids_set
            else:
                curr_distinct_allele = {all_alleles_ids_list[i] for i, val in enumerate(curr_index_nuc_list) if counts[val] == 1}
            if index+gap_length < seq_len:
                next_index_nuc_list = [scaffold_to_seq_dict[allele][index+gap_length] for allele in scaffold_to_seq_dict.keys()]
                while curr_index_nuc_list == next_index_nuc_list and index+gap_length < seq_len:
                    gap_length += 1
                    next_index_nuc_list = [scaffold_to_seq_dict[allele][index + gap_length] for allele in
                                           scaffold_to_seq_dict.keys()]
                snp = SNP_Obj(index, curr_distinct_allele, gap_length)
                if snp.position > primer_length:
                    result_list.append(snp)
                index += gap_length
            else:
                snp = SNP_Obj(index, curr_distinct_allele)
                if snp.position > primer_length:
                    result_list.append(snp)
                index += 1
    return result_list


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
        scaffold_to_seq_dict = {seq_tup[0].split(":")[0][1:]: seq_tup[1] for seq_tup in gene_sequences_dict[exon_region]}
        curr_exon_region_snps_list = create_snps_dict(scaffold_to_seq_dict, distinct_alleles_num, primer_length)
        snps_dict[exon_region] = curr_exon_region_snps_list
    return snps_dict
