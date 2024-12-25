from typing import List, Dict, Tuple
from collections import Counter

from SNP_Obj_2 import SNP_Obj


def create_idx_nuc_dict(scaffold_to_seq_dict, index, gap_length, distinct_alleles_num):
    idx_nuc_lst = [scaffold_to_seq_dict[allele][index + gap_length] for allele in scaffold_to_seq_dict.keys()]
    # create "ID" for current position as the dictionary curr_idx_nuc_dict
    idx_nuc_dict = {"gap_idxs": [i for i in range(distinct_alleles_num) if idx_nuc_lst[i] == "-"],
                    "nuc_idx": [i for i in range(distinct_alleles_num) if idx_nuc_lst[i] != "-"]}
    only_nucs_list = [idx_nuc_lst[i] for i in range(distinct_alleles_num) if idx_nuc_lst[i] != "-"]
    if len(only_nucs_list) > 1:
        idx_nuc_dict.update({"same_nucs": all(only_nucs_list[0] == nuc for nuc in only_nucs_list)})

    return idx_nuc_dict


def get_allele_sets_list(scaffold_to_seq_dict, index):
    nucs_dict = {}
    for allele in scaffold_to_seq_dict.keys():
        curr_nuc = scaffold_to_seq_dict[allele][index]
        if curr_nuc not in nucs_dict:
            nucs_dict[curr_nuc] = {allele}
        else:
            nucs_dict[curr_nuc].add(allele)
    return list(nucs_dict.values())


def create_snps_list(scaffold_to_seq_dict: Dict[str, str], distinct_alleles_num: int, primer_length: int) -> List[SNP_Obj]:
    """

    :param scaffold_to_seq_dict: dictionary of scaffold ID -> sequence of scaffold.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param primer_length: minimum length of the primer sequence.
    :return:
    """
    result_list = []
    seq_len = len(list(scaffold_to_seq_dict.values())[0])
    index = 0
    while index < seq_len:
        curr_idx_nuc_lst = [scaffold_to_seq_dict[allele][index] for allele in scaffold_to_seq_dict.keys()]
        counts = Counter(curr_idx_nuc_lst)
        if all(curr_idx_nuc_lst[0] == nuc for nuc in curr_idx_nuc_lst):  # no SNP in current index
            index += 1
        elif all(curr_idx_nuc_lst[i] != "-" for i in range(distinct_alleles_num)):  # no indel in current index
            alleles_sets_lst = get_allele_sets_list(scaffold_to_seq_dict, index)
            snp = SNP_Obj(index, alleles_sets_lst)
            if snp.position > primer_length:
                result_list.append(snp)
            index += 1
        elif all(counts[val] == 1 for val in
                 counts):  # all nucleotides are different in current index. ex. ('A', '-', 'G')
            alleles_sets_lst = get_allele_sets_list(scaffold_to_seq_dict, index)
            snp = SNP_Obj(index, alleles_sets_lst)
            if snp.position > primer_length:
                result_list.append(snp)
            index += 1
        else:  # indel/s in current position with same nucleotide/s in other alleles. ex. ('A', '-', 'A'), ('G', '-', '-')
            gap_length = 1
            alleles_sets_lst = get_allele_sets_list(scaffold_to_seq_dict, index)
            if index + gap_length < seq_len:
                # create "ID" for current position as the dictionary curr_idx_nuc_dict
                curr_idx_nuc_dict = create_idx_nuc_dict(scaffold_to_seq_dict, index, 0, distinct_alleles_num)
                # create "ID" for next position as the dictionary next_index_nuc_dict
                next_index_nuc_dict = create_idx_nuc_dict(scaffold_to_seq_dict, index, gap_length, distinct_alleles_num)
                # find gap length by comparing the "ID" dictionaries between positions
                while curr_idx_nuc_dict == next_index_nuc_dict and index + gap_length < seq_len:
                    gap_length += 1
                    if index + gap_length < seq_len:
                        # create "ID" for next position as the dictionary next_index_nuc_dict
                        next_index_nuc_dict = create_idx_nuc_dict(scaffold_to_seq_dict, index, gap_length,
                                                                  distinct_alleles_num)
                    else:
                        break
                snp = SNP_Obj(index, alleles_sets_lst, gap_length)
                if snp.position > primer_length:
                    result_list.append(snp)
                index += gap_length
            else:
                snp = SNP_Obj(index, alleles_sets_lst)
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
        scaffold_to_seq_dict = {seq_tup[0].split(":")[0][1:]: seq_tup[1] for seq_tup in
                                gene_sequences_dict[exon_region]}
        curr_exon_region_snps_list = create_snps_list(scaffold_to_seq_dict, distinct_alleles_num, primer_length)
        snps_dict[exon_region] = curr_exon_region_snps_list
    return snps_dict


def valid_amplicon_2(i: int, j: int, exon_snps_list: List[SNP_Obj], distinct_alleles_num: int) -> bool:
    """
    Check if SNPs i to j in exon_snps_list are enough to distinguish between the different alleles. The union of the sets
    of the alleles that each SNP covers should be at least the size of the number of different alleles minus 1 to create
    a set that will cover all the different alleles.

    :param i: index of current SNP.
    :param j: index of SNP downstream of "i".
    :param exon_snps_list: list of the SNPs of the current exon.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :return: True if sub-list of SNPs i to j are valid to distinguish between the alleles of an amplicon sequence.
    """
    singleton_lst = []
    k = i
    while len(singleton_lst) < distinct_alleles_num - 1 and k < j + 1:
        curr_snp = exon_snps_list[k]
        for set1 in curr_snp.alleles_sets_lst:
            if len(set1) == 1:
                if set1 not in singleton_lst:
                    singleton_lst.append(set1)
                    continue
            else:
                for m in range(k+1, j+1):
                    next_snp = exon_snps_list[m]
                    for set2 in next_snp.alleles_sets_lst:
                        intersect = set1.intersection(set2)
                        if len(intersect) == 1:  # intersection is singleton
                            if intersect not in singleton_lst:
                                singleton_lst.append(intersect)
                                if len(singleton_lst) >= distinct_alleles_num - 1:
                                    return True
        k += 1

    if len(singleton_lst) >= distinct_alleles_num - 1:
        return True
    else:
        return False
