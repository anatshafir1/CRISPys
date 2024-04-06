
from typing import List, Dict, Tuple
from collections import Counter

from Amplicon_construction.SNP_Obj import SNP_Obj


def create_snps_dict(fasta_sequences_lst: List, distinct_alleles_num: int) -> Dict[int, Tuple[str]]:
    """

    :param fasta_sequences_lst:
    :param distinct_alleles_num:
    :return:
    """
    snps_dict = {}
    only_sequences_lst = [fasta_sequences_lst[i][1] for i in range(distinct_alleles_num)]
    only_sequences_zip = zip(*only_sequences_lst)
    for index, locus in enumerate(only_sequences_zip):
        if not all(locus[0] == nucleotide for nucleotide in locus):
            snps_dict[index] = locus  # the indices of nucleotides in the tuple are the alleles numbers

    return snps_dict


def snps_dict_to_lst(snps_dict: Dict[int, Tuple[str]]) -> List[SNP_Obj]:
    """

    :param snps_dict:
    :return:
    """
    snps_lst = []
    for index in snps_dict:
        counter = Counter(snps_dict[index])
        most_common_nuc = counter.most_common(1)[0][0]
        if most_common_nuc == "-":
            continue
        else:
            current_tuple = snps_dict[index]
            list_snps_nuc = [current_tuple.index(nuc) for nuc in current_tuple if nuc != most_common_nuc]
            snp = SNP_Obj(int(index), set(list_snps_nuc))
            if snp.position_in_sequence > 18:
                snps_lst.append(snp)
    return snps_lst


def get_snps(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], distinct_alleles_num: int) -> Dict[int, List[SNP_Obj]]:
    """

    :param gene_sequences_dict:
    :param distinct_alleles_num:
    :return:
    """
    snps_dict = {}
    for exon_region in gene_sequences_dict:
        curr_exon_region_snps_dict = create_snps_dict(gene_sequences_dict[exon_region], distinct_alleles_num)
        curr_exon_region_snps_list = snps_dict_to_lst(curr_exon_region_snps_dict)
        snps_dict[exon_region] = curr_exon_region_snps_list
    return snps_dict
