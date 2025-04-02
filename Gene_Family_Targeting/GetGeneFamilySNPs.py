from typing import Dict, Tuple, List

from Amplicon_construction.GetSNPs import get_snps
from Amplicon_construction.SNP_Obj import SNP_Obj


def get_snps_dict(genes_exons_seq_dict: Dict[str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]],
                  distinct_alleles_num: int, primer_length: int) -> Dict[str, Dict[int, List[SNP_Obj]]]:
    print("finding snps for gene family".upper().center(40, "#"))
    genes_snps_dict = {}
    for gene_name in genes_exons_seq_dict:
        gene_seqs_dict = genes_exons_seq_dict[gene_name][0]
        gene_snps_dict = get_snps(gene_seqs_dict, distinct_alleles_num, primer_length)
        genes_snps_dict[gene_name] = gene_snps_dict
    return genes_snps_dict
