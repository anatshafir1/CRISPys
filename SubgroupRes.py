"""SubgroupRes class containing file"""
from typing import List


class SubgroupRes:
    """
    An accessory class to store results of the Algorithm run with gene homology taken in consideration
    """
    def __init__(self, genes_lst: List, candidate_lst: List, name: str, genes_in_node: List = [], family_name: str = "gene_family"):
        self.genes_lst = genes_lst
        self.candidates_list = candidate_lst
        self.name = name
        self.genes_in_node = genes_in_node
        self.family_name = family_name
