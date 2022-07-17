"""SubgroupRes class containing file"""
from typing import List


class SubgroupRes:
    """
    An accessory class to store results of the Algorithm run with gene homology taken in consideration
    """
    def __init__(self, genes_lst: List, candidate_lst: List, name: str):
        self.genes_lst = genes_lst
        self.candidates_list = candidate_lst
        self.name = name
