"""Amplicon class file"""
from typing import List

from Crispys import globals
from Amplicon_construction.SNP_Obj import SNP_Obj


class Amplicon_Obj:
    """A class representing a Constructed Amplicon for gene editing.
        a toy example:


        """

    def __init__(self, chromosome: int, gene: str, size: int, sgrna: str, primer_length: int, snps: List[SNP_Obj]):

        self.chromosome = chromosome
        self.gene = gene
        self.size = size
        self.sgrna = sgrna
        self.safety_padding_around_target = globals.safety_padding_around_target
        self.primer_length = primer_length
        self.snps = snps

    def __str__(self):
        return f"{self.chromosome}, {self.gene}, {self.size}, {self.sgrna}, {self.snps}"

    def __repr__(self):
        return self.__str__()
