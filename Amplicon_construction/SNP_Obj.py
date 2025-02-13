"""SNP class file"""

from typing import List


class SNP_Obj:
    """A class representing a Single Nucleotide Polymorphism site in a gene.
        a toy example:

        10234, {3}

        In this SNP object the polymorphisms occur at position 10234 of the aligned alleles. Allele 1 and 2  are
        similar in their nucleotides in that position, while allele 3 has a different nucleotide in the position.
        """

    def __init__(self, position: int, alleles_sets_lst: List[set], gap_length: int = 1):
        self.position = position
        self.alleles_sets_lst = alleles_sets_lst
        self.gap_length = gap_length

    def __str__(self):
        return f"{self.position}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (self.position == other.position and self.alleles_sets_lst == other.alleles_sets_lst and
                self.gap_length == other.gap_length)

    def __hash__(self):
        return hash(self.__str__())

    def update_snp_index(self, length_to_subtract):
        self.position = self.position - length_to_subtract
