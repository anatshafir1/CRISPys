"""SNP class file"""


class SNP_Obj:
    """A class representing a Single Nucleotide Polymorphism site in a gene.
        a toy example:
        10234, {3}
        In this SNP object the polymorphisms occur at position 10234 of the aligned alleles. Allele 1 and 2  are
        similar in their nucleotides in that position, while allele 3 has a different nucleotide in the position.
        """

    def __init__(self, position_in_sequence: int, different_alleles_set: set):
        self.position_in_sequence = position_in_sequence
        self.different_alleles_set = different_alleles_set

    def __str__(self):
        return f"{self.position_in_sequence}, {self.different_alleles_set}"

    def __repr__(self):
        return self.__str__()
