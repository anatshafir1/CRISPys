"""Amplicon class file"""
from typing import List

from Primers_Obj import Primers_Obj
from Target_Obj import Target_Obj
from SNP_Obj import SNP_Obj


class Amplicon_Obj:
    """A class representing a Constructed Amplicon for gene editing.
        a toy example:


        """

    def __init__(self, exon_num: int, scaffold: str, strand: str, sequence: str, start_idx: int, end_idx: int, snps_median: float, snps_mean: float, target: Target_Obj, snps: List[SNP_Obj], primers: Primers_Obj):

        self.scaffold = scaffold
        self.strand = strand
        self.exon_num = exon_num
        self.sequence = sequence
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.size = end_idx - start_idx + 1
        self.snps_median = snps_median
        self.snps_mean = snps_mean
        self.target = target
        self.snps = snps
        self.primers = primers

    def __str__(self):
        return f"scaffold: {self.scaffold}, size: {self.size}, target: <{self.target}>, SNPs: <{self.snps}>, primers: <{self.primers}>"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.size

    def to_dict(self):
        """Create a dictionary of the Amplicon object"""
        self_dict = self.__dict__.copy()
        self_dict.pop("target")
        self_dict.pop("snps")
        self_dict.pop("primers")
        snps_str = ""
        for snp in self.snps:
            snps_str += f"{snp}"
        self_dict["snps"] = snps_str
        self_dict.update(self.target.to_dict())
        self_dict.update(self.primers.to_dict())
        return self_dict

    def update_snps_indices(self, length_to_subtract):
        for snp in self.snps:
            snp.update_snp_index(length_to_subtract)
