"""Amplicon class file"""
from typing import List

from Amplicon_construction.Primers_Obj import Primers_Obj
from Amplicon_construction.Target_Obj import Target_Obj
from Amplicon_construction.SNP_Obj import SNP_Obj


def calc_amplicon_size(seq):
    nuc_count = 0
    for char in seq:
        if char != "-":
            nuc_count += 1
    return nuc_count


class Amplicon_Obj:
    """A class representing a Constructed Amplicon for gene editing.
        a toy example:


        """

    def __init__(self, exon_num: int, scaffold: str, strand: str, sequence: str, start_idx: int, end_idx: int, snps_median: float, snps_mean: float, target, snps: List[SNP_Obj], primers: Primers_Obj, orig_exon_num: int = 0):

        self.scaffold = scaffold
        self.strand = strand
        self.exon_num = exon_num
        self.orig_exon_num = orig_exon_num
        self.sequence = sequence
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.size = calc_amplicon_size(sequence)
        self.snps_median = snps_median
        self.snps_mean = snps_mean
        self.target = target
        self.snps = snps
        self.primers = primers
        self.off_targets = []

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
        self_dict.pop("off_targets")
        snps_str = ""
        for snp in self.snps:
            snps_str += f"{snp}"
        self_dict["snps"] = snps_str
        self_dict.update(self.target.to_dict())
        self_dict.update(self.primers.to_dict())
        if len(self.off_targets) > 0:
            self_dict.update(self.off_targets[0].to_dict(1))
        if len(self.off_targets) > 1:
            self_dict.update(self.off_targets[1].to_dict(2))
        return self_dict

    def update_snps_indices(self, length_to_subtract):
        for snp in self.snps:
            snp.update_snp_index(length_to_subtract)

    def sort_off_targets(self):
        """
        Sort the off_targets_list by the off-targets scores from highest to lowest
        """
        self.off_targets.sort(key=lambda off: -off.score)


class OffTarget:
    """
    This class contains off-target information for an sgRNA candidate, it is intended to be part of a
    ActivationCandidate object as an item of the off_target_list
    """
    def __init__(self, seq: str, chromosome: str, start_position: int, strand: str, number_of_mismatches: int):
        self.seq = seq
        self.chromosome = chromosome
        self.start_position = start_position
        self.strand = strand
        self.number_of_mismatches = number_of_mismatches
        self.score = -1

    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.start_position == other.start_position

    def __repr__(self):
        return f"{self.seq}, {str(self.chromosome)}, {self.start_position}, {self.strand}, {self.number_of_mismatches}, {round(float(self.score), 4)}"

    def to_dict(self, num: str):
        """Create a dictionary of the ActivationCandidate object"""

        self_dict = {f"off{num} seq": self.seq, f"off{num} chromosome": self.chromosome,
                     f"off{num} position": self.start_position, f"off{num} mms": self.number_of_mismatches,
                     f"off{num} score": self.score}
        return self_dict
