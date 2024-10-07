"""sgRNA target class file"""
from typing import List, Dict, Set


class Target_Obj:
    """A class representing an sgRNA target for gene editing.
        a toy example:

        "ACTGACTGACTGACTGACTGAGG", 100, 122, "+"

        """

    def __init__(self, seq: str, start_idx: int, end_idx: int, strand: str, scaffold=""):
        self.seq = seq
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.strand = strand
        self.scaffold = scaffold

    def __str__(self):
        return f"{self.seq}, {self.start_idx}, {self.end_idx}, {self.strand}"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.end_idx - self.start_idx + 1

    def to_dict(self):
        return {"gRNA+PAM": self.seq, "gRNA_start": self.start_idx, "gRNA_end": self.end_idx,
                "gRNA_strand": self.strand}


class Combined_Target_Obj:

    def __init__(self, start_idx: int, targets_list: List[Target_Obj] = [], snp_dict: Dict[int, Dict[str, Set]] = {}, k_alleles_list: List[Set[str]] = []):

        self.start_idx = start_idx
        self.targets_list = targets_list
        self.snp_dict = snp_dict
        self.k_alleles_list = k_alleles_list

    def __str__(self):
        return f"{self.start_idx}, {self.k_alleles_list}"

    def __repr__(self):
        return self.__str__()
