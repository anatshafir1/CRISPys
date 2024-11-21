"""sgRNA target class file"""
from typing import List, Dict, Set


class Target_Obj:
    """A class representing an sgRNA target for gene editing.
        a toy example:

        "ACTGACTGACTGACTGACTGAGG", 100, 122, "+"

        """

    def __init__(self, seq: str, start_idx: int, end_idx: int, strand: str, scaffold="", ungapped_seq=""):
        self.seq = seq
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.strand = strand
        self.scaffold = scaffold
        self.ungapped_seq = ungapped_seq

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

    def __init__(self, start_idx: int, end_idx: int, targets_list: List[Target_Obj] = [], sg_perm: List[str] = [],
                 offscores_dict: Dict[str, Dict[str, float]] = {}, cut_alleles: Set[str] = set(), chosen_sg: str = "",
                 chosen_sg_score: float = 0.0):

        self.start_idx = start_idx
        self.end_idx = end_idx
        self.targets_list = targets_list
        self.sg_perm = sg_perm
        self.offscores_dict = offscores_dict
        self.cut_alleles = cut_alleles
        self.chosen_sg = chosen_sg
        self.chosen_sg_score = chosen_sg_score

    def to_dict(self, scaffold: str, strand: str):
        pam = ""
        target_scaffold = ""
        for target in self.targets_list:
            if target.scaffold == list(self.cut_alleles)[0]:
                pam = target.seq[20:23]
            if target.scaffold == scaffold:
                target_scaffold = "+" if target.strand == strand else "-"

        score = self.offscores_dict[self.chosen_sg][scaffold]
        sgandpam = self.chosen_sg + pam if scaffold in self.cut_alleles else "NA"
        return {"gRNA+PAM": sgandpam, "MOFF-score": score, "gRNA_start": self.start_idx,
                "gRNA_end": self.end_idx, "gRNA_strand": target_scaffold}

    def __str__(self):
        return f"{self.start_idx}, {self.chosen_sg}, {self.cut_alleles}"

    def __repr__(self):
        return self.__str__()
