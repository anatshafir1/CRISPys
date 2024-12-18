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

    def __eq__(self, other):
        return self.seq == other.seq and self.start_idx == other.start_idx and self.scaffold == other.scaffold

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

    def __str__(self):
        return f"{self.start_idx}, {self.chosen_sg}, {self.cut_alleles}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return self.chosen_sg == other.chosen_sg and self.start_idx == other.start_idx

    def to_dict(self, scaffold: str, strand: str):
        pam = ""
        target_strand = ""
        for target in self.targets_list:
            if target.scaffold == list(self.cut_alleles)[0]:
                pam = target.seq[20:23]
            if target.scaffold == scaffold:
                target_strand = "+" if target.strand == strand else "-"

        score = self.offscores_dict[self.chosen_sg][scaffold]
        sgandpam = self.chosen_sg + pam if scaffold in self.cut_alleles else "NA"
        return {"gRNA+PAM": sgandpam, "MOFF-score": score, "gRNA_start": self.start_idx,
                "gRNA_end": self.end_idx, "gRNA_strand": target_strand}


class MultiplexTarget:

    def __init__(self, up_start: int, up_end: int, up_seq: str, down_start: int, down_end: int, down_seq: str,
                 multiplex_score: float, up_targets_list: List[Target_Obj], down_targets_list: List[Target_Obj]):

        self.up_start = up_start
        self.up_end = up_end
        self.up_seq = up_seq
        self.down_start = down_start
        self.down_end = down_end
        self.down_seq = down_seq
        self.multiplex_score = multiplex_score
        self.up_targets_list = up_targets_list
        self.down_targets_list = down_targets_list

    def __str__(self):
        return f"{self.up_start}, {self.up_seq}, {self.down_start}, {self.down_seq}, {self.multiplex_score}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return self.up_start == other.up_start and self.up_seq == other.up_seq and self.down_start == other.down_start and self.down_seq == other.down_seq

    def to_dict(self, scaffold: str, strand: str):
        up_pam = ""
        down_pam = ""
        up_target_strand = ""
        down_target_strand = ""
        for target in self.up_targets_list:
            if target.scaffold == scaffold:
                up_pam = target.seq[20:23]
                up_target_strand = "+" if target.strand == strand else "-"
        for target in self.down_targets_list:
            if target.scaffold == scaffold:
                down_pam = target.seq[20:23]
                down_target_strand = "+" if target.strand == strand else "-"
        score = self.multiplex_score
        up_sgandpam = self.up_seq + up_pam
        down_sgandpam = self.down_seq + down_pam
        return {"up_gRNA+PAM": up_sgandpam, "up_gRNA_start": self.up_start,
                "up_gRNA_end": self.up_end, "up_gRNA_strand": up_target_strand, "down_gRNA+PAM": down_sgandpam,
                "down_gRNA_start": self.down_start, "down_gRNA_end": self.down_end,
                "down_gRNA_strand": down_target_strand, "multiplex_score": score}
