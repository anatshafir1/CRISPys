"""sgRNA target class file"""
from typing import List, Dict, Set


class Target_Obj:
    """A class representing an sgRNA target for gene editing.
        a toy example:

        "ACTGACTGACTGACTGACTGAGG", 100, 122, "+"

        """

    def __init__(self, seq: str, start_idx: int, end_idx: int, strand: str, scaffold="", ungapped_seq="", rank: float = 0,
                 score: int = 0, exon_num: int = 0):
        self.seq = seq
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.strand = strand
        self.scaffold = scaffold
        self.ungapped_seq = ungapped_seq
        self.rank = rank
        self.score = score
        self.exon_num = exon_num

    def __str__(self):
        return f"{self.rank};{self.start_idx};{self.end_idx};{self.strand}"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.end_idx - self.start_idx + 1

    def __eq__(self, other):
        return self.seq == other.seq and self.start_idx == other.start_idx and self.scaffold == other.scaffold

    def to_dict(self):
        return {"gRNA+PAM": self.seq, "gRNA_start": self.start_idx, "gRNA_end": self.end_idx,
                "gRNA_strand": self.strand}


class Family_Target_Obj(Target_Obj):
    """A class representing an sgRNA target for gene editing.
        a toy example:

        "ACTGACTGACTGACTGACTGAGG", 100, 122, "+"

        """

    def __init__(self, seq: str, start_idx: int, end_idx: int, strand: str, scaffold="", ungapped_seq="", rank: float = 0,
                 score: int = 0, cut_alleles_lst: List = None, exon_num: int = 0, sgrna: str = ""):
        super().__init__(seq, start_idx, end_idx, strand, scaffold, ungapped_seq, rank, score, exon_num)
        self.cut_alleles_lst = cut_alleles_lst
        self.sgrna = sgrna

    def to_family_singleplex_dict(self, family_targeting: int):
        if family_targeting == 1:
            return {"target+PAM": self.seq, "target_start": self.start_idx,
                    "target_end": self.end_idx, "target_strand": self.strand}
        elif family_targeting == 2:
            if self.rank % 0.2 == 0:
                return {"up_target+PAM": "", "up_target_start": "",
                        "up_target_end": "", "up_target_strand": "",
                        "down_target+PAM": self.seq, "down_target_start": self.start_idx,
                        "down_target_end": self.end_idx, "down_target_strand": self.strand}
            else:
                return {"up_target+PAM": self.seq, "up_target_start": self.start_idx,
                        "up_target_end": self.end_idx, "up_target_strand": self.strand,
                        "down_target+PAM": "", "down_target_start": "",
                        "down_target_end": "", "down_target_strand": ""}


class Combined_Target_Obj:

    def __init__(self, start_idx: int, end_idx: int, targets_list: List[Target_Obj] = None, sg_perm: List[str] = None,
                 offscores_dict: Dict[str, Dict[str, float]] = None, cut_alleles: Set[str] = None, chosen_sg: str = "",
                 chosen_sg_score: float = 0.0):

        self.start_idx = start_idx
        self.end_idx = end_idx
        self.targets_list = targets_list
        self.sg_perm = sg_perm
        self.offscores_dict = offscores_dict  # {sgRNA: {Scaffold: score}}
        self.cut_alleles = cut_alleles
        self.chosen_sg = chosen_sg
        self.chosen_sg_score = chosen_sg_score

    def __str__(self):
        return f"{self.start_idx};{self.chosen_sg};{self.cut_alleles}"

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
                 multiplex_score: float, up_targets_list: List[Target_Obj], down_targets_list: List[Target_Obj], exon_num: int, rank: int):

        self.up_start = up_start
        self.up_end = up_end
        self.up_seq = up_seq
        self.down_start = down_start
        self.down_end = down_end
        self.down_seq = down_seq
        self.multiplex_score = multiplex_score
        self.up_targets_list = up_targets_list
        self.down_targets_list = down_targets_list
        self.exon_num = exon_num
        self.start_idx = up_start
        self.end_idx = down_end
        self.rank = rank

    def __str__(self):
        return f"{self.rank};{round(self.multiplex_score, 4)};{self.up_start};{self.up_seq};{self.down_start};{self.down_seq}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (self.up_start == other.up_start and self.up_seq == other.up_seq and self.down_start == other.down_start
                and self.down_seq == other.down_seq)

    def __hash__(self):
        return hash(self.__str__())

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
        return {"up_target+PAM": up_sgandpam, "up_target_start": self.up_start,
                "up_target_end": self.up_end, "up_target_strand": up_target_strand, "down_target+PAM": down_sgandpam,
                "down_target_start": self.down_start, "down_target_end": self.down_end,
                "down_target_strand": down_target_strand, "multiplex_score": score}


class FamilyMultiplexTarget(MultiplexTarget):

    def __init__(self, up_start: int, up_end: int, up_seq: str, down_start: int, down_end: int, down_seq: str,
                 multiplex_score: float, up_targets_list: List, down_targets_list: List, exon_num: int, rank: int,
                 up_target_strand, down_target_strand, up_sgrna, down_sgrna):
        """

        :param up_start:
        :param up_end:
        :param up_seq: upstream target + PAM sequence
        :param down_start:
        :param down_end:
        :param down_seq: downstream target + PAM sequence
        :param multiplex_score:
        :param up_targets_list:
        :param down_targets_list:
        :param exon_num:
        :param rank:
        :param up_target_strand:
        :param down_target_strand:
        :param up_sgrna:
        :param down_sgrna:
        """
        super().__init__(up_start, up_end, up_seq, down_start, down_end, down_seq, multiplex_score, up_targets_list,
                         down_targets_list, exon_num, rank)

        self.up_target_strand = up_target_strand
        self.down_target_strand = down_target_strand
        self.up_sgrna = up_sgrna
        self.down_sgrna = down_sgrna

    def to_family_multiplex_dict(self):

        return {"up_target+PAM": self.up_seq, "up_target_start": self.up_start,
                "up_target_end": self.up_end, "up_target_strand": self.up_target_strand,
                "down_target+PAM": self.down_seq, "down_target_start": self.down_start,
                "down_target_end": self.down_end, "down_target_strand": self.down_target_strand}


class sgRNA:
    def __init__(self, start: int, end: int, seq: str, score_dict: Dict[str, float], targets_list: List[Target_Obj]):

        self.start = start
        self.end = end
        self.seq = seq
        self.score_dict = score_dict
        self.targets_list = targets_list

    def __str__(self):
        return f"{self.start}, {self.seq}"

    def __repr__(self):
        return self.__str__()
