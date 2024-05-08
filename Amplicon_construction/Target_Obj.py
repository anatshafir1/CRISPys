"""sgRNA target class file"""


class Target_Obj:
    """A class representing an sgRNA target for gene editing.
        a toy example:

        "ACTGACTGACTGACTGACTGAGG", 100, 122, "+"

        """

    def __init__(self, seq: str, start_idx: int, end_idx: int, strand: str):
        self.seq = seq
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.strand = strand

    def __str__(self):
        return f"{self.seq}, {self.start_idx}, {self.end_idx}, {self.strand}"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.end_idx - self.start_idx + 1

    def to_dict(self):
        return {"gRNA+PAM": self.seq, "gRNA_start": self.start_idx, "gRNA_end": self.end_idx,
                "gRNA_strand": self.strand}
