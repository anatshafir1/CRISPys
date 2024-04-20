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
        self.length = end_idx - start_idx + 1

    def __str__(self):
        return f"{self.seq}, {self.start_idx}, {self.end_idx}, {self.strand}"

    def __repr__(self):
        return self.__str__()

    def to_dict(self):
        return {"target_seq": self.seq, "target_start_idx": self.start_idx, "target_end_idx": self.end_idx,
                "target_strand": self.strand, "target_len": self.length}
