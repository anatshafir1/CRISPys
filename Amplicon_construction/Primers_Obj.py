class Primers_Obj:
    """A class representing a pair of primers given by primer3 algorithm and their parameters.
        a toy example:



        """

    def __init__(self, primer_penalty: float, left_sequence: str, right_sequence: str, left_start_idx: int,
                 left_tm: float, right_start_idx: int, right_tm: float):
        self.primer_penalty = primer_penalty
        self.left_sequence = left_sequence
        self.right_sequence = right_sequence
        self.left_start_idx = left_start_idx
        self.left_tm = left_tm
        self.right_start_idx = right_start_idx
        self.right_tm = right_tm

    def __str__(self):
        return f"penalty: {self.primer_penalty}, sequences: <{self.left_sequence},{self.right_sequence}>, " \
               f"left idx&tm: <{self.left_start_idx},{self.left_tm}>, right idx&tm: <{self.right_start_idx},{self.right_tm}>"

    def __repr__(self):
        return self.__str__()

    def to_dict(self):
        return {"penalty": self.primer_penalty, "primer_left": self.left_sequence, "tm_left": self.left_tm,
                "primer_right": self.right_sequence, "tm_right": self.right_tm}
