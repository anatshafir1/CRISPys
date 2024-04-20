class Primers_Obj:
    """A class representing a pair of primers given by primer3 algorithm and their parameters.
        a toy example:



        """

    def __init__(self, primer_penalty: float, left_sequence: str, right_sequence: str, left_start_idx: int,
                 left_length: int,
                 right_start_idx: int, right_length: int):
        self.primer_penalty = primer_penalty
        self.left_sequence = left_sequence
        self.right_sequence = right_sequence
        self.left_start_idx = left_start_idx
        self.left_length = left_length
        self.right_start_idx = right_start_idx
        self.right_length = right_length

    def __str__(self):
        return f"penalty: {self.primer_penalty}, sequences: <{self.left_sequence},{self.right_sequence}>, " \
               f"left idx&length: <{self.left_start_idx},{self.left_length}>, right idx&length: <{self.right_start_idx},{self.right_length}>"

    def __repr__(self):
        return self.__str__()

    def to_dict(self):
        return self.__dict__
