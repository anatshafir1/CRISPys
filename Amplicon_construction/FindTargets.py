"""Extracting target sequences from upstream sites"""
import regex
from typing import List, Tuple, Dict
from Crispys import globals


def give_complementary(seq: str) -> str:
    """Given a DNA sequence (5' to 3') this function returns its antisense sequence (also 5' to 3'). This is used to
    find possible cut sites for CRISPR in the antisense strand of a given DNA sequence.
    :param seq: input DNA sequence
    :return: antisense sequence for the input
    """
    complementary_seq_list = []
    for i in range(len(seq)):
        if seq[len(seq) - 1 - i] == 'A':
            complementary_seq_list.append('T')
        elif seq[len(seq) - 1 - i] == 'T':
            complementary_seq_list.append('A')
        elif seq[len(seq) - 1 - i] == 'C':
            complementary_seq_list.append('G')
        elif seq[len(seq) - 1 - i] == 'G':
            complementary_seq_list.append('C')
        elif seq[len(seq) - 1 - i] == 'N':
            complementary_seq_list.append('N')
    return ''.join(complementary_seq_list)


def find_exon_targets(exon: str, pams: Tuple, max_amplicon_len: int, primer_length: int, target_len: int) -> List[Tuple[str, int]]:
    """
    This function is used to find CRISPR target site sequences from an input DNA sequence. Using regex this
    function searches for all the patterns of 23 letters long strings with all the PAM sequences in 'pams' in their end,
    in the sense and the antisense strands of the input DNA sequence. The function then returns a list of all the found
    potential targets.

    :param exon: DNA sequence of a gene TSS upstream site
    :param pams: type of PAM
    :return: a list of targets sequences that CRISPR can target
    """
    pam_len = len(pams[0])
    target_without_pam_len = target_len - pam_len
    found_targets = []
    exon_len = len(exon)
    complementary_strand = give_complementary(exon)  # create complementary strand to search for targets on it
    exon_region_add = max_amplicon_len - primer_length - target_len - globals.safety_padding_around_target * 2
    # loop over different PAM's
    for i in range(len(pams)):
        target_and_pam = "." * target_without_pam_len + pams[i]
        compiled = regex.compile(target_and_pam)
        found_sense_targets = regex.finditer(compiled, exon)
        found_antisense_targets = regex.finditer(compiled, complementary_strand)

        found_targets += [(seq.group(0), exon_region_add + seq.start(), exon_region_add + seq.end(), "+") for seq in found_sense_targets if 'N' not in seq.group(0)]
        found_targets += [(seq.group(0), exon_region_add + exon_len - seq.end(), exon_region_add + exon_len - seq.start(), "-") for seq in found_antisense_targets if 'N' not in seq.group(0)]
    return sorted(found_targets, key=lambda x: x[1])


def get_targets(gene_sequences_dict: Dict, pams: Tuple, max_amplicon_len: int, primer_length: int, target_len: int) -> Dict[int, List[Tuple[str, int]]]:

    targets_dict = {}
    exon_start_idx = max_amplicon_len - primer_length - target_len - globals.safety_padding_around_target * 2
    exon_end_idx = exon_start_idx + max_amplicon_len
    for exon_region in gene_sequences_dict:
        only_exon = gene_sequences_dict[exon_region][1][exon_start_idx:exon_end_idx]
        exon_targets_tuple = find_exon_targets(only_exon, pams, max_amplicon_len, primer_length, target_len)
        targets_dict[exon_region] = exon_targets_tuple
    return targets_dict
