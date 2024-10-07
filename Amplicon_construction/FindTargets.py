"""Extracting target sequences from upstream sites"""
import regex
from typing import List, Tuple, Dict

from Target_Obj import Target_Obj
from FindSNPTargets import get_snp_targets


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


def find_targets_in_sequence(exon_region: str, pams: Tuple, max_amplicon_len: int, primer_length: int, cut_location: int,
                             target_surrounding_region: int, target_len: int) -> List[Target_Obj]:
    """
    This function is used to find CRISPR target site sequences from an input DNA sequence. Using regex this
    function searches for all the patterns of 23 letters long strings with all the PAM sequences in 'pams' in their end,
    in the sense and the antisense strands of the input DNA sequence. The function then returns a list of all the found
    potential targets.

    :param exon_region: DNA sequence of a gene TSS upstream site
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param primer_length:
    :param max_amplicon_len:
    :param pams: type of PAM
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :return: a list of targets sequences that CRISPR can target as tuples of target sequence, target start, target end, strand (relative to given exon_region sequence)
    """

    found_targets = []

    exon_region_len = len(exon_region)
    # calculate range of intron sites that are allowed be used to construct an amplicon, but where sgRNA target sites are now allowed
    intron_region_added = max_amplicon_len - primer_length - cut_location - target_surrounding_region  # 255 by default
    target_allowed_start_idx = intron_region_added - target_len + cut_location  # 239 by default
    target_allowed_end_idx = exon_region_len - intron_region_added
    allowed_exon_region_for_targets = exon_region[target_allowed_start_idx: target_allowed_end_idx]

    complementary_strand = give_complementary(exon_region[intron_region_added: exon_region_len - target_allowed_start_idx])  # create complementary strand to search for targets on it

    pam_len = len(pams[0])
    target_without_pam_len = target_len - pam_len
    # loop over different PAM's
    for i in range(len(pams)):
        target_and_pam = "." * target_without_pam_len + pams[i]
        compiled = regex.compile(target_and_pam)
        found_sense_targets = regex.finditer(compiled, allowed_exon_region_for_targets)
        found_antisense_targets = regex.finditer(compiled, complementary_strand)

        found_targets += [Target_Obj(seq.group(0), target_allowed_start_idx + seq.start(), target_allowed_start_idx + seq.end() - 1, "+") for seq in found_sense_targets if 'N' not in seq.group(0)]
        found_targets += [Target_Obj(seq.group(0), exon_region_len - target_allowed_start_idx - seq.end(), exon_region_len - target_allowed_start_idx - seq.start() - 1, "-") for seq in found_antisense_targets if 'N' not in seq.group(0)]
    return sorted(found_targets, key=lambda target: target.start_idx)


def filter_duplicates(exon_targets):
    target_pos_dict = {target.start_idx: target for target in exon_targets}
    return sorted(list(target_pos_dict.values()), key=lambda target: target.start_idx)


def get_targets(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], pams: Tuple, max_amplicon_len: int, primer_length: int, cut_location: int,
                target_surrounding_region: int, target_len: int, k: int) -> Dict[int, List]:
    """

    :param gene_sequences_dict: dictionary of exon num -> list of tuples representing alleles where tuple[0] is scaffold name 
    (example format: ">scaffold10132:437703-438762(+)") and tuple[1] is allele sequence.
    :param pams:
    :param max_amplicon_len: maximum length of the amplicon
    :param primer_length: minimum length of the primer sequence
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :param k:
    :return:
    """
    targets_dict = {}
    if k > 0:
        targets_dict = get_snp_targets(gene_sequences_dict, pams, max_amplicon_len, primer_length, cut_location, target_surrounding_region, target_len, k)
        return targets_dict
    else:
        for exon_region in gene_sequences_dict:
            exon_region_seq = gene_sequences_dict[exon_region][0][1]
            exon_targets = find_targets_in_sequence(exon_region_seq, pams, max_amplicon_len, primer_length, cut_location,
                                                    target_surrounding_region, target_len)
            unique_exon_targets = filter_duplicates(exon_targets)
            targets_dict[exon_region] = unique_exon_targets
        return targets_dict
