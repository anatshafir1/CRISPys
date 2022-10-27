"""A module for finding potential genomic targets in genes for sgRNA cleaving"""
__author__ = 'ItayM5'

from typing import List, Dict
import regex


def fill_genes_targets_dict(genes_exons_dict: Dict, pam_included, where_in_gene: float, start_with_g: bool,
                            pams: int) -> Dict[str, List[str]]:
    """
    creates a dictionary of genes names and the potential targets sequences for sgRNA, so that
    keys = gene names (strings) and values = lists of target sequences in the gene (list of strs).

    :param genes_exons_dict: dictionary of gene names to list of exon sequences
    :param pam_included: whether PAMs are included in target sequence
    :param where_in_gene: ignore targets sites downstream to the fractional part of the gene
    :param start_with_g: True if the guide sequence should start with G
    :param pams: the pams by which the searching function ("get_sites") finds potential sgRNA target sites
    :return: dictionary of {keys: gene names. values: list of target sequences in gene}
    """
    genes_targets_dict = {}
    genes_targets_pos_dict = {}
    for gene_name in genes_exons_dict.keys():
        genes_targets_dict[gene_name] = get_targets_sites_from_exons_lst(
            genes_exons_dict[gene_name],
            pam_included,
            where_in_gene,
            start_with_g,
            pams
        )[0]

        genes_targets_pos_dict[gene_name] = get_targets_sites_from_exons_lst(
            genes_exons_dict[gene_name],
            pam_included,
            where_in_gene,
            start_with_g,
            pams
        )[1]

    return genes_targets_dict, genes_targets_pos_dict


def get_targets_sites_from_exons_lst(exons_lst: List, pam_included, where_in_gene: float = 1,
                                     start_with_g: bool = False, pams: int = 0) -> List[str]:
    """Given an exons list of a single gene, this function loops through the exons sequences and finds potential target
    sites of sgRNA for each of them. The function then returns a list of all the potential target sites of the single
    gene. The function uses a complementary function ("get_sites") to find the targets in the exons, using the given pams.

    :param exons_lst: list of exons sequences of a single gene
    :param pam_included: whether PAMs are included in target sequence
    :param where_in_gene: ignore targets sites downstream to the fractional part of the gene
    :param start_with_g: True if the guide sequence should start with G
    :param pams: the pams by which the searching function ("get_sites") finds potential sgRNA target sites
    :return: 2 lists, a list of all the potential sgRNA target sites for the gene and a list of all the potential sgRNA target sites for the gene with the position strans and pam
    """
    original_range_in_gene = [0, where_in_gene]
    if original_range_in_gene[1] <= original_range_in_gene[0]:
        print("The range of the targets on the gene is not in the right format")
        exit(-1)
    list_of_targets = []
    list_of_targets_with_pos = []
    exon_lengths = [len(exon) for exon in exons_lst]
    gene_length = sum(exon_lengths)
    # range in gene - used for deciding what parts to consider in the gene (multiply the 'where_in_gene' argument by
    # the length of the gene)
    range_in_gene = [int(r * gene_length) for r in original_range_in_gene]
    # accumulate the length of exons
    for i in range(1, len(exon_lengths)):
        exon_lengths[i] = exon_lengths[i - 1] + exon_lengths[i]
    # loop on each exon
    for i in range(len(exons_lst)):
        # if there is only one exon check target
        if i == 0:
            if range_in_gene[0] < exon_lengths[i]:
                targets_tpl = get_sites(
                    exons_lst[i][range_in_gene[0]: min(exon_lengths[i], range_in_gene[1])],
                    pam_included,
                    exons_lst,
                    start_with_g,
                    pams
                )
                list_of_targets += targets_tpl[0]
                list_of_targets_with_pos += targets_tpl[1]
        elif max(range_in_gene[0], exon_lengths[i - 1]) < min(exon_lengths[i], range_in_gene[1]):
            targets_tpl = get_sites(
                exons_lst[i][
                    max(range_in_gene[0] - exon_lengths[i - 1], 0):
                    min(exon_lengths[i] - exon_lengths[i - 1], range_in_gene[1] - exon_lengths[i - 1])
                    ],
                pam_included,
                exons_lst,
                start_with_g,
                pams
            )
            list_of_targets += targets_tpl[0]
            list_of_targets_with_pos += targets_tpl[1]

    return  list_of_targets, list_of_targets_with_pos


def get_sites(exon: str, pam_included, exons_lst: list, start_with_g: bool = False, pams: int = 0):
    """
    This function is used to find CRISPR cut site sequences (targets) from an input DNA sequence, which is
    usually an exon. the minimum length of an input exon str is 23. Using regex this function searches for all the
    patterns of 23 (depending on the scoring function) letters long strings with a PAM sequence in their end, in
    the sense and the antisense strands of the given exon. The function then returns a list of all the found potential
    targets with length of 20 or 23 (depending on the scoring function).

    :param exon: sequence of a gene or exon
    :param pam_included: whether PAMs are included in target sequence
    :param start_with_g: True if the guide sequence should start with G
    :param pams: type of PAM (can be [NGG] or [NGG, NAG])
    :return: a list of targets sequences that CRISPR can target with or without the PAM
    """
    target_len = 20
    if pams == 0:
        pams = ['GG']
    elif pams == 1:
        pams = ['GG', 'AG']
    list_of_targets = []
    list_of_targets_and_pos = []
    if len(exon) < target_len + 3:
        return list_of_targets
        # loop over different PAM's
    for i in range(len(pams)):
        if start_with_g:
            target_and_pam = "G" + "." * target_len + pams[i]
        else:
            target_and_pam = "." * (target_len + 1) + pams[i]
        compiled = regex.compile(target_and_pam)
        # find all sense target
        found_sense_targets = regex.findall(compiled, exon, overlapped=True)
        # find all anti-sense target
        found_antisense_targets = regex.findall(compiled, give_complementary(exon), overlapped=True)

        # add position, pam and strand to sense results
        sense_targets_and_pos = add_positions_and_strand(found_sense_targets, exons_lst, "+")
        # add position, pam and strand to anti-sense results
        antisense_targets_and_pos = add_positions_and_strand(found_antisense_targets, exons_lst, "-")

        target_with_pos = sense_targets_and_pos
        target_with_pos += antisense_targets_and_pos
        # functions that take targets of length 23 (used for scoring scheme that needs the PAM, e.g. gold_off)
        if pam_included:
            found_targets = [seq for seq in found_sense_targets if 'N' not in seq]
            found_targets += [seq for seq in found_antisense_targets if 'N' not in seq]

        # functions that take targets of length 20
        else:
            found_targets = [seq[:-3] for seq in found_sense_targets if 'N' not in seq[:-3]]
            found_targets += [seq[:-3] for seq in found_antisense_targets if 'N' not in seq[:-3]]
        list_of_targets += found_targets
        list_of_targets_and_pos += target_with_pos
    return list_of_targets, list_of_targets_and_pos

def add_positions_and_strand(targets_lst, exons_lst, strand):
    """
    This function was added by me (Udi) in order to get the position and strand of each target and add
     it to the output of crispys
    Args:
        targets_lst: list of target sequence
        exons_lst: list of exons
        strand: the sense or anti sense you want to search for
    it returns list of tuples with (target, pam_seq, pos, strand)
    """
    target_pos_lst = []

    seq = "".join(exons_lst)
    if strand == '-':
        seq = give_complementary(seq)
    for target in targets_lst:
        pos = seq.find(target)
        target_pos_lst.append((target[0:20], target[20:23], pos, strand))
    return target_pos_lst


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
