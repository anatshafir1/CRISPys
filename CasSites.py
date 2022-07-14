
__author__ = 'ItayM5'

from typing import List, Dict
import regex
import Distance_matrix_and_UPGMA


def fill_genes_targets_dict(genes_exons_dict: Dict, scoring_function, where_in_gene: float, min_length: int,
                            max_length: int, start_with_g: bool, pams: int):
    """
    creates a dictionary of genes names and the potential targets sequences for sgRNA, so that
    keys = gene names (strings) and values = lists of target sequences in the gene (list of strs).

    :param genes_exons_dict: dictionary of gene names to list of exon sequences
    :param scoring_function: the chosen scoring function for the algorithm run
    :param where_in_gene: ignore targets sites downstream to the fractional part of the gene
    :param min_length: defines the minimum length of the target sequence in the gene
    :param max_length: defines the maximum length of the target sequence in the gene
    :param start_with_g: True if the guide sequence should start with G
    :param pams: the pams by which the searching function ("get_sites") finds potential sgRNA target sites
    :return: dictionary of {keys: gene names. values: list of target sequences in gene}
    :rtype: dict
    """
    genes_targets_dict = {}
    for gene_name in genes_exons_dict.keys():
        genes_targets_dict[gene_name] = get_targets_sites_from_exons_lst(
            genes_exons_dict[gene_name],
            scoring_function,
            where_in_gene,
            min_length,
            max_length,
            start_with_g,
            pams
        )
    return genes_targets_dict


def get_targets_sites_from_exons_lst(exons_lst: List, scoring_function, where_in_gene: float = 1, min_length: int = 20,
                                     max_length: int = 20, start_with_g: bool = False, pams: int = 0) -> List:
    """Given an exons list of a single gene, this function loops through the exons sequences and finds potential target
    sites of sgRNA for each of them. The function then returns a list of all the potential target sites of the single
    gene. The function uses a complementary function ("get_sites") to find the targets in the exons, using the given pams.

    :param exons_lst: list of exons sequences of a single gene
    :param scoring_function: the chosen scoring function for the algorithm run
    :param where_in_gene: ignore targets sites downstream to the fractional part of the gene
    :param min_length: defines the minimum length of the target sequence in the gene
    :param max_length: defines the maximum length of the target sequence in the gene
    :param start_with_g: True if the guide sequence should start with G
    :param pams: the pams by which the searching function ("get_sites") finds potential sgRNA target sites
    :return: list of all the potential sgRNA target sites for the gene
    :rtype: list
    """
    original_range_in_gene = [0, where_in_gene]
    if original_range_in_gene[1] <= original_range_in_gene[0]:
        print("The range of the targets on the gene is not in the right format")
        exit(-1)
    if max_length < min_length:
        print("The range of the lengths of the sgRNA is not in the right format")
        exit(-1)
    list_of_targets = []
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
                list_of_targets += get_sites(
                    exons_lst[i][range_in_gene[0]: min(exon_lengths[i], range_in_gene[1])],
                    scoring_function,
                    min_length,
                    max_length,
                    start_with_g,
                    pams
                )
        elif max(range_in_gene[0], exon_lengths[i - 1]) < min(exon_lengths[i], range_in_gene[1]):
            list_of_targets += get_sites(
                exons_lst[i][
                    max(range_in_gene[0] - exon_lengths[i - 1], 0):
                    min(exon_lengths[i] - exon_lengths[i - 1], range_in_gene[1] - exon_lengths[i - 1])
                    ],
                scoring_function,
                min_length,
                max_length,
                start_with_g,
                pams
            )
    return list_of_targets


def get_sites(exon: str, scoring_function, min_length: int = 20, max_length: int = 20, start_with_g: bool = False,
              pams: int = 0):
    """
    This function is used to find CRISPR cut site sequences (targets) from an input DNA sequence, which is
    usually an exon. the minimum length of an input exon str is 23. Using regex this function searches for all the
    patterns of 20 or 23 (depending on the scoring function) letters long strings that a PAM sequence in their end, in
    the sense and the antisense strands of the given exon. The function then returns a list of all the found potential
    targets.

    :param exon: sequence of a gene or exon
    :param scoring_function: the chosen scoring function for the algorithm run
    :param min_length: minimum length of target site
    :param max_length: maximum length of target site
    :param start_with_g: True if the guide sequence should start with G
    :param pams: type of PAM (can be [NGG] or [NGG, NAG])
    :return: a list of targets sequences that CRISPR can target with or without the PAM
    :rtype: list
    """
    if pams == 0:
        pams = ['GG']
    elif pams == 1:
        pams = ['GG', 'AG']
    list_of_targets = []
    found_targets = []
    if len(exon) < max_length + 3:
        return list_of_targets
    for length in range(min_length, max_length + 1):
        # loop over different PAM's
        for i in range(len(pams)):
            if start_with_g:
                target_and_pam = "G" + "." * length + pams[i]
            else:
                target_and_pam = "." * (length + 1) + pams[i]
            compiled = regex.compile(target_and_pam)
            found_sense_targets = regex.findall(compiled, exon, overlapped=True)
            found_antisense_targets = regex.findall(compiled, give_complementary(exon), overlapped=True)
            # functions that take targets of length 23 (used for scoring scheme that needs the PAM, e.g. gold_off)
            if scoring_function == Distance_matrix_and_UPGMA.gold_off_func:
                found_targets = [seq for seq in found_sense_targets if 'N' not in seq]
                found_targets += [seq for seq in found_antisense_targets if 'N' not in seq]
            # functions that take targets of length 20
            else:
                found_targets = [seq[:-3] for seq in found_sense_targets if 'N' not in seq[:-3]]
                found_targets += [seq[:-3] for seq in found_antisense_targets if 'N' not in seq[:-3]]
            list_of_targets += found_targets
    return list_of_targets


def give_complementary(seq: str):
    """Given a DNA sequence (5' to 3') this function returns its antisense sequence (also 5' to 3'). This is used to
    find possible cut sites for CRISPR in the antisense strand of a given DNA sequence.

    :param seq: input DNA sequence
    :return: antisense sequence for the input
    :rtype: str
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
