from itertools import combinations, groupby
from operator import attrgetter
from typing import Dict, List, Tuple

from Amplicon_construction.FindTargets import get_all_targets, create_comb_targets, create_sgrna_permutations, \
    calculate_off_scores
from Amplicon_construction.Target_Obj import MultiplexTarget, Combined_Target_Obj, sgRNA


def calc_multiplex_score(upstream_target: sgRNA, downstream_target: sgRNA, allele_ids_lst: List[str]) -> float:
    """

    :param upstream_target: potential upstream sgRNA
    :param downstream_target: potential upstream sgRNA
    :param allele_ids_lst: list of all allele IDs
    :return: multiplex score of upstream and downstream sgRNAs
    """
    tot_score = 1.0
    for allele in allele_ids_lst:
        upstream_sg_no_cut = 1 - 0.9*upstream_target.score_dict[allele]
        downstream_sg_no_cut = 1 - 0.9*downstream_target.score_dict[allele]
        allele_cut = 1 - (upstream_sg_no_cut*downstream_sg_no_cut)
        tot_score *= allele_cut
    return tot_score


def create_sgrna_pairs(comb_targets_lst: List[Combined_Target_Obj]) -> List[Tuple[sgRNA, sgRNA]]:
    """

    :param comb_targets_lst: List of Combined_Target_Obj objects
    :return: List of all potential pairs of sgRNAs extracted from the comb_targets_lst
    """
    sg_list = []
    for comb_target in comb_targets_lst:
        for sg in comb_target.sg_perm:
            sg_list.append(sgRNA(comb_target.start_idx, comb_target.end_idx, sg, comb_target.offscores_dict[sg], comb_target.targets_list))

    sg_pairs_list = [(sg1, sg2) for sg1, sg2 in combinations(sg_list, 2) if sg1.start != sg2.start]

    return sg_pairs_list


def create_multiplex_targets(targets_dict: Dict[int, List[Combined_Target_Obj]],
                             allele_ids_lst: List[str]) -> List[MultiplexTarget]:
    """
    Create combinations of targets using Combined Targets and their chosen sgRNA sequence and score. Create multiplex
    targets from the combinations.

    :param targets_dict: dictionary of exon number -> List of targets as Combined_Target_Obj
    :param allele_ids_lst: list of all allele IDs
    :return: dictionary of exon number -> List of targets as MultiplexTarget
    """
    all_multiplex_targets_list = []
    multiplex_targets_dict = {}
    for exon in targets_dict:
        multiplex_targets_lst = []
        comb_targets_lst = targets_dict[exon]
        target_pairs = create_sgrna_pairs(comb_targets_lst)
        for target_pair in target_pairs:
            # distinguish which target is upstream and which is downstream
            if target_pair[0].start < target_pair[1].start:
                upstream_target = target_pair[0]
                downstream_target = target_pair[1]
            else:
                upstream_target = target_pair[1]
                downstream_target = target_pair[0]

            multiplex_score = calc_multiplex_score(upstream_target, downstream_target, allele_ids_lst)
            multiplex_target = MultiplexTarget(upstream_target.start, upstream_target.end, upstream_target.seq,
                                               downstream_target.start, downstream_target.end, downstream_target.seq,
                                               multiplex_score, upstream_target.targets_list, downstream_target.targets_list, exon)
            multiplex_targets_lst.append(multiplex_target)

        score_sorted_targets = sorted(multiplex_targets_lst, key=lambda tg: (tg.up_start, -tg.multiplex_score))
        # for every start index take only pair with the highest score:
        max_score_sorted_targets = [next(group) for _, group in groupby(score_sorted_targets, key=attrgetter('up_start'))]

        all_multiplex_targets_list.extend(max_score_sorted_targets)
    # store multiplex targets in exons dict:
    sorted_all_targets_list = sorted(all_multiplex_targets_list, key=lambda tg: -tg.multiplex_score)
    return sorted_all_targets_list


def get_multiplex_targets(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], pams: Tuple, max_amplicon_len: int,
                primer_length: int, cut_location: int, target_surrounding_region: int, target_len: int,
                          distinct_alleles_num: int) -> List[MultiplexTarget]:
    """

    :param gene_sequences_dict: dictionary of exon num -> list of tuples representing alleles where tuple[0] is scaffold name
    (example format: ">scaffold10132:437703-438762(+)") and tuple[1] is allele sequence.
    :param pams: tuple of PAM sequences of the Cas protein in use
    :param max_amplicon_len: maximum length of the amplicon
    :param primer_length: minimum length of the primer sequence
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :return: a dictionary of exon number -> List of targets as Target_Obj or Combined_Target_Obj, depending on the tool in use.
    """
    print("searching for potential targets".upper().center(40, "#"))
    allele_ids_lst = [gene_sequences_dict[1][i][0].split(":")[0][1:] for i in range(distinct_alleles_num)]
    all_targets_dict = get_all_targets(gene_sequences_dict, pams, max_amplicon_len, primer_length, cut_location,
                                       target_surrounding_region, target_len)
    relevant_targets_dict = create_comb_targets(all_targets_dict, target_len)
    create_sgrna_permutations(relevant_targets_dict)
    calculate_off_scores(relevant_targets_dict)
    all_multiplex_targets_list = create_multiplex_targets(relevant_targets_dict, allele_ids_lst)
    return all_multiplex_targets_list


def get_multiplex_targets_dict(multiplex_target: MultiplexTarget,
                               gene_exon_seqs_dict: Dict[int, List[Tuple[str, str]]]) -> Dict[int, List]:
    multiplex_targets_dict = {}
    for exon in gene_exon_seqs_dict:
        if exon == multiplex_target.exon_num:
            multiplex_targets_dict[exon] = [multiplex_target]
        else:
            multiplex_targets_dict = []
    return multiplex_targets_dict
