
from typing import List, Tuple
from SubgroupRes import SubgroupRes
from Candidate import Candidate


def candidate_min_mismatch_sum(candidate: Candidate, omega: float, genes_list: List = None) -> int:
    """
    Calculates the sum of minimum number of mismatches between a candidate (sgRNA) and the genes it cleaves with
    probability higher than omega

    :param candidate: current sgRNA candidate
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param genes_list:
    :return:
    """
    if not genes_list:
        genes_list = [gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] >= omega]
    res = 0
    for gene in genes_list:
        targets_lst = candidate.targets_dict.get(gene)
        res += min(list(map(lambda sub_lst: len(sub_lst[1]), targets_lst)))
    return res


def cover_above_thr(candidate: Candidate, omega: float, genes_list: List = None) -> Tuple[float, int, List[str]]:
    """

    :param candidate: current sgRNA candidate
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param genes_list:
    :return:
    """
    cover_prob, num_of_covered = 1, 0
    if not genes_list:
        genes_list = [gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] >= omega]
    for gene in genes_list:
        if candidate.genes_score_dict.get(gene, 0) >= omega:
            cover_prob *= candidate.genes_score_dict[gene]
            num_of_covered += 1
    return cover_prob, num_of_covered, genes_list


def contained(candidate_i: Candidate, candidate_j: Candidate, use_thr: int, omega: float) -> bool:
    """

    :param candidate_i: current sgRNA candidate
    :param candidate_j: current sgRNA candidate
    :param use_thr:
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :return:
    """
    if use_thr:
        cover_prob_i, num_of_covered_i, genes_list = cover_above_thr(candidate_i, omega)
        # first, check if j has all of i targets
        for key, value in candidate_i.targets_dict.items():
            if candidate_i.genes_score_dict.get(key, 0) >= omega > candidate_j.genes_score_dict.get(key, 0):
                return False  # for testing. should be False

        cover_prob_j, num_of_covered_j, _ = cover_above_thr(candidate_j, omega)
        if num_of_covered_j > num_of_covered_i:
            if candidate_i.seq == 'GGCGAATGAGGAGGCTGAAG' and candidate_j.seq == 'GGCGAATGAGGAGGTTAAAG':
                print(num_of_covered_j, num_of_covered_i)
            return True
        if num_of_covered_j == num_of_covered_i:
            if cover_prob_j > cover_prob_i:
                return True
            if cover_prob_j == cover_prob_i:
                if candidate_min_mismatch_sum(candidate_j, omega, genes_list) > candidate_min_mismatch_sum(candidate_i, omega, genes_list):  # total num of mm over best site in relevant genes
                    return True

        return False  # for testing. should be False
    else:
        for key, value in candidate_i.targets_dict.items():
            if candidate_j.genes_score_dict.get(key, 0) < 0.005:  # epsilon
                return False
        return candidate_i.cut_expectation < candidate_j.cut_expectation


def remove_rep_subgroup(candidates_lst: List[Candidate], use_thr: int, omega: float) -> List[Candidate]:
    """

    :param candidates_lst: list of sgRNA candidates as Candidate objects
    :param use_thr:
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :return:
    """
    res = list()
    for i in range(len(candidates_lst)):
        to_append = 1
        for j in range(len(candidates_lst)):
            if i != j and contained(candidates_lst[i], candidates_lst[j], use_thr, omega):
                to_append = 0
        if to_append:
            res.append(candidates_lst[i])
    return res


def remove_repetitions_in_targets_sites(subgroup_list: List[SubgroupRes], use_thr: int, omega: float) -> List[SubgroupRes]:
    """

    :param subgroup_list:
    :param use_thr:
    :param omega:
    :return:
    """
    # if alg == 'gene_homology':
    res = list()
    for i in range(len(subgroup_list)):
        subgroup_candidates_lst = remove_rep_subgroup(subgroup_list[i].candidates_list, use_thr, omega)
        subgroup_candidates_lst[0].off_targets = True
        res.append(SubgroupRes(subgroup_list[i].genes_lst, subgroup_candidates_lst, subgroup_list[i].name))
    return res
    # else:
    #     res = remove_rep_subgroup(candidates_lst, use_thr, omega)
    #     for i in range(min(5, len(res))):
    #         res[i].off_targets = True
    #     return res
