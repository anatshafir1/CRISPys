
import copy
from typing import List, Dict

from Candidate import Candidate


def prob_cover_genes_lst(candidate: Candidate, genes_lst: List) -> float:
    """

    :param candidate:
    :param genes_lst:
    :return:
    """
    cover_all = 1
    for gene in genes_lst:
        cover_all *= candidate.genes_score_dict[gene]
    return cover_all


def find_set_cover(best_permutations: List, targets_genes_dict: Dict, omega: float) -> List[Candidate]:
    """for now, might not work in a case when there is a gene that isn't covered by any of the permutations in the best_permutations. not finished. can make it more readable"""
    temp_best_perm = copy.copy(best_permutations[0].candidates_list)
    res = list()
    uncovered_genes = set()
    for target in targets_genes_dict:
        for gene in targets_genes_dict[target]:
            uncovered_genes.add(gene)
    while (len(uncovered_genes)) > 0 and len(temp_best_perm) > 0:
        best_current_perm, best_num_of_covered, best_prob_of_covered = None, 0, 0  # best_current_perm is the whole tuple
        i = 0
        while i < (len(temp_best_perm)):
            new_genes_covered = list()  # 0
            for gene, score in temp_best_perm[i].genes_score_dict.items():
                if gene in uncovered_genes and score >= omega:
                    new_genes_covered.append(gene)
            if len(new_genes_covered) == 0:
                i += 1
                continue
            elif len(new_genes_covered) >= best_num_of_covered:
                if len(new_genes_covered) > best_num_of_covered or prob_cover > best_prob_of_covered:  # cover more gene or cover the same amount with greater prob.
                    prob_cover = prob_cover_genes_lst(temp_best_perm[i], new_genes_covered)
                    best_num_of_covered, best_prob_of_covered = len(new_genes_covered), prob_cover
                    best_current_perm = temp_best_perm[i]
            i += 1
        if best_current_perm:
            res.append(best_current_perm)
            for gene, score in best_current_perm.genes_score_dict.items():
                if gene in uncovered_genes and score >= omega:  # there is a probability that this gene had already been covered by a previous sgRNA
                    uncovered_genes.remove(gene)
    return res
