import pickle
import Candidate
import SubgroupRes
from typing import List, Dict

def subgroups2dict(subgroup_lst: List) -> Dict:
    """
    This function take a list of subgroups objects and ouput a dictionary of sequence:candidate
    Args:
        subgroup_lst: list of subgroups objects

    Returns: a dictionary of candidates

    """
    candidates_dict = {}
    for group in subgroup_lst:
        for candidate in group.candidates_list:
            if candidate.seq in candidates_dict:
                if len(candidate.genes_score_dict) > len(candidates_dict[candidate.seq].genes_score_dict):
                    candidates_dict[candidate.seq] = candidate
            else:
                candidates_dict[candidate.seq] = candidate
    return candidates_dict


def get_gene_names(subgroup_lst: List) -> set:
    """
    This function returns a set of all genes
    Args:
        subgroup_lst: list of subgroups objects

    Returns: set of genes names

    """
    gene_names_lst = []
    for group in subgroup_lst:
        gene_names_lst += group.genes_lst
    return set(gene_names_lst)


def get_relative_score(candidate: Candidate, coef_dict: Dict) -> float:
    """
    This function calculate the score of a candidate after recalibration of each gene score with a coefficient,
    the coefficient is coming from a dictionary of gene:coefficient
    Args:
        candidate: A candidate object
        coef_dict: a dictionary of gene:soefficient

    Returns: returns the score of a candidate considering the weight of each gene

    """
    score_total = 0
    for gene in coef_dict:
        try:
            score_total += (coef_dict[gene] * candidate.genes_score_dict[gene])
        # if the gene is not in the candidate dict go to the next gene
        except KeyError:
            continue
    return score_total

def select_candidate(candidates_dict: Dict, genes_coef_dict: Dict) -> Candidate:
    """
    This function take a dictionary of candidates and go over each one and calculates its score using a dictionary of coeffitionet for each gene score
    It return the candidate that got the highest score
    Args:
        candidates_dict: a dictionary of seq:candidate
        genes_coef_dict: a dictionary of gene:coefficient

    Returns:
        The candidate with the best score
    """
    # get the candidate with the best score
    high_score = 0
    best_candidate = None
    for candidate in candidates_dict.values():
        score = get_relative_score(candidate, genes_coef_dict)
        if score > high_score:
            best_candidate = candidate
            high_score = score
    return best_candidate


def recalc_coef_dict(candidate: Candidate, coef_dict: Dict, delta: int = 0.999):
    """
    This function update the coefficients in the gene:coef dictionary based on the scores of the genes in a candidate
    it reduce the existing coefficient value with the multiplication of it with the score of the gene from a given candidate (times some delta factor that prevent it to be zero)
    Args:
        candidate: a 'best candidate' that was selected in previous execution of 'select_candidate'
        coef_dict: the existing coefficients dictionary
        delta: a factor close to 1 that is used to prevent the coefficient to be zero

    Returns: it changes the existing coefficient dictionary

    """
    for gene in coef_dict:
        try:
            coef_dict[gene] = coef_dict[gene] * (1 - (delta * candidate.genes_score_dict[gene]))
        except KeyError:
            continue


def get_can_positions(candidate: Candidate) -> set:
    """
    This function return a set tuples with the position and strnad of each target in a candidate
    Args:
        candidate: A candidate object

    Returns: A set of tuples with the position and strand of each target, for example {(313, '+'), (192, '+')}

    """
    pos_set = set()
    for targets in candidate.targets_dict.values():
        for target in targets:
            pos_set.add((target[3], target[4]))
    return pos_set



def choose_candidates(subgroup_list: List, n_sgrnas: int =2) -> Dict:
    """
    This is the main function that takes CRISPys output (subgroup list) and retruns the n best guides for that results (family)
    Args:
        subgroup_list: an output of crispys for certain family
        n_sgrnas: number of guide to output

    Returns: a list of subgroupRes object (with only one item for competability with main function in Stage0)

    """
    # get gene names for the family
    genes_names = get_gene_names(subgroup_list)
    # make a dictionary of seq:candidate from crispys results
    candidates_dict = subgroups2dict(subgroup_list)
    # create initial coefficient dictionary of gene:coef (with coef = 1)
    genes_coef_dict = {gene: coef for gene, coef in zip(genes_names, [1 for i in genes_names])}
    selected_candidates = {}
    pos_lst = []
    i = 0
    # select candidates for the amount specified in n_sgrnas
    while i < n_sgrnas:
        # print(f"The gene coefficients dictionary: {genes_coef_dict}\n")
        # select best guide according to the 'genes_coef_dict'
        best_candidate = select_candidate(candidates_dict, genes_coef_dict)
        # calculate the guide position
        pos = get_can_positions(best_candidate)
        skip_candidate = False
        # check if a guide with the same position is already selected, if so, ignore the new one and find another
        for p in pos_lst:
            if pos == p:
                skip_candidate = True
        if skip_candidate:
            del(candidates_dict[best_candidate.seq])
            continue
        pos_lst.append(pos)
        # store the selected guide in a dictionary
        selected_candidates[best_candidate.seq] = best_candidate
        # re-calculate the coefficients dictionary according to the guide you found
        recalc_coef_dict(best_candidate, genes_coef_dict)
        i += 1
    # make output to a subgroup list
    cand_list = [can for can in selected_candidates.values()]
    genes = []
    for can in selected_candidates.values():
        genes += can.genes_score_dict.keys()
    genes = [set(genes)]
    genes_lst = genes
    name = "selected_candidates"
    subgroup = SubgroupRes.SubgroupRes(genes_lst, cand_list, name)
    return [subgroup]




# 8 genes
# with open("/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out/res_in_lst.p", 'rb') as f:
#     res = pickle.load(f)
#
# 2 genes
# with open("/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out2/res_in_lst.p", 'rb') as f:
#     res = pickle.load(f)

# 3 genes
# with open("/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out1/res_in_lst.p", 'rb') as f:
#     res = pickle.load(f)
#

# cand_sub = choose_candidates(res, n_sgrnas=7)
# for i in cand_sub.candidates_list:
#     print(f"{i}\n{i.genes_score_dict.keys()}\n\n")

