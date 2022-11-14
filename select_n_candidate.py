import pickle
import Candidate
import SubgroupRes
import copy
from typing import List, Dict



class BestSgGroup:
    """
    This class is made to store a collection of subgroups object, such that each contains group of candidates that where
    selected fo multiplexing.
    The multiplexing group is designed so that each subgroup contain one candidate that is the same
     (referred to as 'best candidate') and another (one or more) that are different in each group this why the class is
     called BestSgGroup and the 'best_candidate' attribute will store the sequence of the best candidate.
    """
    def __init__(self, best_candidate: Candidate.Candidate = None, subgroups: List = list(), all_candidates: List = list()):
        self.best_candidate = best_candidate
        self.subgroups = subgroups
        self.all_candidates = all_candidates

    def __str__(self):
        return f"{self.best_candidate} , {self.subgroups}"



def subgroup2dict(subgroup: SubgroupRes.SubgroupRes) -> Dict:
    """
    This function take a subgroup object and ouput a dictionary of sequence:candidate
    Args:
        subgroup: subgroups objects

    Returns: a dictionary of candidates
    """
    candidates_dict = {}
    for candidate in subgroup.candidates_list:
        if candidate.seq not in candidates_dict:
            candidates_dict[candidate.seq] = candidate
    return candidates_dict


def get_gene_names(subgroup_lst: List) -> set:
    """
    **Used to get gene names from list of subgroups DEPRACTED**
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


def recalc_coef_dict(candidate: Candidate, coef_dict: Dict, delta: int = 0.9):
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


def choose_candidates(subgroup: SubgroupRes.SubgroupRes, n_sgrnas: int = 2, best_candidate: Candidate = None,
                      pos_lst:List=None):
    """
    This is the main function that takes CRISPys output (subgroup) and returns the n best guides that will target
     as many genes as possible.
    Args:
        subgroup: a subgroup obhect conatining a list of candidates
        n_sgrnas: number of guide to output

    Returns: subgroup object containig a list of candidates

    """

    # get gene names for the family
    genes_names = subgroup.genes_lst
    # make a dictionary of seq:candidate from crispys results
    candidates_dict = subgroup2dict(subgroup)
    # create initial coefficient dictionary of gene:coef (with coef = 1)
    genes_coef_dict = {gene: coef for gene, coef in zip(genes_names, [1 for i in genes_names])}
    selected_candidates = {}

    # select best guide according to the initial 'genes_coef_dict'
    if not best_candidate:
        best_candidate = select_candidate(candidates_dict, genes_coef_dict)
        # calculate the guide position
        pos = get_can_positions(best_candidate)
        pos_lst = [pos]
    # store the selected guide in a dictionary
    selected_candidates[best_candidate.seq] = best_candidate
    # re-calculate the coefficients dictionary according to the guide you found
    recalc_coef_dict(best_candidate, genes_coef_dict)

    # select the rest (other than the best) of the candidates for the amount specified in n_sgrnas
    i = 1
    while i < n_sgrnas:
        # check if no candidates left
        if not candidates_dict:
            print(f"No more candidates to choose from in {subgroup.name}")
            return None
        candidate = select_candidate(candidates_dict, genes_coef_dict)
        # calculate the guide position
        pos = get_can_positions(candidate)
        skip_candidate = False
        # check if a guide with the same position is already selected, if so, ignore the new one and find another
        for p in pos_lst:
            if pos == p:
                skip_candidate = True
        if skip_candidate:
            del (candidates_dict[candidate.seq])
            continue
        pos_lst.append(pos)
        # store the selected guide in a dictionary
        selected_candidates[candidate.seq] = candidate
        # re-calculate the coefficients dictionary according to the guide you found
        recalc_coef_dict(candidate, genes_coef_dict)
        i += 1
    # make output to a subgroup list
    cand_list = [can for can in selected_candidates.values()]
    genes = []
    for can in selected_candidates.values():
        genes += can.genes_score_dict.keys()
    genes_lst = genes
    name = "selected_candidates"
    subgroup = SubgroupRes.SubgroupRes(list(set(genes_lst)), cand_list, name)
    return subgroup, best_candidate, pos_lst


def get_candidats_groups(subgroup: SubgroupRes.SubgroupRes, m_groups: int, n_sgrnas: int):
    """
    This function takes a group of genes and candidates as SubgroupRes object and finds m_groups of n_sgrnas each
    It is used for multiplexing while the n_sgrnas is the amount of guides in a single vector and the m_groups is the
    number of groups to return
    The algorithm takes the best sgRNA and match him with different guides in each m_group
    It returns a list of subgroups the length of m_groups, the amount of candidates in each subgroup is n_sgrnas
    Args:
        subgroup: SubgrouRes object
        m_groups: number of groups of guides (the number of BestSgGroups to create)
        n_sgrnas: number of guides in each group (for multiplexing)

    Returns:

    """
    # make a copy of subgroup
    subgroup_temp = copy.deepcopy(subgroup)
    # get the first group of sgRNAs and the best guide in the group
    multiplx_candidates = choose_candidates(subgroup_temp, n_sgrnas)
    current_best = BestSgGroup()
    # store the 'best' candidate
    current_best.subgroups = [multiplx_candidates[0]]
    current_best.best_candidate, pos_lst = multiplx_candidates[1], multiplx_candidates[2]
    # add the candidate to the all list
    current_best.all_candidates = copy.copy(current_best.subgroups[0].candidates_list)
    while m_groups > 1:
        # remove the found sg from the subgroup
        subgroup_temp.candidates_list = [can for can in subgroup_temp.candidates_list if can not in current_best.all_candidates]
        # choose the next group of sgRNA that will be joined with the 'best guide' found above
        try:
            pair_group, best_candidate, pos_lst = choose_candidates(subgroup_temp, n_sgrnas, current_best.best_candidate, pos_lst)
            # add the SubgrouRes object containing the list of guides to the results
            current_best.subgroups.append(pair_group)
            current_best.all_candidates += [can for can in pair_group.candidates_list if can not in current_best.all_candidates]
            m_groups -= 1
        except TypeError:
            print(f"No more candidate in group {subgroup.name}")
            return
    return current_best

def get_n_candidates(subgroup_lst: List, number_of_groups, n_with_best_guide, n_sgrnas: int = 2) -> Dict:
    """
    This function takes CRISPys results (a list of SubgroupRes) and 'number of groups' that specify the number of sgRNA
    groups for each subgroup. the 'n with best guide' specify the number of groups for each 'best guide'
    and the 'n sgrna' is the number of guides inside each group of best guide (the amount for multiplexing)
    The results are collected as follow:
     - each group of guides (the multiplexing) is in a subgroup object.
     - each 'best guide' group is in a list so we get a list of subgroups per best guide
     - each best guide list of subgroups is stored in a dictionary where the key is the best guide sequence.
     - the dictionaries are grouped to a list for each subgroup of the original output so final output is:
     list[dict1(best_guide1:list[SubgroupRes1, SubgroupRes2, SubgroupResn]), dict2(best_guide2), dictn()]

    Args:
        subgroup_lst:
        n_with_best_guide:
        n_sgrnas:

    Returns: a dictionary

    """
    output_dict = {}
    for subgroup in subgroup_lst:
        subgroup_dict = subgroup2dict(subgroup)
        bestsgroup_dict = {}
        pos_lst = []
        n = number_of_groups
        while n > 0:
            # get results for one group of genes
            bestsgroup = get_candidats_groups(subgroup, n_with_best_guide, n_sgrnas)
            # if threr are no more candidate it will return None
            if not bestsgroup:
                break
            # store the group of guides in a dictionary with best_sg_seq: BestSgGroup object
            bestsgroup_dict[bestsgroup.best_candidate.seq] = bestsgroup

            # check if 'best' we have got has the same position as other in the group and if so remove it from the
            # results dictionary and the group list
            best_pos = get_can_positions(bestsgroup.best_candidate)
            skip_candidate = False
            if pos_lst:
                for pos in pos_lst:
                    if pos == best_pos:
                        skip_candidate = True

            if skip_candidate:
                subgroup.candidates_list.remove(subgroup_dict[bestsgroup.best_candidate.seq])
                del bestsgroup_dict[bestsgroup.best_candidate.seq]
                continue
            else:
                pos_lst.append(best_pos)

            # remove the 'best candidate' from the list of candidates
            subgroup.candidates_list.remove(subgroup_dict[bestsgroup.best_candidate.seq])

            n -= 1
        output_dict[subgroup.name] = bestsgroup_dict
    return output_dict


# 8 genes
# with open("/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out/res_in_lst.p", 'rb') as f:
#     res = pickle.load(f)
    # for sub in res1:
    #     if sub.name == "Inner4":
    #         res = [sub]

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

# groups = get_n_candidates(res, number_of_groups=20, n_with_best_guide=5, n_sgrnas=2)
# print(groups)
# print(len(groups))
