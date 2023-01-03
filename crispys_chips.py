import pickle
import Candidate
import SubgroupRes
import copy
from typing import List, Dict
from Bio.Seq import reverse_complement
import random
random.seed(42)
"""
This module will takes CRISPys output (list of SubgroupRes) and will output group of n guides that will target the 
maximum number of the genes in a group of genes (internal node).
It will first select the 'best' guide, that is, the one that will target the most genes with maximum expectation score, 
than, it will select another guide that will complement the target to capture more genes that are not capture with the 
previous guide, this will continue until we get the number of guide we specify (usually the multiplex is of 2 guides).
Next, this peocess will repeat but without selecting the 'additional' guides, meaning, it will keep the first ('best') 
guide and remove from the list of all guides that completed the multiplex in previus step, and will find new ones to add to the 'bset' guide. 
the 'best' one, this argument is controled with the 'n_with_best_guide' parameter.
each multiplex (list of guides for one vector) is stored in a SubgroupRes object and a group of them are stored in a BestSgGroup object.
A representation of the BestSgGroup:
 
                         
                        --------------
                        | multiplx 1 | each multiplex is a SubgroupRes object with n candidates (specify in the n_sgrnas parameter) 
                        | multiplx 2 |
BestSgGroup object ->   |    .       |
                        |    .       |
                        | multiplx n |
                        --------------
each BestSgGroup has 'best candidate' sttribute that show the details of the guide that is common to all multiplex.

for each CRISPys output of gRNAs that target group of genes we can have multiple BestSgGroup objects each one with its 
own 'best' guide, so the next step is to remove the 'best guide' we selected from the CRISPys results and repeat the
first part again to create the next group of multiplex, this parameter is defined in 'number_of_groups'
so we get:

   BestSgGroup object 1    BestSgGroup object 2   .   .   .  BestSgGroup object n
   --------------          --------------                     --------------
   | multiplx 1 |          | multiplx 1 |                     | multiplx 1 |
   | multiplx 2 |          | multiplx 2 |                     | multiplx 2 |
   |    .       |          |    .       |        .    .   .   |    .       |
   |    .       |          |    .       |                     |    .       |
   | multiplx n |          | multiplx n |                     | multiplx n |
   --------------          --------------                     --------------
The described above is applied to each subgroup (represent internal node in the gene tree) of CRISPys results and
the output object is a dictionary that each key is the internal node name and the value is a dictionary that each key is
the sequence of the 'best' gRNA and the value is the BestSgGroup object of that 'best guide' 
(in each object are all multiplexes for that 'best guide') 

"""


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


def mark_duplicates(list_of_candidates):
    """
    This function check if candidate is duplicated sgRNAs with identical sequences.
    if the candidate is duplicate it uses the 'compare_duplicate_candidate' to compare the two candidate and
    mark the one the needs to be filtered as can.dup = True
    :param list_of_candidates: A list of CandidateWithOffTargets objects
    :return:
    """
    sequence_to_candidate_dict = {}
    for candidate in list_of_candidates:
        candidate.dup = False
        if candidate.seq not in sequence_to_candidate_dict:
            sequence_to_candidate_dict[candidate.seq] = candidate
            continue
        compare_duplicate_candidate(sequence_to_candidate_dict[candidate.seq], candidate)


def compare_duplicate_candidate(candidate_1, candidate_2):
    """
    This function mark the candidate that needs to be filtered between two duplicated sgRNA candidates.
    It first check wich one has the highest cut_expectation score and if equal, than how many genes each is covering
    The candidate to remove is marked with candidate.dup=True
    :param candidate_1: a dupicated Candidate object.
    :param candidate_2:a duplicated Candidate  object.
    :return:
    """
    assert candidate_2.seq == candidate_1.seq
    # Pick the candidate with the higher cut expectation.
    if candidate_1.cut_expectation > candidate_2.cut_expectation:
        candidate_2.dup = True
    elif candidate_1.cut_expectation == candidate_2.cut_expectation:
        # Pick the candidate from the largest subgroup. ??? (is that right TODO: check this assumption
        if len(candidate_1.genes_in_node) <= len(candidate_2.genes_in_node):
            candidate_1.dup = True


def filter_dup_from_subgroup(subgroup_lst):
    """
    The main function to remove duplicates candidates, it create a list of all candidates and feed it to 'mark_duplicates'
    The candidate that remain after filtering is the one that either have higher cut expectation score or if the score
    is the same will leave the one targeting more genes.
    Args:
        subgroup_lst: a list of subgroupRes (output of crispys)

    Returns:
        filter out candidates that have duplicate
    """
    # add the names of gene in node to each candidate
    add_node_gens_to_candidate(subgroup_lst)
    # go over each subgroup and delete duplicates
    list_of_all_candidates = [can for sub in subgroup_lst for can in sub.candidates_list]
    mark_duplicates(list_of_all_candidates)
    # go over the groups and remove duplicates
    for sub in subgroup_lst:
        for i, can in enumerate(sub.candidates_list):
            if can.dup:
                del sub.candidates_list[i]
        if not sub.candidates_list:
            subgroup_lst.remove(sub)


def add_node_gens_to_candidate(subgroup_res_list):
    """
    Add the 'genes_in_node' attribute to candidate in order to know how many genes in the candidate internal node
    Args:
        subgroup_res_list: list of subgroups

    Returns:

    """
    for sub in subgroup_res_list:
        for can in sub.candidates_list:
            can.genes_in_node = sub.genes_in_node


def remove_candidates_with_restriction_site(subgroup_lst: List, restriction_site: str = "None"):
    """
    This function returns a list of candidates which don't contain a given restriction site.

    Returns:

    """
    if restriction_site == "None":
        return subgroup_lst
    for subgroup in subgroup_lst:
        for i, cand in enumerate(subgroup.candidates_list):
            if (restriction_site in cand.seq) or (restriction_site in reverse_complement(cand.seq)):
                del(subgroup.candidates_list[i])

def add_singletons_to_subgroup(subgroup_list: list, number_of_singltones: int = 5) -> list:
    """
    This function takes the output of crispys that contain subgroupRes object with singletons (targeting one gene)
    and add the candidates of a singletons results to the object holding candidate for multiple genes (internal node in the gene tree)
    which one of the is the gene targeted by the singleton candidates.
    notice that the function add singltones by looking at the genes targeted in the internal node and not all the gene in the internal node,
    if you want to add singletons targeting all genes in the internal node change 'subgroup.genes_lst' to 'subgroup.genes_in_node'
    Args:
        subgroup_list: list of subgroupRes objects (output of crispys)
        number_of_singltones: the number of singletons from each gene to add to the subgroups of internal node

    Returns:
            A list of subgroupRes objects each one contain candidates for internal node in the gene tree (with singletons)
    """
    singletons_dict = {}
    subgroup_list_no_singles = []
    # create new list with only subgroups targetiong multiple genes and store the singletons in a dictionary of genes_name:candidates_list
    for subgroup in subgroup_list:
        if len(subgroup.genes_lst) == 1:
            singletons_dict[subgroup.genes_lst[0]] = subgroup.candidates_list
        else:
            subgroup_list_no_singles.append(subgroup)
    # sort singletons by on target score
    singletons_dict = {k: sorted(v, key=lambda x: x.on_target_score, reverse=True) for k, v in singletons_dict.items()}

    # go over the list of subgroups with no singletons and add the singleton candidate if it is not exists
    # the number of sinbgletones to add for each gene is define in the 'number_of_singltones' variable
    for subgroup in subgroup_list_no_singles:
        singletons2add = []
        # make a list of candidates sequences from the nun-singletons subgroup
        subgroup_cand_seqs = [can.seq for can in subgroup.candidates_list]
        # go over the singletons dictionary and check if the singleton target a gene of the subgroup
        for single_gene in singletons_dict.keys():
            if single_gene in subgroup.genes_in_node:
                # for 'number_of_singltones' times go over each candidate in the list and if it is not in the subgroup,
                # add it to the list of singltones
                i = 0
                for single_candidate in singletons_dict[single_gene][0: number_of_singltones]:
                        if i == number_of_singltones:
                            break
                        if single_candidate.seq not in subgroup_cand_seqs:
                            singletons2add.append(single_candidate)
                            i += 1
        # shfulle the singletons list and add it to subgroup candidate list
        random.shuffle(singletons2add)
        subgroup.candidates_list.extend(singletons2add)
    return subgroup_list_no_singles


def subgroup2dict(subgroup: SubgroupRes.SubgroupRes) -> Dict:
    """
    This function take a subgroup object and output a dictionary of sequence:candidate
    Args:
        subgroup: subgroups objects

    Returns: a dictionary of candidates
    """
    candidates_dict = {}
    for candidate in subgroup.candidates_list:
        # if the candidate is in the dictionary replace it with the same one that have higher cut expectation (if exist)
        if candidate.seq not in candidates_dict or candidate.cut_expectation > candidates_dict[candidate.seq].cut_expectation:
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
    This function take a dictionary of candidates and go over each one and calculates its score using a dictionary of
     coefficient for each gene score that determine the weight of each gene in the final score.
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


def check_overlap_positions(candidate: Candidate, can_pos_dict=None):
    """
    This function check if candidates overlap, it used in two stages: 1) check for overlap candidates inside multiplex
    and 2) check if 'best candidate' overlap with previous 'best' that been chosen
    Args:
        candidate: A candidate object
        can_pos_dict: positions dictionary of gene:positions (to compare with)
        check_bestgroup_overlap: a flag to know if comapring 'best' candidates
    Returns:

    """
    # if no position dictionary supplied create and return such dictionary
    if not can_pos_dict:
        can_pos_dict = {candidate.seq : dict()}
        for gene, target in candidate.targets_dict.items():
            # take only the positions in the best match of the gene (the first item in the list)
            pos_set = {(target[0][3], target[0][4])}
            can_pos_dict[candidate.seq][gene] = pos_set
        return can_pos_dict

    if can_pos_dict:
        for can in can_pos_dict.keys():
            overlaps = []
            for gene in candidate.targets_dict:
                # get candidate positions
                target = candidate.targets_dict[gene]
                can_pos_set = {(target[0][3], target[0][4])}
                try:
                    if can_pos_dict[can][gene] == can_pos_set:
                        overlaps.append(True)
                    else:
                        overlaps.append(False)
                # if the gene is not in the positional dictionary they are not fully overlap
                except KeyError:
                    overlaps.append(False)

            if all(overlaps):
                return True
        # add the new target to dictionary
        can_pos_dict[candidate.seq] = dict()
        for gene in candidate.targets_dict.keys():
            target = candidate.targets_dict[gene][0]
            pos_set = {(target[3], target[4])}
            can_pos_dict[candidate.seq][gene] = pos_set

        return False



def choose_candidates(subgroup: SubgroupRes.SubgroupRes, n_sgrnas: int = 2, best_candidate: Candidate = None,
                      pos_dict: Dict = None):
    """
    This function takes SubgroupRes object and returns n guides that will target
     as many genes as possible. this is the the function that produce the multiplex with n_sgrnas guides targeting the
     most genes in the subgroup
     when it is run for the first time on subgroup it is run without the a 'best_candidate' argument and it retruns
     the first multiplex with the 'Best sgRNA' as subgroup object (an object used in CRISPys to store results of
     internal node, used here for multiplex).
     When it is run subsequently with the 'Best_candidate' argument it finds the rest of the gRNA for the multiplex and
     retruns the multiplex (as subgroup object)
    Args:
        subgroup: a subgroup obhect conatining a list of candidates
        n_sgrnas: number of guide to output
        best_candidate: a 'best' candidate that was already chosen

    Returns: subgroup object containing a list of candidates, Candidate object containing the 'best candidate'

    """

    # get gene names for the family/node
    genes_names = subgroup.genes_in_node
    # make a dictionary of seq:candidate from crispys results
    candidates_dict = subgroup2dict(subgroup)
    # create initial coefficient dictionary of gene:coef (with coef = 1)
    genes_coef_dict = {gene: coef for gene, coef in zip(genes_names, [1 for i in genes_names])}
    # initate a dict that will store the gRNAs selected
    selected_candidates = {}

    # select best guide according to the initial 'genes_coef_dict'
    if not best_candidate:
        if len(candidates_dict) == 1:
            return None
        best_candidate = select_candidate(candidates_dict, genes_coef_dict)
        # calculate the guide position
        pos_dict = check_overlap_positions(best_candidate)
    # store the selected 'best' guide in a dictionary
    selected_candidates[best_candidate.seq] = best_candidate
    # re-calculate the coefficients dictionary according to the 'best' guide you found
    recalc_coef_dict(best_candidate, genes_coef_dict)

    # select the rest (other than the best) of the candidates for the amount specified in n_sgrnas
    i = 1
    while i < n_sgrnas:
        # check if no candidates left
        if not candidates_dict:
            return None
        # select
        candidate = select_candidate(candidates_dict, genes_coef_dict)
        # calculate the guide position
        skip_candidate = check_overlap_positions(candidate, pos_dict)

        if skip_candidate:
            del (candidates_dict[candidate.seq])
            continue
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
    # make a tuple of candidates sequence and use it as a name for the subgroup
    name = tuple(cand.seq for cand in cand_list)
    subgroup = SubgroupRes.SubgroupRes(list(set(genes_lst)), cand_list, name, subgroup.genes_in_node)
    return subgroup, best_candidate, pos_dict


def get_best_groups(subgroup: SubgroupRes.SubgroupRes, m_groups: int, n_sgrnas: int = 2):
    """
    This function takes crispys output as SubgroupRes object and finds m_groups of n_sgrnas.
    It is used for multiplexing while the n_sgrnas is the amount of guides in a single plasmid (the multiplex) and
    the m_groups is the number of groups of sg with the same 'best candidate' to return
    The algorithm takes the best sgRNA and match it with different guides to create m_group of multiplex each one with n_sgrnas guides.
    It returns a BestSgGroup object with the attribute 'subgroups' that store a list of subgroups the length of m_groups,
     each subgroup is a SubgroupRes object with n_sgrnas as the amount of candidates in its candidates_list
    Args:
        subgroup: SubgrouRes object
        m_groups: number of groups of guides (the number of BestSgGroups to create)
        n_sgrnas: number of guides in each group (for multiplexing)

    Returns:
            BestSgGroup object
    """
    # make a copy of subgroup
    subgroup_temp = copy.deepcopy(subgroup)
    # if no more candidates left stop the search
    if not subgroup_temp.candidates_list:
        return None
    # get the first group of sgRNAs, the best guide in the group and the positions of the candidate
    first_multiplex = choose_candidates(subgroup_temp, n_sgrnas)
    # check if a multiplex is found
    if not first_multiplex:
        return None
    # save the result to a BestSgGroup object, this object is design to hold all multiplex of the same 'best' sgRNA
    current_best = BestSgGroup()
    # store the subgroupres object with the candidates
    current_best.subgroups = [first_multiplex[0]]
    # add the 'best' and the postion
    current_best.best_candidate, pos_dict = first_multiplex[1], first_multiplex[2]
    # store a list of all candidates
    current_best.all_candidates = copy.copy(current_best.subgroups[0].candidates_list)
    while m_groups > 1:
        # remove the found sg from the subgroup (recreate it without them)
        subgroup_temp.candidates_list = [can for can in subgroup_temp.candidates_list if
                                         can not in current_best.all_candidates]
        # check that the are candidates left in the subgroup
        if not subgroup_temp.candidates_list:
            return current_best
        # choose the next group of sgRNA that will be joined with the 'best guide' found above
        try:
            # get the SubgroupRes, best_candidate and the updated positions list
            multiplx_group, best_candidate, pos_dict = choose_candidates(subgroup_temp, n_sgrnas,
                                                                         current_best.best_candidate, pos_dict)
            # add the SubgrouRes object containing the list of guides to the results
            current_best.subgroups.append(multiplx_group)
            current_best.all_candidates += [can for can in multiplx_group.candidates_list if
                                            can not in current_best.all_candidates]
            m_groups -= 1
        except TypeError:
            # print(f"No more candidate in group {current_best.best_candidate.seq} in node {subgroup_temp.name}")
            return current_best
    return current_best


def chips_main(subgroup_lst: List, number_of_groups, n_with_best_guide, n_sgrnas: int = 2,
               restriction_site: str = "None") -> Dict:
    """
    This function takes CRISPys results (a list of SubgroupRes) and 'number of groups' that specify the number of sgRNA
    groups for each subgroup. the 'n with best guide' specify the number of multiplex groups for each 'best guide'
    and the 'n sgrna' is the number of guides inside each group of best guide (the amount for multiplexing)
    The results are collected as follow:
     - each group of guides (the multiplexing) is in a subgroup object.
     - each 'best guide' group is in a list so we get a list of subgroups per best guide
     - each best guide list of subgroups is stored in a dictionary where the key is the best guide sequence.
     - the dictionaries are grouped to a list for each subgroup of the original output so final output is:
     list[dict1(best_guide1:list[SubgroupRes1, SubgroupRes2, SubgroupResn]), dict2(best_guide2), dictn()]

    Args:
        subgroup_lst: list of subgroups (CRISPys output)
        n_with_best_guide: the number of multiplex with the same 'best' guide
        n_sgrnas: number of guide in multiplex

    Returns: a dictionary

    """
    # Remove duplicates (skip this step fo now)
    # filter_dup_from_subgroup(subgroup_lst)
    # remove candidates with restriction site
    remove_candidates_with_restriction_site(subgroup_lst, restriction_site)
    # insert singleton subgroup to subgroups without singleton and create a list of subgroups without sinlgetons
    new_subgroups_lst = add_singletons_to_subgroup(subgroup_lst)
    # initiate result dictionary
    output_dict = {}
    # go over each group of results (for each internal node in gene tree)
    for subgroup in new_subgroups_lst:
        subgroup_dict = subgroup2dict(subgroup)
        # check if no candidate in node result
        if len(subgroup_dict) == 0:
            print(f"No CRISPys results for node {subgroup.name}")
            continue
        # initiate results dictionary
        bestsgroup_dict = {}
        # choose multiplex groups
        n = number_of_groups
        while n > 0:
            # get results for a group of guides with the same 'best' guide
            bestsgroup = get_best_groups(subgroup, n_with_best_guide, n_sgrnas)
            # if there are no more candidate it will return None
            if not bestsgroup:
                break
            # store the group of guides in a dictionary with best_sg_seq: BestSgGroup object
            bestsgroup_dict[bestsgroup.best_candidate.seq] = bestsgroup

            # check if 'best' we have got has the same position as other in the group and if so remove it from the
            # results dictionary and the group list
            # if it is the first best group, create the positional dictionary, otherwise compare with previous positions
            if len(bestsgroup_dict) == 1:
                pos_dict = check_overlap_positions(bestsgroup.best_candidate)
            else:
                skip_candidate = check_overlap_positions(bestsgroup.best_candidate, pos_dict)
                if skip_candidate:
                    subgroup.candidates_list.remove(subgroup_dict[bestsgroup.best_candidate.seq])
                    del bestsgroup_dict[bestsgroup.best_candidate.seq]
                    continue

            # remove the 'best candidate' from the list of candidates
            subgroup.candidates_list.remove(subgroup_dict[bestsgroup.best_candidate.seq])
            if not subgroup.candidates_list:
                break
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

# # 3 genes
# with open("/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out1/res_in_lst.p", 'rb') as f:
#     res = pickle.load(f)
#     print(f"Number of subgroups is: {len(res)}")
#

# cand_sub = choose_candidates(res, n_sgrnas=7)
# for i in cand_sub.candidates_list:
#     print(f"{i}\n{i.genes_score_dict.keys()}\n\n")

# groups = get_n_candidates(res, number_of_groups=5, n_with_best_guide=3, n_sgrnas=3)
# print(groups)
# print(len(groups))


# create_output_multiplex("/groups/itay_mayrose/udiland/crispys_test/test_files_git/for_debug/out", res, groups)
