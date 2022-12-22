from typing import List, Dict, Set

from Candidate import Candidate
from SubgroupRes import SubgroupRes


def create_singletons_subgroup(gene, singletons_on_target_function, genes_targets_dict,
                               number_of_singletons):
    """
    This function creates a SubgroupRes object for a single gene containing singletons.
    Args:
        genes_targets_dict: A dictionary of gene -> targets
        singletons_on_target_function: An on-target function used for evaluating singletons
        gene: A string - the name of the gene
        number_of_singletons: The number of singleton candidates created per gene
    Returns: a SubgroupRes object containing singleton sgRNAs that target the gene

    """
    candidates_list = []
    singletons_scores_list = singletons_on_target_function(genes_targets_dict[gene])
    for singleton_sequence, score in zip(genes_targets_dict[gene], singletons_scores_list):
        singleton_candidate = Candidate(seq=singleton_sequence[:20], cut_expectation=1, genes_score_dict={gene: 1},
                                        mismatch_site_dict={gene: [[singleton_sequence[:20], {}]]})
        singleton_candidate.on_target_score = score
        candidates_list.append(singleton_candidate)
    candidates_list = sorted(candidates_list, key=lambda candidate: candidate.on_target_score, reverse=True)[
                      :number_of_singletons]
    singleton_subgroup = SubgroupRes(genes_lst=[gene], candidate_lst=candidates_list, name=gene, genes_in_node=gene)
    return singleton_subgroup


def singletons_main(genes_targets_dict: Dict, singletons_on_target_function, results: List[SubgroupRes],
                    number_of_singletons: int, genes_of_interest_set: Set):
    """
    This is the main function for adding single-gene sgRNAs - singletons_from_crispys.
    For each gene from the set of genes of interest, it appends a new SubgroupRes object containing
    singleton candidates, alongside their on-target score.
    Args:
        genes_targets_dict: A dictionary of gene -> targets
        singletons_on_target_function: An on-target function used for evaluating singletons
        results: A list of SubgroupRes objects given by CRISPys
        number_of_singletons: The number of singleton candidates created per gene
        genes_of_interest_set: A set of genes of interest. if a list was not defined, this function will find singletons_from_crispys
        for each gene in the family.
    Returns: add sub

    """
    if not genes_of_interest_set:
        genes_of_interest_set = {gene for gene in genes_targets_dict}
    for gene in genes_of_interest_set:
        singleton_subgroup = create_singletons_subgroup(gene, singletons_on_target_function, genes_targets_dict,
                                                        number_of_singletons)
        results.append(singleton_subgroup)
