
from typing import List, Dict
import Distance_matrix_and_UPGMA
import Stage1
import Stage3
from Candidate import Candidate
from TreeConstruction_changed import CladeNew


def targets_tree_top_down(best_permutations: List[Candidate], node: CladeNew, omega: float, targets_genes_dict: Dict[str, List[str]],
                          scoring_function, max_target_polymorphic_sites: int = 12, cfd_dict: Dict = None):
    """
    Given an initial input of genomic targets UPGMA tree root, the function traverses the tree in a top-town (depth first) order.
    For each node with creates a dictionary of genes -> targets (leaves under the node) found in them, if the number of
    polymorphisms between the targets is less than 'max_target_polymorphic_sites'. Then the best candidate sgRNAs for
    the targets are found. The candidate sgRNAs are stored as a list of Candidate object.

    :param best_permutations: stored results as a list of Candidate objects from each function call
    :param node: starting at targets tree root
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param scoring_function: scoring function of the potential targets
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    """
    if node.is_terminal():
        return
    if len(node.polymorphic_sites) < max_target_polymorphic_sites:  # was hardcoded to be 12, change it to be the max_target_polymorphic_sites argument Udi 24/02/22
        # make current_genes_targets_dict
        current_genes_targets_dict = dict()
        for target in node.node_leaves:
            genes_leaf_from = targets_genes_dict[target]  # which gene is this target came from. usually it will be only 1 gene
            for gene_name in genes_leaf_from:
                if gene_name in current_genes_targets_dict:
                    if target not in current_genes_targets_dict[gene_name]:
                        current_genes_targets_dict[gene_name] += [target]
                else:
                    current_genes_targets_dict[gene_name] = [target]
        if len(current_genes_targets_dict) == 1:
            return
        # get candidates for the current node
        list_of_candidates = Stage3.return_candidates(current_genes_targets_dict, omega, scoring_function, node, cfd_dict)
        # current best perm is a tuple with the perm and metadata of this perm.
        if list_of_candidates:
            best_permutations += list_of_candidates
        return
    else:
        targets_tree_top_down(best_permutations, node.clades[0], omega, targets_genes_dict, scoring_function, max_target_polymorphic_sites,
                              cfd_dict)
        targets_tree_top_down(best_permutations, node.clades[1], omega, targets_genes_dict, scoring_function, max_target_polymorphic_sites,
                              cfd_dict)


def stage_two_main(targets_list: List[str], targets_names: List[str], targets_genes_dict: Dict[str, List[str]], omega: float,
                   scoring_function, max_target_polymorphic_sites: int = 12, cfd_dict: Dict = None) -> List[Candidate]:
    """
    the main function of stage 2 of the algorithm. the function creates a UPGMA tree from the potential targets in
    'targets_list' and the given scoring function. Then finds the best sgRNA for the potential targets in the tree.

    :param targets_list: list of all the target sequences found in stage0
    :param targets_names: a deep copy of targets_list
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param scoring_function: scoring function of the potential targets
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return: list of Candidate objects representing the best sgRNAs to cleave the targets in the UPGMA tree
    """
    best_permutations = []
    if len(targets_list) == 1:
        print("only one sgRNA in the group")
        genes = targets_genes_dict[targets_list[0]]
        candidate = Candidate(targets_list[0])
        candidate.fill_default_fields(genes)
        best_permutations.append(candidate)
    else:
        upgma_tree = Distance_matrix_and_UPGMA.return_targets_upgma(targets_list, targets_names, scoring_function, cfd_dict)
        Stage1.fill_nodes_leaves_list(upgma_tree)
        fill_polymorphic_sites(upgma_tree.root)
        targets_tree_top_down(best_permutations, upgma_tree.root, omega, targets_genes_dict, scoring_function,
                              max_target_polymorphic_sites, cfd_dict)
    return best_permutations


def two_seqs_differences_set(seq1: str, seq2: str) -> set:
    """
    Returns a set of polymorphic sites between the given sequences

    :param seq1: first sequence
    :param seq2: second sequence
    :return: set of polymorphic sites
    """
    differences = set()  # key: place of disagreement. value: the suggestions of each side
    seq1 = seq1[:20]
    seq2 = seq2[:20]  # the pam is not considered when computing PS sites
    for i in range(1, len(seq2) - len(seq1)):
        differences.add(len(seq2) - i)
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences.add(i)
    return differences


def fill_polymorphic_site_node(node: CladeNew):
    """
    Accessory function for 'fill_polymorphic_sites' function

    :param node: current node in the function recursion
    """
    # find differences between the representors
    polymorphic_site_set = set()
    if node.clades and len(node.clades) > 1:
        polymorphic_site_set = two_seqs_differences_set(node.clades[0].node_leaves[0],
                                                        (node.clades[1].node_leaves[0]))
        # update the rest of the sites
        for clade in node.clades:
            polymorphic_site_set.update(clade.polymorphic_sites)
    node.fill_polymorphic_sites(polymorphic_site_set)


def fill_polymorphic_sites(node: CladeNew):
    """
    For each node in the given targets' tree find the polymorphic sites between the targets under the node (leaves of
    the node) and fill the node's 'polymorphic_sites' with the set of polymorphic sites.

    :param node: tree's root
    """
    if not node:
        return
    if node.clades and len(node.clades) > 1:
        fill_polymorphic_sites(node.clades[0])
        fill_polymorphic_sites(node.clades[1])
    fill_polymorphic_site_node(node)
