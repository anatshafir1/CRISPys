
from typing import List, Dict
from Bio.Phylo import BaseTree
import Distance_matrix_and_UPGMA
import Stage3
import Candidate
from TreeConstruction_changed import CladeNew


def targets_tree_top_down(best_permutations: List, node: CladeNew, omega: float, targets_genes_dict: Dict,
                          scoring_function, max_target_polymorphic_sites: int = 12, cfd_dict: Dict = None):
    """

    :param best_permutations:
    :param node: starting at targets tree root
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param scoring_function: scoring function of the potential targets
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return:
    :rtype: list
    """
    if len(node.polymorphic_sites) < max_target_polymorphic_sites:  # was hardcoded to be 12, change it to be the max_target_polymorphic_sites argument Udi 24/02/22
        # make current_genes_sg_dict
        current_genes_targets_dict = dict()
        for target in node.node_targets:
            genes_leaf_from = targets_genes_dict[target]  # which gene is this target came from. usually it will be only 1 gene
            for gene_name in genes_leaf_from:
                if gene_name in current_genes_targets_dict:
                    if target not in current_genes_targets_dict[gene_name]:
                        current_genes_targets_dict[gene_name] += [target]
                else:
                    current_genes_targets_dict[gene_name] = [target]
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


def stage_two_main(targets_list: List, targets_names: List, targets_genes_dict: Dict, omega: float,
                   scoring_function, max_target_polymorphic_sites: int = 12, cfd_dict: Dict = None) -> List:
    """
    the main function of stage 2 of the algorithm. the function creates a UPGMA tree from the potential targets in
    'targets_list' and the given scoring function.

    :param targets_list: list of all the target sequences found in stage0
    :param targets_names: a deep copy of targets_list
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param scoring_function: scoring function of the potential targets
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return:
    :rtype: list
    """
    best_permutations = []
    if len(targets_list) == 1:
        print("only one sgRNA in the group")
        genes = targets_genes_dict[targets_list[0]]
        candidate = Candidate.Candidate(targets_list[0])
        candidate.fill_default_fields(genes)
        best_permutations.append(candidate)
    else:
        upgma_tree = Distance_matrix_and_UPGMA.return_targets_upgma(targets_list, targets_names, scoring_function, cfd_dict)
        fill_leaves_sets(upgma_tree)
        fill_polymorphic_sites(upgma_tree.root)
        targets_tree_top_down(best_permutations, upgma_tree.root, omega, targets_genes_dict, scoring_function,
                              max_target_polymorphic_sites, cfd_dict)
    return best_permutations


def fill_leaves_sets(tree: BaseTree):
    """
    For each node in the tree add list of the targets in its clade.

    :param tree: a UPGMA tree of potential targets
    """
    # fill the first line of nodes
    leaves = tree.get_terminals()
    for leaf in tree.leaves:
        leaf_clade = list(filter(lambda clade: (clade.name == leaf), leaves))[0]
        leaf_clade.add_node_target(leaf)
        path_to_leaf = tree.get_path(leaf_clade)
        path_to_leaf += [tree.root]
        for node in path_to_leaf[::-1]:
            if leaf_clade.node_targets[0] not in node.node_targets:
                node.add_node_target(leaf_clade.node_targets[0])


def two_seqs_differences_set(seq1: str, seq2: str) -> set:
    """return a set of where the two sequences are different"""
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
    :param node:
    """
    # find differences between the representors
    polymorphic_site_set = set()
    if node.clades and len(node.clades) > 1:
        polymorphic_site_set = two_seqs_differences_set(node.clades[0].node_targets[0],
                                                        (node.clades[1].node_targets[0]))
        # update the rest of the sites
        for clade in node.clades:
            polymorphic_site_set.update(clade.polymorphic_sites)
    node.fill_polymorphic_sites(polymorphic_site_set)


def fill_polymorphic_sites(node: CladeNew):
    """
    :param node: tree's root
    :return:
    """
    if not node:
        return
    if node.clades and len(node.clades) > 1:
        fill_polymorphic_sites(node.clades[0])
        fill_polymorphic_sites(node.clades[1])
    fill_polymorphic_site_node(node)
