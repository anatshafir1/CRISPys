from typing import List, Dict

import Distance_matrix_and_UPGMA
import Stage3
import Candidate
import Metric


def targets_tree_top_down(best_permutations, node, omega, targets_genes_dict, scoring_function, PS_number=12,
                          cfd_dict=None):
    if len(node.polymorphic_sites_set) < PS_number:  # was hardcoded to be 12, change it to be the PS_number argument Udi 24/02/22
        # make current_genes_sg_dict
        current_genes_sg_dict = dict()
        for target in node.node_targets:
            genes_leaf_from = targets_genes_dict[
                target]  # which gene is this target came from. usually it will be only 1 gene
            for gene_name in genes_leaf_from:
                if gene_name in current_genes_sg_dict:
                    if target not in current_genes_sg_dict[gene_name]:
                        current_genes_sg_dict[gene_name] += [target]
                else:
                    current_genes_sg_dict[gene_name] = [target]
        list_of_candidates = Stage3.find_candidates(current_genes_sg_dict, omega, scoring_function, node, cfd_dict)
        # current best perm is a tuple with the perm and metadata of this perm.
        if list_of_candidates:
            best_permutations += list_of_candidates
        return
    else:
        targets_tree_top_down(best_permutations, node.clades[0], omega, targets_genes_dict, scoring_function, PS_number,
                              cfd_dict)
        targets_tree_top_down(best_permutations, node.clades[1], omega, targets_genes_dict, scoring_function, PS_number,
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
        upgma_tree = return_targets_upgma(targets_list, targets_names, scoring_function, cfd_dict)
        fill_leaves_sets(upgma_tree, targets_genes_dict)
        fill_polymorphic_sites(upgma_tree.root)
        targets_tree_top_down(best_permutations, upgma_tree.root, omega, targets_genes_dict, scoring_function,
                              max_target_polymorphic_sites, cfd_dict)
    return best_permutations


def return_targets_upgma(seq_list, names_list, distance_function, cfd_dict=None):
    '''input:  a list of names and a list of sequences, calibrated
    output: an upgma instance.
    '''
    vectors_list = Metric.pos_in_metric_general(seq_list, distance_function)
    # create the distance matrix
    matrix = Distance_matrix_and_UPGMA.make_initial_matrix(vectors_list)
    m2 = Distance_matrix_and_UPGMA.make_distance_matrix(names_list, matrix)
    # apply UPGMA, return a target tree
    upgma1 = Distance_matrix_and_UPGMA.make_UPGMA(m2)
    return upgma1


def fill_leaves_sets(tree, sg_genes_dict):
    '''this version is not competable to genes tree.
    can be combine with fill_distance_from_leaves_function'''
    # fill the first line of nodes
    for leaf in tree.leaves:
        leaf.add_node_target(leaf.name)
        current_candidate = Candidate.Candidate(leaf.name)
        current_candidate.fill_default_fields(sg_genes_dict[leaf.name])
        leaf.set_candidates_DS()
        leaf.candidates[leaf.name] = current_candidate
        node = leaf
        while node.parent:
            for leaf in node.node_targets:
                if leaf not in node.parent.node_targets:
                    node.parent.add_node_target(leaf)
            node = node.parent


def two_seqs_differences_set(seq1, seq2):
    '''return a list of where the two sequences are different'''
    differences = set()  # key: place of disagreement. value: the suggestions of each side
    seq1 = seq1[:20]
    seq2 = seq2[:20]  # the pam is not considered when computing PS sites
    for i in range(1, len(seq2) - len(seq1)):
        differences.add(len(seq2) - i)
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences.add(i)
    return differences


def fill_polymorphic_site_node(node):
    # find differences between the representors
    polymorphic_site_set = set()
    if node.clades and len(node.clades) > 1:
        polymorphic_site_set = two_seqs_differences_set(node.clades[0].node_targets[0],
                                                        (node.clades[1].node_targets[0]))
        # update the rest of the sites
        for clade in node.clades:
            polymorphic_site_set.update(clade.polymorphic_sites_set)
    node.set_polymorphic_sites_set(polymorphic_site_set)


def fill_polymorphic_sites(node):
    '''
    :param node: tree's root
    :return:
    '''
    if not node:
        return
    if node.clades and len(node.clades) > 1:
        fill_polymorphic_sites(node.clades[0])
        fill_polymorphic_sites(node.clades[1])
    fill_polymorphic_site_node(node)
