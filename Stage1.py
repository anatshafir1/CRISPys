__author__ = 'GH'

import pickle
from os.path import dirname, abspath
from typing import List, Dict
import Distance_matrix_and_UPGMA
import SubgroupRes
from TreeConstruction_changed import CladeNew
from Bio.Phylo import BaseTree
from globals import *
import Stage2
import copy


def default_alg(input_targets_genes_dict: Dict, omega: float, scoring_function,
                max_target_polymorphic_sites: int = 12) -> List:
    """
    Called by the main function when choosing the default algorithm run. Given a dictionary of potential genomic targets
    this function returns a list of candidate sgRNAs (Candidate objects) that are the best suitable to target the input
    genes. The candidates list is then sorted in reverse order according to the cut expectation of each candidate,
    and the sorted list is returned.

    :param input_targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param scoring_function: scoring function of the potential targets
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :return: list of Candidate objects which are the best suitable to target the input genes
    :rtype: list
    """
    potential_targets_list = list(input_targets_genes_dict.keys())
    targets_names = copy.deepcopy(potential_targets_list)
    script_path = dirname(abspath(__file__))
    cfd_dict = pickle.load(open(script_path + "/cfd_dict.p", 'rb'))
    best_permutations = Stage2.stage_two_main(potential_targets_list, targets_names, input_targets_genes_dict, omega,
                                              scoring_function, max_target_polymorphic_sites, cfd_dict)
    best_permutations.sort(key=lambda item: item.cut_expectation, reverse=True)  # sort for the print
    return best_permutations


########################################################################################################################


def gene_homology_alg(genes_list: List, genes_names: List, genes_targets_dict: Dict, targets_genes_dict: Dict,
                      omega: float, output_path: str, scoring_function, internal_node_candidates: int,
                      max_target_polymorphic_sites: int = 12):
    """
    Called by the main function when choosing algorithm with gene homology taken in consideration. Creates a UPGMA tree
    from the input genes by their homology. Writes the tree to a newick format file and a preorder format file. Then
    finds the candidate sgRNAs (Candidate objects) that are the best suitable to target each of the input genes and
    stores them in subgroups for each gene. The function returns a list of the subgroups (SubgroupRes objects).

    :param genes_list: a list of the genes sequences input into the algorithm
    :param genes_names: a list of the names of the genes input into the algorithm
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param genes_targets_dict: a dictionary of gene -> list of potential targets found in the gene
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param output_path: the path to which the results will be stored
    :param scoring_function: scoring function of the potential targets
    :param internal_node_candidates: number of sgRNAs designed for each homology subgroup
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :return:
    :rtype: list
    """
    # make a tree and distance matrix of the genes
    genes_upgma_tree = Distance_matrix_and_UPGMA.return_protdist_upgma(genes_list, genes_names, output_path)
    # store the genes UPGMA tree in a newick format file
    write_newick_to_file(genes_upgma_tree.root, output_path)
    tree_to_file(genes_upgma_tree.root, output_path)
    fill_leaves_sets_genes_tree(genes_upgma_tree)  # tree leaves are genes
    # making the sgList for gene homology algorithm:
    list_of_subgroups = []
    script_path = dirname(abspath(__file__))  # add check if scf_funct
    cfd_dict = pickle.load(open(script_path + "/cfd_dict.p", 'rb'))
    genes_tree_top_down(list_of_subgroups, genes_upgma_tree.root, omega, genes_targets_dict, targets_genes_dict,
                        scoring_function, internal_node_candidates, max_target_polymorphic_sites, cfd_dict)
    return list_of_subgroups


def write_newick_to_file(tree_root: CladeNew, path: str):
    """
    This function stores a UPGMA tree clade in a newick format file.

    :param tree_root: a UPGMA tree clade to store
    :param path: to output path to which the newick format file will be stored
    """
    file = open(path + "/tree.newick", 'w')
    write_newick(tree_root, file, 0)
    file.write(';')
    file.close()


def write_newick(node: CladeNew, file, index: int):
    """
    An accessory function to store a UPGMA tree in a newick format file. Recursively write the tree nodes to a newick
    format file

    :param node: current node in the function call
    :param file: file object in which the tree is stored
    :param index: clade index
    """
    index += 1
    if node.is_terminal():
        file.write(node.name)
    else:
        file.write('(')
        for i in range(len(node.clades)):
            write_newick(node.clades[i], file, index)
            if i < len(node.clades) - 1:
                file.write(',')
        file.write(')n' + str(index))


def tree_to_file(tree_root: CladeNew, path: str):
    """
    This function stores a UPGMA tree clade in a preorder text format file.

    :param tree_root: a UPGMA tree clade to store
    :param path: to output path to which the preorder format file will be stored
    """
    lst = list()
    tree_preorder(tree_root, lst)
    print(lst, file=open(path + "/GenesTree.txt", 'w'))


def tree_preorder(node: CladeNew, lst: List):
    """
    An accessory function to store a UPGMA tree in a preorder format file. Recursively write the tree nodes to a
    preorder format file.

    :param node: current node of the function call
    :param lst: current list of preorder format
    """
    if not node:
        return
    lst.append(node.name)
    if node.clades:
        tree_preorder(node.clades[0], lst)
        if len(node.clades) > 1:
            tree_preorder(node.clades[1], lst)


def fill_leaves_sets_genes_tree(tree: BaseTree):
    """
    Given a UPGMA tree of genes the function fills the trees' nodes 'node_targets' with the gene names of the leaves
    under each node.
    :param tree: a UPGMA tree of the genes input to the algorithm run
    """
    # fill the first line of nodes
    for leaf in tree.leaves:
        leaf.add_node_target(leaf)
        node = leaf
        while node.parent:
            for leaf in node.node_targets:
                if leaf not in node.parent.node_targets:
                    node.parent.add_node_target(leaf)
            node = node.parent


# ############################################# Gene Homology top down ############################################### #


def genes_tree_top_down(res: List, node: CladeNew, omega: float, genes_targets_dict: Dict, targets_genes_dict: Dict,
                        scoring_function: Dict, internal_node_candidates: int = 10,
                        max_target_polymorphic_sites: int = 12, cfd_dict=None):
    """
    Given an initial input of genes UPGMA tree root the function traverses the tree in a top-town (depth first) order.
    For each node

    :param res: the result as a list of SubgroupRes objects
    :param node: the current node in the targets UPGMA tree for which the function is called
    :param omega: Threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param genes_targets_dict: a dictionary of gene -> list of potential targets found in the gene
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param scoring_function: scoring function of the potential targets
    :param internal_node_candidates: number of sgRNAs designed for each homology subgroup
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return: list of SubgroupsRes object which represent the best candidates for the algorithm run
    :rtype: List
    """
    # making the genes_targets dict for this subtree and the targets_genes_dict to send to the intermediate algorithm
    current_targets_genes_dict = dict()
    current_genes_targets_dict = dict()
    targets_list = list()
    targets_names = list()
    for leaf in node.node_targets:  # leaf here is a gene. taking only the relevant genes
        current_genes_targets_dict[leaf.name] = genes_targets_dict[leaf.name]
        # filling the targets to genes dict
        for target in current_genes_targets_dict[leaf.name]:
            if target in current_targets_genes_dict:
                current_targets_genes_dict[target] += [leaf.name]
            else:
                current_targets_genes_dict[target] = [leaf.name]
            if target not in targets_list:  # creating a list of target sequences and a list of target names
                targets_list.append(target)
                targets_names.append(target)
    if N_genes_in_node >= len(node.node_targets) > 1:  # I added the 'N_genes_in_node' from globals.py. Udi 16/03/22
        best_permutations = Stage2.stage_two_main(targets_list, targets_names, current_targets_genes_dict, omega,
                                                  scoring_function, max_target_polymorphic_sites, cfd_dict)
        if not best_permutations:
            return
        best_permutations.sort(key=lambda item: item.cut_expectation, reverse=True)
        current_best_perm = best_permutations[:internal_node_candidates]  # the best sg at the current set cover
        res.append(SubgroupRes.SubgroupRes(get_genes_list(best_permutations), current_best_perm, node.name))
    if not node.clades:
        return  # if the function recursion reached a final node (a leaf) - steps out of the current function call
    if node.clades[0]:
        genes_tree_top_down(res, node.clades[0], omega, genes_targets_dict, targets_genes_dict,
                            scoring_function, internal_node_candidates, max_target_polymorphic_sites,
                            cfd_dict)
    if node.clades[1]:
        genes_tree_top_down(res, node.clades[1], omega, genes_targets_dict, targets_genes_dict,
                            scoring_function, internal_node_candidates, max_target_polymorphic_sites,
                            cfd_dict)


def get_genes_list(candidates_lst: List):
    """
    This function takes a list of Candidate objects (representing sgRNA sequences) and returns a list of genes for which
    the candidates were found.

    :param candidates_lst: list of Candidate objects
    :return: list of gene names in which the candidates are found
    :rtype: list
    """
    res = set()
    for candidate in candidates_lst:
        for gene in candidate.genes_score_dict.keys():
            res.add(gene)
    return list(res)
