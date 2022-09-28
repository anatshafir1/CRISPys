"""Stage one"""
__author__ = 'GH'

from os.path import dirname, abspath
import copy
import pickle
from typing import List, Dict

from Bio.Phylo import BaseTree
from Distance_matrix_and_UPGMA import return_protdist_upgma
from Metric import cfd_funct
from SubgroupRes import SubgroupRes
from TreeConstruction_changed import CladeNew
from globals import *
from Stage2 import stage_two_main


def default_alg(input_targets_genes_dict: Dict[str, List[str]], omega: float, off_scoring_function, on_scoring_function,
                max_target_polymorphic_sites: int = 12, singletons: int = 0) -> List[SubgroupRes]:
    """
    Called by the main function when choosing the default algorithm run. Given a dictionary of potential genomic targets
    this function returns a list of candidate sgRNAs (Candidate objects) that are the best suitable to target the input
    genes. The candidates list is then sorted in reverse order according to the cut expectation of each candidate,
    and the sorted list is returned.

    :param input_targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param singletons: optional choice to include singletons (sgRNAs that target only 1 gene) in the results
    :return: list of Candidate objects which are the best suitable to target the input genes
    """
    potential_targets_list = list(input_targets_genes_dict.keys())
    targets_names = copy.deepcopy(potential_targets_list)
    cfd_dict = None
    if off_scoring_function == cfd_funct:
        script_path = dirname(abspath(__file__))
        cfd_dict = pickle.load(open(script_path + "/cfd_dict.p", 'rb'))
    best_permutations = stage_two_main(potential_targets_list, targets_names, input_targets_genes_dict, omega,
                                       off_scoring_function, on_scoring_function, max_target_polymorphic_sites,
                                       cfd_dict, singletons)
    best_permutations.sort(key=lambda item: item.cut_expectation, reverse=True)
    res = [SubgroupRes(get_genes_list(best_permutations), best_permutations, "total")]
    return res


########################################################################################################################


def gene_homology_alg(genes_list: List, genes_names: List, genes_targets_dict: Dict, targets_genes_dict: Dict,
                      omega: float, output_path: str, off_scoring_function, on_scoring_function, internal_node_candidates: int,
                      max_target_polymorphic_sites: int = 12, singletons: int = 0, slim_output: bool = False) -> List[SubgroupRes]:
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
    :param output_path: the output_path to which the results will be stored
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :param internal_node_candidates: number of sgRNAs designed for each homology subgroup
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param singletons: optional choice to include singletons (sgRNAs that target only 1 gene) in the results
    :param slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
    :return:
    """
    # make a tree and distance matrix of the genes
    genes_upgma_tree = return_protdist_upgma(genes_list, genes_names, output_path)
    # store the genes UPGMA tree in a newick format file
    if not slim_output:
        write_newick_to_file(genes_upgma_tree.root, output_path)
    fill_nodes_leaves_list(genes_upgma_tree)  # tree leaves are genes
    # making the sgList for gene homology algorithm:
    list_of_subgroups = []
    cfd_dict = None
    if off_scoring_function == cfd_funct:
        script_path = dirname(abspath(__file__))
        cfd_dict = pickle.load(open(script_path + "/cfd_dict.p", 'rb'))
    genes_tree_top_down(list_of_subgroups, genes_upgma_tree.root, omega, genes_targets_dict, targets_genes_dict,
                        off_scoring_function, on_scoring_function, internal_node_candidates, max_target_polymorphic_sites, cfd_dict, singletons)
    return list_of_subgroups


def write_newick_to_file(tree_root: CladeNew, path: str):
    """
    This function stores a UPGMA tree clade in a newick format file.

    :param tree_root: a UPGMA tree clade to store
    :param path: to output output_path to which the newick format file will be stored
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


def fill_nodes_leaves_list(tree: BaseTree):
    """
    Given a UPGMA tree of genes the function fills the trees' nodes 'node_leaves' with the gene names of the leaves
    under each node.

    :param tree: a UPGMA tree of the genes input to the algorithm run
    """
    # fill the first line of nodes
    leaves = tree.get_terminals()
    for leaf in tree.leaves:
        leaf_clade = list(filter(lambda clade: (clade.name == leaf), leaves))[0]
        leaf_clade.add_nodes_leaves(leaf)
        path_to_leaf = tree.get_path(leaf_clade)
        path_to_leaf += [tree.root]
        for node in path_to_leaf[::-1]:
            if leaf_clade.node_leaves[0] not in node.node_leaves:
                node.add_nodes_leaves(leaf_clade.node_leaves[0])


# ############################################# Gene Homology top down ############################################### #


def genes_tree_top_down(res: List, node: CladeNew, omega: float, genes_targets_dict: Dict[str, List[str]],
                        targets_genes_dict: Dict[str, List[str]], off_scoring_function, on_scoring_function, internal_node_candidates: int = 10,
                        max_target_polymorphic_sites: int = 12, cfd_dict=None, singletons: int = 0):
    """
    Given an initial input of genes UPGMA tree root the function traverses the tree in a top-town (depth first) order.
    For each node creates a dictionary of node's genes (leaves under the node) -> targets found in them, and then find
    the best candidate sgRNA for the targets. The candidates for each genes subgroup are stored as a SubgroupRes object.

    :param res: the result as a list of SubgroupRes objects
    :param node: the current node in the targets UPGMA tree for which the function is called
    :param omega: Threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param genes_targets_dict: a dictionary of gene -> list of potential targets found in the gene
    :param targets_genes_dict: a dictionary of target -> list of genes in which it was found
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :param internal_node_candidates: number of sgRNAs designed for each homology subgroup
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :param singletons: optional choice to include singletons (sgRNAs that target only 1 gene) in the results
    """
    # making the genes_targets dict for this subtree and the targets_genes_dict to send to the intermediate algorithm
    current_targets_genes_dict = dict()
    current_genes_targets_dict = dict()
    targets_list = list()
    targets_names = list()
    if N_genes_in_node >= len(
            node.node_leaves) > singletons:  # I added the 'N_genes_in_node' from globals.py. Udi 16/03/22. check if the node is one gene
        for leaf in node.node_leaves:  # leaf here is a gene. taking only the relevant genes
            current_genes_targets_dict[leaf] = genes_targets_dict[leaf]
            # filling the targets to genes dict
            for target in current_genes_targets_dict[leaf]:
                if target in current_targets_genes_dict:
                    current_targets_genes_dict[target] += [leaf]
                else:
                    current_targets_genes_dict[target] = [leaf]
                if target not in targets_list:  # creating a list of target sequences and a list of target names
                    targets_list.append(target)
                    targets_names.append(target)

        best_permutations = stage_two_main(targets_list, targets_names, current_targets_genes_dict, omega,
                                           off_scoring_function, on_scoring_function, max_target_polymorphic_sites, cfd_dict, singletons)
        if not best_permutations:
            return
        best_permutations.sort(key=lambda item: item.cut_expectation, reverse=True)
        current_best_perm = best_permutations[:internal_node_candidates]  # the best sg at the current set cover
        res.append(SubgroupRes(get_genes_list(best_permutations), current_best_perm, node.name))
    if not node.clades:
        return  # if the function recursion reached a final node (a leaf) - steps out of the current function call
    if node.clades[0]:
        genes_tree_top_down(res, node.clades[0], omega, genes_targets_dict, targets_genes_dict,
                            off_scoring_function, on_scoring_function, internal_node_candidates, max_target_polymorphic_sites,
                            cfd_dict, singletons)
    if node.clades[1]:
        genes_tree_top_down(res, node.clades[1], omega, genes_targets_dict, targets_genes_dict,
                            off_scoring_function, on_scoring_function, internal_node_candidates, max_target_polymorphic_sites,
                            cfd_dict, singletons)


def get_genes_list(candidates_lst: List) -> List[str]:
    """
    This function takes a list of Candidate objects (representing sgRNA sequences) and returns a list of genes for which
    the candidates were found.

    :param candidates_lst: list of Candidate objects
    :return: list of gene names in which the candidates are found
    """
    res = set()
    for candidate in candidates_lst:
        for gene in candidate.genes_score_dict.keys():
            res.add(gene)
    return list(res)
