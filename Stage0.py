"""Main"""
__author__ = 'Gal Hyams'

from typing import List, Dict
import CasSites
from CRISPR_Net.CrisprNetLoad import load_crispr_net
import Stage1
import Distance_matrix_and_UPGMA
import timeit
import pickle
import Metric
import argparse
import os
import sys
import make_tree_display_CSV
from MOFF.MoffLoad import load_moff
from SubgroupRes import SubgroupRes
from Candidate import Candidate
import glob

# get the output_path of this script file


PATH = os.path.dirname(os.path.realpath(__file__))


def sort_expectation(subgroups: List[SubgroupRes]):
    """
    Given a list of candidates the function sorts them by their cut expectation - the sum of cutting probabilities for
    all the genes the candidate cuts, and then by the number of mismatches between the candidate and its targets.

    :param subgroups: the result of the algorithm run as a list of candidates
    """
    for i in range(len(subgroups)):
        subgroups[i].candidates_list.sort(key=lambda item: (item.cut_expectation, item.total_num_of_mismatches()),
                                          reverse=True)


def sort_subgroup(candidates: List[Candidate], omega: float):
    """
    Accessory function for sorting candidates when sorting with threshold was chosen. For each candidate the function calculates
    the number of genes that the candidate sgRNA cuts with probability higher than omega, and the product of the cleaving
    probability across all the genes the candidate cleaves.

    :param candidates: current sgRNA candidate as a Candidate object
    :param omega: input omega threshold in the algorithm run
    """
    for candidate in candidates:
        num_of_genes_above_thr = 0
        cleave_all = 1
        for gene, score in candidate.genes_score_dict.items():
            if score >= omega:
                cleave_all *= score
                num_of_genes_above_thr += 1
        candidate.cleve_all_above_thr = cleave_all
        candidate.num_of_genes_above_thr = num_of_genes_above_thr
    candidates.sort(key=lambda item: (item.num_of_genes_above_thr, item.cleave_all_above_thr), reverse=True)


def sort_threshold(subgroups: List[SubgroupRes], omega: float):
    """
    Sort the candidates by number of genes with cut probability > omega and then by the probability to cleave all of
    these genes.

    :param subgroups: current sgRNA candidate as a Candidate object or as a SubgroupRes object
    :param omega: input omega threshold in the algorithm run
    """
    for i in range(len(subgroups)):
        sort_subgroup(subgroups[i].candidates_list, omega)


def choose_scoring_function(input_scoring_function: str):
    """
    This function translates the chosen input for the scoring function and returns its pointer in the algorithm files.

    :param input_scoring_function: chosen scoring function by the user
    :return: scoring function pointer to use in the algorithm
    :rtype: function
    """
    if input_scoring_function == "gold_off":
        return Distance_matrix_and_UPGMA.gold_off_func
    elif input_scoring_function == "CrisprMIT":
        return Distance_matrix_and_UPGMA.MITScore
    elif input_scoring_function == "CCTop":
        return Distance_matrix_and_UPGMA.ccTop
    elif input_scoring_function == "cfd_funct" or input_scoring_function == "cfd":
        return Metric.cfd_funct
    elif input_scoring_function == "DeepHF" or input_scoring_function == "deephf":
        return Distance_matrix_and_UPGMA.deephf
    elif input_scoring_function == "ucrispr" or input_scoring_function == "uCRISPR":
        return Distance_matrix_and_UPGMA.ucrispr
    elif input_scoring_function == "crispr_net" or input_scoring_function == "CRISPR_Net" or input_scoring_function == "crisprnet":
        return load_crispr_net()
    elif input_scoring_function == "moff":
        return load_moff()
    else:
        print("Did not specify valid scoring function")


def fill_genes_exons_dict(fasta_file: str) -> Dict[str, str]:
    """
    this function takes a fasta format file of genes and sequences and creates a dictionary where key are gene names
    and values are (one or more) exon sequences. If the input file has multiple exons per gene, the value will be
    a list of the exons in the gene.

    :param fasta_file: input text file output_path of gene names and their sequences (or their exons sequences) as lines
    :return: dictionary where key are genes names and values are sequences
    """
    file = open(fasta_file, 'r')
    gene_name = ""
    gene_seq = ""
    lines = file.readlines()
    i = 0
    genes_exons_dict = {}  # key: gene name. value: list of exons
    while i <= len(lines):
        if i == len(lines) or lines[i][0] == '>':
            if len(gene_seq) > 0 and gene_name != "":  # add the gene
                if gene_name not in genes_exons_dict:
                    genes_exons_dict[gene_name] = [gene_seq]
                else:
                    genes_exons_dict[gene_name] = genes_exons_dict[gene_name] + [gene_seq]
                gene_seq = ""
            if i != len(lines):
                gene_name = lines[i][1:].strip()  # without the '>' and the '\n'
        elif lines[i] != "\n":
            gene_seq += lines[i].strip()
        i += 1
    return genes_exons_dict


def get_genes_list(genes_exons_dict: Dict[str, str]) -> List[str]:
    """
    this function extracts the whole genes sequences from an input dictionary of {gene name: exons of gene} and returns
    a list of the genes.

    :param genes_exons_dict: dictionary of {gene name: exons of gene}
    :return: list of whole genes sequences
    """
    genes_list = []
    for gene_name in genes_exons_dict:
        genes_list.append("".join(genes_exons_dict[gene_name]))
    return genes_list


def inverse_genes_targets_dict(genes_targets_dict: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """
    inverts a dictionary of {gene names: list of target seqs in the gene} to a dictionary of
    {target sequence: list of gene names where the targets are found}.

    :param genes_targets_dict: dictionary of {keys: gene names. values: list of target seqs}
    :return: a dictionary of {target sequence: list of gene names where the targets are found}
    """
    targets_genes_dict = {}  # keys: target sequences. values: list of gene names in which the targets are found.
    for gene_name in genes_targets_dict:
        for target in genes_targets_dict[gene_name]:
            if target in targets_genes_dict:
                targets_genes_dict[target] += [gene_name]
            else:
                targets_genes_dict[target] = [gene_name]
    return targets_genes_dict


def CRISPys_main(fasta_file: str, output_path: str, alg: str = 'default', where_in_gene: float = 1, use_thr: int = 1,
                 omega: float = 1, scoring_function: str = "cfd_funct", min_length: int = 20, max_length: int = 20,
                 start_with_g: bool = False, internal_node_candidates: int = 10, max_target_polymorphic_sites: int = 12,
                 pams: int = 0, singletons: int = 0, slim_output: bool = False) -> List[SubgroupRes]:
    """
    Algorithm main function

    :param fasta_file: input text file output_path of gene names and their sequences (or their exons sequences) as lines
    :param output_path: the output_path to the directory in which the output files will be written
    :param alg: the type of the algorithm run - with gene homology or without
    :param where_in_gene: ignore targets sites downstream to the fractional part of the gene
    :param use_thr:
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param scoring_function: scoring function of the potential targets
    :param min_length: minimal length of the target site
    :param max_length: maximal length of the target site
    :param start_with_g: defines whether target sites are obligated to start with a G codon
    :param internal_node_candidates: number of sgRNAs designed for each homology subgroup
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param pams: the pams by which potential sgRNA target sites will be searched
    :param singletons: optional choice to include singletons (sgRNAs that target only 1 gene) in the results
    :param slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
    :return: List of sgRNA candidates as a SubgroupRes objects or Candidates object, depending on the algorithm run type
    """
    start = timeit.default_timer()
    # set the recursion limit to prevent recursion error
    sys.setrecursionlimit(10 ** 6)
    # choosing the scoring function:
    scoring_function_targets = choose_scoring_function(scoring_function)
    genes_exons_dict = fill_genes_exons_dict(fasta_file)  # gene name -> list of exons
    # find the potential sgRNA target sites for each gene:
    genes_targets_dict = CasSites.fill_genes_targets_dict(genes_exons_dict, scoring_function_targets, where_in_gene,
                                                          min_length, max_length,
                                                          start_with_g,
                                                          pams)  # gene names -> list of target seqs
    targets_genes_dict = inverse_genes_targets_dict(genes_targets_dict)  # target seq -> list of gene names
    genes_names_list = list(genes_targets_dict.keys())
    genes_list = get_genes_list(genes_exons_dict)  # a list of all the input genes in the algorithm
    res = []
    if alg == 'gene_homology':
        res = Stage1.gene_homology_alg(genes_list, genes_names_list, genes_targets_dict, targets_genes_dict, omega,
                                       output_path, scoring_function_targets, internal_node_candidates,
                                       max_target_polymorphic_sites, singletons, slim_output)
    elif alg == 'default':  # alg == "default" (look in article for better name)
        res = Stage1.default_alg(targets_genes_dict, omega, scoring_function_targets, max_target_polymorphic_sites,
                                 singletons)

    if use_thr:
        sort_threshold(res, omega)
    else:
        sort_expectation(res)

    if slim_output:
        for path in glob.glob(os.path.join(output_path,"*")):
            if path != fasta_file and not os.path.isdir(path):
                os.remove(path)
        pickle.dump(res, open(output_path + "/res_in_lst.p", "wb"))
    else:
        pickle.dump(res, open(output_path + "/res_in_lst.p", "wb"))
        pickle.dump(genes_names_list, open(output_path + "/genes_names.p", "wb"))
        # add saving the gene_list in pickle in order to produce the results like in the server version - Udi 28/02/22
        pickle.dump(genes_list, open(output_path + '/genes_list.p', 'wb'))
        pickle.dump(targets_genes_dict, open(output_path + "/sg_genes.p", "wb"))
        # new output function taken from the crispys server code. Udi 13/04/2022
        make_tree_display_CSV.tree_display(output_path, alg == 'gene_homology')

    stop = timeit.default_timer()
    if not slim_output:
        with open("time.txt", 'w') as time_file:
            time_file.write(str(stop - start))
    return res


def parse_arguments(parser_obj: argparse.ArgumentParser):
    """
    using a pars_obj object this function parses command line strings into python objects. the chosen parameters for the
    algorithm run are taken from this function.

    :param parser_obj: object for parsing command line strings into python objects
    :return: parsed parameters for the algorithm run
    :rtype: argparse.Namespace
    """
    parser_obj.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The output_path to the input fasta '
                                                                                 'file')
    parser_obj.add_argument('output_path', type=str, metavar='<output_path>',
                            help='THe output_path to the directory in which the output files will be written')
    parser_obj.add_argument('--alg', type=str, default='default', help='Choose "gene_homology" to considering homology')
    parser_obj.add_argument('--where_in_gene', type=float, default=1,
                            help='input a number between 0 to 1 in order to ignore targets sites downstream to the '
                                 'fractional part of the gene')
    parser_obj.add_argument('--use_thr', '-t', type=int, default=1,
                            help='0 for using sgRNA to gain maximal gaining score among all of the input genes or 1 for'
                                 ' the maximal cleavage likelihood only among genes with score higher than the average.'
                                 ' Default: 0.')
    parser_obj.add_argument('--omega', '-v', type=float, default=0.43,
                            help='the value of the threshold. A number between 0 to 1 (included). Default: 0.43')
    parser_obj.add_argument('--scoring_function', '-s', type=str, default='cfd_funct',
                            help='the scoring function of the targets. Optional scoring systems are: cfd_funct('
                                 'default), gold_off, CrisprMIT and CCtop. Additional scoring function may be added '
                                 'by the user or by request.')
    parser_obj.add_argument('--min_target_len', '-l', type=int, default=20, help='minimal length of the target site. '
                                                                                 'Default:20')
    parser_obj.add_argument('--max_target_len', '-m', type=bool, default=20, help='maximal length of the target site. '
                                                                                  'Default:20')
    parser_obj.add_argument('--start_with_g', '-g', type=bool, default=0,
                            help='1 if the target sites are obligated to start with a G codon or 0 otherwise. '
                                 'Default: 0.')
    parser_obj.add_argument('--internal_node_candidates', '-i', type=int, default=10,
                            help='when choosing the considered homology option, this is the number of sgRNAs designed '
                                 'for each homology sub-group. Default: 10')
    parser_obj.add_argument('--max_target_polymorphic_sites', '-ps', type=int, default=12,
                            help='the maximal number of possible polymorphic sites in a target. Default: 12')
    parser_obj.add_argument('--pams', type=int, default=0,
                            help='0 to search NGG pam or 1 to search for NGG and NAG. Default: 0')
    parser_obj.add_argument('--singletons', choices=[0, 1], type=int, default=0,
                            help='0 to return results with singletons (sgRNAs that target only 1 gene) 1 to exclude'
                                 ' singletons. Default: 0')
    parser_obj.add_argument('--slim_output', type=bool, default=False,
                            help='optional choice to store only "res_in_lst" as the result of the algorithm run.'
                                 'Default: False')

    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    CRISPys_main(fasta_file=args.fasta_file,
                 output_path=args.output_path,
                 alg=args.alg,
                 where_in_gene=args.where_in_gene,
                 use_thr=args.use_thr,
                 omega=args.omega,
                 scoring_function=args.scoring_function,
                 min_length=args.min_target_len,
                 max_length=args.max_target_len,
                 start_with_g=args.start_with_g,
                 internal_node_candidates=args.internal_node_candidates,
                 max_target_polymorphic_sites=args.max_target_polymorphic_sites,
                 pams=args.pams,
                 singletons=args.singletons,
                 slim_output=args.slim_output)
