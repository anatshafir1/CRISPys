"""Main"""
__author__ = 'Gal Hyams'

import timeit
import pickle
import argparse
import os
import sys
from typing import List, Dict
from pathlib import Path

from make_tree_display_CSV import tree_display, create_output_multiplex
from CasSites import fill_genes_targets_dict
from Stage1 import default_alg, gene_homology_alg
from Distance_matrix_and_UPGMA import MITScore, ccTop, gold_off_func, ucrispr, default_on_target
from Metric import cfd_funct
from SubgroupRes import SubgroupRes
from Candidate import Candidate
from CRISPR_Net.CrisprNetLoad import load_crispr_net
from MOFF.MoffLoad import load_moff
from DeepHF.LoadDeepHF import load_deephf
from crispys_chips import chips_main
# from singletons import singletons_main
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices'

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
    This function translates the chosen input for the scoring function and returns its pointer in the algorithm files
    and whether the function takes PAMs into its calculation.

    :param input_scoring_function: chosen scoring function by the user
    :return: scoring function pointer to use in the algorithm
    :rtype: function, boolean
    """
    pam_included = False
    if input_scoring_function == "CrisprMIT":
        return MITScore, pam_included
    elif input_scoring_function == "CCTop":
        return ccTop, pam_included
    elif input_scoring_function == "cfd_funct" or input_scoring_function == "cfd":
        return cfd_funct, pam_included
    elif input_scoring_function == "gold_off":
        pam_included = True
        return gold_off_func, pam_included
    elif input_scoring_function == "ucrispr" or input_scoring_function == "uCRISPR":
        pam_included = True
        return ucrispr, pam_included
    elif input_scoring_function == "crispr_net" or input_scoring_function == "crisprnet":
        pam_included = True
        return load_crispr_net(), pam_included
    elif input_scoring_function == "moff":
        pam_included = True
        return load_moff(), pam_included
    elif input_scoring_function == "DeepHF" or input_scoring_function == "deephf":
        pam_included = True
        return load_deephf(), pam_included
    elif input_scoring_function == "default":
        pam_included = True
        return default_on_target, pam_included
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
            if len(gene_seq) > 0 and gene_name != "":
                if gene_name not in genes_exons_dict:
                    genes_exons_dict[gene_name] = [gene_seq]
                else:
                    genes_exons_dict[gene_name] = genes_exons_dict[gene_name] + [gene_seq]
                gene_seq = ""
            if i != len(lines):
                gene_name = lines[i][1:].strip().upper()  # without the '>' and the '\n'
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

    NOTE: the target list was modified to be a list of tuples, each tuple is of (target, pos, strand),
    so the function here was also edited to take the first element of the tuple as the target (UDI 27/10/22)
    """
    targets_genes_dict = {}  # keys: target sequences. values: list of gene names in which the targets are found.
    for gene_name in genes_targets_dict:
        for target in genes_targets_dict[gene_name]:
            if target in targets_genes_dict:
                targets_genes_dict[target] += [gene_name]
            else:
                targets_genes_dict[target] = [gene_name]
    return targets_genes_dict


def add_coord_pam(res: List[SubgroupRes], genes_target_with_position: Dict) -> List[SubgroupRes]:
    """
    This function will add to targets of each candidate the pam the position in exome and the strand
    Args:
        res: list of SubgroupRes
        genes_target_with_position: a dictionary of gene: list of targets with position

    Returns:
        list of SubgroupRes
    """
    for subgroup in res:
        for candidate in subgroup.candidates_list:
            new_targets_dict = {}
            for gene in candidate.targets_dict:
                new_targets_list = []
                for target in candidate.targets_dict[gene]:
                    for target_list in genes_target_with_position[gene]:
                        if target_list[0] == target[0]:
                            new_target = target + [target_list[1], target_list[2], target_list[3]]
                            if new_target not in new_targets_list:
                                new_targets_list.append(new_target)
                new_targets_dict[gene] = new_targets_list
            candidate.targets_dict = new_targets_dict
    return res


def get_genes_of_interest_set(genes_of_interest_file, genes_exons_dict):
    """
    This function creates a genes of interest set from the genes_of_interest_file.
    :param genes_exons_dict: a dictionary genes -> exons
    :param genes_of_interest_file: a text file containing the genes of interest, seperated by rows.
    :return: a set of all genes of interest that appear in the family.
    """
    with open(genes_of_interest_file, 'r') as f:
        genes_from_file = f.readlines()
        return {gene.strip().upper() for gene in genes_from_file if gene.strip().upper() in genes_exons_dict}


def remove_sgrnas_without_gene_of_interest(res, genes_of_interest_set):
    """
    This function takes the results of crispys and removes any sgRNAs that don't target any genes from the set of genes
    of interest.
    :param res: The results of CRISPys
    :param genes_of_interest_set: A set of all genes of interest.
    :return: a new res containing the filtered subgroup results, without subgroups with no sgRNAs.
    """
    new_res = []
    for subgroup in res:
        new_candidates_list = [candidate for candidate in subgroup.candidates_list if
                               set(candidate.genes_score_dict).intersection(genes_of_interest_set)]
        new_subgroup = SubgroupRes(name=subgroup.name, genes_lst=subgroup.genes_lst, candidate_lst=new_candidates_list)
        if new_subgroup.candidates_list:
            new_res.append(new_subgroup)
    return new_res


def add_family_name(fasta_file_path: str, results: List[SubgroupRes]):
    """
    This function updates the family name attribute to each SubgroupRes object in CRISPys results.
    This function assumes that the name of the family is contained in the fasta, and that the name itself doesn't
    contain underscores.
    Args:
        fasta_file_path: path to the fasta file containing the exons of a gene family
        results: a list of CRISPys results as SubgroupRes objects
    Returns: None

    """
    family_name = Path(fasta_file_path).stem.split("_")[0]
    for subroup_res in results:
        subroup_res.family_name = family_name


def delete_file(file_path):
    """
    This function deletes a file if it exists.
    Args:
        file_path: path to file
    Returns: None
    """
    try:
        os.remove(file_path)
    except OSError:
        return


def CRISPys_main(fasta_file: str, output_path: str, output_name: str = "crispys_output",
                 genes_of_interest_file: str = 'None',
                 alg: str = "default",
                 where_in_gene: float = 1, use_thr: int = 1,
                 omega: float = 1, off_scoring_function: str = "cfd_funct", on_scoring_function: str = "default",
                 start_with_g: int = 0, internal_node_candidates: int = 10, max_target_polymorphic_sites: int = 12,
                 pams: int = 0, singletons_from_crispys: int = 0, slim_output: int = 0, set_cover: int = 0,
                 chips: int = 0, number_of_groups: int = 20, n_with_best_guide: int = 5, n_sgrnas: int = 2,
                 desired_genes_fraction_threshold: float = -1.0, singletons: int = 0,
                 singletons_on_target_function: str = "ucrispr", number_of_singletons: int = 50) -> List[SubgroupRes]:
    """
    Algorithm main function


    :param fasta_file: input text file output_path of gene names and their sequences (or their exons sequences) as lines
    :param output_path: the output_path to the directory in which the output files will be written
    :param output_name: the name that would be given to the crispys output.
    :param genes_of_interest_file: path to a txt file consisting of a "gene" column with genes of interest.
    :param alg: the type of the algorithm run - with gene homology or without
    :param where_in_gene: ignore targets sites downstream to the fractional part of the gene
    :param use_thr:
    :param omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param off_scoring_function: off target scoring function
    :param on_scoring_function: on target scoring function
    :param start_with_g: defines whether target sites are obligated to start with a G codon
    :param internal_node_candidates: number of sgRNAs designed for each homology subgroup
    :param max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
    :param pams: the pams by which potential sgRNA target sites will be searched
    :param singletons_from_crispys: optional choice to include singletons given by CRISPys
    :param slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
    :param set_cover: if 1, will output the minimal amount of guides that will capture all genes
    :param chips: if 1, output n candidates that will cover the most number of genes in the family, default is 0.
    :param number_of_groups: how many groups of 'best guide' to choose
    :param n_with_best_guide: for each group of 'best guide' how many multiplex to return
    :param n_sgrnas: the number of guides in each multiplex
    :param desired_genes_fraction_threshold: If a list of genes of interest was entered: the minimal fraction of genes
           of interest. CRISPys will ignore internal nodes with lower or equal fraction of genes of interest.
    :param singletons: select 1 to create singletons (sgRNAs candidates that target a single gene).
    :param number_of_singletons: the number of singletons that will be included for each gene.
    :param singletons_on_target_function: The on-target scoring function used for evaluating singletons.
    :return: List of sgRNA candidates as a SubgroupRes objects or Candidates object, depending on the algorithm run type

    """
    start = timeit.default_timer()
    # set the recursion limit to prevent recursion error
    sys.setrecursionlimit(10 ** 6)
    # choosing the scoring function:
    genes_exons_dict = fill_genes_exons_dict(fasta_file)  # gene name -> list of exons


    if len(genes_exons_dict) == 1 and not (singletons or singletons_from_crispys):
        print("family contains a single gene")
        return []

    off_scoring_function, pam_included = choose_scoring_function(off_scoring_function)
    on_scoring_function, pam_included = choose_scoring_function(on_scoring_function)

    genes_of_interest_set = {}
    if genes_of_interest_file != "None":
        genes_of_interest_set = get_genes_of_interest_set(genes_of_interest_file, genes_exons_dict)
        if not genes_of_interest_set:
            print(f"Job ended  for {output_name}. This family contains no genes from the list of interest")
            return []
    else:
        desired_genes_fraction_threshold = -1.0
    # find the potential sgRNA target sites for each gene:
    genes_targets_dict, genes_target_with_position = fill_genes_targets_dict(genes_exons_dict, pam_included,
                                                                             where_in_gene, start_with_g, pams)
    targets_genes_dict = inverse_genes_targets_dict(genes_targets_dict)  # target seq -> list of gene names
    genes_names_list = list(genes_targets_dict.keys())
    genes_list = get_genes_list(genes_exons_dict)  # a list of all the input genes in the algorithm
    res = []

    if alg == 'gene_homology':
        res = gene_homology_alg(genes_list, genes_names_list, genes_targets_dict, targets_genes_dict,
                                genes_of_interest_set, omega, output_path, off_scoring_function, on_scoring_function,
                                internal_node_candidates, max_target_polymorphic_sites, singletons_from_crispys,
                                desired_genes_fraction_threshold, slim_output)
    elif alg == 'default':  # alg == "default". automatically used on single-gene families.
        res = default_alg(targets_genes_dict, omega, off_scoring_function, on_scoring_function,
                          max_target_polymorphic_sites, singletons_from_crispys)
        # if singletons:
        #     singletons_main(genes_targets_dict, singletons_on_target_function, res, genes_of_interest_set,
        #                     number_of_singletons)
    if genes_of_interest_set:
        res = remove_sgrnas_without_gene_of_interest(res, genes_of_interest_set)
    if use_thr:
        sort_threshold(res, omega)
    else:
        sort_expectation(res)

    res = add_coord_pam(res, genes_target_with_position)
    add_family_name(fasta_file, res)
    pickle.dump(res, open(os.path.join(output_path, f"{output_name}.p"), "wb"))
    if alg == 'gene_homology':
        tree_display(output_path, res, genes_list, targets_genes_dict, omega, set_cover,
                     consider_homology=True, output_name=output_name)
    if alg == 'default':
        tree_display(output_path, res, genes_list, targets_genes_dict, omega, set_cover,
                     consider_homology=False, output_name=output_name)

    if chips:
        multiplex_dict = chips_main(res, number_of_groups, n_with_best_guide, n_sgrnas)
        # write results to csv
        create_output_multiplex(output_path, res, multiplex_dict, number_of_groups, n_with_best_guide, n_sgrnas,
                                output_name=output_name)
        pickle.dump(multiplex_dict, open(f"{output_path}/{output_name}_chips_dict.p", "wb"))

    if slim_output:
        delete_file(os.path.join(output_path, 'mafft_output_aligned_fasta.fa'))
        delete_file(os.path.join(output_path, 'genes_fasta_for_mafft.fa'))
        delete_file(os.path.join(output_path, 'infile'))
        delete_file(os.path.join(output_path, 'outfile.txt'))
    else:
        pickle.dump(genes_names_list, open(output_path + "/genes_names.p", "wb"))
        # add saving the gene_list in pickle in order to produce the results like in the server version - Udi 28/02/22
        pickle.dump(genes_list, open(output_path + '/genes_list.p', 'wb'))
        pickle.dump(targets_genes_dict, open(output_path + "/sg_genes.p", "wb"))
        # new output function taken from the crispys server code. Udi 13/04/2022

    stop = timeit.default_timer()
    if not slim_output:
        with open(os.path.join(output_path, "time.txt"), 'w') as time_file:
            time_file.write(str(stop - start))
    print(f"Job ended successfully for {output_name}. Output files created")
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
                            help='THe output_path to the directory in which the output files will be written.')
    parser_obj.add_argument('output_name', type=str, default='crispys_output',
                            help="The name that would be given to the crispys output.")
    parser_obj.add_argument('--genes_of_interest_file', type=str, metavar='<gene_list_path>', default='None',
                            help="path to a csv file consisting of a 'gene' column with genes of interest, "
                                 "and a 'family' column containing the family of the gene. the file must include a "
                                 "header row.")
    parser_obj.add_argument('--alg', type=str, default='default', help='Choose "gene_homology" to considering homology')
    parser_obj.add_argument('--where_in_gene', type=float, default=1,
                            help='input a number between 0 to 1 in order to ignore targets sites downstream to the '
                                 'fractional part of the gene')
    parser_obj.add_argument('--use_thr', '-t', type=int, default=1,
                            help='0 for using sgRNA to gain maximal gaining score among all of the input genes or 1 for'
                                 ' the maximal cleavage likelihood only among genes with score higher than the average.'
                                 ' Default: 1.')
    parser_obj.add_argument('--omega', '-v', type=float, default=0.43,
                            help='the value of the threshold. A number between 0 to 1 (included). Default: 0.43')
    parser_obj.add_argument('--off_scoring_function', '-s', type=str, default='cfd_funct',
                            help='the off scoring function of the targets. Optional scoring systems are: cfd_funct('
                                 'default), gold_off, CrisprMIT and CCtop. Additional scoring function may be added '
                                 'by the user or by request.')
    parser_obj.add_argument('--on_scoring_function', '-n', type=str, default='default',
                            help='the on scoring function of the targets. Optional scoring systems are: deephf. '
                                 'Additional scoring function may be added by the user or by request.')
    parser_obj.add_argument('--start_with_g', '-g', choices=[0, 1], type=int, default=0,
                            help='1 if the target sites are obligated to start with a G codon or 0 otherwise. '
                                 'Default: 0.')
    parser_obj.add_argument('--internal_node_candidates', '-i', type=int, default=10,
                            help='when choosing the considered homology option, this is the number of sgRNAs designed '
                                 'for each homology sub-group. Default: 10')
    parser_obj.add_argument('--max_target_polymorphic_sites', '-ps', type=int, default=12,
                            help='the maximal number of possible polymorphic sites in a target. Default: 12')
    parser_obj.add_argument('--pams', type=int, default=0,
                            help='0 to search NGG pam or 1 to search for NGG and NAG. Default: 0')

    parser_obj.add_argument('--singletons_from_crispys', choices=[0, 1], type=int, default=1,
                            help='1 to return results with singletons given by CRISPys (sgRNAs that target only 1 gene),'
                                 ' 0 to exclude'
                                 ' singletons given by CRISPys. Default: 1')

    parser_obj.add_argument('--slim_output', '-slim', choices=[0, 1], type=int, default=0,

                            help='optional choice to store only "res_in_lst" as the result of the algorithm run.'
                                 'Default: 0')
    parser_obj.add_argument('--set_cover', choices=[0, 1], type=int, default=0,
                            help='optional choice to output the minimal amount of guides that will capture all genes.'
                                 'Default: 0')
    parser_obj.add_argument('--chips', '-chips', choices=[0, 1], type=int, default=0,
                            help="optional: use 'chips' output that output the best n candidate that will target the "
                                 "highest number of genes in the family, if 0 output all"
                                 'Default: 0')
    parser_obj.add_argument('--number_of_groups', '-n_groups', type=int, default=0,
                            help="If using 'chips', for each group of 'best guide' how many multiplex to return"
                                 "by a 'best' guide)" 'Default: 20')
    parser_obj.add_argument('--n_with_best_guide', '-n_in_best', type=int, default=0,
                            help="If using 'chips', the number of multiplex results created per 'best' guide"
                                 'Default: 5')
    parser_obj.add_argument('--n_sgrnas', '-n_sgrnas', type=int, default=0,
                            help="If using 'chips', the number of guides in each multiplex" 'Default: 2')

    parser_obj.add_argument('--desired_genes_fraction_threshold', '-desired_genes_thr', type=float, default=-1.0,

                            help="If a list of genes of interest was entered: the minimal fraction of genes of "
                                 "interest. CRISPys will ignore internal nodes with lower or equal fraction of genes "
                                 "of interest. "
                                 'Default: -1.0')
    parser_obj.add_argument('--singletons', '-singletons', choices=[0, 1], type=int, default=0,
                            help="optional: select 1 to create singletons (sgRNAs candidates that target a single gene)"
                                 'Default: 0')
    parser_obj.add_argument('--singletons_on_target_function', type=str, default='ucrispr',
                            help='the on scoring functionthat ')
    parser_obj.add_argument('--number_of_singletons', '-num_singletons', type=int, default=50,
                            help="optional: the number of singleton candidates to include for each gene"
                                 'Default: 50')
    arguments = parser_obj.parse_known_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    CRISPys_main(fasta_file=args.fasta_file,
                 output_path=args.output_path,
                 output_name=args.output_name,
                 genes_of_interest_file=args.genes_of_interest_file,
                 alg=args.alg,
                 where_in_gene=args.where_in_gene,
                 use_thr=args.use_thr,
                 omega=args.omega,
                 off_scoring_function=args.off_scoring_function,
                 on_scoring_function=args.on_scoring_function,
                 start_with_g=args.start_with_g,
                 internal_node_candidates=args.internal_node_candidates,
                 max_target_polymorphic_sites=args.max_target_polymorphic_sites,
                 pams=args.pams,
                 chips=args.chips,
                 number_of_groups=args.number_of_groups,
                 n_with_best_guide=args.n_with_best_guide,
                 n_sgrnas=args.n_sgrnas,
                 singletons_from_crispys=args.singletons_from_crispys,
                 slim_output=args.slim_output,
                 set_cover=args.set_cover,
                 desired_genes_fraction_threshold=args.desired_genes_fraction_threshold,
                 singletons=args.singletons,
                 singletons_on_target_function=args.singletons_on_target_function,
                 number_of_singletons=args.number_of_singletons)
