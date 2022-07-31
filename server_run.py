import Stage0
from Stage0 import *
import make_tree_display
import set_cover_greedy
from remove_rep import *
import argparse


def run_server(fasta_file: str, output_path: str, alg: str, where_in_gene: float, use_thr: int, omega: float, scoring_function: str,
               min_length: int, max_length: int, start_with_g: bool, internal_node_candidates: int,
               max_target_polymorphic_sites: int, pams: int, off_targets='Not_selected'):
    """
    This function is a wrapper for crispys the is used to run the server version. i.e. to get output in html format and
    to be able to search for off targets using crista.
    all the arguments are passed to the main crispys function except for the 'off_target' that is handheld here

    Returns: in addition to the regular crispys output it will create html report of the whole data and cover set
    (all function written by Gal with minor changes if any)
    """

    # run the 'local' crispys
    CRISPys_main(fasta_file, output_path, alg, where_in_gene, use_thr, omega, scoring_function, min_length, max_length,
                 start_with_g, internal_node_candidates, max_target_polymorphic_sites, pams)

    # read the results of crispys
    with open(output_path + "/res_in_lst.p", "rb") as cri_res:
        res = pickle.load(cri_res)

    with open(output_path + "/sg_genes.p", "rb") as sg_genes:
        targets_genes_dict = pickle.load(sg_genes)

    # create set cover html output
    if alg == 'default' and use_thr == 1:
        print('use thr: ', use_thr)
        greedy_cover = set_cover_greedy.find_set_cover(res, targets_genes_dict, omega)
        for candidate in greedy_cover:
            candidate.off_targets = True
        pickle.dump(greedy_cover, open(output_path + '/greedy_cover.p', 'wb'))
        make_tree_display.tree_display(output_path, alg == 'gene_homology', 'greedy_set_cover', genomeAssembly=off_targets)

    # create html of all results
    make_tree_display.tree_display(output_path, alg == 'gene_homology', genomeAssembly=off_targets, use_thr=use_thr)

    # make a removed repetition results.
    removed_rep = remove_repetitions_in_targets_sites(res, use_thr, omega)
    pickle.dump(removed_rep, open(output_path + '/res_in_lst_removed_rep.p', 'wb'))

    make_tree_display.tree_display(output_path, alg == 'gene_homology', data='removed_repetitions', genomeAssembly=off_targets,
                                   use_thr=use_thr)


def parse_arguments(parser_obj: argparse.ArgumentParser):
    """
    using a pars_obj object this function parses command line strings into python objects. the chosen parameters for the
    algorithm run are taken from this function.

    :param parser_obj: object for parsing command line strings into python objects
    :return: parsed parameters for the algorithm run
    """
    parser_obj.add_argument('--off_targets', type=str, default='Not_selected',
                        help='Name of genome to search for off-target with CRISTA')
    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args_main = Stage0.parse_arguments(parser)
    args_add = parse_arguments(parser)
    run_server(fasta_file=args_main.fasta_file,
               output_path=args_main.path,
               alg=args_main.alg,
               where_in_gene=args_main.where_in_gene,
               use_thr=args_main.t,
               omega=args_main.v,
               scoring_function=args_main.s,
               min_length=args_main.l,
               max_length=args_main.m,
               start_with_g=args_main.g,
               internal_node_candidates=args_main.i,
               max_target_polymorphic_sites=args_main.ps,
               pams=args_main.PAMs,
               off_targets=args_add.off_targets)
