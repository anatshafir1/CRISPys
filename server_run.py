
from Stage0 import *
import make_tree_display
import set_cover_greedy
from remove_rep import *
import argparse


def run_server(fasta_file: str, output_path: str, alg: str, where_in_gene: float, use_thr: int, Omega: float, scoring_function: str,
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
    CRISPys_main(fasta_file, output_path, alg, where_in_gene, use_thr, Omega, scoring_function, min_length, max_length,
                 start_with_g, internal_node_candidates, max_target_polymorphic_sites, pams)

    # read the results of crispys
    with open(output_path + "/res_in_lst.p", "rb") as cri_res:
        res = pickle.load(cri_res)

    with open(output_path + "/sg_genes.p", "rb") as sg_genes:
        sg_genes_dict = pickle.load(sg_genes)

    # create set cover html output
    if alg == 'default' and use_thr > 0:
        print('use thr: ', use_thr)
        greedy_cover = set_cover_greedy.find_set_cover(res, sg_genes_dict, Omega)
        for c in greedy_cover:
            c.off_targets = True
        pickle.dump(greedy_cover, open(output_path + '/greedy_cover.p', 'wb'))
        make_tree_display.tree_display(output_path, alg == 'gene_homology', 'greedy_set_cover', genomeAssembly=off_targets)

    # create html of all results
    make_tree_display.tree_display(output_path, alg == 'gene_homology', genomeAssembly=off_targets, use_thr=use_thr)

    # make a removed repetition results.
    removed_rep = remove_repetitions_in_targets_sites(res, alg, use_thr, Omega)
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
    parser_obj.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The output_path to the input fasta '
                                                                                 'file')
    parser_obj.add_argument('output_path', type=str, metavar='<output_path>',
                            help='THe output_path to the directory in which the output files will be written')
    parser_obj.add_argument('--alg', type=str, default='default', help='Choose "gene_homology" to considering homology')
    parser_obj.add_argument('--where_in_gene', type=float, default=1,
                            help='input a number between 0 to 1 in order to ignore targets sites downstream to the '
                                 'fractional part of the gene')
    parser_obj.add_argument('--use_thr', '-t', type=int, default=0,
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
    parser_obj.add_argument('--off_targets', type=str, default='Not_selected',
                        help='Name of genome to search for off-target with CRISTA')
    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    run_server(fasta_file=args.fasta_file,
               output_path=args.path,
               alg=args.alg,
               where_in_gene=args.where_in_gene,
               use_thr=args.t,
               Omega=args.v,
               scoring_function=args.s,
               min_length=args.l,
               max_length=args.m,
               start_with_g=args.g,
               internal_node_candidates=args.i,
               max_target_polymorphic_sites=args.ps,
               pams=args.PAMs,
               off_targets=args.off_targets)
