from Stage0 import *
import make_tree_display # from server
import set_cover_greedy # from server
from remove_rep import *
import argparse

def run_server(fasta_file, path, alg, where_in_gene, use_thr, Omega, df_targets,
                     protdist_outfile, min_length, max_length, start_with_G,
                     internal_node_candidates, PS_number, PAMs, off_targets='Not_selected'):
    """
    This function is a wrapper for crispys the is used to run the server version. i.e. to get output in html format and
    to be able to search for off targets using crista.
    all the arguments are passed to the main crispys function except for the 'off_target' that is handheld here

    Returns: in addition to the regular crispys output it will create html report of the whole data and cover set
    (all function written by Gal with minor changes if any)
    """

    # run the 'local' crispys
    CRISPys_main(fasta_file, path, alg, where_in_gene, use_thr, Omega, df_targets,
                     protdist_outfile, min_length, max_length, start_with_G,
                     internal_node_candidates, PS_number, PAMs)

    #read the results of crispys
    with open(path + "/res_in_lst.p", "rb") as cri_res:
        res = pickle.load(cri_res)

    with open(path + "/sg_genes.p", "rb") as sg_genes:
        sg_genes_dict = pickle.load(sg_genes)

    # create set cover html output
    if alg == 'A' and use_thr > 0:
        print( 'use thr: ', use_thr )
        greedy_cover = set_cover_greedy.find_set_cover(res, sg_genes_dict, Omega)
        for c in greedy_cover:
            c.off_targets = True
        pickle.dump( greedy_cover, open( path + '/greedy_cover.p', 'wb' ) )
        make_tree_display.tree_display( path, alg == 'E', 'greedy_set_cover', genomeAssembly=off_targets )

    # create html of all results
    make_tree_display.tree_display(path, alg == 'E', genomeAssembly=off_targets, use_thr=use_thr)

    # make a removed repetition results.
    removed_rep = remove_repetitions_in_targets_sites(res, alg, use_thr, Omega)
    pickle.dump(removed_rep, open(path + '/res_in_lst_removed_rep.p', 'wb'))

    make_tree_display.tree_display(path, alg == 'E', data='removed_repetitions', genomeAssembly=off_targets,
                                    use_thr=use_thr)


def parse_arguments(parser):
    parser.add_argument('fasta_file', type=str, metavar='<fasta_file>', help='The path to the input fasta file')
    parser.add_argument('path', type=str, metavar='<path>', help='The path to the directory in which the output files will be written')
    parser.add_argument('--alg', type=str, default='A', help='Choose E for considering homology')
    parser.add_argument('--where_in_gene', type=float, default=1, help='input a number between 0 to 1 in order to ignore targets sites downstream to the corresponding gene prefix')
    parser.add_argument('--t', type=int, default= 0, help='0 for using sgRNA to gain maximal gaining score among all of the input genes or 1 for the maximal cleavage likelihood only among genes with score higher than the average. Default: 0.')
    parser.add_argument('--v', type=float, default= 0.45, help='the value of the threshold. A number between 0 to 1 (included). Default: 0.43')
    parser.add_argument('--s', type=str, default='cfd_funct', help='the scoring function of the targets. Optional scoring systems are: cfd_funct (default), CrisprMIT and CCtop. Additinal scoring function may be added by the user or by request.')
    parser.add_argument('--p', type=str, default='outfile', help='protDist output file name. Default: "outfile"')
    parser.add_argument('--l', type=int, default= 20, help='minimal length of the target site. Default:20')
    parser.add_argument('--m', type=int, default= 20, help = 'maximal length of the target site, Default:20')
    parser.add_argument('--g', type=int, default= 0, help='1 if the target sites are obligated to start with a G codon or 0 otherwise. Default: 0.')
    parser.add_argument('--i', type=int, default= 10, help='when choosing the consider homology option, this is the number of sgRNAs designed for each homology sub-group. Default: 10')
    parser.add_argument('--ps', type=int, default= 12, help='the maximal number of possible polymorphic sites in a target. Default: 12')
    parser.add_argument('--PAMs', type=int, default=0, help='type of PAM to consider, 0=NGG, 1=NGG & NAG')
    parser.add_argument('--off_targets', type=str, default='Not_selected', help='Name of genome to search for off-target with CRISTA')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    run_server(fasta_file = args.fasta_file,
                path = args.path,
                alg = args.alg,
                where_in_gene = args.where_in_gene,
                use_thr = args.t,
                Omega=args.v,
                df_targets = args.s,
                protdist_outfile = args.p,
                min_length=args.l,
                max_length=args.m,
                start_with_G = args.g,
                internal_node_candidates=args.i,
                PS_number = args.ps,
                PAMs = args.PAMs,
                off_targets = args.off_targets)

# run_server(fasta_file="/groups/itay_mayrose/udiland/crispys_test/compare_server_git/test/HOM04D000221_5.txt",
#              path="/groups/itay_mayrose/udiland/crispys_test/compare_server_git/test/out_git",
#              alg = 'E',
#              where_in_gene = 0.8,
#              use_thr = 1,
#              Omega = 0.6,
#              df_targets = Metric.cfd_funct,
#              protdist_outfile = "outfile",
#              min_length= 20,
#              max_length = 20,
#              start_with_G = False,
#              internal_node_candidates = 200,
#              PS_number = 12,
#              PAMs=0,
#              off_targets='Not_selected') # example 'hg19' , for none use 'Not_selected'


