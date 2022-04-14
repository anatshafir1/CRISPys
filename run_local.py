from Stage0 import *


CRISPys_main(fasta_file="/groups/itay_mayrose/udiland/crispys_test/compare_server_git/test/HOM04D000221_5.txt",
             path="/groups/itay_mayrose/udiland/crispys_test/compare_server_git/test/out_git",
             alg = 'E',
             where_in_gene = 0.8,
             use_thr = 1,
             Omega = 0.6,
             df_targets = Metric.cfd_funct,
             protodist_outfile = "outfile",
             min_length= 20,
             max_length = 20,
             start_with_G = False,
             internal_node_candidates = 200,
             PS_number = 5,
             PAMs=0)
