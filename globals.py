"""
this script is for storing global variables
"""

# this variable will determine the maximum number of genes that will be considered for finding
# its used in Stage1_h.py with function:
N_genes_in_node = 10


# this variable will determine the number of cores used for the multi_processing used by gold_off.predict
# its used in Distance_matrix_and_UPGMA.gold_off_func:
n_cores_for_gold_off = 10