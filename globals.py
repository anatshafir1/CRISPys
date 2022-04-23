"""
this script is for storing global variables
"""

# this variable will determine the maximum number of genes that will be considered for finding
# it's used in Stage1_h.py with function:
N_genes_in_node = 10


# this variable will determine the number of cores used for the multi_processing used by gold_off.predict
# important: set n_process to 1 when debugging, otherwise the code can get stuck
# it's used in Distance_matrix_and_UPGMA.gold_off_func:
n_cores_for_gold_off = 10


# this variable stores the name of the model used by gold_off_func the file needs to be in the code directory
# it's used in Distance_matrix_and_UPGMA.gold_off_func:
xgb_model_name = "regression_gs_log_max.xgb"