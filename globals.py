import os
"""
this script is for storing global variables
"""
# get th epath of the scripts directory
PATH = os.path.dirname(os.path.realpath(__file__))

# This variable will determine the maximum number of genes that will be considered for finding
# it's used in Stage1_h.py with function:
N_genes_in_node = 10

# This variable defines the seed for CRISPys' random functions
seed = 1234

# This variable stores the maximum vector size that is allowed for the distance transformation.
# This variable used in Distance_matrix_and_UPGMA.pos_in_metric_general
vector_size_cutoff = 1000

# this variable will determine the number of cores used for the multi_processing used by gold_off.predict
# important: currently, gold off only supports 10 n_cores.
n_cores_for_gold_off = 10

# This variable stores the name of the model used by gold_off_func the file needs to be in the code directory.
# It's used in Distance_matrix_and_UPGMA.gold_off_func:
xgb_model_name = "regression_union_log_max.xgb"

# in case you need to activate conda when using remote debug
CONDA = "source /groups/itay_mayrose/udiland/miniconda3/etc/profile.d/conda.sh; conda activate crispys"

#command that is used to connect to server and run crista
ssh_conect = 'ssh bioseq@powerlogin "module load python/python-anaconda3.6.5 && '

def set_crisprnet_model(model):
    """
    A function to make a global variable 0f the model of crispr-net
    Args:
        model: crispr_net trained model

    Returns: global variable
    """
    global crisprnet_loaded_model
    crisprnet_loaded_model = model
