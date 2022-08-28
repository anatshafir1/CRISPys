"""
this script is for storing global variables
"""
import os

# get the path of the scripts directory
CODE_PATH = os.path.dirname(os.path.realpath(__file__))

# this variable will determine the maximum number of genes that will be considered for finding
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

# command that is used to connect to server and run crista
ssh_connect = 'ssh bioseq@powerlogin "module load python/python-anaconda3.6.5 && '

mafft_path = "/groups/itay_mayrose/josefbrook/miniconda3/envs/crispys/bin/mafft"
protdist_path = "/groups/itay_mayrose/josefbrook/miniconda3/envs/crispys/bin/protdist"
# mafft_path = None
# protdist_path = None

crisprnet_loaded_model = None
moff_loaded_model = None
moff_mtx1 = None
moff_mtx2 = None
