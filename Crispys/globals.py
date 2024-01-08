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

# This variable determines the distance metric that would be used for constructing the distance matrix for the
# target tree. There are four metrics available:
"""
1. all_targets: for a list of n targets, this method transforms each target into a vector of length n, where the jth
entry of the vector i is equal to scoring_function(targets[j],targets[i]).
For example, the vector of the 1st target would look like this:
[scoring_function(targets[0],targets[0]),scoring_function(targets[1],targets[0]),...,scoring_function(targets[n],targets[0])]

2. full_perturbation: for a list of n targets, this method transforms each target into a vector of length of vector_size_cutoff,
entry of the vector i is equal to scoring_function(metric_sequences[j],targets[i]).
metric_sequences is a list of sequences produced by repeatedly sampling a target from the list of targets, and inserting
up to three random mutations.
For example, the vector of the 1st target would look like this:
[scoring_function(metric_sequences[0],targets[0]),scoring_function(metric_sequences[1],targets[0]),...,scoring_function(metric_sequences[n],targets[0])]

3. mixed_perturbation: similar to full_perturbation. Here, if the total number of targets is above vector_size_cutoff,
a number equal to vector_size_cutoff of targets would be picked for the list of metric sequences. Otherwise, the
list of targets will be taken, and the list is completed to vector_size_cutoff by adding perturbed targets as described
in the full_perturbation metric.

4. average_score_linkage: this metric computes the distance between targets A and B by taking the average 
   between scoring_function(A,B) and scoring_function(B,A)
"""
target_distance_metric = "average_score_linkage"
# this variable will determine the number of cores used for the multi_processing used by gold_off.predict
# important: currently, gold off only supports 10 n_cores.
n_cores_for_gold_off = 1

# This variable stores the name of the model used by gold_off_func the file needs to be in the code directory.
# It's used in Distance_matrix_and_UPGMA.gold_off_func:
xgb_model_name = "regression_union_log_max.xgb"

# in case you need to activate conda when using remote debug
CONDA = "source /groups/itay_mayrose/udiland/miniconda3/etc/profile.d/conda.sh; conda activate crispys"

# command that is used to connect to server and run crista
ssh_connect = 'ssh bioseq@powerlogin "module load python/python-anaconda3.6.5 && '


# protdist_path = "/groups/itay_mayrose/josefbrook/miniconda3/envs/crispys/bin/protdist"
protdist_path = 0

# Models of scoring functions loadings
deephf_config = None
deephf_loaded_model = None
crisprnet_loaded_model = None
moff_loaded_model = None
moff_mtx1 = None
moff_mtx2 = None

# Amplicon parameters
safety_padding_around_target = 20


def createHeaderJob(path, job_name, queue="itaym", ncpu=1, mem=16):
    """
    A function to create qsub file with activating crispys conda env

    :param path: path to log files
    :param job_name: job name
    :param queue:
    :param ncpu: cpu number default 1
    :param mem: memory to use (in gb) default 16
    :return: a string that can be used to write sh file to run on the cluster (need to add command before running on the cluster)
    """
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += f"#PBS -q {queue}\n"
    text += "#PBS -N " + job_name + "\n"
    text += "#PBS -e " + path + "/" + job_name + ".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name + ".OU" + "\n"
    text += "#PBS -l select=ncpus=" + str(ncpu) + ":mem=" + str(mem) + "gb\n"
    text += "source ~/.bashrc\n"
    text += "export PATH='$CONDA_PREFIX/bin:$PATH'\n"
    text += "conda activate /groups/itay_mayrose/udiland/miniconda3/envs/crispys\n"
    return text
