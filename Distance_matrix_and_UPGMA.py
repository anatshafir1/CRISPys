
from typing import List, Dict
import numpy as np
from Bio.Phylo import TreeConstruction, BaseTree
import math
import Metric
import TreeConstruction_changed
import gold_off
from os.path import dirname, abspath
import os
import globals
import re
import mafft_and_phylip
from numpy import clip
import subprocess
from CRISPR_Net import Encoder_sgRNA_off
from DeepHF.scripts import prediction_util


def gold_off_func(sg_seq_list: List, target_seq_list: List) -> List[float]:  # Omer Caldararu 28/3
    """
    Scoring function based on gold-off regressor. This function uses a model.xgb file.

    :param sg_seq_list: a list of sgRNA sequence or sequences
    :param target_seq_list: a list of target sequences
    :return: A list of scores where list[i] = score between the ith sgRNA and the ith target
    """
    if len(sg_seq_list[0]) == 20:
        for i in range(len(sg_seq_list)):
            sg_seq_list[i] = sg_seq_list[i] + target_seq_list[i][-3:]
    # get the xgb model output_path
    script_path = dirname(abspath(__file__))
    xgb_model_path = script_path + "/" + globals.xgb_model_name
    list_of_scores = gold_off.predict(sg_seq_list, target_seq_list, xgb_model_path, include_distance_feature=True,
                                      n_process=globals.n_cores_for_gold_off, model_type="regression")
    list_of_scores = clip(list_of_scores, 0, 1)  # clipping is done when the score is above 1 or below 0
    return [1 - score for score in list_of_scores]


def MITScore(seq1: str, seq2: str) -> float:
    """from CRISPR-MIT PAM comes at the end of the string"""
    distance, first_mm, last_mm = 0, -1, -1

    first_arg = 1
    M = [0, 0, 0.14, 0, 0, 0.395, 0.317, 0, 0.398, 0.079, 0.445, 0.508, 0.613, 0.851, 0.723, 0.828, 0.615, 0.804, 0.685,
         0.583]

    if len(seq2) < len(seq1):  # putting the longer sequence as seq2
        temp = seq1
        seq1 = seq2
        seq2 = temp
    distance += len(seq2) - len(seq1)
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
            last_mm = i
            first_arg = first_arg * (1 - M[i])
            if first_mm == -1:
                first_mm = i
    if distance == 0:
        original_score = 1
    else:
        d_avg = (last_mm - first_mm) / distance
        original_score = first_arg / ((4 * (19 - d_avg) / 19 + 1) * distance ** 2)

    return 1 - original_score


def ccTop(sgseq: str, target_seq: str) -> float:
    """

    :param sgseq:
    :param target_seq:
    :return:
    """
    assert len(sgseq) == len(target_seq)
    max_score = sum([math.pow(1.2, i + 1) for i in range(len(sgseq))])
    mm = [i + 1 if sgseq[i] != target_seq[i] else 0 for i in range(len(sgseq))]
    curScore = sum(list(map(lambda x: pow(1.2, x) if x != 0 else 0, mm)))
    return curScore / max_score  # a value between 0 and 1


def ucrispr(sg_seq_list: List) -> List[float]:
    """
    This function will run the uCRISPR algorithm for a list of targets and will return a list of the on-target scores
    IMPORTANT: before running you need to give the path to data tables that are part of uCRISPR
    enter this to the .sh file
    'export DATAPATH=<path to folder>/uCRISPR/RNAstructure/data_tables/'
    Also make sure the uCRISPR file inside the uCRISPR folder has exe permission
    Args:
        sg_seq_list: list of sgrnas (with PAM)

    Returns:
        a list of ucrispr on-target scores
    """
    # make a file with guides for inputs to ucrispr
    with open(f"{globals.CODE_PATH}/targets.txt", "w") as f:
        for sg in sg_seq_list:
            f.write(f"{sg}\n")

    # run ucrispr in terminal
    p = subprocess.run([f"{globals.CODE_PATH}/uCRISPR/uCRISPR", "-on", f"{globals.CODE_PATH}/targets.txt"],
                       stdout=subprocess.PIPE)

    # pares the results
    res_lst = p.stdout.decode('utf-8').split("\n")
    # delete input file
    subprocess.run(["rm", f"{globals.CODE_PATH}/targets.txt"])
    # return a list of the results
    return [float(i.split(" ")[1]) for i in res_lst[1:len(res_lst) - 1]]


def crisprnet(candidate_lst: List, target_lst: List) -> List[float]:
    """
    This function take list of sgrnas (candidate) and lisr of targets and returns a list of  1 - crispr_net score
    Args:
        candidate_lst: list of candidates
        target_lst: list of targets

    Returns:
        list of crispr_net scores
    """
    input_codes = []
    for seqs in zip(candidate_lst, target_lst):
        on_seq = seqs[0]
        off_seq = seqs[1]
        en = Encoder_sgRNA_off.Encoder(on_seq=on_seq, off_seq=off_seq)
        input_codes.append(en.on_off_code)
    input_codes = np.array(input_codes)
    input_codes = input_codes.reshape((len(input_codes), 1, 24, 7))
    y_pred = globals.crisprnet_loaded_model.predict(input_codes).flatten()
    return [1 - float(y) for y in y_pred]


def deephf(target_lst: List) -> List:
    """
    This function use the model of deephf that was improved by Yaron Orenstein`s lab
    Args:
        target_lst: list of targets with PAM

    Returns:
        list of on-target scores
    """
    # take 21 nt from targets
    targets = [target[0:21] for target in target_lst]
    # get deephf scores
    scores = prediction_util.get_predictions(os.path.join(
        f"{globals.CODE_PATH}", "DeepHF", "models", "model1", "no_bio", "multi_task", "parallel/"), targets)
    return list(scores)


def make_upgma(distance_matrix: TreeConstruction._DistanceMatrix) -> TreeConstruction_changed.TreeNew:
    """use by the doc in http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceTreeConstructor-class.html"""
    constructor = TreeConstruction_changed.DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)
    return tree


def make_distance_matrix(names: List, vectors_list: List) -> TreeConstruction._DistanceMatrix:
    """
    Given a list of names of the sequences, a list of vectors where the i-th vector is the position of the i-th target in a
    len(vectors_list) dimensional space, this function returns a distance matrix object, to use in the UPGMA function.

    :param names: list of potential targets sequences
    :param vectors_list: list of lists representing a lower triangular matrix
    :return: distance matrix instance for the UPGMA tree function
    """
    matrix = []
    for i in range(len(vectors_list)):
        row = []
        for j in range(i + 1):
            tempDistance = np.linalg.norm(np.array(vectors_list[i]) - np.array(vectors_list[j]))
            row.append(tempDistance)
        matrix += [row]
    distance_matrix = TreeConstruction._DistanceMatrix(names, matrix)
    return distance_matrix


def make_distance_matrix_from_protdist(output_path: str, names_list: List) -> TreeConstruction._DistanceMatrix:
    """
    this function creates a list of list representing a lower triangular distance matrix using distances created
    by PHYLIP's protDist method. this function is used to make a distance matrix object and to build a UPGMA tree.

    :param output_path: output_path of file with
    :param names_list:
    :return: lower triangular distance matrix
    """
    protdist_outfile = output_path + "/outfile"
    os.rename(protdist_outfile, protdist_outfile + ".txt")
    in_file = open(protdist_outfile + ".txt", 'r')
    p = re.compile("[a-zA-Z][a-zA-Z]")
    temp_res = re.split(p, in_file.read())[1:]
    matrix = []
    i = 1
    for line in temp_res:
        line_split = re.split(" +", line)[1:]
        to_append = list(map(lambda x: float(x.rstrip()), line_split[:i]))
        if to_append:
            matrix.append(to_append)
        i += 1
    distance_matrix = TreeConstruction._DistanceMatrix(names_list, matrix)
    return distance_matrix


def return_protdist_upgma(seq_list, names_list, output_path) -> TreeConstruction_changed.TreeNew:
    """
    this function is called by 'gene_homology_alg' to create a UPGMA tree and a distance matrix using PHYLIP's protDist.
    Given a list of DNA sequences (genes) and their names, this function returns a UPGMA tree and a distance matrix
    based on the genes homology, using MAFFT algorithm to align the genes and PHYLIP's protDist algorithm.

    :param output_path: output_path to phylip format file
    :param seq_list: list of sequence for the UPGMA
    :param names_list: the name of the sequences at seq_list, at the same order
    :return: UPGMA tree of genes and a distance matrix
    """
    new_names_list = list()
    for i in range(len(names_list)):
        new_names_list.append("GG" + str(i))  # for the regex
    mafft_and_phylip.create_protdist(new_names_list, seq_list, output_path)
    distance_matrix = make_distance_matrix_from_protdist(output_path, names_list)
    tree = make_upgma(distance_matrix)
    return tree


def return_targets_upgma(targets_seqs_list: List, names_list: List, scoring_function, cfd_dict: Dict) -> TreeConstruction_changed.TreeNew:
    """the function creates a UPGMA tree object from the potential targets in 'targets_seqs_list' using the given
    'scoring_function', in 3 steps:
    1) for each target in targets_seqs_list and using the given scoring function, creates a vectors list representing
    the location of the target in a multidimensional space.
    2) creates a distance matrix object from the targets names in 'names_list' and the initial matrix.
    3) creates a UPGMA tree object using the distance matrix

    :param targets_seqs_list: list of all the target sequences found in stage0
    :param names_list: a deep copy of targets_list
    :param scoring_function: scoring function of the potential targets
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return: potential targets' tree hierarchically clustered by UPGMA.
    """
    # create a list of vectors for the targets, which is then used to create the distance matrix
    vectors_list = Metric.pos_in_metric_general(targets_seqs_list, scoring_function, cfd_dict)
    # create the distance matrix
    distance_matrix = make_distance_matrix(names_list, vectors_list)
    # apply UPGMA, return a target tree
    targets_tree = make_upgma(distance_matrix)
    return targets_tree
