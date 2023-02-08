"""Scoring functions and tree calculations"""
import subprocess
from os.path import dirname, abspath, join
from os import rename
import math
from typing import List, Dict
import itertools
import numpy as np
import re
from Bio.Phylo.TreeConstruction import _DistanceMatrix

import globals
from mafft_and_phylip import create_protdist
from CRISPR_Net import Encoder_sgRNA_off
from DeepHF.scripts import prediction_util
from MOFF.MOFF_prediction import MOFF_score
from Metric import pos_in_metric_general, cfd_funct
from TreeConstruction_changed import TreeNew, DistanceTreeConstructor
from gold_off import predict


# ###################################### off target functions ###################################### #


def MITScore(seq1: str, seq2: str, for_metric: bool = False) -> float:
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
    if for_metric:
        return 1 - original_score
    else:
        return original_score


def ccTop(sgseq: str, target_seq: str, for_metric: bool = False) -> float:
    """

    :param sgseq:
    :param target_seq:
    :param for_metric: True for the score for distance calculations, False for the reciprocal of the score
    (ccTOP, The efficiency score is the reciprocal of curScore/max_score)
    :return:
    """
    assert len(sgseq) == len(target_seq)
    max_score = sum([math.pow(1.2, i + 1) for i in range(len(sgseq))])
    mm = [i + 1 if sgseq[i] != target_seq[i] else 0 for i in range(len(sgseq))]
    curScore = sum(list(map(lambda x: pow(1.2, x) if x != 0 else 0, mm)))
    if for_metric:
        return curScore / max_score
    else:
        return 1 - curScore / max_score


def gold_off_func(sg_seq_list: List[str], target_seq_list: List[str], for_metric: bool = False) -> List[
    float]:  # Omer  28/3
    """
    Scoring function based on gold-off regressor. This function uses a model.xgb file.

    :param sg_seq_list: a list of sgRNA sequence or sequences
    :param target_seq_list: a list of target sequences
    :param for_metric: True for the reciprocal of the score for distance calculations, False for the score itself
    :return: A list of scores where list[i] = score between the ith sgRNA and the ith target
    """
    if len(sg_seq_list[0]) == 20:
        for i in range(len(sg_seq_list)):
            sg_seq_list[i] = sg_seq_list[i] + target_seq_list[i][-3:]
    # get the xgb model output_path
    script_path = dirname(abspath(__file__))
    xgb_model_path = script_path + "/" + globals.xgb_model_name
    list_of_scores = predict(sg_seq_list, target_seq_list, xgb_model_path, include_distance_feature=True,
                             n_process=globals.n_cores_for_gold_off, model_type="regression")
    list_of_scores = np.clip(list_of_scores, 0, 1)  # clipping is done when the score is above 1 or below 0
    if for_metric:
        return [1 - float(score) for score in list_of_scores]
    else:
        return [float(score) for score in list_of_scores]


def crisprnet(candidate_lst: List[str], target_lst: List[str], for_metric: bool = False) -> List[float]:
    """
    This function take list of sgrnas (candidate) and list of targets and returns a list of  1 - crispr_net score

    :param candidate_lst: list of candidates
    :param target_lst: list of targets
    :param for_metric: True for the reciprocal of the score for distance calculations, False for the score itself
    :return: list of crispr_net scores (or their reciprocals if for_metric is false)
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
    if for_metric:
        return [1 - float(y) for y in y_pred]
    else:
        return list(y_pred)


def moff(candidate_lst: List[str], target_lst: List[str], for_metric: bool = False) -> List[float]:
    """
    Calling MOFF algorithm, this function take list of sgrnas (candidate) and list of targets and
    returns a list of  1 - MOFF score

    :param candidate_lst: list of candidates
    :param target_lst: list of targets
    :param for_metric: True for the reciprocal of the score for distance calculations, False for the score itself
    :return: list of scores (or 1-scores for metric calculation)
    """
    scores = MOFF_score(globals.moff_mtx1, globals.moff_mtx2, candidate_lst, target_lst)
    if for_metric:
        return [1 - score for score in scores]
    else:
        return list(scores)


# ###################################### on target functions ###################################### #

def deephf(target_lst: List[str], for_metric: bool = False) -> List[float]:
    """
    This function use the model of deephf that was improved by Yaron Orenstein`s lab
    :param target_lst: list of targets with PAM
    :param for_metric: True for the reciprocal of the score for distance calculations, False for the score itself
    :return: list of on-target scores
    """
    # take 21 nt from targets
    targets = [target[0:21] for target in target_lst]
    # get deephf scores
    scores = prediction_util.get_predictions(globals.deephf_loaded_model, globals.deephf_config, targets)
    if for_metric:
        return [1 - score for score in scores]
    else:
        return list(scores)


def ucrispr(sg_seq_list: List[str], output_path: str) -> List[float]:
    """
    This function will run the uCRISPR algorithm for a list of targets and will return a list of the on-target scores
    IMPORTANT: before running you need to give the path to data tables that are part of uCRISPR
    enter this to the .sh file or to the bashrc file
    'export DATAPATH=<path to folder>/uCRISPR/RNAstructure/data_tables/'
    Also make sure the uCRISPR file inside the uCRISPR folder has exe permission

    :param sg_seq_list: list of sgrnas (with PAM)
    :param output_path: the output file
    :return: a list of ucrispr on-target scores
    """
    # make a file with guides for inputs to ucrispr
    with open(f"{output_path}/targets.txt", "w") as f:
        for sg in sg_seq_list:
            f.write(f"{sg}\n")
    # run ucrispr in terminal
    p = subprocess.run([f"{globals.CODE_PATH}/uCRISPR/uCRISPR", "-on", f"{output_path}/targets.txt", f"{output_path}/"],
                       stdout=subprocess.PIPE)
    # parse the results
    res_lst = p.stdout.decode('utf-8').split("\n")
    with open(f"{output_path}/ucrisper_out.txt", "w") as u_out:
        for whatever in res_lst:
            u_out.write(f"{whatever}\n")
    # delete input file
    subprocess.run(["rm", f"{output_path}/targets.txt"])
    # return a list of the results
    return [float(i.split(" ")[1]) for i in res_lst[1:len(res_lst) - 1]]


def default_on_target(target_lst: List[str], for_metric: bool = False) -> List[int]:
    """
    This function is used in the case when no on-target scoring function is chosen. Returns a list of ones in the length
    of the given target list.
    :param target_lst: list of targets with PAM
    :param for_metric: True for the reciprocal of the score for distance calculations, False for the score itself
    :return: list of ones in the length of the given target list.
    """
    result = [1 for _ in target_lst]
    return result


# ###################################### UPGMA tree functions ###################################### #

def make_upgma(distance_matrix: _DistanceMatrix) -> TreeNew:
    """use by the doc in http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceTreeConstructor-class.html"""
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)
    return tree


def make_distance_matrix(names: List, vectors_list: List) -> _DistanceMatrix:
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
    distance_matrix = _DistanceMatrix(names, matrix)
    return distance_matrix


def make_distance_matrix_from_protdist(output_path: str, names_list: List) -> _DistanceMatrix:
    """
    this function creates a list of list representing a lower triangular distance matrix using distances created
    by PHYLIP's protDist method. this function is used to make a distance matrix object and to build a UPGMA tree.

    :param output_path: output_path of file with
    :param names_list:
    :return: lower triangular distance matrix
    """
    protdist_outfile = output_path + "/outfile"
    rename(protdist_outfile, protdist_outfile + ".txt")
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
    distance_matrix = _DistanceMatrix(names_list, matrix)
    return distance_matrix


def return_protdist_upgma(seq_list: List[str], names_list: List[str], output_path: str) -> TreeNew:
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
    create_protdist(new_names_list, seq_list, output_path)
    distance_matrix = make_distance_matrix_from_protdist(output_path, names_list)
    tree = make_upgma(distance_matrix)
    return tree


def make_distance_matrix_from_average(targets_seqs_list: List[str], names_list: List[str],
                                      off_scoring_function) -> _DistanceMatrix:
    """
    This function creates a pseudo-distance matrix based on the average of the scores between target pairs:
    DistanceMatrix[i][j] = (score(targets[i],targets[j]) + score(targets[j],targets[i]))/2
    if score(A,B) is 0.2 and score(B,A) is 0.4, then the distance between A and B is 0.3
    :param targets_seqs_list: a list of target sequences
    :param names_list: a deep copy of targets_list
    :param off_scoring_function: the off target scoring function
    :return: a _DistanceMatrix object
    """
    targets_list_1 = []
    targets_list_2 = []
    for target_pair in itertools.product(targets_seqs_list, repeat=2):  # prepare the input for the scoring function
        targets_list_1.append(target_pair[0])
        targets_list_2.append(target_pair[1])
    list_of_scores = off_scoring_function(targets_list_1, targets_list_2, for_metric=True)
    n = len(targets_seqs_list)
    chunks_list = [list_of_scores[i:i + n] for i in range(0, len(list_of_scores), n)]  # divide the scores into chunks
    distance_matrix = []
    for i in range(n):
        matrix_row = []
        for j in range(i + 1):
            distance = (chunks_list[i][j] + chunks_list[j][i]) / 2  # calculate the average
            matrix_row.append(distance)
        distance_matrix.append(matrix_row)
    return _DistanceMatrix(names_list, distance_matrix)


def return_targets_upgma(targets_seqs_list: List[str], names_list: List[str], off_scoring_function, on_scoring_function,
                         cfd_dict: Dict) -> TreeNew:
    """the function creates a UPGMA tree object from the potential targets in 'targets_seqs_list' using the given
    'scoring_function', in 3 steps:
    1) for each target in targets_seqs_list and using the given scoring function, creates a vectors list representing
    the location of the target in a multidimensional space.
    2) creates a distance matrix object from the targets names in 'names_list' and the initial matrix.
    3) creates a UPGMA tree object using the distance matrix

    :param targets_seqs_list: list of all the target sequences found in stage0
    :param names_list: a deep copy of targets_list
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return: potential targets' tree hierarchically clustered by UPGMA.
    """
    if off_scoring_function == cfd_funct:
        # create a list of vectors for the targets, which is then used to create the distance matrix
        vectors_list = pos_in_metric_general(targets_seqs_list, off_scoring_function, on_scoring_function, cfd_dict)
        # create the distance matrix
        distance_matrix = make_distance_matrix(names_list, vectors_list)
        # apply UPGMA, return a target tree
    else:
        distance_matrix = make_distance_matrix_from_average(targets_seqs_list, names_list, off_scoring_function)
    targets_tree = make_upgma(distance_matrix)
    return targets_tree
