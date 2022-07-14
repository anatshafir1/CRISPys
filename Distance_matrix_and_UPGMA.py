
from typing import List
import numpy as np
from Bio.Phylo import TreeConstruction
import math
import TreeConstruction_changed
import gold_off
from os.path import dirname, abspath
import os
import globals
import re
import mafft_and_phylip
from numpy import clip


def gold_off_func(sg_seq_list, target_seq_list, dicti=None):  # Omer Caldararu 28/3
    """
    Scoring function based on gold-off regressor. This function uses a model
    .xgb file.
    Args:
        sg_seq_list: a list of sgRNA sequence or sequences
        target_seq_list: a list of target sequences
        dicti: None

    Returns: A list of scores where list[i] = score between the ith sgRNA and the ith target

    """
    if len(sg_seq_list[0]) == 20:
        for i in range(len(sg_seq_list)):
            sg_seq_list[i] = sg_seq_list[i] + target_seq_list[i][-3:]
    # get the xgb model path
    script_path = dirname(abspath(__file__))
    xgb_model_path = script_path + "/" + globals.xgb_model_name
    list_of_scores = gold_off.predict(sg_seq_list, target_seq_list, xgb_model_path, include_distance_feature=True,
                                      n_process=globals.n_cores_for_gold_off, model_type="regression")
    list_of_scores = clip(list_of_scores, 0, 1)  # clipping is done when the score is above 1 or below 0
    return [1 - score for score in list_of_scores]


def MITScore(seq1, seq2, dicti=None):
    '''from CRISPR-MIT
    PAM cames at the end of the string'''
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


def ccTop(sgseq, target_seq, dicti=None):
    assert len(sgseq) == len(target_seq)
    max_score = sum([math.pow(1.2, i + 1) for i in range(len(sgseq))])
    mm = [i + 1 if sgseq[i] != target_seq[i] else 0 for i in range(len(sgseq))]
    curScore = sum(list(map(lambda x: pow(1.2, x) if x != 0 else 0, mm)))
    return curScore / max_score  # a value between 0 and 1


def make_UPGMA(dm):
    '''use by the doc in http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceTreeConstructor-class.html'''
    constructor = TreeConstruction_changed.DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    return tree


def make_distance_matrix(names, initial_matrix):
    '''input: list of names of the sequences, and the output of 'make_initial_matrix'
    output: a distance matrix, in a format adequate to the UPGMA function'''
    m = TreeConstruction._DistanceMatrix(names, initial_matrix)
    return m


def make_initial_matrix(vectors_list):
    """
    calculates a distance matrix for a list of points. The matrix is then used to construct the target tree.
    Args:
        vectors_list: a list of vectors where the ith vector is the position of the ith target
        in a len(vectors_list) dimensional space.
    Returns: a triangular distance matrix where matrix[i][j] <- distance(vector_i,vector_j)
    """
    distance_matrix = []
    for i in range(len(vectors_list)):
        row = []
        for j in range(i + 1):
            tempDistance = np.linalg.norm(np.array(vectors_list[i]) - np.array(vectors_list[j]))
            row.append(tempDistance)
        distance_matrix += [row]
    return distance_matrix


def make_initial_matrix_from_protdist(output_path: str):
    """
    this function creates a list of list representing a lower triangular distance matrix using distances created
    by PHYLIP's protDist method. this function is used to make a distance matrix object and to build a UPGMA tree.

    :param output_path: output_path of file with
    :return: lower triangular distance matrix
    :rtype: list
    """
    protdist_outfile = output_path + "/outfile"
    os.rename(protdist_outfile, protdist_outfile + ".txt")
    in_file = open(protdist_outfile + ".txt", 'r')
    p = re.compile("[a-zA-Z][a-zA-Z]")
    temp_res = re.split(p, in_file.read())[1:]
    res = []
    i = 1
    for line in temp_res:
        line_split = re.split(" +", line)[1:]
        to_append = list(map(lambda x: float(x.rstrip()), line_split[:i]))
        if to_append:
            res.append(to_append)
        i += 1
    return res


def return_protdist_upgma(seq_list, names_list, output_path):
    '''
    with protdist
    :param seq_list: list of sequence for the UPGMA
    :param names_list: the name of the sequences at seq_list, at the same order
    :return: UPGMA tree
    '''
    new_names_list = list()
    for i in range(len(names_list)):
        new_names_list.append("GG" + str(i))  # for the regex
    mafft_and_phylip.create_protdist(new_names_list, seq_list, output_path)  # to uncomment when really running the tool on the server
    matrix = make_initial_matrix_from_protdist(output_path)
    distance_matrix = make_distance_matrix(names_list, matrix)
    return make_UPGMA(distance_matrix), distance_matrix
