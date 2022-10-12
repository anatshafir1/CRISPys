"""Vector calculations"""
import numpy as np
from functools import reduce
import Distance_matrix_and_UPGMA
import random
import globals
from typing import List, Dict


def pos_in_metric_general_single_batch(list_of_targets: List[str], metric_seqs_list: List[str],
                                       off_scoring_function, on_scoring_function) -> List[List[float]]:
    """
    This function is called from pos_in_metric_general for cases where the scoring function takes a list of several
    targets, rather than a single target. This function calls the distance_function on all the targets in a single batch.
    input_target_list = [target1, target1, target1,..., target2...]
    input_other_targets_list = [target1, target2, target3, ...]
    concatenated_vectors_list = the result of the distance function, which is then divided into smaller lists.
    output_format : [scoring_function(target_1,target_1),
    scoring_function(target_2,target_1),...,scoring_function(target_j,target_i)]

    :param list_of_targets: a list of targets
    :param metric_seqs_list: a list of sampled sequences used for the metric calculation.
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :return: a vector of distances between those strings
    """
    input_target_list = []
    input_metric_sequences_list = []
    for target in list_of_targets:
        input_target_list.extend([target] * len(metric_seqs_list))
        # update the PAM in each constant target according to the PAM of the target
        metric_seqs_with_pam = [f"{t[:20]}{target[20:]}" for t in metric_seqs_list]
        input_metric_sequences_list.extend(metric_seqs_with_pam)
    vectors_list_off = off_scoring_function(input_metric_sequences_list, input_target_list)
    output_vectors_list = []
    if on_scoring_function == Distance_matrix_and_UPGMA.default_on_target:  # if no on-target scoring function was given
        for i in range(0, len(vectors_list_off), len(metric_seqs_list)):
            output_vectors_list.append(vectors_list_off[i:i + len(metric_seqs_list)])
        return output_vectors_list
    else:
        vectors_list_on = on_scoring_function(list_of_targets)
        for i in range(0, len(vectors_list_off), len(metric_seqs_list)):
            output_vectors_list.append([x * y for x, y in zip(vectors_list_off[i:i + len(metric_seqs_list)], vectors_list_on)])
        return output_vectors_list


def pos_in_metric_general(list_of_targets: List[str], off_scoring_function, on_scoring_function, cfd_dict: Dict) -> List[List[float]]:
    """
    This function takes a list of targets and creates a new list of vectors,
    where each target is converted into a point in number-of-targets dimensional space.
    That way the properties of distance (e.g. symmetry and the triangle inequality)
    are kept when creating the distance matrix.
    e.g. vectors_list[i] = [scoring_function(1,i),scoring_function(2,i),...]

    :param list_of_targets: a list of all targets found in Stage0
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return: a list of vectors, each representing the location of the target in a multidimensional space
    """

    random.seed(globals.seed)
    if off_scoring_function == cfd_funct:
        list_of_vectors = []
        for target in list_of_targets:
            list_of_vectors.append(pos_in_metric_cfd(target, cfd_dict))
        return list_of_vectors
    elif off_scoring_function == Distance_matrix_and_UPGMA.gold_off_func or off_scoring_function == Distance_matrix_and_UPGMA.crisprnet or off_scoring_function == Distance_matrix_and_UPGMA.moff:
        with_pam = True
        metric_sequences_list = create_list_of_metric_sequences(list_of_targets, with_pam)
        return pos_in_metric_general_single_batch(list_of_targets, metric_sequences_list, off_scoring_function, on_scoring_function)
    elif off_scoring_function == Distance_matrix_and_UPGMA.ccTop or off_scoring_function == Distance_matrix_and_UPGMA.MITScore:
        with_pam = False
        metric_sequences_list = create_list_of_metric_sequences(list_of_targets, with_pam)
        if on_scoring_function == Distance_matrix_and_UPGMA.default_on_target:  # if no on-target scoring function was given
            list_of_vectors = []
            for target in list_of_targets:
                score_vector = []
                for sequence in metric_sequences_list:
                    score_vector.append(off_scoring_function(sequence, target[:20]))
                list_of_vectors.append(score_vector)
            return list_of_vectors
        else:
            on_scores = on_scoring_function(list_of_targets)
            list_of_vectors = []
            for target in list_of_targets:
                score_vector = []
                i = 0
                for sequence in metric_sequences_list:
                    score_vector.append(off_scoring_function(sequence, target[:20]) * on_scores[i])
                    i += 1
                list_of_vectors.append(score_vector)
            return list_of_vectors


def create_list_of_metric_sequences(list_of_targets: List[str], with_pam: bool) -> List[str]:
    """
	This function creates a list of sequences of length vector_size_cutoff, which are then used for the distance
	transformation. If the number of targets is smaller than vector_size_cutoff, this function takes all targets,
	and fills in the list by randomly picking a target from the list,
	and appending the input of create_perturbed_target() . This process is repeated until the cutoff size is reached.
	If the number of targets is larger or equal to this variable, shuffle the list of targets and take only the first
	vector_size_cutoff targets.

	:param list_of_targets: a list of targets
	:param with_pam: depending on the scoring function - to include PAM or not
	:return: a list of sampled sequences of length @globals.vector_size_cutoff.
	"""
    if len(list_of_targets) >= globals.vector_size_cutoff:
        shuffled_targets_list = random.sample(list_of_targets, globals.vector_size_cutoff)[:globals.vector_size_cutoff]
        if with_pam:
            return shuffled_targets_list
        else:
            return [t[:20] for t in shuffled_targets_list]
    elif len(list_of_targets) < globals.vector_size_cutoff:
        if with_pam:
            metric_sequences_list = []+list_of_targets
        else:
            metric_sequences_list = [t[:20] for t in list_of_targets]
        for i in range(globals.vector_size_cutoff - len(metric_sequences_list)):
            perturbed_target = create_perturbed_target(random.choice(metric_sequences_list))
            metric_sequences_list.append(perturbed_target)
        return metric_sequences_list


def create_perturbed_target(target: str) -> str:
    """
    This function takes an input target, and creates a perturbed sequence.
    The function randomly chooses between 1 and 3 positions, and inserts
    a single substitution in each mutation index.

    :param target: a target sequence
    :return: a perturbed target.
    """
    number_of_mutations = random.choice(range(1, 4))
    mutation_indices = sorted(random.sample(range(20), number_of_mutations))
    perturbed_target = list(target)
    for j in mutation_indices:
        nucleotide_choices = [char for char in 'ACGT' if char != perturbed_target[j]]
        perturbed_target[j] = random.choice(nucleotide_choices)
    return ''.join(perturbed_target)


def pos_in_metric_cfd(target_seq: str, cfd_dict: Dict) -> List[float]:
    """
    This function is called by "pos_in_metric_general" in "Metric"
    Given a potential sgRNA target sequence and a dictionary of mismatches and their scores for CFD function this
    function builds a vector of length 80 of the scores for mismatches for each nucleotide in the target sequence.
    position 1-4 in the vector are scores for matches/mismatches of the target's first nucleotide with A,C,G,U
    respectively, and so on. a match will always get a score of 1 in the vector
    e.g.: target seq = 'GTGGCAATCCCCATAGACGA' will return
    vector = [0.857142857, 0.714285714, 1.0, 1.0,..., 1.0, 0.7, 0.3, 1.0] where vector[0] is the score for mismatch
    of G with A, vector[1] is the score for mismatch of G with C, and so on.

    :param target_seq: potential target sequence
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function
    :return: a list representing a vector with CFD scores for each nucleotide in the target (see article "Methods", p.9)
    """
    nucs = ['A', 'C', 'G', 'U']
    point = np.zeros(len(target_seq) * len(nucs))
    i = 0
    for pos in range(len(target_seq)):
        for nuc in nucs:
            key = ('r' + nuc + ':d' + target_seq[pos], pos + 1)
            if key in cfd_dict:
                point[i] = cfd_dict[('r' + nuc + ':d' + target_seq[pos], pos + 1)]
            else:
                point[i] = 1
            i += 1
    return list(point)


def cfd_funct(sgRNA: str, target: str, cfd_dict: Dict) -> float:
    """An implementation of the CFD function"""
    return 1 - reduce(lambda x, y: x * y,
                      map(lambda i: cfd_dict[('r' + sgRNA[i] + ':d' + target[i], i + 1)] if sgRNA[i] != target[i] else 1,
                          [j for j in range(0, 20)]))
