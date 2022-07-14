
import pickle
import numpy as np
from functools import reduce
import Distance_matrix_and_UPGMA
import random
from os.path import dirname, abspath
from globals import vector_size_cutoff
from typing import List

random.seed(1234)


def pos_in_metric_general_single_batch(list_of_targets: List, metric_sequences_list: List, distance_function) -> List:
    """
	This function is called from pos_in_metric_general for cases where the scoring function
	takes a list of several targets, rather than a single target.
	This function calls the distance_function on all the targets in a single batch.
	:param list_of_targets: a list of targets.
	:param metric_sequences_list: a list of sampled sequences used for the metric calculation.
	:param distance_function: the distance function.
	:return: a list of vectors, each of length globals.vector_size_cutoff.
	input_target_list = [target1, target1, target1,..., target2 ...]
	input_metric_sequences_list = [metric_sequence1, metric_sequence2, metric_sequence3, ...,metric_sequencei, ...]
	concatenated_vectors_list = the result of the distance function, which is then divided into smaller lists.
	output_format : [distance_function(metric_sequence1,target_1),distance_function(metric_sequence2,target_1),...]
	"""
    input_target_list = []
    input_metric_sequences_list = []
    for target in list_of_targets:
        input_target_list.extend([target] * len(metric_sequences_list))
        # update the PAM in each constant target according to the PAM of the target
        metric_sequences_with_pam = [f"{t}{target[20:]}" for t in metric_sequences_list]
        input_metric_sequences_list.extend(metric_sequences_with_pam)
    concatenated_vectors_list = distance_function(input_metric_sequences_list, input_target_list)
    output_vectors_list = []
    for i in range(0, len(concatenated_vectors_list), len(metric_sequences_list)):
        output_vectors_list.append(concatenated_vectors_list[i:i + len(metric_sequences_list)])
    return output_vectors_list


def pos_in_metric_general(list_of_targets: List, distance_function) -> List:
    """
	This function is used for transforming the scores given by the scoring into vectors that can be used
	to calculate the distances between the targets. That way the properties of distance are kept
	(e.g. symmetry and the triangle inequality).
	These distances are then used for the construction of the target tree.
	The function takes a list of targets and creates a new list of vectors,
	where the ith vector represents the ith target
	and a list of sampled targets with perturbations.
	Args:
		list_of_targets: a list of all targets
		distance_function: the distance function
	Returns: a list of vectors, each representing the location of the target in a
	multidimensional space.
	"""
    if distance_function == cfd_funct:
        list_of_vectors = []
        for target in list_of_targets:
            list_of_vectors.append(pos_in_metric_cfd(target, dicti=None))
        return list_of_vectors
    elif distance_function == Distance_matrix_and_UPGMA.gold_off_func:
        return pos_in_metric_general_single_batch(list_of_targets, create_list_of_metric_sequences(list_of_targets),
                                                  distance_function)
    elif distance_function == Distance_matrix_and_UPGMA.ccTop or distance_function == Distance_matrix_and_UPGMA.MITScore:
        metric_sequences_list = create_list_of_metric_sequences(list_of_targets)
        list_of_vectors = []
        for target in list_of_targets:
            score_vector = []
            for sequence in metric_sequences_list:
                score_vector.append(distance_function(sequence, target))
            list_of_vectors.append(score_vector)
        return list_of_vectors


def create_list_of_metric_sequences(list_of_targets: List) -> List:
    """
	This function creates a list of sequences of length vector_size_cutoff, which are then used for the distance
	transformation. If the number of targets is smaller than vector_size_cutoff, this function takes all targets,
	and fills in the list by randomly picking a target from the list,
	and appending the input of create_perturbed_target() . This process is repeated until the cutoff size is reached.
	If the number of targets is larger or equal to this variable, shuffle the list of targets and take only the first
	vector_size_cutoff targets.
	:param list_of_targets: a list of targets
	:return: a list of sampled sequences of length @globals.vector_size_cutoff.
	"""
    if len(list_of_targets) >= vector_size_cutoff:
        shuffled_targets_list = random.sample(list_of_targets, vector_size_cutoff)[:vector_size_cutoff]
        return [t[:20] for t in shuffled_targets_list]
    elif len(list_of_targets) < vector_size_cutoff:
        metric_sequences_list = [t[:20] for t in list_of_targets]
        for i in range(vector_size_cutoff - len(list_of_targets)):
            perturbed_target = create_perturbed_target(random.choice(list_of_targets))
            metric_sequences_list.append(perturbed_target)
        return metric_sequences_list


def create_perturbed_target(target: str) -> str:
    """
	This function takes an input target, and creates a perturbed sequence.
	The function randomly chooses between 1 and 3 positions, and inserts
	a single subtitution in each mutation index.
	:param target: a target sequence
	:return: a perturbed target.
	"""
    number_of_mutations = random.choice(range(1, 4))
    mutation_indices = sorted(random.sample(range(20), number_of_mutations))
    perturbed_target = list(target)[:20]
    for j in mutation_indices:
        nucleotide_choices = list({'A', 'C', 'G', 'T'} - {perturbed_target[j]})
        perturbed_target[j] = random.choice(nucleotide_choices)
    return ''.join(perturbed_target)


def pos_in_metric_cfd(t, dicti):
    '''
	:param t: target
	 implement a version of the cfd score, in which
	:return:
	there is a bug here - the code and the dictionary do not fit. -which bug?? omer  2/4
	'''
    if not dicti:
        script_path = dirname(abspath(__file__))
        dicti = pickle.load(open(script_path + "/cfd_dict.p", 'rb'))  # added full path to cfd omer 14/4
    Nucs = ['A', 'C', 'G', 'U']
    point = np.zeros(len(t) * len(Nucs))
    i = 0
    for pos in range(len(t)):
        for Nuc in Nucs:
            key = ('r' + Nuc + ':d' + t[pos], pos + 1)
            if key in dicti:
                point[i] = dicti[('r' + Nuc + ':d' + t[pos], pos + 1)]
            else:
                point[i] = 1
            i += 1
    return list(point)


def cfd_funct(sgRNA, target, dicti=None):
    '''my implementation of this function'''
    if not dicti:
        script_path = dirname(abspath(__file__))
        dicti = pickle.load(open(script_path + "/cfd_dict.p", 'rb'))  # added full path to cfd omer 7/4

    return 1 - reduce(lambda x, y: x * y,
                      map(lambda i: dicti[('r' + sgRNA[i] + ':d' + target[i], i + 1)] if sgRNA[i] != target[i] else 1,
                          [j for j in range(0, 20)]))