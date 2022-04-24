import pickle
import numpy as np
from functools import reduce
import Distance_matrix_and_UPGMA
import random
from os.path import dirname, abspath, isfile
import string

import Metric

random.seed(1234)


def pos_in_metric_general_single_batch(list_of_targets, constant_targets_list, distance_function):
    """
	This function takes a list of targets and creates a new list of vectors,
	where each target is converted into a point in number-of-targets dimensional space.
	This function calls the distance_function a single time.
	:param list_of_targets: a list of targets
	:param distance_function: the distance function
	:param constant_targets_list: a list of constant targets
	:return: a vector of distances between those strings
	input_target_list = [target1, target1, target1,..., target2...]
	input_other_targets_list = [constant_target1, constant_target2, constant_target3, ...]
	concatenated_vectors_list = the result of the distance function, which is then divided into smaller lists.
	output_format : [distance_function(target_1,target_1),distance_function(target_2,target_1),...,distance_function(target_j,target_i)]
	"""
    input_target_list = []
    input_other_targets_list = []
    for target in list_of_targets:
        input_target_list.extend([target] * len(constant_targets_list))
		# update the PAM in each constant target according to the PAM of the target
        constant_targets_with_pam = [f"{t}{target[20:]}" for t in constant_targets_list]
        input_other_targets_list.extend(constant_targets_with_pam)
    concatenated_vectors_list = distance_function(input_other_targets_list, input_target_list)
    output_vectors_list = []
    for i in range(0, len(concatenated_vectors_list), len(constant_targets_list)):
        output_vectors_list.append(concatenated_vectors_list[i:i + len(constant_targets_list)])
    return output_vectors_list


def pos_in_metric_general(list_of_targets, distance_function):
    """
    This function is used for transforming the scores given by the scoring into vectors that can be used
    to calculate the distances between the targets. That way the properties of distance are kept
    (e.g. symmetry and the triangle inequality).
    These distances are then used for the construction of the target tree.
	The function takes a list of targets and creates a new list of vectors,
	where each target is a vector of distances between a target
	 and a list of constant targets.
	are kept when creating the distance matrix.
	Args:
		list_of_targets: a list of all targets
		distance_function: the distance function
	Returns: a list of vectors, each representing the location of the target in a
	multidimensional space of size len(constant_targets_list)
	vectors_list[i] = [distance_function(constant_target1,target_i),distance_function(constant_target2,target_i),...]
	"""
    script_path = dirname(abspath(__file__))
    with open(script_path + "/list_of_constant_targets.txt") as f:
        constant_targets_list = f.read().splitlines()
    if distance_function == Distance_matrix_and_UPGMA.gold_off_func:
        return pos_in_metric_general_single_batch(list_of_targets, constant_targets_list, distance_function)
    elif distance_function == cfd_funct:
        list_of_vectors = []
        for target in list_of_targets:
            list_of_vectors.append(pos_in_metric_cfd_np(target, dicti=None))
        return list_of_vectors
    elif distance_function == Distance_matrix_and_UPGMA.ccTop or distance_function == Distance_matrix_and_UPGMA.MITScore:
        list_of_vectors = []
        for target in list_of_targets:
            target_vector = []
            for constant_target in constant_targets_list:
                target_vector.append(distance_function(constant_target, target))
            list_of_vectors.append(target_vector)
        return list_of_vectors

def pos_in_metric_cfd(t, cfd_dict = None):
	'''
	:param t: target
	 construct a distance vector for a given target, using the cfd score dict
	 example of a value in the cfd dic: ('rT:dA', 17): 0.6 means that the score of
	 a mismatch between T (in Nucs) and A (in the target) at position 17 is 0.6
	 @point example for a target with a G at position 1: [0.857142857,0.714285714,1,0.857142857,(...the scores in the other positions...)]
	 the length of the vector should be 80 if the target length is 20
	:return:
	'''
	if not dicti:
		dicti = pickle.load(open("cfd_dict.p",'rb'))
	Nucs = ['A','C','G', 'U']
	point = [0 for i in range(len(t)*len(Nucs))]
	i=0
	for pos in range(len(t)):
		for Nuc in Nucs:
			key = ('r'+ Nuc +':d'+ t[pos], pos+1)
			if key in dicti:
				point[i] = dicti[('r'+ Nuc +':d'+ t[pos], pos+1)]
			else:
				point[i] = 1
			i += 1
	return point


def pos_in_metric_cfd_np(t, dicti):
	'''
	:param t: target
	 implement a version of the cfd score, in which
	:return:
	there is a bug here - the code and the dictinary do not fit. -which bug?? omer  2/4
	'''
	if not dicti:
		script_path = dirname(abspath(__file__))
		dicti = pickle.load(open(script_path+"/cfd_dict.p",'rb')) #added full path to cfd omer 14/4
	Nucs = ['A','C','G', 'U']
	point = np.zeros(len(t)*len(Nucs))
	#point = [0 for i in range(len(t)*len(Nucs))]
	i=0
	for pos in range(len(t)):
		for Nuc in Nucs:
			key = ('r'+ Nuc +':d'+ t[pos], pos+1)
			if key in dicti:
				point[i] = dicti[('r'+ Nuc +':d'+ t[pos], pos+1)]
			else:
				point[i] = 1
			i += 1
	return list(point)

def cfd_funct(sgRNA, target, dicti = None):
	'''my implementation of this function'''
	if not dicti:
		script_path = dirname(abspath(__file__))
		dicti = pickle.load(open(script_path+"/cfd_dict.p",'rb')) #added full path to cfd omer 7/4

	return 1 - reduce(lambda x, y: x*y, map(lambda i: dicti[('r'+sgRNA[i]+':d'+target[i], i+1)] if sgRNA[i] != target[i] else 1, [j for j in range(0, 20)]))
	#multiply all of this in one frase. didn't did it yet. it is calleed reduce



def find_dist(p1, p2):
	return (sum([(p1[i] - p2[i])**2 for i in range(len(p1))]))**0.5


def find_dist_np(p1, p2):
	"""
	scoring function used when using the cfd scoring function
	to calculate distance.
	Args:
		p1: vector of the 1st seq
		p2: vector of the 2nd seq

	Returns: the distance between p1 and p2
	"""
	return np.linalg.norm(np.array(p1) - np.array(p2))



def find_dist_t(t1, t2, cfd_dict = None):
	p1, p2 = pos_in_metric_cfd(t1, cfd_dict), pos_in_metric_cfd(t2, cfd_dict)
	return find_dist(p1, p2)



def make_pos_dict(inpath = "D:\\Lab\\Cdata\\Relevant articles\\STable 19 FractionActive_dlfc_lookup.txt"):
	'''
	the dictionary manufacter here is sutable for comparing the match when the RNA sequence is represented as it's complementary
	'''
	give_compl = lambda x: 'G' if x == 'C' else 'C' if x == 'G' else 'T' if x == 'A' else 'A'
	infile = open(inpath, 'r')
	next(infile)
	dicti = dict()

	for line in infile:
		line_as_array = line.split('\t')
		line_as_array[0] = 'r' +  give_compl(line_as_array[0][1]) + line_as_array[0][2:]
		type, pos, score = line_as_array[0], line_as_array[1], line_as_array[5]
		dicti[(type, int(pos))] = float(score)
	pickle.dump(dicti, open("cfd_dict.p",'wb'))

