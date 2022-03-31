import pickle
import numpy as np
from functools import reduce
import Distance_matrix_and_UPGMA
import random
import string
random.seed(1234)

def pos_in_metric_general(t, df, base, cfd_dict = None):
	'''
	:param t: target
	:param df: distance function
	:param base: a list of strings, spreading the space. each string is an sgRNA
	:return: a vector of distances between those strings
	'''
	if cfd_dict:
		Vetorize = np.vectorize(lambda sg: df(sg, t, cfd_dict))
	else:
		Vetorize = np.vectorize(lambda sg: df(sg,t))
	return Vetorize(base)

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


	there is a bug here - the code and the dictinary dose not fit.
	'''
	if not dicti:
		dicti = pickle.load(open("cfd_dict.p",'rb'))
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
	return point

def cfd_funct(sgRNA, target, dicti = None):
	'''my implementation of this function'''
	if not dicti:
		dicti = pickle.load(open("cfd_dict.p",'rb'))

	return 1 - reduce(lambda x, y: x*y, map(lambda i: dicti[('r'+sgRNA[i]+':d'+target[i], i+1)] if sgRNA[i] != target[i] else 1, [j for j in range(0, 20)]))
	#multiply all of this in one frase. didn't did it yet. it is calleed reduce



def find_dist(p1, p2):
	return (sum([(p1[i] - p2[i])**2 for i in range(len(p1))]))**0.5

def find_dist_np(p1, p2):
	"""
	scoring function used when useing the cfd scoring function
	to calculate distance.
	Args:
		p1: vector of length 80 of the 1st seq
		p2: vector of length 80 of the 2nd seq

	Returns: the distance between p1 and p2

	"""
	return np.linalg.norm(p1 - p2)


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

