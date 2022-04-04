__author__ = 'ItayM5'
import re ##used to be regex
import regex
import Metric
import Distance_matrix_and_UPGMA

def get_sites(gene, df, min_length=20, max_length=20, start_with_G=False, where_in_gene = 1):
	'''
	:param gene:
	:param min_length:
	:param max_length:
	:param start_with_G:
	:param where_in_gene: forword to this position the sgRNA are ignored
	:return:
	'''
	list_of_targets = []
	if len(gene) < max_length+3:
		return list_of_targets
	for length in range(min_length, max_length +1):
		if (start_with_G):
			SiteAndPAM = "G" + "."*length + "GG" #it is acually NGG
		else:
			SiteAndPAM = "."*(length +1) + "GG" #it is acually NGG
		compiled = regex.compile(SiteAndPAM)
		where_in_gene = int(len(gene)*where_in_gene)
		founds_sense = regex.findall(compiled, gene[:where_in_gene], overlapped=True)
		founds_antisense = regex.findall(compiled, give_complementary(gene)[:where_in_gene], overlapped=True)
		# functions that take targets of length 20
		if df == Metric.cfd_funct or df == Distance_matrix_and_UPGMA.ccTop or df == Distance_matrix_and_UPGMA.MITScore:
			founds = [seq[:-3] for seq in founds_sense if 'N' not in seq[:-3]] + [seq[:-3] for seq in founds_antisense if 'N' not in seq[:-3]]
		# functions that take targets of length 23
		elif df == Distance_matrix_and_UPGMA.gold_off_func:
			founds = [seq for seq in founds_sense if 'N' not in seq] + [seq for seq in founds_antisense if 'N' not in seq]
		list_of_targets += founds
	#print(res)
	return list_of_targets


def give_complementary(seq):
	complementary_seq_list = []
	for i in range(len(seq)) :
		if seq[len(seq)-1-i] == 'A':
			complementary_seq_list.append('T')
		elif seq[len(seq)-1-i] == 'T':
			complementary_seq_list.append('A')
		elif seq[len(seq)-1-i] == 'C':
			complementary_seq_list.append('G')
		elif seq[len(seq)-1-i] == 'G':
			complementary_seq_list.append('C')
		elif seq[len(seq)-1-i] == 'N':
			complementary_seq_list.append('N')
	return ''.join(complementary_seq_list)

def give_complementary_old(seq): #can be removed ?
	complementary_seq_list = []
	for letter in seq:
		if letter == 'A':
			complementary_seq_list.append('T')
		elif letter == 'T':
			complementary_seq_list.append('A')
		elif letter == 'C':
			complementary_seq_list.append('G')
		elif letter == 'G':
			complementary_seq_list.append('C')
		elif letter == 'N':
			complementary_seq_list.append('N')
	return ''.join(complementary_seq_list)

def find_offtagrets(seq, chromo_folder):
	'''
	:param seq:
	:param chromo_folder: a folder in which there are exactly all the chromosomse
	:return:
	'''

def get_targets_sites_from_exons_lst(exons_lst, df, original_range_in_gene = [0,1], min_length= 20, max_length = 20,start_with_G = False):
	if original_range_in_gene[1] <= original_range_in_gene[0]:
		print("The range of the targts on the gene is not in the right format")
		exit(-1)
	if max_length < min_length:
		print("The range of the lengths of the sgRNA is not in the right format")
		exit(-1)
	list_of_targets = []
	lengths = [len(exon) for exon in exons_lst]
	gene_length = sum(lengths)
    #where in gene - used for deciding what parts to consider in the gene
	range_in_gene = [int(r*gene_length) for r in original_range_in_gene]
	##exons_lst = list(map(lambda seq: seq.upper(), exons_lst)) #might be unneccesary

	#accumolate the length of exons
	for i in range(1, len(lengths)):
		lengths[i] = lengths[i-1] + lengths[i]
	for i in range(len(exons_lst)):
		if i == 0:
			if range_in_gene[0] < lengths[i]:
				#if range_in_gene[1]*gene_length > lengths[i]:
				list_of_targets += get_sites(exons_lst[i][range_in_gene[0] : min(lengths[i], range_in_gene[1])], df, min_length, max_length, start_with_G, where_in_gene = 1)
		elif max(range_in_gene[0], lengths[i-1]) < min(lengths[i], range_in_gene[1]):
			list_of_targets += get_sites(exons_lst[i][max(range_in_gene[0]  - lengths[i-1], 0) : min(lengths[i] - lengths[i-1], range_in_gene[1] - lengths[i-1])], df, min_length, max_length, start_with_G, where_in_gene = 1)
	#print(list_of_targets)
	return list_of_targets
