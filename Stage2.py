import UPGMA
import Stage3
import copy
import math
import Candidate
import Metric



def print_res_to_file(res, input_sg_genes_dict, Omega):
	''' input(res) format:	array of permutations_DS'''
	##first, get genes set, for the file name:
	genes_set = set()
	for sg, genes in input_sg_genes_dict.items():  #genes  is a list of genes
		for gene in genes:
			genes_set.add(gene.split("RNA")[0])
	file_name = str(genes_set)[1:-1] + ".txt"
	f = open(file_name, "w")
	for permDS in res:
		f.write("SG:\n")
		f.write(permDS[0])
		f.write("\n")
		f.write("score:\n")
		f.write(str(permDS[2]))
		f.write("\n")
		f.write("cleavage site for each gene:\n")
		next_line = str(permDS[4])
		next_line = next_line.replace("RNAfile","")
		f.write(next_line)
		f.write("\n\n")
	f.close()

def stopping_condition(lowest_widest_prob, Omega, num_of_sites = 11):
	'''should continue going up the tree?'''
	return lowest_widest_prob < (1 - (1 - Omega)**(1/num_of_sites))

	#return lowest_widest_prob < (1 - math.sqrt(1-Omega))
def top_down(best_permutations_DS, node, Omega, sg_genes_dict, targets_df, cfd_dict = None, PS_number = 12):
	'''
	:param node:
	:param current_genes_sg_dict:
	:param Omega: can be removed already??
	:return:
	'''
	if len(node.polymorphic_sites_set) < PS_number: # was hardcoded to be 12, change it to be the PS_number argument Udi 24/02/22
		#make current_genes_sg_dict
		current_genes_sg_dict = dict()
		for target in node.node_targets_DS:
			genes_leaf_from = sg_genes_dict[target]  #which gene is this target came from. usually it will be only 1 gene
			for gene_name in genes_leaf_from:
				if gene_name in current_genes_sg_dict:
					if target not in current_genes_sg_dict[gene_name]:
						current_genes_sg_dict[gene_name] += [target]
				else:
					current_genes_sg_dict[gene_name] = [target]
		current_res = Stage3.find_Uno_sgRNA(current_genes_sg_dict, Omega, targets_df, node, cfd_dict, PS_number) #current best perm is a tuple with the perm and metedata of this perm. in this option, node.candidtes_DS is updated in the Naive
		#print(current_genes_sg_dict)
		if current_res:
			best_permutations_DS  += current_res
		return
	else:
		top_down(best_permutations_DS, node.clades[0], Omega, sg_genes_dict, targets_df, cfd_dict, PS_number)
		top_down(best_permutations_DS, node.clades[1], Omega, sg_genes_dict, targets_df, cfd_dict, PS_number)


def bottem_up(node, current_genes_sg_dict, Omega, targets_df, cfd_dict = None):
	'''check find the key seq for the group here, how good it is, and continue going up the tree if stopping_condition() sais so'''
	if node.colour == 'b':
		return
	if not(current_genes_sg_dict):
		current_genes_sg_dict = {}
	for target in node.node_targets_DS:
		genes_leaf_from = sg_genes_dict[target]  #which gene is this target came from. usually it will be only 1 gene
		for gene_name in genes_leaf_from:
			if gene_name in current_genes_sg_dict:
				if target not in current_genes_sg_dict[gene_name]:
					current_genes_sg_dict[gene_name] += [target]
			else:
				current_genes_sg_dict[gene_name] = [target]

	#reuse these two lines for running time optimization after it will all work well
	if (len(current_genes_sg_dict) == 1):  #only 1 gene.
		current_best_perm = find_best_sg_for_single_gene(gene_name, current_genes_sg_dict[gene_name])  #(max_seq, max_fraction, max_cut_prob, genes_list, match_sites_list])

	else:
		current_best_perm = Stage3.find_Uno_sgRNA(current_genes_sg_dict, Omega, targets_df, node, cfd_dict) #current best perm is a tuple with the perm and metedata of this perm. in this option, node.candidtes_DS is updated in the Naive
	if current_best_perm == None:
		return  #end of the search on this brance
	##continue up##
	for item in current_best_perm:  #can be writen in a more complexity efficient way. maybe later.
		if (item):
			best_permutations_DS.append(item) #(perm, fraction genes being cut among all the genes, probability to cut all the genes in genes list, genes_list)
	node.set_colour('b')
	if (node.parent):
		bottem_up(node.parent, current_genes_sg_dict, Omega, target_df, cfd_dict)



def bottem_up_tree(upgmaTree, Omega, target_df, cfd_dict):
	for i in range(len(upgmaTree.leaves_DS)):
		bottem_up(upgmaTree.leaves_DS[i], None, Omega, target_df, cfd_dict)


def call_it_all(sgList, sgNames, input_sg_genes_dict, Omega, df_targets, cfd_dict = None, PS_number = 12):
	best_permutations_DS = []
	if len(sgList) == 1:
		print("only one sgRNA in the group")
		genes = input_sg_genes_dict[sgList[0]]
		c = Candidate.Candidate(sgList[0])
		c.fill_default_fildes(genes)
		best_permutations_DS.append(c)
	else:
        #create the distance matrix of sgRNAs and perform UPGMA on it
		upgmaTree = return_upgma(sgList,sgNames, df_targets, cfd_dict)
        #
		fill_leaves_sets(upgmaTree, input_sg_genes_dict)
		fill_PS(upgmaTree.root)
		top_down(best_permutations_DS,upgmaTree.root, Omega, input_sg_genes_dict, df_targets, cfd_dict, PS_number)
	return  best_permutations_DS



def find_set_cover():
	'''for now, might won't work in a case when there is a gene that isn't covered by any of the permutations in the best_permutations_DS. not finished. can make it more readble'''
	res = [best_permutations_DS[0]]
	temp_best_perm_DS = copy.copy(best_permutations_DS)
	uncovered_genes = set()
	for sg, genesLst in sg_genes_dict.items():
		for gene in genesLst:
			uncovered_genes.add(gene)
	for gene in res[0][3]:  #the genes_list_of the best perm
		uncovered_genes.remove(gene)  ##still need to verify it's in the same format as were added in the uncovered_genes_set
	while(len(uncovered_genes)) > 0 and len(temp_best_perm_DS) > 0:
	##going over all the permutations, and return the permutation that cover the maximal amount of genes haven't been covered yet, in the highest probability among the maximal covered permutations
		best_current_perm, best_num_of_coverd, best_prob_of_covered = None, 0,0  #best_current_perm is the hole tuple
		i = 0
		while i < (len(temp_best_perm_DS)):
			num_of_coverd = 0
			for gene in temp_best_perm_DS[i][3]:
				if gene in uncovered_genes:
					num_of_coverd += 1
			if num_of_coverd == 0:
				del temp_best_perm_DS[i]
			elif num_of_coverd >= best_num_of_coverd:## and temp_best_perm_DS[i][2] > best_prob_of_covered:  ##need to check if 2 is the right index, and not 1.
				best_num_of_coverd, best_prob_of_covered = num_of_coverd, temp_best_perm_DS[i][2]
				best_current_perm = temp_best_perm_DS[i]
				i+=1
			else:
				i+=1
		if(best_current_perm):
			res.append(best_current_perm)
			for gene in best_current_perm[3]:
				if gene in uncovered_genes: #there is a probability that this gene had already been covered bya prevuis sgRNA
					uncovered_genes.remove(gene)
	return res

def return_upgma(seq_list, names_list, df, cfd_dict = None):
	'''input:  a list of names and a list of sequences, calibrated
	output: an upgma instance.
	'''
	if df == Metric.cfd_funct:
		#base = seq_list
		#metric_seq_list = list(map(lambda t: Metric.pos_in_metric_general(t,df,base, cfd_dict), seq_list))
		
        #create a list of target vectors
		seq_list = list(map(lambda t: Metric.pos_in_metric_cfd_np(t, cfd_dict), seq_list)) #to uncomment

	#	df = Metric.find_dist_t  #if prev line is not  is use #to uncomment
		df = Metric.find_dist_np
    #create the distance matrix
	matrix = UPGMA.make_initiale_matrix(df, seq_list)
	m2 = UPGMA.make_distance_matrix(names_list, matrix)  #shuold be m2 = UPGMA.make_distance_matrix(names_list, matrix)
    #apply UPGMA, return a target tree
	upgma1 = UPGMA.make_UPGMA(m2)
	return upgma1

def find_distance_from_leaf_naive(node):
	if node.is_terminal():
		return 0
	else:
		return node.clades[0].branch_length + find_distance_from_leaf_naive(node.clades[0])

def fill_distance_from_leaves(tree):
	'''dinamic programing'''
	for leaf in tree.leaves_DS:  ##leaves is a python array
		leaf.set_distance_from_leaf(0)
	for leaf in tree.leaves_DS:
		node = leaf.parent
		while node and not (node.distance_from_leaf):
			node.distance_from_leaf = node.clades[0].distance_from_leaf + node.clades[0].branch_length
			node = node.parent

def fill_leaves_sets_Genes_tree_as_well(tree, sg_genes_dict, genes_tree = False):
	'''can be combine with fill_distance_from_leaves_function'''
	##fill the first line of nodes
	for leaf in tree.leaves_DS: ##node_targets_DS is a python array
		leaf.add_node_target(leaf.name)
		if not(genes_tree):
			leaf.set_candidates_DS()  #sg_genes_dict[leaf.name] is a list of genes which this target site is on
			current_candidate = Candidate.Candidate(leaf.name)
			current_candidate.fill_default_fildes(sg_genes_dict[leaf.name])
			leaf.candidates_DS[leaf.name] = current_candidate
		else:
			#'node_targets_DS' will be used to hold the genes; it is set to an empty list when the node is cunstracted. Maybe if this algorithm will be really bottems up, it will changed.
			leaf.add_node_target[leaf.name]

		node = leaf
		while(node.parent):
			for leaf in node.node_targets_DS:
				if leaf not in node.parent.node_targets_DS:
					node.parent.add_node_target(leaf)
			node = node.parent

def fill_leaves_sets(tree, sg_genes_dict):
	'''this version is not competable to genes tree.
	can be combine with fill_distance_from_leaves_function'''
	##fill the first line of nodes
	for leaf in tree.leaves_DS: ##node_targets_DS is a python array
		leaf.add_node_target(leaf.name)
		current_candidate = Candidate.Candidate(leaf.name)
		current_candidate.fill_default_fildes(sg_genes_dict[leaf.name])
		leaf.set_candidates_DS()  #sg_genes_dict[leaf.name] is a list of genes which this target site is on
		leaf.candidates_DS[leaf.name] = current_candidate
		node = leaf
		while(node.parent):
			for leaf in node.node_targets_DS:
				if leaf not in node.parent.node_targets_DS:
					node.parent.add_node_target(leaf)
			node = node.parent

def two_sequs_differeces_int(seq1,seq2):
	'''return a list of where the two sequences are different'''
	differences = 0  ##key: place of disagreement. value: the suggestions of each side
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	for i in range(1,len(seq2) - len(seq1)):
		differences |= math.pow(2,len(seq2) - i)
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			differences |= math.pow(2,i)
	return differences


def two_sequs_differeces_set(seq1,seq2):
	'''return a list of where the two sequences are different'''
	differences = set()  ##key: place of disagreement. value: the suggestions of each side
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	for i in range(1,len(seq2) - len(seq1)):
		differences.add(len(seq2) - i)
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			differences.add(i)
	return differences

def fill_PS_node(node):
	#find diferences between the representors
	P_S_set = set()
	if node.clades and len(node.clades)>1:
		P_S_set = two_sequs_differeces_set(node.clades[0].node_targets_DS[0], (node.clades[1].node_targets_DS[0]))
			#update the rest of the sites
		for clade in node.clades:
			P_S_set.update(clade.polymorphic_sites_set)
	node.set_polymorphic_sites_set(P_S_set)

def fill_PS(node):
	'''
	:param node: tree's root
	:return:
	'''
	if not node:
		return
	if node.clades and len(node.clades)>1:
		fill_PS(node.clades[0])
		fill_PS(node.clades[1])
	fill_PS_node(node)

def fill_leaves_sets_genes_tree(tree):
	'''can be combine with fill_distance_from_leaves_function'''
	##fill the first line of nodes
	for leaf in tree.leaves_DS: ##leaves_DS is a python array
		leaf.add_node_target(leaf)
		node = leaf
		while(node.parent):
			for leaf in node.node_targets_DS:
				if leaf not in node.parent.node_targets_DS:
					node.parent.add_node_target(leaf)
			node = node.parent


def fill_sg_genes_dict(input_sg_genes_dict):
	global sg_genes_dict
	sg_genes_dict = input_sg_genes_dict

def test_call_il_all():
	#sg_genes_dict = {"acgtacgt": ["gene1", "gene2"], "acgtagct": ["gene1"], "acgtattg":["gene2"], "atgcacgt": ["gene2", "gene3"], "atgcatgc":["gene3"]}  # a more complex one
	#genes_sg_dict = {"gene1":["acgtacgt","acgtagct"], "gene2":["acgtacgt", "acgtattg", "atgcacgt"], "gene3":["atgcacgt", "atgcatgc"]}
	genes_sg_dict = {"gene1": ["acgtacgt"], "gene2": ["acgtagct"], "gene3": ["acgtattg"], "gene4": ["atgcacgt"], "gene5": ["atgcatgc"]}
	sg_genes_dict = {"acgtacgt" : ["gene1"], "acgtagct" : ["gene2"] , "acgtattg":["gene3"] , "atgcacgt": ["gene4"] , "atgcatgc": ["gene5"]}
	sgNames = ["acgtacgt", "acgtagct", "acgtattg", "atgcacgt", "atgcatgc"]
	sgList = ["acgtacgt", "acgtagct", "acgtattg", "atgcacgt", "atgcatgc"]
	print("using upgma result", call_it_all(sgList, sgNames, sg_genes_dict, 0.11))
	print("naive algo result", Stage3.find_Uno_sgRNA(genes_sg_dict, 0.11))

def test_call_il_all_length20():
	genes_sg_dict = {"gene1": ["acgtacgtgtacgtacgtgt"], "gene2": ["acgtagctctacgtagctct"], "gene3": ["acgtattgtgacgtattgtg"], "gene4": ["atgcacgtgtatgcacgtgt"], "gene5": ["atatgcatgcatgcatgc"]}
	sg_genes_dict = {"acgtacgtgtacgtacgtgt" : ["gene1"], "acgtagctctacgtagctct" : ["gene2"] , "acgtattgtgacgtattgtg":["gene3"] , "atgcacgtgtatgcacgtgt": ["gene4"] , "atatgcatgcatgcatgc": ["gene5"]}
	sgNames = ["acgtacgtgtacgtacgtgt", "acgtagctctacgtagctct", "acgtattgtgacgtattgtg", "atgcacgtgtatgcacgtgt", "atatgcatgcatgcatgc"]
	sgList = ["acgtacgtgtacgtacgtgt", "acgtagctctacgtagctct", "acgtattgtgacgtattgtg", "atgcacgtgtatgcacgtgt", "atatgcatgcatgcatgc"]
	print("using upgma result", call_it_all(sgList, sgNames, sg_genes_dict))
	print("naive algo result", Naive.find_Uno_sgRNA(genes_sg_dict, Omega))



def find_best_sg_for_single_gene(gene_name,sg_lst):
	'''
	:param current_genes_sg_dict: a dictionary with only on key
	:return: current_best_perm, lowest_widest_prob. current_best_perm is of the form: (max_seq, fraction genes being cut among all the genes, probability to cut all the genes in genes list, genes_list, match_sites_list]), lowest_widest_prob
	'''
	return [Candidate.Candidate(sg_lst[0], 1, {gene_name:1}, {gene_name:[[sg_lst[0],{}]]})]


def find_best_sg_for_single_gene_naiveMC_returns_single(gene_name,sg_lst):
	''' the older version, sutable for when naive didn't make set cover
	:param current_genes_sg_dict: a dictionary with only on key
	:return: current_best_perm, lowest_widest_prob. current_best_perm is of the form: (max_seq, fraction genes being cut among all the genes, probability to cut all the genes in genes list, genes_list, match_sites_list]), lowest_widest_prob
	'''
	#
	return Candidate.Candidate(sg_lst[0], 1, {gene_name:1}, {gene_name:[]})

	#return [sg_lst[0],1,1,[gene_name],[]], 1  ##to change to something more sophisticated and maybe more accurate

###############from the old algorithm###########

def two_sequs_differeces(seq1,seq2):
	'''return a list of where the two sequences are different'''
	differences = {}  ##key: place of disagreement. value: the suggestions of each side
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	for i in range(1,len(seq2) - len(seq1)):
		differences[len(seq2) - i] = seq2[len(seq2) -i]
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			differences[i] = set([seq1[i], seq2[i]])
	return differences

def wheres_the_differences_not_naive(leave_DS):
	'''return a list of places which at least two sequences are different at'''
	checked_group = set() ##all the memmbers at this group have already been checked reagarding each other
	differences_in_checked_group = []
	for i in range(len(leave_DS)): ##node_targets_DS is a python array
		differences_from_checked_group = []
		not_to_check = set()
		##seen_cheaked = False  ##if the sequence had alrady tested against a memmber from the cheacked group
		for j in range(i, len(leave_DS)):
			if j in checked_group:
				current_differences = two_sequs_differeces(leave_DS[i], leave_DS[j], not_to_check)
				differences_from_checked_group += current_differences
				for i in current_differences:
					not_to_check.add(i)


def wheres_the_differences(leave_DS):
	''' return a dict: key is a place in which at least two sequences are different at, and value is a set of each letter that was in a seq in this place '''
	differences = {}  #key: place where there is a difference. value: letter apeared in this place
	for i in range(len(leave_DS)): ##node_targets_DS is a python array
		for j in range(i, len(leave_DS)):
			current_differences = two_sequs_differeces(leave_DS[i], leave_DS[j])
			for t in current_differences:
				if t in differences:
					differences[t] = differences[t] | current_differences[t]
				else:
					differences[t] = current_differences[t]
	return differences

def test_fill_distance():
	a = "aret"
	b = "ardw"
	c = "brdw"
	seq_list = [a,b,c]
	names = ["a", "b", "c"]
	matrix = UPGMA.make_initiale_matrix(UPGMA.p_distance,seq_list)
	m2 = UPGMA.make_distance_matrix(names, matrix)
	print("names")
	print(m2.names)
	m3 = m2.__repr__()
	upgma1 = UPGMA.make_UPGMA(m2)
	fill_leaves_sets(upgma1)
	fill_distance_from_leaves(upgma1)
	node = list(upgma1.root.leaves_DS)[0]
	while(node):
		node = node.parent
