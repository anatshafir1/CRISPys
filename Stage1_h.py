
__author__ = 'GH'

##using UPGMA insted of NJ

import Stage1
import Stage2
import copy
from subgroup_res import Subgroup_res

def print_res_to_file_E(res, input_sg_genes_dict, path=''):
	''' input(res) format:	array of permutations_DS'''
	##first, get genes set, for the file name:

	genes_set = set()
	for sg, genes in input_sg_genes_dict.items():  #genes  is a list of genes
		for gene in genes:
			genes_set.add(gene.split("RNA")[0])
	if path != '':
		file_name = path+ "\output.txt"
	else:
		file_name = str(genes_set)[1:-1] + ".txt"
	print("output path = ", file_name)
	f = open(file_name, "w")
	for permDS in res:
		f.write("SG:\n")
		f.write(permDS[0])
		f.write("\n")
		f.write("score:\n")
		f.write(str(permDS[2]))
		f.write("\n")
		f.write("Cleaved genes:\n")
		for gene_tup in permDS[3]:
			f.write(gene_tup[0] + "  ")
		f.write("\n")
		f.write("cleavage site for each gene:\n")
		next_line = str(permDS[4])
		next_line = next_line.replace("RNAfile","")
		f.write(next_line)
		f.write("\n\n")
	f.close()

def fill_sg_genes_dict(input_sg_genes_dict):
	global sg_genes_dict
	sg_genes_dict = input_sg_genes_dict

def fill_genes_sg_dict(input_genes_sg_dict):
	global genes_sg_dict
	genes_sg_dict = input_genes_sg_dict


#############################################################################################################


def call_it_all(genesList, genesNames, input_sg_genes_dict, input_genes_sg_dict, Omega, protdist_outfile, pylip_temps_path, df_targets, internal_node_candidates, cfd_dict = None, PS_number = 12):
	# make a tree and distance matrix of the genes
	upgmaTree, distance_matrix = Stage1.return_UPGMA(genesList, genesNames, protdist_outfile, pylip_temps_path) #to uncomment when using wighted
	write_newik_to_file(upgmaTree.root, pylip_temps_path)
	tree_to_file(upgmaTree.root, pylip_temps_path)
	Stage2.fill_leaves_sets_genes_tree(upgmaTree)  # as apposed to the intermediate algorithem, here leaves are genes
	fill_sg_genes_dict(input_sg_genes_dict)
	fill_genes_sg_dict(input_genes_sg_dict)
	#making the sgList for Algorithm B:
	sgList = list(input_sg_genes_dict.keys())
	#print('sg list',sgList)
	sgNames = copy.deepcopy(sgList)
	res = []
	#E_bottem_up_tree(res, upgmaTree, Omega, df_targets)
	E_top_down(res, upgmaTree.root, Omega, sg_genes_dict, df_targets, internal_node_candidates, cfd_dict, PS_number)

	return res


def write_newik_to_file(node, path):
	f = open(path + "/tree.newick", 'w')
	write_newik(node, f, 0)
	f.write(';')
	f.close()

def write_newik(node, f, index):
	index += 1
	if node.is_terminal():
		#print(node.name ,f)
		f.write(node.name)
		#print(node.name)
		#print("node clades:", node.clades)
	else:
		f.write('(')
	#print('(', f)
	#print('(')
		for i in range(len(node.clades)):
			write_newik(node.clades[i], f, index)
		#print(',', f)
			if i < len(node.clades) -1:# and node.clades[i].is_terminal() and node.clades[i+1].is_terminal():
				f.write(',')
		#print(',')
		f.write(')n' + str(index))
	#print(')', f)
	#print(')')

def print_tree_to_file(tree, path):
	'''still need to be tested'''
	print(tree, file=open(path + "GenesTree.txt", 'w'))
	return


def stopping_condition(current_best_perm):
	'''rapper version: return True if need to stop'''
	if not current_best_perm:
		return False
	return current_best_perm[1] < 1  #if at least 100% of the genes in the subtree are cut, continue


def E_bottem_up(res, node, current_sg_genes_dict, current_genes_sg_dict, sgList, sgNames, Omega, df_targets, internal_node_candidates = 10):
	'''on the genes tree. Caling the buttoms up algorithem with a sg genes dict sutable for the subtree
	#print(node.leaves_DS)
	if node.colour == 'b':
		return
	##making the genes_sg dict for this subtree and the sg_genes_dict to send to the intermadiate algorithm
	if not(current_sg_genes_dict):
		current_sg_genes_dict = {}
	if not(current_genes_sg_dict):
		current_genes_sg_dict = {}
	if not(sgList):
		sgList = []
	if not(sgNames):
		sgNames = []
	for leaf in node.node_targets_DS: ##leaf here is a gene. taking only the relevant genes
		current_genes_sg_dict[leaf.name] = genes_sg_dict[leaf.name]
		##filling the sg genes dict
		for sg in current_genes_sg_dict[leaf.name]:
			current_sg_genes_dict[sg] = sg_genes_dict[sg] ##the checking if this sg is already in the dict just be more expensive overall
			if sg not in sgList:
				sgList.append(sg)
				sgNames.append(sg)

	if len(node.node_targets_DS) > 10:
		#print(node.node_targets_DS)
		return
	if len(current_genes_sg_dict) > 1 :  #more then one gene

		#current_set_cover = bottemsUpAlgorithm.find_best_sg_for_single_gene(leaf.name, sgList )#lowest_of_widest is not in use in this function
	#else:
		#get the set cover from the bottem up algorithm
		best_permutations_DS = bottemsUpAlgorithm.call_it_all(sgList, sgNames, current_sg_genes_dict, Omega, df_targets)##call_it_all(sgList, sgNames, input_sg_genes_dict, Omega)## Naive.find_Uno_sgRNA(current
		#print("current set cover= ",current_set_cover)
		#if current_set_cover == None:
		#	return
		if not (best_permutations_DS):
			return

		best_permutations_DS.sort(key = lambda item: item.cut_expectation, reverse=True)
		#call_MULTICRISPR_Wrapper.sort_thr(best_permutations_DS, Omega)
		current_best_perm = best_permutations_DS[:internal_node_candidates]  #the best sg at the current set cover
		#print(current_best_perm)
	#print(current_best_perm)
	#if current_best_perm[1] == 1: #all the genes were cought
		#res.append(current_best_perm)
		#res += current_best_perm
		res.append(Subgroup_res(get_genes_list(best_permutations_DS) ,current_best_perm, node.name))
	node.set_colour('b')

	if (node.parent):
		E_bottem_up(res, node.parent, current_sg_genes_dict, current_genes_sg_dict, sgList, sgNames, Omega, df_targets, internal_node_candidates)  ##this line is adopted to the rapper algorithm
	'''



def E_bottem_up_tree(res, upgmaTree, Omega, df_targets):
	for i in range(len(upgmaTree.leaves_DS)):
		E_bottem_up(res, upgmaTree.leaves_DS[i], None,None,None,None, Omega, df_targets)

def tree_to_file(node, path):
	lst = list()
	tree_prerder(node,lst)
	print(lst, file=open(path + "/GenesTree.txt", 'w'))

def tree_prerder(node, lst):
	if not node:
		return
	#print(node.name, file=open(path + "GenesTree.txt", 'w'))
	lst.append(node.name)
	if node.clades:
		tree_prerder(node.clades[0], lst)
		if len(node.clades)>1:
			tree_prerder(node.clades[1], lst)

########################################################################################E top down######################################3
def E_top_down(res, node, Omega, sg_genes_dict, df_targets, internal_node_candidates = 10, cfd_dict = None, PS_number = 12): #(res, node, current_sg_genes_dict, current_genes_sg_dict, sgList, sgNames, Omega, df_targets)
	'''
	:param node:
	:param current_genes_sg_dict:
	:param Omega: can be removed already??
	:return:
	'''
	#if len(node.polymorphic_sites_set) < 11: #change to 12!
			##making the genes_sg dict for this subtree and the sg_genes_dict to send to the intermadiate algorithm
	current_sg_genes_dict = dict()
	current_genes_sg_dict = dict()
	sgList = list()
	sgNames = list()
	for leaf in node.node_targets_DS: ##leaf here is a gene. taking only the relevant genes
		current_genes_sg_dict[leaf.name] = genes_sg_dict[leaf.name]
		##filling the sg genes dict
		for sg in current_genes_sg_dict[leaf.name]:
			#current_sg_genes_dict[sg] = sg_genes_dict[sg] ###here is the abnormality!! ##the checking if this sg is already in the dict just be more expensive overall
			
			#will the folowing be clearer?
			#untab the folowing 3 lines:
			if sg in current_sg_genes_dict:
				current_sg_genes_dict[sg] = current_sg_genes_dict[sg] +  [leaf.name]
			else:
				current_sg_genes_dict[sg] = [leaf.name] ###here is the abnormality!! ##the checking if this sg is already in the dict just be more expensive overall
			
			if sg not in sgList:
				sgList.append(sg)
				sgNames.append(sg)		

	#if len(current_genes_sg_dict) > 1 :  #more then one gene
	#print("num of sg: ", len(current_sg_genes_dict))
	if len(node.node_targets_DS) < 9 and len(node.node_targets_DS) > 1: #should change the first parameter to at least 10
		#current_set_cover = bottemsUpAlgorithm.find_best_sg_for_single_gene(leaf.name, sgList )#lowest_of_widest is not in use in this function
	#else:
		#get the set cover from the bottem up algorithm
		#print("genes: ")
		#print(set(list(itertools.chain.from_iterable(current_sg_genes_dict.values()))))
		best_permutations_DS = Stage2.call_it_all(sgList, sgNames, current_sg_genes_dict, Omega, df_targets, cfd_dict, PS_number)##call_it_all(sgList, sgNames, input_sg_genes_dict, Omega)## Naive.find_Uno_sgRNA(current
		#print("current set cover= ",current_set_cover)
		#if current_set_cover == None:
		#	return
		if not (best_permutations_DS):
			return

		best_permutations_DS.sort(key = lambda item: item.cut_expectation, reverse=True)
		#call_MULTICRISPR_Wrapper.sort_thr(best_permutations_DS, Omega)
		current_best_perm = best_permutations_DS[:internal_node_candidates]  #the best sg at the current set cover
		#print(current_best_perm)

	#print(current_best_perm)
	#if current_best_perm[1] == 1: #all the genes were cought
		#res.append(current_best_perm)
		#res += current_best_perm
		#res += [current_best_perm]
		res.append(Subgroup_res(get_genes_list(best_permutations_DS) ,current_best_perm, node.name))
		# print(res[0].genes_lst)

	#else:

	if not node.clades:
		return

	if node.clades[0]:
		E_top_down(res, node.clades[0], Omega, sg_genes_dict, df_targets, internal_node_candidates, cfd_dict, PS_number)
	if node.clades[1]:
		E_top_down(res, node.clades[1], Omega, sg_genes_dict, df_targets, internal_node_candidates, cfd_dict, PS_number)

def get_genes_list(candidates_lst):
	res = set()
	for c in candidates_lst:
		for gene in c.genes_score_dict.keys():
			res.add(gene)
	return list(res)
