__author__ = 'Gal Hyams'
import sys
import CasSites
import Stage1
import Stage1_h
import UPGMA
import timeit
import pickle
import Metric

def call_using_CasSites_server_new(fasta_file, path , alg = 'A', where_in_gene = 1, use_thr = 0,  Omega = 1, df_targets = Metric.cfd_funct, protodist_outfile = "outfile", min_length= 20, max_length = 20,start_with_G = False, internal_node_candidates = 10, PS_number = 12):
	start = timeit.default_timer()
	cfd_dict = None
	if isinstance(where_in_gene, str):
		where_in_gene = float(where_in_gene.strip())
	if isinstance(Omega, str):
		Omega = float(Omega.strip())
	if isinstance(use_thr, str):
		use_thr = int(use_thr.strip())
	#df_targets = UPGMA.MITScore
	#df_targets = UPGMA.cfd_func
	#choosing the distance function
	if df_targets == "MITScore":
		df_targets = UPGMA.MITScore
	if df_targets == "cfd_funct" or df_targets == Metric.cfd_funct:
		df_targets = Metric.cfd_funct
		cfd_dict = pickle.load(open("cfd_dict.p",'rb'))
	if df_targets == "ccTop":
		df_targets = UPGMA.ccTop
	#df_targets = UPGMA.p_distance
	protodist_outfile = path + "/" + protodist_outfile
	original_range_in_gene = [0, where_in_gene]
	genes_sg_dict = {}
	sg_genes_dict = {}
	genesNames = []
	genesList = []
	f = open(fasta_file,'r')
	gene_name = ""
	gene_seq = ""
	lines = f.readlines()
	i = 0
	genes_exons_dict = {}  #key: gene name. value: list of exons
	while i <= len(lines):
	#stage 1: make  gene: sequence dictionary
		if i == len(lines) or lines[i][0] == '>':
			if len(gene_seq) > 0 and gene_name != "": #add the gene
				if gene_name not in genes_exons_dict:
					genes_exons_dict[gene_name] = [gene_seq]
				else:
					genes_exons_dict[gene_name] = genes_exons_dict[gene_name] + [gene_seq]
				gene_seq = ""
			if i != len(lines): # lines[i-1][0] == '>':
				gene_name = lines[i][1:].strip() #without the '>' and the '\n'
		elif lines[i] != "\n":
			gene_seq += lines[i].strip()
		i+=1
	#stage 2: find the target sites
	for gene_name in genes_exons_dict.keys():
		genes_sg_dict[gene_name] = CasSites.get_targets_sites_from_exons_lst(genes_exons_dict[gene_name], original_range_in_gene, min_length, max_length,start_with_G)
		genesNames.append(gene_name)
		genesList.append("".join(genes_exons_dict[gene_name]))
		#filling up the sg_genes_dict
		for sg in genes_sg_dict[gene_name]:
			if sg in sg_genes_dict:
				sg_genes_dict[sg] = sg_genes_dict[sg] + [gene_name]
			else:
				sg_genes_dict[sg] = [gene_name]
	if alg == 'E':
		res = Stage1_h.call_it_all(genesList, genesNames, sg_genes_dict, genes_sg_dict, Omega, protodist_outfile, path, df_targets, internal_node_candidates, cfd_dict, PS_number)
	else:
		res = Stage1.call_it_all(genesList, genesNames, sg_genes_dict, genes_sg_dict, Omega, protodist_outfile, path, df_targets, cfd_dict, PS_number) #thies line have been change to be sutable for wrapper
	if use_thr:
		sort_thr(res, Omega, alg == 'E')
	else:
		sort_expectation(res, alg == 'E')
	#remove the folowing two lines when using CRISPysCover
	if len(res)>200:
		res = res[:200]
	Stage1.print_res_to_csvV2(res, sg_genes_dict, genesList, genesNames, path, alg == 'E')
	Stage1.print_res_to_csvV3(res, sg_genes_dict, genesList, genesNames, path, alg =='E')
	pickle.dump(res, open(path + "/res_in_lst.p", "wb"))
	pickle.dump(genesNames, open(path + "/genesNames.p", "wb"))
	stop = timeit.default_timer()
	print("time: ", stop - start)
	time_file = open("time.txt", 'w')
	time_file.write(str(stop - start))
	time_file.close
	return res

def sort_expectation(candidates_DS, homology):
	def sort_subgroup(candidates_DS):
		candidates_DS.sort(key = lambda item: (item.cut_expectation, item.total_num_of_mm()), reverse=True)
	if not homology:
		sort_subgroup(candidates_DS)
	else:
		for i in range(len(candidates_DS)):
			sort_subgroup(candidates_DS[i].candidate_lst)

	
	
def sort_thr(candidates_DS, Omega, homology):
	'''dort the candidates DS by num of genes with cut prob> Omega and then by the probobility to cleave all of these genes'''
	def sort_subgroup(candidates_DS, Omega):
		for candidate in candidates_DS:
			num_of_genes_above_thr = 0
			cleave_all = 1
			for gene, score in candidate.genes_score_dict.items():
				if score >= Omega:
					cleave_all *= score
					num_of_genes_above_thr += 1
			candidate.cleve_all_above_thr = cleave_all
			candidate.num_of_genes_above_thr = num_of_genes_above_thr
		candidates_DS.sort(key = lambda item: (item.num_of_genes_above_thr, item.cleave_all_above_thr), reverse = True)
	if not homology:
		sort_subgroup(candidates_DS, Omega)
	else:
		for i in range(len(candidates_DS)):
			sort_subgroup(candidates_DS[i], Omega)


def leave_only_relevant_sgRNA(res):
	if len(res) < 1:
		return
	candidates_to_del = []
	for i in range(len(res) -1, -1,-1):
		if res[i].cut_expectation < 1:
			del res[i]
		elif i < len(res) - 1:
			for j in range(i+1,len(res)):
				if j >= len(res):
					continue
				if res[i].seq == res[j].seq: # there is no need in both of them. Who is sutable for more genes?
					if res[i].cut_expectation <= res[j].cut_expectation:
						del res[i]
					else:
						del res[j]


def find_unrepresented_genes(genesNames, res):
	'''res format: list of list. each sublist: [perm sequence, int, int, list of genes, list of list: differences for each gene'''
	unrepresented = set(genesNames)  ##all the genes in the family
	represented = set()
	for perm in res:
		for gene in perm.genes_score_dict.keys():
			represented.add(gene)
	##now the represented set contains all the represented genes
	unrepresented = unrepresented - represented
	return list(unrepresented)


def add_unrepresented(genes_not_represented, genesNames, genesList):
	if len(genes_not_represented) == 1:  #find the first
		index_in_lst = genesNames.index(genes_not_represented[0])
		gene_seq = genesList[index_in_lst]

def remove_repetitions_in_targets_sites(res):
	'''haven't been tested yet'''
	targets = set()
	for i in range(len(res)):
		removed_flag = 0
		for target_tuple in candidate.targets_dict.values():
			if target_tuple[0] in targets:
				to_remove.append[i]
				removed_flag = 1
		if removed_flag == 0:
			for target_tuple in candidate.targets_dict.values():
				targets.add(target_tuple[0])
	#now, remove
	for index in range(len(to_remove) -1, -1,-1):
		del res[to_remove[index]]

	
			
			


if __name__ == "__main__":

	call_using_CasSites_server_new(*sys.argv[1:])

