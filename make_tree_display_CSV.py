import pickle
import Candidate
import subgroup_res
import CasSites


def sub_tree_display(path, candidates_lst, f, consider_homology, counter, genes_names, genes_lst):

	def chagne_mm_to_lowercase(target_str, mm_lst):
		'''mm is mismatch'''
		target_in_lst = list(target_str)
		for place in mm_lst:
			target_in_lst[place] = target_in_lst[place].lower()
		return ''.join(target_in_lst)

	header_row = "sgRNA index,sgRNA,Score,Genes,Genes score,Target site,#mms,Position\n"#new

	if consider_homology == True:
		#make a set of the genes in the group to write before each table. added by Udi 13/04/22
		genes = []
		for candidate in candidates_lst:
			genes += [g for g in candidate.genes_score_dict.keys()]
		genes = set(genes)
		f.write("Genes in group=,"+str(genes)+"\n")
		# f.write("table id="+str(counter)+"\n")

	f.write(header_row)
	sgRNA_index = 0
	for c in candidates_lst:
		sgRNA_index += 1

		num_of_targets = 0
		for targets in c.targets_dict.values():
			num_of_targets += len(targets)
	#    print("num_of_targets: ",num_of_targets)
		first_target = 1
		first_gene = 1
		l = list(c.targets_dict.items())
		l.sort(key = lambda item: c.genes_score_dict[item[0]], reverse = True)

		for gene, targets in l:
			targets.sort(key = lambda target: len(target[1]))#Galll!!! sort the targets by number of mms 
			gene_seq = genes_lst[genes_names.index(gene)]
			seen_sites = dict()
			first_target = 1
			for target in targets:
				pos = find_pos(target,gene_seq, seen_sites)

				if first_target == 1 and first_gene==1:

					f.write(str(sgRNA_index) + '.,' + c.seq+","+str(c.cut_expectation)[:5])
					f.write(","+gene)
					score = str(c.genes_score_dict[gene])
					if len(score)> 5:
						score = score[:5]

					f.write(","+score)
					f.write(","+chagne_mm_to_lowercase(target[0], target[1].keys()))#+"\n")
					f.write(","+str(len(target[1])))
					f.write(","+ pos)
					f.write("\n")
					first_target = 0
					continue
				if first_target != 1:
					f.write(str(sgRNA_index) + ".,,,,,")

					f.write(chagne_mm_to_lowercase(target[0], target[1].keys()))
					f.write(","+str(len(target[1])))
					f.write(","+pos)
					f.write("\n")
				if first_target == 1 and first_gene != 1:
					f.write(str(sgRNA_index) + ".,,,")
					score = str(c.genes_score_dict[gene])
					if len(score)> 5:
						score = score[:5]

					f.write(gene)
					f.write(","+score)
					f.write(","+chagne_mm_to_lowercase(target[0], target[1].keys()))
					f.write(","+str(len(target[1])))
					f.write(","+pos)
					f.write("\n")

					first_target = 0
			first_gene = 0

def genes_of_sub_candidates_lst(sub_candidates_lst):
	res = set()
	for c in sub_candidates_lst:
		for gene in c.genes_score_dict.keys():
			res.add(gene)
	return list(res)

def test_find_pos():
	target, gene_sequence = ['AAAAAAAAAAAAAAAAAAGG'], 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGG'
	seen_sites = dict()
	position = 0
	for i in range(99):
		#r = find_pos(target, gene_sequence, seen_sites)
		position = gene_sequence.find(target[0], position+1)
		if position == -1:
			break
	r = find_pos(target, gene_sequence, seen_sites)


	
def find_pos(target, gene_sequence, seen_sites):
	'''sgRNA_targets is a list of target sites
	returns'''
	#seen_sites = dict()
	#res = dict() #key:targets val: position in gene
	#for target in sgRNA_targets:
	target_seq = target[0]
	if target_seq in seen_sites:
		directions_lst = seen_sites[target_seq]
	else:
		directions_lst = [0,0]
	position = gene_sequence.find(target_seq, directions_lst[0])#+ '+' #item[0] is the target site seq
	if position != -1:
		update_seen_sites_dict(seen_sites, target_seq, 0, position)
		position = str(position) + '+'
	#update_positions_dict(seen_sites, target ,sg_position)
	else:
		position = CasSites.give_complementary(gene_sequence).find(target_seq, directions_lst[1])# + "-"
		update_seen_sites_dict(seen_sites, target_seq, 1, position)
		position = str(position) + '-'
	if position == -1:
		position = ''
		#update_positions_dict(res, target ,sg_position)
		#update_seen_sites_dict(seen_sites, target_seq, 1)
	return position

def update_positions_dict(d, site_seq, position):
	if site_seq in d:
		d[site_seq] = d[site_seq] + [position]
	else:
		d[site_seq] = [position]


def update_seen_sites_dict(d, site_seq, direction, position):
	'''
	d: dict: key: site_seq; val: [num of occurrences, directions_lst]
	direction: 0 for sense, 1 for antisence
	'''
	if site_seq in d:
		directions_lst = d[site_seq]
		directions_lst[direction] = position + 1
		d[site_seq] = directions_lst
	else:
		directions_lst = [0,0]
		directions_lst[direction] = position + 1
		d[site_seq] = directions_lst

		
def tree_display(path, consider_homology = False, set_cover = False):
	if set_cover:
		candidates_lst = pickle.load(open(path + "/greedy_cover.p", "rb"))
	else:
		candidates_lst = pickle.load(open(path + "/res_in_lst.p", "rb"))
	genes_names = pickle.load(open(path + "/genesNames.p", "rb"))
	genes_list = pickle.load(open(path + '/genesList.p', 'rb'))

	filepath = path + "/CRISPys_output.csv"
	f = open(filepath, 'w')
	counter = 0;
	if consider_homology == False:
		sub_tree_display(path, candidates_lst, f, consider_homology,counter, genes_names, genes_list)
	else:
		f.write("The designed sgRNAs for the genes in your input are listed in the table below. Every section of the table corresponds to a homologous genes subgroup as specified by the internal nodes of the constructed genes tree.<br>The name of the subgroup and the list of genes are given in the header of each section.\n")

		for subgroup_item in candidates_lst:
			# create the main table
			counter +=1
			# genes = str(subgroup_item.genes_lst)[1:-1].replace("'","")

			sub_tree_display(path, subgroup_item.candidate_lst, f, consider_homology, counter, genes_names, genes_list)
		counter = 0;	

	f.close()
#test

if __name__ == "__main__":
#	print("ASDF")
	#test_find_pos()
	tree_display("/groups/itay_mayrose/galhyams/1516893877")
