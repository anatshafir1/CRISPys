import subgroup_res


def remove_repetitions_in_targets_sites(candidates_lst, alg, use_thr, Omega):
	'''updated version 16_2
	keep only with at least 1 unique targets'''

	def sum_num_of_mm_in_lowest_mm_target_in_relevant_genes(candidate, Omega, genes_list = None):
		if not genes_list:
			genes_list = [gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] >= Omega]
		res = 0
		for gene in genes_list:
			targets_lst = candidate.targets_dict.get(gene) #  targets_lst: [[target, sg_site_differents]]
			res += min(list(map(lambda sub_lst: len(sub_lst[1]), targets_lst)))
		#print('sum number of mm', res)
		return res

	def cover_above_thr(candidate, Omega, genes_list = None):
		cover_prob, num_of_covered = 1, 0
		if not genes_list:
			genes_list = [gene for gene in candidate.genes_score_dict if candidate.genes_score_dict[gene] >= Omega]
		for gene in genes_list:
			if candidate.genes_score_dict.get(gene, 0) >= Omega:
				cover_prob *= candidate.genes_score_dict[gene]
				num_of_covered += 1
		return cover_prob, num_of_covered, genes_list


	def contained(candidate_i, candidate_j, use_thr, Omega):
		'''dose candidate_i contained in candidate_j'''
		if use_thr:
			cover_prob_i, num_of_covered_i, genes_list = cover_above_thr(candidate_i, Omega)
		##first, check if j has all of i targets
			for key, value in candidate_i.targets_dict.items():
				if candidate_i.genes_score_dict.get(key,0) >= Omega and candidate_j.genes_score_dict.get(key,0) < Omega:# Omega:
					return False #for testing. shuold be False

			cover_prob_j, num_of_covered_j, _ = cover_above_thr(candidate_j, Omega)
			if num_of_covered_j > num_of_covered_i:
				if candidate_i.seq == 'GGCGAATGAGGAGGCTGAAG' and candidate_j.seq == 'GGCGAATGAGGAGGTTAAAG':
					print(num_of_covered_j, num_of_covered_i)
				return True
			if num_of_covered_j == num_of_covered_i:
				if cover_prob_j > cover_prob_i:
					return True
				if cover_prob_j == cover_prob_i:
					if sum_num_of_mm_in_lowest_mm_target_in_relevant_genes(candidate_j, Omega, genes_list) > sum_num_of_mm_in_lowest_mm_target_in_relevant_genes(candidate_i, Omega, genes_list) : #total num of mm over best site in relevant genes
						return True

			return False #for testing. should be False
		else:
			for key, value in candidate_i.targets_dict.items():
				if candidate_j.genes_score_dict.get(key,0) < 0.005: #epsilon
				#if not key in candidate_j.genes_score_dict:
					return False
			return candidate_i.cut_expectation < candidate_j.cut_expectation

	def remove_rep_subgroup(candidates_lst, use_thr, Omega):
		res = list()
		for i in range(len(candidates_lst)):
			to_append = 1
			for j in range(len(candidates_lst)):
				if (i != j and contained(candidates_lst[i], candidates_lst[j], use_thr, Omega)):
					to_append = 0
					#print('not appended')
			if to_append:
				#candidates_lst[i].off_targets = True
				res.append(candidates_lst[i])
		return res

	if alg == 'E':
		res = list()
		for i in range(len(candidates_lst)):
			subgroup_candidates_lst = remove_rep_subgroup(candidates_lst[i].candidate_lst,use_thr, Omega)
			subgroup_candidates_lst[0].off_targets = True
			res.append(subgroup_res.Subgroup_res(candidates_lst[i].genes_lst, subgroup_candidates_lst, candidates_lst[i].name))
		return res
	else:
		res = remove_rep_subgroup(candidates_lst, use_thr, Omega)
		for i in range(min(5, len(res))):
			res[i].off_targets = True
		return res