
__author__ = 'GH'

# the naive algorithm - for a given family, with its group of sgRNA, find the most promising sgRNAs
import Candidate
import Distance_matrix_and_UPGMA
import Metric
import random

random.seed(1234)


def find_candidates(genes_targets_dict, Omega, scoring_function, node, cfd_dict=None):
    """ not Uno any more....
	genes_targets_dict: keys are genes names, values are lists of sgRNA sequences suitable to this gene. classed = a list of
	the candidates DS of the children of the node """
    # stage one: make a list of all the sgRNAs
    list_of_sg = []
    for key, val in genes_targets_dict.items():
        list_of_sg += val
    # stage two: find the suitable sgRNA:
    for_single_gene = False
    return return_candidates(list_of_sg, list_of_sg[0], genes_targets_dict, Omega, scoring_function, node, for_single_gene, cfd_dict)


def generate_scores(genes_targets_dict, list_of_candidates, scoring_function, cfd_dict=None):  # Omer caldararu 24/3
    """
	generates a data structure that contains the candidates and their off-target scores.
	(in the case of gold off, or any other function that can accept several sgRNA's in a single call)
	Args:
		genes_targets_dict: a dictionary : gene name -> targets within the gene
		list_of_candidates: a list of all possible candidates, given by all_perms()
		scoring_function: the scoring function
		cfd_dict: cfd dictionary used for the cfd function

	Returns: scores_dict = {gene : [(target,candidates_target_scores) for target in the gene]}
	"""
    scores_dict = {}
    if scoring_function == Distance_matrix_and_UPGMA.gold_off_func:
        return generate_scores_one_batch(genes_targets_dict, list_of_candidates, scoring_function, scores_dict)
    for gene in genes_targets_dict.keys():
        scores_dict[gene] = []
        for target in genes_targets_dict[gene]:
            if scoring_function == Distance_matrix_and_UPGMA.ccTop or scoring_function == Distance_matrix_and_UPGMA.MITScore or scoring_function == Metric.cfd_funct:
                candidates_target_scores = list(
                    map(lambda sg: scoring_function(sg, target, cfd_dict), list_of_candidates))
            scores_dict[gene].append((target, candidates_target_scores))
    return scores_dict


def generate_scores_one_batch(genes_targets_dict, list_of_candidates, scoring_function, scores_dict):
    """
	generates a data structure that contains the candidates and their off-target scores,
	using a single call of the scoring function.
	Args:
		genes_targets_dict: a dictionary : gene name -> targets within the gene
		list_of_candidates: a list of all possible candidates, given by all_perms()
		scoring_function: the scoring function
		scores_dict: an empty dictionary
	Returns: scores_dict = {gene : [(target,candidates_target_scores) for target in the gene]},
	"""
    genes_list = list(genes_targets_dict.keys())
    batch_targets_list = []
    batch_candidates_list = []
    for gene in genes_list:
        for target in genes_targets_dict[gene]:
            batch_targets_list += [target] * len(list_of_candidates)
        batch_candidates_list += list_of_candidates * len(genes_targets_dict[gene])
    list_of_all_scores = scoring_function(batch_candidates_list, batch_targets_list)
    i = 0
    for gene in genes_list:
        scores_dict[gene] = []
        for target in genes_targets_dict[gene]:
            candidates_target_scores = list_of_all_scores[i:i + len(list_of_candidates)]
            scores_dict[gene].append((target, candidates_target_scores))
            i += len(list_of_candidates)
    return scores_dict


def return_candidates(list_of_targets, initial_seq, genes_sg_dict, Omega, df, node, for_single_gene=False,
                      cfd_dict=None):
    dict_of_different_places = wheres_the_differences_linear(
        list_of_targets)  # node_targets is a python array. where_the_differences.
    node.polymorphic_sites = dict_of_different_places
    list_of_different_places = list(dict_of_different_places.items())
    list_of_different_places.sort(key=lambda item: item[0])
    # going over all the permutations
    list_of_perms_seqs = all_perms(initial_seq, None, list_of_different_places)
    list_of_candidates = []  # a list of tuples: (candidate_str,fraction_of_cut, cut_expectation, genes_list)
    scores_dict = generate_scores(genes_sg_dict, list_of_perms_seqs, df, cfd_dict)
    for i in range(len(list_of_perms_seqs)):
        targets_dict = {}  # a list of tuples: (gene name, list of target of this gene that might be cut by the candidate_str)
        genes_covering = []  # a list of tuples: (gene name, probability to be cut).
        for gene in scores_dict.keys():
            prob_gene_will_not_cut = 1  # the probability that a gene will not be cut by the candidate
            list_of_targets = []  # for later knowing where the candidate_str might cut in each gene (when writing the output)
            num_of_cuts_per_gene = 0  # in use only in the single gene version

            for target, candidates_target_scores in scores_dict[gene]:
                if candidates_target_scores[i] == 1:  # in case the distance is 1 (it means that the score is 0 and
                    # there isn't attachment of the guide and target nad no cut event) don't consider the target
                    continue
                candidate_cut_prob = 1 - candidates_target_scores[i]
                sg_site_differences = two_seqs_differences(list_of_perms_seqs[i],
                                                           target)  # the differences between the ith candidate and the target

                list_of_targets.append([target[:20], sg_site_differences])
                prob_gene_will_not_cut = prob_gene_will_not_cut * (
                            1 - candidate_cut_prob)  # lowering the not cut prob in each sgRNA
                num_of_cuts_per_gene += candidate_cut_prob
            prob_gene_cut = 1 - prob_gene_will_not_cut
            if for_single_gene:  # is this necessary? omer 11/4
                genes_covering.append((gene, num_of_cuts_per_gene))
            if prob_gene_cut >= Omega:
                targets_dict[gene] = list_of_targets
                genes_covering.append((gene, prob_gene_cut))
        cut_expectation = 0  # the probability the permutated sequence will cut all the genes, that the probability each
        # of them will be cut is greater then omega
        genes_score_dict = {}  # a dict of genes: genes considered cut by this sequence, and cut prob
        for tuple in genes_covering:  # tuple : (gene name, probability to be cut)
            cut_expectation += tuple[1]  # the prob to cut all the genes
            genes_score_dict[tuple[0]] = tuple[1]
        if len(genes_covering) >= 1:  # If the candidate has at least one gene with a score above omega, add it to the result  omer 18/04
            current_candidate = Candidate.Candidate(list_of_perms_seqs[i], cut_expectation, genes_score_dict,
                                                    targets_dict)
            list_of_candidates.append(current_candidate)
    del list_of_perms_seqs
    return list_of_candidates


def all_perms(initial_seq, list_of_seqs, list_of_differences):
    """each recursive call add the next part to the sequences. the result is sequences off each of the parameters
	list of differences : list of tuples: (place, set of letters)"""
    if len(list_of_differences) == 0:  # the stopping condition
        if list_of_seqs:
            return list_of_seqs
        elif initial_seq:
            return [initial_seq[:20]]
        else:
            return []
    else:
        new_list_of_seqs = []
        if not list_of_seqs:  # initialising the list of sequences
            list_of_seqs = []
            list_of_seqs.append(initial_seq[:list_of_differences[0][0]])
        for seq in list_of_seqs:
            for letter in list_of_differences[0][1]:
                if len(list_of_differences) > 1:
                    new_list_of_seqs.append(seq + letter + initial_seq[len(seq) + 1:list_of_differences[1][0]])
                    # the place of the next versatile letter place
                else:
                    new_list_of_seqs.append(seq + letter + initial_seq[len(seq) + 1:20])
        del list_of_seqs
        return all_perms(initial_seq, new_list_of_seqs, list_of_differences[1:])


def two_seqs_differences(seq1, seq2):
    """return a list of where the two sequences are different"""
    differences = {}  # key: place of disagreement. value: the suggestions of each side
    seq1 = seq1[:20]
    seq2 = seq2[:20]  # in cases the PAM is not sliced - don't take PAM into account
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences[i] = [seq1[i], seq2[i]]
    return differences


def wheres_the_differences_linear(leaves):
    differences = dict()  # key: place where there is a difference. value: letter appeared in this place
    if len(leaves) < 2:
        return differences
    ref = leaves[0]
    for i in range(1, len(leaves)):  # node_targets is a python array
        current_differences = two_seqs_differences(ref, leaves[i])
        for t in current_differences:
            if t in differences:
                differences[t] = list(set(differences[t] + current_differences[t]))
            else:
                differences[t] = current_differences[t]
    return differences
