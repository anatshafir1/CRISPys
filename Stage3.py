"""Finding sgRNA candidates"""
__author__ = 'GH'

# the naive algorithm - for a given family, with its group of sgRNA, find the most promising sgRNAs
from typing import List, Dict, Tuple
from Candidate import Candidate
from Distance_matrix_and_UPGMA import gold_off_func, crisprnet, moff, ccTop, ucrispr, MITScore, default_on_target
from Metric import cfd_funct
from TreeConstruction_changed import CladeNew


def generate_scores(genes_targets_dict: Dict[str, List[str]], list_of_candidates: List[str], off_scoring_function,
                    on_scoring_function, cfd_dict=None) -> Dict[str, List[Tuple[str, List[float]]]]:
    """
	generates a data structure that contains the candidates and their off-target scores.
	(in the case of gold off, or any other function that can accept several sgRNA's in a single call)
	
	:param genes_targets_dict: a dictionary : gene name -> targets within the gene
	:param list_of_candidates: a list of all possible candidates, given by all_perms()
	:param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
	:param cfd_dict: cfd dictionary used for the cfd function
	:return: scores_dict = {gene : [(target,candidates_target_scores) for target in the gene]}
	"""
    scores_dict = {}
    if off_scoring_function in (gold_off_func, crisprnet, moff) \
            or getattr(off_scoring_function, "func", "") == ucrispr:
        return generate_scores_one_batch(genes_targets_dict, list_of_candidates, off_scoring_function,
                                         on_scoring_function)
    if on_scoring_function == default_on_target:  # if no on-target scoring function was given
        for gene in genes_targets_dict:
            scores_dict[gene] = []
            for target in genes_targets_dict[gene]:
                candidates_target_scores = []
                for sg in list_of_candidates:
                    if off_scoring_function == cfd_funct:
                        candidates_target_scores += [off_scoring_function(sg, target, cfd_dict)]
                    elif off_scoring_function == ccTop or off_scoring_function == MITScore:
                        candidates_target_scores += [off_scoring_function(sg, target[:20])]
                scores_dict[gene].append((target, candidates_target_scores))
        return scores_dict
    else:
        batch_candidates_list = []
        for gene in genes_targets_dict:
            for target in genes_targets_dict[gene]:
                candidates_with_pam = [f"{c[:20]}{target[20]}" for c in list_of_candidates]
                batch_candidates_list += candidates_with_pam
        on_scores = on_scoring_function(batch_candidates_list)
        i = 0
        for gene in genes_targets_dict:
            scores_dict[gene] = []
            for target in genes_targets_dict[gene]:
                candidates_target_scores = []
                for sg in list_of_candidates:
                    if off_scoring_function == cfd_funct:
                        candidates_target_scores += [off_scoring_function(sg, target, cfd_dict) * on_scores[i]]
                    elif off_scoring_function == ccTop or off_scoring_function == MITScore:
                        candidates_target_scores += [off_scoring_function(sg, target[:20]) * on_scores[i]]
                    i += 1
                scores_dict[gene].append((target, candidates_target_scores))
        return scores_dict


def generate_scores_one_batch(genes_targets_dict: Dict[str, List[str]], list_of_candidates: List[str],
                              off_scoring_function, on_scoring_function) -> Dict[str, List[Tuple[str, List[float]]]]:
    """
    generates a data structure that contains the candidates and their off-target scores,
    using a single call of the scoring function.
    :param genes_targets_dict: a dictionary : gene name -> targets within the gene
    :param list_of_candidates: a list of all possible candidates, given by all_perms()
    :param off_scoring_function: the off target scoring function
    :param on_scoring_function: the on target scoring function
    :return: scores_dict = {gene : [(target,candidates_target_scores) for target in the gene]},
    """
    scores_dict = {}
    batch_targets_list = []
    batch_candidates_list = []
    for gene in genes_targets_dict:
        for target in genes_targets_dict[gene]:
            batch_targets_list.extend([target] * len(list_of_candidates))
            candidates_with_pam = [f"{c[:20]}{target[20:]}" for c in list_of_candidates]
            batch_candidates_list += candidates_with_pam
    off_scores = off_scoring_function(batch_candidates_list, batch_targets_list, for_metric=False)
    if on_scoring_function == default_on_target:  # if no on-target scoring function was given
        i = 0
        for gene in genes_targets_dict:
            scores_dict[gene] = []
            for target in genes_targets_dict[gene]:
                candidates_target_scores = off_scores[i:i + len(list_of_candidates)]
                scores_dict[gene].append((target, candidates_target_scores))
                i += len(list_of_candidates)
        return scores_dict
    else:
        on_scores = on_scoring_function(batch_candidates_list)
        all_scores = [x * y for x, y in zip(off_scores, on_scores)]
        i = 0
        for gene in genes_targets_dict:
            scores_dict[gene] = []
            for target in genes_targets_dict[gene]:
                candidates_target_scores = all_scores[i:i + len(list_of_candidates)]
                scores_dict[gene].append((target, candidates_target_scores))
                i += len(list_of_candidates)
        return scores_dict


def return_candidates(genes_targets_dict: Dict[str, List[str]], omega: float, off_scoring_function, on_scoring_function,
                      node: CladeNew, cfd_dict: Dict = None) -> List[Candidate]:
    """

    :param genes_targets_dict: a dictionary of gene -> list of potential targets found in the gene
    :param omega: Threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
    :param off_scoring_function: the off target scoring function
	:param on_scoring_function: the on target scoring function
    :param node: current node in the targets UPGMA tree that the targets in genes_targets_dict belong to
    :param cfd_dict: a dictionary of mismatches and their scores for the CFD function

    :return:
    """
    # stage one: make a list of all the sgRNAs
    list_of_targets = []
    for gene in genes_targets_dict:
        list_of_targets += genes_targets_dict[gene]
    # stage two: find the suitable sgRNA:
    dict_of_different_places = find_polymorphic_sites(list_of_targets)
    node.polymorphic_sites = dict_of_different_places
    list_of_different_places = list(dict_of_different_places.items())
    list_of_different_places.sort(key=lambda item: item[0])
    # going over all the permutations
    list_of_perms_seqs = all_perms(str(list_of_targets[0]), [], list_of_different_places)
    list_of_candidates = []  # a list of tuples: (candidate_str, fraction_of_cut, cut_expectation, genes_list)
    scores_dict = generate_scores(genes_targets_dict, list_of_perms_seqs, off_scoring_function, on_scoring_function,
                                  cfd_dict)
    for i in range(len(list_of_perms_seqs)):
        targets_dict = {}  # a list of tuples: (gene name, list of target of this gene that might be cut by the candidate_str)
        genes_covering = []  # a list of tuples: (gene name, probability to be cut).
        for gene in scores_dict:
            prob_gene_will_not_cut = 1  # the probability that a gene will not be cut by the candidate
            list_of_targets = []  # for later knowing where the candidate_str might cut in each gene (when writing the output)
            num_of_cuts_per_gene = 0  # in use only in the single gene version
            for target, candidates_target_scores in scores_dict[gene]:
                candidate_cut_prob = candidates_target_scores[i]
                if candidate_cut_prob == 0:  # in case the score is 0 and
                    # there isn't attachment of the guide and target nad no cut event) don't consider the target
                    continue
                sg_site_differences = two_seqs_differences(list_of_perms_seqs[i],
                                                           target)  # the differences between the i-th candidate and the target
                list_of_targets.append([target[:20], sg_site_differences])
                prob_gene_will_not_cut = prob_gene_will_not_cut * (
                        1 - candidate_cut_prob)  # lowering the not cut prob in each sgRNA
                num_of_cuts_per_gene += candidate_cut_prob
            prob_gene_cut = 1 - prob_gene_will_not_cut
            if prob_gene_cut >= omega and len(list_of_targets) > 0:
                targets_dict[gene] = list_of_targets
                genes_covering.append((gene, prob_gene_cut))
        # check if the potential candidate covers less than two genes
        if len(genes_covering) < 2:
            continue
        cut_expectation = 0  # the probability the permutated sequence will cut all the genes, that the probability each
        # of them will be cut is greater then omega
        genes_score_dict = {}  # a dict of genes: genes considered cut by this sequence, and cut prob
        for gene_cut_prob in genes_covering:  # tuple : (gene name, probability to be cut)
            cut_expectation += gene_cut_prob[1]  # the prob to cut all the genes
            genes_score_dict[gene_cut_prob[0]] = gene_cut_prob[1]
        if len(genes_covering) >= 1:  # If the candidate has at least one gene with a score above omega, add it to the result omer 18/04
            current_candidate = Candidate(list_of_perms_seqs[i], cut_expectation, genes_score_dict, targets_dict)
            list_of_candidates.append(current_candidate)
    del list_of_perms_seqs
    return list_of_candidates


def all_perms(initial_seq: str, list_of_seqs: List[str], list_of_differences: List[Tuple[int, List[str]]]) -> List[str]:
    """Given an initial sequence and a list of possible polymorphisms and their indices in that sequence, this function
    creates a list of all the possible permutations of the initial sequence.
    list_of_seqs is initialized on the first call of the function. each recursive call adds to list_of_seqs the
    permutations produced with the next index from list_of_differences, and advances the next call to start from the
    next index in list_of_differences. the recursion stops when len(list_of_differences) = 0.
    e.g. all_perms("ACTG", list(), [(0, [T,G]), (3, [A,T])]) will return ['TCTA', 'TCTT', 'GCTA', 'GCTT']

    :param initial_seq: a sequence to create permutations for
    :param list_of_seqs: list of permutations. Initially None. The function creates it during the recursion.
    :param list_of_differences: list of tuples of polymorphisms and their locations: (index, set of nucleotides)
    :return: list of permutations of the initial sequence
    """
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
            list_of_seqs = [initial_seq[:list_of_differences[0][0]]]
        for seq in list_of_seqs:
            for letter in list_of_differences[0][1]:
                if len(list_of_differences) > 1:
                    new_list_of_seqs.append(seq + letter + initial_seq[len(seq) + 1:list_of_differences[1][0]])
                    # the place of the next versatile letter place
                else:
                    new_list_of_seqs.append(seq + letter + initial_seq[len(seq) + 1:20])
        del list_of_seqs
        return all_perms(initial_seq, new_list_of_seqs, list_of_differences[1:])


def two_seqs_differences(seq1: str, seq2: str) -> Dict[int, List[str]]:
    """Given two sequences of genomic targets this function finds the differences between them returns the index of
    the difference (starting at 0) with the residing nucleotides of each sequence at that location. The differences are
    stored in a dictionary where keys are indices of the differences and the values are lists of length 2 where value[0]
    is the nucleotide of seq1 and value[1] is the nucleotide of seq2 at the given index.
    e.g. two_seqs_differences("ACTG", "GCTA") will return {0: [A, G], 3: [G, A]}.

    :param seq1: first genomic target sequence
    :param seq2: second genomic target sequence
    :return: a dictionary representing the differences between the sequences
    """
    differences = {}  # key: place of disagreement. value: the suggestions of each side
    seq1 = seq1[:20]
    seq2 = seq2[:20]  # in cases the PAM is not sliced - don't take PAM into account
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences[i] = [seq1[i], seq2[i]]
    return differences


def find_polymorphic_sites(leaves: List) -> Dict[int, List[str]]:
    """Given a list of genomic targets sequences this function takes the first sequence (leaves_ds[0]) as a reference
    and finds the polymorphic sites between it and the rest of the sequences, and returns a dictionary representing
    these polymorphisms, where keys are indices of the polymorphisms (starting at 0) and values are lists of possible
    nucleotides at the given index.
    e.g. find_polymorphic_sites(["ACTG", "GCTA", "AGTC"]) will return:
    {0: ['A', 'G'], 1:['C', 'G'] 3: ['G', 'C', 'A']}.

    :param leaves: a list of strings of genomic targets sequences.
    :return: dictionary representing the polymorphisms between first and the other sequences of the input list.
    """
    differences = dict()  # key: place where there is a difference. value: letter appeared in this place
    if len(leaves) < 2:
        return differences
    ref = leaves[0]
    for i in range(1, len(leaves)):  # node_leaves is a python array
        current_differences = two_seqs_differences(ref, leaves[i])
        for polymorphism_site in current_differences:
            if polymorphism_site in differences:
                differences[polymorphism_site] = list(
                    set(differences[polymorphism_site] + current_differences[polymorphism_site]))
            else:
                differences[polymorphism_site] = current_differences[polymorphism_site]
    return differences
