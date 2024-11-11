from itertools import combinations
from typing import List, Dict, Tuple

from Amplicon_construction.FindTargets import find_targets_in_sequence, give_complementary
from Amplicon_construction.globals import max_polymorphic_sites
from Amplicon_construction.Target_Obj import Target_Obj, Combined_Target_Obj
from Amplicon_construction.FindOffTargets import moff
from math import prod


def get_ungapped_target_seq(position_targets_lst, exon_alleles_lst, target_len):
    """
    For every target object in the given position_targets_lst find the target sequence without gaps and add it as
    new attribute to the target object. IMPORTANT: the function considers the protospacer sequence to be located
    upstream to the PAM sequence.

    :param position_targets_lst: list of Target_Obj objects
    :param exon_alleles_lst: list of tuples representing alleles where tuple[0] is scaffold name
    (example format: ">scaffold10132:437703-438762(+)") and tuple[1] is allele sequence
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    """
    for i, target in enumerate(position_targets_lst):
        if "-" in target.seq:
            ungapped_seq_lst = list(target.seq.replace("-", ""))
            seq_to_add = []
            idx = target.start_idx if target.strand == "+" else target.end_idx
            while len(ungapped_seq_lst) + len(seq_to_add) < target_len:
                if target.strand == "+":
                    if exon_alleles_lst[i][1][idx - 1] != "-":
                        seq_to_add.append(exon_alleles_lst[i][1][idx - 1])
                    idx = idx - 1
                else:  # target.strand == "-"
                    if exon_alleles_lst[i][1][idx + 1] != "-":
                        seq_to_add.append(exon_alleles_lst[i][1][idx + 1])
                    idx = idx + 1
            if target.strand == "+":
                new_target_seq = "".join(seq_to_add) + "".join(ungapped_seq_lst)
            else:  # target.strand == "-"
                new_target_seq = give_complementary("".join(seq_to_add)) + "".join(ungapped_seq_lst)
            target.ungapped_seq = new_target_seq
        else:
            target.ungapped_seq = target.seq
    # print("Ungapped sequences updated")


def get_all_targets(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], pams: Tuple, max_amplicon_len: int,
                    primer_length: int, cut_location: int, target_surrounding_region: int,
                    target_len: int) -> Dict[int, Dict[int, List[Target_Obj]]]:
    """
    For every allele sequence find all targets (by searching PAM sequences). For every found target get sequences in the
    same indices from the other alleles. Save the targets in a Target_Obj and save the Target_Obj that start at the same
    indices in dictionary of start_index: list of Target_Obj.

    :param gene_sequences_dict: dictionary of exon num -> list of tuples representing alleles where tuple[0] is scaffold name
    (example format: ">scaffold10132:437703-438762(+)") and tuple[1] is allele sequence.
    :param pams:
    :param max_amplicon_len: maximum length of the amplicon
    :param primer_length: minimum length of the primer sequence
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :return: dictionary of exon num -> dictionary of target index in exon -> list of Target_Obj from every allele of the
    exon
    """
    targets_dict = {}
    for exon_num in gene_sequences_dict:  # iterate over every exon of the gene
        exon_targets = {}
        exon_alleles_lst = gene_sequences_dict[exon_num]
        for allele in exon_alleles_lst:  # iterate over every allele (scaffold) of the exon
            exon_region_seq = allele[1]
            allele_targets_list = find_targets_in_sequence(exon_region_seq, pams, max_amplicon_len, primer_length,
                                                      cut_location, target_surrounding_region, target_len)  # get list of all targets (as Target_Obj) for current allele
            for target in allele_targets_list:  # iterate over all targets found on current allele to save them and their parallels (on the other alleles) to a dictionary
                position_targets_lst = []
                for i in range(len(exon_alleles_lst)):  # iterate over every allele (scaffold) of the exon to get all the current target's parallels
                    scaffold = exon_alleles_lst[i][0].split(":")[0][1:]
                    if target.strand == "+":
                        cur_allele_target_seq = exon_alleles_lst[i][1][target.start_idx:target.end_idx + 1]
                    else:
                        cur_allele_target_seq = give_complementary(exon_alleles_lst[i][1][target.start_idx:target.end_idx + 1])
                    position_targets_lst.append(Target_Obj(cur_allele_target_seq, target.start_idx, target.end_idx, target.strand, scaffold))
                get_ungapped_target_seq(position_targets_lst, exon_alleles_lst, target_len)
                if target.start_idx not in exon_targets:
                    exon_targets[target.start_idx] = position_targets_lst
        targets_dict[exon_num] = exon_targets
    return targets_dict


def filter_relevant_targets(targets_dict: Dict[int, Dict[int, List[Target_Obj]]], target_len: int) -> Dict[int, List[Combined_Target_Obj]]:
    """

    :param targets_dict:
    :param target_len:
    :return:
    """
    relevant_targets_dict = {}
    for exon_num in targets_dict:
        exon_comb_targets = []
        cur_exon_targets_dict = targets_dict[exon_num]
        for target_position in cur_exon_targets_dict:
            position_targets_lst = cur_exon_targets_dict[target_position]
            zipped_seqs = zip(*[target.seq for target in position_targets_lst])
            for i, nucs in enumerate(zipped_seqs):
                if not all(nucs[0] == nuc for nuc in nucs):
                    if i != 20:  # the 'N' in the 'NGG' PAM
                        combined_target = Combined_Target_Obj(target_position, target_position+target_len, position_targets_lst)
                        exon_comb_targets.append(combined_target)
                        break

        sorted_targets = sorted(exon_comb_targets, key=lambda target: target.start_idx)
        relevant_targets_dict[exon_num] = sorted_targets

    return relevant_targets_dict


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


def get_initial_seq(comb_target: Combined_Target_Obj, strange_target_scaffolds: List[str]):
    for target in comb_target.targets_list:
        if target.scaffold not in strange_target_scaffolds:
            return target.seq


def get_list_of_differences(comb_target_snps_dict: Dict[int, Dict[str, List[str]]],
                            strange_target_scaffolds: List[str]) -> List[Tuple[int, List[str]]]:
    list_of_differences = []
    for index in comb_target_snps_dict:
        pos_nucs = [nuc for nuc in comb_target_snps_dict[index] if not (len(comb_target_snps_dict[index][nuc]) == 1 and
                                                                        comb_target_snps_dict[index][nuc][0] in strange_target_scaffolds)]
        list_of_differences.append((index, pos_nucs))

    return list_of_differences


def create_sgrna_permutations(relevant_targets_dict: Dict[int, List[Combined_Target_Obj]], pam_idxs: Tuple[int] = (22, 23)):
    """

    :param relevant_targets_dict:
    :param pam_idxs:
    """
    for exon_num in relevant_targets_dict:
        comb_target_lst = relevant_targets_dict[exon_num]
        for comb_target in comb_target_lst:
            scaffold_to_num_polymorphic_sites = {target.scaffold: 0 for target in comb_target.targets_list}
            number_to_scaffold_dict = {i: target.scaffold for i, target in enumerate(comb_target.targets_list)}
            zipped_seqs = zip(*[target.ungapped_seq for target in comb_target.targets_list])
            comb_target_snps_dict = {}  # {snp_position: {nucleotide: list of scaffold_IDs}}
            for index, nucs in enumerate(zipped_seqs):
                if not all(nucs[0] == nucleotide for nucleotide in nucs):  # SNP at current index
                    pos_target_dict = {}  # {nucleotide: list of scaffold_IDs}
                    for i in range(len(nucs)):
                        if nucs[i].upper() not in pos_target_dict:
                            pos_target_dict[nucs[i].upper()] = [number_to_scaffold_dict[i]]
                        else:
                            pos_target_dict[nucs[i].upper()].append(number_to_scaffold_dict[i])
                    for nuc in pos_target_dict:
                        if len(pos_target_dict[nuc]) == 1:
                            scaffold_to_num_polymorphic_sites[pos_target_dict[nuc][0]] += 1
                    comb_target_snps_dict[index] = pos_target_dict
                if index == pam_idxs[0] - 3:
                    break
            strange_target_scaffolds = []
            if len(comb_target_snps_dict) > max_polymorphic_sites:
                for scaffold in scaffold_to_num_polymorphic_sites:
                    if scaffold_to_num_polymorphic_sites[scaffold] >= len(comb_target_snps_dict)/2:  # check if any of the targets is "accountable" for more than half of the SNPs
                        strange_target_scaffolds.append(scaffold)
                if len(strange_target_scaffolds) < 1:  # none of the targets is "accountable" for more than half of the SNPs. Target will not be used.
                    comb_target.sg_perm = []
                    continue
                else:
                    list_of_differences = get_list_of_differences(comb_target_snps_dict, strange_target_scaffolds)
                    initial_seq = get_initial_seq(comb_target, strange_target_scaffolds)
                    comb_target.sg_perm = all_perms(initial_seq, [], list_of_differences)
            else:
                list_of_differences = get_list_of_differences(comb_target_snps_dict, strange_target_scaffolds)
                initial_seq = get_initial_seq(comb_target, strange_target_scaffolds)
                comb_target.sg_perm = all_perms(initial_seq, [], list_of_differences)
    print("Combined targets sgRNA permutations created")


def create_scores_dict(on_targets_list, off_targets_list, moff_scores):

    scores_dict = {f"{on_targets_list[i]},{off_targets_list[i]}": round(moff_scores[i], 4) for i in range(len(moff_scores))}
    return scores_dict


def calculate_off_scores(relevant_targets_dict: [Dict[int, List[Combined_Target_Obj]]], pam_idxs: Tuple[int] = (22, 23)):
    """
    Calculate the MOFF scores of every potential sgRNA against every target sequence of the Combined target and store
    the scores in Combined target's offscore_dict. Offscore_dict example:
    {"CCGGCTATGACAACCTTCAGAGG": {"scaffold346": 0.4, "scaffold406": 0.3, "scaffold68253": 0.9},
     "CCGACTATGACAACCTTCAGAGG": {"scaffold346": 0.2, "scaffold406": 0.4, "scaffold68253": 0.8}}

    :param relevant_targets_dict:
    :param pam_idxs:
    """
    on_targets_list = []
    off_targets_list = []
    for exon_num in relevant_targets_dict:
        for comb_target in relevant_targets_dict[exon_num]:
            comb_on_targets = []
            comb_off_targets = []
            comb_target.offscores_dict = {}  # dictionary of {sgRNA sequence: {target scaffold ID: MOFF score}}
            for sg in comb_target.sg_perm:
                comb_target.offscores_dict[sg] = {}  # dictionary of {target scaffold ID: MOFF score}
                for target in comb_target.targets_list:
                    faulty_pam = False  # check if targets which are not in k_scaffolds have faulty PAMs (not "NGG")
                    for pam_idx in pam_idxs:
                        if target.seq[pam_idx - 1] != "G":
                            comb_target.offscores_dict[sg][target.scaffold] = 0.0
                            faulty_pam = True
                            break

                    if not faulty_pam:
                        comb_on_targets.append(sg)
                        comb_off_targets.append(target.seq)
            on_targets_list.extend(comb_on_targets)
            off_targets_list.extend(comb_off_targets)
    moff_scores = moff(on_targets_list, off_targets_list)
    scores_dict = create_scores_dict(on_targets_list, off_targets_list, moff_scores)

    for exon_num in relevant_targets_dict:
        for comb_target in relevant_targets_dict[exon_num]:
            for sg in comb_target.sg_perm:
                for target in comb_target.targets_list:
                    if f"{sg},{target.seq}" in scores_dict:
                        comb_target.offscores_dict[sg][target.scaffold] = scores_dict[f"{sg},{target.seq}"]
    print("Combined targets MOFF scores updated")


def calculate_sg_rank_scores(relevant_targets_dict: Dict[int, List[Combined_Target_Obj]], k: int):
    """

    :param relevant_targets_dict:
    :param k:
    :return:
    """
    new_relevant_targets_dict = {}
    for exon_num in relevant_targets_dict:
        new_relevant_targets_dict[exon_num] = []
        comb_target_lst = relevant_targets_dict[exon_num]
        for comb_target in comb_target_lst:
            max_sg_score = 0.0
            chosen_sg = ""
            cut_alleles = set()
            for sg in comb_target.offscores_dict:
                combs_of_k = combinations(list(comb_target.offscores_dict[sg].items()), k)  # all combinations of size k of {target scaffold: score} for current sg
                for combo in combs_of_k:
                    product_k = prod(pair[1] for pair in combo)  # product of k "on" target scores
                    n_minus_k_lst = [pair for pair in list(comb_target.offscores_dict[sg].items()) if pair not in combo]
                    product_n_minus_k = prod((1-pair[1]) for pair in n_minus_k_lst)  # product of n-k 1 minus "off" target scores
                    tot_score = product_k * product_n_minus_k
                    if tot_score > max_sg_score:
                        max_sg_score = tot_score
                        chosen_sg = sg
                        cut_alleles = {pair[0] for pair in combo}
            comb_target.chosen_sg = chosen_sg  # without PAM
            comb_target.chosen_sg_score = max_sg_score
            comb_target.cut_alleles = cut_alleles
            if len(cut_alleles) > 0:
                new_relevant_targets_dict[exon_num].append(comb_target)

    return new_relevant_targets_dict


def get_snp_targets(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], pams: Tuple, max_amplicon_len: int,
                    primer_length: int, cut_location: int, target_surrounding_region: int, target_len: int,
                    k: int) -> Dict[int, List[Combined_Target_Obj]]:
    """

    :param gene_sequences_dict: dictionary of exon num -> list of tuples representing alleles where tuple[0] is scaffold name
    (example format: ">scaffold10132:437703-438762(+)") and tuple[1] is allele sequence.
    :param pams:
    :param max_amplicon_len: maximum length of the amplicon
    :param primer_length: minimum length of the primer sequence
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream)
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer
    :param k: number of alleles to target with a single gRNA
    :return:
    """
    all_targets_dict = get_all_targets(gene_sequences_dict, pams, max_amplicon_len, primer_length, cut_location,
                                       target_surrounding_region, target_len)
    relevant_targets_dict = filter_relevant_targets(all_targets_dict, target_len)
    create_sgrna_permutations(relevant_targets_dict)
    calculate_off_scores(relevant_targets_dict)
    new_relevant_targets_dict = calculate_sg_rank_scores(relevant_targets_dict, k)
    return new_relevant_targets_dict
