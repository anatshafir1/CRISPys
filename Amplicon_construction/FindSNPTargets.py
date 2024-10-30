from typing import List, Dict, Tuple, Set

from Amplicon_construction.FindTargets import find_targets_in_sequence, give_complementary
from Target_Obj import Target_Obj, Combined_Target_Obj
from FindOffTargets import moff


def get_all_targets(gene_sequences_dict: Dict[int, List[Tuple[str, str]]], pams: Tuple, max_amplicon_len: int,
                    primer_length: int, cut_location: int, target_surrounding_region: int,
                    target_len: int) -> Dict[int, Dict[int, List[Target_Obj]]]:
    """

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
                                                      cut_location,
                                                      target_surrounding_region,
                                                      target_len)  # get list of all targets (as Target_Obj) for current allele
            for target in allele_targets_list:  # iterate over all targets found on current allele to save them and their parallels (on the other alleles) to a dictionary
                position_targets_lst = []
                for i in range(len(exon_alleles_lst)):  # iterate over every allele (scaffold) of the exon to get all the current target's parallels
                    scaffold = exon_alleles_lst[i][0].split(":")[0][1:]
                    if target.strand == "+":
                        cur_allele_target_seq = exon_alleles_lst[i][1][target.start_idx:target.end_idx + 1]
                    else:
                        cur_allele_target_seq = give_complementary(
                            exon_alleles_lst[i][1][target.start_idx:target.end_idx + 1])
                    position_targets_lst.append(
                        Target_Obj(cur_allele_target_seq, target.start_idx, target.end_idx, target.strand, scaffold))
                if target.start_idx not in exon_targets:
                    exon_targets[target.start_idx] = position_targets_lst
        targets_dict[exon_num] = exon_targets
    return targets_dict


def get_pos_target_snps_dict(position_targets_lst: List[Target_Obj]) -> Dict[int, Dict[str, Set[str]]]:
    """

    :param position_targets_lst:
    :return:
    """

    target_snps_dict = {}
    target_seq_lst = [target.seq for target in position_targets_lst]
    number_to_scaffold_dict = {i: target.scaffold for i, target in enumerate(position_targets_lst)}
    target_seq_zip = zip(*target_seq_lst)
    for index, locus in enumerate(target_seq_zip):
        if not all(locus[0] == nucleotide for nucleotide in locus):
            pos_target_dict = {"A": set(), "T": set(), "C": set(), "G": set(), "N": set(), "-": set()}
            for i in range(len(locus)):
                pos_target_dict[locus[i].upper()].add(number_to_scaffold_dict[i])
            target_snps_dict[index] = pos_target_dict
    return target_snps_dict


def get_k_scaffolds(position_targets_snps_dict: Dict[int, Dict[str, Set[str]]], k: int, all_scaffolds_set: Set[str],
                    pam_idxs: Tuple[int] = (22, 23)) -> List[Set[str]]:
    """

    :param position_targets_snps_dict:
    :param k:
    :param pam_idxs:
    :param all_scaffolds_set:
    :return: list of unique sets of size k with scaffold ids
    """
    pam_start_idx = pam_idxs[0]
    pam_end_idx = pam_idxs[1]
    k_scaffolds_list = []
    all_values_list = []
    faulty_pam_scaffolds = set()  # set for scaffold ids whose target PAM sequences are faulty. e.g. for Cas9 the PAM is not "NGG"
    for pos in position_targets_snps_dict:  # create a list of all the values (sets of scaffold names) of all the dictionaries
        curr_pos_scaffold_sets = [value for value in position_targets_snps_dict[pos].values() if len(value) > 0]
        all_values_list.extend(curr_pos_scaffold_sets)
        if pos in range(pam_start_idx, pam_end_idx + 1):
            curr_pos_nuc_scaffolds_set = position_targets_snps_dict[pos]["G"]  # current position in PAM must be "G"
            scaffolds_not_in_set = curr_pos_nuc_scaffolds_set.difference(
                all_scaffolds_set)  # set of scaffolds that do not have "G" in current position of PAM
            faulty_pam_scaffolds.update(scaffolds_not_in_set)

    if k == 1:  # trivial case
        all_k_scaffolds_list = [value for value in all_values_list if
                                len(value) == 1 and len(value.intersection(faulty_pam_scaffolds)) == 0]
        unique_list_of_sets = list(set(frozenset(s) for s in all_k_scaffolds_list))
        unique_k_scaffolds_list = [set(s) for s in unique_list_of_sets]
        return unique_k_scaffolds_list
    for scaffold_ids in all_values_list:  # k > 1
        good_set = False
        if len(scaffold_ids) == k and len(scaffold_ids.intersection(faulty_pam_scaffolds)) == 0:
            good_set = True
            for other_scaffold_ids in all_values_list:
                if len(scaffold_ids.intersection(other_scaffold_ids)) < k:
                    good_set = False
                    break
        if good_set:
            k_scaffolds_list.append(scaffold_ids)
    return k_scaffolds_list


def find_relevant_targets(targets_dict: Dict[int, Dict[int, List[Target_Obj]]], k: int, all_scaffolds_set: Set[str], target_len: int) -> \
                        Dict[int, List[Combined_Target_Obj]]:
    """

    :param targets_dict:
    :param k: number of alleles to target with a single gRNA
    :param all_scaffolds_set:
    :param target_len:
    :return:
    """
    relevant_targets_dict = {}
    for exon_num in targets_dict:
        exon_comb_targets = []
        cur_exon_targets_dict = targets_dict[exon_num]
        for pos_target in cur_exon_targets_dict:
            position_targets_lst = cur_exon_targets_dict[pos_target]
            position_targets_snps_dict = get_pos_target_snps_dict(position_targets_lst)
            if len(position_targets_snps_dict) == 0:
                continue
            else:
                k_scaffolds_list = get_k_scaffolds(position_targets_snps_dict, k, all_scaffolds_set)
                if len(k_scaffolds_list) > 0:
                    combined_target = Combined_Target_Obj(pos_target, pos_target+target_len, position_targets_lst,
                                                          position_targets_snps_dict, k_scaffolds_list)
                    exon_comb_targets.append(combined_target)

        sorted_targets = sorted(exon_comb_targets, key=lambda target: target.start_idx)
        relevant_targets_dict[exon_num] = sorted_targets

    return relevant_targets_dict


def create_scores_dict(on_targets_list, off_targets_list, moff_scores):

    scores_dict = {f"{on_targets_list[i]},{off_targets_list[i]}": round(moff_scores[i], 4) for i in range(len(moff_scores))}
    return scores_dict


def calculate_off_scores(relevant_targets_dict, pam_idxs: Tuple[int] = (22, 23)):

    on_targets_list = []
    off_targets_list = []
    for exon_num in relevant_targets_dict:
        for comb_target in relevant_targets_dict[exon_num]:
            comb_on_targets = []
            comb_off_targets = []
            comb_target.offscores_dict = {}
            for k_scaffolds in comb_target.k_alleles_list:  # iterate over every set of k scaffolds
                k_scaffolds_on_targets = []
                k_scaffolds_off_targets = []
                curr_k_scaffolds_set = k_scaffolds
                curr_k_scaffolds_tup = tuple(curr_k_scaffolds_set)
                comb_target.offscores_dict[curr_k_scaffolds_tup] = {}
                for target in comb_target.targets_list:
                    if target.scaffold not in curr_k_scaffolds_set:
                        faulty_pam = False  # check if targets which are not in k_scaffolds have faulty PAMs (not "NGG")
                        for pam_idx in pam_idxs:
                            if pam_idx-1 in comb_target.snp_dict:
                                if target.scaffold not in comb_target.snp_dict[pam_idx-1]["G"]:
                                    comb_target.offscores_dict[curr_k_scaffolds_tup][target.scaffold] = 1
                                    faulty_pam = True
                                    break
                        if not faulty_pam:  # check if targets which are not in k_scaffolds have indels in sequence while k_scaffolds don't
                            indel_in_target = False
                            for pos in comb_target.snp_dict:
                                if target.scaffold in comb_target.snp_dict[pos]["-"] and curr_k_scaffolds_set.intersection(comb_target.snp_dict[pos]["-"]) < 1:
                                    comb_target.offscores_dict[curr_k_scaffolds_tup][target.scaffold] = 2
                                    indel_in_target = True
                                    break
                        else:
                            continue
                        if not indel_in_target:  # not faulty PAM and no indels in sequence - "regular" off-target
                            k_scaffolds_off_targets.append(target.seq)
                    else:
                        k_scaffolds_on_targets.append(target.seq)
                #  equalize lengths of comb_off_targets and comb_on_targets
                if len(k_scaffolds_off_targets) > 0:
                    if len(k_scaffolds_on_targets) < len(k_scaffolds_off_targets):
                        for _ in range(len(k_scaffolds_off_targets)-len(k_scaffolds_on_targets)):
                            k_scaffolds_on_targets.append(k_scaffolds_on_targets[0])
                    else:
                        k_scaffolds_on_targets = k_scaffolds_on_targets[:len(k_scaffolds_off_targets)]
                else:
                    k_scaffolds_on_targets = []
                comb_on_targets.extend(k_scaffolds_on_targets)
                comb_off_targets.extend(k_scaffolds_off_targets)
            on_targets_list.extend(comb_on_targets)
            off_targets_list.extend(comb_off_targets)
    moff_scores = moff(on_targets_list, off_targets_list)
    scores_dict = create_scores_dict(on_targets_list, off_targets_list, moff_scores)

    for exon_num in relevant_targets_dict:
        for comb_target in relevant_targets_dict[exon_num]:
            for k_scaffolds in comb_target.k_alleles_list:  # iterate over every set of k scaffolds
                curr_k_scaffolds_tup = tuple(k_scaffolds)
                for on_target in comb_target.targets_list:
                    if on_target.scaffold in k_scaffolds:
                        for off_target in comb_target.targets_list:
                            if off_target.scaffold not in k_scaffolds:
                                if f"{on_target.seq},{off_target.seq}" in scores_dict:
                                    comb_target.offscores_dict[curr_k_scaffolds_tup][off_target.scaffold] = scores_dict[f"{on_target.seq},{off_target.seq}"]


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
    all_scaffolds_set = {gene_sequences_dict[1][i][0].split(":")[0][1:] for i in range(len(gene_sequences_dict[1]))}
    all_targets_dict = get_all_targets(gene_sequences_dict, pams, max_amplicon_len, primer_length, cut_location,
                                       target_surrounding_region, target_len)
    relevant_targets_dict = find_relevant_targets(all_targets_dict, k, all_scaffolds_set, target_len)
    calculate_off_scores(relevant_targets_dict)
    return relevant_targets_dict
