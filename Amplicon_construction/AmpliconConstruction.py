"""Constructing Amplicons from targets and SNPs"""
import statistics
import sys
from typing import List, Tuple, Dict

import pandas as pd

from Amplicon_Obj import Amplicon_Obj
from FindOffTargets import filt_off_targets
from FindTargets import get_targets
from GetSequences import extract_exons_regions
from GetSNPs import get_snps
from GetPrimers import get_primers
from SNP_Obj import SNP_Obj
from Target_Obj import Target_Obj

import warnings

warnings.filterwarnings("ignore")


def get_relevant_targets(gene_targets_dict: Dict[int, List[Target_Obj]], gene_snps_dict: Dict[int, List[SNP_Obj]]) -> \
        Dict[int, List[Target_Obj]]:
    """
    loop over every target in every exon of the gene, and check if any SNPs of that exon are on the target sequence.

    :param gene_targets_dict: dictionary of exon numbers -> a list potential targets of the exon.
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region.
    :return:  dictionary of exon numbers -> a list potential targets of the exon without SNPs in their sequences.
    """

    def ranges_overlap(range1, range2):
        start1, end1 = range1
        start2, end2 = range2
        return start1 <= end2 and start2 <= end1

    new_targets_dict = {}
    for exon in gene_snps_dict:
        new_exon_targets_lst = []
        for target in gene_targets_dict[exon]:
            valid_target = True
            target_range = target.start_idx, target.end_idx
            for snp in gene_snps_dict[exon]:
                snp_range = snp.position, snp.position + snp.gap_length
                if ranges_overlap(target_range, snp_range):
                    valid_target = False
                    break
            if valid_target:
                new_exon_targets_lst.append(target)
        new_targets_dict[exon] = new_exon_targets_lst

    return new_targets_dict


def valid_amplicon(i: int, j: int, exon_snps_list: List[SNP_Obj], distinct_alleles_num: int) -> bool:
    """
    Check if SNPs i to j in exon_snps_list are enough to distinguish between the different alleles. The union of the sets
    of the alleles that each SNP covers should be at least the size of the number of different alleles minus 1 to create
    a set that will cover all the different alleles.

    :param i: index of current SNP.
    :param j: index of SNP downstream of "i".
    :param exon_snps_list: list of the SNPs of the current exon.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :return: True if sub-list of SNPs i to j are valid to distinguish between the alleles of an amplicon sequence.
    """
    allele_set = exon_snps_list[i].different_alleles_set
    for k in range(i + 1, j + 1):
        allele_set = allele_set.union(exon_snps_list[k].different_alleles_set)
        if len(allele_set) >= distinct_alleles_num - 1:
            return True
    return False


def valid_sgrna_target(i: int, j: int, exon_snps_list: List[SNP_Obj], distinct_alleles_num: int,
                       curr_target, max_amplicon_len: int, primer_length: int, target_surrounding_region: int,
                       k: int) -> List[SNP_Obj]:
    """
    Check if current target can be used to construct an amplicon.

    :param i: index of current SNP.
    :param j: index of SNP downstream of "i".
    :param exon_snps_list: list of the SNPs of the current exon.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param curr_target: potential target to construct an amplicon with.
    :param max_amplicon_len: maximum length of the amplicon.
    :param primer_length: minimum length of the primer sequence.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param k: number of alleles to target.
    :return: list of SNPs that are valid to distinguish between the alleles of an amplicon sequence.
    """
    valid_snps_list = []
    max_range_from_snp = max_amplicon_len - primer_length * 2 - target_surrounding_region - 1
    first_snp = exon_snps_list[i]
    last_snp = exon_snps_list[j]
    # Make sure that target is not too far upstream or downstream from first and last SNPs
    if not (curr_target.end_idx < first_snp.position + max_range_from_snp and curr_target.start_idx >
            last_snp.position - max_range_from_snp):
        return []
    # Create a list of SNPs that can be used to calculate the Amplicon statistics
    for m in range(i, j + 1):
        # Add to valid_snps_list SNPs that do not fall on the target or the target buffer region - for all Tools
        if not curr_target.start_idx - target_surrounding_region <= exon_snps_list[m].position <= \
               curr_target.end_idx + target_surrounding_region:
            valid_snps_list.append(exon_snps_list[m])
        if k > 0:  # Tool 2 in use.
            #  Add to valid_snps_list SNPs that fall on the target and match the expected alleles to be cut
            if curr_target.start_idx <= exon_snps_list[m].position <= curr_target.end_idx:  # SNP on target
                if exon_snps_list[m].different_alleles_set == curr_target.cut_alleles:  # SNP match expected alleles to be cut
                    valid_snps_list.append(exon_snps_list[m])
    if len(valid_snps_list) == 0:  # if no SNPs can be used do discern between the alleles
        return []
    # Make sure that the new SNPs list in enough for a valid amplicon (enough to distinguish between all the alleles)
    if valid_amplicon(0, len(valid_snps_list) - 1, valid_snps_list, distinct_alleles_num):
        return valid_snps_list
    return []


def calculate_snps_statistics(valid_snps_for_target: List[SNP_Obj], all_alleles_set: set) -> Tuple[float, float]:
    """
    given a list of SNPs that are valid to distinguish between the alleles of an amplicon sequence, calculate the
    median and mean number of SNPs on each allele.

    :param valid_snps_for_target: list of SNPs for amplicon statistics calculation.
    :param all_alleles_set:  set of all the scaffold IDs.
    :return: tuple of median and mean number of SNPs on each allele.
    """
    counts_dict = {allele: 0 for allele in all_alleles_set}
    for snp in valid_snps_for_target:
        for allele in snp.different_alleles_set:
            counts_dict[allele] += 1
    snps_median = statistics.median(list(counts_dict.values()))
    snps_mean = statistics.mean(list(counts_dict.values()))

    return round(snps_median, 4), round(snps_mean, 4)


# noinspection PyTypeChecker
def create_candidate_amplicon(snps_median: float, snps_mean: float, candidate_amplicon_snps_lst: List[SNP_Obj],
                              target, max_amplicon_len: int, primer_length: int, exon_num: int,
                              target_surrounding_region: int, min_amplicon_len: int, target_len: int) -> Amplicon_Obj:
    """
    Given the parameters of the potential amplicon, calculate potential range of the amplicon and create and return an
    Amplicon_Obj object, if the range is more than the minimum.

    :param snps_median: median number of SNPs on each allele.
    :param snps_mean: mean number of SNPs on each allele.
    :param candidate_amplicon_snps_lst: list of all SNPs in the potential range of the candidate amplicon
    :param target: potential target to construct an amplicon with.
    :param max_amplicon_len: maximum length of the amplicon.
    :param primer_length: minimum length of the primer sequence.
    :param exon_num: current exon number, for which amplicons are constructed.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param min_amplicon_len: minimum length of the amplicon, defined by user.
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer.
    :return: potential amplicon as an Amplicon_Obj.
    """
    first_snp = candidate_amplicon_snps_lst[0]
    last_snp = candidate_amplicon_snps_lst[-1]
    max_range_from_snp = max_amplicon_len - primer_length - 1
    max_range_from_target = max_amplicon_len - target_len - target_surrounding_region - primer_length
    min_amp_start_index = target.start_idx - max_range_from_target if target.end_idx > last_snp.position else last_snp.position - max_range_from_snp
    max_amp_end_index = target.end_idx + max_range_from_target if target.end_idx < first_snp.position else first_snp.position + max_range_from_snp
    if max_amp_end_index - min_amp_start_index < min_amplicon_len:  # assert amplicon is not too short
        return None
    return Amplicon_Obj(exon_num, min_amp_start_index, max_amp_end_index, snps_median, snps_mean, target,
                        candidate_amplicon_snps_lst)


def get_candidate_amplicons(i: int, j: int, exon_snps_lst: List[SNP_Obj], distinct_alleles_num: int,
                            exon_targets_lst: List, max_amplicon_len: int, primer_length: int, exon_num: int,
                            target_surrounding_region: int, min_amplicon_len: int, all_alleles_set: set, k: int,
                            target_len: int) -> List[Amplicon_Obj]:
    """
    Given SNPs i to j - try to construct amplicon with every target of the current exon. return valid amplicons.

    :param i: index of current SNP.
    :param j: index of SNP downstream of "i".
    :param exon_snps_lst: list of the SNPs of the current exon.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param exon_targets_lst: list of potential targets of the current exon.
    :param max_amplicon_len: maximum length of the amplicon, defined by user.
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run.
    :param exon_num: current exon number, for which amplicons are constructed.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param min_amplicon_len: minimum length of the amplicon, defined by user.
    :param all_alleles_set: set of all the scaffold IDs.
    :param k: number of alleles to target.
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer.
    :return: list of potential amplicons as Amplicon_Obj objects for the current exon.
    """

    candidate_amplicons_list = []
    for target in exon_targets_lst:
        # check if a valid amplicon can be constructed with given SNPs and current target
        valid_snps_for_target = valid_sgrna_target(i, j, exon_snps_lst, distinct_alleles_num, target,
                                                   max_amplicon_len, primer_length, target_surrounding_region,
                                                   k)  # SNPs for Amplicon statistics
        if len(valid_snps_for_target) > 0:
            candidate_amplicon_snps_lst = exon_snps_lst[i:j + 1]  # all the SNPs of the Amplicon Candidate
            snps_median, snps_mean = calculate_snps_statistics(valid_snps_for_target,
                                                               all_alleles_set)  # statistics are counted for SNPs that are not located near the target (target surrounding region)
            if snps_median > 0:
                candidate_amplicon = create_candidate_amplicon(snps_median, snps_mean, candidate_amplicon_snps_lst,
                                                               target, max_amplicon_len, primer_length, exon_num,
                                                               target_surrounding_region, min_amplicon_len, target_len)
                if candidate_amplicon is not None:  # if amplicon wasn't too short
                    candidate_amplicons_list.append(candidate_amplicon)

    return candidate_amplicons_list


def construct_amplicons(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                        gene_snps_dict: Dict[int, List[SNP_Obj]], gene_targets_dict: Dict[int, List[Target_Obj]],
                        max_amplicon_len: int, primer_length: int, distinct_alleles_num: int,
                        target_surrounding_region: int, min_amplicon_len: int, k: int, target_len: int) -> List[Amplicon_Obj]:
    """
    Given a dictionary of sgRNA targets, a dictionary of SNP of a gene and a dictionary of exon region sequences -
    for every SNP of every exon find possible amplicons that can be constructed using that SNP.

    :param gene_exon_regions_seqs_dict: dictionary of exon numbers -> a list of allele sequences and their description.
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region.
    :param gene_targets_dict: dictionary of exon numbers -> a list potential targets of the exon.
    :param max_amplicon_len: maximum length of the amplicon.
    :param primer_length: minimum length of the primer sequence.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param min_amplicon_len: minimum length of the amplicon, defined by user.
    :param k: number of alleles to target.
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer.
    :return: list of potential amplicons as Amplicon_Obj objects.
    """
    all_alleles_set = set(seq_tup[0].split(":")[0][1:] for seq_tup in gene_exon_regions_seqs_dict[1])
    candidate_amplicons_list = []
    max_dist_snp = max_amplicon_len - primer_length * 2 - 1
    for exon in gene_snps_dict:
        exon_targets_lst = gene_targets_dict[exon]
        if len(exon_targets_lst) == 0:  # if no targets were found for current exon - skip over current exon
            continue
        current_exon_region = gene_exon_regions_seqs_dict[exon][0][1]
        exon_snps_lst = gene_snps_dict[exon]
        max_range_from_snp = min(len(current_exon_region), max_dist_snp)
        for i in range(len(exon_snps_lst)):
            j = i + 1
            current_snp = exon_snps_lst[i]
            while j < len(exon_snps_lst) and exon_snps_lst[j].position < current_snp.position + max_range_from_snp:
                if valid_amplicon(i, j, exon_snps_lst,
                                  distinct_alleles_num):  # check if snps i to j are enough do distinct between different alleles
                    current_snp_candidate_amplicons = get_candidate_amplicons(i, j, exon_snps_lst, distinct_alleles_num,
                                                                              exon_targets_lst, max_amplicon_len,
                                                                              primer_length, exon,
                                                                              target_surrounding_region,
                                                                              min_amplicon_len, all_alleles_set, k,
                                                                              target_len)
                    candidate_amplicons_list.extend(current_snp_candidate_amplicons)
                j += 1

    return candidate_amplicons_list


def save_results_to_csv(res_amplicons_lst: List[Amplicon_Obj], out_path: str, k: int):
    """
    create a dataframe of the result amplicons with all their parameters and save it as a CSV file in out_path.

    :param res_amplicons_lst: list of resulted amplicons objects from the algorithm run.
    :param out_path: path to output directory where algorithm results will be saved.
    :param k: number of alleles to target with a single gRNA.
    """
    dicts_list = []
    for rank, amplicon in enumerate(res_amplicons_lst):
        dicts_list.extend([amplicon.scaffold_amplicons[scaffold_amplicon].to_dict(rank+1, k) for scaffold_amplicon in
                           amplicon.scaffold_amplicons])
    df = pd.DataFrame(dicts_list)
    df.to_csv(out_path + "/results.csv", index=False)


def get_candidates_scaffold_positions(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]]) -> Dict[str, Tuple[int, int]]:
    """
    Create a dictionary of allele scaffold: start,end indices.

    :param gene_exon_regions_seqs_dict: dictionary of exon number -> list of tuples of allele IDs and their sequences.
    :return: dictionary of allele scaffold: start,end indices.
    """
    scaffolds = []
    positions = []
    last_exon_num = len(gene_exon_regions_seqs_dict)
    for first_exon, last_exon in zip(gene_exon_regions_seqs_dict[1], gene_exon_regions_seqs_dict[last_exon_num]):
        if first_exon[0][-2] == "+":  # consider strand
            gene_start = int(first_exon[0].split(":")[1][:-3].split("-")[0])
            gene_end = int(last_exon[0].split(":")[1][:-3].split("-")[1])
        else:
            gene_start = int(last_exon[0].split(":")[1][:-3].split("-")[0])
            gene_end = int(first_exon[0].split(":")[1][:-3].split("-")[1])
        scaffolds.append(first_exon[0].split(":")[0][1:])
        positions.append((gene_start, gene_end))
    candidates_scaffold_positions = {scaffolds[i]: positions[i] for i in range(len(scaffolds))}
    return candidates_scaffold_positions


def get_amplicons(max_amplicon_len_category: int, primer_length: int, target_surrounding_region: int, cut_location: int,
                  annotations_file_path: str,
                  out_path: str, genome_fasta_file: str, distinct_alleles_num: int, pams: Tuple[str], target_len: int,
                  primer3_core_path: str, n: int, filter_off_targets: int, k: int) -> List[Amplicon_Obj]:
    """

    :param max_amplicon_len_category: category of maximum length of the amplicon, defined by user.
    :param primer_length: minimum length of the primer sequence, defined by the user in the algorithm run.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param cut_location: number of nucleotides upstream to the PAM sequence where the Cas should cut (negative number if downstream).
    :param annotations_file_path: path to GFF file with annotations of the genome.
    :param out_path: path to output directory where algorithm results will be saved.
    :param genome_fasta_file: path to input FASTA format file of the genome.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param pams: tuple of PAM sequences of the Cas protein in use.
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer.
    :param primer3_core_path: A string of the full path of the primer3 core file.
    :param n: desired maximum number of amplicons to return.
    :param filter_off_targets: choose whether to filter amplicons with 'strong' off-targets for their gRNAs, or return
    them in the results.
    :param k: number of alleles to target with a single gRNA.
    :return: list of amplicons for given gene.
    """
    amplicon_ranges = [(200, 300), (300, 500), (500, 1000)]
    max_amplicon_len = max(amplicon_ranges[max_amplicon_len_category - 1])
    min_amplicon_len = min(amplicon_ranges[max_amplicon_len_category - 1])
    gene_exon_regions_seqs_dict, original_exon_indices_dict = extract_exons_regions(max_amplicon_len, primer_length,
                                                                                    target_surrounding_region,
                                                                                    cut_location,
                                                                                    annotations_file_path, out_path,
                                                                                    genome_fasta_file)
    gene_snps_dict = get_snps(gene_exon_regions_seqs_dict, distinct_alleles_num, primer_length)
    gene_targets_dict = get_targets(gene_exon_regions_seqs_dict, pams, max_amplicon_len, primer_length, cut_location,
                                    target_surrounding_region, target_len, k)
    if k > 0:  # Tool 2 in use, targeting k alleles. SNPs allowed in target sequences.
        relevant_gene_targets_dict = gene_targets_dict
    else:  # Tool 1 in use. No SNPs allowed in target sequences.
        relevant_gene_targets_dict = get_relevant_targets(gene_targets_dict, gene_snps_dict)
    candidate_amplicons_list = construct_amplicons(gene_exon_regions_seqs_dict, gene_snps_dict,
                                                   relevant_gene_targets_dict, max_amplicon_len, primer_length,
                                                   distinct_alleles_num, target_surrounding_region, min_amplicon_len, k,
                                                   target_len)
    sorted_candidate_amplicons = sorted(candidate_amplicons_list,
                                        key=lambda amplicon: (amplicon.snps_median, amplicon.snps_mean), reverse=True)
    # create a dictionary of current gene scaffold:
    candidates_scaffold_positions = get_candidates_scaffold_positions(gene_exon_regions_seqs_dict)
    if filter_off_targets:  # find gRNA off-targets and filter amplicons by scores, then get primers.
        filtered_sorted_candidate_amplicons = filt_off_targets(sorted_candidate_amplicons.copy(), genome_fasta_file,
                                                               out_path, pams, candidates_scaffold_positions, k)
        if len(filtered_sorted_candidate_amplicons) == 0:
            print("Zero Amplicons left after filtering strong Off-Targets")
            sys.exit()
        amplicon_obj_with_primers = get_primers(gene_exon_regions_seqs_dict, filtered_sorted_candidate_amplicons,
                                                out_path,
                                                primer3_core_path, n, amplicon_ranges[max_amplicon_len_category - 1],
                                                distinct_alleles_num, target_surrounding_region, filter_off_targets,
                                                genome_fasta_file, pams, candidates_scaffold_positions,
                                                original_exon_indices_dict, max_amplicon_len, gene_snps_dict, k)
    else:  # Get primers, then find gRNA off-targets.
        amplicon_obj_with_primers = get_primers(gene_exon_regions_seqs_dict, sorted_candidate_amplicons, out_path,
                                                primer3_core_path, n, amplicon_ranges[max_amplicon_len_category - 1],
                                                distinct_alleles_num, target_surrounding_region, filter_off_targets,
                                                genome_fasta_file, pams, candidates_scaffold_positions,
                                                original_exon_indices_dict, max_amplicon_len, gene_snps_dict, k)
    if len(amplicon_obj_with_primers) > 0:
        if k > 0:
            sorted_amplicon_obj_with_primers = sorted(amplicon_obj_with_primers,
                                                      key=lambda amplicon: (amplicon.target.chosen_sg_score, amplicon.snps_median,
                                                      amplicon.snps_mean), reverse=True)
        else:
            sorted_amplicon_obj_with_primers = sorted(amplicon_obj_with_primers,
                                                      key=lambda amplicon: (amplicon.snps_median, amplicon.snps_mean),
                                                      reverse=True)
        save_results_to_csv(sorted_amplicon_obj_with_primers, out_path, k)
        return sorted_amplicon_obj_with_primers
    else:
        if max_amplicon_len_category < 3:
            get_amplicons(max_amplicon_len_category + 1, primer_length, target_surrounding_region, cut_location,
                          annotations_file_path, out_path, genome_fasta_file, distinct_alleles_num, pams, target_len,
                          primer3_core_path, n, filter_off_targets, k)
        else:
            print("No amplicons found")
