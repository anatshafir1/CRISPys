from typing import Tuple, Dict, List

from Amplicon_construction.AmpliconConstruction import construct_amplicons, filter_redundancies_and_sort
from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.GetPrimers import get_gene_family_primers
from Amplicon_construction.MultiplexToSingleplex import split_family_target, get_singleplex_amplicons
from Amplicon_construction.SNP_Obj import SNP_Obj
from Amplicon_construction.Target_Obj import FamilyMultiplexTarget, Family_Target_Obj
from CRISPys_master.Candidate import Candidate
from Gene_Family_Targeting.GetGeneFamilyTargets import get_genes_multiplex_targets_dict, get_genes_single_targets_dict


def create_single_sgrna_amplicons(sorted_sgrnas: List[Candidate], genes_exons_seq_dict: Dict[
    str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]],
                                  genes_snps_dict: Dict[str, Dict[int, List[SNP_Obj]]], max_amplicon_len: int,
                                  primer_length: int, distinct_alleles_num: int, target_surrounding_region: int,
                                  min_amplicon_len: int, target_len: int, out_path: str, primer3_core_path: str,
                                  n: int, amplicon_ranges: List[Tuple[int, int]],
                                  max_amplicon_len_category: int, cut_location: int, family_targeting: int):
    """
    create amplicons for every gene of single sgrnas candidates

    :param sorted_sgrnas:
    :param genes_exons_seq_dict:
    :param genes_snps_dict:
    :param max_amplicon_len:
    :param primer_length:
    :param distinct_alleles_num:
    :param target_surrounding_region:
    :param min_amplicon_len:
    :param target_len:
    :param out_path:
    :param primer3_core_path:
    :param n:
    :param amplicon_ranges:
    :param max_amplicon_len_category:
    :param cut_location:
    :param family_targeting:
    :return:
    """
    sgrna_amplicons_dict = {}
    sgrna_seq_to_object_dict = {}
    for index, sgrna in enumerate(sorted_sgrnas):
        print(f"sgRNA {index+1} in progress")
        if len(sgrna_amplicons_dict) == n:
            break
        sgrna_failed = False
        sgrna_seq = sgrna.seq
        sgrna_targets_dict = get_genes_single_targets_dict(sgrna, genes_exons_seq_dict, max_amplicon_len, primer_length,
                                                           cut_location, target_surrounding_region)
        gene_amplicons_with_primers = {gene: [] for gene in sgrna_targets_dict}
        for gene_name in sgrna_targets_dict:
            amps_for_gene = False
            exons_targets_dicts_list = sgrna_targets_dict[gene_name]
            gene_orig_exons_dict = genes_exons_seq_dict[gene_name.lower()][1]
            for exons_targets_dict in exons_targets_dicts_list:
                gene_candidate_amplicons = construct_amplicons(genes_snps_dict[gene_name.lower()],
                                                               exons_targets_dict, max_amplicon_len, primer_length,
                                                               distinct_alleles_num, target_surrounding_region,
                                                               min_amplicon_len, 0, target_len, 0)
                if gene_candidate_amplicons:  # Successfully created Candidate Amplicons for current targets dict of current gene
                    filt_sorted_gene_candidate_amplicons = filter_redundancies_and_sort(gene_candidate_amplicons, 0)
                    gene_amplicons = get_gene_family_primers(genes_exons_seq_dict[gene_name.lower()][0],
                                                         filt_sorted_gene_candidate_amplicons, out_path,
                                                         primer3_core_path, n,
                                                         amplicon_ranges[max_amplicon_len_category - 1],
                                                         distinct_alleles_num, target_surrounding_region,
                                                         gene_orig_exons_dict, max_amplicon_len,
                                                         genes_snps_dict[gene_name.lower()], 0, min_amplicon_len,
                                                         0, family_targeting)
                    if gene_amplicons:  # Successfully created Amplicons for current targets dict of current gene. store results and continue to next targets dict
                        gene_amplicons_with_primers[gene_name].extend(gene_amplicons)
                        amps_for_gene = True
                    else:  # Failed to create Amplicons with primers for current targets dict of current gene. continue to next targets dict of current gene
                        continue
                else:  # Failed to create Candidate Amplicons for current targets dict of current gene. continue to next targets dict of current gene
                    continue
            if amps_for_gene:  # Successfully created Amplicons for current gene. continue to next gene
                continue
            else:  # Failed to create Amplicons for current gene. break and continue to next sgRNA
                sgrna_failed = True
                break
        if not sgrna_failed:
            sgrna_amplicons_dict[sgrna_seq] = gene_amplicons_with_primers
            sgrna_seq_to_object_dict[sgrna_seq] = sgrna
    return sgrna_amplicons_dict, sgrna_seq_to_object_dict


def get_multiplex_target(gene_targets_dict: Dict[int, List[FamilyMultiplexTarget]]) -> FamilyMultiplexTarget:
    target = None
    for exon in gene_targets_dict:
        if gene_targets_dict[exon]:
            target = gene_targets_dict[exon][0]
            break
    return target


def get_single_targets(gene_targets_dict: Dict[int, List[Family_Target_Obj]]) -> Tuple[Family_Target_Obj, Family_Target_Obj]:
    """
    Get the pair of Family Target objects from current targets dictionary. Current targets dictionary must have exactly
    2 targets.
    :param gene_targets_dict: dictionary of exon number -> list of Family Target objects
    :return:
    """
    targets_list = []
    for exon in gene_targets_dict:
        if gene_targets_dict[exon]:
            for target in gene_targets_dict[exon]:
                targets_list.append(target)

    if targets_list[0].start_idx < targets_list[1].start_idx:
        upstream_target = targets_list[0]
        downstream_target = targets_list[1]
    else:
        upstream_target = targets_list[1]
        downstream_target = targets_list[0]

    return upstream_target, downstream_target


def create_multiplex_sgrna_amplicons(sorted_sgrna_pairs: List[Tuple[Candidate, Candidate]], genes_exons_seq_dict: Dict[
    str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]], genes_snps_dict: Dict[str, Dict[int, List[SNP_Obj]]],
                                     max_amplicon_len: int, primer_length: int, distinct_alleles_num: int,
                                     target_surrounding_region: int,  min_amplicon_len: int, target_len: int, out_path: str,
                                     primer3_core_path: str,  n: int, amplicon_ranges: List[Tuple[int, int]],
                                     max_amplicon_len_category: int, cut_location: int, family_targeting: int) -> Tuple[Dict[str, Dict[str, List[Amplicon_Obj]]], Dict[str, Tuple[Candidate, Candidate]]]:
    sgrna_amplicons_dict = {}
    sgrna_pair_seq_to_object = {}
    failed_targets = {}
    for index, sgrna_pair in enumerate(sorted_sgrna_pairs):
        rank = index + 1
        if len(sgrna_amplicons_dict) == n:
            break
        print(f"pair {index+1} in progress")
        sgrna_pair_str = f"{sgrna_pair[0].seq}_{sgrna_pair[1].seq}"
        sgrna_pair_targets_dict = get_genes_multiplex_targets_dict(sgrna_pair, genes_exons_seq_dict, max_amplicon_len,
                                                                   primer_length, cut_location, target_surrounding_region,
                                                                   failed_targets, rank)
        if not sgrna_pair_targets_dict:
            continue
        sgrna_pair_failed = False
        gene_amplicons_with_primers = {}
        for gene_name in sgrna_pair_targets_dict:
            amps_for_gene = False
            exons_targets_dicts_list = sgrna_pair_targets_dict[gene_name]
            gene_orig_exons_dict = genes_exons_seq_dict[gene_name.lower()][1]
            for exons_targets_dict in exons_targets_dicts_list:
                if exons_targets_dict["multiplex"]:  # Current sgRNA pair is on same exon, attempt multiplex amplicons

                    multiplex_candidate_amplicons = construct_amplicons(genes_snps_dict[gene_name.lower()], exons_targets_dict,
                        max_amplicon_len, primer_length, distinct_alleles_num,
                        target_surrounding_region, min_amplicon_len, 0, target_len,
                        1)
                    if not multiplex_candidate_amplicons:  # Failed to construct Multiplex Candidate amplicons. Try singles
                        multiplex_target = get_multiplex_target(exons_targets_dict)
                        up_target, down_target = split_family_target(multiplex_target)
                        singleplex_amplicons = get_singleplex_amplicons(
                            up_target, down_target, genes_snps_dict[gene_name.lower()], max_amplicon_len, primer_length,
                            distinct_alleles_num, target_surrounding_region, min_amplicon_len, target_len,
                            genes_exons_seq_dict[gene_name.lower()][0], out_path, primer3_core_path, n,
                            amplicon_ranges[max_amplicon_len_category - 1], gene_orig_exons_dict, failed_targets, family_targeting)
                        if not singleplex_amplicons:  # Failed to construct Single Target amplicons. sgRNA pair failed.
                            continue
                        else:  # Successfully Created Single Target amplicons. Store results.
                            gene_amplicons_with_primers[gene_name] = singleplex_amplicons
                            amps_for_gene = True
                            break
                    else:   # Constructed multiplex candidate amplicons. Find primers
                        filt_sorted_multiplex_candidate_amplicons = filter_redundancies_and_sort(
                            multiplex_candidate_amplicons,  0)
                        multiplex_amplicons = get_gene_family_primers(
                            genes_exons_seq_dict[gene_name.lower()][0], filt_sorted_multiplex_candidate_amplicons, out_path,
                            primer3_core_path, n, amplicon_ranges[max_amplicon_len_category - 1], distinct_alleles_num,
                            target_surrounding_region, gene_orig_exons_dict, max_amplicon_len, genes_snps_dict[gene_name.lower()],
                            0, min_amplicon_len, 1, family_targeting)

                        if not multiplex_amplicons:  # Failed to find primers for multiplex candidate amplicons. Try singles
                            multiplex_target = filt_sorted_multiplex_candidate_amplicons[0].target
                            up_target, down_target = split_family_target(multiplex_target)
                            singleplex_amplicons = get_singleplex_amplicons(
                                up_target, down_target, genes_snps_dict[gene_name.lower()], max_amplicon_len, primer_length,
                                distinct_alleles_num, target_surrounding_region, min_amplicon_len, target_len,
                                genes_exons_seq_dict[gene_name.lower()][0], out_path, primer3_core_path, n,
                                amplicon_ranges[max_amplicon_len_category - 1], gene_orig_exons_dict, failed_targets, family_targeting)
                            if not singleplex_amplicons:  # Failed to construct Single Target amplicons.
                                continue
                            else:  # Failed to construct Single Target amplicons. sgRNA pair failed.
                                gene_amplicons_with_primers[gene_name] = singleplex_amplicons
                                amps_for_gene = True
                                break
                        else:  # Successfully Found Primers for Multiplex amplicons. Store results.
                            gene_amplicons_with_primers[gene_name] = multiplex_amplicons
                            amps_for_gene = True
                            break

                else:  # Current sgRNA pair is on different exons or too far apart. Construct Single Target amplicons
                    up_target, down_target = get_single_targets(exons_targets_dict)
                    singleplex_amplicons = get_singleplex_amplicons(
                        up_target, down_target, genes_snps_dict[gene_name.lower()], max_amplicon_len, primer_length,
                        distinct_alleles_num, target_surrounding_region, min_amplicon_len, target_len,
                        genes_exons_seq_dict[gene_name.lower()][0], out_path, primer3_core_path, n,
                        amplicon_ranges[max_amplicon_len_category - 1], gene_orig_exons_dict, failed_targets, family_targeting)
                    if not singleplex_amplicons:  # Failed to construct Single Target amplicons. sgRNA pair failed.
                        continue
                    else:  # Successfully Created Single Target amplicons. Store results.
                        gene_amplicons_with_primers[gene_name] = singleplex_amplicons
                        amps_for_gene = True
                        break
            if amps_for_gene:
                continue
            else:
                sgrna_pair_failed = True
                break
        if sgrna_pair_failed:  # Failed to Construct amplicons for all genes in Family for current sgRNA pair Targets Dict.
            continue
        else:   # Successfully Constructed amplicons for all genes in Family for current sgRNA pair Targets Dict. Store results.
            sgrna_amplicons_dict[sgrna_pair_str] = gene_amplicons_with_primers
            sgrna_pair_seq_to_object[sgrna_pair_str] = sgrna_pair

    return sgrna_amplicons_dict, sgrna_pair_seq_to_object
