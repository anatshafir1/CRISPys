from typing import List, Dict, Tuple

from Amplicon_construction.AmpliconConstruction import construct_amplicons, filter_redundancies_and_sort
from Amplicon_construction.MultiplexToSingleplex import split_target
from Amplicon_construction.SNP_Obj import SNP_Obj
from Amplicon_construction.Target_Obj import MultiplexTarget
from Gene_Family_Targeting.CreatesgRNAamplicons import get_singleplex_amplicons
from Gene_Multiplexing.GetMultiplexTargets import get_multiplex_targets_dict


def create_multiplex_amplicons(multiplex_targets_list: List[MultiplexTarget],
                               gene_exon_seqs_dict: Dict[int, List[Tuple[str, str]]],
                               gene_snps_dict: Dict[int, List[SNP_Obj]], original_exon_indices_dict: Dict[str, Dict[int, int]], max_amplicon_len: int,
                               primer_length: int, distinct_alleles_num: int = 3,
                               target_surrounding_region: int = 20,
                               min_amplicon_len: int = 200, target_len: int = 23, out_path: str = "",
                               primer3_core_path: str = "",
                               n: int = 5, amplicon_ranges: List[Tuple[int, int]] = [(200, 300)],
                               max_amplicon_len_category: int = 0, cut_location: int = 7):
    multiplex_amplicons_dict = {}
    singleplex_amplicons_dict = {}
    multiplex_seq_to_object = {}
    failed_targets = {}
    for index, multiplex_target in enumerate(multiplex_targets_list):
        if len(multiplex_amplicons_dict) == n:
            break
        print(f"pair {index + 1} in progress")
        multiplex_target_str = f"{multiplex_target.up_seq}_{multiplex_target.down_seq}"
        sgrna_pair_targets_dict = get_multiplex_targets_dict(multiplex_target, gene_exon_seqs_dict)

        multiplex_target_failed = False
        gene_amplicons_with_primers = {}
        multiplex_candidate_amplicons, targets_for_singleplex_lst = construct_amplicons(
            gene_snps_dict, sgrna_pair_targets_dict,
            max_amplicon_len, primer_length, distinct_alleles_num,
            target_surrounding_region, min_amplicon_len, 0, target_len,
            1)
        if not multiplex_candidate_amplicons:  # Failed to construct Multiplex Candidate amplicons. Try singles
            up_target, down_target = split_target(multiplex_target)
            singleplex_amplicons = get_singleplex_amplicons(
                up_target, down_target, gene_snps_dict, max_amplicon_len, primer_length,
                distinct_alleles_num, target_surrounding_region, min_amplicon_len, target_len,
                gene_exon_seqs_dict, out_path, primer3_core_path, n,
                amplicon_ranges[max_amplicon_len_category - 1], original_exon_indices_dict, failed_targets)
            if not singleplex_amplicons:  # Failed to construct Single Target amplicons. sgRNA pair failed.
                continue
            else:  # Successfully Created Single Target amplicons. Store results.
                singleplex_amplicons_dict[multiplex_target_str] = singleplex_amplicons
                multiplex_seq_to_object[multiplex_target_str] = multiplex_target
        else:  # Constructed multiplex candidate amplicons. Find primers
            filt_sorted_multiplex_candidate_amplicons = filter_redundancies_and_sort(
                multiplex_candidate_amplicons, [], 0,
                0)
            multiplex_amplicons = get_gene_family_primers(
                genes_exons_seq_dict[gene_name][0], filt_sorted_multiplex_candidate_amplicons, out_path,
                primer3_core_path, n, amplicon_ranges[max_amplicon_len_category - 1], distinct_alleles_num,
                target_surrounding_region, gene_orig_exons_dict, max_amplicon_len, genes_snps_dict[gene_name],
                0, min_amplicon_len, 1)

            if not multiplex_amplicons:  # Failed to find primers for multiplex candidate amplicons. Try singles
                multiplex_target = filt_sorted_multiplex_candidate_amplicons[0].target
                up_target, down_target = split_family_target(multiplex_target)
                singleplex_amplicons = get_singleplex_amplicons(
                    up_target, down_target, genes_snps_dict[gene_name], max_amplicon_len, primer_length,
                    distinct_alleles_num, target_surrounding_region, min_amplicon_len, target_len,
                    genes_exons_seq_dict[gene_name][0], out_path, primer3_core_path, n,
                    amplicon_ranges[max_amplicon_len_category - 1], gene_orig_exons_dict, failed_targets)
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

