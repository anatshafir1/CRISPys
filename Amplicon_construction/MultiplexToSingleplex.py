from typing import List, Dict, Tuple

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.GetPrimers import get_gene_family_primers
from Amplicon_construction.SNP_Obj import SNP_Obj
from Amplicon_construction.Target_Obj import Target_Obj, MultiplexTarget, FamilyMultiplexTarget, Family_Target_Obj


def split_target(multiplex_target: MultiplexTarget) -> Tuple[Target_Obj, Target_Obj]:
    up_target = Target_Obj(multiplex_target.up_seq, multiplex_target.up_start, multiplex_target.up_end,
                           multiplex_target.up_targets_list[0].strand, "", multiplex_target.up_seq,
                           multiplex_target.rank + 0.1, 0, multiplex_target.exon_num)
    down_target = Target_Obj(multiplex_target.down_seq, multiplex_target.down_start, multiplex_target.down_end,
                             multiplex_target.down_targets_list[0].strand, "", multiplex_target.down_seq,
                             multiplex_target.rank + 0.2, 0, multiplex_target.exon_num)

    return up_target, down_target


def split_family_target(multiplex_target: FamilyMultiplexTarget) -> Tuple[Family_Target_Obj, Family_Target_Obj]:
    up_target = Family_Target_Obj(multiplex_target.up_seq, multiplex_target.up_start, multiplex_target.up_end,
                                  multiplex_target.up_target_strand, "", multiplex_target.up_seq, multiplex_target.rank + 0.1, 0,
                                  multiplex_target.up_targets_list, multiplex_target.exon_num)
    down_target = Family_Target_Obj(multiplex_target.down_seq, multiplex_target.down_start, multiplex_target.down_end,
                                    multiplex_target.down_target_strand, "", multiplex_target.down_seq, multiplex_target.rank + 0.2, 0,
                                    multiplex_target.down_targets_list, multiplex_target.exon_num)

    return up_target, down_target


def construct_singleplex_candidates(up_target: Target_Obj, down_target: Target_Obj,
                                    gene_snps_dict: Dict[int, List[SNP_Obj]],
                                    max_amplicon_len: int, primer_length: int, distinct_alleles_num: int,
                                    target_surrounding_region: int, min_amplicon_len: int, k: int,
                                    target_len: int,
                                    multiplex: int, failed_targets: Dict[str, List]):
    from Amplicon_construction.AmpliconConstruction import construct_amplicons

    up_exon_num = up_target.exon_num
    up_target_dict = {up_exon_num: [up_target]}
    up_relevant_snps_dict = {up_exon_num: gene_snps_dict[up_exon_num]}
    up_candidate_single_trg_amplicons = construct_amplicons(up_relevant_snps_dict, up_target_dict, max_amplicon_len,
                                                            primer_length,
                                                            distinct_alleles_num, target_surrounding_region,
                                                            min_amplicon_len,
                                                            k, target_len, multiplex)
    if not up_candidate_single_trg_amplicons:
        if up_target.seq not in failed_targets:
            failed_targets[up_target.seq] = [up_target]
        else:
            if up_target not in failed_targets[up_target.seq]:
                failed_targets[up_target.seq].append(up_target)
        return
    down_exon_num = down_target.exon_num
    down_target_dict = {down_exon_num: [down_target]}
    down_relevant_snps_dict = {down_exon_num: gene_snps_dict[down_exon_num]}
    down_candidate_singleplex_amplicons = construct_amplicons(down_relevant_snps_dict, down_target_dict,
                                                              max_amplicon_len,
                                                              primer_length, distinct_alleles_num,
                                                              target_surrounding_region,
                                                              min_amplicon_len, k, target_len, multiplex)
    if not down_candidate_singleplex_amplicons:
        if down_target.seq not in failed_targets:
            failed_targets[down_target.seq] = [down_target]
        else:
            if down_target not in failed_targets[down_target.seq]:
                failed_targets[down_target.seq].append(down_target)
        return
    else:
        return up_candidate_single_trg_amplicons, down_candidate_singleplex_amplicons


def get_singleplex_amplicons(up_target: Target_Obj, down_target: Target_Obj, gene_snps_dict, max_amplicon_len, primer_length, distinct_alleles_num,
                             target_surrounding_region, min_amplicon_len, target_len, gene_exon_regions_seqs_dict,
                             out_path, primer3_core_path, n, amplicon_range, gene_orig_exons_dict,
                             failed_targets: Dict[str, List], family_targeting: int) -> List[Amplicon_Obj]:
    from Amplicon_construction.AmpliconConstruction import filter_redundancies_and_sort
    singleplex_candidate_amplicons = construct_singleplex_candidates(up_target, down_target, gene_snps_dict,
                                                                     max_amplicon_len, primer_length,
                                                                     distinct_alleles_num,
                                                                     target_surrounding_region,
                                                                     min_amplicon_len, 0, target_len, 0, failed_targets)
    if not singleplex_candidate_amplicons:  # Failed to construct Upstream and Downstream Target Candidate Amplicons
        return
    else:  # Successfully constructed Upstream and Downstream Target Candidate Amplicons. Search for Primers for Upstream Candidate Amplicons.
        filt_sorted_up_trg_candidate_amplicons = filter_redundancies_and_sort(singleplex_candidate_amplicons[0], 0)
        up_trg_amplicons = get_gene_family_primers(gene_exon_regions_seqs_dict,
                                                   filt_sorted_up_trg_candidate_amplicons, out_path,
                                                   primer3_core_path, n,
                                                   amplicon_range,
                                                   distinct_alleles_num, target_surrounding_region,
                                                   gene_orig_exons_dict,
                                                   max_amplicon_len, gene_snps_dict, 0, min_amplicon_len,
                                                   0, family_targeting)
        if not up_trg_amplicons:  # Failed to find primers for Upstream Target Candidate Amplicons.
            if up_target.seq not in failed_targets:
                failed_targets[up_target.seq] = [up_target]
            else:
                if up_target not in failed_targets[up_target.seq]:
                    failed_targets[up_target.seq].append(up_target)
            return
        else:  # Successfully constructed Upstream Target Amplicons. Try Downstream Target Amplicons.
            filt_sorted_down_trg_candidate_amplicons = filter_redundancies_and_sort(
                singleplex_candidate_amplicons[1], 0)
            down_trg_amplicons = get_gene_family_primers(gene_exon_regions_seqs_dict,
                                                         filt_sorted_down_trg_candidate_amplicons,
                                                         out_path,
                                                         primer3_core_path, n,
                                                         amplicon_range,
                                                         distinct_alleles_num,
                                                         target_surrounding_region,
                                                         gene_orig_exons_dict,
                                                         max_amplicon_len,
                                                         gene_snps_dict, 0,
                                                         min_amplicon_len,
                                                         0, family_targeting)
            if not down_trg_amplicons:  # Failed to find primers for Downstream Target Candidate Amplicons.
                if down_target.seq not in failed_targets:
                    failed_targets[down_target.seq] = [down_target]
                else:
                    if down_target not in failed_targets[down_target.seq]:
                        failed_targets[down_target.seq].append(down_target)
                return
            else:
                return up_trg_amplicons + down_trg_amplicons
