from typing import Tuple, Dict, List

import pandas as pd

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.FindPrimerOffTargets import multiplex_primers_off_targets
from Amplicon_construction.GetSNPs import get_snps
from Amplicon_construction.GetSequences import extract_exons_regions
from Amplicon_construction.Target_Obj import MultiplexTarget
from Gene_Multiplexing.CreateMultiplexAmplicons import create_multiplex_amplicons
from Gene_Multiplexing.GetMultiplexOffTargets import filter_multiplex_off_targets, get_off_target_from_df, \
    get_off_targets
from Gene_Multiplexing.GetMultiplexTargets import get_multiplex_targets


def save_sgrnas_results_dict(multiplex_seq_to_object: Dict[str, MultiplexTarget], out_path: str):
    multiplex_results_dict = {"vector_ID": [], "multiplex_score": [], "upstream_sgRNA": [], "downstream_sgRNA": []}
    for multiplex in multiplex_seq_to_object:
        multiplex_results_dict["vector_ID"].append(multiplex)
        multiplex_results_dict["multiplex_score"].append(multiplex_seq_to_object[multiplex].multiplex_score)
        multiplex_results_dict["upstream_sgRNA"].append(multiplex_seq_to_object[multiplex].up_seq)
        multiplex_results_dict["downstream_sgRNA"].append(multiplex_seq_to_object[multiplex].down_seq)
    df = pd.DataFrame(multiplex_results_dict)
    df.to_csv(out_path + "/multiplex_results.csv", index=False)


def save_sgrna_amplicons_dict(multiplex_amplicons_dict: Dict[str, List[Amplicon_Obj]], out_path: str):

    records = []
    for vector_id in multiplex_amplicons_dict:
        for amp in multiplex_amplicons_dict[vector_id]:
            for scaffold_amp in amp.scaffold_amplicons:
                params_dict = {
                    "vector_ID": vector_id
                }
                params_dict.update(amp.scaffold_amplicons[scaffold_amp].to_dict(amp.target.rank, 0, 1, 0))
                records.append(params_dict)
    df = pd.DataFrame(records)
    df.to_csv(out_path + "/multiplex_amplicons_results.csv", index=False)


def save_results(multiplex_amplicons_dict: Dict[str, List[Amplicon_Obj]], multiplex_seq_to_object: Dict[str, MultiplexTarget], out_path: str):
    save_sgrnas_results_dict(multiplex_seq_to_object, out_path)
    save_sgrna_amplicons_dict(multiplex_amplicons_dict, out_path)


def get_gene_multiplex_amps(max_amplicon_len_category: int, primer_length: int, target_surrounding_region: int,
                            cut_location: int, annotations_file_path: str, out_path: str, genome_fasta_file: str,
                            distinct_alleles_num: int, pams: Tuple[str], target_len: int, primer3_core_path: str,
                            n: int, filter_off_targets: int):
    amplicon_ranges = [(200, 300), (300, 500), (500, 1000)]
    max_amplicon_len = max(amplicon_ranges[max_amplicon_len_category - 1])
    min_amplicon_len = min(amplicon_ranges[max_amplicon_len_category - 1])
    gene_exon_regions_seqs_dict, original_exon_indices_dict = extract_exons_regions(max_amplicon_len, primer_length,
                                                                                    target_surrounding_region,
                                                                                    cut_location,
                                                                                    annotations_file_path, out_path,
                                                                                    genome_fasta_file)
    gene_snps_dict = get_snps(gene_exon_regions_seqs_dict, distinct_alleles_num, primer_length)
    multiplex_targets_list = get_multiplex_targets(gene_exon_regions_seqs_dict, pams,
                                                   max_amplicon_len, primer_length, cut_location,
                                                   target_surrounding_region, target_len,
                                                   distinct_alleles_num)
    if filter_off_targets:
        multiplex_targets_list, off_targets_df = filter_multiplex_off_targets(multiplex_targets_list, out_path, genome_fasta_file, pams,
                                                              gene_exon_regions_seqs_dict)
    multiplex_amplicons_dict, multiplex_seq_to_object = create_multiplex_amplicons(
            multiplex_targets_list, gene_exon_regions_seqs_dict, gene_snps_dict, original_exon_indices_dict, max_amplicon_len,
            primer_length, distinct_alleles_num, target_surrounding_region, min_amplicon_len,
            target_len, out_path, primer3_core_path, n, amplicon_ranges,
            max_amplicon_len_category)
    if filter_off_targets:
        get_off_target_from_df(off_targets_df, multiplex_seq_to_object, out_path)
    else:
        get_off_targets(multiplex_seq_to_object, out_path, genome_fasta_file, pams, gene_exon_regions_seqs_dict)
    multiplex_primers_off_targets(multiplex_amplicons_dict, out_path, genome_fasta_file, max_amplicon_len, gene_exon_regions_seqs_dict, 0)
    save_results(multiplex_amplicons_dict, multiplex_seq_to_object, out_path)
    return
