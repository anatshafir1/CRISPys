import pickle
from typing import Tuple, Dict, Any, List

import pandas as pd

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.FindPrimerOffTargets import gene_family_primers_off_targets
from CRISPys_master.Stage0 import CRISPys_main
from Gene_Family_Targeting.CreatesgRNAamplicons import create_single_sgrna_amplicons, create_multiplex_sgrna_amplicons
from Gene_Family_Targeting.FindSGRNAOffTargets import filter_sgrna_off_targets, get_off_targets, get_off_target_from_df
from Gene_Family_Targeting.GetGeneFamilySNPs import get_snps_dict
from Gene_Family_Targeting.GetGeneFamilySequences import create_crispys_input_fasta, get_sequences_dict
from Gene_Family_Targeting.GetGeneFamilyTargets import get_sgrnas_from_crispys_output
from globals import REGION_OF_GENE_TO_CUT, MAX_POLYMORPHIC_SITES


# Parameters for CRISPys
output_name = "crispys_out"
genes_of_interest_file = 'None'
alg = "default"
where_in_gene = REGION_OF_GENE_TO_CUT
omega = 0.5
off_scoring_function = "moff"
on_scoring_function = "default"
start_with_g = 0
internal_node_candidates = 10
max_target_polymorphic_sites = MAX_POLYMORPHIC_SITES
crispys_pams = 0
slim_output = 0
set_cover = 0
min_desired_genes_fraction = -1.0
singletons = 0
singletons_on_target_function = "ucrispr"
number_of_singletons = 50
max_gap_distance = 3
export_tree = 0
run4chips = 0


def save_sgrnas_results_dict(sgrna_seq_to_object_dict: Dict[str, Any], sgrnas_dict: Dict[str, Tuple], family_targeting: int, out_path: str):
    sgrnas_results_dict = {}
    if family_targeting == 1:
        sgrnas_results_dict = {"vector_ID": [], "multiplex_score": [], "sgrna_seq": [], "cut_genes": []}
    elif family_targeting == 2:
        sgrnas_results_dict = {"vector_ID": [], "multiplex_score": [], "sgrna_seq_1": [], "cut_genes_sgrna_1": [],
                               "sgrna_seq_2": [], "cut_genes_sgrna_2": []}
    for sgrna in sgrna_seq_to_object_dict:
        if family_targeting == 1:
            sgrnas_results_dict["vector_ID"].append(sgrna)
            sgrnas_results_dict["multiplex_score"].append(sgrnas_dict[sgrna][2])
            sgrnas_results_dict["sgrna_seq"].append(sgrna)
            sgrnas_results_dict["cut_genes"].append(sgrnas_dict[sgrna][4][sgrna])
        elif family_targeting == 2:
            seq1 = sgrna.split("_")[0]
            seq2 = sgrna.split("_")[1]
            sgrnas_results_dict["vector_ID"].append(sgrna)
            sgrnas_results_dict["multiplex_score"].append(sgrnas_dict[sgrna][2])
            sgrnas_results_dict["sgrna_seq_1"].append(seq1)
            sgrnas_results_dict["cut_genes_sgrna_1"].append(sgrnas_dict[sgrna][4][seq1])
            sgrnas_results_dict["sgrna_seq_2"].append(seq2)
            sgrnas_results_dict["cut_genes_sgrna_2"].append(sgrnas_dict[sgrna][4][seq2])
    df = pd.DataFrame(sgrnas_results_dict)
    df.to_csv(out_path + "/sgrna_results.csv", index=False)


def save_sgrna_amplicons_dict(sgrna_amplicons_dict: Dict[str, Dict[str, List[Amplicon_Obj]]], out_path: str, family_targeting: int):

    records = []
    for sgrna_id, gene_dict in sgrna_amplicons_dict.items():
        for gene_name, amplicons_list in gene_dict.items():
            for amp in amplicons_list:
                for scaffold_amp in amp.scaffold_amplicons:
                    params_dict = {
                        "vector_ID": sgrna_id,
                        "gene_name": gene_name}
                    params_dict.update(amp.scaffold_amplicons[scaffold_amp].to_dict(amp.target.rank, 0, 0, family_targeting))
                    records.append(params_dict)
    df = pd.DataFrame(records)
    df.to_csv(out_path + "/sgrna_amplicons_results.csv", index=False)


def save_results(sgrna_amplicons_dict: Dict[str, Dict[str, List[Amplicon_Obj]]], sgrnas_dict: Dict[str, Tuple],
                 sgrna_seq_to_object_dict: Dict[str, Any], family_targeting: int, out_path: str):
    save_sgrnas_results_dict(sgrna_seq_to_object_dict, sgrnas_dict, family_targeting, out_path)
    save_sgrna_amplicons_dict(sgrna_amplicons_dict, out_path, family_targeting)


def get_gene_family_amps(max_amplicon_len_category: int, primer_length: int, target_surrounding_region: int, cut_location: int,
                  annotations_file_path: str, out_path: str, genome_fasta_file: str, distinct_alleles_num: int,
                  pams: Tuple[str], target_len: int, primer3_core_path: str, n: int, filter_off_targets: int,
                  family_targeting: int):
    # crispys_input_fasta = create_crispys_input_fasta(annotations_file_path, out_path, genome_fasta_file)
    # res_grnas = CRISPys_main(crispys_input_fasta, out_path, output_name, genes_of_interest_file, alg, where_in_gene, omega,
    #                    off_scoring_function, on_scoring_function, start_with_g, internal_node_candidates,
    #                    max_target_polymorphic_sites, crispys_pams, slim_output, set_cover, min_desired_genes_fraction,
    #                    singletons, singletons_on_target_function, number_of_singletons, max_gap_distance, export_tree,
    #                    run4chips)
    with open("/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/output/test/tool4_test/crispys_out.p",
              "rb") as file:
        res_grnas = pickle.load(file)
    amplicon_ranges = [(200, 300), (300, 500), (500, 1000)]
    max_amplicon_len = max(amplicon_ranges[max_amplicon_len_category - 1])
    min_amplicon_len = min(amplicon_ranges[max_amplicon_len_category - 1])
    genes_exons_seq_dict = get_sequences_dict(max_amplicon_len, primer_length, target_surrounding_region, cut_location, annotations_file_path,
                                              out_path, genome_fasta_file)
    genes_snps_dict = get_snps_dict(genes_exons_seq_dict, distinct_alleles_num, primer_length)
    sgrnas_dict, sorted_sgrnas_list = get_sgrnas_from_crispys_output(res_grnas[0], genes_exons_seq_dict, family_targeting)
    if filter_off_targets:
        sgrnas_dict, sorted_sgrnas_list, off_targets_df = filter_sgrna_off_targets(sgrnas_dict, sorted_sgrnas_list, out_path, genome_fasta_file, pams, genes_exons_seq_dict, family_targeting)
    if family_targeting == 1:
        sgrna_amplicons_dict, sgrna_seq_to_object_dict = create_single_sgrna_amplicons(sorted_sgrnas_list, genes_exons_seq_dict, genes_snps_dict, max_amplicon_len,
                                                  primer_length, distinct_alleles_num, target_surrounding_region, min_amplicon_len,
                                                  target_len, out_path, primer3_core_path, n, amplicon_ranges,
                                                  max_amplicon_len_category, cut_location, family_targeting)

    else:  # family_targeting == 2:
        sgrna_amplicons_dict, sgrna_seq_to_object_dict = create_multiplex_sgrna_amplicons(
            sorted_sgrnas_list, genes_exons_seq_dict, genes_snps_dict, max_amplicon_len,
            primer_length, distinct_alleles_num, target_surrounding_region, min_amplicon_len,
            target_len, out_path, primer3_core_path, n, amplicon_ranges,
            max_amplicon_len_category, cut_location, family_targeting)
    if filter_off_targets:
        get_off_target_from_df(off_targets_df, sgrna_seq_to_object_dict, family_targeting, out_path)
    else:
        get_off_targets(sgrna_seq_to_object_dict, out_path, genome_fasta_file, pams, genes_exons_seq_dict, family_targeting)
    gene_family_primers_off_targets(sgrna_amplicons_dict, out_path, genome_fasta_file, max_amplicon_len, genes_exons_seq_dict, family_targeting)
    save_results(sgrna_amplicons_dict, sgrnas_dict, sgrna_seq_to_object_dict, family_targeting, out_path)
    return
