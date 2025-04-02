from typing import Dict, List, Tuple

from Amplicon_construction.AmpliconConstruction import get_candidates_scaffold_positions
from Amplicon_construction.FindOffTargets import run_bwa, extract_off_targets, remove_on_targets_df
from Amplicon_construction.Target_Obj import MultiplexTarget
from Gene_Family_Targeting.FindSGRNAOffTargets import calc_off_scores_for_df
from globals import OFF_TARGET_FILTER_CUTOFF


def create_bwa_multiplex_input_fasta(multiplex_targets_list: List[MultiplexTarget], out_path: str):
    grnas_fasta = out_path + "/multiplex_targets_input.fasta"
    unique_sgrna_list = []
    for multiplex_target in multiplex_targets_list:
        seq1 = multiplex_target.up_seq
        seq2 = multiplex_target.down_seq
        if seq1 not in unique_sgrna_list:
            unique_sgrna_list.append(seq1)
        if seq2 not in unique_sgrna_list:
            unique_sgrna_list.append(seq2)
    out_str = ""
    for grna in unique_sgrna_list:
        out_str += f">{grna}\n{grna}\n"
    with open(grnas_fasta, 'w') as f:
        f.write(out_str)
    return grnas_fasta


def search_multiplex_off_targets(multiplex_targets_list: List[MultiplexTarget], out_path: str, genome_fasta_file: str,
                                 pams: Tuple[str], gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]]):
    print("searching for multiplex off targets".upper().center(40, "#"))
    bwa_input_fasta = create_bwa_multiplex_input_fasta(multiplex_targets_list, out_path)
    off_targets_sam = run_bwa(bwa_input_fasta, genome_fasta_file, out_path)
    off_targets_df = extract_off_targets(off_targets_sam, genome_fasta_file, pams)
    candidates_scaffold_positions = get_candidates_scaffold_positions(gene_exon_regions_seqs_dict)
    filtered_off_targets_df = remove_on_targets_df(off_targets_df, candidates_scaffold_positions)
    calc_off_scores_for_df(filtered_off_targets_df)
    return filtered_off_targets_df


def filter_multiplex_off_targets(multiplex_targets_list: List[MultiplexTarget],
                                 out_path: str, genome_fasta_file: str, pams: Tuple[str],
                                 gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]]) -> List[MultiplexTarget]:

    updated_multiplex_targets_list = []
    off_targets_df = search_multiplex_off_targets(multiplex_targets_list, out_path, genome_fasta_file, pams,
                                           gene_exon_regions_seqs_dict)
    strong_off_sgrnas = off_targets_df.loc[off_targets_df["off_scores"] > OFF_TARGET_FILTER_CUTOFF, 'crRNA'].tolist()
    for multiplex_target in multiplex_targets_list:
        if multiplex_target.up_seq not in strong_off_sgrnas and multiplex_target.down_seq not in strong_off_sgrnas:
            updated_multiplex_targets_list.append(multiplex_target)
    return updated_multiplex_targets_list
