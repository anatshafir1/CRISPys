
from typing import Dict, List, Tuple, Any

from pandas import DataFrame

from Amplicon_construction.FindOffTargets import run_bwa, extract_off_targets, calc_off_scores_for_df
from CRISPys_master.Candidate import Candidate
from globals import OFF_TARGET_FILTER_CUTOFF


def create_bwa_sgrna_input_fasta(sgrna_seq_object_dict: Dict[str, Candidate], out_path: str, family_targeting: int):
    grnas_fasta = out_path + "/gRNA_input.fasta"
    unique_sgrna_list = []
    for seq in sgrna_seq_object_dict:
        if family_targeting == 1:
            if seq not in unique_sgrna_list:
                unique_sgrna_list.append(seq)
        elif family_targeting == 2:
            seq1 = seq.split("_")[0]
            seq2 = seq.split("_")[1]
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


def get_family_scaffolds_positions_dict(
        genes_exons_seq_dict: Dict[str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]]) -> Dict[str, List[Tuple[int, int]]]:
    """
    create dictionary of {scaffold_id names: list of gene region indices for the scaffold_id}
    """
    from Amplicon_construction.AmpliconConstruction import get_gene_alleles_positions
    scaffolds_positions_dict = {}
    for gene in genes_exons_seq_dict:
        gene_scaffold_pos_dict = get_gene_alleles_positions(genes_exons_seq_dict[gene][0])
        for scaffold_id in gene_scaffold_pos_dict:
            if scaffold_id not in scaffolds_positions_dict:
                scaffolds_positions_dict[scaffold_id] = gene_scaffold_pos_dict[scaffold_id]
            else:
                scaffolds_positions_dict[scaffold_id].extend(gene_scaffold_pos_dict[scaffold_id])
    return scaffolds_positions_dict


def remove_on_targets_df(off_targets_df: DataFrame,
                         genes_exons_seq_dict: Dict[str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]]):
    scaffolds_positions_dict = get_family_scaffolds_positions_dict(genes_exons_seq_dict)

    def on_target(row):
        mms = row['Mismatches']
        chrom = row['Chromosome']
        position = row['Position']
        if int(mms) == 0:
            if chrom in scaffolds_positions_dict:
                for start, end in scaffolds_positions_dict[chrom]:
                    if start <= int(position) <= end:
                        return True
        return False

    rows_to_remove = off_targets_df[off_targets_df.apply(on_target, axis=1)].index
    filtered_off_targets_df = off_targets_df.drop(index=rows_to_remove, inplace=False)
    return filtered_off_targets_df


def search_sgrna_off_targets(sgrna_seq_dict: Dict[str, Any], out_path: str, genome_fasta: str,
                             pams, genes_exons_seq_dict: Dict[
            str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]], family_targeting: int):
    print("searching for sgrna off targets".upper().center(40, "#"))
    bwa_input_fasta = create_bwa_sgrna_input_fasta(sgrna_seq_dict, out_path, family_targeting)
    off_targets_sam = run_bwa(bwa_input_fasta, genome_fasta, out_path)
    off_targets_df = extract_off_targets(off_targets_sam, genome_fasta, pams)
    filtered_off_targets_df = remove_on_targets_df(off_targets_df, genes_exons_seq_dict)
    calc_off_scores_for_df(filtered_off_targets_df)
    return filtered_off_targets_df


def filter_sgrna_off_targets(sgrna_seq_dict: Dict[str, Any], sorted_sgrnas_list: List, out_path: str, genome_fasta: str,
                          pams: Tuple, genes_exons_seq_dict: Dict[str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]],
                       family_targeting: int):
    updated_sgrna_seq_dict = {}
    updated_sorted_sgrnas_list = []
    off_targets_df = search_sgrna_off_targets(sgrna_seq_dict, out_path, genome_fasta, pams, genes_exons_seq_dict, family_targeting)
    strong_off_sgrnas = off_targets_df.loc[off_targets_df["off_scores"] > OFF_TARGET_FILTER_CUTOFF, 'crRNA'].tolist()
    for sgrna in sorted_sgrnas_list:
        if family_targeting == 1:
            if sgrna.seq not in strong_off_sgrnas:
                updated_sorted_sgrnas_list.append(sgrna)
                updated_sgrna_seq_dict[sgrna.seq] = sgrna_seq_dict[sgrna.seq]
        elif family_targeting == 2:
            if sgrna[0].seq not in strong_off_sgrnas and sgrna[1].seq not in strong_off_sgrnas:
                updated_sorted_sgrnas_list.append(sgrna)
                updated_sgrna_seq_dict[f"{sgrna[0].seq}_{sgrna[1].seq}"] = sgrna_seq_dict[f"{sgrna[0].seq}_{sgrna[1].seq}"]
    return updated_sgrna_seq_dict, updated_sorted_sgrnas_list, off_targets_df


def get_off_target_from_df(off_target_df: DataFrame, sgrna_seq_to_object_dict, family_targeting: int, out_path: str):
    unique_sgrna_list = []
    for seq in sgrna_seq_to_object_dict:
        if family_targeting == 1:
            if seq not in unique_sgrna_list:
                unique_sgrna_list.append(seq)
        elif family_targeting == 2:
            seq1 = seq.split("_")[0]
            seq2 = seq.split("_")[1]
            if seq1 not in unique_sgrna_list:
                unique_sgrna_list.append(seq1)
            if seq2 not in unique_sgrna_list:
                unique_sgrna_list.append(seq2)
    filtered_off_targets_df = off_target_df[off_target_df['crRNA'].isin(unique_sgrna_list)]
    filtered_off_targets_df.to_csv(out_path + "/off_targets.csv", index=False)


def get_off_targets(sgrna_seq_dict: Dict[str, Any], out_path: str, genome_fasta: str, pams,
                    genes_exons_seq_dict: Dict[str, Tuple[Dict[int, List[Tuple[str, str]]], Dict[str, Dict[int, int]]]],
                    family_targeting: int):
    off_targets_df = search_sgrna_off_targets(sgrna_seq_dict, out_path, genome_fasta, pams, genes_exons_seq_dict, family_targeting)
    filtered_df = off_targets_df.groupby('crRNA', group_keys=False).apply(lambda x: x.nlargest(2, 'off_scores'))
    filtered_df.to_csv(out_path + "/off_targets.csv", index=False)
