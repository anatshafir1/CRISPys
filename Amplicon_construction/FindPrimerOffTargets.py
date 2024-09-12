"""Module for finding Primers off-targets"""
import os
import re
import pandas as pd
from typing import Dict, List, Tuple

import subprocess
from Bio import SeqIO
from pandas import DataFrame

from Amplicon_Obj import Amplicon_Obj

from globals import OFF_PRIMER_CUTOFF


def create_bwa_input(candidate_amplicons_list: List[Amplicon_Obj], off_primers_fasta: str) -> str:
    """
    :param candidate_amplicons_list: A list of Amplicon_Obj objects
    :param off_primers_fasta: A path to the algorithm result folder where a folder for crispritz will be created
    :return: A path to the input for xxx and will write the input file to xxx
    """
    out = ""
    unique_primers = []

    for candidate in candidate_amplicons_list:  # go over each candidate and get the guide sequence
        left_primer_seq = candidate.primers.left_sequence  # get left primer sequence
        right_primer_seq = candidate.primers.right_sequence  # get right primer sequence
        if left_primer_seq not in unique_primers:
            unique_primers.append(left_primer_seq)
        if right_primer_seq not in unique_primers:
            unique_primers.append(right_primer_seq)
    for primer in unique_primers:
        out += f">{primer}\n{primer}\n"
    with open(off_primers_fasta, 'w') as f:
        f.write(out)
    return off_primers_fasta


def check_bwa_index_files(genome_fasta: str) -> bool:
    """
    Check if the BWA index files for the given genome FASTA file exist.

    Parameters:
    genome_fasta (str): Path to the genome FASTA file.

    Returns:
    bool: True if all BWA index files exist, False otherwise.
    """
    index_files = [genome_fasta + ext for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']]

    return all(os.path.isfile(f) for f in index_files)


def index_genome(genome_fasta: str):
    """
    Index the genome FASTA file using BWA.

    Parameters:
    genome_fasta (str): Path to the genome FASTA file.
    """
    try:
        subprocess.run(['bwa', 'index', genome_fasta], check=True)
        print(f"Indexing of {genome_fasta} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while indexing the genome: {e}")


def run_bwa(candidate_amplicons_list: List[Amplicon_Obj], genome_fasta: str, out_path: str) -> str:
    """run off-target search with BWA and return path to result SAM file

    :param candidate_amplicons_list:
    :param genome_fasta:
    :param out_path:
    :return:
    """
    print("Searching for primers off-targets with BWA")
    primers_input_fasta_path = out_path + "/primers_input.fasta"
    # create a primers input file for search
    primers_fasta = create_bwa_input(candidate_amplicons_list, primers_input_fasta_path)
    # check if the genome is indexed. index if not
    if check_bwa_index_files(genome_fasta):
        print(f"BWA index files for {genome_fasta} already exist.")
    else:
        print(f"BWA index files for {genome_fasta} do not exist. Indexing the genome...")
        index_genome(genome_fasta)

    # run bwa alignment for primers off-target search
    alignment_sai = out_path + "/primers_alignment.sai"
    bwa_aln_cmd = ["bwa", "aln", "-N", "-n", "2", "-o", "0", "-l", "20", "-k", "4", genome_fasta, primers_fasta]
    with open(alignment_sai, "w") as out:
        result = subprocess.run(bwa_aln_cmd, stdout=out, stderr=subprocess.PIPE)
        # Check for errors
        if result.returncode != 0:
            print(f"Error: {result.stderr.decode('utf-8')}")
        else:
            print("alignment executed successfully.")

    # create, sort and index alignment SAM file
    alignment_sam = out_path + "/primers_alignment.sam"
    bwa_samse_cmd = ['bwa', 'samse', '-n', '10000000', '-f', alignment_sam, genome_fasta, alignment_sai, primers_fasta]
    with open(alignment_sam, "w") as out:
        result = subprocess.run(bwa_samse_cmd, stdout=out, stderr=subprocess.PIPE)
        # Check for errors
        if result.returncode != 0:
            print(f"Error: {result.stderr.decode('utf-8')}")
        else:
            print("SAM created successfully.")

    return alignment_sam


def extract_off_targets(sam_file: str, genome_fasta: str, out_path: str) -> DataFrame:
    """

    :param sam_file:
    :param genome_fasta:
    :param out_path: path to output directory where algorithm results will be saved
    :return:
    """
    off_targets = {'orig_primer': [], 'off_primer': [], 'Chromosome': [], 'Position': [], 'Direction': [], 'Mismatches': []}
    reference_genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue  # Skip header lines

            columns = line.split('\t')
            orig_primer = columns[0]

            # handle first match:
            scaffold = columns[2]
            position = columns[3]
            cigar = columns[5]
            strand = "-" if int(columns[1]) & 0x10 != 0 else "+"   # inferring strand from the FLAG's 5th bit
            start = int(position) - 1  # Convert to 0-based indexing
            ref_seq = reference_genome[scaffold].seq
            off_target_len = sum(int(length) for length, op in re.findall(r'(\d+)([MID])', cigar) if op != 'D')
            if not off_target_len < len(orig_primer):
                end = start + off_target_len
                # Handle strand orientation
                if strand == "-":
                    sequence = ref_seq[start:end]
                    sequence = sequence.reverse_complement()
                else:
                    sequence = ref_seq[start:end]
                off_targets['orig_primer'].append(orig_primer)
                off_targets['off_primer'].append(str(sequence))
                off_targets['Chromosome'].append(scaffold)
                off_targets['Position'].append(position)
                off_targets['Direction'].append(strand)
                off_targets['Mismatches'].append(columns[12].split(":")[-1])

            # handle rest of the matches:
            xa_tag = next((col for col in columns if col.startswith('XA:Z:')), None)  # skip to the column where the off-targets are
            if xa_tag:
                # Extract the off-target alignments from the XA:Z: tag
                off_target_entries = xa_tag[5:].split(';')[:-1]  # Remove the trailing empty entry

                for entry in off_target_entries:
                    scaffold, strand_position, cigar, mismatch = entry.split(',')
                    strand = strand_position[0]
                    position = strand_position[1:]
                    # Extract sequence from the reference genome
                    ref_seq = reference_genome[scaffold].seq
                    start = int(position) - 1  # Convert to 0-based indexing
                    off_target_len = sum(int(length) for length, op in re.findall(r'(\d+)([MID])', cigar) if op != 'D')
                    if not off_target_len < len(orig_primer):
                        end = start + off_target_len
                        # Handle strand orientation
                        if strand == "-":
                            sequence = ref_seq[start:end]
                            sequence = sequence.reverse_complement()
                        else:
                            sequence = ref_seq[start:end]
                        off_targets['orig_primer'].append(orig_primer)
                        off_targets['off_primer'].append(str(sequence))
                        off_targets['Chromosome'].append(scaffold)
                        off_targets['Position'].append(position)
                        off_targets['Direction'].append(strand)
                        off_targets['Mismatches'].append(mismatch)

        #  create a dataframe from the dictionary
        off_targets_df = pd.DataFrame.from_dict(off_targets)

    off_targets_df.to_csv(out_path+"/off_primers.csv", index=False)
    return off_targets_df


def filter_off_targets(off_targets_df: DataFrame, candidates_scaffold_positions: Dict[str, Tuple[int, int]],
                       out_path: str) -> DataFrame:
    """
    Filter the original primers from the primers off-target dataframe

    :param off_targets_df: DataFrame of found off-targets for the primers, including the primers themselves
    :param candidates_scaffold_positions: Dictionary of scaffold IDs to current gene indices in the scaffold
    :param out_path:
    :return:
    """
    scaffold_lst = list(candidates_scaffold_positions.keys())
    # Convert 'Position' and 'Mismatches' columns to integer
    off_targets_df['Position'] = off_targets_df['Position'].astype(int)
    off_targets_df['Mismatches'] = off_targets_df['Mismatches'].astype(int)

    no_mismatches = off_targets_df['Mismatches'] == 0
    in_scaffold = off_targets_df['Chromosome'].isin(scaffold_lst)
    # Check positional condition in the scaffold
    positional_condition = off_targets_df.apply(
        lambda row: candidates_scaffold_positions[row['Chromosome']][0] < row['Position'] <
                    candidates_scaffold_positions[row['Chromosome']][1]
        if row['Chromosome'] in candidates_scaffold_positions else False, axis=1)

    filter_condition = no_mismatches & in_scaffold & positional_condition
    filtered_off_targets_df = off_targets_df[~filter_condition]

    for index, row in filtered_off_targets_df.iterrows():
        if row["Chromosome"] in scaffold_lst:  # check if any of the filtered primers is on one of the gene's scaffolds (alleles)
            rows_with_scaffold = off_targets_df[(off_targets_df["Chromosome"] == row["Chromosome"]) & (off_targets_df["Direction"] != row["Direction"])]  # filter rows with current scaffold and opposite direction
            rows_to_add = rows_with_scaffold[~rows_with_scaffold.isin(filtered_off_targets_df).all(axis=1)]  # filter out rows from off_targets_df that are already in filtered_off_targets_df
            filtered_off_targets_df = filtered_off_targets_df.append(rows_to_add)

    return filtered_off_targets_df


def get_problematic_primers(filtered_off_targets_pd: DataFrame, max_amplicon_len: int):
    """

    :param filtered_off_targets_pd:
    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :return:
    """

    problematic_primers = {'off_primers_scaffold': [], 'orig_rev_primer_seq': [], 'off_rev_primer_seq': [],
                           'off_rev_primer_pos': [], 'off_rev_primer_mms': [], 'orig_fwd_primer_seq': [],
                           'off_fwd_primer_seq': [], 'off_fwd_primer_pos': [], 'off_fwd_primer_mms': [], 'distance': []}
    scaffold_groups = filtered_off_targets_pd.groupby('Chromosome')
    for _, group in scaffold_groups:
        revs_df = group[group['Direction'] == "-"].reset_index()
        revs = revs_df.values.tolist()
        fwds_df = group[group['Direction'] == "+"].reset_index()
        fwds = fwds_df.values.tolist()
        for i in range(len(revs)):
            for j in range(len(fwds)):
                distance = abs(revs[i][4] - fwds[j][4])
                if distance <= max_amplicon_len + OFF_PRIMER_CUTOFF:
                    problematic_primers['off_primers_scaffold'].append(revs[i][3])
                    problematic_primers['orig_rev_primer_seq'].append(revs[i][1])
                    problematic_primers['off_rev_primer_seq'].append(revs[i][2])
                    problematic_primers['off_rev_primer_pos'].append(revs[i][4])
                    problematic_primers['off_rev_primer_mms'].append(revs[i][6])
                    problematic_primers['orig_fwd_primer_seq'].append(fwds[j][1])
                    problematic_primers['off_fwd_primer_seq'].append(fwds[j][2])
                    problematic_primers['off_fwd_primer_pos'].append(fwds[j][4])
                    problematic_primers['off_fwd_primer_mms'].append(fwds[j][6])
                    problematic_primers['distance'].append(distance)

    problematic_primers_df = pd.DataFrame.from_dict(problematic_primers)
    return problematic_primers_df


def get_primers_off_targets(candidate_amplicons_list: List[Amplicon_Obj], genome_fasta_file: str, out_path: str,
                    candidates_scaffold_positions: Dict[str, Tuple[int, int]], max_amplicon_len: int):
    """
    Find the off-targets for each candidate using BWA, store them in the candidate's off_targets_list attribute
    as a list of OffTarget objects. Then calculates the off-target scores for each off-target of each candidate and
    store the scores in each off-target's score attribute.

    :param candidate_amplicons_list: A list of amplicon candidates
    :param out_path: output path for the algorithm results
    :param genome_fasta_file: path to input FASTA format file of the genome
    :param candidates_scaffold_positions: dictionary of allele scaffold -> gene allele start,end indices
    :param max_amplicon_len: maximum length of the amplicon, defined by user
    :return
    """

    # run off target search
    off_targets_sam = run_bwa(candidate_amplicons_list, genome_fasta_file, out_path)
    # extract off-targets from SAM file to pandas DataFrame
    off_targets_pd = extract_off_targets(off_targets_sam, genome_fasta_file, out_path)
    # filter on-targets
    filtered_off_targets_pd = filter_off_targets(off_targets_pd, candidates_scaffold_positions, out_path)
    # find problematic primers
    problematic_primers_df = get_problematic_primers(filtered_off_targets_pd, max_amplicon_len)
    problematic_primers_df.to_csv(out_path+"/problematic_primers.csv", index=False)
    return problematic_primers_df
