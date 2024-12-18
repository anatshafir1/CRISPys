"""Module for finding off-targets"""
import os
import time
import re
import pandas as pd
from typing import Dict, List, Tuple

import subprocess
from Bio import SeqIO
from pandas import DataFrame

from Amplicon_Obj import Amplicon_Obj, OffTarget

from MOFF.MOFF_prediction import MOFF_score
from MOFF.MoffLoad import mtx1, mtx2, model


# create dictionary of sequence:candidate
def create_sequence_to_candidate_dict(candidate_amplicons_list: List[Amplicon_Obj], k: int, multiplex: int,
                                      sgrna: str = None) -> Dict[str, List[Amplicon_Obj]]:
    """
    This function takes a list of sgRNA candidates, and returns a sequence to candidate dictionary.

    :param candidate_amplicons_list: a list of sgRNA candidates
    :param k: number of alleles to target with a single gRNA.
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param sgrna: define which sgRNA to use in case of multiplex.
    :return: a dictionary: sequence -> an ActivationCandidate object where candidate.seq = sequence.
    """
    sequence_to_candidate_dict = {}
    up_sequence_to_candidate_dict = {}
    down_sequence_to_candidate_dict = {}
    if multiplex:
        for i, candidate in enumerate(candidate_amplicons_list):
            if sgrna == "up":
                target_seq = candidate.target.up_seq
                if target_seq in sequence_to_candidate_dict:
                    up_sequence_to_candidate_dict[target_seq] += [candidate_amplicons_list[i]]
                else:
                    up_sequence_to_candidate_dict[target_seq] = [candidate_amplicons_list[i]]
            elif sgrna == "down":
                target_seq = candidate.target.down_seq
                if target_seq in sequence_to_candidate_dict:
                    down_sequence_to_candidate_dict[target_seq] += [candidate_amplicons_list[i]]
                else:
                    down_sequence_to_candidate_dict[target_seq] = [candidate_amplicons_list[i]]
        if sgrna == "up":
            return up_sequence_to_candidate_dict
        elif sgrna == "down":
            return down_sequence_to_candidate_dict
    else:
        for i, candidate in enumerate(candidate_amplicons_list):
            if k > 0:
                target_seq = candidate.target.chosen_sg
            else:
                target_seq = candidate.target.seq[:-3]
            if target_seq in sequence_to_candidate_dict:
                sequence_to_candidate_dict[target_seq] += [candidate_amplicons_list[i]]
            else:
                sequence_to_candidate_dict[target_seq] = [candidate_amplicons_list[i]]
        return sequence_to_candidate_dict


# add to each "Candidate" its off-targets
def add_off_targets(off_targets_df, sequence_to_candidate_dict: Dict[str, List[Amplicon_Obj]], multiplex: int,
                    sgrna: str = None):
    """
    This function adds all found off-targets to each CandidateWithOffTargets using the crispritz results.

    :param off_targets_df: The output of crispritz as a pd datatable, where each row is a potential offtarget.
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence.
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param sgrna: define which sgRNA to use in case of multiplex.

    """
    # apply the 'get_off_target' function on each row in the crispritz table results
    off_targets_df.apply(get_off_target, args=(sequence_to_candidate_dict, multiplex, sgrna), axis=1)
    return


# This function will be used in an apply command. it reads crispritz results and create an offtarget object out of each one
def get_off_target(x, sequence_to_candidate_dict: Dict[str, List[Amplicon_Obj]], multiplex, sgrna):
    """
    A function to use with apply on crispritz result table
    it takes a row of crispritz results and a dictionary of sequence:candidate, and make an OffTarget
    object from the crispritz results and add it to the candidate offtargets list.

    :param x: a row in crispritz results table
    :param sequence_to_candidate_dict: a dictionary of sequence:candidate
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param sgrna: define which sgRNA to use in case of multiplex.
    """
    candidates_list = sequence_to_candidate_dict[x['crRNA'][:20]]  # get list of candidate amplicons with current on-target (chosen_sg)
    off_target = OffTarget(x['DNA'].upper(), x['Chromosome'].split(" ")[0], int(x['Position']), x['Direction'],
                           int(x['Mismatches']))  # create off-target object
    legit_letters = True
    for char in off_target.seq:  # allow only nucleotide letters in off-target sequence
        if char not in {"A", "C", "T", "G"}:
            legit_letters = False
            break
    if legit_letters:  # add current off-target to all amplicons with matching on target (candidates_list)
        for candidate in candidates_list:
            if multiplex:
                if sgrna == "up":
                    if off_target not in candidate.up_off_targets:
                        candidate.up_off_targets.append(off_target)
                elif sgrna == "down":
                    if off_target not in candidate.down_off_targets:
                        candidate.down_off_targets.append(off_target)
            else:
                if off_target not in candidate.off_targets:
                    candidate.off_targets.append(off_target)
    return


def moff(candidate_lst: List[str], target_lst: List[str]) -> List[float]:
    """
    Calling MOFF algorithm, this function take list of sgrnas (candidate) and list of targets and
    returns a list of MOFF score

    :param candidate_lst: list of candidates
    :param target_lst: list of targets
    :return: list of MOFF score
    """
    scores = MOFF_score(mtx1, mtx2, model, candidate_lst, target_lst)
    return list(scores)


def calculate_scores(candidate_amplicons_list: List[Amplicon_Obj], k: int, multiplex: int, sgrna: str = None):
    """
    Calculate the off-target scores for each off-target of each candidate and store the scores in each off-target's
    score attribute.

    :param candidate_amplicons_list: a list of Amplicon candidates
    :param k: number of alleles to target with a single gRNA
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param sgrna: define which sgRNA to use in case of multiplex.
    """
    if multiplex:
        batch_off_targets_list = []
        batch_candidates_list = []
        if sgrna == "up":
            for candidate in candidate_amplicons_list:
                batch_candidates_list += [candidate.target.up_seq for _ in range(len(candidate.up_off_targets))]
                batch_off_targets_list += [up_off_target.seq for up_off_target in candidate.up_off_targets]
        elif sgrna == "down":
            for candidate in candidate_amplicons_list:
                batch_candidates_list += [candidate.target.down_seq for _ in range(len(candidate.down_off_targets))]
                batch_off_targets_list += [down_off_target.seq for down_off_target in candidate.down_off_targets]
        t0 = time.perf_counter()
        scores = moff(batch_candidates_list, batch_off_targets_list)
        t1 = time.perf_counter()
        print(f"scoring function for {len(batch_candidates_list)} off-targets ran in {t1 - t0} seconds")
        i = 0
        for candidate in candidate_amplicons_list:
            if sgrna == "up":
                if len(candidate.up_off_targets) > 0:
                    for up_off in candidate.up_off_targets:
                        up_off.score = round(scores[i], 4)
                        i += 1
                    # sort the off-targets in the off_targets_list by score from highest to lowest
                    candidate.sort_off_targets()
            elif sgrna == "down":
                if len(candidate.down_off_targets) > 0:
                    for down_off in candidate.down_off_targets:
                        down_off.score = round(scores[i], 4)
                        i += 1
                    # sort the off-targets in the off_targets_list by score from highest to lowest
                    candidate.sort_off_targets()
    else:
        batch_off_targets_list = []
        batch_candidates_list = []
        for candidate in candidate_amplicons_list:
            batch_candidates_list += [candidate.target.chosen_sg if k > 0 else candidate.target.seq for _ in range(len(candidate.off_targets))]
            batch_off_targets_list += [off_target.seq for off_target in candidate.off_targets]
        t0 = time.perf_counter()
        scores = moff(batch_candidates_list, batch_off_targets_list)
        t1 = time.perf_counter()
        print(f"scoring function for {len(batch_candidates_list)} off-targets ran in {t1 - t0} seconds")
        i = 0
        for candidate in candidate_amplicons_list:
            if len(candidate.off_targets) > 0:
                for off in candidate.off_targets:
                    off.score = round(scores[i], 4)
                    i += 1
                # sort the off-targets in the off_targets_list by score from highest to lowest
                candidate.sort_off_targets()


def create_bwa_input(candidate_amplicons_list: List[Amplicon_Obj], grnas_fasta: str, k: int, multiplex: int,
                     sgrna: str) -> str:
    """
    create input fasta file for BWA for all the gRNAs of the candidate amplicons

    :param candidate_amplicons_list: A list of candidate amplicons objects
    :param grnas_fasta: output path where the BWA input gRNAs fasta will be created
    :param k: number of alleles to target with a single gRNA
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param sgrna: define which sgRNA to use in case of multiplex.
    :return: path of the BWA input gRNAs fasta
    """
    out = ""
    unique_grnas = []

    for candidate in candidate_amplicons_list:  # go over each candidate and get the guide sequence
        if multiplex:
            if sgrna == "up":
                grna_seq_no_pam = candidate.target.up_seq
                if grna_seq_no_pam not in unique_grnas:
                    unique_grnas.append(grna_seq_no_pam)
            elif sgrna == "down":
                grna_seq_no_pam = candidate.target.down_seq
                if grna_seq_no_pam not in unique_grnas:
                    unique_grnas.append(grna_seq_no_pam)
            continue
        if k > 0:  # Tool 2 in use
            grna_seq_no_pam = candidate.target.chosen_sg
        else:
            grna_seq_no_pam = candidate.target.seq[:-3]  # get gRNA target sequence
        if grna_seq_no_pam not in unique_grnas:
            unique_grnas.append(grna_seq_no_pam)
    for grna in unique_grnas:
        out += f">{grna}\n{grna}\n"
    with open(grnas_fasta, 'w') as f:
        f.write(out)
    return grnas_fasta


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


def run_bwa(candidate_amplicons_list: List[Amplicon_Obj], genome_fasta: str, out_path: str, k: int, multiplex: int,
            sgrna: str = None) -> str:
    """run off-target search with BWA and return path to result SAM file

    :param candidate_amplicons_list:
    :param genome_fasta:
    :param out_path:
    :param k: number of alleles to target with a single gRNA
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param sgrna: define which sgRNA to use in case of multiplex.
    :return:
    """
    print("Searching for gRNA off-targets with BWA")
    grna_input_fasta_path = out_path + "/gRNA_input.fasta"
    # create a gRNA input file for search
    grnas_fasta = create_bwa_input(candidate_amplicons_list, grna_input_fasta_path, k, multiplex, sgrna)
    # check if the genome is indexed. index if not
    if check_bwa_index_files(genome_fasta):
        print(f"BWA index files for {genome_fasta} already exist.")
    else:
        print(f"BWA index files for {genome_fasta} do not exist. Indexing the genome...")
        index_genome(genome_fasta)

    # run bwa alignment for gNRA off-target search
    alignment_sai = out_path + "/alignment.sai"
    bwa_aln_cmd = ["bwa", "aln", "-N", "-n", "4", "-o", "0", "-l", "20", "-k", "4", genome_fasta, grnas_fasta]
    with open(alignment_sai, "w") as out:
        result = subprocess.run(bwa_aln_cmd, stdout=out, stderr=subprocess.PIPE)
        # Check for errors
        if result.returncode != 0:
            print(f"Error: {result.stderr.decode('utf-8')}")
        else:
            print("alignment executed successfully.")

    # create, sort and index alignment SAM file
    alignment_sam = out_path + "/alignment.sam"
    bwa_samse_cmd = ['bwa', 'samse', '-n', '10000000', '-f', alignment_sam, genome_fasta, alignment_sai, grnas_fasta]
    with open(alignment_sam, "w") as out:
        result = subprocess.run(bwa_samse_cmd, stdout=out, stderr=subprocess.PIPE)
        # Check for errors
        if result.returncode != 0:
            print(f"Error: {result.stderr.decode('utf-8')}")
        else:
            print("SAM created successfully.")

    return alignment_sam


def extract_off_targets(sam_file: str, genome_fasta: str, pams: Tuple) -> DataFrame:
    """

    :param sam_file:
    :param genome_fasta:
    :param pams:
    :return:
    """
    off_targets = {'crRNA': [], 'DNA': [], 'Chromosome': [], 'Position': [], 'Direction': [], 'Mismatches': []}
    reference_genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue  # Skip header lines

            columns = line.split('\t')
            grna = columns[0]

            # handle first match:
            scaffold = columns[2]
            position = columns[3]
            cigar = columns[5]
            strand = "-" if int(columns[1]) & 0x10 != 0 else "+"  # inferring strand from the FLAG's 5th bit
            start = int(position) - 1  # Convert to 0-based indexing
            ref_seq = reference_genome[scaffold].seq  # sequence of scaffold where the off-target was found
            off_target_len = sum(int(length) for length, op in re.findall(r'(\d+)([MID])', cigar) if op != 'D')
            if not off_target_len < len(grna):
                end = start + off_target_len
                if strand == "-":
                    sequence = ref_seq[start - 3:end]
                    sequence = sequence.reverse_complement()
                else:
                    sequence = ref_seq[start:end + 3]
                # check if PAM exists - otherwise it's not an off-target:
                str_seq = str(sequence)
                if len(str_seq) == 23:
                    if sequence[-3:] in pams:
                        off_targets['crRNA'].append(grna)
                        off_targets['DNA'].append(str(sequence))
                        off_targets['Chromosome'].append(scaffold)
                        off_targets['Position'].append(position)
                        off_targets['Direction'].append(strand)
                        off_targets['Mismatches'].append(columns[12].split(":")[-1])

            # handle rest of the matches:
            xa_tag = next((col for col in columns if col.startswith('XA:Z:')),
                          None)  # skip to the column where the off-targets are
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
                    if not off_target_len < len(grna):
                        end = start + off_target_len
                        # Handle strand orientation
                        if strand == "-":
                            sequence = ref_seq[start - 3:end]
                            sequence = sequence.reverse_complement()
                        else:
                            sequence = ref_seq[start:end + 3]
                        # check if PAM exist - otherwise it's not an off-target:
                        str_seq = str(sequence)
                        if len(str_seq) == 23:
                            if sequence[-3:] in pams:
                                off_targets['crRNA'].append(grna)
                                off_targets['DNA'].append(str(sequence))
                                off_targets['Chromosome'].append(scaffold)
                                off_targets['Position'].append(position)
                                off_targets['Direction'].append(strand)
                                off_targets['Mismatches'].append(mismatch)

        #  create a dataframe from the dictionary
        off_targets_df = pd.DataFrame.from_dict(off_targets)

    return off_targets_df


def remove_on_targets(candidate_amplicons_list: List[Amplicon_Obj],
                      candidates_scaffold_positions: Dict[str, Tuple[int, int]], multiplex: int):
    for candidate_amplicon in candidate_amplicons_list:
        if multiplex:
            new_up_off_targets_lst = []
            new_down_off_targets_lst = []
            for up_off in candidate_amplicon.up_off_targets:
                if up_off.number_of_mismatches == 0:
                    off_target_scaffold = up_off.chromosome
                    if off_target_scaffold in candidates_scaffold_positions:
                        cand_amp_start = candidates_scaffold_positions[off_target_scaffold][0]
                        cand_amp_end = candidates_scaffold_positions[off_target_scaffold][1]
                        if up_off.start_position in range(cand_amp_start, cand_amp_end):
                            continue
                    else:
                        if up_off not in new_up_off_targets_lst:
                            new_up_off_targets_lst.append(up_off)
                else:
                    if up_off not in new_up_off_targets_lst:
                        new_up_off_targets_lst.append(up_off)
            for down_off in candidate_amplicon.down_off_targets:
                if down_off.number_of_mismatches == 0:
                    off_target_scaffold = down_off.chromosome
                    if off_target_scaffold in candidates_scaffold_positions:
                        cand_amp_start = candidates_scaffold_positions[off_target_scaffold][0]
                        cand_amp_end = candidates_scaffold_positions[off_target_scaffold][1]
                        if down_off.start_position in range(cand_amp_start, cand_amp_end):
                            continue
                    else:
                        if down_off not in new_down_off_targets_lst:
                            new_down_off_targets_lst.append(down_off)
                else:
                    if down_off not in new_down_off_targets_lst:
                        new_down_off_targets_lst.append(down_off)
            candidate_amplicon.up_off_targets = new_up_off_targets_lst
            candidate_amplicon.down_off_targets = new_down_off_targets_lst
        else:
            new_off_targets_lst = []
            for off in candidate_amplicon.off_targets:
                if off.number_of_mismatches == 0:
                    off_target_scaffold = off.chromosome
                    if off_target_scaffold in candidates_scaffold_positions:
                        cand_amp_start = candidates_scaffold_positions[off_target_scaffold][0]
                        cand_amp_end = candidates_scaffold_positions[off_target_scaffold][1]
                        if off.start_position in range(cand_amp_start, cand_amp_end):
                            continue
                    else:
                        if off not in new_off_targets_lst:
                            new_off_targets_lst.append(off)
                else:
                    if off not in new_off_targets_lst:
                        new_off_targets_lst.append(off)
            candidate_amplicon.off_targets = new_off_targets_lst


def get_off_targets(candidate_amplicons_list: List[Amplicon_Obj], genome_fasta_file: str, out_path: str, pams: Tuple,
                    candidates_scaffold_positions: Dict[str, Tuple[int, int]], k: int, multiplex: int):
    """
    Find the off-targets for each candidate using BWA, store them in the candidate's off_targets_list attribute
    as a list of OffTarget objects. Then calculates the off-target scores for each off-target of each candidate and
    store the scores in each off-target's score attribute.

    :param candidate_amplicons_list: a list of amplicon candidates.
    :param out_path: output path for the algorithm results.
    :param genome_fasta_file: path to input FASTA format file of the genome.
    :param pams: tuple of PAM sequences of the Cas protein in use.
    :param candidates_scaffold_positions: dictionary of allele scaffold -> gene allele start,end indices.
    :param k: number of alleles to target with a single gRNA.
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    """

    if multiplex:
        # run off target search
        up_off_targets_sam = run_bwa(candidate_amplicons_list, genome_fasta_file, out_path, k, multiplex, "up")
        # extract off-targets from SAM file to pandas DataFrame
        up_off_targets_pd = extract_off_targets(up_off_targets_sam, genome_fasta_file, pams)
        # create a dictionary of sequence -> candidate
        up_sequence_to_candidate_dict = create_sequence_to_candidate_dict(candidate_amplicons_list, k, multiplex, "up")
        # add the found off-targets of each candidate to the candidate's off_targets_list
        add_off_targets(up_off_targets_pd, up_sequence_to_candidate_dict, multiplex, "up")

        # run off target search
        down_off_targets_sam = run_bwa(candidate_amplicons_list, genome_fasta_file, out_path, k, multiplex, "down")
        # extract off-targets from SAM file to pandas DataFrame
        down_off_targets_pd = extract_off_targets(down_off_targets_sam, genome_fasta_file, pams)
        # create a dictionary of sequence -> candidate
        down_sequence_to_candidate_dict = create_sequence_to_candidate_dict(candidate_amplicons_list, k, multiplex, "down")
        # add the found off-targets of each candidate to the candidate's off_targets_list
        add_off_targets(down_off_targets_pd, down_sequence_to_candidate_dict, multiplex, "down")
    else:
        # run off target search
        off_targets_sam = run_bwa(candidate_amplicons_list, genome_fasta_file, out_path, k, multiplex)
        # extract off-targets from SAM file to pandas DataFrame
        off_targets_pd = extract_off_targets(off_targets_sam, genome_fasta_file, pams)
        # create a dictionary of sequence -> candidate
        sequence_to_candidate_dict = create_sequence_to_candidate_dict(candidate_amplicons_list, k, multiplex)
        # add the found off-targets of each candidate to the candidate's off_targets_list
        add_off_targets(off_targets_pd, sequence_to_candidate_dict, multiplex)
        # remove on-targets
    remove_on_targets(candidate_amplicons_list, candidates_scaffold_positions, multiplex)
    # calculate the off-target scores for each off_target of each candidate
    if multiplex:
        calculate_scores(candidate_amplicons_list, k, multiplex, "up")
        calculate_scores(candidate_amplicons_list, k, multiplex, "down")
    else:
        calculate_scores(candidate_amplicons_list, k, multiplex)
    for candidate_amplicon in candidate_amplicons_list:
        for scaffold_amplicon in candidate_amplicon.scaffold_amplicons:
            candidate_amplicon.scaffold_amplicons[scaffold_amplicon].off_targets = candidate_amplicon.off_targets
            candidate_amplicon.scaffold_amplicons[scaffold_amplicon].up_off_targets = candidate_amplicon.up_off_targets
            candidate_amplicon.scaffold_amplicons[scaffold_amplicon].down_off_targets = candidate_amplicon.down_off_targets


def filt_off_targets(candidate_amplicons_list: List[Amplicon_Obj], genome_fasta_file: str, out_path: str, pams: Tuple,
                     candidates_scaffold_positions: Dict[str, Tuple[int, int]], k: int, multiplex: int) -> List[Amplicon_Obj]:
    """

    :param candidate_amplicons_list: a list of amplicon candidates
    :param genome_fasta_file: path to input FASTA format file of the genome
    :param out_path: output path for the algorithm results
    :param pams: tuple of PAM sequences of the Cas protein in use
    :param candidates_scaffold_positions: dictionary of allele scaffold -> gene allele start,end indices
    :param k: number of alleles to target with a single gRNA
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :return:
    """

    filtered_sorted_candidate_amplicons = []
    get_off_targets(candidate_amplicons_list, genome_fasta_file, out_path, pams, candidates_scaffold_positions,
                    k, multiplex)

    for candidate in candidate_amplicons_list:
        if multiplex:
            if candidate.up_off_targets[0].score < 0.15 and candidate.down_off_targets[0].score < 0.15:
                filtered_sorted_candidate_amplicons.append(candidate)
        else:
            if candidate.off_targets[0].score < 0.15:
                filtered_sorted_candidate_amplicons.append(candidate)

    return filtered_sorted_candidate_amplicons
