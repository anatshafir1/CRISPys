"""Module for finding off-targets"""
import os
import sys
import time

import pandas as pd
from typing import Dict, List

from Amplicon_Obj import Amplicon_Obj, OffTarget

from MOFF.MOFF_prediction import MOFF_score
from MOFF.MoffLoad import mtx1, mtx2, model


def create_crispritz_input_file(candidate_amplicons_list: List[Amplicon_Obj], crispritz_path: str) -> str:
    """
    :param candidate_amplicons_list: A list of CandidateWithOffTargets objects
    :param crispritz_path: A path to the crispys result folder where a folder for crispritz will be created
    :return: A path to the input for xxx and will write the input file to xxx
    """
    crispritz_infile = os.path.join(crispritz_path, 'crispritz_infile.txt')
    out = ""

    for candidate in candidate_amplicons_list:  # go over each candidate and get the guide sequence
        out += f"{candidate.target.seq[:-3]}NNN\n"
    with open(crispritz_infile, 'w') as f:
        f.write(out)
    return crispritz_infile


def create_pam_file(pam_file_path: str) -> str:
    """
    create the pam file for crispritz (for NGG pam)

    :param pam_file_path: path to folder location
    :return: pam file path
    """
    pam_file = f"{pam_file_path}/pamNGG.txt"
    if not os.path.exists(pam_file_path):
        os.makedirs(pam_file_path)
        with open(pam_file, "w") as f:
            f.write("NNNNNNNNNNNNNNNNNNNNNGG 3")
        return pam_file
    else:
        return pam_file


def run_crispritz(candidate_amplicons_list: List[Amplicon_Obj], output_path: str, genome_by_chr_path: str) -> pd.DataFrame:
    """
    This function runs crispritz and returns its output
    nucleotide codes (e.g. 'NRG' matches both NGG and NAG).

    :param candidate_amplicons_list: A list of CandidateWithOffTargets objects
    :param output_path: A path containing the output of CRISPys
    :param genome_by_chr_path: path to the folder where the files of each chromosome fasta file.
    :return: The output of crispritz as pd.DataFrame, where each row is a potential offtarget.
    """
    crispritz_path = f"{output_path}/crispritz"
    os.makedirs(crispritz_path, exist_ok=True)
    # create crispritz input file
    crispritz_infile = create_crispritz_input_file(candidate_amplicons_list, crispritz_path)
    pam_file_path = output_path + "/pams"
    pam_file = create_pam_file(pam_file_path)
    # run crispritz
    t0 = time.perf_counter()
    os.system(f"python {sys.exec_prefix}/bin/crispritz.py search {genome_by_chr_path}/ {pam_file} {crispritz_infile} {crispritz_path}/crispritz -mm 4 -r > /dev/null 2>&1")
    t1 = time.perf_counter()
    print(f"crispritz ran in {t1-t0} seconds")
    # get results
    crispritz_results = pd.read_csv(f"{crispritz_path}/crispritz.targets.txt", sep="\t")
    print(f"Number of off-targets: {crispritz_results.shape[0]}")
    return crispritz_results


# create dictionary of sequence:candidate
def create_sequence_to_candidate_dict(candidate_amplicons_list: List[Amplicon_Obj]):
    """
    This function takes a list of sgRNA candidates, and returns a sequence to candidate dictionary.

    :param candidate_amplicons_list: a list of sgRNA candidates
    :return: a dictionary: sequence -> an ActivationCandidate object where candidate.seq = sequence.
    """
    sequence_to_candidate_dict = {}
    for i, candidate in enumerate(candidate_amplicons_list):
        if candidate.target.seq[:-3] in sequence_to_candidate_dict:
            sequence_to_candidate_dict[candidate.target.seq[:-3]] += [candidate_amplicons_list[i]]
        else:
            sequence_to_candidate_dict[candidate.target.seq[:-3]] = [candidate_amplicons_list[i]]
    return sequence_to_candidate_dict


# add to each "Candidate" its off-targets
def add_crispritz_off_targets(crispritz_results, sequence_to_candidate_dict: Dict[str, List[Amplicon_Obj]]):
    """
    This function adds all found off-targets to each CandidateWithOffTargets using the crispritz results.

    :param crispritz_results: The output of crispritz as a pd datatable, where each row is a potential offtarget.
    :param sequence_to_candidate_dict: sequence -> a CandidateWithOffTargets object with the proper sequence
    """
    # apply the 'get_off_target' function on each row in the crispritz table results
    crispritz_results.apply(get_off_target, args=(sequence_to_candidate_dict,), axis=1)
    return


# This function will be used in an apply command' it reads crispritz results and create an offtarget object out of each one
def get_off_target(x, sequence_to_candidate_dict: Dict[str, List[Amplicon_Obj]]):
    """
    A function to use with apply on crispritz result table
    it takes a row of crispritz results and a dictionary of sequence:candidate, and make an OffTarget
    object from the crispritz results and add it to the candidate offtargets list.

    :param x: a row in crispritz results table
    :param sequence_to_candidate_dict: a dictionary of sequence:candidate
    """
    candidates_list = sequence_to_candidate_dict[x['crRNA'][:20]]
    off_target = OffTarget(x['DNA'].upper(), x['Chromosome'].split(" ")[0], int(x['Position']), x['Direction'], int(x['Mismatches']))
    legit_letters = True
    for char in off_target.seq:
        if char not in {"A", "C", "T", "G"}:
            legit_letters = False
            break
    if legit_letters:
        for candidate in candidates_list:
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


def calculate_scores(candidate_amplicons_list: List[Amplicon_Obj]):
    """
    Calculate the off-target scores for each off-target of each candidate and store the scores in each off-target's
    score attribute.

    :param candidate_amplicons_list: a list of Amplicon candidates
    """
    batch_off_targets_list = []
    batch_candidates_list = []
    for candidate in candidate_amplicons_list:
        batch_candidates_list += [candidate.target.seq for _ in range(len(candidate.off_targets))]
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


def get_off_targets(candidate_amplicons_list: List[Amplicon_Obj], genome_chroms_path: str, out_path: str):
    """
    Find the off-targets for each candidate using crispritz, store them in the candidate's off_targets_list attribute
    as a list of OffTarget objects. Then calculates the off-target scores for each off-target of each candidate and
    store the scores in each off-target's score attribute.

    :param candidate_amplicons_list: a list of sgRNA candidates
    :param genome_chroms_path: path to directory in which the chromosome FASTA files are stored
    :param out_path: output path for the algorithm results
    """
    # run crispritz
    off_targets_pd = run_crispritz(candidate_amplicons_list, out_path, genome_chroms_path)
    # create a dictionary of sequence -> candidate
    sequence_to_candidate_dict = create_sequence_to_candidate_dict(candidate_amplicons_list)
    # add the found off-targets of each candidate to the candidate's off_targets_list
    add_crispritz_off_targets(off_targets_pd, sequence_to_candidate_dict)
    # remove on-targets
    for candidate_amplicon in candidate_amplicons_list:
        scaffold = candidate_amplicon.scaffold
        amplicon_range_start = candidate_amplicon.start_idx
        amplicon_range_end = candidate_amplicon.end_idx
        new_off_targets_lst = []
        for off in candidate_amplicon.off_targets:
            if off.number_of_mismatches == 0:
                if off.chromosome == scaffold and (off.start_position in range(amplicon_range_start, amplicon_range_end)):
                    continue
            else:
                new_off_targets_lst.append(off)
        candidate_amplicon.off_targets = new_off_targets_lst
    # calculate the off-target scores for each off_target of each candidate
    calculate_scores(candidate_amplicons_list)
