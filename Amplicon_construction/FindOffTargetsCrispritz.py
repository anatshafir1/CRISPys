import os
import sys
import time
from typing import List

import pandas as pd

from Amplicon_Obj import Amplicon_Obj


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
    :param output_path: A path containing the output of SPG
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
    off_targets_df = pd.read_csv(f"{crispritz_path}/crispritz.targets.txt", sep="\t")
    print(f"Number of off-targets: {off_targets_df.shape[0]}")
    return off_targets_df
