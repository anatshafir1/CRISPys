import pandas as pd
import typing
import pickle
import subprocess
import subgroup_res, Candidate
import os
import time

def off_target_search(crispys_output_path: str):
    with open(crispys_output_path, 'rb') as f:
        subgroup_list = pickle.load(f)
    sgRNA_list = []
    for subgroup in subgroup_list:
        sgRNA_list += subgroup.candidate_lst

    return crispys_output_path


def create_cas_offinder_input_path(list_of_sgRNAs, genome_path='/groups/itay_mayrose/caldararu/cas_offinder/osa.fasta',
                                   max_number_of_mismatches=4, pam='NGG'):
    """
    This function creates an output file for the cas-offinder tool
    :param pam:
    :param list_of_sgRNAs:
    :param genome_path:
    :param max_number_of_mismatches:
    :return:
    """
    infile = "/groups/itay_mayrose/caldararu/cas_offinder/input.txt"
    outfile = "/groups/itay_mayrose/caldararu/cas_offinder/output_pycharm.txt"
    out = f"{genome_path}\nNNNNNNNNNNNNNNNNNNNN{pam}\n"
    for sgRNA in list_of_sgRNAs:
        out += f"{sgRNA} {max_number_of_mismatches}\n"
    with open(infile, 'w') as f:
        f.write(out)
    subprocess.Popen(["cas-offinder","/groups/itay_mayrose/caldararu/cas_offinder/input.txt","C","/groups/itay_mayrose/caldararu/cas_offinder/output_pycharm.txt"])
    while not os.path.exists(outfile):
        time.sleep(0.2)
    columns = ['sgRNA','chromosome','start_position','off_target_sequence','strand','number_of_mismatches']
    df = pd.read_csv(outfile, sep="\t", header=None,index_col = False,names = columns)
    return df

def get_genomic_feature():
    """
    this function takes the off target and annotation data as input,
    and returns the feature of the off target
    :return:
    """
    pass
create_cas_offinder_input_path(
    ['GGATCGCGCGCTTGGTATGA', 'TGAGACCGGCACGGTCTGTT', 'GCCGTATACTGCGTCATGTT', 'GGATCGTGCGCTTGGTATGA',
     'GGATCGTGCGTTTGGTATGA', 'TTGGAAGATGCGTGGAGATC', 'GAATCATGCGTTTGGTATGA', 'CATGTACATGCTTCAGATAT',
     'GACTTGGCTTCAGTGGTCTC', 'GACTTGGCTCCAGTGGTCTC'])

"""
Off target search and filtering:
    1. Find all off targets (sequence and coordinates). - Done!
    2. Find the genomic feature of each off target.
    3. Check if the off target falls on a gene not from the family.
    4. filter the sgRMAs based on thresholds for different features.
    5. return the filtered list of sgRNAs under the same format (subgroup_res which contains list of candidates).
"""