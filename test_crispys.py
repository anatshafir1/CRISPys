import argparse
import os
import time
import globals
import pandas as pd

def createHeaderJob(path, job_name, queue, ncpu=1, mem=16):
    """
    A function to create qsub file with activating crispys conda env and ssh to 0-247 machine

    :param path: path to log files
    :param job_name: job name
    :param queue:
    :param ncpu: cpu number default 1
    :param mem: memory to use (in gb) default 16
    :return: a string that can be used to write sh file to run on the cluster (need to add command before running on the cluster)
    """
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += f"#PBS -q {queue}\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N " + job_name + "\n"
    text += "#PBS -e " + path + "/" + job_name + ".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name + ".OU" + "\n"
    text += "#PBS -p  3\n"
    text += "#PBS -l select=ncpus=" + str(ncpu) + ":mem=" + str(mem) + "gb\n"
    text += "source ~/.bashrc\n"
    text += "export PATH='$CONDA_PREFIX/bin:$PATH'\n"
    text += "conda activate crispys\n"
    return text


test_folders = ["gain_score/t_1", "gain_score/t_0", "N_internal_node/10",
                "N_internal_node/200", "scoring/CrisprMIT", "scoring/CCtop", "scoring/gold_off", "scoring/cfd",
                "scoring/crispr_net", "scoring/DeepHF", "where_in_gene/0.8", "where_in_gene/0.4", "algo/E", "algo/A",
                "threshold/th_0.8", "threshold/th_0.45", "N_poly_sites/12", "N_poly_sites/2", "PAM/pams_GG",
                "PAM/pams_GGAG", "singletons_from_crispys/with_singletons", "singletons_from_crispys/without_singletons"]


def run_crispys_test(code_folder, res_folder, queue):
    """
    This function will run crispys with different parameters with provided code
    :param code_folder: The path to the folder where your crispys code is
    :param res_folder: The path to your results folder
    :param queue:
    :return: write a folder for each test with the results of crispys and also write test summary
    """

    # clear folder old content
    os.system(f"rm -r {res_folder}/*")
    # make folders for output if needed
    for folder in test_folders:
        if not os.path.isdir(res_folder + "/" + folder.split("/")[0]):
            os.system("mkdir " + res_folder + "/" + folder.split("/")[0])
        if not os.path.isdir(res_folder + "/" + folder):
            os.system("mkdir " + res_folder + "/" + folder)

    # run crispys for each test
    # gain score t1:
    header = createHeaderJob(path=res_folder + "/gain_score/t_1", job_name="gain_t1", queue=queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_1 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/gain_score/t_1/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/gain_score/t_1/Crispys.sh")

    # gain score t0:
    header = createHeaderJob(res_folder + "/gain_score/t_0", "gain_t0", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_0 --alg gene_homology --use_thr 0 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/gain_score/t_0/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/gain_score/t_0/Crispys.sh")

    # N_poly_sites/12
    header = createHeaderJob(res_folder + "/N_poly_sites/12", "poly_12", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/N_poly_sites/12 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 20 --where_in_gene 0.8 --scoring_function cfd_funct --max_target_polymorphic_sites 12"

    with open(res_folder + "/N_poly_sites/12/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_poly_sites/12/Crispys.sh")

    # N_poly_sites/2
    header = createHeaderJob(res_folder + "/N_poly_sites/2", "poly_2", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/N_poly_sites/2 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct --max_target_polymorphic_sites 2"

    with open(res_folder + "/N_poly_sites/2/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_poly_sites/2/Crispys.sh")

    # N_internal_node/10
    header = createHeaderJob(res_folder + "/N_internal_node/10", "in_10", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_6/HOM04D000221_6.txt {res_folder}/N_internal_node/10 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 10 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/N_internal_node/10/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_internal_node/10/Crispys.sh")

    # N_internal_node/200
    header = createHeaderJob(res_folder + "/N_internal_node/200", "in_200", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_6/HOM04D000221_6.txt {res_folder}/N_internal_node/200 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/N_internal_node/200/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_internal_node/200/Crispys.sh")

    # scoring/CrisprMIT
    header = createHeaderJob(res_folder + "/scoring/CrisprMIT", "fun_MIT", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CrisprMIT --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function CrisprMIT"

    with open(res_folder + "/scoring/CrisprMIT/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/CrisprMIT/Crispys.sh")

    # scoring/CCtop
    header = createHeaderJob(res_folder + "/scoring/CCtop", "fun_CCtop", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CCtop --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function CCTop"

    with open(res_folder + "/scoring/CCtop/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/CCtop/Crispys.sh")

    # scoring/cfd
    header = createHeaderJob(res_folder + "/scoring/cfd", "fun_cfd", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/cfd --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/scoring/cfd/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/cfd/Crispys.sh")

    # scoring/gold_off
    header = createHeaderJob(res_folder + "/scoring/gold_off", "gold_off", queue, ncpu=globals.n_cores_for_gold_off)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/gold_off --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function gold_off"

    with open(res_folder + "/scoring/gold_off/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/gold_off/Crispys.sh")

    # scoring/crispr_net
    header = createHeaderJob(res_folder + "/scoring/crispr_net", "crispr_net", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/crispr_net --alg gene_homology --use_thr 1 --omega 0.5 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function crispr_net"

    with open(res_folder + "/scoring/crispr_net/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/crispr_net/Crispys.sh")

    # scoring/DeepHF
    header = createHeaderJob(res_folder + "/scoring/DeepHF", "DeepHF", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/DeepHF --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function deephf"

    with open(res_folder + "/scoring/DeepHF/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/DeepHF/Crispys.sh")

    # where_in_gene/0.8
    header = createHeaderJob(res_folder + "/where_in_gene/0.8", "gene_x0.8", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.8 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/where_in_gene/0.8/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/where_in_gene/0.8/Crispys.sh")

    # where_in_gene/0.4
    header = createHeaderJob(res_folder + "/where_in_gene/0.4", "gene_x0.4", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.4 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.4 --scoring_function cfd_funct"

    with open(res_folder + "/where_in_gene/0.4/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/where_in_gene/0.4/Crispys.sh")

    # algo/E
    header = createHeaderJob(res_folder + "/algo/E", "algo_E", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/E --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/algo/E/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/algo/E/Crispys.sh")

    # algo/A
    header = createHeaderJob(res_folder + "/algo/A", "algo_A", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/A --alg default --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/algo/A/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/algo/A/Crispys.sh")

    # threshold/th_0.8
    header = createHeaderJob(res_folder + "/threshold/th_0.8", "th_0.8", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.8 --alg gene_homology --use_thr 1 --omega 0.8 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/threshold/th_0.8/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/threshold/th_0.8/Crispys.sh")

    # threshold/th_0.45
    header = createHeaderJob(res_folder + "/threshold/th_0.45", "th_0.45", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.45 --alg gene_homology --use_thr 1 --omega 0.45 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct"

    with open(res_folder + "/threshold/th_0.45/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/threshold/th_0.45/Crispys.sh")

    # PAM/pams_GG
    header = createHeaderJob(res_folder + "/PAM/pams_GG", "pams_GG", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/PAM/pams_GG --alg gene_homology --use_thr 1 --omega 0.45 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct --pams 0"

    with open(res_folder + "/PAM/pams_GG/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/PAM/pams_GG/Crispys.sh")

    # PAM/pams_GG and GA
    header = createHeaderJob(res_folder + "/PAM/pams_GGAG", "pams_GGAG", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/PAM/pams_GGAG --alg gene_homology --use_thr 1 --omega 0.45 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct --pams 1"

    with open(res_folder + "/PAM/pams_GGAG/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/PAM/pams_GGAG/Crispys.sh")

    # singletons_from_crispys/with_singletons
    header = createHeaderJob(res_folder + "/singletons_from_crispys/with_singletons", "with_singletons", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/singletons_from_crispys/with_singletons --alg gene_homology --use_thr 1 --omega 0.45 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct --singletons_from_crispys 1"

    with open(res_folder + "/singletons_from_crispys/with_singletons/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/singletons_from_crispys/with_singletons/Crispys.sh")

    # singletons_from_crispys/without_singletons
    header = createHeaderJob(res_folder + "/singletons_from_crispys/without_singletons", "without_singletons", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/singletons_from_crispys/without_singletons --alg gene_homology --use_thr 1 --omega 0.45 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct --singletons_from_crispys 0"

    with open(res_folder + "/singletons_from_crispys/without_singletons/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/singletons_from_crispys/without_singletons/Crispys.sh")

    # genes_of_interest/genes_of_interest
    header = createHeaderJob(res_folder + "/genes_of_interest/genes_of_interest", "genes_of_interest", queue)
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/caldararu/test_crispys/HOM05D000360_1.fa {res_folder}/singletons_from_crispys/without_singletons -- --alg gene_homology --use_thr 1 --omega 0.45 --internal_node_candidates 200 --where_in_gene 0.8 --scoring_function cfd_funct --singletons_from_crispys 1 --min_desired_genes_fraction  0"

    with open(res_folder + "/singletons_from_crispys/without_singletons/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/singletons_from_crispys/without_singletons/Crispys.sh")


def compare_output(old_res_folder, new_res_folder):
    """
        This function will take the results of a new crispys server version run and check for any difference in the
        output files compared to an existing output of crispys ## This is updated to the latest format (server)
        :param old_res_folder: The path to existing crispys test results (made by 'run_crispys_test()')
        :param new_res_folder: The path to 'new' output of crispys server results on the test data.
        :return: The function will write a tsv file with the differences between the two outputs in each test.
        """
    # open file to store results
    res = open(new_res_folder + "/diff_res" + ".tsv", "w")

    # go over each result file and compare it with reference results
    for folder in test_folders:
        # wait for the file to be written
        while not os.path.exists(new_res_folder + "/" + folder + "/CRISPys_output.csv"):
            time.sleep(5)

        # if there is a file, read it and the old version and compare with new
        if os.path.isfile(new_res_folder + "/" + folder + "/CRISPys_output.csv"):
            time.sleep(15)  # the file is created and then written, so it might be analyzed before the finish. so I added wait time
            # open the crispys outputs
            with open(new_res_folder + "/" + folder + "/CRISPys_output.csv", "r") as new, open(
                    old_res_folder + "/" + folder + "/CRISPys_output.csv", "r") as old:
                # special check for DeepHF output (on-score results)
                if folder == "scoring/DeepHF":
                    # get old sg
                    deephf_old = old.readlines()
                    deephf_old = set([line.split(",")[0] for line in deephf_old])
                    # get new sg
                    deephf_new = new.readlines()
                    deephf_new = set([line.split(",")[0] for line in deephf_new])

                    # report sg added to DeepHF
                    if deephf_new.difference(deephf_old):
                        sg_diff_add = deephf_new.difference(
                            deephf_old)  # give a set of sgrnas that are in the new results only
                        line = [sg for sg in deephf_new if sg in sg_diff_add]  # take a list of new added sgs and their score
                        # write the new added sg and there score
                        res.write(f"in {folder} these sgRNA were added:\t {line}\n")

                    # report sg missing in DeepHF
                    if deephf_old.difference(deephf_new):
                        sg_diff_miss = deephf_old.difference(deephf_new)
                        line = [sg for sg in deephf_old if sg in sg_diff_miss]
                        # write the new added sg and there score
                        res.write(f"in {folder} these sgRNA were missing:\t {line}\n")

                    else:
                        res.write(f"No difference in {folder}\n")

                    # go to next result folder (next iteration)
                    continue

                # read the old results lines as list of lines
                old_cri_code = old.readlines()
                # take only lines of that starts with digit (because this is the line of sgrna)
                old_cri_code = [line.split(",") for line in old_cri_code if (line[0].isdigit())]
                # convert to data frame in order to clean the table
                df = pd.DataFrame(old_cri_code,
                                  columns=['sgRNA index', 'sgRNA', 'Score', 'Genes', 'Genes score', 'Target site',
                                           '#mms', 'Position'])

                # check if the table is empty
                if df.empty:
                    print(f"no results for {folder}")
                    res.write(f"no results for {folder}")
                    continue
                # remove empty lines
                sg_lst = df['sgRNA'].dropna().tolist()
                # remove empty string
                sg_lst = list(filter(None, sg_lst))
                # convert to set
                sg_old = set(sg_lst)

                # make a list of all sgrnas and their scores
                df = df.dropna()

                def get_sg_score(row):
                    if row['sgRNA']:
                        return row['sgRNA'], row['Score']

                sg_and_score_old = df.apply(get_sg_score, axis=1)
                sg_and_score_old = sg_and_score_old.tolist()
                sg_and_score_old = list(filter(None, sg_and_score_old))

                # read the new results lines as list of lines
                new_cri_code = new.readlines()
                # take only lines of that starts with digit
                new_cri_code = [line.split(",") for line in new_cri_code if (line[0].isdigit())]
                # convert to data frame in order to clean the table
                df = pd.DataFrame(new_cri_code,
                                  columns=['sgRNA index', 'sgRNA', 'Score', 'Genes', 'Genes score', 'Target site',
                                           '#mms', 'Position'])
                # remove empty lines
                server_sg_lst = df['sgRNA'].dropna().tolist()
                # remove empty string
                server_sg_lst = list(filter(None, server_sg_lst))
                # convert to set
                sg_new = set(server_sg_lst)

                # make a list of all sgrnas and their scores
                df = df.dropna()

                sg_and_score_new = df.apply(get_sg_score, axis=1)
                sg_and_score_new = sg_and_score_new.tolist()
                sg_and_score_new = list(filter(None, sg_and_score_new))

                # check if sg in old but not in new
                if sg_new.difference(sg_old):
                    sg_diff_add = sg_new.difference(sg_old)  # give a set of sgrnas that are in the new results only
                    line = [sg for sg in sg_and_score_new if
                         sg[0] in sg_diff_add]  # take a list of new added sgs and their score
                    # write the new added sg and there score
                    res.write(f"in {folder} these sgRNA were added:\t {line}\n")

                if sg_old.difference(sg_new):
                    sg_diff_miss = sg_old.difference(
                        sg_new)  # give a set of sgrnas that are in the old results only (missing from new results)
                    line = [sg for sg in sg_and_score_old if
                         sg[0] in sg_diff_miss]  # take a list of missing sgs and their score
                    # write the new added sg and there score
                    res.write(f"in {folder} these sgRNA were missing:\t {line}\n")

                else:
                    res.write(f"No difference in {folder}\n")

    res.close()


def main(ref_folder, queue="itaym", code_folder=None, res_folder_new=None, mode="run_and_compare"):
    """
    Main function that run the tests

    :param ref_folder: path to folder where the 'old' crispys results are (what we want to compare with)
    :param code_folder: path to the crispys code folder
    :param queue:
    :param res_folder_new: path to folder where the 'new' crispys results are (generated by the new code we wrote)
    :param mode: the mode of comparison we want to make. can be 'run_and_compare' to run the new code and compare with existing results,
                                                        or 'run' for just running crispys on all the tests,
                                                        or 'compare' to compare to output of the test folders without running crispys.
    :return: depend on the option. can create new folders with crispys results and a summary of the comparison.
    """

    if mode == "run_and_compare":
        run_crispys_test(code_folder, res_folder_new, queue)
        compare_output(ref_folder, res_folder_new)

    if mode == "run":
        run_crispys_test(code_folder, ref_folder, queue)

    if mode == "compare":
        compare_output(ref_folder, res_folder_new)


def parse_arguments(parser):
    parser.add_argument('--code_folder', '-code', type=str, help='The path to the crispys code folder')
    parser.add_argument('--ref_folder', '-ref', type=str, help='The path to the crispys results (reference)')
    parser.add_argument('--res_folder_new', '-new', type=str, help='The path to the new crispys results')
    parser.add_argument('--mode', '-mode', default="run_and_compare", type=str,
                        help="mode of action, choose between 'run_and_compare', 'run' and 'compare'")
    parser.add_argument('--queue', '-q', default="itaym", help="select the queue to run")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    main(code_folder=args.code_folder,
         ref_folder=args.ref_folder,
         res_folder_new=args.res_folder_new,
         mode=args.mode,
         queue=args.queue)
#
# main(code_folder="/groups/itay_mayrose/udiland/crispys_code/CRISPys",
#      ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git/reference",
#      res_folder_new="/groups/itay_mayrose/udiland/crispys_test/test_files_git/res",
#      mode="run_and_compare", queue="itaymaa")

# main(code_folder="/groups/itay_mayrose/udiland/remote_deb/crispys_git/CRISPys",
#      ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git/reference",
#      res_folder_new="/groups/itay_mayrose/udiland/crispys_test/test_files_git/res",
#      mode="run_and_compare", queue="itaym")

# run crispys code to create reference folders
# main(code_folder="/groups/itay_mayrose/udiland/remote_deb/crispys_git/CRISPys",
#      ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git/reference",
#      mode="run")

# run crispys with server code and compare to existing git results
# main(code_folder="/groups/itay_mayrose/udiland/remote_deb/crispys_git",
#      ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git/reference",
#      res_folder_new="/groups/itay_mayrose/udiland/crispys_test/compare_server_git/compare_res",
#      mode="run_and_compare",
#      compare="server")

# main(ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git",
#      res_folder_new="/groups/itay_mayrose/udiland/crispys_test/compare_server_git/compare_res",
#      mode="compare",
#      compare="server")
