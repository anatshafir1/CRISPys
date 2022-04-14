import argparse
import os
import time

import pandas as pd


def createHeaderJob(path, job_name, ncpu=1, mem=16):
    """
    A function to create qsub file with activating crispys conda env and ssh to 0-247 machine
    :param path: path to log files
    :param job_name: job name
    :param ncpu: cpu number default 1
    :param mem: memory to use (in gb) default 16
    :return: a string the can be use to write sh file to run on the cluster (need to add command before running on the cluster)
    """
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += "#PBS -q itaym\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N " + job_name + "\n"
    text += "#PBS -e " + path + "/" + job_name + ".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name + ".OU" + "\n"
    text += "#PBS -p  3\n"
    text += "#PBS -l select=ncpus=" + str( ncpu ) + ":mem=" + str( mem ) + "gb\n"
    text += "source ~/.bashrc\n"
    text += "export PATH='$CONDA_PREFIX/bin:$PATH'\n"
    text += "conda activate crispys\n"
    return text


test_folders = ["gain_score/t_1", "gain_score/t_0", "N_internal_node/10",
                "N_internal_node/200", "scoring/CrisprMIT", "scoring/CCtop", "scoring/gold_off", "scoring/cfd",
                "where_in_gene/0.8", "where_in_gene/0.4", "algo/E", "algo/A", "threshold/th_0.8",
                "threshold/th_0.45", "N_poly_sites/12", "N_poly_sites/2", "PAM/pams_GG", "PAM/pams_GGAG"]


def run_crispys_test(code_folder, res_folder, code="git"):
    """
    This function will run crispys with different partametrs with provided code
    :param code_folder: The path to the folder where your crispys code is
    :param res_folder: The path to your results folder
    :return: write a folder for each test with the rsults of crispys and also wite test summary
    """

    # make folders for output if needed
    for folder in test_folders:
        if not os.path.isdir(res_folder + "/" + folder.split("/")[0]):
            os.system("mkdir " + res_folder + "/" + folder.split("/")[0])
        if not os.path.isdir(res_folder + "/" + folder):
            os.system("mkdir " + res_folder + "/" + folder)

    #### run cripys for each test
    # gain score t1:
    header = createHeaderJob(res_folder + "/gain_score/t_1", "gain_t1")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_1 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_1 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open(res_folder + "/gain_score/t_1/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/gain_score/t_1/Crispys.sh")

    # gain score t0:
    header = createHeaderJob(res_folder + "/gain_score/t_0", "gain_t0")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_0 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_0 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open( res_folder + "/gain_score/t_0/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/gain_score/t_0/Crispys.sh" )

    # N_poly_sites/12
    header = createHeaderJob( res_folder + "/N_poly_sites/12", "poly_12" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/N_poly_sites/12 --alg E --t 1 --v 0.8 --i 20 --where_in_gene 0.8 --s cfd_funct --ps 12"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/N_poly_sites/12 --alg E --t 1 --v 0.8 --i 20 --where_in_gene 0.8 --s cfd_funct --ps 12"

    with open( res_folder + "/N_poly_sites/12/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system("qsub " + res_folder + "/N_poly_sites/12/Crispys.sh" )

    # N_poly_sites/2
    header = createHeaderJob(res_folder + "/N_poly_sites/2", "poly_2")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/N_poly_sites/2 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct --ps 2"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/N_poly_sites/2 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct --ps 2"

    with open(res_folder + "/N_poly_sites/2/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_poly_sites/2/Crispys.sh")

    # N_internal_node/10
    header = createHeaderJob(res_folder + "/N_internal_node/10", "in_10")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/N_internal_node/10 --alg E --t 1 --v 0.8 --i 10 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/N_internal_node/10 --alg E --t 1 --v 0.8 --i 10 --where_in_gene 0.8 --s cfd_funct"

    with open(res_folder + "/N_internal_node/10/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_internal_node/10/Crispys.sh")

    # N_internal_node/200
    header = createHeaderJob(res_folder + "/N_internal_node/200", "in_200")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/N_internal_node/200 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/N_internal_node/200 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open(res_folder + "/N_internal_node/200/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_internal_node/200/Crispys.sh")

    # scoring/CrisprMIT
    header = createHeaderJob(res_folder + "/scoring/CrisprMIT", "fun_MIT")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CrisprMIT --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s MITScore"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CrisprMIT --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s MITScore"

    with open(res_folder + "/scoring/CrisprMIT/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/CrisprMIT/Crispys.sh")

    # scoring/CCtop
    header = createHeaderJob(res_folder + "/scoring/CCtop", "fun_CCtop")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CCtop --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s ccTop"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CCtop --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s ccTop"

    with open(res_folder + "/scoring/CCtop/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/CCtop/Crispys.sh")

    # scoring/cfd
    header = createHeaderJob(res_folder + "/scoring/cfd", "fun_cfd")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/cfd --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/cfd --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open(res_folder + "/scoring/cfd/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/cfd/Crispys.sh")

    # scoring/gold_off
    header = createHeaderJob(res_folder + "/scoring/gold_off", "gold_off")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/gold_off --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s gold_off"

    with open( res_folder + "/scoring/gold_off/Crispys.sh", "w" ) as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/gold_off/Crispys.sh" )

    # where_in_gene/0.8
    header = createHeaderJob(res_folder + "/where_in_gene/0.8", "gene_x0.8")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.8 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.8 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open(res_folder + "/where_in_gene/0.8/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/where_in_gene/0.8/Crispys.sh")

    # where_in_gene/0.4
    header = createHeaderJob(res_folder + "/where_in_gene/0.4", "gene_x0.4")
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.4 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.4 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.4 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.4 --s cfd_funct"

    with open(res_folder + "/where_in_gene/0.4/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system( "qsub " + res_folder + "/where_in_gene/0.4/Crispys.sh" )

    # algo/E
    header = createHeaderJob( res_folder + "/algo/E", "algo_E" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/E --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/E --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open( res_folder + "/algo/E/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/algo/E/Crispys.sh" )

    # algo/A
    header = createHeaderJob( res_folder + "/algo/A", "algo_A" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/A --alg A --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/A --alg A --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open( res_folder + "/algo/A/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/algo/A/Crispys.sh" )

    # threshold/th_0.8
    header = createHeaderJob( res_folder + "/threshold/th_0.8", "th_0.8" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.8 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.8 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open( res_folder + "/threshold/th_0.8/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/threshold/th_0.8/Crispys.sh" )

    # threshold/th_0.45
    header = createHeaderJob( res_folder + "/threshold/th_0.45", "th_0.45" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.45 --alg E --t 1 --v 0.45 --i 200 --where_in_gene 0.8 --s cfd_funct"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.45 --alg E --t 1 --v 0.45 --i 200 --where_in_gene 0.8 --s cfd_funct"

    with open( res_folder + "/threshold/th_0.45/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/threshold/th_0.45/Crispys.sh" )


    # PAM/pams_GG
    header = createHeaderJob( res_folder + "/PAM/pams_GG", "pams_GG" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/PAM/pams_GG --alg E --t 1 --v 0.45 --i 200 --where_in_gene 0.8 --s cfd_funct --PAMs 0"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/PAM/pams_GG --alg E --t 1 --v 0.45 --i 200 --where_in_gene 0.8, --s cfd_funct --PAMs 0"

    with open( res_folder + "/PAM/pams_GG/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/PAM/pams_GG/Crispys.sh" )


    # PAM/pams_GG and GA
    header = createHeaderJob( res_folder + "/PAM/pams_GGAG", "pams_GGAG" )
    if code == "git":
        cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/PAM/pams_GGAG --alg E --t 1 --v 0.45 --i 200 --where_in_gene 0.8 --s cfd_funct --PAMs 1"
    if code == "server":
        cmd = f"python {code_folder}/call_MULTICRISPR_Wrapper.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D004565/HOM04D004565.txt {res_folder}/PAM/pams_GGAG --alg E --t 1 --v 0.45 --i 200 --where_in_gene 0.8, --s cfd_funct --PAMs 1"

    with open( res_folder + "/PAM/pams_GGAG/Crispys.sh", "w" ) as f:
        f.write( header + "\n" + cmd )
    os.system( "qsub " + res_folder + "/PAM/pams_GGAG/Crispys.sh" )

def compare_output(old_res_folder, new_res_folder):
    """
          This function will take the results of a new crispys run and check for any difference in the
        output files compared to an existing output of crispys
        :param old_res_folder: The path to existing crispys test results (made by 'run_crispys_test()')
        :param new_res_folder: The path to 'new' output of crispys results on the test data.
        :return: The function will write a tsv file with the differences between the two outputs in each test.
        """
    # open file to store results
    res = open(new_res_folder + "/diff_res" + ".tsv", "w")

    # go over each result file and compare it with reference results
    for folder in test_folders:
        # wait for the file to be written
        while not os.path.exists( new_res_folder + "/" + folder + "/output.csv" ):
            time.sleep(15)

        # if there is a file, read it and the old version and compare the them
        if os.path.isfile( new_res_folder + "/" + folder + "/output.csv" ):
            time.sleep( 15 )
            # open the crispys outputs
            with open(new_res_folder + "/" + folder + "/output.csv", "r") as new, open(
                    old_res_folder + "/" + folder + "/output.csv", "r") as old:

                # read the old results lines as list of lines
                old_cri_code = old.readlines()
                old_cri_code = [line for line in old_cri_code if
                                (not line.startswith( "genes" )) and (not line.startswith( "sgRNA" ))]

                # make a set of all sgrnas
                sg_old = set( [line.split( "," )[0] for line in old_cri_code] )
                sg_and_score_old = [(line.split( "," )[0], line.split( "," )[1]) for line in old_cri_code]

                # read the new results lines as list of lines
                new_cri_code = new.readlines()
                # take only lines of sgRNAs
                new_cri_code = [line for line in new_cri_code if
                                (not line.startswith( "genes" )) and (not line.startswith( "sgRNA" ))]

                # make a set of all sgrnas
                sg_new = set( [line.split( "," )[0] for line in new_cri_code] )
                # make a list of all sgrnas and ther score
                sg_and_score_new = [(line.split( "," )[0], line.split( "," )[1]) for line in new_cri_code]

                # check if sg in old but not in new
                if sg_new.difference( sg_old ):
                    sg_diff_add = sg_new.difference( sg_old )  # give a set of sgrnas that are in the new results only
                    l = [sg for sg in sg_and_score_new if
                         sg[0] in sg_diff_add]  # take a list of new added sg's and their score
                    # write the new added sg and there score
                    res.write( f"in {folder} these sgRNA were added:\t {l}\n" )

                if (sg_old.difference( sg_new )):
                    sg_diff_miss = sg_old.difference(
                        sg_new )  # give a set of sgrnas that are in the old results only (missing from new results)
                    l = [sg for sg in sg_and_score_old if
                         sg[0] in sg_diff_miss]  # take a list of missing sg's and their score
                    # write the new added sg and there score
                    res.write( f"in {folder} these sgRNA were missing:\t {l}\n" )

                else:
                    res.write( f"No difference in {folder}\n" )

    res.close()


def compare_output_with_server(old_res_folder, new_res_folder):
    """
        This function will take the results of a new crispys server version run and check for any difference in the
        output files compared to an existing output of crispys git version
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
            time.sleep(15)

        # if there is a file, read it and the old version and compare the them
        if os.path.isfile(new_res_folder + "/" + folder + "/CRISPys_output.csv"):
            time.sleep( 15 )
            # open the crispys outputs
            with open( new_res_folder + "/" + folder + "/CRISPys_output.csv", "r" ) as new, open(
                    old_res_folder + "/" + folder + "/CRISPys_output.csv", "r" ) as old:

                # read the old results lines as list of lines
                old_cri_code = old.readlines()
                old_cri_code = [line for line in old_cri_code if
                                (not line.startswith( "genes" )) and (not line.startswith( "sgRNA" ))]

                # make a set of all sgrnas
                sg_old = set( [line.split( "," )[0] for line in old_cri_code] )
                sg_and_score_old = [(line.split( "," )[0], line.split( "," )[1]) for line in old_cri_code]

                # read the new results lines as list of lines
                new_cri_code = new.readlines()
                # take only lines of that starts with digit
                new_cri_code = [line.split( "," ) for line in new_cri_code if (line[0].isdigit())]
                # convert to data frame in order to clean the table
                df = pd.DataFrame( new_cri_code,
                                   columns=['sgRNA index', 'sgRNA', 'Score', 'Genes', 'Genes score', 'Target site',
                                            '#mms', 'Position'] )
                # rmeove empty lines
                server_sg_lst = df['sgRNA'].dropna().tolist()
                # remove empty string
                server_sg_lst = list( filter( None, server_sg_lst ) )
                # convert to set
                sg_new = set( server_sg_lst )

                # make a list of all sgrnas and ther score
                df = df.dropna()

                def get_sg_score(row):
                    if row['sgRNA']:
                        return (row['sgRNA'], row['Score'])

                sg_and_score_new = df.apply( get_sg_score, axis=1 )
                sg_and_score_new = sg_and_score_new.tolist()
                sg_and_score_new = list( filter( None, sg_and_score_new ) )

                # check if sg in old but not in new
                if sg_new.difference( sg_old ):
                    sg_diff_add = sg_new.difference( sg_old )  # give a set of sgrnas that are in the new results only
                    l = [sg for sg in sg_and_score_new if
                         sg[0] in sg_diff_add]  # take a list of new added sg's and their score
                    # write the new added sg and there score
                    res.write( f"in {folder} these sgRNA were added:\t {l}\n" )

                if (sg_old.difference( sg_new )):
                    sg_diff_miss = sg_old.difference(
                        sg_new )  # give a set of sgrnas that are in the old results only (missing from new results)
                    l = [sg for sg in sg_and_score_old if
                         sg[0] in sg_diff_miss]  # take a list of missing sg's and their score
                    # write the new added sg and there score
                    res.write( f"in {folder} these sgRNA were missing:\t {l}\n" )

                else:
                    res.write( f"No difference in {folder}\n" )

    res.close()


def compare_output_new_format(old_res_folder, new_res_folder):
    """
        This function will take the results of a new crispys server version run and check for any difference in the
        output files compared to an existing output of crispys git version
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
            time.sleep(15)

        # if there is a file, read it and the old version and compare with new
        if os.path.isfile(new_res_folder + "/" + folder + "/CRISPys_output.csv"):
            time.sleep( 15 )
            # open the crispys outputs
            with open( new_res_folder + "/" + folder + "/CRISPys_output.csv", "r" ) as new, open(
                    old_res_folder + "/" + folder + "/CRISPys_output.csv", "r" ) as old:

                # read the old results lines as list of lines
                old_cri_code = old.readlines()
                # take only lines of that starts with digit
                old_cri_code = [line.split( "," ) for line in old_cri_code if (line[0].isdigit())]
                # convert to data frame in order to clean the table
                df = pd.DataFrame(old_cri_code,
                                   columns=['sgRNA index', 'sgRNA', 'Score', 'Genes', 'Genes score', 'Target site',
                                            '#mms', 'Position'])
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
                        return (row['sgRNA'], row['Score'])

                sg_and_score_old = df.apply(get_sg_score, axis=1)
                sg_and_score_old = sg_and_score_old.tolist()
                sg_and_score_old = list(filter(None, sg_and_score_old))

                # read the new results lines as list of lines
                new_cri_code = new.readlines()
                # take only lines of that starts with digit
                new_cri_code = [line.split( "," ) for line in new_cri_code if (line[0].isdigit())]
                # convert to data frame in order to clean the table
                df = pd.DataFrame( new_cri_code,
                                   columns=['sgRNA index', 'sgRNA', 'Score', 'Genes', 'Genes score', 'Target site',
                                            '#mms', 'Position'] )
                # rmeove empty lines
                server_sg_lst = df['sgRNA'].dropna().tolist()
                # remove empty string
                server_sg_lst = list( filter( None, server_sg_lst ) )
                # convert to set
                sg_new = set( server_sg_lst )

                # make a list of all sgrnas and their scores
                df = df.dropna()

                sg_and_score_new = df.apply(get_sg_score, axis=1)
                sg_and_score_new = sg_and_score_new.tolist()
                sg_and_score_new = list(filter(None, sg_and_score_new))

                # check if sg in old but not in new
                if sg_new.difference(sg_old):
                    sg_diff_add = sg_new.difference(sg_old)  # give a set of sgrnas that are in the new results only
                    l = [sg for sg in sg_and_score_new if
                         sg[0] in sg_diff_add]  # take a list of new added sg's and their score
                    # write the new added sg and there score
                    res.write(f"in {folder} these sgRNA were added:\t {l}\n")

                if (sg_old.difference(sg_new)):
                    sg_diff_miss = sg_old.difference(
                        sg_new)  # give a set of sgrnas that are in the old results only (missing from new results)
                    l = [sg for sg in sg_and_score_old if
                         sg[0] in sg_diff_miss]  # take a list of missing sg's and their score
                    # write the new added sg and there score
                    res.write(f"in {folder} these sgRNA were missing:\t {l}\n")

                else:
                    res.write(f"No difference in {folder}\n")

    res.close()

def main(ref_folder, code_folder=None, res_folder_new=None, mode="run_and_compare", compare="git"):
    """
    Main function that run the tests
    :param code_folder: path to the crispys code folder
    :param res_folder: path to folder where the 'old' crispys results are (what we want to compare with)
    :param res_folder_new: path to folder where the 'new' crispys results are (generated by the new code we wrote)
    :param mode: the mode of comparison we want to make. can be 'run_and_compare' to run the new code and compare with existing results,
                                                        or 'run' for just running crispys on all the tests,
                                                        or 'compare' to compare to output of the test folders without running crispys.
    :param compare: type of results we use for comparison (what is the format of the 'new' results). can be 'git' or 'server'.
    :return: depend on the option. can create new folders with crispys results and a summary of the comparison.
    """

    if compare == "git":
        if mode == "run_and_compare":
            run_crispys_test(code_folder, res_folder_new)
            compare_output_new_format(ref_folder, res_folder_new)

        if mode == "run":
            run_crispys_test(code_folder, ref_folder)

        if mode == "compare":
            compare_output_new_format(ref_folder, res_folder_new)


def parse_arguments(parser):
    parser.add_argument( '--code_folder', '-code', type=str, help='The path to the crispys code folder' )
    parser.add_argument( '--ref_folder', '-ref', type=str, help='The path to the crispys results (reference)' )
    parser.add_argument( '--res_folder_new', '-new', type=str, help='The path to the new crispys results' )
    parser.add_argument( '--mode', '-m', default="run_and_compare", type=str,
                         help="mode of action, choose between 'run_and_compare', 'run' and 'compare'" )
    parser.add_argument( '--compare', '-c', default="git", type=str,
                         help="crispys type of output to use as the 'new', choose between 'git' and 'server'" )
    args = parser.parse_args()
    return args


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     args = parse_arguments( parser )
#     main( code_folder=args.code_folder,
#           ref_folder=args.ref_folder,
#           res_folder_new=args.res_folder_new,
#           mode=args.mode,
#           compare=args.compare )

# main(code_folder="/groups/itay_mayrose/udiland/crispys_code/CRISPys",
#      ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git/reference",
#      res_folder_new="/groups/itay_mayrose/udiland/crispys_test/test_files_git/res",
#      mode="run_and_compare")

main(code_folder="/groups/itay_mayrose/udiland/remote_deb/crispys_git/CRISPys",
     ref_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git/reference",
     res_folder_new="/groups/itay_mayrose/udiland/crispys_test/test_files_git/res",
     mode="run_and_compare")

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
