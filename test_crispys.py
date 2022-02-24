import argparse
import os
import time


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
    text += "#PBS -l select=ncpus=" + str(ncpu) + ":mem=" + str(mem) + "gb\n"
    text += "source ~/.bashrc\n"
    text += "export PATH='$CONDA_PREFIX/bin:$PATH'\n"
    text += "conda activate crispys\n"
    return text


test_folders = ["gain_score/t_1", "gain_score/t_0", "N_internal_node/10",
                "N_internal_node/200", "scoring/CrisprMIT", "scoring/CCtop", "scoring/CCtop", "scoring/cfd",
                "where_in_gene/0.8", "where_in_gene/0.4", "algo/E", "algo/A", "threshold/th_0.8",
                "threshold/th_0.45", "N_poly_sites/12", "N_poly_sites/2"]



def run_crispys_test(code_folder, res_folder):
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
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_1 --alg E --t 1 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/gain_score/t_1/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/gain_score/t_1/Crispys.sh")

    # gain score t0:
    header = createHeaderJob(res_folder + "/gain_score/t_0", "gain_t0")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000012_7/HOM04D000012_7.txt {res_folder}/gain_score/t_0 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/gain_score/t_0/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/gain_score/t_0/Crispys.sh")

    # N_poly_sites/12
    header = createHeaderJob(res_folder + "/N_poly_sites/12", "poly_12")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000034/HOM04D000034.txt {res_folder}/N_poly_sites/12 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct --ps 12"
    with open(res_folder + "/N_poly_sites/12/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_poly_sites/12/Crispys.sh")

    # N_poly_sites/2
    header = createHeaderJob(res_folder + "/N_poly_sites/2", "poly_2")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000034/HOM04D000034.txt {res_folder}/N_poly_sites/2 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct --ps 2"
    with open(res_folder + "/N_poly_sites/2/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_poly_sites/2/Crispys.sh")

    # N_internal_node/10
    header = createHeaderJob(res_folder + "/N_internal_node/10", "in_10")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/N_internal_node/10 --alg E --t 0 --v 0.8 --i 10 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/N_internal_node/10/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_internal_node/10/Crispys.sh")

    # N_internal_node/200
    header = createHeaderJob(res_folder + "/N_internal_node/200", "in_200")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/N_internal_node/200 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/N_internal_node/200/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/N_internal_node/200/Crispys.sh")

    # scoring/CrisprMIT
    header = createHeaderJob(res_folder + "/scoring/CrisprMIT", "fun_MIT")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CrisprMIT --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s MITScore"
    with open(res_folder + "/scoring/CrisprMIT/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/CrisprMIT/Crispys.sh")

    # scoring/CCtop
    header = createHeaderJob(res_folder + "/scoring/CCtop", "fun_CCtop")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/CCtop --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s ccTop"
    with open(res_folder + "/scoring/CCtop/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/CCtop/Crispys.sh")

    # scoring/cfd
    header = createHeaderJob(res_folder + "/scoring/cfd", "fun_cfd")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000350/HOM04D000350.txt {res_folder}/scoring/cfd --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/scoring/cfd/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/scoring/cfd/Crispys.sh")

    # where_in_gene/0.8
    header = createHeaderJob(res_folder + "/where_in_gene/0.8", "gene_x0.8")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.8 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/where_in_gene/0.8/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/where_in_gene/0.8/Crispys.sh")

    # where_in_gene/0.4
    header = createHeaderJob(res_folder + "/where_in_gene/0.4", "gene_x0.4")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000054_6/HOM04D000054_6.txt {res_folder}/where_in_gene/0.4 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.4 --s cfd_funct"
    with open(res_folder + "/where_in_gene/0.4/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/where_in_gene/0.4/Crispys.sh")

    # algo/E
    header = createHeaderJob(res_folder + "/algo/E", "algo_E")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/E --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/algo/E/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/algo/E/Crispys.sh")

    # algo/A
    header = createHeaderJob(res_folder + "/algo/A", "algo_A")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000632/HOM04D000632.txt {res_folder}/algo/A --alg A --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/algo/A/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/algo/A/Crispys.sh")

    # threshold/th_0.8
    header = createHeaderJob(res_folder + "/threshold/th_0.8", "th_0.8")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.8 --alg E --t 0 --v 0.8 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/threshold/th_0.8/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/threshold/th_0.8/Crispys.sh")

    # threshold/th_0.45
    header = createHeaderJob(res_folder + "/threshold/th_0.45", "th_0.45")
    cmd = f"python {code_folder}/Stage0.py /groups/itay_mayrose/udiland/crispys_test/test_files_git/HOM04D000221_5/HOM04D000221_5.txt {res_folder}/threshold/th_0.45 --alg E --t 0 --v 0.45 --i 200 --where_in_gene 0.8 --s cfd_funct"
    with open(res_folder + "/threshold/th_0.45/Crispys.sh", "w") as f:
        f.write(header + "\n" + cmd)
    os.system("qsub " + res_folder + "/threshold/th_0.45/Crispys.sh")


def compare_output(old_res_folder, new_res_folder):
    """
        This function will take the results of a new crispys run and check for any difference in the output file
        :param res_folder:
        :return:
        """
    # open file to store results
    res = open(new_res_folder + "/diff_res" + ".tsv", "w")

    # go over each result file and compare it with reference results
    for folder in test_folders:
        # wait for the file to be written
         while not os.path.exists(new_res_folder + "/" + folder + "/output.csv"):
             time.sleep(15)

         # if there is a file, read it and the old version and compare the them
         if os.path.isfile(new_res_folder + "/" + folder + "/output.csv"):
            time.sleep(15)
            # open the crispys outputs
            with open(new_res_folder + "/" + folder + "/output.csv", "r") as new, open(
                    old_res_folder + "/" + folder + "/output.csv", "r") as old:

                # read the old results lines as list of lines
                old_cri_code = old.readlines()
                old_cri_code = [line for line in old_cri_code if (not line.startswith("genes")) and (not line.startswith("sgRNA"))]

                #make a set of all sgrnas
                sg_old = set([line.split(",")[0] for line in old_cri_code])
                sg_and_score_old = [(line.split(",")[0], line.split(",")[1]) for line in old_cri_code]

                # read the new results lines as list of lines
                new_cri_code = new.readlines()
                # take only lines of sgRNAs
                new_cri_code = [line for line in new_cri_code if (not line.startswith("genes")) and (not line.startswith("sgRNA"))]

                #make a set of all sgrnas
                sg_new = set([line.split(",")[0] for line in new_cri_code])
                #make a list of all sgrnas and ther score
                sg_and_score_new = [(line.split(",")[0], line.split(",")[1]) for line in new_cri_code]


                # check if sg in old but not in new
                if sg_new.difference(sg_old):
                    sg_diff_add = sg_new.difference(sg_old) # give a set of sgrnas that are in the new results only
                    l = [sg for sg in sg_and_score_new if sg[0] in sg_diff_add] # take a list of new added sg's and their score
                    # write the new added sg and there score
                    res.write(f"in {folder} these sgRNA were added:\t {l}\n")

                if (sg_old.difference(sg_new)):
                    sg_diff_miss = sg_old.difference(sg_new)  # give a set of sgrnas that are in the old results only (missing from new results)
                    l = [sg for sg in sg_and_score_old if sg[0] in sg_diff_miss]  # take a list of missing sg's and their score
                    # write the new added sg and there score
                    res.write(f"in {folder} these sgRNA were missing:\t {l}\n")

                else:
                    res.write(f"No difference in {folder}\n")

    res.close()


def main(code_folder, res_folder, res_folder_new=None, mode="run_and_compare"):
    if mode == "run_and_compare":
        run_crispys_test(code_folder, res_folder_new)
        compare_output(res_folder, res_folder_new)

    if mode == "run":
        run_crispys_test(code_folder, res_folder)

    if mode == "compare":
        compare_output(res_folder, res_folder_new)


def parse_arguments(parser):
    parser.add_argument('--code_folder', '-code', type=str, help='The path to the crispys code folder')
    parser.add_argument('--res_folder', '-ref', type=str, help='The path to the crispys results (reference)')
    parser.add_argument('--res_folder_new', '-new', type=str, help='The path to the new crispys results')
    parser.add_argument('--mode', '-m', default="run_and_compare", type=str, help="mode of action, choose between 'run_and_compare', 'run' and 'compare'")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    main(code_folder=args.code_folder,
                 res_folder=args.res_folder,
                 res_folder_new=args.res_folder_new,
                 mode=args.mode)

# main(code_folder="/groups/itay_mayrose/udiland/crispys_code/CRISPys",
#      res_folder="/groups/itay_mayrose/udiland/crispys_test/test_files_git",
#      res_folder_new="/groups/itay_mayrose/udiland/crispys_test/test_files_git/res",
#      mode="run_and_compare")



