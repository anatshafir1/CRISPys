import os
from test_crispys import createHeaderJob
import logging
import argparse


def create_crispys_command(code_path: str, fam_fasta_path: str, fam_dir_path: str,
                           genes_of_interest_file: str = "None",
                           output_name: str = "CRISPys", algorithm: str = "default",
                           where_in_gene: float = 0.8, omega: float = 0.43,
                           off_scoring_function: str = "cfd",
                           on_scoring_function: str = "default", start_with_g: int = 0,
                           internal_node_candidates: int = 10,
                           max_target_polymorphic_sites: int = 12, pams: int = 0,
                           slim_output: int = 0,
                           set_cover: int = 0, min_desired_genes_fraction: float = -1.0, singletons: int = 0,
                           singletons_on_target_function: str = "ucrispr", number_of_singletons: int = 50,
                           max_gap_distance: int = 0, export_tree: int = 0, run4chips: int = 0) -> str:
    """
    This function creates a string for the CRISPys command.
    Args:
        code_path: path to CRISPys code folder
        fam_dir_path: path to family directory
        fam_fasta_path: path to family fasta
        fam_fasta_path: input text file output_path of gene names and their sequences (or their exons sequences) as lines
        fam_dir_path: the output_path to the directory in which the output files will be written
        output_name: the name that would be given to the crispys output.
        genes_of_interest_file: path to a txt file consisting of a "gene" column with genes of interest.
        algorithm: the type of the algorithm run - with gene homology or without
        where_in_gene: ignore targets sites downstream to the fractional part of the gene
        omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
        off_scoring_function: off target scoring function
        on_scoring_function: on target scoring function
        start_with_g: defines whether target sites are obligated to start with a G codon
        internal_node_candidates: number of sgRNAs designed for each homology subgroup
        max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
        pams: the pams by which potential sgRNA target sites will be searched
        slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
        set_cover: if True will output the minimal amount of guides that will capture all genes
        min_desired_genes_fraction: If a list of genes of interest was entered: the minimal fraction of genes
        of interest. CRISPys will ignore internal nodes with lower or equal fraction of genes of interest.
        singletons: select 1 to create singletons (sgRNAs candidates that target a single gene).
        number_of_singletons: the number of singletons that will be included for each gene.
        singletons_on_target_function: The on-target scoring function used for evaluating singletons.
        max_gap_distance: max_gap_distance: The maximal distance that is allowed between the genes targeted by the sgRNA
        export_tree: to output a pickle file with the gene tree object
        run4chips: do not filter sgRNAs that target genes that are not in the list

    Returns: The CRISPys command as a string
    """
    command = ''
    if singletons_on_target_function == "ucrispr" or on_scoring_function == "ucrispr":
        command += f"export DATAPATH={code_path}/uCRISPR/RNAstructure/data_tables/\n"
    command += f"python {os.path.join(code_path, 'Stage0.py')} "
    command += f"{fam_fasta_path} "
    command += f"{fam_dir_path} "
    command += f"--output_name {output_name} "
    command += f"--genes_of_interest_file {genes_of_interest_file} "
    command += f"--alg {algorithm} "
    command += f"--where_in_gene {where_in_gene} "
    command += f"--omega {omega} "
    command += f"--off_scoring_function {off_scoring_function} "
    command += f"--on_scoring_function {on_scoring_function} "
    command += f"--start_with_g {start_with_g} "
    command += f"--internal_node_candidates {internal_node_candidates} "
    command += f"--max_target_polymorphic_sites {max_target_polymorphic_sites} "
    command += f"--pams {pams} "
    command += f"--slim_output {slim_output} "
    command += f"--set_cover {set_cover} "
    command += f"--min_desired_genes_fraction {min_desired_genes_fraction} "
    command += f"--singletons {singletons} "
    command += f"--singletons_on_target_function {singletons_on_target_function} "
    command += f"--number_of_singletons {number_of_singletons} "
    command += f"--max_gap_distance {max_gap_distance} "
    command += f"--export_tree {export_tree} "
    command += f"--run4chips {run4chips}"

    return command


def contains_genes_of_interest(fam_fasta_path, set_of_genes_of_interest):
    """
    This function checks if a family has at least one gene from a list of interest.
    :param set_of_genes_of_interest: A set of genes of interest.
    :param fam_fasta_path: Path to the fasta file of a gene family.
    :return: True if at least one gene is in the set of genes of interest. Else False
    """
    if not set_of_genes_of_interest:
        return True
    with open(fam_fasta_path, 'r') as f:
        lines = f.readlines()

    fasta_gene_names_set = {gene.strip(">\n") for gene in lines if gene.startswith(">")}  # Extract genes from the fasta
    intersection_set = fasta_gene_names_set.intersection(set_of_genes_of_interest)  # check if there is an intersection
    if intersection_set:
        return True
    return False


def run(main_folder_path: str, genes_of_interest_file: str = "None", ncpu: int = 1, mem: int = 8,
        queue="itaym",
        output_name: str = "CRISPys", include_family_name_in_output: int = 1, alg: str = "default",
        where_in_gene: float = 0.8, use_thr: int = 1, omega: float = 0.43, off_scoring_function: str = "cfd",
        on_scoring_function: str = "default", start_with_g: int = 0, internal_node_candidates: int = 10,
        max_target_polymorphic_sites: int = 12, pams: int = 0, slim_output: int = 0,
        set_cover: int = 0, min_desired_genes_fraction: float = -1.0, singletons: int = 0,
        singletons_on_target_function: str = "ucrispr", number_of_singletons: int = 50, max_gap_distance: int = 3,
        check_for_genes_of_interest: bool = False, export_tree: int = 0, run4chips: int = 0):
    """
    A wrapper function to run CRISPys on the cluster for multiple folders.
    Args:
        include_family_name_in_output: if set to 1: adds the family name to each output.
        ncpu: The number of cores to use in each job
        mem: amount of memory for each job
        queue: Name of queue
        main_folder_path: The path to the main folder. the main folder contains a list of directories for each family.
        CRISPys args:
        output_name: the name that would be given to the crispys output.
        genes_of_interest_file: path to a txt file consisting of a "gene" column with genes of interest.
        alg: the type of the algorithm run - with gene homology or without
        where_in_gene (float): ignore targets sites downstream to the fractional part of the gene
        use_thr:
        omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
        off_scoring_function: off target scoring function
        on_scoring_function: on target scoring function
        start_with_g: defines whether target sites are obligated to start with a G codon
        internal_node_candidates: number of sgRNAs designed for each homology subgroup
        max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
        pams: the pams by which potential sgRNA target sites will be searched
        slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
        set_cover: if True will output the minimal amount of guides that will capture all genes
        min_desired_genes_fraction: If a list of genes of interest was entered: the minimal fraction of genes
        of interest. CRISPys will ignore internal nodes with lower or equal fraction of genes of interest.
        singletons: select 1 to create singletons (sgRNAs candidates that target a single gene).
        number_of_singletons: the number of singletons that will be included for each gene.
        singletons_on_target_function: The on-target scoring function used for evaluating singletons.
        max_gap_distance: max_gap_distance: The maximal distance that is allowed between the genes targeted by the sgRNA
        export_tree: if to output a pickle file of gene tree
        run4chips:
        check_for_genes_of_interest: optional: checks if the fasta file contains any genes of interest before submitting
        the job for that particular family
    Returns:
        :param
    """
    code_path = os.path.dirname(os.path.abspath(__file__))
    families = os.listdir(main_folder_path)
    family_output_name = output_name
    set_of_genes_of_interest = set()
    if check_for_genes_of_interest and genes_of_interest_file != "None":
        with open(genes_of_interest_file, 'r') as f:
            genes_of_interest_lines = f.readlines()
        set_of_genes_of_interest = {gene.strip() for gene in genes_of_interest_lines}
    for i, family in enumerate(families):
        fam_dir_path = os.path.join(main_folder_path, family)
        if os.path.isdir(fam_dir_path) and not family.startswith("."):
            if include_family_name_in_output:
                family_output_name = f"{family}_{output_name}"
            fam_fasta_path = os.path.join(fam_dir_path, f"{family}.fa")
            if check_for_genes_of_interest and not contains_genes_of_interest(fam_fasta_path, set_of_genes_of_interest):
                continue
            logging.debug(f"Run CRISPys on family: {family}\n")
            header = createHeaderJob(fam_dir_path, job_name=family_output_name, ncpu=ncpu, mem=mem, queue=queue)
            sh_file = os.path.join(fam_dir_path, f"crispys_{family_output_name}.sh")
            command = create_crispys_command(code_path=code_path, fam_fasta_path=fam_fasta_path,
                                             fam_dir_path=fam_dir_path,
                                             genes_of_interest_file=genes_of_interest_file,
                                             output_name=family_output_name,
                                             algorithm=alg, where_in_gene=where_in_gene,
                                             omega=omega, off_scoring_function=off_scoring_function,
                                             on_scoring_function=on_scoring_function, start_with_g=start_with_g,
                                             internal_node_candidates=internal_node_candidates,
                                             max_target_polymorphic_sites=max_target_polymorphic_sites, pams=pams,
                                             slim_output=slim_output,
                                             set_cover=set_cover,
                                             min_desired_genes_fraction=min_desired_genes_fraction,
                                             singletons=singletons, number_of_singletons=number_of_singletons,
                                             singletons_on_target_function=singletons_on_target_function,
                                             max_gap_distance=max_gap_distance,
                                             export_tree=export_tree, run4chips=run4chips)
            with open(sh_file, "w") as f:
                f.write(f"{header}\n{command}")
            os.system(f"qsub {sh_file}")
            logging.debug(f"Job submitted for {family}")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="A wrapper function to run CRISPys on the cluster for multiple folders.")

    parser.add_argument("main_folder_path", type=str,
                        help="The path to the main folder. The main folder contains a list of directories for each family.")
    parser.add_argument("--genes_of_interest_file", type=str, default="None",
                        help="Path to a txt file consisting of a 'gene' column with genes of interest.")
    parser.add_argument("--ncpu", type=int, default=1, help="The number of cores to use in each job")
    parser.add_argument("--mem", type=int, default=8, help="Amount of memory for each job")
    parser.add_argument("--queue", type=str, default="itaym", help="Name of queue")
    parser.add_argument("--output_name", type=str, default="CRISPys",
                        help="The name that would be given to the CRISPys output")
    parser.add_argument("--include_family_name_in_output", type=int, default=1,
                        help="If set to 1, adds the family name to each output")
    parser.add_argument("--alg", type=str, default="default",
                        help="The type of the algorithm run - with gene homology or without")
    parser.add_argument("--where_in_gene", type=float, default=0.8,
                        help="Ignore target sites downstream to the fractional part of the gene")
    parser.add_argument("--use_thr", type=int, default=1)
    parser.add_argument("--omega", type=float, default=0.43,
                        help="Threshold of targeting propensity of a gene by a considered sgRNA")
    parser.add_argument("--off_scoring_function", type=str, default="cfd", help="Off-target scoring function")
    parser.add_argument("--on_scoring_function", type=str, default="default", help="On-target scoring function")
    parser.add_argument("--start_with_g", type=int, default=0,
                        help="Defines whether target sites are obligated to start with a G codon")
    parser.add_argument("--internal_node_candidates", type=int, default=10,
                        help="Number of sgRNAs designed for each homology subgroup")
    parser.add_argument("--max_target_polymorphic_sites", type=int, default=12,
                        help="The maximal number of possible polymorphic sites in a target")
    parser.add_argument("--pams", type=int, default=0,
                        help="The PAMs by which potential sgRNA target sites will be searched")
    parser.add_argument("--slim_output", type=int, default=0,
                        help="Optional choice to store only 'res_in_lst' as the result of the algorithm run")
    parser.add_argument("--set_cover", type=int, default=0,
                        help="If True, will output the minimal amount of guides that will capture all genes")
    parser.add_argument("--min_desired_genes_fraction", type=float, default=-1.0,
                        help="If a list of genes of interest was entered, the minimal fraction of genes of interest. CRISPys will ignore internal nodes with lower or equal fraction of genes of interest.")
    parser.add_argument("--singletons", type=int, default=0,
                        help="Select 1 to create singletons (sgRNAs candidates that target a single gene)")
    parser.add_argument("--number_of_singletons", type=int, default=50,
                        help="The number of singletons that will be included for each gene")
    parser.add_argument("--max_gap_distance", type=int, default=3,
                        help="The maximal distance that is allowed between the genes targeted by the sgRNA")
    parser.add_argument("--check_for_genes_of_interest", action="store_true",
                        help="Optional: checks if the fasta file contains any genes of interest before submitting the job for that particular family")
    parser.add_argument("--export_tree", type=int, default=0)
    parser.add_argument("--run4chips", type=int, default=0)

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    logging.basicConfig(level=logging.DEBUG,  # Set the desired log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                        filename=os.path.join(args.main_folder_path,
                                              f"{args.output_name}_CRISPys.log"))
    run(**vars(args))