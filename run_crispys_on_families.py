import os
import time
from test_crispys import createHeaderJob


def create_crispys_command(code_path: str, fam_fasta_path: str, fam_dir_path: str,
                           genes_of_interest_file: str = "None",
                           output_name: str = "CRISPys", algorithm: str = "default",
                           where_in_gene: int = 0.8, use_thr: int = 1, omega: int = 0.43,
                           off_scoring_function: str = "cfd",
                           on_scoring_function: str = "default", start_with_g: int = 0,
                           internal_node_candidates: int = 10,
                           max_target_polymorphic_sites: int = 12, pams: int = 0, singletons: int = 0,
                           slim_output: int = False,
                           set_cover: int = False, desired_genes_fraction_threshold: float = -1.0) -> str:
    """
    This function creates a string for the CRISPys command.
    Args:
        code_path: path to CRISPys code (Stage0)
        fam_dir_path: path to family directory
        fam_fasta_path: path to family fasta
        fam_fasta_path: input text file output_path of gene names and their sequences (or their exons sequences) as lines
        fam_dir_path: the output_path to the directory in which the output files will be written
        output_name: the name that would be given to the crispys output.
        genes_of_interest_file: path to a txt file consisting of a "gene" column with genes of interest.
        algorithm: the type of the algorithm run - with gene homology or without
        where_in_gene: ignore targets sites downstream to the fractional part of the gene
        use_thr:
        omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
        off_scoring_function: off target scoring function
        on_scoring_function: on target scoring function
        start_with_g: defines whether target sites are obligated to start with a G codon
        internal_node_candidates: number of sgRNAs designed for each homology subgroup
        max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
        pams: the pams by which potential sgRNA target sites will be searched
        singletons: optional choice to include singletons (sgRNAs that target only 1 gene) in the results
        slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
        set_cover: if True will output the minimal amount of guides that will capture all genes
        desired_genes_fraction_threshold: If a list of genes of interest was entered: the minimal fraction of genes
        of interest. CRISPys will ignore internal nodes with lower or equal fraction of genes of interest.
    Returns: The CRISPys command as a string
    """
    command = f"python {code_path} "
    command += f"{fam_fasta_path} "
    command += f"{fam_dir_path} "
    command += f"{output_name} "
    command += f"--genes_of_interest_file {genes_of_interest_file} "
    command += f"--alg {algorithm} "
    command += f"--where_in_gene {where_in_gene} "
    command += f"--use_thr {use_thr} "
    command += f"--omega {omega} "
    command += f"--off_scoring_function {off_scoring_function} "
    command += f"--on_scoring_function {on_scoring_function} "
    command += f"--start_with_g {start_with_g} "
    command += f"--internal_node_candidates {internal_node_candidates} "
    command += f"--max_target_polymorphic_sites {max_target_polymorphic_sites} "
    command += f"--pams {pams} "
    command += f"--singletons {singletons} "
    command += f"--slim_output {slim_output} "
    command += f"--set_cover {set_cover} "
    command += f"--desired_genes_fraction_threshold {desired_genes_fraction_threshold} "
    return command


def run(code_path: str, main_folder_path: str, genes_of_interest_file: str = "None", ncpu: int = 1, mem: int = 8,
        queue="itaym",
        output_name: str = "CRISPys", include_family_name_in_output: int = 1, algorithm: str = "default",
        where_in_gene: int = 0.8, use_thr: int = 1, omega: int = 0.43, off_scoring_function: str = "cfd",
        on_scoring_function: str = "default", start_with_g: int = 0, internal_node_candidates: int = 10,
        max_target_polymorphic_sites: int = 12, pams: int = 0, singletons: int = 0, slim_output: int = 0,
        set_cover: int = 0, desired_genes_fraction_threshold: float = -1.0):
    """
    A wrapper function to run CRISPys on the cluster for multiple folders.
    Args:
        include_family_name_in_output: if set to 1: adds the family name to each output.
        ncpu: The number of cores to use in each job
        mem: amount of memory for each job
        queue: Name of queue
        code_path: The path to the CRISPys code (Stage0)
        main_folder_path: The path to the main folder. the main folder contains a list of directories for each family.
        CRISPys args:
        output_name: the name that would be given to the crispys output.
        genes_of_interest_file: path to a txt file consisting of a "gene" column with genes of interest.
        algorithm: the type of the algorithm run - with gene homology or without
        where_in_gene: ignore targets sites downstream to the fractional part of the gene
        use_thr:
        omega: threshold of targeting propensity of a gene by a considered sgRNA (see article p. 4)
        off_scoring_function: off target scoring function
        on_scoring_function: on target scoring function
        start_with_g: defines whether target sites are obligated to start with a G codon
        internal_node_candidates: number of sgRNAs designed for each homology subgroup
        max_target_polymorphic_sites: the maximal number of possible polymorphic sites in a target
        pams: the pams by which potential sgRNA target sites will be searched
        singletons: optional choice to include singletons (sgRNAs that target only 1 gene) in the results
        slim_output: optional choice to store only 'res_in_lst' as the result of the algorithm run
        set_cover: if True will output the minimal amount of guides that will capture all genes
        desired_genes_fraction_threshold:
    Returns:
        Runs CRISPys-CRUNCH on each folder
    """
    families = os.listdir(main_folder_path)
    family_output_name = output_name
    for i, family in enumerate(families):
        if i % 700 == 0 and i > 0:
            time.sleep(153)
        fam_dir_path = os.path.join(main_folder_path, family)
        if os.path.isdir(fam_dir_path) and not family.startswith("."):
            if include_family_name_in_output:
                family_output_name = f"{family}_{output_name}"
            fam_fasta_path = os.path.join(fam_dir_path, f"{family}.fa")
            header = createHeaderJob(fam_dir_path, job_name=family_output_name, ncpu=ncpu, mem=mem, queue=queue)
            sh_file = os.path.join(fam_dir_path, f"crispys_{family_output_name}.sh")
            command = create_crispys_command(code_path=code_path, fam_fasta_path=fam_fasta_path,
                                             fam_dir_path=fam_dir_path,
                                             genes_of_interest_file=genes_of_interest_file,
                                             output_name=family_output_name,
                                             algorithm=algorithm, where_in_gene=where_in_gene, use_thr=use_thr,
                                             omega=omega, off_scoring_function=off_scoring_function,
                                             on_scoring_function=on_scoring_function, start_with_g=start_with_g,
                                             internal_node_candidates=internal_node_candidates,
                                             max_target_polymorphic_sites=max_target_polymorphic_sites, pams=pams,
                                             singletons=singletons, slim_output=slim_output, set_cover=set_cover,
                                             desired_genes_fraction_threshold=desired_genes_fraction_threshold)
            with open(sh_file, "w") as f:
                f.write(f"{header}\n{command}")
            os.system(f"qsub {sh_file}")


if __name__ == '__main__':
    run(main_folder_path="/groups/itay_mayrose/caldararu/crispys_arabidopsis/families/",
        code_path="/groups/itay_mayrose/caldararu/tmp/crispys/Stage0.py", include_family_name_in_output=True,
        genes_of_interest_file="/groups/itay_mayrose/caldararu/crispys_arabidopsis/genes_of_interest.txt",
        desired_genes_fraction_threshold=0.0, algorithm="gene_homology", slim_output=1, off_scoring_function="moff",
        omega=0.22, output_name="moff_0.22", internal_node_candidates=200, singletons=0, mem=16, ncpu=1,
        max_target_polymorphic_sites=12, set_cover=0)