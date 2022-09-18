from globals import createHeaderJob
import os

def run(crispys_code, main_fam_folder, alg="gene_homology", where_in_gene=0.7, use_thr=1, omega=0.43,
        scoring_function='cfd_funct', internal_node_candidate=200,
        max_target_polymorphic_sites=12, pams=0, singletons=1, slim_output=True, queue="itaym", ncpu=1, mem=16):
    families = os.listdir(main_fam_folder)
    for family in families:
        fam_path = os.path.join(main_fam_folder, family)
        fasta = family + ".txt"
        if os.path.isdir(fam_path) and not family.startswith("."):
            header = createHeaderJob(fam_path, job_name=family, ncpu=ncpu, mem=mem, queue=queue)
            command = f"cd {fam_path}\npython {crispys_code}/Stage0.py {fasta} {fam_path} --alg {alg}" \
                      f" --where_in_gene {where_in_gene} " \
                      f"-t {use_thr} -v {omega} -s {scoring_function} -i {internal_node_candidate}" \
                      f" -ps {max_target_polymorphic_sites} --pams {pams} --singletons {singletons} " \
                      f"--slim_output {slim_output}"
            with open(f"{fam_path}/Crispys.sh", "w") as f:
                f.write(f"{header}\n{command}")
            os.system(f"qsub {fam_path}/Crispys.sh")


