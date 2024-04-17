import subprocess
from typing import Tuple, List, Dict

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj
from Amplicon_construction.Primers_Obj import Primers_Obj


def run_primer3(env_name: str, primer3_core: str, parameters_file: str) -> str:
    try:
        # activate the primer3 environment and run primer3
        result = subprocess.run(f"source activate {env_name} && {primer3_core} {parameters_file}", shell=True,
                                capture_output=True, text=True, check=True)
        output = result.stdout
        return output

    except subprocess.CalledProcessError as e:
        print(f"Error activating Conda environment: {env_name}")
        print(e)


def create_param_file(in_path: str, seq_id: str, seq: str, seg_target: str, product_range: str, excluded: str) -> None:
    params = [f"SEQUENCE_ID={seq_id}",
              f"SEQUENCE_TEMPLATE={seq}",
              f"SEQUENCE_TARGET={seg_target}",
              f"PRIMER_TASK=generic",
              f"PRIMER_PICK_LEFT_PRIMER=1",
              f"PRIMER_PICK_INTERNAL_OLIGO=0",
              f"PRIMER_PICK_RIGHT_PRIMER=1",
              f"PRIMER_PRODUCT_SIZE_RANGE={product_range}",
              f"SEQUENCE_EXCLUDED_REGION={excluded}",
              f"PRIMER_EXPLAIN_FLAG=1",
              "="]

    with open(in_path + "/param_primer1", "w") as param_file:
        params_str = '\n'.join(map(str, params))
        param_file.write(params_str)


def handle_primer3_output(output: str) -> Primers_Obj:
    primer_penalty = 0
    left_sequence = ""
    right_sequence = ""
    left_start_idx = 0
    right_start_idx = 0
    left_length = 0
    right_length = 0
    output_lst = output.split("\n")
    for i in range(5):
        for line in output_lst:
            if line.startswith(f"PRIMER_PAIR_{i}_PENALTY"):
                primer_penalty = float(line.split("=")[1])
            elif line.startswith(f"PRIMER_LEFT_{i}_SEQUENCE"):
                left_sequence = line.split("=")[1]
            elif line.startswith(f"PRIMER_RIGHT_{i}_SEQUENCE"):
                right_sequence = line.split("=")[1]
            elif line.startswith(f"PRIMER_LEFT_{i}="):
                left_idx_len = line.split("=")[1]
                left_start_idx = int(left_idx_len.split(",")[0])
                left_length = int(left_idx_len.split(",")[1])
            elif line.startswith(f"PRIMER_RIGHT_{i}="):
                right_idx_len = line.split("=")[1]
                right_start_idx = int(right_idx_len.split(",")[0])
                right_length = int(right_idx_len.split(",")[1])

    primers = Primers_Obj(primer_penalty, left_sequence, right_sequence, left_start_idx, left_length, right_start_idx, right_length)
    return primers


def modify_primer3_input(exon_region_seq: str, candidate_amplicon: Amplicon_Obj):

    seq_id = "gene_id"
    seq = exon_region_seq[candidate_amplicon.start_idx:candidate_amplicon.end_idx + 1]
    seg_target = f"{candidate_amplicon.target.start_idx},{candidate_amplicon.target.length}"
    product_range = "200-300"
    excluded_ranges = ""
    for snp in candidate_amplicon.snps:
        excluded_ranges += f"{snp.position_in_sequence},1 "
    return seq_id, seq, seg_target, product_range, excluded_ranges[:-1]


def get_primers(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]], sorted_candidate_amplicons: List[Amplicon_Obj],
                in_path: str, primer3_env_path: str, primer3_core_path: str, parameters_file_path: str, n: int) -> \
                List[Amplicon_Obj]:
    """


    :param gene_exon_regions_seqs_dict:
    :param sorted_candidate_amplicons:
    :param in_path:
    :param primer3_env_path:
    :param primer3_core_path:
    :param parameters_file_path:
    :param n:
    :return:
    """
    amplicons = []
    for candidate_amplicon in sorted_candidate_amplicons:
        seq_id, seq, seg_target, product_range, excluded_ranges = modify_primer3_input(gene_exon_regions_seqs_dict[1][0][1], candidate_amplicon)
        create_param_file(in_path, seq_id, seq, seg_target, product_range, excluded_ranges)
        primer3_res = run_primer3(primer3_env_path, primer3_core_path, parameters_file_path)
        if "PRIMER_LEFT_0" in primer3_res:
            primers = handle_primer3_output(primer3_res)
            start_idx = primers.left_start_idx
            end_idx = primers.right_start_idx + primers.right_length
            snps_median = candidate_amplicon.snps_median
            snps_mean = candidate_amplicon.snps_mean
            target = candidate_amplicon.target
            snps = candidate_amplicon.snps
            amplicon = Amplicon_Obj(start_idx, end_idx, snps_median, snps_mean, target, snps, primers)
            if len(amplicons) < n:
                amplicons.append(amplicon)
            else:
                break
    return amplicons
