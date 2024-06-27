import subprocess
from typing import Tuple, List, Dict

from Amplicon_Obj import Amplicon_Obj
from FindOffTargets import get_off_targets
from SNP_Obj import SNP_Obj
from Target_Obj import Target_Obj
from Primers_Obj import Primers_Obj


def run_primer3(primer3_core: str, parameters_file: str) -> str:
    try:
        # run primer3
        result = subprocess.run(f"{primer3_core} {parameters_file}", shell=True,
                                capture_output=True, text=True, check=True)
        output = result.stdout
        return output

    except subprocess.CalledProcessError as e:
        print(f"Error running {primer3_core}")
        print(e)


def create_param_file(parameters_file_path: str, seq_id: str, seq: str, seg_target: str, product_range: str,
                      excluded: str) -> None:
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

    with open(parameters_file_path, "w") as param_file:
        params_str = '\n'.join(map(str, params))
        param_file.write(params_str)


def handle_primer3_output(output: str) -> Primers_Obj:
    """

    :param output:
    :return:
    """
    primer_penalty = 0
    left_sequence = ""
    right_sequence = ""
    left_start_idx = 0
    right_start_idx = 0
    left_tm = 0.0
    right_tm = 0.0
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
            elif line.startswith(f"PRIMER_RIGHT_{i}="):
                right_idx_len = line.split("=")[1]
                right_start_idx = int(right_idx_len.split(",")[0])
            elif line.startswith(f"PRIMER_LEFT_{i}_TM="):
                left_tm = line.split("=")[1]
            elif line.startswith(f"PRIMER_RIGHT_{i}_TM="):
                right_tm = line.split("=")[1]
    primers = Primers_Obj(primer_penalty, left_sequence, right_sequence, left_start_idx, left_tm, right_start_idx,
                          right_tm)
    return primers


def modify_primer3_input(exon_region_seq: str, candidate_amplicon: Amplicon_Obj, amplicon_range,
                         target_surrounding_region: int) -> Tuple[str, str, str, str, str]:
    """

    :param exon_region_seq:
    :param candidate_amplicon:
    :param amplicon_range:
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :return:
    """
    seq_id = "gene_id"
    seq = exon_region_seq[candidate_amplicon.start_idx:candidate_amplicon.end_idx + 1]
    seg_target = f"{candidate_amplicon.target.start_idx - candidate_amplicon.start_idx - target_surrounding_region},{len(candidate_amplicon.target) + 2 * target_surrounding_region}"
    product_range = f"{amplicon_range[0]}-{amplicon_range[1]}"
    excluded_ranges = f"{candidate_amplicon.snps[0].position_in_sequence - candidate_amplicon.start_idx},{candidate_amplicon.snps[-1].position_in_sequence - candidate_amplicon.snps[0].position_in_sequence + 1}"
    return seq_id, seq, seg_target, product_range, excluded_ranges


def build_amplicon(primers, gene_exon_regions_seqs_dict, candidate_amplicon, i: int, exon_num: int) -> Amplicon_Obj:
    """

    :param primers:
    :param gene_exon_regions_seqs_dict:
    :param candidate_amplicon:
    :param i:
    :param exon_num:
    :return:
    """
    exon_region_params = gene_exon_regions_seqs_dict[exon_num][i][0].split(":")
    scaffold = exon_region_params[0][1:]
    scaffold_strand = exon_region_params[1][-2:-1]
    original_exon_region_start_idx = int(exon_region_params[1][:-3].split("-")[0]) + 1
    original_exon_region_end_idx = int(exon_region_params[1][:-3].split("-")[1])

    if scaffold_strand == "+":  # Amplicon's allele on genome forward strand
        amplicon_start_idx = primers.left_start_idx + candidate_amplicon.start_idx + original_exon_region_start_idx
        amplicon_end_idx = primers.right_start_idx + candidate_amplicon.start_idx + original_exon_region_start_idx
        target_start_idx = original_exon_region_start_idx + candidate_amplicon.target.start_idx
        target_end_idx = original_exon_region_start_idx + candidate_amplicon.target.end_idx
    else:  # scaffold_strand == "-". Amplicon's allele on genome reverse strand
        amplicon_start_idx = original_exon_region_end_idx - candidate_amplicon.start_idx - primers.right_start_idx
        amplicon_end_idx = original_exon_region_end_idx - candidate_amplicon.start_idx - primers.left_start_idx
        target_start_idx = original_exon_region_end_idx - candidate_amplicon.target.end_idx
        target_end_idx = original_exon_region_end_idx - candidate_amplicon.target.start_idx

    sequence = gene_exon_regions_seqs_dict[exon_num][i][1][
               candidate_amplicon.start_idx + primers.left_start_idx: candidate_amplicon.start_idx + primers.right_start_idx + 1]
    snps_median = candidate_amplicon.snps_median
    snps_mean = candidate_amplicon.snps_mean
    target_strand = "+" if candidate_amplicon.target.strand == scaffold_strand else "-"
    new_target = Target_Obj(candidate_amplicon.target.seq, target_start_idx, target_end_idx, target_strand)
    snps = [SNP_Obj(snp.position_in_sequence, snp.different_alleles_set) for snp in candidate_amplicon.snps]
    amplicon = Amplicon_Obj(exon_num, scaffold, scaffold_strand, sequence, amplicon_start_idx, amplicon_end_idx,
                            snps_median, snps_mean,
                            new_target, snps, primers)
    amplicon.update_snps_indices(candidate_amplicon.start_idx + primers.left_start_idx)

    return amplicon


def get_primers(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                sorted_candidate_amplicons: List[Amplicon_Obj],
                out_path: str, primer3_core_path: str, n: int, amplicon_range: Tuple[int, int],
                distinct_alleles_num: int,
                target_surrounding_region: int, filter_off_targets: int, genome_chroms_path: str) -> List[Amplicon_Obj]:
    # noinspection GrazieInspection
    """


        :param gene_exon_regions_seqs_dict:
        :param sorted_candidate_amplicons:
        :param out_path: path to output directory where algorithm results will be saved
        :param primer3_core_path:
        :param n: maximum number of Amplicons to return in the results
        :param amplicon_range:
        :param distinct_alleles_num:
        :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
        :param filter_off_targets: choose whether to filter amplicons with 'strong' off-targets for their gRNAs, or return them
        in the results.
        :param genome_chroms_path: path to the directory in which the genome separated by scaffold fasta files are.
        :return:
        """
    amplicons = []

    for candidate_amplicon in sorted_candidate_amplicons:
        search_end = False
        exon_num = candidate_amplicon.exon_num
        parameters_file_path = out_path + "/param_primer"
        seq_id, seq, seg_target, product_range, excluded_ranges = modify_primer3_input(
            gene_exon_regions_seqs_dict[exon_num][0][1], candidate_amplicon, amplicon_range, target_surrounding_region)
        create_param_file(parameters_file_path, seq_id, seq, seg_target, product_range, excluded_ranges)
        primer3_res = run_primer3(primer3_core_path, parameters_file_path)
        if "PRIMER_LEFT_0" in primer3_res:  # PRIMERS FOUND
            primers = handle_primer3_output(primer3_res)

            for i in range(distinct_alleles_num):
                amplicon = build_amplicon(primers, gene_exon_regions_seqs_dict, candidate_amplicon, i, exon_num)
                if len(amplicons) < n * distinct_alleles_num:  # 50:  # up to 'n' amplicons with 'distinct_alleles_num' alleles each
                    amplicons.append(amplicon)
                else:
                    search_end = True
                    break
        else:  # NO PRIMERS FOUND
            continue

        if search_end is True:

            return amplicons
            # get_off_targets(amplicons, genome_chroms_path, out_path)
            #
            # if filter_off_targets:  # filter amplicons with 'strong' off targets
            #     new_amplicons_lst = []
            #     for amplicon in amplicons:
            #         if len(amplicon.off_targets) > 0:
            #             if amplicon.off_targets[0].score < 0.15:
            #                 new_amplicons_lst.append(amplicon)
            #     if len(new_amplicons_lst) >= n * distinct_alleles_num:
            #         return new_amplicons_lst[:n * distinct_alleles_num]
            #     else:
            #         amplicons = new_amplicons_lst
            #
            # else:  # Do not filter amplicons
            #     return amplicons
