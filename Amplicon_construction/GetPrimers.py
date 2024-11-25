import subprocess
from typing import Tuple, List, Dict

from Amplicon_Obj import Amplicon_Obj, ScaffoldAmplicon
from Target_Obj import Combined_Target_Obj, Target_Obj
from FindPrimerOffTargets import get_primers_off_targets
from FindOffTargets import get_off_targets
from SNP_Obj import SNP_Obj
from Primers_Obj import Primers_Obj


def run_primer3(primer3_core: str, parameters_file: str) -> str:
    try:
        # run primer3
        result = subprocess.run(f"{primer3_core} {parameters_file}", shell=True, capture_output=True, text=True, check=True)
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


def modify_primer3_input(exon_region_seq: str, candidate_amplicon: Amplicon_Obj, amplicon_range: Tuple[int, int],
                         target_surrounding_region: int, gene_snps_dict: Dict[int, List[SNP_Obj]]) -> Tuple[str, str, str, str, str]:
    """

    :param exon_region_seq:
    :param candidate_amplicon:
    :param amplicon_range:
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param gene_snps_dict:
    :return:
    """
    # calculate target start and end indices. target here is the area of the sequence around which Primer3 should search for primers.
    target_start_idx = min(candidate_amplicon.snps[0].position,
                           candidate_amplicon.target.start_idx - target_surrounding_region) - candidate_amplicon.start_idx
    target_end_idx = max(candidate_amplicon.snps[-1].position + candidate_amplicon.snps[-1].gap_length,
                         candidate_amplicon.target.end_idx + target_surrounding_region) - candidate_amplicon.start_idx
    target_len = target_end_idx - target_start_idx

    seq_id = "gene_id"

    seq_with_gaps = exon_region_seq[candidate_amplicon.start_idx:candidate_amplicon.end_idx + 1]

    seq = seq_with_gaps.replace("-", "")  # TODO take gaps into account for all parameters
    seg_target = f"{target_start_idx},{target_len}"
    product_range = f"{amplicon_range[0]}-{amplicon_range[1]}"
    excluded_ranges = ""
    relevant_snp = [snp for snp in gene_snps_dict[candidate_amplicon.exon_num]
                    if candidate_amplicon.start_idx <= snp.position
                    and snp.position + snp.gap_length <= candidate_amplicon.end_idx]
    for snp in relevant_snp:
        excluded_ranges += f"{snp.position - candidate_amplicon.start_idx},1 "
    return seq_id, seq, seg_target, product_range, excluded_ranges[:-1]


def build_amplicon(primers: Primers_Obj, allele_seq_tup: Tuple[str, str], candidate_amplicon: Amplicon_Obj,
                   original_exon_indices_dict: Dict[str, Dict[int, int]], k: int) -> ScaffoldAmplicon:
    """
    calculate amplicon parameters of current scaffold and create a ScaffoldAmplicon object.

    :param primers: pair of primers and their parameters of the current candidate amplicons, as Primers_Obj object.
    :param allele_seq_tup: tuple of scaffold_ID, allele sequence of current allele.
    :param candidate_amplicon:
    :param original_exon_indices_dict: Dictionary of scaffold ID -> dictionary of exon numbers after filtering -> exon
    original numbers in the annotations file.
    :param k: number of alleles to target with a single gRNA.
    :return:
    """
    exon_num = candidate_amplicon.exon_num
    # Extract scaffold ID, scaffold strand and exon allele start and end indices (from annotations file)
    exon_region_params = allele_seq_tup[0].split(":")
    scaffold = exon_region_params[0][1:]
    scaffold_strand = exon_region_params[1][-2:-1]
    original_exon_region_start_idx = int(exon_region_params[1][:-3].split("-")[0]) + 1
    original_exon_region_end_idx = int(exon_region_params[1][:-3].split("-")[1])
    # Extract sequence parts for current allele
    allele_seq = allele_seq_tup[1]
    seq_to_amp_start = allele_seq[:candidate_amplicon.start_idx + primers.left_start_idx + 1]
    seq_to_target = allele_seq[:candidate_amplicon.target.start_idx + 1]
    seq_to_amp_end = allele_seq[:candidate_amplicon.start_idx + primers.right_start_idx + 1]
    gaps_to_amp_start = seq_to_amp_start.count("-")
    gaps_to_target = seq_to_target.count("-")
    gaps_to_amp_end = seq_to_amp_end.count("-")
    if scaffold_strand == "+":  # Amplicon's allele on genome forward strand
        amplicon_start_idx = original_exon_region_start_idx + candidate_amplicon.start_idx + primers.left_start_idx - gaps_to_amp_end
        amplicon_end_idx = original_exon_region_start_idx + candidate_amplicon.start_idx + primers.right_start_idx - gaps_to_amp_start
        target_start_idx = original_exon_region_start_idx + candidate_amplicon.target.start_idx - gaps_to_target
        target_end_idx = original_exon_region_start_idx + candidate_amplicon.target.end_idx - gaps_to_target
    else:  # scaffold_strand == "-". Amplicon's allele on genome reverse strand
        amplicon_start_idx = original_exon_region_end_idx - candidate_amplicon.start_idx - primers.right_start_idx + gaps_to_amp_end
        amplicon_end_idx = original_exon_region_end_idx - candidate_amplicon.start_idx - primers.left_start_idx + gaps_to_amp_start
        target_start_idx = original_exon_region_end_idx - candidate_amplicon.target.end_idx + gaps_to_target
        target_end_idx = original_exon_region_end_idx - candidate_amplicon.target.start_idx + gaps_to_target

    sequence = allele_seq[candidate_amplicon.start_idx + primers.left_start_idx: candidate_amplicon.start_idx + primers.right_start_idx + 1]
    snps_median = candidate_amplicon.snps_median
    snps_mean = candidate_amplicon.snps_mean
    if k > 0:  # Tool 2 in use
        new_target = Combined_Target_Obj(target_start_idx, target_end_idx, candidate_amplicon.target.targets_list,
                                         candidate_amplicon.target.sg_perm,
                                         candidate_amplicon.target.offscores_dict, candidate_amplicon.target.cut_alleles,
                                         candidate_amplicon.target.chosen_sg,  candidate_amplicon.target.chosen_sg_score
                                         )
    else:
        target_strand = "+" if candidate_amplicon.target.strand == scaffold_strand else "-"
        new_target = Target_Obj(candidate_amplicon.target.seq, target_start_idx, target_end_idx, target_strand)
    snps = [SNP_Obj(snp.position, snp.different_alleles_set) for snp in candidate_amplicon.snps]
    orig_exon_num = original_exon_indices_dict[scaffold][exon_num]
    scaffold_amplicon = ScaffoldAmplicon(scaffold, scaffold_strand, sequence, exon_num, amplicon_start_idx, amplicon_end_idx,
                            snps_median, snps_mean, new_target, snps, primers, orig_exon_num)
    scaffold_amplicon.off_targets = candidate_amplicon.off_targets
    scaffold_amplicon.update_snps_indices(candidate_amplicon.start_idx + primers.left_start_idx)

    return scaffold_amplicon


def filer_amplicons(amplicons: List[Amplicon_Obj]) -> List[Amplicon_Obj]:
    """
    filter amplicons with the same primers and target, and return amplicons only with longest snps list
    :param amplicons:
    :return:
    """
    filtered_set = set()
    sorted_by_snps = sorted(amplicons, key=lambda amp: len(amp.snps), reverse=True)
    for amplicon in sorted_by_snps:
        if amplicon not in filtered_set:
            filtered_set.add(amplicon)
    return list(filtered_set)


def get_primers(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                sorted_candidate_amplicons: List[Amplicon_Obj], out_path: str, primer3_core_path: str, n: int,
                amplicon_range: Tuple[int, int], distinct_alleles_num: int, target_surrounding_region: int,
                filter_off_targets: int, genome_fasta_path: str, pams: Tuple, candidates_scaffold_positions: Dict[str, Tuple[int, int]],
                original_exon_indices_dict: Dict[str, Dict[int, int]], max_amplicon_len: int, gene_snps_dict: Dict[int, List[SNP_Obj]], k: int) -> List[Amplicon_Obj]:
    # noinspection GrazieInspection
    """


    :param gene_exon_regions_seqs_dict: Dictionary of exon number -> list tuples of (scaffold_ID, allele sequence).
    :param sorted_candidate_amplicons: list of potential amplicons, sorted by median number of SNPs per allele.
    :param out_path: path to output directory where algorithm results will be saved.
    :param primer3_core_path: A string of the full path of the primer3 core file.
    :param n: maximum number of Amplicons to return in the results.
    :param amplicon_range: tuple of minimum amplicon size and maximum amplicon size.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param filter_off_targets: choose whether to filter amplicons with 'strong' off-targets for their gRNAs, or return them
    in the results.
    :param genome_fasta_path: path to the directory in which the genome fasta file.
    :param pams: tuple of PAM sequences of the Cas protein in use.
    :param candidates_scaffold_positions: Dictionary of scaffold ID to tuples of gene start index and gene end index.
    :param original_exon_indices_dict: Dictionary of scaffold ID to dictionary of exon numbers after filtering to their
    original numbers in the annotations.
    :param max_amplicon_len: maximum length of the amplicon.
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region.
    :param k: number of alleles to target with a single gRNA.
    :return: list of results amplicons with primers.
    """
    amplicons = []

    for index, candidate_amplicon in enumerate(sorted_candidate_amplicons):
        exon_num = candidate_amplicon.exon_num
        exon_region_seqs = gene_exon_regions_seqs_dict[exon_num]
        parameters_file_path = out_path + "/param_primer"
        # primers should be same for all alleles therefore any allele sequence is acceptable (exon_region_seqs[0][1])
        seq_id, seq, seg_target, product_range, excluded_ranges = modify_primer3_input(exon_region_seqs[0][1],
                                                                                       candidate_amplicon, amplicon_range, target_surrounding_region, gene_snps_dict)
        create_param_file(parameters_file_path, seq_id, seq, seg_target, product_range, excluded_ranges)
        primer3_res = run_primer3(primer3_core_path, parameters_file_path)

        if "PRIMER_LEFT_0" in primer3_res:  # PRIMERS FOUND
            primers = handle_primer3_output(primer3_res)
            candidate_amplicon.primers = primers
            for i in range(distinct_alleles_num):  # Build amplicon for every allele
                allele_seq_tup = exon_region_seqs[i]
                scaffold_amplicon = build_amplicon(primers, allele_seq_tup, candidate_amplicon, original_exon_indices_dict, k)
                candidate_amplicon.scaffold_amplicons[scaffold_amplicon.scaffold] = scaffold_amplicon
            if len(amplicons) < n:  # up to 'n' amplicons
                amplicons.append(candidate_amplicon)
            else:
                break
        else:  # NO PRIMERS FOUND
            continue
    filtered_amplicons = filer_amplicons(amplicons)
    if not filter_off_targets:  # search for off targets
        get_off_targets(filtered_amplicons, genome_fasta_path, out_path, pams, candidates_scaffold_positions, k)
        get_primers_off_targets(filtered_amplicons, genome_fasta_path, out_path, candidates_scaffold_positions, max_amplicon_len)
        return filtered_amplicons
    else:
        return filtered_amplicons
