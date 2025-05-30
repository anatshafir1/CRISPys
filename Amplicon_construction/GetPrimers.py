import re

import subprocess
from typing import Tuple, List, Dict

from Amplicon_construction.Amplicon_Obj import Amplicon_Obj, ScaffoldAmplicon
from Amplicon_construction.FindTargets import give_complementary
from Amplicon_construction.Target_Obj import Combined_Target_Obj, Target_Obj, MultiplexTarget, FamilyMultiplexTarget, \
    Family_Target_Obj
from Amplicon_construction.FindPrimerOffTargets import get_primers_off_targets
from Amplicon_construction.FindOffTargets import get_off_targets
from Amplicon_construction.SNP_Obj import SNP_Obj
from Amplicon_construction.Primers_Obj import Primers_Obj


def run_primer3(primer3_core: str, parameters_file: str) -> str:
    try:
        # run primer3
        result = subprocess.run(f"{primer3_core} {parameters_file}", shell=True, capture_output=True, text=True,
                                check=True)
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


def handle_primer3_output(output: str, exon_region_seqs: List[Tuple[str, str]]) -> Primers_Obj:
    """

    :param output:
    :param exon_region_seqs:
    :return:
    """
    primer_penalty = 0
    left_sequence = ""
    right_sequence = ""
    left_tm = 0.0
    right_tm = 0.0
    output_lst = output.split("\n")
    for line in output_lst:
        if line.startswith(f"PRIMER_PAIR_0_PENALTY"):
            primer_penalty = float(line.split("=")[1])
        elif line.startswith(f"PRIMER_LEFT_0_SEQUENCE"):
            left_sequence = line.split("=")[1]
        elif line.startswith(f"PRIMER_RIGHT_0_SEQUENCE"):
            right_sequence = line.split("=")[1]
        elif line.startswith(f"PRIMER_LEFT_0_TM="):
            left_tm = line.split("=")[1]
        elif line.startswith(f"PRIMER_RIGHT_0_TM="):
            right_tm = line.split("=")[1]
    left_start_in_aligned = ""
    right_start_in_aligned = ""
    for allele in exon_region_seqs:
        left_matches = list(re.finditer(left_sequence, allele[1]))
        if left_matches:
            seq = left_matches[0]
            left_start_in_aligned = seq.start()
        right_matches = list(re.finditer(give_complementary(right_sequence), allele[1]))
        if right_matches:
            seq = right_matches[0]
            right_start_in_aligned = seq.end() - 1
        if left_start_in_aligned and right_start_in_aligned:
            break
    primers = Primers_Obj(primer_penalty, left_sequence, right_sequence, left_start_in_aligned, left_tm,
                          right_start_in_aligned, right_tm)
    return primers


def get_gaps_dict(exon_region_seq: str, candidate_amplicon: Amplicon_Obj, gene_snps_dict: Dict[int, List[SNP_Obj]],
                  multiplex: int):
    """
    create a dictionary of all relevant objects in the aligned sequence and the number of gaps from the start of the
    sequence up to the position of the object. SNP keys will be "snp" + "_" + "{position in sequence}", e.g. snp_15.
    Target keys will be "target_" + "start"/"end" "_{position in sequence}" e.g. target_start_35.
    """
    gaps_dict = {}
    relevant_exon_num = candidate_amplicon.exon_num
    relevant_snps_lst = gene_snps_dict[relevant_exon_num]
    for snp in relevant_snps_lst:
        gaps_to_snp = exon_region_seq[:snp.position].count("-")
        gaps_dict[f"snp_{snp.position}"] = gaps_to_snp
    if multiplex:
        gaps_to_up_target_start = exon_region_seq[:candidate_amplicon.target.up_start].count("-")
        gaps_to_up_target_end = exon_region_seq[:candidate_amplicon.target.up_end].count("-")
        gaps_to_down_target_start = exon_region_seq[:candidate_amplicon.target.down_start].count("-")
        gaps_to_down_target_end = exon_region_seq[:candidate_amplicon.target.down_end].count("-")
        gaps_dict[f"up_target_start"] = gaps_to_up_target_start
        gaps_dict[f"up_target_end"] = gaps_to_up_target_end
        gaps_dict[f"down_target_start"] = gaps_to_down_target_start
        gaps_dict[f"down_target_end"] = gaps_to_down_target_end
    else:
        gaps_to_target_start = exon_region_seq[:candidate_amplicon.target.start_idx].count("-")
        gaps_to_target_end = exon_region_seq[:candidate_amplicon.target.end_idx].count("-")
        gaps_dict[f"target_start"] = gaps_to_target_start
        gaps_dict[f"target_end"] = gaps_to_target_end
    return gaps_dict


def modify_primer3_input(exon_region_seq: str, candidate_amplicon: Amplicon_Obj, amplicon_range: Tuple[int, int],
                         target_surrounding_region: int, gene_snps_dict: Dict[int, List[SNP_Obj]],
                         multiplex: int) -> Tuple[str, str, str, str, str, Dict]:
    """

    :param exon_region_seq:
    :param candidate_amplicon:
    :param amplicon_range:
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed
    :param gene_snps_dict:
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :return:
    """
    gaps_dict = get_gaps_dict(exon_region_seq, candidate_amplicon, gene_snps_dict, multiplex)
    # calculate target start and end indices. target here is the area of the sequence around which Primer3 should search for primers.
    ungapped_first_snp_pos = candidate_amplicon.snps[0].position - gaps_dict[
        f"snp_{candidate_amplicon.snps[0].position}"]
    ungapped_last_snp_pos = candidate_amplicon.snps[-1].position - gaps_dict[
        f"snp_{candidate_amplicon.snps[-1].position}"]
    if multiplex:
        primer_search_left_boundary_idx = min(ungapped_first_snp_pos,
                                              candidate_amplicon.target.up_start - target_surrounding_region -
                                              gaps_dict["up_target_start"])
        primer_search_right_boundary_idx = max(ungapped_last_snp_pos,
                                               candidate_amplicon.target.down_end + target_surrounding_region -
                                               gaps_dict["down_target_end"])
    else:
        primer_search_left_boundary_idx = min(ungapped_first_snp_pos,
                                              candidate_amplicon.target.start_idx - target_surrounding_region -
                                              gaps_dict[f"target_start"])
        primer_search_right_boundary_idx = max(ungapped_last_snp_pos,
                                               candidate_amplicon.target.end_idx + target_surrounding_region -
                                               gaps_dict[f"target_end"])

    primers_target_seq_len = primer_search_right_boundary_idx - primer_search_left_boundary_idx

    seq_id = "gene_id"
    seq = exon_region_seq.replace("-", "")
    primers_target_seq = f"{primer_search_left_boundary_idx + 1},{primers_target_seq_len}"
    product_range = f"{amplicon_range[0]}-{amplicon_range[1]}"
    excluded_ranges = ""
    for snp in gene_snps_dict[candidate_amplicon.exon_num]:
        if snp.position > 0:
            snp_key = f"snp_{snp.position}"
            if snp.position - gaps_dict[snp_key] < len(seq):
                excluded_ranges += f"{snp.position - gaps_dict[snp_key]},1 "
    return seq_id, seq, primers_target_seq, product_range, excluded_ranges[:-1], gaps_dict


def build_amplicon(primers: Primers_Obj, allele_seq_tup: Tuple[str, str], candidate_amplicon: Amplicon_Obj,
                   original_exon_indices_dict: Dict[str, Dict[int, int]], k: int, multiplex: int,
                   gaps_dict: Dict[str, int], family_targeting: int) -> ScaffoldAmplicon:
    """
    calculate amplicon parameters of current scaffold_id and create a ScaffoldAmplicon object.

    :param primers: pair of primers and their parameters of the current candidate amplicons, as Primers_Obj object.
    :param allele_seq_tup: tuple of scaffold_ID, allele sequence of current allele.
    :param candidate_amplicon:
    :param original_exon_indices_dict: Dictionary of scaffold_id ID -> dictionary of exon numbers after filtering -> exon.
    original numbers in the annotations file.
    :param k: number of alleles to target with a single gRNA.
    :param multiplex: choose whether to plan 2 sgRNA or 1 sgRNA.
    :param gaps_dict: dictionary of object in aligned sequences and the number of gaps from sequences start to their position
    :param family_targeting:
    :return:
    """
    exon_num = candidate_amplicon.exon_num
    # Extract scaffold_id ID, scaffold_id strand and exon allele start and end indices (from annotations file)
    exon_region_params = allele_seq_tup[0].split("::")
    allele_id = exon_region_params[0][1:]
    scaffold_id = allele_id.split(";")[0]
    allele_strand = exon_region_params[1][-2]
    original_exon_region_start_idx = int(exon_region_params[1].split(":")[1][:-3].split("-")[0]) + 1
    original_exon_region_end_idx = int(exon_region_params[1].split(":")[1][:-3].split("-")[1])
    # Extract sequence parts for current allele
    allele_seq = allele_seq_tup[1]
    up_target_start_idx = 0
    up_target_end_idx = 0
    down_target_start_idx = 0
    down_target_end_idx = 0
    target_start_idx = 0
    target_end_idx = 0

    if allele_strand == "+":  # Amplicon's allele on genome forward strand
        amplicon_start_idx = original_exon_region_start_idx + primers.left_start_idx
        amplicon_end_idx = original_exon_region_start_idx + primers.right_start_idx
        if multiplex:
            up_target_start_idx = original_exon_region_start_idx + candidate_amplicon.target.up_start - gaps_dict[
                f"up_target_start"]
            up_target_end_idx = original_exon_region_start_idx + candidate_amplicon.target.up_end - gaps_dict[
                f"up_target_start"]
            down_target_start_idx = original_exon_region_start_idx + candidate_amplicon.target.down_start - gaps_dict[
                f"down_target_start"]
            down_target_end_idx = original_exon_region_start_idx + candidate_amplicon.target.down_end - gaps_dict[
                f"down_target_start"]
        else:
            target_start_idx = original_exon_region_start_idx + candidate_amplicon.target.start_idx - gaps_dict[
                f"target_start"]
            target_end_idx = original_exon_region_start_idx + candidate_amplicon.target.end_idx - gaps_dict[
                f"target_start"]
    else:  # scaffold_strand == "-". Amplicon's allele on genome reverse strand
        amplicon_start_idx = original_exon_region_end_idx - primers.right_start_idx
        amplicon_end_idx = original_exon_region_end_idx - primers.left_start_idx
        if multiplex:
            up_target_start_idx = original_exon_region_end_idx - candidate_amplicon.target.up_end + gaps_dict[
                f"up_target_start"]
            up_target_end_idx = original_exon_region_end_idx - candidate_amplicon.target.up_start + gaps_dict[
                f"up_target_start"]
            down_target_start_idx = original_exon_region_end_idx - candidate_amplicon.target.down_end + gaps_dict[
                f"down_target_start"]
            down_target_end_idx = original_exon_region_end_idx - candidate_amplicon.target.down_start + gaps_dict[
                f"down_target_start"]
        else:
            target_start_idx = original_exon_region_end_idx - candidate_amplicon.target.end_idx + gaps_dict[
                f"target_start"]
            target_end_idx = original_exon_region_end_idx - candidate_amplicon.target.start_idx + gaps_dict[
                f"target_start"]

    sequence = allele_seq[primers.left_start_idx: primers.right_start_idx + 1]
    snps_median = candidate_amplicon.snps_median
    snps_mean = candidate_amplicon.snps_mean
    if family_targeting == 1:
        new_target = Family_Target_Obj(candidate_amplicon.target.seq, target_start_idx, target_end_idx, candidate_amplicon.target.strand,
                                       candidate_amplicon.target.ungapped_seq, candidate_amplicon.target.rank,
                                       candidate_amplicon.target.score, candidate_amplicon.target.cut_alleles_lst,
                                       candidate_amplicon.target.exon_num, candidate_amplicon.target.sgrna)
    elif family_targeting == 2:
        if multiplex:
            new_target = FamilyMultiplexTarget(up_target_start_idx, up_target_end_idx, candidate_amplicon.target.up_seq,
                                         down_target_start_idx, down_target_end_idx, candidate_amplicon.target.down_seq,
                                         round(candidate_amplicon.target.multiplex_score, 4),
                                         candidate_amplicon.target.up_targets_list,
                                         candidate_amplicon.target.down_targets_list, candidate_amplicon.target.exon_num, candidate_amplicon.target.rank,
                                         candidate_amplicon.target.up_target_strand, candidate_amplicon.target.down_target_strand,
                                         candidate_amplicon.target.up_sgrna, candidate_amplicon.target.down_sgrna)
        else:
            new_target = Family_Target_Obj(candidate_amplicon.target.seq, target_start_idx, target_end_idx,
                                           candidate_amplicon.target.strand,
                                           candidate_amplicon.target.ungapped_seq, candidate_amplicon.target.rank,
                                           candidate_amplicon.target.score, candidate_amplicon.target.cut_alleles_lst,
                                           candidate_amplicon.target.exon_num, candidate_amplicon.target.sgrna)
    elif k > 0:  # Tool 2 in use
        new_target = Combined_Target_Obj(target_start_idx, target_end_idx, candidate_amplicon.target.targets_list,
                                         candidate_amplicon.target.sg_perm,
                                         candidate_amplicon.target.offscores_dict,
                                         candidate_amplicon.target.cut_alleles,
                                         candidate_amplicon.target.chosen_sg, candidate_amplicon.target.chosen_sg_score
                                         )
    elif multiplex:
        new_target = MultiplexTarget(up_target_start_idx, up_target_end_idx, candidate_amplicon.target.up_seq,
                                     down_target_start_idx, down_target_end_idx, candidate_amplicon.target.down_seq,
                                     round(candidate_amplicon.target.multiplex_score, 4),
                                     candidate_amplicon.target.up_targets_list,
                                     candidate_amplicon.target.down_targets_list, candidate_amplicon.target.exon_num, candidate_amplicon.target.rank)

    else:
        target_strand = "+" if candidate_amplicon.target.strand == allele_strand else "-"
        new_target = Target_Obj(candidate_amplicon.target.seq, target_start_idx, target_end_idx, target_strand, scaffold_id,
                                candidate_amplicon.target.ungapped_seq, candidate_amplicon.target.rank, candidate_amplicon.target.score,
                                candidate_amplicon.target.exon_num, candidate_amplicon.target.allele_id)

    snps = [SNP_Obj(snp.position, snp.alleles_sets_lst) for snp in candidate_amplicon.snps]
    orig_exon_num = original_exon_indices_dict[allele_id][exon_num]
    scaffold_amplicon = ScaffoldAmplicon(scaffold_id, allele_strand, sequence, exon_num, amplicon_start_idx,
                                         amplicon_end_idx, snps_median, snps_mean, new_target, snps, primers, orig_exon_num, allele_id=allele_id)
    scaffold_amplicon.off_targets = candidate_amplicon.off_targets
    scaffold_amplicon.update_snps_indices(primers.left_start_idx)

    return scaffold_amplicon


def get_candidate_primers(curr_rank_amp: Amplicon_Obj, gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                          out_path: str, amplicon_range: Tuple[int, int], target_surrounding_region: int,
                          gene_snps_dict: Dict[int, List[SNP_Obj]], multiplex: int, primer3_core_path: str,
                          distinct_alleles_num: int, original_exon_indices_dict: Dict[str, Dict[int, int]], k: int,
                          family_targeting: int):
    exon_num = curr_rank_amp.exon_num
    exon_region_seqs = gene_exon_regions_seqs_dict[exon_num]
    parameters_file_path = out_path + "/param_primer"
    # primers should be same for all alleles therefore any allele sequence is acceptable (exon_region_seqs[0][1])
    primer3_input = modify_primer3_input(exon_region_seqs[0][1], curr_rank_amp, amplicon_range,
                                         target_surrounding_region, gene_snps_dict, multiplex)
    seq_id, seq, seg_target, product_range, excluded_ranges, gaps_dict = primer3_input
    create_param_file(parameters_file_path, seq_id, seq, seg_target, product_range, excluded_ranges)
    primer3_res = run_primer3(primer3_core_path, parameters_file_path)

    if "PRIMER_LEFT_0" in primer3_res:  # PRIMERS FOUND
        primers = handle_primer3_output(primer3_res, exon_region_seqs)
        curr_rank_amp.primers = primers
        for i in range(distinct_alleles_num):  # Build amplicon for every allele
            allele_seq_tup = exon_region_seqs[i]
            scaffold_amplicon = build_amplicon(primers, allele_seq_tup, curr_rank_amp,
                                               original_exon_indices_dict, k, multiplex, gaps_dict, family_targeting)
            curr_rank_amp.scaffold_amplicons[scaffold_amplicon.allele_id] = scaffold_amplicon


def get_primers(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                sorted_candidate_amplicons: List[Amplicon_Obj], out_path: str, primer3_core_path: str, n: int,
                amplicon_range: Tuple[int, int], distinct_alleles_num: int, target_surrounding_region: int,
                filter_off_targets: int, genome_fasta_path: str, pams: Tuple,
                candidates_scaffold_positions: Dict[str, List[Tuple[int, int]]],
                original_exon_indices_dict: Dict[str, Dict[int, int]], max_amplicon_len: int,
                gene_snps_dict: Dict[int, List[SNP_Obj]], k: int, min_amplicon_len: int) -> List[Amplicon_Obj]:
    # noinspection GrazieInspection
    """

    :param gene_exon_regions_seqs_dict: Dictionary of exon number -> list tuples of (scaffold_ID, allele sequence).
    :param sorted_candidate_amplicons: list of potential amplicons, sorted number of SNPs.
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
    :param candidates_scaffold_positions: Dictionary of scaffold_id ID to tuples of gene start index and gene end index.
    :param original_exon_indices_dict: Dictionary of scaffold_id ID to dictionary of exon numbers after filtering to their
    original numbers in the annotations.
    :param max_amplicon_len: maximum length of the amplicon.
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region.
    :param k: number of alleles to target with a single gRNA.
    :param min_amplicon_len:
    :return: list of results amplicons with primers.
    """
    # print("finding primers".upper().center(40, "#"))
    amplicons = []
    for candidate_amplicon in sorted_candidate_amplicons:
        get_candidate_primers(candidate_amplicon, gene_exon_regions_seqs_dict, out_path, amplicon_range,
                              target_surrounding_region, gene_snps_dict, 0, primer3_core_path,
                              distinct_alleles_num, original_exon_indices_dict, k, 0)
        if (len(candidate_amplicon.scaffold_amplicons) == distinct_alleles_num and
                all(min_amplicon_len <= candidate_amplicon.scaffold_amplicons[scaf_amp].size <= max_amplicon_len for
                    scaf_amp in candidate_amplicon.scaffold_amplicons)):
            if len(amplicons) < n:  # up to 'n' amplicons
                if candidate_amplicon not in amplicons:  # check needed for some candidates had different start and end indices with potential to get different primers but same primers were found
                    amplicons.append(candidate_amplicon)
            else:
                break

        else:  # NO PRIMERS FOUND
            continue

    if not filter_off_targets:  # search for off targets
        if len(amplicons) > 0:
            get_off_targets(amplicons, genome_fasta_path, out_path, pams, candidates_scaffold_positions, k)
            get_primers_off_targets(amplicons, genome_fasta_path, out_path, candidates_scaffold_positions,
                                    max_amplicon_len, 0)
        return amplicons
    else:
        return amplicons


def get_gene_family_primers(gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                            sorted_candidate_amplicons: List[Amplicon_Obj], out_path: str, primer3_core_path: str,
                            n: int,
                            amplicon_range: Tuple[int, int], distinct_alleles_num: int, target_surrounding_region: int,
                            original_exon_indices_dict: Dict[str, Dict[int, int]], max_amplicon_len: int,
                            gene_snps_dict: Dict[int, List[SNP_Obj]], k: int,
                            min_amplicon_len: int, multiplex: int, family_targeting: int) -> List[Amplicon_Obj]:
    # print("finding primers".upper().center(40, "#"))
    amplicons = []

    for candidate_amplicon in sorted_candidate_amplicons:
        get_candidate_primers(candidate_amplicon, gene_exon_regions_seqs_dict, out_path, amplicon_range,
                              target_surrounding_region, gene_snps_dict, multiplex, primer3_core_path,
                              distinct_alleles_num, original_exon_indices_dict, k, family_targeting)
        if (len(candidate_amplicon.scaffold_amplicons) == distinct_alleles_num and
                all(min_amplicon_len <= candidate_amplicon.scaffold_amplicons[scaf_amp].size <= max_amplicon_len for
                    scaf_amp in candidate_amplicon.scaffold_amplicons)):
            if len(amplicons) < n:  # up to 'n' amplicons
                if candidate_amplicon not in amplicons:  # check needed for some candidates had different start and end indices with potential to get different primers but same primers were found
                    amplicons.append(candidate_amplicon)
            else:
                break

        else:  # NO PRIMERS FOUND
            continue
    return amplicons
