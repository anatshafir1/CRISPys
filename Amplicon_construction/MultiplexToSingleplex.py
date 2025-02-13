
from typing import List, Dict, Tuple

from Amplicon_Obj import Amplicon_Obj
from SNP_Obj import SNP_Obj
from Target_Obj import Target_Obj, MultiplexTarget


def split_target(multiplex_target: MultiplexTarget) -> Tuple[Target_Obj, Target_Obj]:

    up_target = Target_Obj(multiplex_target.up_seq, multiplex_target.up_start, multiplex_target.up_end,
                           multiplex_target.up_targets_list[0].strand, "", multiplex_target.up_seq,
                           multiplex_target.rank + 0.1)
    down_target = Target_Obj(multiplex_target.down_seq, multiplex_target.down_start, multiplex_target.down_end,
                             multiplex_target.down_targets_list[0].strand, "", multiplex_target.down_seq,
                             multiplex_target.rank + 0.2)

    return up_target, down_target


def get_relevant_snps(gene_targets_dict: Dict[int, List[Target_Obj]], gene_snps_dict: Dict[int, List[SNP_Obj]]) -> \
        Dict[int, List[SNP_Obj]]:
    """
    loop over every target in every exon of the gene, and check if any SNPs of that exon are on the target sequence.

    :param gene_targets_dict: dictionary of exon numbers -> a list potential targets of the exon.
    :param gene_snps_dict: dictionary of exon numbers -> a list of SNPs of the exon region.
    :return:  dictionary of exon numbers -> a list potential targets of the exon without SNPs in their sequences.
    """

    def ranges_overlap(range1, range2):
        start1, end1 = range1
        start2, end2 = range2
        return start1 <= end2 and start2 <= end1

    new_snps_dict = {}
    for exon in gene_targets_dict:
        new_exon_snps_lst = []
        for target in gene_targets_dict[exon]:
            target_range = target.start_idx, target.end_idx
            for snp in gene_snps_dict[exon]:
                snp_range = snp.position, snp.position + snp.gap_length
                if ranges_overlap(target_range, snp_range):
                    continue
                else:
                    new_exon_snps_lst.append(snp)
        new_snps_dict[exon] = new_exon_snps_lst

    return new_snps_dict


def create_amplicon(target_dict: Dict[int, List[Target_Obj]], gene_snps_dict: Dict[int, List[SNP_Obj]],
                    max_amplicon_len: int, primer_length: int, distinct_alleles_num: int,
                    target_surrounding_region: int, min_amplicon_len: int, target_len: int, genome_fasta_file: str,
                    out_path: str, pams: Tuple, gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]],
                    primer3_core_path: str, n: int, amplicon_size_range: Tuple[int, int],
                    original_exon_indices_dict: Dict[str, Dict[int, int]]) -> Amplicon_Obj:

    from Amplicon_construction.AmpliconConstruction import construct_amplicons, get_candidates_scaffold_positions
    from Amplicon_construction.GetPrimers import get_primers
    relevant_gene_snps_dict = get_relevant_snps(target_dict, gene_snps_dict)
    candidate_amplicons_list = construct_amplicons(relevant_gene_snps_dict, target_dict, max_amplicon_len,
                                                   primer_length, distinct_alleles_num, target_surrounding_region,
                                                   min_amplicon_len, 0, target_len, 0)
    if len(candidate_amplicons_list[0]) == 0:
        return
    sorted_amplicon_obj = sorted(candidate_amplicons_list[0], key=lambda amp: (amp.snps_median, amp.snps_mean),
                                 reverse=True)
    candidates_scaffold_positions = get_candidates_scaffold_positions(gene_exon_regions_seqs_dict)
    amplicon_obj_with_primers = get_primers(gene_exon_regions_seqs_dict, sorted_amplicon_obj, out_path,
                                            primer3_core_path, n, amplicon_size_range,
                                            distinct_alleles_num, target_surrounding_region, 1,
                                            genome_fasta_file, pams, candidates_scaffold_positions,
                                            original_exon_indices_dict, max_amplicon_len, gene_snps_dict, 0,
                                            0, primer_length, min_amplicon_len, target_len, genome_fasta_file)
    if len(amplicon_obj_with_primers[0]) == 0:
        return
    return amplicon_obj_with_primers[0][0]


def construct_singleplex_candidates(multiplex_targets: List[MultiplexTarget], gene_snps_dict: Dict[int, List[SNP_Obj]],
               max_amplicon_len: int, primer_length: int, distinct_alleles_num: int,
               target_surrounding_region: int, min_amplicon_len: int, k: int, target_len: int,
               multiplex: int):

    from AmpliconConstruction import construct_amplicons

    singleplex_candidate_amplicons = []
    for target in multiplex_targets:
        up_found, down_found = False, False
        exon_num = target.exon_num
        up_target, down_target = split_target(target)
        up_target_dict = {exon_num: [up_target]}
        relevant_snps_dict = {exon_num: gene_snps_dict[exon_num]}
        up_candidate_singleplex_amplicons = construct_amplicons(relevant_snps_dict, up_target_dict, max_amplicon_len, primer_length,
                                                                distinct_alleles_num, target_surrounding_region, min_amplicon_len,
                                                                k, target_len, multiplex)
        if len(up_candidate_singleplex_amplicons[0]) == 0:  # skip to next target since further search is irrelevant
            continue
        else:
            up_found = True
        down_target_dict = {exon_num: [down_target]}
        down_candidate_singleplex_amplicons = construct_amplicons(relevant_snps_dict, down_target_dict, max_amplicon_len,
                                                                primer_length, distinct_alleles_num, target_surrounding_region,
                                                                min_amplicon_len, k, target_len, multiplex)
        if len(down_candidate_singleplex_amplicons[0]) == 0:
            continue
        else:
            down_found = True
        if up_found and down_found:
            singleplex_candidate_amplicons.extend(up_candidate_singleplex_amplicons[0])
            singleplex_candidate_amplicons.extend(down_candidate_singleplex_amplicons[0])

    return singleplex_candidate_amplicons


def singleplex(candidate_multiplex_amplicon: Amplicon_Obj, gene_snps_dict: Dict[int, List[SNP_Obj]],
               max_amplicon_len: int, primer_length: int, distinct_alleles_num: int, target_surrounding_region: int,
               min_amplicon_len: int, target_len: int, genome_fasta_file: str, out_path: str, pams: Tuple,
               gene_exon_regions_seqs_dict: Dict[int, List[Tuple[str, str]]], primer3_core_path: str, n: int,
               amplicon_size_range: Tuple[int, int], original_exon_indices_dict: Dict[str, Dict[int, int]],
               ) -> List[Amplicon_Obj]:
    """
    :param candidate_multiplex_amplicon:
    :param gene_snps_dict:
    :param max_amplicon_len: maximum length of the amplicon.
    :param primer_length: minimum length of the primer sequence.
    :param distinct_alleles_num: number of distinct alleles of the gene.
    :param target_surrounding_region: buffer regions around sgRNA target (upstream and downstream) where primers are not allowed.
    :param min_amplicon_len: minimum length of the amplicon, defined by user.
    :param target_len: number of nucleotides in sgRNA target: PAM + protospacer.
    :param genome_fasta_file: path to input FASTA format file of the genome.
   :param out_path: path to output directory where algorithm results will be saved.
   :param pams: tuple of PAM sequences of the Cas protein in use.
   :param gene_exon_regions_seqs_dict: Dictionary of exon number -> list tuples of (scaffold_ID, allele sequence).
   :param primer3_core_path: A string of the full path of the primer3 core file.
   :param n: desired maximum number of amplicons to return.
   :param amplicon_size_range: tuple of minimum amplicon size and maximum amplicon size.
   :param original_exon_indices_dict: Dictionary of scaffold ID to dictionary of exon numbers after filtering to their
    original numbers in the annotations.
   :return:
   """
    print("attempting single amplicons construction for current sgrna pair".upper().center(40, "#"))
    up_target, down_target = split_target(candidate_multiplex_amplicon.target)
    current_exon_snps_dict = {
        candidate_multiplex_amplicon.exon_num: gene_snps_dict[candidate_multiplex_amplicon.exon_num]}
    up_target_dict = {candidate_multiplex_amplicon.exon_num: [up_target]}
    down_target_dict = {candidate_multiplex_amplicon.exon_num: [down_target]}
    up_target_amplicon = create_amplicon(up_target_dict, current_exon_snps_dict, max_amplicon_len, primer_length,
                                         distinct_alleles_num, target_surrounding_region, min_amplicon_len,
                                         target_len, genome_fasta_file, out_path, pams, gene_exon_regions_seqs_dict,
                                         primer3_core_path, n, amplicon_size_range, original_exon_indices_dict)
    if up_target_amplicon is None:
        return []
    down_target_amplicon = create_amplicon(down_target_dict, current_exon_snps_dict, max_amplicon_len, primer_length,
                                           distinct_alleles_num, target_surrounding_region, min_amplicon_len,
                                           target_len, genome_fasta_file, out_path, pams, gene_exon_regions_seqs_dict,
                                           primer3_core_path, n, amplicon_size_range, original_exon_indices_dict)
    if down_target_amplicon is None:
        return []
    up_target_amplicon.rank += 0.1
    down_target_amplicon.rank += 0.2
    return [up_target_amplicon, down_target_amplicon]
