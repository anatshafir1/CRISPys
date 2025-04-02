from typing import Tuple

from Amplicon_construction.GetSNPs import get_snps
from Amplicon_construction.GetSequences import extract_exons_regions
from Gene_Multiplexing.CreateMultiplexAmplicons import create_multiplex_amplicons
from Gene_Multiplexing.GetMultiplexOffTargets import filter_multiplex_off_targets
from Gene_Multiplexing.GetMultiplexTargets import get_multiplex_targets


def get_gene_multiplex_amps(max_amplicon_len_category: int, primer_length: int, target_surrounding_region: int,
                            cut_location: int,
                            annotations_file_path: str, out_path: str, genome_fasta_file: str,
                            distinct_alleles_num: int,
                            pams: Tuple[str], target_len: int, primer3_core_path: str, n: int, filter_off_targets: int,
                            k: int,
                            multiplex: int):
    amplicon_ranges = [(200, 300), (300, 500), (500, 1000)]
    max_amplicon_len = max(amplicon_ranges[max_amplicon_len_category - 1])
    min_amplicon_len = min(amplicon_ranges[max_amplicon_len_category - 1])
    gene_exon_regions_seqs_dict, original_exon_indices_dict = extract_exons_regions(max_amplicon_len, primer_length,
                                                                                    target_surrounding_region,
                                                                                    cut_location,
                                                                                    annotations_file_path, out_path,
                                                                                    genome_fasta_file)
    gene_snps_dict = get_snps(gene_exon_regions_seqs_dict, distinct_alleles_num, primer_length)
    multiplex_targets_list = get_multiplex_targets(gene_exon_regions_seqs_dict, pams,
                                                   max_amplicon_len, primer_length, cut_location,
                                                   target_surrounding_region, target_len,
                                                   distinct_alleles_num)
    if filter_off_targets:
        multiplex_targets_list = filter_multiplex_off_targets(multiplex_targets_list, out_path, genome_fasta_file, pams,
                                                              gene_exon_regions_seqs_dict)
    multiplex_amplicons_dict = create_multiplex_amplicons(
            multiplex_targets_list, gene_exon_regions_seqs_dict, gene_snps_dict, original_exon_indices_dict, max_amplicon_len,
            primer_length, distinct_alleles_num, target_surrounding_region, min_amplicon_len,
            target_len, out_path, primer3_core_path, n, amplicon_ranges,
            max_amplicon_len_category, cut_location)
    return
