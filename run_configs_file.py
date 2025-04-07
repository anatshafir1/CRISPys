import argparse

import yaml

from Amplicon_construction.AmpliconConstruction import get_amplicons
from Gene_Family_Targeting.RunFamilyTargeting import get_gene_family_amps
from Gene_Multiplexing.RunGeneMultiplexing import get_gene_multiplex_amps


def main(configs_file):
    with open(configs_file) as f:
        data = yaml.safe_load(f)

    if data["family_targeting"]:
        amps = get_gene_family_amps(data["max_amplicon_len"], data["min_primer_len"],
                                    data["target_surrounding_region"], data["cut_location"],
                                    data["annotations_path"], data["out_path"],
                                    data["genome_fasta_path"], data["num_of_alleles"],
                                    tuple(data["pams"].split(",")), data["target_len"],
                                    data["primer3_core_path"], data["max_amplicons"],
                                    data["filter_off_targets"], data["family_targeting"])
    elif data["multiplex"]:
        amps = get_gene_multiplex_amps(data["max_amplicon_len"], data["min_primer_len"],
                                       data["target_surrounding_region"], data["cut_location"],
                                       data["annotations_path"], data["out_path"],
                                       data["genome_fasta_path"], data["num_of_alleles"],
                                       tuple(data["pams"].split(",")), data["target_len"],
                                       data["primer3_core_path"], data["max_amplicons"],
                                       data["filter_off_targets"])
    else:
        amps = get_amplicons(data["max_amplicon_len"], data["min_primer_len"],
                             data["target_surrounding_region"], data["cut_location"],
                             data["annotations_path"], data["out_path"],
                             data["genome_fasta_path"], data["num_of_alleles"],
                             tuple(data["pams"].split(",")), data["target_len"],
                             data["primer3_core_path"], data["max_amplicons"],
                             data["filter_off_targets"], data["k_alleles"],
                             data["multiplex"])
    return amps


def parse_arguments(parser_obj: argparse.ArgumentParser):
    parser_obj.add_argument('--configs_file', '-c', type=str, help='configurations file path')
    arguments = parser_obj.parse_args()
    return arguments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parse_arguments(parser)
    main(configs_file=args.configs_file)
