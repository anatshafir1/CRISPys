import pickle

from Gene_Family_Targeting.CreatesgRNAamplicons import create_single_sgrna_amplicons, create_multiplex_sgrna_amplicons
from Gene_Family_Targeting.GetGeneFamilySNPs import get_snps_dict
# from CRISPys_master.Stage0 import CRISPys_main
from Gene_Family_Targeting.GetGeneFamilySequences import get_sequences_dict, create_crispys_input_fasta
from Gene_Family_Targeting.GetGeneFamilyTargets import get_genes_single_targets_dict, \
    get_sgrnas_from_crispys_output
from globals import OFF_TARGET_FILTER_CUTOFF, MAX_POLYMORPHIC_SITES, REGION_OF_GENE_TO_CUT

primer3_core_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/sgRNA_polyploids_design_code/Amplicon_construction/primer3/src/primer3_core"
annotations_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/output/test/crispys_test/test_MLO_annotations.txt"
genome_fasta = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Banana_GAL.Phased_Scaffolds.fasta"
output_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/output/test/crispys_test/"
output_name = "crispys_out"
amplicon_ranges = [(200, 300), (300, 500), (500, 1000)]

genes_of_interest_file = 'None'
alg = "default"
where_in_gene = REGION_OF_GENE_TO_CUT
omega = 0.5
off_scoring_function = "moff"
on_scoring_function = "default"
start_with_g = 0
internal_node_candidates = 10
max_target_polymorphic_sites = MAX_POLYMORPHIC_SITES
pams = 0
slim_output = 0
set_cover = 0
min_desired_genes_fraction = -1.0
singletons = 0
singletons_on_target_function = "ucrispr"
number_of_singletons = 50
max_gap_distance = 3
export_tree = 0
run4chips = 0

# crispys_input_fasta = create_crispys_input_fasta(annotations_path, output_path, genome_fasta)

# res_grnas = CRISPys_main(crispys_input_fasta, output_path, output_name, genes_of_interest_file, alg, where_in_gene, omega,
#                    off_scoring_function, on_scoring_function, start_with_g, internal_node_candidates,
#                    max_target_polymorphic_sites, pams, slim_output, set_cover, min_desired_genes_fraction,
#                    singletons, singletons_on_target_function, number_of_singletons, max_gap_distance, export_tree,
#                    run4chips)

with open("/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/output/test/crispys_test/crispys_out.p",
          "rb") as file:
    res_grnas = pickle.load(file)

genes_exons_seq_dict = get_sequences_dict(300, 18, 20, 7, annotations_path,
                                          output_path, genome_fasta)

genes_snps_dict = get_snps_dict(genes_exons_seq_dict, 3, 18)

# genes_targets_dict = get_genes_single_targets_dict(res_grnas[0], genes_exons_seq_dict, 300, 18,
#                                                     7, 20)


sgrna_pairs_dict, sorted_sgrna_pairs_list = get_sgrnas_from_crispys_output(res_grnas[0],
                                                                           genes_exons_seq_dict, 1)

# sgrna_amplicons_dict, sgrna_seq_to_sgrna_dict = create_single_sgrna_amplicons(res_grnas[0], genes_exons_seq_dict, genes_snps_dict, genes_targets_dict, 300,
#                                           18, 3, 20, 200,
#                                           23, output_path, primer3_core_path, 5, amplicon_ranges,
#                                           1, 0, genome_fasta, ("CGG", "AGG", "GGG", "TGG"))

sgrna_amplicons_dict, sgrna_seq_to_sgrna_dict = create_multiplex_sgrna_amplicons(
    sorted_sgrna_pairs_list, genes_exons_seq_dict, genes_snps_dict, 300,
    18, 3, 20, 200,
    23, output_path, primer3_core_path, 1, amplicon_ranges,
    1, 7)
