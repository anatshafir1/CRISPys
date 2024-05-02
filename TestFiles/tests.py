from Amplicon_construction.run_configs_file import main

# max_amplicon_len = 1
# min_primer_len = 18
# cut_location = 7  # number of nucleotides upstream to the PAM (depending on the PAM) sequence where the Cas should cut (negative number if downstream)
# annotations_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/CONSTANS1/CONSTANS1_gene_annotations.txt"
# # annotations_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/PDS/PDS_annotations.txt"
# # annotations_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/ADH1/ADH1_corrected_annotations.txt"
# out_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/CONSTANS1"
# # out_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/PDS"
# # out_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Manual_amplicon/ADH1"
# genome_fasta_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/Banana_GAL.Phased_Scaffolds.fasta"
# num_of_alleles = 2
# pams = ("AGG", "TGG", "CGG", "GGG")
# target_len = 23
#
# primer3_core_path = "/groups/itay_mayrose/anatshafir1/CRISPR_projects/primer3/primer3/src/primer3_core"
# primer3_env_path = "/groups/itay_mayrose/anatshafir1/miniconda3/envs/primer3/"
# parameters_file_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/param_primer1"
# in_path = "/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design"
#
# max_amplicons = 5
#
# whatever = get_amplicons(max_amplicon_len, min_primer_len, cut_location, annotations_path, out_path, genome_fasta_path,
#                          num_of_alleles, pams, target_len, primer3_core_path, primer3_env_path, parameters_file_path, in_path, max_amplicons)

main("/groups/itay_mayrose/josefbrook/projects/sgRNA_Polyploids_Design/sgRNA_polyploids_design_code/Amplicon_construction/configs_file.yml")
