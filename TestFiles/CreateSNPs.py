import GetSNPs_2

seq1 = "aAcaGatgaa".lower()
seq2 = "aAcaCatgaa".lower()
seq3 = "aTcaCatgaa".lower()
seq4 = "aTcaCatgaa".lower()
scaf1 = ">scaffold1:212074-212740(-)"
scaf2 = ">scaffold2:18465316-18465982(-)"
scaf3 = ">scaffold3:14465465-14466131(-)"
scaf4 = ">scaffold4:14465465-14466131(-)"

gene_sequences_dict = {1: [(scaf1, seq1), (scaf2, seq2), (scaf3, seq3), (scaf4, seq4)]}

snps_dict = GetSNPs_2.get_snps(gene_sequences_dict, 4, 0)
print(snps_dict)
valid = GetSNPs_2.valid_amplicon_2(0, len(snps_dict[1]) - 1, snps_dict[1], 4)
print(valid)
