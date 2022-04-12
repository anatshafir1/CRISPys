__author__ = 'ItayM5'
import re ##used to be regex
import regex

#1324
#write simhing to test git commit-push-pull :

def get_sites(gene, min_length=20, max_length=20, start_with_G=False, where_in_gene = 1, PAMs= ["GG"]):
	'''
	This function is used to get CRISPR cut site sequence from an input dna sequence, is used to find site in gene/exon
	:param gene: sequence of a gene or exon
	:param min_length: minimum length of site
	:param max_length: maximum length of site
	:param start_with_G: if the guide should start with G
	:param where_in_gene: forword to this position the sgRNA are ignored
	:param PAMs: type of pam (can be [NGG] or [NGG, NAG]
	:return: a list of sites sequences that CRISPR can target
	'''
	res = []
	if len(gene) < max_length+3:
		return res
	for length in range(min_length, max_length +1):
		for i in range(len(PAMs)):
			if (start_with_G):
				SiteAndPAM = "G" + "."*length + PAMs[i]
			else:
				SiteAndPAM = "."*(length +1) + PAMs[i]
			compiled = regex.compile(SiteAndPAM)
			where_in_gene = int(len(gene)*where_in_gene)
			founds_sense = regex.findall(compiled, gene[:where_in_gene], overlapped=True)
			founds_antisense = regex.findall(compiled, give_complementary(gene)[:where_in_gene], overlapped=True)
			founds = [seq[:-3] for seq in founds_sense if 'N' not in seq[:-3]] + [seq[:-3] for seq in founds_antisense if 'N' not in seq[:-3]]
			res += founds
	return res


def give_complementary(seq):
	res = []
	for i in range(len(seq)) :
		if seq[len(seq)-1-i] == 'A':
			res.append('T')
		elif seq[len(seq)-1-i] == 'T':
			res.append('A')
		elif seq[len(seq)-1-i] == 'C':
			res.append('G')
		elif seq[len(seq)-1-i] == 'G':
			res.append('C')
		elif seq[len(seq)-1-i] == 'N':
			res.append('N')
	return ''.join(res)

def give_complementary_old(seq):
	res = []
	for letter in seq:
		if letter == 'A':
			res.append('T')
		elif letter == 'T':
			res.append('A')
		elif letter == 'C':
			res.append('G')
		elif letter == 'G':
			res.append('C')
		elif letter == 'N':
			res.append('N')
	return ''.join(res)

def find_offtagrets(seq, chromo_folder):
	'''
	:param seq:
	:param chromo_folder: a folder in which there are exactly all the chromosomse
	:return:
	'''

def get_targets_sites_from_exons_lst(exons_lst, original_range_in_gene = [0,1], min_length= 20, max_length = 20,start_with_G = False, PAMs = ["GG"]):
	if original_range_in_gene[1] <= original_range_in_gene[0]:
		print("The range of the targts on the gene is not in the right format")
		exit(-1)
	if max_length < min_length:
		print("The range of the lengths of the sgRNA is not in the right format")
		exit(-1)
	res = []
	# make a list with exons length
	lengths = [len(exon) for exon in exons_lst]
	gene_length = sum(lengths)
    #range in gene - used for deciding what parts to consider in the gene (multiply the 'where_in_gene' argument by the length of the gene)
	range_in_gene = [int(r*gene_length) for r in original_range_in_gene]

	where_in_gene = 1
	#accumolate the length of exons
	for i in range(1, len(lengths)):
		lengths[i] = lengths[i-1] + lengths[i]
	# loop on each exon
	for i in range(len(exons_lst)):
		# if there is only one exon check target
		if i == 0:
			if range_in_gene[0] < lengths[i]:
				res += get_sites(exons_lst[i][range_in_gene[0] : min(lengths[i], range_in_gene[1])], min_length, max_length, start_with_G, where_in_gene, PAMs)
		# else go over each exon and get sites
		elif max(range_in_gene[0], lengths[i-1]) < min(lengths[i], range_in_gene[1]):
			res += get_sites(exons_lst[i][max(range_in_gene[0]  - lengths[i-1], 0) : min(lengths[i] - lengths[i-1], range_in_gene[1] - lengths[i-1])], min_length, max_length, start_with_G, where_in_gene, PAMs)
	# return a list of sites for the gene
	return res

def test_2():
	gene = ["TTTATGTCAACTTTTTCAATCTAATAGATCAATGAATTGTAAACTTTTTTCGACCACAAAATGATGCTTCCAAATACAAACAAAACCTGATGCAATCAGTCAATACCTTCCAACTTTAGAACACATATATGTAGCAATGCTCCTACAGTTTACTTTTCTATCTTTTAGCCTAATCATTTACTCTCATATTTTTTCTTTAAACTAGAAAGTTCAGAATCCAAATATAATATCATCTCCTTCTCTCTATTACAGCAATGGTTTTGGTTGATAACCATGCTGGAAAAGATGGTGCAGAAGATGGTAATATGGTTGATTTTCGAGGAAATCCGGTGGATAAGTCTAGGACAGGGGGATGGCTAGCTGCAGGACTTATCCTAGGAACTGAGCTATCAGAAAGGGTATGTGTTATGGGGATTTCGATGAATTTAGTGACGTACTTAGTTGGAGATTTACATCTTCCATCCTCCAAATCTGCCAACATTGTCACCAATTTCATGGGGACACTTAATCTTCTTGGTCTTCTAGGTGGTTTCTTGGCAGATGCTAAACTCGGACGTTATCTGACTGTTGGAATCTTTGCTTCAATTGCTGCTGTGGGGGTTACGCTTTTGACATTGGCGACATCCATTCCAGGCATGAAGCCGCCTGAATGTAACCCAAGAAAAAGTGGTCACTGCATTGAAGCCAGTGGCCAGCAGCTTGCTCTTCTCTATACGGCGCTTTACATCCTAGCTCTTGGTGGTGGTGGAATTAAGTCAAATGTCTCCGGGTTTGGTTCAGACCAATTTGACTCATCAGATCCTAAGGAGAACAAGTCCATGATATACTTCTTCAACAGATTCTATTTCTGCATAAGCCTTGGTTCTCTGTTTGCAGTGACTGTGCTGGTGTACTTACAAGACAATGTAGGAAGAGGATGGGGATATGGGATATCAGCAGGCACAATGGTCCTCGGGGTCGCTGTATTGATTGGTGGAACGACGTTGTATCGATTCAAGAAGCCTCAAGGAAGTCCTTTGACTATCATATGGAGGGTTCTGCTTTTAGCTTGGAGGAAGAGAAAGCTTAGTTACCCTTCTGATACTGGCTTCTTGAATGAATATCACAATGCCAAAGTCCCACATACACATATGTTGAGGTGTCTTGACAAGGCAGCCATTCTTGATGACTCTGCAGCTGCAAATGAGAATAGCAAGAATCGTTGGATAGTTTCAACAGTTACAGAAGTCGAAGAAGTGAAAATGGTGCTCAAATTGATTCCCATATGGTCCACATGCATACTTTTTTGGACAGTATACTCTCAGATGAATACCTTCACCATTGAACAAGCTACCTTCATGAACCGGAATGTTGGAAACTTTGCTGTCCCTGCAGGTTCCTTATCCGTGTTTCTCTTTATTAGCATACTTCTGTTTACTTCCATAAACGAAAGGGTCACAGTTCGTATTGCCAGAAAAATCACTCACAACAGCCAAGGAATCACAAGCCTTCAGAGAGTTGGAATTGGACTACTACTCTCTATTGTTGGTATGGTAGCTTCAGCTCTGGTAGAAAAACGACGAAGGGAACATGCCATCCATCATAACTTCAAGATAAGCGCGTTTTGGTTAGTGCCTCAATTCTTCATTGTAGGTGCTGGGGAAGCTTTTGCCTATGTAGGACAGCTAGAGTTTTTCATCAGGGAGGCACCAGAAGGGATGAAATCTATGAGCACAGGCCTATTTCTCAGCACACTCTCGATGGGATATTTCGTGAGTAGTTTGCTAGTATTCGTTGTACAGAAAGCAACAAAAGGAAGATGGCTTAAAAGCAATTTAAACAAAGGAAAACTGGATTTATTCTACTGGTTGCTAGCAGTTCTCGGAGTAATTAATTTCTTGATTTTCATTGCATTTTCAATGAAACACCAATACAAGGTGCAGAAACTTAGCAGTATTGAGGATTCTGCAGAAGAGCTCGGGAGTTGGAAGGATTTGACCCTCGACAACAAGGAAAAGAAACTCGAAGCAGACGAGAAGGTGGAAGCTTAAATACAGCATATTAGCTTTCAATGAATCATTCATTTCCAGAGTTTGTAATATAGAACCGTATTCAATTATCAAAGACGTCAATACAAATTTGCTACCAGTCTTGAGTTCTGTTTAGATTAAAACCTTGGATATTAGAGTGCAGAAATATGATCAATTCAGAAAGATATTTACACTTCAAATTCTCACTAAA"]
	print(get_targets_sites_from_exons_lst(gene))
	print(get_sites("".join(gene)))

if __name__ == "__main__":
	test_2()