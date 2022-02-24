__author__ = 'GH'
##Unix vession

import os
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO # app to read fasta files (added by Udi)
import time

# get the absolute path of the script
PATH = os.path.dirname(os.path.realpath(__file__))

####################################################################################################################################################################
##the pipeline: make FASTA file, run mafft on it to get the output in a fasta file, convert this FASTA to PHYLIP, and run PHYLIP's tool for making distance matrix##
####################################################################################################################################################################

def make_fasta_file(names_lst, seq_lst, outfile):
	''' should outfile have .fa ending ?'''
	f = open(outfile, 'w')
	for i in range(len(names_lst)):
		f.write('>'+names_lst[i]+'\n')
		f.write(seq_lst[i] + '\n')
	f.close()

def call_mafft_0(in_file, out_file):
	#mafft_exe = "D:\Gal\MultiCrisper\mafft-7.245-win64\mafft-win\mafft.bat"
	#in_file = "../Doc/examples/opuntia.fasta"
	#mafft_cline = MafftCommandline(mafft_exe, input=in_file)
	mafft_cline = MafftCommandline(input=in_file)	
	print(mafft_cline)
	stdout, stderr = mafft_cline()
	with open(out_file, "w") as handle:
		handle.write(stdout)
	##from Bio import AlignIO  ##not in use for now
	##	align = AlignIO.read("aligned.fasta", "fasta")  ##not in use for now
	return out_file

def call_mafft_distout(in_file, out_file):
	#mafft_exe = "D:\Gal\MultiCrisper\mafft-7.245-win64\mafft-win\mafft.bat"
	#in_file = "../Doc/examples/opuntia.fasta"
	cycles = 0
	while not os.path.isfile(in_file):  # wait for the in_file to be created added by Udi 03022022
		if cycles < 11:
			time.sleep(10)
			cycles += 1
			if cycles == 10:
				print(f"wheited for 100 seconds for {in_file} and couldnt find it, check out function 'call_mafft_distout' in mafft_and_phylip.py")
				break
	command = "mafft --distout" + in_file + " > "+  out_file
	os.system(command)


def FASTA_to_PHYLIP_old(in_file, out_file):
	os.system("perl Fasta2Phylip.pl "+ in_file + " " + out_file)
	return out_file

# def FASTA_to_PHYLIP(in_f, out_f):
# 	os.system('perl ' + PATH +  '/convertMsaFormat.pl '+in_f + ' ' +out_f+' fasta phylip')

def FASTA_to_PHYLIP(in_f, out_f):
	'''
	A new function to change format of the alignments from fasta to phylip using biopython
	This function replace the old one that used a perl script (above)
	written by Udi 25/01/22
	'''
	aligned_genes = list(SeqIO.parse(in_f, "fasta"))
	SeqIO.write(aligned_genes, out_f, "phylip")

def call_protdist_using_q(phylip_file, protdist_file, outpath):
	#print("protdist file = ",protdist_file)
	#make sh file:
	ssh_path = outpath +'/qsub.sh'
	ssh_file = open(ssh_path ,'w')
	ssh_file.write('#!/bin/tcsh\n#$ -N MultiCRISPR_1452010925\n#$ -S /bin/tcsh\n#$ -cwd\n#$ -l bioseq\n#$ -e /bioseq/data/results/multicrispr/1452010925/$JOB_NAME.$JOB_ID.ER\n#$ -o /bioseq/data/results/multicrispr/1452010925/$JOB_NAME.$JOB_ID.OU\ncd /bioseq/data/results/multicrispr/1452010925\necho "/bioseq/data/results/multicrispr/1452010925/phylip_file.ph\nF\n/bioseq/data/results/multicrispr/1452010925/protdist_file.ph\nY" | protdist')
	#submit to q:
	os.system('ssh bioseq@lecs2 qsub '+ ssh_path) #/bioseq/data/results/multicrispr/1452010925/qsub.sh')
	#os.system('echo "' + phylip_file+ '\nF\n'+ protdist_file+'\nY" | protdist')  #for my file
	while not (os.path.isfile(protdist_file)):
		time.sleep(1)


	print("is the file exist?", os.path.isfile(protdist_file))

def call_protdist(phylip_file, protdist_file):
	print("protdist file = ",protdist_file)
	os.system('tcsh -c echo "' + phylip_file+ '\nF\n'+ protdist_file+'\nY" | protdist')  #for my file
	print("is the file exist?", os.path.isfile(protdist_file))
	#print('echo ' + phylip_file+ '\\nF\\nCALM_HUMAN_ProtDist\\nY" | protdist')
	#os.system('echo '+ "/groups/pupko/haim/CALM_HUMAN.GoodProteins.Q05_S05_ID_35_ProtPerOrg1.Top25.phylip\nF\nCALM_HUMAN_ProtDist\nY | protdist")  #haim's file

def call_protdist1(phylip_file):
	os.system("protdist "+phylip_file)

def call_phylip_UPGMA(protdist_file):
	'''to change the last command line argiment to be input'''
	os.system('echo "CALM_HUMAN_ProtDist\nF\nCALM_HUMAN_ProtDist.UPGMA\nN\nY" | neighbor ; mv outtree CALM_HUMAN_ProtDist.UPGMA.tree')

def delete_temps(*paths):
	for path in paths:
		os.remove(path)


def runIt(names_lst, seq_lst, protdist_outfile, out_path):
	#fasta_outfile, aligned_file, phylip_file = out_path+"/fasta_outfile.fa", out_path+"/aligned_mattf.fa",  out_path+"/phylip_file.ph"
	fasta_outfile, aligned_file, phylip_file = out_path+"/fasta_outfile.fa", out_path+"/aligned_mattf.fa",  out_path+"/infile"

	make_fasta_file(names_lst, seq_lst, fasta_outfile)
	call_mafft_0(fasta_outfile, aligned_file)
	#call_mafft_distout(fasta_outfile,fasta_dist)
	FASTA_to_PHYLIP(aligned_file, phylip_file)
	#call_protdist_using_q(phylip_file, protdist_outfile, out_path)
	call_protdist1(phylip_file, protdist_outfile, out_path)
	#call_protdist(phylip_file, protdist_outfile)
	#call_phylip_UPGMA(protdist_file)
	#delete_temps(fasta_outfile, aligned_file, phylip_file)

def call_protdist1(phylip_file, protdist_outfile, outpath):
	#os.system("protdist "+phylip_file)
	os.chdir(outpath)
	os.system('echo "Y\r\n" | protdist')  #for my file


def test1():
	names_lst = ['s1','s2','s3']
	seq_lst = ['ACGT', 'ACGG', 'ACGA']
	runIt(names_lst,seq_lst)
def test2():
    genesNames = ["acggacgtgtacgtacgtgtt", "acggagctctacgtagctctt", "acggattgtgacgtattgtgt", "atgcacgtgtatgcacgtgtt", "ataggcatgcatgcatgctgg"]
    genesList = ["acggacgtgtacgtacgtgtt", "acggagctctacgtagctctt", "acggattgtgacgtattgtgt", "atgcacgtgtatgcacgtgtt", "ataggcatgcatgcatgctgg"]
    runIt(genesNames,genesList)

def test_mafft_distout():
	names_lst = ['s1','s2','s3']
	seq_lst = ['ACGT', 'ACGG', 'ACGA']
	fasta_outfile, aligned_file, phylip_file = "fasta_outfile.fa", "aligned_mattf.fa",  "phylip_file.ph"
	make_fasta_file(names_lst, seq_lst, fasta_outfile)
	call_mafft_distout(fasta_outfile, aligned_file)



if __name__ == '__main__':
	test1()
	#test_mafft_distout()
	#out_path = os.path.dirname(os.path.realpath(__file__)) #deafult out_path for windows
