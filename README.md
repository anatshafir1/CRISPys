# CRISPys - A tool for the optimal design of sgRNAs for simultaneous cleavage of multiple genes.

CRISPys is a powerful tool designed to facilitate the optimal design of sgRNAs for the simultaneous cleavage of multiple genes. It clusters potential CRISPR-Cas9 targets within a set of genomic sequences and generates highly promising sgRNA designs for each cluster.

This is the standalone version of CRISPys, an online version is available at http://multicrispr.tau.ac.il/.

### Installation
CRISPys can run on a GNU/Linux-based operating system only. 
After clonning of the repository, we recommend creating a conda environment and installing all the necessary dependencies using the crispys.yml file.
If uCRISPR is being used, please extract the contents of the uCRISPR.zip file and add the global DATAPATH variable to match the extracted directory (see running example below). 

### Input File
CRISPys requires a Fasta file as input, containing genomic sequences from various resources (genes), with each sequence accompanied by a unique header. It is also possible to include different coding sequences of the same gene by assigning them the same header.

Here's an example of an input file with 3 genes:

    >Gene1
    ATGCGCACGCCCTACCTCCATGATCCACTGACGTCCCTGAGGCTGCAATACATG
    >Gene1
    CAACCGGGCAGTCTCCGCGGCAAGTCCTAGTGCAATGGGGCTTTTTTTTAA
    >Gene2
    ATGTGCACAGGCTAAATCCATGATCCACTGACGTCCCTGAGGCACAGTGACCAATACATG
    >Gene3
    ATGTGCACAGGCTAAATACATGATAACACTATCCTATCCGTGGGGCACAGTGACCAATACCAC

### Running CRISPys
To run CRISPys, use the Stage0.py script, use ```python Stage0.py --help``` to view all parameters.
For most users, configuring the following parameters will suffice:

* input_file: Path to the input file.
* output_dir: Absolute path to the output directory.
* algorithm: Choose between "default" or "gene_homology" options. The "default" option finds sgRNAs for all genes in the input file, while the "gene_homology" option creates a sequence homology-based tree and outputs sgRNAs for each subgroup.
* where_in_gene: Determines the percentage of gene length used to search for CRISPR targets. For example, setting this variable to 0.5 will search for targets only in the first half of the gene.

Additional configurable parameters include:

* scoring_function: Determines the function used to score the efficiency of the sgRNA. Refer to the help page for available functions.
* omega: Sets the threshold for the sgRNA score. Guides with scores below the threshold will not be considered. The threshold depends on the scoring function used, and users should be familiar with the score distribution to select an appropriate value. (Note: In the past, we used 0.43 for CFD and 0.15 for MOFF, but this is not a recommendation.)

#### An example of CRISPys run:
If uCRISPR is selected as the scoring function or the "--singletons" option is set to 1 (with the "singletons_on_target_function" option set to ucrispr), please run the following command before executing:
    export DATAPATH=/<path_to_repo>/CRISPys/uCRISPR/RNAstructure/data_tables/

Execute the following command to run CRISPys with the provided parameters:

    python Stage0.py /home/fam_seq/fam1.fa /home/CRISPys_output/fam1 --output_name fam1 --alg gene_homology --where_in_gene 0.8 --omega 0.15 --off_scoring_function moff --singletons 1 --number_of_singletons 5 --slim_output 1

The command above will execute CRISPys using the input fasta file, wrinting output for each subgroup of genes (generated by using the 'gene_homology'), while each sgRNA score will be above moff score of 0.15 and CRISPR targets will be taken from up to 80% of the gene, the output will also contain 5 singletons be for each gene. 

The example output fam1.csv file:
```
Genes in group=	['GENE2'	 'GENE3'	 'GENE1']						
sgRNA index	sgRNA	    Score	Genes	Genes score	Target site	#mms	Position	strand	PAM
1	TCAGGGACGTCAGTGGATCA	2	GENE2	1	TCAGGGACGTCAGTGGATCA	0	20	        -	TGG
1			                GENE1	1	TCAGGGACGTCAGTGGATCA	0	65	        -	TGG
2	ATGATCCACTGACGTCCCTG	2	GENE2	1	ATGATCCACTGACGTCCCTG	0	19	        +	AGG
2			                GENE1	1	ATGATCCACTGACGTCCCTG	0	19	        +	AGG
3	CTGAGCCTCAGGGACGTCAG    0.611	GENE2	0.426	CTGtGCCTCAGGGACGTCAG	1	13  	        -	TGG
3			                GENE1	0.185	tgcAGCCTCAGGGACGTCAG	3	58  	        -	TGG
4	CGGAGCCTCAGGGACGTCAG	0.56	GENE1	0.347	tGcAGCCTCAGGGACGTCAG	2	58  	        -	TGG
4			                GENE2	0.212	CtGtGCCTCAGGGACGTCAG	2	13	        -	TGG
5	TGCTGCCTCAGGGACGTCAG	0.519	GENE1	0.358	TGCaGCCTCAGGGACGTCAG	1	58      	-	TGG
5			                GENE2	0.161	ctgTGCCTCAGGGACGTCAG	3	13	        -	TGG
6	CGCTGCCTCAGGGACGTCAG	0.439	GENE2	0.27	CtgTGCCTCAGGGACGTCAG	2	13	        -	TGG
6			                GENE1	0.169	tGCaGCCTCAGGGACGTCAG	2	58	        -	TGG
7	TTGAGCCTCAGGGACGTCAG	0.416	GENE1	0.236	TgcAGCCTCAGGGACGTCAG	2	58	        -	TGG
7			                GENE2	0.179	cTGtGCCTCAGGGACGTCAG	2	13	        -	TGG
Genes in group=	['GENE2']								
sgRNA index	sgRNA	    Score	Genes	Genes score	Target site	#mms	Position	strand	PAM
1	ATGATCCACTGACGTCCCTG	0.92	GENE2	1	ATGATCCACTGACGTCCCTG	0	19	        +	AGG
2	CTGTGCCTCAGGGACGTCAG	0.908	GENE2	1	CTGTGCCTCAGGGACGTCAG	0	13	        -	TGG
3	TCAGGGACGTCAGTGGATCA	0.9	GENE2	1	TCAGGGACGTCAGTGGATCA	0	20	        -	TGG
Genes in group=	['GENE3']								
sgRNA index	sgRNA	    Score	Genes	Genes score	Target site	#mms	Position	strand	PAM
1	ATAACACTATCCTATCCGTG	0.937	GENE3	1	ATAACACTATCCTATCCGTG	0	22	        +	GGG
2	GATAACACTATCCTATCCGT	0.934	GENE3	1	GATAACACTATCCTATCCGT	0	21	        +	GGG
3	TGATAACACTATCCTATCCG	0.914	GENE3	1	TGATAACACTATCCTATCCG	0	20	        +	TGG
Genes in group=	['GENE1']								
sgRNA index	sgRNA	    Score	Genes	Genes score	Target site	#mms	Position	strand	PAM
1	GGGACGTCAGTGGATCATGG	0.92	GENE1	1	GGGACGTCAGTGGATCATGG	0	68	        -	AGG
2	ATGATCCACTGACGTCCCTG	0.92	GENE1	1	ATGATCCACTGACGTCCCTG	0	19	        +	AGG
3	TGCAGCCTCAGGGACGTCAG	0.905	GENE1	1	TGCAGCCTCAGGGACGTCAG	0	58	        -	TGG
4	TCAGGGACGTCAGTGGATCA	0.9	GENE1	1	TCAGGGACGTCAGTGGATCA	0	65	        -	TGG
5	GACTTGCCGCGGAGACTGCC	0.892	GENE1	1	GACTTGCCGCGGAGACTGCC	0	25	        -	CGG
```

---

For any inquiries, please feel free to contact the author at galhyams@gmail.com or omer.cald123@gmail.com



For commercial use, kindly reach out to the author for further arrangements.




