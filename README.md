#CRISPys - A tool for the optimal design of sgRNAs for simultaneous cleavage of multiple genes.

CRISPys is a powerful tool designed to facilitate the optimal design of sgRNAs for the simultaneous cleavage of multiple genes. It clusters potential CRISPR-Cas9 targets within a set of genomic sequences and generates highly promising sgRNA designs for each cluster.

Please note that this is the standalone version of CRISPys. An online version of the tool is available at http://multicrispr.tau.ac.il/.

### Installation
To run CRISPys on a GNU/Linux-based operating system, clone the repository.

We recommend creating a conda environment and installing all the necessary dependencies using the crispys.yml file.
If you intend to use uCRISPR as the scoring function, extract the contents of the uCRISPR.zip file and update the value of the global DATAPATH variable to match the extracted directory. 
### Input File
CRISPys requires a Fasta file as input, containing genomic sequences from various resources (genes), with each sequence accompanied by a unique header. It is also possible to include different coding sequences of the same gene by assigning them the same header.

Here's an example of an input file:

    >Gene1
    ATGCGCACGCCCTACCTCCATGATCCACTGACGTCCCTGAGGCTGCAATACATG
    >Gene1
    CAACCGGGCAGTCTCCGCGGCAAGTCCTAGTGCAATGGGGCTTTTTTTTAA
    >Gene2
    ATGTGCACAGGCTAAATCCATGATCCACTGACGTCCCTGAGGCACAGTGACCAATACATG
    >Gene3
    ATGTGCACAGGCTAAATACATGATAACACTATCCTATCCGTGGGGCACAGTGACCAATACCAC

### Running CRISPys
To run CRISPys, use the Stage0.py script. Run Stage0.py --help to view all available arguments.
For most users, configuring the following parameters will suffice:

* input_file: Path to the input file.
* output_dir: Path to the output directory.
* algorithm: Choose between "default" or "gene_homology" options. The "default" option finds sgRNAs for all genes in the input file, while the "gene_homology" option creates a sequence homology-based tree and outputs sgRNAs for each subgroup.
* where_in_gene: Determines the percentage of gene length used to search for CRISPR targets. For example, setting this variable to 0.5 will search for targets only in the first half of the gene.

Additional configurable parameters include:

* scoring_function: Determines the function used to score the efficiency of the sgRNA. Refer to the help page for available functions.
* omega: Sets the threshold for the sgRNA score. Guides with scores below the threshold will not be considered. The threshold depends on the scoring function used, and users should be familiar with the score distribution to select an appropriate value. (Note: In the past, we used 0.43 for CFD and 0.15 for MOFF, but this is not a recommendation.)

#### An example of CRISPys run:
If uCRISPR was selected as the scoring function, run the following command before executing:

    export DATAPATH=/<path_to_repo>/CRISPys/uCRISPR/RNAstructure/data_tables/

Execute the following command to run CRISPys with the provided parameters:

    python Stage0.py /home/fam_seq/fam1.fa /home/CRISPys_output/fam1 --output_name fam1 --alg gene_homology --where_in_gene 0.8 --omega 0.15 --off_scoring_function moff --internal_node_candidates 200 --max_target_polymorphic_sites 12

For any inquiries, please feel free to contact the author at galhyams@gmail.com.

For commercial use, kindly reach out to the author for further arrangements.




