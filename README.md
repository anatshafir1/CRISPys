A tool for the optimal design of sgRNAs for simultaneous cleavage of multiple genes.

Given a set of genomic sequences as input, CRISPys first clusters all the potential CRISPR-Cas9 targets that are located within them, and then designs the most promising sgRNAs to each cluster.

This is the standalone version of CRISPys. An online tool can be found at http://multicrispr.tau.ac.il/.

### Installation
This version will run on gnu/linux based OS, after cloning the repository we recommend using conda to create an environment with all the program needed using the crispys.yml file.

If you want CRISPys to output singletons (sgRNAs that target one site) you need to unzip the uCRISPR.zip file and add the DATAPATH variable.
Before running CRISPys execute:
export DATAPATH=/path_to_repo/CRISPys/uCRISPR/RNAstructure/data_tables/

### Input file
CRISPys take as an input a Fasta file of genomic sequences from multiple resources (genes) each one with its uniqe header, it is possible to enter different exons of the same gene by giving them the same header.

### Run CRISPys
To run CRISPys use the Stage0.py, run 'Stage0.py --help' to see all arguments.\
The common user will probably only need to configure the input file, output dir, algorithm, where_in_gene, scoring function and omega (most of the times the other parameters can be set to default).\
* The algorithm (--alg) argument can take two values: 'default' or 'gene_homolog'. The default option will find sgRNAs for all the genes in the input file, and the gene homology option will create a tree based on sequence homology and will output sgRNAs for each subgroup of that tree, this approach described in the paper as strategy III.\The 'where in gene' argument will determine the percentage of the gene length that will be used to look for CRISPR targets, for example: if the user would like to search for targets only in the first half of the gene he should set tish variable to 0.5. 

  
* The scoring function argument determine the function that is used to score the efficiency of the sgRNA, see the help page of available functions. 
  

* The omega argument is the threshold of the sgRNA score, guide with score less than the threshold will not be considered, the threshold is dependent on the scoring function that been used, the user needs to be familiar with the distribution of the score in order to select the threshold, previously we used 0.43 for CFD and 0.15 for MOFF but this is not a recommendation.

#### An example of CRISPys run:
export DATAPATH=/home/CRISPys/uCRISPR/RNAstructure/data_tables/ \
python Stage0.py /home/fam_seq/fam1.fa /home/CRISPys_output/fam1 \
--output_name fam1 \
--alg gene_homology \
--where_in_gene 0.8 \
--omega 0.15 \
--off_scoring_function moff \
--internal_node_candidates 200 \
--max_target_polymorphic_sites 12 \

For any questions, please contact the author at galhyams@gmail.com

For any commercial use, please contact the author. 
