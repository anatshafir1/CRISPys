# CRISPys
•	Background: Due to extensive history of local and large-scale genomic duplications, many eukaryotic genomes harbor homologous gene families of partially overlapping functions. This redundancy often leads to mutational robustness such that the inactivation of one gene often results in no or minimal phenotypic consequence. The CRISPR system has been recently adopted as a genome editing technique of eukaryotic genomes. The system is directed to the genomic site using a programmed single-guide RNA (sgRNA) that base-pairs with the DNA target, subsequently leading to a site-specific double-strand break, frequently leading to inactivation of the encoded protein. Yet, the low specificity of the system could, in fact, be harnessed to enable a rational design of an sgRNA that would target multiple genes simultaneously.  
•	Main advances: Previous works for accomplishing this task is to align the sequences of the given genes and then locate highly similar CRISPR-Cas9 target sites among them using the consensus sequence, while allowing for few mismatches between the consensus and each of the aligned sequences. This approach may miss the most promising sgRNA candidates. We here propose a graph-based method for finding of the most effective sgRNA candidates under alternative sets of assumptions and based the ground for farther research on choosing a group of sgRNA aim to silence large gene families. In addition, we suggested a novel approaches for transforming any topological space to a metric space, which was needed in our research.

This is the standalone version of CRISPys, a tool for designing the most promising sgRNA candidates for cleaving several genes simultaneously. A webwerver can be found at http://multicrispr.tau.ac.il/. 

Running instractions:

dependencies: Python 3.5 oh higher with the Biopython module; For designing the sgRNAs while considering genes homology, please install MAFFT (https://mafft.cbrc.jp/alignment/software/) and keep the Protdist (http://evolution.genetics.washington.edu/phylip/doc/protdist.html) execution file (supplied here for convenience) in the code directory as well.

CMD running command:
Python stage0.py path_of_input_file running_directory_path

Additional arguments that may be hendeled (under the following order; see the webserber and paper for detailed description of all of this arguments):


alg (algorithm): ither A or E. Chooes A for the defluts CRISPys algorithm, or E for considering genes homology when desining the sgRNAs. The considering homology option can run only on Unix. Default: A.

where_in_gene: floating point number between 0 to 1. Search for targets only on the first 'where_in_gene' fraction  of the gene, in order to increase the liklihood that cleavage in the chosen target will be upstream to the active site and results in a loss of function mutation. Default: 1.

use_thr. Boolean argument. Default: 0. dsing sgRNA to gain maximal gaining score among all of the input genes (0) or for the maximal cleavage liklihood only among genes with score higher then the avrage.

Omega: the value of the thrashold. Recomended value: 0.43

df_targets: the scoring fundtion of the targets. Default:  Metric.cfd_funct

protodist_outfile. Default: "outfile".

min_length, of the target site. Default:20

max_length, of the target site, Default:20

start_with_G: rather the target sites are obligated to start with a G codon. Default:False

internal_node_candidates. When choosing the consider homology option, this is the number of sgRNAs desinged to each homology sub-group.
Default: 10

The maximal number of polymorphiic site of a targets cluster which the sgRNA are accordingly desinded. Default: 12
