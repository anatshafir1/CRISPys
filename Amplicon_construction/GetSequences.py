import subprocess
from typing import List

import pandas as pd


def annotations_to_df(annotations_file_path: str):
    annotations_df = pd.read_csv(annotations_file_path, sep="\s+")  # sep='\s+' uses a regular expression to match one or more whitespace characters as the separator
    filtered_by_gene = annotations_df[annotations_df['feature'] == 'gene']
    start_char = '='
    end_char = ';'
    filtered_by_gene['start_position'] = filtered_by_gene['attribute'].str.find(start_char) + 1
    filtered_by_gene['end_position'] = filtered_by_gene['attribute'].str.find(end_char)
    filtered_by_gene['gene_id'] = filtered_by_gene.apply(lambda row: row['attribute'][row['start_position']:row['end_position']], axis=1)

    return filtered_by_gene


def get_gene_sites(out_path: str, fasta_file: str, filtered_gene_df) -> List[str]:
    """Find the sequences upstream to the gene's TSS

    :param out_path: the directory path to which the algorithm will store the results
    :param fasta_file: path to input FASTA format file of the genome
    :param filtered_gene_df:
    :return:
    """
    # create BED format file
    filtered_gene_df.to_csv(out_path + '/genes.bed', sep='\t',
                            columns=['seqname', 'start', 'end', 'gene_id', 'score', 'strand'],
                            header=False, index=False)
    bed_file = out_path + "/genes.bed"
    # run bedtools
    seq = subprocess.run(['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-nameOnly', '-s'],
                         stdout=subprocess.PIPE)
    sites_list = seq.stdout.decode().split()
    return sites_list


def extract_genes(annotations_file_path, out_path: str, genome_fasta_file: str):
    filtered_annotations_df = annotations_to_df(annotations_file_path)
    genes_list = get_gene_sites(out_path, genome_fasta_file, filtered_annotations_df)
    print(genes_list)
