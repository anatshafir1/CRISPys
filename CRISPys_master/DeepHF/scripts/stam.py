import pandas as pd
import argparse
import os
from cutadapt.adapters import AdapterParser
from cutadapt.seqio import Sequence

def get_config():
    print('Receiving parser values')
    parser = argparse.ArgumentParser()

    # Paths
    parser.add_argument('--data_path', type=str, default='../data/')  # [wt,esp,hf]
    parser.add_argument('--sra_path', type=str, default='../fastq_files/original_files/')
    parser.add_argument('--fastq_path', type=str, default='../fastq_files/')


    # Data
    parser.add_argument('-sup2', dest='sup_1', action='store_false')
    parser.add_argument('--data_type', type=str, default='wt')  # [wt,esp,hf]

    config = parser.parse_args()

    return config

def fq2df(config, report_df):
    sort_and_count_list = report_df.loc['sort_and_count', :]

    for file_name, val in sort_and_count_list.iteritems():
        if val == 0:
            file_index = 0
            while os.path.exists(config.fastq_path + 'split_files/' + file_name + f'_{file_index}_masked.fastq'):
                file_path = config.fastq_path + config.fastq_path + 'split_files/' + file_name + f'_{file_index}_masked.fastq'
                file_index += 1
                df = pd.read_csv(file_path, header=None)

                report_df.loc['sort_and_count', file_name] = 1
                report_df.to_csv(config.fastq_path + 'report.csv')
                df.reset_index(inplace=True, drop=True)
                df_read = df[df.index % 2 == 1]
                df_read.columns = ['read']
                df_read.reset_index(inplace=True, drop=True)
                df_read.columns = ['read']
                df_unique = df_read.read.value_counts().reset_index()
                df_unique.columns = ['read', 'counts']
                return df_unique

            a=0

def write_fasta(x, fasta_file):
    print('>read_{}|{}'.format(x['index'], x['counts']), x['read'][1:20], file=fasta_file, sep='\n')

def df2fa(config, report_df):
    create_fasta_list = report_df.loc['create_fasta', :]

    if not os.path.exists(config.fastq_path + 'fasta_files'):
        os.mkdir(config.fastq_path + 'fasta_files')

    for file_name, val in create_fasta_list.iteritems():
        if val == 0:
            print(f'Creating fasta file from  {file_name}_unified_df.pkl')
            file_path = config.fastq_path + f'dataframes/{file_name}_unified_df.pkl'
            unique_df = pd.read_pickle(file_path)
            fasta_file = open(config.fastq_path + f'fasta_files/{file_name}.fa', 'w')

            unique_df.reset_index().apply(lambda x: write_fasta(x, fasta_file), axis=1)
            fasta_file.close()

    a=0
config = get_config()
report_df = pd.read_csv(config.fastq_path + 'report.csv')
report_df.set_index('Unnamed: 0', inplace=True)
# fq2df(config, report_df)

# adapter_parser = AdapterParser(
#     colorspace=None,
#     max_error_rate=0.2,
#     min_overlap=17,
#     read_wildcards=False,
#     adapter_wildcards=False,
#     indels=True)
# seq_95_150 = 'CTTTTTTGCCACCATGGAAAAAAAAACTCCAAAACCCTGGTGAGCAAGGGCGAGG'
# read_seq = ('read_ing',seq_95_150)
#
# read = Sequence( name=read_seq[0], sequence=read_seq[1] )
# # -a 3'
# a = []
# # -b both
# b = []
# # -g 5'
# g = ['gRNA_prefix_suffix={}...{}'.format('gccaccATGG', 'CTGgtgagc')]
#
# adapter = adapter_parser.parse_multi( a, b, g)[0]
# r = adapter.match_to(read)
#
# if r is not None:
#     prefix_range = (r.front_match.rstart, r.front_match.rstop, r.front_match.errors)
#     suffix_ragne = (
#     r.back_match.rstart + r.front_match.rstop, r.back_match.rstop + r.front_match.rstop, r.back_match.errors)
#     target_range = (r.front_match.rstop, r.back_match.rstart + r.front_match.rstop)
#     read_prefix_seq = seq_95_150[prefix_range[0]:prefix_range[1]]
#     read_suffix_seq = seq_95_150[suffix_ragne[0]:suffix_ragne[1]]
#     read_gRNA_PAM = seq_95_150[target_range[0]:target_range[1]]
valid_dir_path = config.fastq_path + 'gRNA_valid_df/'

unique_df = pd.read_pickle(valid_dir_path + f'wt_1_valid_df.pkl')
print(unique_df.head())
