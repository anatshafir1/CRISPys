import pandas as pd
from dotmap import DotMap
from pathlib import Path,PosixPath

opt= DotMap()
opt.path = Path('../data/').expanduser()
opt.df_ref_gRNA = opt.path.joinpath('suplementry1.csv')

df_ref_gRNA_name = pd.read_csv(opt.df_ref_gRNA)

def get_gRNA_seq(row):
    gRNA = row['Seq'][:20]
    return gRNA
# df_ref_gRNA_name['gRNA_seq'] = df_ref_gRNA_name.apply( lambda x:get_gRNA_seq(x), axis=1 )

def get_pos(row):
    pos = row['Seq'].rfind(row['gRNA_seq'])
    return pos
# df_ref_gRNA_name['target_pos'] = df_ref_gRNA_name.apply( lambda x:get_pos(x), axis=1 )

def get_xfix(row):
    suffix_start = row['target_pos'] + 23
    suffix_end = suffix_start + 10
    prefix_start = row['target_pos'] - 10
    prefix_end = row['target_pos']
    suffix = row['Seq'][suffix_start:suffix_end]
    prefix = row['Seq'][prefix_start:prefix_end]
    # print(suffix)
    # print(prefix)
    return prefix, suffix


# df_ref_gRNA_name[['prefix','suffix']] = df_ref_gRNA_name.apply( lambda x: pd.Series( get_xfix(x) ), axis=1 )


def change_fq2df(file_path):
    df = pd.read_csv(file_path,header=None)
    df.reset_index(inplace=True,drop=True)
    #提取每4行中的第2行
    df_read = df[df.index % 4 == 1]
    df_read.columns=['read']
    df_read.reset_index(inplace=True,drop=True)
    #对每个read的个数计数并且生成 dataframe
    df_read.columns=['read']
    df_unique = df_read.read.value_counts().reset_index()
    df_unique.columns=['read','counts']
    return df_unique

sra_path = 'C:/Softwares/sratoolkit.2.10.8-win64/bin/'
file = sra_path + 'SRR8935792.fastq'
df_reads = change_fq2df(file)
stem = file.name
df_reads.to_pickle(f'{stem}.pkl')


a=0