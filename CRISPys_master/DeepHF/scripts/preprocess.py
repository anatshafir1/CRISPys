import numpy as np
import pickle
import pandas as pd
import os
import time
from multiprocessing import Pool
from . import feature_util

def char2int(char):
    if char == 'A':
        return 1
    if char == 'T':
        return 2
    if char == 'C':
        return 3
    if char == 'G':
        return 4
    else:
        print('Recived wrong char {} - exiting'.format(char))
        exit(1)



def split_to_train_test_file(config):
    # Creating a test file that is made of 15% of all gRNAs that have more than 100 reads in all 3 enzymes.
    # This test file will be used in all the tests for comparison.
    # The rest of the data will be used as the training data.
    if config.data_source == 'new':
        dir_path = 'data/preprocessed_data/new_data/'
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        if os.path.exists(dir_path + 'test_seq.csv'):
            return
        eff_df = pd.read_csv('data/preprocessed_data/final_efficiency_with_bio.csv')
        eff_df.drop('index', axis='columns', inplace=True)
        eff_df.rename(columns={'stem': 'epi1', 'dG': 'epi2', 'dG_binding_20': 'epi3', 'dg_binding_7to20': 'epi4',
                               'GC > 10': 'epi5', 'GC < 10': 'epi6', 'GC count': 'epi7', 'Tm global_False': 'epi8',
                                 '5mer_end_False': 'epi9', '8mer_middle_False': 'epi10', '4mer_start_False': 'epi11'}
                      , inplace=True)

        # For the test file we use only reads with threshold >= 100
        wt_cond = eff_df.wt_reads_sum >= 100
        esp_cond = eff_df.esp_reads_sum >= 100
        hf_cond = eff_df.hf_reads_sum >= 100
        test_df = eff_df[wt_cond & esp_cond & hf_cond].sample(frac=0.15).sort_index()
        train_df = eff_df.drop(test_df.index)

        test_df.to_csv(dir_path + 'test.csv', index=False)
        train_df.to_csv(dir_path + 'train_valid.csv', index=False)

    elif config.data_source == 'old':
        # The original data is given in suplementry2 without the biofeatures
        dir_path = 'data/preprocessed_data/original_data/'
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        if os.path.exists(dir_path + 'test.csv'):
            return

        if os.path.exists(dir_path + 'full_data.csv'):
            eff_df_with_bio = pd.read_csv(dir_path + 'full_data.csv')
        else:
            eff_df = pd.read_csv('data/suplementry2.csv')

            feature_options = {
                "testing_non_binary_target_name": 'ranks',
                'include_pi_nuc_feat': True,
                "gc_features": True,
                "nuc_features": False,
                "include_Tm": True,
                "include_structure_features": True,
                "order": 3,
                "num_proc": 20,
                "normalize_features": None
            }

            # Create the biofeatures and add them to the dataframe
            feature_sets = feature_util.featurize_data(eff_df, feature_options)
            gc_above_10 = feature_sets['gc_above_10']
            gc_below_10 = feature_sets['gc_below_10']
            gc_count = feature_sets['gc_count']
            Tm = feature_sets['Tm']
            dG_features = feature_sets['dG_features']
            dG_features.reset_index(inplace=True)
            dG_features.drop(columns=['index'], inplace=True)
            eff_df_with_bio = pd.concat([eff_df, dG_features, gc_above_10, gc_below_10, gc_count, Tm], axis=1)

            eff_df_with_bio.rename(columns={'stem': 'epi1', 'dG': 'epi2', 'dG_binding_20': 'epi3', 'dg_binding_7to20': 'epi4',
                                   'GC > 10': 'epi5', 'GC < 10': 'epi6', 'GC count': 'epi7', 'Tm global_False': 'epi8',
                                     '5mer_end_False': 'epi9', '8mer_middle_False': 'epi10', '4mer_start_False': 'epi11',
                                   'Wt_Efficiency': 'wt_mean_eff', 'eSpCas 9_Efficiency': 'esp_mean_eff',
                                   'SpCas9-HF1_Efficiency': 'hf_mean_eff'}
                          , inplace=True)

            eff_df_with_bio.to_csv(dir_path + 'full_data.csv', index=False)


        test_df = eff_df_with_bio.sample(frac=0.15).sort_index()
        train_df = eff_df_with_bio.drop(test_df.index)

        test_df.to_csv(dir_path + 'test.csv', index=False)
        train_df.to_csv(dir_path + 'train_valid.csv', index=False)

def prepare_test_seq(config):
    # Read the test file and organize it in a pickled Seq object
    if config.data_source == 'old':
        dir_path = 'data/preprocessed_data/original_data/'
    elif config.data_source == 'new':
        dir_path = 'data/preprocessed_data/new_data/'
    else:
        print('UnKnown data source')
        exit(1)

    if os.path.exists(dir_path + 'test_seq.pkl'):
        return
    test_df = pd.read_csv(dir_path + 'test.csv')
    test_seq = Seq()
    for index, row in test_df.iterrows():
        print('line: {}'.format(index))
        seq = row['21mer']
        wt_eff = row['wt_mean_eff']
        esp_eff = row['esp_mean_eff']
        hf_eff = row['hf_mean_eff']

        biofeat = []
        for i in range(1, 12):
            biofeat.append(row['epi{}'.format(i)])

        test_seq.add_seq(seq, wt=wt_eff, esp=esp_eff, hf=hf_eff, biofeat=biofeat)

    with open(dir_path + 'test_seq.pkl', "wb") as fp:
        pickle.dump(test_seq, fp)

def prepare_train_valid_seq(config):
    # Read the train file and organize it in a pickled Seq object
    if config.data_source == 'old':
        dir_path = 'data/preprocessed_data/original_data/'
    elif config.data_source == 'new':
        dir_path = 'data/preprocessed_data/new_data/'
    else:
        print('UnKnown data source')
        exit(1)

    if not os.path.exists(dir_path + 'train.csv'):
        # Split to train and valid
        train_val_df = pd.read_csv(dir_path + 'train_valid.csv')
        valid_df = train_val_df.sample(frac=0.1).sort_index()
        train_df = train_val_df.drop(valid_df.index)

        valid_df.to_csv(dir_path + 'valid.csv', index=False)
        train_df.to_csv(dir_path + 'train.csv', index=False)
    else:
        valid_df = pd.read_csv(dir_path + 'valid.csv')
        train_df = pd.read_csv(dir_path + 'train.csv')



    if not os.path.exists(dir_path + 'valid_seq.pkl'):

    # Prepare valid data in Seq object
        valid_seq = Seq()
        for index, row in valid_df.iterrows():
            print('line: {}'.format(index))
            seq = row['21mer']
            wt_eff, esp_eff, hf_eff = row['wt_mean_eff'], row['esp_mean_eff'], row['hf_mean_eff']
            if config.data_source == 'new':
                wt_reads_sum, esp_reads_sum, hf_reads_sum = row['wt_reads_sum'], row['esp_reads_sum'], row['hf_reads_sum']
                confidence = [wt_reads_sum, esp_reads_sum, hf_reads_sum]
                for ind, con in enumerate(confidence):
                    if np.isnan(con):
                        confidence[ind] = 0
            else:
                confidence = None

            biofeat = []
            for i in range(1, 12):
                biofeat.append(row['epi{}'.format(i)])

            valid_seq.add_seq(seq, wt=wt_eff, esp=esp_eff, hf=hf_eff, biofeat=biofeat, conf=confidence)

        with open(dir_path + 'valid_seq.pkl', "wb") as fp:
            pickle.dump(valid_seq, fp)


    # Prepare train data in Seq object
    if (config.data_source == 'old') or (config.data_source == 'new' and config.weighted_loss):
        train_seq = Seq()
        for index, row in train_df.iterrows():
            print('line: {}'.format(index))
            seq = row['21mer']
            wt_eff, esp_eff, hf_eff = row['wt_mean_eff'], row['esp_mean_eff'], row['hf_mean_eff']
            if config.data_source == 'new':
                wt_reads_sum, esp_reads_sum, hf_reads_sum = row['wt_reads_sum'], row['esp_reads_sum'], row['hf_reads_sum']
                confidence = [wt_reads_sum, esp_reads_sum, hf_reads_sum]
                for ind, con in enumerate(confidence):
                    if np.isnan(con):
                        confidence[ind] = 0
            else:
                confidence = None

            biofeat = []
            for i in range(1, 12):
                biofeat.append(row['epi{}'.format(i)])

            train_seq.add_seq(seq, wt=wt_eff, esp=esp_eff, hf=hf_eff, biofeat=biofeat, conf=confidence)

        if config.data_source == 'new':
            dir_path += 'weighted_loss/'
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)

        with open(dir_path + 'train_seq.pkl', "wb") as fp:
            pickle.dump(train_seq, fp)

    elif (config.data_source == 'new' and not config.weighted_loss):
        dir_path += 'row_reads/'
        if not os.path.exists(dir_path):
            os.mkdir(dir_path)
        crt_row_reads_model_data(config, dir_path)

def prepare_leenay():
    leenay = open('data/Leenay_mean_eff_coordinates.tsv')

    leenay_seq = Seq(has_biofeat=False)
    _ = leenay.readline()
    for line in leenay:
        data = line.split('\t')
        mer = data[7] + data[8][0]
        eff = data[4]
        leenay_seq.add_seq(mer, eff, 0, 0)

    with open("data/preprocessed_data/leenay_seq.pkl", "wb") as fp:
        pickle.dump(leenay_seq, fp)


## Row reads model preperations
def crt_row_reads_model_data(config, dir_path):
    full_train_df = pd.read_csv(dir_path + 'train.csv')

    # Extract only the enzyme columns and biofeatures
    biofeatures_cols = [col for col in full_train_df.columns if 'epi' in col]

    relevant_cols = ['21mer', f'{config.enzyme}_reads_sum', f'{config.enzyme}_edited_read_counts'] + biofeatures_cols
    enzyme_train_df = full_train_df[relevant_cols]
    enzyme_train_df.rename(columns={f'{config.enzyme}_reads_sum': 'reads_sum',
                                    f'{config.enzyme}_edited_read_counts': 'edited_read_counts'}
                           , inplace=True)

    enzyme_train_df.dropna(subset=['reads_sum'], inplace=True)

    # Create train_file
    temp_dir_path = dir_path + 'temp_train_files/'
    if not os.path.exists(temp_dir_path):
        os.mkdir(temp_dir_path)

    length = enzyme_train_df.shape[0]
    parts = length//1000
    for i in range(parts):

        if i == parts-1:
            print(f'from {i*1000} to end')
            train_df = enzyme_train_df[i * 1000:]
        else:
            print(f'from {i*1000} to {(i+1)*1000}')
            train_df = enzyme_train_df[i * 1000: (i + 1) * 1000]
        df_split = np.array_split(train_df, 100)
        pool = Pool(8)
        start = time.time()
        a = pool.map(prepare_seq, df_split)

        train_seq = RowReads(config.enzyme)

        for seq in a:
            train_seq.X = np.concatenate((train_seq.X, seq.X), axis=0)
            train_seq.X_biofeat = np.concatenate((train_seq.X_biofeat, seq.X_biofeat), axis=0)
            train_seq.y = np.concatenate((train_seq.y, seq.y), axis=0)
        pool.close()
        pool.join()
        end = time.time()
        print(f'time: {end - start}')

        perm = np.random.permutation(len(train_seq.X))
        train_seq.X = train_seq.X[perm]
        train_seq.X_biofeat = train_seq.X_biofeat[perm]
        train_seq.y = train_seq.y[perm]
        with open(temp_dir_path + f'row_read_{config.enzyme}_seq_{i}.pkl', "wb") as fp:
            pickle.dump(train_seq, fp)

    # Unite files
    full_train_seq = RowReads(config.enzyme)
    for i in range(parts):
        print(f'Uniting part {i}')
        with open(temp_dir_path + f'row_read_{config.enzyme}_seq_{i}.pkl', "rb") as fp:
            train_seq = pickle.load(fp)
        full_train_seq.X = np.concatenate((full_train_seq.X, train_seq.X), axis=0)
        full_train_seq.X_biofeat = np.concatenate((full_train_seq.X_biofeat, train_seq.X_biofeat), axis=0)
        full_train_seq.y = np.concatenate((full_train_seq.y, train_seq.y), axis=0)

    perm = np.random.permutation(len(full_train_seq.X))
    full_train_seq.X = full_train_seq.X[perm]
    full_train_seq.X_biofeat = full_train_seq.X_biofeat[perm]
    full_train_seq.y = full_train_seq.y[perm]

    with open(dir_path + f'train_{config.enzyme}_seq.pkl', "wb") as fp:
        pickle.dump(full_train_seq, fp, protocol=4)


def prepare_seq(df_split):
    train_seq = RowReads()

    for index, row in df_split.iterrows():
        seq = row['21mer']
        reads_sum = row['reads_sum']
        edited_count = row['edited_read_counts']
        if np.isnan(edited_count):
            continue

        biofeat = []
        for i in range(1, 12):
            biofeat.append(row['epi{}'.format(i)])

        train_seq.add_seq(seq, biofeat, int(reads_sum), int(edited_count))

    return train_seq

def check_results():
    full_train_df = pd.read_csv('data/preprocessed_data/train.csv')
    esp_reads = full_train_df['esp_reads_sum']
    esp_edited_read = full_train_df['esp_edited_read_counts']
    esp_reads.dropna(inplace=True)
    esp_edited_read.dropna(inplace=True)

    esp_reads_sum = esp_reads.sum()
    esp_edited_sum = esp_edited_read.sum()

    with open(f'data/preprocessed_data/row_reads/train_esp_seq.pkl', "rb") as fp:
        train_seq = pickle.load(fp)
    result_reads_sum = train_seq.X.shape[0]
    unique, counts = np.unique(train_seq.y, return_counts=True)


# Data structure objects

class Seq(object):
    # This class is the serialaized data for all 3 cas9 types
    def __init__(self):
        self.X = np.empty((0, 22), np.uint8)
        self.X_biofeat = np.empty((0, 11), np.float16)
        self.y = {'wt': np.empty(0, np.float16), 'esp': np.empty(0, np.float16), 'hf': np.empty(0, np.float16)}
        self.confidence = {'wt': np.empty(0, np.uint16), 'esp': np.empty(0, np.uint16), 'hf': np.empty(0, np.uint16)}


    def add_seq(self, seq, wt=None, esp=None, hf=None, biofeat=None, conf=None):
        # The '0' represent the beginning of the sequence (help for the RNN)
        mer_array = np.array([0], dtype=np.uint8)

        for char in seq:
            num = char2int(char)
            mer_array = np.append(mer_array, np.array([num], dtype=np.uint8), axis=0)

        mer_array = np.expand_dims(mer_array, 0)
        self.X = np.concatenate((self.X, mer_array), axis=0)

        biofeat = np.array([float(epi) for epi in biofeat], dtype=np.float16)
        biofeat = np.expand_dims(biofeat, 0)
        self.X_biofeat = np.concatenate((self.X_biofeat,biofeat), axis=0)
        self.y['wt'] = np.append(self.y['wt'], np.array([wt], dtype=np.float16), axis=0)
        self.y['esp'] = np.append(self.y['esp'], np.array([esp], dtype=np.float16), axis=0)
        self.y['hf'] = np.append(self.y['hf'], np.array([hf], dtype=np.float16), axis=0)
        if conf is not None:
            self.confidence['wt'] = np.append(self.confidence['wt'], np.array([conf[0]], dtype=np.uint16), axis=0)
            self.confidence['esp'] = np.append(self.confidence['esp'], np.array([conf[1]], dtype=np.uint16), axis=0)
            self.confidence['hf'] = np.append(self.confidence['hf'], np.array([conf[2]], dtype=np.uint16), axis=0)



    def get_data(self):
        if self.has_biofeat:
            return self.X, self.X_biofeat, self.y
        else:
            return self.X, self.y

class RowReads(object):
    def __init__(self, enzyme=None):
        self.X = np.empty((0, 22), np.uint8)
        self.X_biofeat = np.empty((0, 11), np.float16)
        self.enzyme = enzyme
        self.y = np.empty(0, np.bool)

    def add_seq(self, seq, biofeat, reads_sum, edited_count):
        mer_array = np.array([0], dtype=np.uint8)

        for char in seq:
            num = char2int(char)
            mer_array = np.append(mer_array, np.array([num], dtype=np.uint8), axis=0)
        mer_array = np.expand_dims(mer_array, 0)
        # print(mer_array.nbytes)
        biofeat = np.array([float(epi) for epi in biofeat], dtype=np.float16)
        biofeat = np.expand_dims(biofeat, 0)

        non_edited = reads_sum - edited_count
        assert non_edited >= 0

        for i in range(edited_count):
            self.X = np.concatenate((self.X, mer_array), axis=0)
            self.X_biofeat = np.concatenate((self.X_biofeat, biofeat), axis=0)
            self.y = np.append(self.y, np.array([1], dtype=np.bool), axis=0)

        for i in range(non_edited):
            self.X = np.concatenate((self.X, mer_array), axis=0)
            self.X_biofeat = np.concatenate((self.X_biofeat, biofeat), axis=0)
            self.y = np.append(self.y, np.array([0], dtype=np.bool), axis=0)

