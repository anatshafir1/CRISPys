import pickle
import numpy as np
from sklearn.model_selection import train_test_split


def get_data(config):
    if config.data_source == 'old':
        dir_path = 'data/preprocessed_data/original_data/'
        with open(dir_path + 'train_seq.pkl', "rb") as fp:
            train_seq = pickle.load(fp)
    elif config.data_source == 'new':
        dir_path = 'data/preprocessed_data/new_data/'
        if config.weighted_loss:
            with open(dir_path + f'weighted_loss/train_seq.pkl', "rb") as fp:
                train_seq = pickle.load(fp)
        else:
            with open(dir_path + f'row_reads/train_{config.enzyme}_seq.pkl', "rb") as fp:
                train_seq = pickle.load(fp)
            # train_seq.X = train_seq.X[:1000]
            # train_seq.X_biofeat = train_seq.X_biofeat[:1000]
            # train_seq.y = train_seq.y[:1000]

    with open(dir_path + 'valid_seq.pkl', "rb") as fp:
        valid_seq = pickle.load(fp)

    with open(dir_path + 'test_seq.pkl', "rb") as fp:
        test_seq = pickle.load(fp)

    DataHandler = {}

    # Loading test and validation data
    if config.enzyme == 'multi_task':
        if config.multi_data == 'parallel':
            print('Loading data for multi_task parallel model')
            parallel_enzyme_loader(config, DataHandler, test_seq, valid_seq, train_seq)

        elif config.multi_data == 'serial':
            print('Loading data for multi_task serial model')
            serial_enzyme_loader(config, DataHandler, test_seq, valid_seq, train_seq)

    else:
        single_enzyme_loader(config, DataHandler, test_seq, valid_seq, train_seq)

    if not config.has_embedding:
        # X_train_val = np.eye(5)[X_train_val]
        DataHandler['X_test'] = np.eye(5)[DataHandler['X_test']]
        DataHandler['X_train'] = np.eye(5)[DataHandler['X_train']]
        DataHandler['X_val'] = np.eye(5)[DataHandler['X_val']]

    if config.eda:
        eda(DataHandler)

    return DataHandler


def single_enzyme_loader(config, DataHandler, test_seq, valid_seq, train_seq):
    DataHandler['y_test'] = test_seq.y[config.enzyme]
    DataHandler['X_test'] = test_seq.X
    DataHandler['X_test_biofeat'] = test_seq.X_biofeat

    valid_indexes = np.logical_not(np.isnan(valid_seq.y[config.enzyme]))
    DataHandler['y_valid'] = valid_seq.y[config.enzyme][valid_indexes]
    DataHandler['X_valid'] = valid_seq.X[valid_indexes]
    DataHandler['X_valid_biofeat'] = valid_seq.X_biofeat[valid_indexes]
    if config.weighted_loss:
        DataHandler['conf_valid'] = valid_seq.confidence[config.enzyme][valid_indexes]

    if config.data_source == 'old' or config.weighted_loss:
        train_indexes = np.logical_not(np.isnan(train_seq.y[config.enzyme]))
        DataHandler['y_train'] = train_seq.y[config.enzyme][train_indexes]
        DataHandler['X_train'] = train_seq.X[train_indexes]
        DataHandler['X_train_biofeat'] = train_seq.X_biofeat[train_indexes]
        if config.weighted_loss:
            DataHandler['conf_train'] = train_seq.confidence[config.enzyme][train_indexes]
    else:
        DataHandler['y_train'] = train_seq.y
        DataHandler['X_train'] = train_seq.X
        DataHandler['X_train_biofeat'] = train_seq.X_biofeat


### Parallel loader for multi_task
def parallel_enzyme_loader(config, DataHandler, test_seq, valid_seq, train_seq):
    # The test data is being tested for every enzyme separately, so we can use nan values
    y_test, _ = parallel_target(test_seq, drop_nans=False)
    DataHandler['y_test'] = y_test
    DataHandler['X_test'] = test_seq.X
    DataHandler['X_test_biofeat'] = test_seq.X_biofeat

    # The validation data cannot hold nan values in the target
    y_valid, valid_indexes = parallel_target(valid_seq, drop_nans=True)
    DataHandler['y_valid'] = y_valid
    DataHandler['X_valid'] = valid_seq.X[valid_indexes]
    DataHandler['X_valid_biofeat'] = valid_seq.X_biofeat[valid_indexes]

    # The train data cannot hold nan values in the target
    if config.data_source == 'old':
        y_train, train_indexes = parallel_target(train_seq, drop_nans=True)
        DataHandler['y_train'] = y_train
        DataHandler['X_train'] = train_seq.X[train_indexes]
        DataHandler['X_train_biofeat'] = train_seq.X_biofeat[train_indexes]

    else:
        DataHandler['X_train'] = train_seq.X
        DataHandler['X_train_biofeat'] = train_seq.X_biofeat
        DataHandler['y_train'] = train_seq.y


def parallel_target(data_seq, drop_nans):
    if drop_nans:
        wt_valid_indexes = np.logical_not(np.isnan(data_seq.y['wt']))
        esp_valid_indexes = np.logical_not(np.isnan(data_seq.y['esp']))
        hf_valid_indexes = np.logical_not(np.isnan(data_seq.y['hf']))
        valid_indexes = wt_valid_indexes & esp_valid_indexes & hf_valid_indexes

    else:
        # This is a manipulation to make all 3 enzymes valid_indexes valid
        valid_indexes = np.ones(len(data_seq.y['wt']), dtype=bool)

    y_data = data_seq.y['wt'][valid_indexes]
    y_data = np.expand_dims(y_data, 1)
    y_data = np.append(y_data, np.expand_dims(data_seq.y['esp'][valid_indexes], axis=1), axis=1)
    y_data = np.append(y_data, np.expand_dims(data_seq.y['hf'][valid_indexes], axis=1), axis=1)

    return y_data, valid_indexes


### Serial loader for multi_task
def serial_enzyme_loader(config, DataHandler, test_seq, valid_seq, train_seq):
    # The test data is being tested for every enzyme separately, so we can use nan values
    y_test, _ = parallel_target(test_seq, drop_nans=False)
    DataHandler['y_test'] = y_test
    DataHandler['X_test'] = test_seq.X
    DataHandler['X_test_biofeat'] = test_seq.X_biofeat

    # The validation data need to be filtered to hold only not nan values. Because the data is being entered serially,
    # we can filter the not nan data for each enzyme separately.
    enzymes = ['wt', 'esp', 'hf']
    y_valid = np.empty(shape=[0, ], dtype=np.float16)
    conf_valid = np.empty(shape=[0, ])
    X_valid = np.empty(shape=[0, 22])
    X_biofeat_valid = np.empty(shape=[0, 14])

    for ind, enzyme in enumerate(enzymes):
        valid_indexes = np.logical_not(np.isnan(valid_seq.y[enzyme]))
        y_valid = np.concatenate((y_valid, valid_seq.y[enzyme][valid_indexes]), axis=0)
        conf_valid = np.concatenate((conf_valid, valid_seq.confidence[enzyme][valid_indexes]), axis=0)
        X_valid = np.concatenate((X_valid, valid_seq.X[valid_indexes]), axis=0)
        x_biofeat_valid_enzyme = valid_seq.X_biofeat[valid_indexes]
        ohe = np.eye(3)[np.zeros((x_biofeat_valid_enzyme.shape[0]), dtype=np.int) + ind]
        x_biofeat_valid_enzyme = np.append(x_biofeat_valid_enzyme, ohe, axis=-1)
        X_biofeat_valid = np.concatenate((X_biofeat_valid, x_biofeat_valid_enzyme), axis=0)

    DataHandler['y_valid'] = y_valid
    DataHandler['conf_valid'] = conf_valid
    DataHandler['X_valid'] = X_valid
    DataHandler['X_valid_biofeat'] = X_biofeat_valid

    # The train data can be loaded same as the validation data.
    # For old data - we need to filter the not nan values
    # For the new data (raw reads) the data is already filtered - we need only to combine the 3 files.
    y_train = np.empty(shape=[0, ], dtype=np.float16)
    conf_train = np.empty(shape=[0, ])
    X_train = np.empty(shape=[0, 22])
    X_biofeat_train = np.empty(shape=[0, 14])
    for ind, enzyme in enumerate(enzymes):
        valid_indexes = np.logical_not(np.isnan(train_seq.y[enzyme]))
        y_train = np.concatenate((y_train, train_seq.y[enzyme][valid_indexes]), axis=0)
        conf_train = np.concatenate((conf_train, train_seq.confidence[enzyme][valid_indexes]), axis=0)
        X_train = np.concatenate((X_train, train_seq.X[valid_indexes]), axis=0)
        x_biofeat_train_enzyme = train_seq.X_biofeat[valid_indexes]
        ohe = np.eye(3)[np.zeros((x_biofeat_train_enzyme.shape[0]), dtype=np.int) + ind]
        x_biofeat_train_enzyme = np.append(x_biofeat_train_enzyme, ohe, axis=-1)
        X_biofeat_train = np.concatenate((X_biofeat_train, x_biofeat_train_enzyme), axis=0)

    DataHandler['y_train'] = y_train
    DataHandler['conf_train'] = conf_train
    DataHandler['X_train'] = X_train
    DataHandler['X_train_biofeat'] = X_biofeat_train


# def train_val_split(X,X_biofeat,y, val_size = 0.1,random_state=40):
#     X_train, X_val, y_train, y_val = train_test_split(
#        X, y, test_size=val_size, random_state=random_state)
#
#     X_train_biofeat, X_val_biofeat, y_train, y_val = train_test_split(
#        X_biofeat, y, test_size=val_size, random_state=random_state)
#
#     return X_train, X_val, X_train_biofeat, X_val_biofeat, y_train, y_val


def eda(DataHandler):
    import matplotlib.pyplot as plt

    from matplotlib import pyplot
    train_conf, valid_conf = DataHandler['conf_train'], DataHandler['conf_valid']

    plt.figure(figsize=(40, 20))
    for i in range(2, 5):
        pyplot.subplot(6, 2, i + 1)
        pyplot.hist(train_conf, range=(1, 10 ** i))
        pyplot.title(f'conf_train range: {1}:{10 ** i}')

    for i in range(2, 5):
        pyplot.subplot(6, 2, i + 4)
        pyplot.hist(valid_conf, range=(1, 10 ** i))
        pyplot.title(f'conf_valid range: {1}:{10 ** i}')

    pyplot.savefig('confidence.png')
    exit(1)

    biofeatures = DataHandler['X_train_biofeat']

    plt.figure(figsize=(40, 30))
    for ind in range(11):
        biofeature = biofeatures[:, ind]
        pyplot.subplot(11, 2, ind + 1)
        pyplot.hist(biofeature)
        pyplot.title('feature' + str(ind))

    pyplot.savefig('features_EDA.png')

    plt.figure()
    pyplot.hist(DataHandler['y_train'], label=['wt', 'esp', 'hf'])
    pyplot.legend()
    pyplot.savefig('efficiency_EDA.png')

    from sklearn.preprocessing import MinMaxScaler, StandardScaler
    scaler = StandardScaler()
    scaler.fit(DataHandler['y_train'])
    y_train = scaler.transform(DataHandler['y_train'])
    plt.figure()
    pyplot.hist(y_train, label=['wt', 'esp', 'hf'])
    pyplot.legend()
    pyplot.savefig('efficiency_EDA_transformed.png')
    exit(1)
