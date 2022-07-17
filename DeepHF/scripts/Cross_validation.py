import scipy as sp
from sklearn.metrics import mean_squared_error, r2_score

import keras
from keras.layers import merge, Embedding, Bidirectional, TimeDistributed
from keras.layers.convolutional import Conv1D
from keras.layers.core import *
from keras.models import *
from keras.layers.pooling import MaxPooling1D
from keras.layers.recurrent import GRU, LSTM
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.optimizers import *

from scripts.feature_util import *
from sklearn.model_selection import KFold, ShuffleSplit, StratifiedKFold, GroupKFold
from sklearn import linear_model

import numpy as np
from keras.callbacks import Callback
import logging
import time
import os
from scripts import network
from scripts import training_util

class GetBest(Callback):
    """Get the best model at the end of training.
	# Arguments
        monitor: quantity to monitor.
        verbose: verbosity mode, 0 or 1.
        mode: one of {auto, min, max}.
            The decision
            to overwrite the current stored weights is made
            based on either the maximization or the
            minimization of the monitored quantity. For `val_acc`,
            this should be `max`, for `val_loss` this should
            be `min`, etc. In `auto` mode, the direction is
            automatically inferred from the name of the monitored quantity.
        period: Interval (number of epochs) between checkpoints.
	# Example
		callbacks = [GetBest(monitor='val_acc', verbose=1, mode='max')]
		mode.fit(X, y, validation_data=(X_eval, Y_eval),
                 callbacks=callbacks)
    """

    def __init__(self, filepath=None, monitor='val_loss', save_best=False, verbose=0,
                 mode='auto', period=1):
        super(GetBest, self).__init__()
        self.monitor = monitor
        self.verbose = verbose
        self.period = period
        self.save_best = save_best
        self.filepath = filepath
        self.best_epochs = 0
        self.epochs_since_last_save = 0

        if mode not in ['auto', 'min', 'max']:
            warnings.warn('GetBest mode %s is unknown, '
                          'fallback to auto mode.' % (mode),
                          RuntimeWarning)
            mode = 'auto'

        if mode == 'min':
            self.monitor_op = np.less
            self.best = np.Inf
        elif mode == 'max':
            self.monitor_op = np.greater
            self.best = -np.Inf
        else:
            if 'acc' in self.monitor or self.monitor.startswith('fmeasure'):
                self.monitor_op = np.greater
                self.best = -np.Inf
            else:
                self.monitor_op = np.less
                self.best = np.Inf

    def on_train_begin(self, logs=None):
        self.best_weights = self.model.get_weights()

    def on_epoch_end(self, epoch, logs=None):
        logs = logs or {}
        self.epochs_since_last_save += 1
        if self.epochs_since_last_save >= self.period:
            self.epochs_since_last_save = 0
            filepath = self.filepath.format(epoch=epoch + 1, **logs)
            current = logs.get(self.monitor)
            if current is None:
                warnings.warn('Can pick best model only with %s available, '
                              'skipping.' % (self.monitor), RuntimeWarning)
            else:
                if self.monitor_op(current, self.best):
                    if self.verbose > 0:
                        print('\nEpoch %05d: %s improved from %0.5f to %0.5f,'
                              ' storing weights.'
                              % (epoch + 1, self.monitor, self.best,
                                 current))
                    self.best = current
                    self.best_epochs = epoch + 1
                    self.best_weights = self.model.get_weights()
                    # self.model.save(filepath, overwrite=True)
                else:
                    if self.verbose > 0:
                        print('\nEpoch %05d: %s did not improve.' %
                              (epoch + 1, self.monitor))

    def on_train_end(self, logs=None):
        if self.verbose > 0:
            print('Using epoch %05d with %s: %0.5f.' % (self.best_epochs, self.monitor,
                                                        self.best))
        self.model.set_weights(self.best_weights)
        # self.model.save(filepath, overwrite=True)


fc_activation_dict = {'1': 'relu', '2': 'tanh', '3': 'sigmoid', '4': 'hard_sigmoid', '0': 'elu'}
initializer_dict = {'1': 'lecun_uniform', '2': 'normal', '3': 'he_normal', '0': 'he_uniform'}
optimizer_dict = {'1': SGD, '2': RMSprop, '3': Adagrad, '4': Adadelta, '5': Adam, '6': Adamax, '0': Nadam}

#
def load_data_kf(X, X_biofeat, y, split=0.15):
    train_valid_data = []
    kf = ShuffleSplit(n_splits=10, test_size=0.1, random_state=33)
    for train_index, valid_index in kf.split(X):
        X_train, X_train_biofeat, X_valid, X_valid_biofeat = X[train_index], X_biofeat[train_index], X[valid_index], \
                                                           X_biofeat[valid_index]
        y_train, y_valid = y[train_index], y[valid_index]

        train_valid_data.append((X_train, X_train_biofeat, X_valid, X_valid_biofeat, y_train, y_valid))
    return train_valid_data




def cross_validation(param, DataHandler):
    cross_v_start = time.time()
    enzime = param['data_type']
    # create logger with 'Model_application'
    logger = logging.getLogger('Model')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    ind = 1
    logfile_name = 'logs/' + enzime + '_CrossValidation_1.log'
    while os.path.exists(logfile_name):
        ind += 1
        logfile_name = 'logs/' + enzime + '_CrossValidation_' + str(ind) + '.log'

    fh = logging.FileHandler(logfile_name)
    fh.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info('###################################################')
    logger.info('Running  10-fold cross-validation for data type {}'.format(enzime))
    logger.info('###################################################')
    logger.info('----------')
    logger.info('params:{}'.format(param))
    logger.info('----------')

    # # load data
    # with open('data/' + param['data_type'] + '_seq_data_array.pkl', 'rb') as handle:
    #     X, X_biofeat, y = pickle.load(handle)
    #     X = X - 1
    # logger.info('Data size: {}'.format(len(X)))
    X_test, X_test_biofeat, y_test = DataHandler['X_test'], DataHandler['X_test_biofeat'], DataHandler['y_test']
    X_train, X_train_biofeat, y_train = DataHandler['X_train'], DataHandler['X_train_biofeat'], DataHandler['y_train']
    data_list = load_data_kf(X_train, X_train_biofeat, y_train)
    # training_r2_scores = []
    # testing_r2_scores = []

    Data = {}
    Data['X_test'] = X_test
    Data['X_test_biofeat'] = X_test_biofeat
    Data['y_test'] = y_test
    if param['data_type'] == 'multi_task':
        spearmanrs = {'wt': [], 'esp': [], 'hf': []}
        loss = {'wt': [], 'esp': [], 'hf': []}
    else:
        spearmanrs = []
        loss = []
    for index, data in enumerate(data_list):
        logger.info('Beginning Cross-Validation num: {}'.format(index))
        train_start = time.time()

        X_train, X_train_biofeat, X_val, X_val_biofeat, y_train, y_val = data[0], data[1], data[2], data[3], data[4], \
                                                                            data[5]


        Data['X_train'] = X_train
        Data['X_train_biofeat'] = X_train_biofeat
        Data['y_train'] = y_train
        Data['X_val'] = X_val
        Data['X_val_biofeat'] = X_val_biofeat
        Data['y_val'] = y_val

        model, output_scale, input_scale = training_util.lstm_model(**param, data=Data)

        # prepare data
        y_test = Data['y_test']
        if output_scale is not None:
            y_test = output_scale.transform(y_test)
        if input_scale is not None:
            X_test_biofeat = input_scale.transform(X_test_biofeat)

        if param['data_type'] == 'multi_task':
            y_test_pred = model.predict([Data['X_test'], Data['X_test_biofeat']])
            data_types = ['wt', 'esp', 'hf']
            for ind, type in enumerate(data_types):
                y_test_true_enzyme = y_test[:, ind]
                y_test_pred_enzyme = y_test_pred[:, ind]
                mse = mean_squared_error(y_test_true_enzyme, y_test_pred_enzyme)
                spearmanr = sp.stats.spearmanr(y_test_true_enzyme, y_test_pred_enzyme)[0]
                logger.info('     type: {}, spearman: {}, mse: {}'.format(type, spearmanr, mse))
                spearmanrs[type].append(spearmanr)
                loss[type].append(mse)
        else:
            y_test_pred = model.predict([Data['X_test'], Data['X_test_biofeat']])
            mse = mean_squared_error(y_test, y_test_pred)
            spearmanr = sp.stats.spearmanr(y_test, y_test_pred)[0]
            logger.info('     spearman: {}'.format(spearmanr))
            logger.info('     mse: {}'.format(mse))

            spearmanrs.append(spearmanr)
            loss.append(mse)


        # y_train_pred = model.predict([X_train, X_train_biofeat])
        # y_test_pred = model.predict([X_test, X_test_biofeat])
        # training_score = r2_score(y_train, y_train_pred)
        # testing_score = r2_score(y_test, y_test_pred)
        # mse = mean_squared_error(y_test, y_test_pred)
        # spearmanr = sp.stats.spearmanr(y_test, y_test_pred)[0]
        train_end = time.time()

        # logger.info('   mse: {}'.format(mse))
        # logger.info('   spearmanr: {}'.format(spearmanr))
        # logger.info('   testing_score: {}'.format(testing_score))
        logger.info('   train time: {}'.format(train_end - train_start))

        # training_r2_scores.append(training_score)
        # testing_r2_scores.append(testing_score)



    if param['data_type'] == 'multi_task':
        r = {'wt': {}, 'esp': {}, 'hf': {}}
        data_types = ['wt', 'esp', 'hf']
        for ind, type in enumerate(data_types):
            r[type]['spearmanr'] = np.mean(spearmanrs[type]), np.std(spearmanrs[type])
            r[type]['loss'] = np.mean(loss[type]), np.std(loss[type])
        logger.info('{}'.format(r))

    else:
        r = {}
        r['spearmanr'] = np.mean(spearmanrs), np.std(spearmanrs)
        r['loss'] = np.mean(loss), np.std(loss)
        logger.info('{}'.format(r))
    r = {}
    # r['training_r2_scores'] = np.mean(training_r2_scores), np.std(training_r2_scores)
    # r['testing_r2_scores'] = np.mean(testing_r2_scores), np.std(testing_r2_scores)
    # r['spearmanr'] = np.mean(spearmanrs), np.std(spearmanrs)
    # r['loss'] = np.mean(loss), np.std(loss)
    # logger.info('{}'.format(r))

    cross_v_end = time.time()
    logger.info('Cross Validation time: {}'.format(cross_v_end - cross_v_start))



