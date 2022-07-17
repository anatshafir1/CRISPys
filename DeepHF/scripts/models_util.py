import keras
from keras.preprocessing import text
from keras.preprocessing import sequence
from keras.layers import merge, Embedding, Bidirectional
from keras.layers.core import *
from keras.models import *
from keras.layers.recurrent import LSTM
from keras.callbacks import ModelCheckpoint, EarlyStopping, LearningRateScheduler
from keras.optimizers import *
import numpy as np
from keras.callbacks import Callback
from sklearn.metrics import mean_squared_error, log_loss
import scipy as sp
import os
import pickle
import tensorflow as tf
import tensorflow.keras.backend as kb

class GetBest(Callback):
    def __init__(self, filepath=None, monitor='val_loss', save_best=False, verbose=0,
                 mode='auto', period=1, val_data=None, checks=['on_epoch_end']):
        super(GetBest, self).__init__()
        self.monitor = monitor
        self.verbose = verbose
        self.period = period
        self.save_best = save_best
        self.filepath = filepath
        self.best_epochs = 0
        self.epochs_since_last_save = 0
        self.val_data = val_data
        self.checks = checks
        self.patience = 5
        self.init_lr = 0.002
        self.lr = self.init_lr

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
        if 'on_epoch_end' in self.checks:
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
                        self.patience = 5
                        # self.model.save(filepath, overwrite=True)
                    else:
                        self.patience -= 1
                        if self.patience == 0:
                            self.patience = 5
                            self.lr /= 2
                            self.model.optimizer.learning_rate.assign(self.lr)
                            print(f'updating learning rate: {self.lr * 2} -> {self.lr}')
                        if self.verbose > 0:
                            print('\nEpoch %05d: %s did not improve.' %
                                  (epoch + 1, self.monitor))
                            print(f'patience: {self.patience}')

    def on_batch_end(self, batch, logs=None):
        if 'on_batch_end' in self.checks:
            if batch%100 == 0:
                loss = logs.get('loss')
                y_test_pred = self.model.predict(self.val_data['val_input'])
                val_loss = mean_squared_error(self.val_data['y_val'], y_test_pred)
                val_spearman = sp.stats.spearmanr(self.val_data['y_val'], y_test_pred)[0]
                if (val_spearman > self.best) or (self.best==np.Inf) :
                    self.best = val_spearman
                    self.best_weights = self.model.get_weights()

                print(f'batch: {batch} , loss: {loss}, val_loss: {val_loss}, val_spearman: {val_spearman}, best: {self.best}')

                a=0


    def on_train_end(self, logs=None):
        if self.verbose > 0:
            print('Using epoch %05d with %s: %0.5f.' % (self.best_epochs, self.monitor,
                                                        self.best))
        self.model.set_weights(self.best_weights)
        if self.save_best:
            self.model.save(self.filepath, overwrite=True)



def scheduler(epoch, lr):
    learning_rate = 0.01  # initial learning rate
    decay_rate = 0.03
    lrate = learning_rate * np.exp(-decay_rate * epoch)
    return lrate


def lstm_model(config, DataHandler):

    X_train, X_train_biofeat, y_train = DataHandler['X_train'], DataHandler['X_train_biofeat'], DataHandler['y_train']
    X_val, X_val_biofeat, y_val = DataHandler['X_valid'], DataHandler['X_valid_biofeat'], DataHandler['y_valid']


    if config.has_embedding:
        sequence_input = Input(name='seq_input', shape=(22,))
        embedding_layer = Embedding(5, config.em_dim, input_length=22)
        embedded = embedding_layer(sequence_input)
        embedded = SpatialDropout1D(config.em_drop)(embedded)
        x = embedded

    else:
        sequence_input = Input(name='seq_input', shape=(22, 5))
        x = sequence_input

    # RNN
    lstm = LSTM(config.rnn_units, dropout=config.rnn_drop,
                kernel_regularizer='l2', recurrent_regularizer='l2',
                recurrent_dropout=config.rnn_rec_drop, return_sequences=True, kernel_initializer=config.initializer)

    x = Bidirectional(lstm)(x)
    x = Flatten()(x)

    # Biological features
    if config.has_biofeatures:
        biological_input = Input(name = 'bio_input', shape = (X_train_biofeat.shape[1],))
        x = keras.layers.concatenate([x, biological_input])


    for l in range(config.fc_num_hidden_layers):
        x = Dense(config.fc_num_units, activation=config.fc_activation, kernel_initializer=config.initializer)(x)
        x = Dropout(config.fc_drop)(x)

    # Finish model
    if config.enzyme == 'multi_task' and config.multi_data == 'parallel':
        mix_output = Dense(3, activation=config.last_activation, name='mix_output', kernel_initializer=config.initializer)(x)
    else:
        mix_output = Dense(1, activation=config.last_activation, name='mix_output', kernel_initializer=config.initializer)(x)


    if config.has_biofeatures:
        model = Model(inputs=[sequence_input, biological_input], outputs=[mix_output])
    else:
        model = Model(inputs=[sequence_input], outputs=[mix_output])


    # model.compile(loss=config.cost_function, optimizer=config.optimizer(lr=config.learning_rate))
    model.compile(loss=config.cost_function, optimizer=config.optimizer(lr=0.002))



    model.summary()
    np.random.seed(1337)

    callback_list = []
    early_stopping = EarlyStopping(monitor='val_loss', patience=15, verbose=1)
    callback_list.append(early_stopping)



    if config.lr_scheduler:
        lr_scheduler = LearningRateScheduler(scheduler)
        callback_list.append(lr_scheduler)

    # Prepare data

    if config.input_scale is not None:
        config.input_scale.fit(X_train_biofeat)
        X_train_biofeat = config.input_scale.transform(X_train_biofeat)
        X_val_biofeat = config.input_scale.transform(X_val_biofeat)

    if config.output_scale is not None:
        config.output_scale.fit(y_train)
        y_train = config.output_scale.transform(y_train)
        y_val = config.output_scale.transform(y_val)



    if config.has_biofeatures:
        train_input = [X_train, X_train_biofeat]
        val_input = [X_val, X_val_biofeat]
    else:
        train_input = [X_train]
        val_input = [X_val]

    save_dir = f'models/{config.model_type}/'
    if config.has_biofeatures:
        save_dir += 'bio/'
    else:
        save_dir += 'no_bio/'
    if config.enzyme == 'multi_task':
        save_dir += f'{config.enzyme}/{config.multi_data}/'
    else:
        save_dir += f'{config.enzyme}/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    get_best_model = GetBest(filepath=f'{save_dir}model' ,val_data={'val_input': val_input, 'y_val': y_val}, verbose=1, mode='min',
                             save_best=config.save_model, checks=config.checks)
    callback_list.append(get_best_model)

    save_scales(config, save_dir)

    loss_weights = DataHandler['conf_train']
    loss_weights_val = DataHandler['conf_valid']

    weight_scale = np.mean(loss_weights)
    loss_weights = loss_weights / weight_scale
    loss_weights_val = loss_weights_val / weight_scale

    return model, train_input, y_train, val_input, y_val, callback_list, config, loss_weights, loss_weights_val

def save_scales(config, save_dir):
    with open(save_dir + 'config.pkl', "wb") as fp:
        pickle.dump(config, fp)


def fully_connected_model(config, DataHandler):

    X_train, X_train_biofeat, y_train = DataHandler['X_train'], DataHandler['X_train_biofeat'], DataHandler['y_train']
    X_val, X_val_biofeat, y_val = DataHandler['X_valid'], DataHandler['X_valid_biofeat'], DataHandler['y_valid']


    if config.has_embedding:
        sequence_input = Input(name='seq_input', shape=(22,))
        embedding_layer = Embedding(5, config.em_dim, input_length=22)
        embedded = embedding_layer(sequence_input)
        embedded = SpatialDropout1D(config.em_drop)(embedded)
        x = embedded

    else:
        sequence_input = Input(name='seq_input', shape=(22, 5))
        x = sequence_input

    # RNN
    # lstm = LSTM(config.rnn_units, dropout=config.rnn_drop,
    #             kernel_regularizer='l2', recurrent_regularizer='l2',
    #             recurrent_dropout=config.rnn_rec_drop, return_sequences=True, kernel_initializer=config.initializer)
    #
    # x = Bidirectional(lstm)(x)
    x = Flatten()(x)

    # Biological features
    if config.has_biofeatures:
        biological_input = Input(name = 'bio_input', shape = (X_train_biofeat.shape[1],))
        x = keras.layers.concatenate([x, biological_input])


    for l in range(config.fc_num_hidden_layers):
        x = Dense(config.fc_num_units, activation=config.fc_activation, kernel_initializer=config.initializer)(x)
        x = Dropout(config.fc_drop)(x)

    # Finish model
    if config.enzyme == 'multi_task' and config.multi_data == 'parallel':
        mix_output = Dense(3, activation=config.last_activation, name='mix_output', kernel_initializer=config.initializer)(x)
    else:
        mix_output = Dense(1, activation=config.last_activation, name='mix_output', kernel_initializer=config.initializer)(x)


    if config.has_biofeatures:
        model = Model(inputs=[sequence_input, biological_input], outputs=[mix_output])
    else:
        model = Model(inputs=[sequence_input], outputs=[mix_output])
    model.compile(loss=config.cost_function, optimizer=config.optimizer(lr=config.learning_rate))


    model.summary()
    np.random.seed(1337)

    callback_list = []
    early_stopping = EarlyStopping(monitor='val_loss', patience=15, verbose=1)
    callback_list.append(early_stopping)



    if config.lr_scheduler:
        lr_scheduler = LearningRateScheduler(scheduler)
        callback_list.append(lr_scheduler)

    # Prepare data

    if config.input_scale is not None:
        config.input_scale.fit(X_train_biofeat)
        X_train_biofeat = config.input_scale.transform(X_train_biofeat)
        X_val_biofeat = config.input_scale.transform(X_val_biofeat)

    if config.output_scale is not None:
        config.output_scale.fit(y_train)
        y_train = config.output_scale.transform(y_train)
        y_val = config.output_scale.transform(y_val)



    if config.has_biofeatures:
        train_input = [X_train, X_train_biofeat]
        val_input = [X_val, X_val_biofeat]
    else:
        train_input = [X_train]
        val_input = [X_val]


    get_best_model = GetBest('models/' + config.enzyme + '_rnn.hd5',val_data={'val_input': val_input, 'y_val': y_val}, verbose=1, mode='min',
                             save_best=config.save_model, checks=config.checks)
    callback_list.append(get_best_model)

    return model, train_input, y_train, val_input, y_val, callback_list, config

