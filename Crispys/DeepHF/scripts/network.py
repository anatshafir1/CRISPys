from keras.optimizers import *
import keras
from keras.layers import merge, Embedding, Bidirectional, TimeDistributed
from keras.models import *
from keras.layers.core import *
from keras.layers.recurrent import GRU, LSTM
from keras.callbacks import Callback
from keras.callbacks import ModelCheckpoint, EarlyStopping
import numpy as np


fc_activation_dict = {'1': 'relu', '2': 'tanh', '3': 'sigmoid', '4': 'hard_sigmoid', '0': 'elu'}
initializer_dict = {'1': 'lecun_uniform', '2': 'normal', '3': 'he_normal', '0': 'he_uniform'}
optimizer_dict = {'1': SGD, '2': RMSprop, '3': Adagrad, '4': Adadelta, '5': Adam, '6': Adamax, '0': Nadam}

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


def lstm_model(batch_size=90, epochs=50, initializer='0', em_dim=44, em_drop=0.2,
               rnn_units=60, rnn_drop=0.6, rnn_rec_drop=0.1, fc_num_hidden_layers=3,
               fc_num_units=320, fc_drop=0.4, fc_activation=None, optimizer=Adam, learning_rate=0.001,
               validation_split=0.1, shuffle=False, data=None, data_type=None):

    initializer = initializer_dict[str(initializer)]
    sequence_input = Input(name='seq_input', shape=(22,))

    embedding_layer = Embedding(5, em_dim, input_length=22)
    embedded = embedding_layer(sequence_input)
    embedded = SpatialDropout1D(em_drop)(embedded)
    x = embedded

    # RNN
    lstm = LSTM(rnn_units, dropout=rnn_drop,
                kernel_regularizer='l2', recurrent_regularizer='l2',
                recurrent_dropout=rnn_rec_drop, return_sequences=True)
    x = Bidirectional(lstm)(x)
    x = Flatten()(x)

    biological_input = Input(name='bio_input', shape=(data['X_train_biofeat'].shape[1],))
    x = keras.layers.concatenate([x, biological_input])

    for l in range(fc_num_hidden_layers):
        x = Dense(fc_num_units, activation=fc_activation)(x)
        x = Dropout(fc_drop)(x)
    # finish model
    mix_output = Dense(1, activation='linear', name='mix_output')(x)

    model = Model(inputs=[sequence_input, biological_input], outputs=[mix_output])
    # model = Model(inputs=[sequence_input], outputs=[mix_output])
    model.compile(loss='mse', optimizer=optimizer(lr=0.001))

    np.random.seed(1337)
    early_stopping = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
    get_best_model = GetBest('Model/hf_rnn.hd5', monitor='val_loss', verbose=1, mode='min')
    print('Beggining training')
    model.fit([data['X_train'], data['X_train_biofeat']],
              # model.fit([X_train],
              data['y_train'],
              batch_size=batch_size,
              epochs=epochs,
              verbose=2,
              validation_split=0.1,
              shuffle=False,
              callbacks=[get_best_model, early_stopping])
    return model
