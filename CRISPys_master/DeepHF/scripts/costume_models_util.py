import tensorflow as tf
import tensorflow.keras.layers as layers
import numpy as np
import scipy as sp
import time
from tqdm import tqdm
import os
import pickle


# Build the model
class lstm(tf.keras.Model):

    def __init__(self, config):
        super(lstm, self).__init__()
        # tf.keras.backend.set_floatx('float16')
        # Configurations:
        self.outputs = 3 if (config.enzyme=='multi_task' and config.multi_data == 'parallel') else 1
        self.best_val_loss = np.Inf
        self.patience = 5
        self.early_stop = False
        self.init_lr = 0.002
        self.lr = self.init_lr
        self.saved_weights = None
        self.save_best = config.save_model

        model_path = f'models/{config.model_type}/bio/{config.enzyme}'
        model_num = 0
        while os.path.exists(f'{model_path}/model_{model_num}'):
            model_num += 1
        self.filepath = f'{model_path}/model_{model_num}'
        # Architecture:
        self.embedding = layers.Embedding(5, config.em_dim, input_length=22)
        self.embedding_drop = layers.SpatialDropout1D(config.em_drop)
        self.lstm = layers.LSTM(config.rnn_units, dropout=config.rnn_drop,
                                kernel_regularizer='l2', recurrent_regularizer='l2',
                                recurrent_dropout=config.rnn_rec_drop, return_sequences=True,
                                kernel_initializer=config.initializer)
        self.bidirectional = layers.Bidirectional(self.lstm)
        self.flatten = layers.Flatten()
        # self.concatenate = layers.concatenate([self.flatten, biological_input])
        self.dense_layers = []
        self.dense_drop_layers = []
        for l in range(config.fc_num_hidden_layers):
            self.dense_layers.append(layers.Dense(config.fc_num_units, activation=config.fc_activation,
                                                     kernel_initializer=config.initializer, kernel_regularizer='l2'))
            self.dense_drop_layers.append(layers.Dropout(config.fc_drop))


        self.out_layer = layers.Dense(self.outputs, activation=config.last_activation, name='mix_output',
                           kernel_initializer=config.initializer)


        self.trainable_layers = [self.embedding, self.bidirectional] + self.dense_layers + [self.out_layer]
        self.train_op = config.optimizer(lr=self.init_lr)



    @tf.function
    def predict(self, inputs, training=False):
        sequence = inputs[0]
        bio_feat = inputs[1]
        x = self.embedding(sequence, training=training)
        x = self.embedding_drop(x, training=training)
        x = self.bidirectional(x, training=training)
        x = self.flatten(x)
        x = layers.concatenate([x, bio_feat], dtype=np.float32)
        for dense, drop in zip(self.dense_layers, self.dense_drop_layers):
            x = dense(x, training=training)
            x = drop(x, training=training)
        y_pred = self.out_layer(x, training=training)
        return y_pred

    # Custom loss fucntion
    # @tf.autograph.experimental.do_not_convert
    @tf.function
    def get_loss(self, inputs, y_true, loss_weights, training=True):

        y_pred = self.predict(inputs, training=training)
        y_pred = tf.dtypes.cast(y_pred, tf.float16)
        if loss_weights is not None:
            loss = y_true - y_pred
            weighted_loss = tf.math.multiply(loss, loss_weights)
            weighted_loss = tf.math.square(weighted_loss)
            weighted_loss_mean = tf.math.reduce_mean(weighted_loss)
            return weighted_loss_mean
        else:
            return tf.keras.losses.MSE(y_true, y_pred)


    # get gradients
    def get_grad(self, X, Y, loss_weights):
        with tf.GradientTape() as tape:
            for trainable_layer in self.trainable_layers:
                tape.watch(trainable_layer.variables)
            L = self.get_loss(X, Y, loss_weights)
            weights = [var for layer in self.trainable_layers for var in layer.variables]
            g = tape.gradient(L, weights)

        return g, weights

    # perform gradient descent
    @tf.function
    def network_learn(self, X, Y, loss_weights):
        g, weights = self.get_grad(X, Y, loss_weights)
        # print(self.var)
        self.train_op.apply_gradients(zip(g, weights))

    def on_epoch_end(self):

        val_loss = self.get_loss(self.X_val, self.Y_val, self.sample_weight_valid, training=False).numpy()
        print(f'val_loss: {val_loss}')
        # for ind in range(self.outputs):

            # spearman.append(sp.stats.spearmanr(self.Y_val[:, ind], y_pred_val[:, ind])[0])
        # val_avg = sum(spearman) / len(spearman)
        # print(f'validation spearman: {spearman} mean: {val_avg}')

        if val_loss < self.best_val_loss:
            self.best_val_loss = val_loss
            self.saved_weights = self.get_weights()
            self.patience = 5
        else:
            self.patience -= 1

        if self.patience == 0:
            if self.lr == self.init_lr / 8:
                self.early_stop = True
                print('early stoping')
            else:
                self.patience = 5
                self.lr /= 2
                self.train_op.learning_rate.assign(self.lr)
                print('changing learning rate')


    def fit(self, X, Y, batch_size, epochs, verbose, validation_data, shuffle, callbacks, sample_weight=None):
        X_val = validation_data[0]
        Y_val = validation_data[1]
        sample_weight_valid = validation_data[2]

        num_batches = X[0].shape[0]//batch_size
        if len(Y.shape) == 1:
            Y = np.expand_dims(Y, axis=1)
            Y_val = np.expand_dims(Y_val, axis=1)
            if sample_weight is not None:
                sample_weight = np.expand_dims(sample_weight, axis=1).astype(np.float16)
                sample_weight_valid = np.expand_dims(sample_weight_valid, axis=1).astype(np.float16)

        sample_weight_batch = None
        self.X, self.Y = X, Y
        self.X_val, self.Y_val = X_val, Y_val
        self.sample_weight_valid = sample_weight_valid

        for epoch in range(epochs):
            if self.early_stop:
                break
            print(f'epoch{epoch}')
            epoch_start = time.time()
            # Shuffle to receive new batches permutation
            if shuffle:
                perm = np.random.permutation(len(X[0]))
            else:
                perm = range(len(X))
            perm_batches = np.array_split(perm, num_batches)
            # progbar = tf.keras.utils.Progbar(len(perm_batches))
            for ind, batch_indexes in zip(tqdm(range(len(perm_batches))), perm_batches):
            # for batch_indexes in perm_batches:

                # progbar.update(ind+1)
                X_batch = [X[0][batch_indexes], X[1][batch_indexes]]
                Y_batch = Y[batch_indexes]
                if sample_weight is not None:
                    sample_weight_batch = sample_weight[batch_indexes]
                self.network_learn(X_batch, Y_batch, sample_weight_batch)


            self.on_epoch_end()
            epoch_end = time.time()
            print(f'time: {int(epoch_end - epoch_start)}')

        self.on_train_end()

    def on_train_end(self):
        if self.saved_weights is not None:
            self.set_weights(self.saved_weights)
        if self.save_best:
            with open(self.filepath, "wb") as fp:
                pickle.dump(self.saved_weights, fp)
            # tf.saved_model.save(self, self.filepath)

    # def set_weights(self, saved_weights):
    #     self.set_weights(self.saved_weights)



def train(config, DataHandler):
    train_input = [DataHandler['X_train'], DataHandler['X_train_biofeat']]
    train_output = DataHandler['y_train']

    valid_input = [DataHandler['X_valid'], DataHandler['X_valid_biofeat']]
    valid_output = DataHandler['y_valid']

    loss_weights = DataHandler['conf_train']
    loss_scaling = np.mean(loss_weights)
    loss_weights = loss_weights / loss_scaling
    loss_weights = np.sqrt(loss_weights)
    loss_weights_valid = DataHandler['conf_valid']
    loss_weights_valid = loss_weights_valid / loss_scaling
    loss_weights_valid = np.sqrt(loss_weights_valid)


    # valid_input = [DataHandler['X_test'], DataHandler['X_test_biofeat']] #Debug
    # valid_output = DataHandler['y_test']

    model = lstm(config)

    return model, train_input, train_output, valid_input, valid_output, loss_weights, loss_weights_valid
    # model.fit(train_input,
    #              train_output,
    #              batch_size=config.batch_size,
    #              epochs=config.epochs,
    #              verbose=2,
    #              validation_data=(valid_input, valid_output),
    #              shuffle=True,
    #              callbacks=None,
    #              sample_weight=loss_weights)
    # y_pred_test = model.predict([DataHandler['X_test'], DataHandler['X_test_biofeat']], training=False)
    # y_true_test = DataHandler['y_test']
    # if len(y_true_test.shape) == 1:
    #     y_true_test = np.expand_dims(y_true_test, axis=1)
    #
    # spearman = []
    # for ind in range(model.outputs):
    #     spearman.append(sp.stats.spearmanr(y_true_test[:, ind], y_pred_test[:, ind])[0])
    # avg = sum(spearman) / len(spearman)
    # print(f'test spearman: {spearman} mean: {avg}')
    #
    # a=0
    #
    #
    # exit(1)
    # a=0