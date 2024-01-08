import numpy as np
from keras.models import *
from keras.optimizers import *
from sklearn.metrics import mean_squared_error, r2_score
import scipy as sp
from scripts import models_util
from scripts import costume_models_util



def scheduler(epoch, lr):
    learning_rate = 0.01  # initial learning rate
    decay_rate = 0.03
    lrate = learning_rate * np.exp(-decay_rate*epoch)
    return lrate


from sklearn.preprocessing import MinMaxScaler, StandardScaler
def debug(config):
    config.em_drop = 0.4
    config.fc_drop = 0.4
    config.batch_size = 110
    config.epochs = 35
    config.em_dim = 68
    config.fc_num_hidden_layers = 1
    config.fc_num_units = 200
    config.fc_activation = 'elu'
    config.optimizer = SGD
    config.cost_function = 'mse'
    config.last_activation = 'linear'
    config.initializer = 'lecun_uniform'
    config.input_scale = None
    config.output_scale = None
    config.rnn_drop = 0.2
    config.rnn_rec_drop = 0.2
    config.rnn_units = 60

    return config



def train_model(config=None, DataHandler=None):

    # config = debug(config)
    if config.model_type == 'model1':
        model, train_input, y_train, val_input, y_val, callback_list, config, loss_weights, loss_weights_val = models_util.lstm_model(config, DataHandler)

    elif config.model_type == 'model2':
        model, train_input, y_train, val_input, y_val, callback_list, config = models_util.fully_connected_model(config, DataHandler,)
    elif config.model_type == 'model3':
        model, train_input, y_train, val_input, y_val, loss_weights, loss_weights_valid = costume_models_util.train(config, DataHandler)
        callback_list = []

    print('Start training')
    model.fit(train_input,
              y_train,
              batch_size=config.batch_size,
              epochs=config.epochs,
              verbose=2,
              # validation_data=(val_input, y_val, loss_weights_val),\
              validation_data=(val_input, y_val, loss_weights_valid),
              shuffle=True,
              callbacks=callback_list,
              sample_weight=loss_weights
              )

    # y_val_pred = model.predict(val_input)
    # sperman = sp.stats.weightedtau(y_val, y_val_pred, weigher=lambda x: 1)[0]
    # sperman_weigted = sp.stats.weightedtau(y_val, y_val_pred, weigher=lambda x: loss_weights_val[x])[0]
    # print(f'spearman: {sperman}')
    # print(f'sperman_weigted: {sperman_weigted}')
    # exit(1)



    if config.simulation_type == 'train':
        spearman = test_model(config, DataHandler, model)
        return
    if config.simulation_type == 'param_search':
        # return model, config.output_scale, config.input_scale
        return model.best_val_loss

def test_model(config, DataHandler, model):
    # prepare data for testing
    y_test = DataHandler['y_test']
    X_test_biofeat = DataHandler['X_test_biofeat']


    if config.output_scale is not None:
        y_test = config.output_scale.transform(y_test)
    if config.input_scale is not None:
        X_test_biofeat = config.input_scale.transform(X_test_biofeat)

    if config.has_biofeatures:
        test_input = [DataHandler['X_test'], X_test_biofeat]
    else:
        test_input = [DataHandler['X_test']]

    if config.enzyme == 'multi_task':
        data_types = ['wt', 'esp', 'hf']
        for ind, type in enumerate(data_types):
            test_indexes = np.logical_not(np.isnan(y_test[:, ind]))
            test_input_enzyme = test_input.copy()
            for i in range(len(test_input_enzyme)):
                test_input_enzyme[i] = test_input[i][test_indexes]
                # test_input_enzyme[1] = test_input[1][test_indexes]
            y_test_true_enzyme = y_test[:, ind][test_indexes]

            if config.multi_data == 'parallel':
                y_test_pred_enzyme = model.predict(test_input_enzyme)
                mse = mean_squared_error(y_test_true_enzyme, y_test_pred_enzyme[:, ind])
                spearman = sp.stats.spearmanr(y_test_true_enzyme, y_test_pred_enzyme[:, ind])[0]
            elif config.multi_data == 'serial':
                x_biofeat_test_enzyme = test_input_enzyme[1]
                ohe = np.eye(3)[np.zeros((x_biofeat_test_enzyme.shape[0]), dtype=np.int) + ind]
                x_biofeat_test_enzyme = np.append(x_biofeat_test_enzyme, ohe, axis=-1)
                test_input_enzyme[1] = x_biofeat_test_enzyme
                y_test_pred_enzyme = model.predict(test_input_enzyme)
                mse = mean_squared_error(y_test_true_enzyme, y_test_pred_enzyme)
                spearman = sp.stats.spearmanr(y_test_true_enzyme, y_test_pred_enzyme)[0]
            print('type: {}, spearman: {}, mse: {}'.format(type, spearman, mse))
    else:
        test_indexes = np.logical_not(np.isnan(y_test))
        test_input[0] = test_input[0][test_indexes]
        test_input[1] = test_input[1][test_indexes]
        y_test = y_test[test_indexes]

        y_test_pred = model.predict(test_input)
        evaluation = mean_squared_error(y_test, y_test_pred)
        spearman = sp.stats.spearmanr(y_test, y_test_pred)[0]
        print('spearman: {}'.format(spearman))
        print('evaluation: {}'.format(evaluation))

    return spearman

