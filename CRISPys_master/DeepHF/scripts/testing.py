import os
from keras.models import load_model
from scripts.Cross_validation import load_data_kf
import pickle
import scipy as sp
from sklearn.metrics import mean_squared_error, r2_score
from scripts import data_handler as dh
import numpy as np
from scripts import costume_models_util
from scripts import data_handler as dh


def test_model(config):
    print('Testing model')
    if not os.path.exists('models/' + config.model_path):
        print('model path {} dosent exists'.format(config.model_path))
        exit(1)
    if config.model_type == 'model3':
        DataHandler = dh.get_data(config)
        config.save_model = False
        model, train_input, y_train, val_input, y_val, loss_weights, loss_weights_valid = costume_models_util.train(config, DataHandler)
        with open('models/' + config.model_path, "rb") as fp:
            saved_weights = pickle.load(fp)
        model.fit(train_input,
                  y_train,
                  batch_size=config.batch_size,
                  epochs=1,  # config.epochs,
                  verbose=2,
                  # validation_data=(val_input, y_val, loss_weights_val),\
                  validation_data=(val_input, y_val, loss_weights_valid),
                  shuffle=True,
                  callbacks=[],
                  sample_weight=loss_weights
                  )
        model.set_weights(saved_weights)
    else:
        model = load_model('models/' +  config.model_path)
    print('Model has been loaded')
    # load data
    if config.enzyme == 'leenay':
        with open('data/preprocessed_data/leenay_seq.pkl', 'rb') as handle:
            leenay_seq = pickle.load(handle)
            X, y = leenay_seq.get_data()
            y = y['wt']
    else:
        DataHandler = dh.get_data(config)
        # with open('data/' + config.data_type + '_seq_data_array.pkl', 'rb') as handle:
        #     X, X_biofeat, y = pickle.load(handle)

    # test data
    test_input = [DataHandler['X_test'], DataHandler['X_test_biofeat']]
    test_indexes = np.logical_not(np.isnan(DataHandler['y_test']))
    test_input[0] = test_input[0][test_indexes]
    test_input[1] = test_input[1][test_indexes]
    y_test = DataHandler['y_test'][test_indexes]

    # Valid data
    valid_input = [DataHandler['X_valid'], DataHandler['X_valid_biofeat']]
    valid_indexes = np.logical_not(np.isnan(DataHandler['y_valid']))
    valid_input[0] = valid_input[0][valid_indexes]
    valid_input[1] = valid_input[1][valid_indexes]
    y_valid = DataHandler['y_valid'][valid_indexes]


    y_pred = model.predict(test_input)
    if config.enzyme == 'multi_task':
        for ind in range(3):
            enzim_pred = y_pred[:,ind]
            spearmanr = sp.stats.spearmanr(DataHandler['y_test'], enzim_pred)[0]
            print('ind: {}, spearman: {}'.format(ind, spearmanr))
    else:
        spearmanr = sp.stats.spearmanr(DataHandler['y_test'], y_pred)
        print('spearman: {}'.format(spearmanr))



    # data_list = load_data_kf(X, X_biofeat, y)
    # training_r2_scores = []
    # testing_r2_scores = []
    # spearmanrs_test = []
    # spearmanrs_train = []
    # loss_train = []
    # loss_test = []
    #
    # for index, data in enumerate(data_list):
    #     X_train, X_train_biofeat, X_test, X_test_biofeat, y_train, y_test = data[0], data[1], data[2], data[3], data[4], \
    #                                                               data[5]
    #     y_train_pred = model.predict([X_train, X_train_biofeat])
    #     training_score = r2_score(y_train, y_train_pred)
    #     mse_train = mean_squared_error(y_train, y_train_pred)
    #     spearmanr_train = sp.stats.spearmanr(y_train, y_train_pred)[0]
    #     training_r2_scores.append(training_score)
    #     spearmanrs_train.append(spearmanr_train)
    #     loss_train.append(mse_train)
    #
    #     y_test_pred = model.predict([X_test, X_test_biofeat])
    #     testing_score = r2_score(y_test, y_test_pred)
    #     mse_test = mean_squared_error(y_test, y_test_pred)
    #     spearmanr_test = sp.stats.spearmanr(y_test, y_test_pred)[0]
    #     testing_r2_scores.append(testing_score)
    #     spearmanrs_test.append(spearmanr_test)
    #     loss_test.append(mse_test)
    #     print('train: {}, test: {}'.format(spearmanr_train, spearmanr_test))

    a=0

def test_means(config):
    if not os.path.exists('models/' + config.model_path):
        print('model path {} dosent exists'.format(config.model_path))
        exit(1)

    DataHandler = dh.get_data(config)
    config.save_model = False
    model, train_input, y_train, val_input, y_val, loss_weights, loss_weights_valid = costume_models_util.train(config,
                                                                                                                DataHandler)
    model.fit(train_input,
              y_train,
              batch_size=config.batch_size,
              epochs=1,  # config.epochs,
              verbose=2,
              # validation_data=(val_input, y_val, loss_weights_val),\
              validation_data=(val_input, y_val, loss_weights_valid),
              shuffle=True,
              callbacks=[],
              sample_weight=loss_weights
              )
    saved_weights_names = next(os.walk('models/' + config.model_path))[2]
    models_preds_arr = []
    test_input = [DataHandler['X_test'], DataHandler['X_test_biofeat']]

    for ind, saved_weights_name in enumerate(saved_weights_names):
        with open('models/' + config.model_path + f'/{saved_weights_name}', "rb") as fp:
            saved_weights = pickle.load(fp)
        model.set_weights(saved_weights)
        y_pred = model.predict(test_input).numpy()
        models_preds_arr.append(y_pred)

    finall_pred = np.zeros((y_pred.shape[0], 1))
    for pred in models_preds_arr:
        finall_pred += pred
    finall_pred /= len(models_preds_arr)
    spearmanr = sp.stats.spearmanr(DataHandler['y_test'], finall_pred)
    print('spearman: {}'.format(spearmanr))

    pred_sum = 0
    for ind, pred in enumerate(models_preds_arr):
        pred_sum += pred
        finall_pred = pred_sum / (ind+1)
        spearmanr = sp.stats.spearmanr(DataHandler['y_test'], finall_pred)
        print('models: {}, spearman: {}'.format(ind+1, spearmanr))

    a=0



