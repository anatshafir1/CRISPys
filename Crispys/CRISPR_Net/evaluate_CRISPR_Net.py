# -*- coding: utf-8 -*-
# @Time    : 7/4/2020 8:17 PM
# @Author  : Jason Lin
# @File    : evaluate_CRISPR_Net.py
# @Software: PyCharm
import os
import pickle as pkl
import encode_data
from keras.models import Model, model_from_json
import keras
import numpy as np
from keras.layers import Input, Dense, Reshape, Flatten, LSTM, Bidirectional
from keras import layers
from keras.layers import Input
from sklearn import metrics
from numpy.random import seed
seed(1)
from tensorflow import set_random_seed
set_random_seed(1)


def conv2d_bn(x,
              filters,
              kernel_size,
              strides=1,
              padding='same',
              activation='relu',
              use_bias=True,
              name=None, trainable=True):
    x = layers.Conv2D(filters,
                      kernel_size,
                      strides=strides,
                      padding=padding,
                      use_bias=use_bias,
                      name=name, trainable=trainable)(x)

    # x = layers.BatchNormalization(axis=-1,scale=True)(x)
    if activation is not None:
        ac_name = None if name is None else name + '_ac'
        x = layers.Activation(activation, name=ac_name)(x)
    return x

def CRISPR_Net_model():
    inputs = Input(shape=(1, 24, 7), name='main_input')
    branch_0 = conv2d_bn(inputs, 10, (1, 1))
    branch_1 = conv2d_bn(inputs, 10, (1, 2))
    branch_2 = conv2d_bn(inputs, 10, (1, 3))
    branch_3 = conv2d_bn(inputs, 10, (1, 5))
    branches = [inputs, branch_0, branch_1, branch_2, branch_3]
    # branches = [branch_0, branch_1, branch_2, branch_3]
    mixed = layers.Concatenate(axis=-1)(branches)
    mixed = Reshape((24, 47))(mixed)
    blstm_out = Bidirectional(LSTM(15, return_sequences=True, input_shape=(24, 47), name="LSTM_out"))(mixed)
    # inputs_rs = Reshape((24, 7))(inputs)
    # blstm_out = layers.Concatenate(axis=-1)([mixed, blstm_out])
    blstm_out = Flatten()(blstm_out)
    x = Dense(80, activation='relu')(blstm_out)
    x = Dense(20, activation='relu')(x)
    x = keras.layers.Dropout(rate=0.35)(x)
    prediction = Dense(1, activation='sigmoid', name='main_output')(x)
    model = Model(inputs, prediction)
    print(model.summary())
    return model


def load_elevation_data():
    if not os.path.exists("./tmp_code"):
        os.mkdir("./tmp_code")
    train_code_file = "./tmp_code/elevation_training_code.pkl"
    if os.path.exists(train_code_file):
        print("Loading Dataset I-1, II-1, II-2, and II-4 ...")
        X, y = pkl.load(open(train_code_file, "rb"))
    else:
        dataset_II1_code, dataset_II1_labels = encode_data.load_elevation_CD33_dataset()
        dataset_II2_code, dataset_II2_labels = encode_data.load_elevation_hmg_dataset()
        dataset_II4_code, dataset_II4_labels = encode_data.load_elevation_guideseq_data()
        X =  np.concatenate([dataset_II4_code, dataset_II2_code, dataset_II1_code], axis=0)
        y  = np.concatenate([dataset_II4_labels, dataset_II2_labels, dataset_II1_labels], axis=0)
        pkl.dump([X, y], open("./tmp_code/elevation_training_code.pkl", "wb"))
    X = X.reshape((len(X), 1, 24, 7))
    return X, y

def load_CIRCLESeq_data():
    if not os.path.exists("./tmp_code"):
        os.mkdir("./tmp_code")
    train_code_file = "./tmp_code/CIRCLE_seq_code.pkl"
    if os.path.exists(train_code_file):
        print("Loading CIRCLE-Seq data (Dataset I-1) ...")
        X, y = pkl.load(open(train_code_file, "rb"))
    else:
        X, y = encode_data.load_CIRCLE_data()
        pkl.dump([X, y], open("./tmp_code/CIRCLE_seq_code.pkl", "wb"))
    X = X.reshape((len(X), 1, 24, 7))
    return X, y

def load_SITESeq_data():
    if not os.path.exists("./tmp_code"):
        os.mkdir("./tmp_code")
    train_code_file = "./tmp_code/SITE_seq_code.pkl"
    if os.path.exists(train_code_file):
        print("Loading CIRCLE-Seq data (Dataset I-1) ...")
        X, y = pkl.load(open(train_code_file, "rb"))
    else:
        X, y = encode_data.load_siteseq_data()
        pkl.dump([X, y], open("./tmp_code/SITE_seq_code.pkl", "wb"))
    X = X.reshape((len(X), 1, 24, 7))
    return X, y

def load_testing_data(ds_type = 'listgarten'):
    if not os.path.exists("./tmp_code"):
        os.mkdir("./tmp_code")
    test_code_file1 = "./tmp_code/testing_code_listgarten.pkl"
    test_code_file2 = "./tmp_code/testing_code_Kleinstiver.pkl"
    test_code_file3 = "./tmp_code/testing_code_listgarten_indel.pkl"
    if os.path.exists(test_code_file1):
        X_test1, y_test1 = pkl.load(open(test_code_file1, "rb"))
    else:
        X_test1, y_test1 = encode_data.load_22sgRNA_data()
        pkl.dump([X_test1, y_test1], open("./tmp_code/testing_code_listgarten.pkl", "wb"))

    if os.path.exists(test_code_file2):
        X_test2, y_test2 = pkl.load(open(test_code_file2, "rb"))
    else:
        X_test2, y_test2 = encode_data.load_Kleinstiver_data()
        pkl.dump([X_test2, y_test2], open("./tmp_code/testing_code_Kleinstiver.pkl", "wb"))

    if os.path.exists(test_code_file3):
        X_test3, y_test3 = pkl.load(open(test_code_file3, "rb"))
    else:
        X_test3, y_test3 = encode_data.load_Listgarten_indel_dataset()
        pkl.dump([X_test3, y_test3], open("./tmp_code/testing_code_listgarten_indel.pkl", "wb"))

    if ds_type == "listgarten":
        print("Loading Listgarten testing dataset (Dataset II-6) ...")
        X_test1 = X_test1.reshape((len(X_test1), 1, 24, 7))
        print("Finished!", "Data_size:", X_test1.shape, "Num of Positive sample:", len(y_test1[y_test1>0]))
        return X_test1, y_test1 
    elif ds_type == "kleinstiver":
        print("Loading Kleinstiver testing dataset (Dataset II-5) ...")
        X_test2 = X_test2.reshape((len(X_test2), 1, 24, 7))
        print("Finished!", "Data_size:", X_test2.shape, "Num of Positive sample:", len(y_test2[y_test2>0]))
        return X_test2, y_test2
    elif ds_type == "listgarten_indel":
        print("Loading Listgarten indel testing dataset (Dataset I-2) ...")
        X_test3 = X_test3.reshape((len(X_test3), 1, 24, 7))
        print("Finished!", "Data_size:", X_test3.shape, "Num of Positive sample:", len(y_test3[y_test3>0]))
        return X_test3, y_test3
    else:
        print("Dataset not found!")
        return 0

def CRISPR_Net_training(X_train, y_train, EPOCHS):
    #  X_train, y_train = load_elevation_data()
    #  X_train.shape (len, 1, 24, 7)
    #  y_train.shape {0, 1}^len
    model = CRISPR_Net_model()
    adam_opt = keras.optimizers.Adam(lr=0.0001)
    model.compile(loss="binary_crossentropy", optimizer=adam_opt)
    model.fit(X_train, y_train, batch_size=10000, epochs=EPOCHS, shuffle=True)
    # SAVE Model
    model_jason = model.to_json()
    model_path = "scoring_models"
    if os.path.isdir(model_path):
        pass
    else:
        os.mkdir(model_path)
    with open(model_path + "/CRISPR_Net_structure_user.json", "w") as jason_file:
        jason_file.write(model_jason)
    model.save_weights(model_path + "/CRISPR_Net_weights_user.h5")
    print("Saved model to disk!")

def CRISPR_Net_testing(weights_file, test_ds='listgarten'):
    X_test, y_test = load_testing_data(ds_type=test_ds)
    # print(X_test[0])
    print(len(y_test[y_test > 0]), len(y_test))
    json_file = open("scoring_models/CRISPR_Net_structure.json", 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    print(loaded_model.summary())
    loaded_model.load_weights(weights_file)
    print("Loaded model from disk!")
    y_pred = loaded_model.predict(X_test).flatten()
    # AUROC and AUPRC
    roc_auc = metrics.roc_auc_score(y_test, y_pred)
    prc, rec, _ = metrics.precision_recall_curve(y_test, y_pred)
    prc[[rec == 0]] = 1.0
    prc_auc = metrics.auc(rec, prc)
    print("Testing results on " + test_ds + " dataset: ")
    print("AUROC: ", roc_auc, "AUPRC: ", prc_auc)
    return roc_auc, prc_auc


# ----------------------------Mismatch----------------------------
# Trained on CIRCLE-Seq and Elevation dataset
weight_file = "scoring_models/CRISPR_Net_CIRCLE_elevation_weights.h5"

# Trained on Elevation dataset
weight_file = "scoring_models/CRISPR_Net_elevation_weights.h5"

# Trained on SITE-Seq dataset 
weight_file = "scoring_models/CRISPR_Net_SITE_weights.h5"

# Trained on CIRCLE-Seq, Elevation and SITE-Seq dataset
weight_file = "scoring_models/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5"
# test on Dataset II/5
roc_auc, prc_auc = CRISPR_Net_testing(weight_file, test_ds='kleinstiver')
# test on Dataset II/6
roc_auc, prc_auc = CRISPR_Net_testing(weight_file, test_ds='listgarten')


# ---------------------------Indel & Mismatch--------------------------
# Trained on CIRCLE-Seq and Elevation indel dataset
weight_file1 = "./scoring_models/CRISPR_Net_CIRCLE_elevation_weights.h5"
# Trained on CIRCLE-Seq indel dataset
weight_file2 = "./scoring_models/CRISPR_Net_CIRCLE_indel_weights.h5"

# test on Dataset I/2
roc_auc, prc_auc = CRISPR_Net_testing(weight_file1, test_ds='listgarten_indel')
roc_auc, prc_auc = CRISPR_Net_testing(weight_file2, test_ds='listgarten_indel')


