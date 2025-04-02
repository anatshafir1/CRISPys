"""DeepHF predictions"""
import pandas as pd
import numpy as np
from . import preprocess
import scipy as sp


def get_predictions(model, config, sequences, decoded=False):
    if not decoded:
        sequences = prepare_arr(sequences)
    predictions = model.predict(sequences)
    # Scale outputs
    scaled_predictions = config.output_scale.inverse_transform(predictions)
    return scaled_predictions[:, 0]


def prepare_arr(sequences):
    decoded_seq = [decode_seq(seq) for seq in sequences]
    decoded_seq = np.concatenate(decoded_seq, axis=0)
    return decoded_seq


def decode_seq(seq):
    mer_array = np.array([0], dtype=np.uint8)
    for char in seq:
        num = preprocess.char2int(char)
        mer_array = np.append(mer_array, np.array([num], dtype=np.uint8), axis=0)
    mer_array = np.expand_dims(mer_array, 0)
    return mer_array


def debug(model_dir, inputs, true_outputs):
    df = pd.read_csv('data/suplementry2.csv')
    sequences = list(df['21mer'])[:1000]
    pred_outputs = get_predictions(model_dir, inputs, decoded=True)

    true_outputs = true_outputs[:, 0]
    test_indexes = np.logical_not(np.isnan(true_outputs))
    true_outputs = true_outputs[test_indexes]
    pred_outputs = pred_outputs[test_indexes]

    spearman = sp.stats.spearmanr(true_outputs, pred_outputs)[0]

    print(spearman)
    a = 0
