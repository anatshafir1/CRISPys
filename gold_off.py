"""gold_off algorithm module"""
import random
import xgboost as xgb
import pandas as pd
import numpy as np
import itertools
from multiprocessing import Pool
from functools import partial

SEED = 10
random.seed(SEED)


######################################################################
def create_nucleotides_to_position_mapping():
    """
    Return the nucleotides to position mapping
    """
    # matrix positions for ('A','A'), ('A','C'),...
    # tuples of ('A','A'), ('A','C'),...
    nucleotides_product = list(itertools.product(*(['ACGT'] * 2)))
    # tuples of (0,0), (0,1), ...
    position_product = [(int(x[0]), int(x[1]))
                        for x in list(itertools.product(*(['0123'] * 2)))]
    nucleotides_to_position_mapping = dict(
        zip(nucleotides_product, position_product))

    # tuples of ('N','A'), ('N','C'),...
    n_mapping_nucleotides_list = [('N', char) for char in ['A', 'C', 'G', 'T']]
    # list of tuples positions corresponding to ('A','A'), ('C','C'), ...
    n_mapping_position_list = [nucleotides_to_position_mapping[(char, char)]
                               for char in ['A', 'C', 'G', 'T']]

    nucleotides_to_position_mapping.update(
        dict(zip(n_mapping_nucleotides_list, n_mapping_position_list)))

    # tuples of ('A','N'), ('C','N'),...
    n_mapping_nucleotides_list = [(char, 'N') for char in ['A', 'C', 'G', 'T']]
    # list of tuples positions corresponding to ('A','A'), ('C','C'), ...
    n_mapping_position_list = [nucleotides_to_position_mapping[(char, char)]
                               for char in ['A', 'C', 'G', 'T']]
    nucleotides_to_position_mapping.update(
        dict(zip(n_mapping_nucleotides_list, n_mapping_position_list)))

    return nucleotides_to_position_mapping


def hamming_distance(s1, s2):
    """

    :param s1:
    :param s2:
    :return:
    :rtype:
    """
    if len(s1) != len(s2):
        raise ValueError("Strand lengths are not equal!")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def pair_encoding(i, target, off_target, distance, nucleotides_to_position_mapping,
                  final_result, include_distance_feature):
    """

    :param i:
    :param target:
    :param off_target:
    :param distance:
    :param nucleotides_to_position_mapping:
    :param final_result:
    :param include_distance_feature:
    """
    intersection_matrices = np.zeros((23, 4, 4), dtype=np.int8)
    for j in range(23):
        matrix_positions = nucleotides_to_position_mapping[(
            target[j], off_target[j])]
        intersection_matrices[j, matrix_positions[0],
                              matrix_positions[1]] = 1

    if include_distance_feature:
        final_result[i, :-1] = intersection_matrices.flatten()
        final_result[i, -1] = distance if distance is not None else \
            hamming_distance(target[:-3] + target[-2:], off_target[:-3] + off_target[-2:])
    else:
        final_result[i, :] = intersection_matrices.flatten()


def build_sequences_encoding(dataset_df, nucleotides_to_position_mapping,
                             include_distance_feature):
    """

    :param dataset_df:
    :param nucleotides_to_position_mapping:
    :param include_distance_feature:
    :return:
    :rtype:
    """
    dataset_df.reset_index(drop=True, inplace=True)
    final_result = np.zeros((len(dataset_df), (23 * 16) + 1),
                            dtype=np.int8) if include_distance_feature else \
        np.zeros((len(dataset_df), 23 * 16), dtype=np.int8)

    if "distance" in dataset_df:
        for (i, (target, off_target, distance)) in enumerate(zip(dataset_df["target"],
                                                                 dataset_df["offtarget_sequence"],
                                                                 dataset_df["distance"])):
            pair_encoding(i, target, off_target, distance, nucleotides_to_position_mapping,
                          final_result, include_distance_feature)
    else:
        for (i, (target, off_target)) in enumerate(zip(dataset_df["target"],
                                                       dataset_df["offtarget_sequence"])):
            pair_encoding(i, target, off_target, None, nucleotides_to_position_mapping,
                          final_result, include_distance_feature)
    return final_result


def build_sequence_features(dataset_df, nucleotides_to_position_mapping,
                            include_distance_feature=False, n_cores=10):
    """
    Build sequence features using the nucleotides to position mapping
    """
    dataset_df.reset_index(drop=True, inplace=True)

    # convert dataset_df["target"] -3 position to 'N'
    dataset_df['target'] = dataset_df['target'].apply(lambda s: s[:-3] + 'N' + s[-2:])
    if n_cores == 1:
        # Avoid multiprocessing when running with a single core
        final_result = build_sequences_encoding(dataset_df,
                                                nucleotides_to_position_mapping=nucleotides_to_position_mapping,
                                                include_distance_feature=include_distance_feature)
        return final_result
    dataset_df_splits = np.array_split(dataset_df, n_cores)
    with Pool(n_cores) as pool:
        final_result = np.concatenate(
            pool.map(partial(build_sequences_encoding,
                             include_distance_feature=include_distance_feature,
                             nucleotides_to_position_mapping=nucleotides_to_position_mapping),
                     dataset_df_splits), axis=0)
    return final_result


##########################################################################


def predict(sg_rnas, off_targets, model_path,
            distances=None, include_distance_feature=False, n_process=10,
            model_type="classification"):
    """

    :param sg_rnas:
    :param off_targets:
    :param model_path:
    :param distances:
    :param include_distance_feature:
    :param n_process:
    :param model_type:
    :return:
    :rtype:
    """
    # load the model
    if model_type == "classification":
        model = xgb.XGBClassifier()
    elif model_type == "regression":
        model = xgb.XGBRegressor()
    else:
        raise ValueError("Invalid model_type. should 'classification' or 'regression'")
    model.load_model(model_path)

    # create features
    dataset_df = pd.DataFrame()
    dataset_df["target"] = pd.Series(sg_rnas)
    dataset_df["offtarget_sequence"] = pd.Series(off_targets)
    if distances is not None:
        dataset_df["distance"] = pd.Series(distances)
    nucleotides_to_position_mapping = create_nucleotides_to_position_mapping()
    features = build_sequence_features(
        dataset_df, nucleotides_to_position_mapping,
        include_distance_feature=include_distance_feature,
        n_cores=n_process)

    predictions = model.predict_proba(features)[:, 1] if model_type == "classification" \
        else model.predict(features)
    return predictions
