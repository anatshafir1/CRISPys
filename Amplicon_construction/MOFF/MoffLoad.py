"""MOFF imports and model loading"""
import os
from keras import models
from json import loads
import keras
import tensorflow.python

PATH = os.path.dirname(os.path.realpath(__file__))


def load_moff():
    """
    Function to load MOFF algorithm model
    """
    moff_mtx1 = loads(open(f"{PATH}/StaticFiles/M1_matrix_dic_D9").read())
    moff_mtx2 = loads(open(f"{PATH}/StaticFiles/M2_matrix_smooth_MLE").read())
    moff_loaded_model = models.load_model(rf"{PATH}/StaticFiles/GOP_model_3.h5")
    return moff_mtx1, moff_mtx2, moff_loaded_model


mtx1, mtx2, model = load_moff()
