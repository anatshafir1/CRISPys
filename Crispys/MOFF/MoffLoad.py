"""MOFF imports and model loading"""

from keras import models
from json import loads
from Crispys import Distance_matrix_and_UPGMA, globals


def load_moff():
    """
    Function to load MOFF algorithm model
    """
    globals.moff_mtx1 = loads(open(f"{globals.CODE_PATH}/MOFF/StaticFiles/M1_matrix_dic_D9").read())
    globals.moff_mtx2 = loads(open(f"{globals.CODE_PATH}/MOFF/StaticFiles/M2_matrix_smooth_MLE").read())
    globals.moff_loaded_model = models.load_model(rf"{globals.CODE_PATH}/MOFF/StaticFiles/GOP_model_3.h5")
    return Distance_matrix_and_UPGMA.moff
