"""DeepHF imports and model loading"""
import pickle
from keras.models import load_model

from Crispys import Distance_matrix_and_UPGMA, globals


def load_deephf():
    """
    Function to load DeepHF algorithm model
    """
    with open("DeepHF/models/parallel/config.pkl", "rb") as fp:
        config = pickle.load(fp)
    model = load_model("DeepHF/models/parallel/model")
    globals.deephf_loaded_model = model
    globals.deephf_config = config
    return Distance_matrix_and_UPGMA.deephf
