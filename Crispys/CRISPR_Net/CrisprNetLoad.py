"""CRISPR_net imports and model loading"""
import warnings
import os
from tensorflow.keras.models import model_from_json
# set tensorflow to use 1 core
from tensorflow.config import threading
from Crispys import Distance_matrix_and_UPGMA, globals

threading.set_inter_op_parallelism_threads(1)

warnings.filterwarnings('ignore')
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = ""


def load_crispr_net():
    # load the model and make it available globally
    json_file = open(f"{globals.CODE_PATH}/CRISPR_Net/scoring_models/CRISPR_Net_structure.json", 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    crisprnet_loaded_model = model_from_json(loaded_model_json)
    crisprnet_loaded_model.load_weights(
        f"{globals.CODE_PATH}/CRISPR_Net/scoring_models/CRISPR_Net_CIRCLE_elevation_SITE_weights.h5")
    globals.crisprnet_loaded_model = crisprnet_loaded_model
    print("Loaded model from disk!")
    return Distance_matrix_and_UPGMA.crisprnet
