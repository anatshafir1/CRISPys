import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from scipy import stats
from scipy import special
import elevation_metrics

def normalise_LFC(lfc_values):
    lfc_values = (lfc_values - min(lfc_values))/(max(lfc_values)-min(lfc_values))
    lfc_values = special.boxcox1p(lfc_values, 0.2)
    return lfc_values

def cal_correlation(pred, actual):
    spear_cor, spear_p = stats.spearmanr(pred, actual)
    weighted_spear_cor = elevation_metrics.spearman_weighted(pred, actual, w=actual)
    return spear_cor, weighted_spear_cor

avana_df = pd.read_csv("../data/Aggregate_dataset/Avana_aggregate_data_codeocean.csv")
avana_label = np.array(avana_df['lentiGuide Avg'])
avana_data = avana_df.values
avana_feature = avana_data[:, 1:24]
avana_label = normalise_LFC(avana_label)

gecko_df = pd.read_csv("../data/Aggregate_dataset/GeCKO_aggregate_data_codeocean.csv")
gecko_label_rep1 = normalise_LFC(np.array(gecko_df['A375RepB_LFC']))
gecko_label_rep2 = normalise_LFC(np.array(gecko_df['A375RepC_LFC']))
gecko_label_rep_avg = normalise_LFC(np.array(gecko_df['A375_LFC_Avg']))
gecko_labels = {"Replicate_1" : gecko_label_rep1,"Replicate_2" : gecko_label_rep2, "Average" : gecko_label_rep_avg}
gecko_data = gecko_df.values
gecko_feature = gecko_data[:, 1:24]

reg = ExtraTreesRegressor(random_state=66, n_estimators=300, max_depth=5, max_features="log2")
reg.fit(avana_feature, avana_label)
# joblib.dump(reg, "./CRISPR-Net-Aggregator.m")

pred_gecko = reg.predict(gecko_feature)
elevation_pred_gecko = gecko_df['Elevation_agg_val']

for rep in gecko_labels.keys():
    print(rep)
    res_dict = {}
    scorr, wscorr = cal_correlation(pred_gecko, gecko_labels[rep])
    res_dict['CRSIPR-Net-Aggregate'] = {'spearman(scipy.stats)':scorr,'weighted spearman(Elevation.spearman)':wscorr}
    scorr, wscorr = cal_correlation(elevation_pred_gecko, gecko_labels[rep])
    res_dict['Elevation-aggregate'] = {'spearman(scipy.stats)':scorr,'weighted spearman(Elevation.spearman)':wscorr}
    print('Elevation-aggregate:\n', res_dict['Elevation-aggregate'])
    print('CRSIPR-Net-Aggregate:\n', res_dict['CRSIPR-Net-Aggregate'])
    print('-'*50)