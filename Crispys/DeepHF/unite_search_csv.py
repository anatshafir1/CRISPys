import argparse
import pandas as pd
import numpy as np
import os
import csv
import GPyOpt
import scripts.HyperParameterSearching as HPS
from scripts import configurations as cfg


def unite_csv(config):

    dir_path = f'HPS/{config.data_source}/{config.model_type}/'
    if config.has_biofeatures:
        dir_path += 'bio/'
    else:
        dir_path += 'no_bio/'

    dir_path = dir_path + f'{config.enzyme}/'


    if config.enzyme == 'multi_task':
        dir_path += f'{config.multi_data}/'

    main_eval_path = dir_path + 'HPS.csv'
    main_evals = pd.read_csv(main_eval_path, index_col=0, delimiter="\t")
    Y = np.array([[x] for x in main_evals["Y"]])
    X = np.array(main_evals.filter(regex="var*"))

    ind = 0
    save_path = dir_path + f'HPS_{ind}.csv' #"logs/HyperParamSearch_{}_{}.csv".format(config.data_type, ind)
    while os.path.exists(save_path):
        new_evals = pd.read_csv(save_path, index_col=0, delimiter="\t")
        Y_new = np.array([[x] for x in new_evals["Y"]])
        X_new = np.array(new_evals.filter(regex="var*"))

        for val, new_params in zip(Y_new, X_new):
            is_new = True
            for old_params in X:
                if np.array_equal(new_params, old_params):
                    is_new = False
                    break
            if is_new:
                X = np.vstack([X, new_params])
                Y = np.vstack([Y, val])

        ind += 1
        save_path = dir_path + f'HPS_{ind}.csv'


    iterations = np.array(range(1, Y.shape[0] + 1))[:, None]
    results = np.hstack((iterations, Y, X))
    header = ['Iteration', 'Y'] + ['var_' + str(k) for k in range(1, X.shape[1] + 1)]

    data = [header] + results.tolist()
    with open(main_eval_path, 'w') as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')
        writer.writerows(data)


def sort_csv(config):
    dir_path = f'HPS/{config.data_source}/{config.model_type}/'
    if config.has_biofeatures:
        dir_path += 'bio/'
    else:
        dir_path += 'no_bio/'

    dir_path = dir_path + f'{config.enzyme}/'


    if config.enzyme == 'multi_task':
        dir_path += f'{config.multi_data}/'

    main_eval_path = dir_path + 'HPS.csv'

    main_evals = pd.read_csv(main_eval_path, index_col=0, delimiter="\t")
    main_evals = main_evals.sort_values(by=['Y'])

    # if config.model_type == 'model1':
    #     main_evals = main_evals.mul([-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    #     main_evals = main_evals.add([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    # elif config.model_type == 'model2':
    #     main_evals = main_evals.mul([-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    #     main_evals = main_evals.add([1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    # elif config.model_type == 'model3':
    #     main_evals = main_evals.mul([-1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    #     main_evals = main_evals.add([1,0,0,0,0,0,0,0,0,0,0,0,0,0])



    bounds = HPS.get_bounds(has_logger=False, model_type=config.model_type)
    titles = ['loss']
    for bound in bounds:
        titles.append(bound['name'])
    dictionary = {}
    for old, new in zip(main_evals.columns, titles):
        dictionary[old] = new
    main_evals.rename(columns = dictionary, inplace=True)
    main_evals.to_csv(dir_path + 'out.csv')

config = cfg.get_parser()
unite_csv(config)
sort_csv(config)