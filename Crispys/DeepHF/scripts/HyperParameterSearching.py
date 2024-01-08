import logging
#import GPy
import GPyOpt
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import pandas as pd
from scripts.training_util import *
from scripts import data_handler as dh
import time
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.compose import ColumnTransformer
import os

# string parameters changed to dictionary which keys are integers
fc_activation_dict = {'1': 'relu', '2': 'tanh', '3': 'sigmoid', '4': 'hard_sigmoid', '0': 'elu'}
optimizer_dict =  {'1': SGD, '2': RMSprop, '3': Adagrad, '4': Adadelta, '5': Adam, '6': Adamax, '0': Nadam}
cost_function_dict = {'0': 'mse', '1': 'binary_crossentropy'}
last_activation_dict = {'0': 'sigmoid', '1': 'linear'}
initializer_dict = {'1': 'lecun_uniform', '2': 'normal', '3': 'he_normal', '0': 'he_uniform'}
input_scale_dict = {'0': 'None', '1': 'norm_all', '2': 'stand_all', '3': 'norm_non_gausian', '4': 'stand_gausian', '5': 'selective'}
output_scale_dict = {'0': 'None', '1': 'stand'}


def get_bounds(enzyme=None, has_logger=True, model_type=None):
    # if config:
    #     glob_config = config
    # parameter space

    basic_bounds = [
        # Discrete
        {'name': 'em_drop', 'type': 'discrete', 'domain': (0.1, 0.2, 0.3, 0.4, 0.5)},
        {'name': 'fc_drop', 'type': 'discrete', 'domain': (0.1, 0.2, 0.4, 0.4, 0.5)},
        # Discrete
        {'name': 'batch_size', 'type': 'discrete', 'domain': (80, 90, 100, 110, 120, 140, 160, 200)},
        {'name': 'em_dim', 'type': 'discrete', 'domain': (16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68)},
        {'name': 'fc_num_hidden_layers', 'type': 'discrete', 'domain': (1, 2, 3, 4)},
        {'name': 'fc_num_units', 'type': 'discrete',
         'domain': (120, 140, 160, 180, 200, 220, 240, 260, 280, 300)},
        # Categorical
        {'name': 'fc_activation', 'type': 'categorical', 'domain': tuple(map(int, tuple(fc_activation_dict.keys())))},
        {'name': 'optimizer', 'type': 'categorical', 'domain': tuple(map(int, tuple(optimizer_dict.keys())))},
        {'name': 'last_activation', 'type': 'categorical', 'domain': tuple(map(int, tuple(last_activation_dict.keys())))},
        {'name': 'initializer', 'type': 'categorical', 'domain': tuple(map(int, tuple(initializer_dict.keys())))}

    ]

    advanced_bounds = [
        {'name': 'cost_function', 'type': 'categorical', 'domain': tuple(map(int, tuple(cost_function_dict.keys())))},
        {'name': 'epochs', 'type': 'discrete', 'domain': (35, 40, 45, 50, 55, 60, 65)},
        {'name': 'input_scale', 'type': 'categorical', 'domain': tuple(map(int, tuple(input_scale_dict.keys())))},
        {'name': 'output_scale', 'type': 'categorical', 'domain': tuple(map(int, tuple(output_scale_dict.keys())))}

    ]
    rnn_bounds = [
        # Discrete
        {'name': 'rnn_drop', 'type': 'discrete', 'domain': (0.1, 0.2, 0.4, 0.4, 0.5)},
        {'name': 'rnn_rec_drop', 'type': 'discrete', 'domain': (0.1, 0.2, 0.4, 0.4, 0.5)},
        {'name': 'rnn_units', 'type': 'discrete',
         'domain': (50, 60, 70, 80, 90, 100, 120, 140, 160, 180)}
    ]

    if model_type == 'model1':
        bounds = basic_bounds + rnn_bounds
    elif model_type == 'model2':
        bounds = basic_bounds + advanced_bounds
    elif model_type == 'model3':
        bounds = basic_bounds + rnn_bounds
    else:
        print(f'No such model: {model_type} ---> exiting')
        exit(1)


    if has_logger:
        for bound in bounds:
            logger.info(bound)
    return bounds



# from dotmap import DotMap


def get_scales(input_scale_str, output_scale_str):
    norm_ix = [2, 3, 6, 7, 8, 9, 10]
    exp_ix = [1]
    if input_scale_str == 'None':
        input_scalar = None
    elif input_scale_str == 'norm_all':
        input_scalar = MinMaxScaler()
    elif input_scale_str == 'stand_all':
        input_scalar = StandardScaler()
    elif input_scale_str == 'norm_non_gausian':
        t = [('e', MinMaxScaler(), exp_ix)]
        input_scalar = ColumnTransformer(transformers=t, remainder='passthrough')
    elif input_scale_str == 'stand_gausian':
        t = [('n', StandardScaler(), norm_ix)]
        input_scalar = ColumnTransformer(transformers=t, remainder='passthrough')
    elif input_scale_str == 'selective':
        t = [('e', MinMaxScaler(), exp_ix), ('n', StandardScaler(), norm_ix)]
        input_scalar = ColumnTransformer(transformers=t, remainder='passthrough')

    if output_scale_str == 'None':
        output_scale = None
    elif output_scale_str == 'stand':
        output_scale = StandardScaler()
    return input_scalar, output_scale

class Opt(object):
    pass

def cat_decode(x):
    # opt = DotMap()
    opt = Opt()
    if glob_config.model_type == 'model1':
        opt.em_drop = float(x[:, 0])
        opt.fc_drop = float(x[:, 1])
        opt.batch_size = int(x[:, 2])
        opt.em_dim = int(x[:, 3])
        opt.fc_num_hidden_layers = int(x[:, 4])
        opt.fc_num_units = int(x[:, 5])
        opt.fc_activation = fc_activation_dict[str(int(x[:, 6]))]
        opt.optimizer = optimizer_dict[str(int(x[:, 7]))]
        opt.last_activation = last_activation_dict[str(int(x[:, 8]))]
        opt.initializer = initializer_dict[str(int(x[:, 9]))]
        opt.cost_function = cost_function_dict[str(int(x[:, 10]))]

        opt.rnn_drop = float(x[:, 11])
        opt.rnn_rec_drop = float(x[:, 12])
        opt.rnn_units = int(x[:, 13])

    elif glob_config.model_type == 'model2':
        print('TODO - sort bounds')
        exit(1)
    elif glob_config.model_type == 'model3':

        opt.em_drop = float(x[:, 0])
        opt.fc_drop = float(x[:, 1])

        opt.batch_size = int(x[:, 2])
        opt.em_dim = int(x[:, 3])
        opt.fc_num_hidden_layers = int(x[:, 4])
        opt.fc_num_units = int(x[:, 5])

        opt.fc_activation = fc_activation_dict[str(int(x[:, 6]))]
        opt.optimizer = optimizer_dict[str(int(x[:, 7]))]
        opt.last_activation = last_activation_dict[str(int(x[:, 8]))]
        opt.initializer = initializer_dict[str(int(x[:, 9]))]

        opt.rnn_drop = float(x[:, 10])
        opt.rnn_rec_drop = float(x[:, 11])
        opt.rnn_units = int(x[:, 12])

    return opt

def global_configurations(opt):
    opt.enzyme = glob_config.enzyme
    opt.model_type = glob_config.model_type
    opt.has_embedding = glob_config.has_embedding
    opt.has_biofeatures = glob_config.has_biofeatures
    opt.learning_rate = glob_config.learning_rate
    opt.lr_scheduler = glob_config.lr_scheduler
    opt.save_model = False
    opt.checks = glob_config.checks
    opt.simulation_type = glob_config.simulation_type
    opt.multi_data = glob_config.multi_data

    if glob_config.model_type == 'model3':
        opt.epochs = 35
        opt.cost_function = 'mse'
        opt.input_scale = None
        opt.output_scale = None

    if glob_config.model_type == 'model1':
        opt.epochs = 65
        opt.input_scale = None
        opt.output_scale = None

    return opt

def f(x):
    opt = cat_decode(x)
    opt = global_configurations(opt)


    param = {'data_type': opt.enzyme,  # mannualy set model_type as enzyme
             'em_drop': opt.em_drop,
             'fc_drop': opt.fc_drop,

             'batch_size': opt.batch_size,
             'epochs': opt.epochs,
             'em_dim': opt.em_dim,
             'fc_num_hidden_layers': opt.fc_num_hidden_layers,
             'fc_num_units': opt.fc_num_units,

             'fc_activation': opt.fc_activation,
             'optimizer': opt.optimizer,
             'cost_function': opt.cost_function,
             'last_activation': opt.last_activation,
             'initializer': opt.initializer,
             'input_scale': opt.input_scale,
             'output_scale': opt.output_scale,
             'verbose': 2}  # disable the training prints

    if glob_config.model_type != 'model2':
        param['rnn_drop'] = opt.rnn_drop
        param['rnn_rec_drop'] = opt.rnn_rec_drop
        param['rnn_units'] = opt.rnn_units



    logger.info('----------')
    global training_num
    logger.info('Starting training {}'.format(training_num))
    training_num += 1
    logger.info('params:{}'.format(param))

    start_clk = time.time()
    # model, output_scale, input_scale = train_model(opt, DataHandler=data)
    val_loss = train_model(opt, DataHandler=data)

    end_clk = time.time()

    # prepare data
    # y_valid = data['y_valid']
    # X_valid_biofeat = data['X_valid_biofeat']
    #
    # if input_scale is not None:
    #     X_valid_biofeat = input_scale.transform(X_valid_biofeat)
    #
    #
    # y_valid_pred = model.predict([data['X_valid'], X_valid_biofeat])
    # if output_scale is not None:
    #     y_valid_pred = output_scale.inverse_transform(y_valid_pred)
    #
    # if glob_config.enzyme == 'multi_task':
    #
    #     data_types = ['wt', 'esp', 'hf']
    #
    #     for ind, type in enumerate(data_types):
    #         y_valid_true_enzyme = y_valid[:, ind]
    #         y_valid_pred_enzyme = y_valid_pred[:,ind]
    #         mse = mean_squared_error(y_valid_true_enzyme, y_valid_pred_enzyme)
    #         spearman = sp.stats.spearmanr(y_valid_true_enzyme, y_valid_pred_enzyme)[0]
    #         logger.info('type: {}, spearman: {}, mse: {}'.format(type, spearman, mse))
    #         if type == 'wt':
    #             evaluation = spearman
    # else:
    #     mse = mean_squared_error(y_valid, y_valid_pred)
    #     spearman = sp.stats.spearmanr(y_valid, y_valid_pred)[0]
    #     evaluation = spearman
    #     logger.info('spearman: {}'.format(spearman))
    #     logger.info('evaluation: {}'.format(mse))
    #
    # logger.info('time: {}'.format(end_clk - start_clk))

    # if np.isnan(evaluation):
    #     return 2.0
    # else:
    #     return 1 - evaluation  # because the Bayesian Optimization will look for minimum and we want maximum.
    return val_loss



def creat_logger(enzyme):
    # create logger with 'Model_application'
    global logger
    logger = logging.getLogger('Model')
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    if not os.path.exists('logs/'):
        os.mkdir('logs/')
    ind = 1
    logfile_name = 'logs/' + enzyme + '_HyperParamSearch_1.log'
    while os.path.exists(logfile_name):
        ind += 1
        logfile_name = 'logs/' + enzyme + '_HyperParamSearch_' + str(ind) + '.log'

    fh = logging.FileHandler(logfile_name)
    fh.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    logger.info('###################################################')
    logger.info('Running hyper parameter search for data type {}'.format(enzyme))
    logger.info('###################################################')


def create_data(DataHandler):
    # X_train, X_train_biofeat, y_train = DataHandler['X_train'], DataHandler['X_train_biofeat'], DataHandler['y_train']
    # X_test, X_test_biofeat, y_test = DataHandler['X_test'], DataHandler['X_test_biofeat'], DataHandler['y_test']
    # X_train, X_val, X_train_biofeat, X_val_biofeat, y_train, y_val = dh.train_val_split(X_train, X_train_biofeat,
                                                                                        # y_train)
    X_train, X_train_biofeat, y_train = DataHandler['X_train'], DataHandler['X_train_biofeat'], DataHandler['y_train']
    X_valid, X_valid_biofeat, y_valid = DataHandler['X_valid'], DataHandler['X_valid_biofeat'], DataHandler['y_valid']

    global data

    if glob_config.data_source == 'new':
        conf_train, conf_valid = DataHandler['conf_train'], DataHandler['conf_valid']

        data = {'X_train': X_train, 'X_train_biofeat': X_train_biofeat, 'y_train': y_train,'conf_train': conf_train, 'X_valid': X_valid,
                'X_valid_biofeat': X_valid_biofeat, 'y_valid': y_valid, 'conf_valid': conf_valid}
    else:
        data = {'X_train': X_train, 'X_train_biofeat': X_train_biofeat, 'y_train': y_train, 'X_valid': X_valid,
                'X_valid_biofeat': X_valid_biofeat, 'y_valid': y_valid}


# time consuming
# set initial_design_numdata and  max_iter as lower integer to reduce searching time

def save_results(enzyme, opt_model, results_dir):
    ind = 0
    save_path = results_dir + "HPS_{}.csv".format(ind)
    while os.path.exists(save_path):
        ind += 1
        save_path = results_dir + "HPS_{}.csv".format(ind)
    opt_model.save_evaluations(save_path)


def param_search(config, DataHandler):

    # create globals
    global glob_config
    glob_config = config
    creat_logger(config.enzyme)
    create_data(DataHandler)
    global training_num
    training_num = 0

    bounds = get_bounds(config.enzyme, model_type=config.model_type)

    logger.info('Starting Bayesian Optimization')

    results_dir = 'HPS/'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    results_dir += config.data_source + '/'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    results_dir += config.model_type + '/'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    if config.has_biofeatures:
        results_dir += 'bio/'
    else:
        results_dir += 'no_bio/'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    results_dir += config.enzyme + '/'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    if config.enzyme == 'multi_task':
        results_dir += config.multi_data + '/'
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)



    results_path = results_dir + 'HPS.csv'
    if os.path.isfile(results_path):
        logger.info('found previeus csv')
        evals = pd.read_csv(results_path, index_col=0, delimiter="\t")
        Y = np.array([[x] for x in evals["Y"]])
        X = np.array(evals.filter(regex="var*"))
        opt_model = GPyOpt.methods.BayesianOptimization(f=f, domain=bounds, initial_design_numdata=0, X=X, Y=Y)
        opt_model.run_optimization(max_iter=5, max_time=72000)
        save_results(config.enzyme, opt_model, results_dir)
    else:
        logger.info('Didnt find previeus csv, beggigning new search with initial_design_numdata = 5')
        opt_model = GPyOpt.methods.BayesianOptimization(f=f, domain=bounds, initial_design_numdata=5)
        save_results(config.enzyme, opt_model, results_dir)

    logger.info("optimized loss: {0}".format(opt_model.fx_opt))
    for i, v in enumerate(bounds):
        name = v['name']
        val = opt_model.x_opt[i]

        if name == 'fc_activation':
            val = fc_activation_dict[str(int(val))]
        if name == 'optimizer':
            val = optimizer_dict[str(int(val))]
        if name == 'cost_function':
            val = cost_function_dict[str(int(val))]
        if name == 'last_activation':
            val = last_activation_dict[str(int(val))]
        if name == 'initializer':
            val = initializer_dict[str(int(val))]

        logger.info('parameter {}:{}'.format(name, val))


















