import argparse
from keras.optimizers import *
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.compose import ColumnTransformer

'''
Model Types:
    model1 - This is the original model, or the mixed model when 'enzyme'='multi_task'
            Inputs: Sequence + biofetures
            Output: Efficiency per enzyme
            Architecture: embeding -> lstm -> FC
              
    model2 - This model is using only the biofeatures as inputs and have only FC layers
            Inputs: biofetures
            Output: Efficiency per enzyme
            Architecture: FC
'''
def get_parser():
    print('Receiving parser values')
    parser = argparse.ArgumentParser(description='Process some integers.')

    # Data
    parser.add_argument('--enzyme', type=str, default=None)  # [wt,esp,hf,multi_task]
    parser.add_argument('--data_file', type=str, default='row_data')  # original file = 'row_data'
    parser.add_argument('-ds', '--data_source', type=str, default='new')  # For original data choose old  ['old', 'new']
    parser.add_argument('-md', '--multi_data', type=str, default='parallel')  # ['parallel', 'serialized] this option is for multi task only


    # Simulation
    parser.add_argument('-s_type', '--simulation_type', type=str, default=None)  # [train, param_search, cross_v, test_model, preprocess]
    parser.add_argument('--model_path', type=str, default='None')  # For the test_model option
    parser.add_argument('--eda', dest='eda', action='store_true')  # For the test_model option


    # Training
    parser.add_argument('--lr_scheduler', dest='lr_scheduler', action='store_true')


    # Model
    parser.add_argument('--dont_save_model', dest='save_model', action='store_false')
    parser.add_argument('-no_em', '--has_embedding', dest='has_embedding', action='store_false')
    parser.add_argument('-no_bio', '--no_biofeatures', dest='has_biofeatures', action='store_false')
    parser.add_argument('-mt', '--model_type', type=str, default='model1')
    parser.add_argument('--weighted_loss', dest='weighted_loss', action='store_true')  # weighted_loss / row_reads model






    config = parser.parse_args()
    # Sanity check
    if config.enzyme is None and config.simulation_type != 'preprocess':
        print('No data type received >>>> exiting')
        exit(1)

    if config.simulation_type is None:
        print('No simulation type received >>>> exiting')
        exit(1)

    config = get_optimized_params(config)

    if config.lr_scheduler:
        config.learning_rate = 0.01
    return config

def get_optimized_params(config):
    # This method will define the model optimized hyper parameters that were found in the Hyper Parameter search

    if config.has_biofeatures:
        if config.enzyme == 'wt':
            config.em_drop = 0.3
            config.rnn_drop = 0.6
            config.rnn_rec_drop = 0.2
            config.fc_drop = 0.1
            config.batch_size = 70
            config.epochs = 40
            config.em_dim = 40
            config.rnn_units = 200
            config.fc_num_hidden_layers = 5
            config.fc_num_units = 50
            config.fc_activation = 'elu'
            config.optimizer = Adamax
            config.learning_rate = 0.001
            config.last_activation = 'linear'
            config.cost_function = 'mse'
            config.initializer = 'glorot_uniform'
            config.input_scale = None
            config.output_scale = None

        elif config.enzyme == 'esp':
            config.em_drop = 0.2
            config.rnn_drop = 0.2
            config.rnn_rec_drop = 0.4
            config.fc_drop = 0.4
            config.batch_size = 70
            config.epochs = 50
            config.em_dim = 44
            config.rnn_units = 70
            config.fc_num_hidden_layers = 2
            config.fc_num_units = 300
            config.fc_activation = 'sigmoid'
            config.optimizer = Adamax
            config.learning_rate = 0.001
            config.last_activation = 'linear'
            config.cost_function = 'mse'
            config.initializer = 'glorot_uniform'
            config.input_scale = None
            config.output_scale = None


        elif config.enzyme == 'hf':
            config.em_drop = 0.2
            config.rnn_drop = 0.5
            config.rnn_rec_drop = 0.4
            config.fc_drop = 0.5
            config.batch_size = 80
            config.epochs = 45
            config.em_dim = 48
            config.rnn_units = 80
            config.fc_num_hidden_layers = 2
            config.fc_num_units = 300
            config.fc_activation = 'sigmoid'
            config.optimizer = Adamax
            config.learning_rate = 0.001
            config.last_activation = 'linear'
            config.cost_function = 'mse'
            config.initializer = 'glorot_uniform'
            config.input_scale = None
            config.output_scale = None


        elif config.enzyme == 'multi_task':
            config.em_drop = 0.1
            config.rnn_drop = 0.2
            config.rnn_rec_drop = 0.1
            config.fc_drop = 0.5
            config.batch_size = 100
            config.epochs = 65
            config.em_dim = 36
            config.rnn_units = 120
            config.fc_num_hidden_layers = 4
            config.fc_num_units = 400
            config.fc_activation = 'sigmoid'
            config.optimizer = RMSprop
            config.cost_function = 'mse'
            config.last_activation = 'sigmoid'
            config.initializer = 'normal'

            norm_ix = [2, 3, 6, 7, 8, 9, 10]
            exp_ix = [1]
            t = [('e', MinMaxScaler(), exp_ix), ('n', StandardScaler(), norm_ix)]
            config.input_scale = StandardScaler()
            config.output_scale = None

            config.learning_rate = 0.001

    # No biofeatures model
    else:
        config.em_drop = 0.2
        config.rnn_drop = 0.4
        config.rnn_rec_drop = 0.4
        config.fc_drop = 0.4
        config.batch_size = 120
        config.epochs = 65
        config.em_dim = 32
        config.rnn_units = 70
        config.fc_num_hidden_layers = 1
        config.fc_num_units = 280
        config.fc_activation = 'hard_sigmoid'
        config.optimizer = RMSprop
        config.cost_function = 'mse'
        config.last_activation = 'linear'
        config.initializer = 'normal'

        norm_ix = [2, 3, 6, 7, 8, 9, 10]
        exp_ix = [1]
        t = [('e', MinMaxScaler(), exp_ix), ('n', StandardScaler(), norm_ix)]
        config.input_scale = ColumnTransformer(transformers=t, remainder='passthrough')
        config.output_scale = StandardScaler()

        config.learning_rate = 0.001

    if config.data_source == 'new':
        config.checks = ['on_batch_end']
        config.batch_size = 1000
        config.epochs = 65
        config.last_activation = 'sigmoid'
        config.cost_function = 'binary_crossentropy'

    else:
        config.checks = ['on_epoch_end']

    if config.data_source == 'new':
        if config.enzyme == 'esp':
            config.em_drop = 0.2
            config.rnn_drop = 0.2
            config.rnn_rec_drop = 0.2
            config.fc_drop = 0.4
            config.batch_size = 110
            config.epochs = 65
            config.em_dim = 64
            config.rnn_units = 120
            config.fc_num_hidden_layers = 2
            config.fc_num_units = 200
            config.fc_activation = 'elu'
            config.optimizer = RMSprop
            config.cost_function = 'mse'
            config.last_activation = 'sigmoid'
            config.initializer = 'he_uniform'

            config.checks = ['on_epoch_end']

        elif config.enzyme == 'multi_task':
            config.em_drop = 0.3
            config.rnn_drop = 0.2
            config.rnn_rec_drop = 0.2
            config.fc_drop = 0.2
            config.batch_size = 90
            config.epochs = 65
            config.em_dim = 32
            config.rnn_units = 140
            config.fc_num_hidden_layers = 2
            config.fc_num_units = 260
            config.fc_activation = 'elu'
            config.optimizer = Adamax
            config.cost_function = 'mse'
            config.last_activation = 'linear'
            config.initializer = 'he_uniform'
            config.input_scale = None
            config.checks = ['on_epoch_end']


    return config