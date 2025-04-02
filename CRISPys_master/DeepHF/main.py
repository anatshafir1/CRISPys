"""standard command line:
     --data_type multi_task -s_type train -em --data_file final_df_1
     --enzyme multi_task -s_type train -ds new --model_type model3 --weighted_loss -md serial
     --enzyme esp -s_type train -ds new --model_type model3 --weighted_loss
     --enzyme esp -s_type param_search -ds new --model_type model1  -md parallel --weighted_loss
"""
if __name__ == '__main__':
    import pandas as pd
    import tensorflow as tf
    from scripts.training_util import *
    # my scripts
    import scripts.configurations as cfg
    import scripts.Cross_validation as Cross_validation
    import scripts.HyperParameterSearching as HPS
    from scripts import testing
    from scripts import preprocess
    from scripts import data_handler as dh
    from scripts import prediction_util as pred

    # Using only needed memory on GPU
    conf = tf.compat.v1.ConfigProto()
    conf.gpu_options.allow_growth = True
    tf.compat.v1.keras.backend.set_session(tf.compat.v1.Session(config=conf))

    # Get configs
    config = cfg.get_parser()

    if config.simulation_type == 'train':
        print('simulation_type = train')
        DataHandler = dh.get_data(config)
        train_model(config, DataHandler)

    if config.simulation_type == 'cross_v':
        param = get_param(config)
        DataHandler = dh.get_data(config)
        Cross_validation.cross_validation(param, DataHandler)

    if config.simulation_type == 'param_search':
        # param = get_param(config)
        DataHandler = dh.get_data(config)
        HPS.param_search(config, DataHandler)

    if config.simulation_type == 'test_model':
        testing.test_model(config)

    if config.simulation_type == 'test_means':
        testing.test_means(config)

    if config.simulation_type == 'preprocess':
        preprocess.split_to_train_test_file(config)
        preprocess.prepare_test_seq(config)
        preprocess.prepare_train_valid_seq(config)

    if config.simulation_type == 'predict':
        model_path = 'models/{}/{}/{}/'.format(config.model_type, 'bio' if config.has_biofeatures else 'no_bio',
                                               config.enzyme)
        if config.enzyme == 'multi_task':
            model_path += config.multi_data + '/'

        DataHandler = dh.get_data(config)
        if config.has_biofeatures:
            inputs = [DataHandler['X_test'], DataHandler['X_test_biofeat']]
        else:
            inputs = [DataHandler['X_test']]
        # pred.get_predictions(model_path, inputs, decoded=True)
        pred.debug(model_path, inputs, DataHandler['y_test'])
        a = 0
