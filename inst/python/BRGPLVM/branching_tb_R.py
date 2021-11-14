import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import gpflow
import BRGP_main as bm 
import BRGP_class as bc
import BRGP_utils as bu
import ChangePoint as CP
import argparse

def initialise (args):
    Y = pd.read_csv (args.root+'/STREAM/STREAM_data.csv', index_col=[0])
    time_prior =pd.read_csv (args.root+'/STREAM/STREAM_pseudotime.csv', index_col=[0]).values 
    branch = pd.read_csv (args.root+'/STREAM/STREAM_branch_labels.csv', index_col=[0])
    Y_data = Y.T.loc[branch.index.values].values 

    branch_num = np.unique (branch.values)
    branch = branch.values.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
    branch_inp = branch.astype (float)

    du_index = branch_inp [:, 0] == 1
    branch_inp, = bu.duplicate_main ([branch_inp], du_index)
    added_len = du_index.astype(int).sum()
    du_branch = np.concatenate ([du_index, np.repeat (False, added_len )] )
    branch_inp [du_branch, 2] = 1
    branch_inp [-added_len:, 1] = 1
    branch_inp [:,0] = 1

    kernel0 = gpflow.kernels.Matern32 (1, variance=0.1)
    kernel1 = gpflow.kernels.Matern32 (1, variance=0.1)
    kernel2 = gpflow.kernels.Matern32 (1, variance=0.1)
    kern_list = gpflow.ParamList([bu.zero_kern(), kernel1, kernel2])

    br_model = bm.BRGP_main (Y_data, kern_list, branch_inp, num_inducing=30,
            X_prior_mean=time_prior, tie_index=du_index)

    br_model.model.compile ()
    return br_model

def training (args):
    br_model = initialise (args)
    br_model.fit (maxiter=args.maxiter) # done on a GPU
    bu.save_model (br_model.model, args.root+'/BRGP/model_dict_matern.hdf5')

def testing (args):
    # the steps below are done on a CPU
    br_model = initialise (args)
    bu.load_model (br_model.model, args.root+'/BRGP/model_dict_matern.hdf5')
    bu.print_kernel_summary (br_model.model)

    br_model.predict (num=500, feature_names=Y.index, with_noise=False,
            min_perc=2, max_perc=98, base_kernel=[0], kernel_pred=[1,2])

    #br_model.plot_prediction ('CGA')
    br_model.prediction.to_csv (args.root+'/BRGP/prediction_matern_500.csv')
    br_model.save_pt (args.root+'/BRGP/infer_pt_matern.csv')

if __name__ == '__main__': 
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=str)
    parser.add_argument('--maxiter', type=int)
    parser.add_argument('--mode', type=str)
    args = parser.parse_args()
    if args.mode == 'training': training (args)
    else: testing (args)
