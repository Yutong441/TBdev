import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
import gpflow
import tensorflow as tf

Y = pd.read_csv ('toy_branch.csv')
Y = Y.iloc [np.random.permutation (Y.shape[0]), :]
#sb.scatterplot (data=Y, x='time', y='y', hue='branch')
branch = Y['branch'].values.reshape ([-1, 1])
time_prior = Y['time'].values.reshape ([-1, 1])
Y_data = Y['y'].values.reshape ([-1, 1])

import BRGP_main as bm 
import BRGP_utils as bu
branch_num = np.unique (branch)
branch = branch.reshape ([-1, 1]) == branch_num.reshape ([1, -1])
branch_inp = np.insert (branch, 0, branch.sum (1, keepdims=True).T, axis=1)
br_model = bm.BRGP_main (Y_data, bu.get_kernel (), branch_inp, num_inducing=50,
        X_prior_mean = time_prior)

br_model.compile ()
br_model.fit (1000)
br_model.predict (min_x=time_prior.min(), max_x=time_prior.max())
br_model.plot_prediction ()
