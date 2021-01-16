import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import GrandPrix

Y_data = pd.read_csv ('data/merged_all.csv', index_col=0)
Y = Y_data.T.values
Y -= Y.mean (0, keepdims=True)
Y = Y [:, Y.std (0) != 0] #remove genes of zero variance
Y /= Y.std (0, keepdims=True)
metadata = pd.read_csv ('data/merged_meta_all.csv', index_col=[0])
new_meta = metadata

(pt, var), model = GrandPrix.fit_model(data=Y, n_latent_dims=3,
                          n_inducing_points=10, 
                          latent_prior_var=1.**2, 
                          kernel={'name':'RBF', 'ls':1.0, 'var':1.0}, 
                          maxitr=10000,
                          disp=True)

GrandPrix.save_model (model.m, 'result/model_no_prior.hdf5')
# save pseudotime
PT_save = pd.DataFrame (np.concatenate ([pt, var], axis=1))
PT_save.columns = ['PT1', 'PT2', 'PT3', 'PT1_var', 'PT2_var', 'PT3_var']
PT_save.index = new_meta.index
PT_save.to_csv ('result/PT_no_prior.csv')

pt_pred_out, var_pred_out = model.m.predict_y (pt)
pt_pred = pd.DataFrame (pt_pred_out.numpy())
var_pred = pd.DataFrame (var_pred_out.numpy())
pt_pred.index = new_meta.index
var_pred.index = new_meta.index
pt_pred.columns = Y_data.index 
var_pred.columns = Y_data.index 

pt_pred.to_csv ('result/PT_pred_mean.csv')
var_pred.to_csv ('result/PT_pred_var.csv')
