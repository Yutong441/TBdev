import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import GrandPrix
import GrandPrixModel

Y_data = pd.read_csv ('data/merged_.csv', index_col=0)
Y = Y_data.T.values
Y -= Y.mean (0, keepdims=True)
Y = Y [:, Y.std (0) != 0] #remove genes of zero variance
Y /= Y.std (0, keepdims=True)
metadata = pd.read_csv ('data/merged_meta_.csv', index_col=[0])

# load model
model = GrandPrixModel.GrandPrixModel (Y, 3, 10)
model.build()
GrandPrix.load_model (model.m, 'result/model_no_prior.hdf5')
GrandPrix.print_kernel_summary (model.m)
model.m.likelihood.variance.read_value ()

PT_save = pd.read_csv ('result/PT_no_prior.csv', index_col=[0])
pt = PT_save [['PT1', 'PT2', 'PT3']].values

pt_pred_out, var_pred_out = model.m.predict_y (pt)
pt_pred = pd.DataFrame (pt_pred_out)
var_pred = pd.DataFrame (var_pred_out)

pt_pred.index = metadata.index
var_pred.index = metadata.index
pt_pred.columns = Y_data.index 
var_pred.columns = Y_data.index 
print ('saving data')
pt_pred.to_csv ('result/PT_pred_mean.csv')
var_pred.to_csv ('result/PT_pred_var.csv')
