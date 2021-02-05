import numpy as np
import pandas as pd
import anndata

def py_mean_log_likelihood (mu, var, x, return_mean=True):
    out = mu**2 + (x**2).mean (axis=0, keepdims=True) -2*mu*x.mean (axis=0,
            keepdims=True)
    out /= -var*2 
    out = out - np.sqrt (var) - 0.5*np.log (2*np.pi)
    if return_mean: return out.mean (axis=1)
    else: return out

def get_likelihood (x, cell_type, feature, mu, var, return_mean=True):
    return np.exp (py_mean_log_likelihood (mu, var, x.X [x.obs [feature] ==
        cell_type], return_mean=return_mean))

def get_mat_likelihood (x, mu, var, return_mean=True):
    return np.exp (py_mean_log_likelihood (mu, var, x, return_mean=return_mean))

def fitting_likelihood (ori_data, new_X, new_meta, group_by, pt_mean, pt_var):
    ori_mean = ori_data.values.mean(1, keepdims=True)
    ori_std  = ori_data.values.std(1, keepdims=True)

    scaled_X = pd.DataFrame ((new_X.values - ori_mean)/ori_std)
    scaled_X.columns = new_meta.index
    scaled_X.index = new_X.index
    all_data = anndata.AnnData (scaled_X.T, obs=new_meta)

    all_types = np.unique (all_data.obs[group_by])
    all_prob = {}
    for one_type in all_types:
        all_prob[one_type] = get_likelihood (all_data, one_type, 'assigned_cluster',
                pt_mean.values, pt_var.values)
    return pd.DataFrame (all_prob)
