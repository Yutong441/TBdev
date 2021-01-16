import sys
import GrandPrixModel
import h5py
import re
import pandas as pd

def fit_model(
        data,
        n_latent_dims=1,
        n_inducing_points=10,
        kernel={'name':'RBF', 'ls':1.0, 'var':1.0},
        mData=None,
        latent_prior_mean=None,
        latent_prior_var=1.,
        latent_mean=None,
        latent_var=0.1,
        inducing_inputs=None,
        fix_parameters = None,
        predict=None,
        jitter=1e-6,
        dtype='float64',
        disp = True,
        maxitr = 1000):
    sys.argv.append(dtype)
    m = GrandPrixModel.GrandPrixModel(data, n_latent_dims, n_inducing_points,
            kernel, mData, latent_prior_mean, latent_prior_var, latent_mean,
            latent_var, inducing_inputs, dtype)
    m.set_jitter_level(jitter)
    m.set_trainable(fix_parameters)
    m.build()

    m.fit(maxiter=maxitr, display=disp)
    posterior = m.get_latent_dims()
    return posterior, m

def save_model(model, save_file):
    '''
    Reference:
    https://gpflow.readthedocs.io/en/master/notebooks/intro_to_gpflow2.html#Saving-and-loading-models
    https://stackoverflow.com/questions/54142917/saving-and-retrieving-the-parameters-of-a-gpflow-model
    '''
    vars = trainable_dict(model)
    with h5py.File(save_file, 'w') as f:
        for name, value in vars.items():
            f[name] = value.read_value ()

def load_model(model, load_file):
    """
    Load a model given by model path
    """

    vars = {}
    def _gather(name, obj):
        if isinstance(obj, h5py.Dataset):
            vars[name] = obj[...]

    with h5py.File(load_file) as f:
        f.visititems(_gather)

    model.assign(vars)

def trainable_dict (model):
    return {param.full_name: param for param in model.parameters if param.trainable }

def print_kernel_summary (model):
    regex = re.compile ('^.*/kern[e]?[l]?/')
    key_list, val_list = [], []
    for val in model.parameters:
        key = val.full_name
        if val.trainable:
            if len (regex.findall (key)) != 0:
                append_key = regex.sub ('', key) 
                key_list.append (append_key)
                val_list.append (val.read_value ())

    return pd.DataFrame (val_list, index=key_list, columns=['value'])

