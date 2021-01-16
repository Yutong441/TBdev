import numpy as np
import pandas as pd
import tensorflow as tf
import gpflow
#import ChangePoint as CP
import BRGP_class as bc
import h5py
import re

def MapTo01(y) -> tf.Tensor:
    '''Scale Y into range [0, 1]'''
    return (y - np.min (y, 0, keepdims=True)) / (np.max(y, 0,
        keepdims=True) - np.min(y, 0, keepdims=True))

def zero_kern ():
    kernel_base = gpflow.ekernels.Linear (1)
    kernel_base.variance = 0.
    kernel_base.variance.trainable = False
    return kernel_base

def get_kernel (): 
    '''
    Hierarchical GP usually requires a base kernel k1 from which other kernels
    branch off. This function includes the change point kernels as the
    branching kernels. The change point kernel switches from a zero kernel into
    a RBF.
    '''
    base_kernel = gpflow.kernels.RBF (variance=tf.cast(1, dtype=tf.float64)) 
    kernel1 = gpflow.kernels.RBF (variance=tf.cast(0.2, dtype=tf.float64))
    kernel2 = gpflow.kernels.RBF (variance=tf.cast(0.2, dtype=tf.float64))

    #cp_kernel1 = gpflow.kernels.ChangePoints ([zero_kernel(), kernel1],
    #        locations=[4.], steepness=1.65)
    #cp_kernel2 = gpflow.kernels.ChangePoints ([zero_kernel(), kernel2],
    #        locations=[4.], steepness=1.65)
    cp_kernel1 = CP.tanh_CP (kernels=[bc.zero_kernel(), kernel1],
            locations=[4.], steepness=[1.65])
    cp_kernel2 = CP.tanh_CP (kernels=[bc.zero_kernel(), kernel2],
            locations=[4.], steepness=[1.65])
    CP_loc = gpflow.params.Parameter ([4.], dtype=tf.float64)
    # tieing the branching time
    cp_kernel1.locations = CP_loc
    cp_kernel2.locations = CP_loc
    return bc.zero_kernel(), kernel1, kernel2 

def get_BP_pred (model, Xnew, branch_num=1, base_kernel=None):
    M = model.num_kernel 
    assert branch_num <= M, 'branch number must be no larger than {}'.format (M)
    if len (Xnew.shape) == 1: Xnew = Xnew.reshape ([-1, 1])
    else: assert Xnew.shape[1] == 1, 'Xnew must be an 1D numpy array'
    BP_new = np.zeros ( [Xnew.shape[0], M])
    BP_new [:, branch_num] = 1
    if base_kernel is not None:
        for i in base_kernel:
            BP_new [:, i] = 1
    return BP_new

def toy_function (z, t=None, sigma=0.1):
    '''Implement the toy example in Penfold 2018
    Args:
        `t`: a 1D array of times
        `z`: branch, either 1 or 2
    Examples:
    >>> t, y_val = toy_function (z=2)
    '''
    if t is None: t = np.random.normal(0, 3, 200)
    y = np.zeros (t.shape) # y is the final output
    y [t <= -np.pi/2] = 0.
    select_t = np.logical_and (t > -np.pi/2, t <= 0.)
    y [select_t] = np.cos ( t  )[select_t] 
    y [t > 0.] = 1.
    if z == 2: y = -y
    return t, y + np.random.normal (0, sigma, y.shape)

def toy_data (t=None, sigma=0.1, toy_func = toy_function):
    '''
    Examples:
    >>> t, y_val, branch = toy_data ()
    >>> plot_data = pd.DataFrame ( [t, y_val, branch] ).T
    >>> plot_data.columns = ['time', 'y', 'branch']
    >>> import seaborn as sb
    >>> sb.scatterplot (data=plot_data, x='time', y='y', hue='branch')
    >>> plt.show()
    >>> plot_data.to_csv ('toy_branch.csv')

    For demo data:
    >>> Y = pd.read_csv (root + 'demo_data.txt', sep='\t')
    >>> # process column names to obtain the grouping information
    >>> # example column names: 'Time 1 Mock BioRep A', 'Time 13 Pathogen 2 BioRep D'
    >>> split_time = [one_name.split ('Time ')[1] for one_name in Y.columns]
    >>> pathogen = ['_'.join (one_name.split (' ')[1:]) for one_name in split_time]
    >>> pathogen = [one_name.split ('_BioRep')[0] for one_name in pathogen]
    >>> # one hot encoding
    >>> pathogen = np.array (pathogen).reshape ([-1, 1])
    >>> one_hot_lab = (pathogen == np.unique (pathogen).reshape ([1, -1])).astype (float)
    '''
    t1, y1 = toy_func (z=1, sigma=sigma)
    t2, y2 = toy_func (z=2, sigma=sigma)
    all_t = np.concatenate ([t1, t2], axis=0)
    all_y = np.concatenate ([y1, y2], axis=0)
    branch = np.concatenate ( [np.ones (t1.shape), np.ones (t2.shape)+1.],
            axis=0 )
    return all_t, all_y, branch

def toy_recombinant (z, t=None, sigma=0.1):
    '''
    >>> import BRGP_utils as bu
    >>> plot_data = bu.toy_data (toy_func = bu.toy_recombinant)
    >>> plot_data = pd.DataFrame (plot_data).T
    >>> plot_data.columns = ['time', 'y', 'branch']
    >>> sb.scatterplot (data=plot_data, x='time', y='y', hue='branch')
    >>> plt.show ()
    >>> plot_data.to_csv ('toy_recombinant.csv')
    '''
    if t is None: t = np.random.normal(0, 3, 200)
    y = np.zeros (t.shape) # y is the final output
    y [t <= -np.pi/2] = 0.
    select_t = np.logical_and (t > -np.pi/2, t <= np.pi/2)
    y [select_t] = np.cos ( t  )[select_t] 
    y [t > np.pi/2] = 0.
    if z == 2: y = -y
    return t, y + np.random.normal (0, sigma, y.shape)

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
        if hasattr (model, 'BP_Z'):
            f['BP_Z'] = model.BP_Z
        if hasattr (model, 'Z'):
            f['Z'] = model.Z.read_value ()

def load_model(model, load_file):
    """
    Load a model given by model path
    """

    var_dict = {}
    def _gather(name, obj):
        if isinstance(obj, h5py.Dataset):
            var_dict [name] = obj[...]

    with h5py.File(load_file) as f:
        f.visititems(_gather)

    model.BP_Z = var_dict ['BP_Z']
    model.Z = var_dict ['Z']
    var_dict.pop ('BP_Z', None)
    var_dict.pop ('Z', None)
    model.assign(var_dict)

def KL_divergence_normal (mu1, mu2, var1, var2, eps=1e-6):
    NQ = 1
    KL = -0.5 * np.log(var1+eps)
    KL += 0.5 * np.log(var2+eps)
    KL -= 0.5 * NQ
    KL += 0.5 * ((mu1 - mu2)**2 + var1) / var2 
    return KL

def DE_by_KL(x, num1, num2):
    '''Calculate the branching score
    However, I prefer to do this in R, which integrates with ggplot better.
    Hence, I only implemented a crude version of this function.
    Example:
    >>> KL_gene = pu.DE_by_KL (br_model, 1, 2)
    '''
    pred_b1 = x.prediction [x.prediction['branch'].isin(['branch'+str(num1)])]
    pred_b2 = x.prediction [x.prediction['branch'].isin(['branch'+str(num2)])]
    pred_mean_b1 = pred_b1 [ 'mean_'+ br_model.features ]
    pred_var_b1  = pred_b1 [ 'var_' + br_model.features ]
    pred_mean_b2 = pred_b2 [ 'mean_'+ br_model.features ]
    pred_var_b2  = pred_b2 [ 'var_' + br_model.features ]
    KL_branch = KL_divergence_normal (pred_mean_b1.values, pred_mean_b2.values,
            pred_var_b1.values, pred_var_b2.values)
    return KL_branch

    #branch_genes = np.percentile (KL_branch, 95, axis=0)
    #branch_genes = pd.DataFrame (branch_genes, index=br_model.features)
    #branch_genes.columns = ['genes']
    #sort_genes = branch_genes.sort_values ('genes', ascending=False)
    #sort_genes.index [:100]

def duplicate_main (x_list, du_index):
    x_du_list = []
    for x in x_list:
        x_du = x [du_index, :]
        x_du_list.append (np.concatenate ([x, x_du], axis=0))
    return x_du_list

def duplicate_main_tf (x_list, du_index):
    x_du_list = []
    for x in x_list:
        x_du = tf.gather (x, du_index)
        x_du_list.append (tf.concat ([x, x_du], axis=0))
    return x_du_list
