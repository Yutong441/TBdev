import numpy as np
import time
import sys, os, re
import pandas as pd
from shutil import copyfile
import tensorflow as tf
import gpflow

def MapTo01(y):
    return (y.copy() - y.min(0)) / (y.max(0) - y.min(0))

class GrandPrixModel(object):
    """
    Parameters
    ----------
    data : `array-like`, shape `N` x `D`
        Observed data, where `N` is the number of samples and `D` is the number of features.
    n_latent_dims : `int`, optional (default: 1)
        Number of latent dimentions to compute.
    n_inducing_points : `int`, optional (default: 10)
        Number of inducing or auxiliary points. 
    kernel : ``gpflow.kernels`` object, optional (default: RBF kernel with lengthscale and variance set to 1.0)
        Kernel functions are used to compute the covariance among datapoints. They impose constraints such as smoothness, periodicity on the function being
        learned that is shared by all datapoints.
        Kernels are parameterize by a set of hyperparameters, i.e. lengthscale, variance, etc which can be optimized during model fitting.
    latent_prior_mean : `array-like`, shape `N` x `n_latent_dims`, optional (default: 0)
        Mean of the prior distribution over the latent dimensions.
    latent_prior_var : `array-like`, shape `N` x `n_latent_dims`, optional (default: 1.)
        Variance of the prior distribution over the latent dimensions.
    latent_mean : `array-like`, shape `N x n_latent_dims`, optional (default: PCA)
        Initial mean values of the distribution over the latent dimensions.
    latent_var : `array-like`, shape `N` x `n_latent_dims`, optional (default: 0.1)
        Initial variance of the distribution over the latent dimensions.
    inducing_inputs : `array-like`, shape `n_inducing_points` x `n_latent_dims`, optional (default: random subset from laten_mean)
        Set of inducing or auxiliary input points.
    dtype : `str`, optional (default: 'float64')
        Floating point data type precision to be used.
    """
    def __init__(self, data, n_latent_dims=1, n_inducing_points=10,
            kernel={'name':'RBF', 'ls':1.0, 'var':1.0}, mData=None,
            latent_prior_mean=None, latent_prior_var=1., latent_mean=None,
            latent_var=0.1, inducing_inputs=None, dtype='float64',
            periodic=False):
        self.Y = None
        self.Q = n_latent_dims
        self.M = n_inducing_points
        self.kern = None
        self.mData = mData
        self.X_prior_mean = None
        self.X_prior_var = None
        self.X_mean = None
        self.X_var = None
        self.Z = None

        self.set_Y(data)
        self.N, self.D = self.Y.shape

        self.set_kern(kernel, periodic=periodic)

        self.set_X_prior_mean(latent_prior_mean)
        self.set_X_prior_var(latent_prior_var)

        self.set_X_mean(latent_mean, latent_prior_mean)
        self.set_X_var(latent_var)

        self.set_inducing_inputs(inducing_inputs)

        self.fitting_time = 0

        with gpflow.defer_build():
            self.m = gpflow.models.BayesianGPLVM(Y=self.Y, kern=self.kern,
                    X_prior_mean=self.X_prior_mean,
                    X_prior_var=self.X_prior_var, X_mean=self.X_mean.copy(),
                    X_var=self.X_var.copy(), Z=self.Z.copy(), M=self.M)
            self.m.likelihood.variance = 0.01

    def __str__(self):
        return "%s"%(self.m)

    def build(self):
        r"""
        Build the model into a Tensorflow graph.
        
        """

        self.m.compile()

    def fit(self, maxiter=1000, display=False):
        r"""
        Fit the BGPLVM model.
        
        Parameters
        ----------
        maxiter : `int`, optional (default: 1000)
            Maximum number of iterations to perform.
        display : *bool*, optional (default: False)
            If set to True, print convergence messages. 
        
        """

        opt = gpflow.train.ScipyOptimizer()

        try:
            t0 = time.time()
            opt.minimize(self.m, maxiter=maxiter, disp=display)
            self.fitting_time = time.time() - t0
        except:
            print('Warning: The model terminates abnormally...')


    def predict(self, Xnew):
        r"""
        Predict posterior mean and variance using the BGPLVM. The prediction can also be done on the unfitted model using the Gaussian Process prior. 
        
        Parameters
        ----------
        Xnew : `array-like`, shape `n_sample` x `n_latent_dims`
            `n_sample` is the number of query points where the prediction will be evaluated.
        
        Returns
        -------
        data_mean : `array-like`, shape `N` x `D`
            Mean values of the predictive distribution at the query points.
        data_var : `array-like`, shape `N` x `D`
            Variance of the predictive distribution at the query points.
        
        """

        if type(Xnew) is int:
            latent_mean = self.get_latent_dims()[0]
            n_points = Xnew
            del Xnew
            Xnew = np.linspace(min(latent_mean), max(latent_mean), n_points)[:, None]
        # assert isinstance(Xnew, np.ndarray)
        return self.m.predict_y(Xnew)

    def get_latent_dims(self):
        r"""
        Get predictive distribuiton over latent dimensions of the Gaussian Process Latent Variable Model.
        
        Returns
        -------
        latent_mean : `array-like`, shape `N` x `n_latent_dims`
            Mean values of the predictive distribution over the latent dimensions.
        latent_var : `array-like`, shape `N` x `n_latent_dims`
            Variance of the predictive distribution over the latent dimensions.
        
        """

        return (self.m.X_mean.read_value()[:, 0:self.Q], self.m.X_var.read_value()[:, 0:self.Q])

    def get_model(self):
        return self.m.as_pandas_table()

    def set_trainable(self, paramlist=None):
        """
        Fix model (auxiliary) parameters and hyperparameters.
        Parameters
        ----------
        paramlist : `list`, optional (default: None)
        :param paramlist:
        :return:
        """
        if paramlist is not None:
            if 'kernel_lengthscales' in paramlist:  self.m.kern.lengthscales.trainable = False
            if 'kernel_variance' in paramlist:  self.m.kern.variance.trainable = False
            if 'likelihood_variance' in paramlist: self.m.likelihood.variance.trainable = False
            if 'inducing_inputs' in paramlist:  self.m.feature.Z.trainable = False
            if 'latent_mean' in paramlist:  self.m.X_mean.trainable = False
            if 'latent_variance' in paramlist: self.m.X_var.trainable = False

    def set_jitter_level(self, jitter_level):
        gpflow.settings.numerics.jitter_level = jitter_level

    def set_Y(self, data):
        self.Y = data

    def set_kern(self, kernel, periodic=False):
        if kernel is not None:
            if 'name' in kernel:
                kernelName = kernel['name']
            if 'ls' in kernel:
                ls = kernel['ls']
            if 'var' in kernel:
                var = kernel['var']
            if 'period' in kernel:
                period = kernel['period']

        if kernelName == 'RBF':
            k = gpflow.ekernels.RBF(self.Q, lengthscales=ls, variance=var, ARD=True)
        elif kernelName == 'Matern32':
            k = gpflow.kernels.Matern32(self.Q, lengthscales=ls, variance=var)
            # k =  k + gpflow.kernels.White(input_dim, variance=0.01)
        elif kernelName == 'Periodic':
            k = gpflow.kernels.PeriodicKernel(self.Q, period=period,
                    lengthscales=ls, variance=var)
        if periodic:
            k = Periodic (base=k)

        self.kern = k

    def set_X_prior_mean(self, X_prior_mean):
        if X_prior_mean is not None:
            if type(X_prior_mean) is str:
                self.X_prior_mean = self.mData[X_prior_mean].values[:, None]
            else:
                #self.X_prior_mean = X_prior_mean
                prior_mean = (X_prior_mean - X_prior_mean.mean())/X_prior_mean.std()
                self.X_prior_mean = prior_mean
        else:
            self.X_prior_mean = np.zeros((self.N, self.Q))

    def set_X_prior_var(self, X_prior_var):
        if X_prior_var is not None:
            if type(X_prior_var) is str:
                self.X_prior_var = self.mData[X_prior_var].values[:, None]
            elif type(X_prior_var) is np.ndarray:
                self.X_prior_var = X_prior_var
            else:
                self.X_prior_var = X_prior_var * np.ones((self.N, self.Q))
        else:
            self.X_prior_var = 1. * np.ones((self.N, self.Q))

    def set_X_mean(self, X_mean, X_prior_mean=None):
        if X_mean is not None:
            if type(X_mean) is str:
                self.X_mean = self.mData[X_mean].values[:, None]
            else:
                self.X_mean = X_mean
        else:
            self.X_mean = MapTo01(gpflow.models.PCA_reduce(self.Y, self.Q))
            # I added the 2 lines below
            if X_prior_mean is not None:
                self.X_mean [:,-1] = ((X_prior_mean - X_prior_mean.mean(
                    keepdims=True))/X_prior_mean.std(keepdims=True))[:,0]

    def set_X_var(self, X_var):
        if type(X_var) is str:
            self.X_var = self.mData[X_var].values[:, None]
        elif type(X_var) is np.ndarray:
            self.X_var = X_var
        else:
            self.X_var = X_var * np.ones((self.N, self.Q))

    def set_inducing_inputs(self, Z):
        if Z is None:
            self.Z = np.random.permutation(self.X_mean.copy())[:self.M]
        else:
            self.Z = np.asarray(Z)

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
        maxitr = 1000,
        return_pred = False,
        pred_value=None,
        pred_noise=True,
        periodic=False):
    sys.argv.append(dtype)

    # zero the mean, set the standard deviation to 1
    m = GrandPrixModel(data, n_latent_dims, n_inducing_points,
            kernel, mData, latent_prior_mean, latent_prior_var, latent_mean,
            latent_var, inducing_inputs, dtype, periodic=periodic)
    m.set_jitter_level(jitter)
    m.set_trainable(fix_parameters)
    m.build()

    m.fit(maxiter=maxitr, display=disp)
    print (print_kernel_summary (m.m))
    if return_pred:
        if pred_value is None:
            posterior = m.get_latent_dims()
        else: posterior = [pred_value]
        if pred_noise:
            pt_pred_out, var_pred_out = m.m.predict_y (posterior[0])
        else:
            pt_pred_out, var_pred_out = m.m.predict_f (posterior[0])
        return m.get_latent_dims(), (pt_pred_out, var_pred_out)
    else:
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

