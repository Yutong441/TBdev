import numpy as np
import pandas as pd
import tensorflow as tf
import gpflow
import plotnine as pn
import time
import BRGP_utils as bu
from BRGP_class import BRGP

class BRGP_main:
    def __init__ (self, data, kernel, BP, X_mean=None, X_var=None,
            X_prior_mean=None, X_prior_var=None, num_inducing=30,
            scale_data=True, tie_index=None):
        '''
        A wrapper class for the `BRGP` class to perform Branch Recombinant
        Gaussian Process.
        Args:
            `data`: a numpy array of shape N and M, with N being the number of
            samples and M features
            `kernel`: a list of kernel with the first one being the base kernel
            `BP`: branch assignment for each sample. The first dimension should
            be N. This function accepts one hot encoded branch data or non
            one encoded version.
            `scale_data`: whether to subtract the mean of the `data` and then
            divide by the standard deviation. This normalizes the data across
            features.
            `tie_index`: which samples in the `BP` are tied with the last
            indices of BP, then the parameters in pseudotime and observed data
            are tied accordingly. It should be a 1D numpy array with the same
            length as the number of samples. Currently, it only supports tieing
            on set of variables
        Attributes:
            `model`: a BRGP instance
        '''
        unique_num= list (np.unique (BP)) 
        if len (unique_num ) == 2 and 1. in unique_num and 0. in unique_num: pass
        else: BP = (BP == np.unique (BP).reshape ([1, -1])).astype (float)

        with gpflow.defer_build ():
            self.model = BRGP (data=data, kernel_list=kernel,
                               X_data_mean=X_mean, X_data_var=X_var, 
                               X_prior_mean = X_prior_mean, X_prior_var=X_prior_var, 
                               num_inducing_variables=num_inducing,
                               BP=BP, tie_index=tie_index, scale_data=scale_data)
            self.model.likelihood.variance = 0.01

    def fit (self, maxiter=10000):
        tic = time.time ()
        opt = gpflow.train.ScipyOptimizer()
        opt.minimize (self.model, disp=True, maxiter=maxiter)
        self.fitting_time = time.time () - tic
        print ('Optimisation takes {:.2f} min'.format (self.fitting_time/60))

    def fit_adam (self, maxiter=10000):
        tic = time.time ()
        opt = tf.optimizers.Adam (0.005)
        for i in range (maxiter):
            opt.minimize (self.model.training_loss,
                          var_list=self.model.trainable_variables,)
        self.fitting_time = time.time () - tic
        print ('Optimisation takes {:.2f} min'.format (self.fitting_time/60))

    def predict (self, num=1000, min_x_val=None, max_x_val=None, feature_names=None,
            with_noise=True, min_perc=0, max_perc=100, base_kernel=None,
            kernel_pred=None):
        '''
        Args:
            `num`: number of data pointed to be predicted
            `max_x` and `min_x`: the range of pseudotimes for which their
            associated features are predicted.
            `with_noise`: whether to include noise in the prediction
            `feature_names`: the names for the features predicted
            `base_kernel`: which kernels are at higher hierarchy. It must be a
            list
            `kernel_pred`: a list of indices stating which kernel to use. If
            None, all kernels/branches will be predicted

        Attributes:
            `prediction`: a dataframe storing the predicted mean and variance
            for each branch. It has the following columns: 'x', 'y_mu',
            'y_var', 'ymin', 'ymax'
        '''
        pred_mean, pred_var  =[], []
        X_mean = self.model.X_mean.read_value()
        if min_x_val is None: 
            min_x = np.percentile (X_mean, min_perc)
        else: min_x = min_x_val
        if max_x_val is None: 
            max_x = np.percentile (X_mean, max_perc)
        else: max_x = max_x_val
        Xnew = np.linspace (min_x, max_x, num=num)

        if kernel_pred is None: kernel_pred = np.arange (self.model.num_kernel)
        for i in kernel_pred:
            print ('analysising branch '+str (i) )
            BP_new = bu.get_BP_pred (self.model, Xnew, i, base_kernel)
            if with_noise:
                mean_i, var_i = self.model.predict_BP (Xnew.reshape([-1, 1]), BP_new)
            else:
                mean_i, var_i = self.model.predict_latent (Xnew.reshape([-1, 1]), BP_new)

            if feature_names is None: feature_names = np.arange ( mean_i.shape[1] )
            mean_col = 'mean_'+ feature_names
            var_col = 'var_'+ feature_names
            mean_i = pd.DataFrame (mean_i, index=np.arange (len(Xnew)),
                    columns=mean_col)
            var_i  = pd.DataFrame (var_i, index=np.arange (len(Xnew)),
                    columns=var_col)

            # no need to do the same for var because they will be cocatenated together
            mean_i ['branch'] = 'branch'+str(i)
            mean_i ['x'] = Xnew
            pred_mean.append (mean_i)
            pred_var.append (var_i)

        pred_mean = pd.concat (pred_mean, axis=0)
        pred_var = pd.concat (pred_var, axis=0)
        self.prediction = pd.concat ([pred_mean, pred_var], axis=1)
        self.features = np.array(feature_names)

    def plot_prediction (self, feature, real_time=None):
        '''
        Plot the predicted features along the pseudotime with uncertainty. The
        actual data points are also plotted in red.
        Args:
            `real_time`: if it is None, the pseudotime estimated for the data
            points will be used.
        '''
        assert feature in self.features, 'cannot find the requested feature'
        assert hasattr (self, 'prediction'), 'need to run self.predict'
        if real_time is None: real_time = self.model.X_mean.read_value()
        feature_index = self.features == feature
        Y = self.model.Y.read_value() [:, feature_index]
        Y = pd.DataFrame (np.concatenate ([real_time, Y.reshape([-1,1])], axis=1))
        Y.columns = ['x', 'y']

        sel_column = ['mean_'+feature, 'var_'+feature, 'branch', 'x']
        plot_data = self.prediction [sel_column]
        plot_data.columns = ['y_mu', 'y_var', 'branch', 'x']
        y_std = np.sqrt (plot_data ['y_var'].values)
        plot_data ['ymin'] = plot_data ['y_mu'].values - 2*y_std
        plot_data ['ymax'] = plot_data ['y_mu'].values + 2*y_std

        Y = Y.iloc [Y['x'].values < max (plot_data['x'].values),:]
        Y = Y.iloc [Y['x'].values > min (plot_data['x'].values),:]

        return (pn.ggplot (plot_data, pn.aes (x='x', y='y_mu') )+
                #pn.geom_line () +
                pn.geom_point (pn.aes (x='x', y='y'), Y, color='red') +
                pn.geom_ribbon (pn.aes (ymin='ymin', ymax='ymax', fill='branch' ), alpha=0.7) +
                pn.theme_minimal ()+
                pn.ggtitle (feature)
                )

    def save_pt (self, save_dir):
        all_pred = np.concatenate ([self.model.X_mean.read_value (), 
            self.model.X_var.read_value () ], axis=1)
        all_pred = pd.DataFrame (all_pred, columns = ['pt_mean', 'pt_var'])
        all_pred.to_csv (save_dir)
