import numpy as np
import tensorflow as tf
import gpflow
from gpflow.params import Parameter
from gpflow.decors import params_as_tensors, autoflow
import BRGP_utils as bu

# need a high jitter to avoid failures of Cholesky decomposition
jitter = 1e-3
# can set to tf.float32 for faster execution, but I get the errors with
# Cholesky decomposition more frequently
float_type = gpflow.settings.tf_float

def PCA_reduce(X, Q):
    """
    A helpful function for linearly reducing the dimensionality of the data X
    to Q.
    :param X: data array of size N (number of points) x D (dimensions)
    :param Q: Number of latent dimensions, Q < D
    :return: PCA projection array of size N x Q.
    """
    assert Q <= X.shape[1], 'Cannot have more latent dimensions than observed'
    evecs, evals = np.linalg.eigh(np.cov(X.T))
    i = np.argsort(evecs)[::-1]
    W = evals[:, i]
    W = W[:, :Q]
    return (X - X.mean(0)).dot(W)

class zero_kernel (gpflow.kernels.Static):
    '''K0 as described in Penfold (2018). It will always output 0.'''
    def K(self, X, X2=None, presliced=False):
        if X2 is None: X2= X
        return tf.zeros ([X.shape[0], X2.shape[0]], dtype=float_type)

class BRGP (gpflow.models.BayesianGPLVM):
    '''
    This script modifies the source code of gpflow for `BayesianGPLVM`, including
    the branch assignment labels in the kernel calculations. The variable names
    and those used in this documentation follow the convention from Penfold
    (2018).
    It is important to run on GPflow2.1.3, not 2.1.4 in which the `len()` does not
    seem to work on inducing variable class.
    I cannot guarantee this script is bug free entirely, but I can at least in
    the compilation stage. The most likely source of error is that Cholesky
    decomposition is not successful.
    Reference: 
    Titsias 2010: http://proceedings.mlr.press/v5/titsias09a/titsias09a.pdf
    Penfold 2018: https://academic.oup.com/bioinformatics/article/34/17/i1005/5093256
    '''
    def __init__ (self, data, kernel_list, BP,
            X_data_mean=None, X_data_var=None, 
            X_prior_mean=None, X_prior_var=None,
            num_inducing_variables=30, 
            tie_index=None, scale_data=False):
        '''
        Most of the arguments here are the same as those for `BayesianGPLVM`,
        except the following two:
        Args:
            `BP`: prior branch assignment probability, should be the same number
            of columns as the number of branches, the same number of rows as
            `X_data_mean`. It should be one hot encoded.
            `kernel_list`: a list of kernels. The first one must be the base
            kernel
            `tie_index`: which samples in the `BP` are tied with the last
            indices of BP, then the parameters in pseudotime and observed data
            are tied accordingly

        Attributes:
            `kernel`: a list of kernels from `kernel_list`
            `num_kernel`: number of kernels used, equal to the number of
            branches plus 1 (base kernel)
            `inducing_variable`: same as `BayesianGPLVM`
            `BP_Z`: branch assignment for the inducing variables
        '''
        if X_prior_mean is None:
            X_prior_mean = np.zeros ([data.shape[0], 1])
        else: X_prior_mean = (X_prior_mean - X_prior_mean.mean())/X_prior_mean.std()

        if X_prior_var is None:
            X_prior_var = np.ones ([data.shape[0], 1])
        if X_data_mean is None:
            #X_data_mean = bu.MapTo01(PCA_reduce(data, 1))
            X_data_mean = X_prior_mean
        if X_data_var is None:
            X_data_var  = np.ones ([data.shape[0], 1])/10.

        N = X_data_mean.shape[0]
        M = len (kernel_list) 
        assert N == data.shape[0], \
                'length of `X_data_mean` should be the same as `data`'

        # if tieing different indices is required
        if tie_index is not None:
            data, = bu.duplicate_main ([data], tie_index)
        # scale the data after duplicating the data points
        if scale_data:
            data = (data - data.mean (0, keepdims=True) )/data.std (0, keepdims=True)

        # return data to the original size because it cannot differ from X_mean
        super (BRGP, self).__init__ (Y=data[:N], X_mean=X_data_mean, X_var=X_data_var,
                kern =kernel_list[1], X_prior_mean=X_prior_mean,
                X_prior_var=X_prior_var, M=num_inducing_variables)
         
        self.kernel = kernel_list
        self.num_kernel = len (kernel_list)

        if tie_index is not None:
            assert N == len (tie_index),  \
                '`num_tie` should have the same length as the number of samples'
            num_tie = tie_index.astype (int).sum()
            assert BP.shape[0] == N + num_tie, \
                    'the number of cell labels should be the same as the number of samples'
            self.X_prior_mean, self.X_prior_var, X_data_mean_Z = \
                    bu.duplicate_main ( [self.X_prior_mean, self.X_prior_var,
                        X_data_mean], tie_index)
            self.tie_index = list (np.arange (N) [ tie_index ])
        else:
            assert BP.shape[0] == N, \
                    'the number of cell labels should be the same as the number of samples'
            self.tie_index = None
            X_data_mean_Z = X_data_mean

        # set up the branch assignment
        self.training_epochs = 0
        # check the shape of the branch matrix
        assert BP.shape[1] == M, \
                'the number of branches should match the kernel number'
        self.BP = BP

        # reset the inducing variables
        rand_index = list(np.random.permutation(len (X_data_mean_Z) ) )
        Z = X_data_mean_Z [rand_index[:num_inducing_variables]]
        self.Z = Parameter (Z)
        self.Z.trainable= False
        # set the assignment probability for the inducing variables
        self.BP_Z = BP [rand_index[:num_inducing_variables]]

    def p_X_p_2D (self, x_tensor, p_n, p_m):
        return tf.reduce_sum (x_tensor*p_n*p_m, axis=1 )

    def p_X_p_3D (self, x_tensor, p_n, p_m) :
        select_mat = tf.expand_dims(p_n, axis=1)*tf.expand_dims(p_m, axis=0) 
        return tf.reduce_sum (x_tensor*select_mat, axis=2 )

    def p_X_p_4D (self, x_tensor, p_n, p_m):
        x_return = x_tensor*p_n [:, tf.newaxis, tf.newaxis, :]
        x_return *= p_m [tf.newaxis, :, tf.newaxis, :]
        x_return *= p_m [tf.newaxis, tf.newaxis, :, :]
        return tf.reduce_sum (x_return, axis=3 )

    def p_K_p (self, X1, p_n, full_cov=True) :
        '''Evaluate K (t, t, z, z)'''
        # stack the kernels along the last axis
        if full_cov:
            K_list = tf.stack ([ self.kernel[i].K (X1) for i in
                range (self.num_kernel)], axis=-1)
        else:
            K_list = tf.stack ([ self.kernel[i].Kdiag (X1) for i in
                range (self.num_kernel)], axis=-1)

        # when full_cov is False, the return tensor is only 2D not 3D
        # need to expand the dimension for consistency
        if len (K_list.shape) != 2:
            return self.p_X_p_3D (K_list, p_n, p_n) 
        else: return tf.reduce_sum (K_list*p_n**2, axis=1)

    @params_as_tensors
    def get_psi_0 (self, pX, pX_var, p_n):
        '''Args: `pX`: the variable to be integrated out'''
        psi0 = tf.stack ([ self.kernel[i].eKdiag (pX, pX_var) for i in range
                (self.num_kernel) ], axis=-1)
        return self.p_X_p_2D (psi0, p_n, p_n) 

    @params_as_tensors
    def get_psi_1 (self, pX, pX_var, pX_prime, p_n, p_m=None):
        '''
        Args: 
            `pX`: the variable to be integrated out
            `pX_prime`: the variable not to be integrated out
        '''
        psi1 = tf.stack ([self.kernel[i].eKxz(pX_prime, pX, pX_var) for i in
            range (self.num_kernel) ], axis=-1)
        return self.p_X_p_3D (psi1, p_n, p_m) 

    @params_as_tensors
    def get_psi_2 (self, pX, pX_var, pX_prime, p_n, p_m=None):
        psi2 = tf.stack ( [self.kernel[i].eKzxKxz(pX_prime, pX, pX_var
            ) for i in range (self.num_kernel)], axis=-1)
        psi2_p = self.p_X_p_4D (psi2, p_n, p_m)
        return tf.reduce_sum (psi2_p, axis=0) 

    @params_as_tensors
    def p_Kuu_p (self, pX, p_n):
        L = tf.stack ([self.kernel[i].K(pX) for i in range (self.num_kernel) ],
                axis=-1)
        L_p = self.p_X_p_3D (L, p_n, p_n)
        return L_p

    @params_as_tensors
    def p_Kuf_p (self, pX, pX_prime, p_n, p_m) :
        L = tf.stack ([self.kernel[i].K(pX, pX_prime) for i in
            range (self.num_kernel) ], axis=-1)
        L_p = self.p_X_p_3D (L, p_n, p_m)
        return L_p

    @params_as_tensors
    def _build_likelihood(self):
        ''' This method is based on `BayesianGPLVM.elbo` 

        NB: in gpflow1, it is necessary to add the decorator @params_as_tensors
        because the kernels can only process tensors not gpflow parameters.
        This does not appear to be necessary in gpflow 2
        '''
        self.training_epochs += 1
        Y_data = tf.cast(self.Y, dtype=float_type)
        BP = tf.cast (self.BP, float_type)
        BP_Z = tf.cast (self.BP_Z, float_type)

        if self.tie_index is not None:
            Y_data, X_mean, X_var = bu.duplicate_main_tf ([Y_data, self.X_mean,
                self.X_var], self.tie_index)
        else:
            X_mean, X_var = self.X_mean, self.X_var

        num_inducing = tf.shape(self.Z)[0]
        psi0 = self.get_psi_0 (X_mean, X_var, BP)
        psi1 = self.get_psi_1 (X_mean, X_var, self.Z, BP, BP_Z)
        psi2 = self.get_psi_2 (X_mean, X_var, self.Z, BP, BP_Z)
        L = self.p_Kuu_p ( self.Z, BP_Z )
        L += jitter*tf.eye (num_inducing, dtype=float_type )
        L = tf.cholesky (L)

        sigma2 = self.likelihood.variance
        sigma = tf.sqrt(sigma2)

        # The below code is almost the same as source code
        # Compute intermediate matrices
        A = tf.matrix_triangular_solve(L, tf.transpose(psi1), lower=True) / sigma
        tmp = tf.matrix_triangular_solve(L, psi2, lower=True)
        AAT = tf.matrix_triangular_solve(L, tf.transpose(tmp), lower=True) / sigma2
        B = AAT + tf.eye(num_inducing, dtype=float_type)
        LB = tf.cholesky(B )
        log_det_B = 2.0 * tf.reduce_sum(tf.log(tf.matrix_diag_part(LB)))
        c = tf.matrix_triangular_solve(LB, tf.matmul(A, Y_data), lower=True) / sigma

        # KL[q(x) || p(x)]
        x_var_shape_dim = len (X_var.get_shape())
        dX_data_var = (
            X_var
            if x_var_shape_dim == 2
            else tf.matrix_diag_part(X_var)
        )
        NQ = tf.cast(tf.size(X_mean), dtype=float_type)
        D = tf.cast(tf.shape(Y_data)[1], dtype=float_type)
        KL = -0.5 * tf.reduce_sum(tf.log(dX_data_var))
        KL += 0.5 * tf.reduce_sum(tf.log(self.X_prior_var))
        KL -= 0.5 * NQ
        KL += 0.5 * tf.reduce_sum(
            (tf.square(X_mean - self.X_prior_mean) + dX_data_var) / self.X_prior_var
        )
        # compute log marginal bound
        ND = tf.cast(tf.size(Y_data), dtype=float_type)
        bound = -0.5 * ND * tf.log(2 * np.pi * sigma2)
        bound += -0.5 * D * log_det_B
        bound += -0.5 * tf.reduce_sum(tf.square(Y_data)) / sigma2
        bound += 0.5 * tf.reduce_sum(tf.square(c))
        bound += -0.5 * D * (tf.reduce_sum(psi0) / sigma2 - tf.reduce_sum(tf.matrix_diag_part(AAT)))
        bound -= KL
        return bound

    @params_as_tensors
    def _build_predict(self, Xnew, BP_new, full_cov=False):
        '''
        This is based on `BayesianGPLVM.predict_f`
        Reference: sparse Gaussian process regression
        http://www.gatsby.ucl.ac.uk/~snelson/SPGP_up.pdf
        Derivation: https://gpflow.readthedocs.io/en/master/notebooks/theory/SGPR_notes.html

        Args:
            `Xnew`: a 1D numpy array
            `branch_num`: predict the Y values given X for one particular
            branch. To use the base kernel, set to 0

        Notations in comments: x is the pseudotime estimated for the original
        data, u is the inducing point, f is the predicted point
        '''
        Y_data = tf.cast (self.Y, dtype=float_type)
        if self.tie_index is not None:
            Y_data, X_mean, X_var = bu.duplicate_main_tf ([Y_data, self.X_mean,
                self.X_var], self.tie_index)
        else:
            X_mean, X_var = self.X_mean, self.X_var
        num_inducing = tf.shape(self.Z)[0]
        psi1 = self.get_psi_1 (X_mean, X_var, self.Z, self.BP, self.BP_Z)
        psi2 = self.get_psi_2 (X_mean, X_var, self.Z, self.BP, self.BP_Z)

        Kus = self.p_Kuf_p ( self.Z, Xnew, self.BP_Z, BP_new)
        L = self.p_Kuu_p ( self.Z, self.BP_Z )
        L += jitter*tf.eye (num_inducing, dtype=float_type)
        L = tf.linalg.cholesky(L)
        sigma2 = self.likelihood.variance
        sigma = tf.sqrt(sigma2)

        A = tf.matrix_triangular_solve(L, tf.transpose(psi1), lower=True) / sigma
        tmp = tf.matrix_triangular_solve(L, psi2, lower=True)
        AAT = tf.matrix_triangular_solve(L, tf.transpose(tmp), lower=True) / sigma2
        B = AAT + tf.eye(num_inducing, dtype=float_type)

        LB = tf.cholesky(B)
        c = tf.matrix_triangular_solve(LB, tf.matmul(A, Y_data), lower=True) / sigma
        tmp1 = tf.matrix_triangular_solve(L, Kus, lower=True)
        tmp2 = tf.matrix_triangular_solve(LB, tmp1, lower=True)
        mean = tf.matmul(tmp2, c, transpose_a=True)

        if full_cov:
            var = (
                self.p_K_p(Xnew, BP_new)
                + tf.matmul(tmp2, tmp2, transpose_a=True)
                - tf.matmul(tmp1, tmp1, transpose_a=True)
            )
            shape = tf.stack([1, 1, tf.shape(Y_data)[1]])
            var = tf.tile(tf.expand_dims(var, 2), shape)
        else:
            # calculate K_{ff} - k_{fu}^T K_{uu}^{-1} K_{uu}^{-1} k_{fu} 
            # - k_{fu}^T K_{uu}^{-1} B^{-1} K_{uu}^{-1} k_{fu} 
            var = (
                self.p_K_p(Xnew, BP_new, full_cov=False)
                + tf.reduce_sum(tmp2**2, axis=0)
                - tf.reduce_sum(tmp1**2, axis=0)
            )
            shape = tf.stack([1, tf.shape(Y_data)[1]])
            var = tf.tile(tf.expand_dims(var, 1), shape)
        return mean + self.mean_function(Xnew), var

    @autoflow((float_type, [None, None]), (float_type, [None, None]))
    def predict_latent (self, Xnew, BP_new):
        return self._build_predict(Xnew, BP_new)

    @autoflow((float_type, [None, None]),(float_type, [None, None]) )
    def predict_BP (self, Xnew, BP_new):
        pred_f_mean, pred_f_var = self._build_predict(Xnew, BP_new)
        return self.likelihood.predict_mean_and_var(pred_f_mean, pred_f_var)

    @autoflow ()
    def elbo (self):
        return self._build_likelihood ()
