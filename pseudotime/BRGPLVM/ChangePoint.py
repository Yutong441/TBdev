from collections.abc import Iterable
import tensorflow as tf
import gpflow
from gpflow.params import Parameter
from typing import List, Optional, Union
from gpflow.decors import params_as_tensors

float_type = gpflow.settings.tf_float
class tanh_CP(gpflow.kernels.Combination):
    '''
    from https://gpflow.readthedocs.io/en/master/_modules/gpflow/kernels/changepoints.html
    except that locations should be gpflow parameters
    '''
    def __init__(
        self, kernels, locations, steepness,
        constraints = None,
        #steepness: Union[float, List[float]] = 0.01,
        name: Optional[str] = None,
    ):
        if len(kernels) != locations._value.shape[0] + 1:
            raise ValueError(
                "Number of kernels ({nk}) must be one more than the number of "
                "changepoint locations ({nl})".format(nk=len(kernels), nl=len(locations))
            )

        if isinstance(steepness, Iterable) and len(steepness) != len(locations):
            raise ValueError(
                "Dimension of steepness ({ns}) does not match number of changepoint "
                "locations ({nl})".format(ns=len(steepness), nl=len(locations))
            )

        super().__init__(kernels, name=name)

        self.locations = locations
        self._set_constraints (constraints)
        #self.steepness = Parameter(steepness, transform=gpflow.transforms.positive)
        self.steepness = steepness

    def _set_constraints (self, constraint):
        if constraint is None: 
            N = len (self.kern_list) - 1
            self.diff = tf.zeros ([N,], dtype=float_type)
        else:
            if len(self.kern_list) != constraint._value.shape[0] + 1:
                raise ValueError(
                    "Number of kernels ({nk}) must be one more than the number of "
                    "constraints locations ({nl})".format(nk=len(kernels), nl=len(constraint))
                )
            self.diff = constraint

    def _set_kernels(self, kernels):
        # it is not clear how to flatten out nested change-points
        self.kernels = kernels

    @params_as_tensors
    def K(self, X, X2 = None, presliced=False):
        sig_X = self._sigmoids(X)  # N1 x 1 x Ncp
        sig_X2 = self._sigmoids(X2) if X2 is not None else sig_X  # N2 x 1 x Ncp

        # `starters` are the sigmoids going from 0 -> 1, whilst `stoppers` go
        # from 1 -> 0, dimensions are N1 x N2 x Ncp
        starters = sig_X * tf.transpose(sig_X2, perm=(1, 0, 2))
        stoppers = (1 - sig_X) * tf.transpose((1 - sig_X2), perm=(1, 0, 2))

        # prepend `starters` with ones and append ones to `stoppers` since the
        # first kernel has no start and the last kernel has no end
        N1 = tf.shape(X)[0]
        N2 = tf.shape(X2)[0] if X2 is not None else N1
        ones = tf.ones((N1, N2, 1), dtype=X.dtype)
        starters = tf.concat([ones, starters], axis=2)
        stoppers = tf.concat([stoppers, ones], axis=2)

        # now combine with the underlying kernels
        kernel_stack = tf.stack([k.K(X, X2) for k in self.kern_list], axis=2)
        return tf.reduce_sum(kernel_stack * starters * stoppers, axis=2)

    @params_as_tensors
    def Kdiag(self, X, presliced=False):
        N = tf.shape(X)[0]
        sig_X = tf.reshape(self._sigmoids(X), (N, -1))  # N x Ncp

        ones = tf.ones((N, 1), dtype=X.dtype)
        starters = tf.concat([ones, sig_X * sig_X], axis=1)  # N x Ncp
        stoppers = tf.concat([(1 - sig_X) * (1 - sig_X), ones], axis=1)

        kernel_stack = tf.stack([k.Kdiag(X) for k in self.kern_list], axis=1)
        return tf.reduce_sum(kernel_stack * starters * stoppers, axis=1)

    @params_as_tensors
    def _sigmoids(self, X: tf.Tensor):
        locations = self.locations + self.diff
        locations = tf.sort(locations)  # ensure locations are ordered
        locations = tf.reshape(locations, (1, 1, -1))
        steepness = tf.reshape(self.steepness, (1, 1, -1))
        #return (1 - tf.tanh((locations -X[:, :, None])/steepness)) /2
        return tf.sigmoid ((X[:, :, None] - locations )*steepness) 
