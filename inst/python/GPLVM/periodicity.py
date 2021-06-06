import tensorflow as tf
import warnings
from gpflow.params import Parameter
from gpflow.decors import params_as_tensors
from gpflow.kernels import Kernel
from gpflow import transforms, settings

class Periodic(Kernel):
    """
    The periodic family of kernels. Can be used to wrap any Stationary kernel
    to transform it into a periodic version. The canonical form (based on the
    SquaredExponential kernel) can be found in Equation (47) of

    D.J.C.MacKay. Introduction to Gaussian processes. In C.M.Bishop, editor,
    Neural Networks and Machine Learning, pages 133--165. Springer, 1998.

    The derivation can be achieved by mapping the original inputs through the
    transformation u = (cos(x), sin(x)).

    For the SquaredExponential base kernel, the result can be expressed as:
        k(r) = σ² exp{ -0.5 sin²(π r / γ) / ℓ² }

    where:
    r is the Euclidean distance between the input points,
    ℓ is the lengthscale parameter,
    σ² is the variance parameter,
    γ is the period parameter.

    (note that usually we have a factor of 4 instead of 0.5 in front but this
    is absorbed into lengthscale hyperparameter).
    """
    def __init__(self, input_dim=None, period=1.0, variance=1.0,
                 lengthscales=1.0, base=None, active_dims=None, name=None):
        """
        The Periodic kernel supports the specification of any Stationary kernel
        as the base through the "base" parameter. If this is provided then
        input_dim, variance, lengthscales and active_dims parameters are
        taken from there.

        If a "base" kernel is not provided, then a SquaredExponential is
        constructed implicitly using the input_dim, variance, lengthscales
        and active_dims parameters.
        """
        if (base is not None) and (input_dim is not None):
            raise ValueError("input_dim should be defined through the base kernel, not explicitly")

        if base is None:
            warnings.warn(
                'Implicit specification of the base kernel for Periodic is '
                'deprecated, use Periodic(base=SquaredExponential({i}, '
                'variance={v}, lengthscales={l}, active_dims={a}), period={p}) '
                'instead.'.format(i=input_dim, v=variance, l=lengthscales,
                                  a=active_dims, p=period),
                DeprecationWarning)
            base = SquaredExponential(
                input_dim,
                variance=variance,
                lengthscales=lengthscales,
                active_dims=active_dims)

        if not isinstance(base, Stationary):
            raise TypeError("Periodic requires a Stationary kernel as the `base`")

        super().__init__(base.input_dim, base.active_dims, name=name)
        self.base = base
        self.period = Parameter(  # No ARD support for period yet
            period, transform=transforms.positive, dtype=settings.float_type)

    @property
    def ARD(self):
        return self.base.ARD

    @params_as_tensors
    def Kdiag(self, X, presliced=False):
        return self.base.Kdiag(X)

    @params_as_tensors
    def K(self, X, X2=None, presliced=False):
        if not presliced:
            X, X2 = self._slice(X, X2)
        if X2 is None:
            X2 = X

        # Introduce dummy dimension so we can use broadcasting
        f = tf.expand_dims(X, -2)  # ... x N x 1 x D
        f2 = tf.expand_dims(X2, -3)  # ... x 1 x M x D

        r = np.pi * (f - f2) / self.period
        scaled_sine = tf.sin(r) / self.base.lengthscales
        try:
            sine_r = tf.reduce_sum(tf.abs(scaled_sine), -1)
            K = self.base.K_r(sine_r)
        except NotImplementedError:
            sine_r2 = tf.reduce_sum(tf.square(scaled_sine), -1)
            K = self.base.K_r2(sine_r2)
        return K
