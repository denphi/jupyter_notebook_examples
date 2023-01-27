import warnings
warnings.simplefilter(action="ignore", category=UserWarning)
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel

class GPR:

    @staticmethod
    def _check_and_regularize_1d_array(x):
        """
        Checks if `x` is a `np.ndarray`. If `x.ndim` is 1,
        then, it turns it into a column vector.
        Otherwise, it leaves it as is.
        """
        x = np.array(x)
        if x.ndim == 1:
            x = x[:, None]
        assert x.ndim == 2
        return x

    def __init__(self, x, y):
        """
        :param x:               The observed inputs (1D `np.ndarray`)
        :param y:               The observed outputs (1D `np.ndarray`).
        """
        self.x = x.to_numpy().reshape(-1, 1) 
        self.y = y.to_numpy().reshape(-1, 1) 

    def run(self, variance=1.,
            length_scale=1.,
            noise_variance=1):
        """
        Perform 1D regression.
        :param variance:        The signal strength of the square exponential
                            covariance function (positive float).
        :param length_scale:    The length scale of the square
                            exponential covariance function (positive float).
        :param noise_variance:  The noise of the model (non-negative float).
        :returns:               A dictionary containing the following elements:
                                    + x:       points on which the predictive
                                                    distribution is actually evaluated
                                                    (1D `np.ndarray`)
                                    + y:         points target on which the predictive
                                                    distribution is actually evaluated
                                                    (1D `np.ndarray`)
                                                    (1D `np.ndarray` of size x_eval.shape[0])
                                    + y_mean:        Mean of predictive distribution a query points.
                                    + y_std:        Standard deviation of predictive distribution at query points
                                    + k:          The kernel
                                    + gpr:        the trained gaussian process model
        """
        variance = float(variance)
        length_scale = float(length_scale)
        noise_variance = float(noise_variance)

        k = 1.0 * RBF(length_scale=length_scale, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(
            noise_level=noise_variance, noise_level_bounds=(1e-5, 1e1)
        )
        gpr = GaussianProcessRegressor(kernel=k, alpha=0.0)
        gpr.fit(self.x, self.y)
        y_mean, y_std = gpr.predict(self.x, return_std=True)
        
        return self.x, self.y, y_mean, y_std, gpr, k

