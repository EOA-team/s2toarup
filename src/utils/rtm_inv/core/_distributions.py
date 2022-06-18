'''
Distributions for sampling RTM parameters
'''

import numpy as np

from scipy.stats import truncnorm
from typing import List, Optional, Union


class Distributions(object):
    """
    Class with statistical distributions for drawning RTM
    samples from a set of input parameters

    For each RTM parameter, min, max and type of distribution
    must be passed
    """

    distributions: List[str] = ['Gaussian', 'Uniform']

    def __init__(
            self,
            min_value: Union[int,float],
            max_value: Union[int,float],
            mean_value: Optional[Union[int,float]] = None,
            std_value: Optional[Union[int,float]] = None
            
        ):
        """
        Creates a new ``Distributions`` class to use for sampling
        RTM parameter values

        :param min_value:
            minimum parameter value
        :param max_value:
            maximum parameter value
        :param mean_value:
            optional mean value to be used for creating a Gaussian
            distribution. If not provided the mean is calculated as
            min_value + 0.5 * (max_value - min_value)
        :param std_value:
            optional standard deviation value to be used for creating
            a Gaussian distribution. If not provided the standard
            deviation is calculated as 0.5 * (max_value - min_value)
        """
        if min_value > max_value:
            raise ValueError('Minimum cannot be greater than maximum')
        if mean_value is None:
            mean_value = min_value + 0.5 * (max_value - min_value)
        if std_value is None:
            std_value = 0.5 * (max_value - min_value)
        self.min_value = min_value
        self.max_value = max_value
        self.mean_value = mean_value
        self.std_value = std_value

    def sample(self,
               distribution: str,
               n_samples: int
            ):
        """
        Returns a ``numpy.ndarray`` with RTM parameter samples drawn from
        a specific distribution

        :param distribution:
            name of the distribution from which to sample. See
            `~core.distributions.Distributions.distributions` for a list
            of distributions currently implemented
        :param n_samples:
            number of samples to draw.
        """
        if distribution == 'Uniform':
            return np.random.uniform(
                low=self.min_value,
                high=self.max_value,
                size=n_samples
            )
        elif distribution == 'Gaussian':
            pass
            
