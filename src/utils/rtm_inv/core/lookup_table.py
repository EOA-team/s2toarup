'''
Module to create lookup-tables (LUT) of synthetic spectra
'''

import lhsmdu
import numpy as np
import pandas as pd

from pathlib import Path
from typing import List, Optional, Union

from rtm_inv._distributions import Distributions

sampling_methods: List[str] = ['LHS']

class LookupTable(object):
    """
    Lookup-table with RTM simulated spectra plus corresponding
    parameterization (leaf and canopy traits)

    :attrib samples:
        RTM trait samples generated using a custom sample strategy
        sampling. RTM-generated spectra are appended as additional
        columns.
    """
    def __init__(
            self,
            params: Union[Path,pd.DataFrame]
        ):
        """
        creates a new ``Lookup Table`` instance

        :param params:
            csv file with RTM parameters (traits), their min and max
            value and selected distribution
        """
        if isinstance(params, Path):
            self._params_df = pd.read_csv(self.params_csv)
        elif isinstance(params, pd.DataFrame):
            self._params_df = params.copy()
        else:
            raise TypeError('Expected Path-object or DataFrame')
        self.samples = None

    @property
    def samples(self) -> Union[pd.DataFrame, None]:
        """
        Trait samples for generating synthetic spectra
        """
        return self._samples

    @samples.setter
    def samples(self, in_df: pd.DataFrame) -> None:
        """
        Trait samples for generating synthetic spectra
        """
        if in_df is not None:
            if not isinstance(in_df, pd.DataFrame):
                raise TypeError(
                    f'Expected a pandas DataFrame instance, got {type(in_df)}'
                )
        self._samples = in_df

    def generate_samples(
            self,
            num_samples: int,
            method: str,
            seed_value: Optional[int] = 0
        ):
        """
        Sample parameter values using a custom sampling scheme.

        Currently supported sampling schemes are:

        - Latin Hypercube Sampling (LHS)
        - ...

        All parameters (traits) are sampled, whose distribution is not set
        as "constant"

        :param num_samples:
            number of samples to draw (equals the size of the resulting
            lookup-table)
        :param method:
            sampling method to apply
        :param seed_value:
            seed value to set to the pseudo-random-number generator. Default
            is zero.
        """
        # set seed to the random number generator
        np.random.seed(seed_value)

        # determine traits to sample (indicated by a distribution different from
        # "Constant"
        traits = self._params_df[
            self._params_df['Distribution'].isin(Distributions.distributions)
        ]
        trait_names = traits['Parameter'].to_list()
        traits = traits.transpose()
        traits.columns = trait_names
        n_traits = len(trait_names)

        # and those traits/ parameters remaining constant
        constant_traits = self._params_df[
            ~self._params_df['Parameter'].isin(trait_names)
        ]
        constant_trait_names = constant_traits['Parameter'].to_list()
        constant_traits = constant_traits.transpose()
        constant_traits.columns = constant_trait_names

        # select method and conduct sampling
        if method.upper() == 'LHS':
            # create LHS matrix
            lhc = lhsmdu.createRandomStandardUniformMatrix(n_traits, num_samples)
            traits_lhc = pd.DataFrame(lhc).transpose()
            traits_lhc.columns = trait_names
            # replace original values in LHS matrix (scaled between 0 and 1) with
            # trait values scaled between their specific min and max
            for trait_name in trait_names:
                traits_lhc[trait_name] = traits_lhc[trait_name] * \
                    traits[trait_name]['Max'] + traits[trait_name]['Min']
        else:
            raise NotImplementedError(f'{method} not found')

        # combine trait samples and constant values into a single DataFrame
        # so that in can be passed to the RTM
        for constant_trait in constant_trait_names:
            # for constant traits the value in the min column is used
            # (actually min and max should be set to the same value)
            traits_lhc[constant_trait] = constant_traits[constant_trait]['Min']

        # set samples to instance variable
        self.samples = traits_lhc
