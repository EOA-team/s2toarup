'''
Fitting of a double logistic curve following the approach proposed by Vrieling et al.
for determining the free curve parameters
'''

import numpy as np
from copy import deepcopy
from scipy.interpolate import interp1d
from typing import Optional
from scipy.signal import savgol_filter

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WHITTAKER-EILERS SMOOTHER in Python 3 using numpy and scipy

based on the work by Eilers [1].
    [1] P. H. C. Eilers, "A perfect smoother", 
        Anal. Chem. 2003, (75), 3631-3636
coded by M. H. V. Werts (CNRS, France)
tested on Anaconda 64-bit (Python 3.6.4, numpy 1.14.0, scipy 1.0.0)

Read the license text at the end of this file before using this software.

Warm thanks go to Simon Bordeyne who pioneered a first (non-sparse) version
of the smoother in Python.
"""

import scipy.sparse as sparse
from scipy.sparse.linalg import splu


def speyediff(N, d, format='csc'):
    """
    (utility function)
    Construct a d-th order sparse difference matrix based on 
    an initial N x N identity matrix
    
    Final matrix (N-d) x N
    """
    
    assert not (d < 0), "d must be non negative"
    shape     = (N-d, N)
    diagonals = np.zeros(2*d + 1)
    diagonals[d] = 1.
    for i in range(d):
        diff = diagonals[:-1] - diagonals[1:]
        diagonals = diff
    offsets = np.arange(d+1)
    spmat = sparse.diags(diagonals, offsets, shape, format=format)
    return spmat


def whittaker_smooth(y, lmbd, d = 2):
    """
    Implementation of the Whittaker smoothing algorithm,
    based on the work by Eilers [1].

    [1] P. H. C. Eilers, "A perfect smoother", Anal. Chem. 2003, (75), 3631-3636
    
    The larger 'lmbd', the smoother the data.
    For smoothing of a complete data series, sampled at equal intervals

    This implementation uses sparse matrices enabling high-speed processing
    of large input vectors
    
    ---------
    
    Arguments :
    
    y       : vector containing raw data
    lmbd    : parameter for the smoothing algorithm (roughness penalty)
    d       : order of the smoothing 
    
    ---------

    Returns :
    
    z       : vector of the smoothed data.
    """

    m = len(y)
    E = sparse.eye(m, format='csc')
    D = speyediff(m, d, format='csc')
    coefmat = E + lmbd * D.conj().T.dot(D)
    z = splu(coefmat).solve(y)
    return z    


#Copyright M. H. V. Werts, 2017
#
#martinus point werts à ens-rennes point fr
#
#This software is a computer program whose purpose is to smooth noisy data.
#
#This software is governed by the CeCILL-B license under French law and
#abiding by the rules of distribution of free software.  You can  use, 
#modify and/ or redistribute the software under the terms of the CeCILL-B
#license as circulated by CEA, CNRS and INRIA at the following URL
#"http://www.cecill.info". 
#
#As a counterpart to the access to the source code and  rights to copy,
#modify and redistribute granted by the license, users are provided only
#with a limited warranty  and the software's author,  the holder of the
#economic rights,  and the successive licensors  have only  limited
#liability. 
#
#In this respect, the user's attention is drawn to the risks associated
#with loading,  using,  modifying and/or developing or reproducing the
#software by the user in light of its specific status of free software,
#that may mean  that it is complicated to manipulate,  and  that  also
#therefore means  that it is reserved for developers  and  experienced
#professionals having in-depth computer knowledge. Users are therefore
#encouraged to load and test the software's suitability as regards their
#requirements in conditions enabling the security of their systems and/or 
#data to be ensured and,  more generally, to use and operate it in the 
#same conditions as regards security. 
#
#The fact that you are presently reading this means that you have had
#knowledge of the CeCILL-B license and that you accept its terms.


class TimeSeriesModel(object):
    """
    Abstract class implementing basic properties of a time
    series model following the "upper-envelope" approach
    by Chen et al. (RSE, 2004).

    :param ts_values:
        values of the observed quantity at equally
        spaced temporal frequency (missing values filled
        with NaN)
    :param ts_inidices:
        optional indices denoting the time dimension.
        If not provided, equals the array indices of
        ts_values
    :param max_iter:
        maximum number of iterations for fitting if the termination
        criterion is not met before.
    """

    def __init__(
            self,
            ts_values: np.ndarray,
            time_indices: np.ndarray,
            max_iter: Optional[int]=100
        ):
        self._ts_modeled = ts_values
        self._time_indices = time_indices
        self._max_iter = max_iter

        # store original values
        self._ts_orig = deepcopy(ts_values)

        # weights for fitting
        self._weights = np.empty(shape=self._ts_orig.shape)

        # fitting-effect indices
        self.fitting_indices = []


    def _calculate_weights(self) -> None:
        """
        Calculates the weights between the original and the
        modeled time series using the methodology proposed by
        Chen et al. (2004).
        """
        # assess difference original vs. modeled
        diff = self._ts_orig - self._ts_modeled
        # inverse maximum absolute difference
        dmax = 1. / np.abs(diff).max()

        def calc_weight(orig: float,
                        modeled: float
                        ) -> float:
            """
            auxiliary function returning the weight per data point
            """
            if orig >= modeled:
                return 1.
            else:
                return 1. - np.abs(orig - modeled) * dmax

        calc_weights = np.vectorize(calc_weight)
        self._weights = calc_weights(self._ts_orig, self._ts_modeled)

    def _update_time_series(self) -> np.array:
        """
        Updates the time series based on Chen et al.'s upper
        envelope approach replacing noisy points in the time
        series with modeled ones while keeping the other points
        """
        def update(
                orig: float,
                modeled: float
            ) -> float:
            """
            auxiliary function performing the update per element
            """
            if orig >= modeled:
                return orig
            else:
                return modeled
        update_data = np.vectorize(update)
        return update_data(self._ts_orig, self._ts_modeled)

    def _calc_fitting_index(self) -> float:
        """
        Calculates the Fitting-Effect index after Chen et al. (2004)
        between modeled and original data and the calculated weights
        per data point
        """
        return np.sum(np.abs(self._ts_modeled - self._ts_orig) * self._weights)


    def _get_time_indices(self) -> None:
        """
        populates time indices if not available from inputs
        """
        if self._time_indices is None:
            self._time_indices = np.array(
                [x for x in range(self._ts_orig.shape[0])]
            )

    def lin_interpol(self
                     ) -> None:
        """
        linear interpolation of the original time series data
        points
        """
        # find indices of original datapoints that are not NaN
        # indices = [idx[0] for idx, val in np.ndenumerate(self._ts_orig) \
        #            if not np.isnan(val)]
        indices = np.arange(self._time_indices[0], self._time_indices[-1]+1)
        y_lin_fit = self._ts_orig
        linear_interpol = interp1d(x=indices,
                                   y=y_lin_fit,
                                   kind='linear',
                                   copy=True,
                                   assume_sorted=True
        )
        self._ts_orig = linear_interpol(self._time_indices)

    def model_time_series(self
                          ) -> np.array:
        """
        Returns the modeled time series using the estimated
        parameters between start and end of the input time
        series
        """
        res = []
        for idx in range(self._ts_modeled.shape[0]):
            res.append(self.model_time_point(t=idx))
        return np.array(res, dtype=np.float32)

class DoubleLogisticModel(TimeSeriesModel):
    """
    Class for defining, storing and fitting
    a double logistic (in this case realized as
    double hyperbolic tangent function) following the
    approach proposed in Vrieling et al. (RSE, 2018)
    """
    def __init__(self, *args, **kwargs):

        # extend __init__ of super class
        super(DoubleLogisticModel, self).__init__(*args, **kwargs)

        # model parameters to fit
        self.a0 = np.nan
        self.a1 = np.nan
        self.a2 = np.nan
        self.a4 = np.nan
        self.a5 = np.nan
        self.a6 = np.nan


    def find_half(self
                  ) -> int:
        """
        Helper function determining the array index
        that splits the time series (equally spaced)
        into two halves.
    
        Returns the index of the "mid" point.
        """
        # return int(self._ts_modeled.shape[0] * 0.5)
        return np.argmax(self._ts_orig)

    def _calc_a0(self):
        """
        Computes the minimum value in the first
        half of the time series denoted as a0.
        """
    
        # determine index splitting the array in two halves
        mid = self.find_half()
        self.a0 = np.nanmin(self._ts_modeled[0:mid])

    def _calc_a1_a4(self):
        """
        Computes the amplitude of the observed quantity
        in the green-up (a1) and senescence (a4) phase
        of the time series. The amplitude is thereby defined
        as the difference between the maximum and the minimum
        value in the first and second half of the time series,
        respectively.
        """
    
        # determine index splitting the array in two halves
        mid = self.find_half()
    
        # calculate a1 as the difference between the max and min
        # value in the first half of the time series
        max_a1 = np.nanmax(self._ts_modeled[0:mid])
        min_a1 = np.nanmin(self._ts_modeled[0:mid])
        self.a1 = max_a1 - min_a1
    
        # calculate a2 according to a1 for the second half of
        # the time series
        max_a4 = np.nanmax(self._ts_modeled[mid:])
        min_a4 = np.nanmin(self._ts_modeled[mid:])
        self.a4 = max_a4 - min_a4

    def _calc_a2_a5(self):
        """
        Computes the inflection point for the first and
        second half of the time series in days. The inflection
        point is defined as the midpoint between the start of
        the time series and the maximum observed value (a2),
        and the midpoint between the maximum observed value and
        the end of the evenly spaced time series (a5).
        """
        # find position of the maximum observed value
        max_idx = np.nanargmax(self._ts_modeled)
    
        self.a2 = int(max_idx * 0.5)
        self.a5 = int(max_idx + 0.5*(self._ts_modeled.shape[0] - max_idx))

    def _calc_a3_a6(self):
        """
        Returns the default values reported by Vrieling et al.
        (RSE, 2018) for a3 and a6 which control the slope at the
        inflection point
        """
        self.a3 = 0.02
        self.a6 = -0.02

    def estimate_parameters(self):
        """
        Computes all six parameters of the model based
        on the current time series
        """
        self._calc_a0()
        self._calc_a1_a4()
        self._calc_a2_a5()
        self._calc_a3_a6()

    def model_time_point(self,
                         t: int
                         ) -> float:
        """
        Returns the estimated value of the modeled
        time series at a given point in time t given
        as array index. To be overwritten be inheriting classes-
        """
        estimate = self.a0 + self.a1 * ((np.tanh((t - self.a2) * self.a3) + 1) * 0.5) \
                    + self.a4 * ((np.tanh((t - self.a5) * self.a6) + 1) * 0.5) - self.a4
        return estimate

    def fit_ts_model(self):
        """
        Fits the parameters of the model based on a weighting scheme
        proposed by Chen et al. (RSE, 2004) making the assumption that
        noise in the satellite data introduces a negative bias.

        The fitting process works iteratively comparing the modeled
        time series with the original time series values by keeping only
        those records for fitting from the original data that are larger
        than the modeled values and replacing smaller original values with
        modeled ones to iteratively remove outliers in negative y-direction.

        After each iteration a "fitting-effect-index" (Chen et al., 2004) is
        computed between the fitted and the original time series data. When
        the fitting index reaches a minimum, the iteration stops.
        """
        # *****************************************************
        # here, Chen et al.'s method starts
        # *****************************************************

        # initialize modeled time series
        self._ts_modeled = deepcopy(self._ts_orig)

        # STEP 2 - FITTING
        # initial model parameter guess
        self.estimate_parameters()
        # and use them to model the time series for the first time
        self._ts_modeled = self.model_time_series()

        # STEP 3 - CALCULATE WEIGHTS
        self._calculate_weights()

        self.fitting_indices.append(self._calc_fitting_index())

        # ITERATIVE PROCEDURE UNITL TERMINATION CRITERION
        for iteration in range(1,self._max_iter):

            # STEP 4 - GENERATION OF A NEW TIME SERIES
            previous_ts_modeled = self._ts_modeled
            self._ts_modeled = self._update_time_series()

            # STEP 5 - FITTING THE UPDATED TIME SERIES
            self.estimate_parameters()
            self._ts_modeled = self.model_time_series()

            # STEP 6 - CALCULATING THE FITTING-EFFECT INDEX
            self.fitting_indices.append(self._calc_fitting_index())

            # CHECK TERMINATION CRITERION
            if iteration >= 2:
                if self.fitting_indices[iteration-2] >= self.fitting_indices[iteration-1] and \
                self.fitting_indices[iteration-1] <= self.fitting_indices[iteration]:
                    self._ts_modeled = previous_ts_modeled
                    return self._ts_modeled
            else:
                if self.fitting_indices[iteration-1] < self.fitting_indices[iteration]:
                    self._ts_modeled = previous_ts_modeled
                    return self._ts_modeled

        return self._ts_modeled

if __name__ == '__main__':

    import pandas as pd
    import matplotlib.pyplot as plt


    data_fpath = '/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_TimeSeries_Analysis/GLAI/GLAI_crops.csv'
    df = pd.read_csv(data_fpath)

    # get dates
    df.date = pd.to_datetime(df.date)
    df['start_date'] = df.date.unique()[0]
    df['time_index'] = (df.date - df.start_date).dt.days
    df.drop('start_date', axis=1, inplace=True)

    # get a single pixel (using its coordinates)
    
    crop = df[df.crop_code == 5].copy()
    coords = crop.geometry.unique()
    pixel = crop[crop.geometry == coords[0]].copy()

    # convert pixel to daily values
    pixel.set_index('date', inplace=True)
    pixel = pixel.asfreq('d')
    pixel = pixel.interpolate(method='linear', axis=0).ffill().bfill()

    # Whittaker smoother
    ws = whittaker_smooth(y=pixel.GLAI.values, lmbd=1000, d=2)

    # double-logistic regression
    dlm = DoubleLogisticModel(
        ts_values=pixel.GLAI.values,
        time_indices=pixel.time_index.values
    )
    dlm.estimate_parameters()
    dlm_fitted = dlm.model_time_series()

    # savitzky golay
    sg = savgol_filter(x=pixel.GLAI.values, window_length=11, polyorder=2, mode='interp')
    
    plt.plot(pixel.index, pixel.GLAI, color='r', label='linear interpolation')
    plt.plot(pixel.index, ws, color='b', label='Whittaker')
    plt.plot(pixel.index, dlm_fitted, color='g', label='DLM')
    plt.plot(pixel.index, sg, color='y', label='Savitzky-Golay')
    plt.legend()
    plt.show()
    




