'''
Actual inversion strategy using a pre-computed Lookup-Table (LUT)
and an image matrix of remotely sensed spectra.
'''

import numpy as np
import pandas as pd

from numba import njit, prange
from typing import List


@njit(parallel=True, cache=True)
def inv_img(
        lut: np.ndarray,
        img: np.ndarray,
        mask: np.ndarray,
        cost_function: str,
        n_solutions: int
    ):
    """
    Lookup-table based inversion on images by minimizing a
    cost function using *n* best solutions to improve numerical
    robustness

    :param lut:
        LUT with synthetic (i.e., RTM-simulated) spectra in the
        spectral resolution of the sensor used. The shape of the LUT
        must equal (num_spectra, num_bands).
    :param img:
        image with sensor spectra. The number of spectral bands must
        match the number  of spectral bands in the LUT. The shape of
        the img must equal (num_bands, num_rows, num_columns).
    :param mask:
        mask of `img.shape[1], img.shape[2]` to skip pixels. If all
        pixels should be processed set all cells in `mask` to False.
    :param cost_function:
        cost function implementing similarity metric between sensor
        synthetic spectra. Currently implemented: 'rmse', 'mae', 'mNSE'
    :param n_solutions:
        number of best solutions to return (where cost function is
        minimal)
    :returns:
        ``np.ndarray`` of shape `(n_solutions, img_rows, img_columns)`
        where for each pixel the `n_solutions` best solutions are returned
        as row indices in the `lut`.
    """
    output_shape = (n_solutions, img.shape[1], img.shape[2])
    lut_idxs = np.zeros(shape=output_shape, dtype='int32')

    for row in prange(img.shape[1]):
        for col in range(img.shape[2]):
            # skip masked pixels
            if mask[row, col]:
                lut_idxs[:,row,col] = -1
                continue
            # get sensor spectrum (single pixel)
            image_ref = img[:,row,col]
            # image_ref_normalized = np.sum(np.abs(image_ref - (np.mean(image_ref))))
            # cost functions (from EnMap box) implemented in a
            # way Numba can handle them (cannot use keywords in numpy functions)
            delta = np.zeros(shape=(lut.shape[0],), dtype='float64')
            for idx in range(lut.shape[0]):
                if cost_function == 'rmse':
                    delta[idx] = np.sqrt(np.mean((image_ref - lut[idx,:])**2))
                elif cost_function == 'mae':
                    delta[idx] = np.sum(np.abs(image_ref - lut[idx,:]))
                elif cost_function == 'contrast_function':
                    delta[idx] = np.sum(
                        -np.log10(lut[idx,:] / image_ref) + lut[idx,:] / image_ref
                    )
            # find the smallest errors between simulated and observed spectra
            # we need the row index of the corresponding entries in the LUT
            delta_sorted = np.argsort(delta)
            lut_idxs[:,row,col] = delta_sorted[0:n_solutions]
    return lut_idxs

@njit(cache=True, parallel=True)
def _retrieve_traits(
        trait_values: np.ndarray,
        lut_idxs: np.ndarray,
    ):
    """
    Uses the indices of the best matching spectra to retrieve the
    corresponding trait values from the LUT

    :param trait_values:
        array with traits entries in the LUT
    :param lut_idxs:
        indices of the best matching entries in the LUT
    :returns:
        3-d image with trait values with shape (n_traits, nrows, ncols)
    """
    n_traits = trait_values.shape[1]
    _, rows, cols = lut_idxs.shape
    # allocate array for storing inversion results
    trait_img_shape = (n_traits, rows, cols)
    trait_img = np.zeros(trait_img_shape, dtype='float64')
    # loop over pixels and write inversion result to trait_img
    for trait in range(n_traits):
        for row in prange(rows):
            for col in range(cols):
                # continue if the pixel was masked before and has therefore no
                # solutions
                if (lut_idxs[:,row,col] == -1).all():
                    trait_img[trait,row,col] = np.nan
                    continue
                trait_img[trait,row,col] = \
                    np.median(trait_values[lut_idxs[:,row,col],trait])
    return trait_img

def retrieve_traits(
        lut: pd.DataFrame,
        lut_idxs: np.ndarray,
        traits: List[str]
    ):
    """
    Extracts traits from a lookup-table on results of `inv_img`

    :param lut:
        complete lookup-table from the RTM forward runs (i.e.,
        spectra + trait values) as ``pd.DataFrame``.
    :param lut_idxs:
        row indices in the `lut` denoting for each image pixel
        the *n* best solutions (smallest value of cost function
        between modelled and observed spectra)
    :param traits:
        name of traits to extract from the `lut`. The output
        array will have as many entries per pixel as traits.
    :param aggregation_function:
        name of the function to aggregate the *n* best solutions
        into a single final one. Calls [np.]median per default.
        Otherwise 'mean' can be passed.
    :returns:
        ``np.ndarray`` of shape `(len(traits), rows, cols)`
        where rows and columns are given by the shape of lut_idx
    """
    trait_values = lut[traits].values
    return _retrieve_traits(trait_values, lut_idxs)
