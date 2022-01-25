#!/usr/bin/python

"""
This scripts implements the error correlation in the
spectral, temporal (along-track), and spatial (across-track)
dimension among the uncertainty contributors available from the
S2 L1C-RUT runs per single uncertainty contributor to get a
user-defined (>=100, ideally) number of L1C TOA scenarios using
Monte Carlo (MC).
In short, each L1C TOA scenario generated reflects the radiometric uncertainty
taking into account that the uncertainty contributors might be
correlated in one or more out of the three dimensions (multiple
combinations as well as strength levels of correlation are possible).

NOTE: The uncertainty scenarios will only cover the selected study area
(in the script termed region of interest, ROI). The ROI MUST have a quadratic
shape and should have a number of pixels in the 10m band that can be divided by
6 (because of the 60m S2 bands).
In our example we use a region of interest of 1200 by 1200 10m pixels (12 by 12 km)
that equals 600 by 600 20m pixels and 200 by 200 60m pixels.

The information about the error correlation is taken from the
paper by Gorroño et al. (2018, https://doi.org/10.1080/22797254.2018.1471739),
Table 1. The script has been developed based on an intensive exchange of
ideas and support by Javier Gorroño (Nov. 2021) which is gratefully
acknowledged.

It took the idea of the MC framework from the Python script written by
Javier Gorroño to support their 2018 publication. The script is available here:
https://github.com/senbox-org/snap-rut/blob/master/src/test/python/s2roiunc_test.py
(link accessed latest: 7th Nov 2021)

"""

import os
import glob
import shutil
import numpy as np
import pandas as pd
import rasterio as rio
import rasterio.mask
import itertools
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional
from typing import Union
from copy import deepcopy
from datetime import datetime
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from logger import get_logger

# setup logger -> will write log file to the /../log directory
logger = get_logger('l1c-monte-carlo')


# define L1C uncertainty contributors available from L1C-RUT
l1c_unc_contributors = [
    'noise',
    'stray_sys',
    'stray_rand',
    # 'x_talk', not considered as suggested by Gorrono et al., 2018
    'ADC',
    'DS',
    'gamma',
    'diff_abs',
    'diff_temp',
    'diff_cos',
    'diff_k',
    'quant'
]

# define S2 bands to loop over
s2_bands = [
    'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'
]
# spatial resolution (meters) of each spectral band
s2_band_res = {
    'B01': 60, 'B02': 10, 'B03': 10, 'B04': 10, 'B05': 20, 'B06': 20, 'B07': 20, 'B08': 10,
    'B8A': 20, 'B09': 60, 'B10': 60, 'B11': 20, 'B12': 20
}

# rut gain factor (RUT outputs are decoded as integers between 0 and 250, 250 means 25%
# relative uncertainty). We rescale the values to fall between 0 and 0.25.
rut_gain = 0.001


def get_roi_geom(
       ulx: Union[int,float],
       uly: Union[int,float],
       min_col: int,
       max_col: int,
       min_row: int,
       max_row: int,
       spatial_res: Union[int,float]
    ) -> Polygon:
    """
    Converts the image coordinates denoting the ROI boundaries
    in terms of row and columns into a shapely Polygon with
    projected coordinates

    :return:
        polygon of the ROI in image projection (i.e., UTM
        coordinates)
    """
    ulx_roi = ulx + min_col * spatial_res
    uly_roi = uly - min_row * spatial_res

    llx_roi = ulx_roi
    lly_roi = uly - max_row * spatial_res

    lrx_roi = ulx + max_col * spatial_res
    lry_roi = lly_roi

    urx_roi = lrx_roi
    ury_roi = uly_roi

    roi_coords = [
        [ulx_roi, uly_roi],
        [urx_roi, ury_roi],
        [lrx_roi, lry_roi],
        [llx_roi, lly_roi]
    ]
    return Polygon(roi_coords)


def upsample_array(
        in_array: np.array,
        scaling_factor: int,
    ) -> np.array:
    """
    takes a 2-dimensional input array (i.e., image matrix) and splits every
    array cell (i.e., pixel) into X smaller ones having all the same value
    as the "super" cell they belong to, where X is the scaling factor (X>=1).
    This way the input image matrix gets a higher spatial resolution without
    changing any of the original pixel values.

    The value of the scaling_factor determines the spatial resolution of the output.
    If scaling_factor = 1 then the input and the output array are the same.
    If scaling_factor = 2 then the output array has a spatial resolution two times
    higher then the input (e.g. from 20 to 10 m), and so on.

    :param array_in:
        2-d array (image matrix)
    :param scaling_factor:
        factor for increasing spatial resolution. Must be greater than/ equal to 1
    :return out_array:
        upsampled array with pixel values in target spatial resolution
    """
    # check inputs
    if scaling_factor < 1:
        raise ValueError('scaling_factor must be greater/equal 1')

    # define output image matrix array bounds
    shape_out = (
        int(in_array.shape[0]*scaling_factor),
        int(in_array.shape[1]*scaling_factor)
    )
    out_array = np.zeros(shape_out, dtype=in_array.dtype)

    # increase resolution using itertools by repeating pixel values
    # scaling_factor times
    counter = 0
    for row in range(in_array.shape[0]):
        column = in_array[row, :]
        out_array[counter:counter+scaling_factor,:] = list(
            itertools.chain.from_iterable(
                itertools.repeat(x, scaling_factor) for x in column))
        counter += scaling_factor
    return out_array



class Band_Data(object):
    """
    Data structure for storing the spectral (i.e., TOA reflectance factors)
    and associated uncertainties contributors per spectral band
    """

    def __init__(self):
        self._r_toa = np.empty(0, dtype=np.uint16)
        self._unc_contrib = {}
    
    @property
    def r_toa(self):
        """
        Array with L1C TOA reflectance values
        """
        return self._r_toa

    @r_toa.setter
    def r_toa(self, arr: np.array):
        self._r_toa = arr

    @property
    def unc_contrib(self):
        """
        Dictionary with uncertainty contributors. The dict
        keys are the single uncertainty contributors and the
        values are the uncertainties obtained from L1C-RUT.
        """
        return self._unc_contrib

    @unc_contrib.setter
    def unc_contrib(self, unc_dict: Dict[str, np.array]):
        self._unc_contrib = unc_dict



def gen_rad_unc_scenarios(
        orig_dataset_path: Path,
        unc_dataset_path: Path,
        scenario_path: Path,
        template_path: Path,
        n_scenarios: int,
        roi_bounds_10m: List[float],
        correlation_table: Optional[Path] = 'uncertainty_contributors_correlation.csv',
        check_contributors_only: Optional[bool] = False
    ) -> None:
    """
    Taking the original Sentinel-2 L1C scene and the radiometric uncertainty
    derived from running the Sentinel-2 radiometric uncertainty toolbox (S2RUT)
    this function creates time n_scenarios possible "alternative" scene realities
    that can be used for further processing, using, e.g., Sen2cor.

    A region of interest (ROI) can speed things up a lot and is therefore strongly
    recommended to use instead of processing the entire image.

    :param orig_dataset_path:
        original S2 scene in .SAFE format
    :param unc_dataset_path:
        outputs of S2RUT per spectral band generated using the provided shell script
    :param scenario_path:
        directory where to store the "alternative" scene realities (scenarios)
    :param template_path:
        template of the corresponding .SAFE folder structure without the
        spectral bands (will be created using this function)
    :param n_scenarios:
        number of scenarios to generate
    :param roi_bounds_10m:
        region of interest (ROI) boundaries in image coordindates in 10m resolution.
        Expected: [col_min, col_max, row_min, row_max]
    :param correlation_table:
        file-path to the CSV with the correlation table obtained from Gorrono et al.
        2018, Table 1. This table describes the correlation of the single uncertainty
        contributors in the three domains (spatial, temporal, spectral).
    """

    # create scenario output folders in .SAFE structure
    scenario_paths = []
    for idx in range(n_scenarios):
        current_scenario_path = scenario_path.joinpath(str(idx+1))
        # copy template
        shutil.copytree(
            template_path,
            current_scenario_path,
            dirs_exist_ok=True
        )
        scenario_paths.append(current_scenario_path)

    # roi bounds in all spatial resolutions for sub-setting the data
    roi_bounds_20m = [int(x/2) for x in roi_bounds_10m]
    roi_bounds_60m = [int(x/6) for x in roi_bounds_10m]
    roi_bounds_all = {
        10: roi_bounds_10m,
        20: roi_bounds_20m,
        60: roi_bounds_60m
    }
    roi_size_all = dict.fromkeys(roi_bounds_all.keys())
    for res in roi_size_all:
        roi_size_all[res] = {
            'n_col': roi_bounds_all[res][1]-roi_bounds_all[res][0],
            'n_row': roi_bounds_all[res][3]-roi_bounds_all[res][2]
        }

    # read the band data and the associated uncertainties
    full_img_size = dict.fromkeys([10, 20, 60])
    mc_input_data = dict.fromkeys(s2_bands)
    band_files = dict.fromkeys(s2_bands)
    band_georeference_info = dict.fromkeys(s2_bands)

    logger.info(
        f'Reading TOA reflectance data and uncertainty contributors for {orig_dataset_path.name}'
    )

    for s2_band in s2_bands:

        # get band name alias (without zero, i.e, B01 -> B1)
        s2_band_alias = s2_band.replace('0','')

        # get spatial resolution and corresponding ROI bounds
        spatial_res = s2_band_res[s2_band]
        roi_bounds = roi_bounds_all[spatial_res]
        min_row, max_row = roi_bounds[2], roi_bounds[3]
        min_col, max_col = roi_bounds[0], roi_bounds[1]

        # instanciate band data struct
        band_data = Band_Data()

        # read TOA reflectance factor data
        r_toa_file = glob.glob(
            orig_dataset_path.joinpath(f'GRANULE/*/IMG_DATA/*_{s2_band}.jp2').as_posix()
        )[0]
        band_files[s2_band] = Path(r_toa_file)

        # using rasterio.mask is much faster than subsetting the array after reading
        # therefore we need to translate the ROI into a shapely polygon using
        # the upper left coordinate information
        with rio.open(r_toa_file, 'r') as src:
            ulx, uly = src.meta['transform'][2], src.meta['transform'][5]

        # generate polygon from ROI bounds for masking
        roi_geom = get_roi_geom(
            ulx=ulx,
            uly=uly,
            min_col=min_col,
            max_col=max_col,
            min_row=min_row,
            max_row=max_row,
            spatial_res=spatial_res
        )

        with rio.open(r_toa_file, 'r') as src:
            n_rows_full = src.height
            n_cols_full = src.width
            # keep the band geo-referencation and related metadata for writing
            meta = src.meta
            meta.update({
                'QUALITY': '100',
                'REVERSIBLE': 'YES'
            })
            band_georeference_info[s2_band] = meta
            # read with roi_geom as mask
            out_band, _ = rio.mask.mask(
                src, 
                [roi_geom], 
                crop=True, 
                all_touched=True
            )
        # remember the original image size and save the band data
        if full_img_size[spatial_res] is None:
            full_img_size[spatial_res] = (n_rows_full, n_cols_full)
        band_data.r_toa = out_band[0,:,:]

        # read the single uncertainty contributors
        unc_contrib_dict = dict.fromkeys(l1c_unc_contributors)
        for l1c_unc_contributor in l1c_unc_contributors:
            
            unc_contrib_file = glob.glob(
                unc_dataset_path.joinpath(
                    f'S2*_rut_{l1c_unc_contributor}_{s2_band_alias}.tif'
                ).as_posix()
            )[0]

            with rio.open(unc_contrib_file, 'r') as src:
                # read with roi_geom as mask
                out_band, _ = rio.mask.mask(
                    src, 
                    [roi_geom], 
                    crop=True, 
                    all_touched=True
                )
            unc_contrib_dict[l1c_unc_contributor] = out_band[0,:,:] * rut_gain

        band_data.unc_contrib = unc_contrib_dict

        # save the extracted reflectance and uncertainty of the current band
        mc_input_data[s2_band] = band_data

    # sample the L1C scenarios using Monte-Carlo simulations. The simulations implement the
    # correlation in the spectral, spatial and temporal domain.
    # The domains spectral, spatial and temporal are defined as
    #
    #    spectral - error correlation between the single spectral bands
    #    spatial  - error correlation along track (North-South direction)
    #    temporal - correlation across track (East-West directon)
    #
    # If a contributor is entirely uncorrelated in all three domains, the contributor
    # is sampled independently for each spectral band and pixel.
    # If a contributor is fully correlated in one or all dimensions, then one distribution
    # is used to sample from all correlated dimensions.
    # If a contributor is partly correlated then two distributions are used; one reflects
    # the independent case and one the fully correlated part. The contributions are then
    # combined based on the strength of the correlation.
    # This approach also allows to account for correlations in a single or two dimensions,
    # only.
    #
    # The knowledge about the uncertainty contributor correlation is taken from
    # Gorrono et al., 2018 (https://doi.org/10.1080/22797254.2018.1471739), Table 1.
    #

    # load the correlation table
    corr_df = pd.read_csv(correlation_table, index_col='uncertainty_contributor')
    dimensions = ['spectral', 'spatial', 'temporal']

    # assess how the single error contributors are correlated
    # fully uncorrelated contributors
    fully_uncorr_contributors = corr_df[dimensions][
        (corr_df.spectral == 0) & (corr_df.temporal == 0) & (corr_df.spatial == 0)
    ].index.tolist()

    # constant error terms
    const_error_terms = corr_df[dimensions][
        (corr_df.spectral == 999) & (corr_df.temporal == 999) & (corr_df.spatial == 999)
    ].index.tolist()

    # fully correlated contributors (correlated, i.e., correlation_coeff == 1 in all dimensions)
    fully_corr_contributors = corr_df[dimensions][
        (corr_df.spectral == 1) & (corr_df.temporal == 1) & (corr_df.spatial == 1)
    ].index.tolist()

    # partly correlated contributors with an correlation_coeff of smaller 1 the spectral
    # dimension
    partly_corr_contributors = corr_df[dimensions][
        (
            (corr_df.spectral.between(0.01,0.99)) | (corr_df.temporal.between(0.01,0.99)) | (corr_df.spatial.between(0.01,0.99))
        ) &
        (
            ((corr_df.temporal == 1) & (corr_df.spatial == 1)) | ((corr_df.spatial == 1) & (corr_df.spectral == 1)) | ((corr_df.spectral == 1) & (corr_df.temporal == 1))
        )
    ].index.tolist()
    
    # contributors correlated in some but not all dimensions
    # these contributors are treated separately because of their more complex correlation implementation
    only_temporally_corr_contributors = corr_df[dimensions][
        (corr_df.spectral == 0) & (corr_df.temporal == 1) & (corr_df.spatial == 0)
    ].index.tolist()

    # contributors correlated across the temporal and spectral domain (but no the spatial)
    spectral_temporally_corr_contributors = corr_df[dimensions][
        (corr_df.spectral == 1) & (corr_df.temporal == 1) & (corr_df.spatial == 0)
    ].index.tolist()

    # empty image matrices for writing the samples to (dtype: uint16)
    img_matrices = dict.fromkeys(full_img_size.keys())
    for res in img_matrices:
        img_matrices[res] = np.zeros(shape=full_img_size[res], dtype=np.uint16)

    # open two arrays for saving the independent (uncorrelated) and dependent (correlated)
    # contributions for later analysis (helps to understand if the uncertainty is mainly
    # caused by random or system effects)
    uncorrelated_part = dict.fromkeys(s2_bands)
    correlated_part = dict.fromkeys(s2_bands)

    # start the iteration process
    for scenario in range(n_scenarios):
        
        logger.info(
            f'Creating scenario {scenario+1}/{n_scenarios} for {orig_dataset_path.name}'
        )

        # empty arrays for storing the errors
        error_band_dict = dict.fromkeys(s2_bands)
        for s2_band in s2_bands:
            error_band_dict[s2_band] = np.zeros_like(mc_input_data[s2_band].r_toa.astype(np.float16))

            if scenario == 0:
                uncorrelated_part[s2_band] = np.zeros_like(error_band_dict[s2_band])

        ######################################################################
        #                                                                    #
        #    ================= THE SAMPLING STARTS HERE =================    #
        #                                                                    #
        ######################################################################
        
        # completely uncorrelated contributors
        for s2_band in s2_bands:

            # get spatial resolution and corresponding array size of the ROI
            cols_and_rows = roi_size_all[s2_band_res[s2_band]]
            num_row = cols_and_rows['n_row']
            num_col = cols_and_rows['n_col']

            # loop over uncorrelated contributors
            for fully_uncorr_contributor in fully_uncorr_contributors:
                uncorr_rut = mc_input_data[s2_band].unc_contrib[fully_uncorr_contributor]

                # check type of distribution
                dist_type = corr_df[corr_df.index == fully_uncorr_contributor]['distribution'].values[0]
                if dist_type == 'normal':
                    uncorr_sample = np.random.normal(
                        loc=0,
                        scale=uncorr_rut,
                        size=(num_row, num_col)
                    )
                elif dist_type == 'uniform':
                    uncorr_sample = np.random.uniform(
                        low=-uncorr_rut * np.sqrt(3),
                        high=uncorr_rut * np.sqrt(3),
                        size=(num_row, num_col)
                    )
                error_band_dict[s2_band] += uncorr_sample

                if scenario == 0:
                    uncorrelated_part[s2_band] += uncorr_sample

            # loop over constant error terms; they are simply added
            for const_error_term in const_error_terms:
                error_band_dict[s2_band] += \
                    mc_input_data[s2_band].unc_contrib[const_error_term]

                if scenario == 0:
                    uncorrelated_part[s2_band] += \
                        mc_input_data[s2_band].unc_contrib[const_error_term]

            # save uncorrelated contributors for the current S2 band
            if scenario == 0:
                fname_uncorrelated = scenario_path.joinpath(
                    f'uncorrelated_contributors_sample_{s2_band}.tif'
                )
                meta = deepcopy(band_georeference_info[s2_band])
                meta.update({'dtype': 'float32', 'driver': 'GTiff'})
                with rio.open(fname_uncorrelated, 'w+', **meta) as dst:
                    dst.write(uncorrelated_part[s2_band], 1)
            

        # fully and partly correlated contributors
        # append these to a list of arrays and concatenate them into a 1d-array
        # this way it is possible to combine spectral bands with different pixel sizes
        # and to sample along all three dimensions.
        # If the spectral dimension has a correlation coefficient smaller 1, than it is
        # necessary to sample for that dimension using a second, independent distribution
        # and combine the samples weighted according to the correlation coefficient
        corr_contributors = fully_corr_contributors + partly_corr_contributors
        for corr_contributor in corr_contributors:

            # we must consider all three dimensions
            band_unc_arr_list = []
            for s2_band in s2_bands:
                band_unc_arr_list.append(mc_input_data[s2_band].unc_contrib[corr_contributor])

            # remember the size of the original arrays so that they can be
            # reshaped into 2d afterwards
            trailing_indices = [0]
            trailing_indices.extend([x.shape[0]*x.shape[1] for x in band_unc_arr_list])

            # concatente all 2d arrays into a single 1d array, axis=None flattens the array
            corr_rut = np.concatenate(band_unc_arr_list, axis=None)

            # sample from the same normal or uniform distribution depending on the contributor
            dist_type = corr_df[corr_df.index == corr_contributor]['distribution'].values[0]
            if dist_type == 'normal':
                corr_sample = np.ones(shape=corr_rut.shape) * corr_rut
                corr_sample = np.random.normal(0, 1, 1)[0] * corr_sample
            elif dist_type == 'uniform':
                corr_sample = np.empty(shape=corr_rut.shape)
                corr_sample = np.random.uniform(-1, 1, 1)[0] * corr_rut * np.sqrt(3)

            # partly contributors with a weaker correlation in one dimension
            # In this case, it is necessary to combine the two samples for that dimensions
            # one fully correlated and one that is independent
            if corr_contributor in partly_corr_contributors:
    
                if corr_df[corr_df.index == corr_contributor]['spectral'].values[0] < 1:
                    # get the correlation coefficient (alpha)
                    alpha = corr_df[corr_df.index == corr_contributor]['spectral'].values[0]

                    # loop over the spectral bands and sample for each band
                    # independently, maintain the correlation in the spatial
                    # and temporal domain
                    indep_band_samples = []
                    for s2_band in s2_bands:
                        corr_spatial_temporal_rut = mc_input_data[s2_band].unc_contrib[corr_contributor]
                        if dist_type == 'normal':
                            corr_spatial_temporal = np.ones(shape=corr_spatial_temporal_rut.shape) * \
                                corr_spatial_temporal_rut
                            corr_spatial_temporal = np.random.normal(0, 1, 1)[0] * corr_spatial_temporal
                        elif dist_type == 'uniform':
                            corr_spatial_temporal = np.empty(shape=corr_spatial_temporal_rut.shape)
                            corr_spatial_temporal = np.random.uniform(-1, 1, 1)[0] * corr_spatial_temporal_rut * np.sqrt(3)
                        indep_band_samples.append(corr_spatial_temporal)
                else:
                    raise Exception('this correlation is not implemented!')

            # undo the flattening of the band arrays
            for idx, s2_band in enumerate(s2_bands):

                band_samples = corr_sample[trailing_indices[idx]:trailing_indices[idx]+trailing_indices[idx+1]].reshape(
                     mc_input_data[s2_band].unc_contrib[corr_contributor].shape
                )
                # add the sample to the error_band_dict if all contributors are correlated
                if scenario == 0:
                    correlated_part[s2_band] = np.zeros_like(uncorrelated_part[s2_band])
                if corr_contributor in fully_corr_contributors:
                    error_band_dict[s2_band] += band_samples
                    if scenario == 0:
                        correlated_part[s2_band] += band_samples
                # or weight it by alpha in case the spectral domain has a correlation coefficient
                # smaller 1
                elif corr_contributor in partly_corr_contributors:
                    error_band_dict[s2_band] += (1 - alpha) * band_samples + \
                        alpha * indep_band_samples[idx]
                    if scenario == 0:
                        correlated_part[s2_band] += (1 - alpha) * band_samples + \
                            alpha * indep_band_samples[idx]

        # correlation in the temporal domain only -> should be u_stray_rand
        for contributor in only_temporally_corr_contributors:

            for s2_band in s2_bands:
                temp_rut = mc_input_data[s2_band].unc_contrib[corr_contributor]
                temp_corr = np.ones(shape=temp_rut.shape) * temp_rut
                num_row, num_col = temp_rut.shape
                for row in range(num_row):
                    temp_corr[row, :] = np.random.normal(
                        0,
                        temp_corr[row, :] + 1e-9, # avoids 0 std (taken from Gorrono et al. 2018 code)
                        num_col
                    )
                # add to other contributors
                error_band_dict[s2_band] += temp_corr
                if scenario == 0:
                    correlated_part[s2_band] += temp_corr

        # implement correlation in the temporal and spatial domain
        # this requires a little tweak so that sampling accross the columns and spectral bands
        # is possible despite the different pixel sizes
        for contributor in spectral_temporally_corr_contributors:
            
            # "resample" all bands to 10m resolution, i.e., repeat the value of a 20m
            # pixel 4 times, and the value of 60m pixel 36 times. Then loop over the
            # rows and sample along the columns of the row in all bands
            resample_bands = [True if x[1] != 10 else False for x in s2_band_res.items()]
            num_row_10m, num_col_10m = mc_input_data['B02'].unc_contrib[contributor].shape

            all_rut = np.empty(shape=(num_row_10m, num_col_10m*len(s2_bands)))
            band_scaling_factors = dict.fromkeys(s2_bands)
            for idx, s2_band in enumerate(s2_bands):
                rut_array = mc_input_data[s2_band].unc_contrib[contributor]
                if resample_bands[idx]:
                    # "resample" 10 and 20m bands
                    scaling_factor = int(s2_band_res[s2_band] / 10) # 10 m is target spatial resolution
                    band_scaling_factors[s2_band] = scaling_factor
                    all_rut[:,idx*num_col_10m:(idx+1)*num_col_10m] = upsample_array(
                        in_array=rut_array,
                        scaling_factor=scaling_factor
                    )
                else:
                    # 10m bands stay as they are
                    all_rut[:,idx*num_col_10m:(idx+1)*num_col_10m] = rut_array
                    band_scaling_factors[s2_band] = 1

            # the sampling process starts here, stack the columns of a row in each band
            # and sample for that array. Afterwards the data has to be brought back into
            # the original pixel size
            unc_samples = np.empty_like(all_rut)
            for row in range(num_row_10m):
                unc_sample_row = np.ones(num_col_10m*len(s2_bands)) * all_rut[row,:]
                unc_samples[row,:] = np.random.normal(
                    loc=0, scale=1, size=1)[0] * unc_sample_row

            for idx, s2_band in enumerate(s2_bands):
                
                band_arr = unc_samples[:,idx*num_col_10m:(idx+1)*num_col_10m]
                # band has 10m pixel size -> nothing to do
                if band_scaling_factors[s2_band] == 1:
                    error_band_dict[s2_band] += band_arr
                    if scenario == 0:
                        correlated_part[s2_band] += band_arr
                # else take if n-th row and column element to obtain the original number
                # of pixels due to the pixel size
                else:
                    error_band_dict[s2_band] += band_arr[
                        0::band_scaling_factors[s2_band],0::band_scaling_factors[s2_band]
                    ]
                    if scenario == 0:
                        correlated_part[s2_band] += band_arr[
                            0::band_scaling_factors[s2_band],0::band_scaling_factors[s2_band]
                        ]

        # save correlated contributors to file
        if scenario == 0:

            for s2_band in s2_bands:
                
                fname_correlated = scenario_path.joinpath(
                    f'correlated_contributors_sample_{s2_band}.tif'
                )
                with rio.open(fname_correlated, 'w+', **meta) as dst:
                    dst.write(correlated_part[s2_band], 1)

        if check_contributors_only:
            return

        ######################################################################
        #                                                                    #
        #    =================     SAVE THE SCENARIO    =================    #
        #                                                                    #
        ######################################################################

        # combine the original L1C band data and the uncertainty (error) term
        # band_data_scenario = band_data_toa * error_data_band + band_data_toa

        for s2_band in s2_bands:

            # define output file location (bit clumbsy due to the .SAFE structure)
            band_fname = band_files[s2_band].name
            # we need to reconstruct the intermediate part of the .SAFE directory
            # to obtain the correct sub-directory for writting the band
            dataset_path = str(band_files[s2_band].parent).split(os.sep)
            dot_safe = [x for x in dataset_path if '.SAFE' in x][0]
            dataset = os.path.sep.join(dataset_path[dataset_path.index(dot_safe)::])
            current_scenario_path = scenario_path.joinpath(str(scenario+1))
            file_dst = current_scenario_path.joinpath(dataset).joinpath(band_fname).as_posix()

            # create the L1C TOA scenario
            l1c_toa_scenario = mc_input_data[s2_band].r_toa + \
                error_band_dict[s2_band] * mc_input_data[s2_band].r_toa

            # plotting for debug/ visualization

            diff = mc_input_data[s2_band].r_toa - l1c_toa_scenario
            # plot difference image for debugging and plausibility checking
            fig, axs = plt.subplots(2)
            divider = make_axes_locatable(axs[0])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            err_img = axs[0].imshow(error_band_dict[s2_band].astype(np.float32)*100, cmap='coolwarm')
            cbar = plt.colorbar(err_img, cax=cax)
            cbar.ax.get_yaxis().labelpad = 20
            cbar.set_label(r'Relative Error $\rho_{TOA}$ [%]', rotation=270)
            axs[0].set_title(s2_band)

            divider = make_axes_locatable(axs[1])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            diff_img = axs[1].imshow(diff, cmap='coolwarm')
            cbar = plt.colorbar(diff_img, cax=cax)
            cbar.ax.get_yaxis().labelpad = 20
            cbar.set_label(r'Absolute Error $\rho_{TOA}$', rotation=270)

            fname_fig_save = current_scenario_path.joinpath(f'{s2_band}_sample.png')
            fig.savefig(fname_fig_save, bbox_inches='tight')
            plt.close(fig)

            # finally, we have to insert the scenario data into the empty image
            # matrix having the full spatial extent of the original S2 scene
            spatial_res_band = s2_band_res[s2_band]
            roi_bounds = roi_bounds_all[spatial_res_band]
            min_row, max_row = roi_bounds[2], roi_bounds[3]
            min_col, max_col = roi_bounds[0], roi_bounds[1]
            img = deepcopy(img_matrices[spatial_res_band])
            img[min_row:max_row, min_col:max_col] = l1c_toa_scenario.astype(np.uint16)

            # write to file (jp2), data type must be np.uint16
            with rio.open(file_dst, 'w+', **band_georeference_info[s2_band]) as dst:
                dst.write(img, 1)

        logger.info(
            f'Created scenario {scenario+1}/{n_scenarios} for {orig_dataset_path.name}'
        )
              

def main(
        orig_datasets_dir: Path,
        unc_datasets_dir: Path,
        scenario_dir: Path,
        roi_bounds_10m: List[int],
        n_scenarios: int,
        check_contributors_only: Optional[bool] = True
    ) -> None:
    """
    main executable function taking care about generating the scenarios and
    the required .SAFE folder structure

    :param orig_datasets_dir:
        directory where the original L1C TOA .SAFE scenes can be found
    :param unc_datasets_dir:
        directory where the associated L1C radiometric uncertainty results
        can be found per contributor and band
    :param scenario_dir:
        directory where to create and store the resulting L1C TOA scenarios
    :param roi_bounds_10m:
        because of the high computational load, the scenarios are created
        for a spatial subset of the original scene, only
    :param n_scenarios:
        number of scenarios to create (>=100).
    """

    # find scenes and their uncertainty
    orig_datasets = glob.glob(orig_datasets_dir.joinpath('S2*MSIL1C*.SAFE').as_posix())
    unc_datasets = glob.glob(unc_datasets_dir.joinpath('*.RUT').as_posix())
    n_datasets = len(orig_datasets)
    
    # loop over the scenes. Before generating the scenarios some
    # preparation is required
    for counter, orig_dataset in enumerate(orig_datasets):

        orig_dataset_path = Path(orig_dataset)
        scene_name = orig_dataset_path.name

        logger.info(f'Working on {scene_name} ({counter+1}/{n_datasets})')

        # find corresponding uncertainty directory
        unc_dataset_path = [
            Path(x) for x in unc_datasets if Path(x).name.split('.')[0] == scene_name.split('.')[0]
        ][0]

        # preparation
        # create a subdirectory in the scenario folder with the scene name without .SAFE
        scene_name_without_safe = scene_name.replace('.SAFE', '')

        scenario_path = scenario_dir.joinpath(scene_name_without_safe)
        if not scenario_path.exists():
            scenario_path.mkdir()
        # create a folder to store the scene template (entire .SAFE structure)
        template_path = scenario_path.joinpath('template')
        if not template_path.exists():
            template_path.mkdir()

        # copy the original L1C scene into the template folder and delete
        # the jp2 files in the GRANULE directory
        shutil.copytree(
            orig_dataset_path,
            template_path.joinpath(scene_name),
            dirs_exist_ok=True
        )

        # delete the jp2 files in the template
        search_expr = '*.SAFE/GRANULE/*/IMG_DATA/*_B*.jp2'
        jp2_files = glob.glob(template_path.joinpath(search_expr).as_posix())
        for jp2_file in jp2_files:
            os.remove(jp2_file)

        logger.info(f'Created template for sampling for {scene_name} ({counter+1}/{n_datasets})')

        # finally, generate the scenarios using Monte Carlo
        gen_rad_unc_scenarios(
            orig_dataset_path=orig_dataset_path,
            unc_dataset_path=unc_dataset_path,
            scenario_path=scenario_path,
            template_path=template_path,
            n_scenarios=n_scenarios,
            roi_bounds_10m=roi_bounds_10m,
            check_contributors_only=check_contributors_only
        )

        logger.info(f'Finished MC simulations for {scene_name} ({counter+1}/{n_datasets})')

if __name__ == '__main__':

    # debug
    # orig_dataset_path = Path(
    #     './../S2A_MSIL1C_orig/done/S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651.SAFE'
    # )
    # unc_dataset_path = Path(
    #     './../S2A_MSIL1C_orig/done/S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651.RUT'
    # )
    # scenario_path = Path(
    #     './../debug'
    # )
    # template_path = Path(
    #     './../debug/template'
    # )
    # n_scenarios = 1
    # roi_bounds_10m = [7200,8400,4200,5400]  # pixel coordinates (i.e., rows and columns)
    #
    # gen_rad_unc_scenarios(orig_dataset_path, unc_dataset_path, scenario_path, template_path, n_scenarios, roi_bounds_10m)
    #

    ### define user inputs

    # TODO: correct paths afterwards
    # directory with L1C data (.SAFE subdirectories)
    orig_datasets_dir = Path('/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL1C_orig')
    
    # directory with radiometric uncertainty outputs (.RUT subdirectories)
    unc_datasets_dir = orig_datasets_dir
    
    # directory where to store the scenarios (a subdirectory will be created for each scene)
    # in which the actual scenarios are placed
    scenario_dir = Path('../S2A_MSIL1C_RUT-Scenarios/contributor_analysis')
    
    # define bounds of the study area (aka region of interest)
    # bounds col_min, col_max, row_min, row_max (image coordinates of the 10m raster)
    roi_bounds_10m = [7200,8400,4200,5400]
    
    # number of scenarios (each scenario is a possible realization of a S2 scene!)
    # TODO: set to 150 afterwards
    n_scenarios = 1 # 150
    
    main(
        orig_datasets_dir,
        unc_datasets_dir,
        scenario_dir,
        roi_bounds_10m,
        n_scenarios,
        check_contributors_only=True # TODO: set to False
    )
