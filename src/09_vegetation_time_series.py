'''
Created on Dec 6, 2021

@author: Lukas Graf
'''

import cv2
import glob
import pandas as pd
import numpy as np

from pathlib import Path
from datetime import datetime
from uncertainties import unumpy
from typing import Tuple
from typing import Optional
from copy import deepcopy
from scipy.interpolate import interp1d

from agrisatpy.io import Sat_Data_Reader
from agrisatpy.io.sentinel2 import S2_Band_Reader
from agrisatpy.utils.constants.sentinel2 import ProcessingLevels


def ts_temporal_compositing(
        ts_values: np.array,
        dates: pd.Series,
        composite_length: Optional[int]=10,
        composition_unit: Optional[str]='d',
        method: Optional[str]='max',
    ) -> pd.DataFrame:
    """
    creates a temporal composite of a time series by aggregating
    (irregularly spaced) time series data at regular temporal units
    (e.g., every 10 days) using one of the following aggregation
    methods:
    
    * mean
    * min
    * max
    * median
    
    If no data points are located within one of the aggregation
    periods, NaN is inserted.

    IMPORTANT: The dataframe must be sorted by date (asc) before
    passing it to this function!

    :param ts_df:
        dataframe containing the temporal information and time
        series data to smooth
    :param colname_values:
        name of the column holding the time series
        values (e.g., NDVI)
    :param composite_length:
        length of the composite in terms if the specified
        temporal unit (see composition_unit)
    :param colname_time:
        name of the column holding the temporal
        dimension (e.g., dates or timestamps) if not already in index
    :param composition_unit:
        temporal unit to be used for creating the composite.
        Default is d(ays).
    :param method:
        aggregation method for creating the composite. Must be
        one out of 'min', 'max', 'mean', 'median'. The default
        is 'max'
    """

    # the creation of the composite starts at the first timestamp available
    local = pd.DataFrame({'date': dates, 'values': ts_values})
    local.index = pd.to_datetime(local.date)
    local.drop('date', inplace=True, axis=1)

    local= local['values'].resample(
        f'{composite_length}{composition_unit}', closed='right'
    ).agg([method])
    # add offset to resampled time stamps since the time series is shifted
    # otherwise to the "left"
    local.index = local.index + pd.DateOffset(int(0.5*composite_length))
    return local


def lin_interpol(
            ts_df: pd.Dataframe
        ) -> None:
        """
        linear interpolation of time series data
        points
        """

        # reindex dataframe between start and end date and interpolate
        rng = pd.date_range(ts_df.index.min(), ts_df.index.max())
        return ts_df.reindex(rng).interpolate()


def moving_average_smooting(
        ts_df: pd.DataFrame,
        colname_values: str,
        window_size: int,
        colname_time: Optional[str]=None,
    ) -> pd.DataFrame:
    """
    Applies a moving average to a time series to smooth
    it using a user-defined window size.

    :param ts_df:
        dataframe containing the temporal information and time
        series data to smooth
    :param colname_values:
        name of the column holding the time series
        values (e.g., NDVI)
    :param window_size:
        size of the moving window in the temporal unit provided
        in the time column (usually days)
    :param colname_time:
        name of the column holding the temporal
        dimension (e.g., dates or timestamps) if not already in index
    """
    local = deepcopy(ts_df)
    if colname_time is not None:
        local.index = pd.to_datetime(local[colname_time])
        local.drop(colname_time, inplace=True, axis=1)
    return local[colname_values].rolling(window=window_size).mean().iloc[window_size-1:].values


def get_data_and_uncertainty_files(
        orig_dataset_dir: Path,
        uncertainty_analysis_dir: Path,
        vi_name: str
    ) -> pd.DataFrame:
    """
    Reads the bandstacked geoTiff (S2 bands) from the original
    Sentinel-2 data and searches for the corresponding relative uncertainty
    information. Stores the files by acquisition date in a pandas dataframe

    :param orig_dataset_dir:
        directory with the original Sentinel-2 data (*.SAFE datasets)
    :param uncertainty_analysis_dir:
        directory where the corresponding relative uncertainty for the index
        or parameter is stored.
    :param vi_name:
        name of the vegetation index or parameter. It is used to find the
        correct files
    :return:
        pandas dataframe with the file location of the bandstacks
        and file with corresponding uncertainty for the selected vegetation
        index/ parameter alongside with the (acquisition) date
    """

    # define search expressions
    # vegetation parameter/ index
    search_expr_bandstack = orig_dataset_dir.joinpath(
        'S2*_MSIL2A_*.SAFE'
    )
    search_expr_unc = uncertainty_analysis_dir.joinpath(
        f'S2*/L3_{vi_name.upper()}_*.tif'
    )

    # get list of files
    bandstack_list = glob.glob(search_expr_bandstack.as_posix())
    unc_file_list = glob.glob(search_expr_unc.as_posix())

    # convert lists to dataframe
    orig_file_df = pd.DataFrame(bandstack_list, columns=['filename_orig'])
    unc_file_df = pd.DataFrame(unc_file_list, columns=['filename_unc'])

    # extract dates (files are matched and ordered by date)
    orig_file_df['date'] = orig_file_df.filename_orig.apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[2][0:8], '%Y%m%d').date()
    )
    unc_file_df['date'] = unc_file_df.filename_unc.apply(
        lambda x: datetime.strptime(Path(x).parent.name.split('_')[2][0:8], '%Y%m%d').date()
    )

    # join dataframes on date
    df = pd.merge(orig_file_df, unc_file_df, on='date')
    # and sort by date
    df.sort_values(by='date', inplace=True)

    return df

def threshold_based_phenology(ts_df, colname_values, colname_time, amplitude_threshold):
    
    res = {}
    pos_idx = np.argmax(ts_df[colname_values])
    res['POS'] = ts_df[colname_time].iloc[pos_idx]

    # define amplitude as min-max spread
    amplitude = np.max(ts_df[colname_values]) - \
        np.min(ts_df[colname_values])

    # define upward branch as located "left" from the POS index
    amplitude_critical = amplitude *  amplitude_threshold / 100.
    sos_idx = np.where(
        ts_df[colname_values].iloc[0:pos_idx] > amplitude_critical)[0][0]
    res['SOS'] = ts_df[colname_time].iloc[sos_idx]
    # and downward right from it
    eos_idx = pos_idx + np.where(
        ts_df[colname_values].iloc[pos_idx+1:] > amplitude_critical)[0][-1]
    res['EOS'] = ts_df[colname_time].iloc[eos_idx]

    res['LENGTH'] = (res['EOS'] - res['SOS']).days

    return res

    
def read_data_and_uncertainty(
        data_df: pd.DataFrame,
        in_file_aoi: Path,
        vi_name: str,
        out_dir_plots: Path
    ) -> Tuple[pd.DataFrame, np.array]:
    """
    This function reads the selected vegetation index (computed on the fly)
    and the associated standard uncertainty in a uarray (from uncertainties)
    and stores the read data as new columns in data_df. In addition, plots
    of the RGB, False-Color Infrared and SCL and the vegetation index are generated.

    :param data_df:
        dataframe returned from ``get_data_and_uncertainty``
    :param in_file_aoi:
        shapefile defining the study area
    :param vi_name:
        name of the vegetation index to compute (e.g., NDVI)
    :param out_dir_plots:
        directory where to save the resulting plots to
        (are stored per image acquisition date)
    :return:
        input dataframe with read data + raster stack (including
        standard uncertainty)
    """

    # loop over datasets (single acquisition dates) and read the data
    # (vegetation index/ parameter + standard uncertainty)
    update_list = []
    ts_stack_list = []

    for _, record in data_df.iterrows():

        # read S2 bandstack data (including scene classification layer, SCL)
        s2_stack = S2_Band_Reader()
        s2_stack.read_from_safe(
            in_dir=Path(record.filename_orig),
            processing_level=ProcessingLevels.L2A,
            in_file_aoi=in_file_aoi
        )

        # mask clouds and cloud shadows (if any) and store the cloudy pixel percentage
        cloud_coverage = s2_stack.get_cloudy_pixel_percentage()

        # plot quicklooks and save them
        plot_dir = out_dir_plots.joinpath(record.date.strftime('%Y%m%d'))
        if not plot_dir.exists():
            plot_dir.mkdir()
        fname_rgb = plot_dir.joinpath('rgb_quicklook.png')
        fname_nir = plot_dir.joinpath('falsecolor-nir_quicklook.png')
        fname_scl = plot_dir.joinpath('scene-classification_quicklook.png')
        fname_vi = plot_dir.joinpath(f'{vi_name.lower()}-nominal_quicklook.png')
        fname_unc = plot_dir.joinpath(f'{vi_name.lower()}-stdunc_quicklook.png')

        fig_rgb = s2_stack.plot_rgb()
        fig_rgb.savefig(fname=fname_rgb, bbox_inches='tight')
        fig_nir = s2_stack.plot_false_color_infrared()
        fig_nir.savefig(fname=fname_nir, bbox_inches='tight')
        fig_scl = s2_stack.plot_scl()
        fig_scl.savefig(fname=fname_scl, bbox_inches='tight')

        # resample spectral bands to 10m (i.e., those bands with 20m resolution)
        s2_stack.resample(
            target_resolution=10,
            resampling_method=cv2.INTER_CUBIC,
            bands_to_exclude=['scl']
        )

        # resample SCL band to 10m
        s2_stack.resample(
            target_resolution=10,
            resampling_method=cv2.INTER_NEAREST_EXACT
        )

        # calculate the vegetation index
        s2_stack.calc_vi(vi=vi_name)
        
        fig_vi = s2_stack.plot_band(band_name=vi_name, colormap='summer')
        fig_vi.savefig(fname=fname_vi, bbox_inches='tight')

        s2_stack.mask_clouds_and_shadows(bands_to_mask=[vi_name])

        # read uncertainty data
        uncertainty_band = Sat_Data_Reader()
        uncertainty_band.read_from_bandstack(fname_bandstack=record.filename_unc)

        fig_unc = uncertainty_band.plot_band(
            band_name='abs_stddev'
        )
        fig_unc.savefig(fname=fname_unc, bbox_inches='tight')

        # apply the cloud mask also to the absolute uncertainty band
        uncertainty_band.add_band(
            band_name='scl',
            band_data=s2_stack.data['scl']
        )

        uncertainty_band.mask(
            name_mask_band='scl',
            mask_values=[3, 8, 9, 10],
            bands_to_mask=['abs_stddev'],
            keep_mask_values=False
        )

        # store data and standard uncertainty
        uarray = unumpy.uarray(
            nominal_values=s2_stack.data[vi_name],
            std_devs=uncertainty_band.data['abs_stddev']
        )
        ts_stack_list.append(uarray)

        update_list.append({
            'date': record.date,
            'cloudy_pixel_percentage': cloud_coverage
        })

    # stack arrays
    stack_ts = np.dstack(ts_stack_list)

    # added update columns to input data frame
    update_df = pd.DataFrame(update_list)
    merged = pd.merge(data_df, update_df, on='date')

    return merged, stack_ts


def phenological_uncertainty(
        dates: pd.Series,
        stack_ts: np.array
    ):
    """
    
    """

    nrows, ncols = stack_ts.shape[0], stack_ts.shape[1]

    for nrow in range(nrows):
        for ncol in range(ncols):
            pixel_ts = stack_ts[nrow,ncol,:].data
            pixel_ts_vals = np.array([x.n for x in pixel_ts])
            pixel_ts_unc = np.array([x.s for x in pixel_ts])

            pixel_phenology_uncertainty(dates, pixel_ts_vals, pixel_ts_unc)


def pixel_phenology_uncertainty(
        dates: pd.Series,
        pixel_ts_vals: np.array,
        pixel_ts_unc: np.array,
        n_scenarios: Optional[int] = 1000
    ):
    """
    Function to calculate the uncertainty of key phenological metrics
    using a set of vegetation index/ parameter observations from different
    points in time and their standard uncertainty to generate `n_scenarios`` 
    of pixel time series.
    """

    n_dates = len(dates)

    for scenario in n_scenarios:

        # generate sample from normal distribution
        error_samples = np.random.normal(
            loc=0,
            scale=pixel_ts_unc,
            size=n_dates
        )
        ts_scenario = pixel_ts_vals + error_samples

        # maximum value 10-day composite
        ts_df = ts_temporal_compositing(
            ts_values=ts_scenario,
            dates=dates
        )

        # linear interpolation
        ts_df_lin = lin_interpol(ts_df=ts_df)

        # moving average smoothing
        ts_df_sm = moving_average_smooting(
            ts_df=ts_df,
            colname_values='values',
            window_size=11,
            colname_time='date'
        )

        # extract phenological metrics
        


if __name__ == '__main__':
    # original Sentinel-2 scenes with vegetation indices
    orig_dataset_dir = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL1C_orig'
    )
    
    # directory with uncertainty analysis results
    uncertainty_analysis_dir = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL2A_Analysis'
    )
    
    # vegetation index to consider
    vi_name = 'NDVI'
    data_df = get_data_and_uncertainty_files(
        orig_dataset_dir=orig_dataset_dir,
        uncertainty_analysis_dir=uncertainty_analysis_dir,
        vi_name=vi_name
    )

    in_file_aoi = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/shp/AOI_Esch_EPSG32632.shp'
    )

    out_dir_plots = uncertainty_analysis_dir.joinpath('scene_quicklooks')
    if not out_dir_plots.exists():
        out_dir_plots.mkdir()

    unc_df, stack_ts = read_data_and_uncertainty(
        data_df=data_df,
        vi_name=vi_name,
        in_file_aoi=in_file_aoi,
        out_dir_plots=out_dir_plots
    )

    dates = unc_df.date

      