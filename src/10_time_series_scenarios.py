'''
Created on Mar 7, 2022

@author: graflu
'''

import geopandas as gpd
import pandas as pd
import numpy as np
import xarray as xr

from copy import deepcopy
from datetime import date
from typing import Dict, Optional
from pathlib import Path

from agrisatpy.core.band import Band
from agrisatpy.core.raster import RasterCollection
from agrisatpy.core.sensors import Sentinel2

from phenolopy import calc_phenometrics, interpolate, remove_outliers, smooth
from _find_datasets import get_data_and_uncertainty_files, read_data_and_uncertainty
from logger import get_logger

logger = get_logger('10_time_series_scenarios')

def _calc_pheno_metrics(
        xds: xr.Dataset,
        debug: Optional[bool] = False,
        point_features: Optional[Path] = None
    ) -> Dict[str, xr.Dataset]:
    """
    Calculation of the phenological metrics using ``Phenolopy``.

    Steps:
        1. Remove outliers using a moving median filter
        2. Interpolate NaNs linearly
        3. Smooth data using Savitzky-Golay
        4. Calculate SOS, POS, EOS using the seasonal amplitude

    :param xds:
        ``xarray`` dataset containing the vegetation parameter/index
        stacked over time
    :param debug:
        if True does extract time series at single points for different
        crop types
    :param point_features:
        point features for extraction if debug == True
    :return:
        dictionary with two items: 'pheno_metrics' is a ``xarray``
        dataset with pixel-based pheno-metrics. 'ds' contains the smoothed
        per-pixel time series as another ``xarray`` dataset.
    """

    # create indexers for extracting point features (optional)
    if debug:
        gdf = gpd.read_file(point_features)
        x_indexer = xr.DataArray(gdf.geometry.x, dims=["point"])
        y_indexer = xr.DataArray(gdf.geometry.y, dims=["point"])

    # remove outliers using median filter
    ds = remove_outliers(xds, method='median', user_factor=2, z_pval=0.05)

    # interpolate nans linearly
    ds = interpolate(ds=ds, method='interpolate_na')

    # interpolate to daily values
    ds = ds.resample(time="1D").interpolate("linear")

    #smooth data using Savitzky-Golay filtering
    ds = smooth(ds=ds, method='savitsky', window_length=11, polyorder=1)

    # calculate the phenometrics using 20% seasonal amplitude
    pheno_ds = calc_phenometrics(
        da=ds['veg_index'],
        peak_metric='pos',
        base_metric='bse',
        method='seasonal_amplitude',
        factor=0.2,
        thresh_sides='two_sided'
    )

    # debug
    if debug:
        # extract data
        veg_index = ds['veg_index'].sel(x=x_indexer, y=y_indexer, method="nearest").to_dataframe()
        sos_times = pheno_ds['sos_times'].sel(x=x_indexer, y=y_indexer, method="nearest").to_dataframe()
        eos_times = pheno_ds['eos_times'].sel(x=x_indexer, y=y_indexer, method="nearest").to_dataframe()
        # convert to GeoDataFrames
        def _to_gdf(df):
            return gpd.GeoDataFrame(
                df, geometry=gpd.points_from_xy(df.x, df.y), crs=32632
            )
        veg_index = _to_gdf(df=veg_index)
        sos_times = _to_gdf(df=sos_times)
        eos_times = _to_gdf(df=eos_times)
        # join on coordinates
        joined = veg_index.sjoin(sos_times)
        joined.drop(columns=['x_left','y_left','x_right','y_right','index_right'], inplace=True)
        ts_gdf = joined.sjoin(eos_times)
        # nearest neighbor spatial join on gdf
        sample_gdf = ts_gdf.sjoin_nearest(right=gdf, lsuffix='l', rsuffix='r')
        sample_gdf.to_file('../S2_TimeSeries_Analysis/sample_pixels.gkpg', driver='GPKG')

    return {'pheno_metrics': pheno_ds, 'ds': ds}

def vegetation_time_series_scenarios(
        ts_stack_dict: Dict[date, Sentinel2],
        n_scenarios: int,
        out_dir_scenarios: Path,
        vi_name: str,
        min_val: float,
        max_val: float,
        sample_points: Path,
        fully_correlated: Optional[bool] = False
    ) -> None:
    """
    Generates the vegetation time series scenarios applying the uncertainty
    in the vegetation parameter/indices. The time series are processed in a
    similar manner TIMESAT does, using the ``Phenolopy`` package.

    To ensure the time series results are meaningful, pixel time series are
    extracted for a selected number of point features specified
    by an ESRI shapefile of point features containing information about the
    underlying crop type. A CSV file containing the extracted pixel time series
    (original data + smoothed time series from scenarios) is written to
    ``out_dir_scenarios``.

    :param ts_stack_dict:
        dictionary with ``SatDataHandler`` instances holding the vegetation parameter
        as one band and its uncertainty as another band per scene date
    :param n_scenarios:
        number of scenarios to generate
    :param out_dir_scenarios:
        directory where to save the results of the scenario runs to
    :param vi_name:
        name of the vegetation parameter/index
    :param dates:
        list of dates. Must equal the length of the ``ts_stack_list``. Is used
        to assign a date to each handler to construct the time series.
    :param min_val:
        minimum allowed value of the index or parameter in order to avoid impossible
        values
    :param max_val:
        maximum allowed value of the index or paramter in order to avoid impossible
        values
    :param sample_points:
        point features for sampling points for debugging
    :param fully_correlated:
        inter-scene correlation (full or None)
    """

    # get coordinates for xarray dataset
    coords = ts_stack_dict[list(ts_stack_dict.keys())[0]][vi_name].coordinates
    dates = ts_stack_dict.keys()
    dates_np = [np.datetime64(x) for x in dates]
    coords.update({'time': dates_np})

    # add stack atttributes for xarray dataset
    attrs = {}
    crs = ts_stack_dict[list(ts_stack_dict.keys())[0]][vi_name].crs
    attrs['crs'] = crs
    attrs['transform'] = tuple(
        ts_stack_dict[list(ts_stack_dict.keys())[0]][vi_name].geo_info.as_affine()
    )

    # loop over the scenarios. In each scenario run generate a xarray containing vegetation
    # index/parameter samples order by date (i.e., the xarray dataset is 3-dimensional: x,y,time)
    # save the calculated Pheno-metrics by the end of each scenario run to a raster file
    for scenario in range(n_scenarios):

        out_dir_scenario = out_dir_scenarios.joinpath(str(scenario+1))
        if not out_dir_scenario.exists():
            out_dir_scenario.mkdir()
        else:
            logger.info(f'Scenario {scenario+1}/{n_scenarios} already exists - skipping')
            continue

        logger.info(f'Running scenario ({scenario+1}/{n_scenarios})')
        # add vegetation samples to 3d numpy array and pass it to xarray dataset
        sample_list = []
        gdf_list = []
        orig_ts_list = []
        for _date in list(dates):
            vi_data = ts_stack_dict[_date][vi_name].values
            vi_unc = ts_stack_dict[_date][f'{vi_name}_unc'].values

            # get samples from uncertainty distribution and add it to the vi_data
            # (can be done directly because we deal with absolute uncertainties)
            # two extreme cases are possible:
            # - full inter-scene correlation
            # - no inter-scene correlation
            # (the "truth" most likely lies somewhere in between)
            if fully_correlated:
                samples = vi_unc
            #     samples = np.random.normal(0, 1, 1)[0] * samples
            else:
                samples = np.random.normal(
                    loc=0,
                    scale=vi_unc,
                    size=vi_data.shape
                )
                samples += vi_data
                # constrain samples to min and max of the parameter
                samples[samples > max_val] = max_val
                samples[samples < min_val] = min_val
            sample_list.append(samples)
            orig_ts_list.append(vi_data)

        stacked = np.ma.stack(sample_list)
        stacked = stacked.filled(np.nan)
        stacked_org = np.ma.stack(orig_ts_list)
        stacked_org  = stacked_org.filled(np.nan)
        # fully correlated case - sample along the temporal axis
        if fully_correlated:
            for row in range(stacked.shape[1]):
                for col in range(stacked.shape[2]):
                    if np.isnan(stacked[:,row,col]).all(): continue
                    stacked[:,row,col] = np.random.normal(0,1,1)[0] * stacked[:,row,col] + stacked_org[:,row,col]
                    # constrain samples to min and max of the parameter
                    stacked[:,row,col][stacked[:,row,col] > max_val] = max_val
                    stacked[:,row,col][stacked[:,row,col] < min_val] = min_val

        # create the 3d numpy array to pass as xarray dataset
        stack = {'veg_index': tuple([('time', 'y','x'), stacked])}

        # create stack for the reference run and the pixel time series of the reference
        # data (do it in the first iteration only to avoid overwritten the data again and
        # again)
        if scenario == 0:
            orig_stack = {'veg_index': tuple([('time','y','x'), stacked_org])}

        # gdf[(gdf.row == 690)&(gdf.col == 478)][['date','GLAI']]
        # construct xarray dataset for the scenario run
        try:
            xds = xr.Dataset(
                stack,
                coords=coords,
                attrs=attrs
            )
            # xds = xds.resample(time="1D").interpolate("linear")
            res = _calc_pheno_metrics(xds)
            pheno_ds = res['pheno_metrics'] # pheno metric results
            ds = res['ds'] # smoothed time series values
        except Exception as e:
            logger.error(f'Error in scenario ({scenario+1}/{n_scenarios}): {e}')

        # do the same for the reference run
        if scenario == 0:
            try:
                xds_ref = xr.Dataset(
                    orig_stack,
                    coords=coords,
                    attrs=attrs
                )
                # xds_ref = xds_ref.resample(time="1D").interpolate("linear")
                res_ref = _calc_pheno_metrics(
                    xds_ref,
                    debug=True,
                    point_features=sample_points
                )
                pheno_ds_ref = res_ref['pheno_metrics'] # pheno metric results
                ds_ref = res_ref['ds'] # smoothed time series values
            except Exception as e:
                logger.error(f'Error in reference run: {e}')

        # save to raster files using the SatDataHandler object since it has the full geoinformation
        # available
        pheno_handler = deepcopy(ts_stack_dict[list(ts_stack_dict.keys())[0]])
        pheno_metrics = [
            'sos_values', 'sos_times', 'pos_values', 'pos_times', 'eos_values', 'eos_times'
        ]

        out_file = out_dir_scenario.joinpath('pheno_metrics.tif')

        for pheno_metric in pheno_metrics:
            pheno_handler.add_band(
                band_constructor=Band,
                band_name=pheno_metric,
                values=eval(f'pheno_ds.{pheno_metric}.data'),
                geo_info=pheno_handler[vi_name].geo_info
            )

        # remove "template" files
        bands_to_drop = [vi_name, f'{vi_name}_unc']
        for band_to_drop in bands_to_drop:
            pheno_handler.drop_band(band_to_drop)

        # save to geoTiff
        pheno_handler.to_rasterio(out_file)

        # write original data
        if scenario == 0:
            pheno_handler = None
            pheno_handler = deepcopy(ts_stack_dict[list(ts_stack_dict.keys())[0]])

            out_dir_scenario = out_dir_scenarios.joinpath('reference')
            if not out_dir_scenario.exists():
                out_dir_scenario.mkdir()
            out_file = out_dir_scenario.joinpath('pheno_metrics.tif')

            for pheno_metric in pheno_metrics:
                pheno_handler.add_band(
                    band_constructor=Band,
                    band_name=pheno_metric,
                    values=eval(f'pheno_ds_ref.{pheno_metric}.data'),
                    geo_info=pheno_handler[vi_name].geo_info
                )

            # remove "template" files
            bands_to_drop = [vi_name, f'{vi_name}_unc']
            for band_to_drop in bands_to_drop:
                pheno_handler.drop_band(band_to_drop)

            # save to geoTiff
            pheno_handler.to_rasterio(out_file)

        logger.info(f'Finished scenario ({scenario+1}/{n_scenarios})')

def main(
        vi_dir: Path,
        uncertainty_analysis_dir: Path,
        out_dir_scenarios: Path,
        vi_name: str,
        sample_polygons: Path,
        ymin: float,
        ymax: float,
        n_scenarios: int,
        crop_periods: Path,
        sample_points: Path,
        fully_correlated: Optional[bool] = False
    ):
    """
    main executable function of this module generating the time series plots

    ATTENTION:
        The functions in this module require some computation time!
    """

    # search files, join and order them by date
    data_df = get_data_and_uncertainty_files(
        vi_dir=vi_dir,
        uncertainty_analysis_dir=uncertainty_analysis_dir,
        vi_name=vi_name
    )

    # read the data and mask clouds, shadows, and unclassified pixels based on SCL information
    file_df, ts_stack_dict = read_data_and_uncertainty(
        data_df=data_df,
        vi_name=vi_name,
        parcels=sample_polygons,
        crop_periods=crop_periods
    )

    # generate scenarios and calculate the pheno-metrics
    if fully_correlated:
        out_dir_scenarios_c = out_dir_scenarios.joinpath('fully_correlated')
    else:
        out_dir_scenarios_c = out_dir_scenarios.joinpath('uncorrelated')
    out_dir_scenarios_c.mkdir(exist_ok=True)

    vegetation_time_series_scenarios(
        ts_stack_dict=ts_stack_dict,
        n_scenarios=n_scenarios,
        out_dir_scenarios=out_dir_scenarios_c,
        vi_name=vi_name,
        min_val=ymin,
        max_val=ymax,
        fully_correlated=fully_correlated,
        sample_points=sample_points
    )

if __name__ == '__main__':

    import sys

    # original Sentinel-2 scenes with vegetation indices
    vi_dir = Path(
        # '../S2A_MSIL1C_orig/*.VIs'
        '/home/graflu/Documents/uncertainty/S2_MSIL1C_orig/*.VIs'
    )
    
    # directory with uncertainty analysis results
    uncertainty_analysis_dir = Path(
        '../S2_MSIL2A_Analysis'
    )

    # define sample polygons (for visualizing the uncertainty per crop type over time)
    sample_polygons = Path('../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')

    # define point sampling locations for visualizing pixel time series
    sample_points = Path('../shp/ZH_Points_2019_EPSG32632_selected-crops.shp')


    # define key crop growth periods
    crop_periods = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/01_Uncertainty/ESCH/scripts_paper_uncertainty/S2_TimeSeries_Analysis/crop_growth_periods-CH.csv'
    )

    # vegetation index to consider
    # arg = sys.argv
    # with open(arg[1], 'r') as src:
    #     vi = src.readlines()
    #     vi_names = [vi[0].replace('\n','')]
    vi_names = ['EVI', 'NDVI', 'GLAI']
    ymins = {'NDVI': -1, 'EVI': -1, 'GLAI': 0}
    ymaxs = {'NDVI': 1, 'EVI': 1, 'GLAI': 7}

    # number of scenarios to generate
    n_scenarios = 1000

    # directory where to save phenological metrics to
    out_dir_scenarios = Path(f'../S2_TimeSeries_Analysis_Test')
    if not out_dir_scenarios.exists():
        out_dir_scenarios.mkdir()

    fully_correlated = [True] # [False, True]

    for idx, vi_name in enumerate(vi_names):

        out_dir_scenarios_vi = out_dir_scenarios.joinpath(vi_name)
        if not out_dir_scenarios_vi.exists():
            out_dir_scenarios_vi.mkdir()

        for corr in fully_correlated:
            main(
                vi_dir=vi_dir,
                uncertainty_analysis_dir=uncertainty_analysis_dir,
                out_dir_scenarios=out_dir_scenarios_vi,
                vi_name=vi_name,
                sample_polygons=sample_polygons,
                ymin=ymins[vi_name],
                ymax=ymaxs[vi_name],
                fully_correlated=corr,
                crop_periods=crop_periods,
                n_scenarios=n_scenarios,
                sample_points=sample_points
            )
