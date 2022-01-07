'''
Similar to the generation of the L1C TOA scenarios, this script generates
a set of vegetation time series scenarios.

For each scenario, the timing and value of the start, peak and end of season
is calculated using the Phenolopy package that is heavily inspired by TIMESAT.

Many thanks @lewistrotter for providing Phenolopy
(see: https://github.com/lewistrotter/Phenolopy) avialable under Apache 2.0 licence.
'''

import glob
import logging
import pandas as pd
import numpy as np
import xarray as xr

from pathlib import Path
from datetime import datetime
from datetime import date
from typing import Tuple
from typing import List
from typing import Dict
from copy import deepcopy
import geopandas as gpd

from _phenolopy import remove_outliers
from _phenolopy import interpolate
from _phenolopy import smooth
from _phenolopy import calc_phenometrics
from agrisatpy.io import SatDataHandler
from agrisatpy.io.sentinel2 import Sentinel2Handler

from logger import get_logger

# setup logger -> will write log file to the /../log directory
logger = get_logger('l4_phenology')


def _calc_pheno_metrics(xds: xr.Dataset) -> Dict[str, xr.Dataset]:
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
    :return:
        dictionary with two items: 'pheno_metrics' is a ``xarray``
        dataset with pixel-based pheno-metrics. 'ds' contains the smoothed
        per-pixel time series as another ``xarray`` dataset.
    """
    
    # remove outliers using median filter
    ds = remove_outliers(ds=xds, method='median', user_factor=2, z_pval=0.05)

    # interpolate nans linearly
    ds = interpolate(ds=ds, method='interpolate_na')

    # smooth data using Savitzky-Golay filter
    ds = smooth(ds=ds, method='savitsky', window_length=3, polyorder=1)

    # calculate the phenometrics
    pheno_ds = calc_phenometrics(
        da=ds['veg_index'],
        peak_metric='pos',
        base_metric='bse',
        method='seasonal_amplitude',
        factor=0.2,
        thresh_sides='two_sided'
    )

    return {'pheno_metrics': pheno_ds, 'ds': ds}


def get_data_and_uncertainty_files(
        vi_dir: Path,
        uncertainty_analysis_dir: Path,
        vi_name: str
    ) -> pd.DataFrame:
    """
    Searches the vegetation parameter/index generated from the original
    Sentinel-2 data alongside with the corresponding relative uncertainty
    information. Stores the files by acquisition date in a pandas dataframe.

    :param vi_dir:
        directory containing the vegetation indices/ parameters derived from the
        original (i.e., non-scenario-based) datasets and the resampled original
        Sentinel-2 data in L2A processing level
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
    search_expression_bandstack = vi_dir.joinpath(
        '*_pixel_division_10m.tiff'
    )
    search_expr_veg_par = vi_dir.joinpath(
        Path('Vegetation_Indices').joinpath(f'VI_*_{vi_name.upper()}.tif')
    )
    search_expr_unc = uncertainty_analysis_dir.joinpath(
        f'S2*/L3_{vi_name.upper()}_*.tif'
    )

    # get list of files
    bandstack_list = glob.glob(search_expression_bandstack.as_posix())
    veg_par_list = glob.glob(search_expr_veg_par.as_posix())
    unc_file_list = glob.glob(search_expr_unc.as_posix())

    # convert lists to dataframe
    bandstack_file_df = pd.DataFrame(bandstack_list, columns=['filename_bandstack'])
    veg_par_file_df = pd.DataFrame(veg_par_list, columns=['filename_veg_par'])
    unc_file_df = pd.DataFrame(unc_file_list, columns=['filename_unc'])

    # extract dates (files are matched and ordered by date)
    veg_par_file_df['date'] = veg_par_file_df.filename_veg_par.apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[1], '%Y%m%d').date()
    )
    unc_file_df['date'] = unc_file_df.filename_unc.apply(
        lambda x: datetime.strptime(Path(x).parent.name.split('_')[2][0:8], '%Y%m%d').date()
    )
    bandstack_file_df['date'] = bandstack_file_df.filename_bandstack.apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[0], '%Y%m%d').date()
    )

    # join dataframes on date
    tmp_df = pd.merge(veg_par_file_df, unc_file_df, on='date')
    df = pd.merge(bandstack_file_df, tmp_df, on='date')
    # and sort by date
    df.sort_values(by='date', inplace=True)

    return df


def read_data_and_uncertainty(
        data_df: pd.DataFrame,
        in_file_aoi: Path,
        vi_name: str,
        out_dir_plots: Path
    ) -> Tuple[pd.DataFrame, List[Sentinel2Handler]]:
    """
    This function reads the selected vegetation index (computed in step 7)
    and the associated standard uncertainty in a uarray (from uncertainties)
    and stores the read data as new columns in data_df.

    :param data_df:
        dataframe returned from ``get_data_and_uncertainty``
    :param in_file_aoi:
        shapefile defining the study area
    :param vi_name:
        name of the vegetation index or parameter to process (e.g., NDVI)
    :param out_dir_plots:
        directory where to store scene quicklooks (might be helpful for
        result interpretation)
    :return:
        input dataframe with read data + list of read vegetation data and
        their absolute uncertainty per image acquisition date
    """

    # loop over datasets (single acquisition dates) and read the data
    # (vegetation index/ parameter + standard uncertainty)
    update_list = []
    ts_stack_dict = {}
    orig_vi_dict = {}
    unc_dict = {}

    for _, record in data_df.iterrows():

        # read S2 banstack plus its SCL file to mask out clouds and shadows
        s2_stack = Sentinel2Handler()
        s2_stack.read_from_bandstack(
            fname_bandstack=Path(record.filename_bandstack),
            in_file_aoi=in_file_aoi,
            band_selection=['B02','B03','B04','B08']
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

        # read Vegetation parameter/index band
        handler = SatDataHandler()
        handler.read_from_bandstack(
            fname_bandstack=Path(record.filename_veg_par),
            in_file_aoi=in_file_aoi
        )
        handler.reset_bandnames(new_bandnames=[vi_name])
        fig_vi = handler.plot_band(band_name=vi_name, colormap='summer')
        fig_vi.savefig(fname=fname_vi, bbox_inches='tight')
        s2_stack.add_band(
            band_name=vi_name,
            band_data=handler.get_band(vi_name),
            snap_band='B02'
        )
        # save the file path to dict
        orig_vi_dict[record.date] = Path(record.filename_veg_par)

        # read vegetation parameter/index uncertainty band
        handler = SatDataHandler()
        handler.read_from_bandstack(
            fname_bandstack=Path(record.filename_unc),
            in_file_aoi=in_file_aoi,
            band_selection=['abs_stddev']
        )

        fig_unc = handler.plot_band(
            band_name='abs_stddev'
        )
        fig_unc.savefig(fname=fname_unc, bbox_inches='tight')

        # apply the cloud mask also to the absolute uncertainty band
        s2_stack.add_band(
            band_name=f'{vi_name}_unc',
            band_data=handler.get_band('abs_stddev'),
            snap_band='B02'
        )
        unc_dict[record.date] = Path(record.filename_unc)

        # mask clouds, shadows, cirrus, dark area and unclassified pixels using the
        # information in the scene classification layer
        s2_stack.mask_clouds_and_shadows(
            bands_to_mask=[vi_name, f'{vi_name}_unc'],
            cloud_classes=[2,3,7,8,9,10]
        )

        update_list.append({
            'date': record.date,
            'cloudy_pixel_percentage': cloud_coverage
        })

        bands_to_drop = ['B02','B03','B04','B08','SCL']
        for band_to_drop in bands_to_drop:
            s2_stack.drop_band(band_to_drop)
        ts_stack_dict[record.date] = s2_stack

    # added update columns to input data frame
    update_df = pd.DataFrame(update_list)
    merged = pd.merge(data_df, update_df, on='date')

    # backup met
    merged.to_csv(out_dir_plots.joinpath('metadata.csv'), index=False)

    return merged, ts_stack_dict, orig_vi_dict, unc_dict


def vegetation_time_series_scenarios(
        ts_stack_dict: Dict[date, Sentinel2Handler],
        orig_vi_dict: Dict[date, Path],
        unc_dict: Dict[date, Path],
        n_scenarios: int,
        out_dir_scenarios: Path,
        vi_name: str,
        sample_points: Path
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

    :param ts_stack_list:
        list of ``SatDataHandler`` instances holding the vegetation parameter
        as one band and its uncertainty as another band.
    :param n_scenarios:
        number of scenarios to generate
    :param out_dir_scenarios:
        directory where to save the results of the scenario runs to
    :param vi_name:
        name of the vegetation parameter/index
    :param dates:
        list of dates. Must equal the length of the ``ts_stack_list``. Is used
        to assign a date to each handler to construct the time series.
    :param sample_points:
        shapefile with points to use for time series pixels samples. The extracted
        time series are plotted and saved to a sub-directory called 'pixel_plots'
    """

    # get coordinates for xarray dataset
    coords = ts_stack_dict[list(ts_stack_dict.keys())[0]].get_coordinates(vi_name)
    dates = ts_stack_dict.keys()
    dates_np = [np.datetime64(x) for x in dates]
    coords.update({'time': dates_np})

    # add stack atttributes for xarray dataset
    attrs = {}
    crs = ts_stack_dict[list(ts_stack_dict.keys())[0]].get_epsg(vi_name)
    attrs['crs'] = crs
    attrs['transform'] = tuple(ts_stack_dict[list(ts_stack_dict.keys())[0]].get_meta(vi_name)['transform'])

    # calculate the x and y array indices required to extract the pixel values
    # at the sample locations from the array holding the smoothed vegetation
    gdf_points = gpd.read_file(sample_points)
    gdf_points['x'] = gdf_points.geometry.x
    gdf_points['y'] = gdf_points.geometry.y

    # map the coordinates to array indices
    def find_nearest_array_index(array, value):
        return np.abs(array - value).argmin()
    # get column (x) indices
    gdf_points['col'] = gdf_points['x'].apply(
        lambda x, coords=coords, find_nearest_array_index=find_nearest_array_index:
        find_nearest_array_index(coords['x'], x)
    )
    # get row (y) indices
    gdf_points['row'] = gdf_points['y'].apply(
        lambda y, coords=coords, find_nearest_array_index=find_nearest_array_index:
        find_nearest_array_index(coords['y'], y)
    )

    # loop over the scenarios. In each scenario run generate a xarray containing vegetation
    # index/parameter samples order by date (i.e., the xarray dataset is 3-dimensional: x,y,time)
    # save the calculated Pheno-metrics by the end of each scenario run to a raster file
    orig_ts_list = []
    for scenario in range(n_scenarios):

        logger.info(f'Running scenario ({scenario+1}/{n_scenarios})')
        # add vegetation samples to 3d numpy array and pass it to xarray dataset
        sample_list = []
        gdf_list = []
        for _date in list(dates):
            vi_data = ts_stack_dict[_date].get_band(vi_name).data
            vi_unc = ts_stack_dict[_date].get_band(f'{vi_name}_unc').data

            # sample the test points
            # sample original pixel values only in first scenario run
            if scenario == 0:
                vegpar_file = orig_vi_dict[_date]
                # TODO: sample from array directly to avoid potential shifts
                gdf_sample_data = SatDataHandler.read_pixels(
                    point_features=sample_points,
                    raster=vegpar_file
                )
                colnames = list(gdf_sample_data.columns)
                if None in colnames:
                    colnames[colnames.index(None)] = vi_name
                gdf_sample_data.columns = colnames
                
                vegpar_unc = unc_dict[_date]
            
                gdf_sample_unc = SatDataHandler.read_pixels(
                    point_features=sample_points,
                    raster=vegpar_unc,
                    band_selection=['abs_stddev']
                )
                gdf_sample_unc.rename({'abs_stddev': 'unc'}, axis=1, inplace=True)

                # join the data frames
                gdf = pd.merge(
                    gdf_sample_data[['id', 'crop_type', 'geometry', vi_name]],
                    gdf_sample_unc[['id', 'unc']],
                    on='id'
                )
                gdf['date'] = _date
                # join row and column
                gdf_joined = pd.merge(
                    gdf,
                    gdf_points[['id', 'row', 'col']],
                    on='id'
                )
                gdf_list.append(gdf_joined)

                # append original raster data to stack once
                orig_ts_list.append(vi_data)

            # get samples from uncertainty distribution and add it to the vi_data
            # (can be done directly because we deal with absolute uncertainties)
            samples = np.random.normal(
                loc=0,
                scale=vi_unc,
                size=vi_data.shape
            )
            samples += vi_data
            sample_list.append(samples)

        # create the 3d numpy array to pass as xarray dataset
        stack = {'veg_index': tuple([('y','x','time'), np.dstack(sample_list)])}

        # create stack for the reference run and the pixel time series of the reference
        # data (do it in the first iteration only to avoid overwritten the data again and
        # again)
        if scenario == 0:
            orig_stack = {'veg_index': tuple([('y','x','time'), np.dstack(orig_ts_list)])}
            gdf = pd.concat(gdf_list)

        # construct xarray dataset for the scenario run
        try:
            xds = xr.Dataset(
                stack,
                coords=coords,
                attrs=attrs
            )
            res = _calc_pheno_metrics(xds)
            pheno_ds = res['pheno_metrics'] # pheno metric results
            ds = res['ds'] # smoothed time series values
        except Exception as e:
            logger.error(f'Error in scenario ({scenario+1}/{n_scenarios}): {e}')

        # get the smoothed time series values from the dataset at the sample locations
        # unfortunately, there seems to be no ready-to-use solution
        unique_points = gdf[gdf.date == gdf.date.unique()[0]][['id', 'row', 'col']]

        scenario_col = f'{vi_name}_{scenario+1}'
        gdf[scenario_col] = np.empty(gdf.shape[0])
        # loop over sample points and add them as new entries to the dataframe
        for _, unique_point in unique_points.iterrows():
            ts_values = ds[dict(x=[unique_point.col], y=[unique_point.row])]['veg_index'].data
            gdf.loc[
                (gdf.row == unique_point.row) & (gdf.col == unique_point.col) & (gdf.id == unique_point.id),
                scenario_col
            ] = ts_values[0,0,:]

        # do the same for the reference run
        if scenario == 0:
            try:
                xds_ref = xr.Dataset(
                    orig_stack,
                    coords=coords,
                    attrs=attrs
                )
                res_ref = _calc_pheno_metrics(xds_ref)
                pheno_ds_ref = res_ref['pheno_metrics'] # pheno metric results
                ds_ref = res_ref['ds'] # smoothed time series values
            except Exception as e:
                logger.error(f'Error in reference run: {e}')

            ref_col = f'{vi_name}_ts_sm'
            gdf[ref_col] = np.empty(gdf.shape[0])
            # loop over sample points and add them as new entries to the dataframe
            for _, unique_point in unique_points.iterrows():
                ts_values = ds_ref[dict(x=[unique_point.col], y=[unique_point.row])]['veg_index'].data
                gdf.loc[
                    (gdf.row == unique_point.row) & (gdf.col == unique_point.col) & (gdf.id == unique_point.id),
                    ref_col
                ] = ts_values[0,0,:]

        # save to raster files using the SatDataHandler object since it has the full geoinformation
        # available
        pheno_handler = deepcopy(ts_stack_dict[list(ts_stack_dict.keys())[0]])
        pheno_metrics = [
            'sos_values', 'sos_times', 'pos_values', 'pos_times', 'eos_values', 'eos_times'
        ]

        out_dir_scenario = out_dir_scenarios.joinpath(str(scenario+1))
        if not out_dir_scenario.exists():
            out_dir_scenario.mkdir()
        out_file = out_dir_scenario.joinpath('pheno_metrics.tif')

        for pheno_metric in pheno_metrics:
            pheno_handler.add_band(
                band_name=pheno_metric,
                band_data=eval(f'pheno_ds.{pheno_metric}.data'),
                snap_band=vi_name
            )

        # remove "template" files
        bands_to_drop = [vi_name, f'{vi_name}_unc']
        for band_to_drop in bands_to_drop:
            pheno_handler.drop_band(band_to_drop)

        # save to geoTiff
        pheno_handler.write_bands(out_file)

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
                    band_name=pheno_metric,
                    band_data=eval(f'pheno_ds_ref.{pheno_metric}.data'),
                    snap_band=vi_name
                )

            # remove "template" files
            bands_to_drop = [vi_name, f'{vi_name}_unc']
            for band_to_drop in bands_to_drop:
                pheno_handler.drop_band(band_to_drop)

            # save to geoTiff
            pheno_handler.write_bands(out_file)

        logger.info(f'Finished scenario ({scenario+1}/{n_scenarios})')

    # save sample points to CSV (can be used for analyzing pixel time series)
    fname_csv = out_dir_scenarios.joinpath(
        f'{vi_name}_{sample_points.name}_pixel_time_series.csv'
    )
    gdf['date'] = gdf['date'].apply(lambda x: x.strftime('%Y-%m-%d'))
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    gdf.drop('geometry', axis=1, inplace=True)
    gdf.to_csv(fname_csv, index=False)


def main(
        vi_dir: Path,
        uncertainty_analysis_dir: Path,
        in_file_aoi: Path,
        out_dir_scenarios: Path,
        out_dir_plots: Path,
        n_scenarios: int,
        vi_name: str,
        sample_points: Path
    ):
    """
    main executable function of this module generating the scenarios for the
    phenological metrics based on the uncertainty in the underlying vegetation
    indices/ parameters calculated in the previous step.
    """

    # search files, join and order them by date
    data_df = get_data_and_uncertainty_files(
        vi_dir=vi_dir,
        uncertainty_analysis_dir=uncertainty_analysis_dir,
        vi_name=vi_name
    )

    # read the data and mask clouds, shadows, and unclassified pixels based on SCL information
    _, ts_stack_dict, orig_vi_dict, unc_dict = read_data_and_uncertainty(
        data_df=data_df,
        vi_name=vi_name,
        in_file_aoi=in_file_aoi,
        out_dir_plots=out_dir_plots
    )

    # actual phenological metrics scenarios
    vegetation_time_series_scenarios(
        ts_stack_dict=ts_stack_dict,
        orig_vi_dict=orig_vi_dict,
        unc_dict=unc_dict,
        n_scenarios=n_scenarios,
        out_dir_scenarios=out_dir_scenarios,
        vi_name=vi_name,
        sample_points=sample_points
    )


if __name__ == '__main__':
    # original Sentinel-2 scenes with vegetation indices
    vi_dir = Path(
        '../S2A_MSIL1C_orig/*.VIs'
    )
    
    # directory with uncertainty analysis results
    uncertainty_analysis_dir = Path(
        '../S2A_MSIL2A_Analysis'
    )

    # extent of study area
    in_file_aoi = Path(
        '../shp/AOI_Esch_EPSG32632.shp'
    )

    # define point sampling locations for visualizing pixel time series
    sample_points = Path('../shp/ZH_Points_2019_EPSG32632_selected-crops.shp')

    # vegetation index to consider
    vi_names = ['NDVI','EVI']

    # directory where to save phenological metrics to
    out_dir_scenarios = Path(f'../S2_TimeSeries_Analysis')
    if not out_dir_scenarios.exists():
        out_dir_scenarios.mkdir()

    # directory where to save plots to
    out_dir_plots = uncertainty_analysis_dir.joinpath('scene_quicklooks')
    if not out_dir_plots.exists():
        out_dir_plots.mkdir()
    
    # number of scenarios to generate
    n_scenarios = 1000

    for vi_name in vi_names:

        out_dir_scenarios_vi = out_dir_scenarios.joinpath(vi_name)
        if not out_dir_scenarios_vi.exists():
            out_dir_scenarios_vi.mkdir()

        main(
            vi_dir=vi_dir,
            uncertainty_analysis_dir=uncertainty_analysis_dir,
            in_file_aoi=in_file_aoi,
            out_dir_scenarios=out_dir_scenarios_vi,
            out_dir_plots=out_dir_plots,
            n_scenarios=n_scenarios,
            vi_name=vi_name,
            sample_points=sample_points
        )
