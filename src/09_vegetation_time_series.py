'''
This script visualizes the temporal development of the vegetation
indices and functional canopy traits together with their uncertainty.
'''

import matplotlib
import pandas as pd
import numpy as np

from pathlib import Path
from datetime import date
from typing import Tuple
from typing import List
from typing import Dict
from copy import deepcopy
import geopandas as gpd
import matplotlib.pyplot as plt

from _find_datasets import get_data_and_uncertainty_files
from agrisatpy.core.sensors import Sentinel2
from agrisatpy.core.raster import RasterCollection

from logger import get_logger
from agrisatpy.core.band import Band


# setup logger -> will write log file to the /../log directory
logger = get_logger('l4_phenology')

plt.style.use('ggplot')
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)

# define SCL classes to mask
# see https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-2-msi/level-2a/algorithm
# for a legend of the class labels
scl_mask_values = [1,2,3,7,8,9,10,11]

def read_data_and_uncertainty(
        data_df: pd.DataFrame,
        parcels: Path,
        vi_name: str
    ) -> Tuple[pd.DataFrame, List[Sentinel2]]:
    """
    This function reads the selected vegetation index (computed in step 7)
    and the associated standard (absolute) uncertainty and stores the read
    data as new columns in data_df.

    :param data_df:
        dataframe returned from ``get_data_and_uncertainty``
    :param parcels:
        shapefile defining the crop parcels
    :param vi_name:
        name of the vegetation index or parameter to process (e.g., NDVI)
    :return:
        input dataframe with read data + list of read vegetation data and
        their absolute uncertainty per image acquisition date
    """

    # loop over datasets (single acquisition dates) and read the data
    # (vegetation index/ parameter + standard uncertainty)
    update_list = []
    ts_stack_dict = {}

    for _, record in data_df.iterrows():

        # read VI data
        vi_file = Path(record.filename_veg_par)
        collection = Sentinel2()
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=vi_file,
            band_idx=1,
            band_name_dst=vi_name,
            vector_features=parcels
        )
        # read vegetation parameter/index absolute uncertainty band
        vi_unc_file = Path(record.filename_unc)
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=vi_unc_file,
            band_idx=4,
            vector_features=parcels,
            band_name_dst=f'{vi_name}_unc'
        )

        # read S2 blue and red bands plus SCL file to mask out clouds and shadows
        scl_file = next(
            Path(record.filename_bandstack).parent.rglob('scene_classification/*.tiff')
        )
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=scl_file,
            band_idx=1,
            band_name_dst='SCL',
            vector_features=parcels
        )
        # collection.add_band(
        #     band_constructor=Band.from_rasterio,
        #     fpath_raster=Path(record.filename_bandstack),
        #     band_idx=1,
        #     band_name_dst='blue',
        #     vector_features=parcels,
        # )
        # collection.add_band(
        #     band_constructor=Band.from_rasterio,
        #     fpath_raster=Path(record.filename_bandstack),
        #     band_idx=3,
        #     band_name_dst='red',
        #     vector_features=parcels,
        # )
        # add crop codes
        collection.add_band(
            band_constructor=Band.from_vector,
            vector_features=parcels,
            geo_info=collection[vi_name].geo_info,
            band_name_src='crop_code',
            band_name_dst='crop_code',
            snap_bounds=collection[vi_name].bounds
        )
        collection.mask(
            mask=collection[vi_name].values.mask,
            bands_to_mask=['crop_code'],
            inplace=True
        )
        # mask clouds, shadows and snow
        collection.mask_clouds_and_shadows(
            bands_to_mask=collection.band_names,
            cloud_classes=[1,2,3,6,7,8,9,10,11],
            inplace=True
        )
        cloud_coverage = collection.get_cloudy_pixel_percentage(cloud_classes=[1,2,3,6,7,8,9,10,11])

        # save masked files to disk (thus, the reference run will also have the correct
        # input)
        vi_masked = record.filename_veg_par
        vi_masked = Path(vi_masked.replace('.tif', '_scl-filtered.tif'))
        collection.to_rasterio(
            fpath_raster=vi_masked,
            band_selection=[vi_name]
        )
        vi_unc_masked = record.filename_unc
        vi_unc_masked = Path(vi_unc_masked.replace('.tif', '_scl-filtered.tif'))
        collection.to_rasterio(
            fpath_raster=vi_unc_masked,
            band_selection=[f'{vi_name}_unc']
        )
        ts_stack_dict[record.date] = collection

        update_list.append({
            'date': record.date,
            'cloudy_pixel_percentage': cloud_coverage,
            'filename_veg_par_scl': vi_masked,
            'filename_unc_scl': vi_unc_masked
        })

    # added update columns to input data frame
    update_df = pd.DataFrame(update_list)
    merged = pd.merge(data_df, update_df, on='date')
    return merged, ts_stack_dict


def percentile(n):
    def percentile_(x):
        return np.percentile(x, n)
    percentile_.__name__ = 'percentile_%s' % n
    return percentile_


def extract_uncertainty_crops(
        ts_stack_dict: Dict[date, RasterCollection],
        file_df: pd.DataFrame,
        vi_name: str,
        out_dir: Path,
    ) -> None:
    """
    Extracts the uncertainty information for the field parcel polygons alongside
    the time series of the vegetation index/parameter. The data is saved to CSV file
    as backup.

    :param ts_stack_dict:
        dictionary with ``RasterCollection`` instances holding the vegetation parameter
        as one band and its uncertainty as another band per scene date
    :param file_df:
        pandas dataframe with links to the original files (parameter values + uncertainty)
        for carrying out the reference run
    :param vi_name:
        name of the vegetation index/ parameter to process
    :param out_dir:
        directory where to save the CSV files to
    :param ymin:
        minimum y value to be used for plotting the vegetation index/parameter values
    :param ymax:
        maximum y value to be used for plotting the vegetation index/parameter values
    """

    # loop over scenes, rasterize the shapefile with the crops and convert it to
    # geodataframe to allow for crop-type specific analysis
    gdf_list = []
    for idx, record in file_df.iterrows():

        # read original (reference) data for the field polygons
        _date = record.date
        scene_handler = deepcopy(ts_stack_dict[_date])

        # convert to geodataframe and save them as CSVs (backup)
        gdf = scene_handler.to_dataframe()
        # drop NaNs (occur because NDVI was filtered)
        gdf = gdf.dropna()
        gdf['date'] = record.date
        gdf_list.append(gdf)
        
        fname_csv = out_dir.joinpath(f'{vi_name}_{idx+1}_crops.csv')
        scene_handler = None
        gdf.to_csv(fname_csv, index=False)

        logger.info(f'Read scene data from {_date} ({idx+1}/{file_df.shape[0]})')

    # concat dataframes
    gdf = pd.concat(gdf_list)
    gdf.to_csv(out_dir.joinpath(f'{vi_name}_crops.csv'), index=False)

def plot_uncertainty_tim_series(
        sample_polygons: Path,
        vi_data_fpath: Path,
        out_dir: Path,
        ymin: float,
        ymax: float
    ):
    """
    Plots the VI and its uncertainty

    :param sample_polygon:
        shapefile with field parcel polygons
    :param vi_data_fpath:
        path to the CSV file with extracted VI data for the different crop types
        and field parcels
    :param out_dir:
        directory where to save the CSV files to
    :param ymin:
        minimum y value to be used for plotting the vegetation index/parameter values
    :param ymax:
        maximum y value to be used for plotting the vegetation index/parameter values
    """

    # add crop name from original shape file
    gdf_polys = gpd.read_file(sample_polygons)

    # loop over crop types and plot their VI curve and its uncertainty over time
    gdf = pd.read_csv(vi_data_fpath)
    crop_codes = np.unique(gdf.crop_code)

    for crop_code in crop_codes:

        # get crop name
        crop_name = gdf_polys[gdf_polys.crop_code == int(crop_code)]['crop_type'].iloc[0]
        crop_gdf = gdf[gdf.crop_code == crop_code].copy()

        # add relative uncertainty by computing the ratio between the absolute uncertainty
        # and the original index value (multiplied by 100 to get % relative uncertainty)
        crop_gdf[f'{vi_name}_rel_unc'] =  crop_gdf[f'{vi_name}_unc'] / crop_gdf[vi_name] * 100

        # aggregate by date
        col_selection = ['date', vi_name, f'{vi_name}_unc', f'{vi_name}_rel_unc']
        crop_gdf_grouped = crop_gdf[col_selection].groupby('date').agg(
            ['mean', 'min', 'max', 'std', percentile(5), percentile(95), 'count']
        )
        crop_gdf_grouped['date'] = pd.to_datetime(crop_gdf_grouped.index)
        pixel_count = crop_gdf_grouped[vi_name , 'count'].max()

        # plot
        fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(15,20))

        # plot original time series showing spread between pixels (not scenarios!)
        ax1.plot(crop_gdf_grouped.date, crop_gdf_grouped[vi_name, 'mean'], color='blue', linewidth=3)
        # central 90% of values
        ax1.fill_between(
            x=crop_gdf_grouped.date,
            y1=crop_gdf_grouped[vi_name, 'percentile_5'],
            y2=crop_gdf_grouped[vi_name, 'percentile_95'],
            color='orange',
            alpha=0.4
        )
        # standard deviation around the mean.
        y1 = crop_gdf_grouped[vi_name, 'mean'] - crop_gdf_grouped[vi_name, 'std']
        # standard deviations smaller ymin are not possible
        y1[y1 < ymin] = ymin
        y2 = crop_gdf_grouped[vi_name, 'mean'] + crop_gdf_grouped[vi_name, 'std']
        # standard deviations larger then ymax are not possible
        y2[y2 > ymax] = ymax
        ax1.fill_between(
            x=crop_gdf_grouped.index,
            y1=y1,
            y2=y2,
            color='red',
            alpha=0.45
        )
        ax1.set_title(vi_name, fontsize=20)
        ax1.set_ylabel(f'{vi_name} [-]', fontsize=24)
        ax1.set_ylim(ymin, ymax)

        # plot uncertainties (absolute)
        unc_name = vi_name + '_unc'
        ax2.plot(crop_gdf_grouped.date, crop_gdf_grouped[unc_name, 'mean'], label='Mean', color='blue', linewidth=3)
        ax2.fill_between(
            x=crop_gdf_grouped.index,
            y1=crop_gdf_grouped[unc_name, 'percentile_5'],
            y2=crop_gdf_grouped[unc_name, 'percentile_95'],
            color='orange',
            alpha=0.4
        )

        # standard deviation around the mean. Values below ymin are not possible
        y1 = crop_gdf_grouped[unc_name, 'mean'] - crop_gdf_grouped[unc_name, 'std']
        y1[y1 < ymin] = ymin
        # the same applies to the upper bound
        y2 = crop_gdf_grouped[unc_name, 'mean'] + crop_gdf_grouped[unc_name, 'std']
        y2[y2 > ymax] = ymax
        ax2.fill_between(
            x=crop_gdf_grouped.index,
            y1=y1,
            y2=y2,
            color='red',
            alpha=0.45
        )
        ax2.set_title(f'Uncertainty in {vi_name} (k=1)', fontsize=20)
        ax2.set_ylabel(f'Absolute Uncertainty [-]', fontsize=24)

        # relative uncertainties
        unc_name = vi_name + '_rel_unc'
        ax3.plot(crop_gdf_grouped.date,crop_gdf_grouped[unc_name, 'mean'], label='Mean', color='blue', linewidth=3)
        ax3.fill_between(
            x=crop_gdf_grouped.index,
            y1=crop_gdf_grouped[unc_name, 'percentile_5'],
            y2=crop_gdf_grouped[unc_name, 'percentile_95'],
            color='orange',
            alpha=0.4,
            label='5-95% Quantile Spread'
        )

        # absolute uncertainties
        # negative values are not possible for relative uncertainties
        y1 = crop_gdf_grouped[unc_name, 'mean']-crop_gdf_grouped[unc_name, 'std']
        y1[y1 < 0] = 0
        y2 = crop_gdf_grouped[unc_name, 'mean']+crop_gdf_grouped[unc_name, 'std']
        ax3.fill_between(
            x=crop_gdf_grouped.index,
            y1=y1,
            y2=y2,
            color='red',
            alpha=0.45,
            label=r'$\pm$ 1 Stddev'
        )
        ax3.set_ylabel(f'Relative Uncertainty [%]', fontsize=24)
       
        ax3.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3), fontsize=20, ncol=3)

        fig.suptitle(f'Time Series of {crop_name} (Pixels: {pixel_count})', fontsize=26)
        fname = f'{vi_name}_{crop_name}_all-pixel-timeseries.png'
        fig.savefig(out_dir.joinpath(fname), dpi=300, bbox_inches='tight')
        plt.close(fig)

        # plot histograms of uncertainty values
        fig = plt.figure(figsize=(7,7))
        ax = fig.add_subplot(111)
        crop_gdf[f'{vi_name}_unc'].hist(ax=ax, bins=100, density=True)
        ax.set_title(f'{crop_name}', fontsize=20)
        ax.set_xlabel(f'{vi_name} Absolute Uncertainty (k=1)', fontsize=24)
        ax.set_ylabel('Relative Frequency [%]', fontsize=24)
        fname = f'{vi_name}_{crop_name}_abs-unc-histogram.png'
        fig.savefig(out_dir.joinpath(fname), dpi=300, bbox_inches='tight')
        plt.close(fig)


def main(
        vi_dir: Path,
        uncertainty_analysis_dir: Path,
        out_dir_scenarios: Path,
        vi_name: str,
        sample_polygons: Path,
        ymin: float,
        ymax: float
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
        parcels=sample_polygons
    )
    
    # check uncertainty in different crop types over time
    extract_uncertainty_crops(
        file_df=file_df,
        ts_stack_dict=ts_stack_dict,
        out_dir=out_dir_scenarios,
        vi_name=vi_name,
    )

    # visualize it
    fname_csv = out_dir_scenarios.joinpath(
        f'{vi_name}_crops.csv'
    )
    plot_uncertainty_tim_series(
        sample_polygons=sample_polygons,
        vi_data_fpath=fname_csv,
        out_dir=out_dir_scenarios,
        ymin=ymin,
        ymax=ymax
    )

if __name__ == '__main__':
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

    # vegetation index to consider
    vi_names = ['GLAI', 'NDVI','EVI']
    ymins = {'NDVI': 0, 'EVI': 0, 'GLAI': 0}
    ymaxs = {'NDVI': 1, 'EVI': 1, 'GLAI': 7}

    # directory where to save phenological metrics to
    out_dir_scenarios = Path(f'../S2_TimeSeries_Analysis')
    if not out_dir_scenarios.exists():
        out_dir_scenarios.mkdir()

    for vi_name in vi_names:

        out_dir_scenarios_vi = out_dir_scenarios.joinpath(vi_name)
        if not out_dir_scenarios_vi.exists():
            out_dir_scenarios_vi.mkdir()

        main(
            vi_dir=vi_dir,
            uncertainty_analysis_dir=uncertainty_analysis_dir,
            out_dir_scenarios=out_dir_scenarios_vi,
            vi_name=vi_name,
            sample_polygons=sample_polygons,
            ymin=ymins[vi_name],
            ymax=ymaxs[vi_name]
        )
