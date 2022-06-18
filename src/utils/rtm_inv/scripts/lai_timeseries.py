'''
Created on Feb 17, 2022

@author: graflu
'''

import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import gaussian_kde

from agrisatpy.core.raster import RasterCollection
from agrisatpy.core.band import Band

search_expr_lai ='S2*_MSIL2A_*.VIs'
# search GPR LAI and ProSAIL LAI
search_expr_lai_prosail = 'VI_*None_10m_GLAI.tif'
search_expr_lai_gpr = 'VI_*None_10m_LAI.tif'
# search NDVI (serves as reference)
search_expr_ndvi = 'VI_*None_10m_NDVI.tif'
# scene classification layer to mask out clouds and shadows
search_expr_scl = '*_SCL.tiff'

plt.style.use('ggplot')
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16

def extract_data(
        lai_dir: Path,
        parcels: Path,
        output_fname: Path
    ):
    """
    Extracts the NDVI, GPR-derived GLAI and ProSAIL-derived GLAI from a series
    of Sentinel-2 scenes. Uses the scene classification layer to mask out all
    pixels not classified as vegetation (4) or bare soil (5)
    """
    df_list = []
    for scene in lai_dir.rglob(str(search_expr_lai)):

        sensing_date = pd.to_datetime(scene.name.split('_')[2][0:8])

        vi_dir = scene.joinpath('Vegetation_Indices')

        ndvi_file = next(vi_dir.glob(search_expr_ndvi))
        prosail_lai_file = next(vi_dir.glob(search_expr_lai_prosail))
        gpr_lai_file = next(vi_dir.glob(search_expr_lai_gpr))

        scl_dir = scene.joinpath('scene_classification')
        scl_file = next(scl_dir.glob(search_expr_scl))

        # read data from raster into RasterCollection
        collection = RasterCollection()
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=ndvi_file,
            band_idx=1,
            band_name_dst='NDVI',
            vector_features=parcels
        )
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=prosail_lai_file,
            band_idx=1,
            band_name_dst='ProSAIL GLAI',
            vector_features=parcels
        )
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=gpr_lai_file,
            band_idx=1,
            band_name_dst='GPR GLAI',
            vector_features=parcels
        )
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=gpr_lai_file,
            band_idx=2,
            band_name_dst='GPR GLAI SD',
            vector_features=parcels
        )
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=scl_file,
            band_idx=1,
            band_name_dst='SCL',
            vector_features=parcels
        )
        collection.add_band(
            band_constructor=Band.from_vector,
            vector_features=parcels,
            geo_info=collection['NDVI'].geo_info,
            band_name_src='crop_code',
            band_name_dst='crop_code',
            snap_bounds=collection['NDVI'].bounds
        )
        collection.mask(
            mask=collection['NDVI'].values.mask,
            bands_to_mask=['crop_code'],
            inplace=True
        )

        # convert to GeoDataFrame
        gdf = collection.to_dataframe()
        # keep SCL classes 4 and 5, only
        clean_pix_gdf = gdf[gdf['SCL'].isin([4,5])].copy()
        clean_pix_gdf['date'] = sensing_date
        clean_pix_gdf['date'] = pd.to_datetime(clean_pix_gdf['date'])
        df_list.append(clean_pix_gdf)

        print(f'Extracted scene {scene.name}')

    # concat dataframe, clean it and save it to csv
    complete_gdf = pd.concat(df_list)
    complete_gdf['x'] = complete_gdf.geometry.x
    complete_gdf['y'] = complete_gdf.geometry.y
    complete_gdf.drop('geometry', axis=1, inplace=True)
    complete_gdf.to_csv(output_fname, index=False)

def plot_crop_timeseries(
        output_fname: Path,
        parcels: Path
    ):
    """
    Plots time series of NDVI, GRP-LAI and ProSAIL LAI for all
    crop types available
    """
    # get mapping of crop code and crop type
    gdf = gpd.read_file(parcels)
    crop_codes = gdf['crop_code'].unique()
    mapping = dict.fromkeys(crop_codes)
    for crop_code in crop_codes:
        mapping[crop_code] = gdf[gdf.crop_code == crop_code]['crop_type'].iloc[0]

    # loop over crops and plot the time series of all pixels of the crop type
    # calculate the mean, stddev and min-max spread for each trait (NDVI, GPR-LAI,
    # ProSAIL-LAI)
    ts_df = pd.read_csv(output_fname)
    col_selection = ['date', 'NDVI', 'ProSAIL GLAI', 'GPR GLAI', 'GPR GLAI SD']

    for crop_code in crop_codes:

        crop_type = mapping[crop_code]
        ts_df_crop = ts_df[ts_df.crop_code == crop_code].copy()
        ts_df_crop_agg = ts_df_crop[col_selection].groupby(
            by='date'
        ).agg(['mean','std','min','max','count'])
        ts_df_crop_agg['date'] = pd.to_datetime(ts_df_crop_agg.index)

        max_pixel_num = ts_df_crop_agg['NDVI', 'count'].max()
        # get number of available pixel per date in percent
        pixel_date_count = ts_df_crop_agg['NDVI', 'count'] / max_pixel_num * 100

        # plot time series
        f, axes = plt.subplots(
            nrows=2,
            figsize=(14,12),
            gridspec_kw={'height_ratios': [3, 1]}
        )

        # NDVI
        ax_second = axes[0].twinx()
        ax_second.plot(
            ts_df_crop_agg.date,
            ts_df_crop_agg['NDVI', 'mean'],
            '-o',
            label='Mean NDVI',
            color='blue'
        )
        ax_second.fill_between(
            ts_df_crop_agg.date,
            ts_df_crop_agg['NDVI', 'mean'] - ts_df_crop_agg['NDVI', 'std'],
            ts_df_crop_agg['NDVI', 'mean'] + ts_df_crop_agg['NDVI', 'std'],
            color='lightblue',
            label=r'$\pm$ 1 Stddev NDVI',
            alpha=0.5
        )
        ax_second.set_ylabel('NDVI [-]', color='blue', fontsize=20)
        ax_second.set_ylim(0,1)

        # both LAI products
        axes[0].plot(
            ts_df_crop_agg.date,
            ts_df_crop_agg['ProSAIL GLAI', 'mean'],
            '-o',
            label='Mean ProSAIL GLAI',
            color='green'
        )
        axes[0].fill_between(
            ts_df_crop_agg.date,
            ts_df_crop_agg['ProSAIL GLAI', 'mean'] - ts_df_crop_agg['ProSAIL GLAI', 'std'],
            ts_df_crop_agg['ProSAIL GLAI', 'mean'] + ts_df_crop_agg['ProSAIL GLAI', 'std'],
            color='green',
            label=r'$\pm$ 1 Stddev ProSAIL GLAI',
            alpha=0.4
        )
        axes[0].plot(
            ts_df_crop_agg.date,
            ts_df_crop_agg['GPR GLAI', 'mean'],
            '-o',
            label='Mean GPR GLAI',
            color='darkgreen'
        )
        axes[0].fill_between(
            ts_df_crop_agg.date,
            ts_df_crop_agg['GPR GLAI', 'mean'] - ts_df_crop_agg['GPR GLAI', 'std'],
            ts_df_crop_agg['GPR GLAI', 'mean'] + ts_df_crop_agg['GPR GLAI', 'std'],
            color='darkgreen',
            label=r'$\pm$ 1 Stddev GPR GLAI',
            alpha=0.5
        )
        axes[0].set_ylabel(r'Green LAI [$m^2$/$m^2$]', fontsize=20, color='g')
        axes[0].set_ylim(0,7)

        axes[1].bar(ts_df_crop_agg.date, pixel_date_count)
        axes[1].set_ylabel('Valid Pixels [%]', fontsize=20)

        f.suptitle(f'{crop_type}\n(#S2 10m Pixels: {max_pixel_num})', size=24)

        lines_labels = [ax.get_legend_handles_labels() for ax in f.axes]
        lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        f.legend(lines, labels, loc='lower center', ncol=3, fontsize=14)

        fname_fig = output_fname.parent.joinpath(f'{crop_type}_NDVI_GLAI.png')
        f.savefig(fname_fig, dpi=200, bbox_inches='tight')
        plt.close(f)

        # scatter plot with kernel density (GPR-LAI vs ProSAIL LAI)
        gpr_lai = ts_df_crop['GPR GLAI'].values
        prosail_lai = ts_df_crop['ProSAIL GLAI'].values
        prosail_lai = prosail_lai[~np.isnan(gpr_lai)]
        gpr_lai = gpr_lai[~np.isnan(gpr_lai)]

        xy = np.vstack([gpr_lai, prosail_lai])
        # get kernel density
        z = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        
        x, y, z = gpr_lai[idx], prosail_lai[idx], z[idx]

        f = plt.figure(figsize=(14,14))
        scatter_axes = plt.subplot2grid((3, 3), (1, 0), rowspan=2, colspan=2)
        x_hist_axes = plt.subplot2grid((3, 3), (0, 0), colspan=2,
                                       sharex=scatter_axes)
        y_hist_axes = plt.subplot2grid((3, 3), (1, 2), rowspan=2,
                                       sharey=scatter_axes)
        
        scatter_axes.scatter(x, y, c=z, s=50)
        # scatter_axes.scatter(gpr_lai, prosail_lai)
        scatter_axes.set_xlim(0,8)
        scatter_axes.set_ylim(0,8)
        scatter_axes.set_xlabel(r'GPR GLAI [$m^2$/$m^2$]')
        scatter_axes.set_ylabel(r'ProSAIL GLAI [$m^2$/$m^2$]')
        
        sns.histplot(
            data=gpr_lai,
            stat='count',
            ax=x_hist_axes,
            label='GPR GLAI [$m^2$/$m^2$]'
        )
        x_hist_axes.legend()
        sns.histplot(
            y=prosail_lai,
            stat='count',
            ax=y_hist_axes,
            label='ProSAIL GLAI [$m^2$/$m^2$]'
        )
        y_hist_axes.legend()

        fname_fig = output_fname.parent.joinpath(f'{crop_type}_GPR-ProSAIL_GLAI.png')
        f.savefig(fname_fig, dpi=200, bbox_inches='tight')
        plt.close(f)

if __name__ == '__main__':

    # set input and output paths
    parcels = Path('/mnt/ides/Lukas/software/scripts_paper_uncertainty/shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    lai_dir = Path('/home/graflu/Documents/uncertainty/S2_MSIL1C_orig')
    output_fname = Path('/home/graflu/Documents/uncertainty/NDVI_GPR-ProSAIL_GLAI_values_crops.csv')
    # output_fname = Path('/mnt/ides/Lukas/04_Work/GPR_LAI/NDVI_GPR-ProSAIL_GLAI_values_crops.csv')

    # extract_data(lai_dir, parcels, output_fname)
    plot_crop_timeseries(output_fname, parcels)
