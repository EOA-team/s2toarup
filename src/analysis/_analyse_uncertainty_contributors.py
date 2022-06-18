'''
This script is a simplified version of 08_analyze_l1c_l2a_l3_scenarios.py
that focuses on quantifying the impact of random and systematic uncertainty
contributors on Sentinel-2 L1C TOA radiometric data.

It takes the (optionally) generated samples from 04_generate_l1c_scenarios.py
that denote the systematic and random uncertainty parts (named as
correlated and uncorrelated contributors, respectively).
'''

import geopandas as gpd
import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
import re

from eodal.core.band import Band
from eodal.core.raster import RasterCollection
from eodal.core.sensors import Sentinel2
from eodal.utils.constants.sentinel2 import s2_gain_factor
from datetime import datetime
from pathlib import Path
from typing import List

from utils.logger import get_logger

s2_gain_factor *= 100

logger = get_logger('_analyze_unc_contrib')

# define Sentinel-2 bands to analyze
s2_band_selection = ['B02', 'B03', 'B04', 'B08']


def calc_uncertainty(
        scene_dir: Path,
        orig_scene: Path,
        out_dir: Path
    ):
    """
    Calculates absolute uncertainty per uncertainty type (random or systematic)
    and spectral band as the standard deviation of all scenario members.

    :param scene_dir:
        directory where the scenes and their MC results are stored
    :param orig_scene
        original S2 scene to use for calculating the relative uncertainty
    :param out_dir:
        directory where to store the results (create sub-directories for the
        scenes analyzed)
    """

    # get scenarios (numbered from 1 to N)
    m = re.compile('\d')
    scenarios = [f for f in os.listdir(scene_dir) if m.search(f)]
    scenarios = [Path(scene_dir).joinpath(x) for x in scenarios]

    out_dir_scene = out_dir.joinpath(Path(scene_dir).name)
    out_dir_scene.mkdir(exist_ok=True)

    logger.info(f'Working on {scene_dir}')

    # loop over spectral bands (uncertainty is always calculated per band)
    # read only the extent of the study area
    fname_rand = Path(scenarios[0]).joinpath(
                f'uncorrelated_contributors_sample_B04.jp2'
            )
    # read data into memory
    rand_handler = RasterCollection()
    rand_handler.add_band(
        band_constructor=Band.from_rasterio,
        fpath_raster=fname_rand,
        band_name_dst='random'
    )
    bdf = gpd.GeoDataFrame(geometry=[rand_handler['random'].bounds])
    bdf.set_crs(crs=rand_handler['random'].crs, inplace=True)
    s2_orig = Sentinel2.from_safe(
        orig_scene,
        vector_features=bdf,
        band_selection=s2_band_selection,
        apply_scaling=False
    )
    # plot FCIR quicklook and save it
    fcir = s2_orig.plot_multiple_bands(band_selection=['B08','B04','B03'])
    fcir.savefig(out_dir.joinpath(orig_scene.name.split('.')[0] + '_fcir.png'))

    for band in s2_band_selection:
        # loop over scenarios
        for idx, scenario in enumerate(scenarios):

            # get random and systematic results
            fname_rand = Path(scenario).joinpath(
                f'uncorrelated_contributors_sample_{band}.jp2'
            )
            fname_sys = Path(scenario).joinpath(
                f'correlated_contributors_sample_{band}.jp2'
            )

            # read data into memory
            rand_handler = RasterCollection()
            rand_handler.add_band(
                band_constructor=Band.from_rasterio,
                fpath_raster=fname_rand,
                band_name_dst='random'
            )
            rand_data = rand_handler.get_band('random').values

            sys_handler = RasterCollection()
            sys_handler.add_band(
                band_constructor=Band.from_rasterio,
                fpath_raster=fname_sys,
                band_name_dst='systematic'
            )
            sys_data = sys_handler.get_band('systematic').values

            if idx == 0:
                rand_array = np.zeros(
                    shape=(len(scenarios), rand_data.shape[0], rand_data.shape[1]),
                    dtype=rand_data.dtype
                )
                sys_array = np.zeros(
                    shape=(len(scenarios), sys_data.shape[0], sys_data.shape[1]),
                    dtype=sys_data.dtype
                )

            rand_array[idx,:,:] = rand_data
            sys_array[idx,:,:] = sys_data

        # calculate standard uncertainty
        rand_unc = np.nanstd(rand_array, axis=0)
        sys_unc = np.nanstd(sys_array, axis=0)
        # convert to relative uncertainties
        rand_unc = rand_unc / s2_orig[band].values * 100.
        sys_unc = sys_unc / s2_orig[band].values * 100.

        fname_sys = out_dir_scene.joinpath(f'systematic_uncertainty_{band}.tif')
        fname_rand = out_dir_scene.joinpath(f'random_uncertainty_{band}.tif')

        sys_handler.add_band(
            band_constructor=Band,
            values=sys_unc,
            band_name='sys_unc',
            geo_info=sys_handler['systematic'].geo_info
        )
        rand_handler.add_band(
            band_constructor=Band,
            values=rand_unc,
            band_name='rand_unc',
            geo_info=rand_handler['random'].geo_info
        )

        # save as raster
        sys_handler.to_rasterio(
            fpath_raster=fname_sys,
            band_selection=['sys_unc']
        )
        rand_handler.to_rasterio(
            fpath_raster=fname_rand,
            band_selection=['rand_unc']
        )

    logger.info(f'Finished {scene_dir}')


def map_uncertainty(unc_results_dir: Path):
    """
    Generates plots of system and random uncertainty components per spectral
    band

    :param unc_results_dir:
        directory where files with derived uncertainty are located
    """

    # find scenes and loop over them to plot uncertainty components
    scenes = glob.glob(unc_results_dir.joinpath('S2*_MSIL1C*').as_posix())

    for scene in scenes:

        for band in s2_band_selection:

            # random results
            fname_random = Path(scene).joinpath(f'random_uncertainty_{band}.tif')
            random_unc_handler = RasterCollection()
            random_unc_handler.add_band(
                band_constructor=Band.from_rasterio,
                fpath_raster=fname_random,
                band_name_dst=f'Random Radiometric Uncertainty {band}'
            )
            random_unc_stats = random_unc_handler.band_summaries()

            # systematic results
            fname_sys = Path(scene).joinpath(f'systematic_uncertainty_{band}.tif')
            sys_unc_handler = RasterCollection()
            sys_unc_handler.add_band(
                band_constructor=Band.from_rasterio,
                fpath_raster=fname_sys,
                band_name_dst=f'Systematic Radiometric Uncertainty {band}'
            )
            sys_unc_stats = sys_unc_handler.band_summaries()

            # determine vmin and vmax for plotting
            minmin = np.min(
                [random_unc_stats['nanmin'].iloc[0], sys_unc_stats['nanmin'].iloc[0]]
            )
            maxmax = np.max(
                [random_unc_stats['nanmax'].iloc[0], sys_unc_stats['nanmax'].iloc[0]]
            )

            fig_random_unc = random_unc_handler.plot_band(
                band_name=f'Random Radiometric Uncertainty {band}',
                colormap='terrain',
                vmin=minmin,
                vmax=maxmax,
                colorbar_label='Absolute Uncertainty (k=1)'
            )
            # save to folder with scene uncertainty results
            fname_fig_random_unc = Path(scene).joinpath(f'random_uncertainty_{band}.png')
            fig_random_unc.savefig(fname_fig_random_unc, dpi=150, bbox_inches='tight')
            plt.close(fig_random_unc)

            fig_sys_unc = sys_unc_handler.plot_band(
                band_name=f'Systematic Radiometric Uncertainty {band}',
                colormap='terrain',
                vmin=minmin,
                vmax=maxmax,
                colorbar_label='Absolute Uncertainty (k=1)'
            )
            fname_fig_sys_unc = Path(scene).joinpath(f'systematic_uncertainty_{band}.png')
            fig_sys_unc.savefig(fname_fig_sys_unc, dpi=150, bbox_inches='tight')
            plt.close(fig_sys_unc)


def analyze_rois(
        unc_results_dir: Path,
        rois: gpd.GeoDataFrame,
        black_list: List[str],
        band_selection: List[str] = ['B04', 'B08'],
    ):
    """
    Extract region of interest (ROI) to display uncertainty values
    for selected land cover types
    """
    # find scenes and loop over them to plot uncertainty components
    scenes = glob.glob(unc_results_dir.joinpath('S2*_MSIL1C*').as_posix())
    scenes = [x for x in scenes if not x.endswith('png')]

    # convert to DataFrame
    scenes = pd.DataFrame({'scenes': scenes})
    # get dates for x axis labels and spacing
    scenes['date'] = scenes['scenes'].apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[2][0:8], '%Y%m%d').date()
    )
    # exclude black-listed dates
    black_list_dates = [datetime.strptime(x, '%Y%m%d').date() for x in black_list]
    scenes = scenes[~scenes.date.isin(black_list_dates)]
    # order by date
    scenes.sort_values(by='date', inplace=True)

    # open figure for plotting uncertainty histograms over time
    crops = rois.crop_type.unique()
    n_crops = len(crops)

    for crop in crops:
        df_list = []
        for _, record in scenes.iterrows():
            for band in band_selection:
                crop_gdf = gdf[gdf.crop_type == crop].copy()
                raster = RasterCollection()
                # random contributors
                raster.add_band(
                    band_constructor=Band.from_rasterio,
                    fpath_raster=Path(record.scenes).joinpath(f'random_uncertainty_{band}.tif'),
                    vector_features=crop_gdf,
                    band_name_dst=f'Random Radiometric Uncertainty {band}'
                )
                df_1 = raster.to_dataframe(
                    band_selection=[f'Random Radiometric Uncertainty {band}']
                )
                df_1 = df_1.rename(
                    columns={f'Random Radiometric Uncertainty {band}': 'uncertainty'}
                )
                df_1['Uncertainty Type'] = 'random'
                df_1['band'] = band
                # systematic contributors
                raster.add_band(
                    band_constructor=Band.from_rasterio,
                    fpath_raster=Path(record.scenes).joinpath(f'systematic_uncertainty_{band}.tif'),
                    vector_features=crop_gdf,
                    band_name_dst=f'Systematic Radiometric Uncertainty {band}'
                )
                df_2 = raster.to_dataframe(
                    band_selection=[f'Systematic Radiometric Uncertainty {band}']
                )
                df_2 = df_2.rename(
                    columns={f'Systematic Radiometric Uncertainty {band}': 'uncertainty'}
                )
                df_2['Uncertainty Type'] = 'systematic'
                df_2['band'] = band
                
                # convert to dataframe
                df = df_1.append(df_2)
                df['date'] = record.date
                df_list.append(df)

        large_df = pd.concat(df_list)

        # violin plots
        palette = ['blue', 'red']
        f, ax = plt.subplots(nrows=1, ncols=2, figsize=(30,10), sharey=True)
        sns.violinplot(
            x='date',
            y='uncertainty',
            hue='Uncertainty Type',
            data=large_df[large_df.band == 'B04'],
            split=True,
            ax=ax[0],
            scale='count',
            scale_hue=False,
            saturation=.75,
            inner=None,
            palette=palette
        )
        ax[0].set_ylabel('Relative Radiometric Uncertainty [%]\n(k=1)', fontsize=16)
        ax[0].set_title('Sentinel-2 B04 (red)')
        labels = ax[0].xaxis.get_ticklabels()
        ax[0].set_xticklabels(labels=labels,rotation=90)
        ax[0].set_ylim(0,5)

        sns.violinplot(
            x='date',
            y='uncertainty',
            hue='Uncertainty Type',
            data=large_df[large_df.band == 'B08'],
            split=True,
            ax=ax[1],
            scale='count',
            scale_hue=False,
            saturation=.75,
            inner=None,
            palette=palette
        )
        ax[1].set_ylabel('')
        ax[1].set_title('Sentinel-2 B08 (near-infrared)')
        ax[1].yaxis.set_label_position('right')
        # ax[1].yaxis.tick_right()
        labels = ax[1].xaxis.get_ticklabels()
        ax[1].set_xticklabels(labels=labels,rotation=90)

        crop_name = crop
        if crop_name == 'Canola': crop_name = 'Rapeseed'
        if crop_name == 'Permament Grasland': crop_name = 'Permanent Grassland'
        if crop_name == 'Extensively Used Grasland': crop_name = 'Extensively Used Grassland'
        if crop_name == 'Corn': crop_name = 'Grain Maize'

        f.suptitle(
            f'{crop_name} (Pixels:  {df.groupby(by="date").agg("count")["uncertainty"].values[0]})',
            fontsize=20
        )
        plt.tight_layout()

        fname = unc_results_dir.joinpath(f'{crop_name}_unc_violin-plots.png')
        f.savefig(fname, bbox_inches='tight', dpi=300)
        logger.info(f'Plotted Uncertainty for {crop}')


if __name__ == '__main__':

    # define regions of interest (different crop types) + forest + settlement
    gdf = gpd.read_file('../../shp/areas_of_interest_uncertainty_contributors_dissolved.gpkg')

    orig_scenes_dir = Path('../../S2_MSIL1C_orig')

    scenario_dir = Path('../../S2_MSIL1C_RUT-Scenarios')
    out_dir = Path('../../S2_MSIL1C_RUT-Scenarios/L1C_Analysis')
    out_dir.mkdir(exist_ok=True)

    plt.style.use('ggplot')
    matplotlib.rc('xtick', labelsize=14) 
    matplotlib.rc('ytick', labelsize=14)

    for scene_path in scenario_dir.glob('*MSIL1C*'):
        # find the corresponding original S2 scene for relative uncertainty calculation
        scene_date = scene_path.name.split('_')[2]
        orig_scene = next(orig_scenes_dir.glob(f'*MSIL1C_{scene_date}*.SAFE'))
        calc_uncertainty(scene_dir=scene_path, orig_scene=orig_scene, out_dir=out_dir)

    # black-list scenes with clouds or snow
    black_list = ['20190214', '20190216', '20190530','20190716']

    analyze_rois(unc_results_dir=out_dir, rois=gdf, black_list=black_list)
