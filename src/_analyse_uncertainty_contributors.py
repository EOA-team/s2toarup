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
import os
import seaborn as sns
import re

from agrisatpy.core.band import Band
from agrisatpy.core.raster import RasterCollection
from agrisatpy.utils.constants.sentinel2 import s2_gain_factor
from pathlib import Path

from logger import get_logger

s2_gain_factor *= 100

logger = get_logger('_analyze_unc_contrib')

# define Sentinel-2 bands to analyze
s2_band_selection = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12']


def calc_uncertainty(
        scene_dir: Path,
        out_dir: Path
    ):
    """
    Calculates absolute uncertainty per uncertainty type (random or systematic)
    and spectral band as the standard deviation of all scenario members.

    :param scene_dir:
        directory where the scenes and their MC results are stored
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
        scene_path: Path,
        roi_file: Path,
        luc_mapping: dict
    ):
    """
    Extract region of interest (ROI) to display uncertainty values
    for selected land cover types
    """

    rois = gpd.read_file(roi_file)

    band_res = dict.fromkeys(s2_band_selection)
    for band in s2_band_selection:

        # random results
        fname_random = scene_path.joinpath(f'random_uncertainty_{band}.tif')
        random_unc_handler = RasterCollection()
        random_unc_handler.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=fname_random,
            vector_features=rois,
            band_name_dst='random uncertainty'
        )
        random_unc_handler.add_band(
            in_file_vector=rois,
            snap_band=random_unc_handler.bandnames[0],
            attribute_selection=['luc_code']
        )
        gdf_rand = random_unc_handler.to_dataframe()

        # systematic results
        fname_sys = scene_path.joinpath(f'systematic_uncertainty_{band}.tif')
        sys_unc_handler = RasterCollection()
        sys_unc_handler.add_band(
            fname_sys,
            in_file_aoi=rois
        )
        sys_unc_handler.add_bands_from_vector(
            in_file_vector=rois,
            snap_band=sys_unc_handler.bandnames[0],
            attribute_selection=['luc_code']
        )
        gdf_sys = sys_unc_handler.to_dataframe()
        gdf_sys = gdf_sys.rename(
            columns={'sys_unc': f'systematic uncertainty'}
        )

        # join into a single dataframe
        gdf_joined = gdf_rand.join(gdf_sys, lsuffix='_l', rsuffix='_r')
        gdf_joined.drop(['geometry_l', 'luc_code_l'], axis=1, inplace=True)

        # apply land code mapping
        gdf_joined['luc'] = gdf_joined.luc_code_r.apply(
            lambda x, luc_mapping=luc_mapping: luc_mapping[x] 
        )
        band_res[band] = gdf_joined

    # make histogram plot of random and systematic components per ROI
    band_names = [f'systematic uncertainty', f'random uncertainty']
    for roi in gdf_joined['luc'].unique():

        roi_code = [k for k,v in luc_mapping.items() if v == roi]
        # roi size in sqm
        area = rois[rois.luc_code == roi_code[0]].iloc[0].geometry.area
        # convert to ha
        area /= (100*100)

        f, axs = plt.subplots(nrows=2, ncols=5, sharex=False, sharey=True, figsize=(20,12))

        nrow = 0
        ncol = 0
        for band in s2_band_selection:

            df = band_res[band].copy()
            unc_values = df[df.luc == roi].copy()
            unc_values[band_names] *= s2_gain_factor

            sns.histplot(
                data=unc_values[band_names],
                stat='density',
                ax=axs[nrow, ncol],
                legend=True
            )
            if ncol == 0:
                axs[nrow, ncol].set_ylabel('Relative Frequency', fontsize=24)
            if ncol == 2 and nrow == 1:
                axs[nrow, ncol].set_xlabel(
                    r'Absolute Uncertainty (k=1) $\rho_{TOA}$ Reflectance Factor [%]',
                    fontsize=24
                )
            axs[nrow, ncol].set_title(band)

            ncol += 1
            if ncol == 5:
                ncol = 0
                nrow += 1

        f.suptitle(f'{roi}\nArea: {np.round(area,1)}ha', fontsize=24)

        # save figure
        fname = scene_path.joinpath(f'{roi.replace(" ","_")}_histogram_plots.png')
        f.savefig(fname, dpi=150, bbox_inches='tight')


if __name__ == '__main__':

    batches = [str(x) for x in range(1,6)]

    # define regions of interest (different crop types) + forest + water
    

    for batch in batches:

        scenario_dir = Path(f'/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_MSIL1C_RUT-Scenarios/batch_{batch}')
        out_dir = Path('/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_MSIL1C_RUT-Scenarios/L1C_Analysis')
        out_dir.mkdir(exist_ok=True)
        # calc_uncertainty(scenario_dir, out_dir)
        # map_uncertainty(unc_results_dir=out_dir)
    
        plt.style.use('ggplot')
        matplotlib.rc('xtick', labelsize=20) 
        matplotlib.rc('ytick', labelsize=20)

        for scene_path in scenario_dir.glob('*MSIL1C*'):
            # calc_uncertainty(scene_dir=scene_path, out_dir=out_dir)
            analyze_rois(scene_path, roi_file, luc_mapping)

        map_uncertainty(unc_results_dir=out_dir)
    
        # scene_path = out_dir.joinpath('S2A_MSIL1C_20190530T103031_N0207_R108_T32TMT_20190530T123429')
        # roi_file = Path('../shp/ZH_Sampling_Locations_Contributor_Analysis.shp')
        # luc_mapping = {
        #     1: 'Water (Lake)',
        #     2: 'Cumulus Cloud',
        #     3: 'Mixed Forest',
        #     4: 'Arable land (green vegetation)'
        # }
        #
        # analyze_rois(scene_path, roi_file, luc_mapping)
    
    