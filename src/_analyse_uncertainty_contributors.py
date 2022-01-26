'''
This script is a simplified version of 08_analyze_l1c_l2a_l3_scenarios.py
that focusses on quantifying the impact of random and systematic uncertainty
contributors on Sentinel-2 L1C TOA radiometric data.

It takes the (optionally) generated samples from 04_generate_l1c_scenarios.py
that denote the systematic and random uncertainty parts (named as
correlated and uncorrelated contributors, respectively).
'''

import glob
import numpy as np
import os
import re

from agrisatpy.io import SatDataHandler
from pathlib import Path

from logger import get_logger

logger = get_logger('_analyze_unc_contrib')

# define Sentinel-2 bands to analyze
s2_band_selection = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12']


def calc_uncertainty(
        scenario_dir: Path,
        out_dir: Path
    ):
    """
    Calculates absolute uncertainty per uncertainty type (random or systematic)
    and spectral band as the standard deviation of all scenario members.

    :param scenario_dir:
        directory where the Monte-Carlo scenarios are stored (organized by scene)
    :param out_dir:
        directory where to store the results (create sub-directories for the
        scenes analyzed)
    """

    # get scenes
    search_expr = 'S2*_MSIL1C_*'
    scenes = glob.glob(scenario_dir.joinpath(search_expr).as_posix())

    for scene in scenes:

        # get scenarios (numbered from 1 to N)
        m = re.compile('\d')
        scenarios = [f for f in os.listdir(scene) if m.search(f)]
        scenarios = [Path(scene).joinpath(x) for x in scenarios]

        out_dir_scene = out_dir.joinpath(Path(scene).name)
        out_dir_scene.mkdir(exist_ok=True)

        logger.info(f'Working on {scene}')

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
                rand_handler = SatDataHandler()
                rand_handler.read_from_bandstack(
                    fname_bandstack=fname_rand
                )
                rand_data = rand_handler.get_band('B1')

                sys_handler = SatDataHandler()
                sys_handler.read_from_bandstack(
                    fname_bandstack=fname_sys
                )
                sys_data = sys_handler.get_band('B1')

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
                band_data=sys_unc,
                band_name='sys_unc',
                snap_band='B1'
            )
            rand_handler.add_band(
                band_data=rand_unc,
                band_name='rand_unc',
                snap_band='B1'
            )

            # save as raster
            sys_handler.write_bands(
                out_file=fname_sys,
                band_names=['sys_unc']
            )
            rand_handler.write_bands(
                out_file=fname_rand,
                band_names=['rand_unc']
            )

        logger.info(f'Finished {scene}')


def map_uncertainty(unc_results_dir: Path):
    """
    Generates plots of system and random uncertainty components per spectral
    band
    """

    pass


if __name__ == '__main__':

    scenario_dir = Path('../S2A_MSIL1C_orig/uncertainty_contributors')
    out_dir = Path('../S2A_MSIL1C_orig/uncertainty_contributors/L1C_Analysis')
    out_dir.mkdir(exist_ok=True)

    calc_uncertainty(scenario_dir, out_dir)
    
    