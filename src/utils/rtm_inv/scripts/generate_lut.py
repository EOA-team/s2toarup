'''
Generates a lookup table (LUT) with RTM forward runs
'''

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from typing import Optional

from agrisatpy.metadata.sentinel2 import parse_s2_scene_metadata
from agrisatpy.utils.sentinel2 import get_S2_platform_from_safe

from rtm_inv.lookup_table import LookupTable
from rtm_inv.rtm_adapter import RTM

plt.style.use('ggplot')

def generate_lut(
        sensor: str,
        lut_params: pd.DataFrame,
        lut_size: Optional[int] = 50000,
        rtm_name: Optional[str] = 'prosail',
        sampling_method: Optional[str] = 'LHS',
        solar_zenith_angle: Optional[float] = None,
        viewing_zenith_angle: Optional[float] = None,
        solar_azimuth_angle: Optional[float] = None,
        viewing_azimuth_angle: Optional[float] = None
    ) -> pd.DataFrame:

    # overwrite angles in LUT DataFrame
    lut_params.loc[lut_params['Parameter'] == 'tts','Min'] = solar_zenith_angle
    lut_params.loc[lut_params['Parameter'] == 'tts','Max'] = solar_zenith_angle
    lut_params.loc[lut_params['Parameter'] == 'tto', 'Min'] = viewing_zenith_angle
    lut_params.loc[lut_params['Parameter'] == 'tto', 'Max'] = viewing_zenith_angle
    # calculate relative azimuth (psi)
    psi = abs(solar_azimuth_angle - viewing_azimuth_angle)
    lut_params.loc[lut_params['Parameter'] == 'psi', 'Min'] = psi
    lut_params.loc[lut_params['Parameter'] == 'psi', 'Max'] = psi

    # get input parameter samples first
    lut = LookupTable(params=lut_params)
    lut.generate_samples(lut_size, sampling_method)

    # and run the RTM in forward mode in the second step
    # outputs get resampled to the spectral resolution of the sensor
    rtm = RTM(lut=lut, rtm=rtm_name)
    lut_simulations = rtm.simulate_spectra(sensor=sensor)
    return lut_simulations

if __name__ == '__main__':

    # define lookup-table parameter ranges and distributions
    path_lut_params = Path('../../parameters/prosail_s2.csv')
    lut_params = pd.read_csv(path_lut_params)

    out_dir = Path('/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_ProSAIL_LUTs')

    # define locations of Sentinel-2 scenes (.SAFE)
    dot_safe_base_dir = Path('/home/graflu/Documents/uncertainty/S2_MSIL1C_orig')

    # loop over Sentinel-2 scenes in L2A processing level
    for in_dir_safe in dot_safe_base_dir.glob('S2*_MSIL2A_*.SAFE'):
        
        # get viewing and illumination angles
        scene_metadata, _ = parse_s2_scene_metadata(in_dir=in_dir_safe)
        solar_zenith_angle = scene_metadata['SUN_ZENITH_ANGLE']
        solar_azimuth_angle = scene_metadata['SUN_AZIMUTH_ANGLE']
        viewing_zenith_angle = scene_metadata['SENSOR_ZENITH_ANGLE']
        viewing_azimuth_angle = scene_metadata['SENSOR_AZIMUTH_ANGLE']
        # get platform
        platform = get_S2_platform_from_safe(dot_safe_name=in_dir_safe)
        # map to full plafform name
        full_names = {'S2A': 'Sentinel2A', 'S2B': 'Sentinel2B'}
        platform = full_names[platform]
    
        # call function to generate lookup-table
        lut = generate_lut(
            sensor=platform,
            lut_params=lut_params,
            solar_zenith_angle=solar_zenith_angle,
            viewing_zenith_angle=viewing_zenith_angle,
            solar_azimuth_angle=solar_azimuth_angle,
            viewing_azimuth_angle=viewing_azimuth_angle
        )

        # box-plot of simulated spectra
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot()
        spectral_boxplot = lut.boxplot(
            column=['B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12'],
            ax=ax
        )
        ax.set_ylabel('Surface Reflectance Factors [-]')
        ax.set_ylim(0,1)
        ax.set_title(f'ProSAIL Simulated {platform} Spectra (N={lut.shape[0]})')
        fname_boxplot = out_dir.joinpath(in_dir_safe.name.replace('.SAFE', '-prosail_spectra.png'))
        fig.savefig(fname_boxplot, dpi=200, bbox_inches='tight')
        plt.close(fig)
        # save lookup-table as pickled object on disk
        fname_out = out_dir.joinpath(in_dir_safe.name.replace('.SAFE', '-prosail_lut.pkl'))
        lut.to_pickle(fname_out)

        print(f'Generated LUT for {in_dir_safe.name}')
    