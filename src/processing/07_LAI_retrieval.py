'''

'''

import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path

from agrisatpy.core.band import Band
from agrisatpy.core.raster import RasterCollection
from agrisatpy.core.sensors import Sentinel2
from rtm_inv.inversion import inv_img, retrieve_traits
from logger import get_logger

logger = get_logger('07_LAI_retrieval')

# settings for inversion
n_solutions = 100
cost_function = 'rmse'
traits = ['lai']

def loop_scenarios(
        scenario_dir: Path,
        field_parcels: Path,
        lut_dir: Path
    ):
    """
    Loops over the scenarios of all S2 scenes found and estimates the green
    leaf area index (GLAI) by inverting ProSAIL using lookup-tables. The
    inversion is carried out on field parcel pixels, only.

    :param scenario_dir:
        directory where the S2 scenes and their scenario runs are stored
    :param field_parcels:
        ESRI shapefile denoting the field parcel pixels
    :param lut_dir:
        directory with the ProSAIL lookup tables for each scene (*.pkl)
    """

    # find scenes for which scenarios are available
    scenes = glob.glob(scenario_dir.joinpath('S2*_MSIL1C*').as_posix())

    # loop over scenes and their scenarios
    for idx, scene in enumerate(scenes):

        logger.info(f'Working on scene {scene} ({idx+1}/{len(scenes)})')

        # find L2A scenes
        scenarios = glob.glob(Path(scene).joinpath('*/S2*_MSIL2A*.SAFE').as_posix())

        # TODO: find ProSAIL LUT for the current scene
        scene = Path(scene)
        platform = scene.name.split('_')[0]
        timestamp = scene.name.split('_')[2]
        
        scene_lut = next(lut_dir.glob(f'{platform}_MSIL2A_{timestamp}*-prosail_lut.pkl'))
        s2_lut = pd.read_pickle(scene_lut)

        # loop over scenarios of the current scene
        for jdx, scenario in enumerate(scenarios):

            logger.info(f'Starting LAI retrieval {jdx+1}/{len(scenarios)} ({scenario})')
            # define input and outputs for the LAI model
            scenario = Path(scenario)
            vi_dir = scenario.parent.joinpath('Vegetation_Indices')
            out_dir = vi_dir

            s2_ds = Sentinel2().from_safe(
                in_dir=scenario,
                read_scl=False,
                vector_features=field_parcels
            )
            # resample to 10m spatial resolution using nearest neigbor interpolation
            s2_ds.resample(target_resolution=10, inplace=True)
            # scale reflectance values between 0 and 1
            s2_ds.scale(inplace=True)
            s2_array = s2_ds.get_values()
            mask = s2_array.mask
            s2_spectra = s2_array.data
            # get simulated spectra (ProSAIL forward runs)
            s2_lut_spectra = s2_lut[s2_ds.band_names].values

            lut_idxs = inv_img(
                lut=s2_lut_spectra,
                img=s2_spectra,
                mask=mask[0,:,:],
                cost_function=cost_function,
                n_solutions=n_solutions,
            )
            trait_img = retrieve_traits(
                lut=s2_lut,
                lut_idxs=lut_idxs,
                traits=traits
            )
            # save LAI to gTiff
            collection = RasterCollection(
                Band,
                geo_info=s2_ds['B12'].geo_info,
                band_name='GLAI',
                values=trait_img[0,:,:]
            )

            date = scenario.name.split('_')[2][0:8]
            tile = scenario.name.split('_')[5]
            sensor = scenario.name.split('_')[0]
            fname_lai = f'VI_{date}_{tile}_MSIL2A_{sensor}_None_10m_GLAI.tif'
            fpath_out = out_dir.joinpath(fname_lai)
            collection.to_rasterio(fpath_raster=fpath_out)

            logger.info(f'Finished LAI retrieval {jdx+1}/{len(scenarios)} ({scenario})')
        logger.info(f'Finished scene {scene} ({idx+1}/{len(scenes)})')
    logger.info('Done')

if __name__ == '__main__':

    # input directories and files
    field_parcels = Path('../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    lut_dir = Path('../S2_ProSAIL_LUTs')

    for batch in range(1,6):
        scenario_dir = Path(f'../S2_MSIL1C_RUT-Scenarios/batch_{batch}')
        loop_scenarios(scenario_dir, field_parcels, lut_dir)
        
    # scenario_dir = Path('/mnt/ides/Lukas/04_Work/GPR_LAI/S2_RUT_Scenarios')
