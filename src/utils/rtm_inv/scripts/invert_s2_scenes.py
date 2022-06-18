'''
Lookup-table based inversion of Sentinel-2 scenes (original S2 scenes in L2A processing
level, no uncertainty applied).

The inversion strategy is based on the median value of the 100 best performing solutions
in terms of the minimum root mean squared error (RMSE) between observed and simulated
spectra.
'''

import numpy as np
import pandas as pd

from pathlib import Path
from agrisatpy.config import get_settings
from agrisatpy.core.band import Band
from agrisatpy.core.raster import RasterCollection
from agrisatpy.core.sensors import Sentinel2

from rtm_inv.inversion import inv_img, retrieve_traits


logger = get_settings().logger

if __name__ == '__main__':

    traits = ['lai']
    msil2a_scenes_dir = Path(
        '/home/graflu/Documents/uncertainty/S2_MSIL1C_orig'
    )
    lut_dir = Path(
        '/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_ProSAIL_LUTs'
    )
    aoi = Path(
        '/mnt/ides/Lukas/software/scripts_paper_uncertainty/shp/AOI_Esch_EPSG32632.shp'
    )

    n_solutions = 100
    cost_function = 'rmse'

    # loop over scenes and invert them using the median of the 100 best fitting solutions
    for scene in msil2a_scenes_dir.rglob('S2*_MSIL2A_*.SAFE'):

        # search for the corresponding lookup-table
        scene_lut = lut_dir.joinpath(scene.name.replace('.SAFE', '-prosail_lut.pkl'))
        s2_lut = pd.read_pickle(scene_lut)

        s2_collection = Sentinel2().from_safe(
            in_dir=scene,
            read_scl=False,
            vector_features=aoi
        )

        # resample to 10m spatial resolution using nearest neighbor interpolation
        s2_collection.resample(target_resolution=10, inplace=True)
        # scale reflectance values between 0 and 1
        s2_collection.scale(inplace=True)

        # invert the S2 scene by comparing ProSAIL simulated to S2 observed spectra
        s2_lut_spectra = s2_lut[s2_collection.band_names].values
        s2_spectra = s2_collection.get_values()
        if isinstance(s2_spectra, np.ma.MaskedArray):
            mask = s2_spectra.mask[0,:,:]
            s2_spectra = s2_spectra.data
        else:
            mask = np.zeros(shape=(s2_spectra.shape[1], s2_spectra.shape[2]), dtype='uint8')
            mask = mask.as_type('bool')

        logger.info(f'Starting inversion of {scene.name}')
        lut_idxs = inv_img(
            lut=s2_lut_spectra,
            img=s2_spectra,
            mask=mask,
            cost_function=cost_function,
            n_solutions=n_solutions,
        )
        trait_img = retrieve_traits(
            lut=s2_lut,
            lut_idxs=lut_idxs,
            traits=traits
        )
        logger.info(f'Finished inversion of {scene.name}')
        # save LAI to gTiff
        collection = RasterCollection(
            Band,
            geo_info=s2_collection['B12'].geo_info,
            band_name='GLAI',
            values=trait_img[0,:,:]
        )
        out_dir_vis = scene.name.replace('.SAFE', '.VIs')
        date = scene.name.split('_')[2][0:8]
        tile = scene.name.split('_')[5]
        sensor = scene.name.split('_')[0]
        fname_lai = f'VI_{date}_{tile}_MSIL2A_{sensor}_None_10m_GLAI.tif'
        fpath_out = msil2a_scenes_dir.joinpath(out_dir_vis).joinpath('Vegetation_Indices').joinpath(fname_lai)
        collection.to_rasterio(fpath_raster=fpath_out)
        logger.info(f'Wrote LAI product to file: {fpath_out}')
    