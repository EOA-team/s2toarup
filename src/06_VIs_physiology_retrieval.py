# -*- coding: utf-8 -*-

'''
This script is used to calculate widely-used vegetation indices (VIs)
for each L2A scenario outcome including the EVI (Enhanced Vegetation Index)
and NDVI (Normalized Difference Vegetation Index).

In addition, it resamples the original Sentinel-2 data in L2A processing level
to 10m spatial resolution using nearest neighbor interpolation and saves the
data clipped to the extent of the study area. For these original datasets NDVI and
EVI are calculated as well; these serve as reference to the scenario runs.

TODO: test this script first on the workstation (autumn and spring)
TODO: then delete the *.VIs folders on the group share and rerun this script also thre
TODO: upload the original datasets from the workstation to the group share
'''

import cv2
import glob

from pathlib import Path
from agrisatpy.operational.resampling.sentinel2.resample_and_stack import _get_output_file_names
from agrisatpy.io.sentinel2 import Sentinel2Handler

from logger import get_logger

# setup logger -> will write log file to the /../log directory
logger = get_logger('l3_vegetation-indices')


def main(
        scenario_dir: Path,
        shapefile_study_area: Path
    ) -> None:
    """
    Main executable method for calculating the NDVI and EVI from Sentinel-2
    data in .SAFE archive structure.

    :param scenario_dir:
        directory where the outcomes of the atmospheric correction are stored as .SAFE
        datasets in L2A processing level
    :param shapefile_study_area:
        shapefile defining the extent of our study area for clipping the data to
    """

    # find scenes for which scenarios are available
    scenes = glob.glob(scenario_dir.joinpath('S2*_MSIL1C*').as_posix())

    # loop over scenes and their scenarios
    for idx, scene in enumerate(scenes):

        logger.info(f'Working on scene {scene} ({idx+1}/{len(scenes)})')

        # find L2A scenes
        scenarios = glob.glob(Path(scene).joinpath('*/S2*_MSIL2A*.SAFE').as_posix())

        # loop over scenarios of the current scene
        for jdx, scenario in enumerate(scenarios):

            logger.info(f'Working on scenario {jdx+1}/{len(scenarios)} ({scenario})')
            # place results in the root of the scenario
            out_dir = Path(scenario).parent

            # read data from SAFE; we only need the bands for the EVI and NDVI
            try:
                handler = Sentinel2Handler()
                handler.read_from_safe(
                    in_dir=Path(scenario), 
                    band_selection=['B02','B04','B08'],
                    polygon_features=shapefile_study_area,
                    read_scl=False
                )
                # calculate the spectral indices
                # define output directory
                vis_dir = out_dir.joinpath('Vegetation_Indices')
                if not vis_dir.exists():
                    vis_dir.mkdir()
    
                # the actual index calculation starts here
                vi_names = ['NDVI', 'EVI']
                fnames = _get_output_file_names(
                    in_dir=Path(scenario),
                    resampling_method='None',
                    target_resolution=10
                )
                in_file = fnames['bandstack']
                for vi_name in vi_names:
                    handler.calc_si(vi_name)
                    vi_fname = vis_dir.joinpath(f'VI_{in_file.split(".")[0]}_{vi_name.upper()}.tif').as_posix()
                    # save to raster
                    handler.write_bands(
                        out_file=vi_fname,
                        band_names=[vi_name]
                    )
                handler = None
            except Exception as e:
                logger.error(f'Scenario {jdx+1} failed ({scenario}): {e}')
                continue
            logger.info(f'Finished scenario {jdx+1}/{len(scenarios)} ({scenario})')

        logger.info(f'Finshed scene {scene} ({idx+1}/{len(scenes)})')


def resample_and_stack_orig_data(
        orig_datasets_dir: Path,
        shapefile_study_area: Path
    ):
    """
    Loops over the original (i.e., reference) Sentinel-2 datasets in L2A processing
    level and resamples the data (spectral bands plus the scene classification layer)
    to 10m spatial resolution. Calculates EVI and NDVI as it is done for single
    scenarios

    :param orig_datasets_dir:
        directory where the S2-datasets in .SAFE format can be found (MSIL2A)
    :param shapefile_study_area:
        shapefile defining the extent of our study area for clipping the data to
    """

    # find Sentinel-2 datasets
    orig_datasets_l2a = glob.glob(Path(orig_datasets_dir).joinpath(
        'S2*_MSIL2A*.SAFE').as_posix())
    
    for idx, orig_dataset in enumerate(orig_datasets_l2a):

        logger.info(f'Processing original dataset {orig_dataset} ({idx+1}/{len(orig_datasets_l2a)})')

        # save results in a folder called as the dataset ending with .VIs instead of
        # .SAFE
        out_dir = Path(orig_dataset).parent
        out_dir = out_dir.joinpath(Path(orig_dataset).name.replace('.SAFE', '.VIs'))
        out_dir.mkdir(exist_ok=True)

        # read dataset from .SAFE for the extent of the study area
        handler = Sentinel2Handler()
        handler.read_from_safe(
            in_dir=Path(orig_dataset), 
            polygon_features=shapefile_study_area
        )

        # resample the bands to 10m using nearest neighbor interpolation
        handler.resample(
            target_resolution=10.,
            resampling_method=cv2.INTER_NEAREST_EXACT
        )

        fnames = _get_output_file_names(
            in_dir=Path(orig_dataset),
            resampling_method='None', # named None to be consistent to what we find in the scenarios
            target_resolution=10
        )

        # write SCL, RGB previews
        fig_rgb = handler.plot_rgb()
        rgb_dir = out_dir.joinpath('rgb_previews')
        rgb_dir.mkdir(exist_ok=True)
        fname_rgb = rgb_dir.joinpath(fnames['rgb_preview'])
        fig_rgb.savefig(fname=fname_rgb, dpi=300, bbox_inches='tight')

        fig_scl = handler.plot_scl()
        scl_dir = out_dir.joinpath('scene_classification')
        scl_dir.mkdir(exist_ok=True)
        fname_scl_preview = scl_dir.joinpath(fnames['scl_preview'])
        fig_scl.savefig(fname=fname_scl_preview, dpi=300, bbox_inches='tight')

        # save rasters (spectral bands, VIs, SCL)
        fname_bandstack = out_dir.joinpath(fnames['bandstack'])
        band_selection = handler.bandnames
        band_selection.remove('scl')
        handler.write_bands(
            out_file=fname_bandstack,
            band_names=band_selection,
            use_band_aliases=True
        )

        fname_scl = scl_dir.joinpath(fnames['scl'])
        handler.write_bands(
            out_file=fname_scl,
            band_names=['scl']
        )

        vi_names = ['NDVI', 'EVI']
        in_file = fnames['bandstack']
        vis_dir = out_dir.joinpath('Vegetation_Indices')
        vis_dir.mkdir(exist_ok=True)

        for vi_name in vi_names:
            handler.calc_si(vi_name)
            vi_fname = vis_dir.joinpath(f'VI_{in_file.split(".")[0]}_{vi_name.upper()}.tif').as_posix()
            # save to raster
            handler.write_bands(
                out_file=vi_fname,
                band_names=[vi_name]
            )

        logger.info(f'Processed original dataset {orig_dataset} ({idx+1}/{len(orig_datasets_l2a)})')


if __name__ == '__main__':

    scenario_dir = Path('../S2A_MSIL1C_RUT-Scenarios')
    
    shapefile_study_area = Path('../shp/AOI_Esch_EPSG32632.shp')
    
    # vegetation indices on scenarios
    main(
        scenario_dir=scenario_dir,
        shapefile_study_area=shapefile_study_area
    )

    # # original Sentinel-2 L2A data
    # orig_l2a_data = Path('../S2A_MSIL1C_orig')
    #
    # # resample bandstacks and SCL for original datasets
    # resample_and_stack_orig_data(
    #     orig_datasets_dir=orig_l2a_data,
    #     shapefile_study_area=shapefile_study_area
    # )
