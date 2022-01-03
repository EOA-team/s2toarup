# -*- coding: utf-8 -*-

'''
This script is used to calculate widely-used
vegetation indices (VIs) for each L2A scenario outcome.
To do so, the S2 L2A outputs after Sen2Cor are
clipped (masked) to the study area and bandstacked
into a single geoTiff file. The spatial resolution of the
20m bands is therefore increased to 10m without
modifying the spectral data since the spatial resampling
procedure comes along with its own uncertainty.
The VIs are stored as geoTiff files in the uncertainty directories
(thus, per scenario).
'''

import glob
from pathlib import Path
from typing import Optional

from agrisatpy.spatial_resampling.sentinel2 import resample_and_stack_s2
from agrisatpy.io.sentinel2 import Sentinel2Handler


def calc_indices(
        in_file: Path,
        out_dir: Path,
        shapefile_study_area: Path
    ) -> None:
    """
    Calculates the NDVI and EVI for a band-stacked,
    resampled Sentinel-2 geoTiff file

    :param in_file:
        file-path to the resampled, band-stacked geoTiff file
    :param out_dir:
        directory where to write the results to. These have the
        same name as the input plus the name of the index calculated:
        <in_file.name>_<index>.tif
    :param shapefile_study_area:
        shapefile defining the extent of our study area for clipping the data to
    """

    # open the file and read those bands required for calculating the indices
    handler = Sentinel2Handler()
    handler.read_from_bandstack(
        fname_bandstack=in_file,
        in_file_aoi=shapefile_study_area,
        band_selection=['B02','B03','B04','B08']
    )

    # define output directory
    vis_dir = out_dir.joinpath('Vegetation_Indices')
    if not vis_dir.exists():
        vis_dir.mkdir()

    # the actual index calculation starts here
    vi_names = ['NDVI', 'EVI']
    for vi_name in vi_names:
        handler.calc_vi(vi_name)
        vi_fname = vis_dir.joinpath(f'VI_{in_file.name.split(".")[0]}_{vi_name.upper()}.tif').as_posix()
        # save to raster
        handler.write_bands(
            out_file=vi_fname,
            band_names=[vi_name]
        )


def main(
        scenario_dir: Path,
        shapefile_study_area: Path,
        is_orig_data: Optional[bool] = False
    ) -> None:
    """
    Main executable method for resampling the Sentinel-2 datasets to 10m spatial
    resolution and calculating the NDVI and EVI.

    There are two processing options:
        a) calculating VIs for the scenario outcomes after Sen2Cor
        b) calculating VIs for the original S2 datasets after Sen2Cor
    Option a) is the default case. To run b) the variable `is_orig_data`
    must be set to True

    :param scenario_dir:
        directory where the outcomes of the atmospheric correction are stored as .SAFE
        datasets in L2A processing level
    :param shapefile_study_area:
        shapefile defining the extent of our study area for clipping the data to
    :param is_orig_data:
        if False (default) assumes working on the scenario outcomes. If True, assumes
        working on the atmospherically corrected original datasets (the only difference
        is a minor modification in the path search pattern)
    """

    # find scenes for which scenarios are available
    if not is_orig_data:
        scenarios = glob.glob(scenario_dir.joinpath('S2*_MSIL1C*').as_posix())
    else:
        scenarios = ['dummy']

    # loop over scenarios
    for scenario in scenarios:

        # find L2A scenes
        if not is_orig_data:
            orig_datasets_l2a = glob.glob(Path(scenario).joinpath('*/S2*_MSIL2A*.SAFE').as_posix())
        else:
            orig_datasets_l2a = glob.glob(Path(scenario_dir).joinpath('S2*_MSIL2A*.SAFE').as_posix())

        # loop over scenes, resample them for the extent of the study area and
        # calculate the spectral indices
        for orig_dataset in orig_datasets_l2a:

            # place results in the root of the scenario
            out_dir = Path(orig_dataset).parent
            if is_orig_data:
                out_dir = out_dir.joinpath(Path(orig_dataset).name.replace('.SAFE', '.VIs'))
                out_dir.mkdir(exist_ok=True)

            # bandstack and mask the data
            # we "resample" the data using "pixel_division" which only increases
            # the pixel resolution without changing the spectral values
            file_dict = resample_and_stack_s2(
                in_dir=Path(orig_dataset),
                out_dir=out_dir,
                pixel_division=True
            )
            
            # calculate the spectral indices using the resampled data
            # and write them to a separate sub-directory
            calc_indices(
                in_file=file_dict['bandstack'],
                out_dir=out_dir,
                shapefile_study_area=shapefile_study_area
            )


if __name__ == '__main__':

    scenario_dir = Path('../S2A_MSIL1C_RUT-Scenarios')
    
    shapefile_study_area = Path('../shp/AOI_Esch_EPSG32632.shp')
    
    # vegetation indices on scenarios
    main(
        scenario_dir=scenario_dir,
        shapefile_study_area=shapefile_study_area
    )

    # original Sentinel-2 L2A data
    orig_l2a_data = Path('../S2A_MSIL1C_orig')
    main(
        scenario_dir=orig_l2a_data,
        shapefile_study_area=shapefile_study_area,
        is_orig_data=True
    )
