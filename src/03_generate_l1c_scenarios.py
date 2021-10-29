'''
Created on Oct 20, 2021

@author:     Lukas Graf (D-USYS, ETHZ)

@purpose:    This script takes the radiometric uncertainty and
             creates L1C scenarios using the original L1C scene.
             To speed things up, it only samples over the study
             area and sets all other parts of the scene to NoData.
             This helps to make the sampling process faster and
             also enables a faster execution of Sen2Cor in the
             next processing step.
'''

import shutil
import glob
from pathlib import Path
from typing import List
import rasterio as rio
import numpy as np
import os
from pip._vendor.pkg_resources import PathMetadata
from patsy import origin


def _sample_fast(
        data_val: float,
        unc_val: float
    ) -> float:
    """
    helper function to sample from random distribution
    for each pixel. It will be called vectorized later
    and save a lot of time...

    :param data_val:
        original data (i.e., reflectance) value
    :param unc_val:
        associated uncertainty (standard uncertainty)
    :return:
        original data + random value sampled from normalized
        Gaussian distribution with standard deviation of the
        uncertainty and mean zero
    """
    rval = np.random.normal(
        loc=0,
        scale=data_val*unc_val,
        size=1
    )
    return data_val + rval


def gen_rad_unc_scenarios(
        orig_dataset_path: Path,
        unc_dataset_path: Path,
        scenario_path: Path,
        template_path: Path,
        n_scenarios: int,
        roi_bounds_10m: List[float]
    ) -> None:
    """
    Taking the original Sentinel-2 L1C scene and the radiometric uncertainty
    derived from running the Sentinel-2 radiometric uncertainty toolbox (S2RUT)
    this function creates n_scenario possible "alternative" scene realities that
    can be used for further processing, using, e.g., Sen2cor.

    A region of interest (ROI) can speed things up a lot and is therefore strongly
    recommended to use instead of processing the entire image.

    :param orig_dataset_path:
        original S2 scene in .SAFE format
    :param unc_dataset_path:
        outputs of S2RUT per spectral band generated using the provided shell script
    :param scenario_path:
        directory where to store the "alternative" scene realities (scenarios)
    :param template_path:
        template of the corresponding .SAFE folder structure without the
        spectral bands (will be created using this function)
    :param n_scenarios:
        number of scenarios to generate
    :param roi_bounds_10m:
        region of interest (ROI) boundaries in image coordindates in 10m resolution.
        Expected: [col_min, col_max, row_min, row_max]
    """

    # create scenario output folders in .SAFE structure
    for idx in range(n_scenarios):
        current_scenario_path = scenario_path.joinpath(str(idx+1))
        # copy template
        shutil.copytree(
            template_path,
            current_scenario_path
        )

    # find spectral bands and uncertainty images
    s2_bands = glob.glob(
        orig_dataset_path.joinpath('GRANULE/*/IMG_DATA/*_B*.jp2').as_posix()
    )
    unc_bands = glob.glob(
        unc_dataset_path.joinpath('*_rut_b*.tif').as_posix()
    )
    unc_band_ids = [x.split('_')[-1].split('.')[0] for x in unc_bands]

    sample_fast_v = np.vectorize(_sample_fast)

    # roi bounds in 20m image coordinates
    roi_bounds_20m = [int(x/2) for x in roi_bounds_10m]
    roi_bounds_60m = [int(x/6) for x in roi_bounds_10m]

    # iterate over each spectral band and generate the scenarios
    bands_10m = ['b02', 'b03', 'b04', 'b08']
    bands_20m = ['b05', 'b06', 'b07', 'b8a', 'b11', 'b12']
    bands_60m = ['b01', 'b09', 'b10']

    for s2_band in s2_bands:

        print(f'Running scenarios for {s2_band}')

        spectral_band_id = Path(s2_band).name.split('_')[2].split('.')[0].lower()

        # check for correct roi bands
        if spectral_band_id in bands_10m:
            roi_bounds = roi_bounds_10m
        elif spectral_band_id in bands_20m:
            roi_bounds = roi_bounds_20m
        else:
            roi_bounds = roi_bounds_60m

        unc_band = unc_bands[unc_band_ids.index(spectral_band_id)]

        # read uncertainty
        with rio.open(unc_band, 'r') as src:
            unc_data = src.read()
        
        # read original spectral values
        with rio.open(s2_band) as src:
            meta = src.meta
            band_data = src.read()
        
        # using the appropriate compression avoids artifacts between original and
        # rasterio derived jp2 files
        meta.update({
            'QUALITY': '100',
            'REVERSIBLE': 'YES',
        })

        # sample from normalized random distribution
        # undo scaling of values
        unc_data = unc_data * 0.001

        # define output directory
        band_fname = Path(s2_band).name
        dataset_path = str(Path(s2_band).parent).split(os.sep)
        dot_safe = [x for x in dataset_path if '.SAFE' in x][0]
        dataset = os.path.sep.join(dataset_path[dataset_path.index(dot_safe)::])
            
        for idx in range(n_scenarios):
            
            unc_realiz = np.zeros(shape=(1, unc_data.shape[1], unc_data.shape[2]))

            # for loop over region of interest (ROI)
            # rows
            for i in range(roi_bounds[2], roi_bounds[3]):
                # cols
                for j in range(roi_bounds[0], roi_bounds[1]):
                    unc_realiz[0,i,j] = np.random.normal(
                        loc=0,
                        scale=unc_data[0,i,j]*band_data[0,i,j],
                        size=1
                    )
                    unc_realiz[0,i,j] = band_data[0,i,j] + unc_realiz[0,i,j]

            # save to scenarios
            current_scenario_path = scenario_path.joinpath(str(idx+1))
            file_dst = current_scenario_path.joinpath(dataset).joinpath(band_fname).as_posix()
            with rio.open(file_dst, 'w', **meta) as dst:
                dst.write(unc_realiz[0,:,:].astype(np.uint16), 1) # must be unsigned int16
                print(f'Wrote scenario {idx+1}/{n_scenarios} for {s2_band}')


def main(
        orig_datasets_dir: Path,
        unc_datasets_dir: Path,
        scenario_dir: Path,
        roi_bounds_10m: List[int],
        n_scenarios: int
    ) -> None:
    """
    main executable function taking care about everything
    """

    # find scenes and their uncertainty
    orig_datasets = glob.glob(orig_datasets_dir.joinpath('*.SAFE').as_posix())
    unc_datasets = glob.glob(unc_datasets_dir.joinpath('*.RUT').as_posix())

    # loop over the scenes. Before generating the scenarios some
    # preparation is required
    for orig_dataset in orig_datasets:

        orig_dataset_path = Path(orig_dataset)
        scene_name = orig_dataset_path.name

        print(f'** Working on {scene_name}')

        # find corresponding uncertainty directory
        unc_dataset_path = [
            Path(x) for x in unc_datasets if Path(x).name.split('.')[0] == scene_name.split('.')[0]
        ][0]

        # preparation
        # create a subdirectory in the scenario folder with the scene name without .SAFE
        scene_name_without_safe = scene_name.replace('.SAFE', '')

        scenario_path = scenario_dir.joinpath(scene_name_without_safe)
        if not scenario_path.exists():
            scenario_path.mkdir()
        # create a folder to store the scene template (entire .SAFE structure)
        template_path = scenario_path.joinpath('template')
        if not template_path.exists():
            template_path.mkdir()

        # copy the original L1C scene into the template folder and delete
        # the jp2 files in the GRANULE directory
        shutil.copytree(orig_dataset_path, template_path.joinpath(scene_name))

        # delete the jp2 files in the template
        search_expr = '*.SAFE/GRANULE/*/IMG_DATA/*_B*.jp2'
        jp2_files = glob.glob(template_path.joinpath(search_expr).as_posix())
        for jp2_file in jp2_files:
            os.remove(jp2_file)

        # finally, generate the scenarios
        gen_rad_unc_scenarios(
            orig_dataset_path=orig_dataset_path,
            unc_dataset_path=unc_dataset_path,
            scenario_path=scenario_path,
            template_path=template_path,
            n_scenarios=n_scenarios,
            roi_bounds_10m=roi_bounds_10m
        )

        # clean up and remove the template directory as it not required any further
        shutil.rmtree(template_path, ignore_errors=True)


if __name__ == '__main__':

    ### define user inputs

    # directory with L1C data (.SAFE subdirectories)
    orig_datasets_dir = Path('./../S2A_MSIL1C_orig')

    # directory with radiometric uncertainty outputs (.RUT subdirectories)
    unc_datasets_dir = orig_datasets_dir

    # directory where to store the scenarios (a subdirectory will be created for each scene)
    # in which the actual scenarios are placed
    scenario_dir = Path('./../S2A_MSIL1C_RUT-Scenarios')

    # define bounds of the study area (encompassing the single regions of interest)
    # bounds col_min, col_max, row_min, row_max (image coordinates of the 10m raster)
    roi_bounds_10m = [7000,8000,4000,5000]

    # number of scenarios (each scenario is a possible realization of a S2 scene!)
    n_scenarios = 100
    
    main(
        orig_datasets_dir,
        unc_datasets_dir,
        scenario_dir,
        roi_bounds_10m,
        n_scenarios
    )
