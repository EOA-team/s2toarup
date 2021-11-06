
import glob
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional
import rasterio as rio


# define L1C uncertainty contributors available from L1C-RUT
l1c_unc_contributors = [
    'noise',
    'stray_sys',
    'stray_rand',
    # 'x_talk', not considered as suggested by Gorrono et al., 2018
    'ADC',
    'DS',
    'gamma',
    'diff_abs',
    'diff_temp',
    'diff_cos',
    'diff_k',
    'quant'
]

# define S2 bands to loop over
s2_bands = [
    'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12'
]
# spatial resolution (meters) of each spectral band
s2_band_res = {
    'B01': 60, 'B02': 10, 'B03': 10, 'B04': 10, 'B05': 20, 'B06': 20, 'B07': 20, 'B08': 10,
    'B8A': 20, 'B09': 60, 'B10': 60, 'B11': 20, 'B12': 20
}


class Band_Data(object):
    """
    Data structure for storing the spectral (i.e., TOA reflectance factors)
    and associated uncertainties contributors per spectral band
    """

    def __init__(self):
        self._r_toa = np.empty(0, dtype=np.uint16)
        self._unc_contrib = {}
    
    @property
    def r_toa(self):
        """
        Array with L1C TOA reflectance values
        """
        return self._r_toa

    @r_toa.setter
    def r_toa(self, arr: np.array):
        self._r_toa = arr

    @property
    def unc_contrib(self):
        """
        Dictionary with uncertainty contributors. The dict
        keys are the single uncertainty contributors and the
        values are the uncertainties obtained from L1C-RUT.
        """
        return self._unc_contrib

    @unc_contrib.setter
    def unc_contrib(self, unc_dict: Dict[str, np.array]):
        self._unc_contrib = unc_dict



def gen_rad_unc_scenarios(
        orig_dataset_path: Path,
        unc_dataset_path: Path,
        scenario_path: Path,
        template_path: Path,
        n_scenarios: int,
        roi_bounds_10m: List[float],
        correlation_table: Optional[Path]='uncertainty_contributors_correlation.csv'
    ) -> None:
    """
    Taking the original Sentinel-2 L1C scene and the radiometric uncertainty
    derived from running the Sentinel-2 radiometric uncertainty toolbox (S2RUT)
    this function creates time n_scenarios possible "alternative" scene realities
    that can be used for further processing, using, e.g., Sen2cor.

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
    :param correlation_table:
        file-path to the CSV with the correlation table obtained from Gorrono et al.
        2018, Table 1. This table describes the correlation of the single uncertainty
        contributors in the three domains (spatial, temporal, spectral).
    """

    # create scenario output folders in .SAFE structure
    for idx in range(n_scenarios):
        current_scenario_path = scenario_path.joinpath(str(idx+1))
        # copy template
        shutil.copytree(
            template_path,
            current_scenario_path,
            dirs_exist_ok=True
        )

    # roi bounds in all spatial resolutions for sub-setting the data
    roi_bounds_20m = [int(x/2) for x in roi_bounds_10m]
    roi_bounds_60m = [int(x/6) for x in roi_bounds_10m]
    roi_bounds_all = {
        10: roi_bounds_10m,
        20: roi_bounds_20m,
        60: roi_bounds_60m
    }
    roi_size_all = dict.fromkeys(roi_bounds_all.keys())
    for res in roi_size_all:
        roi_size_all[res] = {
            'n_col': roi_bounds_all[res][1]-roi_bounds_all[res][0],
            'n_row': roi_bounds_all[res][3]-roi_bounds_all[res][2]
        }

    # read the band data and the associated uncertainties
    full_img_size = dict.fromkeys([10, 20, 60])
    mc_input_data = dict.fromkeys(s2_bands)
    for s2_band in s2_bands:

        # get band name alias (without zero, i.e, B01 -> B1)
        s2_band_alias = s2_band.replace('0','')

        # get spatial resolution and corresponding ROI bounds
        spatial_res = s2_band_res[s2_band]
        roi_bounds = roi_bounds_all[spatial_res]
        min_row, max_row = roi_bounds[2], roi_bounds[3]
        min_col, max_col = roi_bounds[0], roi_bounds[1]

        # instanciate band data struct
        band_data = Band_Data()

        # read TOA reflectance factor data
        r_toa_file = glob.glob(
            orig_dataset_path.joinpath(f'GRANULE/*/IMG_DATA/*_{s2_band}.jp2').as_posix()
        )[0]

        with rio.open(r_toa_file, 'r') as src:
            n_rows_full = src.height
            n_cols_full = src.width
            r_toa = src.read(1)

        # remember the original image size
        if full_img_size[spatial_res] is None:
            full_img_size[spatial_res] = (n_rows_full, n_cols_full)

        band_data.r_toa = r_toa[min_row:max_row,min_col:max_col]

        # read the single uncertainty contributors
        unc_contrib_dict = dict.fromkeys(l1c_unc_contributors)
        for l1c_unc_contributor in l1c_unc_contributors:
            
            unc_contrib_file = glob.glob(
                unc_dataset_path.joinpath(f'S2*_rut_{l1c_unc_contributor}_{s2_band_alias}.tif').as_posix()
            )[0]

            with rio.open(unc_contrib_file, 'r') as src:
                unc = src.read(1)
            unc_contrib_dict[l1c_unc_contributor] = unc[min_row:max_row,min_col:max_col]
        band_data.unc_contrib = unc_contrib_dict

        # save the extracted reflectance and uncertainty of the current band
        mc_input_data[s2_band] = band_data

    # sample the L1C scenarios using Monte-Carlo simulations. The simulations implement the
    # correlation in the spectral, spatial and temporal domain.
    # The domains spectral, spatial and temporal are defined as
    #
    #    spectral - error correlation between the single spectral bands
    #    spatial  - error correlation along track (North-South direction)
    #    temporal - correlation across track (East-West directon)
    #
    # If a contributor is entirely uncorrelated in all three domains, the contributor
    # is sampled independently for each spectral band and pixel.
    # If a contributor is fully correlated in one or all dimensions, then one distribution
    # is used to sample from all correlated dimensions.
    # If a contributor is partly correlated then two distributions are used; one reflects
    # the independent case and one the fully correlated part. The contributions are then
    # combined based on the strength of the correlation.
    # This approach also allows to account for correlations in a single or two dimensions,
    # only.
    #
    # The knowledge about the uncertainty contributor correlation is taken from
    # Gorrono et al., 2018 (https://doi.org/10.1080/22797254.2018.1471739), Table 1.
    #

    # load the correlation table
    corr_df = pd.read_csv(correlation_table, index_col='uncertainty_contributor')
    dimensions = ['spectral', 'spatial', 'temporal']

    # assess how the single error contributors are correlated
    # fully uncorrelated contributors
    fully_uncorr_contributors = corr_df[dimensions][
        (corr_df.spectral == 0) & (corr_df.temporal == 0) & (corr_df.spatial == 0)
    ].index

    # constant error terms
    const_error_terms = corr_df[dimensions][
        (corr_df.spectral == 999) & (corr_df.temporal == 999) & (corr_df.spatial == 999)
    ].index.tolist()

    # fully correlated contributors (correlated in all dimensions)
    fully_corr_contributors = corr_df[dimensions][
        (corr_df.spectral == 1) & (corr_df.temporal == 1) & (corr_df.spatial == 1)
    ].index.tolist()
    
    # partly correlated contributors -> all others

    # empty image matrices for writing the samples to (dtype: uint16)
    img_matrices = dict.fromkeys(full_img_size.keys())
    for res in img_matrices:
        img_matrices[res] = np.zeros(shape=full_img_size[res], dtype=np.uint16)

    # start the iteration process
    for scenario in range(n_scenarios):
        
        print(f'Creating scenario {scenario+1}/{n_scenarios}')

        # empty arrays for storing the errors
        error_band_dict = dict.fromkeys(s2_bands)
        for s2_band in s2_bands:
            error_band_dict[s2_band] = np.zeros_like(mc_input_data[s2_band].r_toa.astype(np.float16))

        for s2_band in s2_bands:

            # get spatial resolution and corresponding array size of the ROI
            cols_and_rows = roi_size_all[s2_band_res[s2_band]]
            num_row = cols_and_rows['n_row']
            num_col = cols_and_rows['n_col']

            # loop over uncorrelated contributors
            for fully_uncorr_contributor in fully_uncorr_contributors:
                uncorr_rut = mc_input_data[s2_band].unc_contrib[fully_uncorr_contributor]

                # check type of distribution
                dist_type = corr_df[corr_df.index == fully_uncorr_contributor]['distribution'].values[0]
                if dist_type == 'normal':
                    uncorr_sample = np.random.normal(
                        loc=0,
                        scale=uncorr_rut,
                        size=(num_row, num_col)
                    ) / 10 # 10 because of scaling of S2-RUT images
                elif dist_type == 'uniform':
                    uncorr_sample = np.random.uniform(
                        low=-uncorr_rut * np.sqrt(3),
                        high=uncorr_rut * np.sqrt(3),
                        size=(numrow, numcol)
                    ) / 10 # 10 because of scaling of S2-RUT images
                error_band_dict[s2_band] += uncorr_sample

            # loop over constant error terms; they are simply added
            for const_error_term in const_error_terms:
                error_band_dict[s2_band] += \
                    mc_input_data[s2_band].unc_contrib[const_error_term] / 10 # 10 because of scaling of S2-RUT images

        # fully correlated contributors
        # append these to a list of arrays and concatenate them into a 1d-array
        # this way it is possible to combine spectral bands with different pixel sizes
        for fully_corr_contributor in fully_corr_contributors:

            band_unc_arr_list = []
            for s2_band in s2_bands:
                band_unc_arr_list.append(mc_input_data[s2_band].unc_contrib[fully_corr_contributor])

            # remember the size of the original arrays so that they can be
            # reshaped into 2d afterwards
            trailing_indices = [0]
            trailing_indices.extend([x.shape[0]*x.shape[1] for x in band_unc_arr_list])

            # concatente all 2d arrays into a single 1d array, axis=None flattens the array
            corr_unc = np.concatenate(band_unc_arr_list, axis=None)

            # sample from the same normal or uniform distribution depending on the contributor
            dist_type = corr_df[corr_df.index == fully_corr_contributor]['distribution'].values[0]
            if dist_type == 'normal':
                corr_rut = np.ones(shape=corr_unc.shape) * corr_unc
                corr_rut = np.random.normal(0, 1, 1)[0] * corr_rut # divide by 10 is not required here because of N(0,1)
            elif dist_type == 'uniform':
                corr_rut = np.empty(shape=corr_unc.shape)
                corr_rut = np.random.uniform(-1, 1, 1)[0] * corr_rut * np.sqrt(3) / 10

            # undo the flattening of the band arrays and add the samples to
            # the error_band_dict
            for idx, s2_band in enumerate(s2_bands):
                band_samples = corr_rut[trailing_indices[idx]:trailing_indices[idx]+trailing_indices[idx+1]].reshape(
                    mc_input_data[s2_band].unc_contrib[fully_corr_contributor].shape
                )
                error_band_dict[s2_band] += uncorr_sample


def main(
        orig_datasets_dir: Path,
        unc_datasets_dir: Path,
        scenario_dir: Path,
        roi_bounds_10m: List[int],
        n_scenarios: int
    ) -> None:
    """
    main executable function taking care about generating the scenarios and
    the required .SAFE folder structure

    :param orig_datasets_dir:
        directory where the original L1C TOA .SAFE scenes can be found
    :param unc_datasets_dir:
        directory where the associated L1C radiometric uncertainty results
        can be found per contributor and band
    :param scenario_dir:
        directory where to create and store the resulting L1C TOA scenarios
    :param roi_bounds_10m:
        because of the high computational load, the scenarios are created
        for a spatial subset of the original scene, only
    :param n_scenarios:
        number of scenarios to create (>=100).
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


if __name__ == '__main__':

    # debug
    orig_dataset_path = Path(
        './../S2A_MSIL1C_orig/done/S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651.SAFE'
    )
    unc_dataset_path = Path(
        './../S2A_MSIL1C_orig/done/S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651.RUT'
    )
    scenario_path = Path(
        './../debug'
    )
    template_path = Path(
        './../debug/template'
    )
    n_scenarios = 1
    roi_bounds_10m = [7000,8000,4000,5000]
    
    gen_rad_unc_scenarios(orig_dataset_path, unc_dataset_path, scenario_path, template_path, n_scenarios, roi_bounds_10m)
    

    # ### define user inputs
    #
    # # directory with L1C data (.SAFE subdirectories)
    # orig_datasets_dir = Path('./../S2A_MSIL1C_orig')
    #
    # # directory with radiometric uncertainty outputs (.RUT subdirectories)
    # unc_datasets_dir = orig_datasets_dir
    #
    # # directory where to store the scenarios (a subdirectory will be created for each scene)
    # # in which the actual scenarios are placed
    # scenario_dir = Path('./../S2A_MSIL1C_RUT-Scenarios')
    #
    # # define bounds of the study area (encompassing the single regions of interest)
    # # bounds col_min, col_max, row_min, row_max (image coordinates of the 10m raster)
    # # TODO: get the bounds from the shapefile of the study area!!
    # roi_bounds_10m = [7000,8000,4000,5000]
    #
    # # number of scenarios (each scenario is a possible realization of a S2 scene!)
    # n_scenarios = 100
    #
    # main(
    #     orig_datasets_dir,
    #     unc_datasets_dir,
    #     scenario_dir,
    #     roi_bounds_10m,
    #     n_scenarios
    # )
