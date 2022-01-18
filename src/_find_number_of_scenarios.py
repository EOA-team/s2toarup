'''
Check how many scenario runs are necessary
'''

import glob
import matplotlib
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from copy import deepcopy
from typing import Optional
from agrisatpy.io import SatDataHandler

from logger import get_logger

plt.style.use('ggplot')
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)

logger = get_logger('_how_many_scenarios')


def plot_uncertainty_number_of_scenarios(
        scenarios_scene_dir: Path,
        sample_polygons: Path,
        vi_name: str,
        out_dir: Path,
        relative_uncertainty: Optional[bool] = True,
        orig_scene_dir: Optional[Path] = None
    ):
    """
    Plots the number of scenarios vs. the corresponding uncertainty
    in one of the vegetation indices/parameters for all pixels
    of a crop type

    :param scenario_scene_dir:
        directory where the outcomes of the S2-RUT scenarios (vegetation
        indices/parameters derived from the L2A data) are located
    :param sample_polygons:
        polygon features with crop type information
    :param vi_name:
        name of the vegetation index/parameter to analyze
    :param relative_uncertainty:
        determines if to calculate relative (default) or absolute
        uncertainty
    :param orig_scene_dir:
        original (reference) vegetation index/parameter data to use
        for relative uncertainty calculation
    """

    # search expression to find all scenario realizations
    search_expr = f'*/Vegetation_Indices/VI*_None_10m_{vi_name}.tif'
    scenario_files = glob.glob(
        scenarios_scene_dir.joinpath(search_expr).as_posix()
    )

    # find original data for relative uncertainty calculation
    if relative_uncertainty:
        orig_file = glob.glob(
            orig_scene_dir.joinpath(search_expr[2:]).as_posix()
        )[0]
        handler = SatDataHandler()
        handler.read_from_bandstack(
            fname_bandstack=orig_file,
            in_file_aoi=sample_polygons,
            full_bounding_box_only=True
        )
        band_name = handler.get_bandnames()[0]
        orig_arr = handler.get_band(band_name)
        orig_arr = orig_arr.data

    # absolute or relative uncertainty
    prefix = 'Absolute'
    unit = '-'
    if relative_uncertainty:
        prefix = 'Relative'
        unit = '%'

    # read all scenario results at the selected point locations into a dataframe
    scenario_array_list = []
    for idx, scenario in enumerate(scenario_files):

        handler = SatDataHandler()
        handler.read_from_bandstack(
            fname_bandstack=scenario,
            in_file_aoi=sample_polygons,
            full_bounding_box_only=True
        )
        band_name = handler.get_bandnames()[0]
        scenario_array_list.append(handler.get_band(band_name))

        logger.info(f'Read scenario {idx+1}/{len(scenario_files)}')

    # map crop type to crop code
    gdf_polys = gpd.read_file(sample_polygons)
    crops = list(gdf_polys.crop_type.unique())

    # stack scenario array into 3d array
    scenarios_stacked = np.stack(scenario_array_list)

    unc_handler = deepcopy(handler)
    handler = None

    # rasterize the field polygons (using their int crop code)
    snap_band = unc_handler.get_bandnames()[0]
    unc_handler.add_bands_from_vector(
        in_file_vector=sample_polygons,
        snap_band=snap_band,
        attribute_selection=['crop_code'],
        blackfill_value=-999.
    )
    crop_code_arr = unc_handler.get_band('crop_code')

    # loop over crop types and calculate the absolute uncertainty in relation to the
    # the number of scenarios considered on the x axis
    for crop in crops:

        logger.info(f'Plotting number of scenarios vs. {prefix.lower()} uncertainty for crop "{crop}"')

        # get crop code
        crop_code = gdf_polys[gdf_polys.crop_type == crop]['crop_code'].iloc[0]
        # and crop mask
        crop_mask = np.zeros_like(crop_code_arr).astype('uint8')
        crop_mask[~np.isin(crop_code_arr, crop_code)] = 1
        crop_pixels = np.unique(crop_mask, return_counts=True)[1][0]

        title_str = f'{crop} - {vi_name}\nNumber of Pixels: {crop_pixels}'

        # loop over the scenarios and save min, mean, std and max uncertainty of all pixels
        # of the crop type for the current number of scenarios
        start_idx = 0
        increment = 1
        counter = start_idx + 2*increment
        n_batches = int(scenarios_stacked.shape[0]/increment)

        res_list = []
        for _ in range(1, n_batches):

            res = {}
            # get uncertainty and apply crop mask
            abs_unc = np.nanstd(scenarios_stacked[start_idx:counter,:,:], axis=0).data

            # convert to relative uncertainty using the original data
            if relative_uncertainty:
                abs_unc = np.divide(abs_unc, orig_arr) * 100.

            abs_unc[crop_mask == 1]= np.nan

            res['scenario_number'] = counter
            res['_min'] = np.nanmin(abs_unc)
            res['_max'] = np.nanmax(abs_unc)
            res['_q05'] = np.nanquantile(abs_unc, 0.05)
            res['_q95'] = np.nanquantile(abs_unc, 0.95)
            res['_mean'] = np.nanmean(abs_unc)
            res['_median'] = np.nanmedian(abs_unc)
            res['_std_plus'] = res['_mean'] + np.nanstd(abs_unc)
            res['_std_minus'] = res['_mean'] - np.nanstd(abs_unc)

            res_list.append(res)
            counter += increment

        df = pd.DataFrame(res_list)
        fname_csv = out_dir.joinpath(
            f'{vi_name}_{crop.replace(" ","-")}_scenarios-{prefix.lower()}-uncertainty.csv'
        )
        df.to_csv(fname_csv, index=False)

        # plot the derived uncertainty values on the y-axis and the number of scenarios
        # used on the x-axis
        fig = plt.figure(num=1, figsize=(15,10))
        ax = fig.add_subplot(111)

        # absolute uncertainty
        ax.plot(df.scenario_number, df._mean, color='blue', label=f'Mean')
        ax.fill_between(
            x=df.scenario_number,
            y1=df._min,
            y2=df._max,
            label='Min-Max Spread',
            color='orange',
            alpha=0.5
        )
        ax.fill_between(
            x=df.scenario_number,
            y1=df._std_minus,
            y2=df._std_plus,
            label=r'$\pm$ 1 Stddev',
            color='red',
            alpha=0.45
        )

        ax.fill_between(
            x=df.scenario_number,
            y1=df._q05,
            y2=df._q95,
            label=r'5-95% Quantile Spread',
            color='red',
            alpha=0.3
        )

        ax.set_xlabel('Number of Scenarios', fontsize=24)
        ax.set_ylabel(f'{prefix} Uncertainty (k=1) [{unit}]', fontsize=24)

        ax.set_title(title_str, fontsize=26)
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=4, fontsize=20)

        fname_out = out_dir.joinpath(f'{vi_name}_{crop.replace(" ","")}_scenarios-{prefix.lower()}-uncertainty.png')
        fig.savefig(fname_out, bbox_inches='tight', dpi=300)
        plt.close(fig)

        logger.info(f'Plotted number of scenarios vs. {prefix.lower()} uncertainty for crop {crop}')
    


if __name__ == '__main__':

    scenes = [
        'S2A_MSIL1C_20190328T102021_N0207_R065_T32TMT_20190328T154025',
        'S2A_MSIL1C_20190530T103031_N0207_R108_T32TMT_20190530T123429',
        'S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651'
    ]

    sample_polygons = Path('../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    vi_names = ['NDVI', 'EVI', 'LAI']
    out_dir = Path('../S2A_MSIL2A_Analysis/how_many_scenarios')
    orig_datasets_dir = Path('../S2A_MSIL1C_orig')

    relative_uncertainty = False

    for scene in scenes:

        scenario_scene_dir = Path(f'../S2A_MSIL1C_RUT-Scenarios/{scene}')
        out_dir_scene = out_dir.joinpath(scene)

        if not out_dir_scene.exists():
            out_dir_scene.mkdir()

        orig_scene_dir = None
        if relative_uncertainty:
            scene_splitted = scene.split('_')
            search_expr = scene_splitted[0] + '_MSIL2A_' + scene_splitted[2]
            orig_scene_dir = glob.glob(
                orig_datasets_dir.joinpath(f'{search_expr}*.VIs').as_posix()
            )[0]
            orig_scene_dir = Path(orig_scene_dir)

        for vi_name in vi_names:

            out_dir_vi = out_dir_scene.joinpath(vi_name)
            if not out_dir_vi.exists():
                out_dir_vi.mkdir()

            plot_uncertainty_number_of_scenarios(
                scenarios_scene_dir=scenario_scene_dir,
                sample_polygons=sample_polygons,
                vi_name=vi_name,
                out_dir=out_dir_vi,
                orig_scene_dir=orig_scene_dir,
                relative_uncertainty=relative_uncertainty
            )
        