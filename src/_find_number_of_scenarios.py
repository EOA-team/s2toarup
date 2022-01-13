'''
Check how many scenario runs are necessary
'''

import glob
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from pathlib import Path
from copy import deepcopy
from agrisatpy.io import SatDataHandler

from logger import get_logger

plt.style.use('ggplot')

logger = get_logger('_how_many_scenarios')


def plot_uncertainty_number_of_scenarios(
        scenarios_scene_dir: Path,
        sample_polygons: Path,
        vi_name: str,
        out_dir: Path
    ):

    # search expression to find all scenario realizations
    search_expr = f'*/Vegetation_Indices/VI*_None_10m_{vi_name}.tif'
    scenario_files = glob.glob(
        scenarios_scene_dir.joinpath(search_expr).as_posix()
    )

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
    crops = list(gdf_polys.crop_code.unique())

    # stack scenario array into 3d array
    scenarios_stacked = np.stack(scenario_array_list)

    unc_handler = deepcopy(handler)
    handler = None

    # rasterize the field polygons (using their int crop code)
    snap_band = handler.get_bandnames()[0]
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

        logger.info(f'Plotting number of scenarios vs. derived uncertainty for crop {crop}')

        # get crop code
        crop_code = gdf_polys[gdf_polys.crop_type == crop].iloc[0]

        title_str = f'{crop} {vi_name}'

        # loop over the scenarios and save min, mean, std and max uncertainty of all pixels
        # of the crop type for the current number of scenarios
        start_idx = 0
        increment = 1
        counter = start_idx + increment
        n_batches = int(scenarios_stacked.shape[0]/increment)

        res_list = []
        for _ in range(n_batches):

            res = {}
            abs_unc = np.nanstd(scenarios_stacked[start_idx:counter,:,:], axis=0)

            # add to handler, mask all but pixels of the current crop type and calculate
            # the statistics
            mask = crop_code_arr[crop_code_arr != crop_code]

            res['scenario_number'] = counter

            res_list.append(res)
            counter += increment


        df = pd.DataFrame(res_list)

        # plot the derived uncertainty values on the y-axis and the number of scenarios
        # used on the x-axis
        fig = plt.figure(num=1, dpi=150)
        ax = fig.add_subplot(111)

        # absolute uncertainty
        ax.plot(scenario_number, abs_unc, label=f'{vi_name} Uncertainty')
        ax.set_xlabel('Number of scenarios', fontsize=16)
        ax.set_ylabel('Absolute Uncertainty (k=1) [-]', fontsize=16)

        ax.set_title(title_str, fontsize=14)
        ax.legend(fontsize=14)

        fname_out = out_dir.joinpath(f'{vi_name}_{crop}_scenarios-uncertainty.png')
        fig.savefig(fname_out, bbox_inches='tight')
        plt.close(fig)

        logger.info(f'Plotted number of scenarios vs. derived uncertainty for crop {crop}')
    


if __name__ == '__main__':

    scenes = [
        'S2A_MSIL1C_20190328T102021_N0207_R065_T32TMT_20190328T154025',
        'S2A_MSIL1C_20190530T103031_N0207_R108_T32TMT_20190530T123429',
        'S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651'
    ]

    sample_polygons = Path('../shp/ZH_Points_2019_EPSG32632_selected-crops_buffered.shp')
    vi_names = ['NDVI', 'EVI']
    out_dir = Path('../S2A_MSIL2A_Analysis/how_many_scenarios')
    orig_datasets_dir = Path('../S2A_MSIL1C_orig')

    for scene in scenes:

        scenario_scene_dir = Path(f'../S2A_MSIL1C_RUT-Scenarios/{scene}')
        out_dir_scene = out_dir.joinpath(scene)

        if not out_dir_scene.exists():
            out_dir_scene.mkdir()

        for vi_name in vi_names:

            out_dir_vi = out_dir_scene.joinpath(vi_name)
            if not out_dir_vi.exists():
                out_dir_vi.mkdir()

            plot_uncertainty_number_of_scenarios(
                scenarios_scene_dir=scenario_scene_dir,
                sample_polygons=sample_polygons,
                vi_name=vi_name,
                out_dir=out_dir_vi
            )
        