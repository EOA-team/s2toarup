'''
Check how many scenario runs are necessary
'''

import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from agrisatpy.io import SatDataHandler

plt.style.use('ggplot')


def plot_uncertainty_number_of_scenarios(
        scenarios_scene_dir: Path,
        sample_points: Path,
        vi_name: str,
        orig_scene_dir: Path,
        out_dir: Path
    ):

    # search expression to find all scenario realizations
    search_expr = f'*/Vegetation_Indices/VI*_None_10m_{vi_name}.tif'
    scenario_files = glob.glob(
        scenarios_scene_dir.joinpath(search_expr).as_posix()
    )

    # read original (reference) data as well to provide relative uncertainties in addition
    orig_file = glob.glob(
        orig_scene_dir.joinpath(search_expr[2:]).as_posix()
    )[0]
    gdf_orig = SatDataHandler.read_pixels(
        point_features=sample_points,
        raster=orig_file
    )

    # read all scenario results at the selected point locations into a dataframe
    for idx, scenario in enumerate(scenario_files):

        if idx == 0:
            gdf = SatDataHandler.read_pixels(
                point_features=sample_points,
                raster=scenario
            )
        else:
            gdf = SatDataHandler.read_pixels(
                point_features=gdf,
                raster=scenario
            )
        gdf.rename(columns={vi_name: f'{vi_name}_{idx+1}'}, inplace=True)

        if idx == 20:
            break

    # loop over points and plot the derived absolute uncertainty in relation to the
    # the number of scenarios considered. Increase the number considered by steps of 10
    # and plot the result
    point_ids = list(gdf.id.unique())

    for point_id in point_ids:
        
        point_gdf = gdf[gdf.id == point_id].copy()
        # get the original pixel value
        reference_value = gdf_orig[gdf_orig.id == point_id][vi_name].iloc[0]

        # get scenario results at this location
        mask = point_gdf.columns.str.contains(f'{vi_name}_*')
        scenario_values = point_gdf.loc[:,mask].values

        # calculate the absolute uncertainties for increasing number of scenarios
        scenario_number = []
        abs_unc = []
        start_idx = 0
        increment = 1
        n_batches = int(scenario_values.shape[1]/increment)

        counter = 10
        for _ in range(n_batches):
            
            scenario_number.append(counter)
            res = np.std(scenario_values[0,start_idx:counter])
            abs_unc.append(res)

            # increment the counter
            counter += increment

        # plot the derived uncertainty values on the y-axis and the number of scenarios
        # used on the x-axis
        fig = plt.figure(num=1, dpi=150)
        ax = fig.add_subplot(111)

        coord_str = f'x={np.round(point_gdf.geometry.x.iloc[0],0)}m, y={np.round(point_gdf.geometry.y.iloc[0],0)}m'
        crop_type = point_gdf.crop_type.iloc[0]
        title_str = f'{crop_type} {vi_name}={np.round(reference_value,4)}\n {coord_str}'
    
        # absolute uncertainty
        ax.plot(scenario_number, abs_unc, label=f'{vi_name} Uncertainty')
        ax.set_xlabel('Number of scenarios', fontsize=16)
        ax.set_ylabel('Absolute Uncertainty (k=1) [-]', fontsize=16)

        # # relative uncertainty
        # ax2 = ax.twinx()
        # rel_unc = [(x/reference_value)*100 for x in abs_unc]
        # ax2.plot(scenario_number, rel_unc, '-x', label=f'{vi_name} rel. Uncertainty')
        # ax2.set_ylabel('Relative Uncertainty (k=1) [%]', fontsize=16, rotation=270)
        
        ax.set_title(title_str, fontsize=14)
        ax.legend(fontsize=14)

        fname_out = out_dir.joinpath(f'{vi_name}_{crop_type}_{coord_str.replace(", ","_").replace("=","")}.png')
        fig.savefig(fname_out, bbox_inches='tight')
        plt.close(fig)

    # save dataframe to CSV
    fname_out_csv = out_dir.joinpath(f'{vi_name}_{sample_points.name.replace(".shp",".csv")}')
    gdf.to_csv(fname_out_csv, index=False)


if __name__ == '__main__':

    scenes = [
        'S2A_MSIL1C_20190328T102021_N0207_R065_T32TMT_20190328T154025',
        'S2A_MSIL1C_20190530T103031_N0207_R108_T32TMT_20190530T123429',
        'S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651'
    ]

    sample_points = Path('../shp/ZH_Points_2019_EPSG32632_selected-crops.shp')
    vi_names = ['NDVI', 'EVI']
    out_dir = Path('../S2A_MSIL2A_Analysis/how_many_scenarios')
    orig_datasets_dir = Path('../S2A_MSIL1C_orig')

    for scene in scenes:

        scenario_scene_dir = Path(f'../S2A_MSIL1C_RUT-Scenarios/{scene}')
        out_dir_scene = out_dir.joinpath(scene)

        if not out_dir_scene.exists():
            out_dir_scene.mkdir()

        scene_splitted = scene.split('_')
        search_expr = scene_splitted[0] + '_MSIL2A_' + scene_splitted[2]
        orig_scene_dir = glob.glob(
            orig_datasets_dir.joinpath(f'{search_expr}*.VIs').as_posix()
        )[0]

        for vi_name in vi_names:

            out_dir_vi = out_dir_scene.joinpath(vi_name)
            if not out_dir_vi.exists():
                out_dir_vi.mkdir()

            plot_uncertainty_number_of_scenarios(
                scenarios_scene_dir=scenario_scene_dir,
                sample_points=sample_points,
                vi_name=vi_name,
                orig_scene_dir=Path(orig_scene_dir),
                out_dir=out_dir_vi
            )
        