'''
Check how many scenario runs are necessary
'''

import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from agrisatpy.io import SatDataHandler


def plot_uncertainty_number_of_scenarios(
        scenarios_scene_dir: Path,
        sample_points: Path,
        vi_name: str,
        out_dir: Path
    ):
    # search expression to find all scenario realizations
    search_expr = f'*/Vegetation_Indices/VI*_{vi_name}.tif'
    scenario_files = glob.glob(
        scenarios_scene_dir.joinpath(search_expr)
    )

    # read all scenario results at the selected point locations into a dataframe
    for idx, scenario in scenario_files:

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
        gdf.rename(columns={vi_name: f'{vi_name}_{scenario+1}'})

    # loop over points and plot the derived absolute uncertainty in relation to the
    # the number of scenarios considered. Increase the number considered by steps of 10
    # and plot the result
    point_ids = list(gdf.id.unique())

    for point_id in point_ids:
        
        point_gdf = gdf[gdf.id == point_id].copy()

        # get scenario results at this location
        mask = point_gdf.columns.str.contains(f'{vi_name}_*')
        scenario_values = point_gdf.loc[:,mask].values()

        # calculate the absolute uncertainties for increasing number of scenarios
        scenario_number = []
        abs_unc = []
        start_idx = 0
        increment = 10
        n_batches = int(scenario_values.shape[0]/increment)

        counter = 10
        for _ in range(n_batches):
            
            scenario_number.append(counter)
            res = np.std(scenario_values[start_idx:counter])
            abs_unc.append(res)

            # increment the counter
            counter += increment

        # plot the derived uncertainty values on the y-axis and the number of scenarios
        # used on the x-axis
        fig = plt.figure(dpi=150)
        ax = fig.add_subplot(111)

        coord_str = f'x={np.round(point_gdf.geometry.x.iloc[0],0)}m, y={np.round(point_gdf.geometry.y.iloc[0],0)}m'
        crop_type = point_gdf.crop_type.iloc[0]
        title_str = f'{crop_type}\n {coord_str}'
    
        ax.plot(scenario_number, abs_unc, 'o', label=f'{vi_name} Uncertainty')
        ax.set_xlabel('Number of scenarios', fontsize=16)
        ax.set_ylabel('Absolute Uncertainty (k=1)', fontsize=16)
        ax.title(title_str, fontsize=14)
        ax.legend(fontsize=14)

        fname_out = out_dir.joinpath(f'{vi_name}_{crop_type}_{coord_str.replace(", "),("_")}.png')
        fig.savefig(fname_out, bbox_inches='tight')

    # save dataframe to CSV
    fname_out_csv = out_dir.joinpath(f'{vi_name}_{sample_points.name.replace(".shp",".csv")}')
    gdf.to_csv(fname_out_csv, index=False)


if __name__ == '__main__':
    
    scenario_scene_dir = Path('../S2A_MSIL1C_RUT-Scenarios/S2A_MSIL1C_20190328T102021_N0207_R065_T32TMT_20190328T154025')
    sample_points = Path('../shp/ZH_Points_2019_EPSG32632_selected-crops.shp')
    vi_names = ['NDVI', 'EVI']
    out_dir = Path('../S2A_MSIL2A_Analysis/how_many_scenarios')
    
    