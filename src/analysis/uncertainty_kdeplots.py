'''
Kernel-density plots of uncertainty distributions for EVI, NDVI and GLAI
'''

import geopandas as gpd
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from pathlib import Path
from typing import Dict, List

plt.style.use('ggplot')
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

def plot_uncertainty_kdeplots(
        res_dir: Path,
        vis: List[str],
        crop_code_mapping: Dict[int,str]
    ) -> None:
    """
    Generates the KDE plots

    :param res_dir:
        directory with uncertainty estimates
    :param vis:
        list of indices and parameters to analyze
    :param crop_code_mapping:
        mapping of crop codes to crop names
    """
    fig, axes = plt.subplots(1, 3, figsize=(15,8))
    for idx, vi in enumerate(vis):
        vi_res_file = res_dir.joinpath(vi).joinpath(f'{vi}_crops.csv')
        vi_df = pd.read_csv(vi_res_file, usecols=['crop_code', vi, f'{vi}_unc'])
        vi_df['rel_unc'] = vi_df[f'{vi}_unc'] / vi_df[vi] * 100
        vi_df
        vi_df['crop'] = vi_df['crop_code'].apply(
            lambda x,
            crop_code_mapping=crop_code_mapping: crop_code_mapping[x]
        )
        vi_df['crop'] = vi_df['crop'].apply(lambda x: 'Rapeseed' if x == 'Canola' else x)
        vi_df['crop'] = vi_df['crop'].apply(lambda x: 'Grain Maize' if x == 'Corn' else x)
        vi_df['crop'] = vi_df['crop'].apply(lambda x: 'Permament Grassland' if x == 'Permament Grasland' else x)
        vi_df['crop'] = vi_df['crop'].apply(lambda x: 'Extensively Used Grassland' if x == 'Extensively Used Grasland' else x)

        # ensure colors are the same for all plots
        sns_palette = sns.color_palette('husl', 11)
        palette = {
            'Permament Grassland': sns_palette[0],
            'Soybean': sns_palette[1],
            'Silage Maize': sns_palette[2],
            'Sugar Beet': sns_palette[3],
            'Extensively Used Grassland': sns_palette[4],
            'Winter Wheat': sns_palette[5],
            'Grain Maize': sns_palette[6],
            'Winter Barley': sns_palette[7],
            'Rapeseed': sns_palette[8],
            'Potato': sns_palette[9],
            'Sunflower': sns_palette[10]
        }

        # get median uncertainties
        median_unc_overall = vi_df['rel_unc'].median()
        median_unc_crops = vi_df.groupby(by='crop').agg('median')['rel_unc']
        median_unc_crops = median_unc_crops.append(pd.Series({'overall': median_unc_overall}))
        # and save them to file
        median_unc_crops.to_csv(res_dir.joinpath(f'{vi}_median_uncertainty.csv'))

        # cut uncertainties larger than 100% for plotting to get also a
        # KDE with more vertices
        vi_df = vi_df[(vi_df['rel_unc'] >= 0) & (vi_df['rel_unc'] <= 100)].copy()

        sns.kdeplot(x='rel_unc', data=vi_df, hue='crop', ax=axes[idx],
                    fill=True, alpha=.5, multiple='stack', palette=palette,
                    hue_order=list(palette.keys()))
        axes[idx].set_xscale('log')
        axes[idx].set_xlabel('Relative Uncertainty (k=1) [%]', fontsize=18)
        axes[idx].title.set_text(vi)
        axes[idx].set_ylim(0,.8)
        axes[idx].set_xlim(0,100)
        if idx < 2:
            axes[idx].get_legend().remove()
        if idx == 2:
            axes[idx].yaxis.set_label_position("right")
            axes[idx].yaxis.tick_right()
        if idx == 1:
            axes[idx].set_ylabel('')
            axes[idx].yaxis.set_ticklabels([])
            axes[idx].yaxis.set_ticks_position('none')
        

    fpath_fig = res_dir.joinpath('vis_rel_unc_kdeplots.png')
    fig.savefig(fpath_fig, dpi=300, bbox_inches='tight')

if __name__ == '__main__':

    res_dir = Path(
        '../../S2_TimeSeries_Analysis'
    )
    vis = ['EVI', 'NDVI', 'GLAI']

    # # median field parcel size without buffer
    # shp = Path('../../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    # gdf = gpd.read_file(shp)
    # median_area_org = gdf.geometry.area.agg(by='median')
    # # to hectar
    # median_area_org *= 1 / (100*100)

    # shapefile with crop type information for the single field parcels
    shapefile_crops = Path('../../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    column_crop_code = 'crop_code'
    column_crop_names = 'crop_type'

    # define mapping of "crop_codes" to crop types
    gdf = gpd.read_file(shapefile_crops)

    # get area per crop in ha
    gdf.dropna(inplace=True)
    # median_area_buf = gdf.geometry.area.agg(by='median')
    # mean_area_buf = gdf.geometry.area.agg(by='mean')
    # # to hectar
    # median_area_buf *= 1 / (100*100)
    dissolved = gdf.dissolve(by='crop_type')
    dissolved['area'] = np.round(dissolved.geometry.area / (100 * 100),1)
    dissolved['area'].to_csv(res_dir.joinpath('crop_types_area.csv'))
    # print(dissolved['area'].to_latex())

    crop_code_mapping = dict(list(gdf.groupby([column_crop_code, column_crop_names]).groups))

    # from agrisatpy.core.sensors import Sentinel2
    #
    # test_scene = Path('/home/graflu/Documents/uncertainty/S2_MSIL1C_orig/S2A_MSIL1C_20190216T102111_N0207_R065_T32TMT_20190216T122039.SAFE')
    # # get number of S2 pixels
    # res = []
    # for crop in crop_code_mapping:
    #     crop_parcels = dissolved[dissolved.crop_code == crop]
    #     s2 = Sentinel2().from_safe(
    #         in_dir=test_scene, band_selection=['B8A','B02'],
    #         vector_features=crop_parcels
    #     )
    #     res.append(
    #         {'crop': crop_code_mapping[crop], 'pixel_count_10m': s2['blue'].values.count()}
    #     )
    #
    # df = pd.DataFrame(res)

    plot_uncertainty_kdeplots(res_dir, vis, crop_code_mapping)
    