'''

'''

from pathlib import Path
import geopandas as gpd
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.style.use('ggplot')
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

def plot_uncertainty_kdeplots(res_dir, vis, crop_code_mapping):

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

    # shapefile with crop type information for the single field parcels
    shapefile_crops = Path('../../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    column_crop_code = 'crop_code'
    column_crop_names = 'crop_type'

    # define mapping of "crop_codes" to crop types
    gdf = gpd.read_file(shapefile_crops)

    # get area per crop in ha
    gdf.dropna(inplace=True)
    dissolved = gdf.dissolve(by='crop_type')
    dissolved['area'] = np.round(dissolved.geometry.area / (100 * 100),1)
    dissolved['area'].to_csv(res_dir.joinpath('crop_types_area.csv'))
    print(dissolved['area'].to_latex())

    crop_code_mapping = dict(list(gdf.groupby([column_crop_code, column_crop_names]).groups))

    plot_uncertainty_kdeplots(res_dir, vis, crop_code_mapping)
    