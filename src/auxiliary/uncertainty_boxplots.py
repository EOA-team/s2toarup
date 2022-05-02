'''

'''

from pathlib import Path
import geopandas as gpd
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('ggplot')
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

def plot_uncertainty_boxplots(res_dir, vis, crop_code_mapping):

    fig, axes = plt.subplots(1, 3, figsize=(15,8))
    for idx, vi in enumerate(vis):
        vi_res_file = res_dir.joinpath(vi).joinpath(f'{vi}_crops.csv')
        vi_df = pd.read_csv(vi_res_file, usecols=['crop_code', vi, f'{vi}_unc'])
        vi_df['rel_unc'] = vi_df[f'{vi}_unc'] / vi_df[vi] * 100
        vi_df
        vi_df['crop_name'] = vi_df['crop_code'].apply(
            lambda x,
            crop_code_mapping=crop_code_mapping: crop_code_mapping[x]
        )
        vi_df['crop_name'] = vi_df['crop_name'].apply(lambda x: 'Rapeseed' if x == 'Canola' else x)
        vi_df['crop_name'] = vi_df['crop_name'].apply(lambda x: 'Grain Maize' if x == 'Corn' else x)
        # crop_count = vi_df.crop_name.value_counts()
        # vi_df['crop_name'] = vi_df['crop_name'].apply(lambda x, crop_count=crop_count:
        #     f'{x} ({crop_count[crop_count.index == x].values[0]})'
        # )
        sns.boxplot(x='crop_name', y='rel_unc', data=vi_df, ax=axes[idx])
        plt.setp(axes[idx].get_xticklabels(), rotation=90)
        axes[idx].set_xlabel('')
        axes[idx].title.set_text(vi)
        axes[idx].set_ylim(-100,500)
        if idx == 0:
            axes[idx].set_ylabel('Relative Uncertainty (k=1) [%]', fontsize=18)
        else:
            axes[idx].set_ylabel('')
            # plt.grid(True)
            axes[idx].yaxis.set_ticklabels([])
            axes[idx].yaxis.set_ticks_position('none')

    fpath_fig = res_dir.joinpath('vis_rel_unc_boxplots.png')
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
    crop_code_mapping = dict(list(gdf.groupby([column_crop_code, column_crop_names]).groups))

    plot_uncertainty_boxplots(res_dir, vis, crop_code_mapping)
    