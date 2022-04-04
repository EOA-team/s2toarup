'''
Created on Apr 4, 2022

@author: graflu
'''

from pathlib import Path
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

from agrisatpy.core.band import Band, GeoInfo

plt.style.use('ggplot')


if __name__ == '__main__':

    lsp_res_dir = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/01_Uncertainty/ESCH/scripts_paper_uncertainty/S2_TimeSeries_Analysis'
    )

    vis = ['EVI', 'NDVI', 'GLAI']
    runs = ['uncorrelated', 'fully_correlated']
    metrics = ['sos_times', 'eos_times']

    # read all data into a large data frame
    for run in runs:
        res = []
        f, axes = plt.subplots(nrows=4, ncols=3)
        for vidx, vi in enumerate(vis):
            search_path = lsp_res_dir.joinpath(f'{vi}/{run}/Uncertainty_Maps/selected_crops')
            counter = 0
            for midx, metric in enumerate(metrics):
                fpath = next(search_path.glob(f'{vi}_{metric}_*_data.csv'))
                df = pd.read_csv(fpath)
                # drop grassland pixels (no meaningful LSP metrics)
                drop_idx = df[
                    df.crop.isin(['Extensively Used Grasland', 'Permament Grasland'])
                ].index
                df.drop(index=drop_idx, inplace=True)

                # cutoff uncertainties larger than 50 days
                unc = df[f'{metric} Uncertainty'].copy()
                unc[unc > 50] = 50

                # plot histogram of LSP metric uncertainty values (all crops)
                unc.hist(
                    bins=100,
                    ax=axes[counter, vidx],
                    density=True
                )
                axes[counter,vidx].set_xlim(0,50)
                metric_label = metric.split('_')[0].upper()
                title = f'{vi} {metric_label}\nUncertainties'
                axes[counter,vidx].set_title(title, fontdict={'fontsize':9})
                if vidx == 0:
                    axes[counter, vidx].set_ylabel('Relative Frequency')

                # plot map
                gdf = df[['x', 'y', f'{metric} Uncertainty', 'crop']].copy()
                gdf = gpd.GeoDataFrame(gdf, geometry=gpd.points_from_xy(gdf.x, gdf.y))
                gdf.set_crs(epsg=32632, inplace=True)
                # convert vector geometry to Band object (rasterize)
                ulx, uly = gdf.total_bounds[0], gdf.total_bounds[-1]
                geo_info = GeoInfo(epsg=32632, ulx=ulx, uly=uly, pixres_x=10, pixres_y=10)
                band = Band.from_vector(
                    vector_features=gdf,
                    geo_info=geo_info,
                    band_name_src=f'{metric} Uncertainty',
                    band_name_dst=f'{metric} Uncertainty',
                    nodata_dst=np.nan,
                    dtype_src='float32'
                )
                
                counter += 1
                band.plot(
                    colormap='Oranges',
                    vmin=0,
                    vmax=50,
                    ax=axes[counter, vidx]
                )
                axes[counter,vidx].set_title('')
                axes[counter,vidx].set_xlabel('')
                axes[counter,vidx].set_ylabel('')

               
                

                df['vi'] = vi
                res.append(df)
        df_large = pd.concat(res)

        df_large
                
