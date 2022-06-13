'''
Created on Jun 13, 2022

@author: graflu
'''

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from datetime import date
from pathlib import Path
from typing import List
from agrisatpy.core.sensors.sentinel2 import Sentinel2
from agrisatpy.utils.sentinel2 import get_S2_acquistion_date_from_safe

plt.style.use('ggplot')

def plot_ts(
        gdf: gpd.GeoDataFrame,
        orig_scenes_dir: Path,
        black_list: List[date],
        out_dir: Path
    ) -> None:
    """
    Plots the median NDVI and central 90% of NDVI values
    calculated from L1C Sentinel-2 data

    :param gdf:
        GeoDataFrame with dissolved field parcel geometries per
        crop type
    :param orig_scenes_dir:
        directory where the original S2 L1C is stored
    :param black_list:
        list of dates to exclude because of high snow and cloud
        cover
    :param out_dir:
        directory where to save the results to
    """

    # loop over crops - one plot per crop type
    for crop in gdf.crop_type.unique():
        # loop over scenes
        stats = []
        for scene in orig_scenes_dir.glob('*MSIL1C*.SAFE'):
            sensing_date = get_S2_acquistion_date_from_safe(scene)
            if sensing_date in black_list:
                continue
            # read B04 (red) and B08 (nir)
            handler = Sentinel2.from_safe(
                in_dir=scene,
                band_selection=['B04', 'B08'],
                vector_features=gdf[gdf.crop_type == crop].copy()
            )
            # calculate the NDVI
            ndvi = handler.calc_si('NDVI')
            ndvi = ndvi.data
            ndvi_stats = {
                'date': sensing_date,
                'median': np.nanmedian(ndvi),
                'q05': np.nanquantile(ndvi, .05),
                'q95': np.nanquantile(ndvi, .95)
            }
            stats.append(ndvi_stats)
        df = pd.DataFrame(stats)
        df = df.sort_values(by='date')
        df.date = pd.to_datetime(df.date)

        # plot the time series
        f, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,7))
        ax.plot(
            df['date'],
            df['median'],
            marker='x',
            markersize=10,
            label='Median',
            color='blue'
        )
        ax.fill_between(
            x=df.date,
            y1=df['q05'],
            y2=df['q95'],
            color='orange',
            alpha=0.6,
            label='Central 90%'
        )
        ax.set_ylim(0,1)
        ax.set_ylabel('NDVI [-]', fontsize=16)
        ax.legend(fontsize=16)
        ax.set_title(crop)
        # save plot
        fname = out_dir.joinpath(f'{crop}_L1C_NDVI_time-series.png')
        f.savefig(fname, dpi=300, bbox_inches='tight')
    

if __name__ == '__main__':

    # define regions of interest (different crop types) + forest + settlement
    gdf = gpd.read_file('../../shp/areas_of_interest_uncertainty_contributors_dissolved.gpkg')

    # directory with original S2 scenes
    orig_scenes_dir = Path('/home/graflu/Documents/uncertainty/S2_MSIL1C_orig')

    # directory where to store plots
    out_dir = Path('../../S2_TimeSeries_Analysis')

    # black-list scenes with clouds or snow
    black_list = [date(2019,2,14), date(2019,2,16), date(2019,5,30), date(2019,7,16)]

    plot_ts(gdf, orig_scenes_dir, black_list, out_dir)
