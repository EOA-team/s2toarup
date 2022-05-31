'''
Plots phenological curves and indicators for sample pixels
'''

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pytz
import random

from datetime import date, datetime, timedelta, timezone
from pathlib import Path

plt.style.use('ggplot')

def plot_samples(
        sample_pixels_dir: Path,
        out_dir: Path
    ):
    """
    Plots randomly selected pixel time series
    """
    # read pixel time series
    vis = ['EVI', 'NDVI', 'GLAI']
    gdf_dict = dict.fromkeys(vis)
    colors = {'EVI': 'blue', 'NDVI': 'orange', 'GLAI': 'green'}
    for vi in vis:
        sample_pixels = sample_pixels_dir.joinpath(f'sample_pixels_{vi}.gpkg')
        gdf_dict[vi] = gpd.read_file(sample_pixels)
    # get available crop types (11)
    crop_types = gdf_dict[vis[0]].crop_type.unique()

    # single pixel plot (sample pixel to explain LSP metrics)
    f2, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,7))
    crop_type = 'Sunflower'
    vi = 'GLAI'
    crop_gdf = gdf_dict[vi][gdf_dict[vi].crop_type == crop_type].copy()
    pixels = crop_gdf.point.unique()
    sel_pidx = 1
    pixel_gdf = crop_gdf[crop_gdf.point == pixels[sel_pidx]].copy()
    # set time zone
    pixel_gdf['time'] = pd.to_datetime(pixel_gdf['time'])
    pixel_gdf['time'] = pixel_gdf['time'].dt.tz_localize(timezone.utc)
    my_timezone = pytz.timezone('Europe/Berlin')
    pixel_gdf['time'] = pixel_gdf['time'].dt.tz_convert(my_timezone)
    # plot time series
    ax.plot(
        pixel_gdf['time'],
        pixel_gdf[vi],
        color='blue',
        label='Smoothed Daily Time Series'
    )
    sos_idx = timedelta(pixel_gdf['sos_times'].iloc[0]-1)
    sos_date = date(2019,1,1) + sos_idx
    ax.vlines(sos_date, 0, 6, color='green', label='SOS')
    eos_idx = timedelta(pixel_gdf['eos_times'].iloc[0]-1)
    eos_date = date(2019,1,1) + eos_idx
    ax.vlines(eos_date, 0, 6, linestyle='dashed', color='green', label='EOS')
    ax.set_ylim(0,6)
    ax.set_ylabel('GLAI $[m^2/m^2]$', fontsize=18)
    ax.set_xlabel('Time [YYYY-MM]', fontsize=18)
    ax.hlines(1.35, sos_date, eos_date, label='LOS')
    pixel_gdf['baseline'] = 1.35
    offset_idx = pixel_gdf.time.min().timetuple().tm_yday
    sos_idx, eos_idx = sos_idx.days, eos_idx.days
    sos_idx -= offset_idx
    eos_idx -= offset_idx 
    ax.fill_between(
        pixel_gdf['time'][sos_idx+1:eos_idx+1],
        pixel_gdf[vi][sos_idx+1:eos_idx+1],
        pixel_gdf.baseline[sos_idx+1:eos_idx+1], color='lightgreen',
        label='Growing Season'
    )
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=14)
    fname = f'lsp_sample_time_series.png'
    fpath = out_dir.joinpath(fname)
    f2.savefig(fpath, bbox_inches='tight')
    
    
    f, axes = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True, figsize=(15,10))
    # loop over crops and select a single pixel
    overall_counter = 0
    row_counter, col_counter = 0, 0
    for crop_type in crop_types:
        # get single pixel from crop type
        for idx, vi in enumerate(vis):
            crop_gdf = gdf_dict[vi][gdf_dict[vi].crop_type == crop_type].copy()
            pixels = crop_gdf.point.unique()
            if idx == 0:
                sel_pidx = random.randint(0, pixels.shape[0]-1)
            try:
                pixel_gdf = crop_gdf[crop_gdf.point == pixels[sel_pidx]].copy()
            except Exception as e:
                print(e)

            # set time zone
            pixel_gdf['time'] = pd.to_datetime(pixel_gdf['time'])
            pixel_gdf['time'] = pixel_gdf['time'].dt.tz_localize(timezone.utc)
            my_timezone = pytz.timezone('Europe/Berlin')
            pixel_gdf['time'] = pixel_gdf['time'].dt.tz_convert(my_timezone)
            # plot time series
            if vi != 'GLAI':
                axes[row_counter, col_counter].plot(
                    pixel_gdf['time'],
                    pixel_gdf[vi],
                    color=colors[vi],
                    label=vi
                )
                axes[row_counter, col_counter].set_ylim(0,1)
            else:
                ax2 = axes[row_counter, col_counter].twinx()
                ax2.plot(
                    pixel_gdf['time'],
                    pixel_gdf[vi],
                    color=colors[vi],
                    label=vi,
                )
                ax2.set_ylim(0,7)
                if col_counter == 3:
                    ax2.set_ylabel(
                        r'GLAI [$m^2/m^2$]',
                        fontsize=14,
                        rotation=270,
                        labelpad=20
                    )
            # plot SOS and EOS
            sos_date = date(2019,1,1) + timedelta(pixel_gdf['sos_times'].iloc[0])
            axes[row_counter, col_counter].vlines(sos_date, 0, 1, color=colors[vi])
            eos_date = date(2019,1,1) + timedelta(pixel_gdf['eos_times'].iloc[0])
            axes[row_counter, col_counter].vlines(eos_date, 0, 1, color=colors[vi], linestyle='dashed')
            plt.setp(axes[row_counter,col_counter].get_xticklabels(), rotation=45)
            crop_name = crop_type
            if crop_name == 'Canola':
                crop_name = 'Rapeseed'
            if crop_name == 'Corn':
                crop_name = 'Grain Maize'
            axes[row_counter, col_counter].set_title(crop_name)
            if col_counter == 0:
                axes[row_counter, col_counter].set_ylabel('Vegetation Index [-]', fontsize=14)
        col_counter += 1
        if row_counter == 2 and col_counter == 3:
            f.delaxes(axes[row_counter,col_counter])
        if overall_counter in [3,7,11]:
            row_counter += 1
            col_counter = 0
        overall_counter += 1
    # set legend
    lines, labels = axes[2,2].get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(
        lines + lines2, labels + labels2,
        loc='center left',
        bbox_to_anchor=(1.3, 0.5),
        fontsize=14
    )
    ax2.set_ylabel(
        r'GLAI [$m^2/m^2$]',
        fontsize=14,
        rotation=270,
        labelpad=20
    )

    # save plot to file
    now = datetime.now().strftime('%d%m%d_%H%M%S')
    fname = f'{now}_sample_time_series.png'
    fpath = out_dir.joinpath(fname)
    f.savefig(fpath, bbox_inches='tight')

        

if __name__ == '__main__':
    sample_pixels_dir = Path('../../S2_TimeSeries_Analysis')
    out_dir = sample_pixels_dir

    plot_samples(sample_pixels_dir, out_dir)
