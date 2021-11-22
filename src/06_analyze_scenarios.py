
'''
@author:    Lukas Graf (D-USYS, ETHZ)

@purpose:   This script is used to analyze the uncertainty
            propagation outcomes before and after Sen2Cor. It generates
            tif files summarizing the scenario spread and hence
            relative standard uncertainty. In addition, it produces
            some maps useful for visually analyzing the results and
            extracts the uncertainty for the single regions of
            interest (ROI) into a handy CSV file format.
'''

import os
import glob
import pandas as pd
import datetime
from datetime import date
from pathlib import Path
import geopandas as gpd
import rasterio as rio
import rasterio.mask
import numpy as np
from scipy.stats import mode
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import colors
from typing import Optional
from typing import Tuple
from typing import Union
from typing import Dict
from rasterio.transform import Affine


# define plotting styles
plt.style.use('ggplot')
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

plt.rcParams['axes.titlesize'] = 15
gain_factor_refl = 0.01

# define S2 band properties
band_dict_l1c = {
    '10': {
        'B02': '*_B02.jp2',
        'B03': '*_B03.jp2',
        'B04': '*_B04.jp2',
        'B08': '*_B08.jp2',
    },
    '20': {
        'B05': '*_B05.jp2',
        'B06': '*_B06.jp2',
        'B07': '*_B07.jp2',
        'B8A': '*_B8A.jp2',
        'B11': '*_B11.jp2',
        'B12': '*_B12.jp2',
    },
    '60': {
        'B01': '*_B01.jp2',
        'B09': '*_B09.jp2'
    }   
}

band_dict_l2a = {
    '10': {
        'AOT': '*_AOT_10m.jp2',
        'B02': '*_B02_10m.jp2',
        'B03': '*_B03_10m.jp2',
        'B04': '*_B04_10m.jp2',
        'B08': '*_B08_10m.jp2',
        'WVP': '*_WVP_10m.jp2'
    },
    '20': {
        'AOT': '*_AOT_20m.jp2',
        'B05': '*_B05_20m.jp2',
        'B06': '*_B06_20m.jp2',
        'B07': '*_B07_20m.jp2',
        'B8A': '*_B8A_20m.jp2',
        'B11': '*_B11_20m.jp2',
        'B12': '*_B12_20m.jp2',
        'WVP': '*_WVP_20m.jp2',
        'SCL': '*_SCL_20m.jp2'
    },

    '60': {
        'AOT': '*_AOT_60m.jp2',
        'B01': '*_B01_60m.jp2',
        'B09': '*_B09_60m.jp2',
        'WVP': '*_WVP_60m.jp2'
    }   
}


def pixel_to_img_coords(
        point_coords: Tuple[Union[int,float]],
        img_transform: Affine
    ) -> Dict[str,int]:
    """
    Takes a tuple of point coordinates and translate
    it into the corresponding row and column in an image.
    Therefore, the point_coords must be in the same projection
    as the image.

    NOTE: The method does not check for out-of-bounds issues
    when the point is not located within the image!

    :param point_coords:
        coordinate tuple in the form (x,y) denoting a point
        in the image
    :param img_transform:from mpl_toolkits.axes_grid1 import make_axes_locatable

        Affine transformation defining the extent and pixel size
        of the image
    :return:
        dictionary with the image row and column.
    """
    # point coordinates
    utm_x = point_coords[0]
    utm_y = point_coords[1]
    # image coordinates
    pix_res_x = img_transform[0]
    pix_res_y = img_transform[4]
    img_ulx = img_transform[2]
    img_uly = img_transform[5]
    # map into pixel coordinates (row, column)
    sel_col = int(np.around((utm_x - img_ulx) / pix_res_x))
    sel_row = int(np.around((utm_y - img_uly) / pix_res_y))

    return {'row': sel_row, 'col': sel_col}


def analyze_scenarios_spatial(
        unc_scenario_dir: Path,
        in_file_shp: Path,
        out_dir: Path,
        processing_level: str,
        select_random_pixels: Optional[bool]=True,
        n_random_pixels: Optional[int]=5
    ):
    """
    Extracts a region of interest (ROI) from a series of
    uncertainty scenarios, stacks them, and compiles the mean and
    standard deviation among all scenarios per pixel. Thus, it is
    possible to obtain spatial uncertainty information.

    :param unc_scenario_dir:
        directory with Sen2Cor runs (L2A) or radiometric uncertainty
        scenarios (L1C)
    :param in_file_shp:
        shape file with a single ROI
    :param out_dir:
        directory where to store dataframes with extracted data
        (also for later analysis)
    :param processing_level:
        either 'L1C' or 'L2A' for identifying the S2 processing level
        (before and after atmospheric correction)
    :param select_random_pixels:
        if set to True (Default) generates histograms of randomly chosen
        pixels showing the histograms obtained from the scenarios
    :param n_random_pixels:
        if select_random_pixels is True a user-defined number of pixels
        is analyzed. The default is 5 pixels
    """
    # check processing level
    if processing_level == 'L1C':
        band_dict_s2 = band_dict_l1c
        abbrev = 'MSIL1C'
    elif processing_level == 'L2A':
        band_dict_s2 = band_dict_l2a
        abbrev = 'MSIL2A'

    # read shapefile
    gdf = gpd.read_file(in_file_shp)

    # extract random pixels if selected within the ROI bounds
    # works only if the shapefile is in the same projection as the sat data
    if select_random_pixels:
        left, bottom, right, top = gdf.iloc[0]['geometry'].bounds
        # sample across rows and columns
        col_samples = np.random.uniform(left, right, n_random_pixels)
        row_samples = np.random.uniform(bottom, top, n_random_pixels)
        coord_samples = list(zip(col_samples, row_samples))

    # loop over single bands and extract data from scenario runs
    for spatial_res in band_dict_s2.keys():

        # construct search expression for finding the files
        if processing_level == 'L1C':
            search_expr = f'*/*_{abbrev}_*/GRANULE/*/IMG_DATA/'
        elif processing_level == 'L2A':
            search_expr = f'*/*_{abbrev}_*/GRANULE/*/IMG_DATA/R{spatial_res}m/'
        band_dict = band_dict_s2[spatial_res]

        # loop over scenarios of the current band and extract ROIs
        for band in band_dict.keys():

            print(
                f'Extracting {band} from {spatial_res}m spatial resolution ({processing_level})'
            )

            # find scenario files
            scenario_files = glob.glob(
                str(unc_scenario_dir.joinpath(search_expr + band_dict[band]))
            )
            n_scenarios = len(scenario_files)

            # check CRS between ROI and raster
            sat_crs = rio.open(scenario_files[0]).crs
            gdf_reprojected = gdf.to_crs(sat_crs)
            feature = gdf_reprojected.iloc[0]

            # loop over scenario members
            n_scenarios = len(scenario_files)
            for idx, scenario_file in enumerate(scenario_files):
                with rio.open(scenario_file, 'r') as src:
                    if idx == 0:
                        meta = src.meta
                    out_band, out_transform = rio.mask.mask(
                            src,
                            [feature["geometry"]],
                            crop=True, 
                            all_touched=True # IMPORTANT!
                        )
                if idx == 0:
                    data_arr = np.empty(shape=(n_scenarios, out_band.shape[1], out_band.shape[2]))
                    meta.update(
                        {
                            "driver": "GTiff",
                            "height": out_band.shape[1],
                            "width": out_band.shape[2], 
                            "transform": out_transform,
                         }
                    )
                data_arr[idx,:,:] = out_band[0,:,:]

            # analysis: calculate min, max, mean and standard deviation per pixel
            # (does not apply to SCL - here the majority vote will be analyzed)
            if band == 'SCL':
                count = 2
            else:
                count = 5
            meta.update(
                {
                    "count": count,
                    "dtype": np.float32
                }
            )

            # if random pixels shall be analyzed then plot their histogram
            if select_random_pixels:
                # calculate the image coordinates for the coordinate tuples
                for coord_tuple in coord_samples:
                    # point coordinates
                    rand_coords = pixel_to_img_coords(
                        point_coords=coord_tuple,
                        img_transform=meta['transform']
                    )
                    # extract pixel value in all scenarios
                    pixel_vals = data_arr[:,rand_coords['row'],rand_coords['col']]
                    # plot histogram using true percentage values of reflectance
                    if band != 'SCL':
                        pixel_vals *= 0.01
                    # plot the histogram of values
                    fig = plt.figure(figsize=(6,8))
                    ax = fig.add_subplot(111)
                    
                    ax.hist(pixel_vals, bins=30, color='cornflowerblue')
                    x = int(np.around(coord_tuple[0]))
                    y = int(np.around(coord_tuple[1]))
                    epsg = sat_crs.to_epsg()

                    if band == 'SCL':
                        title_str = 'Scene Classification Layer Samples '
                        xlabel = 'Scene Classifcation Class'
                    elif band == 'WVP':
                        title_str = 'Atmospheric Water Vapor Column'
                        xlabel = 'Water Vapor Column [cm]'
                    elif band == 'AOT':
                        title_str = 'Aerosol Optical Thickness @550nm'
                        xlabel = 'Aerosol Optical Thichkness [-]'
                    else:
                        if processing_level == 'L1C':
                            title_str = r'$\rho_{TOA}$ Samples '
                            xlabel = r'$\rho_{TOA}$ Reflectance Factor [%]'
                        else:
                            title_str = r'$\rho_{BOA}$ Samples '
                            xlabel = r'$\rho_{BOA}$ Reflectance Factor [%]'

                    ax.set_title(
                        title_str + f'{band} (N={n_scenarios})\nx = {x}m, y = {y}m (EPSG:{epsg})',
                        fontsize=16
                    )
                    ax.set_ylabel('Absolute Frequency [-]', fontsize=14)
                    ax.set_xlabel(xlabel, fontsize=14)
                    if band != 'SCL':
                        avg = np.mean(pixel_vals)
                        std = np.std(pixel_vals)
                        ymax = ax.get_ylim()[1]
                        ax.vlines(
                            x=avg,
                            ymin=0,
                            ymax=ymax,
                            color='firebrick',
                            linewidth=4,
                            label='Mean'
                        )
                        ax.vlines(
                            x=avg-std,
                            ymin=0,
                            ymax=ymax,
                            color='lightcoral',
                            linewidth=4,
                            linestyle='dashed',
                            label=r'$\pm$ 1 Std-Dev'
                        )
                        ax.vlines(
                            x=avg+std,
                            ymin=0,
                            ymax=ymax,
                            color='lightcoral',
                            linewidth=4,
                            linestyle='dashed',
                        )
                        ax.legend(fontsize=14)
                    # save plots
                    histo_dir = out_dir.joinpath('pixel_histograms')
                    if not histo_dir.exists(): histo_dir.mkdir()
                    fname_plot = histo_dir.joinpath(
                        f'{processing_level}_{band}_{spatial_res}m_{x}_{y}_histogram.png'
                    )
                    fig.savefig(fname_plot, bbox_inches='tight')
                    plt.close(fig)

            fname = f'{processing_level}_{band}_{spatial_res}m_{in_file_shp.name.split(".")[0]}.tif'

            if band == 'SCL':
                
                scl_analysis = mode(data_arr, axis=0)

                with rio.open(out_dir.joinpath(fname), 'w', **meta) as dst:

                    # majority vote
                    dst.set_band_description(1, 'majority_vote')
                    dst.write_band(1, scl_analysis.mode[0,:,:])

                    # confidence (percentage of scenario members voting for the decision)
                    dst.set_band_description(2, 'confidence')
                    dst.write_band(2, scl_analysis.count[0,:,:]/n_scenarios*100.)

            else:

                with rio.open(out_dir.joinpath(fname), 'w', **meta) as dst:
                    # min
                    dst.set_band_description(1, 'min')
                    dst.write_band(1, np.nanmin(data_arr, axis=0))
                    
                    # max
                    dst.set_band_description(2, 'max')
                    dst.write_band(2, np.nanmax(data_arr, axis=0))
    
                    # mean
                    dst.set_band_description(3, 'mean')
                    dst.write_band(3, np.nanmean(data_arr, axis=0))
    
                    # absolute standard deviation
                    dst.set_band_description(4, 'abs_stddev')
                    dst.write_band(4, np.nanstd(data_arr, axis=0))
    
                    # standard uncertainty -> normalize stack of scenarios
                    dst.set_band_description(5, 'rel_std_unc')
                    rel_std = np.nanstd(data_arr, axis=0) / np.nanmean(data_arr, axis=0)
                    dst.write_band(5, rel_std * 100)


def unc_maps(
        analysis_results_l1c: str,
        analysis_results_l2a: str,
        analysis_results_aot: Path,
        analysis_results_wvp: Path,
        analysis_results_scl: Path,
        out_dir: Path
    ) -> None:
    """
    Maps raster values of L1C and L2A uncertainty values to reveal
    spatial pattern of uncertainty and their land cover/ use dependency

    :param l1c_scenarios:
        basename of tif files with L1C RUT output
    :param l2a_scenarios:
        basename of the atmospheric correction uncertainty output
    :param analysis_results_aot:
        file with aerosol optical thickness (result of atcorr process)
    :param analysis_results_wvp:
        file with water vapor content (result of atcorr process)
    :param l1c_band_idx:
        L1C RUT band index
    :param l2a_band_idx:
        L2A atcorr uncertainty band index
    """

    # get files
    l1c_scenarios = glob.glob(analysis_results_l1c.as_posix())
    l2a_scenarios = glob.glob(analysis_results_l2a.as_posix())
    aot = glob.glob(analysis_results_aot.as_posix())[0]
    wvp = glob.glob(analysis_results_wvp.as_posix())[0]
    scl = glob.glob(analysis_results_scl.as_posix())[0]

    l1c_band_idx = 5
    l2a_band_idx = 5

    # loop over bands
    band_list = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B11', 'B12']
    atmospheric_parameters = ['AOT', 'WVP']
    preclass = ['SCL']
    band_list.extend(atmospheric_parameters)
    band_list.extend(preclass)

    for band in band_list:

        print(f'Working on band {band}')

        if band not in atmospheric_parameters and band not in preclass:
        
            l2a_raster = [x for x in l2a_scenarios if Path(x).name.split('_')[1] == band][0]
            l1c_raster = [x for x in l1c_scenarios if Path(x).name.split('_')[1] == band][0]
    
            # read rasters into arrays
            with rio.open(l1c_raster, 'r') as src:
                l1c_data = src.read(l1c_band_idx)
                meta = src.meta
                bounds = src.bounds
    
            with rio.open(l2a_raster, 'r') as src:
                l2a_data = src.read(l2a_band_idx)
                # since the S2-RUT has only 1 decimal precision, round
                # the L2A data
                l2a_data = np.round(l2a_data, 1)
    
            single_fig, single_axs = plt.subplots(1, 2, figsize=(10,20))
    
            epsg = meta['crs'].to_epsg()
    
            # for colormap: find minimum of minima & maximum of maxima
            minmin_l1c = np.round(np.nanmin([np.min(l1c_data), np.nanmin(l1c_data)]), 0)
            maxmax_l1c = np.round(np.nanmax([np.max(l1c_data), np.nanmax(l1c_data)]), 0)
            minmin_l2a = np.round(np.nanmin([np.min(l2a_data), np.nanmin(l2a_data)]), 0)
            maxmax_l2a = np.round(np.nanmax([np.max(l2a_data), np.nanmax(l2a_data)]), 0)
    
            # cut values higher than 10%, otherwise there is not much to see in the image
            labelpad = 20
            if maxmax_l2a > 10.:
                maxmax_l2a = 10
    
            # map
            im_l1c = single_axs[0].imshow(
                l1c_data, vmin=minmin_l1c, vmax=maxmax_l1c, cmap='bwr', interpolation='none',
                extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
            )
            single_axs[0].title.set_text(f'L1C TOA {band}')
            im_l2a = single_axs[1].imshow(
                l2a_data, vmin=minmin_l2a, vmax=maxmax_l2a, cmap='bwr', interpolation='none',
                extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
            )
            single_axs[1].title.set_text(f'L2A BOA {band}')
    
            # add colormap: add_axes[left, bottom, width, heigth)
            divider = make_axes_locatable(single_axs[0])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            single_fig.colorbar(im_l1c, cax=cax, orientation='vertical')
            divider = make_axes_locatable(single_axs[1])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            cbar_l2a = single_fig.colorbar(im_l2a, cax=cax, orientation='vertical')
            cbar_l2a.set_label('Relative Uncertainty [%] (k=1)', rotation=270, fontsize=16,
                           labelpad=labelpad, y=0.5)
    
            single_axs[0].set_xlabel(f'X [m] (EPSG:{epsg})', fontsize=14)
            single_axs[1].set_xlabel(f'X [m] (EPSG:{epsg})', fontsize=14)
            single_axs[0].set_ylabel(f'Y [m] (EPSG:{epsg})', fontsize=14)
    
            single_axs[0].xaxis.set_ticks(np.arange(bounds.left, bounds.right, 5000))
            single_axs[1].xaxis.set_ticks(np.arange(bounds.left, bounds.right, 5000))
            single_axs[0].yaxis.set_ticks(np.arange(bounds.bottom, bounds.top, 5000))
            single_axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
            single_axs[1].yaxis.set_ticks(np.arange(bounds.bottom, bounds.top, 5000))
            single_axs[1].set_yticklabels([])
            single_axs[0].grid(False)
            single_axs[1].grid(False)

        elif band in atmospheric_parameters:
            
            if band == 'AOT':
                with rio.open(aot, 'r') as src:
                    atm_data = src.read(l1c_band_idx)
                    meta = src.meta
                    bounds = src.bounds
            elif band == 'WVP':
                with rio.open(wvp, 'r') as src:
                    atm_data = src.read(l1c_band_idx)
                    meta = src.meta
                    bounds = src.bounds

            single_fig, single_axs = plt.subplots(1, 1, figsize=(10,10))
            single_axs.grid(False)
            epsg = meta['crs'].to_epsg()
    
            # for colormap: find minimum of minima & maximum of maxima
            minmin = np.round(np.nanmin(atm_data), 0)
            maxmax = np.round(np.nanmax(atm_data),0)
    
            # cut values higher than 10%, otherwise there is not much to see in the L1C image
            labelpad = 20

            # map
            im_l1c = single_axs.imshow(
                atm_data, vmin=minmin, vmax=maxmax, cmap='bwr', interpolation='none',
                extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
            )
            single_axs.title.set_text(f'L2A Atmospheric {band}')
    
            # add colormap: add_axes[left, bottom, width, heigth)
            cbar_ax = single_fig.add_axes([0.92, 0.39, 0.04, 0.21])
            n_ticks = int((maxmax - minmin) + 1)
            v1 = np.linspace(minmin, maxmax, n_ticks, endpoint=True)
            cbar_ticks_text = [int(i) for i in v1]
            cbar_ticks_text[-1] = f'>{cbar_ticks_text[-1]}'
            cbar = single_fig.colorbar(im_l2a, cax=cbar_ax)
            cbar.ax.set_yticklabels(cbar_ticks_text)
            cbar.set_label('Uncertainty [%] (k=1)', rotation=270, fontsize=16,
                           labelpad=labelpad, y=0.45)
    
            single_axs.set_xlabel(f'X [m] (EPSG:{epsg})', fontsize=14)
            single_axs.set_ylabel(f'Y [m] (EPSG:{epsg})', fontsize=14)
    
            single_axs.xaxis.set_ticks(np.arange(bounds.left, bounds.right, 5000))
            single_axs.yaxis.set_ticks(np.arange(bounds.bottom, bounds.top, 5000))
            single_axs.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))

        elif band in preclass:
            
            with rio.open(scl, 'r') as src:
                scl_data = src.read()
                meta = src.meta
                bounds = src.bounds
            
            single_fig, single_axs = plt.subplots(1, 2, figsize=(10,20))
            single_axs[0].grid(False)
            single_axs[1].grid(False)
            epsg = meta['crs'].to_epsg()

            cmap = colors.ListedColormap(
                ['black', 'red', 'grey', 'brown', 'green', 'yellow', 'blue', 'lightgrey',
                 'lightsteelblue', 'lavender', 'cyan', 'magenta']
            )

            vote = single_axs[0].imshow(
                scl_data[0,:,:].astype(int), cmap=cmap,
                extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
            )
            single_axs[0].title.set_text(f'L2A SCL (Majority Vote)')

            cbar_vote = single_fig.colorbar(
                vote, orientation='horizontal', ax=single_axs[0],
                ticks=np.linspace(0,11,12).astype(int), pad=0.06
            )
            cbar_vote.set_label('SCL', fontsize=16)

            conf = single_axs[1].imshow(
                scl_data[1,:,:], vmin=0, vmax=100, cmap='Greens', interpolation='none',
                extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
            )
            single_axs[1].title.set_text(f'L2A SCL (Confidence)')

            cbar_conf = single_fig.colorbar(
                conf, orientation='horizontal', ax=single_axs[1], pad=0.06
            )
            cbar_conf.set_label('% Scenario Members Voting for Class', fontsize=16)

            single_axs[0].set_xlabel(f'X [m] (EPSG:{epsg})', fontsize=14)
            single_axs[1].set_xlabel(f'X [m] (EPSG:{epsg})', fontsize=14)
            single_axs[0].set_ylabel(f'Y [m] (EPSG:{epsg})', fontsize=14)
    
            single_axs[0].xaxis.set_ticks(np.arange(bounds.left, bounds.right, 5000))
            single_axs[1].xaxis.set_ticks(np.arange(bounds.left, bounds.right, 5000))
            single_axs[0].yaxis.set_ticks(np.arange(bounds.bottom, bounds.top, 5000))
            single_axs[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
            single_axs[1].yaxis.set_ticks(np.arange(bounds.bottom, bounds.top, 5000))
            single_axs[1].set_yticklabels([])


        # save figure
        fname = out_dir.joinpath(f'{band}_l1c-l2a_uncertainty-map.png')
        single_fig.savefig(fname, bbox_inches='tight')
        plt.close(single_fig)


def _get_roi_mean(
        raster_file: str,
        band_idx: int,
        geom: Polygon,
        nodata_value: float
    ) -> float:
    """
    Calculates the mean of the pixel in a ROI
    """
    with rio.open(raster_file, 'r') as src:
        out_band, _ = rio.mask.mask(
            src,
            [geom],
            crop=True, 
            all_touched=True # IMPORTANT!
        )
    out_band = out_band[band_idx-1,:,:]
    out_band[out_band == nodata_value] = np.nan
    return np.nanmean(out_band)


def extract_roi_unc(
        analysis_results_l1c: Path,
        analysis_results_l2a: Path,
        analysis_results_aot: Path,
        analysis_results_wvp: Path,
        shapefile_rois: Path,
        id_column: str,
        img_date: date
    ) -> pd.DataFrame:
    """
    Once the study areas has been extracted and the statistics per pixel
    have been compiled (using analyze_scenarios_spatial) the single regions
    of interest (ROIs) can be analyzed. These are polygons consisting of 1
    up to N pixels and can represent, e.g., different land use/ cover classes.
    The function extracts the relative uncertainty estimates (band 5 in the
    resulting tif-files of analyze_scenarios_spatial) and computes a spatial
    average. The results are then stored in a pandas dataframe and written to
    CSV.

    :param analysis_results_l1c:
        path object with wildcard to find all .tif files in L1C level output
        from analyze_scenarios_spatial
    :param analysis_results_l2a:
        path object with wildcard to find all .tif files in L2A level output
        from analyze_scenarios_spatial
    :param analysis_results_aot:
        file with aerosol optical thickness (result of atcorr process)
    :param analysis_results_wvp:
        file with water vapor content (result of atcorr process)
    :param shapefile_rois:
        shapefile with polygons defining the single ROIs (each ROI is a feature
        in the file).
    :param id_column:
        name of the column with the unique identifier of each ROI
    :param date:
        acquisition date of the original S2 image (helps to track uncertainty
        over time)
    :return:
        pandas dataframe with the results
    """

    # get files
    l1c_scenarios = glob.glob(analysis_results_l1c.as_posix())
    l2a_scenarios = glob.glob(analysis_results_l2a.as_posix())
    aot = glob.glob(analysis_results_aot.as_posix())[0]
    wvp = glob.glob(analysis_results_wvp.as_posix())[0]

    l1c_band_idx = 5

    # read ROIs from file into a GeoDataFrame
    gdf = gpd.read_file(shapefile_rois)

    # loop over bands
    band_list = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B11', 'B12']
    atmospheric_parameters = ['WVP', 'AOT']
    atmospheric_dict = {
        'WVP': wvp,
        'AOT': aot
    }
    band_list.extend(atmospheric_parameters)

    res = []
    for _, roi in gdf.iterrows():

        # store results in dict
        band_res_l1c = {}
        band_res_l2a = {}

        band_res_l1c['date'] = img_date
        band_res_l2a['date'] = img_date
        band_res_l1c['processing_level'] = 'L1C'
        band_res_l2a['processing_level'] = 'L2A'
        band_res_l1c['ROI'] = roi[id_column]
        band_res_l2a['ROI'] = roi[id_column]

        for band in band_list:

            # get raster for the band
            if band not in atmospheric_parameters:
                l2a_raster = [x for x in l2a_scenarios if Path(x).name.split('_')[1] == band][0]
                l1c_raster = [x for x in l1c_scenarios if Path(x).name.split('_')[1] == band][0]
    
                # get ROI mean of the band
                band_res_l1c[band] = _get_roi_mean(
                    raster_file=l1c_raster,
                    band_idx=l1c_band_idx,
                    geom=roi['geometry'],
                    nodata_value=0.
                )
    
                band_res_l2a[band] = _get_roi_mean(
                    raster_file=l2a_raster,
                    band_idx=l1c_band_idx,
                    geom=roi['geometry'],
                    nodata_value=0.
                )
            else:
                band_res_l2a[band] = _get_roi_mean(
                    raster_file=atmospheric_dict[band],
                    band_idx=l1c_band_idx,
                    geom=roi['geometry'],
                    nodata_value=0.
                )

        res.append(band_res_l1c)
        res.append(band_res_l2a)

    return pd.DataFrame(res)


def main(
        unc_scenario_dir_home: Path,
        out_dir_home: Path,
        in_file_shp: Path,
        in_file_shp_rois: Path,
        id_column: str
    ) -> None:
    """
    main executable function calling all other functions
    required to analyze the results of the uncertainty
    propagation study
    """
    # loop over single scenes and extract the uncertainty of the study area
    # and the regions of interest (ROIs)
    dir_list = next(os.walk(unc_scenario_dir_home.as_posix()))[1]
    for scene in dir_list:

        unc_scenario_dir = unc_scenario_dir_home.joinpath(scene)
        out_dir = out_dir_home.joinpath(scene)
        if not out_dir.exists():
            out_dir.mkdir()

        print(f'Analyzing Uncertainty of {unc_scenario_dir.name}')

        # processing levels of the data; we analyze L1C and L2A
        processing_levels = ['L1C', 'L2A']
    
        #    STEP_1      ANALYZE THE SCENARIOS BY READING ALL REALIZATIONS
        #                FOR YOUR STUDY AREA
        #
        #    FOR THE SPECTRAL BANDS, THE WATER VAPOR AND AEROSOL BAND THIS
        #    RESULTS IN A NEW RASTER FILE WITH 5 BANDS CONTAINING THE MIN,
        #    MAX, MEAN, STD AND STANDARD UNCERTAINTY DENOTING THE SPREAD
        #    AMONG THE SCENARIO MEMBERS

        for processing_level in processing_levels:
        
            analyze_scenarios_spatial(
                unc_scenario_dir=unc_scenario_dir,
                in_file_shp=in_file_shp,
                out_dir=out_dir,
                processing_level=processing_level
            )
    
        #    STEP_2        ANALYSIS AND VISUALIZATION OF UNCERTAINTY
    
        # uncertainty maps
        analysis_results_l1c = out_dir.joinpath('L1C_B*.tif')
        analysis_results_l2a = out_dir.joinpath('L2A_B*.tif')
        analysis_results_wvp = out_dir.joinpath('L2A_WVP_60m_*.tif')
        analysis_results_aot = out_dir.joinpath('L2A_AOT_20m_*.tif')
        analysis_results_scl = out_dir.joinpath('L2A_SCL*.tif')
    
        out_dir_maps = out_dir.joinpath('maps')
        if not out_dir_maps.exists():
            out_dir_maps.mkdir()
    
        unc_maps(
            analysis_results_l1c=analysis_results_l1c,
            analysis_results_l2a=analysis_results_l2a,
            analysis_results_aot=analysis_results_aot,
            analysis_results_wvp=analysis_results_wvp,
            analysis_results_scl=analysis_results_scl,
            out_dir=out_dir_maps
        )
    
        # acquisition date of the S2 image
        img_date = datetime.datetime.strptime(
            unc_scenario_dir.name.split('_')[2].split('T')[0],
            '%Y%m%d'
        ).date()
    
        # extract uncertainty for the selected ROIs
        unc_roi_df = extract_roi_unc(
            analysis_results_l1c=analysis_results_l1c,
            analysis_results_l2a=analysis_results_l2a,
            analysis_results_aot=analysis_results_aot,
            analysis_results_wvp=analysis_results_wvp,
            shapefile_rois=in_file_shp_rois,
            id_column=id_column,
            img_date=img_date
        )
    
        # save (backup) dataframe as csv
        csv_dir = out_dir.joinpath('csv')
        if not csv_dir.exists():
            csv_dir.mkdir()
        fname_csv = csv_dir.joinpath(
            f'spectral-band_l1c-l2a_uncertainty_{in_file_shp_rois.name.split(".")[0]}.csv'
        )
        unc_roi_df.to_csv(fname_csv, index=False)


if __name__ == '__main__':

    ### user inputs
    
    # shapefile (or other vector format) defining the extent of the study area
    in_file_shp = Path(
        './../shp/AOI_Esch_EPSG32632.shp'
    )
    in_file_shp_rois = Path(
        './../shp/ZH_Polygons_2019_EPSG32632_selected-crops.shp'
    )
    id_column = 'crop_type'

    # directory containing the raster realizations
    unc_scenario_dir_home = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL1C_RUT-Scenarios/batch_1'
    )
    # directory where to save the resulting files to
    out_dir_home = Path(
        './../S2A_MSIL2A_Analysis'
    )

    main(
        unc_scenario_dir_home=unc_scenario_dir_home,
        out_dir_home=out_dir_home,
        in_file_shp=in_file_shp,
        in_file_shp_rois=in_file_shp_rois,
        id_column=id_column
    )
