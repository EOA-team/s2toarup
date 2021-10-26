
import glob
import math
import pandas as pd
from pathlib import Path
import geopandas as gpd
import rasterio as rio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker

from agrisatpy.processing.extraction.utils import raster2table
from agrisatpy.processing.extraction.sentinel2 import S2singlebands2table
from copy import deepcopy


plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

plt.rcParams['axes.titlesize'] = 15
gain_factor_refl = 0.01

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



def analyze_scenarios_spatial(
        unc_scenario_dir: Path,
        in_file_shp: Path,
        out_dir: Path,
        processing_level: str,
        **kwargs
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

            print(f'Extracting {band} from {spatial_res}m spatial resolution')

            scenario_files = glob.glob(
                str(unc_scenario_dir.joinpath(search_expr + band_dict[band]))
            )

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
            meta.update(
                {
                    "count": 5,
                    "dtype": np.float64
                 }
            )

            fname = f'{processing_level}_{band}_{spatial_res}m_{in_file_shp.name.split(".")[0]}.tif'
            
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
                data_arr_norm = data_arr / np.linalg.norm(data_arr)
                dst.set_band_description(5, 'rel_std_unc')
                rel_std = np.nanstd(data_arr, axis=0) / np.nanmean(data_arr, axis=0)
                dst.write_band(5, rel_std * 100)

        

def extract_scenarios_roi(
        unc_scenario_dir: Path,
        in_file_shp: Path,
        out_dir: Path,
        processing_level: str,
        **kwargs
    ) -> None:
    """
    Function to extract spread among scenario members for selected
    polygons defining regions of interest. Will analyze L1C and
    L2A data to quantify the impact of both processing levels and
    their uncertainty on the reflectance quanities and derived
    products such as the aerosol optical depth (AOT), the water vapor
    (WVP) and the scene classification layer (SCL). This function
    works very well on small for small ROIs. For larger ROIs (10 sqkm and
    more), using a purely raster-based method is faster!

    :param unc_scenario_dir:
        directory with Sen2Cor runs (L2A) or radiometric uncertainty
        scenarios (L1C)
    :param in_file_shp:
        shape file with ROIs
    :param out_dir:
        directory where to store dataframes with extracted data
        (also for later analysis)
    :param processing_level:
        either 'L1C' or 'L2A' for identifying the S2 processing level
        (before and after atmospheric correction)
    """

    # check processing level
    if processing_level == 'L1C':
        band_dict_s2 = band_dict_l1c
        abbrev = 'MSIL1C'
    elif processing_level == 'L2A':
        band_dict_s2 = band_dict_l2a
        abbrev = 'MSIL2A'

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

            print(f'Extracting {band} from {spatial_res}m spatial resolution')

            scenario_files = glob.glob(
                str(unc_scenario_dir.joinpath(search_expr + band_dict[band]))
            )
            res = []
            for idx, scenario_file in enumerate(scenario_files):
                df = raster2table(
                    in_file=Path(scenario_file),
                    in_file_polys=Path(in_file_shp),
                    out_colnames=[band],
                    **kwargs
                )
                df['scenario'] = idx + 1
                res.append(df)
            band_df = pd.concat(res)
            fname_band_df = f'{processing_level}_{band}_{spatial_res}m_{in_file_shp.name.split(".")[0]}.csv'
            band_df.to_csv(out_dir.joinpath(fname_band_df), index=False)


def unc_boxplots(
        l1c_roi_scenarios: str,
        l2a_roi_scenarios: str,
        l1c_original_data: Path,
        l2a_original_data: Path,
        id_column: str,
        land_use_info: Path,
        out_dir: Path
    ) -> None:
    """
    Takes extracted ROI data (CSV files) and makes a box plot for each
    numerical variable (spectral bands, AOT and WV) taking into account
    the L1C uncertainty, the L2A uncertainty and the original values.

    :param l1c_roi_scenario:
        basename of csv files with extracted L1C scenarios for the ROIs
    :param l2a_roi_scenario:
        basename of csv files with extracted L2A scenarios for the ROIs
    :param l1_original_data:
        csv file with original L1C data for the ROIs
    :param l2a_original_data:
        csv file with original L1C data for the ROIs
    :param id_column:
        name of the column identifying each ROI
    :param land_use_info:
        file specifying the land use information for each ROI
    """

    # get files
    l1c_scenarios = glob.glob(l1c_roi_scenarios.as_posix())
    l2a_scenarios = glob.glob(l2a_roi_scenarios.as_posix())
    l1c_original = glob.glob(l1c_original_data.as_posix())
    l2a_original = glob.glob(l2a_original_data.as_posix())
    luc_info = gpd.read_file(land_use_info)

    # save uncertainties per land cover class
    uncertainties = []

    # loop over the different spatial resolutions and over the spectral
    # bands
    for resolution in band_dict_l1c:

        # find corresponding original datasets (contain multiple bands)
        l1c_orig_res = [x for x in l1c_original if x.split('_')[-1].startswith(resolution)]
        l1c_orig_df = pd.read_csv(l1c_orig_res[0])
        l2a_orig_res = [x for x in l2a_original if x.split('_')[-1].startswith(resolution)]
        l2a_orig_df = pd.read_csv(l2a_orig_res[0])

        # number of ROIs
        rois = l1c_orig_df[id_column].unique()
        roi_luc_mapping = dict.fromkeys(rois)
        for roi in roi_luc_mapping:
            roi_luc_mapping[roi] = luc_info[luc_info[id_column] == roi]['crop_type'].iloc[0]

        # join with landuse information
        l1c_orig_df = pd.merge(l1c_orig_df, luc_info, on=id_column)
        l2a_orig_df = pd.merge(l2a_orig_df, luc_info, on=id_column)

        n_rois = len(rois)

        # group original data by ROIs
        l1c_orig_grouped = l1c_orig_df.groupby(id_column).mean()
        l2a_orig_grouped = l2a_orig_df.groupby(id_column).mean()

        # numnber of spectral bands with the current spatial resolution
        band_dict = band_dict_l1c[resolution]
        n_bands = len(band_dict)

        # loop over single bands here
        fig, axs = plt.subplots(n_bands, n_rois, figsize=(35,40))

        for idx, band in enumerate(band_dict):

            inner_band_dict = {'band': band}

            # read scenarios of the current band
            l1c_scenario_band = pd.read_csv(
                [x for x in l1c_scenarios if Path(x).name.split('_')[1] == band][0]
            )
            l2a_scenario_band = pd.read_csv(
                [x for x in l2a_scenarios if Path(x).name.split('_')[1] == band][0]
            )

            # join land use info
            l1c_scenario_band = pd.merge(l1c_scenario_band, luc_info, on='fid')
            l2a_scenario_band = pd.merge(l2a_scenario_band, luc_info, on='fid')

            # group by scenario and ROI to get the ROI mean per scenario
            l1c_grouped = l1c_scenario_band.groupby(by=['scenario', id_column]).mean()
            l2a_grouped = l2a_scenario_band.groupby(by=['scenario', id_column]).mean()

            # boxplot
            for jdx, roi in enumerate(rois):

                l1c_grouped_roi = l1c_grouped.query(f'{id_column} == {roi}')
                l2a_grouped_roi = l2a_grouped.query(f'{id_column} == {roi}')

                # calculate standard uncertainty as variation coefficient
                l1c_roi_unc = (np.std(l1c_grouped_roi[band]) / np.mean(l1c_grouped_roi[band])) * 1000
                l2a_roi_unc = (np.std(l2a_grouped_roi[band]) / np.mean(l2a_grouped_roi[band])) * 1000

                inner_band_dict[f'{roi_luc_mapping[roi]}_roi_L1C'] = l1c_roi_unc
                inner_band_dict[f'{roi_luc_mapping[roi]}_roi_L2A'] = l2a_roi_unc

                axs[idx][jdx].boxplot(
                    [
                        l1c_grouped_roi[band]*gain_factor_refl,
                        l2a_grouped_roi[band]*gain_factor_refl
                    ]
                )
                axs[idx][jdx].set_xticklabels(['L1C', 'L2A'])
                [t.set_color(i) for (i,t) in zip(['r','b'], axs[idx][jdx].xaxis.get_ticklabels())]

                axs[idx][jdx].axhline(
                    l1c_orig_grouped[l1c_orig_grouped.index==roi][band].values * gain_factor_refl,
                    color='r',
                    label='L1C original'
                )
                axs[idx][jdx].axhline(
                    l2a_orig_grouped[l2a_orig_grouped.index==roi][band].values * gain_factor_refl,
                    color='b',
                    label='L2A original'
                )
                axs[idx][jdx].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

                if idx == 0:
                    axs[idx][jdx].title.set_text(roi_luc_mapping[roi])
                if jdx == 0:
                    axs[idx][jdx].set_ylabel(f'Reflectance {band} (%)', fontsize=28)

            uncertainties.append(inner_band_dict)

        # save plot
        fname = f'{resolution}m_bands_l1c-l2a_refl-boxplots.png'
        fig.savefig(out_dir.joinpath(fname), bbox_inches='tight')
        plt.close(fig)
    
    # save the uncertainties per ROI into a dataframe and then to CSV
    unc_df = pd.DataFrame(uncertainties)
    unc_df.sort_values(by='band', ascending=True, inplace=True)
    unc_df.to_csv(out_dir.joinpath('spectral_bands_l1c-l2a_unc.csv'), index=False)


def unc_maps(
        l1c_scenarios: str,
        l2a_scenarios: str,
        l1c_band_idx: int,
        l2a_band_idx: int,
        out_dir: Path
    ) -> None:
    """
    Maps raster values of L1C and L2A uncertainty values to reveal
    spatial pattern of uncertainty and their land cover/ use dependency

    :param l1c_scenarios:
        basename of tif files with L1C RUT output
    :param l2a_scenarios:
        basename of the atmospheric correction uncertainty output
    :param l1c_band_idx:
        L1C RUT band index
    :param l2a_band_idx:
        L2A atcorr uncertainty band index
    """

    # get files
    l1c_scenarios = glob.glob(l1c_scenarios.as_posix())
    l2a_scenarios = glob.glob(l2a_scenarios.as_posix())

    # loop over bands
    band_list = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B11', 'B12']

    for band in band_list:

        print(f'Working on band {band}')
        
        l2a_raster = [x for x in l2a_scenarios if Path(x).name.split('_')[1] == band][0]
        l1c_raster = [
            x for x in l1c_scenarios if Path(x).name.split('_')[-1].split('.')[0] == band.lower()
        ][0]

        # read rasters into arrays
        with rio.open(l1c_raster, 'r') as src:
            l1c_data = src.read(l1c_band_idx) * 0.1
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
        minmin = math.trunc(np.min([np.min(l1c_data), np.min(l2a_data)]))
        maxmax = np.round(np.max([np.max(l1c_data), np.max(l2a_data)]),0)

        # cut values higher than 10%, otherwise there is not much to see in the L1C image
        labelpad = 20
        if maxmax > 10.:
            maxmax = 10

        # map
        im_l1c = single_axs[0].imshow(
            l1c_data, vmin=minmin, vmax=maxmax, cmap='bwr', interpolation='none',
            extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
        )
        single_axs[0].title.set_text(f'L1C TOA {band}')
        im_l2a = single_axs[1].imshow(
            l2a_data, vmin=minmin, vmax=maxmax, cmap='bwr', interpolation='none',
            extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
        )
        single_axs[1].title.set_text(f'L2A TOC {band}')

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

        del(cbar)
        del(single_fig)


def analyze_scl(
        scl_scenarios: Path,
        out_dir: Path
    ):
    """
    Analyzes the differences in the scene classification layer among scenario
    members to assses the impact of uncertainty on the scene classification
    outcomes and hence the selection of "valid" pixels
    """
    



if __name__ == '__main__':

    # uncertainty maps
    unc_res_raster_dir = Path('/mnt/ides/Lukas/owncloud/PhenomEn_Project/PAPER_1/S2A_MSIL2A_Analysis/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621')
    l1c_scenarios = unc_res_raster_dir.joinpath('*_rut_b*.tif')
    l2a_scenarios = unc_res_raster_dir.joinpath('L2A_B*.tif')
    l1c_band_idx = 1
    l2a_band_idx = 5
    out_dir = Path('/mnt/ides/Lukas/04_Work/Uncertainty/Atcorr')

    unc_maps(l1c_scenarios, l2a_scenarios, l1c_band_idx, l2a_band_idx, out_dir)

    # box plots of ROIs (different landuses)
    unc_res_roi_dir = Path('/mnt/ides/Lukas/owncloud/PhenomEn_Project/PAPER_1/S2A_MSIL2A_Analysis/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621/ROIs')
    l1c_roi_scenarios = unc_res_roi_dir.joinpath('L1C_*.csv')
    l2a_roi_scenarios = unc_res_roi_dir.joinpath('L2A_*.csv')
    l1c_original_data = unc_res_roi_dir.joinpath('*_MSIL1C_*.csv')
    l2a_original_data = unc_res_roi_dir.joinpath('*_MSIL2A_*.csv')
    land_use_info = Path('/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/STUDY_AREA/ZH_Polygons_2019_EPSG32632_samples_scenarios.shp')
    id_column = 'fid'

    unc_boxplots(
        l1c_roi_scenarios=l1c_roi_scenarios,
        l2a_roi_scenarios=l2a_roi_scenarios,
        l1c_original_data=l1c_original_data,
        l2a_original_data=l2a_original_data,
        id_column=id_column,
        land_use_info=land_use_info,
        out_dir=out_dir
    )

    # define inputs
    original_scene_l1c = Path(
        '/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL1C_orig/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621.SAFE'
    )
    original_scene_l2a = Path(
        '/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL1C_orig/S2B_MSIL2A_20190830T102029_N9999_R065_T32TMT_20211020T160133.SAFE'
    )
    unc_scenario_dir = Path(
        '/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL1C_RUT-Scenarios/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621'
    )
    in_file_shp = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/STUDY_AREA/ZH_Polygons_2019_EPSG32632_samples_scenarios.shp'
    )
    options = {
        'buffer': 0.,
        'id_column': 'fid'
    }
    processing_levels = ['L2A']
    # define directory where to backup extracted pixel values
    out_dir = Path(
        '/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL2A_Analysis/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621/ROIs'
    )

    # single ROIs (different land use types)
    for processing_level in processing_levels:

        extract_scenarios_roi(
            unc_scenario_dir=unc_scenario_dir,
            in_file_shp=in_file_shp,
            out_dir=out_dir,
            processing_level=processing_level,
            **options
        )
        # extract original data (i.e., without uncertainty scenarios)
        if processing_level == 'L1C':
            in_file = original_scene_l1c
            is_l2a = False
        elif processing_level == 'L2A':
            in_file = original_scene_l2a
            is_l2a = True
        resolutions = [10, 20, 60]
        for resolution in resolutions:
            df, scl = S2singlebands2table(
                in_dir=in_file,
                in_file_polys=in_file_shp,
                filter_clouds=False,
                is_L2A=is_l2a,
                product_date='20190830',
                resolution=resolution,
                **options
            )
            fname = f'{in_file.name.split(".")[0]}_{resolution}m.csv'
            df.to_csv(out_dir.joinpath(fname), index=False)
            if is_l2a and resolution == 20:
                fname_scl = f'{in_file.name.split(".")[0]}_scl.csv'
                scl.to_csv(out_dir.joinpath(fname_scl), index=False)
        
    #################################3

    # entire region
    in_file_shp = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/STUDY_AREA/AOI_Esch_EPSG32632-large.shp'
    )
    analyze_scenarios_spatial(
        unc_scenario_dir=unc_scenario_dir,
        in_file_shp=in_file_shp,
        out_dir=out_dir,
        processing_level=processing_level
    )    