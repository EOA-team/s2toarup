'''
Created on Dec 6, 2021

@author: Lukas Graf
'''

import cv2
import glob
import pandas as pd
import numpy as np

from pathlib import Path
from datetime import datetime
from uncertainties import unumpy
from typing import Tuple

from agrisatpy.io import Sat_Data_Reader
from agrisatpy.io.sentinel2 import S2_Band_Reader
from agrisatpy.utils.constants.sentinel2 import ProcessingLevels


def get_data_and_uncertainty_files(
        orig_dataset_dir: Path,
        uncertainty_analysis_dir: Path,
        vi_name: str
    ) -> pd.DataFrame:
    """
    Reads the bandstacked geoTiff (S2 bands) from the original
    Sentinel-2 data and searches for the corresponding relative uncertainty
    information. Stores the files by acquisition date in a pandas dataframe

    :param orig_dataset_dir:
        directory with the original Sentinel-2 data (*.SAFE datasets)
    :param uncertainty_analysis_dir:
        directory where the corresponding relative uncertainty for the index
        or parameter is stored.
    :param vi_name:
        name of the vegetation index or parameter. It is used to find the
        correct files
    :return:
        pandas dataframe with the file location of the bandstacks
        and file with corresponding uncertainty for the selected vegetation
        index/ parameter alongside with the (acquisition) date
    """

    # define search expressions
    # vegetation parameter/ index
    search_expr_bandstack = orig_dataset_dir.joinpath(
        'S2*_MSIL2A_*.SAFE'
    )
    search_expr_unc = uncertainty_analysis_dir.joinpath(
        f'S2*/L3_{vi_name.upper()}_*.tif'
    )

    # get list of files
    bandstack_list = glob.glob(search_expr_bandstack.as_posix())
    unc_file_list = glob.glob(search_expr_unc.as_posix())

    # convert lists to dataframe
    orig_file_df = pd.DataFrame(bandstack_list, columns=['filename_orig'])
    unc_file_df = pd.DataFrame(unc_file_list, columns=['filename_unc'])

    # extract dates (files are matched and ordered by date)
    orig_file_df['date'] = orig_file_df.filename_orig.apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[2][0:8], '%Y%m%d').date()
    )
    unc_file_df['date'] = unc_file_df.filename_unc.apply(
        lambda x: datetime.strptime(Path(x).parent.name.split('_')[2][0:8], '%Y%m%d').date()
    )

    # join dataframes on date
    df = pd.merge(orig_file_df, unc_file_df, on='date')
    # and sort by date
    df.sort_values(by='date', inplace=True)

    return df


def read_data_and_uncertainty(
        data_df: pd.DataFrame,
        in_file_aoi: Path,
        vi_name: str,
        out_dir_plots: Path
    ) -> Tuple[pd.DataFrame, np.array]:
    """
    This function reads the selected vegetation index (computed on the fly)
    and the associated standard uncertainty in a uarray (from uncertainties)
    and stores the read data as new columns in data_df. In addition, plots
    of the RGB, False-Color Infrared and SCL and the vegetation index are generated.

    :param data_df:
        dataframe returned from ``get_data_and_uncertainty``
    :param in_file_aoi:
        shapefile defining the study area
    :param vi_name:
        name of the vegetation index to compute (e.g., NDVI)
    :param out_dir_plots:
        directory where to save the resulting plots to
        (are stored per image acquisition date)
    :return:
        input dataframe with read data + standard uncertainty
    """

    # loop over datasets (single acquisition dates) and read the data
    # (vegetation index/ parameter + standard uncertainty)
    update_list = []
    ts_stack_list = []

    for _, record in data_df.iterrows():

        # read S2 bandstack data (including scene classification layer, SCL)
        s2_stack = S2_Band_Reader()
        s2_stack.read_from_safe(
            in_dir=Path(record.filename_orig),
            processing_level=ProcessingLevels.L2A,
            in_file_aoi=in_file_aoi
        )

        # mask clouds and cloud shadows (if any) and store the cloudy pixel percentage
        cloud_coverage = s2_stack.get_cloudy_pixel_percentage()

        # plot quicklooks and save them
        plot_dir = out_dir_plots.joinpath(record.date.strftime('%Y%m%d'))
        if not plot_dir.exists():
            plot_dir.mkdir()
        fname_rgb = plot_dir.joinpath('rgb_quicklook.png')
        fname_nir = plot_dir.joinpath('falsecolor-nir_quicklook.png')
        fname_scl = plot_dir.joinpath('scene-classification_quicklook.png')
        fname_vi = plot_dir.joinpath(f'{vi_name.lower()}-nominal_quicklook.png')
        fname_unc = plot_dir.joinpath(f'{vi_name.lower()}-stdunc_quicklook.png')

        fig_rgb = s2_stack.plot_rgb()
        fig_rgb.savefig(fname=fname_rgb, bbox_inches='tight')
        fig_nir = s2_stack.plot_false_color_infrared()
        fig_nir.savefig(fname=fname_nir, bbox_inches='tight')
        fig_scl = s2_stack.plot_scl()
        fig_scl.savefig(fname=fname_scl, bbox_inches='tight')

        # resample spectral bands to 10m (i.e., those bands with 20m resolution)
        s2_stack.resample(
            target_resolution=10,
            resampling_method=cv2.INTER_CUBIC,
            bands_to_exclude=['scl']
        )

        # resample SCL band to 10m
        s2_stack.resample(
            target_resolution=10,
            resampling_method=cv2.INTER_NEAREST_EXACT
        )

        # calculate the vegetation index
        s2_stack.calc_vi(vi=vi_name)
        
        fig_vi = s2_stack.plot_band(band_name=vi_name, colormap='summer')
        fig_vi.savefig(fname=fname_vi, bbox_inches='tight')

        s2_stack.mask_clouds_and_shadows(bands_to_mask=[vi_name])

        # read uncertainty data
        uncertainty_band = Sat_Data_Reader()
        uncertainty_band.read_from_bandstack(fname_bandstack=record.filename_unc)

        fig_unc = uncertainty_band.plot_band(
            band_name='abs_stddev'
        )
        fig_unc.savefig(fname=fname_unc, bbox_inches='tight')

        # apply the cloud mask also to the absolute uncertainty band
        uncertainty_band.add_band(
            band_name='scl',
            band_data=s2_stack.data['scl']
        )

        uncertainty_band.mask(
            name_mask_band='scl',
            mask_values=[3, 8, 9, 10],
            bands_to_mask=['abs_stddev'],
            keep_mask_values=False
        )

        # store data and standard uncertainty
        uarray = unumpy.uarray(
            nominal_values=s2_stack.data[vi_name],
            std_devs=uncertainty_band.data['abs_stddev']
        )
        ts_stack_list.append(uarray)

        update_list.append({
            'date': record.date,
            'cloudy_pixel_percentage': cloud_coverage
        })

    # stack arrays
    stack_ts = np.dstack(ts_stack_list)

    # added update columns to input data frame
    update_df = pd.DataFrame(update_list)
    merged = pd.merge(data_df, update_df, on='date')

    return merged, stack_ts


def get_pixel(stack_df):
    pass

if __name__ == '__main__':
    # original Sentinel-2 scenes with vegetation indices
    orig_dataset_dir = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL1C_orig'
    )
    
    # directory with uncertainty analysis results
    uncertainty_analysis_dir = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL2A_Analysis'
    )
    
    # vegetation index to consider
    vi_name = 'NDVI'
    data_df = get_data_and_uncertainty_files(
        orig_dataset_dir=orig_dataset_dir,
        uncertainty_analysis_dir=uncertainty_analysis_dir,
        vi_name=vi_name
    )

    in_file_aoi = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/shp/AOI_Esch_EPSG32632.shp'
    )

    out_dir_plots = uncertainty_analysis_dir.joinpath('scene_quicklooks')
    if not out_dir_plots.exists():
        out_dir_plots.mkdir()

    unc_df = read_data_and_uncertainty(
        data_df=data_df,
        vi_name=vi_name,
        in_file_aoi=in_file_aoi,
        out_dir_plots=out_dir_plots
    )

      