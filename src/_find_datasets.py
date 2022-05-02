'''
Helper function to find the original vegetation index/parameter data and its
uncertainty information for all scenes available.
'''

import glob
import geopandas as gpd
import pandas as pd

from datetime import datetime
from pathlib import Path
from typing import List, Tuple

from agrisatpy.core.band import Band
from agrisatpy.core.sensors import Sentinel2

from logger import get_logger

logger = get_logger('read_data')


def get_data_and_uncertainty_files(
        vi_dir: Path,
        uncertainty_analysis_dir: Path,
        vi_name: str
    ) -> pd.DataFrame:
    """
    Searches the vegetation parameter/index generated from the original
    Sentinel-2 data alongside with the corresponding relative uncertainty
    information. Stores the files by acquisition date in a pandas dataframe.

    :param vi_dir:
        directory containing the vegetation indices/ parameters derived from the
        original (i.e., non-scenario-based) datasets and the resampled original
        Sentinel-2 data in L2A processing level
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
    search_expression_bandstack = vi_dir.joinpath(
        '*_None_10m.jp2'
    )
    search_expr_veg_par = vi_dir.joinpath(
        Path('Vegetation_Indices').joinpath(f'VI_*None_10m_{vi_name.upper()}.tif')
    )
    if vi_name == 'GLAI':
        vi_name = 'LAI'
    search_expr_unc = uncertainty_analysis_dir.joinpath(
        f'S2*/L3_{vi_name.upper()}_*bounds.tif'
    )

    # get list of files
    bandstack_list = glob.glob(search_expression_bandstack.as_posix())
    veg_par_list = glob.glob(search_expr_veg_par.as_posix())
    unc_file_list = glob.glob(search_expr_unc.as_posix())

    # convert lists to dataframe
    bandstack_file_df = pd.DataFrame(bandstack_list, columns=['filename_bandstack'])
    veg_par_file_df = pd.DataFrame(veg_par_list, columns=['filename_veg_par'])
    unc_file_df = pd.DataFrame(unc_file_list, columns=['filename_unc'])

    # extract dates (files are matched and ordered by date)
    veg_par_file_df['date'] = veg_par_file_df.filename_veg_par.apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[1], '%Y%m%d').date()
    )
    unc_file_df['date'] = unc_file_df.filename_unc.apply(
        lambda x: datetime.strptime(Path(x).parent.name.split('_')[2][0:8], '%Y%m%d').date()
    )
    bandstack_file_df['date'] = bandstack_file_df.filename_bandstack.apply(
        lambda x: datetime.strptime(Path(x).name.split('_')[0], '%Y%m%d').date()
    )

    # join dataframes on date
    tmp_df = pd.merge(veg_par_file_df, unc_file_df, on='date')
    df = pd.merge(bandstack_file_df, tmp_df, on='date')
    # and sort by date
    df.sort_values(by='date', inplace=True)

    return df

def read_data_and_uncertainty(
        data_df: pd.DataFrame,
        parcels: Path,
        vi_name: str,
        crop_periods: Path
    ) -> Tuple[pd.DataFrame, List[Sentinel2]]:
    """
    This function reads the selected vegetation index (computed in step 7)
    and the associated standard (absolute) uncertainty and stores the read
    data as new columns in data_df.

    :param data_df:
        dataframe returned from ``get_data_and_uncertainty``
    :param parcels:
        shapefile defining the crop parcels
    :param vi_name:
        name of the vegetation index or parameter to process (e.g., NDVI)
    :param crop_periods:
        CSV file specifying temporal range of crop growth to avoid errors in the LSP
        metrics due to inter-crops
    :return:
        input dataframe with read data + list of read vegetation data and
        their absolute uncertainty per image acquisition date
    """

    # loop over datasets (single acquisition dates) and read the data
    # (vegetation index/ parameter + standard uncertainty)
    update_list = []
    ts_stack_dict = {}

    # read crop periods into data frame
    crop_periods = pd.read_csv(crop_periods)

    rdx = 1
    for _, record in data_df.iterrows():

        # read field parcels into dataframe and delete those parcels with crops that
        # are not in their main growing period
        parcels_df = gpd.read_file(parcels)
        # drop empty geometries
        drop_idx = parcels_df[parcels_df.geometry == None].index
        parcels_df.drop(index=drop_idx, inplace=True)

        # read VI data
        vi_file = Path(record.filename_veg_par)
        collection = Sentinel2()
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=vi_file,
            band_idx=1,
            band_name_dst=vi_name,
            vector_features=parcels_df
        )
        # read vegetation parameter/index absolute uncertainty band
        vi_unc_file = Path(record.filename_unc)
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=vi_unc_file,
            band_idx=4,
            vector_features=parcels_df,
            band_name_dst=f'{vi_name}_unc'
        )

        # read S2 blue and red bands plus SCL file to mask out clouds and shadows
        scl_file = next(
            Path(record.filename_bandstack).parent.rglob('scene_classification/*.tiff')
        )
        collection.add_band(
            band_constructor=Band.from_rasterio,
            fpath_raster=scl_file,
            band_idx=1,
            band_name_dst='SCL',
            vector_features=parcels_df
        )
        # add crop codes
        collection.add_band(
            band_constructor=Band.from_vector,
            vector_features=parcels_df,
            geo_info=collection[vi_name].geo_info,
            band_name_src='crop_code',
            band_name_dst='crop_code',
            snap_bounds=collection[vi_name].bounds
        )
        collection.mask(
            mask=collection[vi_name].values.mask,
            bands_to_mask=['crop_code'],
            inplace=True
        )
        # check the current sensing date and mask crops not in their growing season
        sensing_date = record.date
        for _,crop in crop_periods.iterrows():
            start = datetime.strptime(crop.Start, '%Y-%m-%d').date()
            end = datetime.strptime(crop.End, '%Y-%m-%d').date()
            if not start <= sensing_date <= end:
                crop_name = crop.Crop
                if crop_name == 'Grain Maize':
                    crop_name = 'Corn'
                if crop_name == 'Rapeseed':
                    crop_name = 'Canola'
                crop_code = parcels_df[parcels_df.crop_type == crop_name]['crop_code'].iloc[0]
                collection.mask(
                    mask='crop_code',
                    mask_values=[crop_code],
                    bands_to_mask=collection.band_names,
                    inplace=True
                )
        # mask clouds, shadows and snow
        collection.mask_clouds_and_shadows(
            bands_to_mask=collection.band_names,
            cloud_classes=[1,2,3,6,7,8,9,10,11],
            inplace=True
        )
        cloud_coverage = collection.get_cloudy_pixel_percentage(cloud_classes=[1,2,3,6,7,8,9,10,11])

        # save masked files to disk (thus, the reference run will also have the correct
        # input)
        vi_masked = record.filename_veg_par
        vi_masked = Path(vi_masked.replace('.tif', '_scl-filtered.tif'))
        collection.to_rasterio(
            fpath_raster=vi_masked,
            band_selection=[vi_name]
        )
        vi_unc_masked = record.filename_unc
        vi_unc_masked = Path(vi_unc_masked.replace('.tif', '_scl-filtered.tif'))
        collection.to_rasterio(
            fpath_raster=vi_unc_masked,
            band_selection=[f'{vi_name}_unc']
        )
        ts_stack_dict[record.date] = collection

        update_list.append({
            'date': record.date,
            'cloudy_pixel_percentage': cloud_coverage,
            'filename_veg_par_scl': vi_masked,
            'filename_unc_scl': vi_unc_masked
        })
        logger.info(f'Read record {rdx}/{data_df.shape[0]}')
        rdx += 1

    # added update columns to input data frame
    update_df = pd.DataFrame(update_list)
    merged = pd.merge(data_df, update_df, on='date')
    return merged, ts_stack_dict
