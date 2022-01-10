'''
Helper function to find the original vegetation index/parameter data and its
uncertainty information for all scenes available.
'''

import glob
import pandas as pd

from datetime import datetime
from pathlib import Path


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
        '*_pixel_division_10m.tiff'
    )
    search_expr_veg_par = vi_dir.joinpath(
        Path('Vegetation_Indices').joinpath(f'VI_*None_10m_{vi_name.upper()}.tif')
    )
    search_expr_unc = uncertainty_analysis_dir.joinpath(
        f'S2*/L3_{vi_name.upper()}_*.tif'
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
