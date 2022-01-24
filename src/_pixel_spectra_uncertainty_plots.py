'''
Plots pixel spectra at user-defined locations and their radiometric uncertainty
(Level 1C and L2A) plus extracts the resulting NDVI, EVI and green LAI values
(Level 3).
'''

import glob
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from agrisatpy.io.sentinel2 import Sentinel2Handler
from agrisatpy.utils.constants.sentinel2 import central_wavelengths
from agrisatpy.utils.constants.sentinel2 import ProcessingLevels
from agrisatpy.utils.constants.sentinel2 import s2_gain_factor
from copy import deepcopy
from pathlib import Path

from logger import get_logger

logger = get_logger('pixel_spectra_unc')


spectral_bands = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12']
veg_params = ['NDVI', 'EVI', 'LAI']
# veg_params = ['NDVI', 'EVI']
processing_levels = ProcessingLevels
s2_gain_factor *= 100

plt.style.use('ggplot')
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)


def plot_pixel_spectra(
        point_features: Path,
        uncertainty_analysis_dir: Path,
        out_dir: Path
    ):
    """
    Plots pixel spectra and their uncertainty alongside the resulting VIs
    and green LAI values and their uncertainty.

    :param point_features:
        ESRI shapefile with point-like locations for selecting pixel spectra
        and the resulting vegetation parameter/indices
    :param uncertainty_analysis_dir:
        directory where the results of the uncertainty analysis are stored
        (L1C, L2A and L3 processing level)
    :param out_dir:
        directory where to save the resulting plots to
    """

    # search available scenes
    search_expr = 'S2*_MSIL1C_*'
    scenes = glob.glob(uncertainty_analysis_dir.joinpath(search_expr).as_posix())

    for idx, scene in enumerate(scenes):

        logger.info(f'Working on {scene} ({idx+1}/{len(scenes)})')

        # save plots to sub-directory in out_dir
        out_dir_scene = out_dir.joinpath(Path(scene).name)
        out_dir_scene.mkdir(exist_ok=True)

        # determine sensor (S2A or S2B) to get central wavelengths of the spectral
        # bands selected (10 and 20m bands)
        sensor = Path(scene).name[0:3]
        wvl = central_wavelengths[sensor]
        wvl = [v for k,v in wvl.items() if k in spectral_bands]
        unit = central_wavelengths['unit']

        # get acquisition date of the scene
        acqui_date = Path(scene).name.split('_')[2][0:8]
        year, month, day = acqui_date[0:4], acqui_date[4:6], acqui_date[6:8]

        # using the point features we can extract the pixel values at each location
        # we need the mean reflectance in the band plus the absolute uncertainty
        gdf_band = None
        for spectral_band in spectral_bands:

            # get L1C and L2A reflectance
            for processing_level in processing_levels:

                proc_level = processing_level.name
                fpath_spectral_band = glob.glob(
                    Path(scene).joinpath(
                        f'{proc_level}_{spectral_band}_*.tif'
                    ).as_posix()
                )[0]

                if gdf_band is None:
                    gdf_band = Sentinel2Handler.read_pixels(
                        point_features=point_features,
                        raster=fpath_spectral_band,
                        band_selection=['mean', 'abs_stddev']
                    )
                else:
                    gdf_band = Sentinel2Handler.read_pixels(
                        point_features=gdf_band,
                        raster=fpath_spectral_band,
                        band_selection=['mean', 'abs_stddev']
                    )
                band_name = f'{spectral_band}_{proc_level}'
                gdf_band = gdf_band.rename(
                    columns={
                        'mean': band_name,
                        'abs_stddev': f'{spectral_band}_unc_{proc_level}'
                    }
                )

        # read vegetation index + parameter values and uncertainty
        for veg_par in veg_params:
            fpath = glob.glob(
                Path(scene).joinpath(f'L3_{veg_par}_*.tif').as_posix()
            )[0]
            gdf_band = Sentinel2Handler.read_pixels(
                    point_features=gdf_band,
                    raster=fpath,
                    band_selection=['mean', 'abs_stddev']
                )
            gdf_band = gdf_band.rename(
                columns={
                    'mean': veg_par,
                    'abs_stddev': f'{veg_par}_unc'
                }
            )

        # define filter expression for L1C, L2A results for plotting
        l1c_cols = gdf_band.columns.str.contains('B*_L1C')
        l1c_unc_cols = gdf_band.columns.str.contains('B*_unc_L1C')
        l1c_cols[l1c_unc_cols] = False

        l2a_cols = gdf_band.columns.str.contains('B*_L2A')
        l2a_unc_cols = gdf_band.columns.str.contains('B*_unc_L2A')
        l2a_cols[l2a_unc_cols] = False

        l3_unc_cols = [f'{x}_unc' for x in veg_params]

        # plot the spectra with uncertainty as error bars, each record in the DataFrame
        # denotes a single pixel which we can plot
        for _, record in gdf_band.iterrows():

            f = plt.figure(num=1, figsize=(10,7))
            ax = f.add_subplot(111)

            title_str = f'{record.crop_type} {year}-{month}-{day}\n' \
                        f'x={np.round(record.geometry.x)}m, y={np.round(record.geometry.y)}m'

            # plot L1C and L2A spectrum of the pixel (scale reflectance values between 0 and 100%)
            ax.errorbar(
                x=wvl,
                y=record[l1c_cols]*s2_gain_factor,
                yerr=record[l1c_unc_cols]*s2_gain_factor,
                label='Level 1C',
                elinewidth=3,
                capsize=3
            )

            ax.errorbar(
                x=wvl,
                y=record[l2a_cols]*s2_gain_factor,
                yerr=record[l2a_unc_cols]*s2_gain_factor,
                label='Level L2A',
                elinewidth=3,
                capsize=3
            )

            ax.set_xlabel(f'Wavelength [{unit}]', fontsize=24)
            ax.set_ylabel('Reflectance [%]', fontsize=24)
            ax.legend(fontsize=20)
            ax.set_title(title_str, fontdict={'fontsize': 24})

            # add table with veg index/ parameter data
            plt.subplots_adjust(bottom=0.3)
            cols_veg_pars = deepcopy(veg_params)
            cols_veg_pars.extend(l3_unc_cols)
            labels_veg_pars = [x.replace('_',' ').replace('unc', r'$\sigma$') for x in cols_veg_pars]
            vi_values = record[cols_veg_pars].values.reshape(1,len(cols_veg_pars)).astype('float')
            vi_values = vi_values.round(5)
            the_table = ax.table(
                cellText=vi_values,
                loc='bottom',
                colLabels=labels_veg_pars,
                bbox=[0., -0.5, 1, 0.25]
            )
            the_table.auto_set_font_size(False)
            the_table.set_fontsize(16)

            fname = out_dir_scene.joinpath(
                f'{record.crop_type}_x{int(np.round(record.geometry.x))}_y{int(np.round(record.geometry.y))}.png'
            )
            f.savefig(fname, dpi=100, bbox_inches='tight')
            plt.close(f)
    
        logger.info(f'Finished {scene} ({idx+1}/{len(scenes)})')

    logger.info('Done')


if __name__ == '__main__':

    point_features = Path('../shp/ZH_Points_2019_EPSG32632_selected-crops.shp')
    uncertainty_analysis_dir = Path('../S2A_MSIL2A_Analysis/1000_scenarios')
    out_dir = Path('../S2A_MSIL2A_Analysis/pixel_spectra')
    out_dir.mkdir(exist_ok=True)

    plot_pixel_spectra(point_features, uncertainty_analysis_dir, out_dir)
    