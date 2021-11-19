# -*- coding: utf-8 -*-

import glob
import rasterio as rio
import numpy as np
from pathlib import Path
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from agrisatpy.processing.resampling.sentinel2 import resample_and_stack_S2

# TODO: Do we want to keep only those indices that do not require resampling??

# map the Sentinel-2 bands to be used for index calculation
s2_band_mapping = {
    'B02': 'blue',
    'B03': 'green',
    'B04': 'red',
    'B05': 'red_edge_1',
    'B06': 'red_edge_2',
    'B07': 'red_edge_3',
    'B08': 'nir'
}


def NDVI(
        red: np.array,
        nir: np.array
    ) -> np.array:
    """
    Calculates the Normalized Difference Vegetation Index
    (NDVI).

    :param red:
        reflectance in the red band (Sentinel-2 B04)
    :param nir:
        reflectance in the near infrared band (Sentinel-2 B08)
        spectrum
    :return:
        NDVI values
    """

    return (nir - red) / (nir + red)


def EVI(
        blue: np.array,
        red: np.array,
        nir: np.array
    ):
    """
    Calculates the Enhanced Vegetation Index (EVI) following the formula
    provided by Huete et al. (2002)

    :param blue:
        reflectance in the blue band (Sentinel-2 B02)
    :param red:
        reflectance in the red band (Sentinel-2 B04)
    :param nir:
        reflectance in the near infrared band (Sentinel-2 B08)
        spectrum
    :return:
        EVI values
    """

    return 2.5 * (nir - red) / (nir + 6*red - 7.5*blue + 1)


def MSAVI(
        red: np.array,
        nir: np.array
    ) -> np.array:
    """
    Calculates the Modified Soil-Adjusted Vegetation Index
    (MSAVI). MSAVI is sensitive to the green leaf area index
    (greenLAI).

    :param red:
        reflectance in the red band (Sentinel-2 B04)
    :param red_edge_3:
       reflectance in the NIR band (Sentinel-2 B08)
    :return:
        MSAVI values
    """
    return 0.5 * (2*nir + 1 - np.sqrt((2*nir + 1)**2 - 8*(nir - red)))


def CI_green(
        green: np.array,
        nir: np.array 
    ) -> np.array:
    """
    Calculates the green chlorophyll index (CI_green) using
    Sentinel-2 bands 3 (green) and 8 (nir) as suggested by Clevers
    et al. (2017, doi:10.3390/rs9050405). It is sensitive to
    canopy chlorophyll concentration (CCC).

    :param green:
        reflectance in the green band (Sentinel-2 B03)
    :param nir:
        reflectance in the NIR band (Sentinel-2 B08)
    """

    return (nir / green) - 1


def NDRE(
        red_edge_1: np.array,
        red_edge_2: np.array
    ) -> np.array:
    """
    Calculates the Normalized Difference Red Edge (NDRE)
    using Sentinel-2 bands 5 and 6

    :param red_edge_1:
        reflectance in the red edge 1 band (Sentinel-2 B05)
    :param red_edge_2:
        reflectance in the red edge 2 band (Sentinel-2 B06)
    """

    return (red_edge_2 - red_edge_1) / (red_edge_2 + red_edge_1)


def TCARI_OSAVI(
        green: np.array,
        red: np.array,
        red_edge_1: np.array,
        red_edge_3: np.array
    ) -> np.array:
    """
    Calculates the ratio of the Transformed Chlorophyll Index (TCARI)
    and the Optimized Soil-Adjusted Vegetation Index (OSAVI). It is sensitive
    to changes in the leaf chlorophyll content (LCC). The Sentinel-2 band
    selection follows the paper by Clevers et al. (2017, doi:10.3390/rs9050405)

    :param green:
        reflectance in the green band (Sentinel-2 B03)
    :param red:
        reflectance in the green band (Sentinel-2 B04)
    :param red_edge_1:
        reflectance in the red edge 1 band (Sentinel-2 B05)
    :param nir:
       reflectance in the red edge 3 band (Sentinel-2 B07)
    :return:
        TCARI/OSAVI values
    """

    TCARI = 3*((red_edge_1 - red) - 0.2*(red_edge_1 - green) * (red_edge_1 / red))
    OSAVI = (1 + 0.16) * (red_edge_3 - red) / (red_edge_3 + red + 0.16)
    tcari_osavi = TCARI/OSAVI
    # clip values to range between 0 and 1 (division by zero might cause inf)
    tcari_osavi[tcari_osavi < 0.] = 0.
    tcari_osavi[tcari_osavi > 1.] = 1.
    return tcari_osavi


def MCARI(
        green: np.array,
        red: np.array,
        red_edge_1: np.array
    ):
    """
    Calculates the Modified Chlorophyll Absorption Ratio Index (MCARI)
    using Sentinel-2 bands 3 (green), 4 (red), and 5 (red edge 1).
    It is sensitive to leaf chlorophyll concentration (LCC).

    :param green:
        reflectance in the green band (Sentinel-2 B03)
    :param red:
        reflectance in the red band (Sentinel-2 B04)
    :param nir:
        reflectance in the NIR band (Sentinel-2 B08)
    """

    return ((red_edge_1 - red) - 0.2 * (red_edge_1 - green)) * (red_edge_1 / red)


def calc_indices(
        in_file: Path,
        out_dir: Path
    ) -> None:
    """
    Calculates the NDVI, TCARI/OSAVI, MCARI and MSAVI for a band-stacked,
    resampled Sentinel-2 geoTiff file

    :param in_file:
        file-path to the resampled, band-stacked geoTiff file
    :param out_dir:
        directory where to write the results to. These have the
        same name as the input plus the name of the index calculated:
        <in_file.name>_<index>.tif
    """

    # open the file and read those bands required for calculating the indices
    # TODO: 10 and 20m band_data dict
    s2_band_data = {}
    with rio.open(in_file, 'r') as src:
        # get geo-referencation information
        meta = src.meta
        bounds = src.bounds
        # read relevant bands and store them in dict
        band_names = src.descriptions
        for idx, band_name in enumerate(band_names):
            if band_name in list(s2_band_mapping.keys()):
                # read as float; otherwise the division in numpy behaves strange
                # apply gain factors to rescale reflectanc between 0 and 1
                s2_band_data[s2_band_mapping[band_name]] = \
                    src.read(idx+1).astype(float) * 0.0001

    # update the file metadata for writing
    meta.update(
        {'count': 1, 'dtype': np.float32}
    )
    epsg = meta['crs'].to_epsg()

    # define output file name
    fname_base = out_dir.joinpath(f'VI_{in_file.name.split(".")[0]}').as_posix()

    # the actual index calculation starts here

    vis_names = ['NDVI', 'EVI', 'TCARI/OSAVI', 'MCARI', 'MSAVI', 'CI_green']
    vis = dict.fromkeys(vis_names)

    ###########    Indices commonly used to describe LSP    #############
    vis['NDVI'] = {
        'data': NDVI(
            red=s2_band_data['red'],
            nir=s2_band_data['nir']
        ),
        'fname': f'{fname_base}_NDVI.tif'
    }

    vis['EVI'] = {
        'data': EVI(
            blue=s2_band_data['blue'],
            red=s2_band_data['red'],
            nir=s2_band_data['nir']
        ),
        'fname': f'{fname_base}_EVI.tif'
    }

    ###########    Indices with physiological sensitivity   #############

    # sensitivity to leaf chlorophyll content
    vis['TCARI/OSAVI'] = {
        'data': TCARI_OSAVI(
            green=s2_band_data['green'],
            red=s2_band_data['red'],
            red_edge_1=s2_band_data['red_edge_1'],
            red_edge_3=s2_band_data['red_edge_3']
        ),
        'fname': f'{fname_base}_TCARI_OSAVI.tif'
    }

    # sensitivity to leaf chlorophyll concentration
    vis['MCARI'] = {
        'data': MCARI(
            green=s2_band_data['green'],
            red=s2_band_data['blue'],
            red_edge_1=s2_band_data['red_edge_1']
        ),
        'fname': f'{fname_base}_MCARI.tif'
    }

    # sensitivity to canopy chlorophyll concentration
    vis['CI_green'] = {
        'data': CI_green(
            green=s2_band_data['green'],
            nir=s2_band_data['nir']
        ),
        'fname': f'{fname_base}_CI-GREEN.tif'
    }
    
    # sensitivity to green leaf area
    vis['MSAVI'] = {
        'data': MSAVI(
            red=s2_band_data['red'],
            nir=s2_band_data['nir']
        ),
        'fname': f'{fname_base}_MSAVI.tif'
    }

    #### plotting ####
    n_rows, n_cols = 2, 3
    single_fig, single_axs = plt.subplots(n_rows, n_cols, figsize=(20,20))
    vi_idx = 0
    for row in range(n_rows):
        for col in range(n_cols):
            current_vi = vis[vis_names[vi_idx]]
            upper = np.quantile(current_vi['data'], 0.95) # 95% percentile
            lower = np.quantile(current_vi['data'], 0.05) # 5% percentile
            im_vi = single_axs[row, col].imshow(
                current_vi['data'], vmin=lower, vmax=upper, cmap='summer',
                extent=[bounds.left,bounds.right,bounds.bottom,bounds.top]
            )
            # add colormap: add_axes[left, bottom, width, heigth)
            divider = make_axes_locatable(single_axs[row,col])
            cax = divider.append_axes('right', size='5%', pad=0.05)
            single_fig.colorbar(im_vi, cax=cax, orientation='vertical')
            single_axs[row,col].title.set_text(vis_names[vi_idx])

            if row == n_rows - 1:
                single_axs[row,col].set_xlabel(f'X [m] (EPSG:{epsg})', fontsize=14)
                single_axs[row,col].xaxis.set_ticks(np.arange(bounds.left, bounds.right, 5000))
            if col == 0:
                single_axs[row,col].set_ylabel(f'Y [m] (EPSG:{epsg})', fontsize=14)
                single_axs[row,col].yaxis.set_ticks(np.arange(bounds.bottom, bounds.top, 5000))
                single_axs[row,col].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0f'))
            if col > 0:
                single_axs[row,col].set_yticklabels([])
            if row < n_rows - 1:
                single_axs[row,col].set_xticklabels([])
            single_axs[row,col].grid(False)
            vi_idx += 1
    single_fig.suptitle(
        'Vegetation Indices \n(lower and upper 5% percentiles exluded)',
        fontsize=20
    )
    single_fig.savefig(f'{fname_base}_maps.png', bbox_inches='tight')
    plt.close(single_fig)

    # write indices to output
    for vi in vis.values():
        with rio.open(vi['fname'], 'w', **meta) as dst:
            dst.write(vi['data'], 1)


def main(
        scenario_dir: Path,
        shapefile_study_area: Path
    ):
    # find scenes for which scenarios are available
    scenarios = glob.glob(scenario_dir.joinpath('S2*_MSIL1C*').as_posix())

    # define spatial resolution to resample (all)
    processing_options = {
        'in_file_aoi': shapefile_study_area,
        'resolution_selection': [10, 20]
    }
    target_resolution = 10

    # loop over scenarios
    for scenario in scenarios:

        # find L2A scenes
        orig_datasets_l2a = glob.glob(Path(scenario).joinpath('*/S2*_MSIL2A*.SAFE').as_posix())

        # loop over scenes, resample them for the extent of the study area and
        # calculate the spectral indices
        for orig_dataset in orig_datasets_l2a:

            # place results in the root of the scenario
            out_dir = Path(orig_dataset).parent

            # TODO: one 10m bandstack, one 20m bandstack???
            # bandstack, mask and resample the data
            path_bandstack = resample_and_stack_S2(
                in_dir=Path(orig_dataset),
                out_dir=out_dir,
                target_resolution=target_resolution,
                masking=True,
                pixel_division=True,
                is_L2A=True,
                **processing_options
            )
            # path_bandstack = Path(orig_dataset).parent.joinpath('20190420_T32TMT_MSIL2A_S2A_pixel_division_10m.tiff')

            # calculate the spectral indices using the resampled data
            calc_indices(
                in_file=path_bandstack,
                out_dir=out_dir
            )


if __name__ == '__main__':

    # scenario_dir = Path('./../S2A_MSIL1C_RUT-Scenarios/batch_*')
    scenario_dir = Path('/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL1C_RUT-Scenarios/batch_1')
    shapefile_study_area = './../shp/AOI_Esch_EPSG32632.shp'

    main(
        scenario_dir=scenario_dir,
        shapefile_study_area=shapefile_study_area
    )
