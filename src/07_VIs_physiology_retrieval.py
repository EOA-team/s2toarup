
import rasterio as rio
import numpy as np
from pathlib import Path

# map the Sentinel-2 bands to be used for index calculation
s2_band_mapping = {
    'green': 'B03',
    'red': 'B04',
    'red_edge_1': 'B05',
    'red_edge_3': 'B07',
    'nir': 'B08'
}

def NDVI(
        red: np.array,
        nir: np.array
    ) -> np.array:
    """
    Calculates the Normalized Difference Vegetation Index
    (NDVI).

    :param red:
        reflectance in the red part of the electro-magnetic spectrum
    :param nir:
        reflectance in the near infrared part of the electro-magnetic
        spectrum
    :return:
        NDVI values
    """

    return (nir - red) / (nir + red)


def TCARI_OSAVI(
        green: np.array,
        red: np.array,
        red_edge_1: np.array,
        red_edge_3: np.array
    ) -> np.array:
    """
    Calculates the ratio of the Transformed Chlorophyll Index (TCARI)
    and the Optimized Soil-Adjusted Vegetation Index (OSAVI).

    :param green:
        reflectance in the green part of the electro-magnetic spectrum
    :param red:
        reflectance in the red part of the electro-magnetic spectrum
    :param red_edge_1:
        reflectance in the red_edge_1 part of the electro-magnetic spectrum
    :param red_edge_3:
       reflectance in the red_edge_3 part of the electro-magnetic spectrum
    :return:
        TCARI/OSAVI values
    """

    TCARI = 3*((red_edge_1 − red) − 0.2*(red_edge_1 − green) * (red_edge_1 / red))
    OSAVI = (1 + 0.16) * (red_edge_3 − red) / (red_edge_3 + red + 0.16)
    return TCARI/OSAVI


def MSAVI(
        red: np.array,
        red_edge_3: np.array
    ) -> np.array:
    """
    Calculates the Modified Soil-Adjusted Vegetation Index
    (MSAVI).

    :param red:
        reflectance in the red part of the electro-magnetic spectrum
    :param red_edge_3:
       reflectance in the red_edge_3 part of the electro-magnetic spectrum
    :return:
        MSAVI values
    """

    return 0.5 * (2*red_edge_3 + 1 - np.sqrt((2*red_edge_3 + 1)**2 - 8*(red_edge_3 - red)))


def MCARI(
        green: np.array,
        red: np.array,
        red_edge_1: np.array
    ) -> np.array:
    """
    Calculates the Modified Chlorophyll Absorption Ratio Index
    (MCARI).

    :param green:
        reflectance in the green part of the electro-magnetic spectrum
    :param red:
        reflectance in the red part of the electro-magnetic spectrum
    :param red_edge_1:
        reflectance in the red_edge_1 part of the electro-magnetic spectrum
    :return:
        MCARI values
    """

    return ((red_edge_1 - red) - 0.2*(red_edge_1 - green)) * (red_edge_1 / red)


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
    s2_band_data = dict.fromkeys(s2_band_mapping.keys())
    with rio.open(in_file, 'r') as src:
        # get geo-referencation information
        meta = src.meta
        # read relevant bands and store them in dict
        band_names = src.descriptions()
        for idx, band_name in enumerate(band_names):
            if band_name in list(s2_band_mapping.keys()):
                s2_band_data[band_name] = src.read(idx+1)

    # update the file metadata for writing
    meta.update(
        'count': 1
    )
    # define output file name
    fname_base = out_dir.joinpath(f'VI_{in_file.name.split(".")[0]}').as_posix()

    # the actual index calculation starts here
    ndvi = NDVI(
        red=s2_band_data['red'],
        nir=s2_band_data['nir']
    )
    fname_ndvi = f'{fname_base}_NDVI.tif'

    tcari_osavi = TCARI_OSAVI(
        green=s2_band_data['green'],
        red=s2_band_data['red'],
        red_edge_1=s2_band_data['red_edge_1'],
        red_edge_3=s2_band_data['red_edge_3']
    )
    fname_tcari_osavi = f'{fname_base}_TCARI_OSAVI.tif'

    mcari = MCARI(
        green=s2_band_data['green'],
        red=s2_band_data['red'],
        red_edge_1=s2_band_data['red_edge_1']
    )
    fname_mcari = f'{fname_base}_MCARI.tif'

    msavi = MSAVI(
        red=s2_band_data['red'],
        red_edge_1=s2_band_data['red_edge_3']
    )
    fname_msavi = f'{fname_base}_MSAVI.tif'

    # write indices to output
    fnames = [fname_ndvi, fname_tcari_osavi, fname_mcari, fname_msavi]
    data_arrays = [ndvi, tcari_osavi, mcari, msavi]
    for fname, array in zip(fnames, data_arrays):
        with rio.open(fname, 'w', **meta) as dst:
            dst.write(array, 1)


def main(
        scenario_dir: Path,
        shapefile_study_area: Path
    ):
    # find scenes for which scenarios are available
    scenarios = glob.glob(scenario_dir.joinpath('S2*_MSIL1C*').as_posix())

    # define spatial resolution to resample (all)
    processing_options = {
        'in_file_aoi': shapefile_study_area,
        'resolution_selection': [10, 20, 60]
    }
    target_resolution = 10

    # loop over scenarios
    for scenario in scenarios:

        # find L2A scenes
        orig_datasets = glob.glob(scenarios.joinpath('*/S2*_MSIL2A*.SAFE').as_posix())

        # loop over scenes, resample them for the extent of the study area and
        # calculate the spectral indices
        for orig_dataset in orig_datasets:

            # place results in the root of the scenario
            out_dir = Path(orig_dataset).parent

            # bandstack, mask and resample the data
            path_bandstack = resample_and_stack_S2(
                in_dir=Path(orig_dataset),
                out_dir=orig_dataset_out_dir,
                target_resolution=target_resolution,
                masking=True,
                pixel_division=True,
                is_L2A=True,
                **processing_options
            )

            # calculate the spectral indices using the resampled data


if __name__ == '__main__':

    scenario_dir = './../S2A_MSIL1C_RUT-Scenarios'
    shapefile_study_area = './../shp/AOI_Esch_EPSG32632.shp'
    
    