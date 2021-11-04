
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio as rio
import matplotlib.pyplot as plt

from agrisatpy.processing.resampling.sentinel2 import resample_and_stack_S2


def calc_band_covariance(
        image_path: Path
    ) -> np.array:
    """
    Calculates the cross-band correlation and covariance
    matrix for all spectral bands in an image.

    :param image_path:
        file path to the multi-band image
    :returns:
        Tuple containing two numpy arrays. First item
        with the correlation matrix, second containing the
        covariance matrix
    """

    # read image data
    with rio.open(image_path, 'r') as src:
        img = src.read()

    # reshape the data
    img_reshaped = img.reshape(img.shape[1]*img.shape[2], -1)

    # return the correlation, covariance matrix
    return (np.corrcoef(img_reshaped.T), np.cov(img_reshaped.T))
    

def main(
        orig_datasets_dir: Path,
        out_dir: Path,
        shapefile_study_area: Path
    ) -> None:
    """
    executable function

    :param orig_datasets_dir:
        directory where the L1C (.SAFE) scenes are stored
    :param out_dir:
        directory where to store the spatially resampled L1C
        scenes clipped to the study area and their spectral
        covariance matrix
    :param shapefile_study_area:
        shapefile defining the study area.
    """

    # find L1C scenes
    orig_datasets = glob.glob(orig_datasets_dir.joinpath('*.SAFE').as_posix())

    # define spatial resolution to resample (all)
    processing_options = {
        'in_file_aoi' : shapefile_study_area,
        'resolution_selection' : [10, 20, 60]
    }
    target_resolution = 10

    # loop over scenes, resample them for the extent of the study area and
    # calculate the covariance matrix between the bands of Sentinel-2
    for orig_dataset in orig_datasets:
        
        orig_dataset_out_dir = out_dir.joinpath(f'{Path(orig_dataset).name.split(".")[0]}.R10m')
        if not orig_dataset_out_dir.exists():
            orig_dataset_out_dir.mkdir()

        # bandstack, mask and resample the data
        path_bandstack = resample_and_stack_S2(
            in_dir=Path(orig_dataset),
            out_dir=orig_dataset_out_dir,
            target_resolution=target_resolution,
            masking=True,
            pixel_division=True,
            is_L2A=False,
            **processing_options
        )

        # use the bandstacked file to calculate the covariance matrix
        corr_ma, cov_ma = calc_band_covariance(image_path=path_bandstack)

        # plot correlation matrix and save plot
        df = pd.DataFrame(corr_ma)
        s2_bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
        df.columns = s2_bands
        df.index = s2_bands
        vals = np.around(df.values,2)

        norm = plt.Normalize(vals.min()-1, vals.max()+1)
        colours = plt.cm.coolwarm(norm(vals))

        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, frameon=False, xticks=[], yticks=[])
        the_table=plt.table(cellText=vals, rowLabels=df.index, colLabels=df.columns, 
                            colWidths = [0.03]*vals.shape[1], loc='center', 
                            cellColours=colours)
        fname = orig_dataset_out_dir.joinpath('correlation-matrix.png')
        fig.savefig(fname)

        # save to disk for the subsequent processing step
        fname = orig_dataset_out_dir.joinpath('covariance-matrix.txt')
        np.savetxt(fname, cov_ma, delimiter=',')
        fname = orig_dataset_out_dir.joinpath('correlation-matrix.txt')
        np.savetxt(fname, corr_ma, delimiter=',')


if __name__ == '__main__':
    
    ### define user inputs

    # directory with L1C data (.SAFE subdirectories)
    orig_datasets_dir = Path('./../S2A_MSIL1C_orig')
    out_dir = orig_datasets_dir

    # shapefile defining the study area bounds
    shapefile_study_area = Path('./../shp/AOI_Esch_EPSG32632.shp')

    main(
        orig_datasets_dir=orig_datasets_dir,
        out_dir=out_dir,
        shapefile_study_area=shapefile_study_area
    )

