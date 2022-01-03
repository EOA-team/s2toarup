#!/usr/bin/python

"""
This script shows how to download Sentinel-2 data from Creodias for a given
region of interest, date range, cloud coverage and processing level.
            
It requires an account at Creodias (https://creodias.eu/) and the account's
username and password set in the environmental variables as
            
    SET CREODIAS_USER=<your-username>
    SET CREODIAS_PASSWORD0<your-password>

It is the FIRST script in the uncertainty processing chain
"""

import geopandas as gpd
from datetime import date
from pathlib import Path
from typing import Optional

from agrisatpy.downloader.sentinel2.creodias import query_creodias
from agrisatpy.downloader.sentinel2.creodias import download_datasets
from agrisatpy.downloader.sentinel2.utils import unzip_datasets
from agrisatpy.utils.constants import ProcessingLevels


def main(
        start_date: date,
        end_date: date,
        processing_level: ProcessingLevels,
        cloud_cover_threshold: int,
        aoi_file: Path,
        download_dir: Path,
        max_records: Optional[int]=200
    ) -> None:
    """
    the main executable function that calls all the required
    processing steps
    """
    bbox_data = gpd.read_file(aoi_file)

    # project to geographic coordinates (required for API query)
    bbox_data.to_crs(4326, inplace=True)
    # use the first feature (all others are ignored)
    bounding_box = bbox_data.geometry.iloc[0]

    # check for available datasets
    datasets = query_creodias(
        start_date=start_date,
        end_date=end_date,
        max_records=max_records,
        processing_level=processing_level,
        bounding_box=bounding_box,
        cloud_cover_threshold=cloud_cover_threshold
    )

    try:
        download_datasets(datasets, download_dir)
    except Exception as e:
        print(e)

    # unzip files
    unzip_datasets(download_dir=download_dir)


if __name__ == '__main__':
    
    ### user inputs
    # processing level
    processing_level = ProcessingLevels.L1C
    # date range
    start_date = date(2019,3,1)
    end_date = date(2019,3,20)
    # max_records defines the maximum number of datasets to download, increase if
    # necessary; however, CREODIAS might impose a limitation...
    max_records = 200

    # filter by cloud cover (all scenes with a cloud cover lower than the threshold
    # will be downloaded)
    cloud_cover_threshold = 20

    # shapefile defining the bounds of your region of interest
    aoi_file = Path(
        './../shp/AOI_Esch_EPSG32632.shp'
    )
    download_dir = Path('./../S2A_MSIL1C_orig')

    main(
        start_date=start_date,
        end_date=end_date,
        processing_level=processing_level,
        cloud_cover_threshold=cloud_cover_threshold,
        aoi_file=aoi_file,
        download_dir=download_dir
    )
