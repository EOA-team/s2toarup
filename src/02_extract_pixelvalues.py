"""
@purpose:    This script extracts and inserts the field parcel spectra into
            the PhenoDB.
            The vector geometries for those parcels used for extracting the
            pixel values must be already in the database (see 00_insert_parcels.py).
            
@author:    Lukas Graf (ETHZ)
"""
import os
import pandas as pd
from datetime import date
from datetime import datetime
from typing import List
from sqlalchemy import create_engine
from sqlalchemy.orm import close_all_sessions
import geopandas as gpd
from typing import Optional
from pathlib import Path

from agrisatpy.processing.extraction import S2bandstack2table
from agrisatpy.utils import reconstruct_path

from phenomen.config import get_settings
from phenomen.pheno_db import insert_pixel_observations
from phenomen.pheno_db import insert_parcel_observations
from phenomen.pheno_db.db_pixels.pixel_observations import update_pixel_observations

# connect to PhenoDB
Settings = get_settings()
engine = create_engine(Settings.DB_URL, echo=Settings.ECHO_DB)

# connect to Sentinel-2 metadata base
META_DB_URL = f'postgresql://{Settings.METADB_USER}:{Settings.METADB_PW}@{Settings.METADB_HOST}:{Settings.METADB_PORT}/{Settings.METADB_NAME}'
metadb_engine = create_engine(META_DB_URL, echo=Settings.ECHO_DB)


def insert_satellite_observations(
        date_start: date,
        date_end: date,
        tile: str,
        buffer: float,
        is_l2a: Optional[bool]=True,
        cloud_coverage_threshold: Optional[float]=100.,
        update: Optional[bool]=False
    ) -> List[dict]:
    """
    Sample script for inserting satellite (i.e., Sentinel-2) spectra
    of field parcel pixels and metadata into the database. Requires the
    Sentinel-2 data available in an AgrisatPy metadata base
    for querying the available Sentinel-2 data.

    :param date_start
        start date for extracting spectral data. It is recommended to place the
        starting date somewhat before the expected begin of the growing period/
        sowing date to make sure to capture the entire season.
    :param date_end:
        end date for extracting spectral data. It is recommended to place the
        end date somewhat after the expected end of the growing season/ harvest
        date to make sure to capture the entire season.
    :param tile:
        name of the Sentinel-2 tile (T<UTM-Zone><3-letter-code>) for which to extract
        data. Field parcels located outside of the tile are ignored. For multiple
        tiles call this function subsequently.
    :param buffer:
        value of the spatial buffer around the field geometries. Suggested is the
        use a buffer of -20 meters.
    :param is_l2a:
        If True (Default) assumes that the data is in L2A processing level
        (Sentinel-2 specific) and has scene classification layer (SCL) information
        from Sen2Cor.
    :param cloud_coverage_threshold:
        optional cloud coverage threshold (0-100%) to filter out scenes too cloudy.
        All scenes above the threshold will be discarded. By default, all scenes are
        used (threshold: 100%).
    :param update:
        if pixel values are inserted for the first time, leave this flag to False.
        If pixel values are should be updated, set the flag to True.
    """

    # translate processing level
    if is_l2a:
        processing_level = 'Level-2A'
    else:
        processing_level = 'Level-1C'

    # query the database to get the parcel geometries and read them into a
    # GeoDataFrame
    id_column = 'parcel_id'

    # first run - no pixels in the DB
    if not update:
        query = f'''
        SELECT
            fid AS {id_column},
            geom
        FROM
            {Settings.DEFAULT_SCHEMA}.{Settings.TABLE_FIELD_PARCELS};'''
    # update - pixels already in database
    else:
        query = f'''
        SELECT
            pixel.fid AS {id_column},
            pixel.pixid,
            pixel.x_orig AS x_coord,
            pixel.y_orig AS y_coord,
            parcel.geom
        FROM
            {Settings.DEFAULT_SCHEMA}.{Settings.TABLE_PARCEL_PIXELS} pixel
        LEFT JOIN {Settings.DEFAULT_SCHEMA}.{Settings.TABLE_FIELD_PARCELS} parcel
        ON pixel.fid=parcel.fid;'''
    gdf = gpd.read_postgis(sql=query, con=engine)
    
    # query the AgriSatPy metadata base
    query = f"""
    select
        proc.scene_id,
        proc.storage_device_ip, 
        proc.storage_device_ip_alias, 
        proc.storage_share, 
        proc.product_uri, 
        proc.bandstack, 
        proc.scl,
        raw.sensing_date,
        raw.spacecraft_name,
        raw.sun_zenith_angle,
        raw.sun_azimuth_angle,
        raw.sensor_zenith_angle,
        raw.sensor_azimuth_angle
    from sentinel2_processed_metadata as proc
    left join sentinel2_raw_metadata as raw
    on proc.product_uri=raw.product_uri
    where
        raw.sensing_date between '{date_start}' and '{date_end}'
    and
        raw.tile_id = '{tile}'
    and
        raw.processing_level = '{processing_level}'
    and
        raw.cloudy_pixel_percentage <= {cloud_coverage_threshold}
    order by sensing_date;
    """
    metadata = pd.read_sql(query, metadb_engine)

    errored = []
    # loop over the sensing dates available
    for idx, record in metadata.iterrows():

        # extract the relevant data to be inserted into the database
        date = record.sensing_date
        spacecraft = record.spacecraft_name
        sun_zenith = record.sun_zenith_angle
        sun_azimuth = record.sun_azimuth_angle
        obs_zenith = record.sensor_zenith_angle
        obs_azimuth = record.sensor_azimuth_angle
        scene_id = record.scene_id
        product_uri = record.product_uri
        bandstack = record.bandstack
        scl = record.scl

        # reconstruct dataset paths
        in_dir = reconstruct_path(
            record=record,
            is_raw_data=False
        )

        print(f'Extracting data from {product_uri} ({date}; {spacecraft}) ({idx+1}/{metadata.shape[0]})')

        if is_l2a:
            # some mess with paths...
            scl = Path(scl).as_posix().replace('\\', os.sep)
            in_file_scl = Path(in_dir).joinpath(scl)
        bandstack = Path(in_dir).joinpath(bandstack)

        # extract the pixel values
        try:
            if is_l2a:
                refl, scl_stats = S2bandstack2table(
                    in_file=bandstack,
                    in_file_scl=in_file_scl,
                    buffer=buffer,
                    id_column=id_column,
                    product_date=date,
                    in_gdf_polys=gdf,
                    filter_clouds=False
                )
                # sum up valid SCL classes
                scl_stats['perc_valid'] = scl_stats.vegetation + scl_stats.non_vegetated
            else:
                refl, _ = S2bandstack2table(
                    in_file=bandstack,
                    buffer=buffer,
                    id_column=id_column,
                    product_date=date,
                    in_gdf_polys=gdf,
                    is_l2a=False
                )
        except Exception as e:
            print(f'Could not extract pixel spectra: {e}')
            errored.append({'bandstack': bandstack, 'date': date,  'msg': e})
            continue

        # join with pixel_ids in case of an update
        if update:
            refl = pd.merge(refl, gdf, on=['x_coord', 'y_coord'])
            refl['parcel_id'] = refl.parcel_id_x

        # insert the Sentinel-2 spectra
        try:
            for _, pixel in refl.iterrows():
                # get the spectral bands (starting with 'b')
                pixel.index = [x.lower() for x in list(pixel.index)]
                value_dict = pixel.filter(regex=('b.*')).to_dict()

                # get the scene classification value (L2A)
                if is_l2a:
                    value_dict['scl'] = pixel['scl_class']
                # L1C processing level; rename reflectance values to bxx_toa
                else:
                    value_dict = {k+'_toa': v for k, v in value_dict.items() if k.startswith('b')}

                # save reference to product_uri and scene id (allows
                # metadata queries also on extracted pixels)
                value_dict['scene_id'] = scene_id
                value_dict['product_uri'] = product_uri

                # insert or update of pixel values
                if not update:
                    _ = insert_pixel_observations(
                        parcel_id=pixel.parcel_id,
                        x_orig=pixel.x_coord,
                        y_orig=pixel.y_coord,
                        epsg_orig=pixel.epsg,
                        date=date,
                        value_dict=value_dict
                    )
                else:
                    _ = update_pixel_observations(
                        parcel_id=pixel.parcel_id,
                        pixel_id=pixel.pixid,
                        date=date,
                        value_dict=value_dict
                    )

        except Exception as e:
            print(f'Could not insert pixel spectra into DB: {e}')
            errored.append({'parcel_id': pixel.parcel_id, 'date': date, 'msg': e})

        # insert parcel observations (i.e., Sentinel-2 metadata)
        try:
            for idx, parcel in gdf.iterrows():
                parcel_id = int(parcel.parcel_id)
                if parcel_id not in refl.parcel_id.unique():
                    continue
                # Level-1C: all pixels are considered valid for the time being
                perc_valid = 100.
                # Level-2A: use SCL information
                if is_l2a:
                    perc_valid = scl_stats[scl_stats.parcel_id == parcel_id]['perc_valid'].values[0]
                value_dict = {
                    'space_craft_name': spacecraft,
                    'perc_valid_pixels': perc_valid,
                    'sun_zenith_angle': sun_zenith,
                    'sun_azimuth_angle': sun_azimuth,
                    'obs_zenith_angle': obs_zenith ,
                    'obs_azimuth_angle': obs_azimuth
                }

                insert_parcel_observations(parcel_id=parcel_id,
                                           date=date,
                                           value_dict=value_dict
                )
        except Exception as e:
            print(f'Could not insert parcel satellite metadata into DB: {e}')
            errored.append({'parcel_id': pixel.parcel_id, 'date': date, 'msg': e})

    # close all db sessions in the end
    close_all_sessions()

    return errored


if __name__ == '__main__':

    date_start = datetime.strptime('2019-01-01', '%Y-%m-%d')
    date_end = datetime.strptime('2019-12-31', '%Y-%m-%d')
    buffer = -20. # meters
    tile = 'T32TMT'
    # processing level is L1C
    # is_l2a = False
    # or L2A
    is_l2a = True
    # first-run or update
    update = True
    # only use scenes with less than 20 percent cloud coverage
    cloud_coverage_percentage = 20.

    err = insert_satellite_observations(
        date_start=date_start,
        date_end=date_end,
        tile=tile,
        buffer=buffer,
        is_l2a=is_l2a,
        cloud_coverage_threshold=cloud_coverage_percentage,
        update=update
    )
