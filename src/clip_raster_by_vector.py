#! /usr/bin/python3
'''
Created on Jul 12, 2021

@author:     Lukas Graf

@purpose:    Clips a raster file (single or multi-band)
             using a vector mask. Uses rasterio and avoids
             therefore direct GDAL calls.
             Intended to be used stand-alone.

Copyright 2021 Lukas Graf

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
or implied. See the License for the specific language governing
permissions and limitations under the License.
'''

import os
import glob
from pathlib import Path
import geopandas as gpd
import rasterio as rio
from rasterio.mask import mask

from pyproj import Transformer
from shapely.geometry import Polygon


def project_polygon(geom: Polygon,
                    src_epsg: int,
                    dst_epsg: int
                    ) -> Polygon:
    """
    projects a polygon shapely geometry from any spatial reference system
    into any target reference system and returns a new geom.

    :param geom:
        input shapely polygon geometry in custom SRID specified by the
        epsg code
    :param src_epsg:
        EPSG code of the input geometry
    :param dst_epsg:
        EPSG code of the target projection
    """
    transformer = Transformer.from_crs(src_epsg,
                                       dst_epsg,
                                       always_xy=True)
    x, y = geom.exterior.coords.xy
    lon, lat = transformer.transform(x, y)
    return Polygon(zip(list(lon), list(lat)))


def clip_raster_by_vector(in_file_vector: str,
                          in_file_raster: str,
                          out_file_raster: str
                          ) -> None:
    """
    clips a raster file by a vector mask assuming that
    the vector only contains one feature that is the mask.

    :param in_file_vector:
        vector mask to be used for clipping. Must be of type Polygon.
        If it is projected in a coordinate system different from the
        input raster, the geometry is projected into the coordinate
        system of the raster before clipping
    :param in_file_raster:
        raster file (preferably geoTiff because rasterio supports this
        file format best) to be clipped by the geometry.
    :param out_file_raster:
        filepath of the clipped raster (always geoTiff because rasterio
        seems to have some problems with other drivers)
    """
    # open vector file using geopandas
    gdf = gpd.read_file(in_file_vector)
    epsg = gdf.crs.to_epsg()
    # check if EPSG codes match and reproject polygon if necessary
    epsg_raster = rio.open(in_file_raster).crs
    geom = gdf['geometry'].iloc[0]
    if epsg != epsg_raster:
        geom = project_polygon(geom=geom,
                               src_epsg=epsg,
                               dst_epsg=epsg_raster)

    # read raster and clip it
    with rio.open(in_file_raster, 'r') as src:
        meta = src.meta
        nodata_raster = src.nodata
        data, transform = rio.mask.mask(src,
                                        [geom],
                                        crop=True,
                                        all_touched=True,
                                        nodata=nodata_raster
                                        )
        
    # update meta for writing
    meta.update({"driver": "GTiff",
                 "height": data.shape[1],
                 "width": data.shape[2], 
                 "transform": transform
                })

    # write clipped data as new raster file
    with rio.open(out_file_raster, 'w', **meta) as dst:
        dst.write(data)


# executable part
if __name__ == '__main__':
    
    in_file_vector = '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/STUDY_AREA/AOI_Esch_EPSG32632-large.shp'
    
    out_dir = '/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL2A_Analysis/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621'
    unc_dir = '/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL1C_orig/S2B_MSIL1C_20190830T102029_N0208_R065_T32TMT_20190830T130621.RUT'

    unc_files = glob.glob(f'{unc_dir}/*rut_b*.tif')

    for unc_file in unc_files:
        
        out_file = Path(unc_file).name
        out_file_raster = out_dir + os.sep + out_file

        try:
            clip_raster_by_vector(
                in_file_vector=in_file_vector,
                in_file_raster=unc_file,
                out_file_raster=out_file_raster
            )
            print('Clipped raster successfully!')
        except Exception as e:
            print(f'Failed to clip raster: {e}')
