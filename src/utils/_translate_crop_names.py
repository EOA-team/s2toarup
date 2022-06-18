'''
Translation of the attributes of the shapefile containing the main crop types
from German into English. Assigns internally used crop codes to allow for
rasterization of crop types.
'''

import geopandas as gpd


fname_shp = '../../shp/ZH_Polygons_2019_EPSG32632_selected-crops.shp'

gdf = gpd.read_file(fname_shp)

translation_dict = {
    'Körnermais': 'Corn',
    'Sonnenblumen zur Speiseölgewinnung': 'Sunflower',
    'Winterweizen ohne Futterweizen swissgranum': 'Winter Wheat',
    'Silo- und Grünmais': 'Silage Maize',
    'Winterraps zur Speiseölgewinnung': 'Canola',
    'Wintergerste': 'Winter Barley',
    'Zuckerrüben': 'Sugar Beet',
    'Kartoffeln': 'Potato',
    'Soja zur Speiseölgewinnung': 'Soybean',
    'Extensiv genutzte Wiesen (ohne Weiden)': 'Extensively Used Grasland',
    'Übrige Dauerwiesen (ohne Weiden)': 'Permament Grasland'
}

# assign custom crop codes
crop_codes = dict.fromkeys(translation_dict.values())
for idx, crop_code in enumerate(crop_codes):
    crop_codes[crop_code] = idx

gdf['crop_type'] = gdf.NUTZUNG.apply(lambda x, translation_dict=translation_dict: translation_dict[x])
gdf['crop_code'] = gdf.crop_type.apply(lambda x, crop_codes=crop_codes: crop_codes[x])

fname_out = '../shp/ZH_Polygons_2019_EPSG32632_selected-crops.shp'

gdf.to_file(fname_out)
