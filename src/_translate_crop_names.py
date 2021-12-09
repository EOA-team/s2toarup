'''
Created on Dec 3, 2021

@author: graflu
'''

import geopandas as gpd


fname_shp = '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/shp/AOI_Esch_CropsNotFinal.shp'

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

gdf['crop_type'] = gdf.NUTZUNG.apply(lambda x, translation_dict=translation_dict: translation_dict[x])

fname_out = '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/shp/ZH_Polygons_2019_EPSG32632_selected-crops.shp'

gdf.to_file(fname_out)