Readme.txt
----------

Field parcel boundaries and crop type information for 2019 (ZH_Polygons_2019_EPSG32632_selected-crops.shp) taken from Canton Zurich.
The polygons can be requested from https://www.geodienste.ch/services/lwb_nutzungsflaechen (Germnan only).

We shrink the polygons in our study area (AOI_Esch_EPSG32632.shp) to the following list of crop types and grasslands (original German name and English translation made by the authors):

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

The sample points (ZH_Points_2019_EPSG32632_selected-crops.shp) were derived from the field polygons:

* First, we randomly generated 5000 points within the spatial extent covered by the polygons.
* Second, we only kept those points laying within the field polygons features (ST_Within).
* To avoid mixed pixel effects close the field boundaries we applied a 20m inward buffering to all polygon features before selecting the points.
* We then assigned the crop type information to the remaining points by spatial join.

All vector pre-processing steps were carried out using QGIS 3.18.3-Zurich on Fedora34.

