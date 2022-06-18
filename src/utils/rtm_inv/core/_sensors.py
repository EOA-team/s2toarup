'''
Sensors currently supported for RTM inversion
'''

class Sensors(object):
    """
    Class of classes of sensors
    """

    class Sentinel2:
        """
        defines properties shared by all members of the Sentinel2 platform
        """
        
        num_bands = 13
        band_names = [
            'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B09', 'B10',
            'B11', 'B12'
        ]
        central_wvls = [
            443, 490, 560, 665, 705, 740, 783, 842, 865, 945, 1375, 1610, 2190
        ]
        band_names_long = [
            'COSTAL-AEROSOL', 'BLUE', 'GREEN', 'RED', 'REDEDGE1',
            'REDEDGE2', 'REDEDGE3', 'NIR', 'NARROW-NIR', 'WATERVAPOUR',
            'SWIR-CIRRUS', 'SWIR1', 'SWIR2'
        ]

    class Sentinel2A(Sentinel2):
        """
        defines Sentinel2A-MSI
        """
        def __init__(self):
            self.name = 'Sentinel2A-MSI'
            self.band_widths = [
                21, 66, 36, 31, 15, 15, 20, 106, 21, 20, 31, 91, 175
            ]

    class Sentinel2B(Sentinel2):
        """
        defines Sentinel2B-MSI
        """
        def __init__(self):
            self.name = 'Sentinel2B-MSI'
            self.band_widths = [
                21, 66, 36, 31, 16, 15, 20, 106, 22, 21, 30, 94, 185
            ]

    class Landsat8:
        """
        defines Landsat8-OLI
        """
        name = 'LANDSAT8-OLI'
        num_bands = 9
        # order of wvl of bands is a bit weired but it's the official ordering
        band_names = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9']
        central_wvls = [440, 480, 560, 655, 865, 1610, 2200, 590, 1370]
