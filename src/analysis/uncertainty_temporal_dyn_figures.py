'''
Created on Jun 3, 2022

@author: graflu
'''

from pathlib import Path
from PIL import Image

def combine_images(fpath_list, fpath_out):
    images = [Image.open(x) for x in fpath_list]
    widths, heights = zip(*(i.size for i in images))
    
    total_width = sum(widths)
    max_height = max(heights)
    
    new_im = Image.new('RGB', (total_width, max_height))
    
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]
    # crop central part, only
    # (left, upper, right, lower)
    box = (0, 350, total_width, max_height)
    cropped_img = new_im.crop(box)
    cropped_img.save(fpath_out)


if __name__ == '__main__':

    res_dir = Path('../../S2_TimeSeries_Analysis')
    vis = ['EVI', 'NDVI', 'GLAI']
    out_dir = Path('../../S2_TimeSeries_Analysis')

    crops = [
        'Canola', 'Corn', 'Extensively Used Grasland', 'Permament Grasland',
        'Silage Maize', 'Soybean', 'Sunflower', 'Sugar Beet', 'Winter Wheat', 'Winter Barley',
        'Potato'
    ]

    # loop over crops
    for crop in crops:
        # get paths to the VIs
        fpaths = []
        for vi in vis:
            fpaths.append(
                next(res_dir.glob(f'{vi}/{vi}_{crop}_all-pixel-timeseries.png'))
             )
        # merge images into a single figure
        fpath_out = out_dir.joinpath(f'{crop}_all-pixels-uncertainty-timeseries.png')
        combine_images(fpath_list=fpaths, fpath_out=fpath_out)
