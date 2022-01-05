'''
Compares uncertainty outcomes with different number of scenario runs (150 vs. 1000)
'''

import glob
import numpy as np
import pandas as pd
from pathlib import Path
from agrisatpy.io import SatDataHandler


def compare(file_1, file_2, out_dir, vi_name):
    """compares uncertainty outcomes from different numbers of scenario runs"""

    handler_1 = SatDataHandler()
    handler_2 = SatDataHandler()
    handler_1.read_from_bandstack(file_1)
    handler_2.read_from_bandstack(file_2)

    # loop over bands and calculate the difference images
    diff_stats = []
    for band_name in handler_1.get_bandnames():
        
        diff = np.subtract(handler_1.get_band(band_name), handler_2.get_band(band_name))
        # add as band
        out_band_name = f'{band_name} DIFFERENCE'
        handler_2.add_band(
            band_name=out_band_name,
            band_data=diff
        )
        # save as plot and geoTiff
        fname_out_base = out_dir.joinpath(f'{vi_name}_{band_name}_difference').as_posix()
        fig_diff = handler_2.plot_band(out_band_name)
        fig_diff.savefig(f'{fname_out_base}.png', dpi=300, bbox_inches='tight')
        handler_2.write_bands(
            out_file=f'{fname_out_base}.tif',
            band_names=[out_band_name]
        )
        # save statistics
        diff_stats.append(
            {
                'band': band_name,
                'mean_difference': np.nanmean(diff),
                'median_difference': np.nanmedian(diff),
                'min_difference': np.nanmin(diff),
                'max_difference': np.nanmax(diff),
                'stddev_difference': np.nanstd(diff),
                'fifth-percentile_difference': np.nanquantile(diff, 0.05),
                'ninety-fifth-percentile_difference': np.nanquantile(diff, 0.95)
             }
        )
    df = pd.DataFrame(diff_stats)
    fname_csv = out_dir.joinpath(f'{vi_name}_difference_statistics.csv')
    df.to_csv(fname_csv, index=False)
    

if __name__ == '__main__':
    
    vi_names = ['EVI', 'NDVI']
    scenes = [
        'S2A_MSIL1C_20190328T102021_N0207_R065_T32TMT_20190328T154025',
        'S2A_MSIL1C_20190530T103031_N0207_R108_T32TMT_20190530T123429',
        'S2A_MSIL1C_20190818T103031_N0208_R108_T32TMT_20190818T124651'
    ]
    path_150 = Path('../S2A_MSIL2A_Analysis/150_scenarios')
    path_1000 = Path('../S2A_MSIL2A_Analysis')
    out_dir = path_1000.joinpath('Comparison')
    if not out_dir.exists():
        out_dir.mkdir()

    # loop over scenes
    for scene in scenes:

        # loop over vegetation indices/ parameters
        for vi_name in vi_names:
            # find the two uncertainty rasters
            search_expr = Path(scene).joinpath(f'L3_{vi_name}_*.tif')
            file_150_scenarios = glob.glob(path_150.joinpath(search_expr).as_posix())[0]
            file_1000_scenarios = glob.glob(path_1000.joinpath(search_expr).as_posix())[0]
            # compare them
            out_dir_scene = out_dir.joinpath(scene)
            if not out_dir_scene.exists():
                out_dir_scene.mkdir()
            compare(
                file_1=file_150_scenarios,
                file_2=file_1000_scenarios,
                out_dir=out_dir_scene,
                vi_name=vi_name
            )
            
