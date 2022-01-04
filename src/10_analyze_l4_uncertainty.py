'''
Calculates the absolute uncertainty in the phenological metrics including
start, peak and end of season (timing and vegetation index/ parameter values at
these phenological stages).
'''

import glob
import numpy as np
from pathlib import Path

from agrisatpy.io import SatDataHandler
from logger import get_logger
from copy import deepcopy

logger = get_logger('l5_uncertainty')


def calc_l4_uncertainty(
        uncertainty_dir: Path,
        out_dir: Path,
        vi_name: str
    ):
    """
    Assesses the uncertainty in the phenological metrics derived in the
    previous step. Reports absolute uncertainties, only.

    :param uncertainty_dir:
        directory containing the phenometric results from the time series
        scenarios of the vegetation indices/ parameters analyzed
    """

    # search the scenarios, organized by vegetation index/ parameter
    vi_uncertainty_dir = uncertainty_dir.joinpath(vi_name)

    # search scenarios
    vi_search_expr = vi_uncertainty_dir.joinpath('*/pheno_metrics.tif').as_posix()
    scenarios = glob.glob(vi_search_expr)

    # pheno-metrics available
    handler_list = []
    for idx, scenario in enumerate(scenarios):
        handler = SatDataHandler()
        handler.read_from_bandstack(
            fname_bandstack=Path(scenario)
        )
        handler_list.append(handler)
        handler = None
        logger.info(f'Reading scenario {idx+1}/{len(scenarios)} ({scenario})')

    # calculate the absolute uncertainty for each phenological metric
    pheno_metrics = dict.fromkeys(handler_list[0].get_bandnames())

    for pheno_metric in pheno_metrics:
        
        # get bandstack of all scenarios of the pheno metric to calculate the standard
        # deviation (=standard uncertainty)
        stack_list = [x.get_band(pheno_metric) for x in handler_list]
        stack_array = np.stack(stack_list)
        standard_unc = np.nanstd(stack_array, axis=0)

        # save to raster files and create a preview plot
        unc_handler = deepcopy(handler_list[0])
        band_name = f'{pheno_metric} Uncertainty'
        unc_handler.add_band(band_name=band_name, band_data=standard_unc)
        fig_unc = unc_handler.plot_band(band_name, colormap='summer')
        fname_out_fig = out_dir.joinpath(f'{vi_name}_{pheno_metric}_abs-uncertainty.png')
        fig_unc.savefig(fname_out_fig, dpi=300, bbox_inches='tight')
        fname_out_raster = fname_out_fig.as_posix().replace('.png','.tif')
        unc_handler.write_bands(
            out_file=fname_out_raster,
            band_names=[band_name]
        )


if __name__ == '__main__':

    uncertainty_dir = Path('../S2_TimeSeries_Analysis')
    out_dir = uncertainty_dir.joinpath('Uncertainty_Maps')
    if not out_dir.exists():
        out_dir.mkdir()

    vi_names = ['NDVI', 'EVI']

    for vi_name in vi_names:
        calc_l4_uncertainty(
            uncertainty_dir=uncertainty_dir,
            out_dir=out_dir,
            vi_name=vi_name
        )
    