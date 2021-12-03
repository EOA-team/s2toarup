'''
Created on Dec 1, 2021

@author: graflu
'''

import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from datetime import datetime
from typing import Optional
from typing import List
from pathlib import Path

plt.style.use('ggplot')

def plot_uncertainty_ts(
        analysis_dir: Path,
        output_dir: Path,
        roi_selection: Optional[List[str]]=None,
        search_expression: Optional[str]='spectral-band_l1c-l2a-l3_uncertainty*.csv'
    ):
    """
    Plots uncertainty extracted for selected land use classes (e.g., different
    crop types) and its development over time.

    :param analysis_dir:
        directory where the uncertainty analysis results are stored
    :param output_dir:
        directory where to store the outputs (resulting graphs)
    :param search_expression:
        string expression with wildcard to identify the CSV files With
        extracted uncertainty
    """

    # search for CSV files with uncertainty values. These are stored per scene
    # (starting with S2*) in a sub-directory named 'csv'
    search_wildcard = analysis_dir.joinpath(f'S2*/csv/{search_expression}').as_posix()
    csv_files = glob.glob(search_wildcard)

    # loop over files and read them into a list of pandas dataframes
    uncertainty_df_list = []
    for csv_file in csv_files:
        uncertainty_df_list.append(pd.read_csv(csv_file))

    # merge them into a single df
    uncertainty_df = pd.concat(uncertainty_df_list)
    # parse date to pandas datetime
    uncertainty_df.date = pd.to_datetime(uncertainty_df.date)

    # get processing levels; one plot for each processing level
    processing_levels = uncertainty_df.processing_level.unique()

    for processing_level in processing_levels:
        uncertainty_df_pl = uncertainty_df[uncertainty_df.processing_level == processing_level].copy()
        # order data by date
        uncertainty_df_pl = uncertainty_df_pl.sort_values(by='date')
        uncertainty_df_pl.dropna(axis=1, how='all', inplace=True)
        # get unique ROIs
        rois = uncertainty_df_pl.ROI.unique()
        # check for user-defined ROI selection (subset)
        if roi_selection is not None:
            roi_selection_set = set(roi_selection)
            rois = [x for x in rois if x in roi_selection_set]

        # identify the quantities for which the uncertainty is available
        uncertainty_quantities = list(uncertainty_df_pl.columns)
        values_to_remove = {'date', 'processing_level', 'ROI', 'area_km2'}
        uncertainty_quantities = [x for x in uncertainty_quantities if x not in values_to_remove]

        # plot ROI uncertainty per date and roi for each spectral band/ VI/ parameter
        for uncertainty_quantity in uncertainty_quantities:
            q_df = uncertainty_df_pl[['date', 'ROI', 'area_km2', uncertainty_quantity]].copy()
            fig, ax = plt.subplots(1, 1, clear=True, num=1, figsize=(15,10))
            for roi in rois:
                q_df_roi = q_df[q_df.ROI == roi].copy()
                label = f'{roi} ({np.round(q_df_roi.area_km2.values[0],2)}' + ' $km^2$)'
                ax.plot(q_df_roi.date, q_df_roi[uncertainty_quantity], label=label)
            ax.title.set_text(f'{processing_level} {uncertainty_quantity} Uncertainty Averaged per Crop Type')
            ax.set_ylabel('Relative Uncertainty (k=1) [%]', fontsize=14)
            # plot legend outside of plot
            ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
            fname = output_dir.joinpath(f'{processing_level}_{uncertainty_quantity}_Rel-Uncertainty-TS.png')
            fig.savefig(fname, bbox_inches='tight')
            plt.close(fig)


if __name__ == '__main__':
    
    analysis_dir = Path(
        '/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL2A_Analysis'
    )

    output_dir = analysis_dir.joinpath('uncertainty_ts_plots')
    if not output_dir.exists(): output_dir.mkdir()

    plot_uncertainty_ts(
        analysis_dir=analysis_dir,
        output_dir=output_dir
    )
        