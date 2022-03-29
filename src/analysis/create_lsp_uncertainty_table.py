'''
Combines LSP uncertainty results into single handy dataframes (tables
for the paper)
'''

import sys
import pandas as pd

from pathlib import Path
from typing import List


def combine_results(
        lsp_res_dir: Path,
        vis: List[str],
        interscene_correlations: List[str],
        lsp_metrics: List[str],
        output_dir: Path
    ) -> None:
    """
    Combines the results of the vegetation indices/ physiological parameters
    and their LSP uncertainty into DataFrames
    """
    for interscene_correlation in interscene_correlations:
        res = []
        for lsp_metric in lsp_metrics:
            df_vis = None
            col_selection = ['crop', 'median', 'mean']
            for vidx, vi in enumerate(vis):
                try:
                    path_res = next(
                        lsp_res_dir.rglob(
                            f'{vi}/{interscene_correlation}/Uncertainty_Maps/' + \
                            f'selected_crops/{vi}_{lsp_metric}_*uncertainty-' + \
                            f'crops-stats.csv'
                        )
                    )
                except Exception as e:
                    print(e)
                    sys.exit()
                df_vi = pd.read_csv(path_res)
                col_mapping = {'median': f'{vi} median', 'mean': f'{vi} mean'}
                try:
                    if vidx == 0:
                        df_vis = df_vi[col_selection].copy()
                        df_vis.rename(
                            columns=col_mapping,
                            inplace=True
                        )
                    else:
                        df_vi.rename(columns=col_mapping, inplace=True)
                        col_selection[1:3] = list(col_mapping.values())
                        df_vis = pd.merge(df_vis, df_vi[col_selection], on='crop', how='left')
                except Exception as e:
                    print(e)
            df_vis['metric'] = df_vi['metric_alias']
            res.append(df_vis)
        df_interscene_corr = pd.concat(res)
        # save to CSV
        fname_out = output_dir.joinpath(
            f'{interscene_correlation}_LSP_Uncertainty.csv'
        )
        df_interscene_corr.to_csv(fname_out, index=False)

if __name__ == '__main__':

    # lsp_res_dir = Path('../../S2_TimeSeries_Analysis')
    lsp_res_dir = Path('/home/graflu/public/Evaluation/Projects/KP0031_lgraf_PhenomEn/01_Uncertainty/ESCH/scripts_paper_uncertainty/S2_TimeSeries_Analysis')
    vis = ['EVI', 'NDVI', 'GLAI']
    interscene_correlations = ['fully_correlated', 'uncorrelated']
    lsp_metrics = ['sos_times', 'eos_times']

    output_dir = lsp_res_dir

    combine_results(lsp_res_dir, vis, interscene_correlations, lsp_metrics, output_dir)

