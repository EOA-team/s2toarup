'''
Created on Mar 30, 2022

@author: graflu
'''

from datetime import date
from pathlib import Path
import pandas as pd


if __name__ == '__main__':

    scenario_dir = Path(
        '/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2_MSIL1C_RUT-Scenarios'
    )

    output_dir = scenario_dir

    # loop over batches, read the SCL relative number of pixels and combine it
    # into a single data frame
    res_list = []
    for batch in range(1,6):
        scl_res_file = next(
            scenario_dir.rglob(
                f'batch_{batch}/SCL_Uncertainty/SCL_relative-number-' \
                'of-pixels-per-class_abs-uncertainty.csv'
            )
        )
        scl_res = pd.read_csv(scl_res_file)
        res_list.append(scl_res)

    scl_df = pd.concat(res_list)

    # create Latex table
    col_selection = [
        'date', 'cloud_shadows', 'vegetation', 'non_vegetated', \
        'cloud_medium_probability', 'cloud_high_probability'
    ]
    # select interesting dates
    sel_dates = [
        date(2019,2,14), # very early in the year
        date(2019,4,20), # spring
        date(2019,5,30), # cumulus clouds
        date(2019,7,16), # close the harvest of cereals
        date(2019,9,19)  # early autumn
    ]
    scl_df['date'] = pd.to_datetime(scl_df['date'])
    scl_df_sel = scl_df[scl_df['date'].isin(sel_dates)].copy()
    scl_df_sel.sort_values(by='date', inplace=True)
    table = scl_df_sel[col_selection].to_latex(index=False)

    fname = output_dir.joinpath('SCL_relative-number-of-pixels-per-class_abs-uncertainty.csv')
    with open(fname, 'w+') as dst:
        dst.write(table)

