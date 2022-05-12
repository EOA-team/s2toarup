'''
Created on Apr 4, 2022

@author: graflu
'''

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
import numpy as np
import seaborn as sns

plt.style.use('ggplot')


if __name__ == '__main__':

    lsp_res_dir = Path(
        '../../S2_TimeSeries_Analysis'
    )

    vis = ['EVI', 'NDVI', 'GLAI']
    runs = ['uncorrelated', 'fully_correlated']
    metrics = ['sos_times', 'eos_times', 'length_of_season']

    # read all data into a large data frame
    for metric in metrics:
        f, axes = plt.subplots(nrows=2, ncols=3, figsize=(15,8))
        for vidx, vi in enumerate(vis):
            for rdx, run in enumerate(runs):
                search_path = lsp_res_dir.joinpath(f'{vi}/{run}/Uncertainty_Maps/selected_crops')
                fpath = next(search_path.glob(f'{vi}_{metric}_*_data.csv'))
                df = pd.read_csv(fpath)
                # drop grassland pixels (no meaningful LSP metrics)
                # drop_idx = df[
                #     df.crop.isin(['Extensively Used Grasland', 'Permament Grasland'])
                # ].index
                # df.drop(index=drop_idx, inplace=True)
                df['crop'] = df['crop'].apply(lambda x: 'Rapeseed' if x == 'Canola' else x)
                df['crop'] = df['crop'].apply(lambda x: 'Grain Maize' if x == 'Corn' else x)
                # crop_count = df.crop.value_counts()
                # df['crop'] = df['crop'].apply(lambda x, crop_count=crop_count:
                #     f'{x} ({crop_count[crop_count.index == x].values[0]})'
                # )

                sns.boxplot(x='crop', y=f'{metric} Uncertainty', data=df, ax=axes[rdx,vidx])
                axes[rdx,vidx].set_ylim(0,160)
                axes[rdx,vidx].set_xlabel('')
                # f.suptitle(r'(a)', fontsize=22, x=0.14, y=0.925)

                if vidx == 0:
                    if run == 'uncorrelated':
                        axes[rdx,vidx].set_ylabel('Zero Scene Correlation', fontsize=14)
                    else:
                        axes[rdx,vidx].set_ylabel('Full Scene Correlation', fontsize=14)
                if vidx == 2:
                    label = metric.split('_')[0].upper() + ' Uncertainty (k=1) [days]'
                    if metric == 'length_of_season':
                        label = 'LOS Uncertainty (k=1) [days]'
                    axes[rdx,vidx].set_ylabel(
                        label,
                        fontsize=14,
                        rotation=270,
                        labelpad=14
                    )
                    axes[rdx,vidx].yaxis.set_label_position("right")
                    axes[rdx,vidx].yaxis.tick_right()

                # if vidx == 0 and rdx == 1:
                    # axes[rdx,vidx].set_title(r'(b)', fontsize=22, loc='left')
                if rdx == 0:
                    axes[rdx,vidx].set_xlabel('')
                    axes[rdx,vidx].xaxis.set_ticklabels([])
                    axes[rdx,vidx].xaxis.set_ticks_position('none')
                    axes[rdx,vidx].title.set_text(vi)
                if rdx == 1:
                    plt.setp(axes[rdx,vidx].get_xticklabels(), rotation=90)
                if 0 < vidx < 2:
                    axes[rdx,vidx].set_ylabel('')
                    axes[rdx,vidx].yaxis.set_ticklabels([])
                    axes[rdx,vidx].yaxis.set_ticks_position('none')


        fpath_fig = lsp_res_dir.joinpath(f'{metric}_uncertainty.png')
        f.savefig(fpath_fig, dpi=300, bbox_inches='tight')
                
                
