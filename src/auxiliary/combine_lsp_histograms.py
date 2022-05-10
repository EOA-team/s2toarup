'''
Created on Apr 4, 2022

@author: graflu
'''

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('ggplot')


if __name__ == '__main__':

    lsp_res_dir = Path(
        '../../S2_TimeSeries_Analysis_Test'
    )

    vis = ['EVI', 'NDVI', 'GLAI']
    runs = ['uncorrelated', 'fully_correlated']
    metrics = ['sos_times', 'eos_times']

    # read all data into a large data frame
    for run in runs:
        res = []
        f, axes = plt.subplots(nrows=2, ncols=3, figsize=(15,8))
        for vidx, vi in enumerate(vis):
            search_path = lsp_res_dir.joinpath(f'{vi}/{run}/Uncertainty_Maps/selected_crops')
            counter = 0
            for midx, metric in enumerate(metrics):
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

                sns.boxplot(x='crop', y=f'{metric} Uncertainty', data=df, ax=axes[midx,vidx])
                axes[midx,vidx].set_ylim(0,160)
                axes[midx,vidx].set_xlabel('')

                if vidx == 0:
                    axes[midx,vidx].set_ylabel(metric.split('_')[0].upper() + ' Uncertainty (k=1) [days]', fontsize=14)
                if midx == 0:
                    axes[midx,vidx].set_xlabel('')
                    axes[midx,vidx].xaxis.set_ticklabels([])
                    axes[midx,vidx].xaxis.set_ticks_position('none')
                    axes[midx,vidx].title.set_text(vi)
                if midx == 1:
                    plt.setp(axes[midx,vidx].get_xticklabels(), rotation=90)
                if vidx > 0:
                    axes[midx,vidx].set_ylabel('')
                    axes[midx,vidx].yaxis.set_ticklabels([])
                    axes[midx,vidx].yaxis.set_ticks_position('none')
        fpath_fig = lsp_res_dir.joinpath(f'{run}_lsp_boxplots.png')
        f.savefig(fpath_fig, dpi=300, bbox_inches='tight')
                
                
