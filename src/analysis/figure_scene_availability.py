'''
Plots the figure showing the availability of S2 scenes in the
study area alongside their cloud coverage.
'''

from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle

plt.style.use('ggplot')
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['legend.fontsize'] = 16

fpath_s2_scenes = Path('../../data/s2_scenes.csv')
df = pd.read_csv(fpath_s2_scenes)
df['sensing_date'] = pd.to_datetime(df['sensing_date'])

f = plt.figure(figsize=(8,5))
ax = f.add_subplot(111)

sns.scatterplot(
    x='sensing_date',
    y='cloudy_pixel_percentage',
    hue='spacecraft_name',
    data=df,
    ax=ax,
    s=60
)
ax.set_xlabel('Sensing Date [YYYY-MM]', fontsize=16)
ax.set_ylabel('Cloudy Pixel Percentage [%]', fontsize=16)

# start = mdates.date2num(d=df.sensing_date.iloc[0])
# end = mdates.date2num(d=df.sensing_date.iloc[18])
# width = end - start
# rect = Rectangle((start, 0), width, 20, alpha=0.2)
# ax.add_patch(rect)

text_x, text_y = mdates.date2num(df.sensing_date.iloc[8]), 20.1
ax.text(text_x, text_y, 'PDGS 02.07', color='blue', fontsize=16)

# draw a vertical line to indicate the transition from one baseline to another
ax.vlines(x=df.sensing_date.iloc[18], ymin=0, ymax=20,
          linestyle='dashed', linewidth=3, color='grey')

# start = mdates.date2num(d=df.sensing_date.iloc[18])
# end = mdates.date2num(d=df.sensing_date.iloc[-1])
# width = end - start
# rect = Rectangle((start, 0), width, 20, alpha=0.2, color='g')
# ax.add_patch(rect)

text_x, text_y = mdates.date2num(df.sensing_date.iloc[21]), 20.1
ax.text(text_x, text_y, 'PDGS 02.08', color='g', fontsize=16)

ax.set_ylim(0,20)
ax.set_xlim(df.sensing_date.iloc[0], df.sensing_date.iloc[-1])

ax.legend(bbox_to_anchor=(0.5, -0.25), loc='lower center', borderaxespad=0., ncol=2)

fname = fpath_s2_scenes.parent.joinpath('Fig_S2-Data-Availability.png')
f.savefig(fname, dpi=200, bbox_inches='tight')

