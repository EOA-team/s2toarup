'''
Calculates the absolute uncertainty in the phenological metrics including
start, peak and end of season (timing and vegetation index/ parameter values at
these phenological stages).

In addition, the module plots the time series of selected pixels extracted in the
previous step.
'''

import glob
import datetime
import matplotlib
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from pathlib import Path
from typing import Dict
from typing import Union
from typing import Any
from typing import Optional

from eodal.core.band import Band
from eodal.core.raster import RasterCollection
from utils.logger import get_logger
from copy import deepcopy

logger = get_logger('l5_uncertainty')

# define start date for converting DOYs (day of year) to date
start_date = datetime.date(2019,1,1)

def get_stats(
        gdf: Union[pd.DataFrame, gpd.GeoDataFrame],
        col_name: str,
        unit: str,
        precision: Optional[Union[int, float]] = 1
    ) -> Dict[str, Any]:
    """
    Calculates base statistics of a dataframe column and returns a summary
    string that can be used for plotting.

    :param gdf:
        dataframe to analyze
    :param col_name:
        name of the column to analyze
    :param unit:
        (phyiscal) unit of the column (used in the string)
    :param precision:
        number of decimal places to round the data to
    :return:
        dict with statistics and summary text string
    """

    stats = {
        'median': gdf[col_name].median(),
        'mean':  gdf[col_name].mean(),
        'q05': gdf[col_name].quantile(0.05),
        'q95': gdf[col_name].quantile(0.95)
    }

    textstr = '\n'.join((
        r'$\mu=%.1f$' % (np.round(stats['mean'], precision), ) + f' {unit}',
        r'$\mathrm{median}=%.0f$' % (stats['median'], ) + f' {unit}',
        r'$5\%$'+r' quantile=%.1f' % (np.round(stats['q05'], precision), ) + f' {unit}',
        r'$95\%$'+r' quantile=%.1f' % (np.round(stats['q95'], precision), ) + f' {unit}'
    ))
    stats['textstr'] = textstr

    return stats

def calc_l4_uncertainty(
        uncertainty_dir: Path,
        shp_crops: Path,
        out_dir: Path,
        vi_name: str
    ) -> None:
    """
    Assesses the uncertainty in the phenological metrics derived in the
    previous step. Reports absolute uncertainties, only.

    :param uncertainty_dir:
        directory containing the phenometric results from the time series
        scenarios of the vegetation indices/ parameters analyzed
    :param out_dir:
        directory where to save the results of the uncertainty analysis to
    :param vi_name:
        name of the vegetation index/parameter to analyze
    """

    plt.style.use('default')
    # search the scenarios, organized by vegetation index/ parameter
    vi_uncertainty_dir = uncertainty_dir

    # search scenarios
    vi_search_expr = vi_uncertainty_dir.joinpath('*/pheno_metrics.tif').as_posix()
    scenarios = glob.glob(vi_search_expr)

    # pheno-metrics available
    handler_list = []
    for idx, scenario in enumerate(scenarios):
        if 'reference' in scenario:
            continue
        handler = RasterCollection.from_multi_band_raster(
            fpath_raster=Path(scenario),
            vector_features=shp_crops
        )
        handler_list.append(handler)
        handler = None
        logger.info(f'Reading scenario {idx+1}/{len(scenarios)} ({scenario})')
        # debug
        # break

    # calculate the absolute uncertainty for each phenological metric
    pheno_metric_keys = handler_list[0].band_names
    pheno_metric_keys.append('length_of_season')
    pheno_metrics = dict.fromkeys(pheno_metric_keys)
    
    for pheno_metric in pheno_metrics:

        if pheno_metric == 'SCL' or pheno_metric == 'crop_code':
            continue
        
        # get bandstack of all scenarios of the pheno metric to calculate the standard
        # deviation (=standard uncertainty)
        if pheno_metric != 'length_of_season':
            stack_list = [x[pheno_metric].values for x in handler_list]
        else:
            # length of season is the difference between EOS and SOS in days
            stack_list = [x['eos_times'].values - x['sos_times'].values for x in handler_list]
        stack_array = np.stack(stack_list)
        standard_unc = np.nanstd(stack_array, axis=0)
        # calculate mean of scenarios
        scenario_mean = np.nanmean(stack_array, axis=0)

        # save to raster files and create a preview plot
        unc_handler = deepcopy(handler_list[0])
        snap_band = unc_handler.band_names[0]
        band_names = [f'{pheno_metric} Uncertainty', f'{pheno_metric} Mean']
        band_data = [standard_unc, scenario_mean]

        for bdx, band_name in enumerate(band_names):
            unc_handler.add_band(
                band_constructor=Band,
                band_name=band_name,
                values=band_data[bdx],
                geo_info=unc_handler[snap_band].geo_info
            )

        # plot as map
        if 'times' in pheno_metric:
            unit = 'days'
            vmin= 0
            vmax = 50
        else:
            unit = '-'
            vmin = 0
            vmax = None
        label = f'Absolute Uncertainty (k=1) [{unit}]'

        fig_unc = unc_handler.plot_band(
            band_names[0],
            colormap='Oranges',
            colorbar_label=label,
            vmin=vmin,
            vmax=vmax,
            fontsize=20
        )

        fname_out_fig = out_dir.joinpath(f'{vi_name}_{pheno_metric}_abs-uncertainty.png')
        fig_unc.savefig(fname_out_fig, dpi=300, bbox_inches='tight')
        plt.close(fig_unc)

        fname_out_raster = fname_out_fig.as_posix().replace('.png','.tif')
        unc_handler.to_rasterio(
            fpath_raster=fname_out_raster,
            band_selection=band_names
        )

def get_uncertainty_maps_and_histograms_by_croptype(
        result_dir: Path,
        vi_name: str,
        pheno_metric: str,
        pheno_metric_alias: str,
        shapefile_crops: Path,
        column_crop_code: str,
        crop_code_mapping: Dict[int, str],
        out_dir: Path
    ):
    """
    Get uncertainty maps of a phenological metric for selected crop types
    (considering all pixels labeled as this crop type), plus histograms of
    the uncertainty values over all pixels of a crop type.

    :param result_dir:
        directory where the results of the uncertainty analysis (in the
        phenological metrics) are stored
    :param vi_name:
        name of the vegetation index/ parameter to process
    :param pheno_metric:
        name of the phenometric to process (e.g., 'sos_times')
    :param pheno_metric_alias:
        alias name of the phenometric (used for adding text to the plots in
        a more readable manner)
    :param shapefile_crops:
        ESRI shapefile with vector features containing crop type information
    :param column_crop_code:
        attribute (column name in the resulting GeoDataFrame) in the shapefile
        denoting the crop codes as integer values
    :param out_dir:
        directory where to save the results to
    """

    # find the file with the phenological uncertainty estimates for the selected
    # crop type
    plt.style.use('seaborn-colorblind')

    search_expr = f'{vi_name}_{pheno_metric}*uncertainty.tif'
    unc_file = glob.glob(result_dir.joinpath(search_expr).as_posix())[0]

    # read data, mask out all pixels not belonging to crop selection
    handler = RasterCollection.from_multi_band_raster(
        fpath_raster=unc_file
    )

    # add shapefile data with crop type codes
    band_names = handler.band_names
    unc_band = band_names[0]
    mean_band = band_names[1]
    handler.add_band(
        band_constructor=Band.from_vector,
        vector_features=shapefile_crops,
        geo_info=handler[unc_band].geo_info,
        band_name_src=column_crop_code,
        band_name_dst=column_crop_code,
        snap_bounds=handler[unc_band].bounds,
        nodata_dst=-9999.
    )

    # mask out all other pixels (not having one of the selected crop types)
    handler.mask(
        mask=column_crop_code,
        mask_values=-9999.,
        bands_to_mask=[unc_band],
        inplace=True
    )

    # plot the uncertainty band now masked to the crop selection
    if 'times' in pheno_metric:
        unit = 'days'
        vmin= 0
        vmax = 50
    else:
        unit = '-'
        vmin = 0
        vmax = None
    label = f'Absolute Uncertainty (k=1) [{unit}]'
    fig_unc = handler.plot_band(
        band_name=unc_band,
        colormap='Oranges',
        colorbar_label=label,
        vmin=vmin,
        vmax=vmax
    )
    fname_out_base = out_dir.joinpath(
        f'{vi_name}_{pheno_metric}_abs-uncertainty_{shapefile_crops.name.split(".")[0]}'
    ).as_posix()
    fig_unc.savefig(f'{fname_out_base}.png', dpi=300, bbox_inches='tight')
    plt.close(fig_unc)

    # convert to dataframe to compute histograms per crop types
    gdf = handler.to_dataframe()

    # drop nan's
    gdf = gdf[~np.isnan(gdf[unc_band])]
    nodata_idxs = gdf[gdf.crop_code == -9999.].index
    gdf.drop(nodata_idxs, inplace=True)

    # plot histograms by crop type (add crop names first)
    gdf['crop'] = gdf.crop_code.apply(
        lambda x, crop_code_mapping=crop_code_mapping: crop_code_mapping[x]
    )

    # histogram of all crops
    crop_unc_list = []
    crop_mean_list = []
    for crop in list(gdf.crop.unique()):

        crop_gdf = gdf[gdf.crop == crop].copy()

        # some statistics
        stats_unc = get_stats(crop_gdf, unc_band, unit=unit)
        stats_unc['crop'] = crop
        stats_unc['vi_name'] = vi_name
        stats_unc['metric'] = pheno_metric
        stats_unc['metric_alias'] = pheno_metric_alias
        crop_unc_list.append(stats_unc)

        textstr_unc = stats_unc['textstr']
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # plot histogram of pixel uncertainties and scenario means
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(16,8))

        ax1.hist(crop_gdf[unc_band], bins=50, density=True)
        ax1.set_ylabel('Relative Frequency', fontsize=24)
        ax1.set_xlabel(f'Absolute Uncertainty [{unit}] (k=1)', fontsize=24)
        ax1.text(0.4, 0.95, textstr_unc, transform=ax1.transAxes, fontsize=16,
                 verticalalignment='top', bbox=props)

        stats_mean = get_stats(crop_gdf, mean_band, unit=unit)
        stats_mean['crop'] = crop
        stats_mean['vi_name'] = vi_name
        stats_mean['metric'] = pheno_metric
        stats_mean['metric_alias'] = pheno_metric_alias
        crop_mean_list.append(stats_mean)
        
        textstr_mean = stats_mean['textstr']
        ax2.hist(crop_gdf[mean_band], bins=50, density=True)
        ax2.set_xlabel(f'Pixel Scenario Mean [{unit}]', fontsize=24)
        ax2.text(0.05, 0.95, textstr_mean, transform=ax2.transAxes, fontsize=16,
                 verticalalignment='top', bbox=props)

        fig.suptitle(
            f'{crop} {pheno_metric_alias.upper()} derived from {vi_name}\n'\
            f'Number of Scenarios=1000; Number of Pixels={crop_gdf.shape[0]}',
            fontsize=26
        )

        fname_out = out_dir.joinpath(
            f'{vi_name}_{crop}_{pheno_metric}_pixel-histograms.png'
        )
        fig.savefig(fname_out, dpi=300, bbox_inches='tight')
        plt.close(fig)

    # save crop statistics to csv
    crop_unc_stats = pd.DataFrame(crop_unc_list)
    crop_unc_stats.drop('textstr', axis=1, inplace=True)
    crop_unc_stats.to_csv(f'{fname_out_base}_uncertainty-crops-stats.csv', index=False)

    crop_mean_stats = pd.DataFrame(crop_mean_list)
    crop_mean_stats.drop('textstr', axis=1, inplace=True)
    crop_mean_stats.to_csv(f'{fname_out_base}_scenario-means-crops-stats.csv', index=False)

    # save dataframe to csv for future analysis
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    gdf.drop('geometry', axis=1, inplace=True)
    gdf.to_csv(f'{fname_out_base}_data.csv', index=False)

def visualize_sample_time_series(
        sample_points_scenarios: Path,
        sample_points_pheno_metrics_reference: Path,
        pheno_metrics_uncertainty_dir: Path,
        vi_name: str,
        ymin: int,
        ymax: int,
        out_dir: Path
    ):
    """
    Plots the randomly extracted pixel time series, their scenario spread
    and the pheno metrics SOS, POS, EOS alongside their uncertainty for the
    selected pixel.

    :param sample_points_scenarios:
        path to the csv file containing the different time series realizations
        based on the underlying uncertainty in the vegetation index/parameter
    :param sample_points_pheno_metrics:
        path to the calculated pheno metrics for the original (reference) time
        series data
    :param pheno_metrics_uncertainty_dir:
        directory where the uncertainty results of the single pheno metrics
        are stored as geoTiff files
    :param vi_name:
        name of the vegetation index/parameter to analyze
    :param ymin:
        minimum y (index/parameter) value to use for displaying in plot
    :param ymax:
        maximum y (index/parameter) value to use for displaying in plot
    :param out_dir:
        directory where to save the plots to
    """

    # read extracted pixel samples as dataframe
    df = pd.read_csv(sample_points_scenarios)
    # drop nans
    df.dropna(inplace=True)
    # convert back to geodataframe
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.x, df.y)
    )
    # set CRS (UTM zone 32N)
    gdf.set_crs(epsg=32632, inplace=True)

    # extract the timing of SOS, POS and EOS from the reference run
    metrics = ['sos_times', 'pos_times', 'eos_times']
    gdf = RasterCollection.read_pixels(
        vector_features=gdf,
        fpath_raster=sample_points_pheno_metrics_reference
    )

    # extract the uncertainty and scenario means of the metrics
    for metric in metrics:
        metric_uncertainty_file = glob.glob(
            pheno_metrics_uncertainty_dir.joinpath(
                f'{vi_name}_{metric}_abs-uncertainty.tif'
            ).as_posix()
        )[0]
        gdf = RasterCollection.read_pixels(
            vector_features=gdf,
            fpath_raster=metric_uncertainty_file
        )

    # loop over the point (using their unique id) and visualize the time series
    unique_points = list(gdf.id.unique())

    for point in unique_points:
        
        # get a copy of the slice of the dataframe denoting the point
        point_gdf = gdf[gdf.id == point].copy()
        point_gdf['date'] = pd.to_datetime(point_gdf['date'])

        # plot the time series
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)

        crop_type = point_gdf.crop_type.iloc[0]
        coordinate_str = f'x={np.round(point_gdf.x.iloc[0])}m, y={np.round(point_gdf.y.iloc[0])}m'
        title_str = f'{crop_type} Sample\n{coordinate_str}'

        # phenological metrics. Need to be convert to dates since pheno metrics are provided
        # as days-of-year (doys)
        metric_color = 'grey'
        for metric in metrics:
            pheno_metric_doy = int(point_gdf[metric].iloc[0])
            pheno_metric_date = start_date + datetime.timedelta(days=pheno_metric_doy)
            ax.vlines(x=pheno_metric_date, ymin=ymin, ymax=ymax,
                      linestyles='dashed', color=metric_color, linewidth=3)
            metric_text = metric.replace('_times', '').upper()
            if metric_text == 'EOS':
                time_delta = 3
            else:
                time_delta = -11
            ax.text(pheno_metric_date+datetime.timedelta(time_delta),
                    ymin+0.1, metric_text, fontsize=20)

            # add uncertainty range around the metric (expressed in days) computed from
            # the mean of the scenarios
            # everything smaller than .5 days will be rounded to the next smaller int
            # and vice versa
            unc = int(np.round(point_gdf[f'{metric} Uncertainty'].iloc[0]))
            scenario_mean = int(np.round(point_gdf[f'{metric} Mean'].iloc[0]))
            scenario_mean_date = start_date + datetime.timedelta(days=scenario_mean)
            ax.vlines(x=pheno_metric_date, ymin=ymin, ymax=ymax,
                      linestyles='dotted', color=metric_color, linewidth=3)

            unc_x1 = scenario_mean_date - datetime.timedelta(unc)

            ax.add_patch(
                patches.Rectangle(
                    (unc_x1, ymin),
                    datetime.timedelta(2*unc),
                    ymax,
                    edgecolor=metric_color,
                    fill=False
                )
            )

        # get a climpse of the time series scenarios and visualize their spread
        mask = point_gdf.columns.str.contains(vi_name +'_\d')
        scenarios = point_gdf.loc[:,mask].copy()
        scenarios['min'] = scenarios.min(axis=1)
        scenarios['max'] = scenarios.max(axis=1)
        scenarios['std_plus'] = scenarios.mean(axis=1) + scenarios.std(axis=1)
        scenarios['std_minus'] = scenarios.mean(axis=1) - scenarios.std(axis=1)

        ax.fill_between(x=point_gdf.date, y1=scenarios['min'], y2=scenarios['max'],
                        color='orange', alpha=0.4, label='Min-Max Scenarios')

        ax.fill_between(x=point_gdf.date, y1=scenarios['std_minus'], y2=scenarios['std_plus'],
                        color='red', alpha=0.45, label=r'$\pm$ 1 Stddev Scenarios')

        # plot time series original data
        ax.plot(point_gdf.date, point_gdf[vi_name], 'bo', label='original')
        ax.plot(point_gdf.date, point_gdf[f'{vi_name}_ts_sm'], color='blue',
                label='original smoothed')

        # legend below the plot
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.3), ncol=4, fontsize=20)
        ax.set_ylabel(f'{vi_name} [-]', fontsize=24)
        ax.set_title(title_str, fontsize=24)
        ax.set_ylim(ymin, ymax)
        plt.xticks(rotation=45)

        # save figure
        fname = out_dir.joinpath(
            f'{vi_name}_{crop_type}_{coordinate_str.replace("=","").replace(", ","_")}.png'
        )
        fig.savefig(fname, dpi=300, bbox_inches='tight')
        plt.close(fig)

if __name__ == '__main__':

    vi_names = ['NDVI', 'EVI', 'GLAI']
    ymins = {'NDVI': 0, 'EVI': 0, 'GLAI': 0}
    ymaxs = {'NDVI': 1, 'EVI': 1, 'GLAI': 7}

    # pheno-metrics to analyze
    pheno_metrics = [
        'length_of_season', 'sos_times', 'pos_times', 'eos_times'
    ]
    pheno_metrics_aliases = [
        'LOS', 'SOS', 'POS', 'EOS'
    ]

    # shapefile with crop type information for the single field parcels
    shapefile_crops = Path('../../shp/ZH_Polygons_2019_EPSG32632_selected-crops_buffered.shp')
    column_crop_code = 'crop_code'
    column_crop_names = 'crop_type'

    # define mapping of "crop_codes" to crop types
    gdf = gpd.read_file(shapefile_crops)
    crop_code_mapping = dict(list(gdf.groupby([column_crop_code, column_crop_names]).groups))

    # two ways of inter-scene correlation
    corr_types = ['uncorrelated', 'fully_correlated']

    for corr_type in corr_types:
        for vi_name in vi_names:
    
            uncertainty_dir = Path(f'../../S2_TimeSeries_Analysis/{vi_name}/{corr_type}')
            out_dir = uncertainty_dir.joinpath('Uncertainty_Maps')
            if not out_dir.exists():
                out_dir.mkdir()
            result_dir = out_dir
    
            out_dir_crops = out_dir.joinpath('selected_crops')
            if not out_dir_crops.exists():
                out_dir_crops.mkdir()
    
            out_dir_ts_plots = out_dir.joinpath('pixel_time_series')
            if not out_dir_ts_plots.exists():
                out_dir_ts_plots.mkdir()
    
    
            calc_l4_uncertainty(
                uncertainty_dir=uncertainty_dir,
                out_dir=out_dir,
                vi_name=vi_name,
                shp_crops=shapefile_crops
            )

    # change plot style here to ggplot (therefore, use two different loops)
    plt.style.use('ggplot')
    matplotlib.rc('xtick', labelsize=18) 
    matplotlib.rc('ytick', labelsize=18) 

    for corr_type in corr_types:
        for vi_name in vi_names:
    
            uncertainty_dir = Path(f'../../S2_TimeSeries_Analysis/{vi_name}/{corr_type}')
            out_dir = uncertainty_dir.joinpath('Uncertainty_Maps')
            result_dir = out_dir
            out_dir_crops = out_dir.joinpath('selected_crops')
            out_dir_ts_plots = out_dir.joinpath('pixel_time_series')
    
            # create maps and histograms of phenometrics
            for idx, pheno_metric in enumerate(pheno_metrics):
                pheno_metric_alias = pheno_metrics_aliases[idx]
                get_uncertainty_maps_and_histograms_by_croptype(
                    result_dir=result_dir,
                    vi_name=vi_name,
                    pheno_metric=pheno_metric,
                    pheno_metric_alias=pheno_metric_alias,
                    shapefile_crops=shapefile_crops,
                    column_crop_code=column_crop_code,
                    crop_code_mapping=crop_code_mapping,
                    out_dir=out_dir_crops
                )
    
            # visualize the randomly selected pixel time series samples
            vi_dir = uncertainty_dir
            # path to pixel samples
            sample_points_scenarios = glob.glob(
                vi_dir.joinpath(f'{vi_name}_*time_series.csv').as_posix()
            )[0]
            
            # path to reference pheno metric results (calculated on original time series data)
            sample_points_pheno_metrics_reference = vi_dir.joinpath(
                'reference'
                ).joinpath(
                    'pheno_metrics.tif'
                )
            
            out_dir_ts_plots_vi = out_dir_ts_plots

            visualize_sample_time_series(
                sample_points_scenarios=sample_points_scenarios,
                sample_points_pheno_metrics_reference=sample_points_pheno_metrics_reference,
                pheno_metrics_uncertainty_dir=result_dir,
                vi_name=vi_name,
                ymin=ymins[vi_name],
                ymax=ymaxs[vi_name],
                out_dir=out_dir_ts_plots_vi
            )
