'''
Calculates the absolute uncertainty in the phenological metrics including
start, peak and end of season (timing and vegetation index/ parameter values at
these phenological stages).
'''

import glob
import datetime
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Dict

from agrisatpy.io import SatDataHandler
from logger import get_logger
from copy import deepcopy

logger = get_logger('l5_uncertainty')

start_date = datetime.date(2019,1,1)


def calc_l4_uncertainty(
        uncertainty_dir: Path,
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

        # plot as map
        if 'times' in pheno_metric:
            unit = 'days'
        else:
            unit = '-'
        label = f'Absolute Uncertainty (k=1) [{unit}]'
        fig_unc = unc_handler.plot_band(
            band_name,
            colormap='coolwarm',
            colorbar_label=label
        )

        fname_out_fig = out_dir.joinpath(f'{vi_name}_{pheno_metric}_abs-uncertainty.png')
        fig_unc.savefig(fname_out_fig, dpi=300, bbox_inches='tight')
        plt.close(fig_unc)
        fname_out_raster = fname_out_fig.as_posix().replace('.png','.tif')
        unc_handler.write_bands(
            out_file=fname_out_raster,
            band_names=[band_name]
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
    (considering all pixels labelled as these crop types), plus histograms of
    the uncertainty values over all pixels of a crop type
    """

    # find the file with the phenological uncertainty estimates for the selected
    # crop type
    plt.style.use('seaborn-colorblind')

    search_expr = f'{vi_name}_{pheno_metric}*uncertainty.tif'
    unc_file = glob.glob(result_dir.joinpath(search_expr).as_posix())[0]

    # read data, mask out all pixels not belonging to crop selection
    handler = SatDataHandler()
    handler.read_from_bandstack(
        fname_bandstack=unc_file,
        in_file_aoi=shapefile_crops
    )

    # add shapefile data with crop type codes
    unc_band = handler.get_bandnames()[0]
    handler.add_bands_from_vector(
        in_file_vector=shapefile_crops,
        snap_band=unc_band,
        attribute_selection=[column_crop_code],
        blackfill_value=-9999.
    )

    # mask out all other pixels (not having one of the selected crop types)
    handler.mask(
        name_mask_band=column_crop_code,
        mask_values=-9999.,
        bands_to_mask=[unc_band]
    )

    # plot the uncertainty band now masked to the crop selection
    if 'times' in pheno_metric:
        unit = 'days'
    else:
        unit = '-'
    label = f'Absolute Uncertainty (k=1) [{unit}]'
    fig_unc = handler.plot_band(
        band_name=unc_band,
        colormap='coolwarm',
        colorbar_label=label
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

    # plot histograms by crop type (add crop names first)
    gdf['crop'] = gdf.crop_code.apply(
        lambda x, crop_code_mapping=crop_code_mapping: crop_code_mapping[x]
    )
    # histogram of all crops
    gdf[unc_band].hist(by=gdf['crop'], bins=50, sharex=True, sharey=True, density=True)
    plt.suptitle(
        f'{pheno_metric_alias.upper()} derived from {vi_name}:\nRelative Frequencies of Absolute Uncertainty (k=1) Values per Crop Type'
    )
    plt.subplots_adjust(top=0.85)
    plt.savefig(
        f'{fname_out_base}_histogram-uncertainties-all-crops.png',
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()

    # save dataframe to csv for future analysis
    gdf['x'] = gdf.geometry.x
    gdf['y'] = gdf.geometry.y
    gdf.drop('geometry', axis=1, inplace=True)
    gdf.to_csv(f'{fname_out_base}_data.csv', index=False)



def visualize_sample_time_series(
        sample_points_scenarios: Path,
        sample_points_pheno_metrics_reference: Path,
        pheno_metrics_uncertainty_dir: Path,
        vi_name: str
    ):
    """
    Plots the randomly extracted pixel time series, their scenario spread
    and the pheno metrics SOS, POS, EOS alongside their uncertainty for the
    selected pixel.

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
    gdf = SatDataHandler.read_pixels(
        point_features=gdf,
        raster=sample_points_pheno_metrics_reference,
        band_selection=metrics
    )

    # extract the uncertainty of the metrics
    for metric in metrics:
        metric_uncertainty_file = glob.glob(
            pheno_metrics_uncertainty_dir.joinpath(
                f'{vi_name}_{metric}_abs-uncertainty.tif'
            ).as_posix()
        )[0]
        gdf = SatDataHandler.read_pixels(
            point_features=gdf,
            raster=metric_uncertainty_file
        )

    # loop over the point (using their unique id) and visualize the time series
    unique_points = list(gdf.id.unique())

    for point in unique_points:
        
        # get a copy of the slice of the dataframe denoting the point
        point_gdf = gdf[gdf.id == point].copy()
        point_gdf['date'] = pd.to_datetime(point_gdf['date'])

        # plot the time series
        fig = plt.figure()
        ax = fig.add_subplot(111)

        crop_type = point_gdf.crop_type.iloc[0]
        coordinate_str = f'x={np.round(point_gdf.x.iloc[0])}m, y={np.round(point_gdf.y.iloc[0])}m'
        title_str = f'{crop_type} Sample\n{coordinate_str}'

        # plot original data
        # time series
        ax.plot(point_gdf.date, point_gdf[vi_name], 'ro', label='original')
        ax.plot(point_gdf.date, point_gdf[f'{vi_name}_ts_sm'], label='original smoothed')
        # phenological metrics. Need to be convert to dates since pheno metrics are provided
        # as days-of-year (doys)
        for metric in metrics:
            pheno_metric_doy = int(point_gdf[metric].iloc[0])
            pheno_metric_date = start_date + datetime.timedelta(days=pheno_metric_doy)
            ax.vlines(x=pheno_metric_date, ymin=0, ymax=1, label='SOS', linestyles='dashed')

        plt.show()
        



if __name__ == '__main__':

    uncertainty_dir = Path('../S2_TimeSeries_Analysis')
    out_dir = uncertainty_dir.joinpath('Uncertainty_Maps')
    if not out_dir.exists():
        out_dir.mkdir()
    result_dir = out_dir

    out_dir_crops = out_dir.joinpath('selected_crops')
    if not out_dir_crops.exists():
        out_dir_crops.mkdir()

    vi_names = ['NDVI', 'EVI']

    # pheno-metrics to analyze
    pheno_metrics = [
        'sos_times', 'pos_times', 'eos_times', 'sos_values', 'pos_values', 'eos_values'
    ]
    pheno_metrics_aliases = [
        'SOS', 'POS', 'EOS', 'SOS Value', 'POS Value', 'EOS Value'
    ]

    # shapefile with crop type information for the single field parcels
    shapefile_crops = Path('../shp/ZH_Polygons_2019_EPSG32632_selected-crops.shp')
    column_crop_code = 'crop_code'
    column_crop_names = 'crop_type'

    # define mapping of "NUTZUNGSCODES" to crop types
    gdf = gpd.read_file(shapefile_crops)
    crop_code_mapping = dict(list(gdf.groupby([column_crop_code, column_crop_names]).groups))

    for vi_name in vi_names:

        # calc_l4_uncertainty(
        #     uncertainty_dir=uncertainty_dir,
        #     out_dir=out_dir,
        #     vi_name=vi_name
        # )
        #
        # # create maps and histograms of phenometrics
        # for idx, pheno_metric in enumerate(pheno_metrics):
        #     pheno_metric_alias = pheno_metrics_aliases[idx]
        #     get_uncertainty_maps_and_histograms_by_croptype(
        #         result_dir=result_dir,
        #         vi_name=vi_name,
        #         pheno_metric=pheno_metric,
        #         pheno_metric_alias=pheno_metric_alias,
        #         shapefile_crops=shapefile_crops,
        #         column_crop_code=column_crop_code,
        #         crop_code_mapping=crop_code_mapping,
        #         out_dir=out_dir_crops
        #     )

        # visualize the randomly selected pixel time series samples
        vi_dir = uncertainty_dir.joinpath(vi_name)
        # path to pixel samples
        sample_points_scenarios = glob.glob(
            vi_dir.joinpath(f'{vi_name}_*time_series.csv').as_posix()
        )[0]

        # path to reference pheno metric results (calculated on original time series data)
        sample_points_pheno_metrics_reference = vi_dir.joinpath('reference').joinpath('pheno_metrics.tif')

        visualize_sample_time_series(
            sample_points_scenarios=sample_points_scenarios,
            sample_points_pheno_metrics_reference=sample_points_pheno_metrics_reference,
            pheno_metrics_uncertainty_dir=result_dir,
            vi_name=vi_name
        )
        

        


        
    