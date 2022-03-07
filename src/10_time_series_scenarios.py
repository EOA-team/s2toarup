'''
Created on Mar 7, 2022

@author: graflu
'''

# def _calc_pheno_metrics(xds: xr.Dataset) -> Dict[str, xr.Dataset]:
#     """
#     Calculation of the phenological metrics using ``Phenolopy``.
#
#     Steps:
#         1. Remove outliers using a moving median filter
#         2. Interpolate NaNs linearly
#         3. Smooth data using Savitzky-Golay
#         4. Calculate SOS, POS, EOS using the seasonal amplitude
#
#     :param xds:
#         ``xarray`` dataset containing the vegetation parameter/index
#         stacked over time
#     :return:
#         dictionary with two items: 'pheno_metrics' is a ``xarray``
#         dataset with pixel-based pheno-metrics. 'ds' contains the smoothed
#         per-pixel time series as another ``xarray`` dataset.
#     """
#
#     # remove outliers using median filter
#     ds = remove_outliers(xds, method='median', user_factor=2, z_pval=0.05)
#
#     # interpolate nans linearly
#     ds = interpolate(ds=ds, method='interpolate_na')
#
#     # TODO: smooth data using Whittaker smoother
#     ds = smooth(ds=ds, method='savitsky', window_length=3, polyorder=1)
#
#     # calculate the phenometrics using 20% seasonal amplitude
#     pheno_ds = calc_phenometrics(
#         da=ds['veg_index'],
#         peak_metric='pos',
#         base_metric='bse',
#         method='seasonal_amplitude',
#         factor=0.2,
#         thresh_sides='two_sided'
#     )
#
#     return {'pheno_metrics': pheno_ds, 'ds': ds}
#


# def vegetation_time_series_scenarios(
#         ts_stack_dict: Dict[date, Sentinel2],
#         file_df: pd.DataFrame,
#         n_scenarios: int,
#         out_dir_scenarios: Path,
#         vi_name: str,
#         sample_points: Path
#     ) -> None:
#     """
#     Generates the vegetation time series scenarios applying the uncertainty
#     in the vegetation parameter/indices. The time series are processed in a
#     similar manner TIMESAT does, using the ``Phenolopy`` package.
#
#     To ensure the time series results are meaningful, pixel time series are
#     extracted for a selected number of point features specified
#     by an ESRI shapefile of point features containing information about the
#     underlying crop type. A CSV file containing the extracted pixel time series
#     (original data + smoothed time series from scenarios) is written to
#     ``out_dir_scenarios``.
#
#     :param ts_stack_dict:
#         dictionary with ``SatDataHandler`` instances holding the vegetation parameter
#         as one band and its uncertainty as another band per scene date
#     :param file_df:
#         pandas dataframe with links to the original files (parameter values + uncertainty)
#         for carrying out the reference run
#     :param n_scenarios:
#         number of scenarios to generate
#     :param out_dir_scenarios:
#         directory where to save the results of the scenario runs to
#     :param vi_name:
#         name of the vegetation parameter/index
#     :param dates:
#         list of dates. Must equal the length of the ``ts_stack_list``. Is used
#         to assign a date to each handler to construct the time series.
#     :param sample_points:
#         shapefile with points to use for time series pixels samples. The extracted
#         time series are plotted and saved to a sub-directory called 'pixel_plots'
#     """
#
#     # get coordinates for xarray dataset
#     coords = ts_stack_dict[list(ts_stack_dict.keys())[0]].get_coordinates(vi_name)
#     dates = ts_stack_dict.keys()
#     dates_np = [np.datetime64(x) for x in dates]
#     coords.update({'time': dates_np})
#
#     # add stack atttributes for xarray dataset
#     attrs = {}
#     crs = ts_stack_dict[list(ts_stack_dict.keys())[0]].geo_info.crs
#     attrs['crs'] = crs
#     attrs['transform'] = tuple(
#         ts_stack_dict[list(ts_stack_dict.keys())[0]].geo_info.as_affine()
#     )
#
#     # calculate the x and y array indices required to extract the pixel values
#     # at the sample locations from the array holding the smoothed vegetation
#     gdf_points = gpd.read_file(sample_points)
#     gdf_points['x'] = gdf_points.geometry.x
#     gdf_points['y'] = gdf_points.geometry.y
#
#     # map the coordinates to array indices
#     def find_nearest_array_index(array, value):
#         return np.abs(array - value).argmin()
#     # get column (x) indices
#     gdf_points['col'] = gdf_points['x'].apply(
#         lambda x, coords=coords, find_nearest_array_index=find_nearest_array_index:
#         find_nearest_array_index(coords['x'], x)
#     )
#     # get row (y) indices
#     gdf_points['row'] = gdf_points['y'].apply(
#         lambda y, coords=coords, find_nearest_array_index=find_nearest_array_index:
#         find_nearest_array_index(coords['y'], y)
#     )
#
#     # loop over the scenarios. In each scenario run generate a xarray containing vegetation
#     # index/parameter samples order by date (i.e., the xarray dataset is 3-dimensional: x,y,time)
#     # save the calculated Pheno-metrics by the end of each scenario run to a raster file
#     orig_ts_list = []
#     for scenario in range(n_scenarios):
#
#         logger.info(f'Running scenario ({scenario+1}/{n_scenarios})')
#         # add vegetation samples to 3d numpy array and pass it to xarray dataset
#         sample_list = []
#         gdf_list = []
#         for _date in list(dates):
#             vi_data = ts_stack_dict[_date].get_values(vi_name)
#             vi_unc = ts_stack_dict[_date].get_values(f'{vi_name}_unc')
#             # sample the test points
#             # sample original pixel (SCL classes have been masked) values only in first scenario run
#             if scenario == 0:
#
#                 vegpar_file = file_df[file_df.date == _date]['filename_veg_par_scl'].iloc[0]
#                 gdf_sample_data = RasterCollection.read_pixels(
#                     vector_features=sample_points,
#                     fpath_raster=vegpar_file
#                 )
#                 colnames = list(gdf_sample_data.columns)
#                 if None in colnames:
#                     colnames[colnames.index(None)] = vi_name
#                 gdf_sample_data.columns = colnames
#
#                 vegpar_unc = file_df[file_df.date == _date]['filename_unc_scl'].iloc[0]
#                 gdf_sample_unc = RasterCollection.read_pixels(
#                     vector_features=sample_points,
#                     fpath_raster=vegpar_unc,
#                     band_names_src=[f'{vi_name}_unc']
#                 )
#                 gdf_sample_unc.rename({f'{vi_name}_unc': 'unc'}, axis=1, inplace=True)
#
#                 # join the data frames
#                 gdf = pd.merge(
#                     gdf_sample_data[['id', 'crop_type', 'geometry', vi_name]],
#                     gdf_sample_unc[['id', 'unc']],
#                     on='id'
#                 )
#                 gdf['date'] = _date
#
#                 # join row and column
#                 gdf_joined = pd.merge(
#                     gdf,
#                     gdf_points[['id', 'row', 'col']],
#                     on='id'
#                 )
#                 gdf_list.append(gdf_joined)
#                 # append original raster data to stack once
#                 orig_ts_list.append(vi_data)
#
#             # get samples from uncertainty distribution and add it to the vi_data
#             # (can be done directly because we deal with absolute uncertainties)
#             # TODO: implement also fully correlated case here
#             samples = np.random.normal(
#                 loc=0,
#                 scale=vi_unc,
#                 size=vi_data.shape
#             )
#             samples += vi_data
#             sample_list.append(samples)
#
#         # create the 3d numpy array to pass as xarray dataset
#         stack = {'veg_index': tuple([('y','x','time'), np.dstack(sample_list)])}
#
#         # create stack for the reference run and the pixel time series of the reference
#         # data (do it in the first iteration only to avoid overwritten the data again and
#         # again)
#         if scenario == 0:
#             orig_stack = {'veg_index': tuple([('y','x','time'), np.dstack(orig_ts_list)])}
#             gdf = pd.concat(gdf_list)
#
#         # construct xarray dataset for the scenario run
#         try:
#             xds = xr.Dataset(
#                 stack,
#                 coords=coords,
#                 attrs=attrs
#             )
#             res = _calc_pheno_metrics(xds)
#             pheno_ds = res['pheno_metrics'] # pheno metric results
#             ds = res['ds'] # smoothed time series values
#         except Exception as e:
#             logger.error(f'Error in scenario ({scenario+1}/{n_scenarios}): {e}')
#
#         # get the smoothed time series values from the dataset at the sample locations
#         # unfortunately, there seems to be no ready-to-use solution
#         unique_points = gdf[gdf.date == gdf.date.unique()[0]][['id', 'row', 'col']]
#
#         scenario_col = f'{vi_name}_{scenario+1}'
#         # time series values
#         gdf[scenario_col] = np.empty(gdf.shape[0])
#
#         # loop over sample points and add them as new entries to the dataframe
#         for _, unique_point in unique_points.iterrows():
#             # get time series values
#             ts_values = ds[dict(x=[unique_point.col], y=[unique_point.row])]['veg_index'].data
#             gdf.loc[
#                 (gdf.row == unique_point.row) & (gdf.col == unique_point.col) & (gdf.id == unique_point.id),
#                 scenario_col
#             ] = ts_values[0,0,:]
#
#         # do the same for the reference run
#         if scenario == 0:
#             try:
#                 xds_ref = xr.Dataset(
#                     orig_stack,
#                     coords=coords,
#                     attrs=attrs
#                 )
#                 res_ref = _calc_pheno_metrics(xds_ref)
#                 pheno_ds_ref = res_ref['pheno_metrics'] # pheno metric results
#                 ds_ref = res_ref['ds'] # smoothed time series values
#             except Exception as e:
#                 logger.error(f'Error in reference run: {e}')
#
#             ref_col = f'{vi_name}_ts_sm'
#             gdf[ref_col] = np.empty(gdf.shape[0])
#             # loop over sample points and add them as new entries to the dataframe
#             for _, unique_point in unique_points.iterrows():
#                 ts_values = ds_ref[dict(x=[unique_point.col], y=[unique_point.row])]['veg_index'].data
#                 gdf.loc[
#                     (gdf.row == unique_point.row) & (gdf.col == unique_point.col) & (gdf.id == unique_point.id),
#                     ref_col
#                 ] = ts_values[0,0,:]
#
#         # save to raster files using the SatDataHandler object since it has the full geoinformation
#         # available
#         pheno_handler = deepcopy(ts_stack_dict[list(ts_stack_dict.keys())[0]])
#         pheno_metrics = [
#             'sos_values', 'sos_times', 'pos_values', 'pos_times', 'eos_values', 'eos_times'
#         ]
#
#         out_dir_scenario = out_dir_scenarios.joinpath(str(scenario+1))
#         if not out_dir_scenario.exists():
#             out_dir_scenario.mkdir()
#         out_file = out_dir_scenario.joinpath('pheno_metrics.tif')
#
#         for pheno_metric in pheno_metrics:
#             pheno_handler.add_band(
#                 band_constructor=Band,
#                 band_name=pheno_metric,
#                 values=eval(f'pheno_ds.{pheno_metric}.data'),
#                 geo_info=pheno_handler[vi_name].geo_info
#             )
#
#         # remove "template" files
#         bands_to_drop = [vi_name, f'{vi_name}_unc']
#         for band_to_drop in bands_to_drop:
#             pheno_handler.drop_band(band_to_drop, inplace=True)
#
#         # save to geoTiff
#         pheno_handler.write_bands(out_file)
#
#         # write original data
#         if scenario == 0:
#             pheno_handler = None
#             pheno_handler = deepcopy(ts_stack_dict[list(ts_stack_dict.keys())[0]])
#
#             out_dir_scenario = out_dir_scenarios.joinpath('reference')
#             if not out_dir_scenario.exists():
#                 out_dir_scenario.mkdir()
#             out_file = out_dir_scenario.joinpath('pheno_metrics.tif')
#
#             for pheno_metric in pheno_metrics:
#                 pheno_handler.add_band(
#                     band_name=pheno_metric,
#                     band_data=eval(f'pheno_ds_ref.{pheno_metric}.data'),
#                     snap_band=vi_name
#                 )
#
#             # remove "template" files
#             bands_to_drop = [vi_name, f'{vi_name}_unc']
#             for band_to_drop in bands_to_drop:
#                 pheno_handler.drop_band(band_to_drop)
#
#             # save to geoTiff
#             pheno_handler.write_bands(out_file)
#
#         logger.info(f'Finished scenario ({scenario+1}/{n_scenarios})')
#
#     # save sample points to CSV (can be used for analyzing pixel time series)
#     fname_csv = out_dir_scenarios.joinpath(
#         f'{vi_name}_{sample_points.name}_pixel_time_series.csv'
#     )
#     gdf['date'] = gdf['date'].apply(lambda x: x.strftime('%Y-%m-%d'))
#     gdf['x'] = gdf.geometry.x
#     gdf['y'] = gdf.geometry.y
#     gdf.drop('geometry', axis=1, inplace=True)
#     gdf.to_csv(fname_csv, index=False)
