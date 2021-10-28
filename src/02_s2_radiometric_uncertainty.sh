#!/usr/bin/env bash

dot_safe_dir="./../S2A_MSIL1C_orig"
cd "$dot_safe_dir"

# find subdirectories and loop over them
IFS=""
mapfile -t dirlist < <( find ${dot_safe_dir} -maxdepth 1 -mindepth 1 -type d -printf '%f\n' )

for dir in ${dirlist[@]}; do

	# check if directory ends with .SAFE
	if [[ "$dir" == *.SAFE ]]
	then

 		# the uncertainty files shall be stored in a directory with the same
 		# name but ending with .RUT (radiometric uncertainty toolbox)
 		replace=".RUT"
 		rut_dir=${dir//.SAFE/$replace}
 
		# create the .RUT directory
 		mkdir "$rut_dir"
 
 		# the resulting files (one uncertainty raster per spectral band) shall be
 		# named like the original dataset but without .SAFE and ending with
 		# _rut_<band>.tif
 		replace="_rut_"
 		unc_files_basename=${dir//.SAFE/$replace}

 		# start the S2-Uncertainty toolbox using SNAP's graph processing CLI
 		# each band is processed separately
 		# empty the cache after each band; otherwise there might be out-of-memory errors
 		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b01.tif -Pband_names=B1
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
 		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b02.tif -Pband_names=B2
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b03.tif -Pband_names=B3
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b04.tif -Pband_names=B4
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b05.tif -Pband_names=B5
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b06.tif -Pband_names=B6
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b07.tif -Pband_names=B7
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b08.tif -Pband_names=B8
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b8a.tif -Pband_names=B8A
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b09.tif -Pband_names=B9
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4

		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b10.tif -Pband_names=B10
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b11.tif -Pband_names=B11
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
 
		gpt S2RutOp -Ssource="${dir}" -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"b12.tif -Pband_names=B12
 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4

    fi

done
