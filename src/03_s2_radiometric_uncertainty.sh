#!/usr/bin/env bash

#
#	@author:	Lukas Graf (D-USYS, ETHZ)
#
#	@purpose:	Calls the Sentinel-2 Radiometric Uncertainty Toolbox
#				(S2-RUT) using the Graph Processing Tool (GPT) of
#				SNAP (Sentinel Application Platform). See the README
#				about how to get it running.
#				The script loops over all L1C scenes found and calculates
#				the radiometric uncertainty for each band and uncertainty 
#				contributor separately.
#				The results are stored in a sub-directory named as the original
#				scene but ending with *.RUT instead of *.SAFE
#

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
 
 		# define band names to loop over the single spectral bands
 		bandList=("B1" "B2" "B3" "B4" "B5" "B6" "B7" "B8" "B8A" "B9" "B10" "B11" "B12")
 
		for band in "${bandList[@]}"; do

	 		# start the S2-Uncertainty toolbox using SNAP's graph processing CLI
	 		# each band and uncertainty contributor is processed separately
	 		# empty the cache after each processing step; otherwise there might be out-of-memory errors
	
	 		# 	************************** uncertainty contributors *************************
	 		# instrument noise
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"noise_"$band".tif -Pband_names="$band" -p Instrument_noise.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# out-of-field straylight systematic part
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"stray_sys_"$band".tif -Pband_names="$band" -p OOF_straylight-systematic.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	 
	 		# out-of-field straylight random part
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"stray_rand_"$band".tif -Pband_names="$band" -p OOF_straylight-random.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# crosstalk (not considered any more, see Gorrono et al., 2018 for a discussion)
	 		# gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"x_talk_"$band".tif -Pband_names="$band" -p Crosstalk.properties
	 		# rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# ADC quantisation
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"ADC_"$band".tif -Pband_names="$band" -p ADC_quantisation.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# dark signal stability
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"DS_"$band".tif -Pband_names="$band" -p DS_stability.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# non-linearty and non-uniformity knowledge
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"gamma_"$band".tif -Pband_names="$band" -p Gamma_knowledge.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# diffuser reflectance absolute knowledge
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"diff_abs_"$band".tif -Pband_names="$band" -p Diffuser-absolute_knowledge.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# diffuse reflectance temporal knowledge
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"diff_temp_"$band".tif -Pband_names="$band" -p Diffuser-temporal_knowledge.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# Angular diffuser knowledge - cosine effect
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"diff_cos_"$band".tif -Pband_names="$band" -p Diffuser-cosine_effect.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# strayligth in calibration mode
	 		gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"diff_k_"$band".tif -Pband_names="$band" -p Diffuser-straylight_residual.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
	
			# L1C image quantisation
			gpt S2RutOp -Ssource="${dir}" -x -f GeoTiff -t "${rut_dir}"/"${unc_files_basename}"quant_"$band".tif -Pband_names="$band" -p L1C_image_quantisation.properties
	 		rm -rf /home/"$USER"/.snap/var/cache/s2tbx/l1c-reader/8.0.4
		
		done

    fi

done
