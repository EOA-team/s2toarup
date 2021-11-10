#!/usr/bin/env bash

#
#	@author		Lukas Graf (D-USYS, ETHZ)
#
#	@purpose	script to run Sen2Cor (Sentinel-2 atmospheric correction)
#				on each L1C realization drawn from the radiometric uncertainty
#				modelling approach.
#				The resulting L2A data is stored in the same directory as the
#				input L1C datasets, i.e., in the single scenario folders.
#
#				NOTE: To speed things up, it might be useful to run several
#				instances of Sen2Cor in parallel. In this case, the variable
#				$unc_scenarios (L#28) might have when the scenarios are located
#				at different locations on the file system.
#
#	@requires	Sen2Cor (atmospheric correction and scene classification software).
#				In our case, version 2.09.00 is used. An executable (Linux)
#				is placed in the /bin directory of this Git-Repository which
#				should run without any additional installation requirements.
#				If another Sen2Cor version should be used, place it either in the
#				/bin directory and replace the version in the $sen2cor_install_dir
#				variable or place at any location you like and adopt the complete
#				path. Make sure that the L2A_Bashrc file can be found and sourced!
#

# define directory where the L1C realizations are located
unc_scenarios="./../S2A_MSIL1C_RUT-Scenarios"

# define directory where Sen2Cor is installed to (we take the binary delivered
# in the /bin directory per default)
# source Sen2Cor bashrc (if not added permanently to PATH)
sen2cor_install_dir="./../bin/Sen2Cor-02.09.00-Linux64"
source "$sen2cor_install_dir"/L2A_Bashrc

cd "$unc_scenarios" 

IFS=""
# find all sub-directories (each scene is a sub-directory without the .SAFE of the original dataset)
mapfile -t scenario_list < <( find ${unc_scenarios} -maxdepth 1 -mindepth 1 -type d -printf '%f\n' )

# loop over single scenes
for unc_dir in ${scenario_list[@]}; do
	# check if directory starts with S2
	if [[ "$unc_dir" == S2* ]]
	then

		# find the scenarios of each scene (each scenario run is a sub-directory labelled starting from
		# 1, 2, ..., n_scenarios)
		mapfile -t dirlist < <( find ${unc_scenarios}/${unc_dir} -maxdepth 2 -mindepth 2 -type d -printf '%f\n' )

		# loop over scenarios
		counter=1
		for dir in ${dirlist[@]}; do
			if [[ "$dir" == S2*_MSIL1C_* ]]
			then
			    scenario="${unc_dir}"/"${counter}"/"${dir}"
			    output_dir="${unc_dir}"/"${counter}"
			    counter=$((counter+1))
			    # call Sen2Cor
			    L2A_Process --resolution 10 --output_dir "${output_dir}" "${scenario}"
			    echo Processed "${scenario}"
			fi
		done
	fi
done
