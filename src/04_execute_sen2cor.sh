#
#	@purpose	script to run Sen2Cor (Sentinel-2 atmospheric correction)
#				on each L1C realization drawn from the radiometric uncertainty
#				modelling approach.
#				The resulting L2A data is stored in the same directory as the
#				input L1C datasets, i.e., in the single scenario folders.
#
#	@requires	Having Sen2Cor installed. In our case, version 2.09.00 was
# 				used. The path to the installation directory of the Sen2cor
#				stand-alone version has to be specified.
#				It is recommended to place it in the /../bin directory of
#				the project root.
#
#	@author		Lukas Graf (D-USYS, ETHZ)
#

# define directory where the L1C realizations are located
unc_scenarios="/run/media/graflu/ETH-KP-SSD6/SAT/S2A_MSIL1C_RUT-Scenarios"
# define directory where Sen2Cor is installed to

# source Sen2Cor bashrc (if not added permanently to PATH)
sen2cor_install_dir="./../bin/Sen2Cor-02.09.00-Linux64"
source "$sen2cor_install_dir"/L2A_Bashrc

cd "$unc_scenarios" 

IFS=""
mapfile -t scenario_list < <( find ${unc_scenarios} -maxdepth 1 -mindepth 1 -type d -printf '%f\n' )

# loop over single scenes
for unc_dir in ${scenario_list[@]}; do
	# check if directory starts with S2
	if [[ "$unc_dir" == S2* ]]
	then

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
