#!/usr/bin/env bash

#
#	@author		Lukas Graf (D-USYS, ETHZ)
#
#	@purpose	script to run Sen2Cor (Sentinel-2 atmospheric correction)
#				on each the original L1C scenes.
#				The resulting L2A data is stored in the same directory as the
#				input L1C datasets, i.e., in the single scenario folders.
#
#				NOTE: To speed things up, it might be useful to run several
#				instances of Sen2Cor in parallel. In this case, the variable
#				$unc_scenarios (L#28) might have when the scenarios are located
#				at different locations on the file system.
#
#	@requires	Sen2Cor (atmospheric correction and scene classification software).
#				In our case, version 2.09.00 is used.
#

# define directory where the L1C realizations are located
orig_datasets="/nfs/nas12.ethz.ch/fs1202/green_groups_kp_public/Evaluation/Projects/KP0031_lgraf_PhenomEn/Uncertainty/ESCH/scripts_paper_uncertainty/S2A_MSIL1C_orig"

# define directory where Sen2Cor is installed to
sen2cor_install_dir="./../bin/Sen2Cor-02.09.00-Linux64"
source "$sen2cor_install_dir"/L2A_Bashrc

cd "$orig_datasets" 

IFS=""
# find all sub-directories (each scene is a sub-directory)
mapfile -t dataset_list < <( find ${orig_datasets} -maxdepth 1 -mindepth 1 -type d -printf '%f\n' )

# loop over single scenes
for orig_dataset in ${dataset_list[@]}; do
	# check if directory starts with S2
	if [[ "$orig_dataset" == *SAFE ]]
	then
		# call Sen2Cor
		echo L2A_Process --resolution 10 "${orig_dataset}"
		echo Processed "${orig_dataset}"

	fi
done
