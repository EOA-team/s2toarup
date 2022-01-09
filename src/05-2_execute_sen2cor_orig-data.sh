#!/usr/bin/env bash

#
#	This shell script has the same purpose as 05_execute_sen2cor.sh but runs
#	on the original data and not on the L1C scenarios generated. The resulting L2A
#	products serve as reference products.
#

# define directory where the L1C realizations are located
orig_datasets="/mnt/ides/Lukas/software/scripts_paper_uncertainty/S2A_MSIL1C_orig/autumn"

# define directory where Sen2Cor is installed to
sen2cor_install_dir="/home/graflu/Downloads/Sen2Cor-02.09.00-Linux64"
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
		L2A_Process --resolution 10 "${orig_dataset}"
		echo Processed "${orig_dataset}"

	fi
done
