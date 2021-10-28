#
#	@purpose	script to run Sen2Cor (Sentinel-2 atmospheric correction)
#				on each L1C realization drawn from the radiometric uncertainty
#				modelling approach.
#
#	@author		Lukas Graf (D-USYS, ETHZ)
#

# source Sen2Cor bashrc (if not add permanently to PATH)
source /mnt/ides/Lukas/software/Sen2Cor-02.09.00-Linux64/L2A_Bashrc

IFS=""
unc_dir="S2A_MSIL1C_20190530T103031_N0207_R108_T32TMT_20190530T123429"

mapfile -t dirlist < <( find ${unc_dir} -maxdepth 2 -mindepth 2 -type d -printf '%f\n' )
# loop over scenarios
counter=1
for dir in ${dirlist[@]}; do
    scenario="${unc_dir}"/"${counter}"/"${dir}"
    output_dir="${unc_dir}"/"${counter}"
    echo "${scenario}"
    counter=$((counter+1))
    # call Sen2Cor
    L2A_Process --resolution 10 --output_dir "${output_dir}" "${scenario}"
done
