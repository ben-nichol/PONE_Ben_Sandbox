#!/usr/bin/env bash

# environment
eval $(/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh)
metaproject='/home/fhenningsen/osc/py2-ext/build'

########################################################################
### INPUT
########################################################################
# directory
directory="/home/fhenningsen/pone/calibration"

# modification scripts
mod_scripts="$directory/scripts"

# default ice tables
ice_default="$directory/water/default"

# simulation script
simulation_scripts="$directory/simulation"
simulation_script="${simulation_scripts}/ppc.py"

# geometry file
gcd=${1}

# number of runs
N=${2}

# output directory
dir_name=${3}

# output file
file_name=${dir_name}/${4}

# run configuration
oversize=${5}
ice_shift_abs=${6}
ice_shift_sca=${7}
dom_eff_corr=${8}
p0_ang_acc=${9}
p1_ang_acc=${10}
phot=${11}
fldr=${12}
fwid=${13}
wid=${14}
wfla=${15}

# get flasher strings
flashers=$(${metaproject}/env-shell.sh python -c "import sys; sys.path.append('${simulation_scripts}'); from helper import flasher_strings as fs; print(' '.join(fs('${dir_name}')))")

echo "Flashers: $flashers"

########################################################################
### APPLY MODIFICATIONS
########################################################################
# copy all the tables files so that a new submission doesnt change things
# and create a new tmp directory for the current submission
tmp_dir="tmp-a${ice_shift_abs}-b${ice_shift_sca}-de${dom_eff_corr}-p0${p0_ang_acc}-p1${p1_ang_acc}"
# make sure that truth is not deleted after normal jobs finished running and parameters were the same
if [[ ${file_name} == *"TRUTH"* ]]; then
    ice_tables="${directory}/water/tmp/truth_${tmp_dir}"
else
    ice_tables="${directory}/water/tmp/${tmp_dir}"
fi

# and copy all default files
mkdir -p $ice_tables
cp ${ice_default}/* ${ice_tables}/.

# Multiplying the offset to the icemodel parameters:
$metaproject/env-shell.sh python ${mod_scripts}/change_water.py --a_corr=$ice_shift_abs --b_corr=$ice_shift_sca --outfile=${ice_tables}/icemodel.dat

# Replace the overall dom oversizing factor and dom efficiency in cfg.txt:
sed -e "2s/.*/$oversize/" -e "3s/.*/$dom_eff_corr/" ${ice_default}/cfg.txt > ${ice_tables}/cfg.txt

# - taking p0 & p1 as input and produce an as.dat file in the ice directory
# format: dima_from_unified.py <p0> <p1> <out_file>
python ${mod_scripts}/dima_from_unified.py "${p0_ang_acc}" "${p1_ang_acc}" "${ice_tables}/as.dat"

########################################################################
### SIMULATION
########################################################################
#Then we run the simulation for all given positions
for flasher in $flashers;
do
    echo "Flashing ${flasher}..."
    $metaproject/env-shell.sh python ${simulation_script} --gcd-file=${gcd}  --output-i3-file="${file_name}" --number-of-photons=${phot} --number-of-runs=${N} --flasher=${flasher} --ice-tables=${ice_tables} --fldr=${fldr} --fwid=${fwid} --wid=${wid} --wfla=${wfla}
done

# remove the current ice dir
#rm -rf $ice_tables
    
