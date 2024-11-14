#Example of a job to be run in a container

source ~/.bashrc #runs bashrc in container

USERNAME=$1
TASKID=$2

python3 GenerateEvents.py -o /scratch/$USERNAME/gen_${TASKID}.i3 -g /cvmfs/software.pacific-neutrino.org/pone_offline/GCD/PONE_5String.i3.gz 
