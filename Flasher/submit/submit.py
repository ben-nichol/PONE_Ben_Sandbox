
#!/usr/bin/env python

import os
import sys
import datetime
from os.path import join, exists

import numpy as np

########################################################################
### USER / EXECUTION PARAMETERS
########################################################################
# user
user = os.getenv('USER')

# output directory 
out_dir = '/data/p-one/%s/isotropic' %(user)

# executable
executable = '/home/users/%s/pone_offline/Flasher/submit/executable.sh' %(user)

# singularity container
singularity = '/data/p-one/icetray_offline.sif'

# memory to request
mem = 8 # GB

########################################################################
### SIMULATION PARAMETERS
########################################################################
# submit 
submit = True

# gcd
gcd = '/home/users/%s/pone_offline/GCD/PONE_Phase1.i3.gz' %(user)

# tag
tag = 'submit-test'

########################################################################
### OUTPUT DIRECTORY GENERATION
########################################################################
# create time-sensitive string to avoid overwriting
t_str = datetime.datetime.now().isoformat('_')[:-7].replace(':', '-')

# create general output folder
out_folder = join(out_dir, + 'sim_' + tag + '_' + t_str)
if not exists(out_folder) and submit:
    print('Creating output directory in {}'.format(out_folder))
    os.makedirs(out_folder)

# create output, log and submit folder for all relevant files
folders = {}
for folder in ['submit', 'job', 'log', 'error', 'out', 'simulation']:
    out_f = join(out_dir, out_folder, folder)
    if not exists(out_f) and submit:
        print('Creating output directory in {}'.format(out_f))
        os.makedirs(out_f)
        folders[folder] = out_f

out_sub   = folders['submit']
out_job   = folders['job']
out_log   = folders['log']
out_error = folders['error']
out_out   = folders['out']
out_data  = folders['simulation']

########################################################################
### SUBMIT
########################################################################
# submit file
out_file    = '%i' %(1)
args        = '%i' %(1)
env         = "HDF5_USE_FILE_LOCKING='FALSE'"
log_str     = 'job_%s' %(out_file)
submit_info = 'executable  = {script} \n\
               +SingularityImage = {singularity} \n\
               universe    = vanilla \n\
               request_gpus = 1 \n\
               request_memory = {mem}GB \n\
               log         = {out_log}/{log_str}.log \n\
               output      = {out_out}/{log_str}.out \n\
               error       = {out_err}/{log_str}.err \n\
               environment = "{env}" \n\
               arguments   = "{args}" \n\
               requirements = HasSingularity \n\
               transfer_executable = True \n\
               queue 1 \n'.format(
                                  script = executable,
                                  singularity = singularity,
                                  mem     = mem,
                                  out_log = out_log,
                                  out_err = out_error,
                                  out_out = out_out,
                                  log_str = log_str,
                                  env     = env,
                                  args    = args,
                                 )

sys.exit()

if submit:
    # write submit file
    sub_file = '%s/%s.sub' %(out_sub, log_str)
    with open(sub_file, 'w') as f:
        f.write(submit_info)

    # submit it
    os.system('condor_submit {}'.format(sub_file))