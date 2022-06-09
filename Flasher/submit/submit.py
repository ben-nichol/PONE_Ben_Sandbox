import os
import sys
import json
import datetime
import numpy as np
from os.path import join, exists, basename


########################################################################
### EXECUTION PARAMETERS
########################################################################
# user
user = os.getenv('USER')

# output directory 
out_dir = '/data/p-one/%s/isotropic' %(user)

# singularity container
singularity = '/data/p-one/icetray_offline_lw.sif'
container = '/home/users/%s/pone_offline/env-shell_Container.sh' %(user)

# memory to request
mem = 8 # GB


########################################################################
### SIMULATION PARAMETERS
########################################################################
# submit 
submit = True

# gcd
gcd = '/home/users/%s/pone_offline/GCD/PONE_Phase1.i3.gz' %(user)

# executable
script = '/home/users/%s/pone_offline/Flasher/submit/clsim_isotropic.py' %(user)

# tag
tag = 'submit-test'

# specify arguments dictionary as used by [script] with --option
# example: --numphotons=100 --> {'numphotons' : 100}
# script must accept --gcd argument (generated automatically)
# script must accept --outfile argument (generated automatically)
args = {'oversize' : 1.0,
        'numevents' : 1,
        'flasherkey' : '5-10',
        'numphotons' : 1e6,
        'fwhm' : 5.0,
        #'wavelength' : 405, # not implemented
        #'optical_medium' : '', # not implemented
       }


########################################################################
### OUTPUT DIRECTORY GENERATION
########################################################################
# create time-sensitive string to avoid overwriting
t_str = datetime.datetime.now().isoformat('_')[:-7].replace(':', '-')

# update log string with arguments
log_str = str(tag) + '_'
for key in sorted(args.keys()):
    log_str += '%s-%s' %(key, args[key])
    log_str += '_'
log_str = log_str[:-1]

# create general output folder
out_folder = join(out_dir, 'sim_' + tag + '_' + t_str)
if not exists(out_folder) and submit:
    print('Creating output directory in {}'.format(out_folder))
    os.makedirs(out_folder)

# create output, log and submit folder for all relevant files
folders = {}
for folder in ['submit', 'job', 'log', 'error', 'out', 'simulation']:
    out_f = join(out_dir, out_folder, folder)
    if not exists(out_f) and submit:
        print('\t-- creating sub-directory in {}'.format(folder))
        os.makedirs(out_f)
        folders[folder] = out_f

out_sub   = folders['submit']
out_job   = folders['job']
out_log   = folders['log']
out_error = folders['error']
out_out   = folders['out']
out_sim   = folders['simulation']

# determine outfile name
out_file = join(out_sim, '%s.i3.bz2' %(tag))
args['outfile'] = out_file

# copy relevant files
os.system('cp %s %s' %(script, join(out_sub, basename(script))) )
os.system('cp %s %s' %(gcd, join(out_sub, basename(gcd))) )

# use copied gcd
args['gcd'] = join(out_sub, basename(gcd))

# save arguments
with open(join(out_sub, 'arguments.txt'), 'w') as f:
     f.write(json.dumps(args))
np.save(join(out_sub, 'arguments.npy'), args)


########################################################################
### GENERATE EXECUTABLE
########################################################################
executable = join(out_sub, 'executable.sh')
with open(executable, 'w') as f:
    # setup environment
    f.write('#!/bin/bash\n')
    f.write('CONTAINER=%s\n' %(container))
    f.write('SCRIPT=%s\n' %(join(out_sub, basename(script))) )
    f.write('\n####\n\n')
    
    # write out arguments
    for key in args:
        f.write('%s="%s"\n' %(key.upper(), args[key]))
    f.write('\n####\n\n')
    
    # python options
    python = ''
    for key in args:
        python += ' --%s ${%s}' %(key, key.upper())
    python += ' '
    
    # execution line
    f.write('bash ${CONTAINER} python ${SCRIPT}%s\n' %(python))


########################################################################
### SUBMIT
########################################################################
# submit file
submit_info = 'executable  = {script} \n\
               +SingularityImage = "{singularity}" \n\
               +TransferOutput = "" \n\
               universe    = vanilla \n\
               request_gpus = 1 \n\
               request_memory = {mem}GB \n\
               log         = {out_log}/{log_str}.log \n\
               output      = {out_out}/{log_str}.out \n\
               error       = {out_err}/{log_str}.err \n\
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
                                  args    = args,
                                 )

if submit:
    # write submit file
    sub_file = '%s.submit' %(join(out_sub, tag))
    with open(sub_file, 'w') as f:
        f.write(submit_info)

    # submit it
    os.system('condor_submit %s' %(sub_file))