import os
import json
import datetime
import itertools
import numpy as np
from os.path import join, exists, basename


###############################################################################
### EXECUTION PARAMETERS
###############################################################################
# user
user = os.getenv('USER')

# output directory 
out_dir = '/data/p-one/%s/isotropic' %(user)

# singularity container
singularity = '/data/p-one/icetray_offline_lw.sif'
container = '/home/users/%s/pone_offline/env-shell_Container.sh' %(user)

# memory to request
mem = 8 # GB


###############################################################################
### SIMULATION PARAMETERS
###############################################################################
# submit 
submit = True

# gcd
gcd = '/home/users/%s/pone_offline/GCD/PONE_10String.i3.gz' %(user)

# executable
script = '/home/users/%s/pone_offline/Flasher/Isotropic/isotropic.py' %(user)

# angular acceptance
angular = '/home/users/%s/pone_offline/Flasher/resources/as.uniform' %(user)

# tag
tag = 'isotropic-angular-acc-test'

# specify arguments dictionary as used by [script] with --option
# example: --numphotons=100 --> {'numphotons' : 100}
# script must accept --gcd argument (generated automatically)
# script must accept --outfile argument (generated automatically)
# use lists; all permutations will be simulated
args = {'oversize'           : [1.0],
        'num-events'         : [100],
        'flasher-key'        : ['1-10', '5-10', '10-10'],
        'num-photons'        : [int(1e10)],
        'fwhm'               : [5.00],
        'detect-emitter'     : [int(True)],
        'wavelength'         : [405],
        #'optical_medium' : [''], # not implemented
       }

# permute values
iter_keys = [k for k in sorted(args.keys())]
iter_list = [args[k] for k in iter_keys]
iters = list(itertools.product(*iter_list))
print('Number of simulations: %i' %(len(iters)))
print('Permutations: %s' %(iters))


###############################################################################
### OUTPUT DIRECTORY GENERATION
###############################################################################
# create time-sensitive string to avoid overwriting
t_str = datetime.datetime.now().isoformat('_')[:-7].replace(':', '-')

# create general output folder
out_folder = join(out_dir, 'sim_' + tag + '_' + t_str)
if not exists(out_folder):
    print('Creating output directory in {}'.format(out_folder))
    os.makedirs(out_folder)

# create output, log and submit folder for all relevant files
folders = {}
for folder in ['submit', 'job', 'log', 'error', 'out', 'simulation']:
    out_f = join(out_dir, out_folder, folder)
    folders[folder] = out_f
    if not exists(out_f):
        print('\t-- creating sub-directory in {}'.format(folder))
        os.makedirs(out_f)

out_sub   = folders['submit']
out_job   = folders['job']
out_log   = folders['log']
out_error = folders['error']
out_out   = folders['out']
out_sim   = folders['simulation']

# copy relevant files
os.system('cp %s %s' %(script, join(out_sub, basename(script))) )
os.system('cp %s %s' %(gcd, join(out_sub, basename(gcd))) )
os.system('cp %s %s' %(angular, join(out_sub, basename(angular))) )


# save arguments
with open(join(out_sub, 'arguments.txt'), 'w') as f:
     f.write(json.dumps(args))
np.save(join(out_sub, 'arguments.npy'), args)


###############################################################################
### ITERATE PERMUTATIONS
###############################################################################
for i, tup in enumerate(iters):
    
    
    ###########################################################################
    ### OUT FILE TAGGING
    # create mock copy of the arguments dict for use at this iteration
    args_temp = {}      
    
    # create individual log string
    log_str = str(tag) + '_'
    for j, key in enumerate(iter_keys):
        val = tup[j]
        log_str += '%s-%s' %(key, val)
        log_str += '_'
        args_temp[key] = val
    log_str = log_str[:-1]
    
    # determine outfile name
    out_file = join(out_sim, '%s.i3.bz2' %(log_str))
    args_temp['out-file'] = out_file
    
    # use copied gcd
    args_temp['gcd'] = join(out_sub, basename(gcd))
    
    # use copied angular acceptance
    args_temp['angular-acceptance'] = join(out_sub, basename(angular))
    
    # save arguments
    with open(join(out_sub, '%s_arguments.txt' %(log_str)), 'w') as f:
         f.write(json.dumps(args_temp))
    np.save(join(out_sub, '%s_arguments.npy' %(log_str)), args_temp)
       
    
    ###########################################################################
    ### GENERATE EXECUTABLE
    
    executable = join(out_sub, '%s_executable.sh' %(log_str))
    with open(executable, 'w') as f:
        # setup environment
        f.write('#!/bin/bash\n')
        f.write('CONTAINER=%s\n' %(container))
        f.write('SCRIPT=%s\n' %(join(out_sub, basename(script))) )
        f.write('\n####\n\n')
        
        # write out arguments
        for key in args_temp:
            var = key.replace('-', '')
            if type(args_temp[key]) == str:
                f.write('%s="%s"\n' %(var.upper(), args_temp[key]))
            if type(args_temp[key]) == int:
                f.write('%s="%i"\n' %(var.upper(), args_temp[key]))
            if type(args_temp[key]) == float:
                f.write('%s="%.5f"\n' %(var.upper(), args_temp[key]))
        f.write('\n####\n\n')
        
        # python options
        python = ''
        for key in args_temp:
            var = key.replace('-', '')
            python += ' --%s ${%s}' %(key, var.upper())
        python += ' '
        
        # execution line
        f.write('bash ${CONTAINER} python ${SCRIPT}%s\n' %(python))
        
    
    ###########################################################################
    ### SUBMIT
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
                                     )
    
    if submit:
        # write submit file
        sub_file = '%s.submit' %(join(out_sub, log_str))
        with open(sub_file, 'w') as f:
            f.write(submit_info)
    
        # submit it
        os.system('condor_submit %s' %(sub_file))
