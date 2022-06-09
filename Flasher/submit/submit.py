
#!/usr/bin/env python

import os
import sys
import time
import datetime
from os.path import join, exists

import numpy as np
from scipy import stats

# quasi-random
import quasi_random_rng

# who
datauser = '/data/user/fhenningsen'

# executable
executable = '/home/fhenningsen/pone/calibration/simulation/executable.sh'
mem        = 2 #GB

# gcd
gcd = "/home/fhenningsen/pone/calibration/gcd/icecube_phyisics-volume-GCD.i3.bz2"

# submit 
submit_data  = True
submit_truth = True

# run parameters
oversize = 16    # DOM oversize factor
N        = 100   # runs
N_truth  = N * 1 # truth runs
N_photon = 1e9   # number of photons per run
dim      = 5     # parameter dimensions
n_params = 1000  # number of (additional) points in n-d grid
n_grid   = 0     # use grid from this index on, if <n_grid> simulated already

# select emission profile
pocam = False
flasher = False if pocam else True

if pocam:
    fldr = 'None'   # default: LED 1 = [0, 359], LED 2 = [360, 719], .. , -1 for cylindrically symmetric
    fwid = '-1'     # default: LED = 9.7, -1 for isotropic emission
    wid  = '5'      # width of the rectangular pulse in ns
    wfla = '405'    # wavelength in nm, if None use wv.dat
    tag = 'pocam'   
    
if flasher:
    fldr = '-1'
    fwid = '9.7'
    wid  = '70'
    wfla = 'None'
    tag = 'flasher'

# flashers
flashers = ['36-40', '36-59',
            '79-47', '80-53',
            '81-51', '82-60',
            '83-45', '84-57', 
            '85-43', '86-50']

# tag for simulation
sim_tag = '%s-5d-nn-n-%i-%i-domR-%i' %(tag, n_grid, n_grid + n_params - 1, oversize)

# truth
tabs  = 1.02
tsca  = 0.95
tdome = 1.05
tp0   = -0.08
tp1   = 0.06

# create time-sensitive folder for output
t_str = datetime.datetime.now().isoformat('_')[:-7].replace(':', '-')
out_folder = join(datauser, 'pone' , 'sim_' + sim_tag + '_' + t_str)
if not exists(out_folder):
    print('Creating output directory in {}'.format(out_folder))
    os.makedirs(out_folder)

# create time-sensitive folder for log files
out_log = join('/scratch/fhenningsen', 'logs', t_str)
if not exists(out_log):
    print('Creating log directory in {}'.format(out_log))
    os.makedirs(out_log)

# create time-sensitive folder for log files
out_sub = join('/scratch/fhenningsen', 'subs', t_str)
if not exists(out_sub):
    print('Creating sub directory in {}'.format(out_sub))
    os.makedirs(out_sub)
    
# save parameters
params = {'gcd'          : gcd,
          'flashers'     : flashers,
          'nph'          : N_photon,
          'scan_n'       : N,
          'truth_n'      : N_truth,
          'truth_abs'    : tabs,
          'truth_sca'    : tsca,
          'truth_domeff' : tdome,
          'truth_p0'     : tp0,
          'truth_p1'     : tp1,
          'fldr'         : fldr,
          'fwid'         : fwid,
          'wid'          : wid,
          'wfla'         : wfla,
          'use_pocam'    : pocam,
         }
np.save(join(out_folder, 'PARAMS.npy'), params)

############## RNG ##############

# create quasi-random input for NN
# quasi-random sampling
z = quasi_random_rng.rnew(dim, n_params + n_grid)

# use normal prior for abs, sca, domeff
# centered at 1, sigma 0.3
# truncate on +- 2 sigma
loc    = 1
scale  = 0.3
clip_a = loc - 2 * scale
clip_b = loc + 2 * scale

# get truncated interval
a, b   = (clip_a - loc) / scale, (clip_b - loc) / scale

# get prior distribution from uniform
# by using inverse CDF
# explanation: CDF transforms distribution to [0,1]
# inverse transforms back
transform = stats.truncnorm(a, b, loc=loc, scale=scale).isf

# replace abs, sca, and domeff priors 
z[:,0] = transform(z[:,0])
z[:,1] = transform(z[:,1])
z[:,2] = transform(z[:,2])

# set p0 range to [-1, 1] uniform
a      = - 1.0
b      =   1.0
z[:,3] = z[:,3] * (b - a) + a

# set p1 range to [-0.2, 0.2] uniform
a      = - 0.2
b      =   0.2
z[:,4] = z[:,4] * (b - a) + a

# optionally start grid from index if other points were simulated before
z = z[n_grid:]

# and save
np.save(join(out_folder, 'GRID.npy'), z)

#################################

print('Number of jobs: %i' %(n_params + 1))

for rng in z:
    
    # get grid pointflashers
    abso, sca, dome, p0, p1, = rng

    # submit file
    out_file    = 'ABS-%.5f_SCA-%.5f_DOME-%.5f_P0-%.5f_P1-%.5f_NPH-%.3e_N-%i' %(abso, sca, dome, p0, p1, N_photon, N)
    args        = "%s %i %s %s %i %.8f %.8f %.8f %.8f %.8f %i %s %s %s %s" %(gcd, N, out_folder, out_file,
                                                                             oversize,
                                                                             abso, sca, dome, p0, p1, N_photon,
                                                                             fldr, fwid, wid, wfla)
    env         = "HDF5_USE_FILE_LOCKING='FALSE'"
    log_str     = 'job_%s' %(out_file)
    submit_info = 'executable  = {script} \n\
        universe    = vanilla \n\
        initialdir = /home/fhenningsen \n\
        request_gpus = 1 \n\
        request_memory = {mem}GB \n\
        log         = {outl}/{logs}.log \n\
        output      = {out}/{logs}.out \n\
        error       = {out}/{logs}.err \n\
        arguments   = "{args}" \n\
        environment = "{env}" \n\
        transfer_executable = True \n\
        queue 1 \n'.format(script = executable,
                    mem    = mem,
                    logs   = log_str,
                    outl   = out_log,
                    out    = out_folder,
                    args   = args,
                    env    = env,
                    )

    if submit_data:

        # write submit file
        sub_file = '%s/%s.sub' %(out_sub, log_str)
        with open(sub_file, 'w') as f:
            f.write(submit_info)

        # submit them
        os.system("condor_submit {}".format(sub_file))

# run truth
out_file    = 'TRUTH_ABS-%.5f_SCA-%.5f_DOME-%.5f_P0-%.5f_P1-%.5f_NPH-%.3e_N-%i' %(tabs, tsca, tdome, tp0, tp1, N_photon, N_truth)
args        = "%s %i %s %s %i %.8f %.8f %.8f %.8f %.8f %i %s %s %s %s" %(gcd, N_truth, out_folder, out_file,
                                                                         oversize,
                                                                         tabs, tsca, tdome, tp0, tp1, N_photon,
                                                                         fldr, fwid, wid, wfla)
env         = "HDF5_USE_FILE_LOCKING='FALSE'"
log_str     = 'job_%s' %(out_file)
submit_info = 'executable  = {script} \n\
    universe    = vanilla \n\
    initialdir = /home/fhenningsen \n\
    request_gpus = 1 \n\
    request_memory = {mem}GB \n\
    log         = {outl}/{logs}.log \n\
    output      = {out}/{logs}.out \n\
    error       = {out}/{logs}.err \n\
    arguments   = "{args}" \n\
    environment = "{env}" \n\
    transfer_executable = True \n\
    queue 1 \n'.format(script = executable,
                mem    = mem,
                logs   = log_str,
                outl   = out_log,
                out    = out_folder,
                args   = args,
                env    = env,
                )

if submit_truth:

    # write submit file
    sub_file = '%s/%s.sub' %(out_sub, log_str)
    with open(sub_file, 'w') as f:
        f.write(submit_info)

    # submit them
    os.system("condor_submit {}".format(sub_file))

print('Number of jobs: %i' %(n_params + 1))
