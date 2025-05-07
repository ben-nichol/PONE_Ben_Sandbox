universe = vanilla
executable = /data/p-one/jstacho/noise/correlatedNoise/noiseGenerator/noiseModule/poneoffline/makeK40Hits.sh

fileNum=$(Process)

arguments = $(fileNum)

logdir=noise_sims
output = /data/p-one/jstacho/condor_logfiles/$(logdir)/out/genModulePONE$(Process).out
log = /data/p-one/jstacho/condor_logfiles/$(logdir)/log/genModulePONE$(Process).log
error = /data/p-one/jstacho/condor_logfiles/$(logdir)/err/genModulePONE$(Process).err

+SingularityImage = "/data/p-one/icetray_offline_lw.sif"

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
requestCpus = 1
requestMemory = 8192
requirements = HasSingularity && CUDACapability && (Machine != "illume-worker-titanxp-09-v2") && (Machine != "illume-worker-titanxp-10-v2")

periodic_remove = (CommittedTime - CommittedSuspensionTime) > 43200

queue 1
