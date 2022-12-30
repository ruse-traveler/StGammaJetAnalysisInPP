#!/bin/bash
# 'TriggerResolutionStudy_Cleanup.sh'
# Derek Anderson
# 09.17.2020
#
# Will delete unnecessary output from
# the trigger resolution study to HPSS.
#
# NOTE: run '*_PushToHPSS.sh'first!

path="/star/data01/pwg/dmawxc/Embedding/ParticleGun"

rm "$path"/recodata/*/*.{event,geant,hist,minimc,runco,starsim,starsim.tags}.root
rm "$path"/simdata/*/*.{fzd,root}
rm "$path"/submitout/*/*.{bfc.log,err,out,sim.log}
rm "$path"/{*.session.xml,sched*.condor,sched*.report,sched*.list,sched*.csh}

# End -------------------------------------------------------------------------
