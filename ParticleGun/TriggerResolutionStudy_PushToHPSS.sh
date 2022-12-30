#!/bin/bash
# 'TriggerResolutionStudy_PushToHPSS.sh'
# Derek Anderson
# 09.17.2020
#
# Will save output from the trigger resolution
# study to HPSS.

name="trigResStudy_subTest2"
date="d17m9y2020"
hpss="/home/dmawxc/TriggerResolutionStudy2020"
path="/star/scratch/DancingWithDerek"

htar -cf "$hpss"/"$name"."$date".tar "$path"/recodata/*/*.{MuDst,event,geant,hist,minimc,runco,starsim,starsim.tags}.root >& hpssRec."$name"."$date".log &
htar -cf "$hpss"/"$name"."$date".tar "$path"/simdata/*/*.{fzd,root} >& hpssSim."$name"."$date".log &
htar -cf "$hpss"/"$name"."$date".tar "$path"/submitout/*/*.{bfc.log,err,out,sim.log} >& hpssOut."$name"."$date".log &
htar -cf "$hpss"/"$name"."$date".tar "$path"/{*.session.xml,sched*.condor,sched*.report,sched*.list,sched*.csh} >& hpssSub."$name"."$date".log &

# End -------------------------------------------------------------------------
