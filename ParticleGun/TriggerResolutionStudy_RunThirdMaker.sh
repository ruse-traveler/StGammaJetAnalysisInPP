#!/bin/bash
# 'TriggerResolutionStudy_RunThirdMaker.sh'
#
# Runs the 'StJetTreeThirdMaker' class over
# the MuDst output of the Particle Gun
# framework.

# i/o parameters
iPath="/star/data01/pwg/dmawxc/Embedding/ParticleGun/output/OneEtaSiteTest"
iSuff="MuDst.root"
iList="pionGridTest_oneEtaSite_pin8pic2.et650h09pi0.d13m5y2022.list"
oPath="./output/OneEtaSiteTest/"
oPref="pionGridTest_oneEtaSite_pin8pic2_file"
oSuff="et650h09.d13m5y2022.root"

# fixed parameters
star=SL14g
nEvt=2000000000



# create list of files
ls $iPath/*.$iSuff > $iList;

# set correct star version
starver $star;

# loop over input
(( nFile=0 ))
for iFile in $(ls -1 $iPath/*$iSuff); do

  # define output file path
  oFile=$oPath$oPref$nFile$oSuff

  # run macro
  root4star -b -q "TriggerResolutionStudy_RunThirdMaker.C($nEvt, \"$oFile\", \"$iFile\")"
  (( nFile ++ ))

done

# End -------------------------------------------------------------------------
