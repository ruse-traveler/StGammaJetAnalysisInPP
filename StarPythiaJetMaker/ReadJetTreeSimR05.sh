#!/bin/bash
# 'ReadJetTreeSimR05.sh'
#
# Use this to run 'ReadJetTreeSim.sh'
# in batch mode

# i/o filepaths
declare -a iFiles
declare -a wFiles
declare -a oFiles

# input paths
iFiles[0]="\"input/NewPi0/pp200py6pt12.goodFiles_waveA.et6100pi0.r05rm1chrg.root\""
iFiles[1]="\"input/NewPi0/pp200py6pt12.goodFiles_waveB.et6100pi0.r05rm1chrg.root\""
iFiles[2]="\"input/NewPi0/pp200py6pt12.goodFiles_waveC.et6100pi0.r05rm1chrg.root\""
iFiles[3]="\"input/NewPi0/pp200py6pt47.goodFiles_wave1.et6100pi0.r05rm1chrg.root\""
iFiles[4]="\"input/NewPi0/pp200py6pt47.goodFiles_wave2.et6100pi0.r05rm1chrg.root\""
iFiles[5]="\"input/NewPi0/pp200py6pt79.goodFiles_wave1.et6100pi0.r05rm1chrg.root\""
iFiles[6]="\"input/NewPi0/pp200py6pt79.goodFiles_wave2a.et6100pi0.r05rm1chrg.root\""
iFiles[7]="\"input/NewPi0/pp200py6pt79.goodFiles_wave2b.et6100pi0.r05rm1chrg.root\""
iFiles[8]="\"input/NewPi0/pp200py6pt79.goodFiles_wave3a.et6100pi0.r05rm1chrg.root\""
iFiles[9]="\"input/NewPi0/pp200py6pt79.goodFiles_wave3b.et6100pi0.r05rm1chrg.root\""
iFiles[10]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4aa.et6100pi0.r05rm1chrg.root\""
iFiles[11]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4ab.et6100pi0.r05rm1chrg.root\""
iFiles[12]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4ba.et6100pi0.r05rm1chrg.root\""
iFiles[13]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4bb.et6100pi0.r05rm1chrg.root\""
iFiles[14]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4ca.et6100pi0.r05rm1chrg.root\""
iFiles[15]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4cb.et6100pi0.r05rm1chrg.root\""
iFiles[16]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4cc.et6100pi0.r05rm1chrg.root\""
iFiles[17]="\"input/NewPi0/pp200py6pt912.goodFiles_all.et6100pi0.r05rm1chrg.root\""

# weight paths
wFiles[0]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[1]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[2]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[3]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[4]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[5]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[6]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[7]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[8]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[9]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[10]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[11]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[12]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[13]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[14]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[15]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[16]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""
wFiles[17]="\"input/Weights/starPythiaWeights_usingQtHistToSmear.et1115vz55tsp008pi0.d21m10y2021.root\""

# output paths
oFiles[0]="\"pp200py6pt12et6.forEtSumCheck_noWeight_pTbinOne_waveA.et8100r05pi0.d15m8y2022.root\""
oFiles[1]="\"pp200py6pt12et6.forEtSumCheck_noWeight_pTbinOne_waveB.et8100r05pi0.d15m8y2022.root\""
oFiles[2]="\"pp200py6pt12et6.forEtSumCheck_noWeight_pTbinOne_waveC.et8100r05pi0.d15m8y2022.root\""
oFiles[3]="\"pp200py6pt47et6.forEtSumCheck_noWeight_pTbinOne_wave1.et8100r05pi0.d15m8y2022.root\""
oFiles[4]="\"pp200py6pt47et6.forEtSumCheck_noWeight_pTbinOne_wave2.et8100r05pi0.d15m8y2022.root\""
oFiles[5]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave1.et8100r05pi0.d15m8y2022.root\""
oFiles[6]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave2a.et8100r05pi0.d15m8y2022.root\""
oFiles[7]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave2b.et8100r05pi0.d15m8y2022.root\""
oFiles[8]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave3a.et8100r05pi0.d15m8y2022.root\""
oFiles[9]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave3b.et8100r05pi0.d15m8y2022.root\""
oFiles[10]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4aa.et8100r05pi0.d15m8y2022.root\""
oFiles[11]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4ab.et8100r05pi0.d15m8y2022.root\""
oFiles[12]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4ba.et8100r05pi0.d15m8y2022.root\""
oFiles[13]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4bb.et8100r05pi0.d15m8y2022.root\""
oFiles[14]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4ca.et8100r05pi0.d15m8y2022.root\""
oFiles[15]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4cb.et8100r05pi0.d15m8y2022.root\""
oFiles[16]="\"pp200py6pt79et6.forEtSumCheck_noWeight_pTbinOne_wave4cc.et8100r05pi0.d15m8y2022.root\""
oFiles[17]="\"pp200py6pt912et6.forEtSumCheck_noWeight_pTbinOne_all.et8100r05pi0.d15m8y2022.root\""

# loop over files
(( nFile=0 ))
for input in ${iFiles[@]}; do

root -b -l <<EOF
  .x ReadJetTreeSimR05.C(true, $input, ${wFiles[$nFile]}, ${oFiles[$nFile]})
EOF
(( nFile++ ))

# end file loop
done

# delete arrays
unset iFiles
unset oFiles
unset wFiles

# End -------------------------------------------------------------------------
