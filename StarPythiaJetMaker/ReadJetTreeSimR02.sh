#!/bin/bash
# 'ReadJetTreeSimR02.sh'
#
# Use this to run 'ReadJetTreeSim.sh'
# in batch mode

# i/o filepaths
declare -a iFiles
declare -a wFiles
declare -a oFiles

# input paths
iFiles[0]="\"input/NewPi0/pp200py6pt12.goodFiles_waveA.et6100pi0.r02rm1chrg.root\""
iFiles[1]="\"input/NewPi0/pp200py6pt12.goodFiles_waveB.et6100pi0.r02rm1chrg.root\""
iFiles[2]="\"input/NewPi0/pp200py6pt12.goodFiles_waveC.et6100pi0.r02rm1chrg.root\""
iFiles[3]="\"input/NewPi0/pp200py6pt47.goodFiles_wave1.et6100pi0.r02rm1chrg.root\""
iFiles[4]="\"input/NewPi0/pp200py6pt47.goodFiles_wave2.et6100pi0.r02rm1chrg.root\""
iFiles[5]="\"input/NewPi0/pp200py6pt79.goodFiles_wave1.et6100pi0.r02rm1chrg.root\""
iFiles[6]="\"input/NewPi0/pp200py6pt79.goodFiles_wave2.et6100pi0.r02rm1chrg.root\""
iFiles[7]="\"input/NewPi0/pp200py6pt79.goodFiles_wave3.et6100pi0.r02rm1chrg.root\""
iFiles[8]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4a.et6100pi0.r02rm1chrg.root\""
iFiles[9]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4b.et6100pi0.r02rm1chrg.root\""
iFiles[10]="\"input/NewPi0/pp200py6pt79.goodFiles_wave4c.et6100pi0.r02rm1chrg.root\""
iFiles[11]="\"input/NewPi0/pp200py6pt912.goodFiles_all.et6100pi0.r02rm1chrg.root\""

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

# output paths
oFiles[0]="\"pp200py6pt12et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_waveA.et6100r02pi0.d21m8y2022.root\""
oFiles[1]="\"pp200py6pt12et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_waveB.et6100r02pi0.d21m8y2022.root\""
oFiles[2]="\"pp200py6pt12et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_waveC.et6100r02pi0.d21m8y2022.root\""
oFiles[3]="\"pp200py6pt47et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave1.et6100r02pi0.d21m8y2022.root\""
oFiles[4]="\"pp200py6pt47et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave2.et6100r02pi0.d21m8y2022.root\""
oFiles[5]="\"pp200py6pt79et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave1.et6100r02pi0.d21m8y2022.root\""
oFiles[6]="\"pp200py6pt79et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave2.et6100r02pi0.d21m8y2022.root\""
oFiles[7]="\"pp200py6pt79et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave3.et6100r02pi0.d21m8y2022.root\""
oFiles[8]="\"pp200py6pt79et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave4a.et6100r02pi0.d21m8y2022.root\""
oFiles[9]="\"pp200py6pt79et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave4b.et6100r02pi0.d21m8y2022.root\""
oFiles[10]="\"pp200py6pt79et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_wave4c.et6100r02pi0.d21m8y2022.root\""
oFiles[11]="\"pp200py6pt912et6.forEtBinCheck_noWeightOrEtNorm_pTbinOne_all.et6100r02pi0.d21m8y2022.root\""

# loop over files
(( nFile=0 ))
for input in ${iFiles[@]}; do

root -b -l <<EOF
  .x ReadJetTreeSimR02.C(true, $input, ${wFiles[$nFile]}, ${oFiles[$nFile]})
EOF
(( nFile++ ))

# end file loop
done

# delete arrays
unset iFiles
unset oFiles
unset wFiles

# End -------------------------------------------------------------------------
