#!/bin/bash
# 'MatchPythiaJets.sh'
#
# Use this to run 'MatchPythiaJets.C' in batch mode

# i/o files
declare -a parFiles
declare -a detFiles
declare -a outFiles

# particle input
parFiles[0]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt510par.default.et920nTrg50Kpt021Kpi0.r05rm1chrg.d8m3y2021.root\""
parFiles[1]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt1015par.default.et920nTrg50Kpt021Kpi0.r05rm1chrg.d8m3y2021.root\""
parFiles[2]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt1520par.default.et920nTrg50Kpt021Kpi0.r05rm1chrg.d8m3y2021.root\""
parFiles[3]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt2025par.default.et920nTrg50Kpt021Kpi0.r05rm1chrg.d8m3y2021.root\""
parFiles[4]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt2530par.default.et920nTrg50Kpt021Kpi0.r05rm1chrg.d8m3y2021.root\""
parFiles[5]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt30par.default.et920nTrg50Kpt021Kpi0.r05rm1chrg.d8m3y2021.root\""

# detector input
detFiles[0]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt510det.default.et920nTrg50Kpt0230pi0.r05rm1chrg.d8m3y2021.root\""
detFiles[1]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt1015det.default.et920nTrg50Kpt0230pi0.r05rm1chrg.d8m3y2021.root\""
detFiles[2]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt1520det.default.et920nTrg50Kpt0230pi0.r05rm1chrg.d8m3y2021.root\""
detFiles[3]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt2025det.default.et920nTrg50Kpt0230pi0.r05rm1chrg.d8m3y2021.root\""
detFiles[4]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt2530det.default.et920nTrg50Kpt0230pi0.r05rm1chrg.d8m3y2021.root\""
detFiles[5]="\"../../../JetReco_pp/PythiaJetMaker/output/March2021/pp200py8pt30det.default.et920nTrg50Kpt0230pi0.r05rm1chrg.d8m3y2021.root\""

# output
outFiles[0]="\"./output/pp200py8pt510.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root\""
outFiles[1]="\"./output/pp200py8pt1015.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root\""
outFiles[2]="\"./output/pp200py8pt1520.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root\""
outFiles[3]="\"./output/pp200py8pt2025.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root\""
outFiles[4]="\"./output/pp200py8pt2530.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root\""
outFiles[5]="\"./output/pp200py8pt30.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root\""

# loop over files
(( iFile=0 ))
for pFile in ${parFiles[@]}; do

# run script
root -b -l <<EOF
  .L MatchPythiaJets.C++
  MatchPythiaJets($pFile, ${detFiles[$iFile]}, ${outFiles[$iFile]}, true)
  .q
EOF
(( iFile++ ))

# end file loop
done

# clean up files
rm *{.d,.so}

# delete arrays
unset parFiles
unset detFiles
unset outFiles
