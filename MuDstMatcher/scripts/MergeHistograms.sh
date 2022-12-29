#!/bin/bash
# 'MergeHistograms.sh'
#
# Use this to merge some histograms.

#for i in 1 2 3 4
#do

root -b <<EOF
.x hadd_files.C("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Embedding/Run12pp/MuDstMatching/pt35.geOut.list","/global/project/projectdirs/star/pwg/starjetc/dmawxc/Embedding/Run12pp/MuDstMatching/output/merged/pt35.geant.root")
EOF

#done
