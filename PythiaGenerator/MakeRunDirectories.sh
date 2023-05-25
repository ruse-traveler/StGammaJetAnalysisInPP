#!/bin/bash
# 'MakeRunDirectories.sh'
#
# Creates output directories for running the Gamma-Jet Analysis
# on the Condor cluster.

cwd=$PWD
out="../pTpartonBins/pTpart5_gam"
run="run"
arg="GammaJetExecute.cmnd"

echo "Creating output directories..." 
mkdir $out
cd $out
for i in `seq 0 99`; do
  echo "  Run $i: $run$i"
  mkdir $run$i
  cp $cwd/$arg $run$i/
done
cd $cwd
echo "Finished!"
