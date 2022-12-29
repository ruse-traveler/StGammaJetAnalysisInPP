#!/bin/bash
# 'MergeFiles.sh'
# Derek Anderson
# 09.25.2020
#
# For merging files using 'hadd_files.C'

# i/o parameters
path="/star/data01/pwg/dmawxc/JetReco_pp/PythiaJetMaker"
pref="pp200py8pt"
suff="det.forUnfolding.et920nTrg50Kpt0230pi0.r05rm1chrg.d28m9y2020.root"
list="pp200py8det.forUnfolding.et920pt0230pi0.r05rm1chrg.d10m10y2020.list"
root="pp200py8det.forUnfolding.et920pt0230pi0.r05rm1chrg.d10m10y2020.root"



# create list of files
nFiles=$(ls -1 $path/$pref*$suff | wc -l)
ls $path/$pref*$suff > $list

# merge files
root -b -q "MergeFiles.C($nFiles, \"$list\", \"$root\")"

# End -------------------------------------------------------------------------
