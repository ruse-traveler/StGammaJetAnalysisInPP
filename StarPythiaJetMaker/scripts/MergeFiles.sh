#!/bin/bash
# 'MergeFiles.sh'
# Derek Anderson
# 09.25.2020
#
# For merging files using 'hadd_files.C'

# i/o parameters
path="/star/data01/pwg/dmawxc/PythiaData/StarPythia/StarPythiaJetMaker/output"
pref="pp200py6"
suff=".root"
list="pp200py6.minBias_allFiles.et6100pi0.r05rm1chrg.d23m2y2022.list"
root="pp200py6.minBias_allFiles.et6100pi0.r05rm1chrg.d23m2y2022.root"



# create list of files
nFiles=$(ls -1 $path/$pref*$suff | wc -l)
ls $path/$pref*$suff > $list

# merge files
root -b -q "MergeFiles.C($nFiles, \"$list\", \"$root\")"

# End -------------------------------------------------------------------------
