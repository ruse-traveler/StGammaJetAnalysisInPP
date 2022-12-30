#!/bin/bash
# 'TriggerResolutionStudy_MergeFiles.sh'
# Derek Anderson
# 09.25.2020
#
# For merging files using 'hadd_files.C'

# i/o parameters
path="/star/data01/pwg/dmawxc/Embedding/ParticleGun/output/OneEtaSiteTest"
pref="pionGridTest_oneEtaSite"
suff=".root"
list="pionGridTest_oneEtaSite_merge.et650h09pi0.d13m5y2022.list"
root="pionGridTest_oneEtaSite_merge.et650h09pi0.d13m5y2022.root"



# create list of files
nFiles=$(ls -1 $path/$pref*$suff | wc -l)
ls $path/$pref*$suff > $list

# merge files
root -b -q "TriggerResolutionStudy_MergeFiles.C($nFiles, \"$list\", \"$root\")"

# End -------------------------------------------------------------------------
