#!/bin/bash
# 'MergeFiles.sh'
# Derek Anderson
# 09.25.2020
#
# For merging files using 'hadd_files.C'

# i/o parameters
path="./output"
pref="pp200r9data"
suff=".root"
list="pp200r9data.forEtCheck_merge.et0100r02.d5m3y2022.list"
root="pp200r9data.forEtCheck_merge.et0100r02.d5m3y2022.root"



# create list of files
nFiles=$(ls -1 $path/$pref*$suff | wc -l)
ls $path/$pref*$suff > $list

# merge files
root -b -q "MergeFiles.C($nFiles, \"$list\", \"$root\")"

# End -------------------------------------------------------------------------
