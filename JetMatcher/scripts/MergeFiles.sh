#!/bin/bash
# 'MergeFiles.sh'
# Derek Anderson
# 09.25.2020
#
# For merging files using 'hadd_files.C'

# i/o parameters
path="/star/data01/pwg/dmawxc/Embedding/Run9pp/JetMatcher/output"
pref="pp200py8pt"
suff=".forJetEffPlotM04_pTbinFine.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root"
list="pp200py8merge.forJetEffPlotM04_pTbinFine.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.list"
root="pp200py8merge.forJetEffPlotM04_pTbinFine.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root"

# create list of files
nFiles=$(ls -1 $path/$pref*$suff | wc -l)
ls $path/$pref*$suff > $list

# merge files
root -b -q "MergeFiles.C($nFiles, \"$list\", \"$root\")"

# End -------------------------------------------------------------------------
