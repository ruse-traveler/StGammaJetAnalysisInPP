#!/bin/bash
# 'CreateInputChain.sh'
# Derek Anderson
# 06.24.2021
#
# For creating a TChain using
# 'CreateInputChain.C'.

# file lists
inList="pp200py6et6pt47.forJetRerun_goodFiles_wave1run14.list"
badList="pp200py6et6pt47.forJetRerun_badFilesForCatchingCorruptFiles_wave1run14.list"
goodList="pp200py6et6pt47.forJetRerun_goodFilesForCatchingCorruptFiles_wave1run14.list"

# output parameters
outFile="pp200py6et6pt47.forJetRerun_goodFilesForCatchingCorruptFiles_wave1run14.root"
outChain="PythiaTree"


# create TChain
root -b -q "CreateInputChain.C(\"$inList\", \"$outFile\", \"$outChain\", \"$badList\", \"$goodList\")"

# End -------------------------------------------------------------------------
