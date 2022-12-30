#!/bin/bash
# 'RenameOuput.sh'
# Derek Anderson
#
# Use this to rename the incredible no.
# of files generated during the job
# submission process.

# old name
oldD="FullJetTree/"
pref=$oldD"tree.87311E8F310B172664D4DAA9A85E46BC_"
suff=".root"
nPre=${#pref}
nSuf=${#suff}
nNam=$(expr $nPre + $nSuf)

# new name
newD="FullJetTree/"
newP=$newD"pp200r9.d13m3y2017."
newS=".root"

printf "\nRenaming scripts...\n"

for script in `ls "$pref"*"$suff"`; do
  nTot=${#script}
  nNum=$(expr $nTot - $nNam)
  iJob=$(echo ${script:$nPre:$nNum})
  rename=$newP"job"$iJob$newS
  printf "  Moving '$script' to '$rename'...\n"
  mv -u $script $rename
done

printf "Finished renaming!\n\n"
