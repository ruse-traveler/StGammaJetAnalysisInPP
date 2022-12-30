#!/bin/bash
# 'RenameScripts.sh'
# Derek Anderson
#
# Use this to rename the incredible no.
# of scripts generated during the job
# submission process.

# old name
pref="sched87311E8F310B172664D4DAA9A85E46BC_"
suff=".csh"
nPre=${#pref}
nSuf=${#suff}
nNam=$(expr $nPre + $nSuf)

# new name
newPre="pp200r9.d13m3y2017."
newSuf=".csh"

printf "\nRenaming scripts...\n"

for script in `ls "$pref"*"$suff"`; do
  nTot=${#script}
  nNum=$(expr $nTot - $nNam)
  iJob=$(echo ${script:$nPre:$nNum})
  rename=$newPre"job"$iJob$newSuf
  printf "  Moving '$script' to '$rename'...\n"
  mv -u $script $rename
done

printf "Finished renaming!\n\n"
