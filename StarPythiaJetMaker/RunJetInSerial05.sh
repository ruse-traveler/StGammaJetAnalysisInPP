#!/bin/bash
# 'RunJetInSerial05.sh'
# Use this to run 'RunJet05.C' and
# change input/output files

loop=$1
if [ -z $loop ]; then
  echo "Please specify 1, 2, 3, or 4!"
  echo "Usage: ./RunJetInSerial02.sh {1,2,3,4}"
  exit
fi

if [[ $loop != "1" && $loop != "2" && $loop != "3" && $loop != "4" ]]; then
  echo "Please specify 1, 2, 3, or 4!"
  echo "Usage: ./RunJetInSerial02.sh {1,2,3,4}"
  exit
fi

# for input lists
declare -a iList1
declare -a iList2
declare -a iList3
declare -a iList4
# for output files
declare -a oFile1
declare -a oFile2
declare -a oFile3
declare -a oFile4

# input lists for loop 1
iList1[0]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run1.list\""
iList1[1]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run2.list\""
iList1[2]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run3.list\""
# input lists for loop 2
iList2[0]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run4.list\""
iList2[1]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run5.list\""
iList2[2]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run6.list\""
iList2[3]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run7.list\""
# input lists for loop 3
iList3[0]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run8.list\""
iList3[1]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run9.list\""
iList3[2]="\"./input/pp200py6et6pt12.forJetRerun_noCorruptTrees_run10.list\""
iList3[3]="\"./input/pp200py6et6pt47.forJetRerun_noCorruptTrees_wave1run14.list\""
# input lists for loop 4
iList4[0]="\"./input/pp200py6et8pt12.forJetRerun_noCorruptTrees_run1a.list\""
iList4[1]="\"./input/pp200py6et8pt12.forJetRerun_noCorruptTrees_run1b.list\""
iList4[2]="\"./input/pp200py6et8pt12.forJetRerun_noCorruptTrees_run1c.list\""

# output path & prefix for loop 1
oFile1[0]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run1\""
oFile1[1]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run2\""
oFile1[2]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run3\""
# output path & prefix for loop 2
oFile2[0]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run4\""
oFile2[1]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run5\""
oFile2[2]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run6\""
oFile2[3]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run7\""
# output path & prefix for loop 3
oFile3[0]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run8\""
oFile3[1]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run9\""
oFile3[2]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt12.forJetRerun_noCorruptTrees_run10\""
oFile3[3]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et6pt47.forJetRerun_noCorruptTrees_wave1run14\""
# output path & prefix for loop 4
oFile4[0]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et8pt12.forJetRerun_run1a\""
oFile4[1]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et8pt12.forJetRerun_run1b\""
oFile4[2]="\"/gpfs01/star/pwg_tasks/hp01/output/derek_round2/pp200py6et8pt12.forJetRerun_run1c\""

# file loop 1
if [ $loop = "1" ]; then
  (( nFile1=0 ))
  for input1 in ${iList1[@]}; do

# run macro
root4star -b -l <<EOF
  .x RunJet05.C("high", $input1, ${oFile1[$nFile1]}, 0.5)
EOF

  (( nFile1++ ))
  done
fi

# file loop 2
if [ $loop = "2" ]; then
  (( nFile2=0 ))
  for input2 in ${iList2[@]}; do

# run macro
root4star -b -l <<EOF
  .x RunJet05.C("high", $input2, ${oFile2[$nFile2]}, 0.5)
EOF

  (( nFile2++ ))
  done
fi

# file loop 3
if [ $loop = "3" ]; then
  (( nFile3=0 ))
  for input3 in ${iList3[@]}; do

# run macro
root4star -b -l <<EOF
  .x RunJet05.C("high", $input3, ${oFile3[$nFile3]}, 0.5)
EOF

  (( nFile3++ ))
  done
fi

# file loop 4
if [ $loop = "4" ]; then
  (( nFile4=0 ))
  for input4 in ${iList4[@]}; do

# run macro
root4star -b -l <<EOF
  .x RunJet05.C("high", $input4, ${oFile4[$nFile4]}, 0.5)
EOF

  (( nFile4++ ))
  done
fi

# delete arrays
unset iList1
unset iList2
unset iList3
unset iList4
unset oFile1
unset oFile2
unset oFile3
unset oFile4

# End -------------------------------------------------------------------------
