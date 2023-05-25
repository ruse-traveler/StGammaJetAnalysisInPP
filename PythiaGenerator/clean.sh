#!/bin/bash
# 'clean.sh'
# Derek Anderson
#
# Clean up the mess you've made!

lvl=$1
cwd=$PWD
dirtydir="../output1"


if [ $lvl == "deepest" ]; then

  echo "Time to get REALLY clean!"
  cd $cwd/$dirtydir
  for i in `seq 0 99`; do
    echo "  Cleaning $dirtydir/run$i..."
    cd run$i
    rm *.log
    rm err.txt
    rm out.txt
    rm *.root
    cd ../
  done
  cd $cwd
  echo "All clean!"

elif [ $lvl == "deep" ]; then

  echo "Time to get pretty clean!"
  cd $cwd/$dirtydir
  for i in `seq 0 99`; do
    echo "  Cleaning $dirtydir/run$i..."
    cd run$i
    rm *.log
    rm err.txt
    rm out.txt
    cd ../
  done
  cd $cwd
  echo "All clean!"

else

  echo "Time to get clean!"
  cd $cwd/$dirtydir
  for i in `seq 0 99`; do
    echo "  Cleaning $dirtydir/run$i..."
    cd run$i
    rm out.txt
    cd ../
  done
  cd $cwd
  echo "All clean!"

fi
