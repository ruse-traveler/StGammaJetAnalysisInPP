#!/bin/bash
# 'MakeFemtoDstParallel.sh'

lvl=$1
typ=$2
if [ -z $lvl ]; then
  echo "HEY! Please specify \"particle\" or \"detector\"."
  exit
fi
if [ -z $typ ]; then
  echo "HEY! Please specify \"charged\" or \"full\"."
  exit
fi

# particle case
if [ $lvl = "particle" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDstParallel.C'(500000,0,-1,true,false,0,"pp200py8par.et920nTrg500Kpt0230pi0.r05rm1chrg.d3m11y2018.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDstParallel.C'(500000,0,-1,true,false,1,"pp200py8par.et920nTrg500Kpt0230pi0.r05rm1full.d3m11y2018.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# detector case
elif [ $lvl = "detector" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDstParallel.C'(500000,0,-1,false,true,0,"pp200py8det.resSys2.et920nTrg500Kpt0230pi0.r05rm1chrg.d4m11y2018.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDstParallel.C'(500000,0,-1,false,true,1,"pp200py8det.resSys2.et920nTrg500Kpt0230pi0.r05rm1full.d4m11y2018.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# other case
else
  echo "Hmmm... check what you typed..."
  exit
fi
