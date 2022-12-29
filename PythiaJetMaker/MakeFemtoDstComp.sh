#!/bin/bash
# 'MakeFemtoDstComp.sh'

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
    root -b -q MakeFemtoDst.C'(300000,0,-1,true,false,0,"./input/Pythia/pp200py8pi0.merged.root","./input/Pythia/pp200py8pi0.merged.root","pp200py8par.forQM22.et6100nTrg300Kpi0.r05rm1chrg.d27m3y2022.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(300000,0,-1,true,false,1,"./input/Pythia/pp200py8pi0.merged.root","./input/Pythia/pp200py8pi0.merged.root","pp200py8par.forQM22.et6100nTrg300Kpi0.r05rm1full.d27m3y2022.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# detector case
elif [ $lvl = "detector" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDst.C'(300000,0,-1,false,true,0,"./input/Pythia/pp200py8pi0.merged.root","./input/Pythia/pp200py8pi0.merged.root","pp200py8det.forQM22.et6100nTrg300Kpi0.r05rm1chrg.d27m3y2022.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(300000,0,-1,false,true,1,"./input/Pythia/pp200py8pi0.merged.root","./input/Pythia/pp200py8pi0.merged.root","pp200py8det.forQM22.et6100nTrg300Kpi0.r05rm1full.d27m3y2022.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# other case
else
  echo "Hmmm... check what you typed..."
  exit
fi
