#!/bin/bash
# 'MakeFemtoDst.sh'

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
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,0,"../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","pp200py8pt510par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,0,"../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","pp200py8pt1015par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,0,"../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","pp200py8pt1520par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,0,"../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","pp200py8pt2025par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,0,"../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","pp200py8pt2530par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,0,"../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","pp200py8pt30par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1chrg.d21m8y2021.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,1,"../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","pp200py8pt510par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,1,"../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","pp200py8pt1015par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,1,"../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","pp200py8pt1520par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,1,"../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","pp200py8pt2025par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,1,"../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","pp200py8pt2530par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,true,false,1,"../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","pp200py8pt30par.resSysV2_detCheck.et920nTrg50Kpt021Kpi0.r05rm1full.d21m8y2021.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# detector case
elif [ $lvl = "detector" ]; then
  if [ $typ = "charged" ]; then
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,0,"../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","pp200py8pt510det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,0,"../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","pp200py8pt1015det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,0,"../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","pp200py8pt1520det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,0,"../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","pp200py8pt2025det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,0,"../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","pp200py8pt2530det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1chrg.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,0,"../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","pp200py8pt30det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1chrg.d21m8y2021.root")'
  elif [ $typ = "full" ]; then
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,1,"../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root","pp200py8pt510det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,1,"../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1015pi0.merged.root","pp200py8pt1015det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,1,"../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt1520pi0.merged.root","pp200py8pt1520det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1full.d8m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,1,"../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2025pi0.merged.root","pp200py8pt2025det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,1,"../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt2530pi0.merged.root","pp200py8pt2530det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1full.d21m8y2021.root")'
    root -b -q MakeFemtoDst.C'(50000,0,-1,false,true,1,"../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","../../PythiaData/PtHatBins/pp200py8pt30pi0.merged.root","pp200py8pt30det.resSysV2_detCheck.et920nTrg50Kpt0230pi0.r05rm1full.d21m8y2021.root")'
  else
    echo "Hmmm... check what you typed..."
    exit
  fi
# other case
else
  echo "Hmmm... check what you typed..."
  exit
fi
