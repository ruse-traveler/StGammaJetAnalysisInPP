// 'MakeFemtoDstParallel.C'
// Derek Anderson
// 08.05.2016
//
// This macro produces a compact tree of jets (the "FemtoDST") via the
// 'StFemtoDstMaker' class.  Use the switches 'particle' and 'detector'
// to turn on/off production of particle-level and detector-level jets
// respectively.
//
// Last updated: 04.17.2018

#include <TSystem>
#include <iostream>
#include "TString.h"

using namespace std;

class StFemtoDstMaker;


// i/o parameters
const TString ipFile("./input/pp200py8pi0.merged.root");
const TString idFile("./input/pp200py8pi0.merged.root");
const TString eFile("./input/pp200r9embed.efficiency.moreVariableBins.d11m4y2018.root");
const TString eHist("hPtEfficiency");
// trigger parameters
const Int_t    tID     = 111;
const Double_t eTmin   = 9.;
const Double_t eTmax   = 20.;
const Double_t hTrkMax = 1.;
const Double_t hTrgMax = 0.9;
// jet parameters
const Int_t    nRM   = 1;
const Double_t rJet  = 0.5;
const Double_t aMin  = 0.65;
const Double_t pTmin = 0.2;
const Double_t pTmax = 30.;
const Double_t qMin  = 0.15;
// systematic parameters
const Bool_t  effSys    = false;
const Float_t effAdjust = -0.05;



void MakeFemtoDstParallel(const Int_t nTrgs=-1, const Int_t StartEvt=0, const Int_t StopEvt=-1, const Bool_t particle, const Bool_t detector, const Int_t type=0, const TString oFile="oops.root") {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("/star/u/dmawxc/PDSF/PythiaJetMaker/.sl73_gcc485/lib/libStFemtoDstMaker.so");

  // lower verbosity
  gErrorIgnoreLevel = kError;


  // create particle-level jets
  if (particle) {
    StFemtoDstMaker p(ipFile, oFile, 0, type);
    p.Init(nRM, tID, rJet, aMin, pTmin, pTmax, qMin, eTmin, eTmax, hTrkMax, hTrgMax);
    p.Make(nTrgs, StartEvt, StopEvt);
    p.Finish();
  }

  // create detector-level jets
  if (detector) {
    StFemtoDstMaker d(idFile, oFile, 1, type);
    d.SetEfficiency(eFile, eHist, effSys, effAdjust);
    d.Init(nRM, tID, rJet, aMin, pTmin, pTmax, qMin, eTmin, eTmax, hTrkMax, hTrgMax);
    d.Make(nTrgs, StartEvt, StopEvt);
    d.Finish();
  }

}

// End ------------------------------------------------------------------------
