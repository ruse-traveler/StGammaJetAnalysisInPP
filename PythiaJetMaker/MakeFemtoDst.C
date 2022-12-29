// 'MakeFemtoDst.C'
// Derek Anderson
// 08.05.2016
//
// This macro produces a compact tree of jets (the "FemtoDST") via the
// 'StFemtoDstMaker' class.  Use the switches 'particle' and 'detector'
// to turn on/off production of particle-level and detector-level jets
// respectively.

#include <TSystem>
#include <iostream>
#include "TString.h"

using namespace std;

class StFemtoDstMaker;


// i/o parameters
const TString ipFile("../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root");
const TString idFile("../../PythiaData/PtHatBins/pp200py8pt510pi0.merged.root");
const TString eFile("./input/Response/smoothedEffWithPlot_forcedToMatchEff_coefficientsFitButExponentsFixed.et920vz55pi0.d6m2y2021.root");
const TString eHist("fEff_forceEffAndPseudoEffToConverge");
// trigger parameters
const Int_t    tID     = 111;
const Double_t eTmin   = 6.;
const Double_t eTmax   = 100.;
const Double_t hTrkMax = 1.;
const Double_t hTrgMax = 0.9;
// jet parameters
const Int_t    nRM   = 1;
const Double_t rJet  = 0.5;
const Double_t aMin  = 0.65;
const Double_t pTmin = 0.2;
const Double_t pTpar = 1000.;
const Double_t pTdet = 30.;
const Double_t qMin  = 0.15;
// efficiency parameters
const Bool_t  isHist    = false;
const Bool_t  isFunc    = true;
const Bool_t  effSys    = false;
const Bool_t  effSmooth = false;
const Float_t effAdjust = 0.04;
const Float_t effAmp    = 0.82;
const Float_t effSig    = 3.29;



void MakeFemtoDst(const Int_t nTrgs=-1, const Int_t StartEvt=0, const Int_t StopEvt=-1, const Bool_t particle, const Bool_t detector, const Int_t type=0, const TString iFileP=ipFile, const TString iFileD=idFile, const TString oFile="oops.root") {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("/star/data01/pwg/dmawxc/JetReco_pp/PythiaJetMaker/.sl73_gcc485/lib/libStFemtoDstMaker.so");

  // lower verbosity
  gErrorIgnoreLevel = kError;


  // create particle-level jets
  if (particle) {
    StFemtoDstMaker p(iFileP, oFile, 0, type);
    p.Init(nRM, tID, rJet, aMin, pTmin, pTpar, qMin, eTmin, eTmax, hTrkMax, hTrgMax);
    p.Make(nTrgs, StartEvt, StopEvt);
    p.Finish();
  }

  // create detector-level jets
  if (detector) {
    StFemtoDstMaker d(iFileD, oFile, 1, type);
    d.SetEfficiency(eFile, eHist, isHist, isFunc, effSmooth, effSys, effAdjust, effAmp, effSig);
    d.Init(nRM, tID, rJet, aMin, pTmin, pTdet, qMin, eTmin, eTmax, hTrkMax, hTrgMax);
    d.Make(nTrgs, StartEvt, StopEvt);
    d.Finish();
  }

}

// End ------------------------------------------------------------------------
