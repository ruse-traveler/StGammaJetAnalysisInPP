// 'MakeMuDstJetTree.C'
// Derek Anderson
// 11.13.2017
//
// This produces a tree of jets from the specified
// input with the specified parameters.
//
// NOTE: type = 0, makes charged jets
//       type = 1, makes full jets
//
// NOTE: trigger = 0, uses pi0 triggers
//       trigger = 1, uses h+- triggers
//       trigger = 2, use  NO trigger


#include "TString.h"
#include "TSystem.h"
#include "TDatime.h"

using namespace std;

class StMuDstJetTreeMaker;


// i/o parameters
static const TString  sInDefault("../../MuDstMatching/output/merged/pt7ff.matchWithMc.root");
static const TString  sOutDefault("pp200r9pt15ff.test.et920vz55pi0.r02rm1chrg.d4m2y2020.root");
static const Double_t pTdefault(7.);
// jet parameters
static const UInt_t   type(0);
static const UInt_t   trigger(2);
static const UInt_t   nRepeat(1);
static const UInt_t   nRemove(1);
static const Double_t rJet(0.5);
static const Double_t aGhost(0.01);
static const Double_t pTjetMin(0.2);
static const Double_t etaGhostMax(1.0 + rJet);
static const Double_t etaJetMax(1.0 - rJet);
static const Double_t etaBkgdMax(1.0);

// systematic parameters
static const Bool_t  adjustEff(false);
static const Float_t effAdjust(0.05);



void MakeMuDstJetTree(const Double_t pTparton=pTdefault, const TString sInput=sInDefault, const TString sOutput=sOutDefault, const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("StMuDstJetTreeMaker");

  // event/trigger parameters
  const Int_t    adcMax(6004);
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);
  const Double_t eEtaMin(0.5);
  const Double_t ePhiMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t etaTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  const Double_t tspPi0Min(0.);
  const Double_t tspPi0Max(0.08);
  const Double_t tspGamMin(0.2);
  const Double_t tspGamMax(0.6);
  // track parameters
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.0);
  const Double_t etaTrkMax(1.0);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(30.);
  // tower parameters
  const Double_t etaTwrMax(0.9);
  const Double_t eTwrMin(0.2);
  const Double_t eTwrMax(100.);
  const Double_t eCorrMin(0.2);
  const Double_t eCorrMax(30.);


  StMuDstJetTreeMaker *jetMaker = new StMuDstJetTreeMaker(isInBatchMode);
  // set parameters
  jetMaker -> SetInputAndOutputFiles(sInput.Data(), sOutput.Data(), pTparton);
  jetMaker -> SetEventParameters(rVtxMax, zVtxMax);
  jetMaker -> SetTriggerParameters(adcMax, eEtaMin, ePhiMin, pProjMax, etaTrgMax, eTtrgMin, eTtrgMax, tspPi0Min, tspPi0Max, tspGamMin, tspGamMax);
  jetMaker -> SetTrackParameters(nFitMin, rFitMin, dcaMax, etaTrkMax, pTtrkMin, pTtrkMax);
  jetMaker -> SetTowerParameters(etaTwrMax, eTwrMin, eTwrMax, eCorrMin, eCorrMax);
  jetMaker -> SetJetParameters(type, nRepeat, nRemove, rJet, aGhost, pTjetMin, etaGhostMax, etaJetMax, etaBkgdMax);
  if (adjustEff) jetMaker -> AdjustTrackEfficiency(adjustEff, effAdjust);
  // find jets
  jetMaker -> Init(trigger);
  jetMaker -> Make();
  jetMaker -> Finish();

}

// End ------------------------------------------------------------------------
