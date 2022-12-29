// 'MakeMcJetTree.C'
// Derek Anderson
// 04.30.2018
//
// This produces a tree of jets from the specified
// input with the specified parameters.
//
// NOTE: type = 0, makes charged jets
//       type = 1, makes full jets.
//
// NOTE: event = 0, do *not* apply event cuts
//       event = 1, apply event cuts
//
// NOTE: trigger = 0, do *not* require trigger
//       trigger = 1, require trigger.


#include "TString.h"
#include "TSystem.h"
#include "TDatime.h"

using namespace std;

class StMcJetTreeMaker;


// i/o parameters
static const TString  sInDefault("../../MuDstMatching/output/merged/pt7ff.matchWithMc.root");
static const TString  sOutDefault("pp200r9pt7ff.forUnfolding_Par.et920pt021Kvz55had.r02rm1chrg.d5m3y2020.root");
static const Double_t pTdefault(7.);
// jet parameters
static const UInt_t   type(0);
static const UInt_t   event(1);
static const UInt_t   trigger(1);
static const UInt_t   nRepeat(1);
static const UInt_t   nRemove(1);
static const Double_t rJet(0.5);
static const Double_t aGhost(0.01);
static const Double_t pTjetMin(0.2);
static const Double_t etaGhostMax(1.0 + rJet);
static const Double_t etaJetMax(1.0 - rJet);
static const Double_t etaBkgdMax(1.0);

// systematic parameters
static const Bool_t  turnOffDecay(false);
static const Bool_t  adjustEff(false);
static const Float_t effAdjust(0.05);



void MakeMcJetTree(const Double_t pTparton=pTdefault, const TString sInput=sInDefault, const TString sOutput=sOutDefault, const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("StMcJetTreeMaker");

  // event/trigger parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);
  const Double_t etaTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  // track parameters
  const Double_t etaTrkMax(1.0);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(1000.);


  StMcJetTreeMaker *mcJetMaker = new StMcJetTreeMaker(isInBatchMode);
  // set parameters
  mcJetMaker -> SetInputAndOutputFiles(sInput.Data(), sOutput.Data(), pTparton);
  mcJetMaker -> SetEventParameters(rVtxMax, zVtxMax);
  mcJetMaker -> SetTriggerParameters(etaTrgMax, eTtrgMin, eTtrgMax);
  mcJetMaker -> SetTrackParameters(etaTrkMax, pTtrkMin, pTtrkMax);
  mcJetMaker -> SetJetParameters(type, nRepeat, nRemove, rJet, aGhost, pTjetMin, etaGhostMax, etaJetMax, etaBkgdMax);
  if (adjustEff)    mcJetMaker -> AdjustTrackEfficiency(adjustEff, effAdjust);
  if (turnOffDecay) mcJetMaker -> TurnOffWeakDecays(turnOffDecay);
  // find jets
  mcJetMaker -> Init();
  mcJetMaker -> Make(event, trigger);
  mcJetMaker -> Finish();

}

// End ------------------------------------------------------------------------
