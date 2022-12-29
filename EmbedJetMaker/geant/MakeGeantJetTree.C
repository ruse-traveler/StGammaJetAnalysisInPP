// 'MakeGeantJetTree.C'
// Derek Anderson
// 10.03.2017
//
// This produces a tree of jets from the specified
// input with the specified parameters.
//
// NOTE: type = 0, makes charged jets
//       type = 1, makes full jets.


#include "TString.h"
#include "TSystem.h"
#include "TDatime.h"

using namespace std;

class StGeantJetTreeMaker;


// i/o parameters
static const TString  sInput("../../MuDstMatching/output/merged/pt35rff.matchWithMc.root");
static const TString  sOuput("pp200r12pt5g.r04rm1full.d6m1y2018.root");
static const Double_t pTparton(5);
// jet parameters
static const UInt_t   type(0);
static const UInt_t   nRepeat(1);
static const UInt_t   nRemove(1);
static const Double_t rJet(0.4);
static const Double_t aGhost(0.01);
static const Double_t pTjetMin(0.2);
static const Double_t etaGhostMax(1.0 + rJet);
static const Double_t etaJetMax(1.0 - rJet);
static const Double_t etaBkgdMax(1.0);


void MakeGeantJetTree(const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/Xsl64_gcc482/lib/libfastjet.so");
  gSystem -> Load("/opt/star/Xsl64_gcc482/lib/libfastjettools.so");
  gSystem -> Load("StGeantJetTreeMaker");

  // event/trigger parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);
  const Double_t etaTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  // track parameters
  const Double_t etaTrkMax(1.0);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(20.);


  StGeantJetTreeMaker *gntJetMaker = new StGeantJetTreeMaker(isInBatchMode);
  // set parameters
  gntJetMaker -> SetInputAndOutputFiles(sInput.Data(), sOuput.Data(), pTparton);
  gntJetMaker -> SetEventParameters(rVtxMax, zVtxMax);
  gntJetMaker -> SetTriggerParameters(etaTrgMax, eTtrgMin, eTtrgMax);
  gntJetMaker -> SetTrackParameters(etaTrkMax, pTtrkMin, pTtrkMax);
  gntJetMaker -> SetJetParameters(type, nRepeat, nRemove, rJet, aGhost, pTjetMin, etaGhostMax, etaJetMax, etaBkgdMax);
  // find jets
  gntJetMaker -> Init();
  gntJetMaker -> Make();
  gntJetMaker -> Finish();

}

// End ------------------------------------------------------------------------
