// 'DoEfficiencyCalculation.C'
// Derek Anderson
// 01.30.2020
//
// Use this to run the class StTrackEfficiencyCalculator,
// which calculates the tracking efficiency and produces
// trees of jets made of generated particles and reconst-
// ructed tracks matched to generated particles.
//
// NOTE: trigger = 0, uses pi0 triggers
//       trigger = 1, uses h+- triggers


#include "TString.h"
#include "TSystem.h"

using namespace std;

class StTrackEfficiencyCalculator;


// i/o parameters
static const TString sInDefault("../../MuDstMatching/output/merged/pt35rff.matchWithMc.root");
//static const TString sOutDefault("pp200r9pt35rff.checkingSigPlusDecays_bigNumEvents.et920pt02100vz55pi0.r02rm1chrg.d2m3y2020.root");
static const TString sOutDefault("test.root");



void DoEfficiencyCalculation(const TString sIn=sInDefault, const TString sOut=sOutDefault, const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("StTrackEfficiencyCalculator");

  // tree parameters
  const TString sParTree("McTracks");
  const TString sDetTree("GfmtoDst_mu");

  // jet parameters
  const UInt_t   nRepeat(1);
  const UInt_t   nRemove(1);
  const Double_t rJet(0.2);
  const Double_t pTmin(0.2);
  const Double_t aGhost(0.01);
  const Double_t hJetMax(1.0 - rJet);
  const Double_t hBkgdMax(1.0);
  const Double_t hGhostMax(1.0 + rJet);

  // calculation parameters
  const UInt_t trigger(0);
  const Bool_t useParTrg(true);
  const Bool_t useDetTrg(false);
  const Bool_t useRecoilTrks(true);
  const Bool_t useSpecificParticles(false);

  StTrackEfficiencyCalculator *calculator = new StTrackEfficiencyCalculator(isInBatchMode);
  calculator -> Init(sIn, sOut, sParTree, sDetTree);
  calculator -> SetJetParameters(nRepeat, nRemove, rJet, pTmin, aGhost, hJetMax, hBkgdMax, hGhostMax);
  calculator -> Make(trigger, useParTrg, useDetTrg, useRecoilTrks, useSpecificParticles);
  calculator -> Finish();

}

// End ------------------------------------------------------------------------
