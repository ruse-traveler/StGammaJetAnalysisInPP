// 'DoTriggerEnergyScaleCalculation.C
// Derek Anderson
// 06.16.2020
//
// Use this to run the class StTriggerEnergyScaleCalculator,
// which calculates the trigger energy scale and produces
// trees of jets made of reconstructed tracks triggered on
// reconstructed pi0s and jets made of generated particles
// triggered on generated pi0s matched to the reconstructed
// ones.
//
// NOTE: TrgType controls what particle parameters to use
//   TrgType = 0, pi0
//   TrgType = 1, photon
//
// NOTE: MatchOrder controls what order the particle-to-
// detector matching is performed
//   MatchOrder = 0, detector-level clusters matched TO
//     particle-level particles
//   MatchOrder = 1, particle-level particles matched TO
//     detector-level clusters

#include "TString.h"
#include "TSystem.h"

using namespace std;

class StTriggerEnergyScaleCalculator;

// global constants
static const UInt_t  NTrgs(2);
static const UInt_t  NTrgBin(3);
static const UInt_t  NTrgTsp(2);
static const UInt_t  NHadIds(1);
static const Float_t TowerSideLength(0.05);

// io parameters
static const TString SOutDefault("test.root");
static const TString SParDefault("output/merged/pt35rff.thirdMakerWithMc.root");
static const TString SDetDefault("output/merged/pt35rff.thirdMakerWithMc.root");
static const TString SParTreeDefault("McTracks");
static const TString SMatTreeDefault("McTracksForMatching");
static const TString SDetTreeDefault("Gfmtodst");

// calculation parameters
static const Bool_t UseParTrigger(true);
static const Bool_t UseDetTrigger(true);
static const Bool_t RequireMatch(true);
static const Bool_t ExcludeNeighbors(true);
static const Bool_t AvoidDoubleCounting(false);

// sets pi0 vs. gamma parameters and match order
static const UInt_t TrgType(1);
static const UInt_t MatchOrder(0);



void DoTriggerEnergyScaleCalculation(const TString sOut=SOutDefault, const TString sPar=SParDefault, const TString sDet=SDetDefault, const TString sParTree=SParTreeDefault, const TString sMatTree=SMatTreeDefault, const TString sDetTree=SDetTreeDefault, const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("StTriggerEnergyScaleCalculator");

  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);

  // detector trigger paramters
  const Int_t    adcMax(6004);
  const Double_t eStrMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t hTrgDetMax(0.9);
  //const Double_t eTtrgDetMin(9.);
  //const Double_t eTtrgDetMax(20.);
  const Double_t eTtrgDetMin(0.);
  const Double_t eTtrgDetMax(100.);
  const Double_t eTbinMin[NTrgBin] = {9., 11., 15.};
  const Double_t eTbinMax[NTrgBin] = {11., 15., 20.};
  const Double_t tspPi0[NTrgTsp]   = {0., 0.08};
  const Double_t tspGam[NTrgTsp]   = {0.2, 0.6};
  //  const Double_t tspPi0[NTrgTsp]   = {0., 100.};
  // const Double_t tspGam[NTrgTsp]   = {0., 100.};

  // particle trigger parameters
  const Double_t cTrgPar(0.);
  const Double_t dRtrgMaxPi0(6. * TMath::Sqrt(2) * TowerSideLength);
  const Double_t dRtrgMaxGam(3. * TMath::Sqrt(2) * TowerSideLength);
  const Double_t hTrgParMax(0.9);
  const Double_t eTtrgMatMin(6.);
  const Double_t eTtrgMatMax(50.);
  const Double_t eTtrgParMin(6.);
  const Double_t eTtrgParMax(50.);
  //  const Double_t dMinSeparationPi0(6. * TMath::Sqrt(2) * TowerSideLength);
  const Double_t dMinSeparationPi0(3. * TMath::Sqrt(2) * TowerSideLength);
  const Double_t dMinSeparationGam(3. * TMath::Sqrt(2) * TowerSideLength);
  //  const Double_t aFitParms[NTrgs]  = {13511, 747053};
  const Double_t aFitParms[NTrgs]  = {747053, 13511};
  //these are switched!  const Double_t bFitParms[NTrgs]  = {5.67, 7.36};
  const Double_t bFitParms[NTrgs]  = {7.36, 5.67};
  const UInt_t   idHadPi0[NHadIds] = {7};
  const UInt_t   idHadGam[NHadIds] = {1};

  // track parameters
  const UInt_t   nFitTrkMin(15);
  const Double_t rFitTrkMin(0.52);
  const Double_t dcaTrkMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkParMin(0.);
  const Double_t pTtrkParMax(100.);
  const Double_t pTtrkDetMin(0.2);
  const Double_t pTtrkDetMax(30.);

  // jet parameters
  const UInt_t   nRepeat(1);
  const UInt_t   nRemove(1);
  const Double_t rJet(0.2);
  const Double_t pTmin(0.2);
  const Double_t aGhost(0.01);
  const Double_t hJetMax(1.0 - rJet);
  const Double_t hBkgdMax(1.0);
  const Double_t hGhostMax(1.0 + rJet);

  StTriggerEnergyScaleCalculator *calculator = new StTriggerEnergyScaleCalculator(isInBatchMode);
  calculator -> Init(sOut, sPar, sDet, sParTree, sMatTree, sDetTree);
  calculator -> SetEventParameters(rVtxMax, zVtxMax);
  calculator -> SetDetectorParameters(adcMax, eStrMin, pProjMax, hTrgDetMax, eTtrgDetMin, eTtrgDetMax, eTbinMin, eTbinMax, tspPi0, tspGam);
  calculator -> SetParticleParameters(cTrgPar, dRtrgMaxPi0, dRtrgMaxGam, hTrgParMax, eTtrgMatMin, eTtrgMatMax, eTtrgParMin, eTtrgParMax, dMinSeparationPi0, dMinSeparationGam, aFitParms, bFitParms, idHadPi0, idHadGam);
  calculator -> SetTrackParameters(nFitTrkMin, rFitTrkMin, dcaTrkMax, hTrkMax, pTtrkParMin, pTtrkParMax, pTtrkDetMin, pTtrkDetMax);
  calculator -> SetJetParameters(nRepeat, nRemove, rJet, pTmin, aGhost, hJetMax, hBkgdMax, hGhostMax);
  calculator -> Make(TrgType, MatchOrder, UseDetTrigger, UseParTrigger, RequireMatch, ExcludeNeighbors, AvoidDoubleCounting);
  calculator -> Finish();

}

// End ------------------------------------------------------------------------
