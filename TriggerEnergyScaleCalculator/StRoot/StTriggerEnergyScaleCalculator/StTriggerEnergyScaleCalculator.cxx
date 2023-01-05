// 'StTriggerEnergyScaleCalculator.cxx'
// Derek Anderson, Saskia Mioduszewski
// 06.15.2020
//
// This class produces trees of jets using
// generated particles and reconstructed
// tracks triggered on generated pi0s
// matched to reconstructed pi0s.  Also
// calculates the trigger energy scale
// and resolution.
//
// NOTE: assumes particle-level and
//   detector-level input trees were
//   filled in same order.

#define StTriggerEnergyScaleCalculator_cxx

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StTriggerEnergyScaleCalculator.h"

using namespace std;
using namespace fastjet;

ClassImp(StTriggerEnergyScaleCalculator);



// general public methods -----------------------------------------------------

void StTriggerEnergyScaleCalculator::Init(const TString sOutput, const TString sParticle, const TString sDetector, const TString sParTree, const TString sMatTree, const TString sDetTree) {

  // announce start
  cout << "\n  Beginning trigger energy resolution script..." << endl;

  // set file and tree names
  _sOut     = sOutput;
  _sPar     = sParticle;
  _sDet     = sDetector;
  _sParTree = sParTree;
  _sMatTree = sMatTree;
  _sDetTree = sDetTree;

  // open files
  _fOut = new TFile(_sOut.Data(), "recreate");
  _fPar = new TFile(_sPar.Data(), "read");
  _fDet = new TFile(_sDet.Data(), "read");
  if (!_fPar || !_fDet) {
    cerr << "PANIC: couldn't open an input file!\n"
         << "  fPar = " << _fPar << "\n"
         << "  fDet = " << _fDet
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input trees
  _fPar -> GetObject(_sParTree.Data(), _tPar);
  _fDet -> GetObject(_sMatTree.Data(), _tMat);
  _fDet -> GetObject(_sDetTree.Data(), _tDet);
  if (!_tPar || !_tDet) {
    cerr << "PANIC: couldn't grab an input tree!\n"
         << "  tPar = " << _tPar << "\n"
         << "  tMat = " << _tMat << "\n"
         << "  tDet = " << _tDet
         << endl;
    return;
  }
  cout << "    Grabbed trees." << endl;

  // initialize trees and histograms
  InitTrees();
  InitHists();

}  // end 'Init(TString, TString, TString, TString, TString)'



void StTriggerEnergyScaleCalculator::Make(const UInt_t trgType, const UInt_t matchOrder, const Bool_t useDetTrigger, const Bool_t useParTrigger, const Bool_t requireMatch, const Bool_t excludeNeighbors, const Bool_t avoidDoubleCounting) {

  // set calculation parameters
  _trgType             = trgType;
  _matchOrder          = matchOrder;
  _useDetTrigger       = useDetTrigger;
  _useParTrigger       = useParTrigger;
  _requireMatch        = requireMatch;
  _excludeNeighbors    = excludeNeighbors;
  _avoidDoubleCounting = avoidDoubleCounting;
  cout << "    Set calculation parameters:\n"
       << "      matchOrder = " << _matchOrder << ", useDetTrigger = " << _useDetTrigger << ", useParTrigger = " << _useParTrigger << "\n"
       << "      requireMatch = " << _requireMatch << ", excludeNeighbors = " << _excludeNeighbors << ", avoidDoubleCounting = " << _avoidDoubleCounting
       << endl;

  // set pi0 and photon fit functions (for weighting)
  _fPi0 = new TF1("fPi0", "[0] * TMath::Power(x, -1. * [1])", _eTtrgParMin, _eTtrgParMax);
  _fGam = new TF1("fGam", "[0] * TMath::Power(x, -1. * [1])", _eTtrgParMin, _eTtrgParMax);
  _fPi0 -> SetParameter(0, _aFitParms[0]);
  _fPi0 -> SetParameter(1, _bFitParms[0]);
  _fGam -> SetParameter(0, _aFitParms[1]);
  _fGam -> SetParameter(1, _bFitParms[1]);
  cout << "    Set fit functions (for weighting).\n"
       << "    Selected trigger species:"
       << endl;

  // select function
  switch (_trgType) {
    case 0:
      _fWeight           = (TF1*) _fPi0 -> Clone();
      _dRtrgMaxUse       = _dRtrgMaxPi0;
      _dMinSeparationUse = _dMinSeparationPi0;
      _tspUse[0]         = _tspPi0[0];
      _tspUse[1]         = _tspPi0[1];
      for (UInt_t iTrg = 0; iTrg < NTrgIds; iTrg++) {
        _idTrgUse[iTrg] = _idTrgPi0[iTrg];
      }
      cout << "      Using pi0 fit, spacing, TSP values, and PIDs." << endl;
      break;
    case 1:
      _fWeight           = (TF1*) _fGam -> Clone();
      _dRtrgMaxUse       = _dRtrgMaxGam;
      _dMinSeparationUse = _dMinSeparationGam;
      _tspUse[0]         = _tspGam[0];
      _tspUse[1]         = _tspGam[1];
      for (UInt_t iTrg = 0; iTrg < NTrgIds; iTrg++) {
        _idTrgUse[iTrg] = _idTrgGam[iTrg];
      }
      cout << "      Using photon fit, spacing, TSP values, and PIDs." << endl;
      break;
    default:
      _fWeight           = (TF1*) _fPi0 -> Clone();
      _dRtrgMaxUse       = _dRtrgMaxPi0;
      _dMinSeparationUse = _dMinSeparationPi0;
      _tspUse[0]         = _tspPi0[0];
      _tspUse[1]         = _tspPi0[1];
      for (UInt_t iTrg = 0; iTrg < NTrgIds; iTrg++) {
        _idTrgUse[iTrg] = _idTrgPi0[iTrg];
      }
      cout << "      Using pi0 fit, spacing, TSP values, and PIDs (check 'TrgType' argument)." << endl;
      break;
  }

  // perform calculation
  switch (_matchOrder) {
    case 0:
      DoDetFirstMatching();
      break;
    case 1:
      DoParFirstMatching();
      break;
    default:
      DoDetFirstMatching();
      break;
  }

}  // end 'Make(UInt_t, UInt_t, Bool_t, Bool_t, Bool_t, Bool_t, Bool_t)'



void StTriggerEnergyScaleCalculator::Finish() {

  // save histograms and trees
  SaveTreesAndHists();

  // close files
  _fOut -> cd();
  _fOut -> Close();
  _fPar -> cd();
  _fPar -> Close();
  _fDet -> cd();
  _fDet -> Close();
  cout << "  Finished trigger energy resolution script!\n" << endl;

}  // end 'Finish()'



// parameter setters ----------------------------------------------------------

void StTriggerEnergyScaleCalculator::SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax) {

  _rVtxMax = rVtxMax;
  _zVtxMax = zVtxMax;
  cout << "    Set event parameters:\n"
       << "      rVtxMax = " << _rVtxMax << ", zVtxMax = " << _zVtxMax
       << endl;

}  // end 'SetEventParameters(Double_t, Double_t)'



void StTriggerEnergyScaleCalculator::SetDetectorParameters(const Int_t adcMax, const Double_t eStrMin, const Double_t pProjMax, const Double_t hTrgDetMax, const Double_t eTtrgDetMin, const Double_t eTtrgDetMax, const Double_t *eTbinMin, const Double_t *eTbinMax, const Double_t *tspPi0, const Double_t *tspGam) {

  _adcMax      = adcMax;
  _eStrMin     = eStrMin;
  _pProjMax    = pProjMax;
  _hTrgDetMax  = hTrgDetMax;
  _eTtrgDetMin = eTtrgDetMin;
  _eTtrgDetMax = eTtrgDetMax;
  _tspPi0[0]   = tspPi0[0];
  _tspPi0[1]   = tspPi0[1];
  _tspGam[0]   = tspGam[0];
  _tspGam[1]   = tspGam[1];
  for (UInt_t iBin = 0; iBin < NTrgBin; iBin++) {
    _eTbinMin[iBin] = eTbinMin[iBin];
    _eTbinMax[iBin] = eTbinMax[iBin];
  }
  cout << "    Set detector trigger parameters:\n"
       << "      adcMax = " << _adcMax << ", eStrMin = " << _eStrMin << ", pProjMax = " << _pProjMax << "\n"
       << "      hTrgDetMax = " << _hTrgDetMax << ", eTtrgDet = (" << _eTtrgDetMin << ", " << _eTtrgDetMax << ")\n"
       << "      tspPi0 = (" << _tspPi0[0] << ", " << _tspPi0[1] << "), tspGam = (" << _tspGam[0] << ", " << _tspGam[1] << ")"
       << endl;

  // announce detector trigger bins
  cout << "    Detector eTtrg bins:" << endl;
  for (UInt_t iBin = 0; iBin < NTrgBin; iBin++) {
    cout << "      Bin[" << iBin << "]: eTtrg = (" << _eTbinMin[iBin] << ", " << _eTbinMax[iBin] << ")" << endl;
  }

}  // end 'SetDetectorParameters(Int_t, Double_t, Double_t, Double_t, Double_t, Double_t, *Double_t, *Double_t, *Double_t, *Double_t)'



void StTriggerEnergyScaleCalculator::SetParticleParameters(const Double_t cTrgPar, const Double_t dRtrgMaxPi0, const Double_t dRtrgMaxGam, const Double_t hTrgParMax, const Double_t eTtrgMatMin, const Double_t eTtrgMatMax, const Double_t eTtrgParMin, const Double_t eTtrgParMax, const Double_t dMinSeparationPi0, const Double_t dMinSeparationGam, const Double_t *aFitParms, const Double_t *bFitParms, const UInt_t *idTrgPi0, const UInt_t *idTrgGam) {

  _cTrgPar           = cTrgPar;
  _dRtrgMaxPi0       = dRtrgMaxPi0;
  _dRtrgMaxGam       = dRtrgMaxGam;
  _hTrgParMax        = hTrgParMax;
  _eTtrgMatMin       = eTtrgMatMin;
  _eTtrgMatMax       = eTtrgMatMax;
  _eTtrgParMin       = eTtrgParMin;
  _eTtrgParMax       = eTtrgParMax;
  _dMinSeparationPi0 = dMinSeparationPi0;
  _dMinSeparationGam = dMinSeparationGam;
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    _aFitParms[iTrg] = aFitParms[iTrg];
    _bFitParms[iTrg] = bFitParms[iTrg];
  }
  for (UInt_t iTrg = 0; iTrg < NTrgIds; iTrg++) {
    _idTrgPi0[iTrg] = idTrgPi0[iTrg];
    _idTrgGam[iTrg] = idTrgGam[iTrg];
  }
  cout << "    Set particle trigger parameters:\n"
       << "      cTrgPar = " << _cTrgPar << ", dRtrgMax(pi0, photon) = (" << _dRtrgMaxPi0 << ", " << _dRtrgMaxGam << "), dMinSeparation(pi0, photon) = (" << _dMinSeparationPi0 << ", " << _dMinSeparationGam << ")\n"
       << "      hTrgParMax = " << _hTrgParMax << ", eTtrgMat = (" << _eTtrgMatMin << ", " << _eTtrgMatMax << "),  eTtrgPar = (" << _eTtrgParMin << ", " << _eTtrgParMax << ")\n"
       << "    Set pi0 and photon fit function parameters (for weighting):\n"
       << "      a(pi0, photon) = (" << _aFitParms[0] << ", " << _aFitParms[1] << ")\n"
       << "      b(pi0, photon) = (" << _bFitParms[0] << ", " << _bFitParms[1] << ")"
       << endl;

  // announce particle trigger id's
  cout << "    Pi0 trigger id's:\n"
       << "      ";
  for (UInt_t iTrg = 0; iTrg < NTrgIds; iTrg++) {
    cout << _idTrgPi0[iTrg];
    if ((iTrg + 1) < NTrgIds) {
      cout << ", ";
    } else {
      cout << endl;
    }
  }
  cout << "    Photon trigger id's:\n"
       << "      ";
  for (UInt_t iTrg = 0; iTrg < NTrgIds; iTrg++) {
    cout << _idTrgGam[iTrg];
    if ((iTrg + 1) < NTrgIds) {
      cout << ", ";
    } else {
      cout << endl;
    }
  }

}  // end 'SetParticleParameters(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, *Double_t, *UInt_t, *UInt_t)'



void StTriggerEnergyScaleCalculator::SetTrackParameters(const UInt_t nFitTrkMin, const Double_t rFitTrkMin, const Double_t dcaTrkMax, const Double_t hTrkMax, const Double_t pTtrkParMin, const Double_t pTtrkParMax, const Double_t pTtrkDetMin, const Double_t pTtrkDetMax) {

  _nFitTrkMin  = nFitTrkMin;
  _rFitTrkMin  = rFitTrkMin;
  _dcaTrkMax   = dcaTrkMax;
  _hTrkMax     = hTrkMax;
  _pTtrkParMin = pTtrkParMin;
  _pTtrkParMax = pTtrkParMax;
  _pTtrkDetMin = pTtrkDetMin;
  _pTtrkDetMax = pTtrkDetMax;

}  // end 'SetTrackParameters()'



void StTriggerEnergyScaleCalculator::SetJetParameters(const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t pTjetMin, const Double_t aGhost, const Double_t hJetMax, const Double_t hBkgdMax, const Double_t hGhostMax) {

  _nRepeat     = nRepeat;
  _nRemove     = nRemove;
  _rJet        = rJet;
  _pTjetMin    = pTjetMin;
  _aGhost      = aGhost;
  _etaJetMax   = hJetMax;
  _etaBkgdMax  = hBkgdMax;
  _etaGhostMax = hGhostMax;
  cout << "    Set jet parameters:\n"
       << "      rJet = " << _rJet << ", pTmin = " << _pTjetMin << ", hMax = " << _etaJetMax << ", nRemove = " << _nRemove << "\n"
       << "      nRepeat = " << _nRepeat << ", aGhost = " << _aGhost << ", hBkgd = " << _etaBkgdMax << ", hGhost = " << _etaGhostMax 
       << endl;

}  // end 'SetJetParameters(UInt_t, UInt_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



// private methods ------------------------------------------------------------

void StTriggerEnergyScaleCalculator::InitTrees() {

  _fOut    -> cd();
  _tJetPar = new TTree("ParticleJetTree", "a tree of particle-level jets");
  _tJetPar -> Branch("eventIndex", &_pEventIndex, "EventIndex/I");
  _tJetPar -> Branch("Refmult", &_pRefmult, "Refmult/D");
  _tJetPar -> Branch("RunId", &_pRunId, "RunId/I");
  _tJetPar -> Branch("NJets", &_pNJets, "NJets/I");
  _tJetPar -> Branch("PartonicPt", &_pPartonicPt, "PartonicPt/D");
  _tJetPar -> Branch("TSP", &_pTSP, "TSP/D");
  _tJetPar -> Branch("TrgEta", &_pTrgEta, "TrgEta/D");
  _tJetPar -> Branch("TrgPhi", &_pTrgPhi, "TrgPhi/D");
  _tJetPar -> Branch("TrgEtPar", &_pTrgEtPar, "TrgEtPar/D");
  _tJetPar -> Branch("TrgEtDet", &_pTrgEtDet, "TrgEtDet/D");
  _tJetPar -> Branch("Rho", &_pRho, "Rho/D");
  _tJetPar -> Branch("Sigma", &_pSigma, "Sigma/D");
  _tJetPar -> Branch("Vz", &_pVz, "Vz/D");
  _tJetPar -> Branch("JetIndex", &_pJetIndex);
  _tJetPar -> Branch("JetNCons", &_pJetNCons);
  _tJetPar -> Branch("JetPt", &_pJetPt);
  _tJetPar -> Branch("JetPtCorr", &_pJetPtCorr);
  _tJetPar -> Branch("JetEta", &_pJetEta);
  _tJetPar -> Branch("JetPhi", &_pJetPhi);
  _tJetPar -> Branch("JetE", &_pJetE);
  _tJetPar -> Branch("JetArea", &_pJetArea);
  _tJetPar -> Branch("JetConsPt", &_pJetConsPt);
  _tJetPar -> Branch("JetConsEta", &_pJetConsEta);
  _tJetPar -> Branch("JetConsPhi", &_pJetConsPhi);
  _tJetPar -> Branch("JetConsE", &_pJetConsE);

  _tJetDet = new TTree("DetectorJetTree", "a tree of detector-level jets");
  _tJetDet -> Branch("eventIndex", &_dEventIndex, "EventIndex/I");
  _tJetDet -> Branch("Refmult", &_dRefmult, "Refmult/D");
  _tJetDet -> Branch("NJets", &_dNJets, "NJets/I");
  _tJetDet -> Branch("RunId", &_dRunId, "RunId/I");
  _tJetDet -> Branch("PartonicPt", &_dPartonicPt, "PartonicPt/D");
  _tJetDet -> Branch("TSP", &_dTSP, "TSP/D");
  _tJetDet -> Branch("TrgEta", &_dTrgEta, "TrgEta/D");
  _tJetDet -> Branch("TrgPhi", &_dTrgPhi, "TrgPhi/D");
  _tJetDet -> Branch("TrgEtPar", &_dTrgEtPar, "TrgEtPar/D");
  _tJetDet -> Branch("TrgEtDet", &_dTrgEtDet, "TrgEtDet/D");
  _tJetDet -> Branch("Rho", &_dRho, "Rho/D");
  _tJetDet -> Branch("Sigma", &_dSigma, "Sigma/D");
  _tJetDet -> Branch("Vz", &_dVz, "Vz/D");
  _tJetDet -> Branch("JetIndex", &_dJetIndex);
  _tJetDet -> Branch("JetNCons", &_dJetNCons);
  _tJetDet -> Branch("JetPt", &_dJetPt);
  _tJetDet -> Branch("JetPtCorr", &_dJetPtCorr);
  _tJetDet -> Branch("JetEta", &_dJetEta);
  _tJetDet -> Branch("JetPhi", &_dJetPhi);
  _tJetDet -> Branch("JetE", &_dJetE);
  _tJetDet -> Branch("JetArea", &_dJetArea);
  _tJetDet -> Branch("JetConsPt", &_dJetConsPt);
  _tJetDet -> Branch("JetConsEta", &_dJetConsEta);
  _tJetDet -> Branch("JetConsPhi", &_dJetConsPhi);
  _tJetDet -> Branch("JetConsE", &_dJetConsE);

  // set autosave for every 5 MB
  _tJetPar -> SetAutoSave(-500000000);
  _tJetDet -> SetAutoSave(-500000000);
  cout << "    Initialized output trees." << endl;

  // set particle branch addresses
  _tPar -> SetMakeClass(1);
  _tPar -> SetBranchAddress("EventId", &_pMcEventId, &_bParMcEventId);
  _tPar -> SetBranchAddress("RunId", &_pMcRunId, &_bParMcRunId);
  _tPar -> SetBranchAddress("NumTrks", &_pMcNumTrks, &_bParMcNumTrks);
  _tPar -> SetBranchAddress("MuVtxX", &_pMuVx, &_bParMuVx);
  _tPar -> SetBranchAddress("MuVtxY", &_pMuVy, &_bParMuVy);
  _tPar -> SetBranchAddress("MuVtxZ", &_pMuVz, &_bParMuVz);
  _tPar -> SetBranchAddress("McVtxX", &_pMcVx, &_bParMcVx);
  _tPar -> SetBranchAddress("McVtxY", &_pMcVy, &_bParMcVy);
  _tPar -> SetBranchAddress("McVtxZ", &_pMcVz, &_bParMcVz);
  _tPar -> SetBranchAddress("IdTrk", &_pMcIdTrk, &_bParMcIdTrk);
  _tPar -> SetBranchAddress("IdGeant", &_pMcIdGeant, &_bParMcIdGeant);
  _tPar -> SetBranchAddress("IdVx", &_pMcIdVx, &_bParMcIdVx);
  _tPar -> SetBranchAddress("IdVxEnd", &_pMcIdVxEnd, &_bParMcIdVxEnd);
  _tPar -> SetBranchAddress("IntrVtx", &_pMcIntrVtx, &_bParMcIntrVtx);
  _tPar -> SetBranchAddress("IsShower", &_pMcIsShower, &_bParMcIsShower);
  _tPar -> SetBranchAddress("Charge", &_pMcCharge, &_bParMcCharge);
  _tPar -> SetBranchAddress("Rapidity", &_pMcRapidity, &_bParMcRapidity);
  _tPar -> SetBranchAddress("Eta", &_pMcEta, &_bParMcEta);
  _tPar -> SetBranchAddress("Phi", &_pMcPhi, &_bParMcPhi);
  _tPar -> SetBranchAddress("Px", &_pMcPx, &_bParMcPx);
  _tPar -> SetBranchAddress("Py", &_pMcPy, &_bParMcPy);
  _tPar -> SetBranchAddress("Pz", &_pMcPz, &_bParMcPz);
  _tPar -> SetBranchAddress("Pt", &_pMcPt, &_bParMcPt);
  _tPar -> SetBranchAddress("Ptot", &_pMcPtot, &_bParMcPtot);
  _tPar -> SetBranchAddress("Energy", &_pMcEnergy, &_bParMcEnergy);

  // set matched particle branch addresses
  _tMat -> SetMakeClass(1);
  _tMat -> SetBranchAddress("EventId", &_mMcEventId, &_bMatMcEventId);
  _tMat -> SetBranchAddress("RunId", &_mMcRunId, &_bMatMcRunId);
  _tMat -> SetBranchAddress("NumTrks", &_mMcNumTrks, &_bMatMcNumTrks);
  _tMat -> SetBranchAddress("MuVtxX", &_mMuVx, &_bMatMuVx);
  _tMat -> SetBranchAddress("MuVtxY", &_mMuVy, &_bMatMuVy);
  _tMat -> SetBranchAddress("MuVtxZ", &_mMuVz, &_bMatMuVz);
  _tMat -> SetBranchAddress("McVtxX", &_mMcVx, &_bMatMcVx);
  _tMat -> SetBranchAddress("McVtxY", &_mMcVy, &_bMatMcVy);
  _tMat -> SetBranchAddress("McVtxZ", &_mMcVz, &_bMatMcVz);
  _tMat -> SetBranchAddress("IdTrk", &_mMcIdTrk, &_bMatMcIdTrk);
  _tMat -> SetBranchAddress("IdGeant", &_mMcIdGeant, &_bMatMcIdGeant);
  _tMat -> SetBranchAddress("IdVx", &_mMcIdVx, &_bMatMcIdVx);
  _tMat -> SetBranchAddress("IdVxEnd", &_mMcIdVxEnd, &_bMatMcIdVxEnd);
  _tMat -> SetBranchAddress("IntrVtx", &_mMcIntrVtx, &_bMatMcIntrVtx);
  _tMat -> SetBranchAddress("IsShower", &_mMcIsShower, &_bMatMcIsShower);
  _tMat -> SetBranchAddress("Charge", &_mMcCharge, &_bMatMcCharge);
  _tMat -> SetBranchAddress("Rapidity", &_mMcRapidity, &_bMatMcRapidity);
  _tMat -> SetBranchAddress("Eta", &_mMcEta, &_bMatMcEta);
  _tMat -> SetBranchAddress("Phi", &_mMcPhi, &_bMatMcPhi);
  _tMat -> SetBranchAddress("Px", &_mMcPx, &_bMatMcPx);
  _tMat -> SetBranchAddress("Py", &_mMcPy, &_bMatMcPy);
  _tMat -> SetBranchAddress("Pz", &_mMcPz, &_bMatMcPz);
  _tMat -> SetBranchAddress("Pt", &_mMcPt, &_bMatMcPt);
  _tMat -> SetBranchAddress("Ptot", &_mMcPtot, &_bMatMcPtot);
  _tMat -> SetBranchAddress("Energy", &_mMcEnergy, &_bMatMcEnergy);

  // set detector branch addresses
  _tDet -> SetMakeClass(1);
  _tDet -> SetBranchAddress("fUniqueID", &fUniqueID, &_bEventList_fUniqueID);
  _tDet -> SetBranchAddress("fBits", &fBits, &_bEventList_fBits);
  _tDet -> SetBranchAddress("runNumber", &runNumber, &_bEventList_runNumber);
  _tDet -> SetBranchAddress("eventNumber", &eventNumber, &_bEventList_eventNumber);
  _tDet -> SetBranchAddress("trigID", &trigID, &_bEventList_trigID);
  _tDet -> SetBranchAddress("nGlobalTracks", &nGlobalTracks, &_bEventList_nGlobalTracks);
  _tDet -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks, &_bEventList_nPrimaryTracks);
  _tDet -> SetBranchAddress("refMult", &refMult, &_bEventList_refMult);
  _tDet -> SetBranchAddress("vpdVz", &vpdVz, &_bEventList_vpdVz);
  _tDet -> SetBranchAddress("xVertex", &xVertex, &_bEventList_xVertex);
  _tDet -> SetBranchAddress("yVertex", &yVertex, &_bEventList_yVertex);
  _tDet -> SetBranchAddress("zVertex", &zVertex, &_bEventList_zVertex);
  _tDet -> SetBranchAddress("bbcZVertex", &bbcZVertex, &_bEventList_bbcZVertex);
  _tDet -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate, &_bEventList_zdcCoincidenceRate);
  _tDet -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate, &_bEventList_bbcCoincidenceRate);
  _tDet -> SetBranchAddress("backgroundRate", &backgroundRate, &_bEventList_backgroundRate);
  _tDet -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate, &_bEventList_bbcBlueBackgroundRate);
  _tDet -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate, &_bEventList_bbcYellowBackgroundRate);
  _tDet -> SetBranchAddress("refMultPos", &refMultPos, &_bEventList_refMultPos);
  _tDet -> SetBranchAddress("refMultNeg", &refMultNeg, &_bEventList_refMultNeg);
  _tDet -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity, &_bEventList_bTOFTrayMultiplicity);
  _tDet -> SetBranchAddress("nVerticies", &nVerticies, &_bEventList_nVerticies);
  _tDet -> SetBranchAddress("MagF", &MagF, &_bEventList_MagF);
  _tDet -> SetBranchAddress("VrtxRank", &VrtxRank, &_bEventList_VrtxRank);
  _tDet -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &FlagEvent_TrgTrkMisMtch, &_bEventList_FlagEvent_TrgTrkMisMtch);
  _tDet -> SetBranchAddress("Etsp", &Etsp, &_bEventList_Etsp);
  _tDet -> SetBranchAddress("ETwrdidT", &ETwrdidT, &_bEventList_ETwrdidT);
  _tDet -> SetBranchAddress("ETwradc11", &ETwradc11, &_bEventList_ETwradc11);
  _tDet -> SetBranchAddress("ETwreneT0", &ETwreneT0, &_bEventList_ETwreneT0);
  _tDet -> SetBranchAddress("ETwreT", &ETwreT, &_bEventList_ETwreT);
  _tDet -> SetBranchAddress("ETwrENET0", &ETwrENET0, &_bEventList_ETwrENET0);
  _tDet -> SetBranchAddress("ETwrphT", &ETwrphT, &_bEventList_ETwrphT);
  _tDet -> SetBranchAddress("ETwrPTower", &ETwrPTower, &_bEventList_ETwrPTower);
  _tDet -> SetBranchAddress("ETwrpidTower", &ETwrpidTower, &_bEventList_ETwrpidTower);
  _tDet -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT, &_bEventList_ETwrmoduleT);
  _tDet -> SetBranchAddress("EClustEneT0", &EClustEneT0, &_bEventList_EClustEneT0);
  _tDet -> SetBranchAddress("EClustetav1", &EClustetav1, &_bEventList_EClustetav1);
  _tDet -> SetBranchAddress("EClustphiv1", &EClustphiv1, &_bEventList_EClustphiv1);
  _tDet -> SetBranchAddress("EEstrpen01", &EEstrpen01, &_bEventList_EEstrpen01);
  _tDet -> SetBranchAddress("EEstrpen02", &EEstrpen02, &_bEventList_EEstrpen02);
  _tDet -> SetBranchAddress("EEstrpen03", &EEstrpen03, &_bEventList_EEstrpen03);
  _tDet -> SetBranchAddress("EEstrpen0", &EEstrpen0, &_bEventList_EEstrpen0);
  _tDet -> SetBranchAddress("EEstrpen1", &EEstrpen1, &_bEventList_EEstrpen1);
  _tDet -> SetBranchAddress("EEstrpen2", &EEstrpen2, &_bEventList_EEstrpen2);
  _tDet -> SetBranchAddress("EEstrpen3", &EEstrpen3, &_bEventList_EEstrpen3);
  _tDet -> SetBranchAddress("EEstrpen4", &EEstrpen4, &_bEventList_EEstrpen4);
  _tDet -> SetBranchAddress("EEstrpen5", &EEstrpen5, &_bEventList_EEstrpen5);
  _tDet -> SetBranchAddress("EEstrpen6", &EEstrpen6, &_bEventList_EEstrpen6);
  _tDet -> SetBranchAddress("EEstrpen7", &EEstrpen7, &_bEventList_EEstrpen7);
  _tDet -> SetBranchAddress("EEstrpen8", &EEstrpen8, &_bEventList_EEstrpen8);
  _tDet -> SetBranchAddress("EEstrpen9", &EEstrpen9, &_bEventList_EEstrpen9);
  _tDet -> SetBranchAddress("EEstrpen10", &EEstrpen10, &_bEventList_EEstrpen10);
  _tDet -> SetBranchAddress("EEstrpen11", &EEstrpen11, &_bEventList_EEstrpen11);
  _tDet -> SetBranchAddress("EEstrpen12", &EEstrpen12, &_bEventList_EEstrpen12);
  _tDet -> SetBranchAddress("EEstrpen13", &EEstrpen13, &_bEventList_EEstrpen13);
  _tDet -> SetBranchAddress("EEstrpen14", &EEstrpen14, &_bEventList_EEstrpen14);
  _tDet -> SetBranchAddress("EEstrpen15", &EEstrpen15, &_bEventList_EEstrpen15);
  _tDet -> SetBranchAddress("ETwrdidE", &ETwrdidE, &_bEventList_ETwrdidE);
  _tDet -> SetBranchAddress("EPstripenp01", &EPstripenp01, &_bEventList_EPstripenp01);
  _tDet -> SetBranchAddress("EPstripenp02", &EPstripenp02, &_bEventList_EPstripenp02);
  _tDet -> SetBranchAddress("EPstripenp03", &EPstripenp03, &_bEventList_EPstripenp03);
  _tDet -> SetBranchAddress("EPstripenp0", &EPstripenp0, &_bEventList_EPstripenp0);
  _tDet -> SetBranchAddress("EPstripenp1", &EPstripenp1, &_bEventList_EPstripenp1);
  _tDet -> SetBranchAddress("EPstripenp2", &EPstripenp2, &_bEventList_EPstripenp2);
  _tDet -> SetBranchAddress("EPstripenp3", &EPstripenp3, &_bEventList_EPstripenp3);
  _tDet -> SetBranchAddress("EPstripenp4", &EPstripenp4, &_bEventList_EPstripenp4);
  _tDet -> SetBranchAddress("EPstripenp5", &EPstripenp5, &_bEventList_EPstripenp5);
  _tDet -> SetBranchAddress("EPstripenp6", &EPstripenp6, &_bEventList_EPstripenp6);
  _tDet -> SetBranchAddress("EPstripenp7", &EPstripenp7, &_bEventList_EPstripenp7);
  _tDet -> SetBranchAddress("EPstripenp8", &EPstripenp8, &_bEventList_EPstripenp8);
  _tDet -> SetBranchAddress("EPstripenp9", &EPstripenp9, &_bEventList_EPstripenp9);
  _tDet -> SetBranchAddress("EPstripenp10", &EPstripenp10, &_bEventList_EPstripenp10);
  _tDet -> SetBranchAddress("EPstripenp11", &EPstripenp11, &_bEventList_EPstripenp11);
  _tDet -> SetBranchAddress("EPstripenp12", &EPstripenp12, &_bEventList_EPstripenp12);
  _tDet -> SetBranchAddress("EPstripenp13", &EPstripenp13, &_bEventList_EPstripenp13);
  _tDet -> SetBranchAddress("EPstripenp14", &EPstripenp14, &_bEventList_EPstripenp14);
  _tDet -> SetBranchAddress("EPstripenp15", &EPstripenp15, &_bEventList_EPstripenp15);
  _tDet -> SetBranchAddress("EclustEnnq1", &EclustEnnq1, &_bEventList_EclustEnnq1);
  _tDet -> SetBranchAddress("EclustEnnq20", &EclustEnnq20, &_bEventList_EclustEnnq20);
  _tDet -> SetBranchAddress("EclustEnnq19", &EclustEnnq19, &_bEventList_EclustEnnq19);
  _tDet -> SetBranchAddress("EclustEnpq1", &EclustEnpq1, &_bEventList_EclustEnpq1);
  _tDet -> SetBranchAddress("EclustEnpq20", &EclustEnpq20, &_bEventList_EclustEnpq20);
  _tDet -> SetBranchAddress("EclustEnpq19", &EclustEnpq19, &_bEventList_EclustEnpq19);
  _tDet -> SetBranchAddress("EclustEnpq21", &EclustEnpq21, &_bEventList_EclustEnpq21);
  _tDet -> SetBranchAddress("EtspAlt", &EtspAlt, &_bEventList_EtspAlt);
  _tDet -> SetBranchAddress("EClustNumTwrIncluded", &EClustNumTwrIncluded, &_bEventList_EClustNumTwrIncluded);
  _tDet -> SetBranchAddress("EClustSecondTwrIndex", &EClustSecondTwrIndex, &_bEventList_EClustSecondTwrIndex);
  _tDet -> SetBranchAddress("EClustSecondTwrEnergy", &EClustSecondTwrEnergy, &_bEventList_EClustSecondTwrEnergy);
  _tDet -> SetBranchAddress("EEstrpModuleE", &EEstrpModuleE, &_bEventList_EEstrpModuleE);
  _tDet -> SetBranchAddress("EPstripModuleP", &EPstripModuleP, &_bEventList_EPstripModuleP);
  _tDet -> SetBranchAddress("EClustSecondTwrModule", &EClustSecondTwrModule, &_bEventList_EClustSecondTwrModule);
  _tDet -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &_bEventList_PrimaryTrackArray_);
  _tDet -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID, &_bPrimaryTrackArray_fUniqueID);
  _tDet -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits, &_bPrimaryTrackArray_fBits);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit, &_bPrimaryTrackArray_nHitsFit);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss, &_bPrimaryTrackArray_nHitsPoss);
  _tDet -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag, &_bPrimaryTrackArray_trackFlag);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ, &_bPrimaryTrackArray_pZ);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX, &_bPrimaryTrackArray_pX);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY, &_bPrimaryTrackArray_pY);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT, &_bPrimaryTrackArray_pT);
  _tDet -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx, &_bPrimaryTrackArray_dEdx);
  _tDet -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge, &_bPrimaryTrackArray_charge);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta, &_bPrimaryTrackArray_tofBeta);
  _tDet -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta, &_bPrimaryTrackArray_eta);
  _tDet -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi, &_bPrimaryTrackArray_phi);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron, &_bPrimaryTrackArray_nSigElectron);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion, &_bPrimaryTrackArray_nSigPion);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon, &_bPrimaryTrackArray_nSigKaon);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton, &_bPrimaryTrackArray_nSigProton);
  _tDet -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag, &_bPrimaryTrackArray_dcag);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits, &_bPrimaryTrackArray_nHits);
  _tDet -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits, &_bPrimaryTrackArray_dEdxHits);
  _tDet -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint, &_bPrimaryTrackArray_firstZPoint);
  _tDet -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint, &_bPrimaryTrackArray_lastZPoint);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron, &_bPrimaryTrackArray_tofSigElectron);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion, &_bPrimaryTrackArray_tofSigPion);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon, &_bPrimaryTrackArray_tofSigKaon);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton, &_bPrimaryTrackArray_tofSigProton);
  _tDet -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight, &_bPrimaryTrackArray_timeOfflight);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength, &_bPrimaryTrackArray_pathLength);
  _tDet -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex, &_bPrimaryTrackArray_trkIndex);
  _tDet -> SetBranchAddress("TowerArray", &TowerArray_, &_bEventList_TowerArray_);
  _tDet -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID, &_bTowerArray_fUniqueID);
  _tDet -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits, &_bTowerArray_fBits);
  _tDet -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId, &_bTowerArray_TwrId);
  _tDet -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng, &_bTowerArray_TwrEng);
  _tDet -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta, &_bTowerArray_TwrEta);
  _tDet -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi, &_bTowerArray_TwrPhi);
  _tDet -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC, &_bTowerArray_TwrADC);
  _tDet -> SetBranchAddress("TowerArray.TwrPed", TowerArray_TwrPed, &_bTowerArray_TwrPed);
  _tDet -> SetBranchAddress("TowerArray.TwrRMS", TowerArray_TwrRMS, &_bTowerArray_TwrRMS);
  _tDet -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex, &_bTowerArray_TwrMatchIdnex);
  _tDet -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk, &_bTowerArray_NoOfmatchedTrk);
  _tDet -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP, &_bTowerArray_TwrMatchP);
  _tDet -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx, &_bTowerArray_TwrPx);
  _tDet -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy, &_bTowerArray_TwrPy);
  _tDet -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz, &_bTowerArray_TwrPz);
  _tDet -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks, &_bTowerArray_fNAssocTracks);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_, &_bTowerArray_fMatchedTracksArray_);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P, &_bTowerArray_fMatchedTracksArray_P);
  cout << "    Initialized input trees." << endl;

}  // end 'InitTrees()'



void StTriggerEnergyScaleCalculator::InitHists() {

  const UInt_t  nNum(10);
  const UInt_t  nTsp(1000);
  const UInt_t  nPhi(360);
  const UInt_t  nEta(8);
  const UInt_t  nEt(200);
  const UInt_t  nEt2D(100);
  const UInt_t  nDr(800);
  const UInt_t  nQt(100);
  const UInt_t  nDq(50);
  const Float_t num[2]  = {0., 10.};
  const Float_t tsp[2]  = {0., 10.};
  const Float_t phi[2]  = {-6.3, 6.3};
  const Float_t eta[2]  = {-2., 2.};
  const Float_t eT[2]   = {0., 100.};
  const Float_t dR[2]   = {0., 20.};
  const Float_t qT[2]   = {0., 10.};
  const Float_t dQ[2]   = {0., 5.};
  const Float_t eT2d[2] = {0., 100.};
  // matching efficiences
  _hNumMatch             = new TH1D("hNumMatch", "", nNum, num[0], num[1]);
  _hNumVsCut             = new TH1D("hNumVsCut", "N_{trg} after successive cuts", 7, 0., 7.);
  _hRatioVsCut           = new TH1D("hRatioVsCut", "ratio of N_{trg} after cuts to total N_{trg}", 6, 0., 6.);
  _hPhiEff               = new TH1D("hPhiEff", "#varphi^{trg} matching efficiency", nPhi, phi[0], phi[1]);
  _hEtaEff               = new TH1D("hEtaEff", "#eta^{trg} matching efficiency", nEta, eta[0], eta[1]);
  _hEtEff                = new TH1D("hEtEff", "E_{T}^{trg} matching efficiency", nEt, eT[0], eT[1]);
  _hPhiForEff[0]         = new TH1D("hPhiForEffDet", "Detector #varphi^{trg} (for #epsilon_{match})", nPhi, phi[0], phi[1]);
  _hPhiForEff[1]         = new TH1D("hPhiForEffPar", "Particle #varphi^{trg} (for #epsilon_{match})", nPhi, phi[0], phi[1]);
  _hEtaForEff[0]         = new TH1D("hEtaForEffDet", "Detector #eta^{trg} (for #epsilon_{match})", nEta, eta[0], eta[1]);
  _hEtaForEff[1]         = new TH1D("hEtaForEffPar", "Particle #eta^{trg} (for #epsilon_{match})", nEta, eta[0], eta[1]);
  _hEtForEff[0]          = new TH1D("hEtForEffDet", "Detector E_{T}^{trg} (for #epsilon_{match})", nEt, eT[0], eT[1]);
  _hEtForEff[1]          = new TH1D("hEtForEffPar", "Particle E_{T}^{trg} (for #epsilon_{match})", nEt, eT[0], eT[1]);
  _hEtVsCut[0]           = new TH1D("hEtVsCut_vtx", "E_{T}^{trg}(reco) after vtx. cuts", nEt, eT[0], eT[1]);
  _hEtVsCut[1]           = new TH1D("hEtVsCut_twr", "E_{T}^{trg}(reco) after bad twr. rejection", nEt, eT[0], eT[1]);
  _hEtVsCut[2]           = new TH1D("hEtVsCut_proj", "E_{T}^{trg}(reco) after p_{proj} cut", nEt, eT[0], eT[1]);
  _hEtVsCut[3]           = new TH1D("hEtVsCut_eta", "E_{T}^{trg}(reco) after #eta^{trg} cut", nEt, eT[0], eT[1]);
  _hEtVsCut[4]           = new TH1D("hEtVsCut_match", "E_{T}^{trg}(reco) after matching", nEt, eT[0], eT[1]);
  _hEtVsCut[5]           = new TH1D("hEtVsCut_et", "E_{T}^{trg}(reco) after E_{T} cut", nEt, eT[0], eT[1]);
  _hEtVsCut[6]           = new TH1D("hEtVsCut_tsp", "E_{T}^{trg}(reco) after TSP cut", nEt, eT[0], eT[1]);
  // trigger tsp
  _hTspTrg[0]            = new TH1D("hTspTrgDetAll", "Detector TSP (all)", nTsp, tsp[0], tsp[1]);
  _hTspTrg[1]            = new TH1D("hTspTrgDetNonMatch", "Detector TSP (not matched)", nTsp, tsp[0], tsp[1]);
  _hTspTrg[2]            = new TH1D("hTspTrgDetMatchWeight", "Detector TSP (matched, weighted)", nTsp, tsp[0], tsp[1]);
  _hTspTrg[3]            = new TH1D("hTspTrgDetMatchNoWeight", "Detector TSP (matched, not weighted)", nTsp, tsp[0], tsp[1]);
  _hTspMatch[0]          = new TH1D("hTspMatch911", "TSP (matched), E_{T}^{trg}(det) = 9 - 11 GeV", nTsp, tsp[0], tsp[1]);
  _hTspMatch[1]          = new TH1D("hTspMatch1115", "TSP (matched), E_{T}^{trg}(det) = 11 - 15 GeV", nTsp, tsp[0], tsp[1]);
  _hTspMatch[2]          = new TH1D("hTspMatch1520", "TSP (matched), E_{T}^{trg}(det) = 15 - 20 GeV", nTsp, tsp[0], tsp[1]);
  _hTspMatchNW[0]        = new TH1D("hTspMatch911_NoWeight", "TSP (matched), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nTsp, tsp[0], tsp[1]);
  _hTspMatchNW[1]        = new TH1D("hTspMatch1115_NoWeight", "TSP (matched), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nTsp, tsp[0], tsp[1]);
  _hTspMatchNW[2]        = new TH1D("hTspMatch1520_NoWeight", "TSP (matched), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nTsp, tsp[0], tsp[1]);
  // trigger phi
  _hPhiTrg[0]            = new TH1D("hPhiTrgDetAll", "Detector #varphi^{trg} (all)", nPhi, phi[0], phi[1]);
  _hPhiTrg[1]            = new TH1D("hPhiTrgDetNonMatch", "Detector #varphi^{trg} (not matched)", nPhi, phi[0], phi[1]);
  _hPhiTrg[2]            = new TH1D("hPhiTrgDetMatchWeight", "Detector #varphi^{trg} (matched, weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[3]            = new TH1D("hPhiTrgDetMatchNoWeight", "Detector #varphi^{trg} (matched, not weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[4]            = new TH1D("hPhiTrgParCandWeight", "Particle #varphi^{trg} (candidate match, weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[5]            = new TH1D("hPhiTrgParCandNoWeight", "Particle #varphi^{trg} (candidate match, not weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[6]            = new TH1D("hPhiTrgParMatchWeight", "Particle #varphi^{trg} (matched, weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[7]            = new TH1D("hPhiTrgParMatchNoWeight", "Particle #varphi^{trg} (matched, not weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[8]            = new TH1D("hPhiTrgParAllWeight", "Particle #varphi^{trg} (all, weighted)", nPhi, phi[0], phi[1]);
  _hPhiTrg[9]            = new TH1D("hPhiTrgParAllNoWeight", "Particle #varphi^{trg} (all, not weighted)", nPhi, phi[0], phi[1]);
  // trigger eta
  _hEtaTrg[0]            = new TH1D("hEtaTrgDetAll", "Detector #eta^{trg} (all)", nEta, eta[0], eta[1]);
  _hEtaTrg[1]            = new TH1D("hEtaTrgDetNonMatch", "Detector #eta^{trg} (not matched)", nEta, eta[0], eta[1]);
  _hEtaTrg[2]            = new TH1D("hEtaTrgDetMatchWeight", "Detector #eta^{trg} (matched, weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[3]            = new TH1D("hEtaTrgDetMatchNoWeight", "Detector #eta^{trg} (matched, not weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[4]            = new TH1D("hEtaTrgParCandWeight", "Particle #eta^{trg} (candidate match, weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[5]            = new TH1D("hEtaTrgParCandNoWeight", "Particle #eta^{trg} (candidate match, not weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[6]            = new TH1D("hEtaTrgParMatchWeight", "Particle #eta^{trg} (matched, weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[7]            = new TH1D("hEtaTrgParMatchNoWeight", "Particle #eta^{trg} (matched, not weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[8]            = new TH1D("hEtaTrgParAllWeight", "Particle #eta^{trg} (all, weighted)", nEta, eta[0], eta[1]);
  _hEtaTrg[9]            = new TH1D("hEtaTrgParAllNoWeight", "Particle #eta^{trg} (all, not weighted)", nEta, eta[0], eta[1]);
  // trigger eT
  _hEtTrg[0]             = new TH1D("hEtTrgDetAll", "Detector E_{T}^{trg} (all)", nEt, eT[0], eT[1]);
  _hEtTrg[1]             = new TH1D("hEtTrgDetNonMatch", "Detector E_{T}^{trg} (not matched)", nEt, eT[0], eT[1]);
  _hEtTrg[2]             = new TH1D("hEtTrgDetMatchWeight", "Detector E_{T}^{trg} (matched, weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[3]             = new TH1D("hEtTrgDetMatchNoWeight", "Detector E_{T}^{trg} (matched, not weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[4]             = new TH1D("hEtTrgParCandWeight", "Particle E_{T}^{trg} (candidate match, weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[5]             = new TH1D("hEtTrgParCandNoWeight", "Particle E_{T}^{trg} (candidate match, not weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[6]             = new TH1D("hEtTrgParMatchWeight", "Particle E_{T}^{trg} (matched, weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[7]             = new TH1D("hEtTrgParMatchNoWeight", "Particle E_{T}^{trg} (matched, not weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[8]             = new TH1D("hEtTrgParAllWeight", "Particle E_{T}^{trg} (all, weighted)", nEt, eT[0], eT[1]);
  _hEtTrg[9]             = new TH1D("hEtTrgParAllNoWeight", "Particle E_{T}^{trg} (all, not weighted)", nEt, eT[0], eT[1]);
  _hEtDet[0]             = new TH1D("hEtDet911", "E_{T}^{trg}(det), E_{T}^{trg}(det) = 9 - 11 GeV", nEt, eT[0], eT[1]);
  _hEtDet[1]             = new TH1D("hEtDet1115", "E_{T}^{trg}(det), E_{T}^{trg}(det) = 11 - 15 GeV", nEt, eT[0], eT[1]);
  _hEtDet[2]             = new TH1D("hEtDet1520", "E_{T}^{trg}(det), E_{T}^{trg}(det) = 15 - 20 GeV", nEt, eT[0], eT[1]);
  _hEtDetNW[0]           = new TH1D("hEtDet911_NoWeight", "E_{T}^{trg}(det), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nEt, eT[0], eT[1]);
  _hEtDetNW[1]           = new TH1D("hEtDet1115_NoWeight", "E_{T}^{trg}(det), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nEt, eT[0], eT[1]);
  _hEtDetNW[2]           = new TH1D("hEtDet1520_NoWeight", "E_{T}^{trg}(det), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nEt, eT[0], eT[1]);
  _hEtMatch[0]           = new TH1D("hEtMatch911", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nEt, eT[0], eT[1]);
  _hEtMatch[1]           = new TH1D("hEtMatch1115", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nEt, eT[0], eT[1]);
  _hEtMatch[2]           = new TH1D("hEtMatch1520", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nEt, eT[0], eT[1]);
  _hEtMatchNW[0]         = new TH1D("hEtMatch911_NoWeight", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nEt, eT[0], eT[1]);
  _hEtMatchNW[1]         = new TH1D("hEtMatch1115_NoWeight", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nEt, eT[0], eT[1]);
  _hEtMatchNW[2]         = new TH1D("hEtMatch1520_NoWeight", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nEt, eT[0], eT[1]);
  // trigger dR
  _hDrTrg[0]             = new TH1D("hDrTrgCandWeight", "Candidate match #DeltaR^{trg} (weighted)", nDr, dR[0], dR[1]);
  _hDrTrg[1]             = new TH1D("hDrTrgCandNoWeight", "Candidate match #DeltaR^{trg} (not weighted)", nDr, dR[0], dR[1]);
  _hDrTrg[2]             = new TH1D("hDrTrgMatchWeight", "Matched particle #DeltaR^{trg} (weighted)", nDr, dR[0], dR[1]);
  _hDrTrg[3]             = new TH1D("hDrTrgMatchNoWeight", "Matched particle #DeltaR^{trg} (not weighted)", nDr, dR[0], dR[1]);
  _hDrMatch[0]           = new TH1D("hDrMatch911", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nDr, dR[0], dR[1]);
  _hDrMatch[1]           = new TH1D("hDrMatch1115", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nDr, dR[0], dR[1]);
  _hDrMatch[2]           = new TH1D("hDrMatch1520", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nDr, dR[0], dR[1]);
  _hDrMatchNW[0]         = new TH1D("hDrMatch911_NoWeight", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nDr, dR[0], dR[1]);
  _hDrMatchNW[1]         = new TH1D("hDrMatch1115_NoWeight", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nDr, dR[0], dR[1]);
  _hDrMatchNW[2]         = new TH1D("hDrMatch1520_NoWeight", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nDr, dR[0], dR[1]);
  // trigger qT
  _hQtTrg[0]             = new TH1D("hQtTrgCandWeight", "Candidate match q_{T}^{trg} (weighted)", nQt, qT[0], qT[1]);
  _hQtTrg[1]             = new TH1D("hQtTrgCandNoWeight", "Candidate match q_{T}^{trg} (not weighted)", nQt, qT[0], qT[1]);
  _hQtTrg[2]             = new TH1D("hQtTrgMatchWeight", "Matched particle q_{T}^{trg} (weighted)", nQt, qT[0], qT[1]);
  _hQtTrg[3]             = new TH1D("hQtTrgMatchNoWeight", "Matched particle q_{T}^{trg} (not weighted)", nQt, qT[0], qT[1]);
  _hQtMatch[0]           = new TH1D("hQtMatch911", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nQt, qT[0], qT[1]);
  _hQtMatch[1]           = new TH1D("hQtMatch1115", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nQt, qT[0], qT[1]);
  _hQtMatch[2]           = new TH1D("hQtMatch1520", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nQt, qT[0], qT[1]);
  _hQtMatchNW[0]         = new TH1D("hQtMatch911_NoWeight", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nQt, qT[0], qT[1]);
  _hQtMatchNW[1]         = new TH1D("hQtMatch1115_NoWeight", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nQt, qT[0], qT[1]);
  _hQtMatchNW[2]         = new TH1D("hQtMatch1520_NoWeight", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nQt, qT[0], qT[1]);
  // trigger dQ
  _hDqTrg[0]             = new TH1D("hDqTrgCandWeight", "Candidate match #Deltaq_{T}^{trg} (weighted)", nDq, dQ[0], dQ[1]);
  _hDqTrg[1]             = new TH1D("hDqTrgCandNoWeight", "Candidate match #Deltaq_{T}^{trg} (not weighted)", nDq, dQ[0], dQ[1]);
  _hDqTrg[2]             = new TH1D("hDqTrgMatchWeight", "Matched particle #Deltaq_{T}^{trg} (weighted)", nDq, dQ[0], dQ[1]);
  _hDqTrg[3]             = new TH1D("hDqTrgMatchNoWeight", "Matched particle #Deltaq_{T}^{trg} (not weighted)", nDq, dQ[0], dQ[1]);
  _hDqMatch[0]           = new TH1D("hDqMatch911", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nQt, qT[0], qT[1]);
  _hDqMatch[1]           = new TH1D("hDqMatch1115", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nQt, qT[0], qT[1]);
  _hDqMatch[2]           = new TH1D("hDqMatch1520", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nQt, qT[0], qT[1]);
  _hDqMatchNW[0]         = new TH1D("hDqMatch911_NoWeight", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nQt, qT[0], qT[1]);
  _hDqMatchNW[1]         = new TH1D("hDqMatch1115_NoWeight", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nQt, qT[0], qT[1]);
  _hDqMatchNW[2]         = new TH1D("hDqMatch1520_NoWeight", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nQt, qT[0], qT[1]);
  // trigger resolutions
  _hQtParVsEtPar         = new TH2D("hQtParVsEtPar", "q_{T}^{trg}(cand) vs. E_{T}^{trg}(cand)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtParVsEtParNW       = new TH2D("hQtParVsEtParNoWeight", "q_{T}^{trg}(cand) vs. E_{T}^{trg}(cand) (not weighted)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtParVsEtDet         = new TH2D("hQtParVsEtDet", "q_{T}^{trg}(cand) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtParVsEtDetNW       = new TH2D("hQtParVsEtDetNoWeight", "q_{T}^{trg}(cand) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hDrParVsEtPar         = new TH2D("hDrParVsEtPar", "#DeltaR^{trg}(cand) vs. E_{T}^{trg}(cand)", nEt2D, eT2d[0], eT2d[1], nDr, dR[0], dR[1]);
  _hDrParVsEtParNW       = new TH2D("hDrParVsEtParNoWeight", "#DeltaR^{trg}(cand) vs. E_{T}^{trg}(cand) (not weighted)", nEt2D, eT2d[0], eT2d[1], nDr, dR[0], dR[1]);
  _hDrParVsEtDet         = new TH2D("hDrParVsEtDet", "#DeltaR^{trg}(cand) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nDr, dR[0], dR[1]);
  _hDrParVsEtDetNW       = new TH2D("hDrParVsEtDetNoWeight", "#DeltaR^{trg}(cand) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nDr, dR[0], dR[1]);
  _hDrParVsQtPar         = new TH2D("hDrParVsQtPar", "#DeltaR^{trg}(cand) vs. q_{T}^{trg}(cand)", nQt, qT[0], qT[1], nDr, dR[0], dR[1]);
  _hDrParVsQtParNW       = new TH2D("hDrParVsQtParNoWeight", "#DeltaR^{trg}(cand) vs. q_{T}^{trg}(cand) (not weighted)", nQt, qT[0], qT[1], nDr, dR[0], dR[1]);
  _hPhiMatchVsDet        = new TH2D("hPhiMatchVsDet", "#varphi^{trg}(match) vs. #varphi^{trg}(det)", nPhi, phi[0], phi[1], nPhi, phi[0], phi[1]);
  _hPhiMatchVsDetNW      = new TH2D("hPhiMatchVsDetNoWeight", "#varphi^{trg}(match) vs. #varphi^{trg}(det) (not weighted)", nPhi, phi[0], phi[1], nPhi, phi[0], phi[1]);
  _hEtaMatchVsDet        = new TH2D("hEtaMatchVsDet", "#eta^{trg}(match) vs. #eta^{trg}(det)", nEta, eta[0], eta[1], nEta, eta[0], eta[1]);
  _hEtaMatchVsDetNW      = new TH2D("hEtaMatchVsDetNoWeight", "#eta^{trg}(match) vs. #eta^{trg}(det) (not weighted)", nEta, eta[0], eta[1], nEta, eta[0], eta[1]);
  _hEtMatchVsDet         = new TH2D("hEtMatchVsDet", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetNW       = new TH2D("hEtMatchVsDetNoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetTrg[0]   = new TH2D("hEtMatchVsDet911", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det), E_{T}^{trg}(det) = 9 - 11 GeV", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetTrg[1]   = new TH2D("hEtMatchVsDet1115", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det), E_{T}^{trg}(det) = 11 - 15 GeV", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetTrg[2]   = new TH2D("hEtMatchVsDet1520", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det), E_{T}^{trg}(det) = 15 - 20 GeV", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetTrgNW[0] = new TH2D("hEtMatchVsDet911_NoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det), E_{T}^{trg}(det) = 9 - 11 GeV (not weighted)", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetTrgNW[1] = new TH2D("hEtMatchVsDet1115_NoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det), E_{T}^{trg}(det) = 11 - 15 GeV (not weighted)", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hEtMatchVsDetTrgNW[2] = new TH2D("hEtMatchVsDet1520_NoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det), E_{T}^{trg}(det) = 15 - 20 GeV (not weighted)", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1]);
  _hDrMatchVsEtDet       = new TH2D("hDrMatchVsEtDet", "#DeltaR^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nDr, dR[0], dR[1]);
  _hDrMatchVsEtDetNW     = new TH2D("hDrMatchVsEtDetNoWeight", "#DeltaR^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nDr, dR[0], dR[1]);
  _hQtMatchVsEtMatch     = new TH2D("hQtMatchVsEtMatch", "q_{T}^{trg}(match) vs. E_{T}^{trg}(match)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtMatchVsEtMatchNW   = new TH2D("hQtMatchVsEtMatchNoWeight", "q_{T}^{trg}(match) vs. E_{T}^{trg}(match) (not weighted)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtMatchVsEtDet       = new TH2D("hQtMatchVsEtDet", "q_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtMatchVsEtDetNW     = new TH2D("hQtMatchVsEtDetNoWeight", "q_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtInvMatchVsEtDet    = new TH2D("hQtInvMatchVsEtDet", "(q_{T}^{trg})^{-1}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hQtInvMatchVsEtDetNW  = new TH2D("hQtInvMatchVsEtDetNoWeight", "(q_{T}^{trg})^{-1}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  _hDqMatchVsEtDet       = new TH2D("hDqMatchVsEtDet", "#Deltaq_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nDq, dQ[0], dQ[1]);
  _hDqMatchVsEtDetNW     = new TH2D("hDqMatchVsEtDetNoWeight", "#Deltaq_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], nDq, dQ[0], dQ[1]);
  _hQtEtMatchVsDet       = new TH3D("hQtMatchVsEtMatchVsEtDet", "q_{T}^{trg}(match) vs. E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], nEt2D, eT2d[0], eT2d[1], nQt, qT[0], qT[1]);
  // errors
  _hNumMatch             -> Sumw2();
  _hPhiEff               -> Sumw2();
  _hEtaEff               -> Sumw2();
  _hEtEff                -> Sumw2();
  _hPhiForEff[0]         -> Sumw2();
  _hPhiForEff[1]         -> Sumw2();
  _hEtaForEff[0]         -> Sumw2();
  _hEtaForEff[1]         -> Sumw2();
  _hEtForEff[0]          -> Sumw2();
  _hEtForEff[1]          -> Sumw2();
  _hEtVsCut[0]           -> Sumw2();
  _hEtVsCut[1]           -> Sumw2();
  _hEtVsCut[2]           -> Sumw2();
  _hEtVsCut[3]           -> Sumw2();
  _hEtVsCut[4]           -> Sumw2();
  _hEtVsCut[5]           -> Sumw2();
  _hEtVsCut[6]           -> Sumw2();
  _hTspTrg[0]            -> Sumw2();
  _hTspTrg[1]            -> Sumw2();
  _hTspTrg[2]            -> Sumw2();
  _hTspTrg[3]            -> Sumw2();
  _hTspMatch[0]          -> Sumw2();
  _hTspMatch[1]          -> Sumw2();
  _hTspMatch[2]          -> Sumw2();
  _hTspMatchNW[0]        -> Sumw2();
  _hTspMatchNW[1]        -> Sumw2();
  _hTspMatchNW[2]        -> Sumw2();
  _hPhiTrg[0]            -> Sumw2();
  _hPhiTrg[1]            -> Sumw2();
  _hPhiTrg[2]            -> Sumw2();
  _hPhiTrg[3]            -> Sumw2();
  _hPhiTrg[4]            -> Sumw2();
  _hPhiTrg[5]            -> Sumw2();
  _hPhiTrg[6]            -> Sumw2();
  _hPhiTrg[7]            -> Sumw2();
  _hPhiTrg[8]            -> Sumw2();
  _hPhiTrg[9]            -> Sumw2();
  _hEtaTrg[0]            -> Sumw2();
  _hEtaTrg[1]            -> Sumw2();
  _hEtaTrg[2]            -> Sumw2();
  _hEtaTrg[3]            -> Sumw2();
  _hEtaTrg[4]            -> Sumw2();
  _hEtaTrg[5]            -> Sumw2();
  _hEtaTrg[6]            -> Sumw2();
  _hEtaTrg[7]            -> Sumw2();
  _hEtaTrg[8]            -> Sumw2();
  _hEtaTrg[9]            -> Sumw2();
  _hEtTrg[0]             -> Sumw2();
  _hEtTrg[1]             -> Sumw2();
  _hEtTrg[2]             -> Sumw2();
  _hEtTrg[3]             -> Sumw2();
  _hEtTrg[4]             -> Sumw2();
  _hEtTrg[5]             -> Sumw2();
  _hEtTrg[6]             -> Sumw2();
  _hEtTrg[7]             -> Sumw2();
  _hEtTrg[8]             -> Sumw2();
  _hEtTrg[9]             -> Sumw2();
  _hEtDet[0]             -> Sumw2();
  _hEtDet[1]             -> Sumw2();
  _hEtDet[2]             -> Sumw2();
  _hEtDetNW[0]           -> Sumw2();
  _hEtDetNW[1]           -> Sumw2();
  _hEtDetNW[2]           -> Sumw2();
  _hEtMatch[0]           -> Sumw2();
  _hEtMatch[1]           -> Sumw2();
  _hEtMatch[2]           -> Sumw2();
  _hEtMatchNW[0]         -> Sumw2();
  _hEtMatchNW[1]         -> Sumw2();
  _hEtMatchNW[2]         -> Sumw2();
  _hDrTrg[0]             -> Sumw2();
  _hDrTrg[1]             -> Sumw2();
  _hDrTrg[2]             -> Sumw2();
  _hDrTrg[3]             -> Sumw2();
  _hDrMatch[0]           -> Sumw2();
  _hDrMatch[1]           -> Sumw2();
  _hDrMatch[2]           -> Sumw2();
  _hDrMatchNW[0]         -> Sumw2();
  _hDrMatchNW[1]         -> Sumw2();
  _hDrMatchNW[2]         -> Sumw2();
  _hQtTrg[0]             -> Sumw2();
  _hQtTrg[1]             -> Sumw2();
  _hQtTrg[2]             -> Sumw2();
  _hQtTrg[3]             -> Sumw2();
  _hQtMatch[0]           -> Sumw2();
  _hQtMatch[1]           -> Sumw2();
  _hQtMatch[2]           -> Sumw2();
  _hQtMatchNW[0]         -> Sumw2();
  _hQtMatchNW[1]         -> Sumw2();
  _hQtMatchNW[2]         -> Sumw2();
  _hDqTrg[0]             -> Sumw2();
  _hDqTrg[1]             -> Sumw2();
  _hDqTrg[2]             -> Sumw2();
  _hDqTrg[3]             -> Sumw2();
  _hDqMatch[0]           -> Sumw2();
  _hDqMatch[1]           -> Sumw2();
  _hDqMatch[2]           -> Sumw2();
  _hDqMatchNW[0]         -> Sumw2();
  _hDqMatchNW[1]         -> Sumw2();
  _hDqMatchNW[2]         -> Sumw2();
  _hQtParVsEtPar         -> Sumw2();
  _hQtParVsEtParNW       -> Sumw2();
  _hQtParVsEtDet         -> Sumw2();
  _hQtParVsEtDetNW       -> Sumw2();
  _hDrParVsEtPar         -> Sumw2();
  _hDrParVsEtParNW       -> Sumw2();
  _hDrParVsEtDet         -> Sumw2();
  _hDrParVsEtDetNW       -> Sumw2();
  _hDrParVsQtPar         -> Sumw2();
  _hDrParVsQtParNW       -> Sumw2();
  _hPhiMatchVsDet        -> Sumw2();
  _hPhiMatchVsDetNW      -> Sumw2();
  _hEtaMatchVsDet        -> Sumw2();
  _hEtaMatchVsDetNW      -> Sumw2();
  _hEtMatchVsDet         -> Sumw2();
  _hEtMatchVsDetNW       -> Sumw2();
  _hEtMatchVsDetTrg[0]   -> Sumw2();
  _hEtMatchVsDetTrg[1]   -> Sumw2();
  _hEtMatchVsDetTrg[2]   -> Sumw2();
  _hEtMatchVsDetTrgNW[0] -> Sumw2();
  _hEtMatchVsDetTrgNW[1] -> Sumw2();
  _hEtMatchVsDetTrgNW[2] -> Sumw2();
  _hDrMatchVsEtDet       -> Sumw2();
  _hDrMatchVsEtDetNW     -> Sumw2();
  _hQtMatchVsEtMatch     -> Sumw2();
  _hQtMatchVsEtMatchNW   -> Sumw2();
  _hQtMatchVsEtDet       -> Sumw2();
  _hQtMatchVsEtDetNW     -> Sumw2();
  _hQtInvMatchVsEtDet    -> Sumw2();
  _hQtInvMatchVsEtDetNW  -> Sumw2();
  _hDqMatchVsEtDet       -> Sumw2();
  _hDqMatchVsEtDetNW     -> Sumw2();
  _hQtEtMatchVsDet       -> Sumw2();

  // trigger profiles
  _pPhiMatchVsDet        = new TProfile("pPhiMatchVsDet", "#varphi^{trg}(match) vs. #varphi^{trg}(det)", nPhi, phi[0], phi[1], "S");
  _pPhiMatchVsDetNW      = new TProfile("pPhiMatchVsDetNoWeight", "#varphi^{trg}(match) vs. #varphi^{trg}(det) (not weighted)", nPhi, phi[0], phi[1], "S");
  _pEtaMatchVsDet        = new TProfile("pEtaMatchVsDet", "#eta^{trg}(match) vs. #eta^{trg}(det)", nEta, eta[0], eta[1], "S");
  _pEtaMatchVsDetNW      = new TProfile("pEtaMatchVsDetNoWeight", "#eta^{trg}(match) vs. #eta^{trg}(det) (not weighted)", nEta, eta[0], eta[1], "S");
  _pEtMatchVsDet         = new TProfile("pEtMatchVsDet", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetNW       = new TProfile("pEtMatchVsDetNoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetTrg[0]   = new TProfile("pEtMatchVsDet911", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetTrg[1]   = new TProfile("pEtMatchVsDet1115", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetTrg[2]   = new TProfile("pEtMatchVsDet1520", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetTrgNW[0] = new TProfile("pEtMatchVsDet911_NoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetTrgNW[1] = new TProfile("pEtMatchVsDet1115_NoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pEtMatchVsDetTrgNW[2] = new TProfile("pEtMatchVsDet1520_NoWeight", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pDrMatchVsEtDet       = new TProfile("pDrMatchVsEtDet", "#DeltaR^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pDrMatchVsEtDetNW     = new TProfile("pDrMatchVsEtDetNoWeight", "#DeltaR^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pQtMatchVsEtMatch     = new TProfile("pQtMatchVsEtMatch", "q_{T}^{trg}(match) vs. E_{T}^{trg}(par)", nEt2D, eT2d[0], eT2d[1], "S");
  _pQtMatchVsEtMatchNW   = new TProfile("pQtMatchVsEtMatchNoWeight", "q_{T}^{trg}(match) vs. E_{T}^{trg}(par) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pQtMatchVsEtDet       = new TProfile("pQtMatchVsEtDet", "q_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pQtMatchVsEtDetNW     = new TProfile("pQtMatchVsEtDetNoWeight", "q_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pQtInvMatchVsEtDet    = new TProfile("pQtInvMatchVsEtDet", "(q_{T}^{trg})^{-1}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pQtInvMatchVsEtDetNW  = new TProfile("pQtInvMatchVsEtDetNoWeight", "(q_{T}^{trg})^{-1}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");
  _pDqMatchVsEtDet       = new TProfile("pDqMatchVsEtDet", "#Deltaq_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt2D, eT2d[0], eT2d[1], "S");
  _pDqMatchVsEtDetNW     = new TProfile("pDqMatchVsEtDetNoWeight", "#Deltaq_{T}^{trg}(match) vs. E_{T}^{trg}(det) (not weighted)", nEt2D, eT2d[0], eT2d[1], "S");

  // trigger checks
  _hDidT_weird                   = new TH1D("hDidT_weird", "didT (qTtrg > 1.5)", 4802, -1., 4801.);
  _hDidT_normal                  = new TH1D("hDidT_normal", "didT (qTtrg < 1.5)", 4802, -1., 4801.);
  _hDidE_weird                   = new TH1D("hDidE_weird", "didE (qTtrg > 1.5)", 18002, -1., 18001.);
  _hDidE_normal                  = new TH1D("hDidE_normal", "didE (qTtrg < 1.5)", 18002, -1., 18001.);
  _hTrigID_weird                 = new TH1D("hTrigID_weird", "trigID (qTtrg > 1.5)", 20, -10., 10.);
  _hTrigID_normal                = new TH1D("hTrigID_normal", "trigID (qTtrg < 1.5)", 20, -10., 10.);
  _hModuleT_weird                = new TH1D("hModuleT_weird", "moduleT (qTtrg > 1.5)", 122, -1., 121.);
  _hModuleT_normal               = new TH1D("hModuleT_normal", "moduleT (qTtrg < 1.5)", 122, -1., 121.);
  _hEeneT0_weird                 = new TH1D("hEeneT0_weird", "eneT0 (twr) (qTtrg > 1.5)", 200, 0., 100.);
  _hEeneT0_normal                = new TH1D("hEeneT0_normal", "eneT0 (twr) (qTtrg < 1.5)", 200, 0., 100.);
  _hEENET0_weird                 = new TH1D("hEENET0_weird", "ENET0 (twr) (qTtrg > 1.5)", 200, 0., 100.);
  _hEENET0_normal                = new TH1D("hEENET0_normal", "ENET0 (twr) (qTtrg < 1.5)", 200, 0., 100.);
  _hEEneT0_weird                 = new TH1D("hEEneT0_weird", "EneT0 (clust) (qTtrg > 1.5)", 200, 0., 100.);
  _hEEneT0_normal                = new TH1D("hEEneT0_normal", "EneT0 (clust) (qTtrg < 1.5)", 200, 0., 100.);
  _hEtReco_weird                 = new TH1D("hEtReco_weird", "eTreco (clust) (qTtrg > 1.5)", 200, 0., 100.);
  _hEtReco_normal                = new TH1D("hEtReco_normal", "eTreco (clust) (qTtrg < 1.5)", 200, 0., 100.);
  _hEtSim_weird                  = new TH1D("hEtSim_weird", "eTsim (par) (qTtrg > 1.5)", 200, 0., 100.);
  _hEtSim_normal                 = new TH1D("hEtSim_normal", "eTsim (par) (qTtrg < 1.5)", 200, 0., 100.);
  _hEenToEEnRatio_weird          = new TH1D("hEenToEEnRatio_weird", "eneT0/EneT0 (qTtrg > 1.5)", 504, -1., 11.);
  _hEenToEEnRatio_normal         = new TH1D("hEenToEEnRatio_normal", "eneT0/EneT0 (qTtrg < 1.5)", 504, -1., 11.);
  _hNumTwrInClust_weird          = new TH1D("hNumTwrInClust_weird", "Num twrs in cluster (qTtrg > 1.5)", 7, -1., 6.);
  _hNumTwrInClust_normal         = new TH1D("hNumTwrInClust_normal", "Num twrs in cluster (qTtrg < 1.5)", 7, -1., 6.);
  _hEtaModule                    = new TH1D("hEtaModule", "Module ID of central eta strip", 122, -1, 121);
  _hPhiModule                    = new TH1D("hPhiModule", "Module ID of central phi strip", 122, -1, 121);
  _hTwrModule                    = new TH1D("hTwrModule", "Module ID of leading tower", 122, -1, 121);
  _hSecModule                    = new TH1D("hSecModule", "Module ID of second tower (only filled if 2 towers were in cluster)", 122, -1, 121);
  _hDiffPhiMod                   = new TH1D("hDiffPhiMod", "idMod(Eta) - idMod(Phi)", 242, -121, 121);
  _hDiffTwrMod                   = new TH1D("hDiffTwrMod", "idMod(Eta) - idMod(Twr)", 242, -121, 121);
  _hDiffSecMod                   = new TH1D("hDiffSecMod", "idMod(Eta) - idMod(2nd) (only filled if 2 towers were in cluster)", 242, -121, 121);
  _hSecondTwrID_weird            = new TH1D("hSecondTwrId_weird", "ID of 2nd twr in cluster (qTtrg > 1.5)", 44, -22., 22.);
  _hSecondTwrID_normal           = new TH1D("hSecondTwrId_normal", "ID of 2nd twr in cluster (qTtrg < 1.5)", 44, -22., 22.);
  _hSecondTwrEnergy_weird        = new TH1D("hSecondTwrEnergy_weird", "Energy of 2nd twr in cluster (qTtrg > 1.5)", 200, 0., 100.);
  _hSecondTwrEnergy_normal       = new TH1D("hSecondTwrEnergy_normal", "Energy of 2nd twr in cluster (qTtrg < 1.5)", 200, 0., 100.);
  _hSecondVsLeadTwrEnergy_weird  = new TH2D("hSecondVsLeadTwrEnergy_weird", "Energy of 2nd vs. lead twr in cluster (qTtrg > 1.5);eneT0;SecondTwrEnergy", 200, 0., 100., 200, 0., 100.);
  _hSecondVsLeadTwrEnergy_normal = new TH2D("hSecondVsLeadTwrEnergy_normal", "Energy of 2nd vs. lead twr in cluster (qTtrg < 1.5);eneT0;SecondTwrEnergy", 200, 0., 100., 200, 0., 100.);
  _hDidT_weird                   -> Sumw2();
  _hDidT_normal                  -> Sumw2();
  _hDidE_weird                   -> Sumw2();
  _hDidE_normal                  -> Sumw2();
  _hTrigID_weird                 -> Sumw2();
  _hTrigID_normal                -> Sumw2();
  _hModuleT_weird                -> Sumw2();
  _hModuleT_normal               -> Sumw2();
  _hEeneT0_weird                 -> Sumw2();
  _hEeneT0_normal                -> Sumw2();
  _hEENET0_weird                 -> Sumw2();
  _hEENET0_normal                -> Sumw2();
  _hEEneT0_weird                 -> Sumw2();
  _hEEneT0_normal                -> Sumw2();
  _hEtReco_weird                 -> Sumw2();
  _hEtReco_normal                -> Sumw2();
  _hEtSim_weird                  -> Sumw2();
  _hEtSim_normal                 -> Sumw2();
  _hEenToEEnRatio_weird          -> Sumw2();
  _hEenToEEnRatio_normal         -> Sumw2();
  _hNumTwrInClust_weird          -> Sumw2();
  _hNumTwrInClust_normal         -> Sumw2();
  _hEtaModule                    -> Sumw2();
  _hPhiModule                    -> Sumw2();
  _hTwrModule                    -> Sumw2();
  _hSecModule                    -> Sumw2();
  _hDiffPhiMod                   -> Sumw2();
  _hDiffTwrMod                   -> Sumw2();
  _hDiffSecMod                   -> Sumw2();
  _hSecondTwrID_weird            -> Sumw2();
  _hSecondTwrID_normal           -> Sumw2();
  _hSecondTwrEnergy_weird        -> Sumw2();
  _hSecondTwrEnergy_normal       -> Sumw2();
  _hSecondVsLeadTwrEnergy_weird  -> Sumw2();
  _hSecondVsLeadTwrEnergy_normal -> Sumw2();
  cout << "    Initialized histograms and profiles." << endl;

}  // end 'InitHists()'



void StTriggerEnergyScaleCalculator::SaveTreesAndHists() {

  // make directories
  TDirectory *dDet = (TDirectory*) _fOut -> mkdir("detector");
  TDirectory *dCan = (TDirectory*) _fOut -> mkdir("candidate");
  TDirectory *dMat = (TDirectory*) _fOut -> mkdir("match");
  TDirectory *dPar = (TDirectory*) _fOut -> mkdir("particle");
  TDirectory *dEff = (TDirectory*) _fOut -> mkdir("efficiency");
  TDirectory *dRes = (TDirectory*) _fOut -> mkdir("resolution");
  TDirectory *dChk = (TDirectory*) _fOut -> mkdir("checks");
  cout << "    Made directories." << endl;

  // save histograms
  dDet                   -> cd();
  _hTspTrg[0]            -> Write();
  _hTspTrg[1]            -> Write();
  _hPhiTrg[0]            -> Write();
  _hPhiTrg[1]            -> Write();
  _hEtaTrg[0]            -> Write();
  _hEtaTrg[1]            -> Write();
  _hEtTrg[0]             -> Write();
  _hEtTrg[1]             -> Write();
  dCan                   -> cd();
  _hPhiTrg[4]            -> Write();
  _hPhiTrg[5]            -> Write();
  _hEtaTrg[4]            -> Write();
  _hEtaTrg[5]            -> Write();
  _hEtTrg[4]             -> Write();
  _hEtTrg[5]             -> Write();
  _hDrTrg[0]             -> Write();
  _hDrTrg[1]             -> Write();
  _hQtTrg[0]             -> Write();
  _hQtTrg[1]             -> Write();
  _hDqTrg[0]             -> Write();
  _hDqTrg[1]             -> Write();
  _hQtParVsEtPar         -> Write();
  _hQtParVsEtDet         -> Write();
  _hDrParVsEtPar         -> Write();
  _hDrParVsEtDet         -> Write();
  _hDrParVsQtPar         -> Write();
  dMat                   -> cd();
  _hNumMatch             -> Write();
  _hTspTrg[2]            -> Write();
  _hTspTrg[3]            -> Write();
  _hTspMatch[0]          -> Write();
  _hTspMatch[1]          -> Write();
  _hTspMatch[2]          -> Write();
  _hTspMatchNW[0]        -> Write();
  _hTspMatchNW[1]        -> Write();
  _hTspMatchNW[2]        -> Write();
  _hPhiTrg[2]            -> Write();
  _hPhiTrg[3]            -> Write();
  _hPhiTrg[6]            -> Write();
  _hPhiTrg[7]            -> Write();
  _hEtaTrg[2]            -> Write();
  _hEtaTrg[3]            -> Write();
  _hEtaTrg[6]            -> Write();
  _hEtaTrg[7]            -> Write();
  _hEtTrg[2]             -> Write();
  _hEtTrg[3]             -> Write();
  _hEtTrg[6]             -> Write();
  _hEtTrg[7]             -> Write();
  _hEtDet[0]             -> Write();
  _hEtDet[1]             -> Write();
  _hEtDet[2]             -> Write();
  _hEtDetNW[0]           -> Write();
  _hEtDetNW[1]           -> Write();
  _hEtDetNW[2]           -> Write();
  _hEtMatch[0]           -> Write();
  _hEtMatch[1]           -> Write();
  _hEtMatch[2]           -> Write();
  _hEtMatchNW[0]         -> Write();
  _hEtMatchNW[1]         -> Write();
  _hEtMatchNW[2]         -> Write();
  _hDrTrg[2]             -> Write();
  _hDrTrg[3]             -> Write();
  _hDrMatch[0]           -> Write();
  _hDrMatch[1]           -> Write();
  _hDrMatch[2]           -> Write();
  _hDrMatchNW[0]         -> Write();
  _hDrMatchNW[1]         -> Write();
  _hDrMatchNW[2]         -> Write();
  _hQtTrg[2]             -> Write();
  _hQtTrg[3]             -> Write();
  _hQtMatch[0]           -> Write();
  _hQtMatch[1]           -> Write();
  _hQtMatch[2]           -> Write();
  _hQtMatchNW[0]         -> Write();
  _hQtMatchNW[1]         -> Write();
  _hQtMatchNW[2]         -> Write();
  _hDqTrg[2]             -> Write();
  _hDqTrg[3]             -> Write();
  _hDqMatch[0]           -> Write();
  _hDqMatch[1]           -> Write();
  _hDqMatch[2]           -> Write();
  _hDqMatchNW[0]         -> Write();
  _hDqMatchNW[1]         -> Write();
  _hDqMatchNW[2]         -> Write();
  dPar                   -> cd();
  _hPhiTrg[8]            -> Write();
  _hPhiTrg[9]            -> Write();
  _hEtaTrg[8]            -> Write();
  _hEtaTrg[9]            -> Write();
  _hEtTrg[8]             -> Write();
  _hEtTrg[9]             -> Write();
  dEff                   -> cd();
  _hPhiForEff[0]         -> Write();
  _hPhiForEff[1]         -> Write();
  _hEtaForEff[0]         -> Write();
  _hEtaForEff[1]         -> Write();
  _hEtForEff[0]          -> Write();
  _hEtForEff[1]          -> Write();
  _hPhiEff               -> Write();
  _hEtaEff               -> Write();
  _hEtEff                -> Write();
  _hNumVsCut             -> Write();
  _hRatioVsCut           -> Write();
  _hEtVsCut[0]           -> Write();
  _hEtVsCut[1]           -> Write();
  _hEtVsCut[2]           -> Write();
  _hEtVsCut[3]           -> Write();
  _hEtVsCut[4]           -> Write();
  _hEtVsCut[5]           -> Write();
  _hEtVsCut[6]           -> Write();
  dRes                   -> cd();
  _hPhiMatchVsDet        -> Write();
  _hPhiMatchVsDetNW      -> Write();
  _pPhiMatchVsDet        -> Write();
  _pPhiMatchVsDetNW      -> Write();
  _hEtaMatchVsDet        -> Write();
  _hEtaMatchVsDetNW      -> Write();
  _pEtaMatchVsDet        -> Write();
  _pEtaMatchVsDetNW      -> Write();
  _hEtMatchVsDet         -> Write();
  _hEtMatchVsDetNW       -> Write();
  _pEtMatchVsDet         -> Write();
  _pEtMatchVsDetNW       -> Write();
  _hEtMatchVsDetTrg[0]   -> Write();
  _pEtMatchVsDetTrg[0]   -> Write();
  _hEtMatchVsDetTrg[1]   -> Write();
  _pEtMatchVsDetTrg[1]   -> Write();
  _hEtMatchVsDetTrg[2]   -> Write();
  _pEtMatchVsDetTrg[2]   -> Write();
  _hEtMatchVsDetTrgNW[0] -> Write();
  _pEtMatchVsDetTrgNW[0] -> Write();
  _hEtMatchVsDetTrgNW[1] -> Write();
  _pEtMatchVsDetTrgNW[1] -> Write();
  _hEtMatchVsDetTrgNW[2] -> Write();
  _pEtMatchVsDetTrgNW[2] -> Write();
  _hDrMatchVsEtDet       -> Write();
  _hDrMatchVsEtDetNW     -> Write();
  _pDrMatchVsEtDet       -> Write();
  _pDrMatchVsEtDetNW     -> Write();
  _hQtMatchVsEtMatch     -> Write();
  _hQtMatchVsEtMatchNW   -> Write();
  _pQtMatchVsEtMatch     -> Write();
  _pQtMatchVsEtMatchNW   -> Write();
  _hQtMatchVsEtDet       -> Write();
  _hQtMatchVsEtDetNW     -> Write();
  _pQtMatchVsEtDet       -> Write();
  _pQtMatchVsEtDetNW     -> Write();
  _hQtInvMatchVsEtDet    -> Write();
  _hQtInvMatchVsEtDetNW  -> Write();
  _pQtInvMatchVsEtDet    -> Write();
  _pQtInvMatchVsEtDetNW  -> Write();
  _pQtMatchVsEtDet       -> Write();
  _pQtMatchVsEtDetNW     -> Write();
  _hDqMatchVsEtDet       -> Write();
  _hDqMatchVsEtDetNW     -> Write();
  _pDqMatchVsEtDet       -> Write();
  _pDqMatchVsEtDetNW     -> Write();
  _hQtEtMatchVsDet       -> Write();


  // trigger checks
  dChk                           -> cd();
  _hDidT_weird                   -> Write();
  _hDidT_normal                  -> Write();
  _hDidE_weird                   -> Write();
  _hDidE_normal                  -> Write();
  _hTrigID_weird                 -> Write();
  _hTrigID_normal                -> Write();
  _hModuleT_weird                -> Write();
  _hModuleT_normal               -> Write();
  _hEeneT0_weird                 -> Write();
  _hEeneT0_normal                -> Write();
  _hEENET0_weird                 -> Write();
  _hEENET0_normal                -> Write();
  _hEEneT0_weird                 -> Write();
  _hEEneT0_normal                -> Write();
  _hEtReco_weird                 -> Write();
  _hEtReco_normal                -> Write();
  _hEtSim_weird                  -> Write();
  _hEtSim_normal                 -> Write();
  _hEenToEEnRatio_weird          -> Write();
  _hEenToEEnRatio_normal         -> Write();
  _hNumTwrInClust_weird          -> Write();
  _hNumTwrInClust_normal         -> Write();
  _hEtaModule                    -> Write();
  _hPhiModule                    -> Write();
  _hTwrModule                    -> Write();
  _hSecModule                    -> Write();
  _hDiffPhiMod                   -> Write();
  _hDiffTwrMod                   -> Write();
  _hDiffSecMod                   -> Write();
  _hSecondTwrID_weird            -> Write();
  _hSecondTwrID_normal           -> Write();
  _hSecondTwrEnergy_weird        -> Write();
  _hSecondTwrEnergy_normal       -> Write();
  _hSecondVsLeadTwrEnergy_weird  -> Write();
  _hSecondVsLeadTwrEnergy_normal -> Write();
  cout << "    Saved histograms." << endl;

  // save trees
  _fOut    -> cd();
  _tJetPar -> Write();
  _tJetDet -> Write();
  cout << "    Saved output trees." << endl;

}  // end 'SaveTreesAndHists()'



void StTriggerEnergyScaleCalculator::DoDetFirstMatching() {

  // for jet-finding
  vector<PseudoJet> parTracks;
  vector<PseudoJet> detTracks;
  vector<PseudoJet> pJetsCS;
  vector<PseudoJet> dJetsCS;
  vector<PseudoJet> pJets;
  vector<PseudoJet> dJets;
  parTracks.clear();
  detTracks.clear();
  pJetsCS.clear();
  dJetsCS.clear();
  pJets.clear();
  dJets.clear();

  // for preventing multi-particle clusters
  Bool_t hasNeighbor[NTrkMax];
  for (UInt_t iTrk = 0; iTrk < NTrkMax; iTrk++) {
    hasNeighbor[iTrk] = false;
  }

  // begin event loops
  const UInt_t nParEvts = _tPar -> GetEntriesFast();
  const UInt_t nMatEvts = _tMat -> GetEntriesFast();
  const UInt_t nDetEvts = _tDet -> GetEntriesFast();
  if (nMatEvts != nDetEvts) {
    cerr << "PANIC: number of matched particle and detector events are different!\n"
         << "  nMatEvts = " << nMatEvts << "\n"
         << "  nDetEvts = " << nDetEvts
         << endl;
    return;
  } else {
    cout << "    Beginning event loops: " << nParEvts << " generated and " << nDetEvts << " reconstructed events to process." << endl;
  }

  // detector event loop
  UInt_t mBytes(0);
  UInt_t dBytes(0);
  UInt_t nMatBytes(0);
  UInt_t nDetBytes(0);
  UInt_t nTrgDet(0);
  UInt_t nTrgMat(0);
  UInt_t nVtxTrg(0);
  UInt_t nTwrTrg(0);
  UInt_t nPproj(0);
  UInt_t nEtaTrg(0);
  UInt_t nMatchTrg(0);
  UInt_t nEtTrg(0);
  UInt_t nTspTrg(0);

  int evt_old = 11111111;
  for (UInt_t iEvt = 0; iEvt < nDetEvts; iEvt++) {

    // load entry
    mBytes     = _tMat -> GetEntry(iEvt);
    dBytes     = _tDet -> GetEntry(iEvt);
    nMatBytes += mBytes;
    nDetBytes += dBytes;
    if ((mBytes < 0) || (dBytes < 0)) {
      cerr << "WARNING: issue with entry " << iEvt << "!\n"
           << "  no. of mat. bytes = " << mBytes << "\n"
           << "  no. of det. bytes = " << dBytes
           << endl;
      break;
    } else {
      if (_isInBatchMode) {
        cout << "      Processing detector event " << iEvt + 1 << "/" << nDetEvts << "..." << endl;
      } else {
        cout << "      Processing detector event " << iEvt + 1 << "/" << nDetEvts << "...\r" << flush;
        if ((iEvt + 1) == nMatEvts) cout << endl;
      }
    }

    // clear vectors and arrays
    parTracks.clear();
    detTracks.clear();
    pJetsCS.clear();
    dJetsCS.clear();
    pJets.clear();
    dJets.clear();

    double weight[NTrkMax];

    /*    for (UInt_t iTrk = 0; iTrk < NTrkMax; iTrk++) {
      hasNeighbor[iTrk] = false;
      }*/

    // particle event info
    const Int_t    runIdMat = _mMcRunId;
    const Int_t    evtIdMat = _mMcEventId;
    const UInt_t   nTrkMat  = _mMcNumTrks;
    const Double_t zVtxMat  = _mMcVz;
    const Double_t rVtxMat  = TMath::Sqrt((_mMcVx * _mMcVx) + (_mMcVy * _mMcVy));

    // detector event info
    const UInt_t   runIdDet = runNumber;
    const UInt_t   evtIdDet = eventNumber;
    const UInt_t   nTrkDet  = nPrimaryTracks;
    const UInt_t   nTwrDet  = TowerArray_;
    const Double_t zVtxDet  = zVertex;
    const Double_t rVtxDet  = TMath::Sqrt((xVertex * xVertex) + (yVertex * yVertex));

    // check if events are the same
    const Bool_t isSameRun = (runIdMat == runIdDet);
    const Bool_t isSameEvt = (evtIdMat == evtIdDet);
    if (!isSameRun || !isSameEvt) {
      cerr << "WARNING: matched particle and detector events have different run or event IDs!\n"
           << "  event number    = " << iEvt << "\n"
           << "  runId(Mat, Det) = (" << runIdMat << ", " << runIdDet << ")\n"
           << "  evtId(Mat, Det) = (" << evtIdMat << ", " << evtIdDet << ")"
           << endl;
      continue;
    }

    // new event 
    if (evtIdDet != evt_old) {

      // initialize weights and neighbor arrays	
      int pionVxId[NTrkMax];
      int indexPion[NTrkMax];
      int gammaVxId[NTrkMax];
      int indexGamma[NTrkMax];
      for (int iMcTrk = 0; iMcTrk < NTrkMax; iMcTrk++) {
        hasNeighbor[iMcTrk] = false;
        weight[iMcTrk]      = 0;
      }

      // particle trigger loop 1
      int Npions  = 0;
      int Ngammas = 0;
      for (UInt_t iTrgMC = 0; iTrgMC < nTrkMat; iTrgMC++) {

        // particle MC info
        const Int_t    pidTrgMc = _mMcIdGeant -> at(iTrgMC);
        const Double_t pTtrgMc  = _mMcPt      -> at(iTrgMC);

        // loop over trigger IDs
        for (UInt_t iTrgId = 0; iTrgId < NTrgIds; iTrgId++) {

          // if pi0, look for decays, otherwise record weights
          if (_idTrgUse[iTrgId] == 7) {

            // if pi0, record index and weight
	    if (pidTrgMc == 7) {
	      weight[iTrgMC]    = _fWeight    -> Eval(pTtrgMc);
	      pionVxId[Npions]  = _mMcIdVxEnd -> at(iTrgMC);
	      indexPion[Npions] = iTrgMC;
              Npions++;
	    }

            // if gamma/e-/e+, check if decay from pi0
	    if (pidTrgMc == 1 || pidTrgMc == 2 || pidTrgMc == 3) {
	      for (int iPion = 0; iPion < Npions; iPion++) {
                if (_mMcIdVx -> at(iTrgMC) == pionVxId[iPion]) {
		  weight[iTrgMC] = weight[(indexPion[iPion])];
                }
		if (pidTrgMc == 1) {
                  gammaVxId[Ngammas]  = _mMcIdVxEnd -> at(iTrgMC);
                  indexGamma[Ngammas] = iTrgMC;
                  Ngammas++;
                }
              }
            }

            // if e-/e+, check if decay from photon
            if (pidTrgMc == 2 || pidTrgMc == 3) {
              for (int iGamma = 0; iGamma < Ngammas; iGamma++) {
                if (_mMcIdVx -> at(iTrgMC) == gammaVxId[iGamma]) {
                  weight[iTrgMC] = weight[(indexGamma[iGamma])];
                }
              }
            }
          } else {
            if (pidTrgMc == _idTrgUse[iTrgId]) {
	      weight[iTrgMC] = _fWeight -> Eval(pTtrgMc);
            }
          }
        }  // end trigger ID loop
      }  // end particle loop
    }  // end if (isNewEvt)

    // filter out bad runs
    Bool_t isGoodRunMat(true);
    Bool_t isGoodRunDet(true);
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (runIdMat == _badRunList[iRun]) isGoodRunMat = false;
      if (runIdDet == _badRunList[iRun]) isGoodRunDet = false;
    }

    // detector event cuts
    const Bool_t isInZcutMat    = (TMath::Abs(zVtxMat) < _zVtxMax);
    const Bool_t isInZcutDet    = (TMath::Abs(zVtxDet) < _zVtxMax);
    const Bool_t isInRcutMat    = (TMath::Abs(rVtxMat) < _rVtxMax);
    const Bool_t isInRcutDet    = (TMath::Abs(rVtxDet) < _rVtxMax);
    const Bool_t isGoodEventMat = (isGoodRunMat && isInRcutMat && isInZcutMat);
    const Bool_t isGoodEventDet = (isGoodRunDet && isInRcutDet && isInZcutDet);
    if (_useParTrigger && !isGoodEventMat) continue;
    if (_useDetTrigger && !isGoodEventDet) continue;

    // detector trigger info
    const Int_t    adcDet       = ETwradc11;
    const UInt_t   idTrgDet     = ETwrdidT;
    const Double_t tspTrgDet    = Etsp;
    const Double_t tspTrgAltDet = EtspAlt;
    const Double_t eStrEta      = EEstrpen4;
    const Double_t eStrPhi      = EPstripenp4;
    const Double_t pProjDet     = ETwrPTower;
    const Double_t hTrgDet      = ETwreT;
    const Double_t hPhysDet     = EClustetav1;
    const Double_t fTrgDet      = EClustphiv1;
    const Double_t eTrgDet      = EClustEneT0;
    const Double_t tTrgDet      = 2. * TMath::ATan(TMath::Exp(-1. * hPhysDet));
    const Double_t eTtrgDet     = eTrgDet * TMath::Sin(tTrgDet);

    // filter out bad towers
    Bool_t isGoodTwrDet = true;
    for (UInt_t iTwr = 0; iTwr < NHotTwrs; iTwr++) {
      if (idTrgDet == _hotTwrList[iTwr]) {
        isGoodTwrDet = false;
        break;
      }
    }

    // detector trigger cuts
    const Bool_t isInAdcCutDet    = (adcDet <= _adcMax);
    const Bool_t isInStrCutDet    = ((eStrEta >= _eStrMin) && (eStrPhi >= _eStrMin));
    const Bool_t isInProjCutDet   = (pProjDet < _pProjMax);
    const Bool_t isInEtaTrgCutDet = ((TMath::Abs(hTrgDet) < _hTrgDetMax) && (TMath::Abs(hPhysDet) < _hTrgDetMax));
    const Bool_t isInEtCutDet     = ((eTtrgDet >= _eTtrgDetMin) && (eTtrgDet < _eTtrgDetMax));
    const Bool_t isInTspCutDet    = ((tspTrgDet > _tspUse[0]) && (tspTrgDet < _tspUse[1]));
    const Bool_t isGoodTrgDet     = (isGoodTwrDet && isInAdcCutDet && isInStrCutDet && isInProjCutDet && isInEtaTrgCutDet);

    // fill efficiency histograms
    _hEtVsCut[0] -> Fill(eTtrgDet);
    ++nVtxTrg;
    if (isGoodTwrDet) {
      _hEtVsCut[1] -> Fill(eTtrgDet);
      ++nTwrTrg;
      if (isInProjCutDet) {
        _hEtVsCut[2] -> Fill(eTtrgDet);
        ++nPproj;
        if (isInEtaTrgCutDet) {
          _hEtVsCut[3] -> Fill(eTtrgDet);
          ++nEtaTrg;
        }  // eta cut
      }  // pProj cut
    }  // bad tower cut

    // apply cuts
    if (_useDetTrigger && !isGoodTrgDet) {
      continue;
    } else {
      // fill unweighted detector histograms
      _hTspTrg[0] -> Fill(tspTrgDet);
      _hPhiTrg[0] -> Fill(fTrgDet);
      _hEtaTrg[0] -> Fill(hTrgDet);
      _hEtTrg[0]  -> Fill(eTtrgDet);
      ++nTrgDet;
    }

    // particle trigger loop 1
    for (UInt_t iTrgMC = 0; iTrgMC < nTrkMat; iTrgMC++) {

      // particle trigger info
      const Int_t    pidTrgMc = _mMcIdGeant -> at(iTrgMC);
      const Double_t cTrgMc   = _mMcCharge  -> at(iTrgMC);
      const Double_t fTrgMc   = _mMcPhi     -> at(iTrgMC);
      const Double_t hTrgMc   = _mMcEta     -> at(iTrgMC);
      const Double_t pTtrgMc  = _mMcPt      -> at(iTrgMC);

      // check pid
      Bool_t isHadronMat(false);
      for (UInt_t iParID = 0; iParID < NTrgIds; iParID++) {
        Bool_t isSameId = (pidTrgMc == _idTrgUse[iParID]);
        if (isSameId) {
          isHadronMat = true;
          break;
        }
      }

      // particle trigger cuts
      const Bool_t isNeutralTrgMat  = (cTrgMc == _cTrgPar);
      const Bool_t isInTrgEtaCutMat = (TMath::Abs(hTrgMc) < _hTrgParMax);
      const Bool_t isInTrgEtCutMat  = ((pTtrgMc > _eTtrgMatMin) && (pTtrgMc < _eTtrgMatMax));
      if (isHadronMat && isNeutralTrgMat && isInTrgEtaCutMat && isInTrgEtCutMat) {

	int index_daughter[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	int numDaughters = 0;
        //initialize
	if (hasNeighbor[iTrgMC]) continue;
	//	hasNeighbor[iTrgMC] = false;
        // check for neighbors
        for (UInt_t iCheck = 0; iCheck < nTrkMat; iCheck++) {

          // check if same particle
          const Bool_t isSameParticle = (iTrgMC == iCheck);
          if (isSameParticle) continue;
	  
	  // check for daughter particles (i.e. skip if from same parent)
	  if (_mMcIdVx -> at(iCheck) == _mMcIdVxEnd -> at(iTrgMC))
	    {
	      //this is a direct daughter of MC particle checking for neighbors
	      index_daughter[numDaughters] = iCheck;
	      numDaughters++;
	      continue;
	    }
	  int daughterOfDaughter = 0;
	  if (numDaughters>0) // at least one daughter particle was found
	    {
	      for (int iDaughter=0; iDaughter<numDaughters; iDaughter++)
		{
		  if (_mMcIdVx -> at(iCheck) == _mMcIdVxEnd -> at(index_daughter[iDaughter]))
		// daughter of daughter , e.g. pi0 -> gamma gamma -> gamma e+ e-
		    daughterOfDaughter = 1;
		}
	      if (daughterOfDaughter > 0) continue;
	    }

          // possible neighbor info
	  const Bool_t CheckIsShower =   _mMcIsShower      -> at(iCheck); // allow for geant particle from interactions with material
	  if (CheckIsShower) continue;  //do not count shower particles as neighbors

          const Double_t pidCheck = _mMcIdGeant -> at(iCheck);
          const Double_t cCheck   = _mMcCharge  -> at(iCheck);
          const Double_t fCheck   = _mMcPhi     -> at(iCheck);
          Double_t hCheck   = _mMcEta     -> at(iCheck);
          Double_t Echeck  = _mMcEnergy      -> at(iCheck);
	  int check_daughter = 0;
	  if (hCheck == -999) { hCheck = TMath::ATanH(_mMcPz -> at(iCheck) / _mMcPtot -> at(iCheck)); /*check_daughter = 1;  cout << "recalculate eta " << hCheck<< endl;*/}

          // calculate delta-r
          Double_t dFcheck = fCheck - fTrgMc;
          if (dFcheck > TMath::Pi()) {
            dFcheck -= TMath::TwoPi();
          }
          if (dFcheck < (-1. * TMath::Pi())) {
            dFcheck += TMath::TwoPi();
          }
          Double_t dHcheck = hCheck - hTrgMc;
          Double_t dRcheck = TMath::Sqrt((dFcheck * dFcheck) + (dHcheck * dHcheck));
	  
	  // const Bool_t isInCheckPtCut    = (pTcheck > 0.);

          const Bool_t isInMinSeparation = (dRcheck < _dMinSeparationUse);
	  //          if (isInCheckPtCut && isInMinSeparation) {
          if (isInMinSeparation) {
            hasNeighbor[iTrgMC] = true;
	    //  cout << "pid (with vxend = " << pidTrgMc << "(" << _mMcIdVxEnd -> at(iTrgMC) << ") has neighbor within " << dRcheck << ", which is a pi = " << pidCheck << "with vx = " <<_mMcIdVx -> at(iCheck)<< endl;
	    //            hasNeighbor[iCheck] = true;
            break;
          }
        }  // end neighbor loop
      }
    }  // end particle trigger loop 1

    // particle trigger loop 2
    Int_t    idTrgMat(-1);
    Bool_t   foundTrgMat(false);
    Bool_t   foundMatch(false);
    UInt_t   iTrgMat(nTrkMat + 1);
    UInt_t   nMatch(0);
    Double_t fTrgMat(0.);
    Double_t hTrgMat(0.);
    Double_t eTtrgMat(0.);
    Double_t dRtrgMat(0.);
    Double_t qTtrgMat(-1.);
    Double_t qTinvMat(-1.);
    Double_t dQtTrgMat(-1.);
    Double_t wTrgMat(-1.);
    for (UInt_t iTrgMC = 0; iTrgMC < nTrkMat; iTrgMC++) {

      // don't consider if particle has neighbor
      if (_excludeNeighbors && hasNeighbor[iTrgMC]) continue;

      // particle trigger info
      const Int_t    pidTrgMc = _mMcIdGeant -> at(iTrgMC);
      const Double_t cTrgMc   = _mMcCharge  -> at(iTrgMC);
      const Double_t fTrgMc   = _mMcPhi     -> at(iTrgMC);
      const Double_t hTrgMc   = _mMcEta     -> at(iTrgMC);
      const Double_t pTtrgMc  = _mMcPt      -> at(iTrgMC);
      //      const Double_t wTrgMc   = _fWeight     -> Eval(pTtrgMc);
      //  const Double_t wTrgMc   = 0; //weight[iTrgMC];
  const Double_t wTrgMc   = weight[iTrgMC];

      // delta-r and qT relative to detector trigger
      Double_t dFtrgMc = fTrgMc - fTrgDet;
      Double_t dHtrgMc = hTrgMc - hPhysDet;
      if (dFtrgMc > TMath::Pi()) {
        dFtrgMc -= TMath::TwoPi();
      }
      if (dFtrgMc < (-1. * TMath::Pi())) {
        dFtrgMc += TMath::TwoPi();
      }

      const Double_t dRtrgMc  = TMath::Sqrt((dFtrgMc * dFtrgMc) + (dHtrgMc * dHtrgMc));
      const Double_t qTtrgMc  = eTtrgDet / pTtrgMc;
      const Double_t qTinvMc  = pTtrgMc / eTtrgDet;
      const Double_t dQtTrgMc = TMath::Abs(eTtrgDet - pTtrgMc) / pTtrgMc;

      // check pid
      Bool_t isHadronMat(false);
      for (UInt_t iParID = 0; iParID < NTrgIds; iParID++) {
        Bool_t isSameId = (pidTrgMc == _idTrgUse[iParID]);
        if (isSameId) {
          isHadronMat = true;
          break;
        }
      }

      // particle trigger cuts
      const Bool_t isNeutralTrgMat  = (cTrgMc == _cTrgPar);
      const Bool_t isInTrgDrCutMat  = (dRtrgMc < _dRtrgMaxUse);
      const Bool_t isInTrgEtaCutMat = (TMath::Abs(hTrgMc) < _hTrgParMax);
      const Bool_t isInTrgEtCutMat  = ((pTtrgMc > _eTtrgMatMin) && (pTtrgMc < _eTtrgMatMax));
      if (isHadronMat && isNeutralTrgMat && isInTrgEtaCutMat && isInTrgEtCutMat) {

        // fill candidate histograms
        _hPhiTrg[4]    -> Fill(fTrgMc, weight[iTrgMC]); 
        _hPhiTrg[5]    -> Fill(fTrgMc);
        _hEtaTrg[4]    -> Fill(hTrgMc, weight[iTrgMC]);
        _hEtaTrg[5]    -> Fill(hTrgMc);
        _hEtTrg[4]     -> Fill(pTtrgMc, weight[iTrgMC]);
        _hEtTrg[5]     -> Fill(pTtrgMc);
        _hDrTrg[0]     -> Fill(dRtrgMc, weight[iTrgMC]);
        _hDrTrg[1]     -> Fill(dRtrgMc);
        _hQtTrg[0]     -> Fill(qTtrgMc, weight[iTrgMC]);
        _hQtTrg[1]     -> Fill(qTtrgMc);
        _hDqTrg[0]     -> Fill(dQtTrgMc, weight[iTrgMC]);
        _hDqTrg[1]     -> Fill(dQtTrgMc);
        _hQtParVsEtPar -> Fill(pTtrgMc, qTtrgMc, weight[iTrgMC]);
        _hQtParVsEtDet -> Fill(eTtrgDet, qTtrgMc, weight[iTrgMC]);
        _hDrParVsEtPar -> Fill(pTtrgMc, dRtrgMc, weight[iTrgMC]);
        _hDrParVsEtDet -> Fill(eTtrgDet, dRtrgMc, weight[iTrgMC]);
        _hDrParVsQtPar -> Fill(qTtrgMc, dRtrgMc, weight[iTrgMC]);

        // record trigger info
        iTrgMat     = iTrgMC;
        fTrgMat     = fTrgMc;
        hTrgMat     = hTrgMc;
        eTtrgMat    = pTtrgMc;
        idTrgMat    = pidTrgMc;
        dRtrgMat    = dRtrgMc;
        qTtrgMat    = qTtrgMc;
        qTinvMat    = qTinvMc;
        dQtTrgMat   = dQtTrgMc;
        wTrgMat     = weight[iTrgMC];
        foundTrgMat = true;

        // apply matching condition
        if (isInTrgDrCutMat) {
          foundMatch = true;
          ++nMatch;
          break;
        }  // end matching condition
      }  // end particle trigger cuts
    }  // end particle trigger loop 2

    // determine if particle trigger is good or not
    Bool_t isGoodMatch  = !(_requireMatch && !foundMatch);
    Bool_t isGoodMatTrg = !(_useParTrigger && !foundTrgMat);
    Bool_t isGoodMatLvl = (isGoodMatTrg && isGoodMatch);

    // fill unmatched histograms
    if (isGoodMatTrg && !foundMatch) {
      _hTspTrg[1] -> Fill(tspTrgDet);
      _hPhiTrg[1] -> Fill(fTrgDet);
      _hEtaTrg[1] -> Fill(hTrgDet);
      _hEtTrg[1]  -> Fill(eTtrgDet);
    }

    // fill matched histograms
    if (isGoodMatTrg && foundMatch) {

      //      wTrgMat   =  weight[iTrgMC]; 

      _hEtVsCut[4] -> Fill(eTtrgDet);
      ++nMatchTrg;
      if (isInEtCutDet) {
        _hEtVsCut[5] -> Fill(eTtrgDet);
        ++nEtTrg;
        if (isInTspCutDet) {
          _hEtVsCut[6] -> Fill(eTtrgDet);
          ++nTspTrg;
        }  // tsp cut
      }  // eTdet cut

      if (isInEtCutDet && isInTspCutDet) {

        // fill weighted detector histograms
        _hTspTrg[2]      -> Fill(tspTrgDet, wTrgMat);
        _hTspTrg[3]      -> Fill(tspTrgDet);
        _hPhiTrg[2]      -> Fill(fTrgDet, wTrgMat);
        _hPhiTrg[3]      -> Fill(fTrgDet);
        _hEtaTrg[2]      -> Fill(hTrgDet, wTrgMat);
        _hEtaTrg[3]      -> Fill(hTrgDet);
        _hEtTrg[2]       -> Fill(eTtrgDet, wTrgMat);
        _hEtTrg[3]       -> Fill(eTtrgDet);
        _hNumMatch       -> Fill(nMatch, wTrgMat);
        // fill matched particle histograms 
        _hPhiTrg[6]      -> Fill(fTrgMat, wTrgMat);
        _hPhiTrg[7]      -> Fill(fTrgMat);
        _hEtaTrg[6]      -> Fill(hTrgMat, wTrgMat);
        _hEtaTrg[7]      -> Fill(hTrgMat);
        _hEtTrg[6]       -> Fill(eTtrgMat, wTrgMat);
        _hEtTrg[7]       -> Fill(eTtrgMat);
        _hDrTrg[2]       -> Fill(dRtrgMat, wTrgMat);
        _hDrTrg[3]       -> Fill(dRtrgMat);
        _hQtTrg[2]       -> Fill(qTtrgMat, wTrgMat);
        _hQtTrg[3]       -> Fill(qTtrgMat);
        _hDqTrg[2]       -> Fill(dQtTrgMat, wTrgMat);
        _hDqTrg[3]       -> Fill(dQtTrgMat);
        _hPhiForEff[0]   -> Fill(fTrgMat, wTrgMat);
        _hEtaForEff[0]   -> Fill(hTrgMat, wTrgMat);
        _hEtForEff[0]    -> Fill(eTtrgMat, wTrgMat);
        // fill 2d/3d histograms
        _hPhiMatchVsDet       -> Fill(fTrgDet, fTrgMat, wTrgMat);
        _hPhiMatchVsDetNW     -> Fill(fTrgDet, fTrgMat);
        _pPhiMatchVsDet       -> Fill(fTrgDet, fTrgMat, wTrgMat);
        _pPhiMatchVsDetNW     -> Fill(fTrgDet, fTrgMat);
        _hEtaMatchVsDet       -> Fill(hTrgDet, hTrgMat, wTrgMat);
        _hEtaMatchVsDetNW     -> Fill(hTrgDet, hTrgMat);
        _pEtaMatchVsDet       -> Fill(hTrgDet, hTrgMat, wTrgMat);
        _pEtaMatchVsDetNW     -> Fill(hTrgDet, hTrgMat);
        _hEtMatchVsDet        -> Fill(eTtrgDet, eTtrgMat, wTrgMat);
        _hEtMatchVsDetNW      -> Fill(eTtrgDet, eTtrgMat);
        _pEtMatchVsDet        -> Fill(eTtrgDet, eTtrgMat, wTrgMat);
        _pEtMatchVsDetNW      -> Fill(eTtrgDet, eTtrgMat);
        _hDrMatchVsEtDet      -> Fill(eTtrgDet, dRtrgMat, wTrgMat);
        _hDrMatchVsEtDetNW    -> Fill(eTtrgDet, dRtrgMat);
        _pDrMatchVsEtDet      -> Fill(eTtrgDet, dRtrgMat, wTrgMat);
        _pDrMatchVsEtDetNW    -> Fill(eTtrgDet, dRtrgMat);
        _hQtMatchVsEtMatch    -> Fill(eTtrgMat, qTtrgMat, wTrgMat);
        _hQtMatchVsEtMatchNW  -> Fill(eTtrgMat, qTtrgMat);
        _pQtMatchVsEtMatch    -> Fill(eTtrgMat, qTtrgMat, wTrgMat);
        _pQtMatchVsEtMatchNW  -> Fill(eTtrgMat, qTtrgMat);
        _hQtMatchVsEtDet      -> Fill(eTtrgDet, qTtrgMat, wTrgMat);
        _hQtMatchVsEtDetNW    -> Fill(eTtrgDet, qTtrgMat);
        _pQtMatchVsEtDet      -> Fill(eTtrgDet, qTtrgMat, wTrgMat);
        _pQtMatchVsEtDetNW    -> Fill(eTtrgDet, qTtrgMat);
        _hQtInvMatchVsEtDet   -> Fill(eTtrgDet, qTinvMat, wTrgMat);
        _hQtInvMatchVsEtDetNW -> Fill(eTtrgDet, qTinvMat);
        _pQtInvMatchVsEtDet   -> Fill(eTtrgDet, qTinvMat, wTrgMat);
        _pQtInvMatchVsEtDetNW -> Fill(eTtrgDet, qTinvMat);
        _hDqMatchVsEtDet      -> Fill(eTtrgDet, dQtTrgMat, wTrgMat);
        _hDqMatchVsEtDetNW    -> Fill(eTtrgDet, dQtTrgMat);
        _pDqMatchVsEtDet      -> Fill(eTtrgDet, dQtTrgMat, wTrgMat);
        _pDqMatchVsEtDetNW    -> Fill(eTtrgDet, dQtTrgMat);
        _hQtEtMatchVsDet      -> Fill(eTtrgDet, eTtrgMat, qTtrgMat);

        // trigger checks
        _hEtaModule  -> Fill(EEstrpModuleE);
        _hPhiModule  -> Fill(EPstripModuleP);
        _hTwrModule  -> Fill(ETwrmoduleT);
        _hDiffPhiMod -> Fill(EEstrpModuleE - EPstripModuleP);
        _hDiffTwrMod -> Fill(EEstrpModuleE - ETwrmoduleT);
        if (EClustNumTwrIncluded > 1) {
          _hSecModule  -> Fill(EClustSecondTwrModule);
          _hDiffSecMod -> Fill(EEstrpModuleE - EClustSecondTwrModule);
        }

        if (qTtrgMat > 1.5) {
          _hDidT_weird                  -> Fill(ETwrdidT);
          _hDidE_weird                  -> Fill(ETwrdidE);
          _hTrigID_weird                -> Fill(trigID);
          _hModuleT_weird               -> Fill(ETwrmoduleT);
          _hEeneT0_weird                -> Fill(ETwreneT0);
          _hEENET0_weird                -> Fill(ETwrENET0);
          _hEEneT0_weird                -> Fill(EClustEneT0);
          _hEtReco_weird                -> Fill(eTtrgDet);
          _hEtSim_weird                 -> Fill(eTtrgMat);
          _hEenToEEnRatio_weird         -> Fill(ETwreneT0 / EClustEneT0);
          _hNumTwrInClust_weird         -> Fill(EClustNumTwrIncluded);
          _hSecondTwrID_weird           -> Fill(EClustSecondTwrIndex);
          _hSecondTwrEnergy_weird       -> Fill(EClustSecondTwrEnergy);
          _hSecondVsLeadTwrEnergy_weird -> Fill(ETwreneT0, EClustSecondTwrEnergy);
        } else {
          _hDidT_normal                  -> Fill(ETwrdidT);
          _hDidE_normal                  -> Fill(ETwrdidE);
          _hTrigID_normal                -> Fill(trigID);
          _hModuleT_normal               -> Fill(ETwrmoduleT);
          _hEeneT0_normal                -> Fill(ETwreneT0);
          _hEENET0_normal                -> Fill(ETwrENET0);
          _hEEneT0_normal                -> Fill(EClustEneT0);
          _hEtReco_normal                -> Fill(eTtrgDet);
          _hEtSim_normal                 -> Fill(eTtrgMat);
          _hEenToEEnRatio_normal         -> Fill(ETwreneT0 / EClustEneT0);
          _hNumTwrInClust_normal         -> Fill(EClustNumTwrIncluded);
          _hSecondTwrID_normal           -> Fill(EClustSecondTwrIndex);
          _hSecondTwrEnergy_normal       -> Fill(EClustSecondTwrEnergy);
          _hSecondVsLeadTwrEnergy_normal -> Fill(ETwreneT0, EClustSecondTwrEnergy);
        }
      }  // end if (isInEtCutDet && isInTspCutDet)

      // fill eTtrg-dependent histograms
      for (UInt_t iEtBin = 0; iEtBin < NTrgBin; iEtBin++) {
        if (isInTspCutDet && ((eTtrgDet > _eTbinMin[iEtBin]) && (eTtrgDet < _eTbinMax[iEtBin]))) {
          _hTspMatch[iEtBin]          -> Fill(tspTrgDet, wTrgMat);
          _hTspMatchNW[iEtBin]        -> Fill(tspTrgDet);
          _hEtDet[iEtBin]             -> Fill(eTtrgDet, wTrgMat);
          _hEtDetNW[iEtBin]           -> Fill(eTtrgDet);
          _hEtMatch[iEtBin]           -> Fill(eTtrgMat, wTrgMat);
          _hEtMatchNW[iEtBin]         -> Fill(eTtrgMat);
          _hDrMatch[iEtBin]           -> Fill(dRtrgMat, wTrgMat);
          _hDrMatchNW[iEtBin]         -> Fill(dRtrgMat);
          _hQtMatch[iEtBin]           -> Fill(qTtrgMat, wTrgMat);
          _hQtMatchNW[iEtBin]         -> Fill(qTtrgMat);
          _hDqMatch[iEtBin]           -> Fill(dQtTrgMat, wTrgMat);
          _hDqMatchNW[iEtBin]         -> Fill(dQtTrgMat);
          _hEtMatchVsDetTrg[iEtBin]   -> Fill(eTtrgDet, eTtrgMat, wTrgMat);
          _hEtMatchVsDetTrgNW[iEtBin] -> Fill(eTtrgDet, eTtrgMat);
          _pEtMatchVsDetTrg[iEtBin]   -> Fill(eTtrgDet, eTtrgMat, wTrgMat);
          _pEtMatchVsDetTrgNW[iEtBin] -> Fill(eTtrgDet, eTtrgMat);
        }
      }  // end eTtrg loop
    }  // end if (isGoodTrg && foundMatch)

    // apply particle-level cuts
    if (!isGoodMatLvl) {
      continue;
    } else {
      ++nTrgMat;
    }

    // detector track loop
    UInt_t nDetQA(0);
    for (UInt_t iTrkDet = 0; iTrkDet < nTrkDet; iTrkDet++) {
      // detector track info
      const UInt_t   nFitDet  = PrimaryTrackArray_nHitsFit[iTrkDet];
      const UInt_t   nPossDet = PrimaryTrackArray_nHitsPoss[iTrkDet];
      const Double_t rFitDet  = (Double_t) nFitDet / (Double_t) nPossDet;
      const Double_t dcaDet   = PrimaryTrackArray_dcag[iTrkDet];
      const Double_t hDet     = PrimaryTrackArray_eta[iTrkDet];
      const Double_t pXdet    = PrimaryTrackArray_pX[iTrkDet];
      const Double_t pYdet    = PrimaryTrackArray_pY[iTrkDet];
      const Double_t pZdet    = PrimaryTrackArray_pZ[iTrkDet];
      const Double_t pTdet    = PrimaryTrackArray_pT[iTrkDet];
      const Double_t eDet     = TMath::Sqrt((pTdet * pTdet) + (pZdet * pZdet) + (MassPi0 * MassPi0));

      // detector track cuts
      const Bool_t isInFitCutDet   = (nFitDet > _nFitTrkMin);
      const Bool_t isInRatioCutDet = (rFitDet > _rFitTrkMin);
      const Bool_t isInDcaCutDet   = (dcaDet < _dcaTrkMax);
      const Bool_t isInEtaCutDet   = (TMath::Abs(hDet) < _hTrkMax);
      const Bool_t isInPtCutDet    = ((pTdet > _pTtrkDetMin) && (pTdet < _pTtrkDetMax));
      if (isInFitCutDet && isInRatioCutDet && isInDcaCutDet && isInEtaCutDet && isInPtCutDet) {
        detTracks.push_back(PseudoJet(pXdet, pYdet, pZdet, eDet));
        ++nDetQA;
      }
    }  // end detector track loop

    // particle track loop
    UInt_t nMatQA(0);
    for (UInt_t iTrkMat = 0; iTrkMat < nTrkMat; iTrkMat++) {
      // particle track info
      const Double_t qMC  = _mMcCharge -> at(iTrkMat);
      const Double_t hMC  = _mMcEta    -> at(iTrkMat);
      const Double_t pTmc = _mMcPt     -> at(iTrkMat);
      const Double_t pXmc = _mMcPx     -> at(iTrkMat);
      const Double_t pYmc = _mMcPy     -> at(iTrkMat);
      const Double_t pZmc = _mMcPz     -> at(iTrkMat);
      const Double_t eMC  = _mMcEnergy -> at(iTrkMat);

      // particle track cuts
      const Bool_t isGoodEta   = (TMath::Abs(hMC) < _hTrkMax);
      const Bool_t isGoodChrg  = (qMC != 0.);
      const Bool_t isGoodPt    = ((pTmc > _pTtrkParMin) && (pTmc < _pTtrkParMax));
      const Bool_t isNotTrgMat = (iTrkMat != iTrgMat);
      if (isGoodChrg && isGoodEta && isGoodPt && isNotTrgMat) {
          parTracks.push_back(PseudoJet(pXmc, pYmc, pZmc, eMC));
          ++nMatQA;
      }
    }  // end particle track loop

    // define jets and jet area
    GhostedAreaSpec areaSpecP(_etaGhostMax, _nRepeat, _aGhost);
    GhostedAreaSpec areaSpecD(_etaGhostMax, _nRepeat, _aGhost);
    AreaDefinition  areaDefP(active_area_explicit_ghosts, areaSpecP);
    AreaDefinition  areaDefD(active_area_explicit_ghosts, areaSpecD);
    JetDefinition   jetDefP(antikt_algorithm, _rJet);
    JetDefinition   jetDefD(antikt_algorithm, _rJet);

    // cluster jets
    ClusterSequenceArea clusterSeqP(parTracks, jetDefP, areaDefP);
    ClusterSequenceArea clusterSeqD(detTracks, jetDefD, areaDefD);
    pJetsCS = sorted_by_pt(clusterSeqP.inclusive_jets(_pTjetMin));
    dJetsCS = sorted_by_pt(clusterSeqD.inclusive_jets(_pTjetMin));

    // fiducial cut
    Selector fidCut = SelectorAbsEtaMax(_etaJetMax);
    pJets = fidCut(pJetsCS);
    dJets = fidCut(dJetsCS);

    // define bkgd jets and jet area
    GhostedAreaSpec areaSpecBkgdP(_etaGhostMax, _nRepeat, _aGhost);
    GhostedAreaSpec areaSpecBkgdD(_etaGhostMax, _nRepeat, _aGhost);
    AreaDefinition  areaDefBkgdP(active_area_explicit_ghosts, areaSpecBkgdP);
    AreaDefinition  areaDefBkgdD(active_area_explicit_ghosts, areaSpecBkgdD);
    JetDefinition   jetDefBkgdP(kt_algorithm, _rJet);
    JetDefinition   jetDefBkgdD(kt_algorithm, _rJet);

    // initialize bkgd estimators
    Selector bkgdCut = SelectorAbsEtaMax(_etaBkgdMax) * (!SelectorNHardest(_nRemove));
    JetMedianBackgroundEstimator bkgdP(bkgdCut, jetDefBkgdP, areaDefBkgdP);
    JetMedianBackgroundEstimator bkgdD(bkgdCut, jetDefBkgdD, areaDefBkgdD);
    Subtractor subP(&bkgdP);
    Subtractor subD(&bkgdD);

#if FASTJET_VERSION_NUMBER >= 30100
  subP.set_use_rho_m(true);
  subD.set_use_rho_m(true);
  subP.set_safe_mass(true);
  subD.set_safe_mass(true);
#endif

    // estimate bkgd
    bkgdP.set_particles(parTracks);
    bkgdD.set_particles(detTracks);
    const Double_t pRhoJet = bkgdP.rho();
    const Double_t dRhoJet = bkgdD.rho();
    const Double_t pSigJet = bkgdP.sigma();
    const Double_t dSigJet = bkgdD.sigma();

    // clear jet info
    _pJetIndex.clear();
    _pJetNCons.clear();
    _pJetPt.clear();
    _pJetPtCorr.clear();
    _pJetEta.clear();
    _pJetPhi.clear();
    _pJetE.clear();
    _pJetArea.clear();
    _pJetConsPt.clear();
    _pJetConsEta.clear();
    _pJetConsPhi.clear();
    _pJetConsE.clear();

    _dJetIndex.clear();
    _dJetNCons.clear();
    _dJetPt.clear();
    _dJetPtCorr.clear();
    _dJetEta.clear();
    _dJetPhi.clear();
    _dJetE.clear();
    _dJetArea.clear();
    _dJetConsPt.clear();
    _dJetConsEta.clear();
    _dJetConsPhi.clear();
    _dJetConsE.clear();

    // jet event info
    _pEventIndex = evtIdMat;
    _pNJets      = (Int_t) pJets.size();
    _pRunId      = runIdMat;
    _pRefmult    = nMatQA;
    _pTSP        = idTrgMat;
    _pTrgEta     = hTrgMat;
    _pTrgPhi     = fTrgMat;
    _pTrgEtPar   = eTtrgMat;
    _pTrgEtDet   = eTtrgDet;
    _pRho        = pRhoJet;
    _pSigma      = pSigJet;
    _pVz         = zVtxMat;

    _dEventIndex = evtIdDet;
    _dNJets      = (Int_t) dJets.size();
    _dRunId      = runIdDet;
    _dRefmult    = nDetQA;
    _dTSP        = Etsp;
    _dTrgEta     = hTrgDet;
    _dTrgPhi     = fTrgDet;
    _dTrgEtPar   = eTtrgMat;
    _dTrgEtDet   = eTtrgDet;
    _dRho        = dRhoJet;
    _dSigma      = dSigJet;
    _dVz         = zVtxDet;

    // particle jet loop
    const UInt_t nParJets = (UInt_t) pJets.size();
    for (UInt_t iJetPar = 0; iJetPar < nParJets; iJetPar++) {
      const Int_t    nCstPar   = (Int_t) pJets[iJetPar].constituents().size();
      const Double_t hJetPar   = pJets[iJetPar].pseudorapidity();
      const Double_t fJetPar   = pJets[iJetPar].phi_std();
      const Double_t aJetPar   = pJets[iJetPar].area();
      const Double_t eJetPar   = pJets[iJetPar].e();
      const Double_t pTjetPar  = pJets[iJetPar].perp();
      const Double_t pTcorrPar = pTjetPar - (aJetPar * pRhoJet);
      _pJetIndex.push_back(iJetPar);
      _pJetNCons.push_back(nCstPar);
      _pJetPt.push_back(pTjetPar);
      _pJetPtCorr.push_back(pTcorrPar);
      _pJetEta.push_back(hJetPar);
      _pJetPhi.push_back(fJetPar);
      _pJetE.push_back(eJetPar);
      _pJetArea.push_back(aJetPar);
    }  // end particle jet loop

    // detector jet loop
    const UInt_t nDetJets = (UInt_t) dJets.size();
    for (UInt_t iJetDet = 0; iJetDet < nDetJets; iJetDet++) {
      const Int_t    nCstDet   = (Int_t) dJets[iJetDet].constituents().size();
      const Double_t hJetDet   = dJets[iJetDet].pseudorapidity();
      const Double_t fJetDet   = dJets[iJetDet].phi_std();
      const Double_t aJetDet   = dJets[iJetDet].area();
      const Double_t eJetDet   = dJets[iJetDet].e();
      const Double_t pTjetDet  = dJets[iJetDet].perp();
      const Double_t pTcorrDet = pTjetDet - (aJetDet * dRhoJet);
      _dJetIndex.push_back(iJetDet);
      _dJetNCons.push_back(nCstDet);
      _dJetPt.push_back(pTjetDet);
      _dJetPtCorr.push_back(pTcorrDet);
      _dJetEta.push_back(hJetDet);
      _dJetPhi.push_back(fJetDet);
      _dJetE.push_back(eJetDet);
      _dJetArea.push_back(aJetDet);
    }  // end particle jet loop

    // fill output trees
    _tJetPar -> Fill();
    _tJetDet -> Fill();

    evt_old = evtIdDet;

  }  // end detector event loop
  cout << "    Detector event loop finished:\n"
       << "      nTrgDet (to match) = " << nTrgDet << "\n"
       << "      nTrgPar (matched)  = " << nTrgMat << "\n"
       << "    Beginning particle event loop..."
       << endl;

  // set nTrg vs. cut values
  _hNumVsCut -> SetBinContent(1, nVtxTrg);
  _hNumVsCut -> SetBinContent(2, nTwrTrg);
  _hNumVsCut -> SetBinContent(3, nPproj);
  _hNumVsCut -> SetBinContent(4, nEtaTrg);
  _hNumVsCut -> SetBinContent(5, nMatchTrg);
  _hNumVsCut -> SetBinContent(6, nEtTrg);
  _hNumVsCut -> SetBinContent(7, nTspTrg);
  _hNumVsCut -> GetXaxis() -> SetBinLabel(1, "vtx. cuts");
  _hNumVsCut -> GetXaxis() -> SetBinLabel(2, "bad twr.");
  _hNumVsCut -> GetXaxis() -> SetBinLabel(3, "p_{proj} cut");
  _hNumVsCut -> GetXaxis() -> SetBinLabel(4, "#eta cut");
  _hNumVsCut -> GetXaxis() -> SetBinLabel(5, "matching");
  _hNumVsCut -> GetXaxis() -> SetBinLabel(6, "E_{T}^{trg} cut");
  _hNumVsCut -> GetXaxis() -> SetBinLabel(7, "TSP cut");

  // set nTrg ratios vs. cut values
  const Double_t rTwrTrg   = (Double_t) nTwrTrg / (Double_t) nVtxTrg;
  const Double_t rPproj    = (Double_t) nPproj / (Double_t) nVtxTrg;
  const Double_t rEtaTrg   = (Double_t) nEtaTrg / (Double_t) nVtxTrg;
  const Double_t rMatchTrg = (Double_t) nMatchTrg / (Double_t) nVtxTrg;
  const Double_t rEtTrg    = (Double_t) nEtTrg / (Double_t) nVtxTrg;
  const Double_t rTspTrg   = (Double_t) nTspTrg / (Double_t) nVtxTrg;
  _hRatioVsCut -> SetBinContent(1, rTwrTrg);
  _hRatioVsCut -> SetBinContent(2, rPproj);
  _hRatioVsCut -> SetBinContent(3, rEtaTrg);
  _hRatioVsCut -> SetBinContent(4, rMatchTrg);
  _hRatioVsCut -> SetBinContent(5, rEtTrg);
  _hRatioVsCut -> SetBinContent(6, rTspTrg);
  _hRatioVsCut -> GetXaxis() -> SetBinLabel(1, "bad twr.");
  _hRatioVsCut -> GetXaxis() -> SetBinLabel(2, "p_{proj} cut");
  _hRatioVsCut -> GetXaxis() -> SetBinLabel(3, "#eta cut");
  _hRatioVsCut -> GetXaxis() -> SetBinLabel(4, "matching");
  _hRatioVsCut -> GetXaxis() -> SetBinLabel(5, "E_{T}^{trg} cut");
  _hRatioVsCut -> GetXaxis() -> SetBinLabel(6, "TSP cut");

  // particle event loop
  Int_t  iPrevEvt(-1);
  UInt_t pBytes(0);
  UInt_t nParBytes(0);
  UInt_t nTrgPar(0);
  for (UInt_t iEvt = 0; iEvt < nParEvts; iEvt++) {

    // load entry
    pBytes     = _tPar -> GetEntry(iEvt);
    nParBytes += pBytes;
    if (pBytes < 0) {
      cerr << "WARNING: issue with entry " << iEvt << "!\n"
           << "  no. of par. bytes = " << pBytes << "\n"
           << endl;
      break;
    } else {
      if (_isInBatchMode) {
        cout << "      Processing particle event " << iEvt + 1 << "/" << nParEvts << "..." << endl;
      } else {
        cout << "      Processing particle event " << iEvt + 1 << "/" << nParEvts << "...\r" << flush;
        if ((iEvt + 1) == nParEvts) cout << endl;
      }
    }

    // particle event info
    const Int_t    runIdPar = _pMcRunId;
    const Int_t    evtIdPar = _pMcEventId;
    const UInt_t   nTrkPar  = _pMcNumTrks;
    const Double_t zVtxPar  = _pMcVz;
    const Double_t rVtxPar  = TMath::Sqrt((_pMcVx * _pMcVx) + (_pMcVy * _pMcVy));

    // don't double count mc triggers
    if (_avoidDoubleCounting) {
      if (iPrevEvt == evtIdPar) {
        continue;
      } else {
        iPrevEvt = evtIdPar;
      }
    }

    // filter out bad runs
    Bool_t isGoodRunPar(true);
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (runIdPar == _badRunList[iRun]) isGoodRunPar = false;
    }

    // particle event cuts
    const Bool_t isInZcutPar    = (TMath::Abs(zVtxPar) < _zVtxMax);
    const Bool_t isInRcutPar    = (TMath::Abs(rVtxPar) < _rVtxMax);
    const Bool_t isGoodEventPar = (isGoodRunPar && isInRcutPar && isInZcutPar);
    if (_useParTrigger && !isGoodEventPar) continue;

    // particle track loop
    for (UInt_t iTrkMC = 0; iTrkMC < nTrkPar; iTrkMC++) {

      // particle trigger info
      const Int_t    pidTrkMc = _pMcIdGeant -> at(iTrkMC);
      const Double_t cTrkMc   = _pMcCharge  -> at(iTrkMC);
      const Double_t fTrkMc   = _pMcPhi     -> at(iTrkMC);
      const Double_t hTrkMc   = _pMcEta     -> at(iTrkMC);
      const Double_t pTtrkMc  = _pMcPt      -> at(iTrkMC);


      // check pid
      Bool_t isHadronPar(false);
      for (UInt_t iParID = 0; iParID < NTrgIds; iParID++) {
        Bool_t isSameId = (pidTrkMc == _idTrgUse[iParID]);
        if (isSameId) {
          isHadronPar = true;
          break;
        }
      }


      // particle trigger cuts
      const Bool_t isNeutralTrkPar  = (cTrkMc == _cTrgPar);
      const Bool_t isInTrkEtaCutPar = (TMath::Abs(hTrkMc) < _hTrgParMax);
      const Bool_t isInTrkEtCutPar  = ((pTtrkMc > _eTtrgParMin) && (pTtrkMc < _eTtrgParMax));
      if (isHadronPar && isNeutralTrkPar && isInTrkEtaCutPar && isInTrkEtCutPar) {

	const Double_t wTrkMc   = _fWeight     -> Eval(pTtrkMc);

	// check for neighbors first
	int index_daughter[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
	int numDaughters = 0;

	hasNeighbor[iTrkMC] = false;
        // check for neighbors
        for (UInt_t iCheck = 0; iCheck < nTrkPar; iCheck++) {

          // check if same particle
          const Bool_t isSameParticle = (iTrkMC == iCheck);
          if (isSameParticle) continue;
	  
	  // check for daughter particles (i.e. skip if from same parent)
	  if (_pMcIdVx -> at(iCheck) == _pMcIdVxEnd -> at(iTrkMC))
	    {
	      //this is a direct daughter of MC particle checking for neighbors
	      index_daughter[numDaughters] = iCheck;
	      numDaughters++;
	      continue;
	    }
	  int daughterOfDaughter = 0;
	  if (numDaughters>0) // at least one daughter particle was found
	    {
	      for (int iDaughter=0; iDaughter<numDaughters; iDaughter++)
		{
		  if (_pMcIdVx -> at(iCheck) == _pMcIdVxEnd -> at(index_daughter[iDaughter]))
		// daughter of daughter , e.g. pi0 -> gamma gamma -> gamma e+ e-
		    daughterOfDaughter = 1;
		}
	      if (daughterOfDaughter > 0) continue;
	    }

          // possible neighbor info
	  const Bool_t CheckIsShower =   _pMcIsShower      -> at(iCheck); // allow for geant particle from interactions with material
	  if (CheckIsShower) continue;  //do not count shower particles as neighbors

          const Double_t pidCheck = _pMcIdGeant -> at(iCheck);
          const Double_t cCheck   = _pMcCharge  -> at(iCheck);
          const Double_t fCheck   = _pMcPhi     -> at(iCheck);
          Double_t hCheck   = _pMcEta     -> at(iCheck);
          Double_t Echeck  = _pMcEnergy      -> at(iCheck);
	  int check_daughter = 0;
	  if (hCheck == -999) { hCheck = TMath::ATanH(_pMcPz -> at(iCheck) / _pMcPtot -> at(iCheck)); /*check_daughter = 1;  cout << "recalculate eta " << hCheck<< endl;*/}

          // calculate delta-r
          Double_t dFcheck = fCheck - fTrkMc;
          if (dFcheck > TMath::Pi()) {
            dFcheck -= TMath::TwoPi();
          }
          if (dFcheck < (-1. * TMath::Pi())) {
            dFcheck += TMath::TwoPi();
          }
          Double_t dHcheck = hCheck - hTrkMc;
          Double_t dRcheck = TMath::Sqrt((dFcheck * dFcheck) + (dHcheck * dHcheck));
	  
	  // const Bool_t isInCheckPtCut    = (pTcheck > 0.);

          const Bool_t isInMinSeparation = (dRcheck < _dMinSeparationUse);
	  //          if (isInCheckPtCut && isInMinSeparation) {
          if (isInMinSeparation) {
            hasNeighbor[iTrkMC] = true;
	    //  cout << "pid (with vxend = " << pidTrgMc << "(" << _mMcIdVxEnd -> at(iTrgMC) << ") has neighbor within " << dRcheck << ", which is a pi = " << pidCheck << "with vx = " <<_mMcIdVx -> at(iCheck)<< endl;
	    //            hasNeighbor[iCheck] = true;
            break;
          }
        }  // end neighbor loop
	// don't consider if particle has neighbor
	if (_excludeNeighbors && hasNeighbor[iTrkMC]) continue;

        _hPhiTrg[8]    -> Fill(fTrkMc, wTrkMc);
        _hPhiTrg[9]    -> Fill(fTrkMc);
        _hEtaTrg[8]    -> Fill(hTrkMc, wTrkMc);
        _hEtaTrg[9]    -> Fill(hTrkMc);
        _hEtTrg[8]     -> Fill(pTtrkMc, wTrkMc);
        _hEtTrg[9]     -> Fill(pTtrkMc);
        _hPhiForEff[1] -> Fill(fTrkMc, wTrkMc);
        _hEtaForEff[1] -> Fill(hTrkMc, wTrkMc);
        _hEtForEff[1]  -> Fill(pTtrkMc, wTrkMc);
        ++nTrgPar;
      }
    }  // end particle track llop
  }  // end particle event loop
  cout << "    Particle event loop finished:\n"
       << "      nTrgPar = " << nTrgPar
       << endl;

  // calculate efficiency
  const Float_t dFbin     = _hPhiForEff[0] -> GetBinWidth(1);
  const Float_t dHbin     = _hEtaForEff[0] -> GetBinWidth(1);
  const Float_t dEtBin    = _hEtForEff[0]  -> GetBinWidth(1);
  const Float_t fNormDet  = dFbin * nTrgMat;
  const Float_t fNormPar  = dFbin * nTrgPar;
  const Float_t hNormDet  = dHbin * nTrgMat;
  const Float_t hNormPar  = dHbin * nTrgPar;
  const Float_t eTnormDet = dEtBin * nTrgMat;
  const Float_t eTnormPar = dEtBin * nTrgPar;
  _hPhiForEff[0] -> Scale(1. / fNormDet);
  _hPhiForEff[1] -> Scale(1. / fNormPar);
  _hEtaForEff[0] -> Scale(1. / hNormDet);
  _hEtaForEff[1] -> Scale(1. / hNormPar);
  _hEtForEff[0]  -> Scale(1. / eTnormDet);
  _hEtForEff[1]  -> Scale(1. / eTnormPar);

  const Float_t parWeight(1.);
  const Float_t detWeight(1.);
  _hPhiEff -> Divide(_hPhiForEff[0], _hPhiForEff[1], detWeight, parWeight);
  _hEtaEff -> Divide(_hEtaForEff[0], _hEtaForEff[1], detWeight, parWeight);
  _hEtEff  -> Divide(_hEtForEff[0], _hEtForEff[1], detWeight, parWeight);
  cout << "    Calculated matching efficiencies." << endl;

}  // end 'DoDetFirstMatching()'



void StTriggerEnergyScaleCalculator::DoParFirstMatching() {

  cout << "    Particle first matching coming soon!" << endl;

}  // end 'DoParFirstMatching()'

// End ------------------------------------------------------------------------
