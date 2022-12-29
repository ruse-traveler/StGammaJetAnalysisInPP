// 'StTrackEfficiencyCalculator.cxx'
// Derek Anderson
// 01.30.2020
//
// This class produces trees of jets using
// generated particles and reconstructed
// tracks matched to generated particles.
// Also calculates tracking efficiency.


#define StTrackEfficiencyCalculator_cxx

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StTrackEfficiencyCalculator.h"

using namespace std;
using namespace fastjet;

ClassImp(StTrackEfficiencyCalculator);



void StTrackEfficiencyCalculator::Init(const TString sInput, const TString sOutput, const TString sParTree, const TString sDetTree) {

  // announce start
  cout << "\n  Beginning efficiency calculation..." << endl;

  // set file and tree names
  _sIn  = sInput;
  _sOut = sOutput;
  _sPar = sParTree;
  _sDet = sDetTree;

  // open files
  _fIn  = new TFile(_sIn.Data(), "read");
  _fOut = new TFile(_sOut.Data(), "recreate");
  if (!_fIn || !_fOut) {
    cerr << "PANIC: couldn't open input file and/or create output file!\n"
         << "       fIn = " << _fIn << ", fOut = " << _fOut
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input trees
  _fIn -> GetObject(_sPar.Data(), _tPar);
  _fIn -> GetObject(_sDet.Data(), _tDet);
  if (!_tPar || !_tDet) {
    cerr << "PANIC: couldn't grab particle and/or detector level tree!\n"
         << "       tPar = " << _tPar << ", tDet = " << _tDet
         << endl;
    return;
  }
  cout << "    Grabbed input trees." << endl;

  // initialize output trees
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
  _tJetPar -> Branch("TrgEt", &_pTrgEt, "TrgEt/D");
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
  _tJetDet -> Branch("TrgEt", &_dTrgEt, "TrgEt/D");
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

}  // end 'Init(TString, TString, TString, TString)'



void StTrackEfficiencyCalculator::Make(const UInt_t trigger, const Bool_t parTrg, const Bool_t detTrg, const Bool_t recoilTrks, const Bool_t trkSpecies) {

  // set global flags
  _fTrig                   = trigger;
  _useParTrg               = parTrg;
  _useDetTrg               = detTrg;
  _useOnlyRecoilTracks     = recoilTrks;
  _useSpecificTrackSpecies = trkSpecies;


  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);
  const UInt_t   badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};

  // trigger parameters
  const UInt_t   idChrg[NChrgIds] = {8, 9, 11, 12, 14, 15, 19, 21, 23, 24, 27, 29, 31, 32, 45, 46};
  const UInt_t   idPion[NPionIds] = {7};
  const Double_t hTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  const Double_t nSigTrgCut(3.);

  // track parameters
  const Int_t    minQaTruth(50);
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(100.);
  const Double_t dFminAway(TMath::PiOver2());
  const Double_t dFmaxAway(3. * TMath::PiOver2());
  const Double_t idTrkSpecies[NTrkSpecies] = {8, 9};


  // set particle branch addresses
  _tPar -> SetMakeClass(1);
  _tPar -> SetBranchAddress("EventId", &mcEventId, &_bMcEventId);
  _tPar -> SetBranchAddress("RunId", &mcRunId, &_bMcRunId);
  _tPar -> SetBranchAddress("NumTrks", &mcNumTrks, &_bMcNumTrks);
  _tPar -> SetBranchAddress("MuVtxX", &muVx, &_bMuVx);
  _tPar -> SetBranchAddress("MuVtxY", &muVy, &_bMuVy);
  _tPar -> SetBranchAddress("MuVtxZ", &muVz, &_bMuVz);
  _tPar -> SetBranchAddress("McVtxX", &mcVx, &_bMcVx);
  _tPar -> SetBranchAddress("McVtxY", &mcVy, &_bMcVy);
  _tPar -> SetBranchAddress("McVtxZ", &mcVz, &_bMcVz);
  _tPar -> SetBranchAddress("IdTrk", &mcIdTrk, &_bMcIdTrk);
  _tPar -> SetBranchAddress("IdGeant", &mcIdGeant, &_bMcIdGeant);
  _tPar -> SetBranchAddress("IdVx", &mcIdVx, &_bMcIdVx);
  _tPar -> SetBranchAddress("IdVxEnd", &mcIdVxEnd, &_bMcIdVxEnd);
  _tPar -> SetBranchAddress("IntrVtx", &mcIntrVtx, &_bMcIntrVtx);
  _tPar -> SetBranchAddress("IsShower", &mcIsShower, &_bMcIsShower);
  _tPar -> SetBranchAddress("Charge", &mcCharge, &_bMcCharge);
  _tPar -> SetBranchAddress("Rapidity", &mcRapidity, &_bMcRapidity);
  _tPar -> SetBranchAddress("Eta", &mcEta, &_bMcEta);
  _tPar -> SetBranchAddress("Phi", &mcPhi, &_bMcPhi);
  _tPar -> SetBranchAddress("Px", &mcPx, &_bMcPx);
  _tPar -> SetBranchAddress("Py", &mcPy, &_bMcPy);
  _tPar -> SetBranchAddress("Pz", &mcPz, &_bMcPz);
  _tPar -> SetBranchAddress("Pt", &mcPt, &_bMcPt);
  _tPar -> SetBranchAddress("Ptot", &mcPtot, &_bMcPtot);
  _tPar -> SetBranchAddress("Energy", &mcEnergy, &_bMcEnergy);

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
  _tDet -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &_bEventList_PrimaryTrackArray_);
  _tDet -> SetBranchAddress("PrimaryTrackArray.fUniqueID", &PrimaryTrackArray_fUniqueID, &_bPrimaryTrackArray_fUniqueID);
  _tDet -> SetBranchAddress("PrimaryTrackArray.fBits", &PrimaryTrackArray_fBits, &_bPrimaryTrackArray_fBits);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nHitsFit", &PrimaryTrackArray_nHitsFit, &_bPrimaryTrackArray_nHitsFit);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", &PrimaryTrackArray_nHitsPoss, &_bPrimaryTrackArray_nHitsPoss);
  _tDet -> SetBranchAddress("PrimaryTrackArray.trackFlag", &PrimaryTrackArray_trackFlag, &_bPrimaryTrackArray_trackFlag);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pdgId", &PrimaryTrackArray_pdgId, &_bPrimaryTrackArray_pdgId);
  _tDet -> SetBranchAddress("PrimaryTrackArray.geantId", &PrimaryTrackArray_geantId, &_bPrimaryTrackArray_geantId);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pZ", &PrimaryTrackArray_pZ, &_bPrimaryTrackArray_pZ);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pX", &PrimaryTrackArray_pX, &_bPrimaryTrackArray_pX);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pY", &PrimaryTrackArray_pY, &_bPrimaryTrackArray_pY);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pT", &PrimaryTrackArray_pT, &_bPrimaryTrackArray_pT);
  _tDet -> SetBranchAddress("PrimaryTrackArray.dEdx", &PrimaryTrackArray_dEdx, &_bPrimaryTrackArray_dEdx);
  _tDet -> SetBranchAddress("PrimaryTrackArray.charge", &PrimaryTrackArray_charge, &_bPrimaryTrackArray_charge);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofBeta", &PrimaryTrackArray_tofBeta, &_bPrimaryTrackArray_tofBeta);
  _tDet -> SetBranchAddress("PrimaryTrackArray.eta", &PrimaryTrackArray_eta, &_bPrimaryTrackArray_eta);
  _tDet -> SetBranchAddress("PrimaryTrackArray.phi", &PrimaryTrackArray_phi, &_bPrimaryTrackArray_phi);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigElectron", &PrimaryTrackArray_nSigElectron, &_bPrimaryTrackArray_nSigElectron);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigPion", &PrimaryTrackArray_nSigPion, &_bPrimaryTrackArray_nSigPion);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigKaon", &PrimaryTrackArray_nSigKaon, &_bPrimaryTrackArray_nSigKaon);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nSigProton", &PrimaryTrackArray_nSigProton, &_bPrimaryTrackArray_nSigProton);
  _tDet -> SetBranchAddress("PrimaryTrackArray.dcag", &PrimaryTrackArray_dcag, &_bPrimaryTrackArray_dcag);
  _tDet -> SetBranchAddress("PrimaryTrackArray.nHits", &PrimaryTrackArray_nHits, &_bPrimaryTrackArray_nHits);
  _tDet -> SetBranchAddress("PrimaryTrackArray.dEdxHits", &PrimaryTrackArray_dEdxHits, &_bPrimaryTrackArray_dEdxHits);
  _tDet -> SetBranchAddress("PrimaryTrackArray.firstZPoint", &PrimaryTrackArray_firstZPoint, &_bPrimaryTrackArray_firstZPoint);
  _tDet -> SetBranchAddress("PrimaryTrackArray.lastZPoint", &PrimaryTrackArray_lastZPoint, &_bPrimaryTrackArray_lastZPoint);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", &PrimaryTrackArray_tofSigElectron, &_bPrimaryTrackArray_tofSigElectron);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigPion", &PrimaryTrackArray_tofSigPion, &_bPrimaryTrackArray_tofSigPion);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", &PrimaryTrackArray_tofSigKaon, &_bPrimaryTrackArray_tofSigKaon);
  _tDet -> SetBranchAddress("PrimaryTrackArray.tofSigProton", &PrimaryTrackArray_tofSigProton, &_bPrimaryTrackArray_tofSigProton);
  _tDet -> SetBranchAddress("PrimaryTrackArray.timeOfflight", &PrimaryTrackArray_timeOfflight, &_bPrimaryTrackArray_timeOfflight);
  _tDet -> SetBranchAddress("PrimaryTrackArray.pathLength", &PrimaryTrackArray_pathLength, &_bPrimaryTrackArray_pathLength);
  _tDet -> SetBranchAddress("PrimaryTrackArray.trkIndex", &PrimaryTrackArray_trkIndex, &_bPrimaryTrackArray_trkIndex);
  _tDet -> SetBranchAddress("TowerArray", &TowerArray_, &_bEventList_TowerArray_);
  _tDet -> SetBranchAddress("TowerArray.fUniqueID", &TowerArray_fUniqueID, &_bTowerArray_fUniqueID);
  _tDet -> SetBranchAddress("TowerArray.fBits", &TowerArray_fBits, &_bTowerArray_fBits);
  _tDet -> SetBranchAddress("TowerArray.TwrId", &TowerArray_TwrId, &_bTowerArray_TwrId);
  _tDet -> SetBranchAddress("TowerArray.TwrEng", &TowerArray_TwrEng, &_bTowerArray_TwrEng);
  _tDet -> SetBranchAddress("TowerArray.TwrEta", &TowerArray_TwrEta, &_bTowerArray_TwrEta);
  _tDet -> SetBranchAddress("TowerArray.TwrPhi", &TowerArray_TwrPhi, &_bTowerArray_TwrPhi);
  _tDet -> SetBranchAddress("TowerArray.TwrADC", &TowerArray_TwrADC, &_bTowerArray_TwrADC);
  _tDet -> SetBranchAddress("TowerArray.TwrPed", &TowerArray_TwrPed, &_bTowerArray_TwrPed);
  _tDet -> SetBranchAddress("TowerArray.TwrRMS", &TowerArray_TwrRMS, &_bTowerArray_TwrRMS);
  _tDet -> SetBranchAddress("TowerArray.TwrMatchIdnex", &TowerArray_TwrMatchIdnex, &_bTowerArray_TwrMatchIdnex);
  _tDet -> SetBranchAddress("TowerArray.NoOfmatchedTrk", &TowerArray_NoOfmatchedTrk, &_bTowerArray_NoOfmatchedTrk);
  _tDet -> SetBranchAddress("TowerArray.TwrMatchP", &TowerArray_TwrMatchP, &_bTowerArray_TwrMatchP);
  _tDet -> SetBranchAddress("TowerArray.TwrPx", &TowerArray_TwrPx, &_bTowerArray_TwrPx);
  _tDet -> SetBranchAddress("TowerArray.TwrPy", &TowerArray_TwrPy, &_bTowerArray_TwrPy);
  _tDet -> SetBranchAddress("TowerArray.TwrPz", &TowerArray_TwrPz, &_bTowerArray_TwrPz);
  _tDet -> SetBranchAddress("TowerArray.fNAssocTracks", &TowerArray_fNAssocTracks, &_bTowerArray_fNAssocTracks);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", &TowerArray_fMatchedTracksArray_, &_bTowerArray_fMatchedTracksArray_);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", &TowerArray_fMatchedTracksArray_P, &_bTowerArray_fMatchedTracksArray_P);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigPi[10]", &TowerArray_fMatchedTracksArray_nSigPi, &_bTowerArray_fMatchedTracksArray_nSigPi);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigK[10]", &TowerArray_fMatchedTracksArray_nSigK, &_bTowerArray_fMatchedTracksArray_nSigK);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigP[10]", &TowerArray_fMatchedTracksArray_nSigP, &_bTowerArray_fMatchedTracksArray_nSigP);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigE[10]", &TowerArray_fMatchedTracksArray_nSigE, &_bTowerArray_fMatchedTracksArray_nSigE);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_dcag[10]", &TowerArray_fMatchedTracksArray_dcag, &_bTowerArray_fMatchedTracksArray_dcag);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_eta[10]", &TowerArray_fMatchedTracksArray_eta, &_bTowerArray_fMatchedTracksArray_eta);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_pT[10]", &TowerArray_fMatchedTracksArray_pT, &_bTowerArray_fMatchedTracksArray_pT);
  _tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nFit[10]", &TowerArray_fMatchedTracksArray_nFit, &_bTowerArray_fMatchedTracksArray_nFit); 
  cout << "    Set input branch addresses." << endl;


  // define particle/detector histograms
  TH1D *hPhiTrk[NLvls][NCuts];
  TH1D *hEtaTrk[NLvls][NCuts];
  TH1D *hPtTrk[NLvls][NCuts];
  // define matching histograms
  TH1D *hPhiForEff[NLvls];
  TH1D *hEtaForEff[NLvls];
  TH1D *hPtForEff[NLvls];
  TH2D *hPhiParVsDet;
  TH2D *hEtaParVsDet;
  TH2D *hPtParVsDet;
  // define efficiency histograms
  TH1D *hPhiEff;
  TH1D *hEtaEff;
  TH1D *hPtEff;
  // define resolution histograms
  TH1D *hPhiRes;
  TH1D *hEtaRes;
  TH1D *hPtRes;
  TH2D *hPhiResVsPhi;
  TH2D *hEtaResVsEta;
  TH2D *hPtResVsPt;
  // define QA histograms
  TH1D *hEtaQA[NHistQA];
  TH1D *hPtQA[NHistQA];
  TH1D *hNfitQA[NHistQA - 1];
  TH1D *hRfitQA[NHistQA - 1];
  TH1D *hDcaQA[NHistQA - 1];
  TH1D *hPtFrac[NHistQA - 1];
  TH1D *hIdTruth[NHistQA - 1];
  TH1D *hQaTruth[NHistQA - 1];
  TH1D *hNcommon[NHistQA - 1];
  // define miscellaneous QA histograms
  TH1D *hNumTrk[NLvls];
  TH1D *hNumJetTrk[NLvls + 2];
  TH1D *hTrkVtxId[2];
  TH1D *hNumMatchVtx[NLvls];
  TH1D *hMatchVtxId[NLvls];
  TH1D *hPtMatchVtx;

  // pT binning
  const UInt_t   nPtHistBins(NPtBins - 1);
  const Double_t pTbinEdges[NPtBins] = {0., 0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 10., 12., 14., 16., 18., 20.};

  // binning
  const UInt_t   nPhi(60);
  const UInt_t   nEta(40);
  const UInt_t   nDf(240);
  const UInt_t   nDh(400);
  const UInt_t   nDpt(1000);
  const UInt_t   nPt(1000);
  const UInt_t   nNum(200);
  const UInt_t   nNfit(50);
  const UInt_t   nRfit(100);
  const UInt_t   nDca(50);
  const UInt_t   nPtFrac(100);
  const UInt_t   nIdTruth(500);
  const UInt_t   nQaTruth(100);
  const UInt_t   nNcommon(50);
  const UInt_t   nNumTrkQA(100);
  const UInt_t   nNumVtxQA(500);
  const Double_t f[2]       = {-3.15, 3.15};
  const Double_t h[2]       = {-1., 1.};
  const Double_t dF[2]      = {-12.56, 12.56};
  const Double_t dH[2]      = {-10., 10.};
  const Double_t dPt[2]     = {-100., 100.};
  const Double_t pT[2]      = {0., 100.};
  const Double_t num[2]     = {0., 200.};
  const Double_t nFit[2]    = {0., 50.};
  const Double_t rFit[2]    = {0., 1.};
  const Double_t dca[2]     = {0., 5.};
  const Double_t hQA[2]     = {-2., 2.};
  const Double_t pTfrac[2]  = {0., 10.};
  const Double_t idTruth[2] = {0., 500.};
  const Double_t qaTruth[2] = {0., 100.};
  const Double_t nCommon[2] = {0., 50.};
  const Double_t nTrkQA[2]  = {0., 100.};
  const Double_t nVtxQA[2]  = {0., 500.};
  // create particle histograms
  hPhiTrk[0][0]   = new TH1D("hPhiBeforeQA_par", "", nPhi, f[0], f[1]);
  hPhiTrk[0][1]   = new TH1D("hPhiAfterQA_par", "", nPhi, f[0], f[1]);
  hEtaTrk[0][0]   = new TH1D("hEtaBeforeQA_par", "", nEta, h[0], h[1]);
  hEtaTrk[0][1]   = new TH1D("hEtaAfterQA_par", "", nEta, h[0], h[1]);
  hPtTrk[0][0]    = new TH1D("hPtBeforeQA_par", "", nPtHistBins, pTbinEdges);
  hPtTrk[0][1]    = new TH1D("hPtAfterQA_par", "", nPtHistBins, pTbinEdges);
  // create detector histograms
  hPhiTrk[1][0]   = new TH1D("hPhiBeforeQA_det", "", nPhi, f[0], f[1]);
  hPhiTrk[1][1]   = new TH1D("hPhiAfterQA_det", "", nPhi, f[0], f[1]);
  hEtaTrk[1][0]   = new TH1D("hEtaBeforeQA_det", "", nEta, h[0], h[1]);
  hEtaTrk[1][1]   = new TH1D("hEtaAfterQA_det", "", nEta, h[0], h[1]);
  hPtTrk[1][0]    = new TH1D("hPtBeforeQA_det", "", nPtHistBins, pTbinEdges);
  hPtTrk[1][1]    = new TH1D("hPtAfterQA_det", "", nPtHistBins, pTbinEdges);
  // create matching histograms
  hPhiForEff[0]   = new TH1D("hPhiForEff_par", "", nPhi, f[0], f[1]);
  hPhiForEff[1]   = new TH1D("hPhiForEff_det", "", nPhi, f[0], f[1]);
  hEtaForEff[0]   = new TH1D("hEtaForEff_par", "", nEta, h[0], h[1]);
  hEtaForEff[1]   = new TH1D("hEtaForEff_det", "", nEta, h[0], h[1]);
  hPtForEff[0]    = new TH1D("hPtForEff_par", "", nPtHistBins, pTbinEdges);
  hPtForEff[1]    = new TH1D("hPtForEff_det", "", nPtHistBins, pTbinEdges);
  hPhiParVsDet    = new TH2D("hPhiParVsDet", "", nPhi, f[0], f[1], nPhi, f[0], f[1]);
  hEtaParVsDet    = new TH2D("hEtaParVsDet", "", nEta, h[0], h[1], nEta, h[0], h[1]);
  hPtParVsDet     = new TH2D("hPtParVsDet", "", nPt, pT[0], pT[1], nPt, pT[0], pT[1]);
  // create efficiency and resolution histograms
  hPhiEff         = new TH1D("hPhiEfficiency", "", nPhi, f[0], f[1]);
  hEtaEff         = new TH1D("hEtaEfficiency", "", nEta, h[0], h[1]);
  hPtEff          = new TH1D("hPtEfficiency", "", nPtHistBins, pTbinEdges);
  hPhiRes         = new TH1D("hPhiRes", "", nDf, dF[0], dF[1]);
  hEtaRes         = new TH1D("hEtaRes", "", nDh, dH[0], dH[1]);
  hPtRes          = new TH1D("hPtRes", "", nDpt, dPt[0], dPt[1]);
  hPhiResVsPhi    = new TH2D("hPhiResVsPhi", "", nPhi, f[0], f[1], nDf, dF[0], dF[1]);
  hEtaResVsEta    = new TH2D("hEtaResVsEta", "", nEta, h[0], h[1], nDh, dH[0], dH[1]);
  hPtResVsPt      = new TH2D("hPtResVsPt", "", nPt, pT[0], pT[1], nDpt, dPt[0], dPt[1]);
  // create QA histograms
  hEtaQA[0]       = new TH1D("hEtaParQA", "", nEta, hQA[0], hQA[1]);
  hEtaQA[1]       = new TH1D("hEtaDetQA", "", nEta, hQA[0], hQA[1]);
  hEtaQA[2]       = new TH1D("hEtaMatQA", "", nEta, hQA[0], hQA[1]);
  hEtaQA[3]       = new TH1D("hEtaUnMatQA", "", nEta, hQA[0], hQA[1]);
  hPtQA[0]        = new TH1D("hPtParQA", "", nPtHistBins, pTbinEdges);
  hPtQA[1]        = new TH1D("hPtDetQA", "", nPtHistBins, pTbinEdges);
  hPtQA[2]        = new TH1D("hPtMatQA", "", nPtHistBins, pTbinEdges);
  hPtQA[3]        = new TH1D("hPtUnMatQA", "", nPtHistBins, pTbinEdges);
  hNfitQA[0]      = new TH1D("hNfitDetQA", "", nNfit, nFit[0], nFit[1]);
  hNfitQA[1]      = new TH1D("hNfitMatQA", "", nNfit, nFit[0], nFit[1]);
  hNfitQA[2]      = new TH1D("hNfitUnMatQA", "", nNfit, nFit[0], nFit[1]);
  hRfitQA[0]      = new TH1D("hRfitDetQA", "", nRfit, rFit[0], rFit[1]);
  hRfitQA[1]      = new TH1D("hRfitMatQA", "", nRfit, rFit[0], rFit[1]);
  hRfitQA[2]      = new TH1D("hRfitUnMatQA", "", nRfit, rFit[0], rFit[1]);
  hDcaQA[0]       = new TH1D("hDcaDetQA", "", nDca, dca[0], dca[1]);
  hDcaQA[1]       = new TH1D("hDcaMatQA", "", nDca, dca[0], dca[1]);
  hDcaQA[2]       = new TH1D("hDcaUnMatQA", "", nDca, dca[0], dca[1]);
  hPtFrac[0]      = new TH1D("hPtFracDetQA", "", nPtFrac, pTfrac[0], pTfrac[1]);
  hPtFrac[1]      = new TH1D("hPtFracMatQA", "", nPtFrac, pTfrac[0], pTfrac[1]);
  hPtFrac[2]      = new TH1D("hPtFracUnMatQA", "", nPtFrac, pTfrac[0], pTfrac[1]);
  hIdTruth[0]     = new TH1D("hIdTruthDetQA", "", nIdTruth, idTruth[0], idTruth[1]);
  hIdTruth[1]     = new TH1D("hIdTruthMatQA", "", nIdTruth, idTruth[0], idTruth[1]);
  hIdTruth[2]     = new TH1D("hIdTruthUnMatQA", "", nIdTruth, idTruth[0], idTruth[1]);
  hQaTruth[0]     = new TH1D("hQaTruthDetQA", "", nQaTruth, qaTruth[0], qaTruth[1]);
  hQaTruth[1]     = new TH1D("hQaTruthMatQA", "", nQaTruth, qaTruth[0], qaTruth[1]);
  hQaTruth[2]     = new TH1D("hQaTruthUnMatQA", "", nQaTruth, qaTruth[0], qaTruth[1]);
  hNcommon[0]     = new TH1D("hNcommonDetQA", "", nNcommon, nCommon[0], nCommon[1]);
  hNcommon[1]     = new TH1D("hNcommonMatQA", "", nNcommon, nCommon[0], nCommon[1]);
  hNcommon[2]     = new TH1D("hNcommonUnMatQA", "", nNcommon, nCommon[0], nCommon[1]);
  // create miscellaneous QA histograms
  hNumTrk[0]      = new TH1D("hNumTrk_par", "After QA", nNum, num[0], num[1]);
  hNumTrk[1]      = new TH1D("hNumTrk_det", "After QA", nNum, num[0], num[1]);
  hNumJetTrk[0]   = new TH1D("hNumJetTrk_par", "", nNumTrkQA, nTrkQA[0], nTrkQA[1]);
  hNumJetTrk[1]   = new TH1D("hNumJetTrk_det", "", nNumTrkQA, nTrkQA[0], nTrkQA[1]);
  hNumJetTrk[2]   = new TH1D("hNumJetTrk_match", "", nNumTrkQA, nTrkQA[0], nTrkQA[1]);
  hNumJetTrk[3]   = new TH1D("hNumJetTrk_notMatch", "", nNumTrkQA, nTrkQA[0], nTrkQA[1]);
  hNumMatchVtx[0] = new TH1D("hNumMatchVtx_par", "", nNumVtxQA, nVtxQA[0], nVtxQA[1]);
  hNumMatchVtx[1] = new TH1D("hNumMatchVtx_det", "", nNumVtxQA, nVtxQA[0], nVtxQA[1]);
  hMatchVtxId[0]  = new TH1D("hMatchVtxId_par", "", nNumVtxQA, nVtxQA[0], nVtxQA[1]);
  hMatchVtxId[1]  = new TH1D("hMatchVtxId_det", "", nNumVtxQA, nVtxQA[0], nVtxQA[1]);
  hTrkVtxId[0]    = new TH1D("hTrkVtxId_match", "", nNumVtxQA, nVtxQA[0], nVtxQA[1]);
  hTrkVtxId[1]    = new TH1D("hTrkVtxId_notMatch", "", nNumVtxQA, nVtxQA[0], nVtxQA[1]);
  hPtMatchVtx     = new TH1D("hPtAfterQA_detOnlyMatchVtxs", "", nPtHistBins, pTbinEdges);
  // errors
  for (UInt_t iLvl = 0; iLvl < NLvls; iLvl++) {
    for (UInt_t iCut = 0; iCut < NCuts; iCut++) {
      hPhiTrk[iLvl][iCut] -> Sumw2();
      hEtaTrk[iLvl][iCut] -> Sumw2();
      hPtTrk[iLvl][iCut]  -> Sumw2();
    }
    hPhiForEff[iLvl] -> Sumw2();
    hEtaForEff[iLvl] -> Sumw2();
    hPtForEff[iLvl]  -> Sumw2();
    hNumTrk[iLvl]    -> Sumw2();
  }
  for (UInt_t iHistQA = 0; iHistQA < NHistQA; iHistQA++) {
    hEtaQA[iHistQA] -> Sumw2();
    hPtQA[iHistQA]  -> Sumw2();
    if (iHistQA < (NHistQA - 1)) {
      hNfitQA[iHistQA]  -> Sumw2();
      hRfitQA[iHistQA]  -> Sumw2();
      hDcaQA[iHistQA]   -> Sumw2();
      hPtFrac[iHistQA]  -> Sumw2();
      hIdTruth[iHistQA] -> Sumw2();
      hQaTruth[iHistQA] -> Sumw2();
      hNcommon[iHistQA] -> Sumw2();
    }
  }
  hPhiParVsDet    -> Sumw2();
  hEtaParVsDet    -> Sumw2();
  hPtParVsDet     -> Sumw2();
  hPhiEff         -> Sumw2();
  hEtaEff         -> Sumw2();
  hPtEff          -> Sumw2();
  hPhiRes         -> Sumw2();
  hEtaRes         -> Sumw2();
  hPtRes          -> Sumw2();
  hPhiResVsPhi    -> Sumw2();
  hEtaResVsEta    -> Sumw2();
  hPtResVsPt      -> Sumw2();
  hNumJetTrk[0]   -> Sumw2();
  hNumJetTrk[1]   -> Sumw2();
  hNumJetTrk[2]   -> Sumw2();
  hNumJetTrk[3]   -> Sumw2();
  hNumMatchVtx[1] -> Sumw2();
  hNumMatchVtx[0] -> Sumw2();
  hMatchVtxId[1]  -> Sumw2();
  hMatchVtxId[0]  -> Sumw2();
  hTrkVtxId[0]    -> Sumw2();
  hTrkVtxId[1]    -> Sumw2();
  hPtMatchVtx     -> Sumw2();
  cout << "    Created histograms." << endl;

  // make profiles
  TProfile *pPhiResVsPhi = new TProfile("pPhiResVsPhi", "", nPhi, f[0], f[1], "S");
  TProfile *pEtaResVsEta = new TProfile("pEtaResVsEta", "", nEta, h[0], h[1], "S");
  TProfile *pPtResVsPt   = new TProfile("pPtResVsPt", "", nPt, pT[0], pT[1], "S");
  cout << "    Created profiles." << endl;


  // for matching tracks
  Bool_t isMatched[NTrkMax];

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

  // for identifying vtxs w/ matches
  vector<Int_t> matchRecoVtxs;
  vector<Int_t> matchMcVtxs;
  matchRecoVtxs.clear();
  matchMcVtxs.clear();


  // begin event loops
  const UInt_t nParEvts    = _tPar -> GetEntriesFast();
  const UInt_t nDetEvts    = _tDet -> GetEntriesFast();
  const Bool_t areNumsSame = (nParEvts == nDetEvts);
  cout << "    Beginning event loops: " << nParEvts << " generated and " << nDetEvts << " reconstructed events to process." << endl;

  UInt_t nEvts(0);
  if (!areNumsSame) {
    cerr << "PANIC: number of generated and reconstructed events are not the same!" << endl;
    return;
  }
  else {
    nEvts = nParEvts;
  }


  Int_t  parBytes(0);
  Int_t  detBytes(0);
  UInt_t nTrgPar(0);
  UInt_t nTrgDet(0);
  UInt_t nTrgTot(0);
  UInt_t nParBytes(0);
  UInt_t nDetBytes(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // load entries
    parBytes   = _tPar -> GetEntry(iEvt);
    detBytes   = _tDet -> GetEntry(iEvt);
    nParBytes += parBytes;
    nDetBytes += detBytes;
    if (parBytes < 0) {
      cerr << "WARNING: issue with particle entry " << iEvt << "!" << endl;
      break;
    }
    if (detBytes < 0) {
      cerr << "WARNING: issue with detector entry " << iEvt << "!" << endl;
      break;
    }

    // announce entry
    if (_isInBatchMode) {
      cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
    }
    else {
      cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
      if ((iEvt + 1) == nEvts) cout << endl;
    }


    // particle event info
    const Int_t    runIdPar = mcRunId;
    const Int_t    evtIdPar = mcEventId;
    const UInt_t   nTrkPar  = mcNumTrks;
    const Double_t rVtxPar  = TMath::Sqrt((mcVx * mcVx) + (mcVy * mcVy));
    const Double_t zVtxPar  = mcVz;

    // detector event info
    const Int_t    runIdDet = runNumber;
    const Int_t    evtIdDet = eventNumber;
    const UInt_t   nTrkDet  = nPrimaryTracks;
    const Double_t rVtxDet  = TMath::Sqrt((xVertex * xVertex) + (yVertex * yVertex));
    const Double_t zVtxDet  = zVertex;

    // check events and runs
    const Bool_t isSameRun = (runIdPar == runIdDet);
    const Bool_t isSameEvt = (evtIdPar == evtIdDet);
    if (!isSameRun || !isSameEvt) {
      cerr << "WARNING: run or event is NOT the same at iEvt = " << iEvt << "\n"
           << "         run(par, det) = (" << runIdPar << ", " << runIdDet <<")\n"
           << "         evt(par, det) = (" << evtIdPar << ", " << evtIdDet
           << endl;
    }

    // check if good run
    Bool_t isGoodRunPar(true);
    Bool_t isGoodRunDet(true);
    for (UInt_t iBadRun = 0; iBadRun < NBadRuns; iBadRun++) {
      Bool_t isSameRunPar = (runIdPar == badRunList[iBadRun]);
      Bool_t isSameRunDet = (runIdDet == badRunList[iBadRun]);
      if (isSameRunPar) isGoodRunPar = false;
      if (isSameRunDet) isGoodRunDet = false;
    }

    // event cuts
    const Bool_t isInRvtxCutPar = (TMath::Abs(rVtxPar) < rVtxMax);
    const Bool_t isInZvtxCutPar = (TMath::Abs(zVtxPar) < zVtxMax);
    const Bool_t isInRvtxCutDet = (TMath::Abs(rVtxDet) < rVtxMax);
    const Bool_t isInZvtxCutDet = (TMath::Abs(zVtxDet) < zVtxMax);
    const Bool_t isGoodEvtPar   = (isGoodRunPar && isInRvtxCutPar && isInZvtxCutPar);
    const Bool_t isGoodEvtDet   = (isGoodRunDet && isInRvtxCutDet && isInZvtxCutDet);


    // particle trigger loop
    Double_t fTrgPar(0.);
    Double_t hTrgPar(0.);
    Double_t eTtrgPar(0.);
    UInt_t   iTrgPar(nTrkPar + 1);
    UInt_t   idTrgPar(0);
    Bool_t   isGoodTrgPar(true);
    Bool_t   foundTrgPar(false);
    if (_useParTrg) {
      for (UInt_t iTrgMC = 0; iTrgMC < nTrkPar; iTrgMC++) {

        // particle trigger info
        const Int_t    pidTrgMc  = mcIdGeant -> at(iTrgMC);
        const Int_t    idVxTrgMc = mcIdVx    -> at(iTrgMC);
        const Double_t qTrgMc    = mcCharge  -> at(iTrgMC);
        const Double_t fTrgMc    = mcPhi     -> at(iTrgMC);
        const Double_t hTrgMc    = mcEta     -> at(iTrgMC);
        const Double_t pTtrgMc   = mcPt      -> at(iTrgMC);

        // check pid
        Bool_t isHadronPar(false);
        switch (_fTrig) {
          case 0:
            for (UInt_t iTrgPion = 0; iTrgPion < NPionIds; iTrgPion++) {
              Bool_t isSameId = (pidTrgMc == idPion[iTrgPion]);
              if (isSameId) {
                isHadronPar = true;
                break;
              }
            }
            break;
          case 1:
            for (UInt_t iTrgChrg = 0; iTrgChrg < NChrgIds; iTrgChrg++) {
              Bool_t isSameId = (pidTrgMc == idChrg[iTrgChrg]);
              if (isSameId) {
                isHadronPar = true;
                break;
              }
            }
            break;
        }

        // determine trigger charge
        Bool_t isGoodChrg(false);
        switch (_fTrig) {
          case 0:
            isGoodChrg = (qTrgMc == 0.);
            break;
          case 1:
            isGoodChrg = (qTrgMc != 0.);
            break;
        }

        // particle trigger cuts
        const Bool_t isGoodStateTrgPar = (idVxTrgMc == 1);
        const Bool_t isChargedTrgPar   = isGoodChrg;
        const Bool_t isInTrgEtaCutPar  = (TMath::Abs(hTrgMc) < hTrgMax);
        const Bool_t isInTrgEtCutPar   = ((pTtrgMc > eTtrgMin) && (pTtrgMc < eTtrgMax));
        if (isHadronPar && isGoodStateTrgPar && isChargedTrgPar && isInTrgEtaCutPar && isInTrgEtCutPar) {
          iTrgPar     = iTrgMC;
          idTrgPar    = pidTrgMc;
          fTrgPar     = fTrgMc;
          hTrgPar     = hTrgMc;
          eTtrgPar    = pTtrgMc;
          foundTrgPar = true;
          break;
        }

      }  // end particle trigger loop
      if (!foundTrgPar)
        isGoodTrgPar = false;
      else
        isGoodTrgPar = true;

    }  // end triggered condition

    // detector trigger loop
    Double_t fTrgDet(0.);
    Double_t hTrgDet(0.);
    Double_t eTtrgDet(0.);
    UInt_t   iTrgDet(nTrkDet + 1);
    Bool_t   isGoodTrgDet(true);
    Bool_t   foundTrgDet(false);
    if (_useDetTrg) {
      for (UInt_t iTrgRE = 0; iTrgRE < nTrkDet; iTrgRE++) {

        // detector trigger info
        const UInt_t   nFitTrgRE  = PrimaryTrackArray_nHitsFit[iTrgRE];
        const UInt_t   nPossTrgRE = PrimaryTrackArray_nHitsPoss[iTrgRE];
        const Double_t rFitTrgRE  = (Double_t) nFitTrgRE / (Double_t) nPossTrgRE;
        const Double_t dcaTrgRE   = PrimaryTrackArray_dcag[iTrgRE];
        const Double_t fTrgRE     = PrimaryTrackArray_phi[iTrgRE];
        const Double_t hTrgRE     = PrimaryTrackArray_eta[iTrgRE];
        const Double_t pTtrgRE    = PrimaryTrackArray_pT[iTrgRE];
        const Double_t nSigElTrg  = TMath::Abs(PrimaryTrackArray_nSigElectron[iTrgRE]);
        const Double_t nSigPiTrg  = TMath::Abs(PrimaryTrackArray_nSigPion[iTrgRE]);
        const Double_t nSigKaTrg  = TMath::Abs(PrimaryTrackArray_nSigKaon[iTrgRE]);
        const Double_t nSigPrTrg  = TMath::Abs(PrimaryTrackArray_nSigProton[iTrgRE]);

        // determine trigger species
        const Bool_t isElectron  = (nSigElTrg < nSigTrgCut);
        const Bool_t isPion      = (nSigPiTrg < nSigTrgCut);
        const Bool_t isKaon      = (nSigKaTrg < nSigTrgCut);
        const Bool_t isProton    = (nSigPrTrg < nSigTrgCut);
        const Bool_t isHadronRE = ((isPion || isKaon || isProton) && !isElectron);

        // detector trigger cuts
        const Bool_t isInTrgFitCutRE   = (nFitTrgRE > nFitMin);
        const Bool_t isInTrgRatioCutRE = (rFitTrgRE > rFitMin);
        const Bool_t isInTrgDcaCutRE   = (dcaTrgRE < dcaMax);
        const Bool_t isInTrgEtaCutRE   = (TMath::Abs(hTrgRE) < hTrgMax);
        const Bool_t isInTrgEtCutRE    = ((pTtrgRE > eTtrgMin) && (pTtrgRE < eTtrgMax));
        if (isHadronRE && isInTrgFitCutRE && isInTrgRatioCutRE && isInTrgDcaCutRE && isInTrgEtaCutRE && isInTrgEtCutRE) {
          iTrgDet     = iTrgRE;
          fTrgDet     = fTrgRE;
          hTrgDet     = hTrgRE;
          eTtrgDet    = pTtrgRE;
          foundTrgDet = true;
        }

      }  // end detector trigger loop
      if (!foundTrgDet)
        isGoodTrgDet = false;
      else
        isGoodTrgDet = true;

    }  // end triggered condition

    // trigger cuts
    Bool_t useEvtPar = (isGoodEvtPar && isGoodTrgPar);
    Bool_t useEvtDet = (isGoodEvtDet && isGoodTrgDet);
    Bool_t useEvt    = (useEvtPar && useEvtDet);
    if (useEvtPar) nTrgPar++;
    if (useEvtDet) nTrgDet++;
    if (useEvt)
      nTrgTot++;
    else
      continue;


    UInt_t nJetTrkP(0);
    UInt_t nJetTrkD(0);
    UInt_t nJetTrkC(0);
    UInt_t nJetTrkM(0);
    UInt_t nJetTrkN(0);
    for (UInt_t iTrkMax = 0; iTrkMax < NTrkMax; iTrkMax++) {
      isMatched[iTrkMax] = false;
    }

    // particle track loop
    UInt_t nPar(0);
    for (UInt_t iTrkMC = 0; iTrkMC < nTrkPar; iTrkMC++) {

      // particle track info
      const Int_t    mcId  = mcIdTrk   -> at(iTrkMC);
      const Int_t    mcPID = mcIdGeant -> at(iTrkMC);
      const Int_t    mcIV  = mcIdVx    -> at(iTrkMC);
      const Double_t qMC   = mcCharge  -> at(iTrkMC);
      const Double_t fMC   = mcPhi     -> at(iTrkMC);
      const Double_t hMC   = mcEta     -> at(iTrkMC);
      const Double_t pTmc  = mcPt      -> at(iTrkMC);
      const Double_t pXmc  = mcPx      -> at(iTrkMC);
      const Double_t pYmc  = mcPy      -> at(iTrkMC);
      const Double_t pZmc  = mcPz      -> at(iTrkMC);
      const Double_t eMC   = mcEnergy  -> at(iTrkMC);

      Double_t dFmc = fMC - fTrgPar;
      if (dFmc < (-1. * TMath::PiOver2())) dFmc += TMath::TwoPi();
      if (dFmc > (3. * TMath::PiOver2()))  dFmc -= TMath::TwoPi();


      // consider only specific species (if need be)
      Bool_t isRelevantSpecies(false);
      for (UInt_t iSpecies = 0; iSpecies < NTrkSpecies; iSpecies++) {
        const Bool_t isSameSpecies = (mcPID == idTrkSpecies[iSpecies]);
        if (isSameSpecies) {
          isRelevantSpecies = true;
          break;
        }
      }
      if (_useSpecificTrackSpecies && !isRelevantSpecies) continue;

      // consider only recoil tracks (if need be)
      const Bool_t isAwaySidePar = ((dFmc > dFminAway) && (dFmc < dFmaxAway));
      if (_useOnlyRecoilTracks && !isAwaySidePar) continue;


      // particle track cuts
      const Bool_t isGoodEta   = (TMath::Abs(hMC) < hTrkMax);
      const Bool_t isGoodChrg  = (qMC != 0.);
      const Bool_t isGoodPt    = (pTmc > 0.);
      const Bool_t isNotTrgPar = (iTrkMC != iTrgPar);
      if (isGoodChrg && isGoodPt && isNotTrgPar) {

        // fill particle QA histograms
        hEtaQA[0] -> Fill(hMC);
        hPtQA[0]  -> Fill(pTmc);

        // fill particle histograms (after QA)
        if (isGoodEta) {
          hPhiTrk[0][1] -> Fill(fMC);
          hEtaTrk[0][1] -> Fill(hMC);
          hPtTrk[0][1]  -> Fill(pTmc);
          hPhiForEff[0] -> Fill(fMC);
          hEtaForEff[0] -> Fill(hMC);
          hPtForEff[0]  -> Fill(pTmc);
          parTracks.push_back(PseudoJet(pXmc, pYmc, pZmc, eMC));
          nJetTrkP++;
        }
        nPar++;


        // match to detector tracks
        for (UInt_t iTrkMat = 0; iTrkMat < nTrkDet; iTrkMat++) {

          // detector track info
          const Int_t    idMat     = (Int_t) PrimaryTrackArray_tofSigElectron[iTrkMat];
          const Int_t    qaMat     = (Int_t) PrimaryTrackArray_tofSigPion[iTrkMat];
          const Int_t    idVxMat   = (Int_t) PrimaryTrackArray_tofSigKaon[iTrkMat];
          const UInt_t   nFitMat   = PrimaryTrackArray_nHitsFit[iTrkMat];
          const UInt_t   nPossMat  = PrimaryTrackArray_nHitsPoss[iTrkMat];
          const UInt_t   nComMat   = (UInt_t) (nFitMat * ((Double_t) qaMat / 100.));
          const Double_t rFitMat   = (Double_t) nFitMat / (Double_t) nPossMat;
          const Double_t dcaMat    = PrimaryTrackArray_dcag[iTrkMat];
          const Double_t fMat      = PrimaryTrackArray_phi[iTrkMat];
          const Double_t hMat      = PrimaryTrackArray_eta[iTrkMat];
          const Double_t pTmat     = PrimaryTrackArray_pT[iTrkMat];
          const Double_t pXmat     = PrimaryTrackArray_pX[iTrkMat];
          const Double_t pYmat     = PrimaryTrackArray_pY[iTrkMat];
          const Double_t pZmat     = PrimaryTrackArray_pZ[iTrkMat];
          const Double_t eMat      = TMath::Sqrt((pTmat *pTmat) + (pZmat * pZmat) + (MassPi0 * MassPi0));
          const Double_t fDifMat   = fMC - fMat;
          const Double_t hDifMat   = hMC - hMat;
          const Double_t pTdifMat  = pTmc - pTmat;
          const Double_t pTfracMat = pTmat / pTmc;

          Double_t dFmat = fMat - fTrgPar;
          if (dFmat < (-1. * TMath::PiOver2())) dFmat += TMath::TwoPi();
          if (dFmat > (3. * TMath::PiOver2()))  dFmat -= TMath::TwoPi();


          // consider only recoil tracks (if need be)
          const Bool_t isAwaySideMat = ((dFmat > dFminAway) && (dFmat < dFmaxAway));
          if (_useOnlyRecoilTracks && !isAwaySideMat) continue;


          // detector track cuts
          const Bool_t isInFitCutMat   = (nFitMat > nFitMin);
          const Bool_t isInRatioCutMat = (rFitMat > rFitMin);
          const Bool_t isInDcaCutMat   = (dcaMat < dcaMax);
          const Bool_t isInEtaCutMat   = (TMath::Abs(hMat) < hTrkMax);
          const Bool_t isInPtCutMat    = ((pTmat > pTtrkMin) && (pTmat < pTtrkMax));
          const Bool_t isNotTrgMat     = (iTrkMat != iTrgDet);
          if (!isInFitCutMat || !isInRatioCutMat || !isInDcaCutMat || !isInEtaCutMat || !isInPtCutMat || !isNotTrgMat) continue;

          // fill detector QA histograms
          hEtaQA[1]   -> Fill(hMat);
          hPtQA[1]    -> Fill(pTmat);
          hNfitQA[0]  -> Fill(nFitMat);
          hRfitQA[0]  -> Fill(rFitMat);
          hDcaQA[0]   -> Fill(dcaMat);
          hPtFrac[0]  -> Fill(pTfracMat);
          hIdTruth[0] -> Fill(idMat);
          hQaTruth[0] -> Fill(qaMat);
          hNcommon[0] -> Fill(nComMat);


          // match track
          const Bool_t isMcMatch   = (idMat == mcId);
          const Bool_t isGoodMatch = (qaMat > minQaTruth);
          if (isMcMatch && isGoodMatch) {
            // fill matched QA histograms
            hEtaQA[2]   -> Fill(hMat);
            hPtQA[2]    -> Fill(pTmat);
            hNfitQA[1]  -> Fill(nFitMat);
            hRfitQA[1]  -> Fill(rFitMat);
            hDcaQA[1]   -> Fill(dcaMat);
            hPtFrac[1]  -> Fill(pTfracMat);
            hIdTruth[1] -> Fill(idMat);
            hQaTruth[1] -> Fill(qaMat);
            hNcommon[1] -> Fill(nComMat);

            // fill efficiency & resolution histograms
            hPhiForEff[1] -> Fill(fMC);
            hEtaForEff[1] -> Fill(hMC);
            hPtForEff[1]  -> Fill(pTmc);
            hPhiRes       -> Fill(fDifMat);
            hEtaRes       -> Fill(hDifMat);
            hPtRes        -> Fill(pTdifMat);
            hPhiResVsPhi  -> Fill(fMC, fDifMat);
            hEtaResVsEta  -> Fill(hMC, hDifMat);
            hPtResVsPt    -> Fill(pTmc, pTdifMat);
            pPhiResVsPhi  -> Fill(fMC, fDifMat);
            pEtaResVsEta  -> Fill(hMC, hDifMat);
            pPtResVsPt    -> Fill(pTmc, pTdifMat);

            // add particle vertex to list of vtxs / a match
            Bool_t mcVtxHasNotBeenAdded(true);
            if (matchMcVtxs.size() > 0) {
              for (UInt_t iMatchMcVtx = 0; iMatchMcVtx < matchMcVtxs.size(); iMatchMcVtx++) {
                Bool_t isSameMcVtx = (mcIV == matchMcVtxs[iMatchMcVtx]);
                if (isSameMcVtx) {
                  mcVtxHasNotBeenAdded = false;
                  break;
                }
              }  // end vertex loop
              if (mcVtxHasNotBeenAdded) {
                matchMcVtxs.push_back(mcIV);
              }
            }  // end if list has entries
            else {
              matchMcVtxs.push_back(mcIV);
            }

            // add detector vertex to list of vtxs w/ a match
            Bool_t recoVtxHasNotBeenAdded(true);
            if (matchRecoVtxs.size() > 0) {
              for (UInt_t iMatchRecoVtx = 0; iMatchRecoVtx < matchRecoVtxs.size(); iMatchRecoVtx++) {
                Bool_t isSameRecoVtx = (idVxMat == matchRecoVtxs[iMatchRecoVtx]);
                if (isSameRecoVtx) {
                  recoVtxHasNotBeenAdded = false;
                  break;
                }
              }  // end vertex loop
              if (recoVtxHasNotBeenAdded) {
                matchRecoVtxs.push_back(idVxMat);
              }
            }  // end if list has entries
            else {
              matchRecoVtxs.push_back(idVxMat);
            }

            // save for jet finding
            isMatched[iTrkMat] = true;
            detTracks.push_back(PseudoJet(pXmat, pYmat, pZmat, eMat));
            nJetTrkM++;
          }  // end matching condition
          else {
            // fill unmatched QA histograms
            hEtaQA[3]   -> Fill(hMat);
            hPtQA[3]    -> Fill(pTmat);
            hNfitQA[2]  -> Fill(nFitMat);
            hRfitQA[2]  -> Fill(rFitMat);
            hDcaQA[2]   -> Fill(dcaMat);
            hPtFrac[2]  -> Fill(pTfracMat);
            hIdTruth[2] -> Fill(idMat);
            hQaTruth[2] -> Fill(qaMat);
            hNcommon[2] -> Fill(nComMat);
          }
        }  // end detector track loop
      }  // end particle track cuts

      // fill particle histograms (before QA)
      if (isGoodChrg && isGoodPt && isNotTrgPar) {
        hPhiTrk[0][0] -> Fill(fMC);
        hEtaTrk[0][0] -> Fill(hMC);
        hPtTrk[0][0]  -> Fill(pTmc);
      }
    }  // end particle track loop

    // fill vertex histograms
    hNumMatchVtx[0] -> Fill(matchMcVtxs.size());
    hNumMatchVtx[1] -> Fill(matchRecoVtxs.size());
    for (UInt_t iMcVtx = 0; iMcVtx < matchMcVtxs.size(); iMcVtx++) {
      hMatchVtxId[0] -> Fill(matchMcVtxs[iMcVtx]);
    }
    for (UInt_t iRecoVtx = 0; iRecoVtx < matchRecoVtxs.size(); iRecoVtx++) {
      hMatchVtxId[1] -> Fill(matchRecoVtxs[iRecoVtx]);
    }


    // detector track loop
    UInt_t nDet(0);
    for (UInt_t iTrkDet = 0; iTrkDet < nTrkDet; iTrkDet++) {

      // detector track info
      const Int_t    idVxDet  = (Int_t) PrimaryTrackArray_tofSigKaon[iTrkDet];
      const UInt_t   nFitDet  = PrimaryTrackArray_nHitsFit[iTrkDet];
      const UInt_t   nPossDet = PrimaryTrackArray_nHitsPoss[iTrkDet];
      const Double_t rFitDet  = (Double_t) nFitDet / (Double_t) nPossDet;
      const Double_t dcaDet   = PrimaryTrackArray_dcag[iTrkDet];
      const Double_t fDet     = PrimaryTrackArray_phi[iTrkDet];
      const Double_t hDet     = PrimaryTrackArray_eta[iTrkDet];
      const Double_t pTdet    = PrimaryTrackArray_pT[iTrkDet];

      Double_t dFdet = fDet - fTrgPar;
      if (dFdet < (-1. * TMath::PiOver2())) dFdet += TMath::TwoPi();
      if (dFdet > (3. * TMath::PiOver2()))  dFdet -= TMath::TwoPi();


      // consider only recoil tracks (if need be)
      const Bool_t isAwaySideDet = ((dFdet > dFminAway) && (dFdet < dFmaxAway));
      if (_useOnlyRecoilTracks && !isAwaySideDet) continue;

      // fill detector histograms (before QA)
      hPhiTrk[1][0] -> Fill(fDet);
      hEtaTrk[1][0] -> Fill(hDet);
      hPtTrk[1][0]  -> Fill(pTdet);


      // detector track cuts
      const Bool_t isInFitCutDet   = (nFitDet > nFitMin);
      const Bool_t isInRatioCutDet = (rFitDet > rFitMin);
      const Bool_t isInDcaCutDet   = (dcaDet < dcaMax);
      const Bool_t isInEtaCutDet   = (TMath::Abs(hDet) < hTrkMax);
      const Bool_t isInPtCutDet    = ((pTdet > pTtrkMin) && (pTdet < pTtrkMax));
      const Bool_t isNotTrgDet     = (iTrkDet != iTrgDet);
      if (isInFitCutDet && isInRatioCutDet && isInDcaCutDet && isInEtaCutDet && isInPtCutDet && isNotTrgDet) {
        // fill vertex histograms
        if (isMatched[iTrkDet]) {
          hTrkVtxId[0] -> Fill(idVxDet);
          nJetTrkC++;
        }
        else {
          hTrkVtxId[1] -> Fill(idVxDet);
          nJetTrkN++;
        }

        // fill vertex pt histogram
        for (UInt_t iMatchVtx = 0; iMatchVtx < matchRecoVtxs.size(); iMatchVtx++) {
          const Bool_t isSameVtx = (idVxDet == matchRecoVtxs[iMatchVtx]);
          if (isSameVtx) {
            hPtMatchVtx -> Fill(pTdet);
          }
        }  // end vertex loop

        // fill detector histograms (after QA)
        hPhiTrk[1][1] -> Fill(fDet);
        hEtaTrk[1][1] -> Fill(hDet);
        hPtTrk[1][1]  -> Fill(pTdet);
        nJetTrkD++;
        nDet++;
      }  // end detector track cuts
    }  // end detector track loop

    // fill track number histograms
    hNumTrk[0]    -> Fill(nPar);
    hNumTrk[1]    -> Fill(nDet);
    hNumJetTrk[0] -> Fill(nJetTrkP);
    hNumJetTrk[1] -> Fill(nJetTrkD);
    hNumJetTrk[2] -> Fill(nJetTrkM);
    hNumJetTrk[3] -> Fill(nJetTrkN);


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
    _pEventIndex = evtIdPar;
    _pNJets      = (Int_t) pJets.size();
    _pRunId      = runIdPar;
    _pRefmult    = nPar;
    _pTSP        = idTrgPar;
    _pTrgEta     = hTrgPar;
    _pTrgPhi     = fTrgPar;
    _pTrgEt      = eTtrgPar;
    _pRho        = pRhoJet;
    _pSigma      = pSigJet;
    _pVz         = zVtxPar;

    _dEventIndex = evtIdDet;
    _dNJets      = (Int_t) dJets.size();
    _dRunId      = runIdDet;
    _dRefmult    = nDet;
    _dTSP        = Etsp;
    _dTrgEta     = hTrgPar;
    _dTrgPhi     = fTrgPar;
    _dTrgEt      = eTtrgPar;
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

    // clear vectors for jet-finding
    parTracks.clear();
    detTracks.clear();
    pJetsCS.clear();
    dJetsCS.clear();
    pJets.clear();
    dJets.clear();

    // clear vectors for vertices
    matchRecoVtxs.clear();
    matchMcVtxs.clear();

  }  // end event loop
  cout << "    Finished event loops!\n"
       << "      " << nTrgPar << " generated and " << nTrgDet << " reconstructed triggers.\n"
       << "      In total, " << nTrgTot << " events were used."
       << endl;


  // normalize by bin width (where relevant)
  const UInt_t nFbins = hPhiTrk[0][0] -> GetNbinsX();
  for (UInt_t iFbin = 1; iFbin < (nFbins + 1); iFbin++) {
    for (UInt_t iLvl = 0; iLvl < NLvls; iLvl++) {
      for (UInt_t iCut = 0; iCut < NCuts; iCut++) {
        const Double_t fBin1    = hPhiTrk[iLvl][iCut] -> GetBinWidth(iFbin);
        const Double_t oldValF1 = hPhiTrk[iLvl][iCut] -> GetBinContent(iFbin);
        const Double_t oldErrF1 = hPhiTrk[iLvl][iCut] -> GetBinError(iFbin);
        hPhiTrk[iLvl][iCut] -> SetBinContent(iFbin, oldValF1 / fBin1);
        hPhiTrk[iLvl][iCut] -> SetBinError(iFbin, oldErrF1 / fBin1);
      }  // end cut loop
      const Double_t fBin2    = hPhiForEff[iLvl] -> GetBinWidth(iFbin);
      const Double_t oldValF2 = hPhiForEff[iLvl] -> GetBinContent(iFbin);
      const Double_t oldErrF2 = hPhiForEff[iLvl] -> GetBinError(iFbin);
      hPhiForEff[iLvl] -> SetBinContent(iFbin, oldValF2 / fBin2);
      hPhiForEff[iLvl] -> SetBinError(iFbin, oldErrF2 / fBin2);
    }  // end level loop
  }  // end phi bin loop

  const UInt_t nHbins = hEtaTrk[0][0] -> GetNbinsX();
  for (UInt_t iHbin = 1; iHbin < (nHbins + 1); iHbin++) {
    for (UInt_t iLvl = 0; iLvl < NLvls; iLvl++) {
      for (UInt_t iCut = 0; iCut < NCuts; iCut++) {
        const Double_t hBin1    = hEtaTrk[iLvl][iCut] -> GetBinWidth(iHbin);
        const Double_t oldValH1 = hEtaTrk[iLvl][iCut] -> GetBinContent(iHbin);
        const Double_t oldErrH1 = hEtaTrk[iLvl][iCut] -> GetBinError(iHbin);
        hEtaTrk[iLvl][iCut] -> SetBinContent(iHbin, oldValH1 / hBin1);
        hEtaTrk[iLvl][iCut] -> SetBinError(iHbin, oldErrH1 / hBin1);
      }  // end cut loop
      const Double_t hBin2    = hEtaForEff[iLvl] -> GetBinWidth(iHbin);
      const Double_t oldValH2 = hEtaForEff[iLvl] -> GetBinContent(iHbin);
      const Double_t oldErrH2 = hEtaForEff[iLvl] -> GetBinError(iHbin);
      hEtaForEff[iLvl] -> SetBinContent(iHbin, oldValH2 / hBin2);
      hEtaForEff[iLvl] -> SetBinError(iHbin, oldErrH2 / hBin2);
    }  // end level loop
  }  // end eta bin loop

  const UInt_t nPtBins = hPtTrk[0][0] -> GetNbinsX();
  for (UInt_t iPtBin = 1; iPtBin < (nPtBins + 1); iPtBin++) {
    for (UInt_t iLvl = 0; iLvl < NLvls; iLvl++) {
      for (UInt_t iCut = 0; iCut < NCuts; iCut++) {
        const Double_t pTbin1   = hPtTrk[iLvl][iCut] -> GetBinWidth(iPtBin);
        const Double_t oldValP1 = hPtTrk[iLvl][iCut] -> GetBinContent(iPtBin);
        const Double_t oldErrP1 = hPtTrk[iLvl][iCut] -> GetBinError(iPtBin);
        hPtTrk[iLvl][iCut] -> SetBinContent(iPtBin, oldValP1 / pTbin1);
        hPtTrk[iLvl][iCut] -> SetBinError(iPtBin, oldErrP1 / pTbin1);
      }  // end cut loop
      const Double_t pTbin2   = hPtForEff[iLvl] -> GetBinWidth(iPtBin);
      const Double_t oldValP2 = hPtForEff[iLvl] -> GetBinContent(iPtBin);
      const Double_t oldErrP2 = hPtForEff[iLvl] -> GetBinError(iPtBin);
      hPtForEff[iLvl] -> SetBinContent(iPtBin, oldValP2 / pTbin2);
      hPtForEff[iLvl] -> SetBinError(iPtBin, oldErrP2 / pTbin2);
    }  // end level loop
    for (UInt_t iHistQA = 0; iHistQA < NHistQA; iHistQA++) {
      const Double_t pTbin3   = hPtQA[iHistQA] -> GetBinWidth(iPtBin);
      const Double_t oldValP3 = hPtQA[iHistQA] -> GetBinContent(iPtBin);
      const Double_t oldErrP3 = hPtQA[iHistQA] -> GetBinError(iPtBin);
      hPtQA[iHistQA] -> SetBinContent(iPtBin, oldValP3 / pTbin3);
      hPtQA[iHistQA] -> SetBinError(iPtBin, oldErrP3 / pTbin3);
    }  // end QA histogram loop
  }  // end pT bin loop

  const UInt_t nPtMatchBins = hPtMatchVtx -> GetNbinsX();
  for (UInt_t iPtMatchBin = 1; iPtMatchBin < (nPtMatchBins + 1); iPtMatchBin++) {
    const Double_t pTbinMatch = hPtMatchVtx -> GetBinWidth(iPtMatchBin);
    const Double_t oldValMatch = hPtMatchVtx -> GetBinContent(iPtMatchBin);
    const Double_t oldErrMatch = hPtMatchVtx -> GetBinError(iPtMatchBin);
    hPtMatchVtx -> SetBinContent(iPtMatchBin, oldValMatch / pTbinMatch);
    hPtMatchVtx -> SetBinError(iPtMatchBin, oldErrMatch / pTbinMatch);
  }  // end pT bin loop
  cout << "    Corrected for bin widths." << endl;


  // calculate efficiencies
  const Double_t parWeight(1.);
  const Double_t detWeight(1.);
  hPhiEff -> Divide(hPhiForEff[1], hPhiForEff[0], detWeight, parWeight);
  hEtaEff -> Divide(hEtaForEff[1], hEtaForEff[0], detWeight, parWeight);
  hPtEff  -> Divide(hPtForEff[1], hPtForEff[0], detWeight, parWeight);
  cout << "    Calculated efficiencies." << endl;


  // make directories and save histograms
  TString     sLvlDir[NLvls + 2] = {"particle", "detector", "resolutions", "QA"};
  TDirectory *dLvls[NLvls + 2];
  for (UInt_t iLvl = 0; iLvl < NLvls; iLvl++) {
    dLvls[iLvl] = (TDirectory*) _fOut -> mkdir(sLvlDir[iLvl].Data());
    dLvls[iLvl] -> cd();
    for(UInt_t iCut = 0; iCut < NCuts; iCut++) {
      hPhiTrk[iLvl][iCut] -> Write();
      hEtaTrk[iLvl][iCut] -> Write();
      hPtTrk[iLvl][iCut]  -> Write();
    }
    hPhiForEff[iLvl] -> Write();
    hEtaForEff[iLvl] -> Write();
    hPtForEff[iLvl]  -> Write();
    hNumTrk[iLvl]    -> Write();
  }
  dLvls[NLvls]     = (TDirectory*) _fOut -> mkdir(sLvlDir[NLvls].Data());
  dLvls[NLvls + 1] = (TDirectory*) _fOut -> mkdir(sLvlDir[NLvls + 1].Data());
  dLvls[NLvls]     -> cd();
  hPhiRes          -> Write();
  hEtaRes          -> Write();
  hPtRes           -> Write();
  hPhiResVsPhi     -> Write();
  hEtaResVsEta     -> Write();
  hPtResVsPt       -> Write();
  pPhiResVsPhi     -> Write();
  pEtaResVsEta     -> Write();
  pPtResVsPt       -> Write();
  dLvls[NLvls + 1] -> cd();
  for (UInt_t iHistQA = 0; iHistQA < NHistQA; iHistQA++) {
    hEtaQA[iHistQA] -> Write();
    hPtQA[iHistQA]  -> Write();
    if (iHistQA < (NHistQA - 1)) {
      hNfitQA[iHistQA]  -> Write();
      hRfitQA[iHistQA]  -> Write();
      hDcaQA[iHistQA]   -> Write();
      hPtFrac[iHistQA]  -> Write();
      hIdTruth[iHistQA] -> Write();
      hQaTruth[iHistQA] -> Write();
      hNcommon[iHistQA] -> Write();
    }
  }
  dLvls[NLvls + 1] -> cd();
  hNumJetTrk[0]    -> Write();
  hNumJetTrk[1]    -> Write();
  hNumJetTrk[2]    -> Write();
  hNumJetTrk[3]    -> Write();
  hNumMatchVtx[0]  -> Write();
  hNumMatchVtx[1]  -> Write();
  hMatchVtxId[0]   -> Write();
  hMatchVtxId[1]   -> Write();
  hTrkVtxId[0]     -> Write();
  hTrkVtxId[1]     -> Write();
  dLvls[1]         -> cd();
  hPtMatchVtx      -> Write();
  _fOut            -> cd();
  hPhiEff          -> Write();
  hEtaEff          -> Write();
  hPtEff           -> Write();
  cout << "    Made directories and saved histograms." << endl;

}  // end 'Make(UInt_t, Bool_t, Bool_t, Bool_t, Bool_t)'



void StTrackEfficiencyCalculator::Finish() {

  // save output trees
  _fOut    -> cd();
  _tJetPar -> Write();
  _tJetDet -> Write();
  cout << "    Saved output trees." << endl;

  // close files
  _fOut -> cd();
  _fOut -> Close();
  _fIn  -> cd();
  _fIn  -> Close();
  cout << "  Finished efficiency calculation!\n" << endl;

}  // end 'Finish()'



void StTrackEfficiencyCalculator::SetJetParameters(const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t pTjetMin, const Double_t aGhost, const Double_t hJetMax, const Double_t hBkgdMax, const Double_t hGhostMax) {

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

// End ------------------------------------------------------------------------
