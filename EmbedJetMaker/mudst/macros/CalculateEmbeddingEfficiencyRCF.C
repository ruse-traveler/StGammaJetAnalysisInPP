// 'CalculateEmbeddingEfficiency.C'
// Derek Anderson
// 03.08.2018
//
// Use this to calculate the tracking
// efficiency using output from the
// 'StJetTreeMcMaker'.  Called by
// 'CalculateEfficiency.C'.


#include <map>
#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t  NTrkMax(5000);
static const UInt_t  NTwrMax(5000);
static const UInt_t  NMatchMax(10);
static const UInt_t  NHotTwr(41);
static const UInt_t  NBadRuns(45);
static const UInt_t  NPhiBins(60);
static const UInt_t  NEtaBins(40);
static const UInt_t  NPtBins(48);
static const UInt_t  NLevel(2);
static const UInt_t  NTrgs(2);
static const UInt_t  NCut(2);
static const UInt_t  NVal(2);
static const TString sTreePar("McTracks");
static const TString sTreeDet("GfmtoDst_mu");
static const TString sInputDefault("../../MuDstMatching/output/merged/pt35rff.matchWithMc.root");
static const TString sOutputDefault("pp200r9pt35.default.root");



void CalculateEmbeddingEfficiency(UInt_t &nPi0Trg, UInt_t &nGamTrg, const Bool_t isInBatchMode=false, const Bool_t isTriggered=false, const TString sInput=sInputDefault, const TString sOutput=sOutputDefault) {

  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);

  // trigger paramters
  const Int_t    adcMax(6004);
  const Double_t eStrMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t hTrgMax(0.9);
  const Double_t eTtrgMin(9.);
  const Double_t eTtrgMax(20.);
  const Double_t tspPi0[2] = {0., 0.08};
  const Double_t tspGam[2] = {0.2, 0.6};

  // track parameters
  const Int_t    minQaTruth(50);
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(100.);
  cout << "\n  Beginning efficiency calculation..." << endl;

  // TEST [01.21.2019]
  const Double_t phiMin(0.);
  const Double_t phiMax(TMath::Pi() / 2.);


  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!" << endl;
    assert(0);
  }
  cout << "    Opened files." << endl;

  // grab input tree
  TTree *tPar;
  TTree *tDet;
  fInput -> GetObject(sTreePar.Data(), tPar);
  fInput -> GetObject(sTreeDet.Data(), tDet);
  if (!tPar) {
    cerr << "PANIC: couldn't grab particle level tree!" << endl;
    assert(0);
  }
  if (!tDet) {
    cerr << "PANIC: couldn't grab detector level tree!" << endl;
    assert(0);
  }
  cout << "    Grabbed trees." << endl;


  // declare particle (event) leaf types
  Int_t    mcEventId;
  Int_t    mcRunId;
  Int_t    mcNumTrks;
  Double_t muVx;
  Double_t muVy;
  Double_t muVz;
  Double_t mcVx;
  Double_t mcVy;
  Double_t mcVz;
  // declare particle (track) leaf types
  vector<Int_t>    *mcIdTrk;
  vector<Int_t>    *mcIdGeant;
  vector<Int_t>    *mcIdVx;
  vector<Int_t>    *mcIdVxEnd;
  vector<Int_t>    *mcIntrVtx;
  vector<Bool_t>   *mcIsShower;
  vector<Double_t> *mcCharge;
  vector<Double_t> *mcRapidity;
  vector<Double_t> *mcEta;
  vector<Double_t> *mcPhi;
  vector<Double_t> *mcPx;
  vector<Double_t> *mcPy;
  vector<Double_t> *mcPz;
  vector<Double_t> *mcPt;
  vector<Double_t> *mcPtot;
  vector<Double_t> *mcEnergy;

  // declare detector leaf types
  UInt_t   fUniqueID;
  UInt_t   fBits;
  Long64_t runNumber;
  Long64_t eventNumber;
  Int_t    trigID;
  Int_t    nGlobalTracks;
  Int_t    nPrimaryTracks;
  Int_t    refMult;
  Double_t vpdVz;
  Double_t xVertex;
  Double_t yVertex;
  Double_t zVertex;
  Double_t bbcZVertex;
  Double_t zdcCoincidenceRate;
  Double_t bbcCoincidenceRate;
  Double_t backgroundRate;
  Double_t bbcBlueBackgroundRate;
  Double_t bbcYellowBackgroundRate;
  Double_t refMultPos;
  Double_t refMultNeg;
  Double_t bTOFTrayMultiplicity;
  Int_t    nVerticies;
  Double_t MagF;
  Double_t VrtxRank;
  Float_t  Etsp;
  Int_t    ETwrdidT;
  Int_t    ETwradc11;
  Float_t  ETwreneT0;
  Float_t  ETwreT;
  Float_t  ETwrENET0;
  Float_t  ETwrphT;
  Float_t  ETwrPTower;
  Float_t  ETwrpidTower;
  Int_t    ETwrmoduleT;
  Float_t  EClustEneT0;
  Float_t  EClustetav1;
  Float_t  EClustphiv1;
  Float_t  EEstrpen01;
  Float_t  EEstrpen02;
  Float_t  EEstrpen03;
  Float_t  EEstrpen0;
  Float_t  EEstrpen1;
  Float_t  EEstrpen2;
  Float_t  EEstrpen3;
  Float_t  EEstrpen4;
  Float_t  EEstrpen5;
  Float_t  EEstrpen6;
  Float_t  EEstrpen7;
  Float_t  EEstrpen8;
  Float_t  EEstrpen9;
  Float_t  EEstrpen10;
  Float_t  EEstrpen11;
  Float_t  EEstrpen12;
  Float_t  EEstrpen13;
  Float_t  EEstrpen14;
  Float_t  EEstrpen15;
  Int_t    ETwrdidE;
  Float_t  EPstripenp01;
  Float_t  EPstripenp02;
  Float_t  EPstripenp03;
  Float_t  EPstripenp0;
  Float_t  EPstripenp1;
  Float_t  EPstripenp2;
  Float_t  EPstripenp3;
  Float_t  EPstripenp4;
  Float_t  EPstripenp5;
  Float_t  EPstripenp6;
  Float_t  EPstripenp7;
  Float_t  EPstripenp8;
  Float_t  EPstripenp9;
  Float_t  EPstripenp10;
  Float_t  EPstripenp11;
  Float_t  EPstripenp12;
  Float_t  EPstripenp13;
  Float_t  EPstripenp14;
  Float_t  EPstripenp15;
  Float_t  EclustEnnq1;
  Float_t  EclustEnnq20;
  Float_t  EclustEnnq19;
  Float_t  EclustEnpq1;
  Float_t  EclustEnpq20;
  Float_t  EclustEnpq19;
  Float_t  EclustEnpq21;
  Int_t    PrimaryTrackArray_;
  UInt_t   PrimaryTrackArray_fUniqueID[NTrkMax];
  UInt_t   PrimaryTrackArray_fBits[NTrkMax];
  Int_t    PrimaryTrackArray_nHitsFit[NTrkMax];
  Int_t    PrimaryTrackArray_nHitsPoss[NTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
  Int_t    PrimaryTrackArray_pdgId[NTrkMax];
  Int_t    PrimaryTrackArray_geantId[NTrkMax];
  Double_t PrimaryTrackArray_pZ[NTrkMax];
  Double_t PrimaryTrackArray_pX[NTrkMax];
  Double_t PrimaryTrackArray_pY[NTrkMax];
  Double_t PrimaryTrackArray_pT[NTrkMax];
  Double_t PrimaryTrackArray_dEdx[NTrkMax];
  Double_t PrimaryTrackArray_charge[NTrkMax];
  Double_t PrimaryTrackArray_tofBeta[NTrkMax];
  Double_t PrimaryTrackArray_eta[NTrkMax];
  Double_t PrimaryTrackArray_phi[NTrkMax];
  Double_t PrimaryTrackArray_nSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_nSigPion[NTrkMax];
  Double_t PrimaryTrackArray_nSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_nSigProton[NTrkMax];
  Double_t PrimaryTrackArray_dcag[NTrkMax];
  Double_t PrimaryTrackArray_nHits[NTrkMax];
  Double_t PrimaryTrackArray_dEdxHits[NTrkMax];
  Double_t PrimaryTrackArray_firstZPoint[NTrkMax];
  Double_t PrimaryTrackArray_lastZPoint[NTrkMax];
  Double_t PrimaryTrackArray_tofSigElectron[NTrkMax];
  Double_t PrimaryTrackArray_tofSigPion[NTrkMax];
  Double_t PrimaryTrackArray_tofSigKaon[NTrkMax];
  Double_t PrimaryTrackArray_tofSigProton[NTrkMax];
  Double_t PrimaryTrackArray_timeOfflight[NTrkMax];
  Double_t PrimaryTrackArray_pathLength[NTrkMax];
  Int_t    PrimaryTrackArray_trkIndex[NTrkMax];
  Int_t    TowerArray_;
  UInt_t   TowerArray_fUniqueID[NTwrMax];
  UInt_t   TowerArray_fBits[NTwrMax];
  Int_t    TowerArray_TwrId[NTwrMax];
  Float_t  TowerArray_TwrEng[NTwrMax];
  Float_t  TowerArray_TwrEta[NTwrMax];
  Float_t  TowerArray_TwrPhi[NTwrMax];
  Float_t  TowerArray_TwrADC[NTwrMax];
  Float_t  TowerArray_TwrPed[NTwrMax];
  Float_t  TowerArray_TwrRMS[NTwrMax];
  Int_t    TowerArray_TwrMatchIdnex[NTwrMax];
  Int_t    TowerArray_NoOfmatchedTrk[NTwrMax];
  Float_t  TowerArray_TwrMatchP[NTwrMax];
  Float_t  TowerArray_TwrPx[NTwrMax];
  Float_t  TowerArray_TwrPy[NTwrMax];
  Float_t  TowerArray_TwrPz[NTwrMax];
  Int_t    TowerArray_fNAssocTracks[NTwrMax];
  Int_t    TowerArray_fMatchedTracksArray_[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_P[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigPi[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigK[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigP[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_nSigE[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_dcag[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_eta[NTwrMax][NMatchMax];
  Float_t  TowerArray_fMatchedTracksArray_pT[NTwrMax][NMatchMax];
  Int_t    TowerArray_fMatchedTracksArray_nFit[NTwrMax][NMatchMax];
  Int_t    TowerArray_fMatchedTracksArray_nPos[NTwrMax][NMatchMax];


  // set particle branch addresses
  tPar -> SetBranchAddress("EventId", &mcEventId);
  tPar -> SetBranchAddress("RunId", &mcRunId);
  tPar -> SetBranchAddress("NumTrks", &mcNumTrks);
  tPar -> SetBranchAddress("MuVtxX", &muVx);
  tPar -> SetBranchAddress("MuVtxY", &muVy);
  tPar -> SetBranchAddress("MuVtxZ", &muVz);
  tPar -> SetBranchAddress("McVtxX", &mcVx);
  tPar -> SetBranchAddress("McVtxY", &mcVy);
  tPar -> SetBranchAddress("McVtxZ", &mcVz);
  tPar -> SetBranchAddress("IdTrk", &mcIdTrk);
  tPar -> SetBranchAddress("IdGeant", &mcIdGeant);
  tPar -> SetBranchAddress("IdVx", &mcIdVx);
  tPar -> SetBranchAddress("IdVxEnd", &mcIdVxEnd);
  tPar -> SetBranchAddress("IntrVtx", &mcIntrVtx);
  tPar -> SetBranchAddress("IsShower", &mcIsShower);
  tPar -> SetBranchAddress("Charge", &mcCharge);
  tPar -> SetBranchAddress("Rapidity", &mcRapidity);
  tPar -> SetBranchAddress("Eta", &mcEta);
  tPar -> SetBranchAddress("Phi", &mcPhi);
  tPar -> SetBranchAddress("Px", &mcPx);
  tPar -> SetBranchAddress("Py", &mcPy);
  tPar -> SetBranchAddress("Pz", &mcPz);
  tPar -> SetBranchAddress("Pt", &mcPt);
  tPar -> SetBranchAddress("Ptot", &mcPtot);
  tPar -> SetBranchAddress("Energy", &mcEnergy);

  // set detector branch addresses
  tDet -> SetMakeClass(1);
  tDet -> SetBranchAddress("fUniqueID", &fUniqueID);
  tDet -> SetBranchAddress("fBits", &fBits);
  tDet -> SetBranchAddress("runNumber", &runNumber);
  tDet -> SetBranchAddress("eventNumber", &eventNumber);
  tDet -> SetBranchAddress("trigID", &trigID);
  tDet -> SetBranchAddress("nGlobalTracks", &nGlobalTracks);
  tDet -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks);
  tDet -> SetBranchAddress("refMult", &refMult);
  tDet -> SetBranchAddress("vpdVz", &vpdVz);
  tDet -> SetBranchAddress("xVertex", &xVertex);
  tDet -> SetBranchAddress("yVertex", &yVertex);
  tDet -> SetBranchAddress("zVertex", &zVertex);
  tDet -> SetBranchAddress("bbcZVertex", &bbcZVertex);
  tDet -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate);
  tDet -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate);
  tDet -> SetBranchAddress("backgroundRate", &backgroundRate);
  tDet -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate);
  tDet -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate);
  tDet -> SetBranchAddress("refMultPos", &refMultPos);
  tDet -> SetBranchAddress("refMultNeg", &refMultNeg);
  tDet -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity);
  tDet -> SetBranchAddress("nVerticies", &nVerticies);
  tDet -> SetBranchAddress("MagF", &MagF);
  tDet -> SetBranchAddress("VrtxRank", &VrtxRank);
  tDet -> SetBranchAddress("Etsp", &Etsp);
  tDet -> SetBranchAddress("ETwrdidT", &ETwrdidT);
  tDet -> SetBranchAddress("ETwradc11", &ETwradc11);
  tDet -> SetBranchAddress("ETwreneT0", &ETwreneT0);
  tDet -> SetBranchAddress("ETwreT", &ETwreT);
  tDet -> SetBranchAddress("ETwrENET0", &ETwrENET0);
  tDet -> SetBranchAddress("ETwrphT", &ETwrphT);
  tDet -> SetBranchAddress("ETwrPTower", &ETwrPTower);
  tDet -> SetBranchAddress("ETwrpidTower", &ETwrpidTower);
  tDet -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT);
  tDet -> SetBranchAddress("EClustEneT0", &EClustEneT0);
  tDet -> SetBranchAddress("EClustetav1", &EClustetav1);
  tDet -> SetBranchAddress("EClustphiv1", &EClustphiv1);
  tDet -> SetBranchAddress("EEstrpen01", &EEstrpen01);
  tDet -> SetBranchAddress("EEstrpen02", &EEstrpen02);
  tDet -> SetBranchAddress("EEstrpen03", &EEstrpen03);
  tDet -> SetBranchAddress("EEstrpen0", &EEstrpen0);
  tDet -> SetBranchAddress("EEstrpen1", &EEstrpen1);
  tDet -> SetBranchAddress("EEstrpen2", &EEstrpen2);
  tDet -> SetBranchAddress("EEstrpen3", &EEstrpen3);
  tDet -> SetBranchAddress("EEstrpen4", &EEstrpen4);
  tDet -> SetBranchAddress("EEstrpen5", &EEstrpen5);
  tDet -> SetBranchAddress("EEstrpen6", &EEstrpen6);
  tDet -> SetBranchAddress("EEstrpen7", &EEstrpen7);
  tDet -> SetBranchAddress("EEstrpen8", &EEstrpen8);
  tDet -> SetBranchAddress("EEstrpen9", &EEstrpen9);
  tDet -> SetBranchAddress("EEstrpen10", &EEstrpen10);
  tDet -> SetBranchAddress("EEstrpen11", &EEstrpen11);
  tDet -> SetBranchAddress("EEstrpen12", &EEstrpen12);
  tDet -> SetBranchAddress("EEstrpen13", &EEstrpen13);
  tDet -> SetBranchAddress("EEstrpen14", &EEstrpen14);
  tDet -> SetBranchAddress("EEstrpen15", &EEstrpen15);
  tDet -> SetBranchAddress("ETwrdidE", &ETwrdidE);
  tDet -> SetBranchAddress("EPstripenp01", &EPstripenp01);
  tDet -> SetBranchAddress("EPstripenp02", &EPstripenp02);
  tDet -> SetBranchAddress("EPstripenp03", &EPstripenp03);
  tDet -> SetBranchAddress("EPstripenp0", &EPstripenp0);
  tDet -> SetBranchAddress("EPstripenp1", &EPstripenp1);
  tDet -> SetBranchAddress("EPstripenp2", &EPstripenp2);
  tDet -> SetBranchAddress("EPstripenp3", &EPstripenp3);
  tDet -> SetBranchAddress("EPstripenp4", &EPstripenp4);
  tDet -> SetBranchAddress("EPstripenp5", &EPstripenp5);
  tDet -> SetBranchAddress("EPstripenp6", &EPstripenp6);
  tDet -> SetBranchAddress("EPstripenp7", &EPstripenp7);
  tDet -> SetBranchAddress("EPstripenp8", &EPstripenp8);
  tDet -> SetBranchAddress("EPstripenp9", &EPstripenp9);
  tDet -> SetBranchAddress("EPstripenp10", &EPstripenp10);
  tDet -> SetBranchAddress("EPstripenp11", &EPstripenp11);
  tDet -> SetBranchAddress("EPstripenp12", &EPstripenp12);
  tDet -> SetBranchAddress("EPstripenp13", &EPstripenp13);
  tDet -> SetBranchAddress("EPstripenp14", &EPstripenp14);
  tDet -> SetBranchAddress("EPstripenp15", &EPstripenp15);
  tDet -> SetBranchAddress("EclustEnnq1", &EclustEnnq1);
  tDet -> SetBranchAddress("EclustEnnq20", &EclustEnnq20);
  tDet -> SetBranchAddress("EclustEnnq19", &EclustEnnq19);
  tDet -> SetBranchAddress("EclustEnpq1", &EclustEnpq1);
  tDet -> SetBranchAddress("EclustEnpq20", &EclustEnpq20);
  tDet -> SetBranchAddress("EclustEnpq19", &EclustEnpq19);
  tDet -> SetBranchAddress("EclustEnpq21", &EclustEnpq21);
  tDet -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_);
  tDet -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID);
  tDet -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits);
  tDet -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit);
  tDet -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss);
  tDet -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag);
  tDet -> SetBranchAddress("PrimaryTrackArray.pdgId", PrimaryTrackArray_pdgId);
  tDet -> SetBranchAddress("PrimaryTrackArray.geantId", PrimaryTrackArray_geantId);
  tDet -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ);
  tDet -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX);
  tDet -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY);
  tDet -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT);
  tDet -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx);
  tDet -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge);
  tDet -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta);
  tDet -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta);
  tDet -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi);
  tDet -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron);
  tDet -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion);
  tDet -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon);
  tDet -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton);
  tDet -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag);
  tDet -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits);
  tDet -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits);
  tDet -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint);
  tDet -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint);
  tDet -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron);
  tDet -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion);
  tDet -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon);
  tDet -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton);
  tDet -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight);
  tDet -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength);
  tDet -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex);
  tDet -> SetBranchAddress("TowerArray", &TowerArray_);
  tDet -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID);
  tDet -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits);
  tDet -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId);
  tDet -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng);
  tDet -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta);
  tDet -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi);
  tDet -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC);
  tDet -> SetBranchAddress("TowerArray.TwrPed", TowerArray_TwrPed);
  tDet -> SetBranchAddress("TowerArray.TwrRMS", TowerArray_TwrRMS);
  tDet -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex);
  tDet -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk);
  tDet -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP);
  tDet -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx);
  tDet -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy);
  tDet -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz);
  tDet -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigPi[10]", TowerArray_fMatchedTracksArray_nSigPi);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigK[10]", TowerArray_fMatchedTracksArray_nSigK);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigP[10]", TowerArray_fMatchedTracksArray_nSigP);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nSigE[10]", TowerArray_fMatchedTracksArray_nSigE);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_dcag[10]", TowerArray_fMatchedTracksArray_dcag);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_eta[10]", TowerArray_fMatchedTracksArray_eta);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_pT[10]", TowerArray_fMatchedTracksArray_pT);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nFit[10]", TowerArray_fMatchedTracksArray_nFit);
  tDet -> SetBranchAddress("TowerArray.fMatchedTracksArray_nPos[10]", TowerArray_fMatchedTracksArray_nPos);
  cout << "    Branches set." << endl;



  // define bad run and hot tower lists
  const UInt_t badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  const UInt_t hotTwrList[NHotTwr] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};
  cout << "    Bad run and hot tower lists defined:\n"
       << "      " << NBadRuns << " bad runs, " << NHotTwr << " hot towers."
       << endl;


  // define particle/detector histograms
  TH1D *hPhiTrk[NLevel][NTrgs][NCut];
  TH1D *hEtaTrk[NLevel][NTrgs][NCut];
  TH1D *hPtTrk[NLevel][NTrgs][NCut];
  // define matching histograms
  TH1D *hPhiForEff[NLevel][NTrgs];
  TH1D *hEtaForEff[NLevel][NTrgs];
  TH1D *hPtForEff[NLevel][NTrgs];
  TH2D *hPhiParVsDet[NTrgs];
  TH2D *hEtaParVsDet[NTrgs];
  TH2D *hPtParVsDet[NTrgs];
  // define efficiency histograms
  TH1D *hPhiEff[NTrgs];
  TH1D *hEtaEff[NTrgs];
  TH1D *hPtEff[NTrgs];
  // define resolution histograms
  TH1D *hPhiRes[NTrgs];
  TH1D *hEtaRes[NTrgs];
  TH1D *hPtRes[NTrgs];
  TH2D *hPhiResVsPhi[NTrgs];
  TH2D *hEtaResVsEta[NTrgs];
  TH2D *hPtResVsPt[NTrgs];

  // pT binning
  const UInt_t   nPtHistBins         = NPtBins - 1;
  const Double_t pTbinEdges[NPtBins] = {0., 0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.26, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 7., 8., 10., 12., 14., 16., 18., 20.};

  // binning
  const UInt_t   nDf     = 240;
  const UInt_t   nDh     = 400;
  const UInt_t   nDpt    = 1000;
  const UInt_t   nPt2d   = 1000;
  //const UInt_t   nPt2d   = 150;
  const Double_t f[2]    = {-3.15, 3.15};
  const Double_t h[2]    = {-1., 1.};
  const Double_t dF[2]   = {-12.56, 12.56};
  const Double_t dH[2]   = {-10., 10.};
  const Double_t dPt[2]  = {-100., 100.};
  const Double_t pT2d[2] = {0., 100.};
  // create particle histograms
  hPhiTrk[0][0][0] = new TH1D("hPhiBeforeQA_piPar", "", NPhiBins, f[0], f[1]);
  hPhiTrk[0][1][0] = new TH1D("hPhiBeforeQA_gaPar", "", NPhiBins, f[0], f[1]);
  hPhiTrk[0][0][1] = new TH1D("hPhiAfterQA_piPar", "", NPhiBins, f[0], f[1]);
  hPhiTrk[0][1][1] = new TH1D("hPhiAfterQA_gaPar", "", NPhiBins, f[0], f[1]);
  hEtaTrk[0][0][0] = new TH1D("hEtaBeforeQA_piPar", "", NEtaBins, h[0], h[1]);
  hEtaTrk[0][1][0] = new TH1D("hEtaBeforeQA_gaPar", "", NEtaBins, h[0], h[1]);
  hEtaTrk[0][0][1] = new TH1D("hEtaAfterQA_piPar", "", NEtaBins, h[0], h[1]);
  hEtaTrk[0][1][1] = new TH1D("hEtaAfterQA_gaPar", "", NEtaBins, h[0], h[1]);
  hPtTrk[0][0][0]  = new TH1D("hPtBeforeQA_piPar", "", nPtHistBins, pTbinEdges);
  hPtTrk[0][1][0]  = new TH1D("hPtBeforeQA_gaPar", "", nPtHistBins, pTbinEdges);
  hPtTrk[0][0][1]  = new TH1D("hPtAfterQA_piPar", "", nPtHistBins, pTbinEdges);
  hPtTrk[0][1][1]  = new TH1D("hPtAfterQA_gaPar", "", nPtHistBins, pTbinEdges);
  // create detector histograms
  hPhiTrk[1][0][0] = new TH1D("hPhiBeforeQA_piDet", "", NPhiBins, f[0], f[1]);
  hPhiTrk[1][1][0] = new TH1D("hPhiBeforeQA_gaDet", "", NPhiBins, f[0], f[1]);
  hPhiTrk[1][0][1] = new TH1D("hPhiAfterQA_piDet", "", NPhiBins, f[0], f[1]);
  hPhiTrk[1][1][1] = new TH1D("hPhiAfterQA_gaDet", "", NPhiBins, f[0], f[1]);
  hEtaTrk[1][0][0] = new TH1D("hEtaBeforeQA_piDet", "", NEtaBins, h[0], h[1]);
  hEtaTrk[1][1][0] = new TH1D("hEtaBeforeQA_gaDet", "", NEtaBins, h[0], h[1]);
  hEtaTrk[1][0][1] = new TH1D("hEtaAfterQA_piDet", "", NEtaBins, h[0], h[1]);
  hEtaTrk[1][1][1] = new TH1D("hEtaAfterQA_gaDet", "", NEtaBins, h[0], h[1]);
  hPtTrk[1][0][0]  = new TH1D("hPtBeforeQA_piDet", "", nPtHistBins, pTbinEdges);
  hPtTrk[1][1][0]  = new TH1D("hPtBeforeQA_gaDet", "", nPtHistBins, pTbinEdges);
  hPtTrk[1][0][1]  = new TH1D("hPtAfterQA_piDet", "", nPtHistBins, pTbinEdges);
  hPtTrk[1][1][1]  = new TH1D("hPtAfterQA_gaDet", "", nPtHistBins, pTbinEdges);
  // matching histograms
  hPhiForEff[0][0] = new TH1D("hPhiForEff_piPar", "", NPhiBins, f[0], f[1]);
  hPhiForEff[0][1] = new TH1D("hPhiForEff_gaPar", "", NPhiBins, f[0], f[1]);
  hPhiForEff[1][0] = new TH1D("hPhiForEff_piDet", "", NPhiBins, f[0], f[1]);
  hPhiForEff[1][1] = new TH1D("hPhiForEff_gaDet", "", NPhiBins, f[0], f[1]);
  hEtaForEff[0][0] = new TH1D("hEtaForEff_piPar", "", NEtaBins, h[0], h[1]);
  hEtaForEff[0][1] = new TH1D("hEtaForEff_gaPar", "", NEtaBins, h[0], h[1]);
  hEtaForEff[1][0] = new TH1D("hEtaForEff_piDet", "", NEtaBins, h[0], h[1]);
  hEtaForEff[1][1] = new TH1D("hEtaForEff_gaDet", "", NEtaBins, h[0], h[1]);
  hPtForEff[0][0]  = new TH1D("hPtForEff_piPar", "", nPtHistBins, pTbinEdges);
  hPtForEff[0][1]  = new TH1D("hPtForEff_gaPar", "", nPtHistBins, pTbinEdges);
  hPtForEff[1][0]  = new TH1D("hPtForEff_piDet", "", nPtHistBins, pTbinEdges);
  hPtForEff[1][1]  = new TH1D("hPtForEff_gaDet", "", nPtHistBins, pTbinEdges);
  hPhiParVsDet[0]  = new TH2D("hPhiParVsDet_pi", "", NPhiBins, f[0], f[1], NPhiBins, f[0], f[1]);
  hPhiParVsDet[1]  = new TH2D("hPhiParVsDet_ga", "", NPhiBins, f[0], f[1], NPhiBins, f[0], f[1]);
  hEtaParVsDet[0]  = new TH2D("hEtaParVsDet_pi", "", NEtaBins, h[0], h[1], NEtaBins, h[0], h[1]);
  hEtaParVsDet[1]  = new TH2D("hEtaParVsDet_ga", "", NEtaBins, h[0], h[1], NEtaBins, h[0], h[1]);
  hPtParVsDet[0]   = new TH2D("hPtParVsDet_pi", "", nPtHistBins, pTbinEdges, nPtHistBins, pTbinEdges);
  hPtParVsDet[1]   = new TH2D("hPtParVsDet_ga", "", nPtHistBins, pTbinEdges, nPtHistBins, pTbinEdges);
  // create efficiency histograms
  hPhiEff[0] = new TH1D("hPhiEfficiency_pi", "", NPhiBins, f[0], f[1]);
  hPhiEff[1] = new TH1D("hPhiEfficiency_ga", "", NPhiBins, f[0], f[1]);
  hEtaEff[0] = new TH1D("hEtaEfficiency_pi", "", NEtaBins, h[0], h[1]);
  hEtaEff[1] = new TH1D("hEtaEfficiency_ga", "", NEtaBins, h[0], h[1]);
  hPtEff[0]  = new TH1D("hPtEfficiency_pi", "", nPtHistBins, pTbinEdges);
  hPtEff[1]  = new TH1D("hPtEfficiency_ga", "", nPtHistBins, pTbinEdges);
  // create resolution histograms
  hPhiRes[0]      = new TH1D("hPhiRes_pi", "", nDf, dF[0], dF[1]);
  hPhiRes[1]      = new TH1D("hPhiRes_ga", "", nDf, dF[0], dF[1]);
  hEtaRes[0]      = new TH1D("hEtaRes_pi", "", nDh, dH[0], dH[1]);
  hEtaRes[1]      = new TH1D("hEtaRes_ga", "", nDh, dH[0], dH[1]);
  hPtRes[0]       = new TH1D("hPtRes_pi", "", nDpt, dPt[0], dPt[1]);
  hPtRes[1]       = new TH1D("hPtRes_ga", "", nDpt, dPt[0], dPt[1]);
  hPhiResVsPhi[0] = new TH2D("hPhiResVsPhi_pi", "", NPhiBins, f[0], f[1], nDf, dF[0], dF[1]);
  hPhiResVsPhi[1] = new TH2D("hPhiResVsPhi_ga", "", NPhiBins, f[0], f[1], nDf, dF[0], dF[1]);
  hEtaResVsEta[0] = new TH2D("hEtaResVsEta_pi", "", NEtaBins, h[0], h[1], nDh, dH[0], dH[1]);
  hEtaResVsEta[1] = new TH2D("hEtaResVsEta_ga", "", NEtaBins, h[0], h[1], nDh, dH[0], dH[1]);
  hPtResVsPt[0]   = new TH2D("hPtResVsPt_pi", "", nPt2d, pT2d[0], pT2d[1], nDpt, dPt[0], dPt[1]);
  hPtResVsPt[1]   = new TH2D("hPtResVsPt_ga", "", nPt2d, pT2d[0], pT2d[1], nDpt, dPt[0], dPt[1]);
  // errors
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    hPhiEff[iTrg]       -> Sumw2();
    hEtaEff[iTrg]       -> Sumw2();
    hPtEff[iTrg]        -> Sumw2();
    hPhiRes[iTrg]      -> Sumw2();
    hEtaRes[iTrg]      -> Sumw2();
    hPtRes[iTrg]       -> Sumw2();
    hPhiParVsDet[iTrg]  -> Sumw2();
    hEtaParVsDet[iTrg]  -> Sumw2();
    hPtParVsDet[iTrg]   -> Sumw2();
    hPhiResVsPhi[iTrg] -> Sumw2();
    hEtaResVsEta[iTrg] -> Sumw2();
    hPtResVsPt[iTrg]   -> Sumw2();
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      hPhiForEff[iLevel][iTrg] -> Sumw2();
      hEtaForEff[iLevel][iTrg] -> Sumw2();
      hPtForEff[iLevel][iTrg]  -> Sumw2();
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        hPhiTrk[iLevel][iTrg][iCut] -> Sumw2();
        hEtaTrk[iLevel][iTrg][iCut] -> Sumw2();
        hPtTrk[iLevel][iTrg][iCut]  -> Sumw2();
      }  // end cut loop
    }  // end level loop
  }  // end trigger loop
  cout << "    Made histograms." << endl;

  // make profiles
  TProfile *pPhiResVsPhi[NTrgs];
  TProfile *pEtaResVsEta[NTrgs];
  TProfile *pPtResVsPt[NTrgs];
  pPhiResVsPhi[0] = new TProfile("pPhiResVsPhi_pi", "", NPhiBins, f[0], f[1], "S");
  pPhiResVsPhi[1] = new TProfile("pPhiResVsPhi_ga", "", NPhiBins, f[0], f[1], "S");
  pEtaResVsEta[0] = new TProfile("pEtaResVsEta_pi", "", NEtaBins, h[0], h[1], "S");
  pEtaResVsEta[1] = new TProfile("pEtaResVsEta_ga", "", NEtaBins, h[0], h[1], "S");
  pPtResVsPt[0]   = new TProfile("pPtResVsPt_pi", "", nPt2d, pT2d[0], pT2d[1], "S");
  pPtResVsPt[1]   = new TProfile("pPtResVsPt_ga", "", nPt2d, pT2d[0], pT2d[1], "S");
  cout << "    Made profiles." << endl;

  // make tuples
  fOutput -> cd();
  TNtuple *nTrkAll   = new TNtuple("nTrkAll", "", "mcIdTruth:mcIdVtx:mcVtxEnd:mcIdGeant:mcEta:mcChrg:mcPt:trkIdTruth:trkQaTruth:trkNumFit:trkNumPoss:trkDca:trkEta:trkPt");
  TNtuple *nTrkMatch = new TNtuple("nTrkMatch", "", "mcIdTruth:mcIdVtx:mcVtxEnd:mcIdGeant:mcEta:mcChrg:mcPt:trkIdTruth:trkQaTruth:trkNumFit:trkNumPoss:trkDca:trkEta:trkPt");
  cout << "    Made NTuples." << endl;


  const UInt_t nParEvts = tPar -> GetEntriesFast();
  const UInt_t nDetEvts = tDet -> GetEntriesFast();
  cout << "    Beginning event loops: " << nDetEvts << " reconstructed and " << nParEvts << " generated events to process." << endl;

  // for selecting events and matched tracks
  map<UInt_t, Bool_t> isGoodEntry;
  map<UInt_t, UInt_t> triggerFlag;


  // no. of triggers
  UInt_t nTrgPar[NTrgs];
  UInt_t nTrgDet[NTrgs];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    nTrgPar[iTrg] = 0;
    nTrgDet[iTrg] = 0;
  }


  // detector event loop
  Int_t  detBytes(0);
  UInt_t nDetBytes(0);
  for (UInt_t iDetEvt = 0; iDetEvt < nDetEvts; iDetEvt++) {

    // load entry
    detBytes   = tDet -> GetEntry(iDetEvt);
    nDetBytes += detBytes;
    if (detBytes < 0) {
      cerr << "WARNING: issue with entry " << iDetEvt << "!" << endl;
      break;
    }
    else {
      if (isInBatchMode) {
        cout << "      Processing detector event " << iDetEvt + 1 << "/" << nDetEvts << "..." << endl;
      }
      else {
        cout << "      Processing detector event " << iDetEvt + 1 << "/" << nDetEvts << "...\r" << flush;
        if ((iDetEvt + 1) == nDetEvts) cout << endl;
      }
    }


    // event info
    const UInt_t   run   = (UInt_t) runNumber;
    const Long64_t nTrks = nPrimaryTracks;
    const Double_t xVtx  = xVertex;
    const Double_t yVtx  = yVertex;
    const Double_t zVtx  = zVertex;
    const Double_t rVtx  = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));

    Bool_t isGoodRun(true);
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (run == badRunList[iRun]) {
        isGoodRun = false;
        break;
      }
    }

    // event cuts
    const Bool_t isInRcut = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut = (TMath::Abs(zVtx) < zVtxMax);

    Bool_t isGoodEvt(true);
    if (!isGoodRun || !isInRcut || !isInZcut)
      isGoodEvt = false;


    // trigger info
    const Int_t    adc   = ETwradc11;
    const UInt_t   idTrg = ETwrdidT;
    const Double_t tsp   = Etsp;
    const Double_t eH4   = EEstrpen4;
    const Double_t eF4   = EPstripenp4;
    const Double_t pProj = ETwrPTower;
    const Double_t hDet  = ETwreT;
    const Double_t hPhys = EClustetav1;
    const Double_t eTrg  = EClustEneT0;
    const Double_t tTrg  = 2. * TMath::ATan(TMath::Exp(-1. * hPhys));
    const Double_t eTtrg = eTrg * TMath::Sin(tTrg);

    Bool_t isGoodTwr(true);
    for (UInt_t iTwr = 0; iTwr < NHotTwr; iTwr++) {
      if (idTrg == hotTwrList[iTwr]) {
        isGoodTwr = false;
        break;
      }
    }

    // trigger cuts
    UInt_t species(2);
    Bool_t isPi0(true);
    Bool_t isGam(true);
    Bool_t isGoodTrg(true);
    if (isTriggered) {
      const Bool_t isInAdcCut    = (adc <= adcMax);
      const Bool_t isInStrCut    = ((eH4 >= eStrMin) && (eF4 >= eStrMin));
      const Bool_t isInProjCut   = (pProj < pProjMax);
      const Bool_t isInEtaTrgCut = (TMath::Abs(hDet) < hTrgMax);
      const Bool_t isInEtCut     = ((eTtrg > eTtrgMin) && (eTtrg < eTtrgMax));
      const Bool_t isInPi0cut    = ((tsp > tspPi0[0]) && (tsp < tspPi0[1]));
      const Bool_t isInGamCut    = ((tsp > tspGam[0]) && (tsp < tspGam[1]));
      const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);
      if (!isGoodEvt || !isGoodTwr || !isInAdcCut || !isInStrCut || !isInProjCut || !isInEtaTrgCut || !isInEtCut || !isInTspCut) {
        isGoodTrg = false;
      }
      else {
        if (isInPi0cut) species = 0;
        if (isInGamCut) species = 1;
        isPi0     = isInPi0cut;
        isGam     = isInGamCut;
        isGoodTrg = true;
      }
    }  // end trigger check

    Bool_t eventCanBeUsed = (isGoodEvt && isGoodTrg);
    if (isGoodTrg) {
      if (isPi0 && eventCanBeUsed) nTrgDet[0]++;
      if (isGam && eventCanBeUsed) nTrgDet[1]++;
      isGoodEntry[iDetEvt] = eventCanBeUsed;
      triggerFlag[iDetEvt] = species;
    }
    else {
      isGoodEntry[iDetEvt] = eventCanBeUsed;
      triggerFlag[iDetEvt] = species;
      continue;
    }

  }  // end detector event loop


  // particle event loop
  Int_t  parBytes(0);
  UInt_t nParBytes(0);
  for (UInt_t iParEvt = 0; iParEvt < nParEvts; iParEvt++) {

    // load entry
    parBytes   = tPar -> GetEntry(iParEvt);
    nParBytes += parBytes;
    if (parBytes < 0) {
      cerr << "WARNING: issue with entry " << iParEvt << "!" << endl;
      break;
    }
    else {
      if (isInBatchMode) {
        cout << "      Processing particle event " << iParEvt + 1 << "/" << nParEvts << "..." << endl;
      }
      else {
        cout << "      Processing particle event " << iParEvt + 1 << "/" << nParEvts << "...\r" << flush;
        if ((iParEvt + 1) == nParEvts) cout << endl;
      }
    }


    // check event
    Bool_t isGoodEvt = isGoodEntry[iParEvt];
    if (!isGoodEvt) continue;

    // determine trigger
    Bool_t isPi0(true);
    Bool_t isGam(true);
    UInt_t trgFlag = triggerFlag[iParEvt];
    switch (trgFlag) {
      case 0:
        isPi0 = true;
        isGam = false;
        nTrgPar[0]++;
        break;
      case 1:
        isPi0 = false;
        isGam = true;
        nTrgPar[1]++;
        break;
      case 2:
        isPi0 = true;
        isGam = true;
        nTrgPar[0]++;
        nTrgPar[1]++;
        break;
    }


    // particle track loop
    const Int_t nMC = mcNumTrks;
    for (UInt_t iMC = 0; iMC < nMC; iMC++) {

      // track info
      const Int_t    mcId  = mcIdTrk   -> at(iMC);
      const Int_t    mcPID = mcIdGeant -> at(iMC);
      const Int_t    mcIV  = mcIdVx    -> at(iMC);
      const Int_t    mcIE  = mcIdVxEnd -> at(iMC);
      const Double_t qMC   = mcCharge  -> at(iMC);
      const Double_t fMC   = mcPhi     -> at(iMC);
      const Double_t hMC   = mcEta     -> at(iMC);
      const Double_t pTmc  = mcPt      -> at(iMC);


      // match to detector track
      const Bool_t isGoodEnd  = (mcIE == 0);
      const Bool_t isGoodEta  = (TMath::Abs(hMC) < 1.);
      const Bool_t isGoodChrg = (qMC != 0.);
      if (isGoodEnd && isGoodEta && isGoodChrg) {


        // fill particle histograms
        if (isPi0) {
          hPhiTrk[0][0][1] -> Fill(fMC);
          hPhiForEff[0][0] -> Fill(fMC);
          hEtaTrk[0][0][1] -> Fill(hMC);
          hEtaForEff[0][0] -> Fill(hMC);
          hPtTrk[0][0][1]  -> Fill(pTmc);
          hPtForEff[0][0]  -> Fill(pTmc);
        }
        if (isGam) {
          hPhiTrk[0][1][1] -> Fill(fMC);
          hPhiForEff[0][1] -> Fill(fMC);
          hEtaTrk[0][1][1] -> Fill(hMC);
          hEtaForEff[0][1] -> Fill(hMC);
          hPtTrk[0][1][1]  -> Fill(pTmc);
          hPtForEff[0][1]  -> Fill(pTmc);
        }

        // load detector entry
        Int_t  detBytes(0);
        UInt_t nDetBytes(0);
        detBytes   = tDet -> GetEntry(iParEvt);
        nDetBytes += detBytes;
        if (detBytes < 0) {
          cerr << "WARNING: issue with entry " << iParEvt << "!" << endl;
          break;
        }

        // detector track loop
        const Long64_t nTrks = nPrimaryTracks;
        for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

          // track info
          const Int_t    idTruth = (Int_t) PrimaryTrackArray_tofSigElectron[iTrk];
          const Int_t    qaTruth = (Int_t) PrimaryTrackArray_tofSigPion[iTrk];
          const UInt_t   nFit    = PrimaryTrackArray_nHitsFit[iTrk];
          const UInt_t   nPoss   = PrimaryTrackArray_nHitsPoss[iTrk];
          const Double_t rFit    = (Double_t) nFit / (Double_t) nPoss;
          const Double_t dca     = PrimaryTrackArray_dcag[iTrk];
          const Double_t chrg    = PrimaryTrackArray_charge[iTrk];
          const Double_t fTrk    = PrimaryTrackArray_phi[iTrk];
          const Double_t hTrk    = PrimaryTrackArray_eta[iTrk];
          const Double_t pTtrk   = PrimaryTrackArray_pT[iTrk];

          // fill detector histograms
          const Bool_t matchesMC   = (idTruth == mcId);
          const Bool_t isGoodMatch = (qaTruth > minQaTruth);
          if (matchesMC && isGoodMatch) {
            if (isPi0) {
              hPhiTrk[1][0][0] -> Fill(fTrk);
              hEtaTrk[1][0][0] -> Fill(hTrk);
              hPtTrk[1][0][0]  -> Fill(pTtrk);
            }
            if (isGam) {
              hPhiTrk[1][1][0] -> Fill(fTrk);
              hEtaTrk[1][1][0] -> Fill(hTrk);
              hPtTrk[1][1][0]  -> Fill(pTtrk);
            }
          }
          nTrkAll -> Fill(mcId, mcIV, mcIE, mcPID, hMC, qMC, pTmc, idTruth, qaTruth, nFit, nPoss, dca, hTrk, pTtrk);


          // track cuts
          const Bool_t isInFitCut    = (nFit > nFitMin);
          const Bool_t isInRatioCut  = (rFit > rFitMin);
          const Bool_t isInDcaCut    = (dca < dcaMax);
          const Bool_t isInChrgCut   = (chrg != 0.);
          const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
          const Bool_t isInPtCut     = ((pTtrk > pTtrkMin) && (pTtrk < pTtrkMax));
          if (!isInFitCut || !isInRatioCut || !isInDcaCut || !isInEtaTrkCut || !isInPtCut) continue;

          // TEST [01.21.2019]
          const Bool_t isInPhiCut = ((fTrk > phiMin) && (fTrk < phiMax));
          if (!isInPhiCut) continue;


          if (matchesMC && isGoodMatch) {
            const Double_t fDiff  = fMC - fTrk;
            const Double_t hDiff  = hMC - hTrk;
            const Double_t pTdiff = pTmc - pTtrk;
            if (isPi0) {
              hPhiTrk[1][0][1] -> Fill(fTrk);
              hPhiForEff[1][0] -> Fill(fMC);
              hPhiRes[0]       -> Fill(fDiff);
              hPhiResVsPhi[0]  -> Fill(fMC, fDiff);
              pPhiResVsPhi[0]  -> Fill(fMC, fDiff);
              hEtaTrk[1][0][1] -> Fill(hTrk);
              hEtaForEff[1][0] -> Fill(hMC);
              hEtaRes[0]       -> Fill(hDiff);
              hEtaResVsEta[0]  -> Fill(hMC, hDiff);
              pEtaResVsEta[0]  -> Fill(hMC, hDiff);
              hPtTrk[1][0][1]  -> Fill(pTtrk);
              hPtForEff[1][0]  -> Fill(pTmc);
              hPtRes[0]        -> Fill(pTdiff);
              hPtResVsPt[0]    -> Fill(pTmc, pTdiff);
              pPtResVsPt[0]    -> Fill(pTmc, pTdiff);
            }
            if (isGam) {
              hPhiTrk[1][1][1] -> Fill(fTrk);
              hPhiForEff[1][1] -> Fill(fMC);
              hPhiRes[1]       -> Fill(fDiff);
              hPhiResVsPhi[1]  -> Fill(fMC, fDiff);
              pPhiResVsPhi[1]  -> Fill(fMC, fDiff);
              hEtaTrk[1][1][1] -> Fill(hTrk);
              hEtaForEff[1][1] -> Fill(hMC);
              hEtaRes[1]       -> Fill(hDiff);
              hEtaResVsEta[1]  -> Fill(hMC, hDiff);
              pEtaResVsEta[1]  -> Fill(hMC, hDiff);
              hPtTrk[1][1][1]  -> Fill(pTtrk);
              hPtForEff[1][1]  -> Fill(pTmc);
              hPtRes[1]        -> Fill(pTdiff);
              hPtResVsPt[1]    -> Fill(pTmc, pTdiff);
              pPtResVsPt[1]    -> Fill(pTmc, pTdiff);
            }
            nTrkMatch -> Fill(mcId, mcIV, mcIE, mcPID, hMC, qMC, pTmc, idTruth, qaTruth, nFit, nPoss, dca, hTrk, pTtrk);
          }  // end matching condition
        }  // end detector track loop
      }  // end matching algorithm

      // fill particle histograms
      if (isPi0) {
        hPhiTrk[0][0][0] -> Fill(fMC);
        hEtaTrk[0][0][0] -> Fill(hMC);
        hPtTrk[0][0][0]  -> Fill(pTmc);
      }
      if (isGam) {
        hPhiTrk[0][1][0] -> Fill(fMC);
        hEtaTrk[0][1][0] -> Fill(hMC);
        hPtTrk[0][1][0]  -> Fill(pTmc);
      }

    }  // end particle track loop 
  }  // end particle event loop

  nPi0Trg = nTrgPar[0];
  nGamTrg = nTrgPar[1];
  cout << "    Event loops finished\n"
       << "      " << nTrgDet[0] << " (reconstructed), " << nTrgPar[0] << " (generated)  pi0 triggers,\n"
       << "      " << nTrgDet[1] << " (reconstructed), " << nTrgPar[1] << " (generated)  gamma triggers."
       << endl;


  // calculate ratios
  const Float_t weight(1.);
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    hPhiEff[iTrg] -> Divide(hPhiForEff[1][iTrg], hPhiForEff[0][iTrg], weight, weight);
    hEtaEff[iTrg] -> Divide(hEtaForEff[1][iTrg], hEtaForEff[0][iTrg], weight, weight);
    hPtEff[iTrg]  -> Divide(hPtForEff[1][iTrg], hPtForEff[0][iTrg], weight, weight);
  }
  cout << "    Ratios calculated." << endl;


  // normalize histograms
  const Double_t nTrg[NLevel][NTrgs] = {{(Double_t) nTrgPar[0], (Double_t) nTrgPar[1]}, {(Double_t) nTrgDet[0], (Double_t) nTrgDet[1]}};
  for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        hPhiTrk[iLevel][iTrg][iCut] -> Scale(1. / nTrg[iLevel][iTrg]);
        hEtaTrk[iLevel][iTrg][iCut] -> Scale(1. / nTrg[iLevel][iTrg]);
        hPtTrk[iLevel][iTrg][iCut]  -> Scale(1. / nTrg[iLevel][iTrg]);
      }  // end cut loop
    }  // end trigger loop
  }  // end level loop
  cout << "    Normalized track distributions." << endl;


  // fit efficiencies
  const UInt_t   fLinF(1);
  const UInt_t   fSizF(2);
  const UInt_t   fColF[NTrgs] = {859, 899};
  const TString  sTrg[NTrgs]  = {"Pi", "Ga"};
  const Double_t fEffGuess(0.5);
  const Double_t hEffGuess(0.5);
  const Double_t pEffGuess(0.87);
  const Double_t pSigGuess(4.0);
  const Double_t fFitRange[2] = {f[0], f[1]};
  const Double_t hFitRange[2] = {-0.7, 0.7};
  const Double_t pFitRange[2] = {1., 7.};

  TF1      *fPhiEff[NTrgs];
  TF1      *fEtaEff[NTrgs];
  TF1      *fPtEff[NTrgs];
  Double_t phiEff[NTrgs][NVal];
  Double_t etaEff[NTrgs][NVal];
  Double_t ptEff[NTrgs][NVal];
  Double_t ptSigE[NTrgs][NVal];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {

    // make names
    TString sPhiE("fPhiEff");
    TString sEtaE("fEtaEff");
    TString sPtE("fPtEff");
    sPhiE += sTrg[iTrg];
    sEtaE += sTrg[iTrg];
    sPtE  += sTrg[iTrg];

    // define fits
    fPhiEff[iTrg] = new TF1(sPhiE.Data(), "[0]", f[0], f[1]);
    fEtaEff[iTrg] = new TF1(sEtaE.Data(), "[0]", h[0], h[1]);
    fPtEff[iTrg]  = new TF1(sPtE.Data(), "[0]*(1-exp(-1.*[1]*x))", pTbinEdges[0], pTbinEdges[nPtHistBins]);
    fPhiEff[iTrg] -> SetParameter(0, fEffGuess);
    fPhiEff[iTrg] -> SetLineColor(fColF[iTrg]);
    fPhiEff[iTrg] -> SetLineStyle(fLinF);
    fPhiEff[iTrg] -> SetLineWidth(fSizF);
    fEtaEff[iTrg] -> SetParameter(0, hEffGuess);
    fEtaEff[iTrg] -> SetLineColor(fColF[iTrg]);
    fEtaEff[iTrg] -> SetLineStyle(fLinF);
    fEtaEff[iTrg] -> SetLineWidth(fSizF);
    fPtEff[iTrg]  -> SetParameters(0, pEffGuess);
    fPtEff[iTrg]  -> SetParameters(1, pSigGuess);
    fPtEff[iTrg]  -> SetLineColor(fColF[iTrg]);
    fPtEff[iTrg]  -> SetLineStyle(fLinF);
    fPtEff[iTrg]  -> SetLineWidth(fSizF);

    // fit histograms
    hPhiEff[iTrg] -> Fit(fPhiEff[iTrg], "BMQ0", "", fFitRange[0], fFitRange[1]);
    hEtaEff[iTrg] -> Fit(fEtaEff[iTrg], "BMQ0", "", hFitRange[0], hFitRange[1]);
    hPtEff[iTrg]  -> Fit(fPtEff[iTrg], "BMQ0", "", pFitRange[0], pFitRange[1]);
    phiEff[iTrg][0] = fPhiEff[iTrg] -> GetParameter(0);
    phiEff[iTrg][1] = fPhiEff[iTrg] -> GetParError(0);
    etaEff[iTrg][0] = fEtaEff[iTrg] -> GetParameter(0);
    etaEff[iTrg][1] = fEtaEff[iTrg] -> GetParError(0);
    ptEff[iTrg][0]  = fPtEff[iTrg]  -> GetParameter(0);
    ptEff[iTrg][1]  = fPtEff[iTrg]  -> GetParError(0);
    ptSigE[iTrg][0] = fPtEff[iTrg]  -> GetParameter(1);
    ptSigE[iTrg][1] = fPtEff[iTrg]  -> GetParError(1);

    // reset visibility
    const Int_t kNotDraw = 1<<9;
    if (hPhiEff[iTrg] -> GetFunction(sPhiE.Data()))
      hPhiEff[iTrg] -> GetFunction(sPhiE.Data()) -> ResetBit(kNotDraw);
    if (hEtaEff[iTrg] -> GetFunction(sEtaE.Data()))
      hEtaEff[iTrg] -> GetFunction(sEtaE.Data()) -> ResetBit(kNotDraw);
    if (hPtEff[iTrg] -> GetFunction(sPtE.Data()))
      hPtEff[iTrg] -> GetFunction(sPtE.Data()) -> ResetBit(kNotDraw);

  }
  cout << "    Fit efficiencies." << endl;


  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColE(1);
  const UInt_t  fMarE(1);
  const UInt_t  fColD[NCut]  = {1, 879};
  const UInt_t  fColM[NTrgs] = {859, 899};
  const UInt_t  fMarD[NCut]  = {7, 4};
  const UInt_t  fMarM[NTrgs] = {7, 4};
  const Float_t fLab(0.03);
  const Float_t fOffsetX(1.);
  const Float_t fPhiRangeX[2]   = {f[0], f[1]};
  const Float_t fPhiRangeY[2]   = {0., 0.47};
  const Float_t fEtaRangeX[2]   = {h[0], h[1]};
  const Float_t fEtaRangeY[2]   = {0., 0.73};
  const Float_t fPtRangeX[2]    = {-0.2, 20.2};
  const Float_t fPtRangeY[2]    = {0.00003, 13.};
  const Float_t fPtEffRangeX[2] = {-0.2, 7.2};
  const Float_t fEffRangeY[2]   = {0., 1.13};
  const TString sTitleXF("#varphi^{trk}");
  const TString sTitleXH("#eta^{trk}");
  const TString sTitleXP("p_{T}^{trk}");
  const TString sTitleYM("counts");
  const TString sTitleYF("(1/N^{trg}) dN^{trk}/d#varphi^{trk}");
  const TString sTitleYFE("#epsilon(#varphi^{trk}) = #epsilon_{#varphi}");
  const TString sTitleYH("(1/N^{trg}) dN^{trk}/d#eta^{trk}");
  const TString sTitleYHE("#epsilon(#eta^{trk}) = #epsilon_{#eta}");
  const TString sTitleYP("(1/N^{trg}) dN^{trk}/dp_{T}^{trk}");
  const TString sTitleYPE("#epsilon(p_{T}^{trk}) = #epsilon_{p} (1 - exp(-#sigma #upoint p_{T}^{trk}))");
  const TString sTitleF[NTrgs] = {"Track #varphi: #pi^{0} trigger", "Track #varphi: #gamma^{rich}"};
  const TString sTitleH[NTrgs] = {"Track #eta: #pi^{0} trigger", "Track #eta: #gamma^{rich}"};
  const TString sTitleP[NTrgs] = {"Track p_{T}: #pi^{0}", "Track p_{T}: #gamma^{rich}"};
  const TString sTitleE[NTrgs] = {"Track efficiency, #epsilon: #pi^{0} trigger", "Track efficiency, #epsilon: #gamma^{rich}"};
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    // phi efficiency histogram
    hPhiEff[iTrg] -> SetLineColor(fColE);
    hPhiEff[iTrg] -> SetMarkerColor(fColE);
    hPhiEff[iTrg] -> SetMarkerStyle(fMarE);
    hPhiEff[iTrg] -> SetTitle(sTitleE[iTrg].Data());
    hPhiEff[iTrg] -> SetTitleFont(fTxt);
    hPhiEff[iTrg] -> GetXaxis() -> SetTitle(sTitleXF.Data());
    hPhiEff[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hPhiEff[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hPhiEff[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hPhiEff[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hPhiEff[iTrg] -> GetXaxis() -> SetRangeUser(fPhiRangeX[0], fPhiRangeX[1]);
    hPhiEff[iTrg] -> GetYaxis() -> SetTitle(sTitleYFE.Data());
    hPhiEff[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hPhiEff[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hPhiEff[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hPhiEff[iTrg] -> GetYaxis() -> SetRangeUser(fEffRangeY[0], fEffRangeY[1]);
    // eta efficiency histogram
    hEtaEff[iTrg] -> SetLineColor(fColE);
    hEtaEff[iTrg] -> SetMarkerColor(fColE);
    hEtaEff[iTrg] -> SetMarkerStyle(fMarE);
    hEtaEff[iTrg] -> SetTitle(sTitleE[iTrg].Data());
    hEtaEff[iTrg] -> SetTitleFont(fTxt);
    hEtaEff[iTrg] -> GetXaxis() -> SetTitle(sTitleXH.Data());
    hEtaEff[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hEtaEff[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hEtaEff[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hEtaEff[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hEtaEff[iTrg] -> GetXaxis() -> SetRangeUser(fEtaRangeX[0], fEtaRangeX[1]);
    hEtaEff[iTrg] -> GetYaxis() -> SetTitle(sTitleYHE.Data());
    hEtaEff[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hEtaEff[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hEtaEff[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hEtaEff[iTrg] -> GetYaxis() -> SetRangeUser(fEffRangeY[0], fEffRangeY[1]);
    // pT efficiency histogram
    hPtEff[iTrg] -> SetLineColor(fColE);
    hPtEff[iTrg] -> SetMarkerColor(fColE);
    hPtEff[iTrg] -> SetMarkerStyle(fMarE);
    hPtEff[iTrg] -> SetTitle(sTitleE[iTrg].Data());
    hPtEff[iTrg] -> SetTitleFont(fTxt);
    hPtEff[iTrg] -> GetXaxis() -> SetTitle(sTitleXP.Data());
    hPtEff[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtEff[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hPtEff[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hPtEff[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hPtEff[iTrg] -> GetXaxis() -> SetRangeUser(fPtEffRangeX[0], fPtEffRangeX[1]);
    hPtEff[iTrg] -> GetYaxis() -> SetTitle(sTitleYPE.Data());
    hPtEff[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtEff[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hPtEff[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hPtEff[iTrg] -> GetYaxis() -> SetRangeUser(fEffRangeY[0], fEffRangeY[1]);
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      // phi histograms
      hPhiForEff[iLevel][iTrg] -> SetLineColor(fColM[iTrg]);
      hPhiForEff[iLevel][iTrg] -> SetMarkerColor(fColM[iTrg]);
      hPhiForEff[iLevel][iTrg] -> SetMarkerStyle(fMarM[iTrg]);
      hPhiForEff[iLevel][iTrg] -> SetTitle(sTitleF[iTrg].Data());
      hPhiForEff[iLevel][iTrg] -> SetTitleFont(fTxt);
      hPhiForEff[iLevel][iTrg] -> GetXaxis() -> SetTitle(sTitleXF.Data());
      hPhiForEff[iLevel][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hPhiForEff[iLevel][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPhiForEff[iLevel][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hPhiForEff[iLevel][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hPhiForEff[iLevel][iTrg] -> GetXaxis() -> SetRangeUser(fPhiRangeX[0], fPhiRangeX[1]);
      hPhiForEff[iLevel][iTrg] -> GetYaxis() -> SetTitle(sTitleYM.Data());
      hPhiForEff[iLevel][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hPhiForEff[iLevel][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hPhiForEff[iLevel][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      // eta histograms
      hEtaForEff[iLevel][iTrg] -> SetLineColor(fColM[iTrg]);
      hEtaForEff[iLevel][iTrg] -> SetMarkerColor(fColM[iTrg]);
      hEtaForEff[iLevel][iTrg] -> SetMarkerStyle(fMarM[iTrg]);
      hEtaForEff[iLevel][iTrg] -> SetTitle(sTitleH[iTrg].Data());
      hEtaForEff[iLevel][iTrg] -> SetTitleFont(fTxt);
      hEtaForEff[iLevel][iTrg] -> GetXaxis() -> SetTitle(sTitleXH.Data());
      hEtaForEff[iLevel][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hEtaForEff[iLevel][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hEtaForEff[iLevel][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hEtaForEff[iLevel][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hEtaForEff[iLevel][iTrg] -> GetXaxis() -> SetRangeUser(fEtaRangeX[0], fEtaRangeX[1]);
      hEtaForEff[iLevel][iTrg] -> GetYaxis() -> SetTitle(sTitleYM.Data());
      hEtaForEff[iLevel][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hEtaForEff[iLevel][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hEtaForEff[iLevel][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      // pT histograms
      hPtForEff[iLevel][iTrg] -> SetLineColor(fColM[iTrg]);
      hPtForEff[iLevel][iTrg] -> SetMarkerColor(fColM[iTrg]);
      hPtForEff[iLevel][iTrg] -> SetMarkerStyle(fMarM[iTrg]);
      hPtForEff[iLevel][iTrg] -> SetTitle(sTitleP[iTrg].Data());
      hPtForEff[iLevel][iTrg] -> SetTitleFont(fTxt);
      hPtForEff[iLevel][iTrg] -> GetXaxis() -> SetTitle(sTitleXP.Data());
      hPtForEff[iLevel][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hPtForEff[iLevel][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hPtForEff[iLevel][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hPtForEff[iLevel][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hPtForEff[iLevel][iTrg] -> GetXaxis() -> SetRangeUser(fPtRangeX[0], fPtRangeX[1]);
      hPtForEff[iLevel][iTrg] -> GetYaxis() -> SetTitle(sTitleYM.Data());
      hPtForEff[iLevel][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hPtForEff[iLevel][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hPtForEff[iLevel][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        // phi histograms
        hPhiTrk[iLevel][iTrg][iCut] -> SetLineColor(fColD[iCut]);
        hPhiTrk[iLevel][iTrg][iCut] -> SetMarkerColor(fColD[iCut]);
        hPhiTrk[iLevel][iTrg][iCut] -> SetMarkerStyle(fMarD[iCut]);
        hPhiTrk[iLevel][iTrg][iCut] -> SetTitle(sTitleF[iTrg].Data());
        hPhiTrk[iLevel][iTrg][iCut] -> SetTitleFont(fTxt);
        hPhiTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitle(sTitleXF.Data());
        hPhiTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
        hPhiTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hPhiTrk[iLevel][iTrg][iCut] -> GetXaxis() -> CenterTitle(fCnt);
        hPhiTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetLabelSize(fLab);
        hPhiTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetRangeUser(fPhiRangeX[0], fPhiRangeX[1]);
        hPhiTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetTitle(sTitleYF.Data());
        hPhiTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
        hPhiTrk[iLevel][iTrg][iCut] -> GetYaxis() -> CenterTitle(fCnt);
        hPhiTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetLabelSize(fLab);
        hPhiTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetRangeUser(fPhiRangeY[0], fPhiRangeY[1]);
        // eta histograms
        hEtaTrk[iLevel][iTrg][iCut] -> SetLineColor(fColD[iCut]);
        hEtaTrk[iLevel][iTrg][iCut] -> SetMarkerColor(fColD[iCut]);
        hEtaTrk[iLevel][iTrg][iCut] -> SetMarkerStyle(fMarD[iCut]);
        hEtaTrk[iLevel][iTrg][iCut] -> SetTitle(sTitleH[iTrg].Data());
        hEtaTrk[iLevel][iTrg][iCut] -> SetTitleFont(fTxt);
        hEtaTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitle(sTitleXH.Data());
        hEtaTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
        hEtaTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hEtaTrk[iLevel][iTrg][iCut] -> GetXaxis() -> CenterTitle(fCnt);
        hEtaTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetLabelSize(fLab);
        hEtaTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetRangeUser(fEtaRangeX[0], fEtaRangeX[1]);
        hEtaTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetTitle(sTitleYH.Data());
        hEtaTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
        hEtaTrk[iLevel][iTrg][iCut] -> GetYaxis() -> CenterTitle(fCnt);
        hEtaTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetLabelSize(fLab);
        hEtaTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetRangeUser(fEtaRangeY[0], fEtaRangeY[1]);
        // pT histograms
        hPtTrk[iLevel][iTrg][iCut] -> SetLineColor(fColD[iCut]);
        hPtTrk[iLevel][iTrg][iCut] -> SetMarkerColor(fColD[iCut]);
        hPtTrk[iLevel][iTrg][iCut] -> SetMarkerStyle(fMarD[iCut]);
        hPtTrk[iLevel][iTrg][iCut] -> SetTitle(sTitleP[iTrg].Data());
        hPtTrk[iLevel][iTrg][iCut] -> SetTitleFont(fTxt);
        hPtTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitle(sTitleXP.Data());
        hPtTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitleFont(fTxt);
        hPtTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetTitleOffset(fOffsetX);
        hPtTrk[iLevel][iTrg][iCut] -> GetXaxis() -> CenterTitle(fCnt);
        hPtTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetLabelSize(fLab);
        hPtTrk[iLevel][iTrg][iCut] -> GetXaxis() -> SetRangeUser(fPtRangeX[0], fPtRangeX[1]);
        hPtTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetTitle(sTitleYP.Data());
        hPtTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetTitleFont(fTxt);
        hPtTrk[iLevel][iTrg][iCut] -> GetYaxis() -> CenterTitle(fCnt);
        hPtTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetLabelSize(fLab);
        hPtTrk[iLevel][iTrg][iCut] -> GetYaxis() -> SetRangeUser(fPtRangeY[0], fPtRangeY[1]);
      }  // end cut loop
    }  // end level loop
  }  // end trigger loop
  cout << "    Set styles." << endl;


  // make labels
  const UInt_t  nDec(3);
  const UInt_t  fAlign(12);
  const UInt_t  fColP(0);
  const UInt_t  fColT[NTrgs] = {859, 899};
  const Float_t xPav[2]      = {0.7, 0.9};
  const Float_t xLeg[2]      = {0.7, 0.9};
  const Float_t yPav[2]      = {0.7, 0.9};
  const Float_t yLeg[2]      = {0.5, 0.7};
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrgKin("E_{T}^{trg} #in (9, 20) GeV, |#eta^{trg}| < 0.9");
  const TString sPhiEff("#epsilon_{#varphi} = ");
  const TString sEtaEff("#epsilon_{#eta} = ");
  const TString sPtEff("#epsilon_{p} = ");
  const TString sPtSigE("#sigma_{p} = ");
  const TString sLegEff[2] = {"particle", "detector"};

  TLegend   *lPhiE[NTrgs];
  TLegend   *lEtaE[NTrgs];
  TLegend   *lPtE[NTrgs];
  TPaveText *pPhiE[NTrgs];
  TPaveText *pEtaE[NTrgs];
  TPaveText *pPtE[NTrgs];
  TString   sFrawE[NVal];
  TString   sHrawE[NVal];
  TString   sPErawE[NVal];
  TString   sPSrawE[NVal];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {

    TString sFtxtE(sPhiEff.Data());
    TString sHtxtE(sEtaEff.Data());
    TString sPEtxtE(sPtEff.Data());
    TString sPStxtE(sPtSigE.Data());
    for (UInt_t iVal = 0; iVal < NVal; iVal++) {
      sFrawE[iVal]   = "";
      sHrawE[iVal]   = "";
      sPErawE[iVal]  = "";
      sPSrawE[iVal]  = "";
      sFrawE[iVal]  += phiEff[iTrg][iVal];
      sHrawE[iVal]  += etaEff[iTrg][iVal];
      sPErawE[iVal] += ptEff[iTrg][iVal];
      sPSrawE[iVal] += ptSigE[iTrg][iVal];

      const UInt_t nFrawE  = sFrawE[iVal].First(".");
      const UInt_t nHrawE  = sHrawE[iVal].First(".");
      const UInt_t nPErawE = sPErawE[iVal].First(".");
      const UInt_t nPSrawE = sPSrawE[iVal].First(".");
      const UInt_t nFtxtE  = (nFrawE + nDec) + 1;
      const UInt_t nHtxtE  = (nHrawE + nDec) + 1;
      const UInt_t nPEtxtE = (nPErawE + nDec) + 1;
      const UInt_t nPStxtE = (nPSrawE + nDec) + 1;
      if (iVal == 0) {
        sFtxtE.Append(sFrawE[iVal].Data(), nFtxtE);
        sHtxtE.Append(sHrawE[iVal].Data(), nHtxtE);
        sPEtxtE.Append(sPErawE[iVal].Data(), nPEtxtE);
        sPStxtE.Append(sPSrawE[iVal].Data(), nPStxtE);
      } 
      else {
        sFtxtE  += " #pm ";
        sHtxtE  += " #pm ";
        sPEtxtE += " #pm ";
        sPStxtE += " #pm ";
        sFtxtE.Append(sFrawE[iVal].Data(), nFtxtE);
        sHtxtE.Append(sHrawE[iVal].Data(), nHtxtE);
        sPEtxtE.Append(sPErawE[iVal].Data(), nPEtxtE);
        sPStxtE.Append(sPSrawE[iVal].Data(), nPStxtE);
      }
    }  // end value loop

    // phi labels
    pPhiE[iTrg] = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pPhiE[iTrg] -> SetFillColor(fColP);
    pPhiE[iTrg] -> SetLineColor(fColP);
    pPhiE[iTrg] -> SetTextColor(fColT[iTrg]);
    pPhiE[iTrg] -> SetTextFont(fTxt);
    pPhiE[iTrg] -> SetTextAlign(fAlign);
    pPhiE[iTrg] -> AddText(sSystem.Data());
    pPhiE[iTrg] -> AddText(sTrgKin.Data());
    pPhiE[iTrg] -> AddText(sFtxtE.Data());
    // eta labels
    pEtaE[iTrg] = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pEtaE[iTrg] -> SetFillColor(fColP);
    pEtaE[iTrg] -> SetLineColor(fColP);
    pEtaE[iTrg] -> SetTextColor(fColT[iTrg]);
    pEtaE[iTrg] -> SetTextFont(fTxt);
    pEtaE[iTrg] -> SetTextAlign(fAlign);
    pEtaE[iTrg] -> AddText(sSystem.Data());
    pEtaE[iTrg] -> AddText(sTrgKin.Data());
    pEtaE[iTrg] -> AddText(sHtxtE.Data());
    // pT labels
    pPtE[iTrg]  = new TPaveText(xPav[0], yPav[0], xPav[1], xPav[1], "NDC NB");
    pPtE[iTrg] -> SetFillColor(fColP);
    pPtE[iTrg] -> SetLineColor(fColP);
    pPtE[iTrg] -> SetTextColor(fColT[iTrg]);
    pPtE[iTrg] -> SetTextFont(fTxt);
    pPtE[iTrg] -> SetTextAlign(fAlign);
    pPtE[iTrg] -> AddText(sSystem.Data());
    pPtE[iTrg] -> AddText(sTrgKin.Data());
    pPtE[iTrg] -> AddText(sPEtxtE.Data());
    pPtE[iTrg] -> AddText(sPStxtE.Data());

    // phi legend
    lPhiE[iTrg] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lPhiE[iTrg] -> SetFillColor(fColP);
    lPhiE[iTrg] -> SetLineColor(fColP);
    lPhiE[iTrg] -> SetTextFont(fTxt);
    lPhiE[iTrg] -> AddEntry(hPhiForEff[0][iTrg], sLegEff[0].Data());
    lPhiE[iTrg] -> AddEntry(hPhiForEff[1][iTrg], sLegEff[1].Data());
    // eta legend
    lEtaE[iTrg] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lEtaE[iTrg] -> SetFillColor(fColP);
    lEtaE[iTrg] -> SetLineColor(fColP);
    lEtaE[iTrg] -> SetTextFont(fTxt);
    lEtaE[iTrg] -> AddEntry(hEtaForEff[0][iTrg], sLegEff[0].Data());
    lEtaE[iTrg] -> AddEntry(hEtaForEff[1][iTrg], sLegEff[1].Data());
    // pT legend
    lPtE[iTrg] = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
    lPtE[iTrg] -> SetFillColor(fColP);
    lPtE[iTrg] -> SetLineColor(fColP);
    lPtE[iTrg] -> SetTextFont(fTxt);
    lPtE[iTrg] -> AddEntry(hPtForEff[0][iTrg], sLegEff[0].Data());
    lPtE[iTrg] -> AddEntry(hPtForEff[1][iTrg], sLegEff[1].Data());

  }  // end trigger loop
  cout << "    Made labels." << endl;



  // make directories and save histograms
  const TString sTrgDir[NTrgs]  = {"pi0", "gamma"};
  const TString sLevDir[NLevel] = {"particle", "detector"};

  TDirectory *dTrgs[NTrgs];
  TDirectory *dLvls[NTrgs][NLevel];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    dTrgs[iTrg] = (TDirectory*) fOutput -> mkdir(sTrgDir[iTrg].Data());
    dTrgs[iTrg] -> cd();
    for (UInt_t iLevel = 0; iLevel < NLevel; iLevel++) {
      dLvls[iTrg][iLevel] = (TDirectory*) dTrgs[iTrg] -> mkdir(sLevDir[iLevel].Data());
      dLvls[iTrg][iLevel] -> cd();
      for (UInt_t iCut = 0; iCut < NCut; iCut++) {
        hPhiTrk[iLevel][iTrg][iCut] -> Write();
        hEtaTrk[iLevel][iTrg][iCut] -> Write();
        hPtTrk[iLevel][iTrg][iCut]  -> Write();
      }  // end cut loop
      hPhiForEff[iLevel][iTrg] -> Write();
      hEtaForEff[iLevel][iTrg] -> Write();
      hPtForEff[iLevel][iTrg]  -> Write();
    }  // end level loop
    dTrgs[iTrg]         -> cd();
    hPhiEff[iTrg]       -> Write();
    hEtaEff[iTrg]       -> Write();
    hPtEff[iTrg]        -> Write();
    hPhiRes[iTrg]      -> Write();
    hEtaRes[iTrg]      -> Write();
    hPtRes[iTrg]       -> Write();
    hPhiParVsDet[iTrg]  -> Write();
    hEtaParVsDet[iTrg]  -> Write();
    hPtParVsDet[iTrg]   -> Write();
    hPhiResVsPhi[iTrg] -> Write();
    hEtaResVsEta[iTrg] -> Write();
    hPtResVsPt[iTrg]   -> Write();
    pPhiResVsPhi[iTrg] -> Write();
    pEtaResVsEta[iTrg] -> Write();
    pPtResVsPt[iTrg]   -> Write();
  }  // end trigger loop
  fOutput   -> cd();
  nTrkAll   -> Write();
  nTrkMatch -> Write();
  cout << "    Made directories." << endl;


  // make plots
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  ticks(1);
  const UInt_t  log(1);
  const Float_t margin(0.);
  const Float_t xPad[NTrgs + 2] = {0., 0.5, 0.5, 1.};
  const TString sPads[NTrgs]    = {"pPi0", "pGamma"};

  TPad    *pPhiTrk[NTrgs];
  TPad    *pPhiForEff[NTrgs];
  TPad    *pEtaTrk[NTrgs];
  TPad    *pEtaForEff[NTrgs];
  TPad    *pPtTrk[NTrgs];
  TPad    *pPtForEff[NTrgs];
  TPad    *pPhiEff[NTrgs];
  TPad    *pEtaEff[NTrgs];
  TPad    *pPtEff[NTrgs];
  TCanvas *cPhiTrk;
  TCanvas *cPhiForEff;
  TCanvas *cEtaTrk;
  TCanvas *cEtaForEff;
  TCanvas *cPtTrk;
  TCanvas *cPtForEff;
  TCanvas *cPhiEff;
  TCanvas *cEtaEff;
  TCanvas *cPtEff;
  fOutput -> cd();
  // phi plots
  cPhiForEff    = new TCanvas("cPhiForEff", "", width, height);
  pPhiForEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPhiForEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPhiForEff[0]    -> SetGrid(grid, grid);
  pPhiForEff[0]    -> SetTicks(ticks, ticks);
  pPhiForEff[0]    -> SetRightMargin(margin);
  pPhiForEff[1]    -> SetGrid(grid, grid);
  pPhiForEff[1]    -> SetTicks(ticks, ticks);
  pPhiForEff[1]    -> SetLeftMargin(margin);
  cPhiForEff       -> cd();
  pPhiForEff[0]    -> Draw();
  pPhiForEff[1]    -> Draw();
  pPhiForEff[0]    -> cd();
  hPhiForEff[0][0] -> Draw();
  hPhiForEff[1][0] -> Draw("same");
  lPhiE[0]         -> Draw();
  pPhiE[0]         -> Draw();
  pPhiForEff[1]    -> cd();
  hPhiForEff[0][1] -> Draw();
  hPhiForEff[1][1] -> Draw("same");
  lPhiE[1]         -> Draw();
  pPhiE[1]         -> Draw();
  cPhiForEff       -> Write();
  cPhiForEff       -> Close();

  // eta plots
  cEtaForEff    = new TCanvas("cEtaForEff", "", width, height);
  pEtaForEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pEtaForEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pEtaForEff[0]    -> SetGrid(grid, grid);
  pEtaForEff[0]    -> SetTicks(ticks, ticks);
  pEtaForEff[0]    -> SetRightMargin(margin);
  pEtaForEff[1]    -> SetGrid(grid, grid);
  pEtaForEff[1]    -> SetTicks(ticks, ticks);
  pEtaForEff[1]    -> SetLeftMargin(margin);
  cEtaForEff       -> cd();
  pEtaForEff[0]    -> Draw();
  pEtaForEff[1]    -> Draw();
  pEtaForEff[0]    -> cd();
  hEtaForEff[0][0] -> Draw();
  hEtaForEff[1][0] -> Draw("same");
  lEtaE[0]         -> Draw();
  pEtaE[0]         -> Draw();
  pEtaForEff[1]    -> cd();
  hEtaForEff[0][1] -> Draw();
  hEtaForEff[1][1] -> Draw("same");
  lEtaE[1]         -> Draw();
  pEtaE[1]         -> Draw();
  cEtaForEff       -> Write();
  cEtaForEff       -> Close();

  // pT plots
  cPtForEff    = new TCanvas("cPtForEff", "", width, height);
  pPtForEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPtForEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPtForEff[0]    -> SetGrid(grid, grid);
  pPtForEff[0]    -> SetTicks(ticks, ticks);
  pPtForEff[0]    -> SetRightMargin(margin);
  pPtForEff[0]    -> SetLogy(log);
  pPtForEff[1]    -> SetGrid(grid, grid);
  pPtForEff[1]    -> SetTicks(ticks, ticks);
  pPtForEff[1]    -> SetLeftMargin(margin);
  pPtForEff[1]    -> SetLogy(log);
  cPtForEff       -> cd();
  pPtForEff[0]    -> Draw();
  pPtForEff[1]    -> Draw();
  pPtForEff[0]    -> cd();
  hPtForEff[0][0] -> Draw();
  hPtForEff[1][0] -> Draw("same");
  lPtE[0]         -> Draw();
  pPtE[0]         -> Draw();
  pPtForEff[1]    -> cd();
  hPtForEff[0][1] -> Draw();
  hPtForEff[1][1] -> Draw("same");
  lPtE[1]         -> Draw();
  pPtE[1]         -> Draw();
  cPtForEff       -> Write();
  cPtForEff       -> Close();

  // phi efficiency plots
  cPhiEff    = new TCanvas("cPhiEff", "", width, height);
  pPhiEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPhiEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPhiEff[0] -> SetGrid(grid, grid);
  pPhiEff[0] -> SetTicks(ticks, ticks);
  pPhiEff[0] -> SetRightMargin(margin);
  pPhiEff[1] -> SetGrid(grid, grid);
  pPhiEff[1] -> SetTicks(ticks, ticks);
  pPhiEff[1] -> SetLeftMargin(margin);
  cPhiEff    -> cd();
  pPhiEff[0] -> Draw();
  pPhiEff[1] -> Draw();
  pPhiEff[0] -> cd();
  hPhiEff[0] -> Draw();
  pPhiE[0]   -> Draw();
  pPhiEff[1] -> cd();
  hPhiEff[1] -> Draw();
  pPhiE[1]   -> Draw();
  cPhiEff    -> Write();
  cPhiEff    -> Close();

  // eta efficiency plots
  cEtaEff    = new TCanvas("cEtaEff", "", width, height);
  pEtaEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pEtaEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pEtaEff[0] -> SetGrid(grid, grid);
  pEtaEff[0] -> SetTicks(ticks, ticks);
  pEtaEff[0] -> SetRightMargin(margin);
  pEtaEff[1] -> SetGrid(grid, grid);
  pEtaEff[1] -> SetTicks(ticks, ticks);
  pEtaEff[1] -> SetLeftMargin(margin);
  cEtaEff    -> cd();
  pEtaEff[0] -> Draw();
  pEtaEff[1] -> Draw();
  pEtaEff[0] -> cd();
  hEtaEff[0] -> Draw();
  pEtaE[0]   -> Draw();
  pEtaEff[1] -> cd();
  hEtaEff[1] -> Draw();
  pEtaE[1]   -> Draw();
  cEtaEff    -> Write();
  cEtaEff    -> Close();

  // pT efficiency plots
  cPtEff    = new TCanvas("cPtEff", "", width, height);
  pPtEff[0] = new TPad(sPads[0], "", xPad[0], 0., xPad[1], 1.);
  pPtEff[1] = new TPad(sPads[1], "", xPad[2], 0., xPad[3], 1.);
  pPtEff[0] -> SetGrid(grid, grid);
  pPtEff[0] -> SetTicks(ticks, ticks);
  pPtEff[0] -> SetRightMargin(margin);
  pPtEff[1] -> SetGrid(grid, grid);
  pPtEff[1] -> SetTicks(ticks, ticks);
  pPtEff[1] -> SetLeftMargin(margin);
  cPtEff    -> cd();
  pPtEff[0] -> Draw();
  pPtEff[1] -> Draw();
  pPtEff[0] -> cd();
  hPtEff[0] -> Draw();
  pPtE[0]   -> Draw();
  pPtEff[1] -> cd();
  hPtEff[1] -> Draw();
  pPtE[1]   -> Draw();
  cPtEff    -> Write();
  cPtEff    -> Close();
  cout << "    Drew plots." << endl;


  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();

  cout << "  Calculation finished!\n" << endl;
  return;

}

// End ------------------------------------------------------------------------
