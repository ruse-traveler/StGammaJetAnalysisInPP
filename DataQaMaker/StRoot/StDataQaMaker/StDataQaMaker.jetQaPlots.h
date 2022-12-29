// 'StDataQaMaker.jetQaPlots.h'
// Derek Anderson
// 09.19.2019
//
// This class produces the data QA
// plots for the neutral-triggered
// pp recoil jet analysis note.

#pragma once

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

using namespace std;
using namespace fastjet;



void StDataQaMaker::MakeJetQaPlots() {

  cout << "    Making jet QA plots..." << endl;

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
  const Double_t eTbinMin[NTrgBins] = {9., 9., 11., 15.};
  const Double_t eTbinMax[NTrgBins] = {20., 11., 15., 20.};
  const Double_t tspPi0[NTrgTsp]    = {0., 0.08};
  const Double_t tspGam[NTrgTsp]    = {0.2, 0.6};

  // track parameters
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(30.);
  const Double_t pTbinMin[NJetBins] = {0.2, 2., 5., 9.};
  const Double_t pTbinMax[NJetBins] = {2., 5., 9., 30.};

  // jet parameters
  const UInt_t   nRepeat(1);
  const UInt_t   nRemove(1);
  const Double_t aGhost(0.01);
  const Double_t hBkgdMax(1.0);
  const Double_t pTjetMin(0.2);
  const Double_t pTjetMax(30.);
  const Double_t rJet[NResVals]      = {0.2, 0.5};
  const Double_t aJetMin[NResVals]   = {0.05, 0.65};
  const Double_t hJetMax[NResVals]   = {0.8, 0.5};
  const Double_t hGhostMax[NResVals] = {1.2, 1.5};

  // misc parameters
  const Double_t pi(TMath::Pi());
  const Double_t tau(TMath::TwoPi());
  const Double_t piOver2(TMath::PiOver2());
  const Double_t piOver4(TMath::PiOver4());
  const Double_t dFjetMin(-1. * piOver2);
  const Double_t dFjetMax(3. * piOver2);
  const Double_t dFrecoilMax(piOver4);
  const TString  sResParams[NResVals] = {"_R02", "_R05"};
  const TString  sPtBins[NJetBins]    = {"pt022", "pt25", "pt59", "pt930"};
  const TString  sDir[NResVals]       = {"r02", "r05"};
  const TString  sResTxt[NResVals]    = {"anti-k_{T} algo., R = 0.2", "anti-k_{T} algo., R = 0.5"};

  // open files
  TFile *fOutput = new TFile(sOutFile.Data(), "recreate");
  TFile *fInput  = new TFile(sInFile.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "      Opened files." << endl;

  // grab input tree
  TTree *tInput;
  fInput -> GetObject(sInTree.Data(), tInput);
  if (!tInput) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    return;
  }
  cout << "      Grabbed tree." << endl;

  // declare input leaf addresses
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
  Int_t    FlagEvent_TrgTrkMisMtch;
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
  Double_t PrimaryTrackArray_nHitsFit[NTrkMax];
  Double_t PrimaryTrackArray_nHitsPoss[NTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
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

  // declare input branches
  TBranch *bEventList_fUniqueID;
  TBranch *bEventList_fBits;
  TBranch *bEventList_runNumber;
  TBranch *bEventList_eventNumber;
  TBranch *bEventList_trigID;
  TBranch *bEventList_nGlobalTracks;
  TBranch *bEventList_nPrimaryTracks;
  TBranch *bEventList_refMult;
  TBranch *bEventList_vpdVz;
  TBranch *bEventList_xVertex;
  TBranch *bEventList_yVertex;
  TBranch *bEventList_zVertex;
  TBranch *bEventList_bbcZVertex;
  TBranch *bEventList_zdcCoincidenceRate;
  TBranch *bEventList_bbcCoincidenceRate;
  TBranch *bEventList_backgroundRate;
  TBranch *bEventList_bbcBlueBackgroundRate;
  TBranch *bEventList_bbcYellowBackgroundRate;
  TBranch *bEventList_refMultPos;
  TBranch *bEventList_refMultNeg;
  TBranch *bEventList_bTOFTrayMultiplicity;
  TBranch *bEventList_nVerticies;
  TBranch *bEventList_MagF;
  TBranch *bEventList_VrtxRank;
  TBranch *bEventList_FlagEvent_TrgTrkMisMtch;
  TBranch *bEventList_Etsp;
  TBranch *bEventList_ETwrdidT;
  TBranch *bEventList_ETwradc11;
  TBranch *bEventList_ETwreneT0;
  TBranch *bEventList_ETwreT;
  TBranch *bEventList_ETwrENET0;
  TBranch *bEventList_ETwrphT;
  TBranch *bEventList_ETwrPTower;
  TBranch *bEventList_ETwrpidTower;
  TBranch *bEventList_ETwrmoduleT;
  TBranch *bEventList_EClustEneT0;
  TBranch *bEventList_EClustetav1;
  TBranch *bEventList_EClustphiv1;
  TBranch *bEventList_EEstrpen01;
  TBranch *bEventList_EEstrpen02;
  TBranch *bEventList_EEstrpen03;
  TBranch *bEventList_EEstrpen0;
  TBranch *bEventList_EEstrpen1;
  TBranch *bEventList_EEstrpen2;
  TBranch *bEventList_EEstrpen3;
  TBranch *bEventList_EEstrpen4;
  TBranch *bEventList_EEstrpen5;
  TBranch *bEventList_EEstrpen6;
  TBranch *bEventList_EEstrpen7;
  TBranch *bEventList_EEstrpen8;
  TBranch *bEventList_EEstrpen9;
  TBranch *bEventList_EEstrpen10;
  TBranch *bEventList_EEstrpen11;
  TBranch *bEventList_EEstrpen12;
  TBranch *bEventList_EEstrpen13;
  TBranch *bEventList_EEstrpen14;
  TBranch *bEventList_EEstrpen15;
  TBranch *bEventList_ETwrdidE;
  TBranch *bEventList_EPstripenp01;
  TBranch *bEventList_EPstripenp02;
  TBranch *bEventList_EPstripenp03;
  TBranch *bEventList_EPstripenp0;
  TBranch *bEventList_EPstripenp1;
  TBranch *bEventList_EPstripenp2;
  TBranch *bEventList_EPstripenp3;
  TBranch *bEventList_EPstripenp4;
  TBranch *bEventList_EPstripenp5;
  TBranch *bEventList_EPstripenp6;
  TBranch *bEventList_EPstripenp7;
  TBranch *bEventList_EPstripenp8;
  TBranch *bEventList_EPstripenp9;
  TBranch *bEventList_EPstripenp10;
  TBranch *bEventList_EPstripenp11;
  TBranch *bEventList_EPstripenp12;
  TBranch *bEventList_EPstripenp13;
  TBranch *bEventList_EPstripenp14;
  TBranch *bEventList_EPstripenp15;
  TBranch *bEventList_EclustEnnq1;
  TBranch *bEventList_EclustEnnq20;
  TBranch *bEventList_EclustEnnq19;
  TBranch *bEventList_EclustEnpq1;
  TBranch *bEventList_EclustEnpq20;
  TBranch *bEventList_EclustEnpq19;
  TBranch *bEventList_EclustEnpq21;
  TBranch *bEventList_PrimaryTrackArray_;
  TBranch *bPrimaryTrackArray_fUniqueID;
  TBranch *bPrimaryTrackArray_fBits;
  TBranch *bPrimaryTrackArray_nHitsFit;
  TBranch *bPrimaryTrackArray_nHitsPoss;
  TBranch *bPrimaryTrackArray_trackFlag;
  TBranch *bPrimaryTrackArray_pZ;
  TBranch *bPrimaryTrackArray_pX;
  TBranch *bPrimaryTrackArray_pY;
  TBranch *bPrimaryTrackArray_pT;
  TBranch *bPrimaryTrackArray_dEdx;
  TBranch *bPrimaryTrackArray_charge;
  TBranch *bPrimaryTrackArray_tofBeta;
  TBranch *bPrimaryTrackArray_eta;
  TBranch *bPrimaryTrackArray_phi;
  TBranch *bPrimaryTrackArray_nSigElectron;
  TBranch *bPrimaryTrackArray_nSigPion;
  TBranch *bPrimaryTrackArray_nSigKaon;
  TBranch *bPrimaryTrackArray_nSigProton;
  TBranch *bPrimaryTrackArray_dcag;
  TBranch *bPrimaryTrackArray_nHits;
  TBranch *bPrimaryTrackArray_dEdxHits;
  TBranch *bPrimaryTrackArray_firstZPoint;
  TBranch *bPrimaryTrackArray_lastZPoint;
  TBranch *bPrimaryTrackArray_tofSigElectron;
  TBranch *bPrimaryTrackArray_tofSigPion;
  TBranch *bPrimaryTrackArray_tofSigKaon;
  TBranch *bPrimaryTrackArray_tofSigProton;
  TBranch *bPrimaryTrackArray_timeOfflight;
  TBranch *bPrimaryTrackArray_pathLength;
  TBranch *bPrimaryTrackArray_trkIndex;
  TBranch *bEventList_TowerArray_;
  TBranch *bTowerArray_fUniqueID;
  TBranch *bTowerArray_fBits;
  TBranch *bTowerArray_TwrId;
  TBranch *bTowerArray_TwrEng;
  TBranch *bTowerArray_TwrEta;
  TBranch *bTowerArray_TwrPhi;
  TBranch *bTowerArray_TwrADC;
  TBranch *bTowerArray_TwrPed;
  TBranch *bTowerArray_TwrRMS;
  TBranch *bTowerArray_TwrMatchIdnex;
  TBranch *bTowerArray_NoOfmatchedTrk;
  TBranch *bTowerArray_TwrMatchP;
  TBranch *bTowerArray_TwrPx;
  TBranch *bTowerArray_TwrPy;
  TBranch *bTowerArray_TwrPz;
  TBranch *bTowerArray_fNAssocTracks;
  TBranch *bTowerArray_fMatchedTracksArray_;
  TBranch *bTowerArray_fMatchedTracksArray_P;

  // set input branches
  tInput -> SetMakeClass(1);
  tInput -> SetBranchAddress("fUniqueID", &fUniqueID, &bEventList_fUniqueID);
  tInput -> SetBranchAddress("fBits", &fBits, &bEventList_fBits);
  tInput -> SetBranchAddress("runNumber", &runNumber, &bEventList_runNumber);
  tInput -> SetBranchAddress("eventNumber", &eventNumber, &bEventList_eventNumber);
  tInput -> SetBranchAddress("trigID", &trigID, &bEventList_trigID);
  tInput -> SetBranchAddress("nGlobalTracks", &nGlobalTracks, &bEventList_nGlobalTracks);
  tInput -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks, &bEventList_nPrimaryTracks);
  tInput -> SetBranchAddress("refMult", &refMult, &bEventList_refMult);
  tInput -> SetBranchAddress("vpdVz", &vpdVz, &bEventList_vpdVz);
  tInput -> SetBranchAddress("xVertex", &xVertex, &bEventList_xVertex);
  tInput -> SetBranchAddress("yVertex", &yVertex, &bEventList_yVertex);
  tInput -> SetBranchAddress("zVertex", &zVertex, &bEventList_zVertex);
  tInput -> SetBranchAddress("bbcZVertex", &bbcZVertex, &bEventList_bbcZVertex);
  tInput -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate, &bEventList_zdcCoincidenceRate);
  tInput -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate, &bEventList_bbcCoincidenceRate);
  tInput -> SetBranchAddress("backgroundRate", &backgroundRate, &bEventList_backgroundRate);
  tInput -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate, &bEventList_bbcBlueBackgroundRate);
  tInput -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate, &bEventList_bbcYellowBackgroundRate);
  tInput -> SetBranchAddress("refMultPos", &refMultPos, &bEventList_refMultPos);
  tInput -> SetBranchAddress("refMultNeg", &refMultNeg, &bEventList_refMultNeg);
  tInput -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity, &bEventList_bTOFTrayMultiplicity);
  tInput -> SetBranchAddress("nVerticies", &nVerticies, &bEventList_nVerticies);
  tInput -> SetBranchAddress("MagF", &MagF, &bEventList_MagF);
  tInput -> SetBranchAddress("VrtxRank", &VrtxRank, &bEventList_VrtxRank);
  tInput -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &FlagEvent_TrgTrkMisMtch, &bEventList_FlagEvent_TrgTrkMisMtch);
  tInput -> SetBranchAddress("Etsp", &Etsp, &bEventList_Etsp);
  tInput -> SetBranchAddress("ETwrdidT", &ETwrdidT, &bEventList_ETwrdidT);
  tInput -> SetBranchAddress("ETwradc11", &ETwradc11, &bEventList_ETwradc11);
  tInput -> SetBranchAddress("ETwreneT0", &ETwreneT0, &bEventList_ETwreneT0);
  tInput -> SetBranchAddress("ETwreT", &ETwreT, &bEventList_ETwreT);
  tInput -> SetBranchAddress("ETwrENET0", &ETwrENET0, &bEventList_ETwrENET0);
  tInput -> SetBranchAddress("ETwrphT", &ETwrphT, &bEventList_ETwrphT);
  tInput -> SetBranchAddress("ETwrPTower", &ETwrPTower, &bEventList_ETwrPTower);
  tInput -> SetBranchAddress("ETwrpidTower", &ETwrpidTower, &bEventList_ETwrpidTower);
  tInput -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT, &bEventList_ETwrmoduleT);
  tInput -> SetBranchAddress("EClustEneT0", &EClustEneT0, &bEventList_EClustEneT0);
  tInput -> SetBranchAddress("EClustetav1", &EClustetav1, &bEventList_EClustetav1);
  tInput -> SetBranchAddress("EClustphiv1", &EClustphiv1, &bEventList_EClustphiv1);
  tInput -> SetBranchAddress("EEstrpen01", &EEstrpen01, &bEventList_EEstrpen01);
  tInput -> SetBranchAddress("EEstrpen02", &EEstrpen02, &bEventList_EEstrpen02);
  tInput -> SetBranchAddress("EEstrpen03", &EEstrpen03, &bEventList_EEstrpen03);
  tInput -> SetBranchAddress("EEstrpen0", &EEstrpen0, &bEventList_EEstrpen0);
  tInput -> SetBranchAddress("EEstrpen1", &EEstrpen1, &bEventList_EEstrpen1);
  tInput -> SetBranchAddress("EEstrpen2", &EEstrpen2, &bEventList_EEstrpen2);
  tInput -> SetBranchAddress("EEstrpen3", &EEstrpen3, &bEventList_EEstrpen3);
  tInput -> SetBranchAddress("EEstrpen4", &EEstrpen4, &bEventList_EEstrpen4);
  tInput -> SetBranchAddress("EEstrpen5", &EEstrpen5, &bEventList_EEstrpen5);
  tInput -> SetBranchAddress("EEstrpen6", &EEstrpen6, &bEventList_EEstrpen6);
  tInput -> SetBranchAddress("EEstrpen7", &EEstrpen7, &bEventList_EEstrpen7);
  tInput -> SetBranchAddress("EEstrpen8", &EEstrpen8, &bEventList_EEstrpen8);
  tInput -> SetBranchAddress("EEstrpen9", &EEstrpen9, &bEventList_EEstrpen9);
  tInput -> SetBranchAddress("EEstrpen10", &EEstrpen10, &bEventList_EEstrpen10);
  tInput -> SetBranchAddress("EEstrpen11", &EEstrpen11, &bEventList_EEstrpen11);
  tInput -> SetBranchAddress("EEstrpen12", &EEstrpen12, &bEventList_EEstrpen12);
  tInput -> SetBranchAddress("EEstrpen13", &EEstrpen13, &bEventList_EEstrpen13);
  tInput -> SetBranchAddress("EEstrpen14", &EEstrpen14, &bEventList_EEstrpen14);
  tInput -> SetBranchAddress("EEstrpen15", &EEstrpen15, &bEventList_EEstrpen15);
  tInput -> SetBranchAddress("ETwrdidE", &ETwrdidE, &bEventList_ETwrdidE);
  tInput -> SetBranchAddress("EPstripenp01", &EPstripenp01, &bEventList_EPstripenp01);
  tInput -> SetBranchAddress("EPstripenp02", &EPstripenp02, &bEventList_EPstripenp02);
  tInput -> SetBranchAddress("EPstripenp03", &EPstripenp03, &bEventList_EPstripenp03);
  tInput -> SetBranchAddress("EPstripenp0", &EPstripenp0, &bEventList_EPstripenp0);
  tInput -> SetBranchAddress("EPstripenp1", &EPstripenp1, &bEventList_EPstripenp1);
  tInput -> SetBranchAddress("EPstripenp2", &EPstripenp2, &bEventList_EPstripenp2);
  tInput -> SetBranchAddress("EPstripenp3", &EPstripenp3, &bEventList_EPstripenp3);
  tInput -> SetBranchAddress("EPstripenp4", &EPstripenp4, &bEventList_EPstripenp4);
  tInput -> SetBranchAddress("EPstripenp5", &EPstripenp5, &bEventList_EPstripenp5);
  tInput -> SetBranchAddress("EPstripenp6", &EPstripenp6, &bEventList_EPstripenp6);
  tInput -> SetBranchAddress("EPstripenp7", &EPstripenp7, &bEventList_EPstripenp7);
  tInput -> SetBranchAddress("EPstripenp8", &EPstripenp8, &bEventList_EPstripenp8);
  tInput -> SetBranchAddress("EPstripenp9", &EPstripenp9, &bEventList_EPstripenp9);
  tInput -> SetBranchAddress("EPstripenp10", &EPstripenp10, &bEventList_EPstripenp10);
  tInput -> SetBranchAddress("EPstripenp11", &EPstripenp11, &bEventList_EPstripenp11);
  tInput -> SetBranchAddress("EPstripenp12", &EPstripenp12, &bEventList_EPstripenp12);
  tInput -> SetBranchAddress("EPstripenp13", &EPstripenp13, &bEventList_EPstripenp13);
  tInput -> SetBranchAddress("EPstripenp14", &EPstripenp14, &bEventList_EPstripenp14);
  tInput -> SetBranchAddress("EPstripenp15", &EPstripenp15, &bEventList_EPstripenp15);
  tInput -> SetBranchAddress("EclustEnnq1", &EclustEnnq1, &bEventList_EclustEnnq1);
  tInput -> SetBranchAddress("EclustEnnq20", &EclustEnnq20, &bEventList_EclustEnnq20);
  tInput -> SetBranchAddress("EclustEnnq19", &EclustEnnq19, &bEventList_EclustEnnq19);
  tInput -> SetBranchAddress("EclustEnpq1", &EclustEnpq1, &bEventList_EclustEnpq1);
  tInput -> SetBranchAddress("EclustEnpq20", &EclustEnpq20, &bEventList_EclustEnpq20);
  tInput -> SetBranchAddress("EclustEnpq19", &EclustEnpq19, &bEventList_EclustEnpq19);
  tInput -> SetBranchAddress("EclustEnpq21", &EclustEnpq21, &bEventList_EclustEnpq21);
  tInput -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &bEventList_PrimaryTrackArray_);
  tInput -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID, &bPrimaryTrackArray_fUniqueID);
  tInput -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits, &bPrimaryTrackArray_fBits);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit, &bPrimaryTrackArray_nHitsFit);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss, &bPrimaryTrackArray_nHitsPoss);
  tInput -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag, &bPrimaryTrackArray_trackFlag);
  tInput -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ, &bPrimaryTrackArray_pZ);
  tInput -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX, &bPrimaryTrackArray_pX);
  tInput -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY, &bPrimaryTrackArray_pY);
  tInput -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT, &bPrimaryTrackArray_pT);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx, &bPrimaryTrackArray_dEdx);
  tInput -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge, &bPrimaryTrackArray_charge);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta, &bPrimaryTrackArray_tofBeta);
  tInput -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta, &bPrimaryTrackArray_eta);
  tInput -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi, &bPrimaryTrackArray_phi);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron, &bPrimaryTrackArray_nSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion, &bPrimaryTrackArray_nSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon, &bPrimaryTrackArray_nSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton, &bPrimaryTrackArray_nSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag, &bPrimaryTrackArray_dcag);
  tInput -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits, &bPrimaryTrackArray_nHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits, &bPrimaryTrackArray_dEdxHits);
  tInput -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint, &bPrimaryTrackArray_firstZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint, &bPrimaryTrackArray_lastZPoint);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron, &bPrimaryTrackArray_tofSigElectron);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion, &bPrimaryTrackArray_tofSigPion);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon, &bPrimaryTrackArray_tofSigKaon);
  tInput -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton, &bPrimaryTrackArray_tofSigProton);
  tInput -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight, &bPrimaryTrackArray_timeOfflight);
  tInput -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength, &bPrimaryTrackArray_pathLength);
  tInput -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex, &bPrimaryTrackArray_trkIndex);
  tInput -> SetBranchAddress("TowerArray", &TowerArray_, &bEventList_TowerArray_);
  tInput -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID, &bTowerArray_fUniqueID);
  tInput -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits, &bTowerArray_fBits);
  tInput -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId, &bTowerArray_TwrId);
  tInput -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng, &bTowerArray_TwrEng);
  tInput -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta, &bTowerArray_TwrEta);
  tInput -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi, &bTowerArray_TwrPhi);
  tInput -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC, &bTowerArray_TwrADC);
  tInput -> SetBranchAddress("TowerArray.TwrPed", TowerArray_TwrPed, &bTowerArray_TwrPed);
  tInput -> SetBranchAddress("TowerArray.TwrRMS", TowerArray_TwrRMS, &bTowerArray_TwrRMS);
  tInput -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex, &bTowerArray_TwrMatchIdnex);
  tInput -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk, &bTowerArray_NoOfmatchedTrk);
  tInput -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP, &bTowerArray_TwrMatchP);
  tInput -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx, &bTowerArray_TwrPx);
  tInput -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy, &bTowerArray_TwrPy);
  tInput -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz, &bTowerArray_TwrPz);
  tInput -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks, &bTowerArray_fNAssocTracks);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_, &bTowerArray_fMatchedTracksArray_);
  tInput -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P, &bTowerArray_fMatchedTracksArray_P);
  cout << "      Set branches." << endl;

  // define bad run and hot tower lists
  const UInt_t badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  const UInt_t hotTwrList[NHotTwr] = {34, 106, 113, 160, 266, 267, 275, 280, 282, 286, 287, 293, 410, 504, 533, 541, 555, 561, 562, 594, 615, 616, 629, 633, 637, 638, 647, 650, 653, 657, 671, 673, 743, 789, 790, 791, 792, 806, 809, 810, 811, 812, 813, 814, 821, 822, 823, 824, 829, 830, 831, 832, 837, 841, 842, 843, 844, 846, 849, 850, 851, 852, 857, 875, 897, 899, 903, 939, 953, 954, 956, 993, 1026, 1046, 1048, 1080, 1081, 1100, 1125, 1130, 1132, 1180, 1197, 1198, 1199, 1200, 1207, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1312, 1348, 1353, 1354, 1388, 1407, 1409, 1434, 1448, 1537, 1567, 1574, 1597, 1612, 1654, 1668, 1713, 1762, 1765, 1766, 1877, 1878, 1984, 2032, 2043, 2054, 2073, 2077, 2092, 2093, 2097, 2107, 2162, 2168, 2214, 2305, 2392, 2409, 2415, 2439, 2459, 2589, 2590, 2633, 2652, 2749, 2834, 2961, 2969, 3005, 3017, 3070, 3071, 3186, 3220, 3289, 3360, 3493, 3494, 3495, 3508, 3588, 3604, 3611, 3668, 3678, 3679, 3690, 3692, 3732, 3738, 3838, 3840, 3927, 3945, 4005, 4006, 4013, 4018, 4019, 4053, 4059, 4124, 4331, 4355, 4357, 4458, 4464, 4500, 4677, 4678, 4684, 4768, 360, 493, 779, 1284, 1306, 1337, 1438, 1709, 2027, 2445, 3407, 3720, 4217, 4288, 95, 96, 296, 316, 443, 479, 555, 562, 637, 671, 709, 740, 743, 796, 857, 897, 899, 915, 953, 1130, 1132, 1294, 1318, 1337, 1348, 1359, 1378, 1427, 1429, 1440, 1537, 1563, 1574, 1709, 1763, 1773, 1819, 1854, 1874, 1936, 1938, 2018, 2043, 2098, 2099, 2256, 2259, 2294, 2514, 2520, 2552, 2589, 2598, 2680, 2706, 2799, 2880, 2897, 2917, 2969, 3020, 3028, 3310, 3319, 3375, 3399, 3504, 3539, 3541, 3679, 3690, 3692, 3718, 3719, 3720, 3738, 3806, 3838, 3840, 3928, 4013, 4017, 4038, 4053, 4057, 4058, 4079, 4097, 4099};  // big list
  //const UInt_t hotTwrList[NHotTwr] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};  // small list
  cout << "      Bad run and hot tower lists defined:\n"
       << "        " << NBadRuns << " bad runs, " << NHotTwr << " hot towers."
       << endl;

  // create jet histograms
  TH1D *hJetN[NResVals][NTrgTsp];
  TH1D *hJetA[NResVals][2];
  TH1D *hJetEta[NResVals][2];
  TH1D *hJetPt[NResVals][2];
  TH1D *hJetPtBin[NResVals][NTrgTsp];
  TH1D *hJetZtBin[NResVals][NTrgTsp];
  TH1D *hJetDfBin[NResVals][NTrgTsp][2];
  TH1D *hJetEtaBin[NResVals][NTrgTsp][2];
  TH1D *hJetDfPtBin[NResVals][NTrgTsp][NJetBins][2];
  TH1D *hJetEtaPtBin[NResVals][NTrgTsp][NJetBins][2];
  TH1D *hJetCorr[NResVals][NTrgTsp];
  TH1D *hJetSig[NResVals][NTrgTsp];
  TH2D *hJetPtVsA[NResVals][2];
  TH2D *hJetPtVsCorr[NResVals][NTrgTsp];
  TH2D *hJetPtVsDf[NResVals][NTrgTsp][2];
  TH2D *hJetPtVsEta[NResVals][NTrgTsp];
  TH2D *hJetDfVsEta[NResVals][NTrgTsp];

  // create jet histogram name bases
  const TString sJetNbase[NTrgTsp]           = {"hJetNumPi0", "hJetNumGam"};
  const TString sJetAbase[2]                 = {"hJetAall", "hJetAcut"};
  const TString sJetEtaBase[2]               = {"hJetEtaAll", "hJetEtaCut"};
  const TString sJetPtJetBase[2]             = {"hJetPtAll", "hJetPtCut"};
  const TString sJetPtBinBase[NTrgTsp]       = {"hJetPtPi0", "hJetPtGam"};
  const TString sJetZtBinBase[NTrgTsp]       = {"hJetZtPi0", "hJetZtGam"};
  const TString sJetDfBinAllBase[NTrgTsp]    = {"hJetDfAllPi0", "hJetDfAllGam"};
  const TString sJetDfBinRecBase[NTrgTsp]    = {"hJetDfRecoilPi0", "hJetDfRecoilGam"};
  const TString sJetEtaBinAllBase[NTrgTsp]   = {"hJetEtaAllPi0", "hJetEtaAllGam"};
  const TString sJetEtaBinRecBase[NTrgTsp]   = {"hJetEtaRecoilPi0", "hJetRecoilGam"};
  const TString sJetDfPtBinAllBase[NTrgTsp]  = {"hJetDfPtAllPi0", "hJetDfPtAllGam"};
  const TString sJetDfPtBinRecBase[NTrgTsp]  = {"hJetDfPtRecoilPi0", "hJetDfPtRecoilGam"};
  const TString sJetEtaPtBinAllBase[NTrgTsp] = {"hJetEtaPtAllPi0", "hJetEtaPtAllGam"};
  const TString sJetEtaPtBinRecBase[NTrgTsp] = {"hJetEtaPtRecoilPi0", "hJetEtaPtRecoilGam"};
  const TString sJetCorrBase[NTrgTsp]        = {"hJetCorrPi0", "hJetCorrGam"};
  const TString sJetSigBase[NTrgTsp]         = {"hJetSigmaPi0", "hJetSigmaGam"};
  const TString sJetPtVsAbase[2]             = {"hJetPtVsAall", "hJetPtVsAcut"};
  const TString sJetPtVsCorrBase[NTrgTsp]    = {"hJetPtVsCorrPi0", "hJetPtVsCorrGam"};
  const TString sJetPtVsDfAllBase[NTrgTsp]   = {"hJetPtVsDfAllPi0", "hJetPtVsDfAllGam"};
  const TString sJetPtVsDfRecBase[NTrgTsp]   = {"hJetPtVsDfRecoilPi0", "hJetPtVsDfRecoilGam"};
  const TString sJetPtVsEtaBase[NTrgTsp]     = {"hJetPtVsEtaPi0", "hJetPtVsEtaGam"};
  const TString sJetDfVsEtaBase[NTrgTsp]     = {"hJetDfVsEtaPi0", "hJetDfVsEtaGam"};

  // create jet histogram names
  TString sNameN[NResVals][NTrgTsp];
  TString sNameA[NResVals][2];
  TString sNameEta[NResVals][2];
  TString sNamePtJet[NResVals][2];
  TString sNamePtBin[NResVals][NTrgTsp];
  TString sNameZtBin[NResVals][NTrgTsp];
  TString sNameDfBin[NResVals][NTrgTsp][2];
  TString sNameEtaBin[NResVals][NTrgTsp][2];
  TString sNameDfPtBin[NResVals][NTrgTsp][NJetBins][2];
  TString sNameEtaPtBin[NResVals][NTrgTsp][NJetBins][2];
  TString sNameCorr[NResVals][NTrgTsp];
  TString sNameSig[NResVals][NTrgTsp];
  TString sNamePtVsA[NResVals][2];
  TString sNamePtVsCorr[NResVals][NTrgTsp];
  TString sNamePtVsDf[NResVals][NTrgTsp][2];
  TString sNamePtVsEta[NResVals][NTrgTsp];
  TString sNameDfVsEta[NResVals][NTrgTsp];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    sNameN[iRes][0]         = sJetNbase[0]         + sResParams[iRes];
    sNameN[iRes][1]         = sJetNbase[1]         + sResParams[iRes];
    sNameA[iRes][0]         = sJetAbase[0]         + sResParams[iRes];
    sNameA[iRes][1]         = sJetAbase[1]         + sResParams[iRes];
    sNameEta[iRes][0]       = sJetEtaBase[0]       + sResParams[iRes];
    sNameEta[iRes][1]       = sJetEtaBase[1]       + sResParams[iRes];
    sNamePtJet[iRes][0]     = sJetPtJetBase[0]     + sResParams[iRes];
    sNamePtJet[iRes][1]     = sJetPtJetBase[1]     + sResParams[iRes];
    sNamePtBin[iRes][0]     = sJetPtBinBase[0]     + sResParams[iRes];
    sNamePtBin[iRes][1]     = sJetPtBinBase[1]     + sResParams[iRes];
    sNameZtBin[iRes][0]     = sJetZtBinBase[0]     + sResParams[iRes];
    sNameZtBin[iRes][1]     = sJetZtBinBase[1]     + sResParams[iRes];
    sNameDfBin[iRes][0][0]  = sJetDfBinAllBase[0]  + sResParams[iRes];
    sNameDfBin[iRes][1][0]  = sJetDfBinAllBase[1]  + sResParams[iRes];
    sNameDfBin[iRes][0][1]  = sJetDfBinRecBase[0]  + sResParams[iRes];
    sNameDfBin[iRes][1][1]  = sJetDfBinRecBase[1]  + sResParams[iRes];
    sNameEtaBin[iRes][0][0] = sJetEtaBinAllBase[0] + sResParams[iRes];
    sNameEtaBin[iRes][1][0] = sJetEtaBinAllBase[1] + sResParams[iRes];
    sNameEtaBin[iRes][0][1] = sJetEtaBinRecBase[0] + sResParams[iRes];
    sNameEtaBin[iRes][1][1] = sJetEtaBinRecBase[1] + sResParams[iRes];
    sNameCorr[iRes][0]      = sJetCorrBase[0]      + sResParams[iRes];
    sNameCorr[iRes][1]      = sJetCorrBase[1]      + sResParams[iRes];
    sNameSig[iRes][0]       = sJetSigBase[0]       + sResParams[iRes];
    sNameSig[iRes][1]       = sJetSigBase[1]       + sResParams[iRes];
    sNamePtVsA[iRes][0]     = sJetPtVsAbase[0]     + sResParams[iRes];
    sNamePtVsA[iRes][1]     = sJetPtVsAbase[1]     + sResParams[iRes];
    sNamePtVsCorr[iRes][0]  = sJetPtVsCorrBase[0]  + sResParams[iRes];
    sNamePtVsCorr[iRes][1]  = sJetPtVsCorrBase[1]  + sResParams[iRes];
    sNamePtVsDf[iRes][0][0] = sJetPtVsDfAllBase[0] + sResParams[iRes];
    sNamePtVsDf[iRes][1][0] = sJetPtVsDfAllBase[1] + sResParams[iRes];
    sNamePtVsDf[iRes][0][1] = sJetPtVsDfRecBase[0] + sResParams[iRes];
    sNamePtVsDf[iRes][1][1] = sJetPtVsDfRecBase[1] + sResParams[iRes];
    sNamePtVsEta[iRes][0]   = sJetPtVsEtaBase[0]   + sResParams[iRes];
    sNamePtVsEta[iRes][1]   = sJetPtVsEtaBase[1]   + sResParams[iRes];
    sNameDfVsEta[iRes][0]   = sJetDfVsEtaBase[0]   + sResParams[iRes];
    sNameDfVsEta[iRes][1]   = sJetDfVsEtaBase[1]   + sResParams[iRes];
    // create jet dF and eta histogram names
    for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
      sNameDfPtBin[iRes][0][iJetBin][0]  = (sJetDfPtBinAllBase[0] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameDfPtBin[iRes][1][iJetBin][0]  = (sJetDfPtBinAllBase[1] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameDfPtBin[iRes][0][iJetBin][1]  = (sJetDfPtBinRecBase[0] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameDfPtBin[iRes][1][iJetBin][1]  = (sJetDfPtBinRecBase[1] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameEtaPtBin[iRes][0][iJetBin][0] = (sJetEtaPtBinAllBase[0] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameEtaPtBin[iRes][1][iJetBin][0] = (sJetEtaPtBinAllBase[1] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameEtaPtBin[iRes][0][iJetBin][1] = (sJetEtaPtBinRecBase[0] + sResParams[iRes]) + sPtBins[iJetBin];
      sNameEtaPtBin[iRes][1][iJetBin][1] = (sJetEtaPtBinRecBase[1] + sResParams[iRes]) + sPtBins[iJetBin];
    }
  }

  // binning
  const UInt_t  nNum(200);
  const UInt_t  nArea(80);
  const UInt_t  nEta(40);
  const UInt_t  nPt(275);
  const UInt_t  nPt2(110);
  const UInt_t  nPt3(55);
  const UInt_t  nZt(40);
  const UInt_t  nDf(30);
  const UInt_t  nCorr(100);
  const UInt_t  nSigma(100);
  const Float_t num[2]   = {0., 200.};
  const Float_t area[2]  = {0., 2.};
  const Float_t eta[2]   = {-2., 2.};
  const Float_t pt[2]    = {-5., 50.};
  const Float_t zt[2]    = {-2., 2.};
  const Float_t df[2]    = {-1.6, 4.7};
  const Float_t corr[2]  = {0., 10.};
  const Float_t sigma[2] = {0., 10.};
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    hJetN[iRes][0]         = new TH1D(sNameN[iRes][0].Data(), "", nNum, num[0], num[1]);
    hJetN[iRes][1]         = new TH1D(sNameN[iRes][1].Data(), "", nNum, num[0], num[1]);
    hJetA[iRes][0]         = new TH1D(sNameA[iRes][0].Data(), "", nArea, area[0], area[1]);
    hJetA[iRes][1]         = new TH1D(sNameA[iRes][1].Data(), "", nArea, area[0], area[1]);
    hJetEta[iRes][0]       = new TH1D(sNameEta[iRes][0].Data(), "", nEta, eta[0], eta[1]);
    hJetEta[iRes][1]       = new TH1D(sNameEta[iRes][1].Data(), "", nEta, eta[0], eta[1]);
    hJetPt[iRes][0]        = new TH1D(sNamePtJet[iRes][0].Data(), "", nPt, pt[0], pt[1]);
    hJetPt[iRes][1]        = new TH1D(sNamePtJet[iRes][1].Data(), "", nPt, pt[0], pt[1]);
    hJetPtBin[iRes][0]     = new TH1D(sNamePtBin[iRes][0].Data(), "", nPt2, pt[0], pt[1]);
    hJetPtBin[iRes][1]     = new TH1D(sNamePtBin[iRes][1].Data(), "", nPt2, pt[0], pt[1]);
    hJetZtBin[iRes][0]     = new TH1D(sNameZtBin[iRes][0].Data(), "", nZt, zt[0], zt[1]);
    hJetZtBin[iRes][1]     = new TH1D(sNameZtBin[iRes][1].Data(), "", nZt, zt[0], zt[1]);
    hJetDfBin[iRes][0][0]  = new TH1D(sNameDfBin[iRes][0][0].Data(), "", nDf, df[0], df[1]);
    hJetDfBin[iRes][1][0]  = new TH1D(sNameDfBin[iRes][1][0].Data(), "", nDf, df[0], df[1]);
    hJetDfBin[iRes][0][1]  = new TH1D(sNameDfBin[iRes][0][1].Data(), "", nDf, df[0], df[1]);
    hJetDfBin[iRes][1][1]  = new TH1D(sNameDfBin[iRes][1][1].Data(), "", nDf, df[0], df[1]);
    hJetEtaBin[iRes][0][0] = new TH1D(sNameEtaBin[iRes][0][0].Data(), "", nEta, eta[0], eta[1]);
    hJetEtaBin[iRes][1][0] = new TH1D(sNameEtaBin[iRes][1][0].Data(), "", nEta, eta[0], eta[1]);
    hJetEtaBin[iRes][0][1] = new TH1D(sNameEtaBin[iRes][0][1].Data(), "", nEta, eta[0], eta[1]);
    hJetEtaBin[iRes][1][1] = new TH1D(sNameEtaBin[iRes][1][1].Data(), "", nEta, eta[0], eta[1]);
    hJetCorr[iRes][0]      = new TH1D(sNameCorr[iRes][0].Data(), "", nCorr, corr[0], corr[1]);
    hJetCorr[iRes][1]      = new TH1D(sNameCorr[iRes][1].Data(), "", nCorr, corr[0], corr[1]);
    hJetSig[iRes][0]       = new TH1D(sNameSig[iRes][0].Data(), "", nSigma, sigma[0], sigma[1]);
    hJetSig[iRes][1]       = new TH1D(sNameSig[iRes][1].Data(), "", nSigma, sigma[0], sigma[1]);
    hJetPtVsA[iRes][0]     = new TH2D(sNamePtVsA[iRes][0].Data(), "", nArea, area[0], area[1], nPt3, pt[0], pt[1]);
    hJetPtVsA[iRes][1]     = new TH2D(sNamePtVsA[iRes][1].Data(), "", nArea, area[0], area[1], nPt3, pt[0], pt[1]);
    hJetPtVsCorr[iRes][0]  = new TH2D(sNamePtVsCorr[iRes][0].Data(), "", nCorr, corr[0], corr[1], nPt3, pt[0], pt[1]);
    hJetPtVsCorr[iRes][1]  = new TH2D(sNamePtVsCorr[iRes][1].Data(), "", nCorr, corr[0], corr[1], nPt3, pt[0], pt[1]);
    hJetPtVsDf[iRes][0][0] = new TH2D(sNamePtVsDf[iRes][0][0].Data(), "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
    hJetPtVsDf[iRes][1][0] = new TH2D(sNamePtVsDf[iRes][1][0].Data(), "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
    hJetPtVsDf[iRes][0][1] = new TH2D(sNamePtVsDf[iRes][0][1].Data(), "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
    hJetPtVsDf[iRes][1][1] = new TH2D(sNamePtVsDf[iRes][1][1].Data(), "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
    hJetPtVsEta[iRes][0]   = new TH2D(sNamePtVsEta[iRes][0].Data(), "", nEta, eta[0], eta[1], nPt3, pt[0], pt[1]);
    hJetPtVsEta[iRes][1]   = new TH2D(sNamePtVsEta[iRes][1].Data(), "", nEta, eta[0], eta[1], nPt3, pt[0], pt[1]);
    hJetDfVsEta[iRes][0]   = new TH2D(sNameDfVsEta[iRes][0].Data(), "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
    hJetDfVsEta[iRes][1]   = new TH2D(sNameDfVsEta[iRes][1].Data(), "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
    for (UInt_t iTrgTsp = 0; iTrgTsp < NTrgTsp; iTrgTsp++) {
      for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
        hJetDfPtBin[iRes][iTrgTsp][iJetBin][0]  = new TH1D(sNameDfPtBin[iRes][iTrgTsp][iJetBin][0].Data(), "", nDf, df[0], df[1]);
        hJetDfPtBin[iRes][iTrgTsp][iJetBin][1]  = new TH1D(sNameDfPtBin[iRes][iTrgTsp][iJetBin][1].Data(), "", nDf, df[0], df[1]);
        hJetEtaPtBin[iRes][iTrgTsp][iJetBin][0] = new TH1D(sNameEtaPtBin[iRes][iTrgTsp][iJetBin][0].Data(), "", nEta, eta[0], eta[1]);
        hJetEtaPtBin[iRes][iTrgTsp][iJetBin][1] = new TH1D(sNameEtaPtBin[iRes][iTrgTsp][iJetBin][1].Data(), "", nEta, eta[0], eta[1]);
      }
    }
    // errors
    hJetN[iRes][0]         -> Sumw2();
    hJetN[iRes][1]         -> Sumw2();
    hJetA[iRes][0]         -> Sumw2();
    hJetA[iRes][1]         -> Sumw2();
    hJetEta[iRes][0]       -> Sumw2();
    hJetEta[iRes][1]       -> Sumw2();
    hJetPt[iRes][0]        -> Sumw2();
    hJetPt[iRes][1]        -> Sumw2();
    hJetPtBin[iRes][0]     -> Sumw2();
    hJetPtBin[iRes][1]     -> Sumw2();
    hJetZtBin[iRes][0]     -> Sumw2();
    hJetZtBin[iRes][1]     -> Sumw2();
    hJetDfBin[iRes][0][0]  -> Sumw2();
    hJetDfBin[iRes][1][0]  -> Sumw2();
    hJetDfBin[iRes][0][1]  -> Sumw2();
    hJetDfBin[iRes][1][1]  -> Sumw2();
    hJetEtaBin[iRes][0][0] -> Sumw2();
    hJetEtaBin[iRes][1][0] -> Sumw2();
    hJetEtaBin[iRes][0][1] -> Sumw2();
    hJetEtaBin[iRes][1][1] -> Sumw2();
    hJetCorr[iRes][0]      -> Sumw2();
    hJetCorr[iRes][1]      -> Sumw2();
    hJetSig[iRes][0]       -> Sumw2();
    hJetSig[iRes][1]       -> Sumw2();
    hJetPtVsA[iRes][0]     -> Sumw2();
    hJetPtVsA[iRes][1]     -> Sumw2();
    hJetPtVsCorr[iRes][0]  -> Sumw2();
    hJetPtVsCorr[iRes][1]  -> Sumw2();
    hJetPtVsDf[iRes][0][0] -> Sumw2();
    hJetPtVsDf[iRes][1][0] -> Sumw2();
    hJetPtVsDf[iRes][0][1] -> Sumw2();
    hJetPtVsDf[iRes][1][1] -> Sumw2();
    hJetPtVsEta[iRes][0]   -> Sumw2();
    hJetPtVsEta[iRes][1]   -> Sumw2();
    hJetDfVsEta[iRes][0]   -> Sumw2();
    hJetDfVsEta[iRes][1]   -> Sumw2();
    for (UInt_t iTrgTsp = 0; iTrgTsp < NTrgTsp; iTrgTsp++) {
      for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
        hJetDfPtBin[iRes][iTrgTsp][iJetBin][0]  -> Sumw2();
        hJetDfPtBin[iRes][iTrgTsp][iJetBin][1]  -> Sumw2();
        hJetEtaPtBin[iRes][iTrgTsp][iJetBin][0] -> Sumw2();
        hJetEtaPtBin[iRes][iTrgTsp][iJetBin][1] -> Sumw2();
      }  // end pT bin loop
    }  // end tsp loop
  }  // end rJet loop

  // no. of evts and triggers
  UInt_t nTrgPi0[NTrgBins];
  UInt_t nTrgGam[NTrgBins];
  for (UInt_t iBin = 0; iBin < NTrgBins; iBin++) {
    nTrgPi0[iBin] = 0;
    nTrgGam[iBin] = 0;
  }

  const UInt_t nEvts = tInput -> GetEntriesFast();
  cout << "      Beginning event loop: " << nEvts << " events to process." << endl;

  // for jet-finding
  vector<PseudoJet> jetCons[NResVals];
  vector<PseudoJet> jetCS[NResVals];
  vector<PseudoJet> jets[NResVals];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    jetCons[iRes].clear();
    jetCS[iRes].clear();
    jets[iRes].clear();
  }

  // event loop
  UInt_t bytes(0);
  UInt_t nBytes(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // load entry
    bytes   = tInput -> GetEntry(iEvt);
    nBytes += bytes;
    if (bytes < 0) {
      cerr << "WARNING: issue with entry " << iEvt << "!" << endl;
      break;
    } else {
      if (isInBatchMode) {
        cout << "        Processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
      } else {
        cout << "        Processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
        if ((iEvt + 1) == nEvts) cout << endl;
      }
    }

    // event info
    const UInt_t   run   = runNumber;
    const Long64_t nTrks = nPrimaryTracks;
    const Double_t xVtx  = xVertex;
    const Double_t yVtx  = yVertex;
    const Double_t zVtx  = zVertex;
    const Double_t rVtx  = TMath::Sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // filter out bad runs
    Bool_t isGoodRun = true;
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (run == badRunList[iRun]) {
        isGoodRun = false;
        break;
      }
    }
    if (!isGoodRun) continue;

    // vertex cuts
    const Bool_t isInRcut  = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut  = (TMath::Abs(zVtx) < zVtxMax);
    const Bool_t isGoodEvt = (isInRcut && isInZcut);
    if (!isGoodEvt) continue;

    // trigger info
    const Int_t    adc    = ETwradc11;
    const UInt_t   idTrg  = ETwrdidT;
    const Double_t tspTrg = Etsp;
    const Double_t eH4    = EEstrpen4;
    const Double_t eF4    = EPstripenp4;
    const Double_t pProj  = ETwrPTower;
    const Double_t hDet   = ETwreT;
    const Double_t hPhys  = EClustetav1;
    const Double_t fTrg   = EClustphiv1;
    const Double_t eTrg   = EClustEneT0;
    const Double_t tTrg   = 2. * TMath::ATan(TMath::Exp(-1. * hPhys));
    const Double_t eTtrg  = eTrg * TMath::Sin(tTrg);

    // filter out bad towers
    Bool_t isGoodTwr = true;
    for (UInt_t iTwr = 0; iTwr < NHotTwr; iTwr++) {
      if (idTrg == hotTwrList[iTwr]) {
        isGoodTwr = false;
        break;
      }
    }
    if (!isGoodTwr) continue;

    // trigger cuts
    const Bool_t isInAdcCut    = (adc <= adcMax);
    const Bool_t isInStrCut    = ((eH4 >= eStrMin) && (eF4 >= eStrMin));
    const Bool_t isInProjCut   = (pProj < pProjMax);
    const Bool_t isInEtaTrgCut = ((TMath::Abs(hDet) < hTrgMax) && (TMath::Abs(hPhys) < hTrgMax));
    const Bool_t isInEtCut     = ((eTtrg >= eTtrgMin) && (eTtrg < eTtrgMax));
    const Bool_t isInPi0cut    = ((tspTrg > tspPi0[0]) && (tspTrg < tspPi0[1]));
    const Bool_t isInGamCut    = ((tspTrg > tspGam[0]) && (tspTrg < tspGam[1]));
    const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);
    const Bool_t isGoodTrigger = (isInAdcCut && isInStrCut && isInProjCut && isInEtaTrgCut && isInEtCut && isInTspCut);
    if (!isGoodTrigger) continue;

    // figure out eT bin
    UInt_t iEtBin(5);
    for (UInt_t iBin = 1; iBin < NTrgBins; iBin++) {
      const Bool_t isInBin = ((eTtrg >= eTbinMin[iBin]) && (eTtrg < eTbinMax[iBin]));
      if (isInBin) {
        iEtBin = iBin;
        break;
      }
    }

    // count triggers
    if (isInPi0cut) {
      nTrgPi0[0]++;
      nTrgPi0[iEtBin]++;
    }
    if (isInGamCut) {
      nTrgGam[0]++;
      nTrgGam[iEtBin]++;
    }

    // reset jet-finding vectors
    for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
      jetCons[iRes].clear();
      jetCS[iRes].clear();
      jets[iRes].clear();
    }

    // track loop
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFitTrk  = PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPossTrk = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFitTrk  = (Double_t) nFitTrk / (Double_t) nPossTrk;
      const Double_t dcaTrk   = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk     = PrimaryTrackArray_eta[iTrk];
      const Double_t pTtrk    = PrimaryTrackArray_pT[iTrk];
      const Double_t pXtrk    = PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk    = PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk    = PrimaryTrackArray_pZ[iTrk];
      const Double_t eTrk     = TMath::Sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (0.140 * 0.140));

      // track cuts
      const Bool_t isInFitCut    = (nFitTrk >= nFitMin);
      const Bool_t isInRatioCut  = (rFitTrk >= rFitMin);
      const Bool_t isInDcaCut    = (dcaTrk < dcaMax);
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtCut     = ((pTtrk > pTtrkMin) && (pTtrk < pTtrkMax));
      const Bool_t isGoodTrk     = (isInFitCut && isInRatioCut && isInDcaCut && isInEtaTrkCut && isInPtCut);
      if (!isGoodTrk) continue;

      // for jet-finding
      for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
        jetCons[iRes].push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));
      }
    }  // end track loop

    // jet clustering
    Double_t rhoJet[NResVals];
    Double_t sigJet[NResVals];

    // loop over resolution parameters
    UInt_t nJetsJet[NResVals];
    for (UInt_t iRes = 0; iRes < NResVals; iRes++) {

      // define jets and jet area
      GhostedAreaSpec areaSpec(hGhostMax[iRes], nRepeat, aGhost);
      AreaDefinition  areaDef(active_area_explicit_ghosts, areaSpec);
      JetDefinition   jetDef(antikt_algorithm, rJet[iRes]);

      // cluster jets
      ClusterSequenceArea clusterSeq(jetCons[iRes], jetDef, areaDef);
      jetCS[iRes] = sorted_by_pt(clusterSeq.inclusive_jets(pTjetMin));

      // fiducial cut
      Selector fidCut = SelectorAbsEtaMax(hJetMax[iRes]);
      jets[iRes]      = fidCut(jetCS[iRes]);

      // define bkgd jets and jet area
      GhostedAreaSpec areaSpecBkgd(hGhostMax[iRes], nRepeat, aGhost);
      AreaDefinition  areaDefBkgd(active_area_explicit_ghosts, areaSpecBkgd);
      JetDefinition   jetDefBkgd(kt_algorithm, rJet[iRes]);

      // initialize bkgd estimators
      Selector bkgdCut = SelectorAbsEtaMax(hBkgdMax) * (!SelectorNHardest(nRemove));
      JetMedianBackgroundEstimator bkgd(bkgdCut, jetDefBkgd, areaDefBkgd);
      Subtractor sub(&bkgd);

#if FASTJET_VERSION_NUMBER >= 30100
      sub.set_use_rho_m(true);
      sub.set_safe_mass(true);
#endif

      // estimate bkgd
      bkgd.set_particles(jetCons[iRes]);
      rhoJet[iRes] = bkgd.rho();
      sigJet[iRes] = bkgd.sigma();
      if (isInPi0cut) {
        hJetSig[iRes][0] -> Fill(sigJet[iRes]);
      }
      if (isInGamCut) {
        hJetSig[iRes][1] -> Fill(sigJet[iRes]);
      }

      // jet loop
      nJetsJet[iRes] = (UInt_t) jets[iRes].size();

       UInt_t nRecJetP(0);
       UInt_t nRecJetG(0);
      for (UInt_t iJet = 0; iJet < nJetsJet[iRes]; iJet++) {

        // jet info
        const Double_t fJet   = jets[iRes][iJet].phi_std();
        const Double_t hJet   = jets[iRes][iJet].pseudorapidity();
        const Double_t aJet   = jets[iRes][iJet].area();
        const Double_t pTjet  = jets[iRes][iJet].perp();
        const Double_t pTcorr = rhoJet[iRes] * aJet;
        const Double_t pTrec  = pTjet - pTcorr;
        const Double_t zTrec  = pTrec / eTtrg;

        Double_t dFjet = fJet - fTrg;
        if (dFjet < dFjetMin) dFjet += tau;
        if (dFjet > dFjetMax) dFjet -= tau;

        Double_t dFre = fJet - fTrg;
        if (dFre < 0.)  dFre += tau;
        if (dFre > tau) dFre -= tau;
        const Double_t dFcut = TMath::Abs(dFre - pi);

        // jet cuts
        const Bool_t isInEtaJetCut = (TMath::Abs(hJet) < hJetMax[iRes]);
        const Bool_t isInPtJetCut  = ((pTjet > pTjetMin) && (pTjet < pTjetMax));
        const Bool_t isInAjetCut   = (aJet > aJetMin[iRes]);
        const Bool_t isRecoilJet   = (dFcut < dFrecoilMax);

        // eta cut
        hJetEta[iRes][0] -> Fill(hJet);
        if (!isInEtaJetCut) {
          hJetEta[iRes][1] -> Fill(hJet);
        }

        // pT cut
        hJetPt[iRes][0] -> Fill(pTjet);
        if (!isInPtJetCut) {
          hJetPt[iRes][1] -> Fill(pTjet);
        }

        // area cut
        hJetA[iRes][0]     -> Fill(aJet);
        hJetPtVsA[iRes][0] -> Fill(aJet, pTjet);
        if (!isInAjetCut) {
          hJetA[iRes][1]     -> Fill(aJet);
          hJetPtVsA[iRes][1] -> Fill(aJet, pTjet);
        }

        // for pT-differential dF histograms
        Bool_t isInJetPtBin[NJetBins];
        for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
          isInJetPtBin[iJetBin] = ((pTjet > pTbinMin[iJetBin]) && (pTjet < pTbinMax[iJetBin]));
        }

        // fill histograms
        if (isInPi0cut && isInAjetCut) {
          hJetPtVsEta[iRes][0] -> Fill(hJet, pTrec);
          hJetDfVsEta[iRes][0] -> Fill(hJet, dFjet);
          if (isInEtaJetCut) {
            hJetDfBin[iRes][0][0]  -> Fill(dFjet);
            hJetEtaBin[iRes][0][0] -> Fill(hJet);
            hJetPtVsDf[iRes][0][0] -> Fill(dFjet, pTrec);
            if (isRecoilJet) {
              hJetDfBin[iRes][0][1]  -> Fill(dFjet);
              hJetEtaBin[iRes][0][1] -> Fill(hJet);
              hJetPtVsDf[iRes][0][1] -> Fill(dFjet, pTrec);
              if (isInPtJetCut) {
                hJetCorr[iRes][0]     -> Fill(pTcorr);
                hJetPtBin[iRes][0]    -> Fill(pTrec);
                hJetZtBin[iRes][0]    -> Fill(zTrec);
                hJetPtVsCorr[iRes][0] -> Fill(pTcorr, pTjet);
                ++nRecJetP;
              }  // end pT cut
            }  // end recoil cut
            for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
              if (isInJetPtBin[iJetBin]) {
                hJetDfPtBin[iRes][0][iJetBin][0]  -> Fill(dFjet);
                hJetEtaPtBin[iRes][0][iJetBin][0] -> Fill(hJet);
                if (isRecoilJet) {
                  hJetDfPtBin[iRes][0][iJetBin][1]  -> Fill(dFjet);
                  hJetEtaPtBin[iRes][0][iJetBin][1] -> Fill(hJet);
                }
              }
            }  // end pT bin loop
          }  // end eta cut
        }  // end pi0 and area cut
        if (isInGamCut && isInAjetCut) {
          hJetPtVsEta[iRes][1] -> Fill(hJet, pTrec);
          hJetDfVsEta[iRes][1] -> Fill(hJet, dFjet);
          if (isInEtaJetCut) {
            hJetDfBin[iRes][1][0]  -> Fill(dFjet);
            hJetEtaBin[iRes][1][0] -> Fill(hJet);
            hJetPtVsDf[iRes][1][0] -> Fill(dFjet, pTrec);
            if (isRecoilJet) {
              hJetDfBin[iRes][1][1]  -> Fill(dFjet);
              hJetEtaBin[iRes][1][1] -> Fill(hJet);
              hJetPtVsDf[iRes][1][1] -> Fill(dFjet, pTrec);
              if (isInPtJetCut) {
                hJetCorr[iRes][1]     -> Fill(pTcorr);
                hJetPtBin[iRes][1]    -> Fill(pTrec);
                hJetZtBin[iRes][1]    -> Fill(zTrec);
                hJetPtVsCorr[iRes][1] -> Fill(pTcorr, pTjet);
                ++nRecJetG;
              }  // end pT cut
            }  // end recoil cut
            for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
              if (isInJetPtBin[iJetBin]) {
                hJetDfPtBin[iRes][1][iJetBin][0]  -> Fill(dFjet);
                hJetEtaPtBin[iRes][1][iJetBin][0] -> Fill(hJet);
                if (isRecoilJet) {
                  hJetDfPtBin[iRes][1][iJetBin][1]  -> Fill(dFjet);
                  hJetEtaPtBin[iRes][1][iJetBin][1] -> Fill(hJet);
                }
              }
            }  // end pT bin loop
          }  // end eta cut
        }  // end gamma and area cut
      }  // end jet loop

      // fill nJet histograms
      if (isInPi0cut) {
        hJetN[iRes][0] -> Fill(nRecJetP);
      }
      if (isInGamCut) {
        hJetN[iRes][1] -> Fill(nRecJetG);
      }

    }  // end rJet loop
  }  // end event loop

  // report no. of pi0's and gamma's
  cout << "      Number of triggers:" << endl;
  for (UInt_t iBin = 0; iBin < NTrgBins; iBin++) {
    cout << "        eTtrg = (" << eTbinMin[iBin] << ", " << eTbinMax[iBin]
         << ") GeV: nPi0 = " << nTrgPi0[iBin] << ", nGam = " << nTrgGam[iBin]
         << endl;
  }

  // normalize relevant histograms
  const Double_t aBin    = (area[1] - area[0]) / (Double_t) nArea;
  const Double_t hBin    = (eta[1] - eta[0]) / (Double_t) nEta;
  const Double_t pTbin2  = (pt[1] - pt[0]) / (Double_t) nPt2;
  const Double_t pTbin3  = (pt[1] - pt[0]) / (Double_t) nPt3;
  const Double_t zTbin   = (zt[1] - zt[0]) / (Double_t) nZt;
  const Double_t dFbin   = (df[1] - df[0]) / (Double_t) nDf;
  const Double_t corrBin = (corr[1] - corr[0]) / (Double_t) nCorr;
  const Double_t sigBin  = (sigma[1] - sigma[0]) / (Double_t) nSigma;
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    for (UInt_t iTsp = 0; iTsp < NTrgTsp; iTsp++) {
      // no. of triggers
      UInt_t nTrg(1);
      if (iTsp == 0) nTrg = nTrgPi0[0];
      if (iTsp == 1) nTrg = nTrgGam[0];

      // norms
      const Double_t nNorm      = 1. / (Double_t) nTrg;
      const Double_t pTnorm2    = 1. / (pTbin2 * nTrg);
      const Double_t zTnorm     = 1. / (zTbin * nTrg);
      const Double_t dFnorm     = 1. / (dFbin * nTrg);
      const Double_t hNorm      = 1. / (hBin * nTrg);
      const Double_t corrNorm   = 1. / (corrBin * nTrg);
      const Double_t sigNorm    = 1. / (sigBin * nTrg);
      const Double_t hPtNorm    = 1. / (hBin * pTbin3 * nTrg);
      const Double_t hDfNorm    = 1. / (hBin * dFbin * nTrg);
      const Double_t dFpTnorm   = 1. / (dFbin * pTbin3 * nTrg);
      const Double_t corrPtNorm = 1. / (corrBin * pTbin3 * nTrg);
      hJetN[iRes][iTsp]         -> Scale(nNorm);
      hJetPtBin[iRes][iTsp]     -> Scale(pTnorm2);
      hJetZtBin[iRes][iTsp]     -> Scale(zTnorm);
      hJetDfBin[iRes][iTsp][0]  -> Scale(dFnorm);
      hJetDfBin[iRes][iTsp][1]  -> Scale(dFnorm);
      hJetEtaBin[iRes][iTsp][0] -> Scale(hNorm);
      hJetEtaBin[iRes][iTsp][1] -> Scale(hNorm);
      hJetCorr[iRes][iTsp]      -> Scale(corrNorm);
      hJetSig[iRes][iTsp]       -> Scale(sigNorm);
      hJetPtVsCorr[iRes][iTsp]  -> Scale(corrPtNorm);
      hJetPtVsDf[iRes][iTsp][0] -> Scale(dFpTnorm);
      hJetPtVsDf[iRes][iTsp][1] -> Scale(dFpTnorm);
      hJetPtVsEta[iRes][iTsp]   -> Scale(hPtNorm);
      hJetDfVsEta[iRes][iTsp]   -> Scale(hDfNorm);
      for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
        hJetDfPtBin[iRes][iTsp][iJetBin][0]  -> Scale(dFnorm);
        hJetDfPtBin[iRes][iTsp][iJetBin][1]  -> Scale(dFnorm);
        hJetEtaPtBin[iRes][iTsp][iJetBin][0] -> Scale(hNorm);
        hJetEtaPtBin[iRes][iTsp][iJetBin][1] -> Scale(hNorm);
      }
    }
    const UInt_t   nTrgTot = nTrgPi0[0] + nTrgGam[0];
    const Double_t aPtNorm = 1. / (aBin * pTbin3 * nTrgTot);
    hJetPtVsA[iRes][0] -> Scale(aPtNorm);
    hJetPtVsA[iRes][1] -> Scale(aPtNorm);
  }
  cout << "      Normalized certain histograms." << endl;

  // find minimum / maximum for 2d comparison plots
  Float_t zCorrVsPt[NResVals][2];
  Float_t zDfVsPt[NResVals][2];
  Float_t zEtaVsPt[NResVals][2];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    const Float_t minCorrP = hJetPtVsCorr[iRes][0] -> GetMinimum(0.);
    const Float_t minCorrG = hJetPtVsCorr[iRes][1] -> GetMinimum(0.);
    const Float_t maxCorrP = hJetPtVsCorr[iRes][0] -> GetMaximum();
    const Float_t maxCorrG = hJetPtVsCorr[iRes][1] -> GetMaximum();
    zCorrVsPt[iRes][0] = TMath::Min(minCorrP, minCorrG);
    zCorrVsPt[iRes][1] = TMath::Max(maxCorrP, maxCorrG);

    const Float_t minDfPA = hJetPtVsDf[iRes][0][0] -> GetMinimum(0.);
    const Float_t minDfPR = hJetPtVsDf[iRes][0][1] -> GetMinimum(0.);
    const Float_t minDfGA = hJetPtVsDf[iRes][1][0] -> GetMinimum(0.);
    const Float_t minDfGR = hJetPtVsDf[iRes][1][1] -> GetMinimum(0.);
    const Float_t minDfP  = TMath::Min(minDfPA, minDfPR);
    const Float_t minDfG  = TMath::Min(minDfGA, minDfGR);
    const Float_t maxDfPA = hJetPtVsDf[iRes][0][0] -> GetMaximum();
    const Float_t maxDfPR = hJetPtVsDf[iRes][0][1] -> GetMaximum();
    const Float_t maxDfGA = hJetPtVsDf[iRes][1][0] -> GetMaximum();
    const Float_t maxDfGR = hJetPtVsDf[iRes][1][1] -> GetMaximum();
    const Float_t maxDfP  = TMath::Max(maxDfPA, maxDfPR);
    const Float_t maxDfG  = TMath::Max(maxDfGA, maxDfGR);
    zDfVsPt[iRes][0] = TMath::Min(minDfP, minDfG);
    zDfVsPt[iRes][1] = TMath::Max(maxDfP, maxDfG);

    const Float_t minEtaP = hJetPtVsEta[iRes][0] -> GetMinimum(0.);
    const Float_t minEtaG = hJetPtVsEta[iRes][1] -> GetMinimum(0.);
    const Float_t maxEtaP = hJetPtVsEta[iRes][0] -> GetMaximum();
    const Float_t maxEtaG = hJetPtVsEta[iRes][1] -> GetMaximum();
    zEtaVsPt[iRes][0] = TMath::Min(minEtaP, minEtaG);
    zEtaVsPt[iRes][1] = TMath::Max(maxEtaP, maxEtaG);
  }

  Float_t zAreaVsPt[2] = {999., 0.};
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    const Float_t minAreaA = hJetPtVsA[iRes][0] -> GetMinimum(0.);
    const Float_t maxAreaA = hJetPtVsA[iRes][0] -> GetMaximum();
    const Float_t minAreaC = hJetPtVsA[iRes][1] -> GetMinimum(0.);
    const Float_t maxAreaC = hJetPtVsA[iRes][1] -> GetMaximum();
    const Float_t minAreaP = zAreaVsPt[0];
    const Float_t maxAreaP = zAreaVsPt[1];

    const Float_t newMin = TMath::Min(minAreaA, minAreaC);
    const Float_t newMax = TMath::Max(maxAreaA, maxAreaC);
    zAreaVsPt[0] = TMath::Min(newMin, minAreaP);
    zAreaVsPt[1] = TMath::Max(newMax, maxAreaP);
  }

  cout << "      Determined z-axis ranges for 2d comparisons." << endl;
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    cout << "        " << sDir[iRes].Data() << " jets:\n"
         << "          pT vs. corr: min = " << zCorrVsPt[iRes][0] << ", max = " << zCorrVsPt[iRes][1] << "\n"
         << "          pT vs. dF:   min = " << zDfVsPt[iRes][0]   << ", max = " << zDfVsPt[iRes][1] << "\n"
         << "          pT vs. eta:  min = " << zEtaVsPt[iRes][0]  << ", max = " << zEtaVsPt[iRes][1]
         << endl;
  }
  cout << "        pT vs. area:\n"
       << "          min = " << zAreaVsPt[0] << ", max = " << zAreaVsPt[1]
       << endl;

  // set styles
  const UInt_t   fColAll(1);
  const UInt_t   fColCut(632);
  const UInt_t   fLinAll(1);
  const UInt_t   fLinCut(1);
  const UInt_t   fFilAll(0);
  const UInt_t   fFilCut(3345);
  const UInt_t   fMarAll(8);
  const UInt_t   fMarCut(8);
  const UInt_t   fTxt(42);
  const UInt_t   fCnt(1);
  const UInt_t   fColTrg[NTrgTsp]        = {859, 899};
  const UInt_t   fColDfAll[NTrgTsp]      = {1, 1};
  const UInt_t   fColDfRec[NTrgTsp]      = {859, 899};
  const UInt_t   fColAllDfPi0[NTrkBins]  = {920, 921, 922, 923};
  const UInt_t   fColAllDfGam[NTrkBins]  = {920, 921, 922, 923};
  const UInt_t   fColRecDfPi0[NTrkBins]  = {590, 591, 593, 596};
  const UInt_t   fColRecDfGam[NTrkBins]  = {622, 623, 625, 628};
  const UInt_t   fColAllEtaPi0[NTrkBins] = {590, 591, 593, 596};
  const UInt_t   fColAllEtaGam[NTrkBins] = {622, 623, 625, 628};
  const UInt_t   fColRecEtaPi0[NTrkBins] = {590, 591, 593, 596};
  const UInt_t   fColRecEtaGam[NTrkBins] = {622, 623, 625, 628};
  const UInt_t   fMarDfAll[NTrgTsp]   = {8, 8};
  const UInt_t   fMarDfRec[NTrgTsp]   = {8, 8};
  const UInt_t   fMarAllPi0[NTrkBins] = {20, 22, 23, 21};
  const UInt_t   fMarAllGam[NTrkBins] = {20, 22, 23, 21};
  const UInt_t   fMarRecPi0[NTrkBins] = {20, 22, 23, 21};
  const UInt_t   fMarRecGam[NTrkBins] = {20, 22, 23, 21};
  const TString  sJetN("N_{recoil}");
  const TString  sJetA("A_{jet}");
  const TString  sJetEta("#eta^{jet}");
  const TString  sJetPt("p_{T}^{jet} [GeV/c]");
  const TString  sJetZt("z_{T}^{jet} = p_{T}^{reco} / E_{T}^{trg}");
  const TString  sJetDf("#Delta#varphi^{jet}");
  const TString  sJetSig("#sigma");
  const TString  sPtReco("p_{T}^{reco} = p_{T}^{jet} - (#rho #upoint A^{jet}) [GeV/c]");
  const TString  sPtCorr("#rho #upoint A^{jet} [GeV/c]");
  const TString  sJetNY("arbitrary units");
  const TString  sJetPtY("(1/N^{trg}) dN^{jet}/dp_{T}^{jet} [GeV/c]^{-1}");
  const TString  sJetZtY("(1/N^{trg}) dN^{jet}/dz_{T}^{jet} [GeV/c]^{-1}");
  const TString  sJetDfY("(1/N^{trg}) dN^{jet}/d#Delta#varphi^{jet}");
  const TString  sJetEtaY("(1/N^{trg}) dN^{jet}/d#eta^{jet} [GeV/c]^{-1}");
  const TString  sPtRecoY("(1/N^{trg}) dN^{jet}/dp_{T}^{reco} [GeV/c]^{-1}");
  const TString  sPtCorrY("(1/N^{trg}) dN^{jet}/d(#rho #upoint A^{jet}) [GeV/c]^{-1}");
  const TString  sJetSigY("arbitrary units");
  const TString  sCount("counts");
  const Double_t fLab(0.04);
  const Double_t fTit(0.04);
  const Double_t fOffX(1.1);
  const Double_t fOffY(1.3);
  const Double_t fOffZ(1.);
  const Double_t fOffL(0.007);
  const Double_t xEtaPlot[2]            = {-1.7, 1.7};
  const Double_t xAreaPlot[NResVals][2] = {{0., 0.625}, {0., 1.325}};
  const Double_t xCorrPlot[2]           = {0., 1.7};
  const Double_t yDfPlot[2]             = {0.00003, 7777.};
  const Double_t yEtaPlot[2]            = {0.00003, 7777.};
  const Double_t yPtPlot[2]             = {-1., 47.};
  const Double_t yCorrPlot[NResVals][2] = {{3., 37.}, {0.0003, 17.}};
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {

    // jet number
    hJetN[iRes][0] -> SetLineColor(fColTrg[0]);
    hJetN[iRes][0] -> SetLineStyle(fLinAll);
    hJetN[iRes][0] -> SetFillColor(fColTrg[0]);
    hJetN[iRes][0] -> SetFillStyle(fFilAll);
    hJetN[iRes][0] -> SetMarkerColor(fColTrg[0]);
    hJetN[iRes][0] -> SetMarkerStyle(fMarAll);
    hJetN[iRes][0] -> SetTitle("");
    hJetN[iRes][0] -> SetTitleFont(fTxt);
    hJetN[iRes][0] -> GetXaxis() -> SetTitle(sJetN.Data());
    hJetN[iRes][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetN[iRes][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetN[iRes][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetN[iRes][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetN[iRes][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetN[iRes][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetN[iRes][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetN[iRes][0] -> GetYaxis() -> SetTitle(sJetNY.Data());
    hJetN[iRes][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetN[iRes][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetN[iRes][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetN[iRes][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetN[iRes][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetN[iRes][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetN[iRes][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetN[iRes][1] -> SetLineColor(fColTrg[1]);
    hJetN[iRes][1] -> SetLineStyle(fLinAll);
    hJetN[iRes][1] -> SetFillColor(fColTrg[1]);
    hJetN[iRes][1] -> SetFillStyle(fFilAll);
    hJetN[iRes][1] -> SetMarkerColor(fColTrg[1]);
    hJetN[iRes][1] -> SetMarkerStyle(fMarAll);
    hJetN[iRes][1] -> SetTitle("");
    hJetN[iRes][1] -> SetTitleFont(fTxt);
    hJetN[iRes][1] -> GetXaxis() -> SetTitle(sJetN.Data());
    hJetN[iRes][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetN[iRes][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetN[iRes][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetN[iRes][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetN[iRes][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetN[iRes][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetN[iRes][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetN[iRes][1] -> GetYaxis() -> SetTitle(sJetNY.Data());
    hJetN[iRes][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetN[iRes][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetN[iRes][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetN[iRes][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetN[iRes][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetN[iRes][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetN[iRes][1] -> GetYaxis() -> CenterTitle(fCnt);

    // jet area
    hJetA[iRes][0]     -> SetLineColor(fColAll);
    hJetA[iRes][0]     -> SetLineStyle(fLinAll);
    hJetA[iRes][0]     -> SetFillColor(fColAll);
    hJetA[iRes][0]     -> SetFillStyle(fFilAll);
    hJetA[iRes][0]     -> SetMarkerColor(fColAll);
    hJetA[iRes][0]     -> SetMarkerStyle(fMarAll);
    hJetA[iRes][0]     -> SetTitle("");
    hJetA[iRes][0]     -> SetTitleFont(fTxt);
    hJetA[iRes][0]     -> GetXaxis() -> SetTitle(sJetA.Data());
    hJetA[iRes][0]     -> GetXaxis() -> SetTitleFont(fTxt);
    hJetA[iRes][0]     -> GetXaxis() -> SetTitleSize(fTit);
    hJetA[iRes][0]     -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetA[iRes][0]     -> GetXaxis() -> SetLabelFont(fTxt);
    hJetA[iRes][0]     -> GetXaxis() -> SetLabelSize(fLab);
    hJetA[iRes][0]     -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetA[iRes][0]     -> GetXaxis() -> CenterTitle(fCnt);
    hJetA[iRes][0]     -> GetXaxis() -> SetRangeUser(xAreaPlot[iRes][0], xAreaPlot[iRes][1]);
    hJetA[iRes][0]     -> GetYaxis() -> SetTitle(sCount.Data());
    hJetA[iRes][0]     -> GetYaxis() -> SetTitleFont(fTxt);
    hJetA[iRes][0]     -> GetYaxis() -> SetTitleSize(fTit);
    hJetA[iRes][0]     -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetA[iRes][0]     -> GetYaxis() -> SetLabelFont(fTxt);
    hJetA[iRes][0]     -> GetYaxis() -> SetLabelSize(fLab);
    hJetA[iRes][0]     -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetA[iRes][0]     -> GetYaxis() -> CenterTitle(fCnt);
    hJetA[iRes][1]     -> SetLineColor(fColCut);
    hJetA[iRes][1]     -> SetLineStyle(fLinCut);
    hJetA[iRes][1]     -> SetFillColor(fColCut);
    hJetA[iRes][1]     -> SetFillStyle(fFilCut);
    hJetA[iRes][1]     -> SetMarkerColor(fColCut);
    hJetA[iRes][1]     -> SetMarkerStyle(fMarCut);
    hJetA[iRes][1]     -> SetTitle("");
    hJetA[iRes][1]     -> SetTitleFont(fTxt);
    hJetA[iRes][1]     -> GetXaxis() -> SetTitle(sJetA.Data());
    hJetA[iRes][1]     -> GetXaxis() -> SetTitleFont(fTxt);
    hJetA[iRes][1]     -> GetXaxis() -> SetTitleSize(fTit);
    hJetA[iRes][1]     -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetA[iRes][1]     -> GetXaxis() -> SetLabelFont(fTxt);
    hJetA[iRes][1]     -> GetXaxis() -> SetLabelSize(fLab);
    hJetA[iRes][1]     -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetA[iRes][1]     -> GetXaxis() -> CenterTitle(fCnt);
    hJetA[iRes][1]     -> GetXaxis() -> SetRangeUser(xAreaPlot[iRes][0], xAreaPlot[iRes][1]);
    hJetA[iRes][1]     -> GetYaxis() -> SetTitle(sCount.Data());
    hJetA[iRes][1]     -> GetYaxis() -> SetTitleFont(fTxt);
    hJetA[iRes][1]     -> GetYaxis() -> SetTitleSize(fTit);
    hJetA[iRes][1]     -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetA[iRes][1]     -> GetYaxis() -> SetLabelFont(fTxt);
    hJetA[iRes][1]     -> GetYaxis() -> SetLabelSize(fLab);
    hJetA[iRes][1]     -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetA[iRes][1]     -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][0] -> SetLineColor(fColAll);
    hJetPtVsA[iRes][0] -> SetLineStyle(fLinAll);
    hJetPtVsA[iRes][0] -> SetFillColor(fColAll);
    hJetPtVsA[iRes][0] -> SetFillStyle(fFilAll);
    hJetPtVsA[iRes][0] -> SetMarkerColor(fColAll);
    hJetPtVsA[iRes][0] -> SetMarkerStyle(fMarAll);
    hJetPtVsA[iRes][0] -> SetTitle("");
    hJetPtVsA[iRes][0] -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetTitle(sJetA.Data());
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsA[iRes][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][0] -> GetXaxis() -> SetRangeUser(xAreaPlot[iRes][0], xAreaPlot[iRes][1]);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetTitle(sJetPt.Data());
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsA[iRes][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][0] -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetTitle("");
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsA[iRes][0] -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][0] -> GetZaxis() -> SetRangeUser(zAreaVsPt[0], zAreaVsPt[1]);
    hJetPtVsA[iRes][1] -> SetLineColor(fColCut);
    hJetPtVsA[iRes][1] -> SetLineStyle(fLinCut);
    hJetPtVsA[iRes][1] -> SetFillColor(fColCut);
    hJetPtVsA[iRes][1] -> SetFillStyle(fFilCut);
    hJetPtVsA[iRes][1] -> SetMarkerColor(fColCut);
    hJetPtVsA[iRes][1] -> SetMarkerStyle(fMarCut);
    hJetPtVsA[iRes][1] -> SetTitle("");
    hJetPtVsA[iRes][1] -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetTitle(sJetA.Data());
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsA[iRes][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][1] -> GetXaxis() -> SetRangeUser(xAreaPlot[iRes][0], xAreaPlot[iRes][1]);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetTitle(sJetPt.Data());
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsA[iRes][1] -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][1] -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetTitle("");
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsA[iRes][1] -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsA[iRes][1] -> GetZaxis() -> SetRangeUser(zAreaVsPt[0], zAreaVsPt[1]);

    // jet eta
    hJetEta[iRes][0]       -> SetLineColor(fColAll);
    hJetEta[iRes][0]       -> SetLineStyle(fLinAll);
    hJetEta[iRes][0]       -> SetFillColor(fColAll);
    hJetEta[iRes][0]       -> SetFillStyle(fFilAll);
    hJetEta[iRes][0]       -> SetMarkerColor(fColAll);
    hJetEta[iRes][0]       -> SetMarkerStyle(fMarAll);
    hJetEta[iRes][0]       -> SetTitle("");
    hJetEta[iRes][0]       -> SetTitleFont(fTxt);
    hJetEta[iRes][0]       -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetEta[iRes][0]       -> GetXaxis() -> SetTitleFont(fTxt);
    hJetEta[iRes][0]       -> GetXaxis() -> SetTitleSize(fTit);
    hJetEta[iRes][0]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetEta[iRes][0]       -> GetXaxis() -> SetLabelFont(fTxt);
    hJetEta[iRes][0]       -> GetXaxis() -> SetLabelSize(fLab);
    hJetEta[iRes][0]       -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetEta[iRes][0]       -> GetXaxis() -> CenterTitle(fCnt);
    hJetEta[iRes][0]       -> GetYaxis() -> SetTitle(sCount.Data());
    hJetEta[iRes][0]       -> GetYaxis() -> SetTitleFont(fTxt);
    hJetEta[iRes][0]       -> GetYaxis() -> SetTitleSize(fTit);
    hJetEta[iRes][0]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetEta[iRes][0]       -> GetYaxis() -> SetLabelFont(fTxt);
    hJetEta[iRes][0]       -> GetYaxis() -> SetLabelSize(fLab);
    hJetEta[iRes][0]       -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetEta[iRes][0]       -> GetYaxis() -> CenterTitle(fCnt);
    hJetEta[iRes][1]       -> SetLineColor(fColCut);
    hJetEta[iRes][1]       -> SetLineStyle(fLinCut);
    hJetEta[iRes][1]       -> SetFillColor(fColCut);
    hJetEta[iRes][1]       -> SetFillStyle(fFilCut);
    hJetEta[iRes][1]       -> SetMarkerColor(fColCut);
    hJetEta[iRes][1]       -> SetMarkerStyle(fMarCut);
    hJetEta[iRes][1]       -> SetTitle("");
    hJetEta[iRes][1]       -> SetTitleFont(fTxt);
    hJetEta[iRes][1]       -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetEta[iRes][1]       -> GetXaxis() -> SetTitleFont(fTxt);
    hJetEta[iRes][1]       -> GetXaxis() -> SetTitleSize(fTit);
    hJetEta[iRes][1]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetEta[iRes][1]       -> GetXaxis() -> SetLabelFont(fTxt);
    hJetEta[iRes][1]       -> GetXaxis() -> SetLabelSize(fLab);
    hJetEta[iRes][1]       -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetEta[iRes][1]       -> GetXaxis() -> CenterTitle(fCnt);
    hJetEta[iRes][1]       -> GetYaxis() -> SetTitle(sCount.Data());
    hJetEta[iRes][1]       -> GetYaxis() -> SetTitleFont(fTxt);
    hJetEta[iRes][1]       -> GetYaxis() -> SetTitleSize(fTit);
    hJetEta[iRes][1]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetEta[iRes][1]       -> GetYaxis() -> SetLabelFont(fTxt);
    hJetEta[iRes][1]       -> GetYaxis() -> SetLabelSize(fLab);
    hJetEta[iRes][1]       -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetEta[iRes][1]       -> GetYaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][0][0] -> SetLineColor(fColTrg[0]);
    hJetEtaBin[iRes][0][0] -> SetLineStyle(fLinAll);
    hJetEtaBin[iRes][0][0] -> SetFillColor(fColTrg[0]);
    hJetEtaBin[iRes][0][0] -> SetFillStyle(fFilAll);
    hJetEtaBin[iRes][0][0] -> SetMarkerColor(fColTrg[0]);
    hJetEtaBin[iRes][0][0] -> SetMarkerStyle(fMarAll);
    hJetEtaBin[iRes][0][0] -> SetTitle("");
    hJetEtaBin[iRes][0][0] -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][0][0] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][0][0] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
    hJetEtaBin[iRes][1][0] -> SetLineColor(fColTrg[1]);
    hJetEtaBin[iRes][1][0] -> SetLineStyle(fLinAll);
    hJetEtaBin[iRes][1][0] -> SetFillColor(fColTrg[1]);
    hJetEtaBin[iRes][1][0] -> SetFillStyle(fFilAll);
    hJetEtaBin[iRes][1][0] -> SetMarkerColor(fColTrg[1]);
    hJetEtaBin[iRes][1][0] -> SetMarkerStyle(fMarAll);
    hJetEtaBin[iRes][1][0] -> SetTitle("");
    hJetEtaBin[iRes][1][0] -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][1][0] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][1][0] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
    hJetEtaBin[iRes][0][1] -> SetLineColor(fColTrg[0]);
    hJetEtaBin[iRes][0][1] -> SetLineStyle(fLinAll);
    hJetEtaBin[iRes][0][1] -> SetFillColor(fColTrg[0]);
    hJetEtaBin[iRes][0][1] -> SetFillStyle(fFilAll);
    hJetEtaBin[iRes][0][1] -> SetMarkerColor(fColTrg[0]);
    hJetEtaBin[iRes][0][1] -> SetMarkerStyle(fMarAll);
    hJetEtaBin[iRes][0][1] -> SetTitle("");
    hJetEtaBin[iRes][0][1] -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][0][1] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][0][1] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
    hJetEtaBin[iRes][1][1] -> SetLineColor(fColTrg[1]);
    hJetEtaBin[iRes][1][1] -> SetLineStyle(fLinAll);
    hJetEtaBin[iRes][1][1] -> SetFillColor(fColTrg[1]);
    hJetEtaBin[iRes][1][1] -> SetFillStyle(fFilAll);
    hJetEtaBin[iRes][1][1] -> SetMarkerColor(fColTrg[1]);
    hJetEtaBin[iRes][1][1] -> SetMarkerStyle(fMarAll);
    hJetEtaBin[iRes][1][1] -> SetTitle("");
    hJetEtaBin[iRes][1][1] -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][1][1] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> CenterTitle(fCnt);
    hJetEtaBin[iRes][1][1] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
    hJetPtVsEta[iRes][0]   -> SetLineColor(fColAll);
    hJetPtVsEta[iRes][0]   -> SetLineStyle(fLinAll);
    hJetPtVsEta[iRes][0]   -> SetFillColor(fColAll);
    hJetPtVsEta[iRes][0]   -> SetFillStyle(fFilAll);
    hJetPtVsEta[iRes][0]   -> SetMarkerColor(fColAll);
    hJetPtVsEta[iRes][0]   -> SetMarkerStyle(fMarAll);
    hJetPtVsEta[iRes][0]   -> SetTitle("");
    hJetPtVsEta[iRes][0]   -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsEta[iRes][0]   -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetTitle(sPtReco.Data());
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsEta[iRes][0]   -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetTitle("");
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsEta[iRes][0]   -> GetZaxis() -> SetRangeUser(zEtaVsPt[iRes][0], zEtaVsPt[iRes][1]);
    hJetPtVsEta[iRes][1]   -> SetLineColor(fColCut);
    hJetPtVsEta[iRes][1]   -> SetLineStyle(fLinCut);
    hJetPtVsEta[iRes][1]   -> SetFillColor(fColCut);
    hJetPtVsEta[iRes][1]   -> SetFillStyle(fFilCut);
    hJetPtVsEta[iRes][1]   -> SetMarkerColor(fColCut);
    hJetPtVsEta[iRes][1]   -> SetMarkerStyle(fMarCut);
    hJetPtVsEta[iRes][1]   -> SetTitle("");
    hJetPtVsEta[iRes][1]   -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsEta[iRes][1]   -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetTitle(sPtReco.Data());
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsEta[iRes][1]   -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetTitle("");
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsEta[iRes][1]   -> GetZaxis() -> SetRangeUser(zEtaVsPt[iRes][0], zEtaVsPt[iRes][1]);
    hJetDfVsEta[iRes][0]   -> SetLineColor(fColAll);
    hJetDfVsEta[iRes][0]   -> SetLineStyle(fLinAll);
    hJetDfVsEta[iRes][0]   -> SetFillColor(fColAll);
    hJetDfVsEta[iRes][0]   -> SetFillStyle(fFilAll);
    hJetDfVsEta[iRes][0]   -> SetMarkerColor(fColAll);
    hJetDfVsEta[iRes][0]   -> SetMarkerStyle(fMarAll);
    hJetDfVsEta[iRes][0]   -> SetTitle("");
    hJetDfVsEta[iRes][0]   -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDfVsEta[iRes][0]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetTitle(sJetDf.Data());
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDfVsEta[iRes][0]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> SetTitle("");
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> SetTitleSize(fTit);
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> SetLabelFont(fTxt);
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> SetLabelSize(fLab);
    hJetDfVsEta[iRes][0]   -> GetZaxis() -> CenterTitle(fCnt);
    hJetDfVsEta[iRes][1]   -> SetLineColor(fColCut);
    hJetDfVsEta[iRes][1]   -> SetLineStyle(fLinCut);
    hJetDfVsEta[iRes][1]   -> SetFillColor(fColCut);
    hJetDfVsEta[iRes][1]   -> SetFillStyle(fFilCut);
    hJetDfVsEta[iRes][1]   -> SetMarkerColor(fColCut);
    hJetDfVsEta[iRes][1]   -> SetMarkerStyle(fMarCut);
    hJetDfVsEta[iRes][1]   -> SetTitle("");
    hJetDfVsEta[iRes][1]   -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetTitle(sJetEta.Data());
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDfVsEta[iRes][1]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetTitle(sJetDf.Data());
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDfVsEta[iRes][1]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> SetTitle("");
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> SetTitleFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> SetTitleSize(fTit);
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> SetLabelFont(fTxt);
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> SetLabelSize(fLab);
    hJetDfVsEta[iRes][1]   -> GetZaxis() -> CenterTitle(fCnt);

    // pT-differential eta histograms
    for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetLineColor(fColAllEtaPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetLineStyle(fLinAll);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetFillColor(fColAllEtaPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetFillStyle(fFilAll);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetMarkerColor(fColAllEtaPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetMarkerStyle(fMarAllPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetTitle("");
      hJetEtaPtBin[iRes][0][iJetBin][0] -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitle(sJetEta.Data());
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetLineColor(fColAllEtaGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetLineStyle(fLinAll);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetFillColor(fColAllEtaGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetFillStyle(fFilAll);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetMarkerColor(fColAllEtaGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetMarkerStyle(fMarAllGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetTitle("");
      hJetEtaPtBin[iRes][1][iJetBin][0] -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitle(sJetEta.Data());
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetLineColor(fColRecEtaPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetLineStyle(fLinAll);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetFillColor(fColRecEtaPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetFillStyle(fFilAll);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetMarkerColor(fColRecEtaPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetMarkerStyle(fMarRecPi0[iJetBin]);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetTitle("");
      hJetEtaPtBin[iRes][0][iJetBin][1] -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitle(sJetEta.Data());
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetLineColor(fColRecEtaGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetLineStyle(fLinAll);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetFillColor(fColRecEtaGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetFillStyle(fFilAll);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetMarkerColor(fColRecEtaGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetMarkerStyle(fMarRecGam[iJetBin]);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetTitle("");
      hJetEtaPtBin[iRes][1][iJetBin][1] -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitle(sJetEta.Data());
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetRangeUser(xEtaPlot[0], xEtaPlot[1]);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitle(sJetEtaY.Data());
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitleSize(fTit);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetLabelSize(fLab);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> CenterTitle(fCnt);
      hJetEtaPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetRangeUser(yEtaPlot[0], yEtaPlot[1]);
    }  // end jet pT bin loop

    // jet dF
    hJetDfBin[iRes][0][0]   -> SetLineColor(fColDfAll[0]);
    hJetDfBin[iRes][0][0]   -> SetLineStyle(fLinAll);
    hJetDfBin[iRes][0][0]   -> SetFillColor(fColDfAll[0]);
    hJetDfBin[iRes][0][0]   -> SetFillStyle(fFilAll);
    hJetDfBin[iRes][0][0]   -> SetMarkerColor(fColDfAll[0]);
    hJetDfBin[iRes][0][0]   -> SetMarkerStyle(fMarDfAll[0]);
    hJetDfBin[iRes][0][0]   -> SetTitle("");
    hJetDfBin[iRes][0][0]   -> SetTitleFont(fTxt);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][0][0]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetTitle(sJetDfY.Data());
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][0][0]   -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
    hJetDfBin[iRes][1][0]   -> SetLineColor(fColDfAll[1]);
    hJetDfBin[iRes][1][0]   -> SetLineStyle(fLinAll);
    hJetDfBin[iRes][1][0]   -> SetFillColor(fColDfAll[1]);
    hJetDfBin[iRes][1][0]   -> SetFillStyle(fFilAll);
    hJetDfBin[iRes][1][0]   -> SetMarkerColor(fColDfAll[1]);
    hJetDfBin[iRes][1][0]   -> SetMarkerStyle(fMarDfAll[1]);
    hJetDfBin[iRes][1][0]   -> SetTitle("");
    hJetDfBin[iRes][1][0]   -> SetTitleFont(fTxt);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][1][0]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetTitle(sJetDfY.Data());
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][1][0]   -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
    hJetDfBin[iRes][0][1]   -> SetLineColor(fColDfRec[0]);
    hJetDfBin[iRes][0][1]   -> SetLineStyle(fLinAll);
    hJetDfBin[iRes][0][1]   -> SetFillColor(fColDfRec[0]);
    hJetDfBin[iRes][0][1]   -> SetFillStyle(fFilAll);
    hJetDfBin[iRes][0][1]   -> SetMarkerColor(fColDfRec[0]);
    hJetDfBin[iRes][0][1]   -> SetMarkerStyle(fMarDfRec[0]);
    hJetDfBin[iRes][0][1]   -> SetTitle("");
    hJetDfBin[iRes][0][1]   -> SetTitleFont(fTxt);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][0][1]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetTitle(sJetDfY.Data());
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][0][1]   -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
    hJetDfBin[iRes][1][1]   -> SetLineColor(fColDfRec[1]);
    hJetDfBin[iRes][1][1]   -> SetLineStyle(fLinAll);
    hJetDfBin[iRes][1][1]   -> SetFillColor(fColDfRec[1]);
    hJetDfBin[iRes][1][1]   -> SetFillStyle(fFilAll);
    hJetDfBin[iRes][1][1]   -> SetMarkerColor(fColDfRec[1]);
    hJetDfBin[iRes][1][1]   -> SetMarkerStyle(fMarDfRec[1]);
    hJetDfBin[iRes][1][1]   -> SetTitle("");
    hJetDfBin[iRes][1][1]   -> SetTitleFont(fTxt);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][1][1]   -> GetXaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetTitle(sJetDfY.Data());
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetTitleSize(fTit);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetLabelSize(fLab);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> CenterTitle(fCnt);
    hJetDfBin[iRes][1][1]   -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
    hJetPtVsDf[iRes][0][0]  -> SetTitle("");
    hJetPtVsDf[iRes][0][0]  -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][0][0]  -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetTitle(sPtReco.Data());
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][0][0]  -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetTitle("");
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][0][0]  -> GetZaxis() -> SetRangeUser(zDfVsPt[iRes][0], zDfVsPt[iRes][1]);
    hJetPtVsDf[iRes][1][0]  -> SetTitle("");
    hJetPtVsDf[iRes][1][0]  -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][1][0]  -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetTitle(sPtReco.Data());
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][1][0]  -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetTitle("");
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][1][0]  -> GetZaxis() -> SetRangeUser(zDfVsPt[iRes][0], zDfVsPt[iRes][1]);
    hJetPtVsDf[iRes][0][1]  -> SetTitle("");
    hJetPtVsDf[iRes][0][1]  -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][0][1]  -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetTitle(sPtReco.Data());
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][0][1]  -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetTitle("");
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][0][1]  -> GetZaxis() -> SetRangeUser(zDfVsPt[iRes][0], zDfVsPt[iRes][1]);
    hJetPtVsDf[iRes][1][1]  -> SetTitle("");
    hJetPtVsDf[iRes][1][1]  -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetTitle(sJetDf.Data());
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][1][1]  -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetTitle(sPtReco.Data());
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][1][1]  -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetTitle("");
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsDf[iRes][1][1]  -> GetZaxis() -> SetRangeUser(zDfVsPt[iRes][0], zDfVsPt[iRes][1]);

    // pT-differential dF histograms
    for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetLineColor(fColAllDfPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetLineStyle(fLinAll);
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetFillColor(fColAllDfPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetFillStyle(fFilAll);
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetMarkerColor(fColAllDfPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetMarkerStyle(fMarAllPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetTitle("");
      hJetDfPtBin[iRes][0][iJetBin][0] -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitle(sJetDf.Data());
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetXaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitle(sJetDfY.Data());
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][0][iJetBin][0] -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetLineColor(fColAllDfGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetLineStyle(fLinAll);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetFillColor(fColAllDfGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetFillStyle(fFilAll);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetMarkerColor(fColAllDfGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetMarkerStyle(fMarAllGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetTitle("");
      hJetDfPtBin[iRes][1][iJetBin][0] -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitle(sJetDf.Data());
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetXaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitle(sJetDfY.Data());
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][1][iJetBin][0] -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetLineColor(fColRecDfPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetLineStyle(fLinAll);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetFillColor(fColRecDfPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetFillStyle(fFilAll);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetMarkerColor(fColRecDfPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetMarkerStyle(fMarRecPi0[iJetBin]);
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetTitle("");
      hJetDfPtBin[iRes][0][iJetBin][1] -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitle(sJetDf.Data());
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetXaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitle(sJetDfY.Data());
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][0][iJetBin][1] -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetLineColor(fColRecDfGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetLineStyle(fLinAll);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetFillColor(fColRecDfGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetFillStyle(fFilAll);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetMarkerColor(fColRecDfGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetMarkerStyle(fMarRecGam[iJetBin]);
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetTitle("");
      hJetDfPtBin[iRes][1][iJetBin][1] -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitle(sJetDf.Data());
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetTitleOffset(fOffX);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetXaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitle(sJetDfY.Data());
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitleFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitleSize(fTit);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetTitleOffset(fOffY);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetLabelFont(fTxt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetLabelSize(fLab);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetLabelOffset(fOffL);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> CenterTitle(fCnt);
      hJetDfPtBin[iRes][1][iJetBin][1] -> GetYaxis() -> SetRangeUser(yDfPlot[0], yDfPlot[1]);
    }  // end jet pT bin loop

    // jet pT and zT
    hJetPt[iRes][0]    -> SetLineColor(fColAll);
    hJetPt[iRes][0]    -> SetLineStyle(fLinAll);
    hJetPt[iRes][0]    -> SetFillColor(fColAll);
    hJetPt[iRes][0]    -> SetFillStyle(fFilAll);
    hJetPt[iRes][0]    -> SetMarkerColor(fColAll);
    hJetPt[iRes][0]    -> SetMarkerStyle(fMarAll);
    hJetPt[iRes][0]    -> SetTitle("");
    hJetPt[iRes][0]    -> SetTitleFont(fTxt);
    hJetPt[iRes][0]    -> GetXaxis() -> SetTitle(sJetPt.Data());
    hJetPt[iRes][0]    -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPt[iRes][0]    -> GetXaxis() -> SetTitleSize(fTit);
    hJetPt[iRes][0]    -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPt[iRes][0]    -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPt[iRes][0]    -> GetXaxis() -> SetLabelSize(fLab);
    hJetPt[iRes][0]    -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPt[iRes][0]    -> GetXaxis() -> CenterTitle(fCnt);
    hJetPt[iRes][0]    -> GetYaxis() -> SetTitle(sCount.Data());
    hJetPt[iRes][0]    -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPt[iRes][0]    -> GetYaxis() -> SetTitleSize(fTit);
    hJetPt[iRes][0]    -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPt[iRes][0]    -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPt[iRes][0]    -> GetYaxis() -> SetLabelSize(fLab);
    hJetPt[iRes][0]    -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPt[iRes][0]    -> GetYaxis() -> CenterTitle(fCnt);
    hJetPt[iRes][1]    -> SetLineColor(fColCut);
    hJetPt[iRes][1]    -> SetLineStyle(fLinCut);
    hJetPt[iRes][1]    -> SetFillColor(fColCut);
    hJetPt[iRes][1]    -> SetFillStyle(fFilCut);
    hJetPt[iRes][1]    -> SetMarkerColor(fColCut);
    hJetPt[iRes][1]    -> SetMarkerStyle(fMarCut);
    hJetPt[iRes][1]    -> SetTitle("");
    hJetPt[iRes][1]    -> SetTitleFont(fTxt);
    hJetPt[iRes][1]    -> GetXaxis() -> SetTitle(sJetPt.Data());
    hJetPt[iRes][1]    -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPt[iRes][1]    -> GetXaxis() -> SetTitleSize(fTit);
    hJetPt[iRes][1]    -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPt[iRes][1]    -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPt[iRes][1]    -> GetXaxis() -> SetLabelSize(fLab);
    hJetPt[iRes][1]    -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPt[iRes][1]    -> GetXaxis() -> CenterTitle(fCnt);
    hJetPt[iRes][1]    -> GetYaxis() -> SetTitle(sCount.Data());
    hJetPt[iRes][1]    -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPt[iRes][1]    -> GetYaxis() -> SetTitleSize(fTit);
    hJetPt[iRes][1]    -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPt[iRes][1]    -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPt[iRes][1]    -> GetYaxis() -> SetLabelSize(fLab);
    hJetPt[iRes][1]    -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPt[iRes][1]    -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtBin[iRes][0] -> SetLineColor(fColTrg[0]);
    hJetPtBin[iRes][0] -> SetLineStyle(fLinAll);
    hJetPtBin[iRes][0] -> SetFillColor(fColTrg[0]);
    hJetPtBin[iRes][0] -> SetFillStyle(fFilAll);
    hJetPtBin[iRes][0] -> SetMarkerColor(fColTrg[0]);
    hJetPtBin[iRes][0] -> SetMarkerStyle(fMarAll);
    hJetPtBin[iRes][0] -> SetTitle("");
    hJetPtBin[iRes][0] -> SetTitleFont(fTxt);
    hJetPtBin[iRes][0] -> GetXaxis() -> SetTitle(sPtReco.Data());
    hJetPtBin[iRes][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtBin[iRes][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtBin[iRes][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtBin[iRes][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtBin[iRes][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtBin[iRes][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtBin[iRes][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtBin[iRes][0] -> GetYaxis() -> SetTitle(sPtRecoY.Data());
    hJetPtBin[iRes][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtBin[iRes][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtBin[iRes][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtBin[iRes][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtBin[iRes][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtBin[iRes][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtBin[iRes][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtBin[iRes][1] -> SetLineColor(fColTrg[1]);
    hJetPtBin[iRes][1] -> SetLineStyle(fLinAll);
    hJetPtBin[iRes][1] -> SetFillColor(fColTrg[1]);
    hJetPtBin[iRes][1] -> SetFillStyle(fFilAll);
    hJetPtBin[iRes][1] -> SetMarkerColor(fColTrg[1]);
    hJetPtBin[iRes][1] -> SetMarkerStyle(fMarAll);
    hJetPtBin[iRes][1] -> SetTitle("");
    hJetPtBin[iRes][1] -> SetTitleFont(fTxt);
    hJetPtBin[iRes][1] -> GetXaxis() -> SetTitle(sPtReco.Data());
    hJetPtBin[iRes][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtBin[iRes][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtBin[iRes][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtBin[iRes][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtBin[iRes][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtBin[iRes][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtBin[iRes][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtBin[iRes][1] -> GetYaxis() -> SetTitle(sPtRecoY.Data());
    hJetPtBin[iRes][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtBin[iRes][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtBin[iRes][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtBin[iRes][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtBin[iRes][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtBin[iRes][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtBin[iRes][1] -> GetYaxis() -> CenterTitle(fCnt);
    hJetZtBin[iRes][0] -> SetLineColor(fColTrg[0]);
    hJetZtBin[iRes][0] -> SetLineStyle(fLinAll);
    hJetZtBin[iRes][0] -> SetFillColor(fColTrg[0]);
    hJetZtBin[iRes][0] -> SetFillStyle(fFilAll);
    hJetZtBin[iRes][0] -> SetMarkerColor(fColTrg[0]);
    hJetZtBin[iRes][0] -> SetMarkerStyle(fMarAll);
    hJetZtBin[iRes][0] -> SetTitle("");
    hJetZtBin[iRes][0] -> SetTitleFont(fTxt);
    hJetZtBin[iRes][0] -> GetXaxis() -> SetTitle(sJetZt.Data());
    hJetZtBin[iRes][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetZtBin[iRes][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetZtBin[iRes][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetZtBin[iRes][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetZtBin[iRes][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetZtBin[iRes][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetZtBin[iRes][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetZtBin[iRes][0] -> GetYaxis() -> SetTitle(sJetZtY.Data());
    hJetZtBin[iRes][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetZtBin[iRes][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetZtBin[iRes][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetZtBin[iRes][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetZtBin[iRes][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetZtBin[iRes][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetZtBin[iRes][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetZtBin[iRes][1] -> SetLineColor(fColTrg[1]);
    hJetZtBin[iRes][1] -> SetLineStyle(fLinAll);
    hJetZtBin[iRes][1] -> SetFillColor(fColTrg[1]);
    hJetZtBin[iRes][1] -> SetFillStyle(fFilAll);
    hJetZtBin[iRes][1] -> SetMarkerColor(fColTrg[1]);
    hJetZtBin[iRes][1] -> SetMarkerStyle(fMarAll);
    hJetZtBin[iRes][1] -> SetTitle("");
    hJetZtBin[iRes][1] -> SetTitleFont(fTxt);
    hJetZtBin[iRes][1] -> GetXaxis() -> SetTitle(sJetZt.Data());
    hJetZtBin[iRes][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetZtBin[iRes][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetZtBin[iRes][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetZtBin[iRes][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetZtBin[iRes][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetZtBin[iRes][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetZtBin[iRes][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetZtBin[iRes][1] -> GetYaxis() -> SetTitle(sJetZtY.Data());
    hJetZtBin[iRes][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetZtBin[iRes][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetZtBin[iRes][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetZtBin[iRes][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetZtBin[iRes][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetZtBin[iRes][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetZtBin[iRes][1] -> GetYaxis() -> CenterTitle(fCnt);

    // jet pT correction
    hJetCorr[iRes][0]     -> SetLineColor(fColTrg[0]);
    hJetCorr[iRes][0]     -> SetLineStyle(fLinAll);
    hJetCorr[iRes][0]     -> SetFillColor(fColTrg[0]);
    hJetCorr[iRes][0]     -> SetFillStyle(fFilAll);
    hJetCorr[iRes][0]     -> SetMarkerColor(fColTrg[0]);
    hJetCorr[iRes][0]     -> SetMarkerStyle(fMarAll);
    hJetCorr[iRes][0]     -> SetTitle("");
    hJetCorr[iRes][0]     -> SetTitleFont(fTxt);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetTitle(sPtCorr.Data());
    hJetCorr[iRes][0]     -> GetXaxis() -> SetTitleFont(fTxt);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetTitleSize(fTit);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetLabelFont(fTxt);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetLabelSize(fLab);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetCorr[iRes][0]     -> GetXaxis() -> CenterTitle(fCnt);
    hJetCorr[iRes][0]     -> GetXaxis() -> SetRangeUser(xCorrPlot[0], xCorrPlot[1]);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetTitle(sPtCorrY.Data());
    hJetCorr[iRes][0]     -> GetYaxis() -> SetTitleFont(fTxt);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetTitleSize(fTit);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetLabelFont(fTxt);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetLabelSize(fLab);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetCorr[iRes][0]     -> GetYaxis() -> CenterTitle(fCnt);
    hJetCorr[iRes][0]     -> GetYaxis() -> SetRangeUser(yCorrPlot[iRes][0], yCorrPlot[iRes][1]);
    hJetCorr[iRes][1]     -> SetLineColor(fColTrg[1]);
    hJetCorr[iRes][1]     -> SetLineStyle(fLinAll);
    hJetCorr[iRes][1]     -> SetFillColor(fColTrg[1]);
    hJetCorr[iRes][1]     -> SetFillStyle(fFilAll);
    hJetCorr[iRes][1]     -> SetMarkerColor(fColTrg[1]);
    hJetCorr[iRes][1]     -> SetMarkerStyle(fMarAll);
    hJetCorr[iRes][1]     -> SetTitle("");
    hJetCorr[iRes][1]     -> SetTitleFont(fTxt);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetTitle(sPtCorr.Data());
    hJetCorr[iRes][1]     -> GetXaxis() -> SetTitleFont(fTxt);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetTitleSize(fTit);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetLabelFont(fTxt);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetLabelSize(fLab);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetCorr[iRes][1]     -> GetXaxis() -> CenterTitle(fCnt);
    hJetCorr[iRes][1]     -> GetXaxis() -> SetRangeUser(xCorrPlot[0], xCorrPlot[1]);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetTitle(sPtCorrY.Data());
    hJetCorr[iRes][1]     -> GetYaxis() -> SetTitleFont(fTxt);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetTitleSize(fTit);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetLabelFont(fTxt);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetLabelSize(fLab);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetCorr[iRes][1]     -> GetYaxis() -> CenterTitle(fCnt);
    hJetCorr[iRes][1]     -> GetYaxis() -> SetRangeUser(yCorrPlot[iRes][0], yCorrPlot[iRes][1]);
    hJetPtVsCorr[iRes][0] -> SetLineColor(fColAll);
    hJetPtVsCorr[iRes][0] -> SetLineStyle(fLinAll);
    hJetPtVsCorr[iRes][0] -> SetFillColor(fColAll);
    hJetPtVsCorr[iRes][0] -> SetFillStyle(fFilAll);
    hJetPtVsCorr[iRes][0] -> SetMarkerColor(fColAll);
    hJetPtVsCorr[iRes][0] -> SetMarkerStyle(fMarAll);
    hJetPtVsCorr[iRes][0] -> SetTitle("");
    hJetPtVsCorr[iRes][0] -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetTitle(sPtCorr.Data());
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsCorr[iRes][0] -> GetXaxis() -> SetRangeUser(xCorrPlot[0], xCorrPlot[1]);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetTitle(sJetPt.Data());
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsCorr[iRes][0] -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetTitle("");
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsCorr[iRes][0] -> GetZaxis() -> SetRangeUser(zCorrVsPt[iRes][0], zCorrVsPt[iRes][1]);
    hJetPtVsCorr[iRes][1] -> SetLineColor(fColCut);
    hJetPtVsCorr[iRes][1] -> SetLineStyle(fLinCut);
    hJetPtVsCorr[iRes][1] -> SetFillColor(fColCut);
    hJetPtVsCorr[iRes][1] -> SetFillStyle(fFilCut);
    hJetPtVsCorr[iRes][1] -> SetMarkerColor(fColCut);
    hJetPtVsCorr[iRes][1] -> SetMarkerStyle(fMarCut);
    hJetPtVsCorr[iRes][1] -> SetTitle("");
    hJetPtVsCorr[iRes][1] -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetTitle(sPtCorr.Data());
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetPtVsCorr[iRes][1] -> GetXaxis() -> SetRangeUser(xCorrPlot[0], xCorrPlot[1]);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetTitle(sJetPt.Data());
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> CenterTitle(fCnt);
    hJetPtVsCorr[iRes][1] -> GetYaxis() -> SetRangeUser(yPtPlot[0], yPtPlot[1]);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetTitle("");
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetTitleFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetTitleSize(fTit);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetTitleOffset(fOffZ);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetLabelFont(fTxt);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetLabelSize(fLab);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> CenterTitle(fCnt);
    hJetPtVsCorr[iRes][1] -> GetZaxis() -> SetRangeUser(zCorrVsPt[iRes][0], zCorrVsPt[iRes][1]);

    // jet sigma
    hJetSig[iRes][0] -> SetLineColor(fColTrg[0]);
    hJetSig[iRes][0] -> SetLineStyle(fLinAll);
    hJetSig[iRes][0] -> SetFillColor(fColTrg[0]);
    hJetSig[iRes][0] -> SetFillStyle(fFilAll);
    hJetSig[iRes][0] -> SetMarkerColor(fColTrg[0]);
    hJetSig[iRes][0] -> SetMarkerStyle(fMarAll);
    hJetSig[iRes][0] -> SetTitle("");
    hJetSig[iRes][0] -> SetTitleFont(fTxt);
    hJetSig[iRes][0] -> GetXaxis() -> SetTitle(sJetSig.Data());
    hJetSig[iRes][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetSig[iRes][0] -> GetXaxis() -> SetTitleSize(fTit);
    hJetSig[iRes][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetSig[iRes][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetSig[iRes][0] -> GetXaxis() -> SetLabelSize(fLab);
    hJetSig[iRes][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetSig[iRes][0] -> GetXaxis() -> CenterTitle(fCnt);
    hJetSig[iRes][0] -> GetYaxis() -> SetTitle(sJetSigY.Data());
    hJetSig[iRes][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetSig[iRes][0] -> GetYaxis() -> SetTitleSize(fTit);
    hJetSig[iRes][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetSig[iRes][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetSig[iRes][0] -> GetYaxis() -> SetLabelSize(fLab);
    hJetSig[iRes][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetSig[iRes][0] -> GetYaxis() -> CenterTitle(fCnt);
    hJetSig[iRes][1] -> SetLineColor(fColTrg[1]);
    hJetSig[iRes][1] -> SetLineStyle(fLinAll);
    hJetSig[iRes][1] -> SetFillColor(fColTrg[1]);
    hJetSig[iRes][1] -> SetFillStyle(fFilAll);
    hJetSig[iRes][1] -> SetMarkerColor(fColTrg[1]);
    hJetSig[iRes][1] -> SetMarkerStyle(fMarAll);
    hJetSig[iRes][1] -> SetTitle("");
    hJetSig[iRes][1] -> SetTitleFont(fTxt);
    hJetSig[iRes][1] -> GetXaxis() -> SetTitle(sJetSig.Data());
    hJetSig[iRes][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetSig[iRes][1] -> GetXaxis() -> SetTitleSize(fTit);
    hJetSig[iRes][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetSig[iRes][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetSig[iRes][1] -> GetXaxis() -> SetLabelSize(fLab);
    hJetSig[iRes][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetSig[iRes][1] -> GetXaxis() -> CenterTitle(fCnt);
    hJetSig[iRes][1] -> GetYaxis() -> SetTitle(sJetSigY.Data());
    hJetSig[iRes][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetSig[iRes][1] -> GetYaxis() -> SetTitleSize(fTit);
    hJetSig[iRes][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetSig[iRes][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetSig[iRes][1] -> GetYaxis() -> SetLabelSize(fLab);
    hJetSig[iRes][1] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetSig[iRes][1] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end resolution parameter loop
  cout << "      Styles set." << endl;

  // determine area lines
  const UInt_t fColLin(1);
  const UInt_t fWidCut(2);
  const UInt_t fWidCirc(2);
  const UInt_t fStyCut(7);
  const UInt_t fStyCirc(1);

  TH1D  *hProjCut;
  TH1D  *hProjCirc;
  TLine *lCut1D[NResVals];
  TLine *lCut2D[NResVals];
  TLine *lCirc1D[NResVals];
  TLine *lCirc2D[NResVals];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {

    // x-coordinates
    const Double_t xCut  = aJetMin[iRes];
    const Double_t xCirc = pi * (rJet[iRes] * rJet[iRes]);

    // 1d y-coordinates
    const UInt_t   iCut1D  = hJetA[iRes][0] -> FindBin(xCut);
    const UInt_t   iCirc1D = hJetA[iRes][0] -> FindBin(xCirc);
    const Double_t yCut1D  = hJetA[iRes][0] -> GetBinContent(iCut1D);
    const Double_t yCirc1D = hJetA[iRes][0] -> GetBinContent(iCirc1D);
    const Double_t yMin1D  = hJetA[iRes][0] -> GetMinimum(0.);

    // 2d y-coordinates
    const UInt_t iCut2D  = hJetPtVsA[iRes][0] -> GetXaxis() -> FindBin(xCut);
    const UInt_t iCirc2D = hJetPtVsA[iRes][0] -> GetXaxis() -> FindBin(xCirc);
    hProjCut  = (TH1D*) hJetPtVsA[iRes][0] -> ProjectionY("hProjCut", iCut2D, iCut2D);
    hProjCirc = (TH1D*) hJetPtVsA[iRes][0] -> ProjectionY("hProjCirc", iCirc2D, iCirc2D);

    const UInt_t   iCutStart2D  = hProjCut  -> FindFirstBinAbove(0.);
    const UInt_t   iCutStop2D   = hProjCut  -> FindLastBinAbove(0.);
    const UInt_t   iCircStart2D = hProjCirc -> FindFirstBinAbove(0.);
    const UInt_t   iCircStop2D  = hProjCirc -> FindLastBinAbove(0.);
    const Double_t yCutStart2D  = hProjCut  -> GetBinLowEdge(iCutStart2D);
    const Double_t yCutStop2D   = hProjCut  -> GetBinLowEdge(iCutStop2D + 1);
    const Double_t yCircStart2D = hProjCirc -> GetBinLowEdge(iCircStart2D);
    const Double_t yCircStop2D  = hProjCirc -> GetBinLowEdge(iCircStop2D + 1);

    // initialize lines
    lCut1D[iRes]  = new TLine(xCut, yMin1D, xCut, yCut1D);
    lCut2D[iRes]  = new TLine(xCut, yCutStart2D, xCut, yCutStop2D);
    lCirc1D[iRes] = new TLine(xCirc, yMin1D, xCirc, yCirc1D);
    lCirc2D[iRes] = new TLine(xCirc, yCircStart2D, xCirc, yCircStop2D);
    lCut1D[iRes]  -> SetLineColor(fColLin);
    lCut1D[iRes]  -> SetLineWidth(fWidCut);
    lCut1D[iRes]  -> SetLineStyle(fStyCut);
    lCut2D[iRes]  -> SetLineColor(fColLin);
    lCut2D[iRes]  -> SetLineWidth(fWidCut);
    lCut2D[iRes]  -> SetLineStyle(fStyCut);
    lCirc1D[iRes] -> SetLineColor(fColLin);
    lCirc1D[iRes] -> SetLineWidth(fWidCirc);
    lCirc1D[iRes] -> SetLineStyle(fStyCirc);
    lCirc2D[iRes] -> SetLineColor(fColLin);
    lCirc2D[iRes] -> SetLineWidth(fWidCirc);
    lCirc2D[iRes] -> SetLineStyle(fStyCirc);
  }
  cout << "      Made lines." << endl;

  // make legends
  const UInt_t  fColLeg(0);
  const UInt_t  fFilLeg(0);
  const UInt_t  fLinLeg(0);
  const UInt_t  fAlnLeg(12);
  const UInt_t  fAlnTxt(32);
  const Float_t xyLegSmall[4]  = {0.1, 0.1, 0.3, 0.2};
  const Float_t xyLegBig[4]    = {0.1, 0.1, 0.3, 0.5};
  const Float_t xyTxtOnePad[4] = {0.3, 0.1, 0.5, 0.3};
  const Float_t xyTxtTwoPad[4] = {0.3, 0.1, 0.5, 0.5};

  TLegend *lJetNum[NResVals];
  TLegend *lJetPt[NResVals];
  TLegend *lJetEta[NResVals];
  TLegend *lJetArea[NResVals];
  TLegend *lJetDfVsPt[NResVals][NTrgTsp];
  TLegend *lJetEtaVsPt[NResVals][NTrgTsp];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {

    lJetNum[iRes]       = new TLegend(xyLegSmall[0], xyLegSmall[1], xyLegSmall[2], xyLegSmall[3]);
    lJetDfVsPt[iRes][0] = new TLegend(xyLegBig[0], xyLegBig[1], xyLegBig[2], xyLegBig[3]);
    lJetNum[iRes]       -> SetFillColor(fColLeg);
    lJetNum[iRes]       -> SetFillStyle(fFilLeg);
    lJetNum[iRes]       -> SetLineColor(fColLeg);
    lJetNum[iRes]       -> SetLineStyle(fLinLeg);
    lJetNum[iRes]       -> SetTextFont(fTxt);
    lJetNum[iRes]       -> SetTextAlign(fAlnLeg);
    lJetDfVsPt[iRes][0] -> SetFillColor(fColLeg);
    lJetDfVsPt[iRes][0] -> SetFillStyle(fFilLeg);
    lJetDfVsPt[iRes][0] -> SetLineColor(fColLeg);
    lJetDfVsPt[iRes][0] -> SetLineStyle(fLinLeg);
    lJetDfVsPt[iRes][0] -> SetTextFont(fTxt);
    lJetDfVsPt[iRes][0] -> SetTextAlign(fAlnLeg);

    // small legends
    lJetPt[iRes]   = (TLegend*) lJetNum[iRes] -> Clone();
    lJetEta[iRes]  = (TLegend*) lJetNum[iRes] -> Clone();

    // big legends
    lJetDfVsPt[iRes][1]  = (TLegend*) lJetDfVsPt[iRes][0] -> Clone();
    lJetEtaVsPt[iRes][0] = (TLegend*) lJetDfVsPt[iRes][0] -> Clone();
    lJetEtaVsPt[iRes][1] = (TLegend*) lJetDfVsPt[iRes][0] -> Clone();
    lJetArea[iRes]       = (TLegend*) lJetDfVsPt[iRes][0] -> Clone();

    // jet number
    lJetNum[iRes] -> AddEntry(hJetN[iRes][0], "number of recoil jets #color[859]{[#pi^{0} trig.]}", "pf");
    lJetNum[iRes] -> AddEntry(hJetN[iRes][1], "number of recoil jets #color[899]{[#gamma_{rich} trig.]}", "pf");

    // jet pT
    lJetPt[iRes] -> AddEntry(hJetPt[iRes][0], "all jets", "pf");
    lJetPt[iRes] -> AddEntry(hJetPt[iRes][1], "excluded jets", "f");

    // jet eta
    lJetEta[iRes] -> AddEntry(hJetEta[iRes][0], "all jets", "pf");
    lJetEta[iRes] -> AddEntry(hJetEta[iRes][1], "excluded jets", "f");

    // jet area
    lJetArea[iRes] -> AddEntry(hJetA[iRes][0], "all jets", "pf");
    lJetArea[iRes] -> AddEntry(hJetA[iRes][1], "excluded jets", "f");
    lJetArea[iRes] -> AddEntry(lCut1D[iRes], "area cut", "l");
    lJetArea[iRes] -> AddEntry(lCirc1D[iRes], "#piR^{2}", "l");

    // jet delta-phi
    lJetDfVsPt[iRes][0] -> AddEntry(hJetDfBin[iRes][0][1], "all accepted jets", "pf");
    lJetDfVsPt[iRes][0] -> AddEntry(hJetDfPtBin[iRes][0][0][1], "p_{T}^{jet} #in (0.2, 2) GeV/c", "pf");
    lJetDfVsPt[iRes][0] -> AddEntry(hJetDfPtBin[iRes][0][1][1], "p_{T}^{jet} #in (2, 5) GeV/c", "pf");
    lJetDfVsPt[iRes][0] -> AddEntry(hJetDfPtBin[iRes][0][2][1], "p_{T}^{jet} #in (5, 9) GeV/c", "pf");
    lJetDfVsPt[iRes][0] -> AddEntry(hJetDfPtBin[iRes][0][3][1], "p_{T}^{jet} #in (9, 30) GeV/c", "pf");
    lJetDfVsPt[iRes][1] -> AddEntry(hJetDfBin[iRes][1][1], "all accepted jets", "pf");
    lJetDfVsPt[iRes][1] -> AddEntry(hJetDfPtBin[iRes][1][0][1], "p_{T}^{jet} #in (0.2, 2) GeV/c", "pf");
    lJetDfVsPt[iRes][1] -> AddEntry(hJetDfPtBin[iRes][1][1][1], "p_{T}^{jet} #in (2, 5) GeV/c", "pf");
    lJetDfVsPt[iRes][1] -> AddEntry(hJetDfPtBin[iRes][1][2][1], "p_{T}^{jet} #in (5, 9) GeV/c", "pf");
    lJetDfVsPt[iRes][1] -> AddEntry(hJetDfPtBin[iRes][1][3][1], "p_{T}^{jet} #in (9, 30) Gev/c", "pf");

    // jet eta
    lJetEtaVsPt[iRes][0] -> AddEntry(hJetEtaBin[iRes][0][0], "all accepted jets", "pf");
    lJetEtaVsPt[iRes][0] -> AddEntry(hJetEtaPtBin[iRes][0][0][0], "p_{T}^{jet} #in (0.2, 2) GeV/c", "pf");
    lJetEtaVsPt[iRes][0] -> AddEntry(hJetEtaPtBin[iRes][0][1][0], "p_{T}^{jet} #in (2, 5) GeV/c", "pf");
    lJetEtaVsPt[iRes][0] -> AddEntry(hJetEtaPtBin[iRes][0][2][0], "p_{T}^{jet} #in (5, 9) GeV/c", "pf");
    lJetEtaVsPt[iRes][0] -> AddEntry(hJetEtaPtBin[iRes][0][3][0], "p_{T}^{jet} #in (9, 30) GeV/c", "pf");
    lJetEtaVsPt[iRes][1] -> AddEntry(hJetEtaBin[iRes][1][0], "all accepted jets", "pf");
    lJetEtaVsPt[iRes][1] -> AddEntry(hJetEtaPtBin[iRes][1][0][0], "p_{T}^{jet} #in (0.2, 2) GeV/c", "pf");
    lJetEtaVsPt[iRes][1] -> AddEntry(hJetEtaPtBin[iRes][1][1][0], "p_{T}^{jet} #in (2, 5) GeV/c", "pf");
    lJetEtaVsPt[iRes][1] -> AddEntry(hJetEtaPtBin[iRes][1][2][0], "p_{T}^{jet} #in (5, 9) GeV/c", "pf");
    lJetEtaVsPt[iRes][1] -> AddEntry(hJetEtaPtBin[iRes][1][3][0], "p_{T}^{jet} #in (9, 30) Gev/c", "pf");

  }  // end resolution parameter loop
  cout << "      Made legends." << endl;

  // make text boxes
  TPaveText *pBoth1P[NResVals];
  TPaveText *pBoth2P[NResVals];
  TPaveText *pTrg1P[NResVals][NTrgTsp];
  TPaveText *pTrg2P[NResVals][NTrgTsp];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {

    pBoth1P[iRes] = new TPaveText(xyTxtOnePad[0], xyTxtOnePad[1], xyTxtOnePad[2], xyTxtOnePad[3], "NDC NB");
    pBoth2P[iRes] = new TPaveText(xyTxtTwoPad[0], xyTxtTwoPad[1], xyTxtTwoPad[2], xyTxtTwoPad[3], "NDC NB");
    pBoth1P[iRes] -> SetFillColor(fColLeg);
    pBoth1P[iRes] -> SetFillStyle(fFilLeg);
    pBoth1P[iRes] -> SetLineColor(fColLeg);
    pBoth1P[iRes] -> SetTextFont(fTxt);
    pBoth1P[iRes] -> SetTextAlign(fAlnTxt);
    pBoth2P[iRes] -> SetFillColor(fColLeg);
    pBoth2P[iRes] -> SetFillStyle(fFilLeg);
    pBoth2P[iRes] -> SetLineColor(fColLeg);
    pBoth2P[iRes] -> SetTextFont(fTxt);
    pBoth2P[iRes] -> SetTextAlign(fAlnTxt);

    pTrg1P[iRes][0] = (TPaveText*) pBoth1P[iRes] -> Clone();
    pTrg1P[iRes][1] = (TPaveText*) pBoth1P[iRes] -> Clone();
    pTrg2P[iRes][0] = (TPaveText*) pBoth2P[iRes] -> Clone();
    pTrg2P[iRes][1] = (TPaveText*) pBoth2P[iRes] -> Clone();
    pBoth1P[iRes]   -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
    pBoth1P[iRes]   -> AddText("#pi^{0}/#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");
    pBoth1P[iRes]   -> AddText(sResTxt[iRes].Data());
    pBoth1P[iRes]   -> AddText("#bf{charged jets}");
    pBoth2P[iRes]   -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
    pBoth2P[iRes]   -> AddText("#pi^{0}/#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");
    pBoth2P[iRes]   -> AddText(sResTxt[iRes].Data());
    pBoth2P[iRes]   -> AddText("#bf{charged jets}");
    pTrg1P[iRes][0] -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
    pTrg1P[iRes][0] -> AddText("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
    pTrg1P[iRes][0] -> AddText(sResTxt[iRes].Data());
    pTrg1P[iRes][0] -> AddText("#bf{charged jets}");
    pTrg1P[iRes][1] -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
    pTrg1P[iRes][1] -> AddText("#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");
    pTrg1P[iRes][1] -> AddText(sResTxt[iRes].Data());
    pTrg1P[iRes][1] -> AddText("#bf{charged jets}");
    pTrg2P[iRes][0] -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
    pTrg2P[iRes][0] -> AddText("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
    pTrg2P[iRes][0] -> AddText(sResTxt[iRes].Data());
    pTrg2P[iRes][0] -> AddText("#bf{charged jets}");
    pTrg2P[iRes][1] -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
    pTrg2P[iRes][1] -> AddText("#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");
    pTrg2P[iRes][1] -> AddText(sResTxt[iRes].Data());
    pTrg2P[iRes][1] -> AddText("#bf{charged jets}");
  }
  cout << "      Made text boxes." << endl;

  // make jet selection boxes
  const UInt_t   fColSel(923);
  const UInt_t   fFilSel(3354);
  const UInt_t   fLinSel(1);
  const Double_t dFrecoil[2]              = {pi - dFrecoilMax, pi + dFrecoilMax};
  const Double_t hRange[NResVals][2]      = {{-1. * hJetMax[0], hJetMax[0]}, {-1. * hJetMax[1], hJetMax[1]}};
  const Double_t pTvsAreaSel[NResVals][4] = {{aJetMin[0], pTjetMin, xAreaPlot[0][1], pTjetMax}, {aJetMin[1], pTjetMin, xAreaPlot[1][1], pTjetMax}};
  const Double_t pTvsCorrSel[NResVals][4] = {{xCorrPlot[0], pTjetMin, xCorrPlot[1], pTjetMax},  {xCorrPlot[0], pTjetMin, xCorrPlot[1], pTjetMax}};
  const Double_t pTvsEtaSel[NResVals][4]  = {{hRange[0][0], yPtPlot[0], hRange[0][1], pTjetMax},  {hRange[1][0], yPtPlot[0], hRange[1][1], pTjetMax}};
  const Double_t pTvsDfSel[NResVals][4]   = {{dFrecoil[0], yPtPlot[0], dFrecoil[1], pTjetMax},    {dFrecoil[0], yPtPlot[0], dFrecoil[1], pTjetMax}};

  TBox *bPtVsAreaSel[NResVals];
  TBox *bPtVsCorrSel[NResVals];
  TBox *bPtVsEtaSel[NResVals];
  TBox *bPtVsDfSel[NResVals];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    bPtVsAreaSel[iRes] = new TBox(pTvsAreaSel[iRes][0], pTvsAreaSel[iRes][1], pTvsAreaSel[iRes][2], pTvsAreaSel[iRes][3]);
    bPtVsCorrSel[iRes] = new TBox(pTvsCorrSel[iRes][0], pTvsCorrSel[iRes][1], pTvsCorrSel[iRes][2], pTvsCorrSel[iRes][3]);
    bPtVsEtaSel[iRes]  = new TBox(pTvsEtaSel[iRes][0], pTvsEtaSel[iRes][1], pTvsEtaSel[iRes][2], pTvsEtaSel[iRes][3]);
    bPtVsDfSel[iRes]   = new TBox(pTvsDfSel[iRes][0], pTvsDfSel[iRes][1], pTvsDfSel[iRes][2], pTvsDfSel[iRes][3]);
    bPtVsAreaSel[iRes] -> SetFillColor(fColSel);
    bPtVsAreaSel[iRes] -> SetFillStyle(fFilSel);
    bPtVsAreaSel[iRes] -> SetLineColor(fColSel);
    bPtVsAreaSel[iRes] -> SetLineStyle(fLinSel);
    bPtVsCorrSel[iRes] -> SetFillColor(fColSel);
    bPtVsCorrSel[iRes] -> SetFillStyle(fFilSel);
    bPtVsCorrSel[iRes] -> SetLineColor(fColSel);
    bPtVsCorrSel[iRes] -> SetLineStyle(fLinSel);
    bPtVsEtaSel[iRes]  -> SetFillColor(fColSel);
    bPtVsEtaSel[iRes]  -> SetFillStyle(fFilSel);
    bPtVsEtaSel[iRes]  -> SetLineColor(fColSel);
    bPtVsEtaSel[iRes]  -> SetLineStyle(fLinSel);
    bPtVsDfSel[iRes]   -> SetFillColor(fColSel);
    bPtVsDfSel[iRes]   -> SetFillStyle(fFilSel);
    bPtVsDfSel[iRes]   -> SetLineColor(fColSel);
    bPtVsDfSel[iRes]   -> SetLineStyle(fLinSel);
  }
  cout << "      Made 2d selection boxes." << endl;

  // draw plots
  const UInt_t  width(750);
  const UInt_t  bigWidth(1500);
  const UInt_t  height(750);
  const UInt_t  bigHeight(1500);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fFrame(0);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fLogY2(0);
  const UInt_t  fLogZ(1);
  const UInt_t  fTick(1);
  const Float_t fMarginBig(0.15);
  const Float_t fMarginBord(0.005);
  const Float_t fMarginSmall(0.02);

  // canvas and pad name bases
  const TString sCanNumBase("cJetNum");
  const TString sCanPtBase("cJetPtCutVsAll");
  const TString sCanEtaBase("cJetEtaCutVsAll");
  const TString sCanAreaBase("cJetArea");
  const TString sCanCorrBase("cJetCorr");
  const TString sCanDfBinBase("cJetDfBin");
  const TString sCanEtaBinBase("cJetEtaBin");
  const TString sPadAreaVsPtBase("pJetAreaVsPt");
  const TString sPadAreaBase("pJetArea");
  const TString sPadCorrVsPtBase[NTrgTsp] = {"pJetCorrVsPtPi0", "pJetCorrVsPtGam"};
  const TString sPadEtaVsPtBase[NTrgTsp]  = {"pJetEtaVsPtPi0", "pJetEtaVsPtGam"};
  const TString sPadDfVsPtBase[NTrgTsp]   = {"pJetDfVsPtPi0", "pJetDfVsPtGam"};
  const TString sPadCorrBase[NTrgTsp]     = {"pJetCorrPi0", "pJetCorrGam"};
  const TString sPadEtaBase[NTrgTsp]      = {"pJetEtaPi0", "pJetEtaGam"};
  const TString sPadDfBase[NTrgTsp]       = {"pJetDfPi0", "pJetDfGam"};

  // loop over resolution parameters
  TCanvas *cJetNum[NResVals];
  TCanvas *cJetPt[NResVals];
  TCanvas *cJetEta[NResVals];
  TCanvas *cJetArea[NResVals];
  TCanvas *cJetCorr[NResVals];
  TCanvas *cJetDfBin[NResVals];
  TCanvas *cJetEtaBin[NResVals];
  TPad    *pJetAreaVsPt[NResVals];
  TPad    *pJetArea[NResVals];
  TPad    *pJetCorrVsPt[NResVals][NTrgTsp];
  TPad    *pJetEtaVsPt[NResVals][NTrgTsp];
  TPad    *pJetDfVsPt[NResVals][NTrgTsp];
  TPad    *pJetCorr[NResVals][NTrgTsp];
  TPad    *pJetEta[NResVals][NTrgTsp];
  TPad    *pJetDf[NResVals][NTrgTsp];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {

    // create canvas names
    TString sCanNum(sCanNumBase.Data());
    TString sCanPt(sCanPtBase.Data());
    TString sCanEta(sCanEtaBase.Data());
    TString sCanArea(sCanAreaBase.Data());
    TString sCanCorr(sCanCorrBase.Data());
    TString sCanDfBin(sCanDfBinBase.Data());
    TString sCanEtaBin(sCanEtaBinBase.Data());
    sCanNum.Append(sResParams[iRes]);
    sCanPt.Append(sResParams[iRes]);
    sCanEta.Append(sResParams[iRes]);
    sCanArea.Append(sResParams[iRes]);
    sCanCorr.Append(sResParams[iRes]);
    sCanDfBin.Append(sResParams[iRes]);
    sCanEtaBin.Append(sResParams[iRes]);

    // create pad names
    TString sPadAreaVsPt(sPadAreaVsPtBase.Data());
    TString sPadCorrVsPtP(sPadCorrVsPtBase[0].Data());
    TString sPadCorrVsPtG(sPadCorrVsPtBase[1].Data());
    TString sPadEtaVsPtP(sPadEtaVsPtBase[0].Data());
    TString sPadEtaVsPtG(sPadEtaVsPtBase[1].Data());
    TString sPadDfVsPtP(sPadDfVsPtBase[0].Data());
    TString sPadDfVsPtG(sPadDfVsPtBase[1].Data());
    TString sPadArea(sPadAreaBase.Data());
    TString sPadCorrP(sPadCorrBase[0].Data());
    TString sPadCorrG(sPadCorrBase[1].Data());
    TString sPadEtaP(sPadEtaBase[0].Data());
    TString sPadEtaG(sPadEtaBase[1].Data());
    TString sPadDfP(sPadDfBase[0].Data());
    TString sPadDfG(sPadDfBase[1].Data());
    sPadAreaVsPt.Append(sResParams[iRes]);
    sPadCorrVsPtP.Append(sResParams[iRes]);
    sPadCorrVsPtG.Append(sResParams[iRes]);
    sPadEtaVsPtP.Append(sResParams[iRes]);
    sPadEtaVsPtG.Append(sResParams[iRes]);
    sPadDfVsPtP.Append(sResParams[iRes]);
    sPadDfVsPtG.Append(sResParams[iRes]);
    sPadArea.Append(sResParams[iRes]);
    sPadCorrP.Append(sResParams[iRes]);
    sPadCorrG.Append(sResParams[iRes]);
    sPadEtaP.Append(sResParams[iRes]);
    sPadEtaG.Append(sResParams[iRes]);
    sPadDfP.Append(sResParams[iRes]);
    sPadDfG.Append(sResParams[iRes]);

    // jet no's
    cJetNum[iRes] = new TCanvas(sCanNum.Data(), "", width, height);
    cJetNum[iRes]  -> SetLogx(fLogX);
    cJetNum[iRes]  -> SetLogy(fLogY);
    cJetNum[iRes]  -> SetGrid(fGrid, fGrid);
    cJetNum[iRes]  -> SetTicks(fTick, fTick);
    cJetNum[iRes]  -> SetBorderMode(fMode);
    cJetNum[iRes]  -> SetBorderSize(fBord);
    cJetNum[iRes]  -> SetFrameBorderMode(fFrame);
    cJetNum[iRes]  -> SetLeftMargin(fMarginBig);
    cJetNum[iRes]  -> SetTopMargin(fMarginSmall);
    cJetNum[iRes]  -> SetRightMargin(fMarginSmall);
    cJetNum[iRes]  -> SetBottomMargin(fMarginBig);
    cJetNum[iRes]  -> cd();
    hJetN[iRes][0] -> Draw();
    hJetN[iRes][1] -> Draw("same");
    lJetNum[iRes]  -> Draw();
    pBoth1P[iRes]  -> Draw();
    fOutput        -> cd();
    cJetNum[iRes]  -> Write();
    cJetNum[iRes]  -> Close();

    // jet pT
    cJetPt[iRes] = new TCanvas(sCanPt.Data(), "", width, height);
    cJetPt[iRes]    -> SetLogx(fLogX);
    cJetPt[iRes]    -> SetLogy(fLogY);
    cJetPt[iRes]    -> SetGrid(fGrid, fGrid);
    cJetPt[iRes]    -> SetTicks(fTick, fTick);
    cJetPt[iRes]    -> SetBorderMode(fMode);
    cJetPt[iRes]    -> SetBorderSize(fBord);
    cJetPt[iRes]    -> SetFrameBorderMode(fFrame);
    cJetPt[iRes]    -> SetLeftMargin(fMarginBig);
    cJetPt[iRes]    -> SetTopMargin(fMarginSmall);
    cJetPt[iRes]    -> SetRightMargin(fMarginSmall);
    cJetPt[iRes]    -> SetBottomMargin(fMarginBig);
    cJetPt[iRes]    -> cd();
    hJetPt[iRes][0] -> Draw();
    hJetPt[iRes][1] -> Draw("same");
    hJetPt[iRes][1] -> Draw("hist same");
    lJetPt[iRes]    -> Draw();
    pBoth1P[iRes]   -> Draw();
    fOutput         -> cd();
    cJetPt[iRes]    -> Write();
    cJetPt[iRes]    -> Close();

    // jet eta
    cJetEta[iRes] = new TCanvas(sCanEta.Data(), "", width, height);
    cJetEta[iRes]    -> SetLogx(fLogX);
    cJetEta[iRes]    -> SetLogy(fLogY);
    cJetEta[iRes]    -> SetGrid(fGrid, fGrid);
    cJetEta[iRes]    -> SetTicks(fTick, fTick);
    cJetEta[iRes]    -> SetBorderMode(fMode);
    cJetEta[iRes]    -> SetBorderSize(fBord);
    cJetEta[iRes]    -> SetFrameBorderMode(fFrame);
    cJetEta[iRes]    -> SetLeftMargin(fMarginBig);
    cJetEta[iRes]    -> SetTopMargin(fMarginSmall);
    cJetEta[iRes]    -> SetRightMargin(fMarginSmall);
    cJetEta[iRes]    -> SetBottomMargin(fMarginBig);
    cJetEta[iRes]    -> cd();
    hJetEta[iRes][0] -> Draw();
    hJetEta[iRes][1] -> Draw("same");
    hJetEta[iRes][1] -> Draw("hist same");
    lJetEta[iRes]    -> Draw();
    pBoth1P[iRes]    -> Draw();
    fOutput          -> cd();
    cJetEta[iRes]    -> Write();
    cJetEta[iRes]    -> Close();

    // jet area
    cJetArea[iRes]     = new TCanvas(sCanArea.Data(), "", width, bigHeight);
    pJetAreaVsPt[iRes] = new TPad(sPadAreaVsPt.Data(), "", 0., 0., 1., 0.5);
    pJetArea[iRes]     = new TPad(sPadArea.Data(),     "", 0., 0.5, 1., 1.);
    pJetAreaVsPt[iRes] -> SetLogx(fLogX);
    pJetAreaVsPt[iRes] -> SetLogy(fLogY2);
    pJetAreaVsPt[iRes] -> SetLogz(fLogZ);
    pJetAreaVsPt[iRes] -> SetGrid(fGrid, fGrid);
    pJetAreaVsPt[iRes] -> SetTicks(fTick, fTick);
    pJetAreaVsPt[iRes] -> SetBorderMode(fMode);
    pJetAreaVsPt[iRes] -> SetBorderSize(fBord);
    pJetAreaVsPt[iRes] -> SetFrameBorderMode(fFrame);
    pJetAreaVsPt[iRes] -> SetLeftMargin(fMarginBig);
    pJetAreaVsPt[iRes] -> SetTopMargin(fMarginBord);
    pJetAreaVsPt[iRes] -> SetRightMargin(fMarginBig);
    pJetAreaVsPt[iRes] -> SetBottomMargin(fMarginBig);
    pJetArea[iRes]     -> SetLogx(fLogX);
    pJetArea[iRes]     -> SetLogy(fLogY2);
    pJetArea[iRes]     -> SetGrid(fGrid, fGrid);
    pJetArea[iRes]     -> SetTicks(fTick, fTick);
    pJetArea[iRes]     -> SetBorderMode(fMode);
    pJetArea[iRes]     -> SetBorderSize(fBord);
    pJetArea[iRes]     -> SetFrameBorderMode(fFrame);
    pJetArea[iRes]     -> SetLeftMargin(fMarginBig);
    pJetArea[iRes]     -> SetTopMargin(fMarginSmall);
    pJetArea[iRes]     -> SetRightMargin(fMarginBig);
    pJetArea[iRes]     -> SetBottomMargin(fMarginBord);
    cJetArea[iRes]     -> cd();
    pJetAreaVsPt[iRes] -> Draw();
    pJetArea[iRes]     -> Draw();
    pJetAreaVsPt[iRes] -> cd();
    hJetPtVsA[iRes][0] -> Draw("colz");
    bPtVsAreaSel[iRes] -> Draw();
    lCut2D[iRes]       -> Draw();
    pBoth2P[iRes]      -> Draw();
    pJetArea[iRes]     -> cd();
    hJetA[iRes][0]     -> Draw("hist ][");
    hJetA[iRes][0]     -> Draw("hist p ][ same");
    hJetA[iRes][1]     -> Draw("hist ][ same");
    hJetA[iRes][1]     -> Draw("hist p ][ same");
    lCut1D[iRes]       -> Draw();
    lCirc1D[iRes]      -> Draw();
    lJetArea[iRes]     -> Draw();
    fOutput            -> cd();
    cJetArea[iRes]     -> Write();
    cJetArea[iRes]     -> Close();

    // jet corr
    cJetCorr[iRes]        = new TCanvas(sCanCorr.Data(), "", bigWidth, bigHeight);
    pJetCorrVsPt[iRes][0] = new TPad(sPadCorrVsPtP.Data(), "", 0., 0., 0.5, 0.5);
    pJetCorrVsPt[iRes][1] = new TPad(sPadCorrVsPtG.Data(), "", 0.5, 0., 1., 0.5);
    pJetCorr[iRes][0]     = new TPad(sPadCorrP.Data(),     "", 0., 0.5, 0.5, 1.);
    pJetCorr[iRes][1]     = new TPad(sPadCorrG.Data(),     "", 0.5, 0.5, 1., 1.);
    pJetCorrVsPt[iRes][0] -> SetLogx(fLogX);
    pJetCorrVsPt[iRes][0] -> SetLogy(fLogY2);
    pJetCorrVsPt[iRes][0] -> SetLogz(fLogZ);
    pJetCorrVsPt[iRes][0] -> SetGrid(fGrid, fGrid);
    pJetCorrVsPt[iRes][0] -> SetTicks(fTick, fTick);
    pJetCorrVsPt[iRes][0] -> SetBorderMode(fMode);
    pJetCorrVsPt[iRes][0] -> SetBorderSize(fBord);
    pJetCorrVsPt[iRes][0] -> SetFrameBorderMode(fFrame);
    pJetCorrVsPt[iRes][0] -> SetLeftMargin(fMarginBig);
    pJetCorrVsPt[iRes][0] -> SetTopMargin(fMarginBord);
    pJetCorrVsPt[iRes][0] -> SetRightMargin(fMarginBord);
    pJetCorrVsPt[iRes][0] -> SetBottomMargin(fMarginBig);
    pJetCorrVsPt[iRes][1] -> SetLogx(fLogX);
    pJetCorrVsPt[iRes][1] -> SetLogy(fLogY2);
    pJetCorrVsPt[iRes][1] -> SetLogz(fLogZ);
    pJetCorrVsPt[iRes][1] -> SetGrid(fGrid, fGrid);
    pJetCorrVsPt[iRes][1] -> SetTicks(fTick, fTick);
    pJetCorrVsPt[iRes][1] -> SetBorderMode(fMode);
    pJetCorrVsPt[iRes][1] -> SetBorderSize(fBord);
    pJetCorrVsPt[iRes][1] -> SetFrameBorderMode(fFrame);
    pJetCorrVsPt[iRes][1] -> SetLeftMargin(fMarginBord);
    pJetCorrVsPt[iRes][1] -> SetTopMargin(fMarginBord);
    pJetCorrVsPt[iRes][1] -> SetRightMargin(fMarginBig);
    pJetCorrVsPt[iRes][1] -> SetBottomMargin(fMarginBig);
    pJetCorr[iRes][0]     -> SetLogx(fLogX);
    pJetCorr[iRes][0]     -> SetLogy(fLogY);
    pJetCorr[iRes][0]     -> SetGrid(fGrid, fGrid);
    pJetCorr[iRes][0]     -> SetTicks(fTick, fTick);
    pJetCorr[iRes][0]     -> SetBorderMode(fMode);
    pJetCorr[iRes][0]     -> SetBorderSize(fBord);
    pJetCorr[iRes][0]     -> SetFrameBorderMode(fFrame);
    pJetCorr[iRes][0]     -> SetLeftMargin(fMarginBig);
    pJetCorr[iRes][0]     -> SetTopMargin(fMarginSmall);
    pJetCorr[iRes][0]     -> SetRightMargin(fMarginBord);
    pJetCorr[iRes][0]     -> SetBottomMargin(fMarginBord);
    pJetCorr[iRes][1]     -> SetLogx(fLogX);
    pJetCorr[iRes][1]     -> SetLogy(fLogY);
    pJetCorr[iRes][1]     -> SetGrid(fGrid, fGrid);
    pJetCorr[iRes][1]     -> SetTicks(fTick, fTick);
    pJetCorr[iRes][1]     -> SetBorderMode(fMode);
    pJetCorr[iRes][1]     -> SetBorderSize(fBord);
    pJetCorr[iRes][1]     -> SetFrameBorderMode(fFrame);
    pJetCorr[iRes][1]     -> SetLeftMargin(fMarginBord);
    pJetCorr[iRes][1]     -> SetTopMargin(fMarginSmall);
    pJetCorr[iRes][1]     -> SetRightMargin(fMarginBig);
    pJetCorr[iRes][1]     -> SetBottomMargin(fMarginBord);
    cJetCorr[iRes]        -> cd();
    pJetCorrVsPt[iRes][0] -> Draw();
    pJetCorrVsPt[iRes][1] -> Draw();
    pJetCorr[iRes][0]     -> Draw();
    pJetCorr[iRes][1]     -> Draw();
    pJetCorrVsPt[iRes][0] -> cd();
    hJetPtVsCorr[iRes][0] -> Draw("col");
    bPtVsCorrSel[iRes]    -> Draw();
    pJetCorrVsPt[iRes][1] -> cd();
    hJetPtVsCorr[iRes][1] -> Draw("colz");
    bPtVsCorrSel[iRes]    -> Draw();
    pJetCorr[iRes][0]     -> cd();
    hJetCorr[iRes][0]     -> Draw();
    pTrg2P[iRes][0]         -> Draw();
    pJetCorr[iRes][1]     -> cd();
    hJetCorr[iRes][1]     -> Draw();
    pTrg2P[iRes][1]         -> Draw();
    fOutput               -> cd();
    cJetCorr[iRes]        -> Write();
    cJetCorr[iRes]        -> Close();

    // delta phi
    cJetDfBin[iRes]     = new TCanvas(sCanDfBin.Data(), "", bigWidth, bigHeight);
    pJetDfVsPt[iRes][0] = new TPad(sPadDfVsPtP.Data(), "", 0., 0., 0.5, 0.5);
    pJetDfVsPt[iRes][1] = new TPad(sPadDfVsPtG.Data(), "", 0.5, 0., 1., 0.5);
    pJetDf[iRes][0]     = new TPad(sPadDfP.Data(),     "", 0., 0.5, 0.5, 1.);
    pJetDf[iRes][1]     = new TPad(sPadDfG.Data(),     "", 0.5, 0.5, 1., 1.);
    pJetDfVsPt[iRes][0]    -> SetLogx(fLogX);
    pJetDfVsPt[iRes][0]    -> SetLogy(fLogY2);
    pJetDfVsPt[iRes][0]    -> SetLogz(fLogZ);
    pJetDfVsPt[iRes][0]    -> SetGrid(fGrid, fGrid);
    pJetDfVsPt[iRes][0]    -> SetTicks(fTick, fTick);
    pJetDfVsPt[iRes][0]    -> SetBorderMode(fMode);
    pJetDfVsPt[iRes][0]    -> SetBorderSize(fBord);
    pJetDfVsPt[iRes][0]    -> SetFrameBorderMode(fFrame);
    pJetDfVsPt[iRes][0]    -> SetLeftMargin(fMarginBig);
    pJetDfVsPt[iRes][0]    -> SetTopMargin(fMarginBord);
    pJetDfVsPt[iRes][0]    -> SetRightMargin(fMarginBord);
    pJetDfVsPt[iRes][0]    -> SetBottomMargin(fMarginBig);
    pJetDfVsPt[iRes][1]    -> SetLogx(fLogX);
    pJetDfVsPt[iRes][1]    -> SetLogy(fLogY2);
    pJetDfVsPt[iRes][1]    -> SetLogz(fLogZ);
    pJetDfVsPt[iRes][1]    -> SetGrid(fGrid, fGrid);
    pJetDfVsPt[iRes][1]    -> SetTicks(fTick, fTick);
    pJetDfVsPt[iRes][1]    -> SetBorderMode(fMode);
    pJetDfVsPt[iRes][1]    -> SetBorderSize(fBord);
    pJetDfVsPt[iRes][1]    -> SetFrameBorderMode(fFrame);
    pJetDfVsPt[iRes][1]    -> SetLeftMargin(fMarginBord);
    pJetDfVsPt[iRes][1]    -> SetTopMargin(fMarginBord);
    pJetDfVsPt[iRes][1]    -> SetRightMargin(fMarginBig);
    pJetDfVsPt[iRes][1]    -> SetBottomMargin(fMarginBig);
    pJetDf[iRes][0]        -> SetLogx(fLogX);
    pJetDf[iRes][0]        -> SetLogy(fLogY);
    pJetDf[iRes][0]        -> SetGrid(fGrid, fGrid);
    pJetDf[iRes][0]        -> SetTicks(fTick, fTick);
    pJetDf[iRes][0]        -> SetBorderMode(fMode);
    pJetDf[iRes][0]        -> SetBorderSize(fBord);
    pJetDf[iRes][0]        -> SetFrameBorderMode(fFrame);
    pJetDf[iRes][0]        -> SetLeftMargin(fMarginBig);
    pJetDf[iRes][0]        -> SetTopMargin(fMarginSmall);
    pJetDf[iRes][0]        -> SetRightMargin(fMarginBord);
    pJetDf[iRes][0]        -> SetBottomMargin(fMarginBord);
    pJetDf[iRes][1]        -> SetLogx(fLogX);
    pJetDf[iRes][1]        -> SetLogy(fLogY);
    pJetDf[iRes][1]        -> SetGrid(fGrid, fGrid);
    pJetDf[iRes][1]        -> SetTicks(fTick, fTick);
    pJetDf[iRes][1]        -> SetBorderMode(fMode);
    pJetDf[iRes][1]        -> SetBorderSize(fBord);
    pJetDf[iRes][1]        -> SetFrameBorderMode(fFrame);
    pJetDf[iRes][1]        -> SetLeftMargin(fMarginBord);
    pJetDf[iRes][1]        -> SetTopMargin(fMarginSmall);
    pJetDf[iRes][1]        -> SetRightMargin(fMarginBig);
    pJetDf[iRes][1]        -> SetBottomMargin(fMarginBord);
    cJetDfBin[iRes]        -> cd();
    pJetDfVsPt[iRes][0]    -> Draw();
    pJetDfVsPt[iRes][1]    -> Draw();
    pJetDf[iRes][0]        -> Draw();
    pJetDf[iRes][1]        -> Draw();
    pJetDfVsPt[iRes][0]    -> cd();
    hJetPtVsDf[iRes][0][0] -> Draw("col");
    bPtVsDfSel[iRes]       -> Draw();
    pTrg2P[iRes][0]        -> Draw();
    pJetDfVsPt[iRes][1]    -> cd();
    hJetPtVsDf[iRes][1][0] -> Draw("colz");
    bPtVsDfSel[iRes]       -> Draw();
    pTrg2P[iRes][1]        -> Draw();
    pJetDf[iRes][0]        -> cd();
    hJetDfBin[iRes][0][0]  -> Draw();
    hJetDfBin[iRes][0][1]  -> Draw("same");
    for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
      hJetDfPtBin[iRes][0][iTrkBin][0] -> Draw("same");
      hJetDfPtBin[iRes][0][iTrkBin][1] -> Draw("same");
    }
    lJetDfVsPt[iRes][0]   -> Draw();
    pJetDf[iRes][1]       -> cd();
    hJetDfBin[iRes][1][0] -> Draw();
    hJetDfBin[iRes][1][1] -> Draw("same");
    for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
      hJetDfPtBin[iRes][1][iTrkBin][0] -> Draw("same");
      hJetDfPtBin[iRes][1][iTrkBin][1] -> Draw("same");
    }
    lJetDfVsPt[iRes][1] -> Draw();
    fOutput             -> cd();
    cJetDfBin[iRes]     -> Write();
    cJetDfBin[iRes]     -> Close();

    // jet eta vs pT
    cJetEtaBin[iRes]     = new TCanvas(sCanEtaBin.Data(), "", bigWidth, bigHeight);
    pJetEtaVsPt[iRes][0] = new TPad(sPadEtaVsPtP.Data(), "", 0., 0., 0.5, 0.5);
    pJetEtaVsPt[iRes][1] = new TPad(sPadEtaVsPtG.Data(), "", 0.5, 0., 1., 0.5);
    pJetEta[iRes][0]     = new TPad(sPadEtaP.Data(),     "", 0., 0.5, 0.5, 1.);
    pJetEta[iRes][1]     = new TPad(sPadEtaG.Data(),     "", 0.5, 0.5, 1., 1.);
    pJetEtaVsPt[iRes][0]   -> SetLogx(fLogX);
    pJetEtaVsPt[iRes][0]   -> SetLogy(fLogY2);
    pJetEtaVsPt[iRes][0]   -> SetLogz(fLogZ);
    pJetEtaVsPt[iRes][0]   -> SetGrid(fGrid, fGrid);
    pJetEtaVsPt[iRes][0]   -> SetTicks(fTick, fTick);
    pJetEtaVsPt[iRes][0]   -> SetBorderMode(fMode);
    pJetEtaVsPt[iRes][0]   -> SetBorderSize(fBord);
    pJetEtaVsPt[iRes][0]   -> SetFrameBorderMode(fFrame);
    pJetEtaVsPt[iRes][0]   -> SetLeftMargin(fMarginBig);
    pJetEtaVsPt[iRes][0]   -> SetTopMargin(fMarginBord);
    pJetEtaVsPt[iRes][0]   -> SetRightMargin(fMarginBord);
    pJetEtaVsPt[iRes][0]   -> SetBottomMargin(fMarginBig);
    pJetEtaVsPt[iRes][1]   -> SetLogx(fLogX);
    pJetEtaVsPt[iRes][1]   -> SetLogy(fLogY2);
    pJetEtaVsPt[iRes][1]   -> SetLogz(fLogZ);
    pJetEtaVsPt[iRes][1]   -> SetGrid(fGrid, fGrid);
    pJetEtaVsPt[iRes][1]   -> SetTicks(fTick, fTick);
    pJetEtaVsPt[iRes][1]   -> SetBorderMode(fMode);
    pJetEtaVsPt[iRes][1]   -> SetBorderSize(fBord);
    pJetEtaVsPt[iRes][1]   -> SetFrameBorderMode(fFrame);
    pJetEtaVsPt[iRes][1]   -> SetLeftMargin(fMarginBord);
    pJetEtaVsPt[iRes][1]   -> SetTopMargin(fMarginBord);
    pJetEtaVsPt[iRes][1]   -> SetRightMargin(fMarginBig);
    pJetEtaVsPt[iRes][1]   -> SetBottomMargin(fMarginBig);
    pJetEta[iRes][0]       -> SetLogx(fLogX);
    pJetEta[iRes][0]       -> SetLogy(fLogY);
    pJetEta[iRes][0]       -> SetGrid(fGrid, fGrid);
    pJetEta[iRes][0]       -> SetTicks(fTick, fTick);
    pJetEta[iRes][0]       -> SetBorderMode(fMode);
    pJetEta[iRes][0]       -> SetBorderSize(fBord);
    pJetEta[iRes][0]       -> SetFrameBorderMode(fFrame);
    pJetEta[iRes][0]       -> SetLeftMargin(fMarginBig);
    pJetEta[iRes][0]       -> SetTopMargin(fMarginSmall);
    pJetEta[iRes][0]       -> SetRightMargin(fMarginBord);
    pJetEta[iRes][0]       -> SetBottomMargin(fMarginBord);
    pJetEta[iRes][1]       -> SetLogx(fLogX);
    pJetEta[iRes][1]       -> SetLogy(fLogY);
    pJetEta[iRes][1]       -> SetGrid(fGrid, fGrid);
    pJetEta[iRes][1]       -> SetTicks(fTick, fTick);
    pJetEta[iRes][1]       -> SetBorderMode(fMode);
    pJetEta[iRes][1]       -> SetBorderSize(fBord);
    pJetEta[iRes][1]       -> SetFrameBorderMode(fFrame);
    pJetEta[iRes][1]       -> SetLeftMargin(fMarginBord);
    pJetEta[iRes][1]       -> SetTopMargin(fMarginSmall);
    pJetEta[iRes][1]       -> SetRightMargin(fMarginBig);
    pJetEta[iRes][1]       -> SetBottomMargin(fMarginBord);
    cJetEtaBin[iRes]       -> cd();
    pJetEtaVsPt[iRes][0]   -> Draw();
    pJetEtaVsPt[iRes][1]   -> Draw();
    pJetEta[iRes][0]       -> Draw();
    pJetEta[iRes][1]       -> Draw();
    pJetEtaVsPt[iRes][0]   -> cd();
    hJetPtVsEta[iRes][0]   -> Draw("col");
    bPtVsEtaSel[iRes]      -> Draw();
    pTrg2P[iRes][0]         -> Draw();
    pJetEtaVsPt[iRes][1]   -> cd();
    hJetPtVsEta[iRes][1]   -> Draw("colz");
    bPtVsEtaSel[iRes]      -> Draw();
    pTrg2P[iRes][1]         -> Draw();
    pJetEta[iRes][0]       -> cd();
    hJetEtaBin[iRes][0][0] -> Draw();
    for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
      hJetEtaPtBin[iRes][0][iTrkBin][0] -> Draw("same");
    }
    lJetEtaVsPt[iRes][0]   -> Draw();
    pJetEta[iRes][1]       -> cd();
    hJetEtaBin[iRes][1][0] -> Draw();
    for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
      hJetEtaPtBin[iRes][1][iTrkBin][0] -> Draw("same");
    }
    lJetEtaVsPt[iRes][1] -> Draw();
    fOutput              -> cd();
    cJetEtaBin[iRes]     -> Write();
    cJetEtaBin[iRes]     -> Close();
  }  // end resolution parameter loop

  // save histograms
  TDirectory *dRes[NResVals];
  for (UInt_t iRes = 0; iRes < NResVals; iRes++) {
    dRes[iRes] = (TDirectory*) fOutput -> mkdir(sDir[iRes].Data());
    dRes[iRes]             -> cd();
    hJetN[iRes][0]         -> Write();
    hJetN[iRes][1]         -> Write();
    hJetA[iRes][0]         -> Write();
    hJetA[iRes][1]         -> Write();
    hJetEta[iRes][0]       -> Write();
    hJetEta[iRes][1]       -> Write();
    hJetPt[iRes][0]        -> Write();
    hJetPt[iRes][1]        -> Write();
    hJetPtBin[iRes][0]     -> Write();
    hJetPtBin[iRes][1]     -> Write();
    hJetZtBin[iRes][0]     -> Write();
    hJetZtBin[iRes][1]     -> Write();
    hJetDfBin[iRes][0][0]  -> Write();
    hJetDfBin[iRes][1][0]  -> Write();
    hJetDfBin[iRes][0][1]  -> Write();
    hJetDfBin[iRes][1][1]  -> Write();
    hJetEtaBin[iRes][0][0] -> Write();
    hJetEtaBin[iRes][1][0] -> Write();
    hJetEtaBin[iRes][0][1] -> Write();
    hJetEtaBin[iRes][1][1] -> Write();
    hJetCorr[iRes][0]      -> Write();
    hJetCorr[iRes][1]      -> Write();
    hJetSig[iRes][0]       -> Write();
    hJetSig[iRes][1]       -> Write();
    hJetPtVsA[iRes][0]     -> Write();
    hJetPtVsA[iRes][1]     -> Write();
    hJetPtVsCorr[iRes][0]  -> Write();
    hJetPtVsCorr[iRes][1]  -> Write();
    hJetPtVsDf[iRes][0][0] -> Write();
    hJetPtVsDf[iRes][1][0] -> Write();
    hJetPtVsDf[iRes][0][1] -> Write();
    hJetPtVsDf[iRes][1][1] -> Write();
    hJetPtVsEta[iRes][0]   -> Write();
    hJetPtVsEta[iRes][1]   -> Write();
    hJetDfVsEta[iRes][0]   -> Write();
    hJetDfVsEta[iRes][1]   -> Write();
    for (UInt_t iTrgTsp = 0; iTrgTsp < NTrgTsp; iTrgTsp++) {
      for (UInt_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
        hJetDfPtBin[iRes][iTrgTsp][iJetBin][0]  -> Write();
        hJetDfPtBin[iRes][iTrgTsp][iJetBin][1]  -> Write();
        hJetEtaPtBin[iRes][iTrgTsp][iJetBin][0] -> Write();
        hJetEtaPtBin[iRes][iTrgTsp][iJetBin][1] -> Write();
      }  // end pT bin loop
    }  // end tsp loop
  }  // end rJet loop

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "    Made jet QA plots!" << endl;

}

// End ------------------------------------------------------------------------
