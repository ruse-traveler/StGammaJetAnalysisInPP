// 'CreateWeirdList.C'
// Derek Anderson
// 04.10.2017
//
// Use this to extract the file name, run ID,
// and event ID of events where 'StMuEvent'
// and 'StEvent' differ.

#include <cassert>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

using namespace std;


// global constants
static const Int_t nTrkMax = 500;
static const Int_t nTwrMax = 4800;
static const Int_t nMatMax = 100;
// i/o parameters
static const TString sTree("Gfmtodst");



void CreateWeirdList(const TString sIn, const TString sOut, const Bool_t inBatchMode=false) {

  cout << "\n  Beginning weird list creator..." << endl;
  gErrorIgnoreLevel = kError;


  TFile *fIn  = new TFile(sIn.Data(), "read");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    assert(fIn);
  }
  TTree *tTree = (TTree*) fIn -> Get(sTree.Data());
  if (!tTree) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    assert(tTree);
  }

  cout << "    Tree grabbed..." << endl;


  // declaration of leaves types
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
  UInt_t   PrimaryTrackArray_;
  UInt_t   PrimaryTrackArray_fUniqueID[nTrkMax];
  UInt_t   PrimaryTrackArray_fBits[nTrkMax];
  Int_t    PrimaryTrackArray_nHitsFit[nTrkMax];
  Int_t    PrimaryTrackArray_nHitsPoss[nTrkMax];
  Int_t    PrimaryTrackArray_trackFlag[nTrkMax];
  Double_t PrimaryTrackArray_pZ[nTrkMax];
  Double_t PrimaryTrackArray_pX[nTrkMax];
  Double_t PrimaryTrackArray_pY[nTrkMax];
  Double_t PrimaryTrackArray_pT[nTrkMax];
  Double_t PrimaryTrackArray_dEdx[nTrkMax];
  Double_t PrimaryTrackArray_charge[nTrkMax];
  Double_t PrimaryTrackArray_tofBeta[nTrkMax];
  Double_t PrimaryTrackArray_eta[nTrkMax];
  Double_t PrimaryTrackArray_phi[nTrkMax];
  Double_t PrimaryTrackArray_nSigElectron[nTrkMax];
  Double_t PrimaryTrackArray_nSigPion[nTrkMax];
  Double_t PrimaryTrackArray_nSigKaon[nTrkMax];
  Double_t PrimaryTrackArray_nSigProton[nTrkMax];
  Double_t PrimaryTrackArray_dcag[nTrkMax];
  Double_t PrimaryTrackArray_nHits[nTrkMax];
  Double_t PrimaryTrackArray_dEdxHits[nTrkMax];
  Double_t PrimaryTrackArray_firstZPoint[nTrkMax];
  Double_t PrimaryTrackArray_lastZPoint[nTrkMax];
  Double_t PrimaryTrackArray_tofSigElectron[nTrkMax];
  Double_t PrimaryTrackArray_tofSigPion[nTrkMax];
  Double_t PrimaryTrackArray_tofSigKaon[nTrkMax];
  Double_t PrimaryTrackArray_tofSigProton[nTrkMax];
  Double_t PrimaryTrackArray_timeOfflight[nTrkMax];
  Double_t PrimaryTrackArray_pathLength[nTrkMax];
  Int_t    PrimaryTrackArray_trkIndex[nTrkMax];
  UInt_t   TowerArray_;
  UInt_t   TowerArray_fUniqueID[nTwrMax];
  UInt_t   TowerArray_fBits[nTwrMax];
  Int_t    TowerArray_TwrId[nTwrMax];
  Float_t  TowerArray_TwrEng[nTwrMax];
  Float_t  TowerArray_TwrEta[nTwrMax];
  Float_t  TowerArray_TwrPhi[nTwrMax];
  Float_t  TowerArray_TwrADC[nTwrMax];
  Float_t  TowerArray_TwrPed[nTwrMax];
  Float_t  TowerArray_TwrRMS[nTwrMax];
  Int_t    TowerArray_TwrMatchIdnex[nTwrMax];
  Int_t    TowerArray_NoOfmatchedTrk[nTwrMax];
  Float_t  TowerArray_TwrMatchP[nTwrMax];
  Float_t  TowerArray_TwrPx[nTwrMax];
  Float_t  TowerArray_TwrPy[nTwrMax];
  Float_t  TowerArray_TwrPz[nTwrMax];
  Int_t    TowerArray_fNAssocTracks[nTwrMax];
  Int_t    TowerArray_fMatchedTracksArray_[nTwrMax][nMatMax];
  Float_t  TowerArray_fMatchedTracksArray_P[nTwrMax][nMatMax];

  // list of branches
  TBranch *b_EventList_fUniqueID;
  TBranch *b_EventList_fBits;
  TBranch *b_EventList_runNumber;
  TBranch *b_EventList_eventNumber;
  TBranch *b_EventList_trigID;
  TBranch *b_EventList_nGlobalTracks;
  TBranch *b_EventList_nPrimaryTracks;
  TBranch *b_EventList_refMult;
  TBranch *b_EventList_vpdVz;
  TBranch *b_EventList_xVertex;
  TBranch *b_EventList_yVertex;
  TBranch *b_EventList_zVertex;
  TBranch *b_EventList_bbcZVertex;
  TBranch *b_EventList_zdcCoincidenceRate;
  TBranch *b_EventList_bbcCoincidenceRate;
  TBranch *b_EventList_backgroundRate;
  TBranch *b_EventList_bbcBlueBackgroundRate;
  TBranch *b_EventList_bbcYellowBackgroundRate;
  TBranch *b_EventList_refMultPos;
  TBranch *b_EventList_refMultNeg;
  TBranch *b_EventList_bTOFTrayMultiplicity;
  TBranch *b_EventList_nVerticies;
  TBranch *b_EventList_MagF;
  TBranch *b_EventList_VrtxRank;
  TBranch *b_EventList_Etsp;
  TBranch *b_EventList_ETwrdidT;
  TBranch *b_EventList_ETwradc11;
  TBranch *b_EventList_ETwreneT0;
  TBranch *b_EventList_ETwreT;
  TBranch *b_EventList_ETwrENET0;
  TBranch *b_EventList_ETwrphT;
  TBranch *b_EventList_ETwrPTower;
  TBranch *b_EventList_ETwrpidTower;
  TBranch *b_EventList_ETwrmoduleT;
  TBranch *b_EventList_EClustEneT0;
  TBranch *b_EventList_EClustetav1;
  TBranch *b_EventList_EClustphiv1;
  TBranch *b_EventList_EEstrpen01;
  TBranch *b_EventList_EEstrpen02;
  TBranch *b_EventList_EEstrpen03;
  TBranch *b_EventList_EEstrpen0;
  TBranch *b_EventList_EEstrpen1;
  TBranch *b_EventList_EEstrpen2;
  TBranch *b_EventList_EEstrpen3;
  TBranch *b_EventList_EEstrpen4;
  TBranch *b_EventList_EEstrpen5;
  TBranch *b_EventList_EEstrpen6;
  TBranch *b_EventList_EEstrpen7;
  TBranch *b_EventList_EEstrpen8;
  TBranch *b_EventList_EEstrpen9;
  TBranch *b_EventList_EEstrpen10;
  TBranch *b_EventList_EEstrpen11;
  TBranch *b_EventList_EEstrpen12;
  TBranch *b_EventList_EEstrpen13;
  TBranch *b_EventList_EEstrpen14;
  TBranch *b_EventList_EEstrpen15;
  TBranch *b_EventList_ETwrdidE;
  TBranch *b_EventList_EPstripenp01;
  TBranch *b_EventList_EPstripenp02;
  TBranch *b_EventList_EPstripenp03;
  TBranch *b_EventList_EPstripenp0;
  TBranch *b_EventList_EPstripenp1;
  TBranch *b_EventList_EPstripenp2;
  TBranch *b_EventList_EPstripenp3;
  TBranch *b_EventList_EPstripenp4;
  TBranch *b_EventList_EPstripenp5;
  TBranch *b_EventList_EPstripenp6;
  TBranch *b_EventList_EPstripenp7;
  TBranch *b_EventList_EPstripenp8;
  TBranch *b_EventList_EPstripenp9;
  TBranch *b_EventList_EPstripenp10;
  TBranch *b_EventList_EPstripenp11;
  TBranch *b_EventList_EPstripenp12;
  TBranch *b_EventList_EPstripenp13;
  TBranch *b_EventList_EPstripenp14;
  TBranch *b_EventList_EPstripenp15;
  TBranch *b_EventList_EclustEnnq1;
  TBranch *b_EventList_EclustEnnq20;
  TBranch *b_EventList_EclustEnnq19;
  TBranch *b_EventList_EclustEnpq1;
  TBranch *b_EventList_EclustEnpq20;
  TBranch *b_EventList_EclustEnpq19;
  TBranch *b_EventList_EclustEnpq21;
  TBranch *b_EventList_PrimaryTrackArray_;
  TBranch *b_PrimaryTrackArray_fUniqueID;
  TBranch *b_PrimaryTrackArray_fBits;
  TBranch *b_PrimaryTrackArray_nHitsFit;
  TBranch *b_PrimaryTrackArray_nHitsPoss;
  TBranch *b_PrimaryTrackArray_trackFlag;
  TBranch *b_PrimaryTrackArray_pZ;
  TBranch *b_PrimaryTrackArray_pX;
  TBranch *b_PrimaryTrackArray_pY;
  TBranch *b_PrimaryTrackArray_pT;
  TBranch *b_PrimaryTrackArray_dEdx;
  TBranch *b_PrimaryTrackArray_charge;
  TBranch *b_PrimaryTrackArray_tofBeta;
  TBranch *b_PrimaryTrackArray_eta;
  TBranch *b_PrimaryTrackArray_phi;
  TBranch *b_PrimaryTrackArray_nSigElectron;
  TBranch *b_PrimaryTrackArray_nSigPion;
  TBranch *b_PrimaryTrackArray_nSigKaon;
  TBranch *b_PrimaryTrackArray_nSigProton;
  TBranch *b_PrimaryTrackArray_dcag;
  TBranch *b_PrimaryTrackArray_nHits;
  TBranch *b_PrimaryTrackArray_dEdxHits;
  TBranch *b_PrimaryTrackArray_firstZPoint;
  TBranch *b_PrimaryTrackArray_lastZPoint;
  TBranch *b_PrimaryTrackArray_tofSigElectron;
  TBranch *b_PrimaryTrackArray_tofSigPion;
  TBranch *b_PrimaryTrackArray_tofSigKaon;
  TBranch *b_PrimaryTrackArray_tofSigProton;
  TBranch *b_PrimaryTrackArray_timeOfflight;
  TBranch *b_PrimaryTrackArray_pathLength;
  TBranch *b_PrimaryTrackArray_trkIndex;
  TBranch *b_EventList_TowerArray_;
  TBranch *b_TowerArray_fUniqueID;
  TBranch *b_TowerArray_fBits;
  TBranch *b_TowerArray_TwrId;
  TBranch *b_TowerArray_TwrEng;
  TBranch *b_TowerArray_TwrEta;
  TBranch *b_TowerArray_TwrPhi;
  TBranch *b_TowerArray_TwrADC;
  TBranch *b_TowerArray_TwrMatchIdnex;
  TBranch *b_TowerArray_NoOfmatchedTrk;
  TBranch *b_TowerArray_TwrMatchP;
  TBranch *b_TowerArray_TwrPx;
  TBranch *b_TowerArray_TwrPy;
  TBranch *b_TowerArray_TwrPz;
  TBranch *b_TowerArray_fNAssocTracks;
  TBranch *b_TowerArray_fMatchedTracksArray_;
  TBranch *b_TowerArray_fMatchedTracksArray_P;

  // set branch addresses
  tTree -> SetMakeClass(1);
  tTree -> SetBranchAddress("fUniqueID", &fUniqueID, &b_EventList_fUniqueID);
  tTree -> SetBranchAddress("fBits", &fBits, &b_EventList_fBits);
  tTree -> SetBranchAddress("runNumber", &runNumber, &b_EventList_runNumber);
  tTree -> SetBranchAddress("eventNumber", &eventNumber, &b_EventList_eventNumber);
  tTree -> SetBranchAddress("trigID", &trigID, &b_EventList_trigID);
  tTree -> SetBranchAddress("nGlobalTracks", &nGlobalTracks, &b_EventList_nGlobalTracks);
  tTree -> SetBranchAddress("nPrimaryTracks", &nPrimaryTracks, &b_EventList_nPrimaryTracks);
  tTree -> SetBranchAddress("refMult", &refMult, &b_EventList_refMult);
  tTree -> SetBranchAddress("vpdVz", &vpdVz, &b_EventList_vpdVz);
  tTree -> SetBranchAddress("xVertex", &xVertex, &b_EventList_xVertex);
  tTree -> SetBranchAddress("yVertex", &yVertex, &b_EventList_yVertex);
  tTree -> SetBranchAddress("zVertex", &zVertex, &b_EventList_zVertex);
  tTree -> SetBranchAddress("bbcZVertex", &bbcZVertex, &b_EventList_bbcZVertex);
  tTree -> SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate, &b_EventList_zdcCoincidenceRate);
  tTree -> SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate, &b_EventList_bbcCoincidenceRate);
  tTree -> SetBranchAddress("backgroundRate", &backgroundRate, &b_EventList_backgroundRate);
  tTree -> SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate, &b_EventList_bbcBlueBackgroundRate);
  tTree -> SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate, &b_EventList_bbcYellowBackgroundRate);
  tTree -> SetBranchAddress("refMultPos", &refMultPos, &b_EventList_refMultPos);
  tTree -> SetBranchAddress("refMultNeg", &refMultNeg, &b_EventList_refMultNeg);
  tTree -> SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity, &b_EventList_bTOFTrayMultiplicity);
  tTree -> SetBranchAddress("nVerticies", &nVerticies, &b_EventList_nVerticies);
  tTree -> SetBranchAddress("MagF", &MagF, &b_EventList_MagF);
  tTree -> SetBranchAddress("VrtxRank", &VrtxRank, &b_EventList_VrtxRank);
  tTree -> SetBranchAddress("Etsp", &Etsp, &b_EventList_Etsp);
  tTree -> SetBranchAddress("ETwrdidT", &ETwrdidT, &b_EventList_ETwrdidT);
  tTree -> SetBranchAddress("ETwradc11", &ETwradc11, &b_EventList_ETwradc11);
  tTree -> SetBranchAddress("ETwreneT0", &ETwreneT0, &b_EventList_ETwreneT0);
  tTree -> SetBranchAddress("ETwreT", &ETwreT, &b_EventList_ETwreT);
  tTree -> SetBranchAddress("ETwrENET0", &ETwrENET0, &b_EventList_ETwrENET0);
  tTree -> SetBranchAddress("ETwrphT", &ETwrphT, &b_EventList_ETwrphT);
  tTree -> SetBranchAddress("ETwrPTower", &ETwrPTower, &b_EventList_ETwrPTower);
  tTree -> SetBranchAddress("ETwrpidTower", &ETwrpidTower, &b_EventList_ETwrpidTower);
  tTree -> SetBranchAddress("ETwrmoduleT", &ETwrmoduleT, &b_EventList_ETwrmoduleT);
  tTree -> SetBranchAddress("EClustEneT0", &EClustEneT0, &b_EventList_EClustEneT0);
  tTree -> SetBranchAddress("EClustetav1", &EClustetav1, &b_EventList_EClustetav1);
  tTree -> SetBranchAddress("EClustphiv1", &EClustphiv1, &b_EventList_EClustphiv1);
  tTree -> SetBranchAddress("EEstrpen01", &EEstrpen01, &b_EventList_EEstrpen01);
  tTree -> SetBranchAddress("EEstrpen02", &EEstrpen02, &b_EventList_EEstrpen02);
  tTree -> SetBranchAddress("EEstrpen03", &EEstrpen03, &b_EventList_EEstrpen03);
  tTree -> SetBranchAddress("EEstrpen0", &EEstrpen0, &b_EventList_EEstrpen0);
  tTree -> SetBranchAddress("EEstrpen1", &EEstrpen1, &b_EventList_EEstrpen1);
  tTree -> SetBranchAddress("EEstrpen2", &EEstrpen2, &b_EventList_EEstrpen2);
  tTree -> SetBranchAddress("EEstrpen3", &EEstrpen3, &b_EventList_EEstrpen3);
  tTree -> SetBranchAddress("EEstrpen4", &EEstrpen4, &b_EventList_EEstrpen4);
  tTree -> SetBranchAddress("EEstrpen5", &EEstrpen5, &b_EventList_EEstrpen5);
  tTree -> SetBranchAddress("EEstrpen6", &EEstrpen6, &b_EventList_EEstrpen6);
  tTree -> SetBranchAddress("EEstrpen7", &EEstrpen7, &b_EventList_EEstrpen7);
  tTree -> SetBranchAddress("EEstrpen8", &EEstrpen8, &b_EventList_EEstrpen8);
  tTree -> SetBranchAddress("EEstrpen9", &EEstrpen9, &b_EventList_EEstrpen9);
  tTree -> SetBranchAddress("EEstrpen10", &EEstrpen10, &b_EventList_EEstrpen10);
  tTree -> SetBranchAddress("EEstrpen11", &EEstrpen11, &b_EventList_EEstrpen11);
  tTree -> SetBranchAddress("EEstrpen12", &EEstrpen12, &b_EventList_EEstrpen12);
  tTree -> SetBranchAddress("EEstrpen13", &EEstrpen13, &b_EventList_EEstrpen13);
  tTree -> SetBranchAddress("EEstrpen14", &EEstrpen14, &b_EventList_EEstrpen14);
  tTree -> SetBranchAddress("EEstrpen15", &EEstrpen15, &b_EventList_EEstrpen15);
  tTree -> SetBranchAddress("ETwrdidE", &ETwrdidE, &b_EventList_ETwrdidE);
  tTree -> SetBranchAddress("EPstripenp01", &EPstripenp01, &b_EventList_EPstripenp01);
  tTree -> SetBranchAddress("EPstripenp02", &EPstripenp02, &b_EventList_EPstripenp02);
  tTree -> SetBranchAddress("EPstripenp03", &EPstripenp03, &b_EventList_EPstripenp03);
  tTree -> SetBranchAddress("EPstripenp0", &EPstripenp0, &b_EventList_EPstripenp0);
  tTree -> SetBranchAddress("EPstripenp1", &EPstripenp1, &b_EventList_EPstripenp1);
  tTree -> SetBranchAddress("EPstripenp2", &EPstripenp2, &b_EventList_EPstripenp2);
  tTree -> SetBranchAddress("EPstripenp3", &EPstripenp3, &b_EventList_EPstripenp3);
  tTree -> SetBranchAddress("EPstripenp4", &EPstripenp4, &b_EventList_EPstripenp4);
  tTree -> SetBranchAddress("EPstripenp5", &EPstripenp5, &b_EventList_EPstripenp5);
  tTree -> SetBranchAddress("EPstripenp6", &EPstripenp6, &b_EventList_EPstripenp6);
  tTree -> SetBranchAddress("EPstripenp7", &EPstripenp7, &b_EventList_EPstripenp7);
  tTree -> SetBranchAddress("EPstripenp8", &EPstripenp8, &b_EventList_EPstripenp8);
  tTree -> SetBranchAddress("EPstripenp9", &EPstripenp9, &b_EventList_EPstripenp9);
  tTree -> SetBranchAddress("EPstripenp10", &EPstripenp10, &b_EventList_EPstripenp10);
  tTree -> SetBranchAddress("EPstripenp11", &EPstripenp11, &b_EventList_EPstripenp11);
  tTree -> SetBranchAddress("EPstripenp12", &EPstripenp12, &b_EventList_EPstripenp12);
  tTree -> SetBranchAddress("EPstripenp13", &EPstripenp13, &b_EventList_EPstripenp13);
  tTree -> SetBranchAddress("EPstripenp14", &EPstripenp14, &b_EventList_EPstripenp14);
  tTree -> SetBranchAddress("EPstripenp15", &EPstripenp15, &b_EventList_EPstripenp15);
  tTree -> SetBranchAddress("EclustEnnq1", &EclustEnnq1, &b_EventList_EclustEnnq1);
  tTree -> SetBranchAddress("EclustEnnq20", &EclustEnnq20, &b_EventList_EclustEnnq20);
  tTree -> SetBranchAddress("EclustEnnq19", &EclustEnnq19, &b_EventList_EclustEnnq19);
  tTree -> SetBranchAddress("EclustEnpq1", &EclustEnpq1, &b_EventList_EclustEnpq1);
  tTree -> SetBranchAddress("EclustEnpq20", &EclustEnpq20, &b_EventList_EclustEnpq20);
  tTree -> SetBranchAddress("EclustEnpq19", &EclustEnpq19, &b_EventList_EclustEnpq19);
  tTree -> SetBranchAddress("EclustEnpq21", &EclustEnpq21, &b_EventList_EclustEnpq21);
  tTree -> SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &b_EventList_PrimaryTrackArray_);
  tTree -> SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID, &b_PrimaryTrackArray_fUniqueID);
  tTree -> SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits, &b_PrimaryTrackArray_fBits);
  tTree -> SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit, &b_PrimaryTrackArray_nHitsFit);
  tTree -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss, &b_PrimaryTrackArray_nHitsPoss);
  tTree -> SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag, &b_PrimaryTrackArray_trackFlag);
  tTree -> SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ, &b_PrimaryTrackArray_pZ);
  tTree -> SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX, &b_PrimaryTrackArray_pX);
  tTree -> SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY, &b_PrimaryTrackArray_pY);
  tTree -> SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT, &b_PrimaryTrackArray_pT);
  tTree -> SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx, &b_PrimaryTrackArray_dEdx);
  tTree -> SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge, &b_PrimaryTrackArray_charge);
  tTree -> SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta, &b_PrimaryTrackArray_tofBeta);
  tTree -> SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta, &b_PrimaryTrackArray_eta);
  tTree -> SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi, &b_PrimaryTrackArray_phi);
  tTree -> SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron, &b_PrimaryTrackArray_nSigElectron);
  tTree -> SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion, &b_PrimaryTrackArray_nSigPion);
  tTree -> SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon, &b_PrimaryTrackArray_nSigKaon);
  tTree -> SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton, &b_PrimaryTrackArray_nSigProton);
  tTree -> SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag, &b_PrimaryTrackArray_dcag);
  tTree -> SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits, &b_PrimaryTrackArray_nHits);
  tTree -> SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits, &b_PrimaryTrackArray_dEdxHits);
  tTree -> SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint, &b_PrimaryTrackArray_firstZPoint);
  tTree -> SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint, &b_PrimaryTrackArray_lastZPoint);
  tTree -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron, &b_PrimaryTrackArray_tofSigElectron);
  tTree -> SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion, &b_PrimaryTrackArray_tofSigPion);
  tTree -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon, &b_PrimaryTrackArray_tofSigKaon);
  tTree -> SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton, &b_PrimaryTrackArray_tofSigProton);
  tTree -> SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight, &b_PrimaryTrackArray_timeOfflight);
  tTree -> SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength, &b_PrimaryTrackArray_pathLength);
  tTree -> SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex, &b_PrimaryTrackArray_trkIndex);
  tTree -> SetBranchAddress("TowerArray", &TowerArray_, &b_EventList_TowerArray_);
  tTree -> SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID, &b_TowerArray_fUniqueID);
  tTree -> SetBranchAddress("TowerArray.fBits", TowerArray_fBits, &b_TowerArray_fBits);
  tTree -> SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId, &b_TowerArray_TwrId);
  tTree -> SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng, &b_TowerArray_TwrEng);
  tTree -> SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta, &b_TowerArray_TwrEta);
  tTree -> SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi, &b_TowerArray_TwrPhi);
  tTree -> SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC, &b_TowerArray_TwrADC);
  tTree -> SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex, &b_TowerArray_TwrMatchIdnex);
  tTree -> SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk, &b_TowerArray_NoOfmatchedTrk);
  tTree -> SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP, &b_TowerArray_TwrMatchP);
  tTree -> SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx, &b_TowerArray_TwrPx);
  tTree -> SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy, &b_TowerArray_TwrPy);
  tTree -> SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz, &b_TowerArray_TwrPz);
  tTree -> SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks, &b_TowerArray_fNAssocTracks);
  tTree -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_, &b_TowerArray_fMatchedTracksArray_);
  tTree -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P, &b_TowerArray_fMatchedTracksArray_P);

  cout << "    Branches set..." << endl;


  // open stream
  ofstream out(sOut.Data());
  if (!out) {
    cerr << "PANIC: couldn't open output stream!" << endl;
    assert(out);
  }


  Int_t nEvts = tTree -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " to process..." << endl;

  // event loop
  Int_t nByte = 0;
  Int_t nTrgs = 0;
  for (Int_t i = 0; i < nEvts; i++) {

    nByte += tTree -> GetEntry(i);
    if (inBatchMode) {
      cout << "      Processing event " << (i + 1) << "/" << nEvts << "..." << endl;
    }
    else {
      cout << "      Processing event " << (i + 1) << "/" << nEvts << "...\r" << flush;
      if ((i + 1) == nEvts) cout << endl;
    }

    TString sFile("");
    TString sRun("");
    TString sEvt("");
    sFile += tTree -> GetCurrentFile() -> GetName();
    sRun  += runNumber;
    sEvt  += eventNumber;
    out << sFile.Data() << " " << sRun.Data() << " " << sEvt.Data() << endl;


  }  // end event loop

  cout << "    Event loop finished!" << endl;


  fIn -> cd();
  fIn -> Close();
  cout << "  Weird list creator finished!\n" << endl;

}

// End ------------------------------------------------------------------------
