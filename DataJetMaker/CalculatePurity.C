// 'CalculatePurity.C'
// Derek Anderson
// 01.18.2018
//
// Use this script to calculate the purity
// of a trigger sample from the output
// of the 'StThirdJetMaker'
//
// NOTE: the 1st eTtrg 'bin' is defined
//       to be the whole range of eTtrg.


#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t NTrkMax(5000);
static const UInt_t NTwrMax(5000);
static const UInt_t NMatchMax(10);
//static const UInt_t NHotTwr(41);
static const UInt_t NHotTwr(302);
static const UInt_t NBadRuns(45);
static const UInt_t NBins(4);
static const UInt_t NHist(NBins - 1);
static const UInt_t NTrgs(2);



void CalculatePurity(const Bool_t isInBatchMode=false) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning purity calculation..." << endl;


  // io parameters
  const TString sOutput("pp200r9purity.forThesis_betterStyle.et920vz55pt1220dca1zt0203df063.d24m9y2021.root");
  const TString sInput("../FullJetTree/merged/pp200r9.merge.root");

  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);

  // trigger paramters
  const Int_t    adcMax(6004);
  const Double_t eStrMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t hTrgMax(0.9);

  // track parameters
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(1.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(1.2);
  const Double_t pTtrkMax(20.);
  const Double_t zTtrkMin(0.2);
  const Double_t zTtrkMax(0.3);
  const TString  sZtTrkMin("0.2");
  const TString  sZtTrkMax("0.3");

  // calculation parameters
  const Double_t dFnear(0.);
  const Double_t dFaway(TMath::Pi());
  const Double_t dFsize(0.63);
  const Double_t sigNear(0.02);
  const Double_t sigAway(0.035);



  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input tree
  TTree *tInput;
  fInput -> GetObject("Gfmtodst", tInput);
  if (!tInput) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    return;
  }
  cout << "    Grabbed tree." << endl;


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
  cout << "    Set branches." << endl;


  // define bad run and hot tower lists
  const UInt_t badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  //const UInt_t hotTwrList[NHotTwr] = {1, 35, 141, 187, 224, 341, 424, 594, 814, 899, 900, 1046, 1128, 1132, 1244, 1382, 1388, 1405, 1588, 1766, 1773, 2066, 2160, 2253, 2281, 2284, 2301, 2303, 2306, 2590, 3007, 3495, 3840, 4043, 4047, 4053, 4057, 4121, 4442, 4569, 4617};
  const UInt_t hotTwrList[NHotTwr] = {34, 106, 113, 160, 266, 267, 275, 280, 282, 286, 287, 293, 410, 504, 533, 541, 555, 561, 562, 594, 615, 616, 629, 633, 637, 638, 647, 650, 653, 657, 671, 673, 743, 789, 790, 791, 792, 806, 809, 810, 811, 812, 813, 814, 821, 822, 823, 824, 829, 830, 831, 832, 837, 841, 842, 843, 844, 846, 849, 850, 851, 852, 857, 875, 897, 899, 903, 939, 953, 954, 956, 993, 1026, 1046, 1048, 1080, 1081, 1100, 1125, 1130, 1132, 1180, 1197, 1198, 1199, 1200, 1207, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1312, 1348, 1353, 1354, 1388, 1407, 1409, 1434, 1448, 1537, 1567, 1574, 1597, 1612, 1654, 1668, 1713, 1762, 1765, 1766, 1877, 1878, 1984, 2032, 2043, 2054, 2073, 2077, 2092, 2093, 2097, 2107, 2162, 2168, 2214, 2305, 2392, 2409, 2415, 2439, 2459, 2589, 2590, 2633, 2652, 2749, 2834, 2961, 2969, 3005, 3017, 3070, 3071, 3186, 3220, 3289, 3360, 3493, 3494, 3495, 3508, 3588, 3604, 3611, 3668, 3678, 3679, 3690, 3692, 3732, 3738, 3838, 3840, 3927, 3945, 4005, 4006, 4013, 4018, 4019, 4053, 4059, 4124, 4331, 4355, 4357, 4458, 4464, 4500, 4677, 4678, 4684, 4768, 360, 493, 779, 1284, 1306, 1337, 1438, 1709, 2027, 2445, 3407, 3720, 4217, 4288, 95, 96, 296, 316, 443, 479, 555, 562, 637, 671, 709, 740, 743, 796, 857, 897, 899, 915, 953, 1130, 1132, 1294, 1318, 1337, 1348, 1359, 1378, 1427, 1429, 1440, 1537, 1563, 1574, 1709, 1763, 1773, 1819, 1854, 1874, 1936, 1938, 2018, 2043, 2098, 2099, 2256, 2259, 2294, 2514, 2520, 2552, 2589, 2598, 2680, 2706, 2799, 2880, 2897, 2917, 2969, 3020, 3028, 3310, 3319, 3375, 3399, 3504, 3539, 3541, 3679, 3690, 3692, 3718, 3719, 3720, 3738, 3806, 3838, 3840, 3928, 4013, 4017, 4038, 4053, 4057, 4058, 4079, 4097, 4099};
  cout << "    Bad run and hot tower lists defined:\n"
       << "      " << NBadRuns << " bad runs, " << NHotTwr << " hot towers."
       << endl;


  // define trigger bins
  const Double_t eTtrgMin[NBins] = {9., 9., 11., 15.};
  const Double_t eTtrgMax[NBins] = {20., 11., 15., 20.};
  const Double_t eTtrgBin[NBins] = {9., 11., 15., 20.};
  const Double_t tspMin[NTrgs]   = {0., 0.2};
  const Double_t tspMax[NTrgs]   = {0.08, 0.6};

  cout << "    Trigger bins:" << endl;
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    cout << "      Bin [" << iBin << "] = (" << eTtrgMin[iBin] << ", " << eTtrgMax[iBin] << ") GeV" << endl;
  }


  // define histograms
  TH1D *hDeltaPhi[NBins][NTrgs];
  TH1D *hDeltaPhiNS[NBins][NTrgs];
  TH1D *hDeltaPhiAS[NBins][NTrgs];
  TH1D *hEtaTrk[NBins][NTrgs];
  TH1D *hPtTrk[NBins][NTrgs];
  TH1D *hZtTrk[NBins][NTrgs];
  TH1D *hZtCalc[NBins][NTrgs];

  // binning
  const UInt_t   nDf   = 60;
  const UInt_t   nH    = 40;
  const UInt_t   nPt   = 40;
  const UInt_t   nZt   = 10;
  const Double_t dF[2] = {-1.57, 4.71};
  const Double_t h[2]  = {-1., 1.};
  const Double_t pT[2] = {0., 40.};
  const Double_t zT[2] = {0., 1.};

  // names
  const TString sDeltaPhi("hDeltaPhi_");
  const TString sDeltaPhiNS("hDeltaPhiNS_");
  const TString sDeltaPhiAS("hDeltaPhiAS_");
  const TString sEtaTrk("hEtaTrk_");
  const TString sPtTrk("hPtTrk_");
  const TString sZtTrk("hZtTrk_");
  const TString sZtCalc("hZtCalc_");
  const TString sTrgNames[NTrgs] = {"pi", "ga"};
  const TString sEtNames[NBins]  = {"920", "911", "1115", "1520"};

  fOutput -> cd();
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {

      // add bin names
      TString sDeltaPhiName(sDeltaPhi.Data());
      TString sDeltaPhiNsName(sDeltaPhiNS.Data());
      TString sDeltaPhiAsName(sDeltaPhiAS.Data());
      TString sEtaTrkName(sEtaTrk.Data());
      TString sPtTrkName(sPtTrk.Data());
      TString sZtTrkName(sZtTrk.Data());
      TString sZtCalcName(sZtCalc.Data());
      sDeltaPhiName   += sTrgNames[iTrg].Data();
      sDeltaPhiName   += sEtNames[iBin].Data();
      sDeltaPhiNsName += sTrgNames[iTrg].Data();
      sDeltaPhiNsName += sEtNames[iBin].Data();
      sDeltaPhiAsName += sTrgNames[iTrg].Data();
      sDeltaPhiAsName += sEtNames[iBin].Data();
      sEtaTrkName     += sTrgNames[iTrg].Data();
      sEtaTrkName     += sEtNames[iBin].Data();
      sPtTrkName      += sTrgNames[iTrg].Data();
      sPtTrkName      += sEtNames[iBin].Data();
      sZtTrkName      += sTrgNames[iTrg].Data();
      sZtTrkName      += sEtNames[iBin].Data();
      sZtCalcName     += sTrgNames[iTrg].Data();
      sZtCalcName     += sEtNames[iBin].Data();

      // create histograms
      hDeltaPhi[iBin][iTrg]   = new TH1D(sDeltaPhiName.Data(), "", nDf, dF[0], dF[1]);
      hDeltaPhiNS[iBin][iTrg] = new TH1D(sDeltaPhiNsName.Data(), "", nDf, dF[0], dF[1]);
      hDeltaPhiAS[iBin][iTrg] = new TH1D(sDeltaPhiAsName.Data(), "", nDf, dF[0], dF[1]);
      hEtaTrk[iBin][iTrg]     = new TH1D(sEtaTrkName.Data(), "", nH, h[0], h[1]);
      hPtTrk[iBin][iTrg]      = new TH1D(sPtTrkName.Data(), "", nPt, pT[0], pT[1]);
      hZtTrk[iBin][iTrg]      = new TH1D(sZtTrkName.Data(), "", nZt, zT[0], zT[1]);
      hZtCalc[iBin][iTrg]     = new TH1D(sZtCalcName.Data(), "", nZt, zT[0], zT[1]);
      hDeltaPhi[iBin][iTrg]   -> Sumw2();
      hDeltaPhiNS[iBin][iTrg] -> Sumw2();
      hDeltaPhiAS[iBin][iTrg] -> Sumw2();
      hEtaTrk[iBin][iTrg]     -> Sumw2();
      hPtTrk[iBin][iTrg]      -> Sumw2();
      hZtTrk[iBin][iTrg]      -> Sumw2();
      hZtCalc[iBin][iTrg]     -> Sumw2();

    }  // end trigger id loop
  }  // end trigger eT loop
  cout << "    Made histograms." << endl;


  const UInt_t nEvts = tInput -> GetEntriesFast();
  cout << "    Beginning event loop: " << nEvts << " events to process." << endl;

  // no. of triggers
  UInt_t nTrgBin[NBins][NTrgs];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      nTrgBin[iBin][iTrg] = 0;
    }
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
    }
    else {
      if (isInBatchMode) {
        cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
      }
      else {
        cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
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

    Bool_t isGoodRun = true;
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (run == badRunList[iRun]) {
        isGoodRun = false;
        break;
      }
    }

    // event cuts
    const Bool_t isInRcut = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut = (TMath::Abs(zVtx) < zVtxMax);
    if (!isGoodRun || !isInRcut || !isInZcut) continue;


    // trigger info
    const Int_t    adc   = ETwradc11;
    const UInt_t   idTrg = ETwrdidT;
    const Double_t tsp   = Etsp;
    const Double_t eH4   = EEstrpen4;
    const Double_t eF4   = EPstripenp4;
    const Double_t pProj = ETwrPTower;
    const Double_t hDet  = ETwreT;
    const Double_t hPhys = EClustetav1;
    const Double_t fTrg  = EClustphiv1;
    const Double_t eTrg  = EClustEneT0;
    const Double_t tTrg  = 2. * TMath::ATan(TMath::Exp(-1. * hPhys));
    const Double_t eTtrg = eTrg * TMath::Sin(tTrg);

    Bool_t isGoodTwr = true;
    for (UInt_t iTwr = 0; iTwr < NHotTwr; iTwr++) {
      if (idTrg == hotTwrList[iTwr]) {
        isGoodTwr = false;
        break;
      }
    }

    // trigger cuts
    const Bool_t isInAdcCut    = (adc <= adcMax);
    const Bool_t isInStrCut    = ((eH4 >= eStrMin) && (eF4 >= eStrMin));
    const Bool_t isInProjCut   = (pProj < pProjMax);
    const Bool_t isInEtaTrgCut = (TMath::Abs(hDet) < hTrgMax);
    const Bool_t isInEtCut     = ((eTtrg > eTtrgMin[0]) && (eTtrg < eTtrgMax[0]));
    const Bool_t isInPi0cut    = ((tsp > tspMin[0]) && (tsp < tspMax[0]));
    const Bool_t isInGamCut    = ((tsp > tspMin[1]) && (tsp < tspMax[1]));
    const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);
    if (!isGoodTwr || !isInAdcCut || !isInStrCut || !isInProjCut || !isInEtaTrgCut || !isInEtCut || !isInTspCut) continue;


    // determine bin
    UInt_t eTbin = 0;
    for (UInt_t iBin = 1; iBin < NBins; iBin++) {
      const Bool_t isInBin = ((eTtrg >= eTtrgMin[iBin]) && (eTtrg < eTtrgMax[iBin]));
      if (isInBin) {
        eTbin = iBin;
        break;
      }
    }
    if (isInPi0cut) {
      nTrgBin[0][0]++;
      nTrgBin[eTbin][0]++;
    }
    if (isInGamCut) {
      nTrgBin[0][1]++;
      nTrgBin[eTbin][1]++;
    }


    // track loop
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFit  = PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPoss = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFit  = (Double_t) nFit / (Double_t) nPoss;
      const Double_t dca   = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk  = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk  = PrimaryTrackArray_phi[iTrk];
      const Double_t pTtrk = PrimaryTrackArray_pT[iTrk];
      const Double_t zTtrk = pTtrk / eTtrg;

      Double_t dFtrk = fTrk - fTrg;
      if (dFtrk < (-1. * TMath::PiOver2())) dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))  dFtrk -= TMath::TwoPi();

      // track cuts
      const Bool_t isInFitCut    = (nFit > nFitMin);
      const Bool_t isInRatioCut  = (rFit > rFitMin);
      const Bool_t isInDcaCut    = (dca < dcaMax);
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtCut     = ((pTtrk > pTtrkMin) && (pTtrk < pTtrkMax));
      if (!isInFitCut || !isInRatioCut || !isInDcaCut || !isInEtaTrkCut || !isInPtCut) continue;

      // TEST [05.07.2020]
      const Double_t dEdxTrk     = PrimaryTrackArray_dEdx[iTrk];
      const Bool_t   isInDeDxCut = (dEdxTrk != 0.);
      if (!isInDeDxCut) continue;


      // fill histograms
      if (isInPi0cut) {
        hEtaTrk[0][0]     -> Fill(hTrk);
        hPtTrk[0][0]      -> Fill(pTtrk);
        hZtTrk[0][0]      -> Fill(zTtrk);
        hEtaTrk[eTbin][0] -> Fill(hTrk);
        hPtTrk[eTbin][0]  -> Fill(pTtrk);
        hZtTrk[eTbin][0]  -> Fill(zTtrk);
      }
      if (isInGamCut) {
        hEtaTrk[0][1]     -> Fill(hTrk);
        hPtTrk[0][1]      -> Fill(pTtrk);
        hZtTrk[0][1]      -> Fill(zTtrk);
        hEtaTrk[eTbin][1] -> Fill(hTrk);
        hPtTrk[eTbin][1]  -> Fill(pTtrk);
        hZtTrk[eTbin][1]  -> Fill(zTtrk);
      }

      const Bool_t isInZtCut  = ((zTtrk > zTtrkMin) && (zTtrk < zTtrkMax));
      const Bool_t isNearSide = (TMath::Abs(dFtrk - dFnear) < dFsize);
      const Bool_t isAwaySide = (TMath::Abs(dFtrk - dFaway) < dFsize);
      if (isInZtCut) {
        if (isInPi0cut) {
          hDeltaPhi[0][0]     -> Fill(dFtrk);
          hZtCalc[0][0]       -> Fill(zTtrk);
          hDeltaPhi[eTbin][0] -> Fill(dFtrk);
          hZtCalc[eTbin][0]   -> Fill(zTtrk);
          if (isNearSide) {
            hDeltaPhiNS[0][0]     -> Fill(dFtrk);
            hDeltaPhiNS[eTbin][0] -> Fill(dFtrk);
          }  // end NS cut
          if (isAwaySide) {
            hDeltaPhiAS[0][0]     -> Fill(dFtrk);
            hDeltaPhiAS[eTbin][0] -> Fill(dFtrk);
          }  // end AS cut
        }  // end pi0 cut
        if (isInGamCut) {
          hDeltaPhi[0][1]     -> Fill(dFtrk);
          hZtCalc[0][1]       -> Fill(zTtrk);
          hDeltaPhi[eTbin][1] -> Fill(dFtrk);
          hZtCalc[eTbin][1]   -> Fill(zTtrk);
          if (isNearSide) {
            hDeltaPhiNS[0][1]     -> Fill(dFtrk);
            hDeltaPhiNS[eTbin][1] -> Fill(dFtrk);
          }  // end NS cut
          if (isAwaySide) {
            hDeltaPhiAS[0][1]     -> Fill(dFtrk);
            hDeltaPhiAS[eTbin][1] -> Fill(dFtrk);
          }  // end AS cut
        }  // end pi0 cut
      }  // end zT cut

    }  // end track loop

  }  // end event loop

  cout << "    Event loop finished:\n"
       << "      " << nTrgBin[0][0] << " pi0 triggers in total,\n"
       << "      " << nTrgBin[0][1] << " gamma triggers in total."
       << endl;


   // normalize histograms
   for (UInt_t iBin = 0; iBin < NBins; iBin++) {
     for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
       const Double_t dFbin  = hDeltaPhi[iBin][iTrg] -> GetBinWidth(17);
       const Double_t hBin   = hEtaTrk[iBin][iTrg]   -> GetBinWidth(17);
       const Double_t pTbin  = hPtTrk[iBin][iTrg]    -> GetBinWidth(17);
       const Double_t zTbin  = hZtTrk[iBin][iTrg]    -> GetBinWidth(17);
       const Double_t dFnorm = dFbin * nTrgBin[iBin][iTrg];
       const Double_t hNorm  = hBin * nTrgBin[iBin][iTrg];
       const Double_t pTnorm = pTbin * nTrgBin[iBin][iTrg];
       const Double_t zTnorm = zTbin * nTrgBin[iBin][iTrg];
       hDeltaPhi[iBin][iTrg]   -> Scale(1. / dFnorm);
       hDeltaPhiNS[iBin][iTrg] -> Scale(1. / dFnorm);
       hDeltaPhiAS[iBin][iTrg] -> Scale(1. / dFnorm);
       hEtaTrk[iBin][iTrg]     -> Scale(1. / hNorm);
       hPtTrk[iBin][iTrg]      -> Scale(1. / pTnorm);
       hZtTrk[iBin][iTrg]      -> Scale(1. / zTnorm);
       hZtCalc[iBin][iTrg]     -> Scale(1. / zTnorm);
     }  // end trigger id loop
   }  // end trigger eT loop
   cout << "    Normalized histograms." << endl;


  // fit correlations
  const UInt_t   fLinL(1);
  const UInt_t   fSizL(4);
  const UInt_t   fColL[NTrgs] = {859, 899};
  const Double_t fGausMin(0.);
  const Double_t fSigma[2]    = {0., dFsize};
  const Double_t fMuNear[2]   = {dFnear - dFsize, dFnear + dFsize};
  const Double_t fMuAway[2]   = {dFaway - dFsize, dFaway + dFaway};
  const Double_t fPedRange[2] = {0.2, 10.};

  TF1     *fCorrs[NBins][NTrgs];
  Double_t pedestal[NBins][NTrgs];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {

      // function name
      TString sFunc("fDeltaPhi_et");
      sFunc += iBin;
      sFunc += "id";
      sFunc += iTrg;

      // get max and estimate pedestal
      const UInt_t   iPedLow  = hDeltaPhi[iBin][iTrg] -> FindBin(dFnear + dFsize);
      const UInt_t   iPedHigh = hDeltaPhi[iBin][iTrg] -> FindBin(dFaway - dFsize);
      const Double_t maximum  = hDeltaPhi[iBin][iTrg] -> GetMaximum();
      const Double_t integral = hDeltaPhi[iBin][iTrg] -> Integral(iPedLow, iPedHigh);
      const Double_t estimate = integral / (iPedHigh - iPedLow);
      const Double_t pedMin   = fPedRange[0] * estimate;
      const Double_t pedMax   = fPedRange[1] * estimate;

      // define fit
      fCorrs[iBin][iTrg] = new TF1(sFunc.Data(), "gaus(0)+gaus(3)+pol0(6)", dF[0], dF[1]);
      fCorrs[iBin][iTrg] -> SetParameters(maximum, dFnear, sigNear, maximum, dFaway, sigAway, estimate);
      fCorrs[iBin][iTrg] -> SetParLimits(0, fGausMin, maximum);
      fCorrs[iBin][iTrg] -> SetParLimits(1, fMuNear[0], fMuNear[1]);
      fCorrs[iBin][iTrg] -> SetParLimits(2, fSigma[0], fSigma[1]);
      fCorrs[iBin][iTrg] -> SetParLimits(3, fGausMin, maximum);
      fCorrs[iBin][iTrg] -> SetParLimits(4, fMuAway[0], fMuAway[1]);
      fCorrs[iBin][iTrg] -> SetParLimits(6, pedMin, pedMax);
      fCorrs[iBin][iTrg] -> SetLineColor(fColL[iTrg]);
      fCorrs[iBin][iTrg] -> SetLineStyle(fLinL);
      fCorrs[iBin][iTrg] -> SetLineWidth(fSizL);

      // fit histograms
      hDeltaPhi[iBin][iTrg] -> Fit(fCorrs[iBin][iTrg], "RBMQ0");
      pedestal[iBin][iTrg] = fCorrs[iBin][iTrg] -> GetParameter(6);

      // reset visibility
      const Int_t kNotDraw = 1<<9;
      if (hDeltaPhi[iBin][iTrg] -> GetFunction(sFunc.Data())) {
        hDeltaPhi[iBin][iTrg] -> GetFunction(sFunc.Data()) -> ResetBit(kNotDraw);
      }

    }  // end trigger id loop
  }  // end trigger eT loop
  cout << "    Fit correlations." << endl;
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    cout << "      Bin [" << iBin << "]: pi0 pedestal = " << pedestal[iBin][0]
         << ", gamma pedestal = " << pedestal[iBin][1]
         << endl;
  }


  // calculate yields
  Double_t content(0.);
  Double_t corrected(0.);
  Double_t error(0.);
  Double_t yieldNS[NBins][NTrgs];
  Double_t errorNS[NBins][NTrgs];
  Double_t yieldAS[NBins][NTrgs];
  Double_t errorAS[NBins][NTrgs];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      yieldNS[iBin][iTrg] = 0.;
      errorNS[iBin][iTrg] = 0.;
      yieldAS[iBin][iTrg] = 0.;
      errorAS[iBin][iTrg] = 0.;
    }
  }

  // near-side yields
  content   = 0.;
  corrected = 0.;
  error     = 0.;
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      const Double_t dFlowNS   = dFnear - dFsize;
      const Double_t dFhighNS  = dFnear + dFsize;
      const UInt_t   lowBinNS  = hDeltaPhi[iBin][iTrg] -> FindBin(dFlowNS);
      const UInt_t   highBinNS = hDeltaPhi[iBin][iTrg] -> FindBin(dFhighNS);
      for (UInt_t iPhiNS = lowBinNS; iPhiNS <= highBinNS; iPhiNS++) {
        content   = hDeltaPhi[iBin][iTrg] -> GetBinContent(iPhiNS);
        error     = hDeltaPhi[iBin][iTrg] -> GetBinError(iPhiNS);
        corrected = content - pedestal[iBin][iTrg];
        if (content > 0.) {
          yieldNS[iBin][iTrg] += corrected;
          errorNS[iBin][iTrg] += (error * error);
        }
      }  // end phi loop
      errorNS[iBin][iTrg] = TMath::Sqrt(errorNS[iBin][iTrg]);
    }  // end trigger id loop
  }  // end trigger eT loop

  // away-side yields
  content   = 0.;
  corrected = 0.;
  error     = 0.;
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      const Double_t dFlowAS   = dFaway - dFsize;
      const Double_t dFhighAS  = dFaway + dFsize;
      const UInt_t   lowBinAS  = hDeltaPhi[iBin][iTrg] -> FindBin(dFlowAS);
      const UInt_t   highBinAS = hDeltaPhi[iBin][iTrg] -> FindBin(dFhighAS);
      for (UInt_t iPhiAS = lowBinAS; iPhiAS <= highBinAS; iPhiAS++) {
        content   = hDeltaPhi[iBin][iTrg] -> GetBinContent(iPhiAS);
        error     = hDeltaPhi[iBin][iTrg] -> GetBinError(iPhiAS);
        corrected = content - pedestal[iBin][iTrg];
        if (content > 0.) {
          yieldAS[iBin][iTrg] += corrected;
          errorAS[iBin][iTrg] += (error * error);
        }
      }  // end phi loop
      errorAS[iBin][iTrg] = TMath::Sqrt(errorAS[iBin][iTrg]);
    }  // end trigger id loop
  }  // end trigger eT loop
  cout << "    Calculated yields." << endl;
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    cout << "      Bin[" << iBin << "]: pi0(NS) = " << yieldNS[iBin][0] << " +- " << errorNS[iBin][0] << ",\tgam(NS) = " << yieldNS[iBin][1] << " +- " << errorNS[iBin][1] << "\n"
         << "              pi0(AS) = " << yieldAS[iBin][0] << " +- " << errorAS[iBin][0] << ",\tgam(AS) = " << yieldAS[iBin][1] << " +- " << errorAS[iBin][1]
         << endl;
  }

  // make yield histograms
  const UInt_t fColY[NTrgs]  = {859, 899};
  const UInt_t fStyYV[NTrgs] = {10, 10};
  const UInt_t fStyYE[NTrgs] = {1, 1};
  const UInt_t fSizYV[NTrgs] = {1, 1};
  const UInt_t fSizYE[NTrgs] = {2, 2};

  TH1D *hYieldNS[NTrgs];
  TH1D *hYieldAS[NTrgs];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    TString sYieldNsName("sYieldNS_");
    TString sYieldAsName("sYieldAS_");
    sYieldNsName += sTrgNames[iTrg].Data();
    sYieldAsName += sTrgNames[iTrg].Data();

    hYieldNS[iTrg] = new TH1D(sYieldNsName.Data(), "", NHist, eTtrgBin);
    hYieldAS[iTrg] = new TH1D(sYieldAsName.Data(), "", NHist, eTtrgBin);
    hYieldNS[iTrg] -> Sumw2();
    hYieldAS[iTrg] -> Sumw2();
    for (UInt_t iHist = 1; iHist < NHist + 1; iHist++) {
      hYieldNS[iTrg] -> SetBinContent(iHist, yieldNS[iHist][iTrg]);
      hYieldAS[iTrg] -> SetBinContent(iHist, yieldAS[iHist][iTrg]);
      hYieldNS[iTrg] -> SetBinError(iHist, errorNS[iHist][iTrg]);
      hYieldAS[iTrg] -> SetBinError(iHist, errorAS[iHist][iTrg]);
    }  // end trigger eT loop
  }  // end trigger id loop

  TLine *lYieldNS[NTrgs];
  TLine *lYieldAS[NTrgs];
  TLine *lYieldNsHi[NTrgs];
  TLine *lYieldAsHi[NTrgs];
  TLine *lYieldNsLo[NTrgs];
  TLine *lYieldAsLo[NTrgs];
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    lYieldNS[iTrg]   = new TLine(eTtrgBin[0], yieldNS[0][iTrg], eTtrgBin[NHist], yieldNS[0][iTrg]);
    lYieldAS[iTrg]   = new TLine(eTtrgBin[0], yieldAS[0][iTrg], eTtrgBin[NHist], yieldAS[0][iTrg]);
    lYieldNsHi[iTrg] = new TLine(eTtrgBin[0], yieldNS[0][iTrg] + errorNS[0][iTrg], eTtrgBin[NHist], yieldNS[0][iTrg] + errorNS[0][iTrg]);
    lYieldAsHi[iTrg] = new TLine(eTtrgBin[0], yieldAS[0][iTrg] + errorAS[0][iTrg], eTtrgBin[NHist], yieldAS[0][iTrg] + errorAS[0][iTrg]);
    lYieldNsLo[iTrg] = new TLine(eTtrgBin[0], yieldNS[0][iTrg] - errorNS[0][iTrg], eTtrgBin[NHist], yieldNS[0][iTrg] - errorNS[0][iTrg]);
    lYieldAsLo[iTrg] = new TLine(eTtrgBin[0], yieldAS[0][iTrg] - errorAS[0][iTrg], eTtrgBin[NHist], yieldAS[0][iTrg] - errorAS[0][iTrg]);
    lYieldNS[iTrg]   -> SetLineColor(fColY[iTrg]);
    lYieldNS[iTrg]   -> SetLineStyle(fStyYV[iTrg]);
    lYieldNS[iTrg]   -> SetLineWidth(fSizYV[iTrg]);
    lYieldAS[iTrg]   -> SetLineColor(fColY[iTrg]);
    lYieldAS[iTrg]   -> SetLineStyle(fStyYV[iTrg]);
    lYieldAS[iTrg]   -> SetLineWidth(fSizYV[iTrg]);
    lYieldNsHi[iTrg] -> SetLineColor(fColY[iTrg]);
    lYieldNsHi[iTrg] -> SetLineStyle(fStyYE[iTrg]);
    lYieldNsHi[iTrg] -> SetLineWidth(fSizYE[iTrg]);
    lYieldAsHi[iTrg] -> SetLineColor(fColY[iTrg]);
    lYieldAsHi[iTrg] -> SetLineStyle(fStyYE[iTrg]);
    lYieldAsHi[iTrg] -> SetLineWidth(fSizYE[iTrg]);
    lYieldNsLo[iTrg] -> SetLineColor(fColY[iTrg]);
    lYieldNsLo[iTrg] -> SetLineStyle(fStyYE[iTrg]);
    lYieldNsLo[iTrg] -> SetLineWidth(fSizYE[iTrg]);
    lYieldAsLo[iTrg] -> SetLineColor(fColY[iTrg]);
    lYieldAsLo[iTrg] -> SetLineStyle(fStyYE[iTrg]);
    lYieldAsLo[iTrg] -> SetLineWidth(fSizYE[iTrg]);
  }
  cout << "    Made yield histograms." << endl;


  // calculate R
  Double_t value(0.);
  Double_t errorSum(0.);
  Double_t errorVal(0.);
  Double_t errorTrg[NTrgs];
  Double_t valueR[NBins];
  Double_t errorR[NBins];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    if ((yieldNS[iBin][0] > 0.) && (yieldNS[iBin][1] > 0.)) {
      value        = yieldNS[iBin][1] / yieldNS[iBin][0];
      errorTrg[0]  = errorNS[iBin][0] / yieldNS[iBin][0];
      errorTrg[1]  = errorNS[iBin][1] / yieldNS[iBin][1];
      errorSum     = (errorTrg[0] * errorTrg[0]) + (errorTrg[1] * errorTrg[1]);
      errorVal     = value * TMath::Sqrt(errorSum);
      valueR[iBin] = value;
      errorR[iBin] = errorVal;
    }
    else {
      valueR[iBin] = 0.;
      errorR[iBin] = 0.;
    }
  }  // end trigger eT loop
  cout << "    Calculated R." << endl;
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    cout << "      Bin [" << iBin << "]: R = " << valueR[iBin] << endl;
  }

  // make R histogram
  const UInt_t fColR(879);
  const UInt_t fStyRV(10);
  const UInt_t fStyRE(1);
  const UInt_t fSizRV(1);
  const UInt_t fSizRE(2);

  TH1D *hPurity = new TH1D("hPurity", "", NHist, eTtrgBin);
  hPurity -> Sumw2();
  for (UInt_t iHist = 1; iHist < NHist + 1; iHist++) {
    hPurity -> SetBinContent(iHist, valueR[iHist]);
    hPurity -> SetBinError(iHist, errorR[iHist]);
  }

  TLine *lAverage = new TLine(eTtrgBin[0], valueR[0], eTtrgBin[NHist], valueR[0]);
  TLine *lAvgUp   = new TLine(eTtrgBin[0], valueR[0] + errorR[0], eTtrgBin[NHist], valueR[0] + errorR[0]);
  TLine *lAvgDown = new TLine(eTtrgBin[0], valueR[0] - errorR[0], eTtrgBin[NHist], valueR[0] - errorR[0]);
  lAverage -> SetLineColor(fColR);
  lAverage -> SetLineStyle(fStyRV);
  lAverage -> SetLineWidth(fSizRV);
  lAvgUp   -> SetLineColor(fColR);
  lAvgUp   -> SetLineStyle(fStyRE);
  lAvgUp   -> SetLineWidth(fSizRE);
  lAvgDown -> SetLineColor(fColR);
  lAvgDown -> SetLineStyle(fStyRE);
  lAvgDown -> SetLineWidth(fSizRE);
  cout << "    Made R histogram." << endl;


  // set styles
  const UInt_t  fMar(20);
  const UInt_t  fFil(0);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fFilNS(3345);
  const UInt_t  fLinNS(2);
  const UInt_t  fMarNS(1);
  const UInt_t  fFilAS(3354);
  const UInt_t  fLinAS(2);
  const UInt_t  fMarAS(2);
  const UInt_t  fLinZ(1);
  const UInt_t  fLinY(1);
  const UInt_t  fCol[NTrgs]  = {859, 899};
  const UInt_t  fMarZ[NTrgs] = {27, 24};
  const UInt_t  fMarN[NTrgs] = {27, 24};
  const UInt_t  fMarA[NTrgs] = {27, 24};
  const Float_t fLab(0.04);
  const Float_t fOffsetL(0.007);
  const Float_t fOffsetX(1.);
  const Float_t fOffsetY(1.5);
  const Float_t dFrange[2] = {0., 1.03};
  const Float_t rRange[2]  = {0., 1.};
  const Float_t yRange[2]  = {0., 10.};
  const TString sTitleX("#Delta#varphi^{trk} = #varphi^{trk} - #varphi^{trg}");
  const TString sTitleY("D_{pp}");
  const TString sTitleXR("E_{T}^{trg} [GeV]");
  const TString sTitleYR("B = D^{NS}_{pp}(#gamma^{rich}) / D^{NS}_{pp}(#pi^{0})");
  const TString sTitleYZ("(1/N^{trg}) dN^{trk}/dz_{T}^{trk}");
  const TString sTitleXN("E_{T}^{trg} [GeV]");
  const TString sTitleXA("E_{T}^{trg} [GeV]");
  const TString sTitleYN("D^{NS}_{pp}");
  const TString sTitleYA("D^{AS}_{pp}");
  const TString sTitleN("Near Side, |#Delta#varphi| < 1.3");
  const TString sTitleA("Away Side, |#Delta#varphi - #pi| < 1.3");
  const TString sTitleXZ[NTrgs] = {"z_{T}^{trk}(#pi^{0})", "z_{T}^{trk}(#gamma^{rich})"};
  const TString sTitle[NTrgs]   = {"#pi^{0} trigger", "#gamma_{rich} trigger"};

  // create zT x-axis title
  TString sXaxis("#color[");
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    sXaxis += fCol[iTrg];
    sXaxis += "]{";
    sXaxis += sTitleXZ[iTrg].Data();
    if (iTrg + 1 == NTrgs)
      sXaxis += "}";
    else
      sXaxis += "}, #color[";
  }

  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      // delta phi histograms
      hDeltaPhi[iBin][iTrg] -> SetFillStyle(fFil);
      hDeltaPhi[iBin][iTrg] -> SetMarkerStyle(fMar);
      hDeltaPhi[iBin][iTrg] -> SetTitle(sTitle[iTrg].Data());
      hDeltaPhi[iBin][iTrg] -> SetTitleFont(fTxt);
      hDeltaPhi[iBin][iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hDeltaPhi[iBin][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hDeltaPhi[iBin][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hDeltaPhi[iBin][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hDeltaPhi[iBin][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hDeltaPhi[iBin][iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetLabelOffset(fOffsetL);
      hDeltaPhi[iBin][iTrg] -> GetYaxis() -> SetRangeUser(dFrange[0], dFrange[1]);

      // near-side delta phi histograms
      hDeltaPhiNS[iBin][iTrg] -> SetFillStyle(fFilNS);
      hDeltaPhiNS[iBin][iTrg] -> SetFillColor(fCol[iTrg]);
      hDeltaPhiNS[iBin][iTrg] -> SetLineStyle(fLinNS);
      hDeltaPhiNS[iBin][iTrg] -> SetLineColor(fCol[iTrg]);
      hDeltaPhiNS[iBin][iTrg] -> SetMarkerStyle(fMarNS);
      hDeltaPhiNS[iBin][iTrg] -> SetMarkerColor(fCol[iTrg]);
      hDeltaPhiNS[iBin][iTrg] -> SetTitle(sTitle[iTrg].Data());
      hDeltaPhiNS[iBin][iTrg] -> SetTitleFont(fTxt);
      hDeltaPhiNS[iBin][iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hDeltaPhiNS[iBin][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hDeltaPhiNS[iBin][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hDeltaPhiNS[iBin][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hDeltaPhiNS[iBin][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hDeltaPhiNS[iBin][iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetLabelOffset(fOffsetL);
      hDeltaPhiNS[iBin][iTrg] -> GetYaxis() -> SetRangeUser(dFrange[0], dFrange[1]);

      // away-side delta phi histograms
      hDeltaPhiAS[iBin][iTrg] -> SetFillStyle(fFilAS);
      hDeltaPhiAS[iBin][iTrg] -> SetFillColor(fCol[iTrg]);
      hDeltaPhiAS[iBin][iTrg] -> SetLineStyle(fLinAS);
      hDeltaPhiAS[iBin][iTrg] -> SetLineColor(fCol[iTrg]);
      hDeltaPhiAS[iBin][iTrg] -> SetMarkerStyle(fMarAS);
      hDeltaPhiAS[iBin][iTrg] -> SetMarkerColor(fCol[iTrg]);
      hDeltaPhiAS[iBin][iTrg] -> SetTitle(sTitle[iTrg].Data());
      hDeltaPhiAS[iBin][iTrg] -> SetTitleFont(fTxt);
      hDeltaPhiAS[iBin][iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
      hDeltaPhiAS[iBin][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hDeltaPhiAS[iBin][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hDeltaPhiAS[iBin][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hDeltaPhiAS[iBin][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hDeltaPhiAS[iBin][iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetLabelOffset(fOffsetL);
      hDeltaPhiAS[iBin][iTrg] -> GetYaxis() -> SetRangeUser(dFrange[0], dFrange[1]);

      // zT histograms
      hZtTrk[iBin][iTrg] -> SetFillStyle(fFil);
      hZtTrk[iBin][iTrg] -> SetFillColor(fCol[iTrg]);
      hZtTrk[iBin][iTrg] -> SetLineStyle(fLinZ);
      hZtTrk[iBin][iTrg] -> SetLineColor(fCol[iTrg]);
      hZtTrk[iBin][iTrg] -> SetMarkerStyle(fMarZ[iTrg]);
      hZtTrk[iBin][iTrg] -> SetMarkerColor(fCol[iTrg]);
      hZtTrk[iBin][iTrg] -> GetXaxis() -> SetTitle(sXaxis.Data());
      hZtTrk[iBin][iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
      hZtTrk[iBin][iTrg] -> GetXaxis() -> CenterTitle(fCnt);
      hZtTrk[iBin][iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
      hZtTrk[iBin][iTrg] -> GetXaxis() -> SetLabelSize(fLab);
      hZtTrk[iBin][iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
      hZtTrk[iBin][iTrg] -> GetYaxis() -> SetTitle(sTitleYZ.Data());
      hZtTrk[iBin][iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
      hZtTrk[iBin][iTrg] -> GetYaxis() -> SetTitleOffset(fOffsetY);
      hZtTrk[iBin][iTrg] -> GetYaxis() -> CenterTitle(fCnt);
      hZtTrk[iBin][iTrg] -> GetYaxis() -> SetLabelSize(fLab);
      hZtTrk[iBin][iTrg] -> GetYaxis() -> SetLabelFont(fTxt);

    }  // end trigger id loop
  }  // end trigger eT loop

  // yield histograms
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    hYieldNS[iTrg] -> SetTitle(sTitleN.Data());
    hYieldNS[iTrg] -> SetFillStyle(fFil);
    hYieldNS[iTrg] -> SetFillColor(fCol[iTrg]);
    hYieldNS[iTrg] -> SetLineStyle(fLinY);
    hYieldNS[iTrg] -> SetLineColor(fCol[iTrg]);
    hYieldNS[iTrg] -> SetMarkerStyle(fMarN[iTrg]);
    hYieldNS[iTrg] -> SetMarkerColor(fCol[iTrg]);
    hYieldNS[iTrg] -> GetXaxis() -> SetTitle(sTitleXN.Data());
    hYieldNS[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hYieldNS[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hYieldNS[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hYieldNS[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hYieldNS[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hYieldNS[iTrg] -> GetYaxis() -> SetTitle(sTitleYN.Data());
    hYieldNS[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hYieldNS[iTrg] -> GetYaxis() -> SetTitleOffset(fOffsetY);
    hYieldNS[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hYieldNS[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hYieldNS[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hYieldNS[iTrg] -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);

    hYieldAS[iTrg] -> SetTitle(sTitleA.Data());
    hYieldAS[iTrg] -> SetFillStyle(fFil);
    hYieldAS[iTrg] -> SetFillColor(fCol[iTrg]);
    hYieldAS[iTrg] -> SetLineStyle(fLinY);
    hYieldAS[iTrg] -> SetLineColor(fCol[iTrg]);
    hYieldAS[iTrg] -> SetMarkerStyle(fMarA[iTrg]);
    hYieldAS[iTrg] -> SetMarkerColor(fCol[iTrg]);
    hYieldAS[iTrg] -> GetXaxis() -> SetTitle(sTitleXA.Data());
    hYieldAS[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hYieldAS[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hYieldAS[iTrg] -> GetXaxis() -> SetTitleOffset(fOffsetX);
    hYieldAS[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hYieldAS[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hYieldAS[iTrg] -> GetYaxis() -> SetTitle(sTitleYA.Data());
    hYieldAS[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hYieldAS[iTrg] -> GetYaxis() -> SetTitleOffset(fOffsetY);
    hYieldAS[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hYieldAS[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hYieldAS[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hYieldAS[iTrg] -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
  }  // end trigger id loop

  // R histogram
  hPurity -> SetFillStyle(fFil);
  hPurity -> SetMarkerStyle(fMar);
  hPurity -> GetXaxis() -> SetTitle(sTitleXR.Data());
  hPurity -> GetXaxis() -> SetTitleFont(fTxt);
  hPurity -> GetXaxis() -> CenterTitle(fCnt);
  hPurity -> GetXaxis() -> SetTitleOffset(fOffsetX);
  hPurity -> GetXaxis() -> SetLabelSize(fLab);
  hPurity -> GetXaxis() -> SetLabelFont(fTxt);
  hPurity -> GetYaxis() -> SetTitle(sTitleYR.Data());
  hPurity -> GetYaxis() -> SetTitleFont(fTxt);
  hPurity -> GetYaxis() -> SetTitleOffset(fOffsetY);
  hPurity -> GetYaxis() -> CenterTitle(fCnt);
  hPurity -> GetYaxis() -> SetLabelSize(fLab);
  hPurity -> GetYaxis() -> SetLabelFont(fTxt);
  hPurity -> GetYaxis() -> SetRangeUser(rRange[0], rRange[1]);
  cout << "    Set styles." << endl;


  // make label
  const UInt_t  nDec(3);
  const UInt_t  fColP(0);
  const UInt_t  fColT(1);
  const UInt_t  fStyP(0);
  const UInt_t  fAlign(12);
  const Float_t xyLabel[4]  = {0.1, 0.1, 0.3, 0.3};
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg[NTrgs] = {"#pi^{0}", "#gamma^{rich}"};

  TPaveText *pInfoZt[NBins];
  TPaveText *pInfoDf[NBins][NTrgs];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {

    // zT labels
    pInfoZt[iBin] = new TPaveText(xyLabel[0], xyLabel[1], xyLabel[2], xyLabel[3], "NDC NB");
    pInfoZt[iBin] -> SetFillColor(fColP);
    pInfoZt[iBin] -> SetFillStyle(fStyP);
    pInfoZt[iBin] -> SetLineColor(fColP);
    pInfoZt[iBin] -> SetLineStyle(fStyP);
    pInfoZt[iBin] -> SetTextFont(fTxt);
    pInfoZt[iBin] -> SetTextAlign(fAlign);

    // create system and kinematic info
    TString sEtTrg("E_{T}^{trg} #in (");
    TString sZtTrk("z_{T}^{trk} #in (");
    TString sKinetic("");
    sEtTrg   += eTtrgMin[iBin];
    sEtTrg   += ", ";
    sEtTrg   += eTtrgMax[iBin];
    sEtTrg   += ")";
    sZtTrk   += sZtTrkMin.Data();
    sZtTrk   += ", ";
    sZtTrk   += sZtTrkMax.Data();
    sZtTrk   += ")";
    sKinetic += sEtTrg.Data();
    sKinetic += " #otimes ";
    sKinetic += sZtTrk.Data();

    // add text
    pInfoZt[iBin] -> AddText(sSystem.Data());
    pInfoZt[iBin] -> AddText(sKinetic.Data());

    // delta-phi plots
    for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
      pInfoDf[iBin][iTrg] = new TPaveText(xyLabel[0], xyLabel[1], xyLabel[2], xyLabel[3], "NDC NB");
      pInfoDf[iBin][iTrg] -> SetFillColor(fColP);
      pInfoDf[iBin][iTrg] -> SetFillStyle(fStyP);
      pInfoDf[iBin][iTrg] -> SetLineColor(fColP);
      pInfoDf[iBin][iTrg] -> SetLineStyle(fStyP);
      pInfoDf[iBin][iTrg] -> SetTextFont(fTxt);
      pInfoDf[iBin][iTrg] -> SetTextColor(fCol[iTrg]);
      pInfoDf[iBin][iTrg] -> SetTextAlign(fAlign);

      // convert into strings
      TString sPedRaw("");
      TString sPedTxt("");
      TString sYieldRaw("");
      TString sYieldTxt("");
      sPedRaw   += pedestal[iBin][iTrg];
      sYieldRaw += yieldNS[iBin][iTrg];

      // display only 3 decimal places
      const UInt_t nPed      = sPedRaw.First(".");
      const UInt_t nYield    = sYieldRaw.First(".");
      const UInt_t nPedTxt   = (nPed + nDec) + 1;
      const UInt_t nYieldTxt = (nYield + nDec) + 1;
      sPedTxt.Append(sPedRaw, nPedTxt);
      sYieldTxt.Append(sYieldRaw, nYieldTxt);

      TString sPed("ped. = ");
      TString sYield("yield = ");
      TString sNumbers("");
      sPed     += sPedTxt;
      sYield   += sYieldTxt;
      sNumbers += sPed.Data();
      sNumbers += ", ";
      sNumbers += sYield.Data();

      // add text
      pInfoDf[iBin][iTrg] -> AddText(sSystem.Data());
      pInfoDf[iBin][iTrg] -> AddText(sKinetic.Data());
      pInfoDf[iBin][iTrg] -> AddText(sNumbers.Data());
    }  // end trigger id loop
  }  // end trigger eT loop

  // make labels for yield plots
  TPaveText *ptNear = new TPaveText(xyLabel[0], xyLabel[1], xyLabel[2], xyLabel[3], "NDC NB");
  TPaveText *ptAway = new TPaveText(xyLabel[0], xyLabel[1], xyLabel[2], xyLabel[3], "NDC NB");
  ptNear -> SetFillColor(fColP);
  ptNear -> SetFillStyle(fStyP);
  ptNear -> SetLineColor(fColP);
  ptNear -> SetLineStyle(fStyP);
  ptNear -> SetTextFont(fTxt);
  ptNear -> SetTextColor(fColT);
  ptNear -> SetTextAlign(fAlign);
  ptAway -> SetFillColor(fColP);
  ptAway -> SetFillStyle(fStyP);
  ptAway -> SetLineColor(fColP);
  ptAway -> SetLineStyle(fStyP);
  ptAway -> SetTextFont(fTxt);
  ptAway -> SetTextColor(fColT);
  ptAway -> SetTextAlign(fAlign);

  TString sEtTrg("E_{T}^{trg} #in (");
  TString sZtTrk("z_{T}^{trk} #in (");
  TString sKinetic("");
  sEtTrg   += eTtrgBin[0];
  sEtTrg   += ", ";
  sEtTrg   += eTtrgBin[NHist];
  sEtTrg   += ")";
  sZtTrk   += sZtTrkMin.Data();
  sZtTrk   += ", ";
  sZtTrk   += sZtTrkMax.Data();
  sZtTrk   += ")";
  sKinetic += sEtTrg.Data();
  sKinetic += " #otimes ";
  sKinetic += sZtTrk.Data();

  ptNear -> AddText(sSystem.Data());
  ptNear -> AddText(sKinetic.Data());
  ptAway -> AddText(sSystem.Data());
  ptAway -> AddText(sKinetic.Data());
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    // convert into strings
    TString sNearValRaw("");
    TString sNearValTxt("");
    TString sNearErrRaw("");
    TString sNearErrTxt("");
    TString sAwayValRaw("");
    TString sAwayValTxt("");
    TString sAwayErrRaw("");
    TString sAwayErrTxt("");
    sNearValRaw += yieldNS[0][iTrg];
    sNearErrRaw += errorNS[0][iTrg];
    sAwayValRaw += yieldAS[0][iTrg];
    sAwayErrRaw += errorAS[0][iTrg];

    // display only 3 decimal places
    const UInt_t nNearVal    = sNearValRaw.First(".");
    const UInt_t nNearErr    = sNearErrRaw.First(".");
    const UInt_t nAwayVal    = sAwayValRaw.First(".");
    const UInt_t nAwayErr    = sAwayErrRaw.First(".");
    const UInt_t nNearValTxt = (nNearVal + nDec) + 1;
    const UInt_t nNearErrTxt = (nNearErr + nDec) + 1;
    const UInt_t nAwayValTxt = (nAwayVal + nDec) + 1;
    const UInt_t nAwayErrTxt = (nAwayErr + nDec) + 1;
    sNearValTxt.Append(sNearValRaw, nNearValTxt);
    sNearErrTxt.Append(sNearErrRaw, nNearErrTxt);
    sAwayValTxt.Append(sAwayValRaw, nAwayValTxt);
    sAwayErrTxt.Append(sAwayErrRaw, nAwayErrTxt);

    TString sNearVal("#color[");
    TString sAwayVal("#color[");
    sNearVal += fCol[iTrg];
    sNearVal += "]{D^{NS}_{pp}(";
    sNearVal += sTrg[iTrg].Data();
    sNearVal += " | ";
    sNearVal += eTtrgBin[0];
    sNearVal += " - ";
    sNearVal += eTtrgBin[NHist];
    sNearVal += " GeV) = ";
    sNearVal += sNearValTxt.Data();
    sNearVal += " #pm ";
    sNearVal += sNearErrTxt.Data();
    sNearVal += "}";
    sAwayVal += fCol[iTrg];
    sAwayVal += "]{D^{AS}_{pp}(";
    sAwayVal += sTrg[iTrg].Data();
    sAwayVal += " | ";
    sAwayVal += eTtrgBin[0];
    sAwayVal += " - ";
    sAwayVal += eTtrgBin[NHist];
    sAwayVal += " GeV) = ";
    sAwayVal += sAwayValTxt.Data();
    sAwayVal += " #pm ";
    sAwayVal += sAwayErrTxt.Data();
    sAwayVal += "}";
    ptNear   -> AddText(sNearVal.Data());
    ptAway   -> AddText(sAwayVal.Data());
  }  // end trigger id loop

  // make label for R plot
  TString sEtAvg("E_{T}^{trg} #in (");
  sEtAvg += eTtrgMin[0];
  sEtAvg += ", ";
  sEtAvg += eTtrgMax[0];
  sEtAvg += ") #otimes z_{T}^{trk} #in (";
  sEtAvg += sZtTrkMin;
  sEtAvg += ", ";
  sEtAvg += sZtTrkMax;
  sEtAvg += ")";

  TString sRraw("");
  TString sEraw("");
  TString sRtxt("");
  TString sEtxt("");
  sRraw += valueR[0];
  sEraw += errorR[0];
  sRtxt += "#LTB#GT = ";
  sEtxt += " #pm ";

  const UInt_t nRraw = sRraw.First(".");
  const UInt_t nEraw = sEraw.First(".");
  const UInt_t nRtxt = (nRraw + nDec) + 1;
  const UInt_t nEtxt = (nEraw + nDec) + 1;
  sRtxt.Append(sRraw, nRtxt);
  sEtxt.Append(sEraw, nEtxt);
  sRtxt.Append(sEtxt);

  TPaveText *pPurity = new TPaveText(xyLabel[0], xyLabel[1], xyLabel[2], xyLabel[3], "NDC NB");
  pPurity -> SetFillColor(fColP);
  pPurity -> SetFillStyle(fStyP);
  pPurity -> SetLineColor(fColP);
  pPurity -> SetLineStyle(fStyP);
  pPurity -> SetTextFont(fTxt);
  pPurity -> SetTextColor(fColR);
  pPurity -> SetTextAlign(fAlign);
  pPurity -> AddText(sSystem.Data());
  pPurity -> AddText("Line indicates average over E_{T}^{trg}");
  pPurity -> AddText(sEtAvg.Data());
  pPurity -> AddText(sRtxt.Data());
  cout << "    Made labels." << endl;


  // make legends
  const Float_t xyLeg[4] = {0.3, 0.1, 0.5, 0.3};

  TLegend *lZtTrk = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  lZtTrk -> SetFillColor(fColP);
  lZtTrk -> SetFillStyle(fStyP);
  lZtTrk -> SetLineColor(fColP);
  lZtTrk -> SetLineStyle(fStyP);
  lZtTrk -> SetTextFont(fTxt);
  lZtTrk -> SetTextAlign(fAlign);
  lZtTrk -> AddEntry(hZtTrk[0][0], "#pi^{0}+h^{#pm}", "fp");
  lZtTrk -> AddEntry(hZtTrk[0][1], "#gamma^{rich}+h^{#pm}", "fp");

  TLegend *lYield = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  lYield -> SetFillColor(fColP);
  lYield -> SetFillStyle(fStyP);
  lYield -> SetLineColor(fColP);
  lYield -> SetLineStyle(fStyP);
  lYield -> SetTextFont(fTxt);
  lYield -> SetTextAlign(fAlign);
  lYield -> AddEntry(hYieldAS[0], "#pi^{0}+h^{#pm}", "fp");
  lYield -> AddEntry(hYieldAS[1], "#gamma^{rich}+h^{#pm}", "fp");
  cout << "    Made legend." << endl;


  // create directories and save histograms
  const TString sDir("eTtrg");

  TDirectory *dBin[NBins];
  for (UInt_t iBin = 0; iBin < NBins; iBin++) {
    // make name
    TString sDirName(sDir.Data());
    sDirName += sEtNames[iBin].Data();

    // make directory
    dBin[iBin] = (TDirectory*) fOutput -> mkdir(sDirName.Data());
    dBin[iBin] -> cd();
    for (iTrg = 0; iTrg < NTrgs; iTrg++) {
      hDeltaPhi[iBin][iTrg]   -> Write();
      hDeltaPhiNS[iBin][iTrg] -> Write();
      hDeltaPhiAS[iBin][iTrg] -> Write();
      hEtaTrk[iBin][iTrg]     -> Write();
      hPtTrk[iBin][iTrg]      -> Write();
      hZtTrk[iBin][iTrg]      -> Write();
    }  // end trigger id loop
  }  // end trigger eT loop
  cout << "    Made directories." << endl;


  // draw plots
  const UInt_t  widthR(750);
  const UInt_t  heightR(750);
  const UInt_t  widthDf(1500);
  const UInt_t  heightDf(750);
  const UInt_t  widthZt(750);
  const UInt_t  heightZt(750);
  const UInt_t  widthY(1500);
  const UInt_t  heightY(750);
  const UInt_t  grid(0);
  const UInt_t  log(1);
  const UInt_t  tick(1);
  const Float_t marginBig(0.15);
  const Float_t marginMed(0.1);
  const Float_t marginBord(0.005);
  const Float_t marginSmall(0.02);
  const Float_t xPadN(0.5);
  const Float_t xPadA(1.);
  const Float_t xPadDf[NTrgs]  = {0.5, 1.};
  const TString sRcanvas("cPurity");
  const TString sDfCanvas("cDeltaPhi_");
  const TString sZtCanvas("cZt_");
  const TString sYcanvas("cYield");
  const TString sYpadN("pNearSide");
  const TString sYpadA("pAwaySide");
  const TString sDfPads[NTrgs] = {"pPi0", "pGamma"};

  TPad    *pDeltaPhi[NBins][NTrgs];
  TCanvas *cDeltaPhi[NBins];
  TCanvas *cZtTrk[NBins];
  fOutput -> cd();

  for (UInt_t iBin = 0; iBin < NBins; iBin++) {

    // make canvas names
    TString sDfCanvasName(sDfCanvas.Data());
    TString sZtCanvasName(sZtCanvas.Data());
    sDfCanvasName += "et";
    sDfCanvasName += sEtNames[iBin].Data();
    sZtCanvasName += "et";
    sZtCanvasName += sEtNames[iBin].Data();

    // delta-phi plots
    dBin[iBin] -> cd();
    cDeltaPhi[iBin]    = new TCanvas(sDfCanvasName.Data(), "", widthDf, heightDf);
    pDeltaPhi[iBin][0] = new TPad(sDfPads[0].Data(), "", 0., 0., xPadDf[0], 1.);
    pDeltaPhi[iBin][1] = new TPad(sDfPads[1].Data(), "", xPadDf[0], 0., xPadDf[1], 1.);
    pDeltaPhi[iBin][0]   -> SetGrid(grid, grid);
    pDeltaPhi[iBin][0]   -> SetTicks(tick, tick);
    pDeltaPhi[iBin][0]   -> SetRightMargin(marginBord);
    pDeltaPhi[iBin][0]   -> SetTopMargin(marginMed);
    pDeltaPhi[iBin][0]   -> SetLeftMargin(marginBig);
    pDeltaPhi[iBin][0]   -> SetBottomMargin(marginBig);
    pDeltaPhi[iBin][1]   -> SetGrid(grid, grid);
    pDeltaPhi[iBin][1]   -> SetTicks(tick, tick);
    pDeltaPhi[iBin][1]   -> SetRightMargin(marginBig);
    pDeltaPhi[iBin][1]   -> SetTopMargin(marginMed);
    pDeltaPhi[iBin][1]   -> SetLeftMargin(marginBord);
    pDeltaPhi[iBin][1]   -> SetBottomMargin(marginBig);
    cDeltaPhi[iBin]      -> cd();
    pDeltaPhi[iBin][0]   -> Draw();
    pDeltaPhi[iBin][1]   -> Draw();
    pDeltaPhi[iBin][0]   -> cd();
    hDeltaPhi[iBin][0]   -> Draw();
    hDeltaPhiNS[iBin][0] -> Draw("SAME LF HIST");
    hDeltaPhiAS[iBin][0] -> Draw("SAME LF HIST");
    pInfoDf[iBin][0]     -> Draw();
    pDeltaPhi[iBin][1]   -> cd();
    hDeltaPhi[iBin][1]   -> Draw();
    hDeltaPhiNS[iBin][1] -> Draw("SAME LF HIST");
    hDeltaPhiAS[iBin][1] -> Draw("SAME LF HIST");
    pInfoDf[iBin][1]     -> Draw();
    cDeltaPhi[iBin]      -> Write();
    cDeltaPhi[iBin]      -> Close();

    // zT plots
    dBin[iBin] -> cd();
    cZtTrk[iBin] = new TCanvas(sZtCanvasName.Data(), "", widthZt, heightZt);
    cZtTrk[iBin]    -> SetGrid(grid, grid);
    cZtTrk[iBin]    -> SetTicks(tick, tick);
    cZtTrk[iBin]    -> SetTopMargin(marginSmall);
    cZtTrk[iBin]    -> SetRightMargin(marginSmall);
    cZtTrk[iBin]    -> SetBottomMargin(marginBig);
    cZtTrk[iBin]    -> SetLeftMargin(marginBig);
    cZtTrk[iBin]    -> SetLogy(log);
    cZtTrk[iBin]    -> cd();
    hZtTrk[iBin][0] -> Draw();
    hZtTrk[iBin][1] -> Draw("same");
    pInfoZt[iBin]   -> Draw();
    lZtTrk          -> Draw();
    cZtTrk[iBin]    -> Write();
    cZtTrk[iBin]    -> Close();

  }  // end trigger eT loop

  // yield plots
  fOutput -> cd();
  TCanvas *cYield = new TCanvas(sYcanvas.Data(), "", widthY, heightY);
  TPad    *pNear  = new TPad(sYpadN.Data(), "", 0., 0., xPadN, 1.);
  TPad    *pAway  = new TPad(sYpadA.Data(), "", xPadN, 0., xPadA, 1.);
  pNear         -> SetGrid(grid, grid);
  pNear         -> SetTicks(tick, tick);
  pNear         -> SetTopMargin(marginMed);
  pNear         -> SetRightMargin(marginSmall);
  pNear         -> SetBottomMargin(marginBig);
  pNear         -> SetLeftMargin(marginBig);
  pAway         -> SetGrid(grid, grid);
  pAway         -> SetTicks(tick, tick);
  pAway         -> SetTopMargin(marginMed);
  pAway         -> SetRightMargin(marginSmall);
  pAway         -> SetBottomMargin(marginBig);
  pAway         -> SetLeftMargin(marginBig);
  cYield        -> cd();
  pNear         -> Draw();
  pAway         -> Draw();
  pNear         -> cd();
  hYieldNS[0]   -> Draw();
  hYieldNS[1]   -> Draw("same");
  lYieldNS[0]   -> Draw();
  lYieldNS[1]   -> Draw();
  lYieldNsHi[0] -> Draw();
  lYieldNsHi[1] -> Draw();
  lYieldNsLo[0] -> Draw();
  lYieldNsLo[1] -> Draw();
  ptNear        -> Draw();
  pAway         -> cd();
  hYieldAS[0]   -> Draw();
  hYieldAS[1]   -> Draw("same");
  lYieldAS[0]   -> Draw();
  lYieldAS[1]   -> Draw();
  lYieldAsHi[0] -> Draw();
  lYieldAsHi[1] -> Draw();
  lYieldAsLo[0] -> Draw();
  lYieldAsLo[1] -> Draw();
  ptAway        -> Draw();
  lYield        -> Draw();
  cYield        -> Write();
  cYield        -> Close();

  // purity plots
  fOutput -> cd();
  TCanvas *cPurity = new TCanvas(sRcanvas.Data(), "", widthR, heightR);
  cPurity  -> SetGrid(grid, grid);
  cPurity  -> SetTicks(tick, tick);
  cPurity  -> SetTopMargin(marginSmall);
  cPurity  -> SetRightMargin(marginSmall);
  cPurity  -> SetBottomMargin(marginBig);
  cPurity  -> SetLeftMargin(marginBig);
  cPurity  -> cd();
  hPurity  -> Draw();
  lAverage -> Draw();
  lAvgUp   -> Draw();
  lAvgDown -> Draw();
  pPurity  -> Draw();
  cPurity  -> Write();
  cPurity  -> Close();
  cout << "    Drew plots." << endl;


  // save yield histograms
  fOutput -> cd();
  for (UInt_t iTrg = 0; iTrg < NTrgs; iTrg++) {
    hYieldNS[iTrg] -> Write();
    hYieldAS[iTrg] -> Write();
  }

  // close files
  fOutput -> cd();
  hPurity -> Write();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Calculation finished!\n" << endl;

}

// End ------------------------------------------------------------------------
