// 'StDataQaMaker.eventQaPlots.h'
// Derek Anderson
// 09.19.2019
//
// This class produces the data QA
// plots for the neutral-triggered
// pp recoil jet analysis note.

#pragma once

using namespace std;



void StDataQaMaker::MakeEventQaPlots() {

  cout << "    Making event QA plots..." << endl;

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
  const Double_t pTbinMin[NTrkBins] = {0.2, 2., 5., 9.};
  const Double_t pTbinMax[NTrkBins] = {2., 5., 9., 20.};

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

  // create event histograms
  TH1D *hEvtNum;
  TH1D *hEvtPrim;
  TH1D *hEvtVz[2];
  TH1D *hEvtVr[2];
  TH2D *hEvtVyVsVx;
  TH2D *hEvtVzVsVx;
  TH2D *hEvtVzVsVy;
  // create trigger histograms
  TH1D *hTrgEt[2];
  TH1D *hTrgEta;
  TH1D *hTrgPhi;
  TH1D *hTrgTsp[NTrgTsp + 1];
  TH1D *hTrgEtBin[NTrgBins][NTrgTsp];
  TH2D *hTrgEtVsEta;
  TH2D *hTrgEtVsPhi;
  TH2D *hTrgEtaVsPhi;
  // create track histograms
  TH1D *hTrkNfit[2];
  TH1D *hTrkRfit[2];
  TH1D *hTrkDca[2];
  TH1D *hTrkEta[2];
  TH1D *hTrkPt[2];
  TH1D *hTrkPtBin[NTrgTsp];
  TH1D *hTrkDfBin[NTrgTsp];
  TH1D *hTrkEtaBin[NTrgTsp];
  TH1D *hTrkDfPtBin[NTrgTsp][NTrkBins];
  TH1D *hTrkEtaPtBin[NTrgTsp][NTrkBins];
  TH2D *hTrkNfitVsPt;
  TH2D *hTrkDcaVsPt;
  TH2D *hTrkPtVsDf[NTrgTsp];
  TH2D *hTrkPtVsEta[NTrgTsp];
  TH2D *hTrkDfVsEta[NTrgTsp];

  const UInt_t  nNum(NTrgCuts);
  const UInt_t  nVz(400);
  const UInt_t  nVr(400);
  const UInt_t  nPrim(100);
  const UInt_t  nEt(50);
  const UInt_t  nEta(40);
  const UInt_t  nPhi(30);
  const UInt_t  nTsp(100);
  const UInt_t  nFit(50);
  const UInt_t  nRat(100);
  const UInt_t  nDca(50);
  const UInt_t  nPt(275);
  const UInt_t  nPt2(110);
  const UInt_t  nPt3(55);
  const UInt_t  nDf(30);
  const Float_t num[2]  = {0., (Float_t) NTrgCuts};
  const Float_t vz[2]   = {-200., 200.};
  const Float_t vr[2]   = {-2., 2.};
  const Float_t prim[2] = {0., 100.};
  const Float_t et[2]   = {0., 50.};
  const Float_t eta[2]  = {-2., 2.};
  const Float_t phi[2]  = {-3.15, 3.15};
  const Float_t tsp[2]  = {0., 1.};
  const Float_t fit[2]  = {0., 50.};
  const Float_t rat[2]  = {0., 1.};
  const Float_t dca[2]  = {0., 5.};
  const Float_t pt[2]   = {-5., 50.};
  const Float_t df[2]   = {-1.6, 4.7};

  // define event histograms
  hEvtNum            = new TH1D("hEvtNum", "", nNum, num[0], num[1]);
  hEvtVz[0]          = new TH1D("hEvtVzAll", "", nVz, vz[0], vz[1]);
  hEvtVz[1]          = new TH1D("hEvtVzCut", "", nVz, vz[0], vz[1]);
  hEvtVr[0]          = new TH1D("hEvtVrAll", "", nVr, vr[0], vr[1]);
  hEvtVr[1]          = new TH1D("hEvtVrCut", "", nVr, vr[0], vr[1]);
  hEvtPrim           = new TH1D("hEvtPrim", "", nPrim, prim[0], prim[1]);
  hEvtVyVsVx         = new TH2D("hEvtVyVsVx", "", nVr, vr[0], vr[1], nVr, vr[0], vr[1]);
  hEvtVzVsVx         = new TH2D("hEvtVzVsVx", "", nVr, vr[0], vr[1], nVz, vz[0], vz[1]);
  hEvtVzVsVy         = new TH2D("hEvtVzVsVy", "", nVr, vr[0], vr[1], nVz, vz[0], vz[1]);
  // define trigger histograms
  hTrgEt[0]          = new TH1D("hTrgEtAll", "", nEt, et[0], et[1]);
  hTrgEt[1]          = new TH1D("hTrgEtCut", "", nEt, et[0], et[1]);
  hTrgEta            = new TH1D("hTrgEta", "", nEta, eta[0], eta[1]);
  hTrgPhi            = new TH1D("hTrgPhi", "", nPhi, phi[0], phi[1]);
  hTrgTsp[0]         = new TH1D("hTrgTspPi0", "", nTsp, tsp[0], tsp[1]);
  hTrgTsp[1]         = new TH1D("hTrgTspGam", "", nTsp, tsp[0], tsp[1]);
  hTrgTsp[2]         = new TH1D("hTrgTspAll", "", nTsp, tsp[0], tsp[1]); 
  hTrgEtBin[0][0]    = new TH1D("hTrgEtPi920", "", nEt, et[0], et[1]);
  hTrgEtBin[0][1]    = new TH1D("hTrgEtGa920", "", nEt, et[0], et[1]);
  hTrgEtBin[1][0]    = new TH1D("hTrgEtPi911", "", nEt, et[0], et[1]);
  hTrgEtBin[1][1]    = new TH1D("hTrgEtGa911", "", nEt, et[0], et[1]);
  hTrgEtBin[2][0]    = new TH1D("hTrgEtPi1115", "", nEt, et[0], et[1]);
  hTrgEtBin[2][1]    = new TH1D("hTrgEtGa1115", "", nEt, et[0], et[1]);
  hTrgEtBin[3][0]    = new TH1D("hTrgEtPi1520", "", nEt, et[0], et[1]);
  hTrgEtBin[3][1]    = new TH1D("hTrgEtGa1520", "", nEt, et[0], et[1]);
  hTrgEtVsEta        = new TH2D("hTrgEtVsEta", "", nEta, eta[0], eta[1], nEt, et[0], et[1]);
  hTrgEtVsPhi        = new TH2D("hTrgEtVsPhi", "", nPhi, phi[0], phi[1], nEt, et[0], et[1]);
  hTrgEtaVsPhi       = new TH2D("hTrgEtaVsPhi", "", nPhi, phi[0], phi[1], nEta, eta[0], eta[1]);
  // define track histograms
  hTrkNfit[0]        = new TH1D("hTrkNfitAll", "", nFit, fit[0], fit[1]);
  hTrkNfit[1]        = new TH1D("hTrkNfitCut", "", nFit, fit[0], fit[1]);
  hTrkRfit[0]        = new TH1D("hTrkRfitAll", "", nRat, rat[0], rat[1]);
  hTrkRfit[1]        = new TH1D("hTrkRfitCut", "", nRat, rat[0], rat[1]);
  hTrkDca[0]         = new TH1D("hTrkDcaAll", "", nDca, dca[0], dca[1]);
  hTrkDca[1]         = new TH1D("hTrkDcaCut", "", nDca, dca[0], dca[1]);
  hTrkEta[0]         = new TH1D("hTrkEtaAll", "", nEta, eta[0], eta[1]);
  hTrkEta[1]         = new TH1D("hTrkEtaCut", "", nEta, eta[0], eta[1]);
  hTrkPt[0]          = new TH1D("hTrkPtAll", "", nPt, pt[0], pt[1]);
  hTrkPt[1]          = new TH1D("hTrkPtCut", "", nPt, pt[0], pt[1]);
  hTrkPtBin[0]       = new TH1D("hTrkPtPi0", "", nPt2, pt[0], pt[1]);
  hTrkPtBin[1]       = new TH1D("hTrkPtGam", "", nPt2, pt[0], pt[1]);
  hTrkDfBin[0]       = new TH1D("hTrkDfPi0", "", nDf, df[0], df[1]);
  hTrkDfBin[1]       = new TH1D("hTrkDfGam", "", nDf, df[0], df[1]);
  hTrkEtaBin[0]      = new TH1D("hTrkEtaPi0", "", nEta, eta[0], eta[1]);
  hTrkEtaBin[1]      = new TH1D("hTrkEtaGam", "", nEta, eta[0], eta[1]);
  hTrkDfPtBin[0][0]  = new TH1D("hTrkDfPtBinPi0pt022", "", nDf, df[0], df[1]);
  hTrkDfPtBin[0][1]  = new TH1D("hTrkDfPtBinPi0pt25", "", nDf, df[0], df[1]);
  hTrkDfPtBin[0][2]  = new TH1D("hTrkDfPtBinPi0pt59", "", nDf, df[0], df[1]);
  hTrkDfPtBin[0][3]  = new TH1D("hTrkDfPtBinPi0pt920", "", nDf, df[0], df[1]);
  hTrkDfPtBin[1][0]  = new TH1D("hTrkDfPtBinGamPt022", "", nDf, df[0], df[1]);
  hTrkDfPtBin[1][1]  = new TH1D("hTrkDfPtBinGamPt25", "", nDf, df[0], df[1]);
  hTrkDfPtBin[1][2]  = new TH1D("hTrkDfPtBinGamPt59", "", nDf, df[0], df[1]);
  hTrkDfPtBin[1][3]  = new TH1D("hTrkDfPtBinGamPt920", "", nDf, df[0], df[1]);
  hTrkEtaPtBin[0][0] = new TH1D("hTrkEtaPtBinPi0pt022", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[0][1] = new TH1D("hTrkEtaPtBinPi0pt25", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[0][2] = new TH1D("hTrkEtaPtBinPi0pt59", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[0][3] = new TH1D("hTrkEtaPtBinPi0pt920", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[1][0] = new TH1D("hTrkEtaPtBinGamPt022", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[1][1] = new TH1D("hTrkEtaPtBinGamPt25", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[1][2] = new TH1D("hTrkEtaPtBinGamPt59", "", nEta, eta[0], eta[1]);
  hTrkEtaPtBin[1][3] = new TH1D("hTrkEtaPtBinGamPt920", "", nEta, eta[0], eta[1]);
  hTrkNfitVsPt       = new TH2D("hTrkNfitVsPt", "", nPt3, pt[0], pt[1], nFit, fit[0], fit[1]);
  hTrkDcaVsPt        = new TH2D("hTrkDcaVsPt", "", nPt3, pt[0], pt[1], nDca, dca[0], dca[1]);
  hTrkPtVsDf[0]      = new TH2D("hTrkPtVsDfPi0", "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
  hTrkPtVsDf[1]      = new TH2D("hTrkPtVsDfGam", "", nDf, df[0], df[1], nPt3, pt[0], pt[1]);
  hTrkPtVsEta[0]     = new TH2D("hTrkPtVsEtaPi0", "", nEta, eta[0], eta[1], nPt3, pt[0], pt[1]);
  hTrkPtVsEta[1]     = new TH2D("hTrkPtVsEtaGam", "", nEta, eta[0], eta[1], nPt3, pt[0], pt[1]);
  hTrkDfVsEta[0]     = new TH2D("hTrkDfVsEtaPi0", "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
  hTrkDfVsEta[1]     = new TH2D("hTrkDfVsEtaGam", "", nEta, eta[0], eta[1], nDf, df[0], df[1]);
  // errors
  hEvtVz[0]          -> Sumw2();
  hEvtVz[1]          -> Sumw2();
  hEvtVr[0]          -> Sumw2();
  hEvtVr[1]          -> Sumw2();
  hEvtPrim           -> Sumw2();
  hEvtVyVsVx         -> Sumw2();
  hEvtVzVsVx         -> Sumw2();
  hEvtVzVsVy         -> Sumw2();
  hTrgEtVsEta        -> Sumw2();
  hTrgEtVsPhi        -> Sumw2();
  hTrgEtaVsPhi       -> Sumw2();
  hTrgEt[0]          -> Sumw2();
  hTrgEt[1]          -> Sumw2();
  hTrgEta            -> Sumw2();
  hTrgPhi            -> Sumw2();
  hTrgTsp[0]         -> Sumw2();
  hTrgTsp[1]         -> Sumw2();
  hTrgTsp[2]         -> Sumw2();
  hTrgEtBin[0][0]    -> Sumw2();
  hTrgEtBin[0][1]    -> Sumw2();
  hTrgEtBin[1][0]    -> Sumw2();
  hTrgEtBin[1][1]    -> Sumw2();
  hTrgEtBin[2][0]    -> Sumw2();
  hTrgEtBin[2][1]    -> Sumw2();
  hTrgEtBin[3][0]    -> Sumw2();
  hTrgEtBin[3][1]    -> Sumw2();
  hTrkNfit[0]        -> Sumw2();
  hTrkNfit[1]        -> Sumw2();
  hTrkRfit[0]        -> Sumw2();
  hTrkRfit[1]        -> Sumw2();
  hTrkDca[0]         -> Sumw2();
  hTrkDca[1]         -> Sumw2();
  hTrkEta[0]         -> Sumw2();
  hTrkEta[1]         -> Sumw2();
  hTrkPt[0]          -> Sumw2();
  hTrkPt[1]          -> Sumw2();
  hTrkPtBin[0]       -> Sumw2();
  hTrkPtBin[1]       -> Sumw2();
  hTrkDfBin[0]       -> Sumw2();
  hTrkDfBin[1]       -> Sumw2();
  hTrkEtaBin[0]      -> Sumw2();
  hTrkEtaBin[1]      -> Sumw2();
  hTrkDfPtBin[0][0]  -> Sumw2();
  hTrkDfPtBin[0][1]  -> Sumw2();
  hTrkDfPtBin[0][2]  -> Sumw2();
  hTrkDfPtBin[0][3]  -> Sumw2();
  hTrkDfPtBin[1][0]  -> Sumw2();
  hTrkDfPtBin[1][1]  -> Sumw2();
  hTrkDfPtBin[1][2]  -> Sumw2();
  hTrkDfPtBin[1][3]  -> Sumw2();
  hTrkEtaPtBin[0][0] -> Sumw2();
  hTrkEtaPtBin[0][1] -> Sumw2();
  hTrkEtaPtBin[0][2] -> Sumw2();
  hTrkEtaPtBin[0][3] -> Sumw2();
  hTrkEtaPtBin[1][0] -> Sumw2();
  hTrkEtaPtBin[1][1] -> Sumw2();
  hTrkEtaPtBin[1][2] -> Sumw2();
  hTrkEtaPtBin[1][3] -> Sumw2();
  hTrkNfitVsPt       -> Sumw2();
  hTrkDcaVsPt        -> Sumw2();
  hTrkPtVsDf[0]      -> Sumw2();
  hTrkPtVsDf[1]      -> Sumw2();
  hTrkPtVsEta[0]     -> Sumw2();
  hTrkPtVsEta[1]     -> Sumw2();
  hTrkDfVsEta[0]     -> Sumw2();
  hTrkDfVsEta[1]     -> Sumw2();

  // no. of evts and triggers
  UInt_t nTrgCut[NTrgCuts];
  UInt_t nTrgPi0[NTrgBins];
  UInt_t nTrgGam[NTrgBins];
  for (UInt_t iCut = 0; iCut < NTrgCuts; iCut++) {
    nTrgCut[iCut] = 0;
  }
  for (UInt_t iBin = 0; iBin < NTrgBins; iBin++) {
    nTrgPi0[iBin] = 0;
    nTrgGam[iBin] = 0;
  }

  const UInt_t nEvts = tInput -> GetEntriesFast();
  cout << "      Beginning event loop: " << nEvts << " events to process." << endl;

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
    nTrgCut[0]++;
    if (isGoodRun) {
      nTrgCut[1]++;
    } else {
      continue;
    }

    // fill event histograms
    hEvtVz[0]  -> Fill(zVtx);
    hEvtVr[0]  -> Fill(rVtx);
    hEvtPrim   -> Fill(nTrks);
    hEvtVyVsVx -> Fill(xVtx, yVtx);
    hEvtVzVsVx -> Fill(xVtx, zVtx);
    hEvtVzVsVy -> Fill(yVtx, zVtx);

    // vertex cuts
    const Bool_t isInRcut = (TMath::Abs(rVtx) < rVtxMax);
    const Bool_t isInZcut = (TMath::Abs(zVtx) < zVtxMax);
    if (isInRcut && isInZcut) {
      nTrgCut[2]++;
    } else {
      hEvtVz[1] -> Fill(zVtx);
      hEvtVr[1] -> Fill(rVtx);
      continue;
    }

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
    if (isGoodTwr) {
      nTrgCut[3]++;
    } else {
      continue;
    }

    // trigger cuts
    const Bool_t isInAdcCut    = (adc <= adcMax);
    const Bool_t isInStrCut    = ((eH4 >= eStrMin) && (eF4 >= eStrMin));
    const Bool_t isInProjCut   = (pProj < pProjMax);
    const Bool_t isInEtaTrgCut = ((TMath::Abs(hDet) < hTrgMax) && (TMath::Abs(hPhys) < hTrgMax));
    const Bool_t isInEtCut     = ((eTtrg >= eTtrgMin) && (eTtrg < eTtrgMax));
    const Bool_t isInPi0cut    = ((tspTrg > tspPi0[0]) && (tspTrg < tspPi0[1]));
    const Bool_t isInGamCut    = ((tspTrg > tspGam[0]) && (tspTrg < tspGam[1]));
    const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);

    // central strip cut
    if (isInAdcCut && isInStrCut) {
      nTrgCut[4]++;
    } else {
      continue;
    }

    // charged momentum cut
    if (isInProjCut) {
      nTrgCut[5]++;
    } else {
      continue;
    }

    // eta cut
    if (isInEtaTrgCut) {
      nTrgCut[6]++;
    } else {
      continue;
    }

    // eT and tsp cuts
    if (isInEtCut) {
      nTrgCut[7]++;
    }
    if (isInEtCut && isInTspCut) {
      nTrgCut[8]++;
    }

    // figure out eT bin
    UInt_t iEtBin(5);
    for (UInt_t iBin = 1; iBin < NTrgBins; iBin++) {
      const Bool_t isInBin = ((eTtrg >= eTbinMin[iBin]) && (eTtrg < eTbinMax[iBin]));
      if (isInBin) {
        iEtBin = iBin;
        break;
      }
    }

    // determine species
    if (isInEtCut && isInPi0cut) {
      hTrgTsp[0]           -> Fill(tspTrg);
      hTrgEtBin[0][0]      -> Fill(eTtrg);
      hTrgEtBin[iEtBin][0] -> Fill(eTtrg);
      nTrgPi0[0]++;
      nTrgPi0[iEtBin]++;
    }
    if (isInEtCut && isInGamCut) {
      hTrgTsp[1]           -> Fill(tspTrg);
      hTrgEtBin[0][1]      -> Fill(eTtrg);
      hTrgEtBin[iEtBin][1] -> Fill(eTtrg);
      nTrgGam[0]++;
      nTrgGam[iEtBin]++;
    }

    // fill other histograms
    hTrgEt[0] -> Fill(eTtrg);
    if (isInTspCut) {
      hTrgEta      -> Fill(hPhys);
      hTrgPhi      -> Fill(fTrg);
      hTrgEtVsEta  -> Fill(hPhys, eTtrg);
      hTrgEtVsPhi  -> Fill(fTrg, eTtrg);
      hTrgEtaVsPhi -> Fill(fTrg, hPhys);
    }

    if (isInEtCut) {
      hTrgTsp[2] -> Fill(tspTrg);
    } else {
      hTrgEt[1] -> Fill(eTtrg);
    }

    // leave only pi0 and gamma candidates
    if (!isInEtCut || !isInTspCut) continue;

    // track loop
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFitTrk  = PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPossTrk = PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFitTrk  = (Double_t) nFitTrk / (Double_t) nPossTrk;
      const Double_t dcaTrk   = PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk     = PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk     = PrimaryTrackArray_phi[iTrk];
      const Double_t pTtrk    = PrimaryTrackArray_pT[iTrk];

      Double_t dFtrk = fTrk - fTrg;
      if (dFtrk < (-1. * TMath::PiOver2())) dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))  dFtrk -= TMath::TwoPi();

      // track cuts
      const Bool_t isInFitCut    = (nFitTrk >= nFitMin);
      const Bool_t isInRatioCut  = (rFitTrk >= rFitMin);
      const Bool_t isInDcaCut    = (dcaTrk < dcaMax);
      const Bool_t isInEtaTrkCut = (TMath::Abs(hTrk) < hTrkMax);
      const Bool_t isInPtCut     = ((pTtrk > pTtrkMin) && (pTtrk < pTtrkMax));

      // nFit cut
      hTrkNfit[0]  -> Fill(nFitTrk);
      hTrkNfitVsPt -> Fill(pTtrk, nFitTrk);
      if (!isInFitCut) {
        hTrkNfit[1] -> Fill(nFitTrk);
      }

      // nFit / nPoss cut
      hTrkRfit[0] -> Fill(rFitTrk);
      if (!isInRatioCut) {
        hTrkRfit[1] -> Fill(rFitTrk);
      }

      // dca cut
      hTrkDca[0]  -> Fill(dcaTrk);
      hTrkDcaVsPt -> Fill(pTtrk, dcaTrk);
      if (!isInDcaCut) {
        hTrkDca[1] -> Fill(dcaTrk);
      }

      // eta cut
      hTrkEta[0] -> Fill(hTrk);
      if (!isInEtaTrkCut) {
        hTrkEta[1] -> Fill(hTrk);
      }

      // pT cut
      hTrkPt[0] -> Fill(pTtrk);
      if (!isInPtCut) {
        hTrkPt[1] -> Fill(pTtrk);
      }

      // apply cuts
      if (!isInFitCut || !isInRatioCut || !isInDcaCut || !isInEtaTrkCut) continue;


      // for pT-differential dF histograms
      Bool_t isInTrkPtBin[NTrkBins];
      for (UInt_t iTrkBin = 0; iTrkBin < NTrkBins; iTrkBin++) {
        isInTrkPtBin[iTrkBin] = ((pTtrk > pTbinMin[iTrkBin]) && (pTtrk < pTbinMax[iTrkBin]));
      }

      // fill histograms
      if (isInPi0cut) {
        hTrkPtBin[0]   -> Fill(pTtrk);
        hTrkPtVsDf[0]  -> Fill(dFtrk, pTtrk);
        hTrkPtVsEta[0] -> Fill(hTrk, pTtrk);
        if (isInPtCut) {
          hTrkDfBin[0]   -> Fill(dFtrk);
          hTrkEtaBin[0]  -> Fill(hTrk);
          hTrkDfVsEta[0] -> Fill(hTrk, dFtrk);
        }
        for (UInt_t iTrkBin = 0; iTrkBin < NTrkBins; iTrkBin++) {
          if (isInTrkPtBin[iTrkBin]) {
            hTrkDfPtBin[0][iTrkBin]  -> Fill(dFtrk);
            hTrkEtaPtBin[0][iTrkBin] -> Fill(hTrk);
          }
        }
      }
      if (isInGamCut) {
        hTrkPtBin[1]   -> Fill(pTtrk);
        hTrkPtVsDf[1]  -> Fill(dFtrk, pTtrk);
        hTrkPtVsEta[1] -> Fill(hTrk, pTtrk);
        if (isInPtCut) {
          hTrkDfBin[1]   -> Fill(dFtrk);
          hTrkEtaBin[1]  -> Fill(hTrk);
          hTrkDfVsEta[1] -> Fill(hTrk, dFtrk);
        }
        for (UInt_t iTrkBin = 0; iTrkBin < NTrkBins; iTrkBin++) {
          if (isInTrkPtBin[iTrkBin]) {
            hTrkDfPtBin[1][iTrkBin]  -> Fill(dFtrk);
            hTrkEtaPtBin[1][iTrkBin] -> Fill(hTrk);
          }
        }
      }
    }  // end track loop
  }  // end event loop

  // report no. of triggers
  cout << "      Event loop finished:" << endl;
  for (UInt_t iCut = 0; iCut < NTrgCuts; iCut++) {
    hEvtNum -> SetBinContent(iCut + 1, nTrgCut[iCut]);
    cout << "        " << nTrgCut[iCut] << " evts. left after cut " << iCut << endl;
  }

  // report no. of pi0's and gamma's
  cout << "      Number of triggers:" << endl;
  for (UInt_t iBin = 0; iBin < NTrgBins; iBin++) {
    cout << "        eTtrg = (" << eTbinMin[iBin] << ", " << eTbinMax[iBin]
         << ") GeV: nPi0 = " << nTrgPi0[iBin] << ", nGam = " << nTrgGam[iBin]
         << endl;
  }

  // sum trigger histograms
  TH1D *hTrgEtSum = (TH1D*) hTrgEtBin[0][0] -> Clone();
  hTrgEtSum -> SetName("hTrgEtSum");
  hTrgEtSum -> Add(hTrgEtBin[0][1]);

  // create profiles
  TProfile *pTrkNfitVsPt = hTrkNfitVsPt -> ProfileX("pTrkNfitVsPt", 1, -1, "S");
  TProfile *pTrkDcaVsPt  = hTrkDcaVsPt  -> ProfileX("pTrkDcaVsPt", 1, -1, "S");

  // normalize relevant histograms
  const Double_t pTbin2 = (pt[1] - pt[0]) / (Double_t) nPt2;
  const Double_t pTbin3 = (pt[1] - pt[0]) / (Double_t) nPt3;
  const Double_t dFbin  = (df[1] - df[0]) / (Double_t) nDf;
  const Double_t etaBin = (eta[1] - eta[0]) / (Double_t) nEta;
  for (UInt_t iTsp = 0; iTsp < NTrgTsp; iTsp++) {

    // no. of triggers
    UInt_t nTrg(1);
    if (iTsp == 0) nTrg = nTrgPi0[0];
    if (iTsp == 1) nTrg = nTrgGam[0];

    // norms
    const Double_t trgNorm   = 1. / (Double_t) nTrg;
    const Double_t pTnorm    = 1. / (pTbin2 * nTrg);
    const Double_t dFnorm    = 1. / (dFbin * nTrg);
    const Double_t etaNorm   = 1. / (etaBin * nTrg);
    const Double_t pTdFnorm  = 1. / (pTbin3 * dFbin);
    const Double_t pTetaNorm = 1. / (pTbin3 * etaBin);
    const Double_t etaDfNorm = 1. / (etaBin * dFbin);

    // normalize histograms
    hTrkPtBin[iTsp]   -> Scale(pTnorm);
    hTrkDfBin[iTsp]   -> Scale(dFnorm);
    hTrkEtaBin[iTsp]  -> Scale(etaNorm);
    hTrkPtVsDf[iTsp]  -> Scale(trgNorm);
    hTrkPtVsDf[iTsp]  -> Scale(pTdFnorm);
    hTrkPtVsEta[iTsp] -> Scale(trgNorm);
    hTrkPtVsEta[iTsp] -> Scale(pTetaNorm);
    hTrkDfVsEta[iTsp] -> Scale(trgNorm);
    hTrkDfVsEta[iTsp] -> Scale(etaDfNorm);
    for (UInt_t iTrkBin = 0; iTrkBin < NTrkBins; iTrkBin++) {
      hTrkDfPtBin[iTsp][iTrkBin]  -> Scale(dFnorm);
      hTrkEtaPtBin[iTsp][iTrkBin] -> Scale(etaNorm);
    }
  }
  cout << "      Normalized certain histograms." << endl;

  // find minimum / maximum for 2d comparison plots
  Float_t zVxyVsVzPlot[2];
  Float_t zEtVsEtaPhiPlot[2];
  Float_t zPtVsDfPlot[2];
  Float_t zPtVsEtaPlot[2];
  Float_t zEtaVsPhiPlot[2];

  // Vxy vs. Vz
  const Double_t minVx = hEvtVzVsVx -> GetMinimum(0.);
  const Double_t minVy = hEvtVzVsVy -> GetMinimum(0.);
  const Double_t maxVx = hEvtVzVsVx -> GetMaximum();
  const Double_t maxVy = hEvtVzVsVy -> GetMaximum();
  zVxyVsVzPlot[0] = TMath::Min(minVx, minVy);
  zVxyVsVzPlot[1] = TMath::Max(maxVx, maxVy);

  // eTtrg vs. eta, phi
  const Double_t minEta = hTrgEtVsEta -> GetMinimum(0.);
  const Double_t minPhi = hTrgEtVsPhi -> GetMinimum(0.);
  const Double_t maxEta = hTrgEtVsEta -> GetMaximum();
  const Double_t maxPhi = hTrgEtVsPhi -> GetMaximum();
  zEtVsEtaPhiPlot[0] = TMath::Min(minEta, minPhi);
  zEtVsEtaPhiPlot[1] = TMath::Max(maxEta, maxPhi);

  // pTtrk vs. dFtrk
  const Double_t minDfP = hTrkPtVsDf[0] -> GetMinimum(0.);
  const Double_t minDfG = hTrkPtVsDf[1] -> GetMinimum(0.);
  const Double_t maxDfP = hTrkPtVsDf[0] -> GetMaximum();
  const Double_t maxDfG = hTrkPtVsDf[1] -> GetMaximum();
  zPtVsDfPlot[0]  = TMath::Min(minDfP, minDfG);
  zPtVsDfPlot[1]  = TMath::Max(maxDfP, maxDfG);

  // pTtrk vs. hTrk
  const Double_t minEtaP = hTrkPtVsEta[0] -> GetMinimum(0.);
  const Double_t minEtaG = hTrkPtVsEta[1] -> GetMinimum(0.);
  const Double_t maxEtaP = hTrkPtVsEta[0] -> GetMaximum();
  const Double_t maxEtaG = hTrkPtVsEta[1] -> GetMaximum();
  zPtVsEtaPlot[0]  = TMath::Min(minEtaP, minEtaG);
  zPtVsEtaPlot[1]  = TMath::Max(maxEtaP, maxEtaG);

  // hTrk vs. dFtrk
  const Double_t minAngP = hTrkDfVsEta[0] -> GetMinimum(0.);
  const Double_t minAngG = hTrkDfVsEta[1] -> GetMinimum(0.);
  const Double_t maxAngP = hTrkDfVsEta[0] -> GetMaximum();
  const Double_t maxAngG = hTrkDfVsEta[1] -> GetMaximum();
  zEtaVsPhiPlot[0] = TMath::Min(minAngP, minAngG);
  zEtaVsPhiPlot[1] = TMath::Max(maxAngP, maxAngG);
  cout << "      Determined z-axis ranges for 2d comparisons:\n"
       << "        Vz vs. Vxy:     min = " << zVxyVsVzPlot[0]    << ", max = " << zVxyVsVzPlot[1]    << "\n"
       << "        eT vs. eta/phi: min = " << zEtVsEtaPhiPlot[0] << ", max = " << zEtVsEtaPhiPlot[1] << "\n"
       << "        pT vs. dF:      min = " << zPtVsDfPlot[0]     << ", max = " << zPtVsDfPlot[1]     << "\n"
       << "        pT vs. eta:     min = " << zPtVsEtaPlot[0]    << ", max = " << zPtVsEtaPlot[1]    << "\n"
       << "        dF vs. eta:     min = " << zEtaVsPhiPlot[0]   << ", max = " << zEtaVsPhiPlot[1]
       << endl;

  // set styles
  const UInt_t  fColAll(923);
  const UInt_t  fColCut(899);
  const UInt_t  fLinAll(1);
  const UInt_t  fLinCut(1);
  const UInt_t  fFilAll(0);
  const UInt_t  fFilCut(3345);
  const UInt_t  fMarAll(8);
  const UInt_t  fMarCut(8);
  const UInt_t  fMarPi0(20);
  const UInt_t  fMarGam(21);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColTrg[NTrgTsp]     = {859, 899};
  const UInt_t  fColPi0[NTrgBins]    = {923, 859, 839, 819};
  const UInt_t  fColGam[NTrgBins]    = {923, 899, 879, 859};
  const UInt_t  fColTrkPi0[NTrkBins] = {590, 591, 593, 596};
  const UInt_t  fColTrkGam[NTrkBins] = {622, 623, 625, 628};
  const UInt_t  fMarTrkPi0[NTrkBins] = {20, 22, 23, 21};
  const UInt_t  fMarTrkGam[NTrkBins] = {20, 22, 23, 21};
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.3);
  const Float_t fOffZ(1.);
  const Float_t fOffL(0.007);
  const Float_t xVz2dPlot[2] = {-55., 55.};
  const Float_t xEtPlot[2]   = {7., 23.};
  const Float_t yEtPlot[2]   = {7., 37777.};
  const TString sEvt("events");
  const TString sCount("counts");
  const TString sPrim("N_{primary}");
  const TString sTrgEt("E_{T}^{trg} [GeV]");
  const TString sTrgEta("#eta^{trg}");
  const TString sTrgPhi("#varphi^{trg}");
  const TString sTrgTsp("TSP");
  const TString sTrkNfit("N_{fit}");
  const TString sTrkRfit("N_{fit} / N_{poss}");
  const TString sTrkDca("DCA [cm]");
  const TString sTrkEta("#eta^{trk}");
  const TString sTrkPt("p_{T}^{trk} [GeV/c]");
  const TString sTrkDf("#Delta#varphi^{trk}");
  const TString sTrkPtY("(1/N^{trg}) dN^{trk}/dp_{T}^{trk} [GeV/c]^{-1}");
  const TString sTrkDfY("(1/N^{trg}) dN^{trk}/d#Delta#varphi^{trk}");
  const TString sTrkEtaY("(1/N^{trg}) dN^{trk}/d#eta^{trk}");
  const TString sVtx[NVtx]     = {"V_{x} [cm]", "V_{y} [cm]", "V_{z} [cm]", "V_{r} [cm]"};
  const TString sCut[NTrgCuts] = {"no cuts", "bad runs removed", "V_{z}, V_{r} cuts", "bad towers removed", "e_{#eta}, e_{#varphi} cuts", "P_{proj} cut", "#eta^{trg} cut", "E_{T}^{trg} cut", "TSP cut"};

  // event no's
  const UInt_t  fNumCol(923);
  const UInt_t  fNumLin(1);
  const UInt_t  fNumFil(1001);
  const Float_t fBar(0.6);
  const Float_t fOffB(0.2);
  hEvtNum -> SetLineColor(fNumCol);
  hEvtNum -> SetLineStyle(fNumLin);
  hEvtNum -> SetFillColor(fNumCol);
  hEvtNum -> SetFillStyle(fNumFil);
  hEvtNum -> SetMarkerColor(fNumCol);
  hEvtNum -> SetMarkerStyle(fNumLin);
  hEvtNum -> SetBarWidth(fBar);
  hEvtNum -> SetBarOffset(fOffB);
  hEvtNum -> SetTitle("");
  hEvtNum -> SetTitleFont(fTxt);
  hEvtNum -> GetXaxis() -> SetTitle("");
  hEvtNum -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtNum -> GetXaxis() -> SetTitleSize(fTit);
  hEvtNum -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtNum -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtNum -> GetXaxis() -> SetLabelSize(fLab);
  hEvtNum -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtNum -> GetXaxis() -> CenterTitle(fCnt);
  hEvtNum -> GetYaxis() -> SetTitle(sEvt.Data());
  hEvtNum -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtNum -> GetYaxis() -> SetTitleSize(fTit);
  hEvtNum -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtNum -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtNum -> GetYaxis() -> SetLabelSize(fLab);
  hEvtNum -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtNum -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iBin = 1; iBin < (NTrgCuts + 1); iBin++) {
    hEvtNum -> GetXaxis() -> SetBinLabel(iBin, sCut[iBin - 1].Data());
  }

  // vertices
  hEvtVz[0]  -> SetLineColor(fColAll);
  hEvtVz[0]  -> SetLineStyle(fLinAll);
  hEvtVz[0]  -> SetFillColor(fColAll);
  hEvtVz[0]  -> SetFillStyle(fFilAll);
  hEvtVz[0]  -> SetMarkerColor(fColAll);
  hEvtVz[0]  -> SetMarkerStyle(fMarAll);
  hEvtVz[0]  -> SetTitle("");
  hEvtVz[0]  -> SetTitleFont(fTxt);
  hEvtVz[0]  -> GetXaxis() -> SetTitle(sVtx[2].Data());
  hEvtVz[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVz[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVz[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVz[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVz[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVz[0]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVz[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVz[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVz[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVz[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVz[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVz[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVz[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVz[0]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVz[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVz[1]  -> SetLineColor(fColCut);
  hEvtVz[1]  -> SetLineStyle(fLinCut);
  hEvtVz[1]  -> SetFillColor(fColCut);
  hEvtVz[1]  -> SetFillStyle(fFilCut);
  hEvtVz[1]  -> SetMarkerColor(fColCut);
  hEvtVz[1]  -> SetMarkerStyle(fMarCut);
  hEvtVz[1]  -> SetTitle("");
  hEvtVz[1]  -> SetTitleFont(fTxt);
  hEvtVz[1]  -> GetXaxis() -> SetTitle(sVtx[2].Data());
  hEvtVz[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVz[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVz[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVz[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVz[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVz[1]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVz[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVz[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVz[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVz[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVz[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVz[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVz[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVz[1]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVz[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVr[0]  -> SetLineColor(fColAll);
  hEvtVr[0]  -> SetLineStyle(fLinAll);
  hEvtVr[0]  -> SetFillColor(fColAll);
  hEvtVr[0]  -> SetFillStyle(fFilAll);
  hEvtVr[0]  -> SetMarkerColor(fColAll);
  hEvtVr[0]  -> SetMarkerStyle(fMarAll);
  hEvtVr[0]  -> SetTitle("");
  hEvtVr[0]  -> SetTitleFont(fTxt);
  hEvtVr[0]  -> GetXaxis() -> SetTitle(sVtx[3].Data());
  hEvtVr[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVr[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVr[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVr[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVr[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVr[0]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVr[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVr[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVr[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVr[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVr[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVr[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVr[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVr[0]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVr[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVr[1]  -> SetLineColor(fColCut);
  hEvtVr[1]  -> SetLineStyle(fLinCut);
  hEvtVr[1]  -> SetFillColor(fColCut);
  hEvtVr[1]  -> SetFillStyle(fFilCut);
  hEvtVr[1]  -> SetMarkerColor(fColCut);
  hEvtVr[1]  -> SetMarkerStyle(fMarCut);
  hEvtVr[1]  -> SetTitle("");
  hEvtVr[1]  -> SetTitleFont(fTxt);
  hEvtVr[1]  -> GetXaxis() -> SetTitle(sVtx[3].Data());
  hEvtVr[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVr[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVr[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVr[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVr[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVr[1]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVr[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVr[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtVr[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVr[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVr[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVr[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVr[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVr[1]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVr[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVyVsVx -> SetTitle("");
  hEvtVyVsVx -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetXaxis() -> SetTitle(sVtx[0].Data());
  hEvtVyVsVx -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVyVsVx -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVyVsVx -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVx -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVyVsVx -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVyVsVx -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVyVsVx -> GetYaxis() -> SetTitle(sVtx[1].Data());
  hEvtVyVsVx -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVyVsVx -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVyVsVx -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVx -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVyVsVx -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVyVsVx -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVyVsVx -> GetZaxis() -> SetTitle("");
  hEvtVyVsVx -> GetZaxis() -> SetTitleFont(fTxt);
  hEvtVyVsVx -> GetZaxis() -> SetTitleSize(fTit);
  hEvtVyVsVx -> GetZaxis() -> SetTitleOffset(fOffZ);
  hEvtVyVsVx -> GetZaxis() -> SetLabelFont(fTxt);
  hEvtVyVsVx -> GetZaxis() -> SetLabelSize(fLab);
  hEvtVyVsVx -> GetZaxis() -> CenterTitle(fCnt);
  hEvtVzVsVx -> SetTitle("");
  hEvtVzVsVx -> SetTitleFont(fTxt);
  hEvtVzVsVx -> GetXaxis() -> SetTitle(sVtx[0].Data());
  hEvtVzVsVx -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVzVsVx -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVzVsVx -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVzVsVx -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVzVsVx -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVzVsVx -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVzVsVx -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVzVsVx -> GetYaxis() -> SetTitle(sVtx[2].Data());
  hEvtVzVsVx -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVzVsVx -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVzVsVx -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVzVsVx -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVzVsVx -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVzVsVx -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVzVsVx -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVzVsVx -> GetYaxis() -> SetRangeUser(xVz2dPlot[0], xVz2dPlot[1]);
  hEvtVzVsVx -> GetZaxis() -> SetTitle("");
  hEvtVzVsVx -> GetZaxis() -> SetTitleFont(fTxt);
  hEvtVzVsVx -> GetZaxis() -> SetTitleSize(fTit);
  hEvtVzVsVx -> GetZaxis() -> SetTitleOffset(fOffZ);
  hEvtVzVsVx -> GetZaxis() -> SetLabelFont(fTxt);
  hEvtVzVsVx -> GetZaxis() -> SetLabelSize(fLab);
  hEvtVzVsVx -> GetZaxis() -> CenterTitle(fCnt);
  hEvtVzVsVx -> GetZaxis() -> SetRangeUser(zVxyVsVzPlot[0], zVxyVsVzPlot[1]);
  hEvtVzVsVy -> SetTitle("");
  hEvtVzVsVy -> SetTitleFont(fTxt);
  hEvtVzVsVy -> GetXaxis() -> SetTitle(sVtx[1].Data());
  hEvtVzVsVy -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtVzVsVy -> GetXaxis() -> SetTitleSize(fTit);
  hEvtVzVsVy -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtVzVsVy -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtVzVsVy -> GetXaxis() -> SetLabelSize(fLab);
  hEvtVzVsVy -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtVzVsVy -> GetXaxis() -> CenterTitle(fCnt);
  hEvtVzVsVy -> GetYaxis() -> SetTitle(sVtx[2].Data());
  hEvtVzVsVy -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtVzVsVy -> GetYaxis() -> SetTitleSize(fTit);
  hEvtVzVsVy -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtVzVsVy -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtVzVsVy -> GetYaxis() -> SetLabelSize(fLab);
  hEvtVzVsVy -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtVzVsVy -> GetYaxis() -> CenterTitle(fCnt);
  hEvtVzVsVy -> GetYaxis() -> SetRangeUser(xVz2dPlot[0], xVz2dPlot[1]);
  hEvtVzVsVy -> GetZaxis() -> SetTitle("");
  hEvtVzVsVy -> GetZaxis() -> SetTitleFont(fTxt);
  hEvtVzVsVy -> GetZaxis() -> SetTitleSize(fTit);
  hEvtVzVsVy -> GetZaxis() -> SetTitleOffset(fOffZ);
  hEvtVzVsVy -> GetZaxis() -> SetLabelFont(fTxt);
  hEvtVzVsVy -> GetZaxis() -> SetLabelSize(fLab);
  hEvtVzVsVy -> GetZaxis() -> CenterTitle(fCnt);
  hEvtVzVsVx -> GetZaxis() -> SetRangeUser(zVxyVsVzPlot[0], zVxyVsVzPlot[1]);

  // no. of tracks
  hEvtPrim -> SetLineColor(fColAll);
  hEvtPrim -> SetLineStyle(fLinAll);
  hEvtPrim -> SetFillColor(fColAll);
  hEvtPrim -> SetFillStyle(fFilAll);
  hEvtPrim -> SetMarkerColor(fColAll);
  hEvtPrim -> SetMarkerStyle(fMarAll);
  hEvtPrim -> SetTitle("");
  hEvtPrim -> SetTitleFont(fTxt);
  hEvtPrim -> GetXaxis() -> SetTitle(sPrim.Data());
  hEvtPrim -> GetXaxis() -> SetTitleFont(fTxt);
  hEvtPrim -> GetXaxis() -> SetTitleSize(fTit);
  hEvtPrim -> GetXaxis() -> SetTitleOffset(fOffX);
  hEvtPrim -> GetXaxis() -> SetLabelFont(fTxt);
  hEvtPrim -> GetXaxis() -> SetLabelSize(fLab);
  hEvtPrim -> GetXaxis() -> SetLabelOffset(fOffL);
  hEvtPrim -> GetXaxis() -> CenterTitle(fCnt);
  hEvtPrim -> GetYaxis() -> SetTitle(sCount.Data());
  hEvtPrim -> GetYaxis() -> SetTitleFont(fTxt);
  hEvtPrim -> GetYaxis() -> SetTitleSize(fTit);
  hEvtPrim -> GetYaxis() -> SetTitleOffset(fOffY);
  hEvtPrim -> GetYaxis() -> SetLabelFont(fTxt);
  hEvtPrim -> GetYaxis() -> SetLabelSize(fLab);
  hEvtPrim -> GetYaxis() -> SetLabelOffset(fOffL);
  hEvtPrim -> GetYaxis() -> CenterTitle(fCnt);

  // eTtrg (no tsp cuts)
  hTrgEt[0]   -> SetLineColor(fColAll);
  hTrgEt[0]   -> SetLineStyle(fLinAll);
  hTrgEt[0]   -> SetFillColor(fColAll);
  hTrgEt[0]   -> SetFillStyle(fFilAll);
  hTrgEt[0]   -> SetMarkerColor(fColAll);
  hTrgEt[0]   -> SetMarkerStyle(fMarAll);
  hTrgEt[0]   -> SetTitle("");
  hTrgEt[0]   -> SetTitleFont(fTxt);
  hTrgEt[0]   -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEt[0]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEt[0]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEt[0]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEt[0]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEt[0]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEt[0]   -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEt[0]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEt[0]   -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEt[0]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEt[0]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEt[0]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEt[0]   -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEt[0]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEt[0]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEt[0]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEt[1]   -> SetLineColor(fColCut);
  hTrgEt[1]   -> SetLineStyle(fLinCut);
  hTrgEt[1]   -> SetFillColor(fColCut);
  hTrgEt[1]   -> SetFillStyle(fFilCut);
  hTrgEt[1]   -> SetMarkerColor(fColCut);
  hTrgEt[1]   -> SetMarkerStyle(fMarCut);
  hTrgEt[1]   -> SetTitle("");
  hTrgEt[1]   -> SetTitleFont(fTxt);
  hTrgEt[1]   -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEt[1]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEt[1]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEt[1]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEt[1]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEt[1]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEt[1]   -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEt[1]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEt[1]   -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEt[1]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEt[1]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEt[1]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEt[1]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEt[1]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEt[1]   -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEt[1]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> SetTitle("");
  hTrgEtVsEta -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetXaxis() -> SetTitle(sTrgEta.Data());
  hTrgEtVsEta -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtVsEta -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtVsEta -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtVsEta -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtVsEta -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEtVsEta -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> GetYaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtVsEta -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtVsEta -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtVsEta -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtVsEta -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtVsEta -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEtVsEta -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> GetZaxis() -> SetTitle("");
  hTrgEtVsEta -> GetZaxis() -> SetTitleFont(fTxt);
  hTrgEtVsEta -> GetZaxis() -> SetTitleSize(fTit);
  hTrgEtVsEta -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrgEtVsEta -> GetZaxis() -> SetLabelFont(fTxt);
  hTrgEtVsEta -> GetZaxis() -> SetLabelSize(fLab);
  hTrgEtVsEta -> GetZaxis() -> CenterTitle(fCnt);
  hTrgEtVsEta -> GetZaxis() -> SetRangeUser(zEtVsEtaPhiPlot[0], zEtVsEtaPhiPlot[1]);
  hTrgEtVsPhi -> SetTitle("");
  hTrgEtVsPhi -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetXaxis() -> SetTitle(sTrgPhi.Data());
  hTrgEtVsPhi -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtVsPhi -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtVsPhi -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtVsPhi -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtVsPhi -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEtVsPhi -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtVsPhi -> GetYaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtVsPhi -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtVsPhi -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtVsPhi -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtVsPhi -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtVsPhi -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEtVsPhi -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtVsPhi -> GetZaxis() -> SetTitle("");
  hTrgEtVsPhi -> GetZaxis() -> SetTitleFont(fTxt);
  hTrgEtVsPhi -> GetZaxis() -> SetTitleSize(fTit);
  hTrgEtVsPhi -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrgEtVsPhi -> GetZaxis() -> SetLabelFont(fTxt);
  hTrgEtVsPhi -> GetZaxis() -> SetLabelSize(fLab);
  hTrgEtVsPhi -> GetZaxis() -> CenterTitle(fCnt);
  hTrgEtVsPhi -> GetZaxis() -> SetRangeUser(zEtVsEtaPhiPlot[0], zEtVsEtaPhiPlot[1]);

  // eTtrg (tsp cuts)
  for (UInt_t iBin = 1; iBin < NTrgBins; iBin++) {
    hTrgEtBin[iBin][0] -> SetLineColor(fColPi0[iBin]);
    hTrgEtBin[iBin][0] -> SetLineStyle(fLinAll);
    hTrgEtBin[iBin][0] -> SetFillColor(fColPi0[iBin]);
    hTrgEtBin[iBin][0] -> SetFillStyle(fFilAll);
    hTrgEtBin[iBin][0] -> SetMarkerColor(fColPi0[iBin]);
    hTrgEtBin[iBin][0] -> SetMarkerStyle(fMarPi0);
    hTrgEtBin[iBin][0] -> SetTitle("");
    hTrgEtBin[iBin][0] -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitle(sTrgEt.Data());
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetLabelOffset(fOffL);
    hTrgEtBin[iBin][0] -> GetXaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][0] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitle(sCount.Data());
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hTrgEtBin[iBin][0] -> GetYaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][0] -> GetYaxis() -> SetRangeUser(yEtPlot[0], yEtPlot[1]);
    hTrgEtBin[iBin][1] -> SetLineColor(fColGam[iBin]);
    hTrgEtBin[iBin][1] -> SetLineStyle(fLinAll);
    hTrgEtBin[iBin][1] -> SetFillColor(fColGam[iBin]);
    hTrgEtBin[iBin][1] -> SetFillStyle(fFilAll);
    hTrgEtBin[iBin][1] -> SetMarkerColor(fColGam[iBin]);
    hTrgEtBin[iBin][1] -> SetMarkerStyle(fMarGam);
    hTrgEtBin[iBin][1] -> SetTitle("");
    hTrgEtBin[iBin][1] -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitle(sTrgEt.Data());
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetLabelOffset(fOffL);
    hTrgEtBin[iBin][1] -> GetXaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][1] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitle(sCount.Data());
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitleSize(fTit);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetLabelSize(fLab);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetLabelOffset(fOffY);
    hTrgEtBin[iBin][1] -> GetYaxis() -> CenterTitle(fCnt);
    hTrgEtBin[iBin][1] -> GetYaxis() -> SetRangeUser(yEtPlot[0], yEtPlot[1]);
  }
  hTrgEtBin[0][0] -> SetLineColor(fColTrg[0]);
  hTrgEtBin[0][0] -> SetLineStyle(fLinAll);
  hTrgEtBin[0][0] -> SetFillColor(fColTrg[0]);
  hTrgEtBin[0][0] -> SetFillStyle(fFilAll);
  hTrgEtBin[0][0] -> SetMarkerColor(fColTrg[0]);
  hTrgEtBin[0][0] -> SetMarkerStyle(fMarPi0);
  hTrgEtBin[0][0] -> SetTitle("");
  hTrgEtBin[0][0] -> SetTitleFont(fTxt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtBin[0][0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][0] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEtBin[0][0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtBin[0][0] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtBin[0][0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][0] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEtBin[0][0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtBin[0][1] -> SetLineColor(fColTrg[1]);
  hTrgEtBin[0][1] -> SetLineStyle(fLinAll);
  hTrgEtBin[0][1] -> SetFillColor(fColTrg[1]);
  hTrgEtBin[0][1] -> SetFillStyle(fFilAll);
  hTrgEtBin[0][1] -> SetMarkerColor(fColTrg[1]);
  hTrgEtBin[0][1] -> SetMarkerStyle(fMarGam);
  hTrgEtBin[0][1] -> SetTitle("");
  hTrgEtBin[0][1] -> SetTitleFont(fTxt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtBin[0][1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][1] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEtBin[0][1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtBin[0][1] -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtBin[0][1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtBin[0][1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtBin[0][1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtBin[0][1] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEtBin[0][1] -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtSum       -> SetLineColor(fColAll);
  hTrgEtSum       -> SetLineStyle(fLinAll);
  hTrgEtSum       -> SetFillColor(fColAll);
  hTrgEtSum       -> SetFillStyle(fFilAll);
  hTrgEtSum       -> SetMarkerColor(fColAll);
  hTrgEtSum       -> SetMarkerStyle(fMarAll);
  hTrgEtSum       -> SetTitle("");
  hTrgEtSum       -> SetTitleFont(fTxt);
  hTrgEtSum       -> GetXaxis() -> SetTitle(sTrgEt.Data());
  hTrgEtSum       -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtSum       -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtSum       -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtSum       -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtSum       -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtSum       -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEtSum       -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtSum       -> GetXaxis() -> SetRangeUser(xEtPlot[0], xEtPlot[1]);
  hTrgEtSum       -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgEtSum       -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtSum       -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtSum       -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtSum       -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtSum       -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtSum       -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEtSum       -> GetYaxis() -> CenterTitle(fCnt);

  // tsp
  hTrgTsp[0] -> SetLineColor(fColTrg[0]);
  hTrgTsp[0] -> SetLineStyle(fLinCut);
  hTrgTsp[0] -> SetFillColor(fColTrg[0]);
  hTrgTsp[0] -> SetFillStyle(fFilCut);
  hTrgTsp[0] -> SetMarkerColor(fColTrg[0]);
  hTrgTsp[0] -> SetMarkerStyle(fMarCut);
  hTrgTsp[0] -> SetTitle("");
  hTrgTsp[0] -> SetTitleFont(fTxt);
  hTrgTsp[0] -> GetXaxis() -> SetTitle(sTrgTsp.Data());
  hTrgTsp[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgTsp[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgTsp[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgTsp[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgTsp[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgTsp[0] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgTsp[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgTsp[0] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgTsp[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgTsp[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgTsp[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgTsp[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgTsp[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgTsp[0] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgTsp[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrgTsp[1] -> SetLineColor(fColTrg[1]);
  hTrgTsp[1] -> SetLineStyle(fLinCut);
  hTrgTsp[1] -> SetFillColor(fColTrg[1]);
  hTrgTsp[1] -> SetFillStyle(fFilCut);
  hTrgTsp[1] -> SetMarkerColor(fColTrg[1]);
  hTrgTsp[1] -> SetMarkerStyle(fMarCut);
  hTrgTsp[1] -> SetTitle("");
  hTrgTsp[1] -> SetTitleFont(fTxt);
  hTrgTsp[1] -> GetXaxis() -> SetTitle(sTrgTsp.Data());
  hTrgTsp[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgTsp[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgTsp[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgTsp[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgTsp[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgTsp[1] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgTsp[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgTsp[1] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgTsp[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgTsp[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgTsp[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgTsp[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgTsp[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgTsp[1] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgTsp[1] -> GetYaxis() -> CenterTitle(fCnt);
  hTrgTsp[2] -> SetLineColor(fColAll);
  hTrgTsp[2] -> SetLineStyle(fLinAll);
  hTrgTsp[2] -> SetFillColor(fColAll);
  hTrgTsp[2] -> SetFillStyle(fFilAll);
  hTrgTsp[2] -> SetMarkerColor(fColAll);
  hTrgTsp[2] -> SetMarkerStyle(fMarAll);
  hTrgTsp[2] -> SetTitle("");
  hTrgTsp[2] -> SetTitleFont(fTxt);
  hTrgTsp[2] -> GetXaxis() -> SetTitle(sTrgTsp.Data());
  hTrgTsp[2] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgTsp[2] -> GetXaxis() -> SetTitleSize(fTit);
  hTrgTsp[2] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgTsp[2] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgTsp[2] -> GetXaxis() -> SetLabelSize(fLab);
  hTrgTsp[2] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgTsp[2] -> GetXaxis() -> CenterTitle(fCnt);
  hTrgTsp[2] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrgTsp[2] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgTsp[2] -> GetYaxis() -> SetTitleSize(fTit);
  hTrgTsp[2] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgTsp[2] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgTsp[2] -> GetYaxis() -> SetLabelSize(fLab);
  hTrgTsp[2] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgTsp[2] -> GetYaxis() -> CenterTitle(fCnt);

  // fTrg vs hTrg
  hTrgEtaVsPhi -> SetTitle("");
  hTrgEtaVsPhi -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetXaxis() -> SetTitle(sTrgPhi.Data());
  hTrgEtaVsPhi -> GetXaxis() -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetXaxis() -> SetTitleSize(fTit);
  hTrgEtaVsPhi -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrgEtaVsPhi -> GetXaxis() -> SetLabelFont(fTxt);
  hTrgEtaVsPhi -> GetXaxis() -> SetLabelSize(fLab);
  hTrgEtaVsPhi -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrgEtaVsPhi -> GetXaxis() -> CenterTitle(fCnt);
  hTrgEtaVsPhi -> GetYaxis() -> SetTitle(sTrgEta.Data());
  hTrgEtaVsPhi -> GetYaxis() -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetYaxis() -> SetTitleSize(fTit);
  hTrgEtaVsPhi -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrgEtaVsPhi -> GetYaxis() -> SetLabelFont(fTxt);
  hTrgEtaVsPhi -> GetYaxis() -> SetLabelSize(fLab);
  hTrgEtaVsPhi -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrgEtaVsPhi -> GetYaxis() -> CenterTitle(fCnt);
  hTrgEtaVsPhi -> GetZaxis() -> SetTitle("");
  hTrgEtaVsPhi -> GetZaxis() -> SetTitleFont(fTxt);
  hTrgEtaVsPhi -> GetZaxis() -> SetTitleSize(fTit);
  hTrgEtaVsPhi -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrgEtaVsPhi -> GetZaxis() -> SetLabelFont(fTxt);
  hTrgEtaVsPhi -> GetZaxis() -> SetLabelSize(fLab);
  hTrgEtaVsPhi -> GetZaxis() -> CenterTitle(fCnt);

  // nFit and nFit/nPoss
  hTrkNfit[0]  -> SetLineColor(fColAll);
  hTrkNfit[0]  -> SetLineStyle(fLinAll);
  hTrkNfit[0]  -> SetFillColor(fColAll);
  hTrkNfit[0]  -> SetFillStyle(fFilAll);
  hTrkNfit[0]  -> SetMarkerColor(fColAll);
  hTrkNfit[0]  -> SetMarkerStyle(fMarAll);
  hTrkNfit[0]  -> SetTitle("");
  hTrkNfit[0]  -> SetTitleFont(fTxt);
  hTrkNfit[0]  -> GetXaxis() -> SetTitle(sTrkNfit.Data());
  hTrkNfit[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkNfit[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkNfit[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkNfit[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkNfit[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkNfit[0]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkNfit[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkNfit[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkNfit[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkNfit[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkNfit[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkNfit[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkNfit[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkNfit[0]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkNfit[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkNfit[1]  -> SetLineColor(fColCut);
  hTrkNfit[1]  -> SetLineStyle(fLinCut);
  hTrkNfit[1]  -> SetFillColor(fColCut);
  hTrkNfit[1]  -> SetFillStyle(fFilCut);
  hTrkNfit[1]  -> SetMarkerColor(fColCut);
  hTrkNfit[1]  -> SetMarkerStyle(fMarCut);
  hTrkNfit[1]  -> SetTitle("");
  hTrkNfit[1]  -> SetTitleFont(fTxt);
  hTrkNfit[1]  -> GetXaxis() -> SetTitle(sTrkNfit.Data());
  hTrkNfit[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkNfit[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkNfit[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkNfit[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkNfit[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkNfit[1]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkNfit[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkNfit[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkNfit[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkNfit[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkNfit[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkNfit[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkNfit[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkNfit[1]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkNfit[1]  -> GetYaxis() -> CenterTitle(fCnt);
  pTrkNfitVsPt -> SetLineColor(fColAll);
  pTrkNfitVsPt -> SetLineStyle(fLinAll);
  pTrkNfitVsPt -> SetFillColor(fColAll);
  pTrkNfitVsPt -> SetFillStyle(fFilAll);
  pTrkNfitVsPt -> SetMarkerColor(fColAll);
  pTrkNfitVsPt -> SetMarkerStyle(fMarAll);
  pTrkNfitVsPt -> SetTitle("");
  pTrkNfitVsPt -> SetTitleFont(fTxt);
  pTrkNfitVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  pTrkNfitVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  pTrkNfitVsPt -> GetXaxis() -> SetTitleSize(fTit);
  pTrkNfitVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  pTrkNfitVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  pTrkNfitVsPt -> GetXaxis() -> SetLabelSize(fLab);
  pTrkNfitVsPt -> GetXaxis() -> SetLabelOffset(fOffL);
  pTrkNfitVsPt -> GetXaxis() -> CenterTitle(fCnt);
  pTrkNfitVsPt -> GetYaxis() -> SetTitle(sTrkNfit.Data());
  pTrkNfitVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  pTrkNfitVsPt -> GetYaxis() -> SetTitleSize(fTit);
  pTrkNfitVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  pTrkNfitVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  pTrkNfitVsPt -> GetYaxis() -> SetLabelSize(fLab);
  pTrkNfitVsPt -> GetYaxis() -> SetLabelOffset(fOffL);
  pTrkNfitVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkNfitVsPt -> SetTitle("");
  hTrkNfitVsPt -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkNfitVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetXaxis() -> SetTitleSize(fTit);
  hTrkNfitVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkNfitVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkNfitVsPt -> GetXaxis() -> SetLabelSize(fLab);
  hTrkNfitVsPt -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkNfitVsPt -> GetXaxis() -> CenterTitle(fCnt);
  hTrkNfitVsPt -> GetYaxis() -> SetTitle(sTrkNfit.Data());
  hTrkNfitVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetYaxis() -> SetTitleSize(fTit);
  hTrkNfitVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkNfitVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkNfitVsPt -> GetYaxis() -> SetLabelSize(fLab);
  hTrkNfitVsPt -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkNfitVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkNfitVsPt -> GetZaxis() -> SetTitle("");
  hTrkNfitVsPt -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkNfitVsPt -> GetZaxis() -> SetTitleSize(fTit);
  hTrkNfitVsPt -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkNfitVsPt -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkNfitVsPt -> GetZaxis() -> SetLabelSize(fLab);
  hTrkNfitVsPt -> GetZaxis() -> CenterTitle(fCnt);

  hTrkRfit[0] -> SetLineColor(fColAll);
  hTrkRfit[0] -> SetLineStyle(fLinAll);
  hTrkRfit[0] -> SetFillColor(fColAll);
  hTrkRfit[0] -> SetFillStyle(fFilAll);
  hTrkRfit[0] -> SetMarkerColor(fColAll);
  hTrkRfit[0] -> SetMarkerStyle(fMarAll);
  hTrkRfit[0] -> SetTitle("");
  hTrkRfit[0] -> SetTitleFont(fTxt);
  hTrkRfit[0] -> GetXaxis() -> SetTitle(sTrkRfit.Data());
  hTrkRfit[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkRfit[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkRfit[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkRfit[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkRfit[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkRfit[0] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkRfit[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkRfit[0] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkRfit[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkRfit[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkRfit[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkRfit[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkRfit[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkRfit[0] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkRfit[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkRfit[1] -> SetLineColor(fColCut);
  hTrkRfit[1] -> SetLineStyle(fLinCut);
  hTrkRfit[1] -> SetFillColor(fColCut);
  hTrkRfit[1] -> SetFillStyle(fFilCut);
  hTrkRfit[1] -> SetMarkerColor(fColCut);
  hTrkRfit[1] -> SetMarkerStyle(fMarCut);
  hTrkRfit[1] -> SetTitle("");
  hTrkRfit[1] -> SetTitleFont(fTxt);
  hTrkRfit[1] -> GetXaxis() -> SetTitle(sTrkRfit.Data());
  hTrkRfit[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkRfit[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkRfit[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkRfit[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkRfit[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkRfit[1] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkRfit[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkRfit[1] -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkRfit[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkRfit[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkRfit[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkRfit[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkRfit[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkRfit[1] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkRfit[1] -> GetYaxis() -> CenterTitle(fCnt);

  // dca
  hTrkDca[0]  -> SetLineColor(fColAll);
  hTrkDca[0]  -> SetLineStyle(fLinAll);
  hTrkDca[0]  -> SetFillColor(fColAll);
  hTrkDca[0]  -> SetFillStyle(fFilAll);
  hTrkDca[0]  -> SetMarkerColor(fColAll);
  hTrkDca[0]  -> SetMarkerStyle(fMarAll);
  hTrkDca[0]  -> SetTitle("");
  hTrkDca[0]  -> SetTitleFont(fTxt);
  hTrkDca[0]  -> GetXaxis() -> SetTitle(sTrkDca.Data());
  hTrkDca[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDca[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDca[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDca[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDca[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDca[0]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDca[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDca[0]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkDca[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDca[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDca[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDca[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDca[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDca[0]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDca[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDca[1]  -> SetLineColor(fColCut);
  hTrkDca[1]  -> SetLineStyle(fLinCut);
  hTrkDca[1]  -> SetFillColor(fColCut);
  hTrkDca[1]  -> SetFillStyle(fFilCut);
  hTrkDca[1]  -> SetMarkerColor(fColCut);
  hTrkDca[1]  -> SetMarkerStyle(fMarCut);
  hTrkDca[1]  -> SetTitle("");
  hTrkDca[1]  -> SetTitleFont(fTxt);
  hTrkDca[1]  -> GetXaxis() -> SetTitle(sTrkDca.Data());
  hTrkDca[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDca[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDca[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDca[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDca[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDca[1]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDca[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDca[1]  -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkDca[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDca[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDca[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDca[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDca[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDca[1]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDca[1]  -> GetYaxis() -> CenterTitle(fCnt);
  pTrkDcaVsPt -> SetLineColor(fColAll);
  pTrkDcaVsPt -> SetLineStyle(fLinAll);
  pTrkDcaVsPt -> SetFillColor(fColAll);
  pTrkDcaVsPt -> SetFillStyle(fFilAll);
  pTrkDcaVsPt -> SetMarkerColor(fColAll);
  pTrkDcaVsPt -> SetMarkerStyle(fMarAll);
  pTrkDcaVsPt -> SetTitle("");
  pTrkDcaVsPt -> SetTitleFont(fTxt);
  pTrkDcaVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  pTrkDcaVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  pTrkDcaVsPt -> GetXaxis() -> SetTitleSize(fTit);
  pTrkDcaVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  pTrkDcaVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  pTrkDcaVsPt -> GetXaxis() -> SetLabelSize(fLab);
  pTrkDcaVsPt -> GetXaxis() -> SetLabelOffset(fOffL);
  pTrkDcaVsPt -> GetXaxis() -> CenterTitle(fCnt);
  pTrkDcaVsPt -> GetYaxis() -> SetTitle(sTrkDca.Data());
  pTrkDcaVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  pTrkDcaVsPt -> GetYaxis() -> SetTitleSize(fTit);
  pTrkDcaVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  pTrkDcaVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  pTrkDcaVsPt -> GetYaxis() -> SetLabelSize(fLab);
  pTrkDcaVsPt -> GetYaxis() -> SetLabelOffset(fOffL);
  pTrkDcaVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDcaVsPt -> SetTitle("");
  hTrkDcaVsPt -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkDcaVsPt -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDcaVsPt -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDcaVsPt -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDcaVsPt -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDcaVsPt -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDcaVsPt -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDcaVsPt -> GetYaxis() -> SetTitle(sTrkDca.Data());
  hTrkDcaVsPt -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDcaVsPt -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDcaVsPt -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDcaVsPt -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDcaVsPt -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDcaVsPt -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDcaVsPt -> GetZaxis() -> SetTitle("");
  hTrkDcaVsPt -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkDcaVsPt -> GetZaxis() -> SetTitleSize(fTit);
  hTrkDcaVsPt -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkDcaVsPt -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkDcaVsPt -> GetZaxis() -> SetLabelSize(fLab);
  hTrkDcaVsPt -> GetZaxis() -> CenterTitle(fCnt);
  hTrkEta[1]     -> GetYaxis() -> SetTitleOffset(fOffY);

  // hTrk
  hTrkEta[0]     -> SetLineColor(fColAll);
  hTrkEta[0]     -> SetLineStyle(fLinAll);
  hTrkEta[0]     -> SetFillColor(fColAll);
  hTrkEta[0]     -> SetFillStyle(fFilAll);
  hTrkEta[0]     -> SetMarkerColor(fColAll);
  hTrkEta[0]     -> SetMarkerStyle(fMarAll);
  hTrkEta[0]     -> SetTitle("");
  hTrkEta[0]     -> SetTitleFont(fTxt);
  hTrkEta[0]     -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkEta[0]     -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkEta[0]     -> GetXaxis() -> SetTitleSize(fTit);
  hTrkEta[0]     -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkEta[0]     -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkEta[0]     -> GetXaxis() -> SetLabelSize(fLab);
  hTrkEta[0]     -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkEta[0]     -> GetXaxis() -> CenterTitle(fCnt);
  hTrkEta[0]     -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkEta[0]     -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkEta[0]     -> GetYaxis() -> SetTitleSize(fTit);
  hTrkEta[0]     -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkEta[0]     -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkEta[0]     -> GetYaxis() -> SetLabelSize(fLab);
  hTrkEta[0]     -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkEta[0]     -> GetYaxis() -> CenterTitle(fCnt);
  hTrkEta[1]     -> SetLineColor(fColCut);
  hTrkEta[1]     -> SetLineStyle(fLinCut);
  hTrkEta[1]     -> SetFillColor(fColCut);
  hTrkEta[1]     -> SetFillStyle(fFilCut);
  hTrkEta[1]     -> SetMarkerColor(fColCut);
  hTrkEta[1]     -> SetMarkerStyle(fMarCut);
  hTrkEta[1]     -> SetTitle("");
  hTrkEta[1]     -> SetTitleFont(fTxt);
  hTrkEta[1]     -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkEta[1]     -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkEta[1]     -> GetXaxis() -> SetTitleSize(fTit);
  hTrkEta[1]     -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkEta[1]     -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkEta[1]     -> GetXaxis() -> SetLabelSize(fLab);
  hTrkEta[1]     -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkEta[1]     -> GetXaxis() -> CenterTitle(fCnt);
  hTrkEta[1]     -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkEta[1]     -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkEta[1]     -> GetYaxis() -> SetTitleSize(fTit);
  hTrkEta[1]     -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkEta[1]     -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkEta[1]     -> GetYaxis() -> SetLabelSize(fLab);
  hTrkEta[1]     -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkEta[1]     -> GetYaxis() -> CenterTitle(fCnt);
  hTrkEtaBin[0]  -> SetLineColor(fColAll);
  hTrkEtaBin[0]  -> SetLineStyle(fLinAll);
  hTrkEtaBin[0]  -> SetFillColor(fColAll);
  hTrkEtaBin[0]  -> SetFillStyle(fFilAll);
  hTrkEtaBin[0]  -> SetMarkerColor(fColAll);
  hTrkEtaBin[0]  -> SetMarkerStyle(fMarAll);
  hTrkEtaBin[0]  -> SetTitle("");
  hTrkEtaBin[0]  -> SetTitleFont(fTxt);
  hTrkEtaBin[0]  -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkEtaBin[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkEtaBin[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkEtaBin[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkEtaBin[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkEtaBin[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkEtaBin[0]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkEtaBin[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkEtaBin[0]  -> GetYaxis() -> SetTitle(sTrkEtaY.Data());
  hTrkEtaBin[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkEtaBin[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkEtaBin[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkEtaBin[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkEtaBin[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkEtaBin[0]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkEtaBin[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkEtaBin[1]  -> SetLineColor(fColAll);
  hTrkEtaBin[1]  -> SetLineStyle(fLinAll);
  hTrkEtaBin[1]  -> SetFillColor(fColAll);
  hTrkEtaBin[1]  -> SetFillStyle(fFilAll);
  hTrkEtaBin[1]  -> SetMarkerColor(fColAll);
  hTrkEtaBin[1]  -> SetMarkerStyle(fMarAll);
  hTrkEtaBin[1]  -> SetTitle("");
  hTrkEtaBin[1]  -> SetTitleFont(fTxt);
  hTrkEtaBin[1]  -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkEtaBin[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkEtaBin[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkEtaBin[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkEtaBin[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkEtaBin[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkEtaBin[1]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkEtaBin[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkEtaBin[1]  -> GetYaxis() -> SetTitle(sTrkEtaY.Data());
  hTrkEtaBin[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkEtaBin[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkEtaBin[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkEtaBin[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkEtaBin[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkEtaBin[1]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkEtaBin[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[0] -> SetTitle("");
  hTrkPtVsEta[0] -> SetTitleFont(fTxt);
  hTrkPtVsEta[0] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkPtVsEta[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtVsEta[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtVsEta[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtVsEta[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtVsEta[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtVsEta[0] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsEta[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[0] -> GetYaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtVsEta[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtVsEta[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtVsEta[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtVsEta[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtVsEta[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtVsEta[0] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsEta[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[0] -> GetZaxis() -> SetTitle("");
  hTrkPtVsEta[0] -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkPtVsEta[0] -> GetZaxis() -> SetTitleSize(fTit);
  hTrkPtVsEta[0] -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkPtVsEta[0] -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkPtVsEta[0] -> GetZaxis() -> SetLabelSize(fLab);
  hTrkPtVsEta[0] -> GetZaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[0] -> GetZaxis() -> SetRangeUser(zPtVsEtaPlot[0], zPtVsEtaPlot[1]);
  hTrkPtVsEta[1] -> SetTitle("");
  hTrkPtVsEta[1] -> SetTitleFont(fTxt);
  hTrkPtVsEta[1] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkPtVsEta[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtVsEta[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtVsEta[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtVsEta[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtVsEta[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtVsEta[1] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsEta[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[1] -> GetYaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtVsEta[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtVsEta[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtVsEta[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtVsEta[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtVsEta[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtVsEta[1] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsEta[1] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[1] -> GetZaxis() -> SetTitle("");
  hTrkPtVsEta[1] -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkPtVsEta[1] -> GetZaxis() -> SetTitleSize(fTit);
  hTrkPtVsEta[1] -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkPtVsEta[1] -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkPtVsEta[1] -> GetZaxis() -> SetLabelSize(fLab);
  hTrkPtVsEta[1] -> GetZaxis() -> CenterTitle(fCnt);
  hTrkPtVsEta[1] -> GetZaxis() -> SetRangeUser(zPtVsEtaPlot[0], zPtVsEtaPlot[1]);

  // pT-differential hTrk histograms
  for (UInt_t iTrkBin = 0; iTrkBin < NTrkBins; iTrkBin++) {
    hTrkEtaPtBin[0][iTrkBin] -> SetLineColor(fColTrkPi0[iTrkBin]);
    hTrkEtaPtBin[0][iTrkBin] -> SetLineStyle(fLinAll);
    hTrkEtaPtBin[0][iTrkBin] -> SetFillColor(fColTrkPi0[iTrkBin]);
    hTrkEtaPtBin[0][iTrkBin] -> SetFillStyle(fFilAll);
    hTrkEtaPtBin[0][iTrkBin] -> SetMarkerColor(fColTrkPi0[iTrkBin]);
    hTrkEtaPtBin[0][iTrkBin] -> SetMarkerStyle(fMarTrkPi0[iTrkBin]);
    hTrkEtaPtBin[0][iTrkBin] -> SetTitle("");
    hTrkEtaPtBin[0][iTrkBin] -> SetTitleFont(fTxt);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetTitle(sTrkEta.Data());
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetTitleSize(fTit);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetLabelSize(fLab);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> SetLabelOffset(fOffL);
    hTrkEtaPtBin[0][iTrkBin] -> GetXaxis() -> CenterTitle(fCnt);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetTitle(sTrkEtaY.Data());
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetTitleSize(fTit);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetLabelSize(fLab);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> SetLabelOffset(fOffL);
    hTrkEtaPtBin[0][iTrkBin] -> GetYaxis() -> CenterTitle(fCnt);
    hTrkEtaPtBin[1][iTrkBin] -> SetLineColor(fColTrkGam[iTrkBin]);
    hTrkEtaPtBin[1][iTrkBin] -> SetLineStyle(fLinAll);
    hTrkEtaPtBin[1][iTrkBin] -> SetFillColor(fColTrkGam[iTrkBin]);
    hTrkEtaPtBin[1][iTrkBin] -> SetFillStyle(fFilAll);
    hTrkEtaPtBin[1][iTrkBin] -> SetMarkerColor(fColTrkGam[iTrkBin]);
    hTrkEtaPtBin[1][iTrkBin] -> SetMarkerStyle(fMarTrkGam[iTrkBin]);
    hTrkEtaPtBin[1][iTrkBin] -> SetTitle("");
    hTrkEtaPtBin[1][iTrkBin] -> SetTitleFont(fTxt);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetTitle(sTrkEta.Data());
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetTitleSize(fTit);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetLabelSize(fLab);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> SetLabelOffset(fOffL);
    hTrkEtaPtBin[1][iTrkBin] -> GetXaxis() -> CenterTitle(fCnt);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetTitle(sTrkEtaY.Data());
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetTitleSize(fTit);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetLabelSize(fLab);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> SetLabelOffset(fOffL);
    hTrkEtaPtBin[1][iTrkBin] -> GetYaxis() -> CenterTitle(fCnt);
  }

  // pTtrk
  hTrkPt[0]    -> SetLineColor(fColAll);
  hTrkPt[0]    -> SetLineStyle(fLinAll);
  hTrkPt[0]    -> SetFillColor(fColAll);
  hTrkPt[0]    -> SetFillStyle(fFilAll);
  hTrkPt[0]    -> SetMarkerColor(fColAll);
  hTrkPt[0]    -> SetMarkerStyle(fMarAll);
  hTrkPt[0]    -> SetTitle("");
  hTrkPt[0]    -> SetTitleFont(fTxt);
  hTrkPt[0]    -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPt[0]    -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPt[0]    -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPt[0]    -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPt[0]    -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPt[0]    -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPt[0]    -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPt[0]    -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPt[0]    -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkPt[0]    -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPt[0]    -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPt[0]    -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPt[0]    -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPt[0]    -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPt[0]    -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPt[0]    -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPt[1]    -> SetLineColor(fColCut);
  hTrkPt[1]    -> SetLineStyle(fLinCut);
  hTrkPt[1]    -> SetFillColor(fColCut);
  hTrkPt[1]    -> SetFillStyle(fFilCut);
  hTrkPt[1]    -> SetMarkerColor(fColCut);
  hTrkPt[1]    -> SetMarkerStyle(fMarCut);
  hTrkPt[1]    -> SetTitle("");
  hTrkPt[1]    -> SetTitleFont(fTxt);
  hTrkPt[1]    -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPt[1]    -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPt[1]    -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPt[1]    -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPt[1]    -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPt[1]    -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPt[1]    -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPt[1]    -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPt[1]    -> GetYaxis() -> SetTitle(sCount.Data());
  hTrkPt[1]    -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPt[1]    -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPt[1]    -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPt[1]    -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPt[1]    -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPt[1]    -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPt[1]    -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtBin[0] -> SetLineColor(fColTrg[0]);
  hTrkPtBin[0] -> SetLineStyle(fLinAll);
  hTrkPtBin[0] -> SetFillColor(fColTrg[0]);
  hTrkPtBin[0] -> SetFillStyle(fFilAll);
  hTrkPtBin[0] -> SetMarkerColor(fColTrg[0]);
  hTrkPtBin[0] -> SetMarkerStyle(fMarAll);
  hTrkPtBin[0] -> SetTitle("");
  hTrkPtBin[0] -> SetTitleFont(fTxt);
  hTrkPtBin[0] -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtBin[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtBin[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtBin[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtBin[0] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPtBin[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtBin[0] -> GetYaxis() -> SetTitle(sTrkPtY.Data());
  hTrkPtBin[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtBin[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtBin[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtBin[0] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPtBin[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtBin[1] -> SetLineColor(fColTrg[1]);
  hTrkPtBin[1] -> SetLineStyle(fLinAll);
  hTrkPtBin[1] -> SetFillColor(fColTrg[1]);
  hTrkPtBin[1] -> SetFillStyle(fFilAll);
  hTrkPtBin[1] -> SetMarkerColor(fColTrg[1]);
  hTrkPtBin[1] -> SetMarkerStyle(fMarAll);
  hTrkPtBin[1] -> SetTitle("");
  hTrkPtBin[1] -> SetTitleFont(fTxt);
  hTrkPtBin[1] -> GetXaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtBin[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtBin[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtBin[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtBin[1] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPtBin[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtBin[1] -> GetYaxis() -> SetTitle(sTrkPtY.Data());
  hTrkPtBin[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtBin[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtBin[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtBin[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtBin[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtBin[1] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPtBin[1] -> GetYaxis() -> CenterTitle(fCnt);

  // dFtrk
  hTrkDfBin[0]   -> SetLineColor(fColAll);
  hTrkDfBin[0]   -> SetLineStyle(fLinAll);
  hTrkDfBin[0]   -> SetFillColor(fColAll);
  hTrkDfBin[0]   -> SetFillStyle(fFilAll);
  hTrkDfBin[0]   -> SetMarkerColor(fColAll);
  hTrkDfBin[0]   -> SetMarkerStyle(fMarAll);
  hTrkDfBin[0]   -> SetTitle("");
  hTrkDfBin[0]   -> SetTitleFont(fTxt);
  hTrkDfBin[0]   -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfBin[0]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[0]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfBin[0]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfBin[0]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[0]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfBin[0]   -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDfBin[0]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfBin[0]   -> GetYaxis() -> SetTitle(sTrkDfY.Data());
  hTrkDfBin[0]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[0]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfBin[0]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfBin[0]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[0]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfBin[0]   -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDfBin[0]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDfBin[1]   -> SetLineColor(fColAll);
  hTrkDfBin[1]   -> SetLineStyle(fLinAll);
  hTrkDfBin[1]   -> SetFillColor(fColAll);
  hTrkDfBin[1]   -> SetFillStyle(fFilAll);
  hTrkDfBin[1]   -> SetMarkerColor(fColAll);
  hTrkDfBin[1]   -> SetMarkerStyle(fMarAll);
  hTrkDfBin[1]   -> SetTitle("");
  hTrkDfBin[1]   -> SetTitleFont(fTxt);
  hTrkDfBin[1]   -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfBin[1]   -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[1]   -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfBin[1]   -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfBin[1]   -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[1]   -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfBin[1]   -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDfBin[1]   -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfBin[1]   -> GetYaxis() -> SetTitle(sTrkDfY.Data());
  hTrkDfBin[1]   -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfBin[1]   -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfBin[1]   -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfBin[1]   -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfBin[1]   -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfBin[1]   -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDfBin[1]   -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[0]  -> SetTitle("");
  hTrkPtVsDf[0]  -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[0]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsDf[0]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[0]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsDf[0]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[0]  -> GetZaxis() -> SetTitle("");
  hTrkPtVsDf[0]  -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[0]  -> GetZaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[0]  -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkPtVsDf[0]  -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[0]  -> GetZaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[0]  -> GetZaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[0]  -> GetZaxis() -> SetRangeUser(zPtVsDfPlot[0], zPtVsDfPlot[1]);
  hTrkPtVsDf[1]  -> SetTitle("");
  hTrkPtVsDf[1]  -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitle(sTrkDf.Data());
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[1]  -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsDf[1]  -> GetXaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitle(sTrkPt.Data());
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[1]  -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkPtVsDf[1]  -> GetYaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[1]  -> GetZaxis() -> SetTitle("");
  hTrkPtVsDf[1]  -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkPtVsDf[1]  -> GetZaxis() -> SetTitleSize(fTit);
  hTrkPtVsDf[1]  -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkPtVsDf[1]  -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkPtVsDf[1]  -> GetZaxis() -> SetLabelSize(fLab);
  hTrkPtVsDf[1]  -> GetZaxis() -> CenterTitle(fCnt);
  hTrkPtVsDf[1]  -> GetZaxis() -> SetRangeUser(zPtVsDfPlot[0], zPtVsDfPlot[1]);
  hTrkDfVsEta[0] -> SetTitle("");
  hTrkDfVsEta[0] -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[0] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfVsEta[0] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[0] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[0] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDfVsEta[0] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[0] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfVsEta[0] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[0] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[0] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDfVsEta[0] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[0] -> GetZaxis() -> SetTitle("");
  hTrkDfVsEta[0] -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[0] -> GetZaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[0] -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkDfVsEta[0] -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[0] -> GetZaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[0] -> GetZaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[0] -> GetZaxis() -> SetRangeUser(zEtaVsPhiPlot[0], zEtaVsPhiPlot[1]);
  hTrkDfVsEta[1] -> SetTitle("");
  hTrkDfVsEta[1] -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitle(sTrkEta.Data());
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[1] -> GetXaxis() -> SetTitleOffset(fOffX);
  hTrkDfVsEta[1] -> GetXaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[1] -> GetXaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[1] -> GetXaxis() -> SetLabelOffset(fOffL);
  hTrkDfVsEta[1] -> GetXaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitle(sTrkDf.Data());
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[1] -> GetYaxis() -> SetTitleOffset(fOffY);
  hTrkDfVsEta[1] -> GetYaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[1] -> GetYaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[1] -> GetYaxis() -> SetLabelOffset(fOffL);
  hTrkDfVsEta[1] -> GetYaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[1] -> GetZaxis() -> SetTitle("");
  hTrkDfVsEta[1] -> GetZaxis() -> SetTitleFont(fTxt);
  hTrkDfVsEta[1] -> GetZaxis() -> SetTitleSize(fTit);
  hTrkDfVsEta[1] -> GetZaxis() -> SetTitleOffset(fOffZ);
  hTrkDfVsEta[1] -> GetZaxis() -> SetLabelFont(fTxt);
  hTrkDfVsEta[1] -> GetZaxis() -> SetLabelSize(fLab);
  hTrkDfVsEta[1] -> GetZaxis() -> CenterTitle(fCnt);
  hTrkDfVsEta[1] -> GetZaxis() -> SetRangeUser(zEtaVsPhiPlot[0], zEtaVsPhiPlot[1]);

  // pT-differential dFtrk histograms
  for (UInt_t iTrkBin = 0; iTrkBin < NTrkBins; iTrkBin++) {
    hTrkDfPtBin[0][iTrkBin] -> SetLineColor(fColTrkPi0[iTrkBin]);
    hTrkDfPtBin[0][iTrkBin] -> SetLineStyle(fLinAll);
    hTrkDfPtBin[0][iTrkBin] -> SetFillColor(fColTrkPi0[iTrkBin]);
    hTrkDfPtBin[0][iTrkBin] -> SetFillStyle(fFilAll);
    hTrkDfPtBin[0][iTrkBin] -> SetMarkerColor(fColTrkPi0[iTrkBin]);
    hTrkDfPtBin[0][iTrkBin] -> SetMarkerStyle(fMarTrkPi0[iTrkBin]);
    hTrkDfPtBin[0][iTrkBin] -> SetTitle("");
    hTrkDfPtBin[0][iTrkBin] -> SetTitleFont(fTxt);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetTitle(sTrkDf.Data());
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetTitleSize(fTit);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetLabelSize(fLab);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> SetLabelOffset(fOffL);
    hTrkDfPtBin[0][iTrkBin] -> GetXaxis() -> CenterTitle(fCnt);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetTitle(sTrkDfY.Data());
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetTitleSize(fTit);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetLabelSize(fLab);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> SetLabelOffset(fOffL);
    hTrkDfPtBin[0][iTrkBin] -> GetYaxis() -> CenterTitle(fCnt);
    hTrkDfPtBin[1][iTrkBin] -> SetLineColor(fColTrkGam[iTrkBin]);
    hTrkDfPtBin[1][iTrkBin] -> SetLineStyle(fLinAll);
    hTrkDfPtBin[1][iTrkBin] -> SetFillColor(fColTrkGam[iTrkBin]);
    hTrkDfPtBin[1][iTrkBin] -> SetFillStyle(fFilAll);
    hTrkDfPtBin[1][iTrkBin] -> SetMarkerColor(fColTrkGam[iTrkBin]);
    hTrkDfPtBin[1][iTrkBin] -> SetMarkerStyle(fMarTrkGam[iTrkBin]);
    hTrkDfPtBin[1][iTrkBin] -> SetTitle("");
    hTrkDfPtBin[1][iTrkBin] -> SetTitleFont(fTxt);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetTitle(sTrkDf.Data());
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetTitleSize(fTit);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetLabelSize(fLab);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> SetLabelOffset(fOffL);
    hTrkDfPtBin[1][iTrkBin] -> GetXaxis() -> CenterTitle(fCnt);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetTitle(sTrkDfY.Data());
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetTitleSize(fTit);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetLabelSize(fLab);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> SetLabelOffset(fOffL);
    hTrkDfPtBin[1][iTrkBin] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "      Styles set." << endl;

  // make legends
  const UInt_t  fColLeg(0);
  const UInt_t  fFilLeg(0);
  const UInt_t  fLinLeg(0);
  const UInt_t  fAlnLeg(12);
  const UInt_t  fAlnTxt(32);
  const Float_t xyLegSmall[4] = {0.1, 0.1, 0.3, 0.2};
  const Float_t xyLegMed[4]   = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyLegBig[4]   = {0.1, 0.1, 0.3, 0.35};
  const Float_t xyLegAlt[4]   = {0.5, 0.1, 0.7, 0.3};
  const Float_t xyTxtSmall[4] = {0.3, 0.1, 0.5, 0.2};
  const Float_t xyTxtBig[4]   = {0.3, 0.1, 0.5, 0.25};

  TLegend *lEvtVz      = new TLegend(xyLegSmall[0], xyLegSmall[1], xyLegSmall[2], xyLegSmall[3]);
  TLegend *lTrgEtPi0   = new TLegend(xyLegMed[0], xyLegMed[1], xyLegMed[2], xyLegMed[3]);
  TLegend *lTrgEtAllG  = new TLegend(xyLegAlt[0], xyLegAlt[1], xyLegAlt[2], xyLegAlt[3]);
  TLegend *lTrkDfVsPtP = new TLegend(xyLegBig[0], xyLegBig[1], xyLegBig[2], xyLegBig[3]);
  lEvtVz      -> SetFillColor(fColLeg);
  lEvtVz      -> SetFillStyle(fFilLeg);
  lEvtVz      -> SetLineColor(fColLeg);
  lEvtVz      -> SetLineStyle(fLinLeg);
  lEvtVz      -> SetTextFont(fTxt);
  lEvtVz      -> SetTextAlign(fAlnLeg);
  lTrgEtPi0   -> SetFillColor(fColLeg);
  lTrgEtPi0   -> SetFillStyle(fFilLeg);
  lTrgEtPi0   -> SetLineColor(fColLeg);
  lTrgEtPi0   -> SetLineStyle(fLinLeg);
  lTrgEtPi0   -> SetTextFont(fTxt);
  lTrgEtPi0   -> SetTextAlign(fAlnLeg);
  lTrgEtAllG  -> SetFillColor(fColLeg);
  lTrgEtAllG  -> SetFillStyle(fFilLeg);
  lTrgEtAllG  -> SetLineColor(fColLeg);
  lTrgEtAllG  -> SetLineStyle(fLinLeg);
  lTrgEtAllG  -> SetTextFont(fTxt);
  lTrgEtAllG  -> SetTextAlign(fAlnLeg);
  lTrkDfVsPtP -> SetFillColor(fColLeg);
  lTrkDfVsPtP -> SetFillStyle(fFilLeg);
  lTrkDfVsPtP -> SetLineColor(fColLeg);
  lTrkDfVsPtP -> SetLineStyle(fLinLeg);
  lTrkDfVsPtP -> SetTextFont(fTxt);
  lTrkDfVsPtP -> SetTextAlign(fAlnLeg);

  // small legends
  TLegend *lTrgEt    = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrgEtSum = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrkNfit  = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrkRfit  = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrkDca   = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrkEta   = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrkPt    = (TLegend*) lEvtVz -> Clone();
  TLegend *lTrkPtBin = (TLegend*) lEvtVz -> Clone();

  // medium legends
  TLegend *lTrgEtGam  = (TLegend*) lTrgEtPi0 -> Clone();
  TLegend *lTrgEtAllP = (TLegend*) lTrgEtPi0 -> Clone();
  TLegend *lTrgTsp    = (TLegend*) lTrgEtPi0 -> Clone();

  // big legends
  TLegend *lTrkDfVsPtG  = (TLegend*) lTrkDfVsPtP -> Clone();
  TLegend *lTrkEtaVsPtP = (TLegend*) lTrkDfVsPtP -> Clone();
  TLegend *lTrkEtaVsPtG = (TLegend*) lTrkDfVsPtP -> Clone();

  // vertices
  lEvtVz -> AddEntry(hEvtVz[0], "all events", "pf");
  lEvtVz -> AddEntry(hEvtVz[1], "excluded events", "f");

  // eTtrg
  lTrgEt     -> AddEntry(hTrgEt[0], "all triggers", "pf");
  lTrgEt     -> AddEntry(hTrgEt[1], "excluded triggers", "f");
  lTrgEtPi0  -> AddEntry(hTrgEtBin[1][0], "E_{T}^{trg} #in (9, 11) GeV", "pf");
  lTrgEtPi0  -> AddEntry(hTrgEtBin[2][0], "E_{T}^{trg} #in (11, 15) GeV", "pf");
  lTrgEtPi0  -> AddEntry(hTrgEtBin[3][0], "E_{T}^{trg} #in (15, 20) GeV", "pf");
  lTrgEtGam  -> AddEntry(hTrgEtBin[1][1], "E_{T}^{trg} #in (9, 11) GeV", "pf");
  lTrgEtGam  -> AddEntry(hTrgEtBin[2][1], "E_{T}^{trg} #in (11, 15) GeV", "pf");
  lTrgEtGam  -> AddEntry(hTrgEtBin[3][1], "E_{T}^{trg} #in (15, 20) GeV", "pf");
  lTrgEtSum  -> AddEntry(hTrgEtBin[0][0], "#pi^{0} triggers", "pf");
  lTrgEtSum  -> AddEntry(hTrgEtBin[0][1], "#gamma_{rich} triggers", "pf");
  lTrgEtSum  -> AddEntry(hTrgEtSum, "both triggers", "pf");
  lTrgEtAllP -> AddEntry(hTrgEtSum, "all accepted triggers", "pf");
  lTrgEtAllP -> AddEntry(hTrgEtBin[1][0], "9 - 11 GeV #pi^{0}", "pf");
  lTrgEtAllP -> AddEntry(hTrgEtBin[2][0], "11 - 15 GeV #pi^{0}", "pf");
  lTrgEtAllP -> AddEntry(hTrgEtBin[3][0], "15 - 20 GeV #pi^{0}", "pf");
  lTrgEtAllG -> AddEntry((TObject*) 0, "", "");
  lTrgEtAllG -> AddEntry(hTrgEtBin[1][1], "9 - 11 GeV #gamma_{rich}", "pf");
  lTrgEtAllG -> AddEntry(hTrgEtBin[2][1], "11 - 15 GeV #gamma_{rich}", "pf");
  lTrgEtAllG -> AddEntry(hTrgEtBin[3][1], "15 - 20 GeV #gamma_{rich}", "pf");

  // tsp
  lTrgTsp -> AddEntry(hTrgTsp[2], "all accepted triggers", "pf");
  lTrgTsp -> AddEntry(hTrgTsp[1], "#gamma_{rich} triggers", "f");
  lTrgTsp -> AddEntry(hTrgTsp[0], "#pi^{0} triggers", "f");

  // nFit and nFit/nPoss
  lTrkNfit -> AddEntry(hTrkNfit[0], "all tracks", "pf");
  lTrkNfit -> AddEntry(hTrkNfit[1], "excluded tracks", "f");
  lTrkRfit -> AddEntry(hTrkRfit[0], "all tracks", "pf");
  lTrkRfit -> AddEntry(hTrkRfit[1], "excluded tracks", "f");

  // dca
  lTrkDca -> AddEntry(hTrkDca[0], "all tracks", "pf");
  lTrkDca -> AddEntry(hTrkDca[1], "excluded tracks", "f");

  // eta
  lTrkEta -> AddEntry(hTrkEta[0], "all tracks", "pf");
  lTrkEta -> AddEntry(hTrkEta[1], "excluded tracks", "f");

  // pT
  lTrkPt    -> AddEntry(hTrkPt[0], "all tracks", "pf");
  lTrkPt    -> AddEntry(hTrkPt[1], "excluded tracks", "f");
  lTrkPtBin -> AddEntry(hTrkPtBin[0], "#pi^{0} trigger", "pf");
  lTrkPtBin -> AddEntry(hTrkPtBin[1], "#gamma_{rich} trigger", "pf");

  // dF vs. pT
  lTrkDfVsPtP -> AddEntry(hTrkDfBin[0], "all accepted tracks", "pf");
  lTrkDfVsPtP -> AddEntry(hTrkDfPtBin[0][0], "p_{T}^{trk} #in (0.2, 2) GeV/c", "pf");
  lTrkDfVsPtP -> AddEntry(hTrkDfPtBin[0][1], "p_{T}^{trk} #in (2, 5) GeV/c", "pf");
  lTrkDfVsPtP -> AddEntry(hTrkDfPtBin[0][2], "p_{T}^{trk} #in (5, 9) GeV/c", "pf");
  lTrkDfVsPtP -> AddEntry(hTrkDfPtBin[0][3], "p_{T}^{trk} #in (9, 20) GeV/c", "pf");
  lTrkDfVsPtG -> AddEntry(hTrkDfBin[0], "all accepted tracks", "pf");
  lTrkDfVsPtG -> AddEntry(hTrkDfPtBin[1][0], "p_{T}^{trk} #in (0.2, 2) GeV/c", "pf");
  lTrkDfVsPtG -> AddEntry(hTrkDfPtBin[1][1], "p_{T}^{trk} #in (2, 5) GeV/c", "pf");
  lTrkDfVsPtG -> AddEntry(hTrkDfPtBin[1][2], "p_{T}^{trk} #in (5, 9) GeV/c", "pf");
  lTrkDfVsPtG -> AddEntry(hTrkDfPtBin[1][3], "p_{T}^{trk} #in (9, 20) GeV/c", "pf");

  // eta vs. pT
  lTrkEtaVsPtP -> AddEntry(hTrkEtaBin[0], "all accepted tracks", "pf");
  lTrkEtaVsPtP -> AddEntry(hTrkEtaPtBin[0][0], "p_{T}^{trk} #in (0.2, 2) GeV/c", "pf");
  lTrkEtaVsPtP -> AddEntry(hTrkEtaPtBin[0][1], "p_{T}^{trk} #in (2, 5) GeV/c", "pf");
  lTrkEtaVsPtP -> AddEntry(hTrkEtaPtBin[0][2], "p_{T}^{trk} #in (5, 9) GeV/c", "pf");
  lTrkEtaVsPtP -> AddEntry(hTrkEtaPtBin[0][3], "p_{T}^{trk} #in (9, 20) GeV/c", "pf");
  lTrkEtaVsPtG -> AddEntry(hTrkEtaBin[0], "all accepted tracks", "pf");
  lTrkEtaVsPtG -> AddEntry(hTrkEtaPtBin[1][0], "p_{T}^{trk} #in (0.2, 2) GeV/c", "pf");
  lTrkEtaVsPtG -> AddEntry(hTrkEtaPtBin[1][1], "p_{T}^{trk} #in (2, 5) GeV/c", "pf");
  lTrkEtaVsPtG -> AddEntry(hTrkEtaPtBin[1][2], "p_{T}^{trk} #in (5, 9) GeV/c", "pf");
  lTrkEtaVsPtG -> AddEntry(hTrkEtaPtBin[1][3], "p_{T}^{trk} #in (9, 20) GeV/c", "pf");
  cout << "      Made legends." << endl;

  // make labels
  TPaveText *pEvt = new TPaveText(xyTxtBig[0], xyTxtBig[1], xyTxtBig[2], xyTxtBig[3], "NDC NB");
  TPaveText *pPi0 = new TPaveText(xyTxtSmall[0], xyTxtSmall[1], xyTxtSmall[2], xyTxtSmall[3], "NDC NB");
  pEvt -> SetFillColor(fColLeg);
  pEvt -> SetFillStyle(fFilLeg);
  pEvt -> SetLineColor(fColLeg);
  pEvt -> SetLineStyle(fLinLeg);
  pEvt -> SetTextFont(fTxt);
  pEvt -> SetTextAlign(fAlnTxt);
  pPi0 -> SetFillColor(fColLeg);
  pPi0 -> SetFillStyle(fFilLeg);
  pPi0 -> SetLineColor(fColLeg);
  pPi0 -> SetLineStyle(fLinLeg);
  pPi0 -> SetTextFont(fTxt);
  pPi0 -> SetTextAlign(fAlnTxt);

  // big labels
  TPaveText *pTrg = (TPaveText*) pEvt -> Clone();

  // small labels
  TPaveText *pGam = (TPaveText*) pPi0 -> Clone();
  TPaveText *pSum = (TPaveText*) pPi0 -> Clone();

  pEvt -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
  pEvt -> AddText("Run9 data, L2gamma stream");
  pEvt -> AddText("Bad runs removed");
  pTrg -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
  pTrg -> AddText("Run9 data, L2gamma stream");
  pTrg -> AddText("Bad runs and towers removed");
  pPi0 -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
  pPi0 -> AddText("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  pGam -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
  pGam -> AddText("#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");
  pSum -> AddText("pp-collisions, #sqrt{s} = 200 GeV");
  pSum -> AddText("#pi^{0}/#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");
  cout << "      Made text boxes." << endl;

  // make trigger and track selection boxes
  const UInt_t  fColSel(923);
  const UInt_t  fFilSel(3354);
  const UInt_t  fLinSel(1);
  const Float_t eTvsEtaSel[4]  = {-0.9, 9., 0.9, 20.};
  const Float_t eTvsPhiSel[4]  = {-3.14, 9., 3.14, 20.};
  const Float_t pTvsDcaSel[4]  = {0.2, 0., 30., 1.};
  const Float_t pTvsNfitSel[4] = {0.2, 15., 30., 47.};
  const Float_t pTvsDfSel[4]   = {-1.6, 0.2, 4.7, 30.};
  const Float_t pTvsEtaSel[4]  = {-1., 0.2, 1., 30.};

  TBox *bEtVsEtaSel  = new TBox(eTvsEtaSel[0], eTvsEtaSel[1], eTvsEtaSel[2], eTvsEtaSel[3]);
  TBox *bEtVsPhiSel  = new TBox(eTvsPhiSel[0], eTvsPhiSel[1], eTvsPhiSel[2], eTvsPhiSel[3]);
  TBox *bPtVsDcaSel  = new TBox(pTvsDcaSel[0], pTvsDcaSel[1], pTvsDcaSel[2], pTvsDcaSel[3]);
  TBox *bPtVsNfitSel = new TBox(pTvsNfitSel[0], pTvsNfitSel[1], pTvsNfitSel[2], pTvsNfitSel[3]);
  TBox *bPtVsDfSel   = new TBox(pTvsDfSel[0], pTvsDfSel[1], pTvsDfSel[2], pTvsDfSel[3]);
  TBox *bPtVsEtaSel  = new TBox(pTvsEtaSel[0], pTvsEtaSel[1], pTvsEtaSel[2], pTvsEtaSel[3]);
  bEtVsEtaSel  -> SetFillColor(fColSel);
  bEtVsEtaSel  -> SetFillStyle(fFilSel);
  bEtVsEtaSel  -> SetLineColor(fColSel);
  bEtVsEtaSel  -> SetLineStyle(fLinSel);
  bEtVsPhiSel  -> SetFillColor(fColSel);
  bEtVsPhiSel  -> SetFillStyle(fFilSel);
  bEtVsPhiSel  -> SetLineColor(fColSel);
  bEtVsPhiSel  -> SetLineStyle(fLinSel);
  bPtVsDcaSel  -> SetFillColor(fColSel);
  bPtVsDcaSel  -> SetFillStyle(fFilSel);
  bPtVsDcaSel  -> SetLineColor(fColSel);
  bPtVsDcaSel  -> SetLineStyle(fLinSel);
  bPtVsNfitSel -> SetFillColor(fColSel);
  bPtVsNfitSel -> SetFillStyle(fFilSel);
  bPtVsNfitSel -> SetLineColor(fColSel);
  bPtVsNfitSel -> SetLineStyle(fLinSel);
  bPtVsDfSel   -> SetFillColor(fColSel);
  bPtVsDfSel   -> SetFillStyle(fFilSel);
  bPtVsDfSel   -> SetLineColor(fColSel);
  bPtVsDfSel   -> SetLineStyle(fLinSel);
  bPtVsEtaSel  -> SetFillColor(fColSel);
  bPtVsEtaSel  -> SetFillStyle(fFilSel);
  bPtVsEtaSel  -> SetLineColor(fColSel);
  bPtVsEtaSel  -> SetLineStyle(fLinSel);
  cout << "      Made 2d selection boxes." << endl;

  // make cut (and TSP) histograms
  const UInt_t tspFillPi0(3345);
  const UInt_t tspFillGam(3354);
  const UInt_t tspFillCut(3345);

  TH1D *hTrgTspPi0     = (TH1D*) hTrgTsp[0]     -> Clone();
  TH1D *hTrgTspGam     = (TH1D*) hTrgTsp[1]     -> Clone();
  TH1D *hTrkNfitCut    = (TH1D*) hTrkNfit[1]    -> Clone();
  TH1D *hTrkRfitCut    = (TH1D*) hTrkRfit[1]    -> Clone();
  TH1D *hTrkDcaCut     = (TH1D*) hTrkDca[1]     -> Clone();
  TH1D *hTrkEtaCut     = (TH1D*) hTrkEta[1]     -> Clone();
  TH1D *hTrkPtCut      = (TH1D*) hTrkPt[1]      -> Clone();
  hTrgTspPi0  -> SetNameTitle("hTrgTspPi0", "");
  hTrgTspPi0  -> SetFillStyle(tspFillPi0);
  hTrgTspGam  -> SetNameTitle("hTrgTspGam", "");
  hTrgTspGam  -> SetFillStyle(tspFillGam);
  hTrkNfitCut -> SetNameTitle("hTrkNfitCut", "");
  hTrkNfitCut -> SetFillStyle(tspFillCut);
  hTrkRfitCut -> SetNameTitle("hTrkRfitCut", "");
  hTrkRfitCut -> SetFillStyle(tspFillCut);
  hTrkDcaCut  -> SetNameTitle("hTrkDcaCut", "");
  hTrkDcaCut  -> SetFillStyle(tspFillCut);
  hTrkEtaCut  -> SetNameTitle("hTrkEtaCut", "");
  hTrkEtaCut  -> SetFillStyle(tspFillCut);
  hTrkPtCut   -> SetNameTitle("hTrkPtCut", "");
  hTrkPtCut   -> SetFillStyle(tspFillCut);

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

  // evt no's
  TCanvas *cEvtNum = new TCanvas("cEvtNum", "", width, height);
  cEvtNum -> SetLogx(fLogX);
  cEvtNum -> SetLogy(fLogY);
  cEvtNum -> SetGrid(fGrid, fGrid);
  cEvtNum -> SetTicks(fTick, fTick);
  cEvtNum -> SetBorderMode(fMode);
  cEvtNum -> SetBorderSize(fBord);
  cEvtNum -> SetFrameBorderMode(fFrame);
  cEvtNum -> SetLeftMargin(fMarginBig);
  cEvtNum -> SetTopMargin(fMarginSmall);
  cEvtNum -> SetRightMargin(fMarginSmall);
  cEvtNum -> SetBottomMargin(fMarginBig);
  cEvtNum -> cd();
  hEvtNum -> Draw("B");
  fOutput -> cd();
  cEvtNum -> Write();
  cEvtNum -> Close();

  // vertices
  TCanvas *cEvtVz = new TCanvas("cEvtVz", "", width, height);
  cEvtVz    -> SetLogx(fLogX);
  cEvtVz    -> SetLogy(fLogY);
  cEvtVz    -> SetGrid(fGrid, fGrid);
  cEvtVz    -> SetTicks(fTick, fTick);
  cEvtVz    -> SetBorderMode(fMode);
  cEvtVz    -> SetBorderSize(fBord);
  cEvtVz    -> SetFrameBorderMode(fFrame);
  cEvtVz    -> SetLeftMargin(fMarginBig);
  cEvtVz    -> SetTopMargin(fMarginSmall);
  cEvtVz    -> SetRightMargin(fMarginSmall);
  cEvtVz    -> SetBottomMargin(fMarginBig);
  cEvtVz    -> cd();
  hEvtVz[0] -> Draw("");
  hEvtVz[1] -> Draw("same");
  hEvtVz[1] -> Draw("hist same");
  lEvtVz    -> Draw();
  pEvt      -> Draw();
  fOutput   -> cd();
  cEvtVz    -> Write();
  cEvtVz    -> Close();

  TCanvas *cEvtVzVsVxVy = new TCanvas("cEvtVzVsVxVy", "", bigWidth, height);
  TPad    *pEvtVzVsVx   = new TPad("pEvtVzVsVx", "", 0., 0., 0.5, 1.);
  TPad    *pEvtVzVsVy   = new TPad("pEvtVzVsVy", "", 0.5, 0., 1., 1.);
  pEvtVzVsVx   -> SetLogx(fLogX);
  pEvtVzVsVx   -> SetLogy(fLogY2);
  pEvtVzVsVx   -> SetLogz(fLogZ);
  pEvtVzVsVx   -> SetGrid(fGrid, fGrid);
  pEvtVzVsVx   -> SetTicks(fTick, fTick);
  pEvtVzVsVx   -> SetBorderMode(fMode);
  pEvtVzVsVx   -> SetBorderSize(fBord);
  pEvtVzVsVx   -> SetFrameBorderMode(fFrame);
  pEvtVzVsVx   -> SetLeftMargin(fMarginBig);
  pEvtVzVsVx   -> SetTopMargin(fMarginSmall);
  pEvtVzVsVx   -> SetRightMargin(fMarginBord);
  pEvtVzVsVx   -> SetBottomMargin(fMarginBig);
  pEvtVzVsVy   -> SetLogx(fLogX);
  pEvtVzVsVy   -> SetLogy(fLogY2);
  pEvtVzVsVy   -> SetLogz(fLogZ);
  pEvtVzVsVy   -> SetGrid(fGrid, fGrid);
  pEvtVzVsVy   -> SetTicks(fTick, fTick);
  pEvtVzVsVy   -> SetBorderMode(fMode);
  pEvtVzVsVy   -> SetBorderSize(fBord);
  pEvtVzVsVy   -> SetFrameBorderMode(fFrame);
  pEvtVzVsVy   -> SetLeftMargin(fMarginBord);
  pEvtVzVsVy   -> SetTopMargin(fMarginSmall);
  pEvtVzVsVy   -> SetRightMargin(fMarginBig);
  pEvtVzVsVy   -> SetBottomMargin(fMarginBig);
  cEvtVzVsVxVy -> cd();
  pEvtVzVsVx   -> Draw();
  pEvtVzVsVy   -> Draw();
  pEvtVzVsVx   -> cd();
  hEvtVzVsVx   -> Draw("col");
  pEvtVzVsVy   -> cd();
  hEvtVzVsVy   -> Draw("colz");
  pEvt         -> Draw();
  fOutput      -> cd();
  cEvtVzVsVxVy -> Write();
  cEvtVzVsVxVy -> Close();

  TCanvas *cEvtVyVsVx = new TCanvas("cEvtVyVsVx", "", width, height);
  cEvtVyVsVx -> SetLogx(fLogX);
  cEvtVyVsVx -> SetLogy(fLogY2);
  cEvtVyVsVx -> SetLogz(fLogZ);
  cEvtVyVsVx -> SetGrid(fGrid, fGrid);
  cEvtVyVsVx -> SetTicks(fTick, fTick);
  cEvtVyVsVx -> SetBorderMode(fMode);
  cEvtVyVsVx -> SetBorderSize(fBord);
  cEvtVyVsVx -> SetFrameBorderMode(fFrame);
  cEvtVyVsVx -> SetLeftMargin(fMarginBig);
  cEvtVyVsVx -> SetTopMargin(fMarginSmall);
  cEvtVyVsVx -> SetRightMargin(fMarginBig);
  cEvtVyVsVx -> SetBottomMargin(fMarginBig);
  cEvtVyVsVx -> cd();
  hEvtVyVsVx -> Draw("colz");
  pEvt       -> Draw();
  fOutput    -> cd();
  cEvtVyVsVx -> Write();
  cEvtVyVsVx -> Close();

  // no. of tracks
  TCanvas *cEvtPrim = new TCanvas("cEvtPrim", "", width, height);
  cEvtPrim -> SetLogx(fLogX);
  cEvtPrim -> SetLogy(fLogY);
  cEvtPrim -> SetGrid(fGrid, fGrid);
  cEvtPrim -> SetTicks(fTick, fTick);
  cEvtPrim -> SetBorderMode(fMode);
  cEvtPrim -> SetBorderSize(fBord);
  cEvtPrim -> SetFrameBorderMode(fFrame);
  cEvtPrim -> SetLeftMargin(fMarginBig);
  cEvtPrim -> SetTopMargin(fMarginSmall);
  cEvtPrim -> SetRightMargin(fMarginSmall);
  cEvtPrim -> SetBottomMargin(fMarginBig);
  cEvtPrim -> cd();
  hEvtPrim -> Draw();
  pEvt     -> Draw();
  fOutput  -> cd();
  cEvtPrim -> Write();
  cEvtPrim -> Close();

  // eTtrg (no tsp cuts)
  TCanvas *cTrgEt = new TCanvas("cTrgEt", "", width, height);
  cTrgEt    -> SetLogx(fLogX);
  cTrgEt    -> SetLogy(fLogY);
  cTrgEt    -> SetGrid(fGrid, fGrid);
  cTrgEt    -> SetTicks(fTick, fTick);
  cTrgEt    -> SetBorderMode(fMode);
  cTrgEt    -> SetBorderSize(fBord);
  cTrgEt    -> SetFrameBorderMode(fFrame);
  cTrgEt    -> SetLeftMargin(fMarginBig);
  cTrgEt    -> SetTopMargin(fMarginSmall);
  cTrgEt    -> SetRightMargin(fMarginSmall);
  cTrgEt    -> SetBottomMargin(fMarginBig);
  cTrgEt    -> cd();
  hTrgEt[0] -> Draw();
  hTrgEt[1] -> Draw("same");
  hTrgEt[1] -> Draw("hist same");
  lTrgEt    -> Draw();
  pTrg      -> Draw();
  fOutput   -> cd();
  cTrgEt    -> Write();
  cTrgEt    -> Close();

  TCanvas *cTrgEtVsEtaPhi = new TCanvas("cTrgEtVsEtaPhi", "", bigWidth, height);
  TPad    *pTrgEtVsEta    = new TPad("pTrgEtVsEta", "", 0., 0., 0.5, 1.);
  TPad    *pTrgEtVsPhi    = new TPad("pTrgEtVsPhi", "", 0.5, 0., 1., 1.);
  pTrgEtVsEta    -> SetLogx(fLogX);
  pTrgEtVsEta    -> SetLogy(fLogY2);
  pTrgEtVsEta    -> SetLogz(fLogZ);
  pTrgEtVsEta    -> SetGrid(fGrid, fGrid);
  pTrgEtVsEta    -> SetTicks(fTick, fTick);
  pTrgEtVsEta    -> SetBorderMode(fMode);
  pTrgEtVsEta    -> SetBorderSize(fBord);
  pTrgEtVsEta    -> SetFrameBorderMode(fFrame);
  pTrgEtVsEta    -> SetLeftMargin(fMarginBig);
  pTrgEtVsEta    -> SetTopMargin(fMarginSmall);
  pTrgEtVsEta    -> SetRightMargin(fMarginBord);
  pTrgEtVsEta    -> SetBottomMargin(fMarginBig);
  pTrgEtVsPhi    -> SetLogx(fLogX);
  pTrgEtVsPhi    -> SetLogy(fLogY2);
  pTrgEtVsPhi    -> SetLogz(fLogZ);
  pTrgEtVsPhi    -> SetGrid(fGrid, fGrid);
  pTrgEtVsPhi    -> SetTicks(fTick, fTick);
  pTrgEtVsPhi    -> SetBorderMode(fMode);
  pTrgEtVsPhi    -> SetBorderSize(fBord);
  pTrgEtVsPhi    -> SetFrameBorderMode(fFrame);
  pTrgEtVsPhi    -> SetLeftMargin(fMarginBord);
  pTrgEtVsPhi    -> SetTopMargin(fMarginSmall);
  pTrgEtVsPhi    -> SetRightMargin(fMarginBig);
  pTrgEtVsPhi    -> SetBottomMargin(fMarginBig);
  cTrgEtVsEtaPhi -> cd();
  pTrgEtVsEta    -> Draw();
  pTrgEtVsPhi    -> Draw();
  pTrgEtVsEta    -> cd();
  hTrgEtVsEta    -> Draw("col");
  bEtVsEtaSel    -> Draw();
  pTrgEtVsPhi    -> cd();
  hTrgEtVsPhi    -> Draw("colz");
  bEtVsPhiSel    -> Draw();
  pSum           -> Draw();
  fOutput        -> cd();
  cTrgEtVsEtaPhi -> Write();
  cTrgEtVsEtaPhi -> Close();

  // eTtrg (tsp cuts)
  TCanvas *cTrgEtBin = new TCanvas("cTrgEtBin", "", bigWidth, height);
  TPad    *pTrgEtPi0 = new TPad("pTrgEtPi0", "", 0., 0., 0.5, 1.);
  TPad    *pTrgEtGam = new TPad("pTrgEtGam", "", 0.5, 0., 1., 1.);
  pTrgEtPi0       -> SetLogx(fLogX);
  pTrgEtPi0       -> SetLogy(fLogY);
  pTrgEtPi0       -> SetGrid(fGrid, fGrid);
  pTrgEtPi0       -> SetTicks(fTick, fTick);
  pTrgEtPi0       -> SetBorderMode(fMode);
  pTrgEtPi0       -> SetBorderSize(fBord);
  pTrgEtPi0       -> SetFrameBorderMode(fFrame);
  pTrgEtPi0       -> SetLeftMargin(fMarginBig);
  pTrgEtPi0       -> SetTopMargin(fMarginSmall);
  pTrgEtPi0       -> SetRightMargin(fMarginBord);
  pTrgEtPi0       -> SetBottomMargin(fMarginBig);
  pTrgEtGam       -> SetLogx(fLogX);
  pTrgEtGam       -> SetLogy(fLogY);
  pTrgEtGam       -> SetGrid(fGrid, fGrid);
  pTrgEtGam       -> SetTicks(fTick, fTick);
  pTrgEtGam       -> SetBorderMode(fMode);
  pTrgEtGam       -> SetBorderSize(fBord);
  pTrgEtGam       -> SetFrameBorderMode(fFrame);
  pTrgEtGam       -> SetLeftMargin(fMarginBord);
  pTrgEtGam       -> SetTopMargin(fMarginSmall);
  pTrgEtGam       -> SetRightMargin(fMarginSmall);
  pTrgEtGam       -> SetBottomMargin(fMarginBig);
  cTrgEtBin       -> cd();
  pTrgEtPi0       -> Draw();
  pTrgEtGam       -> Draw();
  pTrgEtPi0       -> cd();
  hTrgEtBin[1][0] -> Draw();
  hTrgEtBin[2][0] -> Draw("same");
  hTrgEtBin[3][0] -> Draw("same");
  lTrgEtPi0       -> Draw();
  pPi0            -> Draw();
  pTrgEtGam       -> cd();
  hTrgEtBin[1][1] -> Draw();
  hTrgEtBin[2][1] -> Draw("same");
  hTrgEtBin[3][1] -> Draw("same");
  lTrgEtGam       -> Draw();
  pGam            -> Draw();
  fOutput         -> cd();
  cTrgEtBin       -> Write();
  cTrgEtBin       -> Close();

  TCanvas *cTrgEtSum = new TCanvas("cTrgEtSum", "", width, height);
  cTrgEtSum       -> SetLogx(fLogX);
  cTrgEtSum       -> SetLogy(fLogY);
  cTrgEtSum       -> SetGrid(fGrid, fGrid);
  cTrgEtSum       -> SetTicks(fTick, fTick);
  cTrgEtSum       -> SetBorderMode(fMode);
  cTrgEtSum       -> SetBorderSize(fBord);
  cTrgEtSum       -> SetFrameBorderMode(fFrame);
  cTrgEtSum       -> SetLeftMargin(fMarginBig);
  cTrgEtSum       -> SetTopMargin(fMarginSmall);
  cTrgEtSum       -> SetRightMargin(fMarginSmall);
  cTrgEtSum       -> SetBottomMargin(fMarginBig);
  hTrgEtBin[0][0] -> Draw();
  hTrgEtBin[0][1] -> Draw("same");
  hTrgEtSum       -> Draw("same");
  lTrgEtSum       -> Draw();
  pSum            -> Draw();
  fOutput         -> cd();
  cTrgEtSum       -> Write();
  cTrgEtSum       -> Close(); 

  TCanvas *cTrgEtAll = new TCanvas("cTrgEtAll", "", width, height);
  cTrgEtAll       -> SetLogx(fLogX);
  cTrgEtAll       -> SetLogy(fLogY);
  cTrgEtAll       -> SetGrid(fGrid, fGrid);
  cTrgEtAll       -> SetTicks(fTick, fTick);
  cTrgEtAll       -> SetBorderMode(fMode);
  cTrgEtAll       -> SetBorderSize(fBord);
  cTrgEtAll       -> SetFrameBorderMode(fFrame);
  cTrgEtAll       -> SetLeftMargin(fMarginBig);
  cTrgEtAll       -> SetTopMargin(fMarginSmall);
  cTrgEtAll       -> SetRightMargin(fMarginSmall);
  cTrgEtAll       -> SetBottomMargin(fMarginBig);
  cTrgEtAll       -> cd();
  hTrgEtSum       -> Draw();
  hTrgEtBin[1][0] -> Draw("same");
  hTrgEtBin[2][0] -> Draw("same");
  hTrgEtBin[3][0] -> Draw("same");
  hTrgEtBin[1][1] -> Draw("same");
  hTrgEtBin[2][1] -> Draw("same");
  hTrgEtBin[3][1] -> Draw("same");
  lTrgEtAllP      -> Draw();
  lTrgEtAllG      -> Draw();
  pSum            -> Draw();
  fOutput         -> cd();
  cTrgEtAll       -> Write();
  cTrgEtAll       -> Close();

  // tsp
  TCanvas *cTrgTsp = new TCanvas("cTrgTsp", "", width, height);
  cTrgTsp    -> SetLogx(fLogX);
  cTrgTsp    -> SetLogy(fLogY2);
  cTrgTsp    -> SetGrid(fGrid, fGrid);
  cTrgTsp    -> SetTicks(fTick, fTick);
  cTrgTsp    -> SetBorderMode(fMode);
  cTrgTsp    -> SetBorderSize(fBord);
  cTrgTsp    -> SetFrameBorderMode(fFrame);
  cTrgTsp    -> SetLeftMargin(fMarginBig);
  cTrgTsp    -> SetTopMargin(fMarginSmall);
  cTrgTsp    -> SetRightMargin(fMarginSmall);
  cTrgTsp    -> SetBottomMargin(fMarginBig);
  cTrgTsp    -> cd();
  hTrgTsp[2] -> Draw();
  hTrgTsp[1] -> Draw("same");
  hTrgTspGam -> Draw("hist same");
  hTrgTsp[0] -> Draw("same");
  hTrgTspPi0 -> Draw("hist same");
  lTrgTsp    -> Draw();
  pSum       -> Draw();
  fOutput    -> cd();
  cTrgTsp    -> Write();
  cTrgTsp    -> Close();

  // fTrg vs. hTrg
  TCanvas *cTrgEtaVsPhi = new TCanvas("cTrgEtaVsPhi", "", width, height);
  cTrgEtaVsPhi -> SetLogx(fLogX);
  cTrgEtaVsPhi -> SetLogy(fLogY2);
  cTrgEtaVsPhi -> SetLogz(fLogZ);
  cTrgEtaVsPhi -> SetGrid(fGrid, fGrid);
  cTrgEtaVsPhi -> SetTicks(fTick, fTick);
  cTrgEtaVsPhi -> SetBorderMode(fMode);
  cTrgEtaVsPhi -> SetBorderSize(fBord);
  cTrgEtaVsPhi -> SetFrameBorderMode(fFrame);
  cTrgEtaVsPhi -> SetLeftMargin(fMarginBig);
  cTrgEtaVsPhi -> SetTopMargin(fMarginSmall);
  cTrgEtaVsPhi -> SetRightMargin(fMarginBig);
  cTrgEtaVsPhi -> SetBottomMargin(fMarginBig);
  hTrgEtaVsPhi -> Draw("colz");
  pSum         -> Draw();
  fOutput      -> cd();
  cTrgEtaVsPhi -> Write();
  cTrgEtaVsPhi -> Close();

  // nFit and nFit/nPoss
  TCanvas *cTrkNfit = new TCanvas("cTrkNfit", "", width, height);
  cTrkNfit    -> SetLogx(fLogX);
  cTrkNfit    -> SetLogy(fLogY);
  cTrkNfit    -> SetGrid(fGrid, fGrid);
  cTrkNfit    -> SetTicks(fTick, fTick);
  cTrkNfit    -> SetBorderMode(fMode);
  cTrkNfit    -> SetBorderSize(fBord);
  cTrkNfit    -> SetFrameBorderMode(fFrame);
  cTrkNfit    -> SetLeftMargin(fMarginBig);
  cTrkNfit    -> SetTopMargin(fMarginSmall);
  cTrkNfit    -> SetRightMargin(fMarginSmall);
  cTrkNfit    -> SetBottomMargin(fMarginBig);
  cTrkNfit    -> cd();
  hTrkNfit[0] -> Draw();
  hTrkNfit[1] -> Draw("same");
  hTrkNfitCut -> Draw("same hist");
  lTrkNfit    -> Draw();
  pSum        -> Draw();
  fOutput     -> cd();
  cTrkNfit    -> Write();
  cTrkNfit    -> Close();

  TCanvas *cTrkNfitVsPt = new TCanvas("cTrkNfitVsPt", "", width, height);
  cTrkNfitVsPt -> SetLogx(fLogX);
  cTrkNfitVsPt -> SetLogy(fLogY2);
  cTrkNfitVsPt -> SetLogz(fLogZ);
  cTrkNfitVsPt -> SetGrid(fGrid, fGrid);
  cTrkNfitVsPt -> SetTicks(fTick, fTick);
  cTrkNfitVsPt -> SetBorderMode(fMode);
  cTrkNfitVsPt -> SetBorderSize(fBord);
  cTrkNfitVsPt -> SetFrameBorderMode(fFrame);
  cTrkNfitVsPt -> SetLeftMargin(fMarginBig);
  cTrkNfitVsPt -> SetTopMargin(fMarginSmall);
  cTrkNfitVsPt -> SetRightMargin(fMarginBig);
  cTrkNfitVsPt -> SetBottomMargin(fMarginBig);
  cTrkNfitVsPt -> cd();
  hTrkNfitVsPt -> Draw("colz");
  bPtVsNfitSel -> Draw();
  pSum         -> Draw();
  fOutput      -> cd();
  cTrkNfitVsPt -> Write();
  cTrkNfitVsPt -> Close();

  TCanvas *cTrkRfit = new TCanvas("cTrkRfit", "", width, height);
  cTrkRfit    -> SetLogx(fLogX);
  cTrkRfit    -> SetLogy(fLogY);
  cTrkRfit    -> SetGrid(fGrid, fGrid);
  cTrkRfit    -> SetTicks(fTick, fTick);
  cTrkRfit    -> SetBorderMode(fMode);
  cTrkRfit    -> SetBorderSize(fBord);
  cTrkRfit    -> SetFrameBorderMode(fFrame);
  cTrkRfit    -> SetLeftMargin(fMarginBig);
  cTrkRfit    -> SetTopMargin(fMarginSmall);
  cTrkRfit    -> SetRightMargin(fMarginSmall);
  cTrkRfit    -> SetBottomMargin(fMarginBig);
  cTrkRfit    -> cd();
  hTrkRfit[0] -> Draw();
  hTrkRfit[1] -> Draw("same");
  hTrkRfitCut -> Draw("hist same");
  lTrkRfit    -> Draw();
  pSum        -> Draw();
  fOutput     -> cd();
  cTrkRfit    -> Write();
  cTrkRfit    -> Close();

  // dca
  TCanvas *cTrkDca = new TCanvas("cTrkDca", "", width, height);
  cTrkDca    -> SetLogx(fLogX);
  cTrkDca    -> SetLogy(fLogY);
  cTrkDca    -> SetGrid(fGrid, fGrid);
  cTrkDca    -> SetTicks(fTick, fTick);
  cTrkDca    -> SetBorderMode(fMode);
  cTrkDca    -> SetBorderSize(fBord);
  cTrkDca    -> SetFrameBorderMode(fFrame);
  cTrkDca    -> SetLeftMargin(fMarginBig);
  cTrkDca    -> SetTopMargin(fMarginSmall);
  cTrkDca    -> SetRightMargin(fMarginSmall);
  cTrkDca    -> SetBottomMargin(fMarginBig);
  cTrkDca    -> cd();
  hTrkDca[0] -> Draw();
  hTrkDca[1] -> Draw("same");
  hTrkDcaCut -> Draw("hist same");
  lTrkDca    -> Draw();
  pSum       -> Draw();
  fOutput    -> cd();
  cTrkDca    -> Write();
  cTrkDca    -> Close();

  TCanvas *cTrkDcaVsPt = new TCanvas("cTrkDcaVsPt", "", width, height);
  cTrkDcaVsPt -> SetLogx(fLogX);
  cTrkDcaVsPt -> SetLogy(fLogY2);
  cTrkDcaVsPt -> SetLogz(fLogZ);
  cTrkDcaVsPt -> SetGrid(fGrid, fGrid);
  cTrkDcaVsPt -> SetTicks(fTick, fTick);
  cTrkDcaVsPt -> SetBorderMode(fMode);
  cTrkDcaVsPt -> SetBorderSize(fBord);
  cTrkDcaVsPt -> SetFrameBorderMode(fFrame);
  cTrkDcaVsPt -> SetLeftMargin(fMarginBig);
  cTrkDcaVsPt -> SetTopMargin(fMarginSmall);
  cTrkDcaVsPt -> SetRightMargin(fMarginBig);
  cTrkDcaVsPt -> SetBottomMargin(fMarginBig);
  cTrkDcaVsPt -> cd();
  hTrkDcaVsPt -> Draw("colz");
  bPtVsDcaSel -> Draw();
  pSum        -> Draw();
  fOutput     -> cd();
  cTrkDcaVsPt -> Write();
  cTrkDcaVsPt -> Close();

  // hTrk
  TCanvas *cTrkEta = new TCanvas("cTrkEta", "", width, height);
  cTrkEta    -> SetLogx(fLogX);
  cTrkEta    -> SetLogy(fLogY);
  cTrkEta    -> SetGrid(fGrid, fGrid);
  cTrkEta    -> SetTicks(fTick, fTick);
  cTrkEta    -> SetBorderMode(fMode);
  cTrkEta    -> SetBorderSize(fBord);
  cTrkEta    -> SetFrameBorderMode(fFrame);
  cTrkEta    -> SetLeftMargin(fMarginBig);
  cTrkEta    -> SetTopMargin(fMarginSmall);
  cTrkEta    -> SetRightMargin(fMarginSmall);
  cTrkEta    -> SetBottomMargin(fMarginBig);
  cTrkEta    -> cd();
  hTrkEta[0] -> Draw();
  hTrkEta[1] -> Draw("same");
  hTrkEtaCut -> Draw("hist same");
  lTrkEta    -> Draw();
  pSum       -> Draw();
  fOutput    -> cd();
  cTrkEta    -> Write();
  cTrkEta    -> Close();

  TCanvas *cTrkEtaBin        = new TCanvas("cTrkEtaBin", "", bigWidth, bigHeight);
  TPad    *pTrkEtaVsPtPi0 = new TPad("pTrkEtaVsPtPi0", "", 0., 0., 0.5, 0.5);
  TPad    *pTrkEtaVsPtGam = new TPad("pTrkEtaVsPtGam", "", 0.5, 0., 1., 0.5);
  TPad    *pTrkEtaPi0     = new TPad("pTrkEtaPi0",     "", 0., 0.5, 0.5, 1.);
  TPad    *pTrkEtaGam     = new TPad("pTrkEtaGam",     "", 0.5, 0.5, 1., 1.);
  pTrkEtaVsPtPi0 -> SetLogx(fLogX);
  pTrkEtaVsPtPi0 -> SetLogy(fLogY2);
  pTrkEtaVsPtPi0 -> SetLogz(fLogZ);
  pTrkEtaVsPtPi0 -> SetGrid(fGrid, fGrid);
  pTrkEtaVsPtPi0 -> SetTicks(fTick, fTick);
  pTrkEtaVsPtPi0 -> SetBorderMode(fMode);
  pTrkEtaVsPtPi0 -> SetBorderSize(fBord);
  pTrkEtaVsPtPi0 -> SetFrameBorderMode(fFrame);
  pTrkEtaVsPtPi0 -> SetLeftMargin(fMarginBig);
  pTrkEtaVsPtPi0 -> SetTopMargin(fMarginBord);
  pTrkEtaVsPtPi0 -> SetRightMargin(fMarginBord);
  pTrkEtaVsPtPi0 -> SetBottomMargin(fMarginBig);
  pTrkEtaVsPtGam -> SetLogx(fLogX);
  pTrkEtaVsPtGam -> SetLogy(fLogY2);
  pTrkEtaVsPtGam -> SetLogz(fLogZ);
  pTrkEtaVsPtGam -> SetGrid(fGrid, fGrid);
  pTrkEtaVsPtGam -> SetTicks(fTick, fTick);
  pTrkEtaVsPtGam -> SetBorderMode(fMode);
  pTrkEtaVsPtGam -> SetBorderSize(fBord);
  pTrkEtaVsPtGam -> SetFrameBorderMode(fFrame);
  pTrkEtaVsPtGam -> SetLeftMargin(fMarginBord);
  pTrkEtaVsPtGam -> SetTopMargin(fMarginBord);
  pTrkEtaVsPtGam -> SetRightMargin(fMarginBig);
  pTrkEtaVsPtGam -> SetBottomMargin(fMarginBig);
  pTrkEtaPi0     -> SetLogx(fLogX);
  pTrkEtaPi0     -> SetLogy(fLogY2);
  pTrkEtaPi0     -> SetGrid(fGrid, fGrid);
  pTrkEtaPi0     -> SetTicks(fTick, fTick);
  pTrkEtaPi0     -> SetBorderMode(fMode);
  pTrkEtaPi0     -> SetBorderSize(fBord);
  pTrkEtaPi0     -> SetFrameBorderMode(fFrame);
  pTrkEtaPi0     -> SetLeftMargin(fMarginBig);
  pTrkEtaPi0     -> SetTopMargin(fMarginSmall);
  pTrkEtaPi0     -> SetRightMargin(fMarginBord);
  pTrkEtaPi0     -> SetBottomMargin(fMarginBord);
  pTrkEtaGam     -> SetLogx(fLogX);
  pTrkEtaGam     -> SetLogy(fLogY2);
  pTrkEtaGam     -> SetGrid(fGrid, fGrid);
  pTrkEtaGam     -> SetTicks(fTick, fTick);
  pTrkEtaGam     -> SetBorderMode(fMode);
  pTrkEtaGam     -> SetBorderSize(fBord);
  pTrkEtaGam     -> SetFrameBorderMode(fFrame);
  pTrkEtaGam     -> SetLeftMargin(fMarginBord);
  pTrkEtaGam     -> SetTopMargin(fMarginSmall);
  pTrkEtaGam     -> SetRightMargin(fMarginBig);
  pTrkEtaGam     -> SetBottomMargin(fMarginBord);
  cTrkEtaBin     -> cd();
  pTrkEtaVsPtPi0 -> Draw();
  pTrkEtaVsPtGam -> Draw();
  pTrkEtaPi0     -> Draw();
  pTrkEtaGam     -> Draw();
  pTrkEtaVsPtPi0 -> cd();
  hTrkPtVsEta[0] -> Draw("col");
  bPtVsEtaSel    -> Draw();
  pPi0           -> Draw();
  pTrkEtaVsPtGam -> cd();
  hTrkPtVsEta[1] -> Draw("colz");
  bPtVsEtaSel    -> Draw();
  pGam           -> Draw();
  pTrkEtaPi0     -> cd();
  hTrkEtaBin[0]  -> Draw();
  for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
    hTrkEtaPtBin[0][iTrkBin] -> Draw("same");
  }
  lTrkEtaVsPtP  -> Draw();
  pTrkEtaGam    -> cd();
  hTrkEtaBin[1] -> Draw();
  for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
    hTrkEtaPtBin[1][iTrkBin] -> Draw("same");
  }
  lTrkEtaVsPtG -> Draw();
  fOutput      -> cd();
  cTrkEtaBin   -> Write();
  cTrkEtaBin   -> Close();

  // pTtrk
  TCanvas *cTrkPt = new TCanvas("cTrkPt", "", width, height);
  cTrkPt    -> SetLogx(fLogX);
  cTrkPt    -> SetLogy(fLogY);
  cTrkPt    -> SetGrid(fGrid, fGrid);
  cTrkPt    -> SetTicks(fTick, fTick);
  cTrkPt    -> SetBorderMode(fMode);
  cTrkPt    -> SetBorderSize(fBord);
  cTrkPt    -> SetFrameBorderMode(fFrame);
  cTrkPt    -> SetLeftMargin(fMarginBig);
  cTrkPt    -> SetTopMargin(fMarginSmall);
  cTrkPt    -> SetRightMargin(fMarginSmall);
  cTrkPt    -> SetBottomMargin(fMarginBig);
  cTrkPt    -> cd();
  hTrkPt[0] -> Draw();
  hTrkPt[1] -> Draw("same");
  hTrkPtCut -> Draw("hist same");
  lTrkPt    -> Draw();
  pSum      -> Draw();
  fOutput   -> cd();
  cTrkPt    -> Write();
  cTrkPt    -> Close();

  TCanvas *cTrkPtBin = new TCanvas("cTrkPtBin", "", width, height);
  cTrkPtBin    -> SetLogx(fLogX);
  cTrkPtBin    -> SetLogy(fLogY);
  cTrkPtBin    -> SetGrid(fGrid, fGrid);
  cTrkPtBin    -> SetTicks(fTick, fTick);
  cTrkPtBin    -> SetBorderMode(fMode);
  cTrkPtBin    -> SetBorderSize(fBord);
  cTrkPtBin    -> SetFrameBorderMode(fFrame);
  cTrkPtBin    -> SetLeftMargin(fMarginBig);
  cTrkPtBin    -> SetTopMargin(fMarginSmall);
  cTrkPtBin    -> SetRightMargin(fMarginSmall);
  cTrkPtBin    -> SetBottomMargin(fMarginBig);
  cTrkPtBin    -> cd();
  hTrkPtBin[0] -> Draw();
  hTrkPtBin[1] -> Draw("same");
  lTrkPtBin    -> Draw();
  pSum         -> Draw();
  fOutput      -> cd();
  cTrkPtBin    -> Write();
  cTrkPtBin    -> Close();

  // dFtrk
  TCanvas *cTrkDfBin     = new TCanvas("cTrkDfBin", "", bigWidth, bigHeight);
  TPad    *pTrkDfVsPtPi0 = new TPad("pTrkDfVsPtPi0", "", 0., 0., 0.5, 0.5);
  TPad    *pTrkDfVsPtGam = new TPad("pTrkDfVsPtGam", "", 0.5, 0., 1., 0.5);
  TPad    *pTrkDfPi0     = new TPad("pTrkDfPi0",     "", 0., 0.5, 0.5, 1.);
  TPad    *pTrkDfGam     = new TPad("pTrkDfGam",     "", 0.5, 0.5, 1., 1.);
  pTrkDfVsPtPi0 -> SetLogx(fLogX);
  pTrkDfVsPtPi0 -> SetLogy(fLogY2);
  pTrkDfVsPtPi0 -> SetLogz(fLogZ);
  pTrkDfVsPtPi0 -> SetGrid(fGrid, fGrid);
  pTrkDfVsPtPi0 -> SetTicks(fTick, fTick);
  pTrkDfVsPtPi0 -> SetBorderMode(fMode);
  pTrkDfVsPtPi0 -> SetBorderSize(fBord);
  pTrkDfVsPtPi0 -> SetFrameBorderMode(fFrame);
  pTrkDfVsPtPi0 -> SetLeftMargin(fMarginBig);
  pTrkDfVsPtPi0 -> SetTopMargin(fMarginBord);
  pTrkDfVsPtPi0 -> SetRightMargin(fMarginBord);
  pTrkDfVsPtPi0 -> SetBottomMargin(fMarginBig);
  pTrkDfVsPtGam -> SetLogx(fLogX);
  pTrkDfVsPtGam -> SetLogy(fLogY2);
  pTrkDfVsPtGam -> SetLogz(fLogZ);
  pTrkDfVsPtGam -> SetGrid(fGrid, fGrid);
  pTrkDfVsPtGam -> SetTicks(fTick, fTick);
  pTrkDfVsPtGam -> SetBorderMode(fMode);
  pTrkDfVsPtGam -> SetBorderSize(fBord);
  pTrkDfVsPtGam -> SetFrameBorderMode(fFrame);
  pTrkDfVsPtGam -> SetLeftMargin(fMarginBord);
  pTrkDfVsPtGam -> SetTopMargin(fMarginBord);
  pTrkDfVsPtGam -> SetRightMargin(fMarginBig);
  pTrkDfVsPtGam -> SetBottomMargin(fMarginBig);
  pTrkDfPi0     -> SetLogx(fLogX);
  pTrkDfPi0     -> SetLogy(fLogY2);
  pTrkDfPi0     -> SetGrid(fGrid, fGrid);
  pTrkDfPi0     -> SetTicks(fTick, fTick);
  pTrkDfPi0     -> SetBorderMode(fMode);
  pTrkDfPi0     -> SetBorderSize(fBord);
  pTrkDfPi0     -> SetFrameBorderMode(fFrame);
  pTrkDfPi0     -> SetLeftMargin(fMarginBig);
  pTrkDfPi0     -> SetTopMargin(fMarginSmall);
  pTrkDfPi0     -> SetRightMargin(fMarginBord);
  pTrkDfPi0     -> SetBottomMargin(fMarginBord);
  pTrkDfGam     -> SetLogx(fLogX);
  pTrkDfGam     -> SetLogy(fLogY2);
  pTrkDfGam     -> SetGrid(fGrid, fGrid);
  pTrkDfGam     -> SetTicks(fTick, fTick);
  pTrkDfGam     -> SetBorderMode(fMode);
  pTrkDfGam     -> SetBorderSize(fBord);
  pTrkDfGam     -> SetFrameBorderMode(fFrame);
  pTrkDfGam     -> SetLeftMargin(fMarginBord);
  pTrkDfGam     -> SetTopMargin(fMarginSmall);
  pTrkDfGam     -> SetRightMargin(fMarginBig);
  pTrkDfGam     -> SetBottomMargin(fMarginBord);
  cTrkDfBin     -> cd();
  pTrkDfVsPtPi0 -> Draw();
  pTrkDfVsPtGam -> Draw();
  pTrkDfPi0     -> Draw();
  pTrkDfGam     -> Draw();
  pTrkDfVsPtPi0 -> cd();
  hTrkPtVsDf[0] -> Draw("col");
  bPtVsDfSel    -> Draw();
  pPi0          -> Draw();
  pTrkDfVsPtGam -> cd();
  hTrkPtVsDf[1] -> Draw("colz");
  bPtVsDfSel    -> Draw();
  pGam          -> Draw();
  pTrkDfPi0     -> cd();
  hTrkDfBin[0]  -> Draw();
  for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
    hTrkDfPtBin[0][iTrkBin] -> Draw("same");
  }
  lTrkDfVsPtP  -> Draw();
  pTrkDfGam    -> cd();
  hTrkDfBin[1] -> Draw();
  for (Int_t iTrkBin = (NTrkBins - 1); iTrkBin > -1; iTrkBin--) {
    hTrkDfPtBin[1][iTrkBin] -> Draw("same");
  }
  lTrkDfVsPtG -> Draw();
  fOutput     -> cd();
  cTrkDfBin   -> Write();
  cTrkDfBin   -> Close();

  TCanvas *cTrkDfVsEta    = new TCanvas("cTrkDfVsEta", "", bigWidth, height);
  TPad    *pTrkDfVsEtaPi0 = new TPad("pTrkDfVsEtaPi0", "", 0., 0., 0.5, 1.);
  TPad    *pTrkDfVsEtaGam = new TPad("pTrkDfVsEtaGam", "", 0.5, 0., 1., 1.);
  pTrkDfVsEtaPi0 -> SetLogx(fLogX);
  pTrkDfVsEtaPi0 -> SetLogy(fLogY2);
  pTrkDfVsEtaPi0 -> SetLogz(fLogZ);
  pTrkDfVsEtaPi0 -> SetGrid(fGrid, fGrid);
  pTrkDfVsEtaPi0 -> SetTicks(fTick, fTick);
  pTrkDfVsEtaPi0 -> SetBorderMode(fMode);
  pTrkDfVsEtaPi0 -> SetBorderSize(fBord);
  pTrkDfVsEtaPi0 -> SetFrameBorderMode(fFrame);
  pTrkDfVsEtaPi0 -> SetLeftMargin(fMarginBig);
  pTrkDfVsEtaPi0 -> SetTopMargin(fMarginSmall);
  pTrkDfVsEtaPi0 -> SetRightMargin(fMarginBord);
  pTrkDfVsEtaPi0 -> SetBottomMargin(fMarginBig);
  pTrkDfVsEtaGam -> SetLogx(fLogX);
  pTrkDfVsEtaGam -> SetLogy(fLogY2);
  pTrkDfVsEtaGam -> SetLogz(fLogZ);
  pTrkDfVsEtaGam -> SetGrid(fGrid, fGrid);
  pTrkDfVsEtaGam -> SetTicks(fTick, fTick);
  pTrkDfVsEtaGam -> SetBorderMode(fMode);
  pTrkDfVsEtaGam -> SetBorderSize(fBord);
  pTrkDfVsEtaGam -> SetFrameBorderMode(fFrame);
  pTrkDfVsEtaGam -> SetLeftMargin(fMarginBord);
  pTrkDfVsEtaGam -> SetTopMargin(fMarginSmall);
  pTrkDfVsEtaGam -> SetRightMargin(fMarginBig);
  pTrkDfVsEtaGam -> SetBottomMargin(fMarginBig);
  cTrkDfVsEta    -> cd();
  pTrkDfVsEtaPi0 -> Draw();
  pTrkDfVsEtaGam -> Draw();
  pTrkDfVsEtaPi0 -> cd();
  hTrkDfVsEta[0] -> Draw("col");
  pPi0           -> Draw();
  pTrkDfVsEtaGam -> cd();
  hTrkDfVsEta[1] -> Draw("colz");
  pGam           -> Draw();
  fOutput        -> cd();
  cTrkDfVsEta    -> Write();
  cTrkDfVsEta    -> Close();
  cout << "      Made plots." << endl;

  // make directories
  TDirectory *dEvt = (TDirectory*) fOutput -> mkdir("EventQA");
  TDirectory *dTrg = (TDirectory*) fOutput -> mkdir("TriggerQA");
  TDirectory *dTrk = (TDirectory*) fOutput -> mkdir("TrackQA");

  // close files
  fOutput            -> cd();
  dEvt               -> cd();
  hEvtNum            -> Write();
  hEvtVz[0]          -> Write();
  hEvtVz[1]          -> Write();
  hEvtVr[0]          -> Write();
  hEvtVr[1]          -> Write();
  hEvtPrim           -> Write();
  hEvtVyVsVx         -> Write();
  hEvtVzVsVx         -> Write();
  hEvtVzVsVy         -> Write();
  dTrg               -> cd();
  hTrgEt[0]          -> Write();
  hTrgEt[1]          -> Write();
  hTrgEta            -> Write();
  hTrgPhi            -> Write();
  hTrgTsp[0]         -> Write();
  hTrgTsp[1]         -> Write();
  hTrgTsp[2]         -> Write();
  hTrgEtBin[0][0]    -> Write();
  hTrgEtBin[0][1]    -> Write();
  hTrgEtBin[1][0]    -> Write();
  hTrgEtBin[1][1]    -> Write();
  hTrgEtBin[2][0]    -> Write();
  hTrgEtBin[2][1]    -> Write();
  hTrgEtBin[3][0]    -> Write();
  hTrgEtBin[3][1]    -> Write();
  hTrgEtSum          -> Write();
  hTrgEtVsEta        -> Write();
  hTrgEtVsPhi        -> Write();
  hTrgEtaVsPhi       -> Write();
  dTrk               -> cd();
  hTrkNfit[0]        -> Write();
  hTrkNfit[1]        -> Write();
  hTrkRfit[0]        -> Write();
  hTrkRfit[1]        -> Write();
  hTrkDca[0]         -> Write();
  hTrkDca[1]         -> Write();
  hTrkEta[0]         -> Write();
  hTrkEta[1]         -> Write();
  hTrkPt[0]          -> Write();
  hTrkPt[1]          -> Write();
  hTrkPtBin[0]       -> Write();
  hTrkPtBin[1]       -> Write();
  hTrkDfBin[0]       -> Write();
  hTrkDfBin[1]       -> Write();
  hTrkEtaBin[0]      -> Write();
  hTrkEtaBin[1]      -> Write();
  hTrkDfPtBin[0][0]  -> Write();
  hTrkDfPtBin[0][1]  -> Write();
  hTrkDfPtBin[0][2]  -> Write();
  hTrkDfPtBin[0][3]  -> Write();
  hTrkDfPtBin[1][0]  -> Write();
  hTrkDfPtBin[1][1]  -> Write();
  hTrkDfPtBin[1][2]  -> Write();
  hTrkDfPtBin[1][3]  -> Write();
  hTrkEtaPtBin[0][0] -> Write();
  hTrkEtaPtBin[0][1] -> Write();
  hTrkEtaPtBin[0][2] -> Write();
  hTrkEtaPtBin[0][3] -> Write();
  hTrkEtaPtBin[1][0] -> Write();
  hTrkEtaPtBin[1][1] -> Write();
  hTrkEtaPtBin[1][2] -> Write();
  hTrkEtaPtBin[1][3] -> Write();
  hTrkNfitVsPt       -> Write();
  hTrkDcaVsPt        -> Write();
  hTrkPtVsDf[0]      -> Write();
  hTrkPtVsDf[1]      -> Write();
  hTrkPtVsEta[0]     -> Write();
  hTrkPtVsEta[1]     -> Write();
  hTrkDfVsEta[0]     -> Write();
  hTrkDfVsEta[1]     -> Write();
  pTrkNfitVsPt       -> Write();
  pTrkDcaVsPt        -> Write();
  fOutput            -> cd();
  fOutput            -> Close();
  fInput             -> cd();
  fInput             -> Close();
  cout << "    Made event QA plots!" << endl;

}

// End ------------------------------------------------------------------------

