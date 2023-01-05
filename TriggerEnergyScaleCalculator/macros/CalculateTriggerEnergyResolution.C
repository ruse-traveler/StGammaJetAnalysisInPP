// 'CalculateTriggerEnergyResolution.C'
// Derek Anderson
// 05.18.2020
//
// Calculates the energy resolution of
// identified pi0 triggers.
//
// NOTE: assumes particle-level and
//   detector-level input trees were
//   filled in same order.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t NEvtMax(100000);
static const UInt_t NTrkMax(5000);
static const UInt_t NTwrMax(5000);
static const UInt_t NMatchMax(10);
static const UInt_t NHotTwrs(302);
static const UInt_t NBadRuns(45);
static const UInt_t NTrgBin(3);
static const UInt_t NTrgTsp(2);
static const UInt_t NHadIds(1);



void CalculateTriggerEnergyResolution(const Bool_t isInBatchMode=false) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning trigger energy resolution script..." << endl;


  // io parameters
  //const TString sOut("pp200r9pt35rff.eTtrgResolutionStudy_forEtDetVsEtParPlot.et650x650vz55pi0.d11m6y2020.root");
  const TString sOut("test920x650.root");
  const TString sPar("output/merged/pt35rff.thirdMakerWithMc.root");
  const TString sDet("output/merged/pt35rff.thirdMakerWithMc.root");
  const TString sParTree("McTracks");
  const TString sDetTree("Gfmtodst");

  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);

  // detector trigger paramters
  const Int_t    adcMax(6004);
  const Double_t eStrMin(0.5);
  const Double_t pProjMax(3.);
  const Double_t hTrgDetMax(0.9);
  const Double_t eTtrgDetMin(9.);
  const Double_t eTtrgDetMax(20.);
  const Double_t eTbinMin[NTrgBin] = {9., 11., 15.};
  const Double_t eTbinMax[NTrgBin] = {11., 15., 20.};
  const Double_t tspPi0[NTrgTsp]   = {0., 0.08};

  // particle trigger parameters
  const Double_t cTrgPar(0.);
  const Double_t dRtrgMax(0.1);
  const Double_t hTrgParMax(0.9);
  const Double_t eTtrgParMin(6.);
  const Double_t eTtrgParMax(50.);
  const UInt_t   idHadTrg[NHadIds] = {7};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fPar = new TFile(sPar.Data(), "read");
  TFile *fDet = new TFile(sDet.Data(), "read");
  if (!fPar || !fDet) {
    cerr << "PANIC: couldn't open an input file!\n"
         << "  fPar = " << fPar << "\n"
         << "  fDet = " << fDet
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input trees
  TTree *tPar;
  TTree *tDet;
  fPar -> GetObject(sParTree.Data(), tPar);
  fDet -> GetObject(sDetTree.Data(), tDet);
  if (!tPar || !tDet) {
    cerr << "PANIC: couldn't grab an input tree!\n"
         << "  tPar = " << tPar << "\n"
         << "  tDet = " << tDet
         << endl;
    return;
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

  // declare detector (event) leaf types
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
  // declare detector (track) leaf types
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
  // declare detector (tower) leaf types
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
  cout << "    Declared input branches." << endl;


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
  tDet -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &FlagEvent_TrgTrkMisMtch);
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
  cout << "    Set input branch addresses." << endl;


  // define bad run and hot tower lists
  const UInt_t badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  const UInt_t hotTwrList[NHotTwrs] = {34, 106, 113, 160, 266, 267, 275, 280, 282, 286, 287, 293, 410, 504, 533, 541, 555, 561, 562, 594, 615, 616, 629, 633, 637, 638, 647, 650, 653, 657, 671, 673, 743, 789, 790, 791, 792, 806, 809, 810, 811, 812, 813, 814, 821, 822, 823, 824, 829, 830, 831, 832, 837, 841, 842, 843, 844, 846, 849, 850, 851, 852, 857, 875, 897, 899, 903, 939, 953, 954, 956, 993, 1026, 1046, 1048, 1080, 1081, 1100, 1125, 1130, 1132, 1180, 1197, 1198, 1199, 1200, 1207, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1312, 1348, 1353, 1354, 1388, 1407, 1409, 1434, 1448, 1537, 1567, 1574, 1597, 1612, 1654, 1668, 1713, 1762, 1765, 1766, 1877, 1878, 1984, 2032, 2043, 2054, 2073, 2077, 2092, 2093, 2097, 2107, 2162, 2168, 2214, 2305, 2392, 2409, 2415, 2439, 2459, 2589, 2590, 2633, 2652, 2749, 2834, 2961, 2969, 3005, 3017, 3070, 3071, 3186, 3220, 3289, 3360, 3493, 3494, 3495, 3508, 3588, 3604, 3611, 3668, 3678, 3679, 3690, 3692, 3732, 3738, 3838, 3840, 3927, 3945, 4005, 4006, 4013, 4018, 4019, 4053, 4059, 4124, 4331, 4355, 4357, 4458, 4464, 4500, 4677, 4678, 4684, 4768, 360, 493, 779, 1284, 1306, 1337, 1438, 1709, 2027, 2445, 3407, 3720, 4217, 4288, 95, 96, 296, 316, 443, 479, 555, 562, 637, 671, 709, 740, 743, 796, 857, 897, 899, 915, 953, 1130, 1132, 1294, 1318, 1337, 1348, 1359, 1378, 1427, 1429, 1440, 1537, 1563, 1574, 1709, 1763, 1773, 1819, 1854, 1874, 1936, 1938, 2018, 2043, 2098, 2099, 2256, 2259, 2294, 2514, 2520, 2552, 2589, 2598, 2680, 2706, 2799, 2880, 2897, 2917, 2969, 3020, 3028, 3310, 3319, 3375, 3399, 3504, 3539, 3541, 3679, 3690, 3692, 3718, 3719, 3720, 3738, 3806, 3838, 3840, 3928, 4013, 4017, 4038, 4053, 4057, 4058, 4079, 4097, 4099};
  cout << "    Bad run and hot tower lists defined:\n"
       << "      " << NBadRuns << " bad runs, " << NHotTwrs << " hot towers."
       << endl;


  // create histograms
  TH1D *hPhiEff;
  TH1D *hEtaEff;
  TH1D *hEtEff;
  TH1D *hPhiForEff[2];
  TH1D *hEtaForEff[2];
  TH1D *hEtForEff[2];
  TH1D *hPhiTrg[3];
  TH1D *hEtaTrg[3];
  TH1D *hEtTrg[3];
  TH1D *hDrTrg[2];
  TH1D *hQtTrg[2];
  TH1D *hDqTrg[2];
  TH1D *hEtMatch[NTrgBin];
  TH1D *hDrMatch[NTrgBin];
  TH1D *hQtMatch[NTrgBin];
  TH1D *hDqMatch[NTrgBin];
  TH2D *hPhiMatchVsDet;
  TH2D *hEtaMatchVsDet;
  TH2D *hEtMatchVsDet;
  TH2D *hDrMatchVsEtDet;
  TH2D *hQtMatchVsEtDet;
  TH2D *hDqMatchVsEtDet;

  const UInt_t  nPhi(360);
  const UInt_t  nEta(8);
  const UInt_t  nEt(200);
  const UInt_t  nDr(50);
  const UInt_t  nQt(100);
  const UInt_t  nDq(50);
  const Float_t phi[2] = {-6.3, 6.3};
  const Float_t eta[2] = {-2., 2.};
  const Float_t eT[2]  = {0., 100.};
  const Float_t dR[2]  = {0., 5.};
  const Float_t qT[2]  = {0., 10.};
  const Float_t dQ[2]  = {0., 5.};
  // matching efficiences
  hPhiEff         = new TH1D("hPhiEff", "#varphi^{trg} matching efficiency", nPhi, phi[0], phi[1]);
  hEtaEff         = new TH1D("hEtaEff", "#eta^{trg} matching efficiency", nEta, eta[0], eta[1]);
  hEtEff          = new TH1D("hEtEff", "E_{T}^{trg} matching efficiency", nEt, eT[0], eT[1]);
  hPhiForEff[0]   = new TH1D("hPhiForEffDet", "Detector #varphi^{trg} (for #epsilon_{match})", nPhi, phi[0], phi[1]);
  hPhiForEff[1]   = new TH1D("hPhiForEffMatch", "Matched #varphi^{trg} (for #epsilon_{match})", nPhi, phi[0], phi[1]);
  hEtaForEff[0]   = new TH1D("hEtaForEffDet", "Detector #eta^{trg} (for #epsilon_{match})", nEta, eta[0], eta[1]);
  hEtaForEff[1]   = new TH1D("hEtaForEffMatch", "Matched #eta^{trg} (for #epsilon_{match})", nEta, eta[0], eta[1]);
  hEtForEff[0]    = new TH1D("hEtForEffDet", "Detector E_{T}^{trg} (for #epsilon_{match})", nEt, eT[0], eT[1]);
  hEtForEff[1]    = new TH1D("hEtForEffMatch", "Matched E_{T}^{trg} (for #epsilon_{match})", nEt, eT[0], eT[1]); 
  // trigger phi
  hPhiTrg[0]      = new TH1D("hPhiTrgDet", "Detector #varphi^{trg}", nPhi, phi[0], phi[1]);
  hPhiTrg[1]      = new TH1D("hPhiTrgCand", "Candidate match #varphi^{trg}", nPhi, phi[0], phi[1]);
  hPhiTrg[2]      = new TH1D("hPhiTrgMatch", "Matched particle #varphi^{trg}", nPhi, phi[0], phi[1]);
  // trigger eta
  hEtaTrg[0]      = new TH1D("hEtaTrgDet", "Detector #eta^{trg}", nEta, eta[0], eta[1]);
  hEtaTrg[1]      = new TH1D("hEtaTrgCand", "Candidate match #eta^{trg}", nEta, eta[0], eta[1]);
  hEtaTrg[2]      = new TH1D("hEtaTrgMatch", "Matched particle #eta^{trg}", nEta, eta[0], eta[1]);
  // trigger eT
  hEtTrg[0]       = new TH1D("hEtTrgDet", "Detector E_{T}^{trg}", nEt, eT[0], eT[1]);
  hEtTrg[1]       = new TH1D("hEtTrgCand", "Candidate match E_{T}^{trg}", nEt, eT[0], eT[1]);
  hEtTrg[2]       = new TH1D("hEtTrgMatch", "Matched particle E_{T}^{trg}", nEt, eT[0], eT[1]);
  hEtMatch[0]     = new TH1D("hEtMatch911", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nEt, eT[0], eT[1]);
  hEtMatch[1]     = new TH1D("hEtMatch1115", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nEt, eT[0], eT[1]);
  hEtMatch[2]     = new TH1D("hEtMatch1520", "E_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nEt, eT[0], eT[1]);
  // trigger dR
  hDrTrg[0]       = new TH1D("hDrTrgCand", "Candidate match #DeltaR^{trg}", nDr, dR[0], dR[1]);
  hDrTrg[1]       = new TH1D("hDrTrgMatch", "Matched particle #DeltaR^{trg}", nDr, dR[0], dR[1]);
  hDrMatch[0]     = new TH1D("hDrMatch911", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nDr, dR[0], dR[1]);
  hDrMatch[1]     = new TH1D("hDrMatch1115", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nDr, dR[0], dR[1]);
  hDrMatch[2]     = new TH1D("hDrMatch1520", "#DeltaR^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nDr, dR[0], dR[1]);
  // trigger qT
  hQtTrg[0]       = new TH1D("hQtTrgCand", "Candidate match q_{T}^{trg}", nQt, qT[0], qT[1]);
  hQtTrg[1]       = new TH1D("hQtTrgMatch", "Matched particle q_{T}^{trg}", nQt, qT[0], qT[1]);
  hQtMatch[0]     = new TH1D("hQtMatch911", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nQt, qT[0], qT[1]);
  hQtMatch[1]     = new TH1D("hQtMatch1115", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nQt, qT[0], qT[1]);
  hQtMatch[2]     = new TH1D("hQtMatch1520", "q_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nQt, qT[0], qT[1]);
  // trigger dQ
  hDqTrg[0]       = new TH1D("hDqTrgCand", "Candidate match #Deltaq_{T}^{trg}", nDq, dQ[0], dQ[1]);
  hDqTrg[1]       = new TH1D("hDqTrgMatch", "Matched particle #Deltaq_{T}^{trg}", nDq, dQ[0], dQ[1]);
  hDqMatch[0]     = new TH1D("hDqMatch911", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 9 - 11 GeV", nQt, qT[0], qT[1]);
  hDqMatch[1]     = new TH1D("hDqMatch1115", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 11 - 15 GeV", nQt, qT[0], qT[1]);
  hDqMatch[2]     = new TH1D("hDqMatch1520", "#Deltaq_{T}^{trg}(match), E_{T}^{trg}(det) = 15 - 20 GeV", nQt, qT[0], qT[1]);
  // trigger resolutions
  hPhiMatchVsDet  = new TH2D("hPhiMatchVsDet", "#varphi^{trg}(match) vs. #varphi^{trg}(det)", nPhi, phi[0], phi[1], nPhi, phi[0], phi[1]);
  hEtaMatchVsDet  = new TH2D("hEtaMatchVsDet", "#eta^{trg}(match) vs. #eta^{trg}(det)", nEta, eta[0], eta[1], nEta, eta[0], eta[1]);
  hEtMatchVsDet   = new TH2D("hEtMatchVsDet", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], nEt, eT[0], eT[1]);
  hDrMatchVsEtDet = new TH2D("hDrMatchVsEtDet", "#DeltaR^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], nDr, dR[0], dR[1]);
  hQtMatchVsEtDet = new TH2D("hQtMatchVsEtDet", "q_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], nQt, qT[0], qT[1]);
  hDqMatchVsEtDet = new TH2D("hDqMatchVsEtDet", "#Deltaq_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], nDq, dQ[0], dQ[1]);
  // errors
  hPhiEff         -> Sumw2();
  hEtaEff         -> Sumw2();
  hEtEff          -> Sumw2();
  hPhiForEff[0]   -> Sumw2();
  hPhiForEff[1]   -> Sumw2();
  hEtaForEff[0]   -> Sumw2();
  hEtaForEff[1]   -> Sumw2();
  hEtForEff[0]    -> Sumw2();
  hEtForEff[1]    -> Sumw2();
  hPhiTrg[0]      -> Sumw2();
  hPhiTrg[1]      -> Sumw2();
  hPhiTrg[2]      -> Sumw2();
  hEtaTrg[0]      -> Sumw2();
  hEtaTrg[1]      -> Sumw2();
  hEtaTrg[2]      -> Sumw2();
  hEtTrg[0]       -> Sumw2();
  hEtTrg[1]       -> Sumw2();
  hEtTrg[2]       -> Sumw2();
  hEtMatch[0]     -> Sumw2();
  hEtMatch[1]     -> Sumw2();
  hEtMatch[2]     -> Sumw2();
  hDrTrg[0]       -> Sumw2();
  hDrTrg[1]       -> Sumw2();
  hDrMatch[0]     -> Sumw2();
  hDrMatch[1]     -> Sumw2();
  hDrMatch[2]     -> Sumw2();
  hQtTrg[0]       -> Sumw2();
  hQtTrg[1]       -> Sumw2();
  hQtMatch[0]     -> Sumw2();
  hQtMatch[1]     -> Sumw2();
  hQtMatch[2]     -> Sumw2();
  hDqTrg[0]       -> Sumw2();
  hDqTrg[1]       -> Sumw2();
  hDqMatch[0]     -> Sumw2();
  hDqMatch[1]     -> Sumw2();
  hDqMatch[2]     -> Sumw2();
  hPhiMatchVsDet  -> Sumw2();
  hEtaMatchVsDet  -> Sumw2();
  hEtMatchVsDet   -> Sumw2();
  hDrMatchVsEtDet -> Sumw2();
  hQtMatchVsEtDet -> Sumw2();
  hDqMatchVsEtDet -> Sumw2();

  // TEST [06.12.2020]
  TH1D *hNumMatch = new TH1D("hNumMatch", "", 10, 0, 10);
  hNumMatch -> Sumw2();

  // trigger profiles
  TProfile *pPhiMatchVsDet  = new TProfile("pPhiMatchVsDet", "#varphi^{trg}(match) vs. #varphi^{trg}(det)", nPhi, phi[0], phi[1], "S");
  TProfile *pEtaMatchVsDet  = new TProfile("pEtaMatchVsDet", "#eta^{trg}(match) vs. #eta^{trg}(det)", nEta, eta[0], eta[1], "S");
  TProfile *pEtMatchVsDet   = new TProfile("pEtMatchVsDet", "E_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], "S");
  TProfile *pDrMatchVsEtDet = new TProfile("pDrMatchVsEtDet", "#DeltaR^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], "S");
  TProfile *pQtMatchVsEtDet = new TProfile("pQtMatchVsEtDet", "q_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], "S");
  TProfile *pDqMatchVsEtDet = new TProfile("pDqMatchVsEtDet", "#Deltaq_{T}^{trg}(match) vs. E_{T}^{trg}(det)", nEt, eT[0], eT[1], "S");
  cout << "    Made histograms and profiles." << endl;


  // begin event loop
  const UInt_t nParEvts = tPar -> GetEntriesFast();
  const UInt_t nDetEvts = tDet -> GetEntriesFast();
  if (nParEvts != nDetEvts) {
    cerr << "PANIC: number of particle and detector events are different!\n"
         << "  nParEvts = " << nParEvts << "\n"
         << "  nDetEvts = " << nDetEvts
         << endl;
    return;
  }
  else {
    cout << "    Beginning event loop: " << nParEvts << " generated and " << nDetEvts << " reconstructed events to process." << endl;
  }

  UInt_t pBytes(0);
  UInt_t dBytes(0);
  UInt_t nParBytes(0);
  UInt_t nDetBytes(0);
  UInt_t nEvtMatch(0);
  UInt_t nTrgPar(0);
  UInt_t nTrgDet(0);
  for (UInt_t iEvt = 0; iEvt < nParEvts; iEvt++) {

    // load entry
    pBytes     = tPar -> GetEntry(iEvt);
    dBytes     = tDet -> GetEntry(iEvt);
    nParBytes += pBytes;
    nDetBytes += dBytes;
    if ((pBytes < 0) || (dBytes < 0)) {
      cerr << "WARNING: issue with entry " << iEvt << "!\n"
           << "  no. of par. bytes = " << pBytes << "\n"
           << "  no. of det. bytes = " << dBytes
           << endl;
      break;
    }
    else {
      if (isInBatchMode) {
        //cout << "      Processing event " << iEvt + 1 << "/" << nParEvts << "..." << endl;
      }
      else {
        //cout << "      Processing event " << iEvt + 1 << "/" << nParEvts << "...\r" << flush;
        //if ((iEvt + 1) == nParEvts) cout << endl;
      }
    }


    // particle event info
    const Int_t    runIdPar = mcRunId;
    const Int_t    evtIdPar = mcEventId;
    const UInt_t   nTrkPar  = mcNumTrks;
    const Double_t zVtxPar  = mcVz;
    const Double_t rVtxPar  = TMath::Sqrt((mcVx * mcVx) + (mcVy * mcVy));

    // detector event info
    const UInt_t   runIdDet = runNumber;
    const UInt_t   evtIdDet = eventNumber;
    const Double_t zVtxDet  = zVertex;
    const Double_t rVtxDet  = TMath::Sqrt((xVertex * xVertex) + (yVertex * yVertex));

    // check if evnets are the same
    const Bool_t isSameRun = (runIdPar == runIdDet);
    const Bool_t isSameEvt = (evtIdPar == evtIdDet);
    if (!isSameRun || !isSameEvt) {
      cerr << "WARNING: particle and detector events have different run or event IDs!\n"
           << "  event number    = " << iEvt << "\n"
           << "  runId(Par, Det) = (" << runIdPar << ", " << runIdDet << ")\n"
           << "  evtId(Par, Det) = (" << evtIdPar << ", " << evtIdDet << ")"
           << endl;
      continue;
    }

    // filter out bad runs
    Bool_t isGoodRunPar(true);
    Bool_t isGoodRunDet(true);
    for (UInt_t iRun = 0; iRun < NBadRuns; iRun++) {
      if (runIdPar == badRunList[iRun]) isGoodRunPar = false;
      if (runIdDet == badRunList[iRun]) isGoodRunDet = false;
    }

    // detector event cuts
    const Bool_t isInZcutPar    = (TMath::Abs(zVtxPar) < zVtxMax);
    const Bool_t isInZcutDet    = (TMath::Abs(zVtxDet) < zVtxMax);
    const Bool_t isInRcutPar    = (TMath::Abs(rVtxPar) < rVtxMax);
    const Bool_t isInRcutDet    = (TMath::Abs(rVtxDet) < rVtxMax);
    const Bool_t isGoodEventPar = (isGoodRunPar && isInRcutPar && isInZcutPar);
    const Bool_t isGoodEventDet = (isGoodRunDet && isInRcutDet && isInZcutDet);
    const Bool_t isGoodEvent    = (isGoodEventPar && isGoodEventDet);
    if (!isGoodEvent)
      continue;
    else
      ++nEvtMatch;


    // detector trigger info
    const Int_t    adcDet    = ETwradc11;
    const UInt_t   idTrgDet  = ETwrdidT;
    const Double_t tspTrgDet = Etsp;
    const Double_t eStrEta   = EEstrpen4;
    const Double_t eStrPhi   = EPstripenp4;
    const Double_t pProjDet  = ETwrPTower;
    const Double_t hTrgDet   = ETwreT;
    const Double_t hPhysDet  = EClustetav1;
    const Double_t fTrgDet   = EClustphiv1;
    const Double_t eTrgDet   = EClustEneT0;
    const Double_t tTrgDet   = 2. * TMath::ATan(TMath::Exp(-1. * hPhysDet));
    const Double_t eTtrgDet  = eTrgDet * TMath::Sin(tTrgDet);

    // filter out bad towers
    Bool_t isGoodTwrDet = true;
    for (UInt_t iTwr = 0; iTwr < NHotTwrs; iTwr++) {
      if (idTrgDet == hotTwrList[iTwr]) {
        isGoodTwrDet = false;
        break;
      }
    }

    // detector trigger cuts
    const Bool_t isInAdcCutDet    = (adcDet <= adcMax);
    const Bool_t isInStrCutDet    = ((eStrEta >= eStrMin) && (eStrPhi >= eStrMin));
    const Bool_t isInProjCutDet   = (pProjDet < pProjMax);
    const Bool_t isInEtaTrgCutDet = ((TMath::Abs(hTrgDet) < hTrgDetMax) && (TMath::Abs(hPhysDet) < hTrgDetMax));
    const Bool_t isInEtCutDet     = ((eTtrgDet >= eTtrgDetMin) && (eTtrgDet < eTtrgDetMax));
    const Bool_t isInPi0cutDet    = ((tspTrgDet > tspPi0[0]) && (tspTrgDet < tspPi0[1]));
    const Bool_t isGoodTrgDet     = (isGoodTwrDet && isInAdcCutDet && isInStrCutDet && isInProjCutDet && isInEtaTrgCutDet && isInEtCutDet && isInPi0cutDet);
    if (!isGoodTrgDet) {
      continue;
    }
    else {
      hPhiForEff[0] -> Fill(fTrgDet);
      hEtaForEff[0] -> Fill(hTrgDet);
      hEtForEff[0]  -> Fill(eTtrgDet);
      hPhiTrg[0]    -> Fill(fTrgDet);
      hEtaTrg[0]    -> Fill(hTrgDet);
      hEtTrg[0]     -> Fill(eTtrgDet);
      ++nTrgDet;
    }


    // TEST [06.12.2020
    UInt_t nMatch(0);

    // particle trigger loop
    Bool_t   foundTrgPar(false);
    UInt_t   iTrgPar(nTrkPar + 1);
    Double_t fTrgPar(0.);
    Double_t hTrgPar(0.);
    Double_t eTtrgPar(0.);
    Double_t dRtrgPar(0.);
    Double_t qTtrgPar(-1.);
    Double_t dQtTrgPar(-1.);
    for (UInt_t iTrgMC = 0; iTrgMC < nTrkPar; iTrgMC++) {

      // particle trigger info
      const Int_t    pidTrgMc  = mcIdGeant -> at(iTrgMC);
      const Double_t cTrgMc    = mcCharge  -> at(iTrgMC);
      const Double_t fTrgMc    = mcPhi     -> at(iTrgMC);
      const Double_t hTrgMc    = mcEta     -> at(iTrgMC);
      const Double_t pTtrgMc   = mcPt      -> at(iTrgMC);

      // delta-r and qT relative to detector trigger
      const Double_t dFtrgMc  = fTrgMc - fTrgDet;
      const Double_t dHtrgMc  = hTrgMc - hPhysDet;
      const Double_t dRtrgMc  = TMath::Sqrt((dFtrgMc * dFtrgMc) + (dHtrgMc * dHtrgMc));
      const Double_t qTtrgMc  = eTtrgDet / pTtrgMc;
      const Double_t dQtTrgMc = TMath::Abs(eTtrgDet - pTtrgMc) / pTtrgMc;

      // check pid
      Bool_t isHadronPar(false);
      for (UInt_t iTrgHad = 0; iTrgHad < NHadIds; iTrgHad++) {
        Bool_t isSameId = (pidTrgMc == idHadTrg[iTrgHad]);
        if (isSameId) {
          isHadronPar = true;
          break;
        }
      }

      // particle trigger cuts
      const Bool_t isNeutralTrgPar  = (cTrgMc == cTrgPar);
      const Bool_t isInTrgDrCutPar  = (dRtrgMc < dRtrgMax);
      const Bool_t isInTrgEtaCutPar = (TMath::Abs(hTrgMc) < hTrgParMax);
      const Bool_t isInTrgEtCutPar  = ((pTtrgMc > eTtrgParMin) && (pTtrgMc < eTtrgParMax));
      if (isHadronPar && isNeutralTrgPar && isInTrgDrCutPar && isInTrgEtaCutPar && isInTrgEtCutPar) {

        // fill candidate histograms
        hPhiTrg[1] -> Fill(fTrgMc);
        hEtaTrg[1] -> Fill(hTrgMc);
        hEtTrg[1]  -> Fill(pTtrgMc);
        hDrTrg[0]  -> Fill(dRtrgMc);
        hQtTrg[0]  -> Fill(qTtrgMc);
        hDqTrg[0]  -> Fill(dQtTrgMc);

        // apply matching condition
        if (isInTrgDrCutPar) {

          // TEST [06.12.2020]
          if (qTtrgMc > 1.5) {
            cout << "  High qTtrg match!\n"
                 << "    (id, vx, end) = (" << pidTrgMc << ", " << (mcIdVx -> at(iTrgMC)) << ", " << (mcIdVxEnd -> at(iTrgMC)) << ")\n"
                 << "    (pT, h, f)    = (" << pTtrgMc << ", " << hTrgMc << ", " << fTrgMc << ")\n"
                 << "    (qT, dR)      = (" << qTtrgMc << ", " << dRtrgMc << ")"
                 << endl;
          }

          iTrgPar     = iTrgMC;
          fTrgPar     = fTrgMc;
          hTrgPar     = hTrgMc;
          eTtrgPar    = pTtrgMc;
          dRtrgPar    = dRtrgMc;
          qTtrgPar    = qTtrgMc;
          dQtTrgPar   = dQtTrgMc;
          foundTrgPar = true;

          // TEST [06.12.2020]
          //break;
          ++nMatch;
        }  // end matching condition
      }  // end particle trigger cuts
    }  // end particle trigger loop
    if (!foundTrgPar) {
      continue;
    }
    else {
      hPhiForEff[1]   -> Fill(fTrgDet);
      hEtaForEff[1]   -> Fill(hTrgDet);
      hEtForEff[1]    -> Fill(eTtrgDet);
      hPhiTrg[2]      -> Fill(fTrgDet);
      hEtaTrg[2]      -> Fill(hTrgDet);
      hEtTrg[2]       -> Fill(eTtrgDet);
      hDrTrg[1]       -> Fill(dRtrgPar);
      hQtTrg[1]       -> Fill(qTtrgPar);
      hDqTrg[1]       -> Fill(dQtTrgPar);
      hPhiMatchVsDet  -> Fill(fTrgDet, fTrgPar);
      pPhiMatchVsDet  -> Fill(fTrgDet, fTrgPar);
      hEtaMatchVsDet  -> Fill(hTrgDet, hTrgPar);
      pEtaMatchVsDet  -> Fill(hTrgDet, hTrgPar);
      hEtMatchVsDet   -> Fill(eTtrgDet, eTtrgPar);
      pEtMatchVsDet   -> Fill(eTtrgDet, eTtrgPar);
      hDrMatchVsEtDet -> Fill(eTtrgDet, dRtrgPar);
      pDrMatchVsEtDet -> Fill(eTtrgDet, dRtrgPar);
      hQtMatchVsEtDet -> Fill(eTtrgDet, qTtrgPar);
      pQtMatchVsEtDet -> Fill(eTtrgDet, qTtrgPar);
      hDqMatchVsEtDet -> Fill(eTtrgDet, dQtTrgPar);
      pDqMatchVsEtDet -> Fill(eTtrgDet, dQtTrgPar);
      ++nTrgPar;
    }

    // TEST [06.12.2020]
    cout << "Num Match = " << nMatch << endl;
    hNumMatch -> Fill(nMatch);

    // fill eTtrg-dependent histograms
    for (UInt_t iEtBin = 0; iEtBin < NTrgBin; iEtBin++) {
      if ((eTtrgDet > eTbinMin[iEtBin]) && (eTtrgDet < eTbinMax[iEtBin])) {
        hEtMatch[iEtBin] -> Fill(eTtrgPar);
        hDrMatch[iEtBin] -> Fill(dRtrgPar);
        hQtMatch[iEtBin] -> Fill(qTtrgPar);
        hDqMatch[iEtBin] -> Fill(dQtTrgPar);
      }
    }  // end eTtrg loop
  }  // end event loop
  cout << "    Event loop finished:\n"
       << "      nEvtMatch = " << nEvtMatch << "\n"
       << "      nTrgDet   = " << nTrgDet << "\n"
       << "      nTrgPar   = " << nTrgPar
       << endl;


  // calculate efficiency
  const Float_t parWeight(1.);
  const Float_t detWeight(1.);
  hPhiEff -> Divide(hPhiForEff[1], hPhiForEff[0], parWeight, detWeight);
  hEtaEff -> Divide(hEtaForEff[1], hEtaForEff[0], parWeight, detWeight);
  hEtEff  -> Divide(hEtForEff[1], hEtForEff[0], parWeight, detWeight);
  cout << "    Calculated matching efficiencies." << endl;


  // make directories
  TDirectory *dDet = (TDirectory*) fOut -> mkdir("detector");
  TDirectory *dCan = (TDirectory*) fOut -> mkdir("candidate");
  TDirectory *dPar = (TDirectory*) fOut -> mkdir("particle");
  TDirectory *dEff = (TDirectory*) fOut -> mkdir("efficiency");
  TDirectory *dRes = (TDirectory*) fOut -> mkdir("resolution");
  cout << "    Made directories." << endl;

  // save histograms
  dDet            -> cd();
  hPhiTrg[0]      -> Write();
  hEtaTrg[0]      -> Write();
  hEtTrg[0]       -> Write();
  dCan            -> cd();
  hPhiTrg[1]      -> Write();
  hEtaTrg[1]      -> Write();
  hEtTrg[1]       -> Write();
  hDrTrg[0]       -> Write();
  hQtTrg[0]       -> Write();
  hDqTrg[0]       -> Write();
  dPar            -> cd();
  hPhiTrg[2]      -> Write();
  hEtaTrg[2]      -> Write();
  hEtTrg[2]       -> Write();
  hEtMatch[0]     -> Write();
  hEtMatch[1]     -> Write();
  hEtMatch[2]     -> Write();
  hDrTrg[1]       -> Write();
  hDrMatch[0]     -> Write();
  hDrMatch[1]     -> Write();
  hDrMatch[2]     -> Write();
  hQtTrg[1]       -> Write();
  hQtMatch[0]     -> Write();
  hQtMatch[1]     -> Write();
  hQtMatch[2]     -> Write();
  hDqTrg[1]       -> Write();
  hDqMatch[0]     -> Write();
  hDqMatch[1]     -> Write();
  hDqMatch[2]     -> Write();
  dEff            -> cd();
  hPhiForEff[0]   -> Write();
  hPhiForEff[1]   -> Write();
  hEtaForEff[0]   -> Write();
  hEtaForEff[1]   -> Write();
  hEtForEff[0]    -> Write();
  hEtForEff[1]    -> Write();
  hPhiEff         -> Write();
  hEtaEff         -> Write();
  hEtEff          -> Write();
  dRes            -> cd();
  hPhiMatchVsDet  -> Write();
  pPhiMatchVsDet  -> Write();
  hEtaMatchVsDet  -> Write();
  pEtaMatchVsDet  -> Write();
  hEtMatchVsDet   -> Write();
  pEtMatchVsDet   -> Write();
  hDrMatchVsEtDet -> Write();
  pDrMatchVsEtDet -> Write();
  hQtMatchVsEtDet -> Write();
  pQtMatchVsEtDet -> Write();
  hDqMatchVsEtDet -> Write();
  pDqMatchVsEtDet -> Write();
  cout << "    Saved histograms." << endl;

  // TEST [06.12.2020]
  fOut      -> cd();
  hNumMatch -> Write();


  // close files
  fOut -> cd();
  fOut -> Close();
  fPar -> cd();
  fPar -> Close();
  fDet -> cd();
  fDet -> Close();
  cout << "  Finished trigger energy resolution script!\n" << endl;

}

// End ------------------------------------------------------------------------
