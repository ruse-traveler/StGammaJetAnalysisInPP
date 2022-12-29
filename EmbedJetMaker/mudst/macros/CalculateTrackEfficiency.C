// 'CalculateTrackEfficiency.C'
// Derek Anderson
// 03.01.2019
//
// Use this to calculate the track
// efficiency as a function of pT,
// eta, and phi from the Run 9 dijet
// embedding sample.
//
// NOTE: this was created to be run
// on RCF (my old code was having
// trouble running on RCF).


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t  NLvls(2);
static const UInt_t  NCuts(2);
static const UInt_t  NHistQA(4);
static const UInt_t  NPtBins(48);
//static const UInt_t  NHadIds(16);
static const UInt_t  NHadIds(1);
static const UInt_t  NBadRuns(45);
static const UInt_t  NTrkMax(5000);
static const UInt_t  NTwrMax(5000);
static const UInt_t  NMatchMax(10);
static const UInt_t  NTrkSpecies(2);
static const TString SParTree("McTracks");
static const TString SDetTree("GfmtoDst_mu");
static const TString SInDefault("../../MuDstMatching/output/merged/pt25rff.matchWithMc.root");
static const TString SOutDefault("pp200r9pt25rff.idEnd0andIdVtxNot1butPtLessThan0andOnlyMatched.et920vz55pi0.d14m1y2020.root");

// global flags
static const Bool_t UseParTrg(true);
static const Bool_t UseDetTrg(false);
static const Bool_t UseOnlyRecoilTracks(true);
static const Bool_t UseSpecificTrackSpecies(false);



void CalculateTrackEfficiency(const TString sIn=SInDefault, const TString sOut=SOutDefault, const Bool_t isInBatchMode=false) {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning efficiency calculation..." << endl;

  // event parameters
  const Double_t rVtxMax(2.);
  const Double_t zVtxMax(55.);
  const UInt_t   badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};
  // trigger parameters
  const Double_t hTrgMax(0.9);
  const Double_t eTtrgMin(12.);
  const Double_t eTtrgMax(20.);
  const Double_t nSigTrgCut(3.);
  //const UInt_t   idHadTrg[NHadIds] = {8, 9, 11, 12, 14, 15, 19, 21, 23, 24, 27, 29, 31, 32, 45, 46};
  const UInt_t   idHadTrg[NHadIds] = {7};
  // track parameters
  const Int_t    minQaTruth(50);
  const UInt_t   nFitMin(15);
  const Double_t rFitMin(0.52);
  const Double_t dcaMax(3.);
  const Double_t hTrkMax(1.);
  const Double_t pTtrkMin(0.2);
  const Double_t pTtrkMax(30.);
  //const Double_t dFminAway(TMath::PiOver2());
  const Double_t dFminAway(TMath::Pi() - 1.3);
  //const Double_t dFmaxAway(3. * TMath::PiOver2());
  const Double_t dFmaxAway(TMath::Pi() + 1.3);
  const Double_t idTrkSpecies[NTrkSpecies] = {8, 9};


  // open files
  TFile *fInput  = new TFile(sIn.Data(), "read");
  TFile *fOutput = new TFile(sOut.Data(), "recreate");
  if (!fInput || !fOutput) {
    cerr << "PANIC: couldn't open input file and/or create output file!\n"
         << "       fInput = " << fInput << ", fOutput = " << fOutput
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input trees
  TTree *tPar;
  TTree *tDet;
  fInput -> GetObject(SParTree.Data(), tPar);
  fInput -> GetObject(SDetTree.Data(), tDet);
  if (!tPar || !tDet) {
    cerr << "PANIC: couldn't grab particle and/or detector level tree!\n"
         << "       tPar = " << tPar << ", tDet = " << tDet
         << endl;
    return;
  }
  cout << "    Grabbed input trees." << endl;


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
  // define miscellaneous histograms
  TH1D *hNumTrk[NLvls];
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
  // create particle histograms
  hPhiTrk[0][0] = new TH1D("hPhiBeforeQA_par", "", nPhi, f[0], f[1]);
  hPhiTrk[0][1] = new TH1D("hPhiAfterQA_par", "", nPhi, f[0], f[1]);
  hEtaTrk[0][0] = new TH1D("hEtaBeforeQA_par", "", nEta, h[0], h[1]);
  hEtaTrk[0][1] = new TH1D("hEtaAfterQA_par", "", nEta, h[0], h[1]);
  hPtTrk[0][0]  = new TH1D("hPtBeforeQA_par", "", nPtHistBins, pTbinEdges);
  hPtTrk[0][1]  = new TH1D("hPtAfterQA_par", "", nPtHistBins, pTbinEdges);
  // create detector histograms
  hPhiTrk[1][0] = new TH1D("hPhiBeforeQA_det", "", nPhi, f[0], f[1]);
  hPhiTrk[1][1] = new TH1D("hPhiAfterQA_det", "", nPhi, f[0], f[1]);
  hEtaTrk[1][0] = new TH1D("hEtaBeforeQA_det", "", nEta, h[0], h[1]);
  hEtaTrk[1][1] = new TH1D("hEtaAfterQA_det", "", nEta, h[0], h[1]);
  hPtTrk[1][0]  = new TH1D("hPtBeforeQA_det", "", nPtHistBins, pTbinEdges);
  hPtTrk[1][1]  = new TH1D("hPtAfterQA_det", "", nPtHistBins, pTbinEdges);
  // create matching histograms
  hPhiForEff[0] = new TH1D("hPhiForEff_par", "", nPhi, f[0], f[1]);
  hPhiForEff[1] = new TH1D("hPhiForEff_det", "", nPhi, f[0], f[1]);
  hEtaForEff[0] = new TH1D("hEtaForEff_par", "", nEta, h[0], h[1]);
  hEtaForEff[1] = new TH1D("hEtaForEff_det", "", nEta, h[0], h[1]);
  hPtForEff[0]  = new TH1D("hPtForEff_par", "", nPtHistBins, pTbinEdges);
  hPtForEff[1]  = new TH1D("hPtForEff_det", "", nPtHistBins, pTbinEdges);
  hPhiParVsDet  = new TH2D("hPhiParVsDet", "", nPhi, f[0], f[1], nPhi, f[0], f[1]);
  hEtaParVsDet  = new TH2D("hEtaParVsDet", "", nEta, h[0], h[1], nEta, h[0], h[1]);
  hPtParVsDet   = new TH2D("hPtParVsDet", "", nPt, pT[0], pT[1], nPt, pT[0], pT[1]);
  // create efficiency and resolution histograms
  hPhiEff       = new TH1D("hPhiEfficiency", "", nPhi, f[0], f[1]);
  hEtaEff       = new TH1D("hEtaEfficiency", "", nEta, h[0], h[1]);
  hPtEff        = new TH1D("hPtEfficiency", "", nPtHistBins, pTbinEdges);
  hPhiRes       = new TH1D("hPhiRes", "", nDf, dF[0], dF[1]);
  hEtaRes       = new TH1D("hEtaRes", "", nDh, dH[0], dH[1]);
  hPtRes        = new TH1D("hPtRes", "", nDpt, dPt[0], dPt[1]);
  hPhiResVsPhi  = new TH2D("hPhiResVsPhi", "", nPhi, f[0], f[1], nDf, dF[0], dF[1]);
  hEtaResVsEta  = new TH2D("hEtaResVsEta", "", nEta, h[0], h[1], nDh, dH[0], dH[1]);
  hPtResVsPt    = new TH2D("hPtResVsPt", "", nPt, pT[0], pT[1], nDpt, dPt[0], dPt[1]);
  // create miscellaneous histograms
  hNumTrk[0]    = new TH1D("hNumTrk_par", "After QA", nNum, num[0], num[1]);
  hNumTrk[1]    = new TH1D("hNumTrk_det", "After QA", nNum, num[0], num[1]);
  // create QA histograms
  hEtaQA[0]     = new TH1D("hEtaParQA", "", nEta, hQA[0], hQA[1]);
  hEtaQA[1]     = new TH1D("hEtaDetQA", "", nEta, hQA[0], hQA[1]);
  hEtaQA[2]     = new TH1D("hEtaMatQA", "", nEta, hQA[0], hQA[1]);
  hEtaQA[3]     = new TH1D("hEtaUnMatQA", "", nEta, hQA[0], hQA[1]);
  hPtQA[0]      = new TH1D("hPtParQA", "", nPtHistBins, pTbinEdges);
  hPtQA[1]      = new TH1D("hPtDetQA", "", nPtHistBins, pTbinEdges);
  hPtQA[2]      = new TH1D("hPtMatQA", "", nPtHistBins, pTbinEdges);
  hPtQA[3]      = new TH1D("hPtUnMatQA", "", nPtHistBins, pTbinEdges);
  hNfitQA[0]    = new TH1D("hNfitDetQA", "", nNfit, nFit[0], nFit[1]);
  hNfitQA[1]    = new TH1D("hNfitMatQA", "", nNfit, nFit[0], nFit[1]);
  hNfitQA[2]    = new TH1D("hNfitUnMatQA", "", nNfit, nFit[0], nFit[1]);
  hRfitQA[0]    = new TH1D("hRfitDetQA", "", nRfit, rFit[0], rFit[1]);
  hRfitQA[1]    = new TH1D("hRfitMatQA", "", nRfit, rFit[0], rFit[1]);
  hRfitQA[2]    = new TH1D("hRfitUnMatQA", "", nRfit, rFit[0], rFit[1]);
  hDcaQA[0]     = new TH1D("hDcaDetQA", "", nDca, dca[0], dca[1]);
  hDcaQA[1]     = new TH1D("hDcaMatQA", "", nDca, dca[0], dca[1]);
  hDcaQA[2]     = new TH1D("hDcaUnMatQA", "", nDca, dca[0], dca[1]);
  hPtFrac[0]    = new TH1D("hPtFracDetQA", "", nPtFrac, pTfrac[0], pTfrac[1]);
  hPtFrac[1]    = new TH1D("hPtFracMatQA", "", nPtFrac, pTfrac[0], pTfrac[1]);
  hPtFrac[2]    = new TH1D("hPtFracUnMatQA", "", nPtFrac, pTfrac[0], pTfrac[1]);
  hIdTruth[0]   = new TH1D("hIdTruthDetQA", "", nIdTruth, idTruth[0], idTruth[1]);
  hIdTruth[1]   = new TH1D("hIdTruthMatQA", "", nIdTruth, idTruth[0], idTruth[1]);
  hIdTruth[2]   = new TH1D("hIdTruthUnMatQA", "", nIdTruth, idTruth[0], idTruth[1]);
  hQaTruth[0]   = new TH1D("hQaTruthDetQA", "", nQaTruth, qaTruth[0], qaTruth[1]);
  hQaTruth[1]   = new TH1D("hQaTruthMatQA", "", nQaTruth, qaTruth[0], qaTruth[1]);
  hQaTruth[2]   = new TH1D("hQaTruthUnMatQA", "", nQaTruth, qaTruth[0], qaTruth[1]);
  hNcommon[0]   = new TH1D("hNcommonDetQA", "", nNcommon, nCommon[0], nCommon[1]);
  hNcommon[1]   = new TH1D("hNcommonMatQA", "", nNcommon, nCommon[0], nCommon[1]);
  hNcommon[2]   = new TH1D("hNcommonUnMatQA", "", nNcommon, nCommon[0], nCommon[1]);
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
  hPhiParVsDet -> Sumw2();
  hEtaParVsDet -> Sumw2();
  hPtParVsDet  -> Sumw2();
  hPhiEff      -> Sumw2();
  hEtaEff      -> Sumw2();
  hPtEff       -> Sumw2();
  hPhiRes      -> Sumw2();
  hEtaRes      -> Sumw2();
  hPtRes       -> Sumw2();
  hPhiResVsPhi -> Sumw2();
  hEtaResVsEta -> Sumw2();
  hPtResVsPt   -> Sumw2();
  cout << "    Created histograms." << endl;

  // make profiles
  TProfile *pPhiParVsDet = new TProfile("pPhiParVsDet", "", nPhi, f[0], f[1], "S");
  TProfile *pEtaParVsDet = new TProfile("pEtaParVsDet", "", nEta, h[0], h[1], "S");
  TProfile *pPtParVsDet  = new TProfile("pPtParVsDet", "", nPt, pT[0], pT[1], "S");
  TProfile *pPhiResVsPhi = new TProfile("pPhiResVsPhi", "", nPhi, f[0], f[1], "S");
  TProfile *pEtaResVsEta = new TProfile("pEtaResVsEta", "", nEta, h[0], h[1], "S");
  TProfile *pPtResVsPt   = new TProfile("pPtResVsPt", "", nPt, pT[0], pT[1], "S");
  cout << "    Created profiles." << endl;


  // begin event loops
  const UInt_t nParEvts    = tPar -> GetEntriesFast();
  const UInt_t nDetEvts    = tDet -> GetEntriesFast();
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
    parBytes   = tPar -> GetEntry(iEvt);
    detBytes   = tDet -> GetEntry(iEvt);
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
    if (isInBatchMode) {
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
    Bool_t   isGoodTrgPar(true);
    Bool_t   foundTrgPar(false);
    if (UseParTrg) {
      for (UInt_t iTrgMC = 0; iTrgMC < nTrkPar; iTrgMC++) {

        // particle trigger info
        const Int_t    pidTrgMc  = mcIdGeant -> at(iTrgMC);
        const Int_t    idVxTrgMc = mcIdVx    -> at(iTrgMC);
        const Int_t    idVeTrgMc = mcIdVxEnd -> at(iTrgMC);
        const Double_t qTrgMc    = mcCharge  -> at(iTrgMC);
        const Double_t fTrgMc    = mcPhi     -> at(iTrgMC);
        const Double_t hTrgMc    = mcEta     -> at(iTrgMC);
        const Double_t pTtrgMc   = mcPt      -> at(iTrgMC);

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
        const Bool_t isGoodStateTrgPar = (idVxTrgMc == 1);
        //const Bool_t isChargedTrgPar   = (qTrgMc != 0.);
        const Bool_t isChargedTrgPar   = (qTrgMc == 0.);
        const Bool_t isInTrgEtaCutPar  = (TMath::Abs(hTrgMc) < hTrgMax);
        const Bool_t isInTrgEtCutPar   = ((pTtrgMc > eTtrgMin) && (pTtrgMc < eTtrgMax));
        if (isHadronPar && isGoodStateTrgPar && isChargedTrgPar && isInTrgEtaCutPar && isInTrgEtCutPar) {
          iTrgPar     = iTrgMC;
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
    if (UseDetTrg) {
      for (UInt_t iTrgRE = 0; iTrgRE < nTrkDet; iTrgRE++) {

        // detector trigger info
        const UInt_t   nFitTrgRE  = PrimaryTrackArray_nHitsFit[iTrgRE];
        const UInt_t   nPossTrgRE = PrimaryTrackArray_nHitsPoss[iTrgRE];
        const Double_t rFitTrgRE  = (Double_t) nFitTrgRE / (Double_t) nPossTrgRE;
        const Double_t dcaTrgRE   = PrimaryTrackArray_dcag[iTrgRE];
        const Double_t chrgTrgRE  = PrimaryTrackArray_charge[iTrgRE];
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


    // particle track loop
    UInt_t nPar(0);
    for (UInt_t iTrkMC = 0; iTrkMC < nTrkPar; iTrkMC++) {

      // particle track info
      const Int_t    mcId  = mcIdTrk   -> at(iTrkMC);
      const Int_t    mcPID = mcIdGeant -> at(iTrkMC);
      const Int_t    mcIV  = mcIdVx    -> at(iTrkMC);
      const Int_t    mcIE  = mcIdVxEnd -> at(iTrkMC);
      const Double_t qMC   = mcCharge  -> at(iTrkMC);
      const Double_t fMC   = mcPhi     -> at(iTrkMC);
      const Double_t hMC   = mcEta     -> at(iTrkMC);
      const Double_t pTmc  = mcPt      -> at(iTrkMC);

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
      if (UseSpecificTrackSpecies && !isRelevantSpecies) continue;

      // consider only recoil tracks (if need be)
      const Bool_t isAwaySidePar = ((dFmc > dFminAway) && (dFmc < dFmaxAway));
      if (UseOnlyRecoilTracks && !isAwaySidePar) continue;

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
        }
        nPar++;

        // match to detector tracks
        for (UInt_t iTrkMat = 0; iTrkMat < nTrkDet; iTrkMat++) {

          // detector track info
          const Int_t    idMat     = (Int_t) PrimaryTrackArray_tofSigElectron[iTrkMat];
          const Int_t    qaMat     = (Int_t) PrimaryTrackArray_tofSigPion[iTrkMat];
          const UInt_t   nFitMat   = PrimaryTrackArray_nHitsFit[iTrkMat];
          const UInt_t   nPossMat  = PrimaryTrackArray_nHitsPoss[iTrkMat];
          const UInt_t   nComMat   = (UInt_t) (nFitMat * ((Double_t) qaMat / 100.));
          const Double_t rFitMat   = (Double_t) nFitMat / (Double_t) nPossMat;
          const Double_t dcaMat    = PrimaryTrackArray_dcag[iTrkMat];
          const Double_t qMat      = PrimaryTrackArray_charge[iTrkMat];
          const Double_t fMat      = PrimaryTrackArray_phi[iTrkMat];
          const Double_t hMat      = PrimaryTrackArray_eta[iTrkMat];
          const Double_t pTmat     = PrimaryTrackArray_pT[iTrkMat];
          const Double_t fDifMat   = fMC - fMat;
          const Double_t hDifMat   = hMC - hMat;
          const Double_t pTdifMat  = pTmc - pTmat;
          const Double_t pTfracMat = pTmat / pTmc;

          Double_t dFmat = fMat - fTrgPar;
          if (dFmat < (-1. * TMath::PiOver2())) dFmat += TMath::TwoPi();
          if (dFmat > (3. * TMath::PiOver2()))  dFmat -= TMath::TwoPi();

          // consider only recoil tracks (if need be)
          const Bool_t isAwaySideMat = ((dFmat > dFminAway) && (dFmat < dFmaxAway));
          if (UseOnlyRecoilTracks && !isAwaySideMat) continue;

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

            // efficiency & resolution histograms
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
          }
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
      if (isGoodChrg && isNotTrgPar) {
        hPhiTrk[0][0] -> Fill(fMC);
        hEtaTrk[0][0] -> Fill(hMC);
        hPtTrk[0][0]  -> Fill(pTmc);
      }
    }  // end particle track loop


    // detector track loop
    UInt_t nDet(0);
    for (UInt_t iTrkDet = 0; iTrkDet < nTrkDet; iTrkDet++) {

      // detector track info
      const UInt_t   nFitDet  = PrimaryTrackArray_nHitsFit[iTrkDet];
      const UInt_t   nPossDet = PrimaryTrackArray_nHitsPoss[iTrkDet];
      const Double_t rFitDet  = (Double_t) nFitDet / (Double_t) nPossDet;
      const Double_t dcaDet   = PrimaryTrackArray_dcag[iTrkDet];
      const Double_t qDet     = PrimaryTrackArray_charge[iTrkDet];
      const Double_t fDet     = PrimaryTrackArray_phi[iTrkDet];
      const Double_t hDet     = PrimaryTrackArray_eta[iTrkDet];
      const Double_t pTdet    = PrimaryTrackArray_pT[iTrkDet];

      Double_t dFdet = fDet - fTrgPar;
      if (dFdet < (-1. * TMath::PiOver2())) dFdet += TMath::TwoPi();
      if (dFdet > (3. * TMath::PiOver2()))  dFdet -= TMath::TwoPi();

      // consider only recoil tracks (if need be)
      const Bool_t isAwaySideDet = ((dFdet > dFminAway) && (dFdet < dFmaxAway));
      if (UseOnlyRecoilTracks && !isAwaySideDet) continue;

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

        // fill detector histograms (after QA)
        hPhiTrk[1][1] -> Fill(fDet);
        hEtaTrk[1][1] -> Fill(hDet);
        hPtTrk[1][1]  -> Fill(pTdet);
        nDet++;

      }  // end detector track cuts
    }  // end detector track loop

    // fill misc. histograms
    hNumTrk[0] -> Fill(nPar);
    hNumTrk[1] -> Fill(nDet);

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
    dLvls[iLvl] = (TDirectory*) fOutput -> mkdir(sLvlDir[iLvl].Data());
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
  dLvls[NLvls]     = (TDirectory*) fOutput -> mkdir(sLvlDir[NLvls].Data());
  dLvls[NLvls + 1] = (TDirectory*) fOutput -> mkdir(sLvlDir[NLvls + 1].Data());
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
  fOutput -> cd();
  hPhiEff -> Write();
  hEtaEff -> Write();
  hPtEff  -> Write();
  cout << "    Made directories and saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Finished efficiency calculation!" << endl;

}

// End ------------------------------------------------------------------------
