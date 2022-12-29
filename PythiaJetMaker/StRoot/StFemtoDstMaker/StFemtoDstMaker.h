// 'StFemtoDstMaker.h'
// Derek Anderson
// 08.04.2016
//
// This class takes output from Pythia, finds jets, and stores in them in a
// compact tree.
//
// Last updated: 04.17.2018

#ifndef StFemtoDstMaker_h
#define StFemtoDstMaker_h

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include "StFemtoDst.h"
#include "StFemtoDstLinkDef.h"

const Int_t    NTmax    = 1000;
const Double_t PionMass = 0.140;
const Double_t HighPt1  = 3.;
const Double_t HighPt2  = 7.;
const Double_t EffPtMin = 0.;
const Double_t EffPtMax = 30.;


class StFemtoDstMaker {

public:

  Int_t      fCurrent;  // current Tree number in a TChain
  TChain     *fChain;   // pointer to the analyzed TTree or TChain
  //TTree      *fChain;   // quick fix [Derek, 11.28.2016]
  TTree      *fDst;     // output tree
  TFile      *oFile;    // output file
  StFemtoDst femto;     // FemtoDst Object

  // Declaration of leaf types
  Int_t    Events_num;
  Int_t    Events_Process;
  Int_t    Events_refmult;
  Int_t    Events_refPos;
  Int_t    Events_refNeg;
  Int_t    Events_runId;
  Double_t Events_MagF;
  Int_t    Events_nVertex;
  Float_t  Events_rankV;
  Double_t Events_pvpdz;
  Double_t Events_BBCz;
  Double_t Events_ZDCz;
  Double_t Events_primVx;
  Double_t Events_primVy;
  Double_t Events_primVz;
  Int_t    Events_TrigId;
  Int_t    Events_TrigStat;
  Int_t    Events_TrigIndex;
  Int_t    Events_PartonId_NS;
  Int_t    Events_PartonId_AS;
  Int_t    Events_PartonStat_NS;
  Int_t    Events_PartonStat_AS;
  Float_t  Events_PartonEta_NS;
  Float_t  Events_PartonEta_AS;
  Float_t  Events_PartonPhi_NS;
  Float_t  Events_PartonPhi_AS;
  Float_t  Events_PartonE_NS;
  Float_t  Events_PartonE_AS;
  Float_t  Events_PartonEt_NS;
  Float_t  Events_PartonEt_AS;
  Float_t  Events_tsp;
  Float_t  Events_iso;
  Int_t    Events_Twr_didT;
  Int_t    Events_Twr_adc11;
  Float_t  Events_Twr_eneT0;
  Float_t  Events_Twr_eT;
  Float_t  Events_Twr_ENET0;
  Float_t  Events_Twr_phT;
  Float_t  Events_Twr_PTower;
  Float_t  Events_Twr_pidTower;
  Float_t  Events_Twr_moduleT;
  Float_t  Events_Clust_EneT0;
  Float_t  Events_Clust_EneR;
  Float_t  Events_Clust_EneH;
  Float_t  Events_Clust_rapv1;
  Float_t  Events_Clust_etav1;
  Float_t  Events_Clust_phiv1;
  Float_t  Events_Estrp_en01;
  Float_t  Events_Estrp_en02;
  Float_t  Events_Estrp_en03;
  Float_t  Events_Estrp_en0;
  Float_t  Events_Estrp_en1;
  Float_t  Events_Estrp_en2;
  Float_t  Events_Estrp_en3;
  Float_t  Events_Estrp_en4;
  Float_t  Events_Estrp_en5;
  Float_t  Events_Estrp_en6;
  Float_t  Events_Estrp_en7;
  Float_t  Events_Estrp_en8;
  Float_t  Events_Estrp_en9;
  Float_t  Events_Estrp_en10;
  Float_t  Events_Estrp_en11;
  Float_t  Events_Estrp_en12;
  Float_t  Events_Estrp_en13;
  Float_t  Events_Estrp_en14;
  Float_t  Events_Estrp_en15;
  Int_t    Events_Twr_didE;
  Float_t  Events_Pstrip_enp01;
  Float_t  Events_Pstrip_enp02;
  Float_t  Events_Pstrip_enp03;
  Float_t  Events_Pstrip_enp0;
  Float_t  Events_Pstrip_enp1;
  Float_t  Events_Pstrip_enp2;
  Float_t  Events_Pstrip_enp3;
  Float_t  Events_Pstrip_enp4;
  Float_t  Events_Pstrip_enp5;
  Float_t  Events_Pstrip_enp6;
  Float_t  Events_Pstrip_enp7;
  Float_t  Events_Pstrip_enp8;
  Float_t  Events_Pstrip_enp9;
  Float_t  Events_Pstrip_enp10;
  Float_t  Events_Pstrip_enp11;
  Float_t  Events_Pstrip_enp12;
  Float_t  Events_Pstrip_enp13;
  Float_t  Events_Pstrip_enp14;
  Float_t  Events_Pstrip_enp15;
  Float_t  Events_clust_Ennq1;
  Float_t  Events_clust_Ennq20;
  Float_t  Events_clust_Ennq19;
  Float_t  Events_clust_Enpq1;
  Float_t  Events_clust_Enpq20;
  Float_t  Events_clust_Enpq19;
  Float_t  Events_clust_Enpq21;
  Int_t    Events_noOfprimaryTrks;
  Float_t  pTracks_pT[NTmax];
  Float_t  pTracks_px[NTmax];
  Float_t  pTracks_py[NTmax];
  Float_t  pTracks_pz[NTmax];
  Float_t  pTracks_Eta[NTmax];
  Float_t  pTracks_Phi[NTmax];
  Float_t  pTracks_dEdx[NTmax];
  Float_t  pTracks_chrg[NTmax];
  Float_t  pTracks_gdca[NTmax];
  Int_t    pTracks_Fp[NTmax];
  Int_t    pTracks_Ppo[NTmax];
  Float_t  pTracks_nSigPi[NTmax];
  Float_t  pTracks_nSigK[NTmax];
  Float_t  pTracks_nSigP[NTmax];
  Float_t  pTracks_nSigE[NTmax];

  // List of branches
  TBranch *b_Events_num;
  TBranch *b_Events_Process;
  TBranch *b_Events_refmult;
  TBranch *b_Events_refPos;
  TBranch *b_Events_refNeg;
  TBranch *b_Events_runId;
  TBranch *b_Events_MagF;
  TBranch *b_Events_nVertex;
  TBranch *b_Events_rankV;
  TBranch *b_Events_pvpdz;
  TBranch *b_Events_BBCz;
  TBranch *b_Events_ZDCz;
  TBranch *b_Events_primVx;
  TBranch *b_Events_primVy;
  TBranch *b_Events_primVz;
  TBranch *b_Events_TrigId;
  TBranch *b_Events_TrigStat;
  TBranch *b_Events_TrigIndex;
  TBranch *b_Events_PartonId_NS;
  TBranch *b_Events_PartonId_AS;
  TBranch *b_Events_PartonStat_NS;
  TBranch *b_Events_PartonStat_AS;
  TBranch *b_Events_PartonEta_NS;
  TBranch *b_Events_PartonEta_AS;
  TBranch *b_Events_PartonPhi_NS;
  TBranch *b_Events_PartonPhi_AS;
  TBranch *b_Events_PartonE_NS;
  TBranch *b_Events_PartonE_AS;
  TBranch *b_Events_PartonEt_NS;
  TBranch *b_Events_PartonEt_AS;
  TBranch *b_Events_tsp;
  TBranch *b_Events_iso;
  TBranch *b_Events_Twr_didT;
  TBranch *b_Events_Twr_adc11;
  TBranch *b_Events_Twr_eneT0;
  TBranch *b_Events_Twr_eT;
  TBranch *b_Events_Twr_ENET0;
  TBranch *b_Events_Twr_phT;
  TBranch *b_Events_Twr_PTower;
  TBranch *b_Events_Twr_pidTower;
  TBranch *b_Events_Twr_moduleT;
  TBranch *b_Events_Clust_EneT0;
  TBranch *b_Events_Clust_EneR;
  TBranch *b_Events_Clust_EneH;
  TBranch *b_Events_Clust_rapv1;
  TBranch *b_Events_Clust_etav1;
  TBranch *b_Events_Clust_phiv1;
  TBranch *b_Events_Estrp_en01;
  TBranch *b_Events_Estrp_en02;
  TBranch *b_Events_Estrp_en03;
  TBranch *b_Events_Estrp_en0;
  TBranch *b_Events_Estrp_en1;
  TBranch *b_Events_Estrp_en2;
  TBranch *b_Events_Estrp_en3;
  TBranch *b_Events_Estrp_en4;
  TBranch *b_Events_Estrp_en5;
  TBranch *b_Events_Estrp_en6;
  TBranch *b_Events_Estrp_en7;
  TBranch *b_Events_Estrp_en8;
  TBranch *b_Events_Estrp_en9;
  TBranch *b_Events_Estrp_en10;
  TBranch *b_Events_Estrp_en11;
  TBranch *b_Events_Estrp_en12;
  TBranch *b_Events_Estrp_en13;
  TBranch *b_Events_Estrp_en14;
  TBranch *b_Events_Estrp_en15;
  TBranch *b_Events_Twr_didE;
  TBranch *b_Events_Pstrip_enp01;
  TBranch *b_Events_Pstrip_enp02;
  TBranch *b_Events_Pstrip_enp03;
  TBranch *b_Events_Pstrip_enp0;
  TBranch *b_Events_Pstrip_enp1;
  TBranch *b_Events_Pstrip_enp2;
  TBranch *b_Events_Pstrip_enp3;
  TBranch *b_Events_Pstrip_enp4;
  TBranch *b_Events_Pstrip_enp5;
  TBranch *b_Events_Pstrip_enp6;
  TBranch *b_Events_Pstrip_enp7;
  TBranch *b_Events_Pstrip_enp8;
  TBranch *b_Events_Pstrip_enp9;
  TBranch *b_Events_Pstrip_enp10;
  TBranch *b_Events_Pstrip_enp11;
  TBranch *b_Events_Pstrip_enp12;
  TBranch *b_Events_Pstrip_enp13;
  TBranch *b_Events_Pstrip_enp14;
  TBranch *b_Events_Pstrip_enp15;
  TBranch *b_Events_clust_Ennq1;
  TBranch *b_Events_clust_Ennq20;
  TBranch *b_Events_clust_Ennq19;
  TBranch *b_Events_clust_Enpq1;
  TBranch *b_Events_clust_Enpq20;
  TBranch *b_Events_clust_Enpq19;
  TBranch *b_Events_clust_Enpq21;
  TBranch *b_Events_noOfprimaryTrks;
  TBranch *b_pTracks_pT;
  TBranch *b_pTracks_px;
  TBranch *b_pTracks_py;
  TBranch *b_pTracks_pz;
  TBranch *b_pTracks_Eta;
  TBranch *b_pTracks_Phi;
  TBranch *b_pTracks_dEdx;
  TBranch *b_pTracks_chrg;
  TBranch *b_pTracks_gdca;
  TBranch *b_pTracks_Fp;
  TBranch *b_pTracks_Ppo;
  TBranch *b_pTracks_nSigPi;
  TBranch *b_pTracks_nSigK;
  TBranch *b_pTracks_nSigP;
  TBranch *b_pTracks_nSigE;

  // QA event histograms
  TH1D *hRefmult;
  TH1D *hTSP;
  TH1D *hEtTrg;
  TH1D *hFtrg;
  TH1D *hHtrg;
  TH1D *hRho;
  TH1D *hNjet;
  TH1D *hNjetG17;
  TH1D *hNjetL17;
  TH1D *hNtrk;
  TH1D *hNtrkG1;
  TH1D *hNtrkG2;
  // QA track histograms
  TH1D *hPtTrk;
  TH1D *hPxTrk;
  TH1D *hPyTrk;
  TH1D *hPzTrk;
  TH1D *hDfTrk;
  TH1D *hDfTrkG1;
  TH1D *hDfTrkG2;
  TH1D *hHtrk;
  TH1D *hPtEffPar;
  TH1D *hPtEffDet;
  TH1D *hPtEffCheck;
  TH1D *hPtRes1D;
  TH2D *hPtRes2D;
  TH2D *hPtVsDfTrk;
  // QA jet histograms
  TH1D *hPtReco;
  TH1D *hPtCorr;
  TH1D *hPtRE;
  TH1D *hPtUE;
  TH1D *hPtSub;
  TH1D *hDfJet;
  TH1D *hHjet;
  TH1D *hAjet;
  // QA cst. histograms
  TH1D *hPtCst;
  TH1D *hDfCst;
  TH1D *hHcst;

  // TEST [04.29.2020]
  TH1D *hDfTrk_pt021;
  TH1D *hDfTrk_pt12;
  TH1D *hDfTrk_pt25;
  TH1D *hDfTrk_pt530;

  StFemtoDstMaker(const TString iName="", const TString oName="", const Int_t lvl=0, const Int_t typ=0, TChain *chain=0);
  //StFemtoDstMaker(const TString iName="", const TString oName="", const Int_t lvl=0, const Int_t typ=0, TTree *chain=0);  // quick fix [Derek, 11.28.2016]
  virtual ~StFemtoDstMaker();

  // public methods
  virtual void SetEfficiency(const TString sEffFile="", const TString sEffHist="", const Bool_t isHist=false, const Bool_t isFunc=false, const Bool_t doHighPtSmooth=false, const Bool_t doEffAdjust=false, const Float_t effNudge=0., const Float_t amplitude=0.94, const Float_t sigma=18.5);
  virtual void Init(const Int_t nRM=1, const Int_t tID=22, const Double_t r=0.5, const Double_t a=0., const Double_t pMin=0.2, const Double_t pMax=20., const Double_t q=0., const Double_t eMin=8., const Double_t eMax=20., const Double_t hTrkMax=1., const Double_t hTrgMax=1.);
  virtual void Make(const Int_t nTrgs=-1, const Int_t StartEvt=0, const Int_t StopEvt=-1);
  virtual void Finish();


private:

  Int_t    level;
  Int_t    type;
  Int_t    Nrm;
  Int_t    trigID;
  Bool_t   doEffSys;
  Bool_t   doEffSmooth;
  Bool_t   effIsHist;
  Bool_t   effIsFunc;
  Float_t  effAdjust;
  Float_t  effAmplitude;
  Float_t  effSigma;
  Double_t Rjet;
  Double_t Amin;
  Double_t pTmin;
  Double_t pTmax;
  Double_t qMin;
  Double_t eTmin;
  Double_t eTmax;
  Double_t EtaTrkMax;
  Double_t EtaTrgMax;

  // eff. and res. histograms
  TH1D *hPtEff;
  TF1  *fPtEff;

  // private methods
  virtual void     InitTree(TChain *chain);
  //virtual void     InitTree(TTree *chain);  // quick fix [Derek, 11.28.2016]
  virtual void     Show(const Long64_t entry = -1);
  virtual Int_t    Cut(const Long64_t entry);
  virtual Int_t    GetEntry(const Long64_t entry);
  virtual Bool_t   Notify();
  virtual Long64_t LoadTree(const Long64_t entry);
  virtual Double_t ApplyDetectorResponse(const Double_t pTpar);


  ClassDef(StFemtoDstMaker, 1)

};


#endif
#ifdef StFemtoDstMaker_cxx


StFemtoDstMaker::StFemtoDstMaker(const TString iName, const TString oName, const Int_t lvl, const Int_t typ, TChain *chain) : fChain(0) {
//StFemtoDstMaker::StFemtoDstMaker(const TString iName, const TString oName, const Int_t lvl, const Int_t typ, TTree *chain) : fChain(0) {


  level = lvl;
  type  = typ;

  if (chain == 0) {
    TFile *f = (TFile*) gROOT -> GetListOfFiles() -> FindObject(iName);
    if (!f || !f->IsOpen())
      f = new TFile(iName);

    // detector response simulated at jet-finding stage
    if (level == 0)
      f -> GetObject("ParTree", chain);
    else
      f -> GetObject("ParTree", chain);
  }
  InitTree(chain);

  TString jetType;
  if (type == 0)
    jetType = "charged";
  else
    jetType = "full";


  // set output file, announce files
  oFile = new TFile(oName, "recreate");
  if (level == 0) {
    cout << "\n  Particle-level maker created! Reading from / writing to:\n"
         << "    in  -- '" << iName << "'\n"
         << "    out -- '" << oName << "'\n"
         << "  Making " << jetType << " jets..."
         << endl;
  }
  else{
    cout << "\n  Detector-level maker created! Reading from / writing to:\n"
         << "    in  -- '" << iName << "'\n"
         << "    out -- '" << oName << "'\n"
         << "  Making " << jetType << " jets..."
         << endl;
  }

}  // end 'StFemtoDstMaker(TTree*, TString, Int_t)'



StFemtoDstMaker::~StFemtoDstMaker() {

  if (!fChain) return;
  delete fChain -> GetCurrentFile();

}  // end '~StFemtoDstMaker()'



void StFemtoDstMaker::InitTree(TChain *chain) {
//void StFemtoDstMaker::InitTree(TTree *chain) {

   // Set branch addresses and branch pointers
  if (!chain) return;
  fChain   = chain;
  fCurrent = -1;
  fChain -> SetMakeClass(1);
  fChain -> SetBranchAddress("Events_num", &Events_num, &b_Events_num);
  fChain -> SetBranchAddress("Events_Process", &Events_Process, &b_Events_Process);
  fChain -> SetBranchAddress("Events_refmult", &Events_refmult, &b_Events_refmult);
  fChain -> SetBranchAddress("Events_refPos", &Events_refPos, &b_Events_refPos);
  fChain -> SetBranchAddress("Events_refNeg", &Events_refNeg, &b_Events_refNeg);
  fChain -> SetBranchAddress("Events_runId", &Events_runId, &b_Events_runId);
  fChain -> SetBranchAddress("Events_MagF", &Events_MagF, &b_Events_MagF);
  fChain -> SetBranchAddress("Events_nVertex", &Events_nVertex, &b_Events_nVertex);
  fChain -> SetBranchAddress("Events_rankV", &Events_rankV, &b_Events_rankV);
  fChain -> SetBranchAddress("Events_pvpdz", &Events_pvpdz, &b_Events_pvpdz);
  fChain -> SetBranchAddress("Events_BBCz", &Events_BBCz, &b_Events_BBCz);
  fChain -> SetBranchAddress("Events_ZDCz", &Events_ZDCz, &b_Events_ZDCz);
  fChain -> SetBranchAddress("Events_primVx", &Events_primVx, &b_Events_primVx);
  fChain -> SetBranchAddress("Events_primVy", &Events_primVy, &b_Events_primVy);
  fChain -> SetBranchAddress("Events_primVz", &Events_primVz, &b_Events_primVz);
  fChain -> SetBranchAddress("Events_TrigId", &Events_TrigId, &b_Events_TrigId);
  fChain -> SetBranchAddress("Events_TrigStat", &Events_TrigStat, &b_Events_TrigStat);
  fChain -> SetBranchAddress("Events_TrigIndex", &Events_TrigIndex, &b_Events_TrigIndex);
  fChain -> SetBranchAddress("Events_PartonId_NS", &Events_PartonId_NS, &b_Events_PartonId_NS);
  fChain -> SetBranchAddress("Events_PartonId_AS", &Events_PartonId_AS, &b_Events_PartonId_AS);
  fChain -> SetBranchAddress("Events_PartonStat_NS", &Events_PartonStat_NS, &b_Events_PartonStat_NS);
  fChain -> SetBranchAddress("Events_PartonStat_AS", &Events_PartonStat_AS, &b_Events_PartonStat_AS);
  fChain -> SetBranchAddress("Events_PartonEta_NS", &Events_PartonEta_NS, &b_Events_PartonEta_NS);
  fChain -> SetBranchAddress("Events_PartonEta_AS", &Events_PartonEta_AS, &b_Events_PartonEta_AS);
  fChain -> SetBranchAddress("Events_PartonPhi_NS", &Events_PartonPhi_NS, &b_Events_PartonPhi_NS);
  fChain -> SetBranchAddress("Events_PartonPhi_AS", &Events_PartonPhi_AS, &b_Events_PartonPhi_AS);
  fChain -> SetBranchAddress("Events_PartonE_NS", &Events_PartonE_NS, &b_Events_PartonE_NS);
  fChain -> SetBranchAddress("Events_PartonE_AS", &Events_PartonE_AS, &b_Events_PartonE_AS);
  fChain -> SetBranchAddress("Events_PartonEt_NS", &Events_PartonEt_NS, &b_Events_PartonEt_NS);
  fChain -> SetBranchAddress("Events_PartonEt_AS", &Events_PartonEt_AS, &b_Events_PartonEt_AS);
  fChain -> SetBranchAddress("Events_tsp", &Events_tsp, &b_Events_tsp);
  fChain -> SetBranchAddress("Events_iso", &Events_iso, &b_Events_iso);
  fChain -> SetBranchAddress("Events_Twr_didT", &Events_Twr_didT, &b_Events_Twr_didT);
  fChain -> SetBranchAddress("Events_Twr_adc11", &Events_Twr_adc11, &b_Events_Twr_adc11);
  fChain -> SetBranchAddress("Events_Twr_eneT0", &Events_Twr_eneT0, &b_Events_Twr_eneT0);
  fChain -> SetBranchAddress("Events_Twr_eT", &Events_Twr_eT, &b_Events_Twr_eT);
  fChain -> SetBranchAddress("Events_Twr_ENET0", &Events_Twr_ENET0, &b_Events_Twr_ENET0);
  fChain -> SetBranchAddress("Events_Twr_phT", &Events_Twr_phT, &b_Events_Twr_phT);
  fChain -> SetBranchAddress("Events_Twr_PTower", &Events_Twr_PTower, &b_Events_Twr_PTower);
  fChain -> SetBranchAddress("Events_Twr_pidTower", &Events_Twr_pidTower, &b_Events_Twr_pidTower);
  fChain -> SetBranchAddress("Events_Twr_moduleT", &Events_Twr_moduleT, &b_Events_Twr_moduleT);
  fChain -> SetBranchAddress("Events_Clust_EneT0", &Events_Clust_EneT0, &b_Events_Clust_EneT0);
  fChain -> SetBranchAddress("Events_Clust_EneR", &Events_Clust_EneR, &b_Events_Clust_EneR);
  fChain -> SetBranchAddress("Events_Clust_EneH", &Events_Clust_EneH, &b_Events_Clust_EneH);
  fChain -> SetBranchAddress("Events_Clust_rapv1", &Events_Clust_rapv1, &b_Events_Clust_rapv1);
  fChain -> SetBranchAddress("Events_Clust_etav1", &Events_Clust_etav1, &b_Events_Clust_etav1);
  fChain -> SetBranchAddress("Events_Clust_phiv1", &Events_Clust_phiv1, &b_Events_Clust_phiv1);
  fChain -> SetBranchAddress("Events_Estrp_en01", &Events_Estrp_en01, &b_Events_Estrp_en01);
  fChain -> SetBranchAddress("Events_Estrp_en02", &Events_Estrp_en02, &b_Events_Estrp_en02);
  fChain -> SetBranchAddress("Events_Estrp_en03", &Events_Estrp_en03, &b_Events_Estrp_en03);
  fChain -> SetBranchAddress("Events_Estrp_en0", &Events_Estrp_en0, &b_Events_Estrp_en0);
  fChain -> SetBranchAddress("Events_Estrp_en1", &Events_Estrp_en1, &b_Events_Estrp_en1);
  fChain -> SetBranchAddress("Events_Estrp_en2", &Events_Estrp_en2, &b_Events_Estrp_en2);
  fChain -> SetBranchAddress("Events_Estrp_en3", &Events_Estrp_en3, &b_Events_Estrp_en3);
  fChain -> SetBranchAddress("Events_Estrp_en4", &Events_Estrp_en4, &b_Events_Estrp_en4);
  fChain -> SetBranchAddress("Events_Estrp_en5", &Events_Estrp_en5, &b_Events_Estrp_en5);
  fChain -> SetBranchAddress("Events_Estrp_en6", &Events_Estrp_en6, &b_Events_Estrp_en6);
  fChain -> SetBranchAddress("Events_Estrp_en7", &Events_Estrp_en7, &b_Events_Estrp_en7);
  fChain -> SetBranchAddress("Events_Estrp_en8", &Events_Estrp_en8, &b_Events_Estrp_en8);
  fChain -> SetBranchAddress("Events_Estrp_en9", &Events_Estrp_en9, &b_Events_Estrp_en9);
  fChain -> SetBranchAddress("Events_Estrp_en10", &Events_Estrp_en10, &b_Events_Estrp_en10);
  fChain -> SetBranchAddress("Events_Estrp_en11", &Events_Estrp_en11, &b_Events_Estrp_en11);
  fChain -> SetBranchAddress("Events_Estrp_en12", &Events_Estrp_en12, &b_Events_Estrp_en12);
  fChain -> SetBranchAddress("Events_Estrp_en13", &Events_Estrp_en13, &b_Events_Estrp_en13);
  fChain -> SetBranchAddress("Events_Estrp_en14", &Events_Estrp_en14, &b_Events_Estrp_en14);
  fChain -> SetBranchAddress("Events_Estrp_en15", &Events_Estrp_en15, &b_Events_Estrp_en15);
  fChain -> SetBranchAddress("Events_Twr_didE", &Events_Twr_didE, &b_Events_Twr_didE);
  fChain -> SetBranchAddress("Events_Pstrip_enp01", &Events_Pstrip_enp01, &b_Events_Pstrip_enp01);
  fChain -> SetBranchAddress("Events_Pstrip_enp02", &Events_Pstrip_enp02, &b_Events_Pstrip_enp02);
  fChain -> SetBranchAddress("Events_Pstrip_enp03", &Events_Pstrip_enp03, &b_Events_Pstrip_enp03);
  fChain -> SetBranchAddress("Events_Pstrip_enp0", &Events_Pstrip_enp0, &b_Events_Pstrip_enp0);
  fChain -> SetBranchAddress("Events_Pstrip_enp1", &Events_Pstrip_enp1, &b_Events_Pstrip_enp1);
  fChain -> SetBranchAddress("Events_Pstrip_enp2", &Events_Pstrip_enp2, &b_Events_Pstrip_enp2);
  fChain -> SetBranchAddress("Events_Pstrip_enp3", &Events_Pstrip_enp3, &b_Events_Pstrip_enp3);
  fChain -> SetBranchAddress("Events_Pstrip_enp4", &Events_Pstrip_enp4, &b_Events_Pstrip_enp4);
  fChain -> SetBranchAddress("Events_Pstrip_enp5", &Events_Pstrip_enp5, &b_Events_Pstrip_enp5);
  fChain -> SetBranchAddress("Events_Pstrip_enp6", &Events_Pstrip_enp6, &b_Events_Pstrip_enp6);
  fChain -> SetBranchAddress("Events_Pstrip_enp7", &Events_Pstrip_enp7, &b_Events_Pstrip_enp7);
  fChain -> SetBranchAddress("Events_Pstrip_enp8", &Events_Pstrip_enp8, &b_Events_Pstrip_enp8);
  fChain -> SetBranchAddress("Events_Pstrip_enp9", &Events_Pstrip_enp9, &b_Events_Pstrip_enp9);
  fChain -> SetBranchAddress("Events_Pstrip_enp10", &Events_Pstrip_enp10, &b_Events_Pstrip_enp10);
  fChain -> SetBranchAddress("Events_Pstrip_enp11", &Events_Pstrip_enp11, &b_Events_Pstrip_enp11);
  fChain -> SetBranchAddress("Events_Pstrip_enp12", &Events_Pstrip_enp12, &b_Events_Pstrip_enp12);
  fChain -> SetBranchAddress("Events_Pstrip_enp13", &Events_Pstrip_enp13, &b_Events_Pstrip_enp13);
  fChain -> SetBranchAddress("Events_Pstrip_enp14", &Events_Pstrip_enp14, &b_Events_Pstrip_enp14);
  fChain -> SetBranchAddress("Events_Pstrip_enp15", &Events_Pstrip_enp15, &b_Events_Pstrip_enp15);
  fChain -> SetBranchAddress("Events_clust_Ennq1", &Events_clust_Ennq1, &b_Events_clust_Ennq1);
  fChain -> SetBranchAddress("Events_clust_Ennq20", &Events_clust_Ennq20, &b_Events_clust_Ennq20);
  fChain -> SetBranchAddress("Events_clust_Ennq19", &Events_clust_Ennq19, &b_Events_clust_Ennq19);
  fChain -> SetBranchAddress("Events_clust_Enpq1", &Events_clust_Enpq1, &b_Events_clust_Enpq1);
  fChain -> SetBranchAddress("Events_clust_Enpq20", &Events_clust_Enpq20, &b_Events_clust_Enpq20);
  fChain -> SetBranchAddress("Events_clust_Enpq19", &Events_clust_Enpq19, &b_Events_clust_Enpq19);
  fChain -> SetBranchAddress("Events_clust_Enpq21", &Events_clust_Enpq21, &b_Events_clust_Enpq21);
  fChain -> SetBranchAddress("Events_noOfprimaryTrks", &Events_noOfprimaryTrks, &b_Events_noOfprimaryTrks);
  fChain -> SetBranchAddress("pTracks_pT", pTracks_pT, &b_pTracks_pT);
  fChain -> SetBranchAddress("pTracks_px", pTracks_px, &b_pTracks_px);
  fChain -> SetBranchAddress("pTracks_py", pTracks_py, &b_pTracks_py);
  fChain -> SetBranchAddress("pTracks_pz", pTracks_pz, &b_pTracks_pz);
  fChain -> SetBranchAddress("pTracks_Eta", pTracks_Eta, &b_pTracks_Eta);
  fChain -> SetBranchAddress("pTracks_Phi", pTracks_Phi, &b_pTracks_Phi);
  fChain -> SetBranchAddress("pTracks_dEdx", pTracks_dEdx, &b_pTracks_dEdx);
  fChain -> SetBranchAddress("pTracks_chrg", pTracks_chrg, &b_pTracks_chrg);
  fChain -> SetBranchAddress("pTracks_gdca", pTracks_gdca, &b_pTracks_gdca);
  fChain -> SetBranchAddress("pTracks_Fp", pTracks_Fp, &b_pTracks_Fp);
  fChain -> SetBranchAddress("pTracks_Ppo", pTracks_Ppo, &b_pTracks_Ppo);
  fChain -> SetBranchAddress("pTracks_nSigPi", pTracks_nSigPi, &b_pTracks_nSigPi);
  fChain -> SetBranchAddress("pTracks_nSigK", pTracks_nSigK, &b_pTracks_nSigK);
  fChain -> SetBranchAddress("pTracks_nSigP", pTracks_nSigP, &b_pTracks_nSigP);
  fChain -> SetBranchAddress("pTracks_nSigE", pTracks_nSigE, &b_pTracks_nSigE);
  Notify();

}  // end 'Init(TTree*)'



void StFemtoDstMaker::Show(const Long64_t entry) {

  if (!fChain) return;
  fChain -> Show(entry);

}  // end 'Show(Long64_t)'



Int_t StFemtoDstMaker::GetEntry(const Long64_t entry) {

  if (!fChain) return 0;
  return fChain -> GetEntry(entry);

}  // end 'GetEntry(Long64_t)'



Int_t StFemtoDstMaker::Cut(const Long64_t entry) {

  return 1;

}  // end 'Cut(Long64_t)'



Bool_t StFemtoDstMaker::Notify() {

  return kTRUE;

}  // end 'Notify()'



Long64_t StFemtoDstMaker::LoadTree(const Long64_t entry) {

  if (!fChain) return -5;
  Long64_t centry = fChain -> LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain -> GetTreeNumber();
    Notify();
  }
  return centry;

}  // end 'LoadTree(Long64_t)'



Double_t StFemtoDstMaker::ApplyDetectorResponse(const Double_t pTpar) {

  // resolution paramters
  Double_t sigCstUse = 0.;
  Double_t sigPtUse  = 0.;
  Double_t sigPt2use = 0.;

  const UInt_t systematicVariation = 0;
  switch (systematicVariation) {
    case 0:
      sigCstUse = 0.014;
      sigPtUse  = 0.01;
      sigPt2use = 0.001;
      break;
    case 1:
      sigCstUse = 0.0445;
      sigPtUse  = 0.0070;
      sigPt2use = 0.0013;
      break;
    case 2:
      sigCstUse = -0.026;
      sigPtUse  = 0.020;
      sigPt2use = 0.003;
      break;
  }
  const Double_t sigCst = sigCstUse;
  const Double_t sigPt  = sigPtUse;
  const Double_t sigPt2 = sigPt2use;


  // resolution calculation
  const Double_t res   = sigCst + (sigPt * pTpar) + (sigPt2 * (pTpar * pTpar));
  const Double_t smear = gRandom -> Gaus(0., res);
  const Double_t pTdet = pTpar + smear;
  const Double_t pTres = (pTpar - pTdet) / pTpar;

  // efficiency calculation
  const Bool_t isHighPt = (pTpar > HighPt1);

  UInt_t effMethod = 0;
  if (effIsHist)
    effMethod = 1;
  else if (effIsFunc)
    effMethod = 2;
  else
    effMethod = 3;

  UInt_t   iEff(0);
  Double_t eff(0);
  switch (effMethod) {
    case 1:
      if (doEffSmooth && isHighPt) {
        iEff = hPtEff -> FindBin(HighPt1);
        eff  = hPtEff -> GetBinContent(iEff);
      }
      else {
        iEff = hPtEff -> FindBin(pTdet);
        eff  = hPtEff -> GetBinContent(iEff);
      }
      break;
    case 2:
      if (doEffSmooth && isHighPt)
        eff = fPtEff -> Eval(HighPt1);
      else
        eff = fPtEff -> Eval(pTdet);
      break;
    case 3:
      eff = fPtEff -> Eval(pTdet);
      break;
  }
  Double_t pass = gRandom -> Uniform(0., 1.);

  // adjust eff (if need be)
  Double_t effUse(eff);
  if (doEffSys) {
    const Double_t effNew = eff * (1. + effAdjust);
    effUse = effNew;
  }


  // for QA
  if (pass <= effUse) {
    hPtRes1D -> Fill(pTres);
    hPtRes2D -> Fill(pTpar, pTres);
  }

  // determine return value
  Double_t pTreturn = 0.;
  if (pass > effUse) {
    pTreturn = -1000.;
  }
  else {
    hPtEffDet -> Fill(pTpar);
    pTreturn  = pTdet;
  }
  hPtEffPar -> Fill(pTpar);

  return pTreturn;

}  // end 'ApplyDetectorResponse(Double_t)'

#endif

// End ------------------------------------------------------------------------
