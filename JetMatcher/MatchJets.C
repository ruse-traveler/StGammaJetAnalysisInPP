// 'MatchJets.C'
// Nihar Sahoo, Derek Anderson
// 12.01.2017

#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TProfile.h"
#include "TDirectory.h"

using namespace std;

// filepaths
static const TString SOutDefault("test.root");
static const TString SParDefault("../JetMaker/mc/output/pp200r9pt15ff.forSmallQtCheck_withImplicitGhosts_Par.et920pt021Kvz55had.r02rm1chrg.d30m3y2020.root");
static const TString SDetDefault("../JetMaker/mudst/output/pp200r9pt15ff.forSmallQtCheck_withImplicitGhosts_Det.noTrigger.r02rm1chrg.d30m3y2020.root");

// jet parameters
static const Double_t Rcut(0.5);       // Rcut = Rjet
static const Double_t MinArea(0.65);   // R02: 0.05, R03: 0.2, R04: 0.35, R05: 0.65, R07: 1.2
static const Double_t MinJetPt(0.2);
static const Double_t MaxJetPt(30.);   // max detector jet pT

// misc parameters
static const Bool_t DoNorm(false);
static const Bool_t UseVariablePtBins(true);
static const Bool_t UseParticleLevelTrigger(true);
static const Bool_t UseDetectorLevelTrigger(false);
static const UInt_t NHadIds(1);
static const UInt_t NJetTypes(7);
static const UInt_t NMatchTypes(5);
static const UInt_t NDirectories(NJetTypes + 1);



void MatchJets(const TString pPath=SParDefault, const TString dPath=SDetDefault, const TString oPath=SOutDefault, Bool_t inBatchMode=false) {

  cout << "\n  Beginning match script!" << endl;

  // event constants
  const Double_t MaxVz(55.);
  const Double_t MaxTrgEta(0.9);
  const Double_t MinTrgEt(9.);
  const Double_t MaxTrgEt(11.);
  const Double_t MinTrgTsp(3.);
  const Double_t MaxTrgTsp(100.);
  const UInt_t   idHadTrg[NHadIds] = {7};

  // matching constants
  const Double_t Qmin(0.5);     // fraction of reconstructed jet pT must be above Qmin
  const Double_t Qmax(1.3);     // fraction of reconstructed jet pT must be below Qmax
  const Double_t HardCut(10.);  // jets w/ pT > HardCut are considered 'hard'
  const Double_t Pcut(-11.);    // pTpar (corrected) must be above this
  const Double_t Dcut(0.2);     // pTdet (corrected) must be above this
  const Double_t FmatchMin(0.);
  const Double_t FmatchMax(TMath::TwoPi());

  // misc. constants
  const Double_t pi(TMath::Pi());
  const Double_t RecoilDf(TMath::PiOver4());

  // "ultrafine" variable binning scheme
  //const UInt_t   nPtBinsX(47);
  //const UInt_t   nPtBinsY(47);
  //const Double_t pTbinsX[48] = {-12., -8., -5., -3., -2.5, -2., -1.5, -1.25, -1., -0.8, -0.6, -0.4, -0.2, -0.1, 0., 0.1, 0.2, 0.4, 0.6, 0.8, 1., 1.25, 1.5, 2., 2.5, 3., 4., 5., 6.5, 8., 10., 12., 14.5, 17., 20., 23., 26.5, 30., 34., 38., 42.5, 47., 52., 57., 62.5, 68., 74.5, 81.};
  //const Double_t pTbinsY[48] = {-12., -8., -5., -3., -2.5, -2., -1.5, -1.25, -1., -0.8, -0.6, -0.4, -0.2, -0.1, 0., 0.1, 0.2, 0.4, 0.6, 0.8, 1., 1.25, 1.5, 2., 2.5, 3., 4., 5., 6.5, 8., 10., 12., 14.5, 17., 20., 23., 26.5, 30., 34., 38., 42.5, 47., 52., 57., 62.5, 68., 74.5, 81.};
  // "fine" variable binning scheme
  //const UInt_t   nPtBinsX(37);
  //const UInt_t   nPtBinsY(37);
  //const Double_t pTbinsX[38] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 14.5, 17., 20., 23., 26.5, 30., 34., 38., 42.5, 47., 52., 57., 62.5, 68., 74.5, 81.};
  //const Double_t pTbinsY[38] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 14.5, 17., 20., 23., 26.5, 30., 34., 38., 42.5, 47., 52., 57., 62.5, 68., 74.5, 81.};
  // "big" variable binning scheme
  //const UInt_t   nPtBinsX(26);
  //const UInt_t   nPtBinsY(26);
  //const Double_t pTbinsX[27] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 5., 8., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  //const Double_t pTbinsY[27] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 5., 8., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  // 9 - 11 GeV gamma binning scheme
  //const UInt_t   nPtBinsX(28);
  //const UInt_t   nPtBinsY(28);
  //const Double_t pTbinsX[29] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 12., 17., 23., 30., 38., 47., 57, 68., 81.};
  //const Double_t pTbinsY[29] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 12., 17., 23., 30., 38., 47., 57, 68., 81.};
  // 11 - 15 GeV gamma binning scheme
  //const UInt_t   nPtBinsX(29);
  //const UInt_t   nPtBinsY(29);
  //const Double_t pTbinsX[30] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  //const Double_t pTbinsY[30] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  // "huge" variable binning scheme
  //const UInt_t   nPtBinsX(21);
  //const UInt_t   nPtBinsY(21);
  //const Double_t pTbinsX[22] = {-11., -8., -6., -4., -2., -1., 0., 1., 2., 4., 6., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60.};
  //const Double_t pTbinsY[22] = {-11., -8., -6., -4., -2., -1., 0., 1., 2., 4., 6., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60.};
  // "giant" variable binning scheme
  const UInt_t   nPtBinsX(21);
  const UInt_t   nPtBinsY(21);
  const Double_t pTbinsX[22] = {-12., -9., -7., -5., -3., -2., 1., 2., 3., 5., 7., 9., 11., 14., 19., 24., 29., 34., 39., 44., 49., 61.};
  const Double_t pTbinsY[22] = {-12., -9., -7., -5., -3., -2., 1., 2., 3., 5., 7., 9., 11., 14., 19., 24., 29., 34., 39., 44., 49., 61.};

  // variable qT binning
  const UInt_t   nQtVar(10);
  const Double_t qTvarBin[11] = {0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0};
  cout << "    Opening files and grabbing trees..." << endl;

  // open files
  TFile *fPar = new TFile(pPath, "read");
  TFile *fDet = new TFile(dPath, "read");
  TFile *fOut = new TFile(oPath, "recreate");
  if (!fPar) {
    cerr << "PANIC: particle input file could not be opened!" << endl;
    assert(fPar);
  }
  if (!fDet) {
    cerr << "PANIC: detector input file could not be opened!" << endl;
    assert(fDet);
  }

  // get trees
  TTree *tPar;
  TTree *tDet;
  if (fPar && fDet) {
    fPar -> GetObject("JetTree", tPar);
    fDet -> GetObject("JetTree", tDet);
  }

  // declare particle event leaves
  cout << "    Setting branch addresses..." << endl;
  Int_t    pEventIndex    = 0;
  Int_t    pRunId         = 0;
  Int_t    pNJets         = 0;
  Double_t pPartonicPt    = 0.;
  Double_t pRefmult       = 0.;
  Double_t pTSP           = 0.;
  Double_t pTrgEta        = 0.;
  Double_t pTrgPhi        = 0.;
  Double_t pTrgEt         = 0.;
  Double_t pRho           = 0.;
  Double_t pSigma         = 0.;
  Double_t pVz            = 0.;
  // declare particle jet leaves
  vector<Double_t> *pJetPt     = 0;
  vector<Double_t> *pJetNCons  = 0;
  vector<Double_t> *pJetIndex  = 0;
  vector<Double_t> *pJetPtCorr = 0;
  vector<Double_t> *pJetPhi    = 0;
  vector<Double_t> *pJetEta    = 0;
  vector<Double_t> *pJetE      = 0;
  vector<Double_t> *pJetArea   = 0;
  // declare particle constituent leaves
  vector<vector<Double_t> > *pJetConsPt  = 0;
  vector<vector<Double_t> > *pJetConsEta = 0;
  vector<vector<Double_t> > *pJetConsPhi = 0;
  vector<vector<Double_t> > *pJetConsE   = 0;
  
  // declare detector event leaves
  Int_t    dEventIndex    = 0;
  Int_t    dRunId         = 0;
  Int_t    dNJets         = 0;
  Double_t dPartonicPt    = 0.;
  Double_t dRefmult       = 0.;
  Double_t dTSP           = 0.;
  Double_t dTrgEta        = 0.;
  Double_t dTrgPhi        = 0.;
  Double_t dTrgEt         = 0.;
  Double_t dRho           = 0.;
  Double_t dSigma         = 0.;
  Double_t dVz            = 0.;
  // declare detector jet leaves
  vector<Double_t> *dJetPt     = 0;
  vector<Double_t> *dJetNCons  = 0;
  vector<Double_t> *dJetIndex  = 0;
  vector<Double_t> *dJetPtCorr = 0;
  vector<Double_t> *dJetEta    = 0;
  vector<Double_t> *dJetPhi    = 0;
  vector<Double_t> *dJetE      = 0;
  vector<Double_t> *dJetArea   = 0;
  // declare detector constituent leaves
  vector<vector<Double_t> > *dJetConsPt  = 0;
  vector<vector<Double_t> > *dJetConsEta = 0;
  vector<vector<Double_t> > *dJetConsPhi = 0;
  vector<vector<Double_t> > *dJetConsE   = 0;

  // declare particle branches
  TBranch *bEventIndexP    = 0;
  TBranch *bRunIdP         = 0;
  TBranch *bNJetsP         = 0;
  TBranch *bRefmultP       = 0;
  TBranch *bPartonicPtP    = 0;
  TBranch *bTspP           = 0;
  TBranch *bTrgEtaP        = 0;
  TBranch *bTrgPhiP        = 0;
  TBranch *bTrgEtP         = 0;
  TBranch *bRhoP           = 0;
  TBranch *bSigmaP         = 0;
  TBranch *bVzP            = 0;
  TBranch *bJetPtP         = 0;
  TBranch *bJetNConsP      = 0;
  TBranch *bJetIndexP      = 0;
  TBranch *bJetPtCorrP     = 0;
  TBranch *bJetEtaP        = 0;
  TBranch *bJetPhiP        = 0;
  TBranch *bJetEP          = 0;
  TBranch *bJetAreaP       = 0;
  TBranch *bJetConsPtP     = 0;
  TBranch *bJetConsEtaP    = 0;
  TBranch *bJetConsPhiP    = 0;
  TBranch *bJetConsEP      = 0;

  // declare detector branches
  TBranch *bEventIndexD    = 0;
  TBranch *bRunIdD         = 0;
  TBranch *bNJetsD         = 0;
  TBranch *bPartonicPtD    = 0;
  TBranch *bRefmultD       = 0;
  TBranch *bTspD           = 0;
  TBranch *bTrgEtaD        = 0;
  TBranch *bTrgPhiD        = 0;
  TBranch *bTrgEtD         = 0;
  TBranch *bRhoD           = 0;
  TBranch *bSigmaD         = 0;
  TBranch *bVzD            = 0;
  TBranch *bJetPtD         = 0;
  TBranch *bJetNConsD      = 0;
  TBranch *bJetIndexD      = 0;
  TBranch *bJetPtCorrD     = 0;
  TBranch *bJetEtaD        = 0;
  TBranch *bJetPhiD        = 0;
  TBranch *bJetED          = 0;
  TBranch *bJetAreaD       = 0;
  TBranch *bJetConsPtD     = 0;
  TBranch *bJetConsEtaD    = 0;
  TBranch *bJetConsPhiD    = 0;
  TBranch *bJetConsED      = 0;

  // set particle branches
  tPar -> SetBranchAddress("eventIndex", &pEventIndex, &bEventIndexP);
  tPar -> SetBranchAddress("RunId", &pRunId, &bRunIdP);
  tPar -> SetBranchAddress("Refmult", &pRefmult, &bRefmultP);
  tPar -> SetBranchAddress("NJets", &pNJets, &bNJetsP);
  tPar -> SetBranchAddress("PartonicPt", &pPartonicPt, &bPartonicPtP);
  tPar -> SetBranchAddress("TSP", &pTSP, &bTspP);
  tPar -> SetBranchAddress("TrgEta", &pTrgEta, &bTrgEtaP);
  tPar -> SetBranchAddress("TrgPhi", &pTrgPhi, &bTrgPhiP);
  tPar -> SetBranchAddress("TrgEt", &pTrgEt, &bTrgEtP);
  tPar -> SetBranchAddress("Rho", &pRho, &bRhoP);
  tPar -> SetBranchAddress("Sigma", &pSigma, &bSigmaP);
  tPar -> SetBranchAddress("Vz", &pVz, &bVzP);
  tPar -> SetBranchAddress("JetIndex", &pJetIndex, &bJetIndexP);
  tPar -> SetBranchAddress("JetPt", &pJetPt, &bJetPtP);
  tPar -> SetBranchAddress("JetNCons", &pJetNCons, &bJetNConsP);
  tPar -> SetBranchAddress("JetPtCorr", &pJetPtCorr, &bJetPtCorrP);
  tPar -> SetBranchAddress("JetEta", &pJetEta, &bJetEtaP);
  tPar -> SetBranchAddress("JetPhi",&pJetPhi, &bJetPhiP); 
  tPar -> SetBranchAddress("JetE", &pJetE, &bJetEP); 
  tPar -> SetBranchAddress("JetArea",&pJetArea, &bJetAreaP);
  tPar -> SetBranchAddress("JetConsPt", &pJetConsPt, &bJetConsPtP);
  tPar -> SetBranchAddress("JetConsEta", &pJetConsEta, &bJetConsEtaP);
  tPar -> SetBranchAddress("JetConsPhi", &pJetConsPhi, &bJetConsPhiP);
  tPar -> SetBranchAddress("JetConsE", &pJetConsE, &bJetConsEP);

  // set detector branches
  tDet -> SetBranchAddress("eventIndex", &dEventIndex, &bEventIndexD);
  tDet -> SetBranchAddress("RunId", &dRunId, &bRunIdD);
  tDet -> SetBranchAddress("Refmult", &dRefmult, &bRefmultD);
  tDet -> SetBranchAddress("NJets", &dNJets, &bNJetsD);
  tDet -> SetBranchAddress("PartonicPt", &dPartonicPt, &bPartonicPtD);
  tDet -> SetBranchAddress("TSP", &dTSP, &bTspD);
  tDet -> SetBranchAddress("TrgEta", &dTrgEta, &bTrgEtaD);
  tDet -> SetBranchAddress("TrgPhi", &dTrgPhi, &bTrgPhiD);
  tDet -> SetBranchAddress("TrgEt", &dTrgEt, &bTrgEtD);
  tDet -> SetBranchAddress("Rho", &dRho, &bRhoD);
  tDet -> SetBranchAddress("Sigma", &dSigma, &bSigmaD);
  tDet -> SetBranchAddress("Vz", &dVz, &bVzD);
  tDet -> SetBranchAddress("JetIndex", &dJetIndex, &bJetIndexD);
  tDet -> SetBranchAddress("JetPt", &dJetPt, &bJetPtD);
  tDet -> SetBranchAddress("JetNCons", &dJetNCons, &bJetNConsD);
  tDet -> SetBranchAddress("JetPtCorr", &dJetPtCorr, &bJetPtCorrD);
  tDet -> SetBranchAddress("JetEta", &dJetEta, &bJetEtaD);
  tDet -> SetBranchAddress("JetPhi",&dJetPhi, &bJetPhiD); 
  tDet -> SetBranchAddress("JetE", &dJetE, &bJetED); 
  tDet -> SetBranchAddress("JetArea",&dJetArea, &bJetAreaD);
  tDet -> SetBranchAddress("JetConsPt", &dJetConsPt, &bJetConsPtD);
  tDet -> SetBranchAddress("JetConsEta", &dJetConsEta, &bJetConsEtaD);
  tDet -> SetBranchAddress("JetConsPhi", &dJetConsPhi, &bJetConsPhiD);
  tDet -> SetBranchAddress("JetConsE", &dJetConsE, &bJetConsED);
  cout << "    Creating histograms..." << endl;

  // results
  TH1D *hEfficiencyA;
  TH1D *hEfficiencyPt;
  TH2D *hResponseA;
  TH2D *hResponseAn;
  TH2D *hResponsePt;
  TH2D *hResponsePtN;
  TH2D *hResponsePtc;
  TH2D *hResponsePtcN;
  // event histograms
  TH1D *hRefmultP;
  TH1D *hRefmultD;
  TH1D *hNumJetsP;
  TH1D *hNumJetsD;
  TH1D *hNumToMatch;
  TH1D *hNumMatched;
  TH1D *hNumHardP;
  TH1D *hNumHardD;
  TH1D *hParArea;
  TH1D *hDetArea;
  TH1D *hParPtCorr;
  TH1D *hDetPtCorr;
  // jet histograms  [0='P',1='D',2='C',3='M',4='J',5='Y',6='N']
  TH1D *hJetArea[NJetTypes];
  TH1D *hJetEta[NJetTypes];
  TH1D *hJetPhi[NJetTypes];
  TH1D *hJetPt[NJetTypes];
  TH1D *hJetPtCorr[NJetTypes];
  TH2D *hJetPhiVsEta[NJetTypes];
  // matching histograms [0='P',1='D',2='M',3='J',4='Y']
  TH1D *hJetQt[NMatchTypes];
  TH1D *hJetDr[NMatchTypes];
  TH1D *hJetS[NMatchTypes];
  TH1D *hJetDp[NMatchTypes];
  TH1D *hJetDq[NMatchTypes];
  TH2D *hJetDrVsPt[NMatchTypes];
  TH2D *hJetQtVsPt[NMatchTypes];
  TH2D *hJetQtVsDr[NMatchTypes];
  TH2D *hJetSvsDr[NMatchTypes];

  // uniform binnings
  const UInt_t   mNum(200);
  const UInt_t   nNum(100);
  const UInt_t   aNum(500);
  const UInt_t   hNum(100);
  const UInt_t   fNum(360);
  const UInt_t   pNum(200);
  const UInt_t   qNum(500);
  const UInt_t   rNum(30);
  const UInt_t   sNum(600);
  const UInt_t   dNum(100);
  const UInt_t   jNum(200);
  const UInt_t   nCst(1000);
  const Double_t mBin[2]   = {0., 200.};
  const Double_t nBin[2]   = {0., 100};
  const Double_t aBin[2]   = {0., 5.};
  const Double_t hBin[2]   = {-5., 5.};
  const Double_t fBin[2]   = {-2.*pi, 2.*pi};
  const Double_t pBin[2]   = {-100., 100.};
  const Double_t qBin[2]   = {0., 50.};
  const Double_t rBin[2]   = {0., 3.};
  const Double_t sBin[2]   = {0., 3.};
  const Double_t dBin[2]   = {-50., 50.};
  const Double_t jBin[2]   = {-10., 10.};
  const Double_t cstBin[2] = {0., 1000.};

  // result histograms
  hEfficiencyA    = new TH1D("hEfficiencyA", "Efficiency, #epsilon(A_{jet}) = N_{det}(A_{jet})/N_{par}(A_{jet})", aNum, aBin[0], aBin[1]);
  if (UseVariablePtBins) {
    hEfficiencyPt = new TH1D("hEfficiencyPt", "Efficiency, #epsilon(p_{T}^{jet}) = N_{det}(p_{T}^{jet})/N_{par}(p_{T}^{jet})", nPtBinsX, pTbinsX);
  } else {
    hEfficiencyPt = new TH1D("hEfficiencyPt", "Efficiency, #epsilon(p_{T}^{jet}) = N_{det}(p_{T}^{jet})/N_{par}(p_{T}^{jet})", pNum, pBin[0], pBin[1]);
  }
  hResponseA      = new TH2D("hResponseA", "Response matrix, jet area; detector; particle", aNum, aBin[0], aBin[1], aNum, aBin[0], aBin[1]);
  hResponseAn     = new TH2D("hResponseAn", "Response matrix, jet area (normalized); detector; particle", aNum, aBin[0], aBin[1], aNum, aBin[0], aBin[1]);
  if (UseVariablePtBins) {
    hResponsePt   = new TH2D("hResponsePt", "Response matrix, jet p_{T}; detector; particle", nPtBinsX, pTbinsX, nPtBinsY, pTbinsY);
    hResponsePtN  = new TH2D("hResponsePtN", "Response matrix, jet p_{T} (normalized); detector; particle", nPtBinsX, pTbinsX, nPtBinsY, pTbinsY);
    hResponsePtc  = new TH2D("hResponsePtc", "Response matrix, jet p_{T}^{corr}; detector; particle", nPtBinsX, pTbinsX, nPtBinsY, pTbinsY);
    hResponsePtcN = new TH2D("hResponsePtcN", "Response matrix, jet p_{T}^{corr} (normalized); detector; particle", nPtBinsX, pTbinsX, nPtBinsY, pTbinsY);
  } else {
    hResponsePt   = new TH2D("hResponsePt", "Response matrix, jet p_{T}; detector; particle", pNum, pBin[0], pBin[1], pNum, pBin[0], pBin[1]);
    hResponsePtN  = new TH2D("hResponsePtN", "Response matrix, jet p_{T} (normalized); detector; particle", pNum, pBin[0], pBin[1], pNum, pBin[0], pBin[1]);
    hResponsePtc  = new TH2D("hResponsePtc", "Response matrix, jet p_{T}^{corr}; detector; particle", pNum, pBin[0], pBin[1], pNum, pBin[0], pBin[1]);
    hResponsePtcN = new TH2D("hResponsePtcN", "Response matrix, jet p_{T}^{corr} (normalized); detector; particle", pNum, pBin[0], pBin[1], pNum, pBin[0], pBin[1]);
  }
  // event histograms
  hRefmultP       = new TH1D("hRefmultP", "Particle Refmult", mNum, mBin[0], mBin[1]);
  hRefmultD       = new TH1D("hRefmultD", "Detector Refmult", mNum, mBin[0], mBin[1]);
  hNumJetsP       = new TH1D("hNumJetsP", "no. of jets, particle", nNum, nBin[0], nBin[1]);
  hNumJetsD       = new TH1D("hNumJetsD", "No. of jets, detector", nNum, nBin[0], nBin[1]);
  hNumToMatch     = new TH1D("hNumToMatch", "No. of jets to match (ie. no. of jets w/ pT > pTcut)", nNum, nBin[0], nBin[1]);
  hNumMatched     = new TH1D("hNumMatched", "No. of jets matched", nNum, nBin[0], nBin[1]);
  hNumHardP       = new TH1D("hNumHardP", "No. of jets w/ p_{T} above a threshold, particle", nNum, nBin[0], nBin[1]);
  hNumHardD       = new TH1D("hNumHardD", "No. of jets w/ p_{T} above a threshold, detector", nNum, nBin[0], nBin[1]);
  hParArea        = new TH1D("hParArea", "Total no. of jets to match per A_{jet} bin (for efficiency)", aNum, aBin[0], aBin[1]);
  hDetArea        = new TH1D("hDetArea", "Total no. of jets matched per A_{jet} bin (for efficiency, no smearing)", aNum, aBin[0], aBin[1]);
  if (UseVariablePtBins) {
    hParPtCorr    = new TH1D("hParPtCorr", "Total no. of jets to match per p_{T}^{corr} bin (for efficiency)", nPtBinsX, pTbinsX);
    hDetPtCorr    = new TH1D("hDetPtCorr", "Total no. of jets matched per p_{T}^{corr} bin (for efficiency, no smearing)", nPtBinsX, pTbinsX);
  } else {
    hParPtCorr    = new TH1D("hParPtCorr", "Total no. of jets to match per p_{T}^{corr} bin (for efficiency)", pNum, pBin[0], pBin[1]);
    hDetPtCorr    = new TH1D("hDetPtCorr", "Total no. of jets matched per p_{T}^{corr} bin (for efficiency, no smearing)", pNum, pBin[0], pBin[1]);
  }
  // particle jets
  hJetArea[0]     = new TH1D("hJetAreaP", "Jet area, particle", aNum, aBin[0], aBin[1]);
  hJetEta[0]      = new TH1D("hJetEtaP", "Jet eta, particle", hNum, hBin[0], hBin[1]);
  hJetPhi[0]      = new TH1D("hJetPhiP", "Jet phi, particle", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[0]     = new TH1D("hJetPtP", "Jet p_{T}, particle", nPtBinsX, pTbinsX);
    hJetPtCorr[0] = new TH1D("hJetPtCorrP", "Jet p_{T}^{corr}, particle", nPtBinsX, pTbinsX);
  } else {
    hJetPt[0]     = new TH1D("hJetPtP", "Jet p_{T}, particle", pNum, pBin[0], pBin[1]);
    hJetPtCorr[0] = new TH1D("hJetPtCorrP", "Jet p_{T}^{corr}, particle", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[0] = new TH2D("hJetPhiVsEtaP", "Jet #varphi vs. #eta, particle", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  // detector jets
  hJetArea[1]     = new TH1D("hJetAreaD", "Jet area, detector", aNum, aBin[0], aBin[1]);
  hJetEta[1]      = new TH1D("hJetEtaD", "Jet eta, detector", hNum, hBin[0], hBin[1]);
  hJetPhi[1]      = new TH1D("hJetPhiD", "Jet phi, detector", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[1]     = new TH1D("hJetPtD", "Jet p_{T}, detector", nPtBinsX, pTbinsX);
    hJetPtCorr[1] = new TH1D("hJetPtCorrD", "Jet p_{T}^{corr}, detector", nPtBinsX, pTbinsX);
  } else {
    hJetPt[1]     = new TH1D("hJetPtD", "Jet p_{T}, detector", pNum, pBin[0], pBin[1]);
    hJetPtCorr[1] = new TH1D("hJetPtCorrD", "Jet p_{T}^{corr}, detector", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[1] = new TH2D("hJetPhiVsEtaD", "Jet #varphi vs. #eta, detector", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  hJetQt[0]       = new TH1D("hJetQtD", "Jet q_{T}, detector (normalization different!)", qNum, qBin[0], qBin[1]);
  hJetDr[0]       = new TH1D("hJetDrD", "Jet #Deltar, detector (normalization different!)", rNum, rBin[0], rBin[1]);
  hJetS[0]        = new TH1D("hJetSd", "Jet s=A_{det}/A_{par}, detector (normalization different!)", sNum, sBin[0], sBin[1]);
  hJetDp[0]       = new TH1D("hJetDpD", "Jet #Deltap_{T}=p_{T}^{det}-p_{T}^{par}, detector (normalization different!)", dNum, dBin[0], dBin[1]);
  hJetDq[0]       = new TH1D("hJetDqD", "Jet #Deltaq_{T}=#Deltap_{T}/p_{T}=q_{T}-1, detector (normalization different!)", jNum, jBin[0], jBin[1]);
  if (UseVariablePtBins) {
    hJetDrVsPt[0] = new TH2D("hJetDrVsPtD", "Jet #Deltar vs. p_{T}^{jet}(detector), detector (normalization different!); p_{T}^{jet}; #Deltar", nPtBinsX, pTbinsX, rNum, rBin[0], rBin[1]);
    hJetQtVsPt[0] = new TH2D("hJetQtVsPtD", "Jet q_{T} vs. p_{T}^{jet}(detector), detector (normalization different!); p_{T}^{jet}; q_{T}", nPtBinsX, pTbinsX, qNum, qBin[0], qBin[1]);
  } else {
    hJetDrVsPt[0] = new TH2D("hJetDrVsPtD", "Jet #Deltar vs. p_{T}^{jet}(detector), detector (normalization different!); p_{T}^{jet}; #Deltar", pNum, pBin[0], pBin[1], rNum, rBin[0], rBin[1]);
    hJetQtVsPt[0] = new TH2D("hJetQtVsPtD", "Jet q_{T} vs. p_{T}^{jet}(detector), detector (normalization different!); p_{T}^{jet}; q_{T}", pNum, pBin[0], pBin[1], qNum, qBin[0], qBin[1]);
  }
  hJetQtVsDr[0]   = new TH2D("hJetQtVsDrD", "Jet q_{T} vs. #Deltar, detector (normalization different!); #Deltar; q_{T}", rNum, rBin[0], rBin[1], qNum, qBin[0], qBin[1]);
  hJetSvsDr[0]    = new TH2D("hJetSvsDrD", "Jet s vs. #Deltar, detector (normalization different!); #Deltar; s", rNum, rBin[0], rBin[1], sNum, sBin[0], sBin[1]);
  // candidate matches
  hJetArea[2]     = new TH1D("hJetAreaC", "Jet area, candidates", aNum, aBin[0], aBin[1]);
  hJetEta[2]      = new TH1D("hJetEtaC", "Jet eta, candidates", hNum, hBin[0], hBin[1]);
  hJetPhi[2]      = new TH1D("hJetPhiC", "Jet phi, candidates", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[2]     = new TH1D("hJetPtC", "Jet p_{T}, candidates", nPtBinsX, pTbinsX);
    hJetPtCorr[2] = new TH1D("hJetPtCorrC", "Jet p_{T}^{corr}, candidates", nPtBinsX, pTbinsX);
  } else {
    hJetPt[2]     = new TH1D("hJetPtC", "Jet p_{T}, candidates", pNum, pBin[0], pBin[1]);
    hJetPtCorr[2] = new TH1D("hJetPtCorrC", "Jet p_{T}^{corr}, candidates", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[2] = new TH2D("hJetPhiVsEtaC", "Jet #varphi vs. #eta, candidates", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  hJetQt[1]       = new TH1D("hJetQtC", "Jet q_{T}, candidates", qNum, qBin[0], qBin[1]);
  hJetDr[1]       = new TH1D("hJetDrC", "Jet #Deltar, candidates", rNum, rBin[0], rBin[1]);
  hJetS[1]        = new TH1D("hJetSc", "Jet s=A_{cand}/A_{par}, candidates", sNum, sBin[0], sBin[1]);
  hJetDp[1]       = new TH1D("hJetDpC", "Jet #Deltap_{T}=p_{T}^{cand}-p_{T}^{par}, candidates (normalization different!)", dNum, dBin[0], dBin[1]);
  hJetDq[1]       = new TH1D("hJetDqC", "Jet #Deltaq_{T}=#Deltap_{T}/p_{T}=q_{T}-1, candidates (normalization different!)", jNum, jBin[0], jBin[1]);
  if (UseVariablePtBins) {
    hJetDrVsPt[1] = new TH2D("hJetDrVsPtC", "Jet #Deltar vs. p_{T}^{jet}(detector), candidates (normalization different!); p_{T}^{jet}; #Deltar", nPtBinsX, pTbinsX, rNum, rBin[0], rBin[1]);
    hJetQtVsPt[1] = new TH2D("hJetQtVsPtC", "Jet q_{T} vs. p_{T}^{jet}(detector), candidates (normalization different!); p_{T}^{jet}; q_{T}", nPtBinsX, pTbinsX, qNum, qBin[0], qBin[1]);
  } else {
    hJetDrVsPt[1] = new TH2D("hJetDrVsPtC", "Jet #Deltar vs. p_{T}^{jet}(detector), candidates (normalization different!); p_{T}^{jet}; #Deltar", pNum, pBin[0], pBin[1], rNum, rBin[0], rBin[1]);
    hJetQtVsPt[1] = new TH2D("hJetQtVsPtC", "Jet q_{T} vs. p_{T}^{jet}(detector), candidates (normalization different!); p_{T}^{jet}; q_{T}", pNum, pBin[0], pBin[1], qNum, qBin[0], qBin[1]);
  }
  hJetQtVsDr[1]   = new TH2D("hJetQtVsDrC", "Jet q_{T} vs. #Deltar, candidates; #Deltar; q_{T}", rNum, rBin[0], rBin[1], qNum, qBin[0], qBin[1]);
  hJetSvsDr[1]    = new TH2D("hJetSvsDrC", "Jet s vs. #Deltar, candidates; #Deltar; s", rNum, rBin[0], rBin[1], sNum, sBin[0], sBin[1]);
  // matches
  hJetArea[3]     = new TH1D("hJetAreaM", "Jet area, matches", aNum, aBin[0], aBin[1]);
  hJetEta[3]      = new TH1D("hJetEtaM", "Jet eta, matches", hNum, hBin[0], hBin[1]);
  hJetPhi[3]      = new TH1D("hJetPhiM", "Jet phi, matches", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[3]     = new TH1D("hJetPtM", "Jet p_{T}, matches", nPtBinsX, pTbinsX);
    hJetPtCorr[3] = new TH1D("hJetPtCorrM", "Jet p_{T}^{corr}, matches", nPtBinsX, pTbinsX);
  } else {
    hJetPt[3]     = new TH1D("hJetPtM", "Jet p_{T}, matches", pNum, pBin[0], pBin[1]);
    hJetPtCorr[3] = new TH1D("hJetPtCorrM", "Jet p_{T}^{corr}, matches", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[3] = new TH2D("hJetPhiVsEtaM", "Jet #varphi vs. #eta, matches", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  hJetQt[2]       = new TH1D("hJetQtM", "Jet q_{T}, matches", qNum, qBin[0], qBin[1]);
  hJetDr[2]       = new TH1D("hJetDrM", "Jet #Deltar, matches", rNum, rBin[0], rBin[1]);
  hJetS[2]        = new TH1D("hJetSm", "Jet s=A_{match}/A_{par}, matches", sNum, sBin[0], sBin[1]);
  hJetDp[2]       = new TH1D("hJetDpM", "Jet #Deltap_{T}=p_{T}^{match}-p_{T}^{par}, matches", dNum, dBin[0], dBin[1]);
  hJetDq[2]       = new TH1D("hJetDqM", "Jet #Deltaq_{T}=#Deltap_{T}/p_{T}=q_{T}-1, matches", jNum, jBin[0], jBin[1]);
  if (UseVariablePtBins) {
    hJetDrVsPt[2] = new TH2D("hJetDrVsPtM", "Jet #Deltar vs. p_{T}^{jet}(detector), matches; p_{T}^{jet}; #Deltar", nPtBinsX, pTbinsX, rNum, rBin[0], rBin[1]);
    hJetQtVsPt[2] = new TH2D("hJetQtVsPtM", "Jet q_{T} vs. p_{T}^{jet}(detector), matches; p_{T}^{jet}; q_{T}", nPtBinsX, pTbinsX, qNum, qBin[0], qBin[1]);
  } else {
    hJetDrVsPt[2] = new TH2D("hJetDrVsPtM", "Jet #Deltar vs. p_{T}^{jet}(detector), matches; #Deltar", pNum, pBin[0], pBin[1], rNum, rBin[0], rBin[1]);
    hJetQtVsPt[2] = new TH2D("hJetQtVsPtM", "Jet q_{T} vs. p_{T}^{jet}(detector), matches; q_{T}", pNum, pBin[0], pBin[1], qNum, qBin[0], qBin[1]);
  }
  hJetQtVsDr[2]   = new TH2D("hJetQtVsDrM", "Jet q_{T} vs. #Deltar, matches", rNum, rBin[0], rBin[1], qNum, qBin[0], qBin[1]);
  hJetSvsDr[2]    = new TH2D("hJetSvsDrM", "Jet s vs. #Deltar, matches", rNum, rBin[0], rBin[1], sNum, sBin[0], sBin[1]);
  // junk (detector jets that weren't matched)
  hJetArea[4]     = new TH1D("hJetAreaJ", "Jet area, junk", aNum, aBin[0], aBin[1]);
  hJetEta[4]      = new TH1D("hJetEtaJ", "Jet eta, junk", hNum, hBin[0], hBin[1]);
  hJetPhi[4]      = new TH1D("hJetPhiJ", "Jet phi, junk", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[4]     = new TH1D("hJetPtJ", "Jet p_{T}, junk", nPtBinsX, pTbinsX);
    hJetPtCorr[4] = new TH1D("hJetPtCorrJ", "Jet p_{T}^{corr}, junk", nPtBinsX, pTbinsX);
  } else {
    hJetPt[4]     = new TH1D("hJetPtJ", "Jet p_{T}, junk", pNum, pBin[0], pBin[1]);
    hJetPtCorr[4] = new TH1D("hJetPtCorrJ", "Jet p_{T}^{corr}, junk", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[4] = new TH2D("hJetPhiVsEtaJ", "Jet #varphi vs. #eta, junk", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  hJetQt[3]       = new TH1D("hJetQtJ", "Jet q_{T}, junk (normalization different!)", qNum, qBin[0], qBin[1]);
  hJetDr[3]       = new TH1D("hJetDrJ", "Jet #Deltar, junk (normalization different!)", rNum, rBin[0], rBin[1]);
  hJetS[3]        = new TH1D("hJetSj", "Jet s=A_{junk}/A_{par}, junk (normalization different!)", sNum, sBin[0], sBin[1]);
  hJetDp[3]       = new TH1D("hJetDpJ", "Jet #Deltap_{T}=p_{T}^{junk}-p_{T}^{par}, junk (normalization different!)", dNum, dBin[0], dBin[1]);
  hJetDq[3]       = new TH1D("hJetDqJ", "Jet #Deltaq_{T}=#Deltap_{T}/p_{T}=q_{T}-1, junk (normalization different!)", jNum, jBin[0], jBin[1]);
  if (UseVariablePtBins) {
    hJetDrVsPt[3] = new TH2D("hJetDrVsPtJ", "Jet #Deltar vs. p_{T}^{jet}(detector), junk (normalization different!); p_{T}^{jet}; #Deltar", nPtBinsX, pTbinsX, rNum, rBin[0], rBin[1]);
    hJetQtVsPt[3] = new TH2D("hJetQtVsPtJ", "Jet q_{T} vs. p_{T}^{jet}(detector), junk (normalization different!); p_{T}^{jet}; q_{T}", nPtBinsX, pTbinsX, qNum, qBin[0], qBin[1]);
  } else {
    hJetDrVsPt[3] = new TH2D("hJetDrVsPtJ", "Jet #Deltar vs. p_{T}^{jet}(detector), junk (normalization different!); p_{T}^{jet}; #Deltar", pNum, pBin[0], pBin[1], rNum, rBin[0], rBin[1]);
    hJetQtVsPt[3] = new TH2D("hJetQtVsPtJ", "Jet q_{T} vs. p_{T}^{jet}(detector), junk (normalization different!); p_{T}^{jet}; q_{T}", pNum, pBin[0], pBin[1], qNum, qBin[0], qBin[1]);
  }
  hJetQtVsDr[3]   = new TH2D("hJetQtVsDrJ", "Jet q_{T} vs. #Deltar, junk (normalization different!)", rNum, rBin[0], rBin[1], qNum, qBin[0], qBin[1]);
  hJetSvsDr[3]    = new TH2D("hJetSvsDrJ", "Jet s vs. #Deltar, junk (normalization different!)", rNum, rBin[0], rBin[1], sNum, sBin[0], sBin[1]);
  // mystery (jets w/ dR > Rjet and |qT-1|<.1)
  hJetArea[5]     = new TH1D("hJetAreaY", "Jet area, mystery", aNum, aBin[0], aBin[1]);
  hJetEta[5]      = new TH1D("hJetEtaY", "Jet eta, mystery", hNum, hBin[0], hBin[1]);
  hJetPhi[5]      = new TH1D("hJetPhiY", "Jet phi, mystery", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[5]     = new TH1D("hJetPtY", "Jet p_{T}, mystery", nPtBinsX, pTbinsX);
    hJetPtCorr[5] = new TH1D("hJetPtCorrY", "Jet p_{T}^{corr}, mystery", nPtBinsX, pTbinsX);
  } else {
    hJetPt[5]     = new TH1D("hJetPtY", "Jet p_{T}, mystery", pNum, pBin[0], pBin[1]);
    hJetPtCorr[5] = new TH1D("hJetPtCorrY", "Jet p_{T}^{corr}, mystery", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[5] = new TH2D("hJetPhiVsEtaY", "Jet #varphi vs. #eta, mystery", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  hJetQt[4]       = new TH1D("hJetQtY", "Jet q_{T}, mystery", qNum, qBin[0], qBin[1]);
  hJetDr[4]       = new TH1D("hJetDrY", "Jet #Deltar, mystery", rNum, rBin[0], rBin[1]);
  hJetS[4]        = new TH1D("hJetSy", "Jet s=A_{?}/A_{par}, mystery", sNum, sBin[0], sBin[1]);
  hJetDp[4]       = new TH1D("hJetDpY", "Jet #Deltap_{T}=p_{T}^{?}-p_{T}^{par}, mystery (normalization different!)", dNum, dBin[0], dBin[1]);
  hJetDq[4]       = new TH1D("hJetDqY", "Jet #Deltaq_{T}=#Deltap_{T}/p_{T}=q_{T}-1, mystery (normalization different!)", jNum, jBin[0], jBin[1]);
  if (UseVariablePtBins) {
    hJetDrVsPt[4] = new TH2D("hJetDrVsPtY", "Jet #Deltar vs. p_{T}^{jet}(detector), mystery (normalization different!); p_{T}^{jet}; #Deltar", nPtBinsX, pTbinsX, rNum, rBin[0], rBin[1]);
    hJetQtVsPt[4] = new TH2D("hJetQtVsPtY", "Jet q_{T} vs. p_{T}^{jet}(detector), mystery (normalization different!); p_{T}^{jet}; q_{T}", nPtBinsX, pTbinsX, qNum, qBin[0], qBin[1]);
  } else {
    hJetDrVsPt[4] = new TH2D("hJetDrVsPtY", "Jet #Deltar vs. p_{T}^{jet}(detector), mystery (normalization different!); p_{T}^{jet}; #Deltar", pNum, pBin[0], pBin[1], rNum, rBin[0], rBin[1]);
    hJetQtVsPt[4] = new TH2D("hJetQtVsPtY", "Jet q_{T} vs. p_{T}^{jet}(detector), mystery (normalization different!); p_{T}^{jet}; q_{T}", pNum, pBin[0], pBin[1], qNum, qBin[0], qBin[1]);
  }
  hJetQtVsDr[4]   = new TH2D("hJetQtVsDrY", "Jet q_{T} vs. #Deltar, mystery", rNum, rBin[0], rBin[1], qNum, qBin[0], qBin[1]);
  hJetSvsDr[4]    = new TH2D("hJetSvsDrY", "Jet s vs. #Deltar, mystery", rNum, rBin[0], rBin[1], sNum, sBin[0], sBin[1]);
  // not matches (particle jets that weren't matched)
  hJetArea[6]     = new TH1D("hJetAreaN", "Jet area, not matches (particle)", aNum, aBin[0], aBin[1]);
  hJetEta[6]      = new TH1D("hJetEtaN", "Jet eta, not matches (particle)", hNum, hBin[0], hBin[1]);
  hJetPhi[6]      = new TH1D("hJetPhiN", "Jet phi, not matches (particle)", fNum, fBin[0], fBin[1]);
  if (UseVariablePtBins) {
    hJetPt[6]     = new TH1D("hJetPtN", "Jet p_{T}, not matches (particle)", nPtBinsX, pTbinsX);
    hJetPtCorr[6] = new TH1D("hJetPtCorrN", "Jet p_{T}^{corr}, not matches (particle)", nPtBinsX, pTbinsX);
  } else {
    hJetPt[6]     = new TH1D("hJetPtN", "Jet p_{T}, not matches (particle)", pNum, pBin[0], pBin[1]);
    hJetPtCorr[6] = new TH1D("hJetPtCorrN", "Jet p_{T}^{corr}, not matches (particle)", pNum, pBin[0], pBin[1]);
  }
  hJetPhiVsEta[6] = new TH2D("hJetPhiVsEtaN", "Jet #varphi vs. #eta, not matches (particle)", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);

  // errors
  hEfficiencyA   -> Sumw2();
  hEfficiencyPt  -> Sumw2();
  hResponseA     -> Sumw2();
  hResponseAn    -> Sumw2();
  hResponsePt    -> Sumw2();
  hResponsePtN   -> Sumw2();
  hResponsePtc   -> Sumw2();
  hResponsePtcN  -> Sumw2();
  hRefmultP      -> Sumw2();
  hRefmultD      -> Sumw2();
  hNumJetsP      -> Sumw2();
  hNumJetsD      -> Sumw2();
  hNumToMatch    -> Sumw2();
  hNumMatched    -> Sumw2();
  hNumHardP      -> Sumw2();
  hNumHardD      -> Sumw2();
  hParArea       -> Sumw2();
  hDetArea       -> Sumw2();
  hParPtCorr     -> Sumw2();
  hDetPtCorr     -> Sumw2();
  for (UInt_t i = 0; i < NJetTypes; i++) {
    hJetArea[i]   -> Sumw2();
    hJetEta[i]    -> Sumw2();
    hJetPhi[i]    -> Sumw2();
    hJetPt[i]     -> Sumw2();
    hJetPtCorr[i] -> Sumw2();
  }
  for (UInt_t i = 0; i < NMatchTypes; i++) {
    hJetQt[i] -> Sumw2();
    hJetDr[i] -> Sumw2();
    hJetS[i]  -> Sumw2();
    hJetDp[i] -> Sumw2();
    hJetDq[i] -> Sumw2();
  }

  // make test and heuristic histograms
  TH1D *hJetNcst_qTtest         = new TH1D("hJetNcst_qTtest", ";N^{cst}", nCst, cstBin[0], cstBin[1]);
  TH1D *hJetArea_qTtest         = new TH1D("hJetArea_qTtest", ";A^{jet}", aNum, aBin[0], aBin[1]);
  TH1D *hJetEta_qTtest          = new TH1D("hJetEta_qTtest", ";#eta^{jet}", hNum, hBin[0], hBin[1]);
  TH1D *hJetPhi_qTtest          = new TH1D("hJetPhi_qTtest", ";#varphi^{jet}", fNum, fBin[0], fBin[1]);
  TH1D *hJetPt_qTtest           = new TH1D("hJetPt_qTtest", ";p_{T}^{raw}(det)", pNum, pBin[0], pBin[1]);
  TH1D *hJetPtCorr_qTtest       = new TH1D("hJetPtCorr_qTtest", "p_{T}^{reco}(det)", pNum, pBin[0], pBin[1]);
  TH2D *hJetPhiVsEta_qTtest     = new TH2D("hJetPhiVsEta_qTtest", ";#eta^{jet};#varphi^{jet}", hNum, hBin[0], hBin[1], fNum, fBin[0], fBin[1]);
  TH2D *hJetParVsDetPt_qTtest   = new TH2D("hJetParVsDetPt_qTtest", ";p_{T}^{raw}(det);p_{T}^{raw}(par)", pNum, pBin[0], pBin[1], pNum, pBin[0], pBin[1]);
  TH2D *hJetParVsDetNcst_qTtest = new TH2D("hJetParVsDetNcst_qTtest", ";N^{cst}(det);N^{cst}(par)", nCst, cstBin[0], cstBin[1], nCst, cstBin[0], cstBin[1]);
  TH1D *hJetQt_qTtest           = new TH1D("hJetQt_qTtest", ";q_{T}^{jet} = p_{T}^{det} / p_{T}^{par}", qNum, qBin[0], qBin[1]);
  TH1D *hJetDr_qTtest           = new TH1D("hJetDr_qTtest", ";#Deltar = #sqrt{#Delta#eta^{2} + #Delta#varphi^{2}}", rNum, rBin[0], rBin[1]);
  TH1D *hJetS_qTtest            = new TH1D("hJetS_qTtest", "s^{jet} = A^{det} / A^{par}", sNum, sBin[0], sBin[1]);
  TH1D *hJetDp_qTtest           = new TH1D("hJetDp_qTtest", "#Deltap_{T}^{jet} = p_{T}^{det} - p_{T}^{par}", dNum, dBin[0], dBin[1]);
  TH1D *hJetDq_qTtest           = new TH1D("hJetDq_qTtest", "#Deltaq_{T}^{jet} = q_{T}^{jet} - 1", jNum, jBin[0], hBin[1]);
  TH1D *hQtVsPtDet_pTdet021     = new TH1D("hQtVsPtDet_pTdet021", "Matched jet q_{T}, p_{T}^{det} #in (0.2, 1) GeV/c", qNum, qBin[0], qBin[1]);
  TH1D *hQtVsPtDet_pTdet12      = new TH1D("hQtVsPtDet_pTdet12", "Matched jet q_{T}, p_{T}^{det} #in (1, 2) GeV/c", qNum, qBin[0], qBin[1]);
  TH1D *hQtVsPtDet_pTdet25      = new TH1D("hQtVsPtDet_pTdet25", "Matched jet q_{T}, p_{T}^{det} #in (2, 5) GeV/c", qNum, qBin[0], qBin[1]);
  TH1D *hQtVsPtDet_pTdet510     = new TH1D("hQtVsPtDet_pTdet510", "Matched jet q_{T}, p_{T}^{det} #in (5, 10) GeV/c", qNum, qBin[0], qBin[1]);
  TH1D *hQtVsPtDet_pTdet10      = new TH1D("hQtVsPtDet_pTdet10", "Matched jet q_{T}, p_{T}^{det} > 10 GeV/c", qNum, qBin[0], qBin[1]);
  TH1D *hQtVsPtPar_pTpar022     = new TH1D("hQtVsPtPar_pTpar022", "Matched q_{T}^{jet}, p_{T}^{par} #in (0.2, 2) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar28      = new TH1D("hQtVsPtPar_pTpar28", "Matched q_{T}^{jet}, p_{T}^{par} #in (2, 8) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar817     = new TH1D("hQtVsPtPar_pTpar817", "Matched q_{T}^{jet}, p_{T}^{par} #in (8, 17) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar0206    = new TH1D("hQtVsPtPar_pTpar0206", "Matched q_{T}^{jet}, p_{T}^{par} #in (0.2, 0.6) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar061     = new TH1D("hQtVsPtPar_pTpar061", "Matched q_{T}^{jet}, p_{T}^{par} #in (0.6, 1) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar115     = new TH1D("hQtVsPtPar_pTpar115", "Matched q_{T}^{jet}, p_{T}^{par} #in (1, 1.5) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar152     = new TH1D("hQtVsPtPar_pTpar152", "Matched q_{T}^{jet}, p_{T}^{par} #in (1.5, 2) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar1730    = new TH1D("hQtVsPtPar_pTpar1730", "Matched p_{T}^{par} vs. q_{T}^{jet}", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar3057    = new TH1D("hQtVsPtPar_pTpar3057", "Matched p_{T}^{par} vs. q_{T}^{jet}", nQtVar, qTvarBin);
  TH1D *hQtVsPtPar_pTpar5781    = new TH1D("hQtVsPtPar_pTpar5781", "Matched p_{T}^{par} vs. q_{T}^{jet}", nQtVar, qTvarBin);
  TH2D *hJetDrVsPt_qTtest       = new TH2D("hJetDrVsPt_qTtest", ";p_{T}^{raw}(det);#Deltar", pNum, pBin[0], pBin[1], rNum, rBin[0], rBin[1]);
  TH2D *hJetQtVsPt_qTtest       = new TH2D("hJetQtVsPt_qTtest", ";p_{T}^{raw}(det);q_{T}^{jet}", pNum, pBin[0], pBin[1], qNum, qBin[0], qBin[1]);
  TH2D *hJetQtVsDr_qTtest       = new TH2D("hJetQtVsDr_qTtest", ";#Deltar;q_{T}^{jet}", rNum, rBin[0], rBin[1], qNum, qBin[0], qBin[1]);
  TH2D *hJetSvsDr_qTtest        = new TH2D("hJetSvsDr_qTtest", ";#Deltar;s^{jet}", rNum, rBin[0], rBin[1], sNum, sBin[0], sBin[1]);
  hJetNcst_qTtest         -> Sumw2();
  hJetArea_qTtest         -> Sumw2();
  hJetEta_qTtest          -> Sumw2();
  hJetPhi_qTtest          -> Sumw2();
  hJetPt_qTtest           -> Sumw2();
  hJetPtCorr_qTtest       -> Sumw2();
  hJetPhiVsEta_qTtest     -> Sumw2();
  hJetParVsDetPt_qTtest   -> Sumw2();
  hJetParVsDetNcst_qTtest -> Sumw2();
  hJetQt_qTtest           -> Sumw2();
  hJetDr_qTtest           -> Sumw2();
  hJetS_qTtest            -> Sumw2();
  hJetDp_qTtest           -> Sumw2();
  hJetDq_qTtest           -> Sumw2();
  hQtVsPtDet_pTdet021     -> Sumw2();
  hQtVsPtDet_pTdet12      -> Sumw2();
  hQtVsPtDet_pTdet25      -> Sumw2();
  hQtVsPtDet_pTdet510     -> Sumw2();
  hQtVsPtDet_pTdet10      -> Sumw2();
  hQtVsPtPar_pTpar022     -> Sumw2();
  hQtVsPtPar_pTpar28      -> Sumw2();
  hQtVsPtPar_pTpar817     -> Sumw2();
  hQtVsPtPar_pTpar0206    -> Sumw2();
  hQtVsPtPar_pTpar061     -> Sumw2();
  hQtVsPtPar_pTpar115     -> Sumw2();
  hQtVsPtPar_pTpar152     -> Sumw2();
  hQtVsPtPar_pTpar1730    -> Sumw2();
  hQtVsPtPar_pTpar3057    -> Sumw2();
  hQtVsPtPar_pTpar5781    -> Sumw2();
  hJetDrVsPt_qTtest       -> Sumw2();
  hJetQtVsPt_qTtest       -> Sumw2();
  hJetQtVsDr_qTtest       -> Sumw2();
  hJetSvsDr_qTtest        -> Sumw2();

  // make histograms for smoothing calculation [05.27.2021, 06.07.2021]
  TH1D *hPtParSmooth_pTpar0206 = new TH1D("hPtParForSmooth_pTpar0206", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (0.2, 0.6) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar061  = new TH1D("hPtParForSmooth_pTpar061", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (0.6, 1) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar12   = new TH1D("hPtParForSmooth_pTpar12", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (1, 2) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar257  = new TH1D("hPtParForSmooth_pTpar257", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (2, 57) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar210  = new TH1D("hPtParForSmooth_pTpar210", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (2, 10) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar220  = new TH1D("hPtParForSmooth_pTpar220", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (2, 20) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar1057 = new TH1D("hPtParForSmooth_pTpar1057", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (10, 57) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtParSmooth_pTpar2057 = new TH1D("hPtParForSmooth_pTpar2057", "Matched p_{T}^{par} (for cross-check, p_{T}^{par} #in (20, 57) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar0206 = new TH1D("hPtDetForSmooth_pTpar0206", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (0.2, 0.6) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar061  = new TH1D("hPtDetForSmooth_pTpar061", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (0.6, 1) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar12   = new TH1D("hPtDetForSmooth_pTpar12", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (1, 2) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar257  = new TH1D("hPtDetForSmooth_pTpar257", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (2, 57) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar210  = new TH1D("hPtDetForSmooth_pTpar210", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (2, 10) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar220  = new TH1D("hPtDetForSmooth_pTpar220", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (2, 20) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar1057 = new TH1D("hPtDetForSmooth_pTpar1057", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (10, 57) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hPtDetSmooth_pTpar2057 = new TH1D("hPtDetForSmooth_pTpar2057", "Matched p_{T}^{det} (for cross-check, p_{T}^{par} #in (20, 57) GeV/c", nPtBinsX, pTbinsX);
  TH1D *hQtForSmooth_pTpar0206 = new TH1D("hQtForSmooth_pTpar0206", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (0.2, 0.6) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar061  = new TH1D("hQtForSmooth_pTpar061", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (0.6, 1) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar12   = new TH1D("hQtForSmooth_pTpar12", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (1, 2) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar257  = new TH1D("hQtForSmooth_pTpar257", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (2, 57) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar210  = new TH1D("hQtForSmooth_pTpar210", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (2, 10) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar220  = new TH1D("hQtForSmooth_pTpar220", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (2, 20) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar1057 = new TH1D("hQtForSmooth_pTpar1057", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (10, 57) GeV/c", nQtVar, qTvarBin);
  TH1D *hQtForSmooth_pTpar2057 = new TH1D("hQtForSmooth_pTpar2057", "Matched q_{T}^{jet} for smoothing, p_{T}^{par} #in (20, 57) GeV/c", nQtVar, qTvarBin);
  TH2D *hJetPtParVsQt          = new TH2D("hJetPtParVsQt", "Matched p_{T}^{par} vs. q_{T}^{jet}", nQtVar, qTvarBin, nPtBinsY, pTbinsY);
  hPtParSmooth_pTpar0206 -> Sumw2();
  hPtParSmooth_pTpar061  -> Sumw2();
  hPtParSmooth_pTpar12   -> Sumw2();
  hPtParSmooth_pTpar257  -> Sumw2();
  hPtParSmooth_pTpar210  -> Sumw2();
  hPtParSmooth_pTpar220  -> Sumw2();
  hPtParSmooth_pTpar1057 -> Sumw2();
  hPtParSmooth_pTpar2057 -> Sumw2();
  hPtDetSmooth_pTpar0206 -> Sumw2();
  hPtDetSmooth_pTpar061  -> Sumw2();
  hPtDetSmooth_pTpar12   -> Sumw2();
  hPtDetSmooth_pTpar257  -> Sumw2();
  hPtDetSmooth_pTpar210  -> Sumw2();
  hPtDetSmooth_pTpar220  -> Sumw2();
  hPtDetSmooth_pTpar1057 -> Sumw2();
  hPtDetSmooth_pTpar2057 -> Sumw2();
  hQtForSmooth_pTpar0206 -> Sumw2();
  hQtForSmooth_pTpar061  -> Sumw2();
  hQtForSmooth_pTpar12   -> Sumw2();
  hQtForSmooth_pTpar257  -> Sumw2();
  hQtForSmooth_pTpar210  -> Sumw2();
  hQtForSmooth_pTpar220  -> Sumw2();
  hQtForSmooth_pTpar1057 -> Sumw2();
  hQtForSmooth_pTpar2057 -> Sumw2();
  hJetPtParVsQt          -> Sumw2();

  // make qt and dr selection plots [10.09.2021]
  TH1D *hDrAll = new TH1D("hJetDrAll", "#Deltar of all accepted detector jets", rNum, rBin[0], rBin[1]);
  TH1D *hDrSel = new TH1D("hJetDrSel", "#Deltar of jets selected for matching", rNum, rBin[0], rBin[1]);
  TH1D *hQtAll = new TH1D("hJetQtAll", "q_{T}^{trg} of all accepted detector jets", qNum, qBin[0], qBin[1]);
  TH1D *hQtSel = new TH1D("hJetQtSel", "q_{T}^{trg} of jets selected for matching", qNum, qBin[0], qBin[1]);
  hDrAll -> Sumw2();
  hDrSel -> Sumw2();
  hQtAll -> Sumw2();
  hQtSel -> Sumw2();

  // make profiles
  TProfile *pResponseA;
  TProfile *pResponseAn;
  TProfile *pResponsePt;
  TProfile *pResponsePtN;
  TProfile *pResponsePtc;
  TProfile *pResponsePtcN;
  pResponseA      = new TProfile("pResponseA", "Response matrix, jet area; detector; particle", aNum, aBin[0], aBin[1], "S");
  pResponseAn     = new TProfile("pResponseAn", "Response matrix, jet area (normalized); detector; particle", aNum, aBin[0], aBin[1], "S");
  if (UseVariablePtBins) {
    pResponsePt   = new TProfile("pResponsePt", "Response matrix, jet p_{T}; detector; particle", nPtBinsX, pTbinsX, "S");
    pResponsePtN  = new TProfile("pResponsePtN", "Response matrix, jet p_{T} (normalized); detector; particle", nPtBinsX, pTbinsX, "S");
    pResponsePtc  = new TProfile("pResponsePtc", "Response matrix, jet p_{T}^{corr}; detector; particle", nPtBinsX, pTbinsX, "S");
    pResponsePtcN = new TProfile("pResponsePtcN", "Response matrix, jet p_{T}^{corr} (normalized); detector; particle", nPtBinsX, pTbinsX, "S");
  } else {
    pResponsePt   = new TProfile("pResponsePt", "Response matrix, jet p_{T}; detector; particle", pNum, pBin[0], pBin[1], "S");
    pResponsePtN  = new TProfile("pResponsePtN", "Response matrix, jet p_{T} (normalized); detector; particle", pNum, pBin[0], pBin[1], "S");
    pResponsePtc  = new TProfile("pResponsePtc", "Response matrix, jet p_{T}^{corr}; detector; particle", pNum, pBin[0], pBin[1], "S");
    pResponsePtcN = new TProfile("pResponsePtcN", "Response matrix, jet p_{T}^{corr} (normalized); detector; particle", pNum, pBin[0], pBin[1], "S");
  }

  // check to make sure there are a reasonable no. of events
  UInt_t fEvt  = 0;
  UInt_t nEvts = 0;
  UInt_t pEvts = (UInt_t) tPar -> GetEntries();
  UInt_t dEvts = (UInt_t) tDet -> GetEntries();
  if (pEvts < dEvts) {
    cerr << "WARNING: There are less particle-level events than detector-level!\n"
         << "         Please double-check that everything is in order...\n"
         << "         nParticle = " << pEvts << ", nDetector = " << dEvts
         << endl;
    fEvt  = 2;
    nEvts = pEvts;
  } else {
    fEvt  = 1;
    nEvts = dEvts;
  }

  // create map of tree-index to event / run no.
  vector<Int_t>          parIndices;
  vector<vector<Int_t> > pMap;
  for (UInt_t iPar = 0; iPar < pEvts; iPar++) {

    // grab entry
    tPar -> GetEntry(iPar);

    // add to map
    parIndices.clear();
    parIndices.push_back((Int_t) iPar);
    parIndices.push_back(pEventIndex);
    parIndices.push_back(pRunId);
    pMap.push_back(parIndices);
  }

  vector<Int_t>          detIndices;
  vector<vector<Int_t> > dMap;
  for (UInt_t iDet = 0; iDet < dEvts; iDet++) {

    // grab entry
    tDet -> GetEntry(iDet);

    // add to map
    detIndices.clear();
    detIndices.push_back((Int_t) iDet);
    detIndices.push_back(dEventIndex);
    detIndices.push_back(dRunId);
    dMap.push_back(detIndices);
  }

  // vectors for matching
  vector<Int_t> jetMatchIndices;
  cout << "    Beginning event loop..." << endl;

  // event loop
  Int_t  breakVal(0);
  Int_t  nBytesP(0);
  Int_t  nBytesD(0);
  UInt_t nShift(0);
  UInt_t nFound(0);
  UInt_t nTrig(0);
  for (UInt_t i = 0; i < nEvts; i++) {

    // locate event in particle or detector tree
    Int_t iParTree = -1;
    Int_t iDetTree = -1;
    switch (fEvt) {
      case 1:
        iDetTree = i;
        tDet -> GetEntry(iDetTree);
        for (UInt_t iPar = 0; iPar < pEvts; iPar++) {
          const Bool_t evtMatch = (pMap[iPar][1] == dEventIndex);
          const Bool_t runMatch = (pMap[iPar][2] == dRunId);
          if (evtMatch && runMatch) {
            iParTree = pMap[iPar][0];
            break;
          }
        }
        break;
      case 2:
        iParTree = i;
        tPar -> GetEntry(iParTree);
        for (UInt_t iDet = 0; iDet < dEvts; iDet++) {
          const Bool_t evtMatch = (dMap[iDet][1] == pEventIndex);
          const Bool_t runMatch = (dMap[iDet][2] == pRunId);
          if (evtMatch && runMatch) {
            iDetTree = dMap[iDet][0];
            break;
          }
        }
        break;
    }  // end switch case
    const Bool_t didNotFindPar = (iParTree == -1);
    const Bool_t didNotFindDet = (iDetTree == -1);
    if (didNotFindPar || didNotFindDet) continue;

    // load entries
    Int_t pBytes = tPar -> GetEntry(iParTree);
    Int_t dBytes = tDet -> GetEntry(iDetTree);
    if (pBytes < 0) {
      cerr << "ERROR: problem with particle event " << i << "...\n" << endl;
      breakVal = 1;
      break;
    }
    if (dBytes < 0) {
      cerr << "ERROR: problem with detector event " << i + nShift << "..." << endl;
      breakVal = 1;
      break;
    }

    // should be same run and event
    Bool_t isSameEvent = (pEventIndex == dEventIndex);
    Bool_t isSameRun   = (pRunId == dRunId);
    if (!isSameEvent || !isSameRun) {
      cerr << "PANIC: event index and run ID are NOT the same! Stopped at i = " << i << "\n"
           << "       ParEvt = " << pEventIndex << ", DetEvt = " << dEventIndex << "\n"
           << "       ParRun = " << pRunId << ", DetRun = " << dRunId
           << endl;
      breakVal = 1;
      break;
    }
    nFound++;

    nBytesP += pBytes;
    nBytesD += dBytes;
    if (!inBatchMode) {
      cout << "      Processing event " << i+1 << "/" << nEvts << "...\r" << flush;
      if (i+1 == nEvts) cout << endl;
    } else  {
      cout << "      Processing event " << i+1 << "/" << nEvts << "..." << endl;
    }

    // particle-level trigger info
    const Double_t eTtrgP  = pTrgEt;
    const Double_t hTrgP   = pTrgEta;
    const UInt_t   tspTrgP = (UInt_t) pTSP;

    Bool_t isHadron(false);
    for (UInt_t iHad = 0; iHad < NHadIds; iHad++) {
      if (tspTrgP == idHadTrg[iHad]) {
        isHadron = true;
        break;
      }
    }

    // particle-level trigger cuts
    const Bool_t isInEtaParCut = (TMath::Abs(hTrgP) < MaxTrgEta);
    const Bool_t isInEtParCut  = ((eTtrgP > MinTrgEt) && (eTtrgP < MaxTrgEt));
    const Bool_t isGoodParTrg  = (isInEtParCut && isInEtaParCut && isHadron);
    if (UseParticleLevelTrigger && !isGoodParTrg) continue;

    // detector-level trigger info
    const Double_t vZtrgD  = dVz;
    const Double_t eTtrgD  = dTrgEt;
    const Double_t hTrgD   = dTrgEta;
    const Double_t tspTrgD = TMath::Abs(dTSP);

    // detector-level trigger cuts
    const Bool_t isInVzDetCut  = (TMath::Abs(vZtrgD) < MaxVz);
    const Bool_t isInEtaDetCut = (TMath::Abs(hTrgD) < MaxTrgEta);
    const Bool_t isInEtDetCut  = ((eTtrgD > MinTrgEt) && (eTtrgD < MaxTrgEt));
    const Bool_t isInTspDetCut = ((tspTrgD >= MinTrgTsp) && (tspTrgD <= MaxTrgTsp));
    const Bool_t isGoodDetTrg  = (isInVzDetCut && isInEtaDetCut && isInEtDetCut && isInTspDetCut);
    if (UseDetectorLevelTrigger && !isGoodDetTrg) continue;
    nTrig++;

    // particle jet loop
    UInt_t nHardP   = 0;
    UInt_t nToMatch = 0;
    UInt_t nMatched = 0;
    UInt_t nPjets   = (UInt_t) pJetEta -> size();
    UInt_t nDjets   = (UInt_t) dJetEta -> size();
    for (UInt_t j = 0; j < nPjets; j++) {

      const Double_t pA   = pJetArea -> at(j);
      const Double_t pH   = pJetEta  -> at(j);
      const Double_t pF   = pJetPhi  -> at(j);
      const Double_t pPt  = pJetPt   -> at(j);
      const Double_t pPtc = pPt - (pRho * pA);
      if (pPt > HardCut) nHardP++;

      // calculate delta phi
      Double_t pDf = pF - pTrgPhi;
      if (pDf < ((-1. * pi) / 2.)) pDf += (2. * pi);
      if (pDf > ((3. * pi) / 2.))  pDf -= (2. * pi);
      const Double_t pDfCut    = TMath::Abs(pDf - pi);
      const Bool_t   isRecoilP = (pDfCut < RecoilDf);

      // for matching
      Double_t hMatchP(pH);
      Double_t fMatchP(pF - pTrgPhi);
      Bool_t   isInParPhiRange(false);
      while (!isInParPhiRange) {
        const Bool_t isAbovePhiMin = (fMatchP >= FmatchMin);
        const Bool_t isUnderPhiMax = (fMatchP <= FmatchMax);

        isInParPhiRange = (isAbovePhiMin && isUnderPhiMax);
        if (!isInParPhiRange) {
          if (!isAbovePhiMin) fMatchP += (2. * pi);
          if (!isUnderPhiMax) fMatchP -= (2. * pi);
        }
      }

      // consider only recoil jets
      if (!isRecoilP) continue;

      // fill particle histograms
      hJetArea[0]     -> Fill(pA);
      hJetEta[0]      -> Fill(pH);
      hJetPhi[0]      -> Fill(pF);
      hJetPt[0]       -> Fill(pPt);
      hJetPtCorr[0]   -> Fill(pPtc);
      hJetPhiVsEta[0] -> Fill(pH, pF);

      if (pPtc < Pcut) {
        continue;
      } else {
        hParArea   -> Fill(pA);    
        hParPtCorr -> Fill(pPtc);
        nToMatch++;
      }

      // matched jets ['b' for best]
      Int_t    bIndex = -999;
      Double_t bH     = -999.;
      Double_t bPt    = -999.;
      Double_t bPtc   = -999.;
      Double_t bF     = -999.;
      Double_t bA     = -999.;
      Double_t bQt    = -999.;
      Double_t bS     = -999.;
      Double_t bDp    = -999.;
      Double_t bDq    = -999.;
      Double_t bDr    = -999.;

      // detector jet loop
      Int_t  bNcst      = -999;
      UInt_t nCandidate = 0;
      Bool_t isMatched  = false;
      for (UInt_t k = 0; k < nDjets; k++) {

        const Double_t dA   = dJetArea -> at(k);
        const Double_t dH   = dJetEta  -> at(k);
        const Double_t dF   = dJetPhi  -> at(k);
        const Double_t dPt  = dJetPt   -> at(k);
        const Double_t dPtc = dPt - (dRho * dA);

        // calculate delta phi
        Double_t dDf = dF - pTrgPhi;
        if (dDf < ((-1. * pi) / 2.)) dDf += (2. * pi);
        if (dDf > ((3. * pi) / 2.))  dDf -= (2. * pi);
        const Double_t dDfCut    = TMath::Abs(dDf - pi);
        const Bool_t   isRecoilD = (dDfCut < RecoilDf);

        // for matching
        Double_t hMatchD(dH);
        Double_t fMatchD(dF - pTrgPhi);
        Bool_t   isInDetPhiRange(false);
        while (!isInDetPhiRange) {
          const Bool_t isAbovePhiMin = (fMatchD >= FmatchMin);
          const Bool_t isUnderPhiMax = (fMatchD <= FmatchMax);

          isInDetPhiRange = (isAbovePhiMin && isUnderPhiMax);
          if (!isInDetPhiRange) {
            if (!isAbovePhiMin) fMatchD += (2. * pi);
            if (!isUnderPhiMax) fMatchD -= (2. * pi);
          }
        }

        // detector cuts
        const Bool_t isInJetAcut    = (dA > MinArea);
        const Bool_t isInJetPtCut   = (dPt > MinJetPt);
        const Bool_t isInAcceptance = (isInJetAcut && isInJetPtCut && isRecoilD);

        // match jets
        Double_t qT   = dPt / pPt;
        Double_t s    = dA / pA;
        Double_t dP   = dPt - pPt;
        Double_t dQ   = dP / pPt;
        Double_t dHpd = TMath::Abs(hMatchD - hMatchP);
        Double_t dFpd = TMath::Abs(fMatchD - fMatchP);
        if (dPtc < Dcut) {
          continue;
        }
        if (dFpd > pi) {
          dFpd = (2. * pi) - dFpd;
        }

        // calculate delta-r
        Double_t dR2 = (dHpd * dHpd) + (dFpd * dFpd);
        Double_t dR  = TMath::Sqrt(dR2);

        // fill detector histograms
        hJetQt[0]     -> Fill(qT);
        hJetDr[0]     -> Fill(dR);
        hJetS[0]      -> Fill(s);
        hJetDp[0]     -> Fill(dP);
        hJetDq[0]     -> Fill(dQ);
        hJetDrVsPt[0] -> Fill(dPt, dR);
        hJetQtVsPt[0] -> Fill(dPt, qT);
        hJetQtVsDr[0] -> Fill(dR, qT);
        hJetSvsDr[0]  -> Fill(dR, s);

        // fill test and heuristic histograms
        if ((dPt >= 0.2) && (dPt < 1.)) {
          hQtVsPtDet_pTdet021 -> Fill(qT);
        }
        if ((dPt >= 1.) && (dPt < 2.)) {
          hQtVsPtDet_pTdet12  -> Fill(qT);
        }
        if ((dPt >= 2.) && (dPt < 5.)) {
          hQtVsPtDet_pTdet25  -> Fill(qT);
        }
        if ((dPt >= 5.) && (dPt < 10.)) {
          hQtVsPtDet_pTdet510 -> Fill(qT);
        }
        if (dPt >= 10.) {
          hQtVsPtDet_pTdet10  -> Fill(qT);
        }


        Bool_t   isBetter   = false;
        Bool_t   isInRcut   = (dR < Rcut);
        Bool_t   isInQcut   = ((qT > Qmin) && (qT < Qmax));
        Double_t matchValue = TMath::Abs(qT - 1.);
        Double_t bestValue  = TMath::Abs(bQt - 1.);
        if (isInRcut && isInQcut && isInAcceptance) {
          isMatched = true;
          isBetter  = (matchValue < bestValue);

          // fill candidate histograms
          hJetArea[2]     -> Fill(dA);
          hJetEta[2]      -> Fill(dH);
          hJetPhi[2]      -> Fill(dF);
          hJetPt[2]       -> Fill(dPt);
          hJetPtCorr[2]   -> Fill(dPtc);
          hJetPhiVsEta[2] -> Fill(dH, dF);
          hJetQt[1]       -> Fill(qT);
          hJetDr[1]       -> Fill(dR);
          hJetS[1]        -> Fill(s);
          hJetDp[1]       -> Fill(dP);
          hJetDq[1]       -> Fill(dQ);
          hJetDrVsPt[1]   -> Fill(dPt, dR);
          hJetQtVsPt[1]   -> Fill(dPt, qT);
          hJetQtVsDr[1]   -> Fill(dR, qT);
          hJetSvsDr[1]    -> Fill(dR, s);
          nCandidate++;
        } else {
          // fill junk histograms
          hJetQt[3]     -> Fill(qT);
          hJetDr[3]     -> Fill(dR);
          hJetS[3]      -> Fill(s);
          hJetDp[3]     -> Fill(dP);
          hJetDq[3]     -> Fill(dQ);
          hJetDrVsPt[3] -> Fill(dPt, dR);
          hJetQtVsPt[3] -> Fill(dPt, qT);
          hJetQtVsDr[3] -> Fill(dR, qT);
          hJetSvsDr[3]  -> Fill(dR, s);
        }

        // fill mystery histograms
        Double_t qCut      = TMath::Abs(qT - 1);
        Bool_t   isNearOne = (qCut < 0.1);
        if (!isInRcut && isNearOne) {
          hJetArea[5]     -> Fill(dA);
          hJetEta[5]      -> Fill(dH);
          hJetPhi[5]      -> Fill(dF);
          hJetPt[5]       -> Fill(dPt);
          hJetPtCorr[5]   -> Fill(dPtc);
          hJetPhiVsEta[5] -> Fill(dH, dF);
          hJetQt[4]       -> Fill(qT);
          hJetDr[4]       -> Fill(dR);
          hJetS[4]        -> Fill(s);
          hJetDp[4]       -> Fill(dP);
          hJetDq[4]       -> Fill(dQ);
          hJetDrVsPt[4]   -> Fill(dPt, dR);
          hJetQtVsPt[4]   -> Fill(dPt, qT);
          hJetQtVsDr[4]   -> Fill(dR, qT);
          hJetSvsDr[4]    -> Fill(dR, s);
        }

        // check if candidate is best match
        if (isMatched && isBetter) {
          bIndex = k;
          bA     = dA;
          bH     = dH;
          bF     = dF;
          bPt    = dPt;
          bPtc   = dPtc;
          bQt    = qT;
          bS     = s;
          bDp    = dP;
          bDq    = dQ;
          bDr    = dR;
          bNcst  = dJetNCons -> at(k);
        }

        // fill qt and dr selection plots [10.09.2021]
        if (isInAcceptance) {
          hDrAll -> Fill(dR);
          hQtAll -> Fill(qT);
          if (isInRcut) hDrSel -> Fill(dR);
          if (isInQcut) hQtSel -> Fill(qT);
        }
      }  // end detector jet loop

      // fill match histograms
      if (isMatched) {
        hResponseA      -> Fill(bA, pA);
        pResponseA      -> Fill(bA, pA);
        hResponseAn     -> Fill(bA, pA);
        pResponseAn     -> Fill(bA, pA);
        hResponsePt     -> Fill(bPt, pPt);
        pResponsePt     -> Fill(bPt, pPt);
        hResponsePtN    -> Fill(bPt, pPt);
        pResponsePtN    -> Fill(bPt, pPt);
        hResponsePtc    -> Fill(bPtc, pPtc);
        pResponsePtc    -> Fill(bPtc, pPtc);
        hResponsePtcN   -> Fill(bPtc, pPtc);
        pResponsePtcN   -> Fill(bPtc, pPtc);
        hDetArea        -> Fill(pA);
        hDetPtCorr      -> Fill(pPtc);
        hJetArea[3]     -> Fill(bA);
        hJetEta[3]      -> Fill(bH);
        hJetPhi[3]      -> Fill(bF);
        hJetPt[3]       -> Fill(bPt);
        hJetPtCorr[3]   -> Fill(bPtc);
        hJetPhiVsEta[3] -> Fill(bH, bF);
        hJetQt[2]       -> Fill(bQt);
        hJetDr[2]       -> Fill(bDr);
        hJetS[2]        -> Fill(bS);
        hJetDp[2]       -> Fill(bDp);
        hJetDq[2]       -> Fill(bDq);
        hJetDrVsPt[2]   -> Fill(bPt, bDr);
        hJetQtVsPt[2]   -> Fill(bPt, bQt);
        hJetQtVsDr[2]   -> Fill(bDr, bQt);
        hJetSvsDr[2]    -> Fill(bDr, bS);
        jetMatchIndices.push_back(bIndex);
        nMatched++;

        // fill test and heuristic algorithms
        Bool_t qTisSmall = ((bQt > 0.15) && (bQt < 0.4));
        if (qTisSmall) {
          hJetNcst_qTtest         -> Fill(bNcst);
          hJetArea_qTtest         -> Fill(bA);
          hJetEta_qTtest          -> Fill(bH);
          hJetPhi_qTtest          -> Fill(bF);
          hJetPt_qTtest           -> Fill(bPt);
          hJetPtCorr_qTtest       -> Fill(bPtc);
          hJetPhiVsEta_qTtest     -> Fill(bH, bF);
          hJetParVsDetPt_qTtest   -> Fill(bPt, pPt);
          hJetParVsDetNcst_qTtest -> Fill(bNcst, pJetNCons -> at(j));
          hJetQt_qTtest           -> Fill(bQt);
          hJetDr_qTtest           -> Fill(bDr);
          hJetS_qTtest            -> Fill(bS);
          hJetDp_qTtest           -> Fill(bDp);
          hJetDq_qTtest           -> Fill(bDq);
          hJetDrVsPt_qTtest       -> Fill(bPt, bDr);
          hJetQtVsPt_qTtest       -> Fill(bPt, bQt);
          hJetQtVsDr_qTtest       -> Fill(bDr, bQt);
          hJetSvsDr_qTtest        -> Fill(bDr, bS);
        }

        // fill test histograms
        if ((pPtc > 0.2) && (pPtc < 0.6)) {
          hQtVsPtPar_pTpar0206 -> Fill(bQt);
        }
        if ((pPtc > 0.6) && (pPtc < 1.)) {
          hQtVsPtPar_pTpar061  -> Fill(bQt);
        }
        if ((pPtc > 1.) && (pPtc < 1.5)) {
          hQtVsPtPar_pTpar115  -> Fill(bQt);
        }
        if ((pPtc > 1.5) && (pPtc < 2.)) {
          hQtVsPtPar_pTpar152  -> Fill(bQt);
        }
        if ((pPtc > 0.2) && (pPtc < 2.)) {
          hQtVsPtPar_pTpar022 -> Fill(bQt);
        }
        if ((pPtc > 2.) && (pPtc < 8.)) {
          hQtVsPtPar_pTpar28  -> Fill(bQt);
        }
        if ((pPtc > 8.) && (pPtc < 17.)) {
          hQtVsPtPar_pTpar817 -> Fill(bQt);
        }
        if ((pPtc > 17.) && (pPtc < 30.)) {
          hQtVsPtPar_pTpar1730 -> Fill(bQt);
        }
        if ((pPtc > 30.) && (pPtc < 57.)) {
          hQtVsPtPar_pTpar3057 -> Fill(bQt);
        }
        if ((pPtc > 57.) && (pPtc < 81.)) {
          hQtVsPtPar_pTpar5781 -> Fill(bQt);
        }

        // fill histograms for smoothing calculation
        if ((pPtc > 0.2) && (pPtc < 0.6)) {
          hPtParSmooth_pTpar0206 -> Fill(pPtc);
          hPtDetSmooth_pTpar0206 -> Fill(bPtc);
          hQtForSmooth_pTpar0206 -> Fill(bQt);
        }
        if ((pPtc > 0.6) && (pPtc < 1.)) {
          hPtParSmooth_pTpar061  -> Fill(pPtc);
          hPtDetSmooth_pTpar061  -> Fill(bPtc);
          hQtForSmooth_pTpar061  -> Fill(bQt);
        }
        if ((pPtc > 1.) && (pPtc < 2.)) {
          hPtParSmooth_pTpar12   -> Fill(pPtc);
          hPtDetSmooth_pTpar12   -> Fill(bPtc);
          hQtForSmooth_pTpar12   -> Fill(bQt);
        }
        if ((pPtc > 2.) && (pPtc < 57.)) {
          hPtParSmooth_pTpar257  -> Fill(pPtc);
          hPtDetSmooth_pTpar257  -> Fill(bPtc);
          hQtForSmooth_pTpar257  -> Fill(bQt);
        }
        if ((pPtc > 2.) && (pPtc < 10.)) {
          hPtParSmooth_pTpar210  -> Fill(pPtc);
          hPtDetSmooth_pTpar210  -> Fill(bPtc);
          hQtForSmooth_pTpar210  -> Fill(bQt);
        }
        if ((pPtc > 2.) && (pPtc < 20.)) {
          hPtParSmooth_pTpar220  -> Fill(pPtc);
          hPtDetSmooth_pTpar220  -> Fill(bPtc);
          hQtForSmooth_pTpar220  -> Fill(bQt);
        }
        if ((pPtc > 10.) && (pPtc < 57.)) {
          hPtParSmooth_pTpar1057 -> Fill(pPtc);
          hPtDetSmooth_pTpar1057 -> Fill(bPtc);
          hQtForSmooth_pTpar1057 -> Fill(bQt);
        }
        if ((pPtc > 20.) && (pPtc < 57.)) {
          hPtParSmooth_pTpar2057 -> Fill(pPtc);
          hPtDetSmooth_pTpar2057 -> Fill(bPtc);
          hQtForSmooth_pTpar2057 -> Fill(bQt);
        }
        hJetPtParVsQt -> Fill(bQt, pPtc);
      } else {
        hJetArea[6]     -> Fill(pA);
        hJetEta[6]      -> Fill(pH);
        hJetPhi[6]      -> Fill(pF);
        hJetPt[6]       -> Fill(pPt);
        hJetPtCorr[6]   -> Fill(pPtc);
        hJetPhiVsEta[6] -> Fill(pH, pF);
      }
    }  // end particle jet loop

    UInt_t matchSize = (UInt_t) jetMatchIndices.size();
    if (nMatched != matchSize) {
      cerr << "ERROR: jetMatchIndices did something weird in event " << i << "..." << endl;
      breakVal = 1;
      break;
    }

    // detector jet loop
    UInt_t nHardD = 0;
    for (UInt_t j = 0; j < nDjets; j++) {

      Double_t dA   = dJetArea -> at(j);
      Double_t dH   = dJetEta  -> at(j);
      Double_t dF   = dJetPhi  -> at(j);
      Double_t dPt  = dJetPt   -> at(j);
      Double_t dPtc = dPt - (dRho * dA);
      if (dPt > HardCut) nHardD++;

      // calculate delta phi
      Double_t dDf = dF - dTrgPhi;
      if (dDf < ((-1. * pi) / 2.)) dDf += (2. * pi);
      if (dDf > ((3. * pi) / 2.))  dDf -= (2. * pi);
      const Double_t dDfCut    = TMath::Abs(dDf - pi);
      const Bool_t   isRecoilD = (dDfCut < RecoilDf);

      if (dPt < MinJetPt) {
        continue;
      }
      if (dPt > MaxJetPt) {
        continue;
      }
      if (dA < MinArea) {
        continue;
      }
      if (!isRecoilD) {
        continue;
      }

      // fill detector histograms
      hJetArea[1]     -> Fill(dA);
      hJetEta[1]      -> Fill(dH);
      hJetPhi[1]      -> Fill(dF);
      hJetPt[1]       -> Fill(dPt);
      hJetPtCorr[1]   -> Fill(dPtc);
      hJetPhiVsEta[1] -> Fill(dH, dF);

      // check if detector jet matches particle jet
      Bool_t isMatch = false;
      for (UInt_t k = 0; k < nMatched; k++) {
        UInt_t m = (UInt_t) jetMatchIndices.at(k);
        if (j == m) {
          isMatch = true;
          break;
        }
      }

      // fill junk histograms
      if (!isMatch) {
        hJetArea[4]     -> Fill(dA);
        hJetEta[4]      -> Fill(dH);
        hJetPhi[4]      -> Fill(dF);
        hJetPt[4]       -> Fill(dPt);
        hJetPtCorr[4]   -> Fill(dPtc);
        hJetPhiVsEta[4] -> Fill(dH, dF);
      }
    }  // end detector jet loop

    // fill event histograms
    hRefmultP   -> Fill(pRefmult);
    hRefmultD   -> Fill(dRefmult);
    hNumJetsP   -> Fill(pNJets);
    hNumJetsD   -> Fill(dNJets);
    hNumToMatch -> Fill(nToMatch);
    hNumMatched -> Fill(nMatched);
    hNumHardP   -> Fill(nHardP);
    hNumHardD   -> Fill(nHardD);

    // clear vector
    jetMatchIndices.clear();

  }  // end event loop

  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    assert(0);
  } else {
    cout << "    Event loop finished!\n"
         << "      nFound = " << nFound << "\n"
         << "      nTrig  = " << nTrig
         << endl;
  }
  cout << "    Calculating efficiency..." << endl;

  // calculate efficiency
  hEfficiencyA  -> Divide(hDetArea, hParArea, 1., 1.);
  hEfficiencyPt -> Divide(hDetPtCorr, hParPtCorr, 1., 1.);
  cout << "    Normalizing histograms..." << endl;

  // bin widths
  const Double_t mWidth = (mBin[1] - mBin[0]) / mNum;
  const Double_t nWidth = (nBin[1] - nBin[0]) / nNum;
  const Double_t aWidth = (aBin[1] - aBin[0]) / aNum;
  const Double_t hWidth = (hBin[1] - hBin[0]) / hNum;
  const Double_t fWidth = (fBin[1] - fBin[0]) / fNum;
  const Double_t pWidth = (pBin[1] - pBin[0]) / pNum;
  // overall normalizations
  const Double_t mNorm  = nEvts * mWidth;
  const Double_t nNorm  = nEvts * nWidth;
  const Double_t aNorm  = nEvts * aWidth;
  const Double_t hNorm  = nEvts * hWidth;
  const Double_t fNorm  = nEvts * fWidth;
  const Double_t pNorm  = nEvts * pWidth;
  const Double_t hfNorm = hNorm * fWidth;
  // normalize histograms
  if (DoNorm) {
    hRefmultP   -> Scale(1. / mNorm);
    hRefmultD   -> Scale(1. / mNorm);
    hNumJetsP   -> Scale(1. / nNorm);
    hNumJetsD   -> Scale(1. / nNorm);
    hNumToMatch -> Scale(1. / nNorm);
    hNumMatched -> Scale(1. / nNorm);
    hNumHardP   -> Scale(1. / nNorm);
    hNumHardD   -> Scale(1. / nNorm);
    for (UInt_t i = 0; i < NJetTypes; i++) {
      hJetArea[i]     -> Scale(1. / aNorm);
      hJetEta[i]      -> Scale(1. / hNorm);
      hJetPhi[i]      -> Scale(1. / fNorm);
      hJetPt[i]       -> Scale(1. / pNorm);
      hJetPtCorr[i]   -> Scale(1. / pNorm);
      hJetPhiVsEta[i] -> Scale(1. / hfNorm);
    }
  }

  // normalize response matrices
  const UInt_t nAbinsX = hResponseAn   -> GetNbinsX();
  const UInt_t nAbinsY = hResponseAn   -> GetNbinsY();
  const UInt_t nPbinsX = hResponsePtN  -> GetNbinsX();
  const UInt_t nPbinsY = hResponsePtN  -> GetNbinsY();
  const UInt_t nCbinsX = hResponsePtcN -> GetNbinsX();
  const UInt_t nCbinsY = hResponsePtcN -> GetNbinsY();
  for (UInt_t i = 1; i < nAbinsY+1; i++) {
    const Double_t aNorm = hResponseAn -> Integral(1, nAbinsX, i, i);
    if (aNorm == 0.) continue;

    for (UInt_t j = 1; j < nAbinsX+1; j++) {
      const Double_t oldCnt = hResponseAn -> GetBinContent(j, i);
      const Double_t oldErr = hResponseAn -> GetBinError(j, i);
      const Double_t xWidth = hResponseAn -> GetXaxis() -> GetBinWidth(j);
      const Double_t yWidth = hResponseAn -> GetYaxis() -> GetBinWidth(i);
      const Double_t dArea  = xWidth * yWidth;
      const Double_t newCnt = (oldCnt / aNorm) * dArea;
      const Double_t newErr = (oldErr / aNorm) * dArea;
      hResponseAn -> SetBinContent(j, i, newCnt);
      hResponseAn -> SetBinError(j, i, newErr);
    }
  }
  for (UInt_t i = 1; i < nPbinsY+1; i++) {
    const Double_t pNorm = hResponsePtN -> Integral(1, nPbinsX, i, i);
    if (pNorm == 0.) continue;

    for (UInt_t j = 1; j < nPbinsX+1; j++) {
      const Double_t oldCnt = hResponsePtN -> GetBinContent(j, i);
      const Double_t oldErr = hResponsePtN -> GetBinError(j, i);
      const Double_t xWidth = hResponsePtN -> GetXaxis() -> GetBinWidth(j);
      const Double_t yWidth = hResponsePtN -> GetYaxis() -> GetBinWidth(i);
      const Double_t dArea  = xWidth * yWidth;
      const Double_t newCnt = (oldCnt / pNorm) * dArea;
      const Double_t newErr = (oldErr / pNorm) * dArea;
      hResponsePtN -> SetBinContent(j, i, newCnt);
      hResponsePtN -> SetBinError(j, i, newErr);
    }
  }
  for (UInt_t i = 1; i < nCbinsY+1; i++) {
    const Double_t cNorm = hResponsePtcN -> Integral(1, nCbinsX, i, i);
    if (cNorm == 0.) continue;

    for (UInt_t j = 1; j < nCbinsX+1; j++) {
      const Double_t oldCnt = hResponsePtcN -> GetBinContent(j, i);
      const Double_t oldErr = hResponsePtcN -> GetBinError(j, i);
      const Double_t xWidth = hResponsePtcN -> GetXaxis() -> GetBinWidth(j);
      const Double_t yWidth = hResponsePtcN -> GetYaxis() -> GetBinWidth(i);
      const Double_t dArea  = xWidth * yWidth;
      const Double_t newCnt = (oldCnt / cNorm) * dArea;
      const Double_t newErr = (oldErr / cNorm) * dArea;
      hResponsePtcN -> SetBinContent(j, i, newCnt);
      hResponsePtcN -> SetBinError(j, i, newErr);
    }
  }

  // create directory structure
  TDirectory *dTest;
  TDirectory *dSmooth;
  TDirectory *dir[NDirectories];
  dir[7]  = (TDirectory*) fOut -> mkdir("EventInfo");
  dir[0]  = (TDirectory*) fOut -> mkdir("ParticleJets");
  dir[1]  = (TDirectory*) fOut -> mkdir("DetectorJets");
  dir[2]  = (TDirectory*) fOut -> mkdir("Candidates");
  dir[3]  = (TDirectory*) fOut -> mkdir("MatchJets");
  dir[4]  = (TDirectory*) fOut -> mkdir("JunkJets");
  dir[5]  = (TDirectory*) fOut -> mkdir("Mystery");
  dir[6]  = (TDirectory*) fOut -> mkdir("NotMatches");
  dTest   = (TDirectory*) fOut -> mkdir("TestAndQA");
  dSmooth = (TDirectory*) fOut -> mkdir("ForSmoothing");

  // save test and heuristic histograms
  dTest                   -> cd();
  hJetNcst_qTtest         -> Write();
  hJetArea_qTtest         -> Write();
  hJetEta_qTtest          -> Write();
  hJetPhi_qTtest          -> Write();
  hJetPt_qTtest           -> Write();
  hJetPtCorr_qTtest       -> Write();
  hJetPhiVsEta_qTtest     -> Write();
  hJetParVsDetPt_qTtest   -> Write();
  hJetParVsDetNcst_qTtest -> Write();
  hJetQt_qTtest           -> Write();
  hJetDr_qTtest           -> Write();
  hJetS_qTtest            -> Write();
  hJetDp_qTtest           -> Write();
  hJetDq_qTtest           -> Write();
  hQtVsPtDet_pTdet021     -> Write();
  hQtVsPtDet_pTdet12      -> Write();
  hQtVsPtDet_pTdet25      -> Write();
  hQtVsPtDet_pTdet510     -> Write();
  hQtVsPtDet_pTdet10      -> Write();
  hQtVsPtPar_pTpar022     -> Write();
  hQtVsPtPar_pTpar28      -> Write();
  hQtVsPtPar_pTpar817     -> Write();
  hQtVsPtPar_pTpar0206    -> Write();
  hQtVsPtPar_pTpar061     -> Write();
  hQtVsPtPar_pTpar115     -> Write();
  hQtVsPtPar_pTpar152     -> Write();
  hQtVsPtPar_pTpar1730    -> Write();
  hQtVsPtPar_pTpar3057    -> Write();
  hQtVsPtPar_pTpar5781    -> Write();
  hJetDrVsPt_qTtest       -> Write();
  hJetQtVsPt_qTtest       -> Write();
  hJetQtVsDr_qTtest       -> Write();
  hJetSvsDr_qTtest        -> Write();
  hDrAll                  -> Write();
  hDrSel                  -> Write();
  hQtAll                  -> Write();
  hQtSel                  -> Write();

  // save histograms for smoothing calculation
  dSmooth                -> cd();
  hPtParSmooth_pTpar0206 -> Write();
  hPtParSmooth_pTpar061  -> Write();
  hPtParSmooth_pTpar12   -> Write();
  hPtParSmooth_pTpar257  -> Write();
  hPtParSmooth_pTpar210  -> Write();
  hPtParSmooth_pTpar220  -> Write();
  hPtParSmooth_pTpar1057 -> Write();
  hPtParSmooth_pTpar2057 -> Write();
  hPtDetSmooth_pTpar0206 -> Write();
  hPtDetSmooth_pTpar061  -> Write();
  hPtDetSmooth_pTpar12   -> Write();
  hPtDetSmooth_pTpar257  -> Write();
  hPtDetSmooth_pTpar210  -> Write();
  hPtDetSmooth_pTpar220  -> Write();
  hPtDetSmooth_pTpar1057 -> Write();
  hPtDetSmooth_pTpar2057 -> Write();
  hQtForSmooth_pTpar0206 -> Write();
  hQtForSmooth_pTpar061  -> Write();
  hQtForSmooth_pTpar12   -> Write();
  hQtForSmooth_pTpar257  -> Write();
  hQtForSmooth_pTpar210  -> Write();
  hQtForSmooth_pTpar220  -> Write();
  hQtForSmooth_pTpar1057 -> Write();
  hQtForSmooth_pTpar2057 -> Write();
  hJetPtParVsQt          -> Write();

  // write all other histograms and close output
  fOut          -> cd();
  hEfficiencyA  -> Write();
  hEfficiencyPt -> Write();
  hResponseA    -> Write();
  pResponseA    -> Write();
  hResponseAn   -> Write();
  pResponseAn   -> Write();
  hResponsePt   -> Write();
  pResponsePt   -> Write();
  hResponsePtN  -> Write();
  pResponsePtN  -> Write();
  hResponsePtc  -> Write();
  pResponsePtc  -> Write();
  hResponsePtcN -> Write();
  pResponsePtcN -> Write();
  dir[7]        -> cd();
  hRefmultP     -> Write();
  hRefmultD     -> Write();
  hNumJetsP     -> Write();
  hNumJetsD     -> Write();
  hNumToMatch   -> Write();
  hNumMatched   -> Write();
  hNumHardP     -> Write();
  hNumHardD     -> Write();
  hParArea      -> Write();
  hDetArea      -> Write();
  hParPtCorr    -> Write();
  hDetPtCorr    -> Write();
  for (UInt_t i = 0; i < NJetTypes; i++) {
    dir[i]          -> cd();
    hJetArea[i]     -> Write();
    hJetEta[i]      -> Write();
    hJetPhi[i]      -> Write();
    hJetPt[i]       -> Write();
    hJetPtCorr[i]   -> Write();
    hJetPhiVsEta[i] -> Write();
  }
  for (UInt_t i = 0; i < NMatchTypes; i++) {
    Int_t iDir = i + 1;
    dir[iDir]     -> cd();
    hJetQt[i]     -> Write();
    hJetDr[i]     -> Write();
    hJetS[i]      -> Write();
    hJetDp[i]     -> Write();
    hJetDq[i]     -> Write();
    hJetDrVsPt[i] -> Write();
    hJetQtVsPt[i] -> Write();
    hJetQtVsDr[i] -> Write();
    hJetSvsDr[i]  -> Write();
  }
  fOut -> cd();
  fOut -> Close();

  // close input
  fPar -> cd();
  fPar -> Close();
  fDet -> cd();
  fDet -> Close();
  cout << "  Matching script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
