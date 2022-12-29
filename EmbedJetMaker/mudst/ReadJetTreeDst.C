// 'ReadJetTreeDst.C'
// Nihar Sahoo, Derek Anderson
// 11.13.2017

#include <vector>
#include <cassert>
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TDirectory.h"

using namespace std;

static const Bool_t   UseVarBins = true;
static const Bool_t   DoTrgNorm  = false;
static const Double_t PiValue    = TMath::Pi();

// filepaths
static const TString iPath("pp200r9pt5ff.chargedJetComp.et920pt0230vz55pi0.r02rm1chrg.d26m3y2019.root");
static const TString oPath("pp200r9pt5ff.chargedJetComp.et920pt0230vz55pi0.r02a005rm1chrg.d27m3y2019.root");

// parameters
static const Double_t MaxVz     = 55.;
static const Double_t MinTrgEt  = 9.;
static const Double_t MaxTrgEt  = 20.;
static const Double_t MaxTrgEta = 0.9;
static const Double_t MinPi0Tsp = 0.;
static const Double_t MaxPi0Tsp = 0.08;
static const Double_t MinGamTsp = 0.2;
static const Double_t MaxGamTsp = 0.6;
static const Double_t MinHadTsp = 3.;

// jet parameters
static const Double_t Rjet     = 0.2;
static const Double_t MinArea  = 0.05;  // R02: 0.05, R03: 0.2, R04: 0.35, R05: 0.65, R07: 1.2
static const Double_t MinJetPt = 0.2;
static const Double_t RecoilDf = PiValue / 4.;



void ReadJetTreeDst(Bool_t inBatchMode=false) {

  gErrorIgnoreLevel = kFatal;
  cout << "\nBeginning reading script!" << endl;

  TFile *fIn  = new TFile(iPath, "read");
  TFile *fOut = new TFile(oPath, "recreate");
  if (!fIn) {
    cerr << "PANIC: input file could not be opened!" << endl;
    assert(fIn);
  }

  // get trees
  TTree *JetTree = 0;
  if (fIn) {
    fIn -> GetObject("JetTree", JetTree);
  }
  if (!JetTree) {
    cerr << "PANIC: 'JetTree' not grabbed!" << endl;
    assert(JetTree);
  }

  // declare event leaves
  cout << "  Setting branch addresses..." << endl;
  Int_t    EventIndex    = 0;
  Int_t    RunId         = 0;
  Int_t    NJets         = 0;
  Double_t PartonicPt    = 0.;
  Double_t Refmult       = 0.;
  Double_t TSP           = 0.;
  Double_t TrgEta        = 0.;
  Double_t TrgPhi        = 0.;
  Double_t TrgEt         = 0.;
  Double_t Rho           = 0.;
  Double_t Sigma         = 0.;
  Double_t Vz            = 0.;
  // declare jet leaves
  vector<Double_t> *JetNCons  = 0;
  vector<Double_t> *JetIndex  = 0;
  vector<Double_t> *JetPt     = 0;
  vector<Double_t> *JetPtCorr = 0;
  vector<Double_t> *JetEta    = 0;
  vector<Double_t> *JetPhi    = 0;
  vector<Double_t> *JetE      = 0;
  vector<Double_t> *JetArea   = 0;
  // declare constituent leaves
  vector<vector<Double_t>> *JetConsPt  = 0;
  vector<vector<Double_t>> *JetConsEta = 0;
  vector<vector<Double_t>> *JetConsPhi = 0;
  vector<vector<Double_t>> *JetConsE   = 0;

  // declare branches
  TBranch *bEventIndex = 0;
  TBranch *bRunId      = 0;
  TBranch *bNJets      = 0;
  TBranch *bPartonicPt = 0;
  TBranch *bRefmult    = 0;
  TBranch *bTsp        = 0;
  TBranch *bTrgEta     = 0;
  TBranch *bTrgPhi     = 0;
  TBranch *bTrgEt      = 0;
  TBranch *bRho        = 0;
  TBranch *bSigma      = 0;
  TBranch *bVz         = 0;
  TBranch *bJetNCons   = 0;
  TBranch *bJetIndex   = 0;
  TBranch *bJetPt      = 0;
  TBranch *bJetPtCorr  = 0;
  TBranch *bJetEta     = 0;
  TBranch *bJetPhi     = 0;
  TBranch *bJetE       = 0;
  TBranch *bJetArea    = 0;
  TBranch *bJetConsPt  = 0;
  TBranch *bJetConsEta = 0;
  TBranch *bJetConsPhi = 0;
  TBranch *bJetConsE   = 0;

  // set branches
  JetTree -> SetBranchAddress("eventIndex", &EventIndex, &bEventIndex);
  JetTree -> SetBranchAddress("RunId", &RunId, &bRunId);
  JetTree -> SetBranchAddress("PartonicPt", &PartonicPt, &bPartonicPt);
  JetTree -> SetBranchAddress("Refmult", &Refmult, &bRefmult);
  JetTree -> SetBranchAddress("NJets", &NJets, &bNJets);
  JetTree -> SetBranchAddress("TSP", &TSP, &bTsp);
  JetTree -> SetBranchAddress("TrgEta", &TrgEta, &bTrgEta);
  JetTree -> SetBranchAddress("TrgPhi", &TrgPhi, &bTrgPhi);
  JetTree -> SetBranchAddress("TrgEt", &TrgEt, &bTrgEt);
  JetTree -> SetBranchAddress("Rho", &Rho, &bRho);
  JetTree -> SetBranchAddress("Sigma", &Sigma, &bSigma);
  JetTree -> SetBranchAddress("Vz", &Vz,&bVz);
  JetTree -> SetBranchAddress("JetIndex", &JetIndex, &bJetIndex);
  JetTree -> SetBranchAddress("JetNCons", &JetNCons, &bJetNCons);
  JetTree -> SetBranchAddress("JetPt", &JetPt, &bJetPt);
  JetTree -> SetBranchAddress("JetPtCorr", &JetPtCorr, &bJetPtCorr);
  JetTree -> SetBranchAddress("JetEta", &JetEta, &bJetEta);
  JetTree -> SetBranchAddress("JetPhi",&JetPhi, &bJetPhi); 
  JetTree -> SetBranchAddress("JetE", &JetE, &bJetE); 
  JetTree -> SetBranchAddress("JetArea",&JetArea, &bJetArea);
  JetTree -> SetBranchAddress("JetConsPt", &JetConsPt, &bJetConsPt);
  JetTree -> SetBranchAddress("JetConsEta", &JetConsEta, &bJetConsEta);
  JetTree -> SetBranchAddress("JetConsPhi", &JetConsPhi, &bJetConsPhi);
  JetTree -> SetBranchAddress("JetConsE", &JetConsE, &bJetConsE);

  // declare histograms [0 for pi0 trigger, 1 for gamma trigger, 2 for hadron]
  cout << "  Creating histograms..." << endl;
  TH1D *hRefmult[3];
  TH1D *hNumJets[3];
  TH1D *hRho[3];
  TH1D *hSigma[3];
  TH1D *hJetArea[3];
  TH1D *hJetEta[3];
  TH1D *hAllPtRaw[3];
  TH1D *hJetPtRaw[3];
  TH1D *hAllDeltaPhi[3];
  TH1D *hBinDeltaPhi[3][3];
  TH1D *hJetDeltaPhi[3];
  TH1D *hJetPtCorr[3];
  TH2D *hJetPtVsDf[3];
  TH2D *hJetPtVsA[3];

  // "fine" variable binning scheme
  const UInt_t   pTnumVar(37);
  const Double_t pTbinVar[38] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 14.5, 17., 20., 23., 26.5, 30., 34., 38., 42.5, 47., 52., 57., 62.5, 68., 74.5, 81.};
  // "big" variable binning scheme
  //const UInt_t   pTnumVar(26);
  //const Double_t pTbinVar[27] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 5., 8., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  // 9 - 11 GeV gamma binning scheme
  //const UInt_t   pTnumVar(28);
  //const Double_t pTbinVar[29] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 12., 17., 23., 30., 38., 47., 57, 68., 81.};
  // 11 - 15 GeV gamma binning scheme
  //const UInt_t   pTnumVar(29);
  //const Double_t pTbinVar[30] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 17., 23., 30., 38., 47., 57., 68., 81.};

  // non-variable binning
  const Int_t    nM  = 200;
  const Int_t    nN  = 50;
  const Int_t    nP  = 50;
  const Int_t    nA  = 500;
  const Int_t    nH  = 100;
  const Int_t    nDf = 360;
  const Int_t    nPt = 240;
  const Double_t m1  = 0.;
  const Double_t m2  = 200.;
  const Double_t n1  = 0.;
  const Double_t n2  = 50.;
  const Double_t p1  = 0.;
  const Double_t p2  = 5.;
  const Double_t a1  = 0.;
  const Double_t a2  = 5.;
  const Double_t h1  = -5.;
  const Double_t h2  = 5.;
  const Double_t dF1 = -1. * (PiValue / 2.);
  const Double_t dF2 = 3. * (PiValue / 2.);
  const Double_t pT1 = -20.;
  const Double_t pT2 = 100.;
  // pi0 histograms
  hRefmult[0]        = new TH1D("hRefmultP", "Refmult, #pi^{0} trigger; refmult; counts", nM, m1, m2);
  hNumJets[0]        = new TH1D("hNumJetsP", "No. of jets, #pi^{0} trigger; N_{jet}; counts", nN, n1, n2);
  hRho[0]            = new TH1D("hRhoP", "Rho, #pi^{0} trigger; #rho; counts", nP, p1, p2);
  hSigma[0]          = new TH1D("hSigmaP", "Sigma, #pi^{0} trigger; #sigma; counts", nP, p1, p2);
  hJetArea[0]        = new TH1D("hJetAreaP", "Recoil jet area, #pi^{0} trigger; A_{jet}; (1/N_{trg}) dN_{jet}/dA_{jet}", nA, a1, a2);
  hJetEta[0]         = new TH1D("hJetEtaP", "Recoil jet #eta, #pi^{0} trigger; #eta; (1/N_{trg}) dN_{jet}/d#eta", nH, h1, h2);
  hAllDeltaPhi[0]    = new TH1D("hAllDeltaPhiP", "All jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[0][0] = new TH1D("hBinDeltaPhiP_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[0][1] = new TH1D("hBinDeltaPhiP_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[0][2] = new TH1D("hBinDeltaPhiP_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hJetDeltaPhi[0]    = new TH1D("hJetDeltaPhiP", "Recoil jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  if (UseVarBins) {
    hAllPtRaw[0]       = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[0]       = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[0]      = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[0]       = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
    hJetPtRaw[0]       = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
    hJetPtCorr[0]      = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", nPt, pT1, pT2);
  }
  hJetPtVsDf[0]      = new TH2D("hJetPtVsDfP", "Jet p_{T}^{raw} vs. #Delta#varphi, #pi^{0} trigger; #Delta#varphi; p_{T}^{raw}", nDf, dF1, dF2, nPt, pT1, pT2);
  hJetPtVsA[0]       = new TH2D("hJetPtVsAp", "Jet p_{T}^{raw} vs. Area, #pi^{0} trigger; A_{jet}; p_{T}^{raw}", nA, a1, a2, nPt, pT1, pT2);
  // gamma histograms
  hRefmult[1]        = new TH1D("hRefmultG", "Refmult, #gamma^{rich} trigger; refmult; counts", nM, m1, m2);
  hNumJets[1]        = new TH1D("hNumJetsG", "No. of jets, #gamma^{rich} trigger; N_{jet}; counts", nN, n1, n2);
  hRho[1]            = new TH1D("hRhoG", "Rho, #gamma^{rich}; #rho; counts", nP, p1, p2);
  hSigma[1]          = new TH1D("hSigmaG", "Sigma, #gamma^{rich} trigger; #sigma; counts", nP, p1, p2);
  hJetArea[1]        = new TH1D("hJetAreaG", "Recoil jet area, #gamma^{rich} trigger; A_{jet}; (1/N_{trg}) dN_{jet}/dA_{jet}", nA, a1, a2);
  hJetEta[1]         = new TH1D("hJetEtaG", "Recoil jet #eta, #gamma^{rich} trigger; #eta; (1/N_{trg}) dN_{jet}/d#eta", nH, h1, h2);
  hAllDeltaPhi[1]    = new TH1D("hAllDeltaPhiG", "All jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[1][0] = new TH1D("hBinDeltaPhiG_pT02", "#Delta#varphi of #gamma^{rich} jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[1][1] = new TH1D("hBinDeltaPhiG_pT1", "#Delta#varphi of #gamma^{rich} jets with p_{T}^{jet}#in(1, 5) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[1][2] = new TH1D("hBinDeltaPhiG_pT5", "#Delta#varphi of #gamma^{rich} jets with p_{T}^{jet}>5 GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hJetDeltaPhi[1]    = new TH1D("hJetDeltaPhiG", "Recoil jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  if (UseVarBins) {
    hAllPtRaw[1]       = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[1]       = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[1]      = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[1]       = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
    hJetPtRaw[1]       = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
    hJetPtCorr[1]      = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", nPt, pT1, pT2);
  }
  hJetPtVsDf[1]      = new TH2D("hJetPtVsDfG", "Jet p_{T}^{raw} vs. #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; p_{T}^{raw}", nDf, dF1, dF2, nPt, pT1, pT2);
  hJetPtVsA[1]       = new TH2D("hJetPtVsAg", "Jet p_{T}^{raw} vs. Area, #gamma^{rich} trigger; A_{jet}; p_{T}^{raw}", nA, a1, a2, nPt, pT1, pT2);
  // hadron histograms
  hRefmult[2]        = new TH1D("hRefmultH", "Refmult, h^{#pm}; refmult; counts", nM, m1, m2);
  hNumJets[2]        = new TH1D("hNumJetsH", "No. of jets, h^{#pm}; N_{jet}; counts", nN, n1, n2);
  hRho[2]            = new TH1D("hRhoH", "Rho, h^{#pm}; #rho; counts", nP, p1, p2);
  hSigma[2]          = new TH1D("hSigmaH", "Sigma, h^{#pm}; #sigma; counts", nP, p1, p2);
  hJetArea[2]        = new TH1D("hJetAreaH", "Recoil jet area, h^{#pm}; A_{jet}; (1/N_{trg}) dN_{jet}/dA_{jet}", nA, a1, a2);
  hJetEta[2]         = new TH1D("hJetEtaH", "Recoil jet #eta, h^{#pm}; #eta; (1/N_{trg}) dN_{jet}/d#eta", nH, h1, h2);
  hAllDeltaPhi[2]    = new TH1D("hAllDeltaPhiH", "All jet #Delta#varphi, h^{#pm}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[2][0] = new TH1D("hBinDeltaPhiH_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, h^{#pm}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[2][1] = new TH1D("hBinDeltaPhiH_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, h^{#pm}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hBinDeltaPhi[2][2] = new TH1D("hBinDeltaPhiH_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, h^{#pm}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  hJetDeltaPhi[2]    = new TH1D("hJetDeltaPhiH", "Recoil jet #Delta#varphi, h^{#pm}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", nDf, dF1, dF2);
  if (UseVarBins) {
    hAllPtRaw[2]       = new TH1D("hAllPtRawH", "All jet p_{T}^{raw}, h^{#pm}; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[2]       = new TH1D("hJetPtRawH", "Recoil jet p_{T}^{raw}, h^{#pm}; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[2]      = new TH1D("hJetPtCorrH", "Recoil jet p_{T}^{corr}, h^{#pm}; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[2]       = new TH1D("hAllPtRawH", "All jet p_{T}^{raw}, h^{#pm}; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
    hJetPtRaw[2]       = new TH1D("hJetPtRawH", "Recoil jet p_{T}^{raw}, h^{#pm}; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", nPt, pT1, pT2);
    hJetPtCorr[2]      = new TH1D("hJetPtCorrH", "Recoil jet p_{T}^{corr}, h^{#pm}; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", nPt, pT1, pT2);
  }
  hJetPtVsDf[2]      = new TH2D("hJetPtVsDfH", "Jet p_{T}^{raw} vs. #Delta#varphi,h^{#pm}; #Delta#varphi; p_{T}^{raw}", nDf, dF1, dF2, nPt, pT1, pT2);
  hJetPtVsA[2]       = new TH2D("hJetPtVsAh", "Jet p_{T}^{raw} vs. Area, h^{#pm}; A_{jet}; p_{T}^{raw}", nA, a1, a2, nPt, pT1, pT2);
  // errors
  for (Int_t i = 0; i < 3; i++) {
    hRefmult[i]        -> Sumw2();
    hNumJets[i]        -> Sumw2();
    hRho[i]            -> Sumw2();
    hSigma[i]          -> Sumw2();
    hJetArea[i]        -> Sumw2();
    hJetEta[i]         -> Sumw2();
    hAllDeltaPhi[i]    -> Sumw2();
    hBinDeltaPhi[i][0] -> Sumw2();
    hBinDeltaPhi[i][1] -> Sumw2();
    hBinDeltaPhi[i][2] -> Sumw2();
    hJetDeltaPhi[i]    -> Sumw2();
    hAllPtRaw[i]       -> Sumw2();
    hJetPtRaw[i]       -> Sumw2();
    hJetPtCorr[i]      -> Sumw2();
    hJetPtVsDf[i]      -> Sumw2();
    hJetPtVsA[i]       -> Sumw2();
  }

  cout << "  Beginning event loop..." << endl;
  Int_t nEvts = (Int_t) JetTree -> GetEntries();

  // event loop
  Int_t nTrgPi0  = 0;
  Int_t nTrgGam  = 0;
  Int_t nTrgHad  = 0;
  Int_t nTrgTot  = 0;
  Int_t nBytes   = 0;
  Int_t breakVal = 0;
  for (Int_t i = 0; i < nEvts; i++) {

    Int_t bytes = JetTree -> GetEntry(i);
    if (bytes < 0) {
      cerr << "ERROR: problem with event " << i << "...\n" << endl;
      breakVal = 1;
      break;
    }

    nBytes += bytes;
    if (!inBatchMode) {
      cout << "    Processing event " << i+1 << "/" << nEvts << "...\r" << flush;
      if (i+1 == nEvts) cout << endl;
    } else {
      cout << "    Processing event " << i+1 << "/" << nEvts << "..." << endl;
    }

    // event info
    const Double_t vZtrg = Vz;
    const Double_t hTrg  = TrgEta;
    const Double_t fTrg  = TrgPhi;
    const Double_t eTtrg = TrgEt;

    // event cuts
    const Bool_t isInVzTrgCut  = (TMath::Abs(vZtrg) < MaxVz);
    const Bool_t isInEtaTrgCut = (TMath::Abs(hTrg) < MaxTrgEta);
    const Bool_t isInEtTrgCut  = ((eTtrg > MinTrgEt) && (eTtrg < MaxTrgEt));
    if (!isInVzTrgCut || !isInEtaTrgCut || !isInEtTrgCut) continue;

    // determine triggers [0 for pi0, 1 for gamma, 2 for hadron]
    const Bool_t inPi0Range = ((TSP >= MinPi0Tsp) && (TSP <= MaxPi0Tsp));
    const Bool_t inGamRange = ((TSP >= MinGamTsp) && (TSP <= MaxGamTsp));
    const Bool_t inHadRange = (TMath::Abs(TSP) >= MinHadTsp);

    Int_t tspFlag = -1;
    if (inPi0Range) {
      tspFlag = 0;
    } else if (inGamRange) {
      tspFlag = 1;
    } else if (inHadRange) {
      tspFlag = 2;
    }

    // fill event histograms
    if (tspFlag == 0) {
      hRefmult[0] -> Fill(Refmult);
      hNumJets[0] -> Fill(NJets);
      hRho[0]     -> Fill(Rho);
      hSigma[0]   -> Fill(Sigma);
      nTrgPi0++;
    } else if (tspFlag == 1) {
      hRefmult[1] -> Fill(Refmult);
      hNumJets[1] -> Fill(NJets);
      hRho[1]     -> Fill(Rho);
      hSigma[1]   -> Fill(Sigma);
      nTrgGam++;
    } else if (tspFlag == 2) {
      hRefmult[2] -> Fill(Refmult);
      hNumJets[2] -> Fill(NJets);
      hRho[2]     -> Fill(Rho);
      hSigma[2]   -> Fill(Sigma);
      nTrgHad++;
    }
    nTrgTot++;

    Int_t nJets = (Int_t) JetEta -> size();
    if (nJets != NJets) {
      cerr << "WARNING: nJets != NJets..." << endl;
    }

    // jet loop
    for (Int_t j = 0; j < nJets; j++) {

      // jet info
      const Double_t hJet  = JetEta    -> at(j);
      const Double_t pT    = JetPt     -> at(j);
      const Double_t pTc   = JetPtCorr -> at(j);
      const Double_t fJet  = JetPhi    -> at(j);
      const Double_t aJet  = JetArea   -> at(j);

      Double_t dF = fJet - fTrg;
      if (dF < (-1. * (PiValue / 2.))) dF += 2.*PiValue;
      if (dF > (3. * (PiValue / 2.)))  dF -= 2.*PiValue;
      Double_t dFre = fJet - fTrg;
      if (dFre < (0. * PiValue)) dFre += 2.*PiValue;
      if (dFre > (2. * PiValue)) dFre -= 2.*PiValue;
      const Double_t dFcut = abs(dFre - PiValue);

      // fill histograms for all jets
      const Bool_t isInBin1 = ((pT > 0.2) && (pT < 1.));
      const Bool_t isInBin2 = ((pT > 1.) && (pT < 5.));
      const Bool_t isInBin3 = (pT > 5.);
      if (tspFlag == 0) {
        hAllPtRaw[0]    -> Fill(pT);
        hAllDeltaPhi[0] -> Fill(dF);
        hJetPtVsDf[0]   -> Fill(dF, pT);
        hJetPtVsA[0]    -> Fill(aJet, pT);
        if (isInBin1) {
          hBinDeltaPhi[0][0] -> Fill(dF);
        }
        if (isInBin2) {
          hBinDeltaPhi[0][1] -> Fill(dF);
        }
        if (isInBin3) {
          hBinDeltaPhi[0][2] -> Fill(dF);
        }
      } else if (tspFlag == 1) {
        hAllPtRaw[1]    -> Fill(pT);
        hAllDeltaPhi[1] -> Fill(dF);
        hJetPtVsDf[1]   -> Fill(dF, pT);
        hJetPtVsA[1]    -> Fill(aJet, pT);
        if (isInBin1) {
          hBinDeltaPhi[1][0] -> Fill(dF);
        }
        if (isInBin2) {
          hBinDeltaPhi[1][1] -> Fill(dF);
        }
        if (isInBin3) {
          hBinDeltaPhi[1][2] -> Fill(dF);
        }
      } else if (tspFlag == 2) {
        hAllPtRaw[2]    -> Fill(pT);
        hAllDeltaPhi[2] -> Fill(dF);
        hJetPtVsDf[2]   -> Fill(dF, pT);
        hJetPtVsA[2]    -> Fill(aJet, pT);
        if (isInBin1) {
          hBinDeltaPhi[2][0] -> Fill(dF);
        }
        if (isInBin2) {
          hBinDeltaPhi[2][1] -> Fill(dF);
        }
        if (isInBin3) {
          hBinDeltaPhi[2][2] -> Fill(dF);
        }
      }

      // jet cuts
      const Bool_t isInPtJetCut = (pT > MinJetPt);
      const Bool_t isInAjetCut  = (aJet > MinArea);
      const Bool_t isRecoil     = (dFcut < RecoilDf);
      if (!isInPtJetCut || !isInAjetCut || !isRecoil) continue;

      // fill histograms
      if (tspFlag == 0) {
        hJetArea[0]     -> Fill(aJet);
        hJetEta[0]      -> Fill(hJet);
        hJetDeltaPhi[0] -> Fill(dF);
        hJetPtRaw[0]    -> Fill(pT);
        hJetPtCorr[0]   -> Fill(pTc);
      } else if (tspFlag == 1) {
        hJetArea[1]     -> Fill(aJet);
        hJetEta[1]      -> Fill(hJet);
        hJetDeltaPhi[1] -> Fill(dF);
        hJetPtRaw[1]    -> Fill(pT);
        hJetPtCorr[1]   -> Fill(pTc);
      } else if (tspFlag == 2) {
        hJetArea[2]     -> Fill(aJet);
        hJetEta[2]      -> Fill(hJet);
        hJetDeltaPhi[2] -> Fill(dF);
        hJetPtRaw[2]    -> Fill(pT);
        hJetPtCorr[2]   -> Fill(pTc);
      }
    }  // end jet loop
  }  // end event loop

  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    cassert(0);
  }
  else {
    cout << "  Event loop finished!\n"
         << "    nPi0 = " << nTrgPi0 << "\n"
         << "    nGam = " << nTrgGam << "\n"
         << "    nHad = " << nTrgHad
         << endl;
  }

  // normalize histograms
  cout << "  Normalizing histograms..." << endl;
  const Double_t mBin    = (m2 - m1) / nM;
  const Double_t nBin    = (n2 - n1) / nN;
  const Double_t pBin    = (p2 - p1) / nP;
  const Double_t aBin    = (a2 - a1) / nA;
  const Double_t hBin    = (h2 - h1) / nH;
  const Double_t dFbin   = (dF2 - dF1) / nDf;
  const Double_t pTbin   = (pT2 - pT1) / nPt;
  const Double_t hJetBin = 2. * (1. - Rjet);
  for (Int_t i = 0; i < 3; i++) {
    Double_t norm = 1;
    if (DoTrgNorm) {
      if (i == 0) {
        norm = nTrgPi0;
      } else if (i == 1) {
        norm = nTrgGam;
      } else if (i == 2) {
        norm = nTrgHad;
      }
    }
    const Double_t mNorm  = mBin * norm;
    const Double_t nNorm  = nBin * norm;
    const Double_t pNorm  = pBin * norm;
    const Double_t aNorm  = aBin * norm * hJetBin;
    const Double_t hNorm  = hBin * norm * hJetBin;
    const Double_t dFnorm = dFbin * norm * hJetBin;
    const Double_t pTnorm = norm * hJetBin;
    hRefmult[i]        -> Scale(1. / mBin);
    hNumJets[i]        -> Scale(1. / nBin);
    hRho[i]            -> Scale(1. / pBin);
    hSigma[i]          -> Scale(1. / pBin);
    hJetArea[i]        -> Scale(1. / aNorm);
    hJetEta[i]         -> Scale(1. / hNorm);
    hAllDeltaPhi[i]    -> Scale(1. / dFnorm);
    hBinDeltaPhi[i][0] -> Scale(1. / dFnorm);
    hBinDeltaPhi[i][1] -> Scale(1. / dFnorm);
    hBinDeltaPhi[i][2] -> Scale(1. / dFnorm);
    hJetDeltaPhi[i]    -> Scale(1. / dFnorm);
    hAllPtRaw[i]       -> Scale(1. / pTnorm);
    hJetPtRaw[i]       -> Scale(1. / pTnorm);
    hJetPtCorr[i]      -> Scale(1. / pTnorm);
    hJetPtVsDf[i]      -> Scale(1. / (norm * dFbin * pTbin * hJetBin));
    hJetPtVsA[i]       -> Scale(1. / (norm * aBin * pTbin * hJetBin));

    // variable bins
    const UInt_t nPtBins = hAllPtRaw[i] -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nPtBins + 1); iBin++) {
      const Double_t pTwidth  = hAllPtRaw[i]  -> GetBinWidth(iBin);
      const Double_t pTallVal = hAllPtRaw[i]  -> GetBinContent(iBin);
      const Double_t pTallErr = hAllPtRaw[i]  -> GetBinError(iBin);
      const Double_t pTrawVal = hJetPtRaw[i]  -> GetBinContent(iBin);
      const Double_t pTrawErr = hJetPtRaw[i]  -> GetBinError(iBin);
      const Double_t pTcorVal = hJetPtCorr[i] -> GetBinContent(iBin);
      const Double_t pTcorErr = hJetPtCorr[i] -> GetBinError(iBin);
      hAllPtRaw[i]  -> SetBinContent(iBin, pTallVal / pTwidth);
      hAllPtRaw[i]  -> SetBinError(iBin, pTallErr / pTwidth);
      hJetPtRaw[i]  -> SetBinContent(iBin, pTrawVal / pTwidth);
      hJetPtRaw[i]  -> SetBinError(iBin, pTrawErr / pTwidth);
      hJetPtCorr[i] -> SetBinContent(iBin, pTcorVal / pTwidth);
      hJetPtCorr[i] -> SetBinError(iBin, pTcorErr / pTwidth);
    }
  }

  // subtract pi0 background
  const Double_t pFrac = 0.4;
  const Double_t gFrac = 0.6;
  TH1D *hPi0 = (TH1D*) hJetPtCorr[0] -> Clone();
  TH1D *hGam = (TH1D*) hJetPtCorr[1] -> Clone();
  TH1D *hDir = (TH1D*) hJetPtCorr[1] -> Clone();
  hPi0 -> Scale(pFrac);
  hDir -> SetNameTitle("hJetPtCorrD", "Recoil jet p_{T}^{corr}, #gamma^{dir} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}");
  hDir -> Add(hGam, hPi0, 1., -1.);
  hDir -> Scale(1. / gFrac);

  // create directories
  TDirectory *dPi0 = (TDirectory*) fOut -> mkdir("Pi0");
  TDirectory *dGam = (TDirectory*) fOut -> mkdir("Gam");
  TDirectory *dHad = (TDirectory*) fOut -> mkdir("Had");

  // write and close output
  fOut -> cd();
  for (Int_t i = 0; i < 3; i++) {
    if (i == 0) {
      dPi0 -> cd();
    } else if (i == 1) {
      dGam -> cd();
    } else if (i == 2) {
      dHad -> cd();
    }
    hRefmult[i]        -> Write();
    hNumJets[i]        -> Write();
    hRho[i]            -> Write();
    hSigma[i]          -> Write();
    hJetArea[i]        -> Write();
    hJetEta[i]         -> Write();
    hAllDeltaPhi[i]    -> Write();
    hBinDeltaPhi[i][0] -> Write();
    hBinDeltaPhi[i][1] -> Write();
    hBinDeltaPhi[i][2] -> Write();
    hJetDeltaPhi[i]    -> Write();
    hAllPtRaw[i]       -> Write();
    hJetPtRaw[i]       -> Write();
    hJetPtCorr[i]      -> Write();
    if (i == 1) {
      hDir             -> Write();
    }
    hJetPtVsDf[i]      -> Write();
    hJetPtVsA[i]       -> Write();
    fOut -> cd();
  }
  fOut -> Close();

  // close input
  fIn -> cd();
  fIn -> Close();
  cout << "Reading script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
