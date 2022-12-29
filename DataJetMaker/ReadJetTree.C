// 'ReadJetTree.C'
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

static const Bool_t   UseVarBins  = false;
static const Bool_t   DoWeighting = false;
static const Bool_t   DoTrgNorm   = true;
static const Double_t PiValue     = TMath::Pi();

// filepaths
static const TString iPath("output/2021/December2021/pp200r9data.forOffAxisCheck_correctCalc.et920pt0230dca1vz55.r05rm1chrg.d22m12y2021.root");
//static const TString iPath("output/2021/December2021/pp200r9data.forOffAxisCheck_correctCalc.et920pt0230dca1vz55.r02rm1chrg.d22m12y2021.root");
static const TString wPath("input/testingWeights_withUpdatedWeightCalculation.mediumSample_allClusterCuts.et620x1115vz55tsp008pi0.d22m1y2021.root");
static const TString oPath("pp200r9data.minBiasPi0check_pTbinOne.et911r05.d5m3y2022.root");

// parameters
static const Double_t MaxVz     = 55.;
static const Double_t MinTrgEt  = 9.;
static const Double_t MaxTrgEt  = 11.;
static const Double_t MaxTrgEta = 0.9;
static const Double_t MinPi0Tsp = 0.;
static const Double_t MaxPi0Tsp = 0.08;
static const Double_t MinGamTsp = 0.2;
static const Double_t MaxGamTsp = 0.6;

// jet parameters
static const Double_t Rjet     = 0.5;
static const Double_t MinArea  = 0.65;  // R02: 0.05, R03: 0.2, R04: 0.35, R05: 0.65, R07: 1.2
static const Double_t MinJetPt = 0.2;
static const Double_t MaxJetPt = 30.;
static const Double_t RecoilDf = PiValue / 4.;



void ReadJetTree(Bool_t inBatchMode=false) {

  gErrorIgnoreLevel = kFatal;
  cout << "\nBeginning reading script!" << endl;

  TFile *fIn  = new TFile(iPath.Data(), "read");
  TFile *fOut = new TFile(oPath.Data(), "recreate");
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

  // get eTtrg weights (if needed)
  TH1D  *hEtWeights;
  TFile *fWeights;
  if (DoWeighting) {
    fWeights = new TFile(wPath.Data(), "read");
    if (!fWeights) {
      cerr << "PANIC: input file for eTtrg weights couldn't be opened!" << endl;
      assert(fWeight);
    }
    hEtWeights = fWeights -> GetObject("hEtWeights", hEtWeights);
    if (!hEtWeights) {
      cerr << "PANIC: 'hEtWeights' not grabbed!" << endl;
      assert(hEtWeights);
    }
  }

  // declare event leaves
  cout << "  Setting branch addresses..." << endl;
  Int_t    EventIndex    = 0;
  Int_t    RunId         = 0;
  Int_t    NJets         = 0;
  Double_t Refmult       = 0.;
  Double_t TSP           = 0.;
  Double_t TrgEta        = 0.;
  Double_t TrgPhi        = 0.;
  Double_t TrgEt         = 0.;
  Double_t Rho           = 0.;
  Double_t Sigma         = 0.;
  Double_t Vz            = 0.;
  // declare jet leaves
  vector<Double_t> *JetNCons         = 0;
  vector<Double_t> *JetIndex         = 0;
  vector<Double_t> *JetPt            = 0;
  vector<Double_t> *JetPtCorr        = 0;
  vector<Double_t> *JetEta           = 0;
  vector<Double_t> *JetPhi           = 0;
  vector<Double_t> *JetE             = 0;
  vector<Double_t> *JetArea          = 0;
  vector<Double_t> *JetPtOffAxisUp   = 0;
  vector<Double_t> *JetPtOffAxisDown = 0;
  // declare constituent leaves
  vector<vector<Double_t>> *JetConsPt  = 0;
  vector<vector<Double_t>> *JetConsEta = 0;
  vector<vector<Double_t>> *JetConsPhi = 0;
  vector<vector<Double_t>> *JetConsE   = 0;

  // declare branches
  TBranch *bEventIndex       = 0;
  TBranch *bRunId            = 0;
  TBranch *bNJets            = 0;
  TBranch *bRefmult          = 0;
  TBranch *bTsp              = 0;
  TBranch *bTrgEta           = 0;
  TBranch *bTrgPhi           = 0;
  TBranch *bTrgEt            = 0;
  TBranch *bRho              = 0;
  TBranch *bSigma            = 0;
  TBranch *bVz               = 0;
  TBranch *bJetNCons         = 0;
  TBranch *bJetIndex         = 0;
  TBranch *bJetPt            = 0;
  TBranch *bJetPtCorr        = 0;
  TBranch *bJetEta           = 0;
  TBranch *bJetPhi           = 0;
  TBranch *bJetE             = 0;
  TBranch *bJetArea          = 0;
  TBranch *bJetPtOffAxisUp   = 0;
  TBranch *bJetPtOffAxisDown = 0;
  TBranch *bJetConsPt        = 0;
  TBranch *bJetConsEta       = 0;
  TBranch *bJetConsPhi       = 0;
  TBranch *bJetConsE         = 0;

  // set branches
  JetTree -> SetBranchAddress("eventIndex", &EventIndex, &bEventIndex);
  JetTree -> SetBranchAddress("RunId", &RunId, &bRunId);
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
  JetTree -> SetBranchAddress("JetPtOffAxisUp", &JetPtOffAxisUp, &bJetPtOffAxisUp);
  JetTree -> SetBranchAddress("JetPtOffAxisDown", &JetPtOffAxisDown, &bJetPtOffAxisDown);
  JetTree -> SetBranchAddress("JetConsPt", &JetConsPt, &bJetConsPt);
  JetTree -> SetBranchAddress("JetConsEta", &JetConsEta, &bJetConsEta);
  JetTree -> SetBranchAddress("JetConsPhi", &JetConsPhi, &bJetConsPhi);
  JetTree -> SetBranchAddress("JetConsE", &JetConsE, &bJetConsE);

  // declare histograms [0 for pi0 trigger, 1 for gamma trigger, 2 for both]
  cout << "  Creating histograms..." << endl;
  TH1D *hTrgEt[3];
  TH1D *hTrgTsp[3];
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
  //const UInt_t   pTnumVar(37);
  //const Double_t pTbinVar[38] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 14.5, 17., 20., 23., 26.5, 30., 34., 38., 42.5, 47., 52., 57., 62.5, 68., 74.5, 81.};
  // "big" variable binning scheme
  //const UInt_t   pTnumVar(26);
  //const Double_t pTbinVar[27] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 5., 8., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  // 9 - 11 GeV gamma binning scheme
  //const UInt_t   pTnumVar(28);
  //const Double_t pTbinVar[29] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 12., 17., 23., 30., 38., 47., 57, 68., 81.};
  // 11 - 15 GeV gamma binning scheme
  //const UInt_t   pTnumVar(29);
  //const Double_t pTbinVar[30] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 4., 5., 6.5, 8., 10., 12., 17., 23., 30., 38., 47., 57., 68., 81.};
  // "huge" variable binning scheme
  const UInt_t   pTnumVar(21);
  const Double_t pTbinVar[22] = {-11., -8., -6., -4., -2., -1., 0., 1., 2., 4., 6., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60.};

  // non-variable binning
  const UInt_t   eTnum(200);
  const UInt_t   tNum(1000);
  const UInt_t   mNum(200);
  const UInt_t   nNum(50);
  const UInt_t   pNum(50);
  const UInt_t   aNum(500);
  const UInt_t   hNum(100);
  const UInt_t   dFnum(360);
  //const UInt_t   pTnum(140);
  const Int_t    pTnum(100);
  const Double_t eTbins[2] = {0., 100.};
  const Double_t tBins[2]  = {0., 10.};
  const Double_t mBins[2]  = {0., 200.};
  const Double_t nBins[2]  = {0., 50.};
  const Double_t pBins[2]  = {0., 5.};
  const Double_t aBins[2]  = {0., 5.};
  const Double_t hBins[2]  = {-5., 5.};
  const Double_t dFbins[2] = {-1. * (PiValue / 2.), 3. * (PiValue / 2.)};
  //const Double_t pTbins[2] = {-20., 50.};
  const Double_t pTbins[2] = {-50., 50.};
  // pi0 histograms
  hTrgEt[0]          = new TH1D("hTrgEtP", "Trigger E_{T}, #pi^{0} trigger; E_{T}^{trg}; counts", eTnum, eTbins[0], eTbins[1]);
  hTrgTsp[0]         = new TH1D("hTrgTspP", "Trigger TSP, #pi^{0} trigger; TSP; counts", tNum, tBins[0], tBins[1]);
  hRefmult[0]        = new TH1D("hRefmultP", "Refmult, #pi^{0} trigger; refmult; counts", mNum, mBins[0], mBins[1]);
  hNumJets[0]        = new TH1D("hNumJetsP", "No. of jets, #pi^{0} trigger; N_{jet}; counts", nNum, nBins[0], nBins[1]);
  hRho[0]            = new TH1D("hRhoP", "Rho, #pi^{0} trigger; #rho; counts", pNum, pBins[0], pBins[1]);
  hSigma[0]          = new TH1D("hSigmaP", "Sigma, #pi^{0} trigger; #sigma; counts", pNum, pBins[0], pBins[1]);
  hJetArea[0]        = new TH1D("hJetAreaP", "Recoil jet area, #pi^{0} trigger; A_{jet}; (1/N_{trg}) dN_{jet}/dA_{jet}", aNum, aBins[0], aBins[1]);
  hJetEta[0]         = new TH1D("hJetEtaP", "Recoil jet #eta, #pi^{0} trigger; #eta; (1/N_{trg}) dN_{jet}/d#eta", hNum, hBins[0], hBins[1]);
  hAllDeltaPhi[0]    = new TH1D("hAllDeltaPhiP", "All jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[0][0] = new TH1D("hBinDeltaPhiP_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[0][1] = new TH1D("hBinDeltaPhiP_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[0][2] = new TH1D("hBinDeltaPhiP_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hJetDeltaPhi[0]    = new TH1D("hJetDeltaPhiP", "Recoil jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  if (UseVarBins) {
    hAllPtRaw[0]       = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[0]       = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[0]      = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[0]       = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtRaw[0]       = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtCorr[0]      = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnum, pTbins[0], pTbins[1]);
  }
  hJetPtVsDf[0]      = new TH2D("hJetPtVsDfP", "Jet p_{T}^{raw} vs. #Delta#varphi, #pi^{0} trigger; #Delta#varphi; p_{T}^{raw}", dFnum, dFbins[0], dFbins[1], pTnum, pTbins[0], pTbins[1]);
  hJetPtVsA[0]       = new TH2D("hJetPtVsAp", "Jet p_{T}^{raw} vs. Area, #pi^{0} trigger; A_{jet}; p_{T}^{raw}", aNum, aBins[0], aBins[1], pTnum, pTbins[0], pTbins[1]);
  // gamma histograms
  hTrgEt[1]          = new TH1D("hTrgEtG", "Trigger E_{T}, #gamma^{rich} trigger; E_{T}^{trg}; counts", eTnum, eTbins[0], eTbins[1]);
  hTrgTsp[1]         = new TH1D("hTrgTspG", "Trigger TSP, #gamma^{rich} trigger; TSP; counts", tNum, tBins[0], tBins[1]);
  hRefmult[1]        = new TH1D("hRefmultG", "Refmult, #gamma^{rich} trigger; refmult; counts", mNum, mBins[0], mBins[1]);
  hNumJets[1]        = new TH1D("hNumJetsG", "No. of jets, #gamma^{rich} trigger; N_{jet}; counts", nNum, nBins[0], nBins[1]);
  hRho[1]            = new TH1D("hRhoG", "Rho, #gamma^{rich}; #rho; counts", pNum, pBins[0], pBins[1]);
  hSigma[1]          = new TH1D("hSigmaG", "Sigma, #gamma^{rich} trigger; #sigma; counts", pNum, pBins[0], pBins[1]);
  hJetArea[1]        = new TH1D("hJetAreaG", "Recoil jet area, #gamma^{rich} trigger; A_{jet}; (1/N_{trg}) dN_{jet}/dA_{jet}", aNum, aBins[0], aBins[1]);
  hJetEta[1]         = new TH1D("hJetEtaG", "Recoil jet #eta, #gamma^{rich} trigger; #eta; (1/N_{trg}) dN_{jet}/d#eta", hNum, hBins[0], hBins[1]);
  hAllDeltaPhi[1]    = new TH1D("hAllDeltaPhiG", "All jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[1][0] = new TH1D("hBinDeltaPhiG_pT02", "#Delta#varphi of #gamma^{rich} jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[1][1] = new TH1D("hBinDeltaPhiG_pT1", "#Delta#varphi of #gamma^{rich} jets with p_{T}^{jet}#in(1, 5) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[1][2] = new TH1D("hBinDeltaPhiG_pT5", "#Delta#varphi of #gamma^{rich} jets with p_{T}^{jet}>5 GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hJetDeltaPhi[1]    = new TH1D("hJetDeltaPhiG", "Recoil jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  if (UseVarBins) {
    hAllPtRaw[1]       = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[1]       = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[1]      = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[1]       = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtRaw[1]       = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtCorr[1]      = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnum, pTbins[0], pTbins[1]);
  }
  hJetPtVsDf[1]      = new TH2D("hJetPtVsDfG", "Jet p_{T}^{raw} vs. #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; p_{T}^{raw}", dFnum, dFbins[0], dFbins[1], pTnum, pTbins[0], pTbins[1]);
  hJetPtVsA[1]       = new TH2D("hJetPtVsAg", "Jet p_{T}^{raw} vs. Area, #gamma^{rich} trigger; A_{jet}; p_{T}^{raw}", aNum, aBins[0], aBins[1], pTnum, pTbins[0], pTbins[1]);
  // all histograms
  hTrgEt[2]          = new TH1D("hTrgEtA", "Trigger E_{T}, all; E_{T}^{trg}; counts", eTnum, eTbins[0], eTbins[1]);
  hTrgTsp[2]         = new TH1D("hTrgTspA", "Trigger TSP, all; TSP; counts", tNum, tBins[0], tBins[1]);
  hRefmult[2]        = new TH1D("hRefmultA", "Refmult, all; refmult; counts", mNum, mBins[0], mBins[1]);
  hNumJets[2]        = new TH1D("hNumJetsA", "No. of jets, all; N_{jet}; counts", nNum, nBins[0], nBins[1]);
  hRho[2]            = new TH1D("hRhoA", "Rho, all; #rho; counts", pNum, pBins[0], pBins[1]);
  hSigma[2]          = new TH1D("hSigmaA", "Sigma, all; #sigma; counts", pNum, pBins[0], pBins[1]);
  hJetArea[2]        = new TH1D("hJetAreaA", "Recoil jet area, all; A_{jet}; (1/N_{trg}) dN_{jet}/dA_{jet}", aNum, aBins[0], aBins[1]);
  hJetEta[2]         = new TH1D("hJetEtaA", "Recoil jet #eta, all; #eta; (1/N_{trg}) dN_{jet}/d#eta", hNum, hBins[0], hBins[1]);
  hAllDeltaPhi[2]    = new TH1D("hAllDeltaPhiA", "All jet #Delta#varphi, all; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[2][0] = new TH1D("hBinDeltaPhiA_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, all; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[2][1] = new TH1D("hBinDeltaPhiA_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, all; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[2][2] = new TH1D("hBinDeltaPhiA_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, all; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hJetDeltaPhi[2]    = new TH1D("hJetDeltaPhiA", "Recoil jet #Delta#varphi, all; #Delta#varphi; (1/N_{trg}) dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  if (UseVarBins) {
    hAllPtRaw[2]       = new TH1D("hAllPtRawA", "All jet p_{T}^{raw}, all; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[2]       = new TH1D("hJetPtRawA", "Recoil jet p_{T}^{raw}, all; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[2]      = new TH1D("hJetPtCorrA", "Recoil jet p_{T}^{corr}, all; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[2]       = new TH1D("hAllPtRawA", "All jet p_{T}^{raw}, all; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtRaw[2]       = new TH1D("hJetPtRawA", "Recoil jet p_{T}^{raw}, all; p_{T}^{raw}; (1/N_{trg}) dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtCorr[2]      = new TH1D("hJetPtCorrA", "Recoil jet p_{T}^{corr}, all; p_{T}^{corr}; (1/N_{trg}) dN_{jet}/dp_{T}^{corr}", pTnum, pTbins[0], pTbins[1]);
  }
  hJetPtVsDf[2]      = new TH2D("hJetPtVsDfA", "Jet p_{T}^{raw} vs. #Delta#varphi,all; #Delta#varphi; p_{T}^{raw}", dFnum, dFbins[0], dFbins[1], pTnum, pTbins[0], pTbins[1]);
  hJetPtVsA[2]       = new TH2D("hJetPtVsAa", "Jet p_{T}^{raw} vs. Area, all; A_{jet}; p_{T}^{raw}", aNum, aBins[0], aBins[1], pTnum, pTbins[0], pTbins[1]);
  // errors
  for (Int_t i = 0; i < 3; i++) {
    hTrgEt[i]          -> Sumw2();
    hTrgTsp[i]         -> Sumw2();
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
  const UInt_t   pTnumVar(21);
  const Double_t pTbinVar[22] = {-11., -8., -6., -4., -2., -1., 0., 1., 2., 4., 6., 8., 11., 15., 20., 25., 30., 35., 40., 45., 50., 60.};

  // TEST [12.08.2021, 12.12.2021, 12.13.2021, 12.31.2021]
  TH1D *hJetUeP;
  TH1D *hJetUeG;
  TH1D *hJetOaP;
  TH1D *hJetOaG;
  TH1D *hJetRhoP;
  TH1D *hJetRhoG;
  TH1D *hJetRhoAP;
  TH1D *hJetRhoAG;
  TH1D *hJetUeAvgP;
  TH1D *hJetUeAvgG;
  TH1D *hJetPtSubP;
  TH1D *hJetPtSubG;
  TH1D *hJetPtUeSubP;
  TH1D *hJetPtUeSubG;
  TH1D *hJetPtOaSubP;
  TH1D *hJetPtOaSubG;
  TH1D *hJetOaDensityP;
  TH1D *hJetOaDensityG;
  TH1D *hJetUeAvgDensityP;
  TH1D *hJetUeAvgDensityG;
  TH2D *hJetPtVsRhoAP;
  TH2D *hJetPtVsRhoAG;

  const Int_t    NumPt = 140;
  const Double_t PtMin = -20.;
  const Double_t PtMax = 50.;
  if (UseVarBins) {
    hJetUeP           = new TH1D("hJetPtUE_pi0", "UE jet pT [pi0]", pTnumVar, pTbinVar);
    hJetUeG           = new TH1D("hJetPtUE_gam", "UE jet pT [gam]", pTnumVar, pTbinVar);
    hJetOaP           = new TH1D("hJetPtOA_pi0", "recoil jet dPtOA = aJet * (sigUp + sigDown) / 2 [pi0]", pTnumVar, pTbinVar);
    hJetOaG           = new TH1D("hJetPtOA_gam", "recoil jet dPtOA = aJet * (sigUp + sigDown) / 2 [gam]", pTnumVar, pTbinVar);
    hJetRhoP          = new TH1D("hJetRho_pi0", "rho [pi0]", pTnumVar, pTbinVar);
    hJetRhoG          = new TH1D("hJetRo_gam", "rho [gam]", pTnumVar, pTbinVar);
    hJetRhoAP         = new TH1D("hJetRhoA_pi0", "recoil jet (rho * A) [pi0]", pTnumVar, pTbinVar);
    hJetRhoAG         = new TH1D("hJetRhoA_gam", "recoil jet (rho * A) [gam]", pTnumVar, pTbinVar);
    hJetUeAvgP        = new TH1D("hJetUeAvg_pi0", "average UE jet pT [pi0]", pTnumVar, pTbinVar);
    hJetUeAvgG        = new TH1D("hJetUeAvg_gam", "average UE jet pT [gam]", pTnumVar, pTbinVar);
    hJetPtSubP        = new TH1D("hJetPtRhoASub_pi0", "recoil jet pTcorr = pTraw - (rho * A) [pi0]", pTnumVar, pTbinVar);
    hJetPtSubG        = new TH1D("hJetPtRhoASub_gam", "recoil jet pTcorr = pTraw - (rho * A) [gam]", pTnumVar, pTbinVar);
    hJetPtUeSubP      = new TH1D("hJetPtUeSub_pi0", "recoil jet pTcorr = pTraw - pTue [pi0]", pTnumVar, pTbinVar);
    hJetPtUeSubG      = new TH1D("hJetPtUeSub_gam", "recoil jet pTcorr = pTraw - pTue [gam]", pTnumVar, pTbinVar);
    hJetPtOaSubP      = new TH1D("hJetPtOaSub_pi0", "recoil jet pTcorr = pTraw - dPtOA [pi0]", pTnumVar, pTbinVar);
    hJetPtOaSubG      = new TH1D("hJetPtOaSub_gam", "recoil jet pTcorr = pTraw - dPtOA [gam]", pTnumVar, pTbinVar);
    hJetOaDensityP    = new TH1D("hJetOaDensity_pi0", "recoil jet OA density = (sigUp + sigDown) / 2 [pi0]", pTnumVar, pTbinVar);
    hJetOaDensityG    = new TH1D("hJetOaDensity_gam", "recoil jet OA density = (sigUp + sigDown) / 2 [gam]", pTnumVar, pTbinVar);
    hJetUeAvgDensityP = new TH1D("hJetUeAvgDensity_pi0", "average UE jet pT / A [pi0]", pTnumVar, pTbinVar);
    hJetUeAvgDensityG = new TH1D("hJetUeAvgDensity_gam", "average UE jet pT / A [gam]", pTnumVar, pTbinVar);
    hJetPtVsRhoAP     = new TH2D("hJetPtVsRhoA_pi0", "recoil jet pTraw vs. (rho * A) [pi0]; pTraw; (rho * A)", pTnumVar, pTbinVar, pTnumVar, pTbinVar);
    hJetPtVsRhoAG     = new TH2D("hJetPtVsRhoA_gam", "recoil jet pTraw vs. (rho * A) [gam]; pTraw; (rho * A)", pTnumVar, pTbinVar, pTnumVar, pTbinVar);
  } else {
    hJetUeP           = new TH1D("hJetPtUE_pi0", "UE jet pT [pi0]", NumPt, PtMin, PtMax);
    hJetUeG           = new TH1D("hJetPtUE_gam", "UE jet pT [gam]", NumPt, PtMin, PtMax);
    hJetOaP           = new TH1D("hJetPtOA_pi0", "recoil jet dPtOA = aJet * (sigUp + sigDown) / 2 [pi0]", NumPt, PtMin, PtMax);
    hJetOaG           = new TH1D("hJetPtOA_gam", "recoil jet dPtOA = aJet * (sigUp + sigDown) / 2 [gam]", NumPt, PtMin, PtMax);
    hJetRhoP          = new TH1D("hJetRho_pi0", "rho [pi0]", NumPt, PtMin, PtMax);
    hJetRhoG          = new TH1D("hJetRho_gam", "rho [gamma]", NumPt, PtMin, PtMax);
    hJetRhoAP         = new TH1D("hJetRhoA_pi0", "recoil jet (rho * A) [pi0]", NumPt, PtMin, PtMax);
    hJetRhoAG         = new TH1D("hJetRhoA_gam", "recoil jet (rho * A) [gam]", NumPt, PtMin, PtMax);
    hJetUeAvgP        = new TH1D("hJetUeAvg_pi0", "average UE jet pT [pi0]", NumPt, PtMin, PtMax);
    hJetUeAvgG        = new TH1D("hJetUeAvg_gam", "average UE jet pT [gam]", NumPt, PtMin, PtMax);
    hJetPtSubP        = new TH1D("hJetPtRhoASub_pi0", "recoil jet pTcorr = pTraw - (rho * A) [pi0]", NumPt, PtMin, PtMax);
    hJetPtSubG        = new TH1D("hJetPtRhoASub_gam", "recoil jet pTcorr = pTraw - (rho * A) [gam]", NumPt, PtMin, PtMax);
    hJetPtUeSubP      = new TH1D("hJetPtUeSub_pi0", "recoil jet pTcorr = pTraw - pTue [pi0]", NumPt, PtMin, PtMax);
    hJetPtUeSubG      = new TH1D("hJetPtUeSub_gam", "recoil jet pTcorr = pTraw - pTue [gam]", NumPt, PtMin, PtMax);
    hJetPtOaSubP      = new TH1D("hJetPtOaSub_pi0", "recoil jet pTcorr = pTraw - dPtOA [pi0]", NumPt, PtMin, PtMax);
    hJetPtOaSubG      = new TH1D("hJetPtOaSub_gam", "recoil jet pTcorr = pTraw - dPtOA [gam]", NumPt, PtMin, PtMax);
    hJetOaDensityP    = new TH1D("hJetOaDensity_pi0", "recoil jet OA density = (sigUp + sigDown) / 2 [pi0]", NumPt, PtMin, PtMax);
    hJetOaDensityG    = new TH1D("hJetOaDensity_gam", "recoil jet OA density = (sigUp + sigDown) / 2 [gam]", NumPt, PtMin, PtMax);
    hJetUeAvgDensityP = new TH1D("hJetUeAvgDensity_pi0", "average UE jet pT / A [pi0]", NumPt, PtMin, PtMax);
    hJetUeAvgDensityG = new TH1D("hJetUeAvgDensity_gam", "average UE jet pT / A [gam]", NumPt, PtMin, PtMax);
    hJetPtVsRhoAP     = new TH2D("hJetPtVsRhoA_pi0", "recoil jet pTraw vs. (rho * A) [pi0]; pTraw; (rho * A)", NumPt, PtMin, PtMax, NumPt, PtMin, PtMax);
    hJetPtVsRhoAG     = new TH2D("hJetPtVsRhoA_gam", "recoil jet pTraw vs. (rho * A) [gam]; pTraw; (rho * A)", NumPt, PtMin, PtMax, NumPt, PtMin, PtMax);
  }
  hJetUeP           -> Sumw2();
  hJetUeG           -> Sumw2();
  hJetOaP           -> Sumw2();
  hJetOaG           -> Sumw2();
  hJetRhoP          -> Sumw2();
  hJetRhoG          -> Sumw2();
  hJetRhoAP         -> Sumw2();
  hJetRhoAG         -> Sumw2();
  hJetUeAvgP        -> Sumw2();
  hJetUeAvgG        -> Sumw2();
  hJetPtSubP        -> Sumw2();
  hJetPtSubG        -> Sumw2();
  hJetPtUeSubP      -> Sumw2();
  hJetPtUeSubG      -> Sumw2();
  hJetPtOaSubP      -> Sumw2();
  hJetPtOaSubG      -> Sumw2();
  hJetOaDensityP    -> Sumw2();
  hJetOaDensityG    -> Sumw2();
  hJetUeAvgDensityP -> Sumw2();
  hJetUeAvgDensityG -> Sumw2();
  hJetPtVsRhoAP     -> Sumw2();
  hJetPtVsRhoAG     -> Sumw2();

  Int_t nEvts = (Int_t) JetTree -> GetEntries();
  cout << "  Beginning event loop..." << endl;

  // event loop
  Int_t nTrgPi0  = 0;
  Int_t nTrgGam  = 0;
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

    // determine triggers [0 for pi0, 1 for gamma, 2 for both]
    const Bool_t inPi0Range = ((TSP >= MinPi0Tsp) && (TSP <= MaxPi0Tsp));
    const Bool_t inGamRange = ((TSP >= MinGamTsp) && (TSP <= MaxGamTsp));
    const Bool_t isInTspCut = (inPi0Range || inGamRange);
    if (!isInTspCut) continue;

    Int_t tspFlag = -1;
    if (inPi0Range) {
      tspFlag = 0;
    } else if (inGamRange) {
      tspFlag = 1;
    }

    // fill event histograms
    if (tspFlag == 0) {
      hTrgEt[0]   -> Fill(eTtrg);
      hTrgTsp[0]  -> Fill(TSP);
      hRefmult[0] -> Fill(Refmult);
      hNumJets[0] -> Fill(NJets);
      hRho[0]     -> Fill(Rho);
      hSigma[0]   -> Fill(Sigma);
      nTrgPi0++;
    } else if (tspFlag == 1) {
      hTrgEt[1]   -> Fill(eTtrg);
      hTrgTsp[1]  -> Fill(TSP);
      hRefmult[1] -> Fill(Refmult);
      hNumJets[1] -> Fill(NJets);
      hRho[1]     -> Fill(Rho);
      hSigma[1]   -> Fill(Sigma);
      nTrgGam++;
    }
    hTrgEt[2]   -> Fill(eTtrg);
    hTrgTsp[2]  -> Fill(TSP);
    hRefmult[2] -> Fill(Refmult);
    hNumJets[2] -> Fill(NJets);
    hRho[2]     -> Fill(Rho);
    hSigma[2]   -> Fill(Sigma);
    nTrgTot++;

    // grab eT weights
    UInt_t   iEtWeight(0);
    Double_t eTweight(1.);
    if (DoWeighting) {
      iEtWeight = hEtWeights -> FindBin(eTtrg);
      eTweight  = hEtWeights -> GetBinContent(iEtWeight);
    }

    Int_t nJets = (Int_t) JetEta -> size();
    if (nJets != NJets) {
      cerr << "WARNING: nJets != NJets..." << endl;
    }

    // UE jet loop
    UInt_t   nUeJet(0);
    Double_t avgUePt(0.);
    Double_t avgUeRho(0.);
    for (Int_t j = 0; j < nJets; j++) {

      // jet info
      const Double_t hJet = JetEta    -> at(j);
      const Double_t pT   = JetPt     -> at(j);
      const Double_t pTc  = JetPtCorr -> at(j);
      const Double_t fJet = JetPhi    -> at(j);
      const Double_t aJet = JetArea   -> at(j);
      const Double_t rJet = pT / aJet;
      const Double_t dJet = Rho * aJet;
      const Double_t pTs  = pT - dJet;

      Double_t dFue = fJet - fTrg;
      if (dFue < (0. * PiValue)) dFue += 2.*PiValue;
      if (dFue > (2. * PiValue)) dFue -= 2.*PiValue;

      // jet cuts
      const Bool_t isInAjetCut = (aJet > MinArea);
      const Bool_t isInLowUE   = ((dFue > (PiValue / 4.)) && (dFue < (PiValue / 2.)));
      const Bool_t isInHighUE  = ((dFue > (3. * PiValue / 2.)) && (dFue < (7. * PiValue / 4.)));
      const Bool_t isInUE      = (isInLowUE || isInHighUE);
      if (!isInAjetCut || !isInUE) continue;

      // histograms and increment sums
      if (tspFlag == 0) {
        hJetUeP -> Fill(pT);
      } else if (tspFlag == 1) {
        hJetUeG -> Fill(pT);
      }
      avgUePt  += pT;
      avgUeRho += rJet;
      ++nUeJet;
    }  // end UE jet loop

    // calculate average UE energy
    if (nUeJet > 0) {
      avgUePt  = avgUePt / nUeJet;
      avgUeRho = avgUeRho / nUeJet;
    } else {
      avgUePt  = 0.;
      avgUeRho = 0.;
    }

    // fill histograms
    if (tspFlag == 0) {
      hJetRhoP          -> Fill(Rho);
      hJetUeAvgP        -> Fill(avgUePt);
      hJetUeAvgDensityP -> Fill(avgUeRho);
    } else if (tspFlag == 1) {
      hJetRhoG          -> Fill(Rho);
      hJetUeAvgG        -> Fill(avgUePt);
      hJetUeAvgDensityG -> Fill(avgUeRho);
    }

    // recoil jet loop
    for (Int_t j = 0; j < nJets; j++) {

      // jet info
      const Double_t hJet    = JetEta           -> at(j);
      const Double_t pT      = JetPt            -> at(j);
      const Double_t pTc     = JetPtCorr        -> at(j);
      const Double_t fJet    = JetPhi           -> at(j);
      const Double_t aJet    = JetArea          -> at(j);
      const Double_t sigUp   = JetPtOffAxisUp   -> at(j);
      const Double_t sigDown = JetPtOffAxisDown -> at(j);
      const Double_t dJet    = Rho * aJet;
      const Double_t sigAvg  = (sigUp + sigDown) / 2.;
      const Double_t dJetOA  = sigAvg * aJet;
      const Double_t pTs     = pT - dJet;
      const Double_t pTue    = pT - avgUePt;
      const Double_t pToa    = pT - dJetOA;

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
      }
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

      // jet cuts
      const Bool_t isInPtJetCut = ((pT > MinJetPt) && (pT < MaxJetPt));
      const Bool_t isInAjetCut  = (aJet > MinArea);
      const Bool_t isRecoil     = (dFcut < RecoilDf);
      if (!isInPtJetCut || !isInAjetCut || !isRecoil) continue;

      // fill pi0/gamma histograms
      if (tspFlag == 0) {
        if (DoWeighting) {
          hJetArea[0]     -> Fill(aJet, eTweight);
          hJetEta[0]      -> Fill(hJet, eTweight);
          hJetDeltaPhi[0] -> Fill(dF, eTweight);
          hJetPtRaw[0]    -> Fill(pT, eTweight);
          hJetPtCorr[0]   -> Fill(pTc, eTweight);
        } else {
          hJetArea[0]     -> Fill(aJet);
          hJetEta[0]      -> Fill(hJet);
          hJetDeltaPhi[0] -> Fill(dF);
          hJetPtRaw[0]    -> Fill(pT);
          hJetPtCorr[0]   -> Fill(pTc);
        }
        hJetOaP        -> Fill(dJetOA);
        hJetRhoAP      -> Fill(dJet);
        hJetPtSubP     -> Fill(pTs);
        hJetPtUeSubP   -> Fill(pTue);
        hJetPtOaSubP   -> Fill(pToa);
        hJetOaDensityP -> Fill(sigAvg);
        hJetPtVsRhoAP  -> Fill(pT, dJet);
      } else if (tspFlag == 1) {
        if (DoWeighting) {
          hJetArea[1]     -> Fill(aJet, eTweight);
          hJetEta[1]      -> Fill(hJet, eTweight);
          hJetDeltaPhi[1] -> Fill(dF, eTweight);
          hJetPtRaw[1]    -> Fill(pT, eTweight);
          hJetPtCorr[1]   -> Fill(pTc, eTweight);
        } else {
          hJetArea[1]     -> Fill(aJet);
          hJetEta[1]      -> Fill(hJet);
          hJetDeltaPhi[1] -> Fill(dF);
          hJetPtRaw[1]    -> Fill(pT);
          hJetPtCorr[1]   -> Fill(pTc);
        }
        hJetOaG        -> Fill(dJetOA);
        hJetRhoAG      -> Fill(dJet);
        hJetPtSubG     -> Fill(pTs);
        hJetPtUeSubG   -> Fill(pTue);
        hJetPtOaSubG   -> Fill(pToa);
        hJetOaDensityG -> Fill(sigAvg);
        hJetPtVsRhoAG  -> Fill(pT, dJet);
      }

      // fill all histograms
      if (DoWeighting) {
        hJetArea[2]     -> Fill(aJet, eTweight);
        hJetEta[2]      -> Fill(hJet, eTweight);
        hJetDeltaPhi[2] -> Fill(dF, eTweight);
        hJetPtRaw[2]    -> Fill(pT, eTweight);
        hJetPtCorr[2]   -> Fill(pTc, eTweight);
      } else {
        hJetArea[2]     -> Fill(aJet);
        hJetEta[2]      -> Fill(hJet);
        hJetDeltaPhi[2] -> Fill(dF);
        hJetPtRaw[2]    -> Fill(pT);
        hJetPtCorr[2]   -> Fill(pTc);
      }
    }  // end recoil jet loop
  }  // end event loop

  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    cassert(0);
  } else {
    cout << "  Event loop finished!" << endl;
  }

  // normalize histograms
  cout << "  Normalizing histograms..." << endl;
  const Double_t eTbin   = (eTbins[1] - eTbins[0]) / eTnum;
  const Double_t tBin    = (tBins[1] - tBins[0]) / tNum;
  const Double_t mBin    = (mBins[1] - mBins[0]) / mNum;
  const Double_t nBin    = (nBins[1] - nBins[0]) / nNum;
  const Double_t pBin    = (pBins[1] - pBins[0]) / pNum;
  const Double_t aBin    = (aBins[1] - aBins[0]) / aNum;
  const Double_t hBin    = (hBins[1] - hBins[0]) / hNum;
  const Double_t dFbin   = (dFbins[1] - dFbins[0]) / dFnum;
  const Double_t pTbin   = (pTbins[1] - pTbins[0]) / pTnum;
  const Double_t hJetBin = 2. * (1. - Rjet);
  for (Int_t i = 0; i < 3; i++) {
    Double_t norm = 1;
    if (DoTrgNorm) {
      if (i == 0) {
        norm = nTrgPi0;
      } else if (i == 1) {
        norm = nTrgGam;
      } else if (i == 2) {
        norm = nTrgTot;
      }
    }
    const Double_t eTnorm = eTbin * norm;
    const Double_t tNorm  = tBin * norm;
    const Double_t mNorm  = mBin * norm;
    const Double_t nNorm  = nBin * norm;
    const Double_t pNorm  = pBin * norm;
    const Double_t aNorm  = aBin * norm * hJetBin;
    const Double_t hNorm  = hBin * norm * hJetBin;
    const Double_t dFnorm = dFbin * norm * hJetBin;
    const Double_t pTnorm = norm * hJetBin;
    hTrgEt[i]          -> Scale(1. / eTnorm);
    hTrgTsp[i]         -> Scale(1. / tNorm);
    hRefmult[i]        -> Scale(1. / mNorm);
    hNumJets[i]        -> Scale(1. / nNorm);
    hRho[i]            -> Scale(1. / pNorm);
    hSigma[i]          -> Scale(1. / pNorm);
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

  // TEST [12.08.2021, 12.12.2021, 12.13.2021, 12.31.2021]
  const Double_t corrPtBin = (PtMax - PtMin) / NumPt;
  const Double_t rhoNormP  = (Double_t) nTrgPi0;
  const Double_t rhoNormG  = (Double_t) nTrgGam;
  const Double_t etaNormP  = nTrgPi0 * hJetBin;
  const Double_t etaNormG  = nTrgGam * hJetBin;
  const Double_t avgNormP  = nTrgPi0 * corrPtBin;
  const Double_t avgNormG  = nTrgGam * corrPtBin;
  const Double_t corrNormP = nTrgPi0 * corrPtBin * hJetBin;
  const Double_t corrNormG = nTrgGam * corrPtBin * hJetBin;
  if (UseVarBins) {
    const UInt_t nPtVar = hJetUeAvgP -> GetNbinsX();
    for (UInt_t iVarBin = 1; iVarBin < (nPtVar + 1); iVarBin++) {
      const Double_t binWidth      = hJetUeAvgP        -> GetBinWidth(iVarBin);
      const Double_t valUeAvgP     = hJetUeAvgP        -> GetBinContent(iVarBin);
      const Double_t valUeAvgG     = hJetUeAvgG        -> GetBinContent(iVarBin);
      const Double_t valUeP        = hJetUeP           -> GetBinContent(iVarBin);
      const Double_t valUeG        = hJetUeG           -> GetBinContent(iVarBin);
      const Double_t valOaP        = hJetOaP           -> GetBinContent(iVarBin);
      const Double_t valOaG        = hJetOaG           -> GetBinContent(iVarBin);
      const Double_t valRhoP       = hJetRhoP          -> GetBinContent(iVarBin);
      const Double_t valRhoG       = hJetRhoG          -> GetBinContent(iVarBin);
      const Double_t valRhoAP      = hJetRhoAP         -> GetBinContent(iVarBin);
      const Double_t valRhoAG      = hJetRhoAG         -> GetBinContent(iVarBin);
      const Double_t valPtSubP     = hJetPtSubP        -> GetBinContent(iVarBin);
      const Double_t valPtSubG     = hJetPtSubG        -> GetBinContent(iVarBin);
      const Double_t valPtUeSubP   = hJetPtUeSubP      -> GetBinContent(iVarBin);
      const Double_t valPtUeSubG   = hJetPtUeSubG      -> GetBinContent(iVarBin);
      const Double_t valPtOaSubP   = hJetPtOaSubP      -> GetBinContent(iVarBin);
      const Double_t valPtOaSubG   = hJetPtOaSubG      -> GetBinContent(iVarBin);
      const Double_t valOaDensityP = hJetOaDensityP    -> GetBinContent(iVarBin);
      const Double_t valOaDensityG = hJetOaDensityG    -> GetBinContent(iVarBin);
      const Double_t valUeDensityP = hJetUeAvgDensityP -> GetBinContent(iVarBin);
      const Double_t valUeDensityG = hJetUeAvgDensityG -> GetBinContent(iVarBin);
      const Double_t errUeAvgP     = hJetUeAvgP        -> GetBinError(iVarBin);
      const Double_t errUeAvgG     = hJetUeAvgG        -> GetBinError(iVarBin);
      const Double_t errUeP        = hJetUeP           -> GetBinError(iVarBin);
      const Double_t errUeG        = hJetUeG           -> GetBinError(iVarBin);
      const Double_t errOaP        = hJetOaP           -> GetBinError(iVarBin);
      const Double_t errOaG        = hJetOaG           -> GetBinError(iVarBin);
      const Double_t errRhoP       = hJetRhoP          -> GetBinError(iVarBin);
      const Double_t errRhoG       = hJetRhoG          -> GetBinError(iVarBin);
      const Double_t errRhoAP      = hJetRhoAP         -> GetBinError(iVarBin);
      const Double_t errRhoAG      = hJetRhoAG         -> GetBinError(iVarBin);
      const Double_t errPtSubP     = hJetPtSubP        -> GetBinError(iVarBin);
      const Double_t errPtSubG     = hJetPtSubG        -> GetBinError(iVarBin);
      const Double_t errPtUeSubP   = hJetPtUeSubP      -> GetBinError(iVarBin);
      const Double_t errPtUeSubG   = hJetPtUeSubG      -> GetBinError(iVarBin);
      const Double_t errPtOaSubP   = hJetPtOaSubP      -> GetBinError(iVarBin);
      const Double_t errPtOaSubG   = hJetPtOaSubG      -> GetBinError(iVarBin);
      const Double_t errOaDensityP = hJetOaDensityP    -> GetBinContent(iVarBin);
      const Double_t errOaDensityG = hJetOaDensityG    -> GetBinContent(iVarBin);
      const Double_t errUeDensityP = hJetUeAvgDensityP -> GetBinContent(iVarBin);
      const Double_t errUeDensityG = hJetUeAvgDensityG -> GetBinContent(iVarBin);
      hJetUeAvgP        -> SetBinContent(iVarBin, valUeAvgP / binWidth);
      hJetUeAvgG        -> SetBinContent(iVarBin, valUeAvgG / binWidth);
      hJetUeP           -> SetBinContent(iVarBin, valUeP / binWidth);
      hJetUeG           -> SetBinContent(iVarBin, valUeG / binWidth);
      hJetOaP           -> SetBinContent(iVarBin, valOaP / binWidth);
      hJetOaG           -> SetBinContent(iVarBin, valOaG / binWidth);
      hJetRhoP          -> SetBinContent(iVarBin, valRhoP / binWidth);
      hJetRhoG          -> SetBinContent(iVarBin, valRhoG / binWidth);
      hJetRhoAP         -> SetBinContent(iVarBin, valRhoAP / binWidth);
      hJetRhoAG         -> SetBinContent(iVarBin, valRhoAG / binWidth);
      hJetPtSubP        -> SetBinContent(iVarBin, valPtSubP / binWidth);
      hJetPtSubG        -> SetBinContent(iVarBin, valPtSubG / binWidth);
      hJetPtUeSubP      -> SetBinContent(iVarBin, valPtUeSubP / binWidth);
      hJetPtUeSubG      -> SetBinContent(iVarBin, valPtUeSubG / binWidth);
      hJetPtOaSubP      -> SetBinContent(iVarBin, valPtOaSubP / binWidth);
      hJetPtOaSubG      -> SetBinContent(iVarBin, valPtOaSubG / binWidth);
      hJetOaDensityP    -> SetBinContent(iVarBin, valOaDensityP / binWidth);
      hJetOaDensityG    -> SetBinContent(iVarBin, valOaDensityG / binWidth);
      hJetUeAvgDensityP -> SetBinContent(iVarBin, valUeDensityP / binWidth);
      hJetUeAvgDensityG -> SetBinContent(iVarBin, valUeDensityG / binWidth);
      hJetUeAvgP        -> SetBinError(iVarBin, errUeAvgP / binWidth);
      hJetUeAvgG        -> SetBinError(iVarBin, errUeAvgG / binWidth);
      hJetUeP           -> SetBinError(iVarBin, errUeP / binWidth);
      hJetUeG           -> SetBinError(iVarBin, errUeG / binWidth);
      hJetOaP           -> SetBinError(iVarBin, errOaP / binWidth);
      hJetOaG           -> SetBinError(iVarBin, errOaG / binWidth);
      hJetRhoP          -> SetBinError(iVarBin, errRhoP / binWidth);
      hJetRhoG          -> SetBinError(iVarBin, errRhoG / binWidth);
      hJetRhoAP         -> SetBinError(iVarBin, errRhoAP / binWidth);
      hJetRhoAG         -> SetBinError(iVarBin, errRhoAG / binWidth);
      hJetPtSubP        -> SetBinError(iVarBin, errPtSubP / binWidth);
      hJetPtSubG        -> SetBinError(iVarBin, errPtSubG / binWidth);
      hJetPtUeSubP      -> SetBinError(iVarBin, errPtUeSubP / binWidth);
      hJetPtUeSubG      -> SetBinError(iVarBin, errPtUeSubG / binWidth);
      hJetPtOaSubP      -> SetBinError(iVarBin, errPtOaSubP / binWidth);
      hJetPtOaSubG      -> SetBinError(iVarBin, errPtOaSubG / binWidth);
      hJetOaDensityP    -> SetBinError(iVarBin, errOaDensityP / binWidth);
      hJetOaDensityG    -> SetBinError(iVarBin, errOaDensityG / binWidth);
      hJetUeAvgDensityP -> SetBinError(iVarBin, errUeDensityP / binWidth);
      hJetUeAvgDensityG -> SetBinError(iVarBin, errUeDensityG / binWidth);
    }
    hJetUeAvgP        -> Scale(1. / etaNormP);
    hJetUeAvgG        -> Scale(1. / etaNormG);
    hJetUeP           -> Scale(1. / etaNormP);
    hJetUeG           -> Scale(1. / etaNormG);
    hJetOaP           -> Scale(1. / etaNormP);
    hJetOaG           -> Scale(1. / etaNormG);
    hJetRhoP          -> Scale(1. / rhoNormP);
    hJetRhoG          -> Scale(1. / rhoNormG);
    hJetRhoAP         -> Scale(1. / etaNormP);
    hJetRhoAG         -> Scale(1. / etaNormG);
    hJetPtSubP        -> Scale(1. / etaNormP);
    hJetPtSubG        -> Scale(1. / etaNormG);
    hJetPtUeSubP      -> Scale(1. / etaNormP);
    hJetPtUeSubG      -> Scale(1. / etaNormG);
    hJetPtOaSubP      -> Scale(1. / etaNormP);
    hJetPtOaSubG      -> Scale(1. / etaNormG);
    hJetOaDensityP    -> Scale(1. / etaNormP);
    hJetOaDensityG    -> Scale(1. / etaNormG);
    hJetUeAvgDensityP -> Scale(1. / etaNormP);
    hJetUeAvgDensityG -> Scale(1. / etaNormG);
  } else {
    hJetUeAvgP        -> Scale(1. / corrNormP);
    hJetUeAvgG        -> Scale(1. / corrNormG);
    hJetUeP           -> Scale(1. / corrNormP);
    hJetUeG           -> Scale(1. / corrNormG);
    hJetOaP           -> Scale(1. / corrNormP);
    hJetOaG           -> Scale(1. / corrNormG);
    hJetRhoP          -> Scale(1. / avgNormP);
    hJetRhoG          -> Scale(1. / avgNormG);
    hJetRhoAP         -> Scale(1. / corrNormP);
    hJetRhoAG         -> Scale(1. / corrNormG);
    hJetPtSubP        -> Scale(1. / corrNormP);
    hJetPtSubG        -> Scale(1. / corrNormG);
    hJetPtUeSubP      -> Scale(1. / corrNormP);
    hJetPtUeSubG      -> Scale(1. / corrNormG);
    hJetPtOaSubP      -> Scale(1. / corrNormP);
    hJetPtOaSubG      -> Scale(1. / corrNormG);
    hJetOaDensityP    -> Scale(1. / corrNormP);
    hJetOaDensityG    -> Scale(1. / corrNormG);
    hJetUeAvgDensityP -> Scale(1. / corrNormP);
    hJetUeAvgDensityG -> Scale(1. / corrNormG);
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
  TDirectory *dAll = (TDirectory*) fOut -> mkdir("All");

  // write and close output
  fOut -> cd();
  for (Int_t i = 0; i < 3; i++) {
    if (i == 0) {
      dPi0 -> cd();
    } else if (i == 1) {
      dGam -> cd();
    } else if (i == 2) {
      dAll -> cd();
    }
    hTrgEt[i]           -> Write();
    hTrgTsp[i]          -> Write();
    hRefmult[i]         -> Write();
    hNumJets[i]         -> Write();
    hRho[i]             -> Write();
    hSigma[i]           -> Write();
    hJetArea[i]         -> Write();
    hJetEta[i]          -> Write();
    hAllDeltaPhi[i]     -> Write();
    hBinDeltaPhi[i][0]  -> Write();
    hBinDeltaPhi[i][1]  -> Write();
    hBinDeltaPhi[i][2]  -> Write();
    hJetDeltaPhi[i]     -> Write();
    hAllPtRaw[i]        -> Write();
    hJetPtRaw[i]        -> Write();
    hJetPtCorr[i]       -> Write();
    if (i == 0) {
      hJetUeP           -> Write();
      hJetOaP           -> Write();
      hJetRhoP          -> Write();
      hJetRhoAP         -> Write();
      hJetUeAvgP        -> Write();
      hJetPtSubP        -> Write();
      hJetPtUeSubP      -> Write();
      hJetPtOaSubP      -> Write();
      hJetOaDensityP    -> Write();
      hJetUeAvgDensityP -> Write();
      hJetPtVsRhoAP     -> Write();
    } else if (i == 1) {
      hDir              -> Write();
      hJetUeG           -> Write();
      hJetOaG           -> Write();
      hJetRhoG          -> Write();
      hJetRhoAG         -> Write();
      hJetUeAvgG        -> Write();
      hJetPtSubG        -> Write();
      hJetPtUeSubG      -> Write();
      hJetPtOaSubG      -> Write();
      hJetOaDensityG    -> Write();
      hJetUeAvgDensityG -> Write();
      hJetPtVsRhoAG     -> Write();
    }
    hJetPtVsDf[i]       -> Write();
    hJetPtVsA[i]        -> Write();
    fOut                -> cd();
  }
  fOut -> Close();

  // close input
  fIn -> cd();
  fIn -> Close();
  cout << "Reading script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
