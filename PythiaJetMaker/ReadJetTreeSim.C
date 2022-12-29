// 'ReadJetTreeSim.C'
// Nihar Sahoo, Derek Anderson
// 11.09.2016

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
static const Bool_t   DoTrgNorm   = false;
static const Double_t PiValue     = TMath::Pi();

// filepaths
static const TString iPath("output/March2022/pp200py8par.forQM22.et6100nTrg300Kpi0.r02rm1chrg.d26m3y2022.root");
static const TString wPath("input/Weights/pythiaWeights_usingQtHistToSmear.et911vz55tsp008pi0.d24m3y2021.root");
static const TString oPath("pp200py8.forEtBinCheck_etBinHalf_noBinCorr.et6100r02pi0.d19m8y2022.root");

// parameters
static const Double_t GamTsp = 22;
static const Double_t Pi0Tsp = 111;
static const Double_t MaxEta = 0.9;
static const Double_t MinPt  = 6.;
static const Double_t MaxPt  = 100.;

// jet parameters
static const Double_t Rjet     = 0.2;
static const Double_t MinArea  = 0.05;  // R02: 0.05, R03: 0.2, R04: 0.35, R05: 0.65, R07: 1.2
static const Double_t MinJetPt = 0.2;
static const Double_t MaxJetPt = 30.;
static const Double_t RecoilDf = PiValue / 4.;



void ReadJetTreeSim(Bool_t inBatchMode=false) {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\nBeginning reading script!" << endl;

  // open files
  TFile *fIn  = new TFile(iPath, "read");
  TFile *fOut = new TFile(oPath, "recreate");
  if (!fIn) {
    cerr << "PANIC: input file could not be opened!" << endl;
    assert(fIn);
  }

  // grab weights
  TFile *fWeight = new TFile(wPath.Data(), "read");
  if (!fWeight && DoWeighting) {
    cerr << "PANIC: weight file couldn't be opened!" << endl;
    assert(fWeight);
  }

  // get trees
  TTree *JetTree = 0;
  if (fIn) {
    fIn -> GetObject("femtoDST", JetTree);
  }
  if (!JetTree) {
    cerr << "PANIC: 'femtoDST' not grabbed!" << endl;
    assert(JetTree);
  }

  // get weights
  TH1D *hWeight;
  if (DoWeighting) {
    hWeight = (TH1D*) fWeight -> Get("hEtWeight");
    if (!hWeight) {
      cerr << "PANIC: couldn't grab weight histogram!" << endl;
      assert(hWeight);
    }
    cout << "  Doing weighting..." << endl;
  }

  // declare event leaves
  cout << "  Setting branch addresses..." << endl;
  Int_t    EventIndex = 0;
  Double_t Refmult    = 0.;
  Double_t TSP        = 0.;
  Double_t TrgEta     = 0.;
  Double_t TrgPhi     = 0.;
  Double_t TrgEt      = 0.;
  Double_t Rho        = 0.;
  Double_t Sigma      = 0.;
  Double_t Vz         = 0.;
  // declare jet leaves
  vector<Double_t> *JetEta    = 0;
  vector<Double_t> *JetPt     = 0;
  vector<Double_t> *JetNCons  = 0;
  vector<Double_t> *JetIndex  = 0;
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
  TBranch *bRefmult    = 0;
  TBranch *bTsp        = 0;
  TBranch *bTrgEta     = 0;
  TBranch *bTrgPhi     = 0;
  TBranch *bTrgEt      = 0;
  TBranch *bRho        = 0;
  TBranch *bSigma      = 0;
  TBranch *bVz         = 0;
  TBranch *bJetEta     = 0;
  TBranch *bJetPt      = 0;
  TBranch *bJetNCons   = 0;
  TBranch *bJetIndex   = 0;
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
  JetTree -> SetBranchAddress("Refmult", &Refmult, &bRefmult);
  JetTree -> SetBranchAddress("TSP", &TSP, &bTsp);
  JetTree -> SetBranchAddress("TrgEta", &TrgEta, &bTrgEta);
  JetTree -> SetBranchAddress("TrgPhi", &TrgPhi, &bTrgPhi);
  JetTree -> SetBranchAddress("TrgEt", &TrgEt, &bTrgEt);
  JetTree -> SetBranchAddress("Rho", &Rho, &bRho);
  JetTree -> SetBranchAddress("Sigma", &Sigma, &bSigma);
  JetTree -> SetBranchAddress("Vz", &Vz,&bVz);
  JetTree -> SetBranchAddress("JetIndex", &JetIndex, &bJetIndex);
  JetTree -> SetBranchAddress("JetPt", &JetPt, &bJetPt);
  JetTree -> SetBranchAddress("JetNCons", &JetNCons, &bJetNCons);
  JetTree -> SetBranchAddress("JetEta", &JetEta, &bJetEta);
  JetTree -> SetBranchAddress("JetPhi",&JetPhi, &bJetPhi); 
  JetTree -> SetBranchAddress("JetE", &JetE, &bJetE); 
  JetTree -> SetBranchAddress("JetArea",&JetArea, &bJetArea);
  JetTree -> SetBranchAddress("JetConsPt", &JetConsPt, &bJetConsPt);
  JetTree -> SetBranchAddress("JetConsEta", &JetConsEta, &bJetConsEta);
  JetTree -> SetBranchAddress("JetConsPhi", &JetConsPhi, &bJetConsPhi);
  JetTree -> SetBranchAddress("JetConsE", &JetConsE, &bJetConsE);

  // declare histograms [0 for pi0 trigger, 1 for gamma trigger]
  cout << "  Creating histograms..." << endl;
  TH1D *hTrgEt[2];
  TH1D *hTrgEtWeight[2];
  TH1D *hTrgEtRatio[2];
  TH1D *hRefmult[2];
  TH1D *hNumJets[2];
  TH1D *hRho[2];
  TH1D *hSigma[2];
  TH1D *hJetArea[2];
  TH1D *hJetEta[2];
  TH1D *hAllPtRaw[2];
  TH1D *hJetPtRaw[2];
  TH1D *hAllDeltaPhi[2];
  TH1D *hBinDeltaPhi[2][3];
  TH1D *hJetDeltaPhi[2];
  TH1D *hJetPtCorr[2];
  TH2D *hJetPtVsDf[2];
  TH2D *hJetPtVsA[2];

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
  const Int_t    eNum  = 200;
  //const Int_t    eNum  = 100;
  const Int_t    mNum  = 200;
  const Int_t    nNum  = 50;
  const Int_t    pNum  = 50;
  const Int_t    aNum  = 500;
  const Int_t    hNum  = 100;
  const Int_t    dFnum = 360;
  //const Int_t    pTnum = 65;
  const Int_t    pTnum = 100;
  const Double_t eBins[2]  = {0., 100.};
  const Double_t mBins[2]  = {0., 200.};
  const Double_t nBins[2]  = {0., 50.};
  const Double_t pBins[2]  = {0., 5.};
  const Double_t aBins[2]  = {0., 5.};
  const Double_t hBins[2]  = {-5., 5.};
  const Double_t dFbins[2] = {-1.* (PiValue / 2.), 3.* (PiValue / 2.)};
  //const Double_t pTbins[2] = {-8., 57.};
  const Double_t pTbins[2] = {-50., 50.};
  // pi0 histograms
  hTrgEt[0]          = new TH1D("hTrgEtP", "Trigger E_{T}, #pi^{0} trigger;E_{T}^{trg};counts", eNum, eBins[0], eBins[1]);
  hTrgEtWeight[0]    = new TH1D("hTrgEtWeightP", "Reweighted trigger E_{T}, #pi^{0};E_{T}^{trg};counts", eNum, eBins[0], eBins[1]);
  hTrgEtRatio[0]     = new TH1D("hTrgEtRatioP", "Ratio of raw to reweighted E_{T}^{trg}, #pi^{0};E_{T}^{trg};reweighted / raw", eNum, eBins[0], eBins[1]);
  hRefmult[0]        = new TH1D("hRefmultP", "Refmult, #pi^{0} trigger; refmult; counts", mNum, mBins[0], mBins[1]);
  hNumJets[0]        = new TH1D("hNumJetsP", "No. of jets, #pi^{0} trigger; N_{jet}; counts", nNum, nBins[0], nBins[1]);
  hSigma[0]          = new TH1D("hSigmaP", "Sigma, #pi^{0} trigger; #sigma; counts", pNum, pBins[0], pBins[1]);
  hRho[0]            = new TH1D("hRhoP", "Rho, #pi^{0} trigger; #rho; counts", pNum, pBins[0], pBins[1]);
  hJetArea[0]        = new TH1D("hJetAreaP", "Recoil jet area, #pi^{0} trigger; A_{jet}; (1/N_{trg})dN_{jet}/dA_{jet}", aNum, aBins[0], aBins[1]);
  hJetEta[0]         = new TH1D("hJetEtaP", "Recoil jet #eta, #pi^{0} trigger; #eta; (1/N_{trg})dN_{jet}/d#eta", hNum, hBins[0], hBins[1]);
hAllDeltaPhi[0]    = new TH1D("hAllDeltaPhiP", "All jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[0][0] = new TH1D("hBinDeltaPhiP_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[0][1] = new TH1D("hBinDeltaPhiP_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[0][2] = new TH1D("hBinDeltaPhiP_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, #pi^{0}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hJetDeltaPhi[0]    = new TH1D("hJetDeltaPhiP", "Recoil jet #Delta#varphi, #pi^{0} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  if (UseVarBins) {
    hAllPtRaw[0]     = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[0]     = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[0]    = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[0]     = new TH1D("hAllPtRawP", "All jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtRaw[0]     = new TH1D("hJetPtRawP", "Recoil jet p_{T}^{raw}, #pi^{0} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtCorr[0]    = new TH1D("hJetPtCorrP", "Recoil jet p_{T}^{corr}, #pi^{0} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}", pTnum, pTbins[0], pTbins[1]);
  }
  hJetPtVsDf[0]      = new TH2D("hJetPtVsDfP", "Jet p_{T}^{raw} vs. #Delta#varphi, #pi^{0} trigger; #Delta#varphi; p_{T}^{raw}", dFnum, dFbins[0], dFbins[1], pTnum, pTbins[0], pTbins[1]);
  hJetPtVsA[0]       = new TH2D("hJetPtVsAp", "Jet p_{T}^{raw} vs. Area, #pi^{0} trigger; A_{jet}; p_{T}^{raw}", aNum, aBins[0], aBins[1], pTnum, pTbins[0], pTbins[1]);
  // gamma histograms
  hTrgEt[1]          = new TH1D("hTrgEtG", "Trigger E_{T}, #gamma^{rich} trigger;E_{T}^{trg};counts", eNum, eBins[0], eBins[1]);
  hTrgEtWeight[1]    = new TH1D("hTrgEtWeightG", "Reweighted trigger E_{T}, #gamma^{rich};E_{T}^{trg};counts", eNum, eBins[0], eBins[1]);
  hTrgEtRatio[1]     = new TH1D("hTrgEtRatioG", "Ratio of raw to reweighted E_{T}^{trg}, #gamma^{rich};E_{T}^{trg};reweighted / raw", eNum, eBins[0], eBins[1]);
  hRefmult[1]        = new TH1D("hRefmultG", "Refmult, #gamma^{rich} trigger; refmult; counts", mNum, mBins[0], mBins[1]);
  hNumJets[1]        = new TH1D("hNumJetsG", "No. of jets, #gamma^{rich} trigger; N_{jet}; counts", nNum, nBins[0], nBins[1]);
  hRho[1]            = new TH1D("hRhoG", "Rho, #gamma^{rich} trigger; #rho; counts", pNum, pBins[0], pBins[1]);
  hSigma[1]          = new TH1D("hSigmaG", "Sigma, #gamma^{rich} trigger; #sigma; counts)", pNum, pBins[0], pBins[1]);
  hJetArea[1]        = new TH1D("hJetAreaG", "Recoil jet area, #gamma^{rich} trigger; A_{jet}; (1/N_{trg})dN_{jet}/dA_{jet}", aNum, aBins[0], aBins[1]);
  hJetEta[1]         = new TH1D("hJetEtaG", "Recoil jet #eta, #gamma^{rich} trigger; #eta; (1/N_{trg})dN_{jet}/d#eta", hNum, hBins[0], hBins[1]);
  hAllDeltaPhi[1]    = new TH1D("hAllDeltaPhiG", "All jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[1][0] = new TH1D("hBinDeltaPhiG_pT02", "#Delta#varphi of all jets with p_{T}^{jet}#in(0.2, 1) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[1][1] = new TH1D("hBinDeltaPhiG_pT1", "#Delta#varphi of all jets with p_{T}^{jet}#in(1, 5) GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hBinDeltaPhi[1][2] = new TH1D("hBinDeltaPhiG_pT5", "#Delta#varphi of all jets with p_{T}^{jet}>5 GeV/c, #gamma^{rich}; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  hJetDeltaPhi[1]    = new TH1D("hJetDeltaPhiG", "Recoil jet #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; (1/N_{trg})dN_{jet}/d#Delta#varphi", dFnum, dFbins[0], dFbins[1]);
  if (UseVarBins) {
    hAllPtRaw[1]     = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtRaw[1]     = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnumVar, pTbinVar);
    hJetPtCorr[1]    = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}", pTnumVar, pTbinVar);
  } else {
    hAllPtRaw[1]     = new TH1D("hAllPtRawG", "All jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtRaw[1]     = new TH1D("hJetPtRawG", "Recoil jet p_{T}^{raw}, #gamma^{rich} trigger; p_{T}^{raw}; (1/N_{trg})dN_{jet}/dp_{T}^{raw}", pTnum, pTbins[0], pTbins[1]);
    hJetPtCorr[1]    = new TH1D("hJetPtCorrG", "Recoil jet p_{T}^{corr}, #gamma^{rich} trigger; p_{T}^{corr}; (1/N_{trg})dN_{jet}/dp_{T}^{corr}", pTnum, pTbins[0], pTbins[1]);
  }
  hJetPtVsDf[1]      = new TH2D("hJetPtVsDfG", "Jet p_{T}^{raw} vs. #Delta#varphi, #gamma^{rich} trigger; #Delta#varphi; p_{T}^{raw}", dFnum, dFbins[0], dFbins[1], pTnum, pTbins[0], pTbins[1]);
  hJetPtVsA[1]       = new TH2D("hJetPtVsAg", "Jet p_{T}^{raw} vs. Area, #gamma^{rich} trigger; A_{jet}; p_{T}^{raw}", aNum, aBins[0], aBins[1], pTnum, pTbins[0], pTbins[1]);
  // errors
  for (Int_t i = 0; i < 2; i++) {
    hTrgEt[i]          -> Sumw2();
    hTrgEtWeight[i]    -> Sumw2();
    hTrgEtRatio[i]     -> Sumw2();
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
  Int_t   nBytes(0);
  Int_t   breakVal(0);
  Int_t   nTrgPi0(0);
  Int_t   nTrgGam(0);
  Float_t nWeightPi0(0.);
  Float_t nWeightGam(0.);
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
    // trigger info
    const Double_t hTrg  = TrgEta;
    const Double_t eTtrg = TrgEt;

     // trigger cuts
    const Bool_t isInTrgEtaCut = (TMath::Abs(hTrg) < MaxEta);
    const Bool_t isInTrgEtCut  = ((eTtrg > MinPt) && (eTtrg < MaxPt));
    if (!isInTrgEtaCut || !isInTrgEtCut) continue;

    // determine triggers [0 for pi0, 1 for gamma-rich]
    Int_t  tspFlag = -1;
    if (TSP == 111) {
      tspFlag = 0;
    } else if (TSP == 22) {
      tspFlag = 1;
    }

    // fill event histograms
    if (tspFlag == 0) {
      hTrgEt[0]   -> Fill(TrgEt);
      hRefmult[0] -> Fill(Refmult);
      hRho[0]     -> Fill(Rho);
      hSigma[0]   -> Fill(Sigma);
    } else if (tspFlag == 1) {
      hTrgEt[1]   -> Fill(TrgEt);
      hRefmult[1] -> Fill(Refmult);
      hRho[1]     -> Fill(Rho);
      hSigma[1]   -> Fill(Sigma);
    }

    // get weights
    Double_t eTweight(1.);
    if (DoWeighting) {
      const UInt_t  iEtWeight = hWeight -> FindBin(TrgEt);
      const Float_t eTbinVal  = hWeight -> GetBinContent(iEtWeight);
      if (eTbinVal > 0.) {
        eTweight = eTbinVal;
      } else {
        eTweight = 0.;
      }
    }

    // increment counters
    if (DoWeighting) {
      if (tspFlag == 0) {
        hTrgEtWeight[0] -> Fill(TrgEt, eTweight);
        nWeightPi0 += eTweight;
        ++nTrgPi0;
      }
      if (tspFlag == 1) {
        hTrgEtWeight[1] -> Fill(TrgEt, eTweight);
        nWeightGam += eTweight;
        ++nTrgGam;
      }
    } else {
      if (tspFlag == 0) {
        hTrgEtWeight[0] -> Fill(TrgEt);
        ++nWeightPi0;
        ++nTrgPi0;
      }
      if (tspFlag == 1) {
        hTrgEtWeight[1] -> Fill(TrgEt);
        ++nWeightGam;
        ++nTrgGam;
      }
    }

    // jet loop
    Int_t nJets = (Int_t) JetEta -> size();
    for (Int_t j = 0; j < nJets; j++) {

      Double_t h     = JetEta  -> at(j);
      Double_t a     = JetArea -> at(j);
      Double_t pT    = JetPt   -> at(j);
      Double_t f     = JetPhi  -> at(j);
      Double_t pTc   = pT - (Rho * a);
      Double_t dF    = f - TrgPhi;
      if (dF < -1.*PiValue/2.) dF += 2.*PiValue;
      if (dF > 3.*PiValue/2.)  dF -= 2.*PiValue;
      Double_t dFre  = f - TrgPhi;
      if (dFre < 0.)         dFre += 2.*PiValue;
      if (dFre > 2.*PiValue) dFre -= 2.*PiValue;
      Double_t dFcut = abs(dFre - PiValue);

      // fill histograms for all jets
      Bool_t isInBin1 = ((pT > 0.2) && (pT < 1.));
      Bool_t isInBin2 = ((pT > 1.) && (pT < 5.));
      Bool_t isInBin3 = (pT > 5.);
      if (tspFlag == 0) {
        hAllPtRaw[0]    -> Fill(pT);
        hAllDeltaPhi[0] -> Fill(dF);
        hJetPtVsDf[0]   -> Fill(dF, pT);
        hJetPtVsA[0]    -> Fill(a, pT);
        if (isInBin1) {
          hBinDeltaPhi[0][0] -> Fill(dF);
        }
        if (isInBin2) {
          hBinDeltaPhi[0][1] -> Fill(dF);
        }
        if (isInBin3) {
          hBinDeltaPhi[0][2] -> Fill(dF);
        }
      }
      if (tspFlag == 1) {
        hAllPtRaw[1]    -> Fill(pT);
        hAllDeltaPhi[1] -> Fill(dF);
        hJetPtVsDf[1]   -> Fill(dF, pT);
        hJetPtVsA[1]    -> Fill(a, pT);
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

      // jet cuts
      const Bool_t isInPtJetCut = ((pT > MinJetPt) && (pT < MaxJetPt));
      const Bool_t isInAjetCut  = (a > MinArea);
      const Bool_t isRecoil     = (dFcut < RecoilDf);
      if (!isInPtJetCut || !isInAjetCut || !isRecoil) continue;

      // fill histograms
      if (DoWeighting) {
        if (tspFlag == 0) {
          hJetArea[0]     -> Fill(a, eTweight);
          hJetEta[0]      -> Fill(h, eTweight);
          hJetDeltaPhi[0] -> Fill(dF, eTweight);
          hJetPtRaw[0]    -> Fill(pT, eTweight);
          hJetPtCorr[0]   -> Fill(pTc, eTweight);
        } else if (tspFlag == 1) {
          hJetArea[1]     -> Fill(a, eTweight);
          hJetEta[1]      -> Fill(h, eTweight);
          hJetDeltaPhi[1] -> Fill(dF, eTweight);
          hJetPtRaw[1]    -> Fill(pT, eTweight);
          hJetPtCorr[1]   -> Fill(pTc, eTweight);
        }
      }
      else {
        if (tspFlag == 0) {
          hJetArea[0]     -> Fill(a);
          hJetEta[0]      -> Fill(h);
          hJetDeltaPhi[0] -> Fill(dF);
          hJetPtRaw[0]    -> Fill(pT);
          hJetPtCorr[0]   -> Fill(pTc);
        } else if (tspFlag == 1) {
          hJetArea[1]     -> Fill(a);
          hJetEta[1]      -> Fill(h);
          hJetDeltaPhi[1] -> Fill(dF);
          hJetPtRaw[1]    -> Fill(pT);
          hJetPtCorr[1]   -> Fill(pTc);
        }
      }
    }  // end jet loop

    if (tspFlag == 0) {
      hNumJets[0] -> Fill(nJets);
    } else if (tspFlag == 1) {
      hNumJets[1] -> Fill(nJets);
    }
  }  // end event loop

  if (breakVal == 1) {
    cerr << "ERROR: Occured during event loop!\n"
         << "       Aborting program!"
         << endl;
    cassert(0);
  } else {
    cout << "  Event loop finished!\n"
         << "    nTrgPi0 = " << nTrgPi0 << ", nTrgPi0(weighted) = " << nWeightPi0 << "\n"
         << "    nTrgGam = " << nTrgGam << ", nTrgGam(weighted) = " << nWeightGam
         << endl;
  }

  // normalize histograms
  cout << "  Normalizing histograms..." << endl;
  const Double_t eBin    = (eBins[1] - eBins[0]) / eNum;
  const Double_t mBin    = (mBins[1] - mBins[0]) / mNum;
  const Double_t nBin    = (nBins[1] - nBins[0]) / nNum;
  const Double_t pBin    = (pBins[1] - pBins[0]) / pNum;
  const Double_t aBin    = (aBins[1] - aBins[0]) / aNum;
  const Double_t hBin    = (hBins[1] - hBins[0]) / hNum;
  const Double_t dFbin   = (dFbins[1] - dFbins[0]) / dFnum;
  const Double_t pTbin   = (pTbins[1] - pTbins[0]) / pTnum;
  const Double_t hJetBin = 2. * (1. - Rjet);
  for (Int_t i = 0; i < 2; i++) {
    Double_t norm  = 1.;
    Double_t normW = 1.;
    if (DoTrgNorm) {
      if (DoWeighting) {
        if (i == 0) {
          norm  = (Double_t) nTrgPi0;
          normW = nWeightPi0;
        } else if (i == 1) {
          norm  = (Double_t) nTrgGam;
          normW = nWeightGam;
        }
      } else {
        if (i == 0)  {
          norm  = (Double_t) nTrgPi0;
          normW = (Double_t) nTrgPi0;
        } else if (i == 1) {
          norm  = (Double_t) nTrgGam;
          normW = (Double_t) nTrgGam;
        }
      }
    }
    const Double_t eNorm   = eBin * norm;
    const Double_t mNorm   = mBin * norm;
    const Double_t nNorm   = nBin * norm;
    const Double_t pNorm   = pBin * norm;
    const Double_t aNorm   = aBin * norm * hJetBin;
    const Double_t aNormW  = aBin * normW * hJetBin;
    const Double_t hNorm   = hBin * norm * hJetBin;
    const Double_t hNormW  = hBin * normW * hJetBin;
    const Double_t dFnorm  = dFbin * norm * hJetBin;
    const Double_t dFnormW = dFbin * normW * hJetBin;
    const Double_t pTnorm  = normW * hJetBin;
    const Double_t pTnormW = normW * hJetBin;
    //hTrgEt[i]       -> Scale(1. / eNorm);
    //hTrgEtWeight[i] -> Scale(1. / eNorm);
    hRefmult[i]     -> Scale(1. / mNorm);
    hNumJets[i]     -> Scale(1. / nNorm);
    hRho[i]         -> Scale(1. / pNorm);
    hSigma[i]       -> Scale(1. / pNorm);
    if (DoWeighting ) {
      hJetArea[i]        -> Scale(1. / aNormW);
      hJetEta[i]         -> Scale(1. / hNormW);
      hAllDeltaPhi[i]    -> Scale(1. / dFnormW);
      hBinDeltaPhi[i][0] -> Scale(1. / dFnormW);
      hBinDeltaPhi[i][1] -> Scale(1. / dFnormW);
      hBinDeltaPhi[i][2] -> Scale(1. / dFnormW);
      hJetDeltaPhi[i]    -> Scale(1. / dFnormW);
      hAllPtRaw[i]       -> Scale(1. / pTnormW);
      hJetPtRaw[i]       -> Scale(1. / pTnormW);
      hJetPtCorr[i]      -> Scale(1. / pTnormW);
      hJetPtVsDf[i]      -> Scale(1. / (normW * dFbin * pTbin * hJetBin));
      hJetPtVsA[i]       -> Scale(1. / (normW * aBin * pTbin * hJetBin));
    } else {
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
    }

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

  // take ratio of eTtrg
  for (UInt_t iTrg = 0 ; iTrg < 2; iTrg++) {
    hTrgEtRatio[iTrg] -> Divide(hTrgEtWeight[iTrg], hTrgEt[iTrg], 1., 1.);
  }

  // create directories
  TDirectory *dPi0 = (TDirectory*) fOut -> mkdir("Pi0");
  TDirectory *dGam = (TDirectory*) fOut -> mkdir("Gam");

  // write and close output
  fOut -> cd();
  for (Int_t i = 0; i < 2; i++) {
    if (i == 0) {
      dPi0 -> cd();
    } else if (i == 1) {
      dGam -> cd();
    }
    hTrgEt[i]          -> Write();
    hTrgEtWeight[i]    -> Write();
    hTrgEtRatio[i]     -> Write();
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
