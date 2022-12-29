// 'ReweightJetSpectra.C'
// Derek Anderson
// 02.03.2021
//
// Takes the JetTree output from the
// 'StJetTreeMaker' class, calculates
// a jet spectrum, and weights it
// according to weights computed by
// the 'CorrectTriggerSpectrum.C'
// code.

#include <vector>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TDirectory.h"

using namespace std;

// global constants
static const UInt_t NTrg(2);
static const UInt_t NVarBins(27);



void ReweightJetSpectrum() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning jet reweighting..." << endl;

  // io parameters
  const TString sOut("triggerReweighting_noReusingEvts_eTsmearNormedTo911.et911pt0230dca1vz55.d25m2y2021.root");
  const TString sJetIn("");
  const TString sWeightIn("input/triggerWeights_et911gamWeights_full650gamSample.et650x911vz55tsp0206gam.d10m3y2021.root");
  const TString sJetTreeIn("JetTree");
  const TString sWeightHistIn("hEtWeights");
  const TString sQtFuncIn("fQtFit");

  // calculation parameters
  const Bool_t   reuseEvents(false);
  const Double_t eTtrgCorr[2] = {6., 50.};

  // trigger parameters
  const Double_t zVtxMax(55.);
  const Double_t hTrgMax(0.9);
  const Double_t eTtrgRaw[2] = {9., 11.};
  const Double_t tspPi0[2]   = {0., 0.08};
  const Double_t tspGam[2]   = {0.2, 0.6};

  // jet parameters
  const Bool_t   doVariableBins(true);
  const Double_t rJet(0.2);
  const Double_t hJetMax(1.);
  const Double_t aJetMin(0.05);  // R02: 0.05, R03: 0.2, R04: 0.35, R05: 0.65, R07: 1.2
  const Double_t dFrecoil(TMath::PiOver4());
  const Double_t pTjetCut[2] = {0.2, 30.};


  // open files
  TFile *fOut      = new TFile(sOut.Data(), "recreate");
  TFile *fJetIn    = new TFile(sJetIn.Data(), "read");
  TFile *fWeightIn = new TFile(sWeightIn.Data(), "read");
  if (!fOut || !fJetIn || !fWeightIn) {
    cerr <<"PANIC: couldn't open file!\n"
         << "      fOut      = " << fOut << "\n"
         << "      fJetIn    = " << fJetIn << "\n"
         << "      fWeightIn = " << fWeightIn << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab jet tree
  TTree *tJets = (TTree*) fJetIn -> Get(sJetTreeIn.Data());
  if (!tJets) {
    cerr << "PANIC: couldn't grab jet tree!\n" << endl;
    return;
  }
  cout << "    Grabbed jet tree." << endl;

  // grab weight input
  TH1D *hWeight = (TH1D*) fWeightIn -> Get(sWeightHistIn.Data());
  TF1  *fQtTrg  = (TF1*)  fWeightIn -> Get(sQtFuncIn.Data());
  if (!hWeight || !fQtTrg) {
    cerr << "PANIC: couldn't grab weight histogram or qT function!\n"
         << "       hWeight = " << hWeight << "\n"
         << "       fQtTrg  = " << fQtTrg
         << endl;
    return;
  }
  hWeight -> SetName("hEtWeightInput");
  fQtTrg  -> SetName("fQtTrgInput");
  cout << "    Grabbed weight histogram and qT function." << endl;


  // variable jet bins
  const UInt_t   pTnumVar(NVarBins - 1);
  const Double_t pTbinVar[NVarBins] = {-12., -8., -5., -3., -2., -1.5, -1., -0.6, -0.2, 0., 0.2, 0.6, 1., 1.5, 2., 3., 5., 8., 12., 17., 23., 30., 38., 47., 57., 68., 81.};

  // declare histograms
  TH1D *hNumTrgRaw[NTrg];
  TH1D *hNumTrgCorr[NTrg];
  TH1D *hNumTrgIn[NTrg];
  TH1D *hNumTrgOut[NTrg];
  TH1D *hEtTrgRaw[NTrg];
  TH1D *hEtTrgCorr[NTrg];
  TH1D *hEtTrgCalc[NTrg];
  TH1D *hPtJetRaw[NTrg];
  TH1D *hPtJetRawN[NTrg];
  TH1D *hPtJetWeight[NTrg];
  TH1D *hPtJetWeightN[NTrg];
  TH1D *hPtJetRatio[NTrg];
  TH1D *hPtJetRatioN[NTrg];

  // binning
  const UInt_t nNum(3);
  const UInt_t eTnum(200);
  const UInt_t pTnum(65);
  const Float_t nBins[2]  = {0., 3.};
  const Float_t eTbins[2] = {0., 100.};
  const Float_t pTbins[2] = {-8., 57.};
  // nTrg histograms
  hNumTrgRaw[0]    = new TH1D("hNumTrgRawP", "N_{#pi^{0}} uncorrected", nNum, nBins[0], nBins[1]);
  hNumTrgRaw[1]    = new TH1D("hNumTrgRawG", "N_{#gamma^{rich}} uncorrected", nNum, nBins[0], nBins[1]);
  hNumTrgCorr[0]   = new TH1D("hNumTrgCorrP", "N_{#pi^{0}} corrected", nNum, nBins[0], nBins[1]);
  hNumTrgCorr[1]   = new TH1D("hNumTrgCorrG", "N_{#gamma^{rich}} corrected", nNum, nBins[0], nBins[1]);
  hNumTrgIn[0]     = new TH1D("hNumTrgInP", "N_{#pi^{0}} in range after back-smear", nNum, nBins[0], nBins[1]);
  hNumTrgIn[1]     = new TH1D("hNumTrgInG", "N_{#gamma^{rich}} in range after back-smear", nNum, nBins[0], nBins[1]);
  hNumTrgOut[0]    = new TH1D("hNumTrgOutP", "N_{#pi^{0}} out of range after back-smear", nNum, nBins[0], nBins[1]);
  hNumTrgOut[1]    = new TH1D("hNumTrgOutG", "N_{#gamma^{rich}} out of range after back-smear", nNum, nBins[0], nBins[1]);
  // eTtrg histograms
  hEtTrgRaw[0]     = new TH1D("hEtTrgRawP", "E_{T}^{trg}(#pi^{0}) uncorrected", eTnum, eTbins[0], eTbins[1]);
  hEtTrgRaw[1]     = new TH1D("hEtTrgRawG", "E_{T}^{trg}(#gamma^{rich}) uncorrected", eTnum, eTbins[0], eTbins[1]);
  hEtTrgCorr[0]    = new TH1D("hEtTrgCorrP", "E_{T}^{trg}(#pi^{0}) corrected", eTnum, eTbins[0], eTbins[1]);
  hEtTrgCorr[1]    = new TH1D("hEtTrgCorrG", "E_{T}^{trg}(#gamma^{rich}) corrected", eTnum, eTbins[0], eTbins[1]);
  hEtTrgCalc[0]    = new TH1D("hEtTrgCorrReuseP", "E_{T}^{trg}(#pi^{0}) corr. (no weight, during qT sampling, regardless of range)", eTnum, eTbins[0], eTbins[1]);
  hEtTrgCalc[1]    = new TH1D("hEtTrgCorrReuseG", "E_{T}^{trg}(#gamma^{rich}) corr. (no weight, during qT smapling, regardless of range)", eTnum, eTbins[0], eTbins[1]);
  // pTjet histograms
  if (doVariableBins) {
    hPtJetRaw[0]     = new TH1D("hPtJetRawP", "Recoil jet p_{T}^{reco}, #pi^{0} trigger", pTnumVar, pTbinVar);
    hPtJetRaw[1]     = new TH1D("hPtJetRawG", "Recoil jet p_{T}^{reco}, #gamma^{rich} trigger", pTnumVar, pTbinVar);
    hPtJetRawN[0]    = new TH1D("hPtJetRawNormP", "Recoil jet p_{T}^{reco} normed by N_{trg}^{raw}, #pi^{0} trigger", pTnumVar, pTbinVar);
    hPtJetRawN[1]    = new TH1D("hPtJetRawNormG", "Recoil jet p_{T}^{reco} normed by N_{trg}^{raw}, #gamma^{rich} trigger", pTnumVar, pTbinVar);
    hPtJetWeight[0]  = new TH1D("hPtJetWeightP", "Recoil jet p_{T}^{reco} (reweighted), #pi^{0} trigger", pTnumVar, pTbinVar);
    hPtJetWeight[1]  = new TH1D("hPtJetWeightG", "Recoil jet p_{T}^{reco} (reweighted), #gamma^{rich} trigger", pTnumVar, pTbinVar);
    hPtJetWeightN[0] = new TH1D("hPtJetWeightNormP", "Recoil jet p_{T}^{reco} (reweighted) normed by N_{trg}^{corr}, #pi^{0} trigger", pTnumVar, pTbinVar);
    hPtJetWeightN[1] = new TH1D("hPtJetWeightNormG", "Recoil jet p_{T}^{reco} (reweighted) normed by N_{trg}^{corr}, #gamma^{rich} trigger", pTnumVar, pTbinVar);
    hPtJetRatio[0]   = new TH1D("hPtJetRatioP", "Ratio between raw and weighted p_{T}^{reco}, #pi^{0} trigger", pTnumVar, pTbinVar);
    hPtJetRatio[1]   = new TH1D("hPtJetRatioG", "Ratio between raw and weighted p_{T}^{reco}, #gamma^{rich} trigger", pTnumVar, pTbinVar);
    hPtJetRatioN[0]  = new TH1D("hPtJetRatioNormP", "Ratio between normalized raw and weighted p_{T}^{reco}, #pi^{0} trigger", pTnumVar, pTbinVar);
    hPtJetRatioN[1]  = new TH1D("hPtJetRatioNormG", "Ratio between normalized raw and weighted p_{T}^{reco}, #gamma^{rich} trigger", pTnumVar, pTbinVar);
  }
  else {
    hPtJetRaw[0]     = new TH1D("hPtJetRawP", "Recoil jet p_{T}^{reco}, #pi^{0} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRaw[1]     = new TH1D("hPtJetRawG", "Recoil jet p_{T}^{reco}, #gamma^{rich} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRawN[0]    = new TH1D("hPtJetRawNormP", "Recoil jet p_{T}^{reco} normed by N_{trg}^{raw}, #pi^{0} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRawN[1]    = new TH1D("hPtJetRawNormG", "Recoil jet p_{T}^{reco} normed by N_{trg}^{raw}, #gamma^{rich} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetWeight[0]  = new TH1D("hPtJetWeightP", "Recoil jet p_{T}^{reco} (reweighted), #pi^{0} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetWeight[1]  = new TH1D("hPtJetWeightG", "Recoil jet p_{T}^{reco} (reweighted), #gamma^{rich} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetWeightN[0] = new TH1D("hPtJetWeightNormP", "Recoil jet p_{T}^{reco} (reweighted) normed by N_{trg}^{corr}, #pi^{0} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetWeightN[1] = new TH1D("hPtJetWeightNormG", "Recoil jet p_{T}^{reco} (reweighted) normed by N_{trg}^{corr}, #gamma^{rich} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRatio[0]   = new TH1D("hPtJetRatioP", "Ratio between raw and weighted p_{T}^{reco}, #pi^{0} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRatio[1]   = new TH1D("hPtJetRatioG", "Ratio between raw and weighted p_{T}^{reco}, #gamma^{rich} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRatioN[0]  = new TH1D("hPtJetRatioNormP", "Ratio between normalized raw and weighted p_{T}^{reco}, #pi^{0} trigger", pTnum, pTbins[0], pTbins[1]);
    hPtJetRatioN[1]  = new TH1D("hPtJetRatioNormG", "Ratio between normalized raw and weighted p_{T}^{reco}, #gamma^{rich} trigger", pTnum, pTbins[0], pTbins[1]);
  }
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hEtTrgRaw[iTrg]     -> Sumw2();
    hEtTrgCorr[iTrg]    -> Sumw2();
    hEtTrgCalc[iTrg]    -> Sumw2();
    hPtJetRaw[iTrg]     -> Sumw2();
    hPtJetRawN[iTrg]    -> Sumw2();
    hPtJetWeight[iTrg]  -> Sumw2();
    hPtJetWeightN[iTrg] -> Sumw2();
    hPtJetRatio[iTrg]   -> Sumw2();
    hPtJetRatioN[iTrg]  -> Sumw2();
  }
  cout << "    Declared histograms." << endl;


  // declare event leaves
  Int_t    EventIndex    = 0;
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
  tJets -> SetBranchAddress("eventIndex", &EventIndex, &bEventIndex);
  tJets -> SetBranchAddress("Refmult", &Refmult, &bRefmult);
  tJets -> SetBranchAddress("NJets", &NJets, &bNJets);
  tJets -> SetBranchAddress("TSP", &TSP, &bTsp);
  tJets -> SetBranchAddress("TrgEta", &TrgEta, &bTrgEta);
  tJets -> SetBranchAddress("TrgPhi", &TrgPhi, &bTrgPhi);
  tJets -> SetBranchAddress("TrgEt", &TrgEt, &bTrgEt);
  tJets -> SetBranchAddress("Rho", &Rho, &bRho);
  tJets -> SetBranchAddress("Sigma", &Sigma, &bSigma);
  tJets -> SetBranchAddress("Vz", &Vz,&bVz);
  tJets -> SetBranchAddress("JetIndex", &JetIndex, &bJetIndex);
  tJets -> SetBranchAddress("JetNCons", &JetNCons, &bJetNCons);
  tJets -> SetBranchAddress("JetPt", &JetPt, &bJetPt);
  tJets -> SetBranchAddress("JetPtCorr", &JetPtCorr, &bJetPtCorr);
  tJets -> SetBranchAddress("JetEta", &JetEta, &bJetEta);
  tJets -> SetBranchAddress("JetPhi",&JetPhi, &bJetPhi); 
  tJets -> SetBranchAddress("JetE", &JetE, &bJetE); 
  tJets -> SetBranchAddress("JetArea",&JetArea, &bJetArea);
  tJets -> SetBranchAddress("JetConsPt", &JetConsPt, &bJetConsPt);
  tJets -> SetBranchAddress("JetConsEta", &JetConsEta, &bJetConsEta);
  tJets -> SetBranchAddress("JetConsPhi", &JetConsPhi, &bJetConsPhi);
  tJets -> SetBranchAddress("JetConsE", &JetConsE, &bJetConsE);
  cout << "    Set jet tree branches." << endl;


  // set random seed
  TDatime *dateTime = new TDatime();
  Int_t    date     = dateTime -> GetDate();
  Int_t    time     = dateTime -> GetTime();
  Int_t    seed     = date * time;
  gRandom -> SetSeed(seed);
  cout << "    Set random seed:\n"
       << "      seed = " << seed
       << endl;

  // qT parameters
  const Double_t ampQt = fQtTrg -> GetParameter(0);
  const Double_t muQt  = fQtTrg -> GetParameter(1);
  const Double_t sigQt = fQtTrg -> GetParameter(2);
  cout << "    Grabbed qT parameters:\n"
       << "      ampl. = " << ampQt << "\n"
       << "      mu    = " << muQt << "\n"
       << "      sigma = " << sigQt
       << endl;


  const UInt_t nEvts = (UInt_t) tJets -> GetEntries();
  cout << "    Beginning event loop:" << endl;


  // event loop
  UInt_t  bytes(0);
  UInt_t  nBytes(0);
  UInt_t  nTrgRaw[NTrg]   = {0, 0};
  UInt_t  nInRange[NTrg]  = {0, 0};
  UInt_t  nOutRange[NTrg] = {0., 0.};
  Float_t nTrgCorr[NTrg]  = {0., 0.};
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    bytes   = tJets -> GetEntry(iEvt);
    nBytes += bytes;
    if (bytes < 0) {
      cerr << "WARNING: issue with event " << iEvt << "!" << endl;
      break;
    }
    else {
      cout << "      Processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
      if (iEvt + 1 == nEvts) cout << endl;
    }


    // event info
    const Double_t zVtx  = Vz;
    const Double_t hTrg  = TrgEta;
    const Double_t fTrg  = TrgPhi;
    const Double_t eTtrg = TrgEt;
    const Double_t tsp   = TSP;

    // event cuts
    const Bool_t isInVzTrgCut  = (TMath::Abs(zVtx) < zVtxMax);
    const Bool_t isInEtaTrgCut = (TMath::Abs(hTrg) < hTrgMax);
    const Bool_t isInEtTrgCut  = ((eTtrg > eTtrgRaw[0]) && (eTtrg < eTtrgRaw[1]));
    const Bool_t isInPi0cut    = ((tsp > tspPi0[0]) && (tsp < tspPi0[1]));
    const Bool_t isInGamCut    = ((tsp > tspGam[0]) && (tsp < tspGam[1]));
    const Bool_t isInTspCut    = (isInPi0cut || isInGamCut);
    if (!isInVzTrgCut || !isInEtaTrgCut || !isInEtTrgCut || !isInTspCut) continue;

    // fill raw trig. histograms
    if (isInPi0cut) {
      hEtTrgRaw[0] -> Fill(eTtrg);
      nTrgRaw[0]++;
    }
    if (isInGamCut) {
      hEtTrgRaw[1] -> Fill(eTtrg);
      nTrgRaw[1]++;
    }


    // back-smear eTtrg
    Double_t qTtrg(0.);
    Double_t eTcorr(0.);
    Bool_t   isInRange(false);
    if (reuseEvents) {
      do {
        qTtrg     = gRandom -> Gaus(muQt, sigQt);
        eTcorr    = eTtrg / qTtrg;
        isInRange = ((eTcorr > eTtrgCorr[0]) && (eTcorr < eTtrgCorr[1]));
        if (isInPi0cut) hEtTrgCalc[0] -> Fill(eTcorr);
        if (isInGamCut) hEtTrgCalc[1] -> Fill(eTcorr);
      } while (!isInRange);
    }
    else {
      qTtrg     = gRandom -> Gaus(muQt, sigQt);
      eTcorr    = eTtrg / qTtrg;
      isInRange = ((eTcorr > eTtrgCorr[0]) && (eTcorr < eTtrgCorr[1]));
      if (isInPi0cut) hEtTrgCalc[0] -> Fill(eTcorr);
      if (isInGamCut) hEtTrgCalc[1] -> Fill(eTcorr);
    }

    // grab weights and fill corrected histograms
    const UInt_t  iEtWeight = hWeight -> FindBin(eTcorr);
    const Float_t eTweight  = hWeight -> GetBinContent(iEtWeight);
    if (isInRange) {
      if (isInPi0cut) {
        hEtTrgCorr[0] -> Fill(eTcorr, eTweight);
        nTrgCorr[0]   += eTweight;
        nInRange[0]++;
      }
      if (isInGamCut) {
        hEtTrgCorr[1] -> Fill(eTcorr, eTweight);
        nTrgCorr[1]   += eTweight;
        nInRange[1]++;
      }
    }
    else {
      if (isInPi0cut) nOutRange[0]++;
      if (isInGamCut) nOutRange[1]++;
      continue;
    }


    // jet loop
    UInt_t nJets = (UInt_t) NJets;
    for (Int_t iJet = 0; iJet < nJets; iJet++) {

      // jet info
      const Double_t hJet   = JetEta    -> at(iJet);
      const Double_t pTjet  = JetPt     -> at(iJet);
      const Double_t pTreco = JetPtCorr -> at(iJet);
      const Double_t fJet   = JetPhi    -> at(iJet);
      const Double_t aJet   = JetArea   -> at(iJet);

      // calculate delta-phi
      Double_t dFjet = fJet - fTrg;
      if (dFjet < (-1. * TMath::PiOver2())) dFjet += TMath::TwoPi();
      if (dFjet > (3. * TMath::PiOver2()))  dFjet -= TMath::TwoPi();
      const Double_t dFcut = TMath::Abs(dFjet - TMath::Pi());


      // jet cuts
      const Bool_t isInPtJetCut = ((pTjet > pTjetCut[0]) && (pTjet < pTjetCut[1]));
      const Bool_t isInAjetCut  = (aJet > aJetMin);
      const Bool_t isRecoil     = (dFcut < dFrecoil);
      if (!isInPtJetCut || !isInAjetCut || !isRecoil) continue;

      // fill pi0/gamma histograms
      if (isInPi0cut) {
        hPtJetRaw[0]     -> Fill(pTreco);
        hPtJetRawN[0]    -> Fill(pTreco);
        hPtJetWeight[0]  -> Fill(pTreco, eTweight);
        hPtJetWeightN[0] -> Fill(pTreco, eTweight);
      }
      if (isInGamCut) {
        hPtJetRaw[1]     -> Fill(pTreco);
        hPtJetRawN[1]    -> Fill(pTreco);
        hPtJetWeight[1]  -> Fill(pTreco, eTweight);
        hPtJetWeightN[1] -> Fill(pTreco, eTweight);
      }
    }  // end jet loop
  }  // end event loop
  cout << "    Event loop finished:\n"
       << "      nTrgRaw(pi0, gam)   = (" << nTrgRaw[0] << ", " << nTrgRaw[1] << ")\n"
       << "      nTrgCorr(pi0, gam)  = (" << nTrgCorr[0] << ", " << nTrgCorr[1] << ")\n"
       << "      nInRange(pi0, gam)  = (" << nInRange[0] << ", " << nInRange[1] << ")\n"
       << "      nOutRange(pi0, gam) = (" << nOutRange[0] << ", " << nOutRange[1] << ")"
       << endl;

  // record trigger no.s
  for (iTrg = 0; iTrg < NTrg; iTrg++) {
    const UInt_t iBin = hNumTrgRaw[iTrg] -> FindBin(1.5);
    hNumTrgRaw[iTrg]  -> SetBinContent(iBin, nTrgRaw[iTrg]);
    hNumTrgCorr[iTrg] -> SetBinContent(iBin, nTrgCorr[iTrg]);
    hNumTrgIn[iTrg]   -> SetBinContent(iBin, nInRange[iTrg]);
    hNumTrgOut[iTrg]  -> SetBinContent(iBin, nOutRange[iTrg]);
  }
  cout << "    Recorded trigger numbers." << endl;


  // normalize jet histograms
  const UInt_t  nJetBins = hPtJetRaw[0] -> GetNbinsX();
  const Float_t hJetBin  = 2. * TMath::Abs(hJetMax - rJet);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    if (doVariableBins) {
      for (UInt_t iBin = 1; iBin < (nJetBins + 1); iBin++) {
        const Float_t binWidth = hPtJetRaw[iTrg]     -> GetBinWidth(iBin);
        const Float_t binValR  = hPtJetRaw[iTrg]     -> GetBinContent(iBin);
        const Float_t binValRN = hPtJetRawN[iTrg]    -> GetBinContent(iBin);
        const Float_t binValW  = hPtJetWeight[iTrg]  -> GetBinContent(iBin);
        const Float_t binValWN = hPtJetWeightN[iTrg] -> GetBinContent(iBin);
        const Float_t binErrR  = hPtJetRaw[iTrg]     -> GetBinError(iBin);
        const Float_t binErrRN = hPtJetRawN[iTrg]    -> GetBinError(iBin);
        const Float_t binErrW  = hPtJetWeight[iTrg]  -> GetBinError(iBin);
        const Float_t binErrWN = hPtJetWeightN[iTrg] -> GetBinError(iBin);
        hPtJetRaw[iTrg]     -> SetBinContent(iBin, binValR / binWidth);
        hPtJetRawN[iTrg]    -> SetBinContent(iBin, binValRN / binWidth);
        hPtJetWeight[iTrg]  -> SetBinContent(iBin, binValW / binWidth);
        hPtJetWeightN[iTrg] -> SetBinContent(iBin, binValWN / binWidth);
        hPtJetRaw[iTrg]     -> SetBinError(iBin, binErrR / binWidth);
        hPtJetRawN[iTrg]    -> SetBinError(iBin, binErrRN / binWidth);
        hPtJetWeight[iTrg]  -> SetBinError(iBin, binErrW / binWidth);
        hPtJetWeightN[iTrg] -> SetBinError(iBin, binErrWN / binWidth);
      }  // end bin loop
      const Float_t rawNorm    = nTrgRaw[iTrg] * hJetBin;
      const Float_t weightNorm = nTrgCorr[iTrg] * hJetBin;
      hPtJetRawN[iTrg]    -> Scale(1. / rawNorm);
      hPtJetWeightN[iTrg] -> Scale(1. / weightNorm);
    }
    else {
      const Float_t binWidth   = hPtJetRaw[iTrg] -> GetBinWidth(2);
      const Float_t rawNorm    = nTrgRaw[iTrg] * hJetBin * binWidth;
      const Float_t weightNorm = nTrgCorr[iTrg] * hJetBin * binWidth;
      hPtJetRawN[iTrg]    -> Scale(1. / rawNorm);
      hPtJetWeightN[iTrg] -> Scale(1. / weightNorm);
    }
  }  // end trigger loop
  cout << "    Normalized jet distributions." << endl;

  // take ratios
  const Float_t rawWeight(1.);
  const Float_t weightWeight(1.);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hPtJetRatio[iTrg]  -> Divide(hPtJetWeight[iTrg], hPtJetRaw[iTrg], weightWeight, rawWeight);
    hPtJetRatioN[iTrg] -> Divide(hPtJetWeightN[iTrg], hPtJetRawN[iTrg], weightWeight, rawWeight);
  }
  cout << "    Took ratios." << endl;


  // create directories
  TDirectory *dPi0 = (TDirectory*) fOut -> mkdir("pi0");
  TDirectory *dGam = (TDirectory*) fOut -> mkdir("gam");
  cout << "    Made directories." << endl;

  // save histograms
  fOut             -> cd();
  hWeight          -> Write();
  fQtTrg           -> Write();
  dPi0             -> cd();
  hNumTrgRaw[0]    -> Write();
  hNumTrgCorr[0]   -> Write();
  hNumTrgIn[0]     -> Write();
  hNumTrgOut[0]    -> Write();
  hEtTrgRaw[0]     -> Write();
  hEtTrgCorr[0]    -> Write();
  hEtTrgCalc[0]    -> Write();
  hPtJetRaw[0]     -> Write();
  hPtJetRawN[0]    -> Write();
  hPtJetWeight[0]  -> Write();
  hPtJetWeightN[0] -> Write();
  hPtJetRatio[0]   -> Write();
  hPtJetRatioN[0]  -> Write();
  dGam             -> cd();
  hNumTrgRaw[1]    -> Write();
  hNumTrgCorr[1]   -> Write();
  hNumTrgIn[1]     -> Write();
  hNumTrgOut[1]    -> Write();
  hEtTrgRaw[1]     -> Write();
  hEtTrgCorr[1]    -> Write();
  hEtTrgCalc[1]    -> Write();
  hPtJetRaw[1]     -> Write();
  hPtJetRawN[1]    -> Write();
  hPtJetWeight[1]  -> Write();
  hPtJetWeightN[1] -> Write();
  hPtJetRatio[1]   -> Write();
  hPtJetRatioN[1]  -> Write();
  cout << "    Saved histograms." << endl;

  // save and close files
  fOut      -> cd();
  fOut      -> Close();
  fJetIn    -> cd();
  fJetIn    -> Close();
  fWeightIn -> cd();
  fWeightIn -> Close();
  cout << "  Jet reweighting finished!\n" << endl;

}

// End ------------------------------------------------------------------------
