// 'MakeSummaryPlot_WithHist.C'
// Derek Anderson
// 06.02.2022
//
// Creates a plot summarizing
// the fully corrected recoil
// jet distributions.
//
// NOTE: use the 'JetRes' variable
//       to change how certain spectra
//       are handled.
//         JetRes = 0 -- R = 0.2
//         JetRes = 1 -- R = 0.5

#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TLine.h"
#include "TError.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"

// global constants
static const UInt_t NVtx(4);
static const UInt_t NErr(4);
static const UInt_t NPlot(2);
static const UInt_t NTrgPi(2);
static const UInt_t NTrgGa(3);
static const UInt_t NPtMax(40);
static const UInt_t NJetRes(2);
static const UInt_t NTrgTot(NTrgPi + NTrgGa);

// jet parameters
static const UInt_t JetRes(0);



void MakeSummaryPlot_WithHist() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distributions..." << endl;

  // io parameters
  const TString sOut("ppSummary.py8vs6_withRawHist_noScales_pTbinHuge.et920r02pi0vsGam.d2m6y2022.root");
  const TString sGraphStat("gIter1");
  const TString sGraphUnf("gIter1");
  const TString sGraphPyF("gIter1");
  const TString sGraphPyH("gIter1");
  const TString sHistPyP("hPtJetCorr_part_Recoil");
  const TString sHistPyG("Gam/hJetPtCorrG");
  const TString sInStat[NTrgTot] = {"et911r02pi0_py8vs6/stats.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et911r02pi0.d7m11y2021.root",
                                    "et1115r02pi0_py8vs6/stats.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et1115r02pi0.d7m11y2021.root",
                                    "et911r02gam_py8vs6/stats.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et911r02gam.d7m11y2021.root",
                                    "et1115r02gam_py8vs6/stats.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et1115r02gam.d7m11y2021.root",
                                    "et1520r02gam_py8vs6/stats.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et1520r02gam.d7m11y2021.root"};
  const TString sInUnf[NTrgTot]  = {"et911r02pi0_py8vs6/totalSys.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et911r02pi0.d7m11y2021.root",
                                    "et1115r02pi0_py8vs6/totalSys.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et1115r02pi0.d7m11y2021.root",
                                    "et911r02gam_py8vs6/totalSys.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et911r02gam.d7m11y2021.root",
                                    "et1115r02gam_py8vs6/totalSys.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et1115r02gam.d7m11y2021.root",
                                    "et1520r02gam_py8vs6/totalSys.shiftedDataPoints_smoothUnfoldRawDet_pTbinHuge.et1520r02gam.d7m11y2021.root"};
  const TString sInPyF[NTrgTot]  = {"et911r02pi0_py8vs6/pythia8.shiftedDataPoints_histWeight_smoothUnfoldRawDet_pTbinHuge.et911r02pi0.d8m11y2021.root",
                                    "et1115r02pi0_py8vs6/pythia8.shiftedDataPoints_histWeight_smoothUnfoldRawDet_pTbinHuge.et1115r02pi0.d8m11y2021.root",
                                    "et911r02gam_py8vs6/pythia8.shiftedDataPoints_histWeight_smoothUnfoldRawDet_pTbinHuge.et911r02gam.d8m11y2021.root",
                                    "et1115r02gam_py8vs6/pythia8.shiftedDataPoints_histWeight_smoothUnfoldRawDet_pTbinHuge.et1115r02gam.d8m11y2021.root",
                                    "et1520r02gam_py8vs6/pythia8.shiftedDataPoints_histWeight_smoothUnfoldRawDet_pTbinHuge.et1520r02gam.d8m11y2021.root"};
  const TString sInPyH[NTrgTot]  = {"et911r02pi0_py8vs6/pythia6.shiftedDataPoints_forPy8vsPy6_pTbinOne.et6100x911r02pi0.d19m5y2022.root",
                                    "et1115r02pi0_py8vs6/pythia6.shiftedDataPoints_forPy8vsPy6_pTbinOne.et6100x1115r02pi0.d19m5y2022.root",
                                    "et911r02gam_py8vs6/pythia6.shiftedDataPoints_forPy8vsPy6_pTbinOne.et6100x911r02gam.d19m5y2022.root",
                                    "et1115r02gam_py8vs6/pythia6.shiftedDataPoints_forPy8vsPy6_pTbinOne.et6100x1115r02gam.d19m5y2022.root",
                                    "et1520r02gam_py8vs6/pythia6.shiftedDataPoints_forPy8vsPy6_pTbinOne.et6100x1520r02gam.d19m5y2022.root"};
  const TString sInHistP[NTrgPi] = {"et911r02pi0_py8vs6/PYTHIA6STAR_AllPthat_Pi0Jet_R02_9to11_28Mar2022.root",
                                     "et1115r02pi0_py8vs6/PYTHIA6STAR_AllPthat_Pi0Jet_R02_11to15_28Mar2022.root"};
  const TString sInHistG[NTrgGa] = {"et911r02gam_py8vs6/pp200py6par.forWeightCheck_withWeightsDefinitelyHist_pTbinOne.et6100x911r02gam.d11m2y2022.root",
                                     "et1115r02gam_py8vs6/pp200py6par.forWeightCheck_withWeightsDefinitelyHist_pTbinOne.et6100x1115r02gam.d11m2y2022.root",
                                     "et1520r02gam_py8vs6/pp200py6par.forWeightCheck_withWeightsDefinitelyHist_pTbinOne.et6100x1520r02gam.d11m2y2022.root"};

  // text parameters
  const TString sSys("p+p #sqrt{s} = 200 GeV");
  const TString sTyp("#pi^{0}/#gamma_{dir}+jet, anti-k_{T}");
  const TString sEtc("R = 0.2");

  // plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T,jet}^{reco,ch} [GeV/c]");
  const TString  sTitleY("(1/N^{trg}) d^{2}N_{jets}/(dp_{T,jet}^{reco,ch} d#eta_{jet}) [GeV/c]^{-1}");
  const Double_t xPlotRange[NPlot] = {0., 30.};
  const Double_t yPlotRange[NPlot] = {0.000003, 3333.};

  // graph parameters
  const Float_t xMinJet(0.);
  const TString sNameStat[NTrgTot] = {"gStat911pi",     "gStat1115pi",    "gStat911ga",     "gStat1115ga",     "gStat1520ga"};
  const TString sNameUnf[NTrgTot]  = {"gSys911pi",      "gSys1115pi",     "gSys911ga",      "gSys1115ga",      "gStat1520ga"};
  const TString sPlotStat[NTrgTot] = {"gPlotStat911pi", "gPlotSys1115pi", "gPlotStat911ga", "gPlotStat1115ga", "gPlotStat1520ga"};
  const TString sPlotUnf[NTrgTot]  = {"gPlotSys911pi",  "gPlotSys1115pi", "gPlotSys911ga",  "gPlotSys1115ga",  "gPlotStat1520ga"};
  const TString sLabels[NTrgTot]   = {"#pi^{0}: E_{T}^{trig} #in (9, 11) GeV",
                                      "#pi^{0}: E_{T}^{trig} #in (11, 15) GeV",
                                      "#gamma_{dir}: E_{T}^{trig} #in (9, 11) GeV",
                                      "#gamma_{dir}: E_{T}^{trig} #in (11, 15) GeV",
                                      "#gamma_{dir}: E_{T}^{trig} #in (15, 20) GeV"};

  // trigger and jet parameters
  const UInt_t  fColSt[NJetRes][NTrgTot]      = {{863, 843, 803, 893, 883}, {863, 843, 803, 893, 883}};
  const UInt_t  fColUn[NJetRes][NTrgTot]      = {{856, 846, 806, 896, 886}, {856, 846, 806, 896, 886}};
  const UInt_t  fColOu[NJetRes][NTrgTot]      = {{604, 436, 404, 636, 620}, {604, 436, 404, 636, 620}};
  const UInt_t  fMarSt[NJetRes][NTrgTot]      = {{20, 22, 21, 23, 29}, {20, 22, 21, 23, 29}};
  const UInt_t  fMarUn[NJetRes][NTrgTot]      = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}};
  const UInt_t  fMarOu[NJetRes][NTrgTot]      = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}};
  const UInt_t  nJetPtBins[NJetRes][NTrgTot]  = {{10, 10, 6, 7, 8}, {10, 10, 6, 7, 8}};
  const Bool_t  goesToZero[NJetRes][NTrgTot]  = {{false, false, true, true, true}, {false, false, true, true, true}};
  const Float_t xMaxJetBins[NJetRes][NTrgTot] = {{30., 30., 11., 15., 20.}, {30., 30., 11., 15., 20.}};
  const Float_t fTrgScale[NJetRes][NTrgTot]   = {{1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}};

  // pythia parameters
  const TString sNameHistP[NTrgPi]            = {"hPy6_911pi",     "hPy6_1115pi"};
  const TString sNameHistG[NTrgGa]            = {"hPy6_911ga",     "hPy6_1115ga",     "hPy6_1520ga"};
  const TString sNamePyF[NTrgTot]             = {"gPy8_911pi",     "gPy8_1115pi",     "gPy8_911ga",     "gPy8_1115ga",     "gPy8_1520ga"};
  const TString sNamePyH[NTrgTot]             = {"gPy6_911pi",     "gPy6_1115pi",     "gPy6_911ga",     "gPy6_1115ga",     "gPy6_1520ga"};
  const TString sPlotPyF[NTrgTot]             = {"gPlotPy8_911pi", "gPlotPy8_1115pi", "gPlotPy8_911ga", "gPlotPy8_1115ga", "gPlotPy8_1520ga"};
  const TString sPlotPyH[NTrgTot]             = {"gPlotPy6_911pi", "gPlotPy6_1115pi", "gPlotPy6_911ga", "gPlotPy6_1115ga", "gPlotPy6_1520ga"};
  const UInt_t  fColPyF[NJetRes][NTrgTot]     = {{601, 433, 401, 633, 617}, {601, 433, 401, 633, 617}};
  const UInt_t  fColPyH[NJetRes][NTrgTot]     = {{593, 425, 393, 625, 609}, {593, 425, 393, 625, 609}};
  const UInt_t  fColHistP[NJetRes][NTrgPi]    = {{592, 424}, {592, 424}};
  const UInt_t  fColHistG[NJetRes][NTrgGa]    = {{392, 624, 608}, {392, 624, 608}};
  const UInt_t  fMarPyF[NJetRes][NTrgTot]     = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}};
  const UInt_t  fMarPyH[NJetRes][NTrgTot]     = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}};
  const UInt_t  fMarHistP[NJetRes][NTrgGa]    = {{24, 26}, {24, 26}};
  const UInt_t  fMarHistG[NJetRes][NTrgGa]    = {{25, 32, 30}, {25, 32, 30}};
  const UInt_t  nJetPyBinsF[NJetRes][NTrgTot] = {{10, 10, 10, 10, 10}, {10, 10, 8, 10, 10}};
  const UInt_t  nJetPyBinsH[NJetRes][NTrgTot] = {{30, 30, 30, 30, 30}, {30, 30, 20, 30, 30}};
  const UInt_t  xMaxPyBinsF[NJetRes][NTrgTot] = {{30., 30., 30., 30., 30.}, {30., 30., 20., 30., 30.}};
  const UInt_t  xMaxPyBinsH[NJetRes][NTrgTot] = {{30., 30., 30., 30., 30.}, {30., 30., 20., 30., 30.}};

  // style parameters
  const UInt_t  fLinStat(1);
  const UInt_t  fLinUnf(1);
  const UInt_t  fLinDet(1);
  const UInt_t  fLinOut(1);
  const UInt_t  fLinArr(1);
  const UInt_t  fLinPyF(1);
  const UInt_t  fLinPyH(9);
  const UInt_t  fLinHistP(1);
  const UInt_t  fLinHistG(1);
  const UInt_t  fWidStat(2);
  const UInt_t  fWidUnf(2);
  const UInt_t  fWidDet(2);
  const UInt_t  fWidOut(2);
  const UInt_t  fWidArr(2);
  const UInt_t  fWidPyF(2);
  const UInt_t  fWidPyH(2);
  const UInt_t  fWidHistP(1);
  const UInt_t  fWidHistG(1);
  const UInt_t  fFilStat(1001);
  const UInt_t  fFilUnf(1001);
  const UInt_t  fFilDet(1001);
  const UInt_t  fFilOut(0);
  const UInt_t  fFilPy(0);
  const UInt_t  fFilHistP(0);
  const UInt_t  fFilHistG(0);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  const Float_t fSizArr(0.05);
  const Float_t fSizMar(1.25);
  const Float_t fSizHat(1.5);

  // select trigger and jets
  UInt_t  fColStat[NTrgTot];
  UInt_t  fColUnf[NTrgTot];
  UInt_t  fColPyWF[NTrgTot];
  UInt_t  fColPyWH[NTrgTot];
  UInt_t  fColOut[NTrgTot];
  UInt_t  fColHistPyP[NTrgPi];
  UInt_t  fColHistPyG[NTrgGa];
  UInt_t  fMarStat[NTrgTot];
  UInt_t  fMarUnf[NTrgTot];
  UInt_t  fMarPyWF[NTrgTot];
  UInt_t  fMarPyWH[NTrgTot];
  UInt_t  fMarOut[NTrgTot];
  UInt_t  fMarHistPyP[NTrgPi];
  UInt_t  fMarHistPyG[NTrgGa];
  UInt_t  nJetPoints[NTrgTot];
  UInt_t  nPyPointsF[NTrgTot];
  UInt_t  nPyPointsH[NTrgTot];
  Bool_t  endIsArrow[NTrgTot];
  Float_t xMaxJet[NTrgTot];
  Float_t xMaxPyF[NTrgTot];
  Float_t xMaxPyH[NTrgTot];
  Float_t fScalers[NTrgTot];
  TString sTrg("");
  TString sJet("");
  switch (JetRes) {
    case 0:
      sTrg = "#pi^{0}/#gamma_{dir} trig.";
      sJet = "anti-k_{T}, R = 0.2";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColStat[iTrg]   = fColSt[0][iTrg];
        fColUnf[iTrg]    = fColUn[0][iTrg];
        fColPyWF[iTrg]   = fColPyF[0][iTrg];
        fColPyWH[iTrg]   = fColPyH[0][iTrg];
        fColOut[iTrg]    = fColOu[0][iTrg];
        fMarStat[iTrg]   = fMarSt[0][iTrg];
        fMarUnf[iTrg]    = fMarUn[0][iTrg];
        fMarPyWF[iTrg]   = fMarPyF[0][iTrg];
        fMarPyWH[iTrg]   = fMarPyH[0][iTrg];
        fMarOut[iTrg]    = fMarOu[0][iTrg];
        nJetPoints[iTrg] = nJetPtBins[0][iTrg];
        nPyPointsF[iTrg] = nJetPyBinsF[0][iTrg];
        nPyPointsH[iTrg] = nJetPyBinsH[0][iTrg];
        endIsArrow[iTrg] = goesToZero[0][iTrg];
        xMaxJet[iTrg]    = xMaxJetBins[0][iTrg];
        xMaxPyF[iTrg]    = xMaxPyBinsF[0][iTrg];
        xMaxPyH[iTrg]    = xMaxPyBinsH[0][iTrg];
        fScalers[iTrg]   = fTrgScale[0][iTrg];
      }
      for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
        fColHistPyP[iTrgP] = fColHistP[0][iTrgP];
        fMarHistPyP[iTrgP] = fMarHistP[0][iTrgP];
      }
      for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
        fColHistPyG[iTrgG] = fColHistG[0][iTrgG];
        fMarHistPyG[iTrgG] = fMarHistG[0][iTrgG];
      }
      break;
    case 1:
      sTrg = "#pi^{0}/#gamma_{dir} trig.";
      sJet = "anti-k_{T}, R = 0.5";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColStat[iTrg]   = fColSt[1][iTrg];
        fColUnf[iTrg]    = fColUn[1][iTrg];
        fColPyWF[iTrg]   = fColPyF[1][iTrg];
        fColPyWH[iTrg]   = fColPyH[1][iTrg];
        fColOut[iTrg]    = fColOu[1][iTrg];
        fMarStat[iTrg]   = fMarSt[1][iTrg];
        fMarUnf[iTrg]    = fMarUn[1][iTrg];
        fMarPyWF[iTrg]   = fMarPyF[1][iTrg];
        fMarPyWH[iTrg]   = fMarPyH[1][iTrg];
        fMarOut[iTrg]    = fMarOu[1][iTrg];
        nJetPoints[iTrg] = nJetPtBins[1][iTrg];
        nPyPointsF[iTrg] = nJetPyBinsF[1][iTrg];
        nPyPointsH[iTrg] = nJetPyBinsH[1][iTrg];
        endIsArrow[iTrg] = goesToZero[1][iTrg];
        xMaxJet[iTrg]    = xMaxJetBins[1][iTrg];
        xMaxPyF[iTrg]    = xMaxPyBinsF[1][iTrg];
        xMaxPyH[iTrg]    = xMaxPyBinsH[1][iTrg];
        fScalers[iTrg]   = fTrgScale[1][iTrg];
      }
      for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
        fColHistPyP[iTrgP] = fColHistP[1][iTrgP];
        fMarHistPyP[iTrgP] = fMarHistP[1][iTrgP];
      }
      for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
        fColHistPyG[iTrgG] = fColHistG[1][iTrgG];
        fMarHistPyG[iTrgG] = fMarHistG[1][iTrgG];
      }
      break;
    default:
      sTrg = "#pi^{0}/#gamma_{dir} trig.";
      sJet = "anti-k_{T}, R = 0.2";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColStat[iTrg]   = fColSt[0][iTrg];
        fColUnf[iTrg]    = fColUn[0][iTrg];
        fColPyWF[iTrg]   = fColPyF[0][iTrg];
        fColPyWH[iTrg]   = fColPyH[0][iTrg];
        fColOut[iTrg]    = fColOu[0][iTrg];
        fMarStat[iTrg]   = fMarSt[0][iTrg];
        fMarUnf[iTrg]    = fMarUn[0][iTrg];
        fMarPyWF[iTrg]   = fMarPyF[0][iTrg];
        fMarPyWH[iTrg]   = fMarPyH[0][iTrg];
        fMarOut[iTrg]    = fMarOu[0][iTrg];
        nJetPoints[iTrg] = nJetPtBins[0][iTrg];
        nPyPointsF[iTrg] = nJetPyBinsF[0][iTrg];
        nPyPointsH[iTrg] = nJetPyBinsH[0][iTrg];
        endIsArrow[iTrg] = goesToZero[0][iTrg];
        xMaxJet[iTrg]    = xMaxJetBins[0][iTrg];
        xMaxPyF[iTrg]    = xMaxPyBinsF[0][iTrg];
        xMaxPyH[iTrg]    = xMaxPyBinsH[0][iTrg];
        fScalers[iTrg]   = fTrgScale[0][iTrg];
      }
      for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
        fColHistPyP[iTrgP] = fColHistP[0][iTrgP];
        fMarHistPyP[iTrgP] = fMarHistP[0][iTrgP];
      }
      for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
        fColHistPyG[iTrgG] = fColHistG[0][iTrgG];
        fMarHistPyG[iTrgG] = fMarHistG[0][iTrgG];
      }
      break;
  }  // end switch-case

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }

  TFile *fInStat[NTrgTot];
  TFile *fInUnf[NTrgTot];
  TFile *fInPyF[NTrgTot];
  TFile *fInPyH[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    fInStat[iTrg] = new TFile(sInStat[iTrg].Data(), "read");
    fInUnf[iTrg]  = new TFile(sInUnf[iTrg].Data(), "read");
    fInPyF[iTrg]  = new TFile(sInPyF[iTrg].Data(), "read");
    fInPyH[iTrg]  = new TFile(sInPyH[iTrg].Data(), "read");
    if (!fInStat[iTrg] || !fInUnf[iTrg] || !fInPyF[iTrg] || !fInPyH[iTrg]) {
      cerr << "PANIC: couldn't open an input file!\n"
           << "       fInStat[" << iTrg << "] = " << fInStat[iTrg] << ", fInUnf[" << iTrg << "] = " << fInUnf[iTrg] << ", fInPyF[" << iTrg << "] = " << fInPyF[iTrg] << ", fInPyH[" << iTrg << "] = " << fInPyH[iTrg] << "\n"
           << endl;
      return;
    }
  }

  TFile *fInHistP[NTrgPi];
  TFile *fInHistG[NTrgGa];
  for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
    fInHistP[iTrgP] = new TFile(sInHistP[iTrgP].Data(), "read");
    if (!fInHistP[iTrgP]) {
      cerr << "PANIC: couldn't open an input pi0 histogram file!\n"
           << "       fInHistP[" << iTrgP << "] = " << fInHistP[iTrgP] << "\n"
           << endl;
      return;
    }
  }
  for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
    fInHistG[iTrgG] = new TFile(sInHistG[iTrgG].Data(), "read");
    if (!fInHistG[iTrgG]) {
      cerr << "PANIC: couldn't open an input gamma histogram file!\n"
           << "       fInHistG[" << iTrgG << "] = " << fInHistG[iTrgG] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input graphs
  TGraphAsymmErrors *gStat[NTrgTot];
  TGraphAsymmErrors *gUnf[NTrgTot];
  TGraphAsymmErrors *gPyF[NTrgTot];
  TGraphAsymmErrors *gPyH[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gStat[iTrg] = (TGraphAsymmErrors*) fInStat[iTrg] -> Get(sGraphStat.Data());
    gUnf[iTrg]  = (TGraphAsymmErrors*) fInUnf[iTrg]  -> Get(sGraphUnf.Data());
    gPyF[iTrg]  = (TGraphAsymmErrors*) fInPyF[iTrg]  -> Get(sGraphPyF.Data());
    gPyH[iTrg]  = (TGraphAsymmErrors*) fInPyH[iTrg]  -> Get(sGraphPyH.Data());
    if (!gStat[iTrg] || !gUnf[iTrg] || !gPyF[iTrg] || !gPyH[iTrg]) {
      cerr << "PANIC: couldn't open an input graph!\n"
           << "       gStat[" << iTrg << "] = " << gStat[iTrg] << ", gUnf[" << iTrg << "] = " << gUnf[iTrg] << ", gPyF[" << iTrg << "] = " << gPyF[iTrg] << ", gPyH[" << iTrg << "] = " << gPyH[iTrg] << "\n"
           << endl;
      return;
    }
    gStat[iTrg] -> SetName(sNameStat[iTrg].Data());
    gUnf[iTrg]  -> SetName(sNameUnf[iTrg].Data());
    gPyF[iTrg]  -> SetName(sNamePyF[iTrg].Data());
    gPyH[iTrg]  -> SetName(sNamePyH[iTrg].Data());
  }
  cout << "    Grabbed graphs." << endl;

  // grab input histograms
  TH1D *hPythiaP[NTrgPi];
  TH1D *hPythiaG[NTrgGa];
  for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
    hPythiaP[iTrgP] = (TH1D*) fInHistP[iTrgP] -> Get(sHistPyP.Data());
    if (!hPythiaP[iTrgP]) {
      cerr << "PANIC: couldn't grab an input pi0 histogram!\n"
           << "       hPythiaP[" << iTrgP << "] = " << hPythiaP[iTrgP] << "\n"
           << endl;
      return;
    }
    hPythiaP[iTrgP] -> SetName(sNameHistP[iTrgP].Data());
  }
  for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
    hPythiaG[iTrgG] = (TH1D*) fInHistG[iTrgG] -> Get(sHistPyG.Data());
    if (!hPythiaG[iTrgG]) {
      cerr << "PANIC: couldn't grab an input gamma histogram!\n"
           << "       hPythiaG[" << iTrgG << "] = " << hPythiaG[iTrgG] << "\n"
           << endl;
      return;
    }
    hPythiaG[iTrgG] -> SetName(sNameHistG[iTrgG].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // create graphs for plotting
  Double_t xPlot[NTrgTot][NPtMax];
  Double_t yPlot[NTrgTot][NPtMax];
  Double_t yPlotLoU[NTrgTot][NPtMax];
  Double_t yPlotHiU[NTrgTot][NPtMax];
  Double_t xPlotPyF[NTrgTot][NPtMax];
  Double_t xPlotPyH[NTrgTot][NPtMax];
  Double_t yPlotPyF[NTrgTot][NPtMax];
  Double_t yPlotPyH[NTrgTot][NPtMax];
  Double_t ePlotZe[NTrgTot][NErr][NPtMax];
  Double_t ePlotSt[NTrgTot][NErr][NPtMax];
  Double_t ePlotUn[NTrgTot][NErr][NPtMax];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    for (UInt_t iPt = 0; iPt < NPtMax; iPt++) {
      xPlot[iTrg][iPt]    = 0.;
      yPlot[iTrg][iPt]    = 0.;
      yPlotLoU[iTrg][iPt] = 0.;
      yPlotHiU[iTrg][iPt] = 0.;
      xPlotPyF[iTrg][iPt] = 0.;
      yPlotPyH[iTrg][iPt] = 0.;
      for (UInt_t iErr = 0; iErr < NErr; iErr++) {
        ePlotZe[iTrg][iErr][iPt] = 0.;
        ePlotSt[iTrg][iErr][iPt] = 0.;
        ePlotUn[iTrg][iErr][iPt] = 0.;
      }
    }  // end bin loop
  }  // end trigger loop

  //TLine             *lEndPoint[NTrgTot];
  //TArrow            *aEndArrow[NTrgTot];
  TGraph            *gPlotLoU[NTrgTot];
  TGraph            *gPlotHiU[NTrgTot];
  TGraph            *gPlotPyF[NTrgTot];
  TGraph            *gPlotPyH[NTrgTot];
  TGraphAsymmErrors *gPlotSt[NTrgTot];
  TGraphAsymmErrors *gPlotUn[NTrgTot];
  cout << "    Reading in data points..." << endl;

  UInt_t iJetPoint(0);
  UInt_t iPyPointF(0);
  UInt_t iPyPointH(0);
  UInt_t nJetCheck(0);
  UInt_t nPyCheckF(0);
  UInt_t nPyCheckH(0);
  UInt_t nPoints(0);
  UInt_t nSimPointsF(0);
  UInt_t nSimPointsH(0);
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {

    // loop over data points
    iJetPoint = 0;
    nJetCheck = 0;
    nPoints   = gStat[iTrg] -> GetN();
    for (UInt_t iPoint = 0; iPoint < nPoints; iPoint++) {

      // get input points
      Double_t xStat(0.);
      Double_t xUnf(0.);
      Double_t xDet(0.);
      Double_t yStat(0.);
      Double_t yUnf(0.);
      Double_t yDet(0.);
      gStat[iTrg] -> GetPoint(iPoint, xStat, yStat);
      gUnf[iTrg]  -> GetPoint(iPoint, xUnf, yUnf);

      // check if data point is in range
      const Bool_t isInRange = ((xStat > xMinJet) && (xStat < xMaxJet[iTrg]));
      if (!isInRange) {
        continue;
      } else {
        nJetCheck++;
      }

      // check if stat and sys yields match
      const Bool_t unfYieldsMatch = (yUnf == yStat);
      if (!unfYieldsMatch) {
        cerr << "      WARNING: stat. and sys. yields don't match at bin #" << iPoint << " for trigger bin #" << iTrg << "!\n"
             << "               yStat = " << yStat << ", yUnf = " << yUnf
             << endl;
      }
      xPlot[iTrg][iJetPoint] = xStat;
      yPlot[iTrg][iJetPoint] = yStat * fScalers[iTrg];

      // get errors
      ePlotSt[iTrg][0][iJetPoint] = gStat[iTrg] -> GetErrorXlow(iPoint);
      ePlotUn[iTrg][0][iJetPoint] = gUnf[iTrg]  -> GetErrorXlow(iPoint);
      ePlotSt[iTrg][1][iJetPoint] = gStat[iTrg] -> GetErrorXhigh(iPoint);
      ePlotUn[iTrg][1][iJetPoint] = gUnf[iTrg]  -> GetErrorXhigh(iPoint);
      ePlotSt[iTrg][2][iJetPoint] = fScalers[iTrg] * (gStat[iTrg] -> GetErrorYlow(iPoint));
      ePlotUn[iTrg][2][iJetPoint] = fScalers[iTrg] * (gUnf[iTrg]  -> GetErrorYlow(iPoint));
      ePlotSt[iTrg][3][iJetPoint] = fScalers[iTrg] * (gStat[iTrg] -> GetErrorYhigh(iPoint));
      ePlotUn[iTrg][3][iJetPoint] = fScalers[iTrg] * (gUnf[iTrg]  -> GetErrorYhigh(iPoint));

      // get upper and lower bands
      yPlotLoU[iTrg][iJetPoint] = yPlot[iTrg][iJetPoint] - ePlotUn[iTrg][2][iJetPoint];
      yPlotHiU[iTrg][iJetPoint] = yPlot[iTrg][iJetPoint] + ePlotUn[iTrg][3][iJetPoint];

      // check if arrow needs to be created
      const Bool_t isLastPoint = ((iJetPoint + 1) == nJetPoints[iTrg]);
      if (endIsArrow[iTrg] && isLastPoint) {

        // get start and stop points
        const Double_t xStart = xPlot[iTrg][iJetPoint] - ePlotSt[iTrg][0][iJetPoint];
        const Double_t xStop  = xPlot[iTrg][iJetPoint] + ePlotSt[iTrg][1][iJetPoint];
        const Double_t yLoSt  = yPlot[iTrg][iJetPoint] - ePlotSt[iTrg][2][iJetPoint];
        const Double_t yLoUn  = yPlot[iTrg][iJetPoint] - ePlotUn[iTrg][2][iJetPoint];
        const Double_t yHiSt  = yPlot[iTrg][iJetPoint] + ePlotSt[iTrg][3][iJetPoint];
        const Double_t yHiUn  = yPlot[iTrg][iJetPoint] + ePlotUn[iTrg][3][iJetPoint];
        const Double_t yLoSy  = yLoUn;
        const Double_t yHiSy  = yHiUn;
        const Double_t yStart = TMath::Max(yHiSt, yHiSy);
        const Double_t yStop  = TMath::Max(yLoSt, yLoSy);

/* will fix later
        // make line
        lEndPoint[iTrg] = new TLine(xStart, yStart, xStop, yStart);
        lEndPoint[iTrg] -> SetLineColor(fColStat[iTrg]);
        lEndPoint[iTrg] -> SetLineWidth(fWidArr);
        lEndPoint[iTrg] -> SetLineStyle(fLinArr);

        // make arrow
        aEndArrow[iTrg] = new TArrow(xPlot[iTrg][iJetPoint], yStart, xPlot[iTrg][iJetPoint], yStop, fSizArr, ">");
        aEndArrow[iTrg] -> SetLineColor(fColStat[iTrg]);
        aEndArrow[iTrg] -> SetLineWidth(fWidArr);
        aEndArrow[iTrg] -> SetLineStyle(fLinArr);
*/
      }
      iJetPoint++;
    }  // end data point loop

    // loop over fit pythia points
    iPyPointF   = 0;
    nPyCheckF   = 0;
    nSimPointsF = gPyF[iTrg] -> GetN();
    for (UInt_t iSimPointF = 0; iSimPointF < nSimPointsF; iSimPointF++) {

      // get fit pythia points
      Double_t xPyF(0.);
      Double_t yPyF(0.);
      gPyF[iTrg] -> GetPoint(iSimPointF, xPyF, yPyF);

      // check if fit pythia point is in range
      const Bool_t isInPyRangeF = ((xPyF > xMinJet) && (xPyF < xMaxPyF[iTrg]));
      if (!isInPyRangeF) {
        continue;
      } else {
        nPyCheckF++;
      }

      // set fit pythia points
      xPlotPyF[iTrg][iPyPointF] = xPyF;
      yPlotPyF[iTrg][iPyPointF] = yPyF * fScalers[iTrg];
      iPyPointF++;
    }  // end fit pythia point loop

    // loop over hist pythia points
    iPyPointH   = 0;
    nPyCheckH   = 0;
    nSimPointsH = gPyH[iTrg] -> GetN();
    for (UInt_t iSimPointH = 0; iSimPointH < nSimPointsH; iSimPointH++) {

      // get hist pythia points
      Double_t xPyH(0.);
      Double_t yPyH(0.);
      gPyH[iTrg] -> GetPoint(iSimPointH, xPyH, yPyH);

      // check if hist pythia point is in range
      const Bool_t isInPyRangeH = ((xPyH > xMinJet) && (xPyH < xMaxPyH[iTrg]));
      if (!isInPyRangeH) {
        continue;
      } else {
        nPyCheckH++;
      }

      // set hist pythia points
      xPlotPyH[iTrg][iPyPointH] = xPyH;
      yPlotPyH[iTrg][iPyPointH] = yPyH * fScalers[iTrg];
      iPyPointH++;
    }  // end hist pythia point loop

    // create graphs
    const Int_t nGraph  = (Int_t) nJetPoints[iTrg];
    const Int_t nGraphF = (Int_t) nPyPointsF[iTrg];
    const Int_t nGraphH = (Int_t) nPyPointsH[iTrg];
    gPlotSt[iTrg]  = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotZe[iTrg][0], ePlotZe[iTrg][1], ePlotSt[iTrg][2], ePlotSt[iTrg][3]);
    gPlotUn[iTrg]  = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotZe[iTrg][0], ePlotZe[iTrg][1], ePlotUn[iTrg][2], ePlotUn[iTrg][3]);
    gPlotLoU[iTrg] = new TGraph(nGraph, xPlot[iTrg], yPlotLoU[iTrg]);
    gPlotHiU[iTrg] = new TGraph(nGraph, xPlot[iTrg], yPlotHiU[iTrg]);
    gPlotPyF[iTrg] = new TGraph(nGraphF, xPlotPyF[iTrg], yPlotPyF[iTrg]);
    gPlotPyH[iTrg] = new TGraph(nGraphH, xPlotPyH[iTrg], yPlotPyH[iTrg]);

/* will fix later
    if (!endIsArrow[iTrg]) {
      gPlotOu[iTrg] = new TGraph(nGraph, xPlot[iTrg], yPlot[iTrg]);
      gPlotSt[iTrg] = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotSt[iTrg][2], ePlotSt[iTrg][3]);
      gPlotUn[iTrg] = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotUn[iTrg][2], ePlotUn[iTrg][3]);
    } else {
      gPlotOu[iTrg] = new TGraphAsymmErrors(nGraph - 1, xPlot[iTrg], yPlot[iTrg]);
      gPlotSt[iTrg] = new TGraphAsymmErrors(nGraph - 1, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotSt[iTrg][2], ePlotSt[iTrg][3]);
      gPlotUn[iTrg] = new TGraphAsymmErrors(nGraph - 1, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotUn[iTrg][2], ePlotUn[iTrg][3]);
    }
*/

    // set data and pythia names
    gPlotSt[iTrg]  -> SetName(sPlotStat[iTrg].Data());
    gPlotUn[iTrg]  -> SetName(sPlotUnf[iTrg].Data());
    gPlotPyF[iTrg] -> SetName(sPlotPyF[iTrg].Data());
    gPlotPyH[iTrg] -> SetName(sPlotPyH[iTrg].Data());

    // create bound names
    TString sNameLoU = gPlotUn[iTrg] -> GetName();
    TString sNameHiU = gPlotUn[iTrg] -> GetName();
    sNameLoU.Append("_LoErr");
    sNameHiU.Append("_HiErr");

    // set bound names
    gPlotLoU[iTrg] -> SetName(sNameLoU.Data());
    gPlotHiU[iTrg] -> SetName(sNameHiU.Data());

    // check if no. of points is reasonable
    if (nJetCheck != nJetPoints[iTrg]) {
      cerr << "      WARNING: no. of accepted data points is off in trigger bin #" << iTrg << "!" << endl;
    }
    if (nPyCheckF != nPyPointsF[iTrg]) {
      cerr << "      WARNING: no. of accepted fit-weighted pythia points is off in trigger bin #" << iTrg << "!" << endl;
    }
    if (nPyCheckH != nPyPointsH[iTrg]) {
      cerr << "      WARNING: no. of accepted hist-weighted pythia points is off in trigger bin #" << iTrg << "!" << endl;
    }
    cout << "      Created plotting graphs for bin #" << iTrg << "..." << endl;

  }  // end trigger loop
  cout << "    Made graphs for plotting." << endl;

  // set styles
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gPlotSt[iTrg]  -> SetLineColor(fColStat[iTrg]);
    gPlotSt[iTrg]  -> SetLineStyle(fLinStat);
    gPlotSt[iTrg]  -> SetLineWidth(fWidStat);
    gPlotSt[iTrg]  -> SetFillColor(fColStat[iTrg]);
    gPlotSt[iTrg]  -> SetFillStyle(fFilStat);
    gPlotSt[iTrg]  -> SetMarkerColor(fColStat[iTrg]);
    gPlotSt[iTrg]  -> SetMarkerStyle(fMarStat[iTrg]);
    gPlotSt[iTrg]  -> SetMarkerSize(fSizMar);
    gPlotSt[iTrg]  -> SetTitle(sTitle.Data());
    gPlotSt[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotSt[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotSt[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotSt[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotSt[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotSt[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    gPlotSt[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    gPlotSt[iTrg]  -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotSt[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotSt[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotSt[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotSt[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotSt[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    gPlotSt[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    gPlotUn[iTrg]  -> SetLineColor(fColUnf[iTrg]);
    gPlotUn[iTrg]  -> SetLineStyle(fLinUnf);
    gPlotUn[iTrg]  -> SetLineWidth(fWidUnf);
    gPlotUn[iTrg]  -> SetFillColor(fColUnf[iTrg]);
    gPlotUn[iTrg]  -> SetFillStyle(fFilUnf);
    gPlotUn[iTrg]  -> SetMarkerColor(fColUnf[iTrg]);
    gPlotUn[iTrg]  -> SetMarkerStyle(fMarUnf[iTrg]);
    gPlotUn[iTrg]  -> SetMarkerSize(fSizMar);
    gPlotUn[iTrg]  -> SetTitle(sTitle.Data());
    gPlotUn[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotUn[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotUn[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotUn[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotUn[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotUn[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    gPlotUn[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    gPlotUn[iTrg]  -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotUn[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotUn[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotUn[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotUn[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotUn[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    gPlotUn[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    gPlotPyF[iTrg] -> SetLineColor(fColPyWF[iTrg]);
    gPlotPyF[iTrg] -> SetLineStyle(fLinPyF);
    gPlotPyF[iTrg] -> SetLineWidth(fWidPyF);
    gPlotPyF[iTrg] -> SetFillColor(fColPyWF[iTrg]);
    gPlotPyF[iTrg] -> SetFillStyle(fFilPy);
    gPlotPyF[iTrg] -> SetMarkerColor(fColPyWF[iTrg]);
    gPlotPyF[iTrg] -> SetMarkerStyle(fMarPyWF[iTrg]);
    gPlotPyF[iTrg] -> SetMarkerSize(fSizMar);
    gPlotPyF[iTrg] -> SetTitle(sTitle.Data());
    gPlotPyF[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotPyF[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotPyF[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotPyF[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotPyF[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotPyF[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    gPlotPyF[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    gPlotPyF[iTrg] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotPyF[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotPyF[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotPyF[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotPyF[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotPyF[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    gPlotPyF[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    gPlotPyH[iTrg] -> SetLineColor(fColPyWH[iTrg]);
    gPlotPyH[iTrg] -> SetLineStyle(fLinPyH);
    gPlotPyH[iTrg] -> SetLineWidth(fWidPyH);
    gPlotPyH[iTrg] -> SetFillColor(fColPyWH[iTrg]);
    gPlotPyH[iTrg] -> SetFillStyle(fFilPy);
    gPlotPyH[iTrg] -> SetMarkerColor(fColPyWH[iTrg]);
    gPlotPyH[iTrg] -> SetMarkerStyle(fMarPyWH[iTrg]);
    gPlotPyH[iTrg] -> SetMarkerSize(fSizMar);
    gPlotPyH[iTrg] -> SetTitle(sTitle.Data());
    gPlotPyH[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotPyH[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotPyH[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotPyH[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotPyH[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotPyH[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    gPlotPyH[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    gPlotPyH[iTrg] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotPyH[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotPyH[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotPyH[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotPyH[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotPyH[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    gPlotPyH[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    gPlotLoU[iTrg] -> SetLineColor(fColOut[iTrg]);
    gPlotLoU[iTrg] -> SetLineStyle(fLinOut);
    gPlotLoU[iTrg] -> SetLineWidth(fWidOut);
    gPlotLoU[iTrg] -> SetFillColor(fColOut[iTrg]);
    gPlotLoU[iTrg] -> SetFillStyle(fFilOut);
    gPlotLoU[iTrg] -> SetMarkerColor(fColOut[iTrg]);
    gPlotLoU[iTrg] -> SetMarkerStyle(fMarOut[iTrg]);
    gPlotLoU[iTrg] -> SetMarkerSize(fSizMar);
    gPlotLoU[iTrg] -> SetTitle(sTitle.Data());
    gPlotLoU[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotLoU[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotLoU[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotLoU[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotLoU[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotLoU[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    gPlotLoU[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    gPlotLoU[iTrg] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotLoU[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotLoU[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotLoU[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotLoU[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotLoU[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    gPlotLoU[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    gPlotHiU[iTrg] -> SetLineColor(fColOut[iTrg]);
    gPlotHiU[iTrg] -> SetLineStyle(fLinOut);
    gPlotHiU[iTrg] -> SetLineWidth(fWidOut);
    gPlotHiU[iTrg] -> SetFillColor(fColOut[iTrg]);
    gPlotHiU[iTrg] -> SetFillStyle(fFilOut);
    gPlotHiU[iTrg] -> SetMarkerColor(fColOut[iTrg]);
    gPlotHiU[iTrg] -> SetMarkerStyle(fMarOut[iTrg]);
    gPlotHiU[iTrg] -> SetMarkerSize(fSizMar);
    gPlotHiU[iTrg] -> SetTitle(sTitle.Data());
    gPlotHiU[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotHiU[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotHiU[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotHiU[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotHiU[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotHiU[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    gPlotHiU[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    gPlotHiU[iTrg] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotHiU[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotHiU[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotHiU[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotHiU[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotHiU[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    gPlotHiU[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end trigger loop
  for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
    hPythiaP[iTrgP] -> SetLineColor(fColHistPyP[iTrgP]);
    hPythiaP[iTrgP] -> SetLineStyle(fLinHistP);
    hPythiaP[iTrgP] -> SetLineWidth(fWidHistP);
    hPythiaP[iTrgP] -> SetFillColor(fColHistPyP[iTrgP]);
    hPythiaP[iTrgP] -> SetFillStyle(fFilHistP);
    hPythiaP[iTrgP] -> SetMarkerColor(fColHistPyP[iTrgP]);
    hPythiaP[iTrgP] -> SetMarkerStyle(fMarHistPyP[iTrgP]);
    hPythiaP[iTrgP] -> SetMarkerSize(fSizMar);
    hPythiaP[iTrgP] -> SetTitle(sTitle.Data());
    hPythiaP[iTrgP] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPythiaP[iTrgP] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPythiaP[iTrgP] -> GetXaxis() -> SetTitleFont(fTxt);
    hPythiaP[iTrgP] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPythiaP[iTrgP] -> GetXaxis() -> SetLabelFont(fTxt);
    hPythiaP[iTrgP] -> GetXaxis() -> SetLabelSize(fLab);
    hPythiaP[iTrgP] -> GetXaxis() -> CenterTitle(fCnt);
    hPythiaP[iTrgP] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    hPythiaP[iTrgP] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPythiaP[iTrgP] -> GetYaxis() -> SetTitleFont(fTxt);
    hPythiaP[iTrgP] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPythiaP[iTrgP] -> GetYaxis() -> SetLabelFont(fTxt);
    hPythiaP[iTrgP] -> GetYaxis() -> SetLabelSize(fLab);
    hPythiaP[iTrgP] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
    hPythiaG[iTrgG] -> SetLineColor(fColHistPyG[iTrgG]);
    hPythiaG[iTrgG] -> SetLineStyle(fLinHistG);
    hPythiaG[iTrgG] -> SetLineWidth(fWidHistG);
    hPythiaG[iTrgG] -> SetFillColor(fColHistPyG[iTrgG]);
    hPythiaG[iTrgG] -> SetFillStyle(fFilHistG);
    hPythiaG[iTrgG] -> SetMarkerColor(fColHistPyG[iTrgG]);
    hPythiaG[iTrgG] -> SetMarkerStyle(fMarHistPyG[iTrgG]);
    hPythiaG[iTrgG] -> SetMarkerSize(fSizMar);
    hPythiaG[iTrgG] -> SetTitle(sTitle.Data());
    hPythiaG[iTrgG] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPythiaG[iTrgG] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPythiaG[iTrgG] -> GetXaxis() -> SetTitleFont(fTxt);
    hPythiaG[iTrgG] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPythiaG[iTrgG] -> GetXaxis() -> SetLabelFont(fTxt);
    hPythiaG[iTrgG] -> GetXaxis() -> SetLabelSize(fLab);
    hPythiaG[iTrgG] -> GetXaxis() -> CenterTitle(fCnt);
    hPythiaG[iTrgG] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    hPythiaG[iTrgG] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPythiaG[iTrgG] -> GetYaxis() -> SetTitleFont(fTxt);
    hPythiaG[iTrgG] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPythiaG[iTrgG] -> GetYaxis() -> SetLabelFont(fTxt);
    hPythiaG[iTrgG] -> GetYaxis() -> SetLabelSize(fLab);
    hPythiaG[iTrgG] -> GetYaxis() -> CenterTitle(fCnt);
  }

  // for legend
  const UInt_t fColPyWfLeg(1);
  const UInt_t fColPyWhLeg(922);
  const UInt_t fColHistLeg(920);
  const UInt_t fLinPlotLeg(1);
  const UInt_t fLinHistLeg(1);
  const UInt_t fWidPlotLeg(2);
  const UInt_t fWidHistLeg(1);
  const UInt_t fFilPlotLeg(0);
  const UInt_t fFilHistLeg(0);
  const UInt_t fMarHistLeg(24);

  TH1D   *hPyLeg   = (TH1D*)   hPythiaP[0] -> Clone();
  TGraph *gPyWfLeg = (TGraph*) gPlotPyF[0] -> Clone();
  TGraph *gPyWhLeg = (TGraph*) gPlotPyH[0] -> Clone();
  hPyLeg   -> SetFillColor(fColHistLeg);
  hPyLeg   -> SetFillStyle(fFilHistLeg);
  hPyLeg   -> SetLineColor(fColHistLeg);
  hPyLeg   -> SetLineStyle(fLinHistLeg);
  hPyLeg   -> SetLineWidth(fWidHistLeg);
  hPyLeg   -> SetMarkerColor(fColHistLeg);
  hPyLeg   -> SetMarkerStyle(fMarHistLeg);
  gPyWfLeg -> SetFillColor(fColPyWfLeg);
  gPyWfLeg -> SetLineColor(fColPyWfLeg);
  gPyWhLeg -> SetFillColor(fColPyWhLeg);
  gPyWhLeg -> SetLineColor(fColPyWhLeg);

  TGraphAsymmErrors *gPlotLeg[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gPlotLeg[iTrg] = (TGraphAsymmErrors*) gPlotSt[iTrg] -> Clone();
    gPlotLeg[iTrg] -> SetLineColor(fColUnf[iTrg]);
    gPlotLeg[iTrg] -> SetLineStyle(fLinPlotLeg);
    gPlotLeg[iTrg] -> SetLineWidth(fWidPlotLeg);
    gPlotLeg[iTrg] -> SetFillColor(fColUnf[iTrg]);
    gPlotLeg[iTrg] -> SetFillStyle(fFilPlotLeg);
  }
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  fAlnLe(12);
  const UInt_t  nColLe(3);
  const Float_t hObj(0.05);
  const Float_t hLeg(((NTrgTot - NTrgPi) + 1) * hObj);
  const Float_t yLeg(0.1 + hLeg);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.5, yLeg};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnLe);
  leg -> SetNColumns(nColLe);
  leg -> AddEntry(gPyWfLeg, "PYTHIA-8 (at Barycenter)", "fl");
  leg -> AddEntry(gPyWhLeg, "PYTHIA-6 (at Barycenter)", "fl");
  leg -> AddEntry(hPyLeg, "PYTHIA-6 (at bin center)", "pf");
  leg -> AddEntry(gPlotLeg[2], sLabels[2].Data(), "pf");
  leg -> AddEntry(gPlotLeg[3], sLabels[3].Data(), "pf");
  leg -> AddEntry(gPlotLeg[4], sLabels[4].Data(), "pf");
  leg -> AddEntry(gPlotLeg[0], sLabels[0].Data(), "pf");
  leg -> AddEntry(gPlotLeg[1], sLabels[1].Data(), "pf");
/* will fix later
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    leg -> AddEntry(gPlotUn[iTrg], sLabels[iTrg].Data(), "fl");
  }
*/
  cout << "    Made legend." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const UInt_t  fAlnTx(32);
  const UInt_t  fAlnLb(22);
  const UInt_t  fFonLb(62);
  const Float_t fTxtTx[NVtx] = {0.5, 0.5, 0.7, 0.65};
  TPaveText *txt = new TPaveText(fTxtTx[0], fTxtTx[1], fTxtTx[2], fTxtTx[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlnTx);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTyp.Data());
  txt -> AddText(sJet.Data());
  cout << "    Made text box." << endl;

  // create frame
  const UInt_t  nFrameX(200);
  const Float_t xFrame[NPlot] = {-100., 100.};
  fOut -> cd();

  TH1D *hFrame = new TH1D("hFrame", "", nFrameX, xFrame[0], xFrame[1]);
  hFrame -> SetDirectory(0);
  hFrame -> SetStats(0);
  hFrame -> SetTitle(sTitle.Data());
  hFrame -> SetTitleFont(fTxt);
  hFrame -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hFrame -> GetXaxis() -> SetTitle(sTitleX.Data());
  hFrame -> GetXaxis() -> SetTitleFont(fTxt);
  hFrame -> GetXaxis() -> SetTitleOffset(fOffX);
  hFrame -> GetXaxis() -> SetLabelFont(fTxt);
  hFrame -> GetXaxis() -> SetLabelSize(fLab);
  hFrame -> GetXaxis() -> CenterTitle(fCnt);
  hFrame -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
  hFrame -> GetYaxis() -> SetTitle(sTitleY.Data());
  hFrame -> GetYaxis() -> SetTitleFont(fTxt);
  hFrame -> GetYaxis() -> SetTitleOffset(fOffY);
  hFrame -> GetYaxis() -> SetLabelFont(fTxt);
  hFrame -> GetYaxis() -> SetLabelSize(fLab);
  hFrame -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Created frame." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.15);
  const Float_t xyFrame[NVtx] = {0., 0., 1., 1.};

  TCanvas *cSummary = new TCanvas("cCorrectedData", "", width, height);
  TPad    *pFrame   = new TPad("pFrame", "", xyFrame[0], xyFrame[1], xyFrame[2], xyFrame[3]); 
  cSummary    -> SetGrid(fGrid, fGrid);
  cSummary    -> SetTicks(fTick, fTick);
  cSummary    -> SetLogx(fLogX);
  cSummary    -> SetLogy(fLogY);
  cSummary    -> SetBorderMode(fMode);
  cSummary    -> SetBorderSize(fBord);
  cSummary    -> SetFrameBorderMode(fFrame);
  cSummary    -> SetLeftMargin(fMarginL);
  cSummary    -> SetRightMargin(fMarginR);
  cSummary    -> SetTopMargin(fMarginT);
  cSummary    -> SetBottomMargin(fMarginB);
  pFrame      -> SetGrid(fGrid, fGrid);
  pFrame      -> SetTicks(fTick, fTick);
  pFrame      -> SetLogx(fLogX);
  pFrame      -> SetLogy(fLogY);
  pFrame      -> SetBorderMode(fMode);
  pFrame      -> SetBorderSize(fBord);
  pFrame      -> SetFrameBorderMode(fFrame);
  pFrame      -> SetLeftMargin(fMarginL);
  pFrame      -> SetRightMargin(fMarginR);
  pFrame      -> SetTopMargin(fMarginT);
  pFrame      -> SetBottomMargin(fMarginB);
  cSummary    -> cd();
  pFrame      -> Draw();
  pFrame      -> cd();
  hFrame      -> Draw("9");
  // draw gamma curves
  gPlotUn[4]  -> Draw("Z3");
  gPlotLoU[4] -> Draw("L");
  gPlotHiU[4] -> Draw("L");
  gPlotSt[4]  -> Draw("P");
  gPlotUn[3]  -> Draw("Z3");
  gPlotLoU[3] -> Draw("L");
  gPlotHiU[3] -> Draw("L");
  gPlotSt[3]  -> Draw("P");
  gPlotUn[2]  -> Draw("Z3");
  gPlotLoU[2] -> Draw("L");
  gPlotHiU[2] -> Draw("L");
  gPlotSt[2]  -> Draw("P");
  // draw pi0 curves
  gPlotUn[1]  -> Draw("Z3");
  gPlotLoU[1] -> Draw("L");
  gPlotHiU[1] -> Draw("L");
  gPlotSt[1]  -> Draw("P");
  gPlotUn[0]  -> Draw("Z3");
  gPlotLoU[0] -> Draw("L");
  gPlotHiU[0] -> Draw("L");
  gPlotSt[0]  -> Draw("P");
  // draw pythia curves
  gPlotPyF[4] -> Draw("C");
  gPlotPyH[4] -> Draw("C");
  gPlotPyF[3] -> Draw("C");
  gPlotPyH[3] -> Draw("C");
  gPlotPyF[2] -> Draw("C");
  gPlotPyH[2] -> Draw("C");
  gPlotPyF[1] -> Draw("C");
  gPlotPyH[1] -> Draw("C");
  gPlotPyF[0] -> Draw("C");
  gPlotPyH[0] -> Draw("C");
  hPythiaG[2] -> Scale(fScalers[4]);
  hPythiaG[2] -> Draw("same");
  hPythiaG[1] -> Scale(fScalers[3]);
  hPythiaG[1] -> Draw("same");
  hPythiaG[0] -> Scale(fScalers[2]);
  hPythiaG[0] -> Draw("same");
  hPythiaP[1] -> Scale(fScalers[1]);
  hPythiaP[1] -> Draw("same");
  hPythiaP[1] -> Scale(fScalers[0]);
  hPythiaP[0] -> Draw("same");
/* will fix later
  for (Int_t iTrgPlot = (NTrgTot - 1); iTrgPlot > -1; iTrgPlot--) {
    if (iTrgPlot == (NTrgTot - 1)) {
      gPlotUn[iTrgPlot]  -> Draw("AZ3");
      gPlotLoU[iTrgPlot] -> Draw("L");
      gPlotHiU[iTrgPlot] -> Draw("L");
      gPlotPyF[iTrgPlot] -> Draw("L");
      gPlotPyH[iTrgPlot] -> Draw("L");
    } else {
      gPlotUn[iTrgPlot]  -> Draw("Z3");
      gPlotLoU[iTrgPlot] -> Draw("L");
      gPlotHiU[iTrgPlot] -> Draw("L");
      gPlotPyF[iTrgPlot] -> Draw("L");
      gPlotPyH[iTrgPlot] -> Draw("L");
    }
    gPlotSt[iTrgPlot] -> Draw("P");
    if (endIsArrow[iTrgPlot]) {
      lEndPoint[iTrgPlot] -> Draw();
      aEndArrow[iTrgPlot] -> Draw();
    }
  }
*/
  txt      -> Draw();
  leg      -> Draw();
  fOut     -> cd();
  cSummary -> Write();
  cSummary -> Close();
  cout << "    Made summary plot." << endl;

  // save graphs
  fOut -> cd();
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gStat[iTrg]    -> Write();
    gUnf[iTrg]     -> Write();
    gPyF[iTrg]     -> Write();
    gPyH[iTrg]     -> Write();
    gPlotSt[iTrg]  -> Write();
    gPlotUn[iTrg]  -> Write();
    gPlotLoU[iTrg] -> Write();
    gPlotHiU[iTrg] -> Write();
    gPlotPyF[iTrg] -> Write();
    gPlotPyH[iTrg] -> Write();
  }
  for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
    hPythiaP[iTrgP] -> Write();
  }
  for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
    hPythiaG[iTrgG] -> Write();
  }
  hFrame -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    fInStat[iTrg] -> cd();
    fInStat[iTrg] -> Close();
    fInUnf[iTrg]  -> cd();
    fInUnf[iTrg]  -> Close();
    fInPyF[iTrg]  -> cd();
    fInPyF[iTrg]  -> Close();
    fInPyH[iTrg]  -> cd();
    fInPyH[iTrg]  -> Close();
  }
  for (UInt_t iTrgP = 0; iTrgP < NTrgPi; iTrgP++) {
    fInHistP[iTrgP] -> cd();
    fInHistP[iTrgP] -> Close();
  }
  for (UInt_t iTrgG = 0; iTrgG < NTrgGa; iTrgG++) {
    fInHistG[iTrgG] -> cd();
    fInHistG[iTrgG] -> Close();
  }
  cout << "  Finished making summary plot!" << endl;

}

// End ------------------------------------------------------------------------
