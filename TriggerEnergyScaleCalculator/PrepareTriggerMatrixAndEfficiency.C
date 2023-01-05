// 'PrepareTriggerMatrixAndEfficiency.C'
// Derek Anderson
// 04.28.2022
//
// Use this to normalize and
// plot a TES/R matrix
//
// NOTE: TrgType controls
// trigger species
//   TrgType = 0: pi0 (apply
//     shape weights)
//   TrgType = 1: gamma (do
//     not apply shape weights)

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TBox.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NParEtCut(1);
static const UInt_t NTrgBins(3);
static const UInt_t NEffHist(3);
static const UInt_t NEffDiv(2);
static const UInt_t NEtCuts(NTrgBins + 1);
static const UInt_t NRange(2);
static const UInt_t NProj(11);
static const UInt_t NVtx(4);
static const UInt_t NPad(2);
static const UInt_t NTxt(2);

// smoothing constants
static const UInt_t NParSmooth(3);
static const UInt_t NSmoothMax(4);
static const UInt_t NSmoothPi0[NTrgBins] = {4, 3, 2};
static const UInt_t NSmoothGam[NTrgBins] = {3, 2, 1};

// trigger species
static const UInt_t TrgType(1);


void PrepareTriggerMatrixAndEfficiency() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning TES/R matrix preparation..." << endl;

  // input parameters
  //const TString sIn("./output/November2022/particleGun.forPaper_withNoWeightMatrix_withTspNoEt.et650x650x0100vz55tsp008pi0.d14m11y2022.root");
  //const TString sIn("./output/November2022/particleGun.forPaper_withNoWeightMatrix_withTspNoEt.et650x650x0100vz55tsp0206gam.d14m11y2022.root");
  const TString sIn("./particleGun.forPwgComments_multiGamPerEvtNoVtx.et650gam.d18m12y2022.root");
  const TString sInMat("resolution/hEtMatchVsDetNoWeight");
  const TString sInPar("particle/hEtTrgParAllNoWeight");
  const TString sInDet[NTrgBins]   = {"match/hEtDet911_NoWeight",   "match/hEtDet1115_NoWeight",   "match/hEtDet1520_NoWeight"};
  const TString sInSmear[NTrgBins] = {"match/hEtMatch911_NoWeight", "match/hEtMatch1115_NoWeight", "match/hEtMatch1520_NoWeight"};
  const TString sInEff[NEffDiv]    = {"match/hEtTrgDetMatchNoWeight", "particle/hEtTrgParAllNoWeight"};

  // output parameters
  const TString sOut("triggerMatrix.forPwgComments_multiGamPerEvtNoVtx.et650x650vz55tsp0206gam.d19m12y2022.root");
  const TString sValTxt("triggerMatrix_values.eTdetVsEtPar_withSaskiaChanges.et650x920vz55tsp0206gam.d26m10y2022.txt");
  const TString sErrTxt("triggerMatrix_errors.eTdetVsEtPar_withSaskiaChanges.et650x920vz55tsp0206gam.d26m10y2022.txt");
  const TString sWeiTxt("triggerMatrix_weights.eTdetVsEtPar_withSaskiaChanges.et650x920vz55tsp0206gam.d26m10y2022.txt");
  const TString sOutMat("hTriggerMatrix");
  const TString sOutPar("hEtPar");
  const TString sOutEff("hEtEff");
  const TString sOutDet[NTrgBins]     = {"hEtDet911",       "hEtDet1115",       "hEtDet1520"};
  const TString sOutBin[NTrgBins]     = {"hEtProj911",      "hEtProj1115",      "hEtProj1520"};
  const TString sOutSmear[NTrgBins]   = {"hEtSmear911",     "hEtSmear1115",     "hEtSmear1520"};
  const TString sOutParNorm[NTrgBins] = {"hEtParNorm911",   "hEtParNorm1115",   "hEtParNorm1520"};
  const TString sBias[NTrgBins]       = {"hTspBias911",     "hTspBias1115",     "hTspBias1520"};
  const TString sWeight[NTrgBins]     = {"hSmearWeight911", "hSmearWeight1115", "hSmearWeight1520"};
  const TString sEffDiv[NEffDiv]      = {"hEtEff_detector", "hEtEff_input"};

  // calculation parameters
  const UInt_t   nRebin(2);
  const Bool_t   rebinHist(false);
  const Bool_t   zeroOutBins(true);
  const Double_t xEtDetLo[NTrgBins]  = {9.,  11., 15.};
  const Double_t xEtDetHi[NTrgBins]  = {11., 15., 20.};
  const Double_t yEtParRange[NRange] = {6., 30.};

  // cut parameters
  const Double_t xEtCutRange[NRange] = {9., 20.};
  const UInt_t   fEtCutCol(1);
  const UInt_t   fEtCutFil(3345);
  const UInt_t   fEtCutLin(1);
  const UInt_t   fEtCutWid(2);

  // projection and weight parameters
  const TString  sLabPrefP("E_{T}^{#pi^{0},det-clust} = ");
  const TString  sLabPrefG("E_{T}^{#gamma,det-clust} = ");
  const TString  sLabSuff(".5 GeV");
  const TString  sProjPref("hProj");
  const TString  sProjPrefW("hProjForWeightEt");
  const TString  sWeightPref("hWeightEt");
  const TString  sWeightPar("hParForWeight");
  const TString  sWeightMat("hWeightMatrix");
  const TString  sProjSuff[NProj] = {"6", "8", "10", "12", "14", "16", "18", "20", "22", "24", "26"};
  const Double_t xProject[NProj]  = {6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5};
  const UInt_t   fColProj[NProj]  = {799, 809, 899,  909,  879,  889,  859,  869,  839,  849,  819};
  const UInt_t   fMarProj[NProj]  = {24,  26,  32,   25,   27,   28,   30,   24,   26,   32,   25};

  // smearing weight parameters
  const UInt_t  fColSmearPi0[NTrgBins] = {859, 899, 819};
  const UInt_t  fColSmearGam[NTrgBins] = {799, 899, 859}; 
  const UInt_t  fMarSmear[NTrgBins]    = {20,  22,  21};
  const TString sLegSmearP[NTrgBins]   = {"9 < E_{T}^{#pi^{0},det-clust} < 11 GeV", "11 < E_{T}^{#pi^{0},det-clust} < 15 GeV", "15 < E_{T}^{#pi^{0},det-clust} < 20 GeV"};
  const TString sLegSmearG[NTrgBins]   = {"9 < E_{T}^{#gamma,det-clust} < 11 GeV",  "11 < E_{T}^{#gamma,det-clust} < 15 GeV",  "15 < E_{T}^{#gamma,det-clust} < 20 GeV"};
  const TString sLegParSP("generated #pi^{0}");
  const TString sLegParSG("generated #gamma");

  // plot parameters
  const TString sTitleN("");
  const TString sTitleB("");
  const TString sTitleW("");
  const TString sTrg("E_{T} [GeV]");
  const TString sParP("E_{T}^{#pi^{0},part-match} [GeV]");
  const TString sParG("E_{T}^{#gamma,part-match} [GeV]");
  const TString sDetP("E_{T}^{#pi^{0},det-clust} [GeV]");
  const TString sDetG("E_{T}^{#gamma,det-clust} [GeV]");
  const TString sNum("counts");
  const TString sRat("TSP bias");
  const TString sEff("#epsilon_{trg} = reco. / thrown");
  const TString sWei("smearing weights");
  const Float_t fTtlE(0.083);
  const Float_t fLbl(0.03);
  const Float_t fLblE(0.063);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.2);
  const Float_t fOffZ(1.);
  const Float_t fOffXE(1.1);
  const Float_t fOffYE(0.58);
  const Float_t xyPlot[NVtx] = {6., 6., 30., 30.};
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fFil(0);
  const UInt_t  fMarS(29);
  const UInt_t  fColS(923);
  const UInt_t  fMarR[NTrgBins] = {24,  26,  25};
  const UInt_t  fMarW[NTrgBins] = {20,  22,  21};
  const UInt_t  fMarE[NEffHist] = {22,  23,  20};
  const UInt_t  fColE[NEffHist] = {899, 859, 923};

  // legend parameters
  const UInt_t  nDec(3);
  const TString sNormP("#LT #int T(E_{T}^{#pi^{0},part-match}, E_{T}^{#pi^{0},det-clust}) dE_{T}^{#pi^{0},part-match} #GT");
  const TString sNormG("#LT #int T(E_{T}^{#gamma,part-match}, E_{T}^{#gamma,det-clust}) dE_{T}^{#gamma,part-match} #GT");
  const TString sHeaderI("#bf{average integral}");
  const TString sHeaderT("select / true");
  const TString sLabelS("sim. particles");
  const TString sEnergy[NTrgBins] = {"9 - 11 GeV", "11 - 15 GeV", "15 - 20 GeV"};
  const TString sLabelRP[NTrgBins] = {"reco. particles, E_{T}^{#pi^{0},det-clust} #in (9, 11) GeV",
                                      "reco. particles, E_{T}^{#pi^{0},det-clust} #in (11, 15) GeV",
                                      "reco. particles, E_{T}^{#pi^{0},det-clust} #in (15, 20) GeV"};
  const TString sLabelRG[NTrgBins] = {"reco. particles, E_{T}^{#gamma,det-clust} #in (9, 11) GeV",
                                      "reco. particles, E_{T}^{#gamma,det-clust} #in (11, 15) GeV",
                                      "reco. particles, E_{T}^{#gamma,det-clust} #in (15, 20) GeV"};
  const TString sLabelTP[NTrgBins] = {"E_{T}^{#pi^{0},det-clust} #in (9, 11) GeV",
                                      "E_{T}^{#pi^{0},det-clust} #in (11, 15) GeV",
                                      "E_{T}^{#pi^{0},det-clust} #in (15, 20) GeV"};
  const TString sLabelTG[NTrgBins] = {"E_{T}^{#gamma,det-clust} #in (9, 11) GeV",
                                      "E_{T}^{#gamma,det-clust} #in (11, 15) GeV",
                                      "E_{T}^{#gamma,det-clust} #in (15, 20) GeV"};
  const TString sLabelE[NEffHist] = {"Reconstructed", "Thrown", "Efficiency, #epsilon_{trg}"};

  // trigger parameters
  const TString sTitlePi0("simulated #pi^{0}");
  const TString sTitleGam("simulated #gamma");
  const TString sLabelPi0("#bf{simulated #pi^{0}}");
  const TString sLabelGam("#bf{simulated #gamma}");
  const TString sWeightPi0("#pi^{0} weights");
  const TString sWeightGam("#gamma weights");
  const UInt_t  fColPi0[NTrgBins] = {859, 839, 819};
  const UInt_t  fColGam[NTrgBins] = {799, 899, 879};

  // text parameters
  const TString sTxtSimPi0("STAR simulation #pi^{0}, 6 < E_{T}^{#pi^{0},part,det-clust} < 50 GeV"); 
  const TString sTxtSimGam("STAR simulation #gamma, 6 < E_{T}^{#gamma,part,det-clust} < 50 GeV");
  const TString sTxtCutPi0("|#eta^{#pi^{0}}| < 0.9, TSP #in (0, 0.08)");
  const TString sTxtCutGam("|#eta^{#gamma}| < 0.9, TSP #in (0.2, 0.6)");

  // generic smoothing parameters
  const Bool_t  smoothWeights(false);
  const TString sSmoothFunc("pol2");
  const TString sSmoothFuncBase[NTrgBins]   = {"fSmooth911", "fSmooth1115", "fSmooth1520"};
  const Float_t pSmoothParGuess[NParSmooth] = {1., 1., 0.5};

  // pi0 smoothing parameters
  const Float_t rSmoothRangesLoP[NTrgBins][NSmoothMax]  = {{6.,  11.,   15.,   22.},  {7.,   12.,   18,  1.}, {6.,  20., 1., 1.}};
  const Float_t rSmoothRangesHiP[NTrgBins][NSmoothMax]  = {{8.5, 15.,   22.,   30.},  {11.,  16.,   30., 1.}, {14., 30., 1., 1.}};
  const Float_t rCorrectRangesLoP[NTrgBins][NSmoothMax] = {{6.,  11.,   14.5,  21.5}, {6.,   12.,   19., 1.}, {6.,  21., 1., 1.}};
  const Float_t rCorrectRangesHiP[NTrgBins][NSmoothMax] = {{8.5, 14.5., 21.5,  30.},  {10.5, 15.5., 30., 1.}, {14., 30., 1., 1.}};

  // gamma smoothing parameters
  const Float_t rSmoothRangesLoG[NTrgBins][NSmoothMax]  = {{6.,  15.,  20., 1.}, {6.,   23., 1., 1.}, {7.,   1., 1., 1.}};
  const Float_t rSmoothRangesHiG[NTrgBins][NSmoothMax]  = {{8.,  20.,  30., 1.}, {10.,  30., 1., 1.}, {14.,  1., 1., 1.}};
  const Float_t rCorrectRangesLoG[NTrgBins][NSmoothMax] = {{6.,  15.,  20., 1.}, {6.,   23., 1., 1.}, {7.,   1., 1., 1.}};
  const Float_t rCorrectRangesHiG[NTrgBins][NSmoothMax] = {{7.5, 19.5, 30., 1.}, {9.5., 30., 1., 1.}, {13.5, 1., 1., 1.}};

  // parse trigger selection
  UInt_t  fColR[NTrgBins];
  UInt_t  fColW[NTrgBins];
  UInt_t  fColSmear[NTrgBins];
  UInt_t  nSmooth[NTrgBins];
  Bool_t  applyBias(false);
  TString sPar("");
  TString sDet("");
  TString sNorm("");
  TString sLabel("");
  TString sLegPar("");
  TString sTitleM("");
  TString sTitleWM("");
  TString sTxtSim("");
  TString sTxtCut("");
  TString sLabPref("");
  TString sLabelR[NTrgBins];
  TString sLabelT[NTrgBins];
  TString sLegSmear[NTrgBins];
  Float_t rSmoothRangesLo[NTrgBins][NSmoothMax];
  Float_t rSmoothRangesHi[NTrgBins][NSmoothMax];
  Float_t rCorrectRangesLo[NTrgBins][NSmoothMax];
  Float_t rCorrectRangesHi[NTrgBins][NSmoothMax];
  switch (TrgType) {
    case 0:
      applyBias = false;
      sPar      = sParP;
      sDet      = sDetP;
      sNorm     = sNormP;
      sLabel    = sLabelPi0;
      sLegPar   = sLegParSP;
      sTitleM   = sTitlePi0;
      sTitleWM  = sWeightPi0;
      sTxtSim   = sTxtSimPi0;
      sTxtCut   = sTxtCutPi0;
      sLabPref  = sLabPrefP;
      for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
        fColR[iTrg]     = fColPi0[iTrg];
        fColSmear[iTrg] = fColSmearPi0[iTrg];
        sLabelR[iTrg]   = sLabelRP[iTrg];
        sLabelT[iTrg]   = sLabelTP[iTrg];
        sLegSmear[iTrg] = sLegSmearP[iTrg];
        nSmooth[iTrg]   = NSmoothPi0[iTrg];
        for (UInt_t iSmooth = 0; iSmooth < NSmoothMax; iSmooth++) {
          rSmoothRangesLo[iTrg][iSmooth]  = rSmoothRangesLoP[iTrg][iSmooth];
          rSmoothRangesHi[iTrg][iSmooth]  = rSmoothRangesHiP[iTrg][iSmooth];
          rCorrectRangesLo[iTrg][iSmooth] = rCorrectRangesLoP[iTrg][iSmooth];
          rCorrectRangesHi[iTrg][iSmooth] = rCorrectRangesHiP[iTrg][iSmooth];
        }
      }
      break;
    case 1:
      applyBias = false;
      sPar      = sParG;
      sDet      = sDetG;
      sNorm     = sNormG;
      sLabel    = sLabelGam;
      sLegPar   = sLegParSG;
      sTitleM   = sTitleGam;
      sTitleWM  = sWeightGam;
      sTxtSim   = sTxtSimGam;
      sTxtCut   = sTxtCutGam;
      sLabPref  = sLabPrefG;
      for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
        fColR[iTrg]     = fColGam[iTrg];
        fColSmear[iTrg] = fColSmearGam[iTrg];
        sLabelR[iTrg]   = sLabelRG[iTrg];
        sLabelT[iTrg]   = sLabelTG[iTrg];
        sLegSmear[iTrg] = sLegSmearG[iTrg];
        nSmooth[iTrg]   = NSmoothGam[iTrg];
        for (UInt_t iSmooth = 0; iSmooth < NSmoothMax; iSmooth++) {
          rSmoothRangesLo[iTrg][iSmooth]  = rSmoothRangesLoG[iTrg][iSmooth];
          rSmoothRangesHi[iTrg][iSmooth]  = rSmoothRangesHiG[iTrg][iSmooth];
          rCorrectRangesLo[iTrg][iSmooth] = rCorrectRangesLoG[iTrg][iSmooth];
          rCorrectRangesHi[iTrg][iSmooth] = rCorrectRangesHiG[iTrg][iSmooth];
        }
      }
      break;
    default:
      applyBias = false;
      sPar      = sParP;
      sDet      = sDetP;
      sNorm     = sNormP;
      sLabel    = sLabelPi0;
      sLegPar   = sLegParSP;
      sTitleM   = sTitlePi0;
      sTitleWM  = sWeightPi0;
      sTxtSim   = sTxtSimPi0;
      sTxtCut   = sTxtCutPi0;
      sLabPref  = sLabPrefP;
      for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
        fColR[iTrg]     = fColPi0[iTrg];
        fColSmear[iTrg] = fColSmearPi0[iTrg];
        sLabelR[iTrg]   = sLabelRP[iTrg];
        sLabelT[iTrg]   = sLabelTP[iTrg];
        sLegSmear[iTrg] = sLegSmearP[iTrg];
        nSmooth[iTrg]   = NSmoothPi0[iTrg];
        for (UInt_t iSmooth = 0; iSmooth < NSmoothMax; iSmooth++) {
          rSmoothRangesLo[iTrg][iSmooth]  = rSmoothRangesLoP[iTrg][iSmooth];
          rSmoothRangesHi[iTrg][iSmooth]  = rSmoothRangesHiP[iTrg][iSmooth];
          rCorrectRangesLo[iTrg][iSmooth] = rCorrectRangesLoP[iTrg][iSmooth];
          rCorrectRangesHi[iTrg][iSmooth] = rCorrectRangesHiP[iTrg][iSmooth];
        }
      }
      break;
  }

  // open files
  TFile *fIn  = new TFile(sIn.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn || !fOut) {
    cerr << "PANIC: couldn't open input or output file!\n"
         << "       fIn = " << fIn << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH2D *hMatrix   = (TH2D*) fIn -> Get(sInMat.Data());
  TH1D *hParticle = (TH1D*) fIn -> Get(sInPar.Data());
  if (!hMatrix || !hParticle) {
    cerr << "PANIC: couldn't grab input matrix or particle histogram!\n"
         << "       hMatrix = " << hMatrix << ", hParticle = " << hParticle << "\n"
         << endl;
    return;
  }
  hMatrix   -> SetName(sOutMat.Data());
  hParticle -> SetName(sOutPar.Data());

  TH1D *hEffDiv[NEffDiv];
  TH1D *hSmeared[NTrgBins];
  TH1D *hDetector[NTrgBins];
  for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
    hEffDiv[iDiv] = (TH1D*) fIn -> Get(sInEff[iDiv].Data());
    if (!hEffDiv[iDiv]) {
      cerr << "PANIC: couldn't grab input efficiency histogram!\n"
           << "       hEffDiv[" << iDiv << "] = " << hEffDiv[iDiv] << "\n"
           << endl;
    }
    hEffDiv[iDiv] -> SetName(sEffDiv[iDiv].Data());
  }
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hSmeared[iTrg]  = (TH1D*) fIn -> Get(sInSmear[iTrg].Data());
    hDetector[iTrg] = (TH1D*) fIn -> Get(sInDet[iTrg].Data());
    if (!hSmeared[iTrg] || !hDetector[iTrg]) {
      cerr << "PANIC: couldn't grab input smeared or detector histogram!\n"
           << "       hsmeared[" << iTrg << "] = " << hSmeared[iTrg] << ", hDetector[" << iTrg << "] = " << hDetector[iTrg] << "\n"
           << endl;
      return;
    }
    hSmeared[iTrg]  -> SetName(sOutSmear[iTrg].Data());
    hDetector[iTrg] -> SetName(sOutDet[iTrg].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  if (rebinHist) {
    for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
      hEffDiv[iDiv] -> Rebin(nRebin);
    }
    for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
      hDetector[iTrg] -> Rebin(nRebin);
    }
    hParticle -> Rebin(nRebin);
    cout << "    Rebinned histograms." << endl;
  }

  // zero out bins outside par range
  const UInt_t nBinsX = hMatrix   -> GetNbinsX();
  const UInt_t nBinsY = hMatrix   -> GetNbinsY();
  const UInt_t nBinsP = hParticle -> GetNbinsX();
  if (zeroOutBins) {

    // zero out matrix bins
    for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {
      for (UInt_t iBinY = 1; iBinY < (nBinsY + 1); iBinY++) {
        const Double_t yCenterM      = hMatrix -> GetYaxis() -> GetBinCenter(iBinY);
        const Bool_t   isInParRangeM = ((yCenterM > yEtParRange[0]) && (yCenterM < yEtParRange[1]));
        if (!isInParRangeM) {
          hMatrix -> SetBinContent(iBinX, iBinY, 0.);
          hMatrix -> SetBinError(iBinX, iBinY, 0.);
        }
      }  // end eTpar loop
    }  // end eTdet loop

    // zero out particle bins
    for (UInt_t iBinP = 1; iBinP < (nBinsP + 1); iBinP++) {
      const Double_t xCenterP      = hParticle -> GetBinCenter(iBinP);
      const Bool_t   isInParRangeP = ((xCenterP > yEtParRange[0]) && (xCenterP < yEtParRange[1]));
      if (!isInParRangeP) {
        hParticle -> SetBinContent(iBinP, 0.);
        hParticle -> SetBinError(iBinP, 0.);
      }
    }

    // zero out smear bins
    for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
      for (UInt_t iBinS = 1; iBinS < (nBinsP + 1); iBinS++) {
        const Double_t xCenterS      = hSmeared[iTrg] -> GetBinCenter(iBinS);
        const Bool_t   isInParRangeS = ((xCenterS > yEtParRange[0]) && (xCenterS < yEtParRange[1]));
        if (!isInParRangeS) {
          hSmeared[iTrg] -> SetBinContent(iBinS, 0.);
          hSmeared[iTrg] -> SetBinError(iBinS, 0.);
        }
      }  // end smear bin loop
    }  // end trigger bin loop
    cout << "    Zeroed out bins outside particle range." << endl;
  }

  // normalize 1d histograms
  const Double_t xBinHist = hParticle -> GetBinWidth(1);
  const Double_t xCheck   = xBinHist / 4.;
  for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
    hEffDiv[iDiv] -> Scale(1. / xBinHist);
  }
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hSmeared[iTrg]  -> Scale(1. / xBinHist);
    hDetector[iTrg] -> Scale(1. / xBinHist);
  }
  hParticle -> Scale(1. / xBinHist);

  const Double_t intPar  = hParticle  -> Integral();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    const UInt_t   iDetLo  = hParticle       -> FindBin(xEtDetLo[iTrg] + xCheck);
    const UInt_t   iDetHi  = hParticle       -> FindBin(xEtDetHi[iTrg] - xCheck);
    const Double_t intDet  = hDetector[iTrg] -> Integral();
    const Double_t intRel  = hParticle       -> Integral(iDetLo, iDetHi);
    const Double_t normDet = intRel / intDet;
    hDetector[iTrg] -> Scale(normDet);
  }

  TH1D *hParForNorm[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {

    // for weight normalization
    hParForNorm[iTrg] = (TH1D*) hParticle -> Clone();
    hParForNorm[iTrg] -> SetName(sOutParNorm[iTrg].Data());

    // TEST [11.29.2022]
    const Double_t intDetNorm = hDetector[iTrg] -> Integral();

    // do normalization
    const UInt_t   iParLo     = hParForNorm[iTrg] -> FindBin(xEtDetLo[iTrg] + xCheck);
    const UInt_t   iParHi     = hParForNorm[iTrg] -> FindBin(xEtDetHi[iTrg] - xCheck);
    const Double_t intParNorm = hParForNorm[iTrg] -> Integral(iParLo, iParHi);
    const Double_t intSmear   = hSmeared[iTrg]    -> Integral();
    const Double_t normPar    = intDetNorm / intParNorm;
    const Double_t normSmear  = intDetNorm / intSmear;
    hParForNorm[iTrg] -> Scale(normPar);
    hSmeared[iTrg]    -> Scale(normSmear);
  }
  cout << "    Normalized 1d histograms." << endl;

  // calculate efficiency
  TH1D *hEff = (TH1D*) hEffDiv[0] -> Clone();
  hEff -> SetName(sOutEff.Data());
  hEff -> Reset("ICES");
  hEff -> Divide(hEffDiv[0], hEffDiv[1], 1., 1.);
  cout << "    Calculated efficiency." << endl;

  // calculate TSP bias
  TH1D *hBias[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hBias[iTrg] = (TH1D*) hDetector[iTrg] -> Clone();
    hBias[iTrg] -> SetName(sBias[iTrg].Data());
    hBias[iTrg] -> Reset("ICES");
    hBias[iTrg] -> Divide(hDetector[iTrg], hParticle, 1., 1.);
  }
  cout << "    Calculated TSP bias." << endl;

  // calculate smearing weights
  TH1D *hSmearWeights[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hSmearWeights[iTrg] = (TH1D*) hParticle -> Clone();
    hSmearWeights[iTrg] -> SetName(sWeight[iTrg].Data());
    hSmearWeights[iTrg] -> Reset("ICES");
    hSmearWeights[iTrg] -> Divide(hSmeared[iTrg], hParForNorm[iTrg], 1., 1.);
  }
  cout << "    Calculated smear weights." << endl;

  // smooth weights (if need be)
  TF1  *fSmoothWeight[NTrgBins][NSmoothMax];
  if (smoothWeights) {
    for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
      for (UInt_t iSmooth = 0; iSmooth < nSmooth[iTrg]; iSmooth++) {

        // create name
        TString sSmoothFuncName(sSmoothFuncBase[iTrg].Data());
        sSmoothFuncName += "_";
        sSmoothFuncName += iSmooth;

        // initialize smoothing functions
        fSmoothWeight[iTrg][iSmooth] = new TF1(sSmoothFuncName.Data(), sSmoothFunc.Data(), yEtParRange[0], yEtParRange[1]);
        for (UInt_t iSmoothPar = 0; iSmoothPar < NParSmooth; iSmoothPar++) {
          fSmoothWeight[iTrg][iSmooth] -> SetParameter(iSmoothPar, pSmoothParGuess[iSmoothPar]);
        }

        // fit weights
        Int_t fitStatus(-1);
        if (iSmooth == 0) {
          fitStatus = (Int_t) hSmearWeights[iTrg] -> Fit(sSmoothFuncName.Data(), "", "", rSmoothRangesLo[iTrg][iSmooth], rSmoothRangesHi[iTrg][iSmooth]);
        } else {
          fitStatus = (Int_t) hSmearWeights[iTrg] -> Fit(sSmoothFuncName.Data(), "+", "", rSmoothRangesLo[iTrg][iSmooth], rSmoothRangesHi[iTrg][iSmooth]);
        }
        const Bool_t isGoodFit = (fitStatus == 0);

        // apply smoothing
        for (UInt_t iBinS = 1; iBinS < (nBinsP + 1); iBinS++) {

          // get bin info
          const Double_t xCenterS = hSmearWeights[iTrg]          -> GetBinCenter(iBinS);
          const Double_t yRawS    = hSmearWeights[iTrg]          -> GetBinContent(iBinS);
          const Double_t ySmooth  = fSmoothWeight[iTrg][iSmooth] -> Eval(xCenterS);
          const Double_t eRawS    = hSmearWeights[iTrg]          -> GetBinError(iBinS);
          const Double_t ePerRawS = eRawS / yRawS;
          const Double_t eSmooth  = ePerRawS * ySmooth;

          // if in range and fit is good, apply smoothing
          const Bool_t isBinNonzero     = ((yRawS > 0.) && (ySmooth > 0.));
          const Bool_t isInCorrectRange = ((xCenterS > rCorrectRangesLo[iTrg][iSmooth]) && (xCenterS < rCorrectRangesHi[iTrg][iSmooth]));
          if (isBinNonzero && isGoodFit && isInCorrectRange) {
            hSmearWeights[iTrg] -> SetBinContent(iBinS, ySmooth);
            hSmearWeights[iTrg] -> SetBinError(iBinS, eSmooth);
          }
        }  // end smeared bin loop
      }  // end smooth loop
    }  // end trigger bin loop
    cout << "    Smoothed smear weights." << endl;
  }

  // remove any degree of non-flatness from matrix
  const UInt_t iEtCutStart = hMatrix -> GetXaxis() -> FindBin(xEtCutRange[0] - xCheck);
  const UInt_t iEtCutStop  = hMatrix -> GetXaxis() -> FindBin(xEtCutRange[1] + xCheck);
  for (UInt_t iBinY = 1; iBinY < (nBinsY + 1); iBinY++) {

    // get scale factor for eTpar slice
    const Double_t totParIntegral   = hMatrix -> Integral(1, nBinsX, iBinY, iBinY);
    const Double_t outParIntegralLo = hMatrix -> Integral(1, iEtCutStart, iBinY, iBinY);
    const Double_t outParIntegralHi = hMatrix -> Integral(iEtCutStop, 1, iBinY, iBinY);
    const Double_t flatScaleFactor  = (outParIntegralLo + outParIntegralHi) / totParIntegral;

    // skip zero slices
    const Bool_t isParIntNonzero = (totParIntegral > 0.);
    if (!isParIntNonzero) continue;

    // scale bins
    for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {
      const Double_t binVal       = hMatrix -> GetBinContent(iBinX, iBinY);
      const Double_t binErr       = hMatrix -> GetBinError(iBinX, iBinY);
      const Double_t scaleVal     = binVal / flatScaleFactor;
      const Double_t scaleErr     = binErr / flatScaleFactor;
      const Bool_t   isBinNonzero = (binVal > 0.);
      if (isBinNonzero) {
        hMatrix -> SetBinContent(iBinX, iBinY, scaleVal);
        hMatrix -> SetBinError(iBinX, iBinY, scaleErr);
      } else {
        hMatrix -> SetBinContent(iBinX, iBinY, 0.);
        hMatrix -> SetBinError(iBinX, iBinY, 0.);
      }
    }  // end 1st x-axis loop
  }  // end 1st y-axis loop

  // normalize eTdet slices of matrix
  for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {

    // get eTdet slice integral
    const Double_t detIntegral = hMatrix -> Integral(iBinX, iBinX, 1, nBinsY);
    const Double_t sliceNorm   = 1. / detIntegral;
    
    // skip zero slices
    const Bool_t isDetIntNonzero = (detIntegral > 0.);
    if (!isDetIntNonzero) continue;

    // scale bins
    for (UInt_t iBinY = 1; iBinY < (nBinsY + 1); iBinY++) {
      const Double_t binVal       = hMatrix -> GetBinContent(iBinX, iBinY);
      const Double_t binErr       = hMatrix -> GetBinError(iBinX, iBinY);
      const Double_t normVal      = binVal * sliceNorm;
      const Double_t normErr      = binErr * sliceNorm;
      const Bool_t   isBinNonzero = (binVal > 0.);
      if (isBinNonzero) {
        hMatrix -> SetBinContent(iBinX, iBinY, normVal);
        hMatrix -> SetBinError(iBinX, iBinY, normErr);
      } else {
        hMatrix -> SetBinContent(iBinX, iBinY, 0.);
        hMatrix -> SetBinError(iBinX, iBinY, 0.);
      }
    }  // end 2nd y-axis loop
  }  // end 2nd x-axis loop

  // check normalization
  UInt_t   nBinAvg[NTrgBins] = {0, 0, 0};
  Double_t avgNorm[NTrgBins] = {0., 0., 0.};
  for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {
    for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
      const Double_t binCenter    = hMatrix -> GetXaxis() -> GetBinCenter(iBinX);
      const Double_t binIntegral  = hMatrix -> Integral(iBinX, iBinX, 1, nBinsY);
      const Bool_t   isInDetRange = ((binCenter > xEtDetLo[iTrg]) && (binCenter < xEtDetHi[iTrg]));
      if (isInDetRange) {
        avgNorm[iTrg] += binIntegral;
        nBinAvg[iTrg]++;
      }
    }  // end trigger loop
  }  // end x-axis loop

  cout << "    Normalized matrix:" << endl;
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    avgNorm[iTrg] = avgNorm[iTrg] / nBinAvg[iTrg];
    cout << "      [" << xEtDetLo[iTrg] << ", " << xEtDetHi[iTrg] << "]: average integral = " << avgNorm[iTrg] << endl;
  }

  // get weight histograms
  TH1D *hWeightPar = (TH1D*) hParticle -> Clone();
  hWeightPar -> SetName(sWeightPar.Data());

  const Double_t intParW  = hWeightPar -> Integral();
  const Double_t normParW = 1. / intParW;
  hWeightPar -> Scale(normParW);

  TH1D *hProjW[NProj];
  TH1D *hWeight[NProj];
  for (UInt_t iProj = 0; iProj < NProj; iProj++) {

    // create names
    TString sNameProjW(sProjPrefW.Data());
    TString sNameWeight(sWeightPref.Data());
    sNameProjW.Append(sProjSuff[iProj].Data());
    sNameWeight.Append(sProjSuff[iProj].Data());

    // do projection
    const UInt_t iProject = hMatrix -> GetXaxis() -> FindBin(xProject[iProj]);
    hProjW[iProj] = (TH1D*) hMatrix -> ProjectionY("", iProject, iProject, "") -> Clone();
    hProjW[iProj] -> SetName(sNameProjW.Data());

    // calculate weights
    hWeight[iProj] = (TH1D*) hProjW[iProj] -> Clone();
    hWeight[iProj] -> SetName(sNameWeight.Data());
    hWeight[iProj] -> Reset("ICES");
    hWeight[iProj] -> Divide(hProjW[iProj], hWeightPar, 1., 1.);
  }
  cout << "    Calculated weights." << endl;

  // create weight matrix
  TH1D *hMatrixW = (TH1D*) hMatrix -> Clone();
  hMatrixW -> SetName(sWeightMat.Data());
  hMatrixW -> Reset("ICES");
  for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {

    // determine projection
    UInt_t iProjUse(0);
    UInt_t iProjBin(0);
    Bool_t isProject(false);
    for (UInt_t iProj = 0; iProj < NProj; iProj++) {
      iProjBin = hMatrixW -> GetXaxis() -> FindBin(xProject[iProj]);
      if (iProjBin == iBinX) {
        iProjUse  = iProj;
        isProject = true;
        break;
      }
    }

    // set contents of matrix
    if (isProject) {
      for (UInt_t iBinY = 1; iBinY < (nBinsY + 1); iBinY++) {
        const Double_t yValW = hWeight[iProjUse] -> GetBinContent(iBinY);
        const Double_t yErrW = hWeight[iProjUse] -> GetBinError(iBinY);
        hMatrixW -> SetBinContent(iBinX, iBinY, yValW);
        hMatrixW -> SetBinError(iBinX, iBinY, yErrW);
      }  // end y-axis loop
    }
  }  // end x-axis loop
  cout << "    Created weight matrix." << endl;

/*
  // read matrix contents to a matrix
  ofstream oValTxt(sValTxt.Data());
  ofstream oErrTxt(sErrTxt.Data());
  ofstream oWeiTxt(sWeiTxt.Data());
  for (UInt_t iBinY = 0; iBinY < (nBinsY + 1); iBinY++) {
    for (UInt_t iBinX = 0; iBinX < (nBinsX + 1); iBinX++) {

      // determine location
      const Bool_t isBeforeBinsY = (iBinY == 0);
      const Bool_t isBeforeBinsX = (iBinX == 0);
      const Bool_t isAtLastBinX  = (iBinX == nBinsX);
      const Bool_t isInBinsY     = (iBinY > 0);
      const Bool_t isInBinsX     = (iBinX > 0);

      // put nothing in top left corner
      if (isBeforeBinsY && isBeforeBinsX) {
        oValTxt << " ";
        oErrTxt << " ";
        oWeiTxt << " ";
      }

      // get bin center and content
      const Double_t xCenter = hMatrix  -> GetXaxis() -> GetBinCenter(iBinX);
      const Double_t yCenter = hMatrix  -> GetYaxis() -> GetBinCenter(iBinY);
      const Double_t content = hMatrix  -> GetBinContent(iBinX, iBinY);
      const Double_t weight  = hMatrixW -> GetBinContent(iBinX, iBinY);
      const Double_t error   = hMatrix  -> GetBinError(iBinX, iBinY);

      // print axes
      if (isBeforeBinsY && isInBinsX) {
        oValTxt << xCenter;
        oValTxt << " ";
        oErrTxt << xCenter;
        oErrTxt << " ";
        oWeiTxt << xCenter;
        oWeiTxt << " ";
      }
      if (isBeforeBinsX && isInBinsY) {
        oValTxt << yCenter;
        oValTxt << " ";
        oErrTxt << yCenter;
        oErrTxt << " ";
        oWeiTxt << yCenter;
        oWeiTxt << " ";
      }

      // print content and error
      if (isInBinsY && isInBinsX) {
        oValTxt << content;
        oValTxt << " ";
        oErrTxt << error;
        oErrTxt << " ";
        oWeiTxt << weight;
        oWeiTxt << " ";
      }

      // print out endlines
      if (isAtLastBinX) {
        oValTxt << endl;
        oErrTxt << endl;
        oWeiTxt << endl;
      }
    }  // end x-axis loop
  }  // end y-axis loop
  cout << "    Printed out matrix." << endl;
*/

  // create projections
  TH1D *hProj[NProj];
  for (UInt_t iProj = 0; iProj < NProj; iProj++) {

    // create name
    TString sNameProj(sProjPref.Data());
    sNameProj.Append(sProjSuff[iProj].Data());

    // do projection
    const UInt_t iProject = hMatrix -> GetXaxis() -> FindBin(xProject[iProj]);
    hProj[iProj] = (TH1D*) hMatrix -> ProjectionY("", iProject, iProject, "") -> Clone();
    hProj[iProj] -> SetName(sNameProj.Data());
  }

  Double_t iEtBin[NTrgBins][NRange] = {{0., 0.}, {0., 0.}, {0., 0.}};
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {
      const Double_t xCenter      = hMatrix -> GetXaxis() -> GetBinCenter(iBinX);
      const Bool_t   isInDetRange = ((xCenter > xEtDetLo[iTrg]) && (xCenter < xEtDetHi[iTrg]));
      const Bool_t   isFirstStart = (iEtBin[iTrg][0] == 0.);
      const Bool_t   isFirstStop  = (iEtBin[iTrg][1] == 0.);
      if (isInDetRange) {

        // set start point
        if (isFirstStart) {
          iEtBin[iTrg][0] = iBinX;
        }

        // set stop point
        if (isFirstStop) {
          iEtBin[iTrg][1] = iBinX;
        } else {
          iEtBin[iTrg][1]++;
        }
      }
    }  // end eTdet loop
  }  // end trigger loop

  TH1D *hEtProj[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hEtProj[iTrg] = (TH1D*) hMatrix -> ProjectionY("", iEtBin[iTrg][0], iEtBin[iTrg][1], "") -> Clone();
    hEtProj[iTrg] -> SetName(sOutBin[iTrg].Data());
  }
  cout << "    Projected distributions." << endl;

  // get matrix minima and maxima
  const Double_t minMatrix = hMatrix  -> GetMinimum(0.);
  const Double_t minWeight = hMatrixW -> GetMinimum(0.);
  const Double_t maxMatrix = hMatrix  -> GetMaximum();
  const Double_t maxWeight = hMatrixW -> GetMaximum();

  // set styles
  hEff       -> SetMarkerColor(fColE[NEffDiv]);
  hEff       -> SetMarkerStyle(fMarE[NEffDiv]);
  hEff       -> SetLineColor(fColE[NEffDiv]);
  hEff       -> SetLineStyle(fLin);
  hEff       -> SetLineWidth(fWid);
  hEff       -> SetFillColor(fColE[NEffDiv]);
  hEff       -> SetFillStyle(fFil);
  hEff       -> SetTitle(sTitleN.Data());
  hEff       -> SetTitleFont(fTxt);
  hEff       -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
  hEff       -> GetXaxis() -> SetLabelFont(fTxt);
  hEff       -> GetXaxis() -> SetLabelSize(fLblE);
  hEff       -> GetXaxis() -> SetTitle(sTrg.Data());
  hEff       -> GetXaxis() -> SetTitleSize(fTtlE);
  hEff       -> GetXaxis() -> SetTitleFont(fTxt);
  hEff       -> GetXaxis() -> SetTitleOffset(fOffXE);
  hEff       -> GetXaxis() -> CenterTitle(fCnt);
  hEff       -> GetYaxis() -> SetLabelFont(fTxt);
  hEff       -> GetYaxis() -> SetLabelSize(fLblE);
  hEff       -> GetYaxis() -> SetTitle(sEff.Data());
  hEff       -> GetYaxis() -> SetTitleSize(fTtlE);
  hEff       -> GetYaxis() -> SetTitleFont(fTxt);
  hEff       -> GetYaxis() -> SetTitleOffset(fOffYE);
  hEff       -> GetYaxis() -> CenterTitle(fCnt);
  hMatrix    -> SetTitle(sTitleM.Data());
  hMatrix    -> SetTitleFont(fTxt);
  hMatrix    -> GetXaxis() -> SetRangeUser(xyPlot[0], xyPlot[2]);
  hMatrix    -> GetXaxis() -> SetLabelFont(fTxt);
  hMatrix    -> GetXaxis() -> SetLabelSize(fLbl);
  hMatrix    -> GetXaxis() -> SetTitle(sDet.Data());
  hMatrix    -> GetXaxis() -> SetTitleFont(fTxt);
  hMatrix    -> GetXaxis() -> SetTitleOffset(fOffX);
  hMatrix    -> GetXaxis() -> CenterTitle(fCnt);
  hMatrix    -> GetYaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
  hMatrix    -> GetYaxis() -> SetLabelFont(fTxt);
  hMatrix    -> GetYaxis() -> SetLabelSize(fLbl);
  hMatrix    -> GetYaxis() -> SetTitle(sPar.Data());
  hMatrix    -> GetYaxis() -> SetTitleFont(fTxt);
  hMatrix    -> GetYaxis() -> SetTitleOffset(fOffY);
  hMatrix    -> GetYaxis() -> CenterTitle(fCnt);
  hMatrix    -> GetZaxis() -> SetRangeUser(minMatrix, maxMatrix);
  hMatrix    -> GetZaxis() -> SetLabelFont(fTxt);
  hMatrix    -> GetZaxis() -> SetLabelSize(fLbl);
  hMatrix    -> GetZaxis() -> SetTitle(sNum.Data());
  hMatrix    -> GetZaxis() -> SetTitleFont(fTxt);
  hMatrix    -> GetZaxis() -> SetTitleOffset(fOffZ);
  hMatrix    -> GetZaxis() -> CenterTitle(fCnt);
  hMatrixW   -> SetTitle(sTitleWM.Data());
  hMatrixW   -> SetTitleFont(fTxt);
  hMatrixW   -> GetXaxis() -> SetRangeUser(xyPlot[0], xyPlot[2]);
  hMatrixW   -> GetXaxis() -> SetLabelFont(fTxt);
  hMatrixW   -> GetXaxis() -> SetLabelSize(fLbl);
  hMatrixW   -> GetXaxis() -> SetTitle(sDet.Data());
  hMatrixW   -> GetXaxis() -> SetTitleFont(fTxt);
  hMatrixW   -> GetXaxis() -> SetTitleOffset(fOffX);
  hMatrixW   -> GetXaxis() -> CenterTitle(fCnt);
  hMatrixW   -> GetYaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
  hMatrixW   -> GetYaxis() -> SetLabelFont(fTxt);
  hMatrixW   -> GetYaxis() -> SetLabelSize(fLbl);
  hMatrixW   -> GetYaxis() -> SetTitle(sPar.Data());
  hMatrixW   -> GetYaxis() -> SetTitleFont(fTxt);
  hMatrixW   -> GetYaxis() -> SetTitleOffset(fOffY);
  hMatrixW   -> GetYaxis() -> CenterTitle(fCnt);
  hMatrixW   -> GetZaxis() -> SetLabelFont(fTxt);
  hMatrixW   -> GetZaxis() -> SetLabelSize(fLbl);
  hMatrixW   -> GetZaxis() -> SetRangeUser(minWeight, maxWeight);
  hParticle  -> SetMarkerColor(fColS);
  hParticle  -> SetMarkerStyle(fMarS);
  hParticle  -> SetLineColor(fColS);
  hParticle  -> SetLineStyle(fLin);
  hParticle  -> SetLineWidth(fWid);
  hParticle  -> SetFillColor(fColS);
  hParticle  -> SetFillStyle(fFil);
  hParticle  -> SetTitle(sTitleN.Data());
  hParticle  -> SetTitleFont(fTxt);
  hParticle  -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
  hParticle  -> GetXaxis() -> SetLabelFont(fTxt);
  hParticle  -> GetXaxis() -> SetLabelSize(fLbl);
  hParticle  -> GetXaxis() -> SetTitle(sTrg.Data());
  hParticle  -> GetXaxis() -> SetTitleFont(fTxt);
  hParticle  -> GetXaxis() -> SetTitleOffset(fOffX);
  hParticle  -> GetXaxis() -> CenterTitle(fCnt);
  hParticle  -> GetYaxis() -> SetLabelFont(fTxt);
  hParticle  -> GetYaxis() -> SetLabelSize(fLbl);
  hParticle  -> GetYaxis() -> SetTitle(sNum.Data());
  hParticle  -> GetYaxis() -> SetTitleFont(fTxt);
  hParticle  -> GetYaxis() -> SetTitleOffset(fOffY);
  hParticle  -> GetYaxis() -> CenterTitle(fCnt);
  hWeightPar -> SetMarkerColor(fColS);
  hWeightPar -> SetMarkerStyle(fMarS);
  hWeightPar -> SetLineColor(fColS);
  hWeightPar -> SetLineStyle(fLin);
  hWeightPar -> SetLineWidth(fWid);
  hWeightPar -> SetFillColor(fColS);
  hWeightPar -> SetFillStyle(fFil);
  hWeightPar -> SetTitle(sTitleW.Data());
  hWeightPar -> SetTitleFont(fTxt);
  hWeightPar -> GetXaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
  hWeightPar -> GetXaxis() -> SetLabelFont(fTxt);
  hWeightPar -> GetXaxis() -> SetLabelSize(fLbl);
  hWeightPar -> GetXaxis() -> SetTitle(sTrg.Data());
  hWeightPar -> GetXaxis() -> SetTitleFont(fTxt);
  hWeightPar -> GetXaxis() -> SetTitleOffset(fOffX);
  hWeightPar -> GetXaxis() -> CenterTitle(fCnt);
  hWeightPar -> GetYaxis() -> SetLabelFont(fTxt);
  hWeightPar -> GetYaxis() -> SetLabelSize(fLbl);
  hWeightPar -> GetYaxis() -> SetTitle(sNum.Data());
  hWeightPar -> GetYaxis() -> SetTitleFont(fTxt);
  hWeightPar -> GetYaxis() -> SetTitleOffset(fOffY);
  hWeightPar -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
    hEffDiv[iDiv] -> SetMarkerColor(fColE[iDiv]);
    hEffDiv[iDiv] -> SetMarkerStyle(fMarE[iDiv]);
    hEffDiv[iDiv] -> SetLineColor(fColE[iDiv]);
    hEffDiv[iDiv] -> SetLineStyle(fLin);
    hEffDiv[iDiv] -> SetLineWidth(fWid);
    hEffDiv[iDiv] -> SetFillColor(fColE[iDiv]);
    hEffDiv[iDiv] -> SetFillStyle(fFil);
    hEffDiv[iDiv] -> SetTitle(sTitleN.Data());
    hEffDiv[iDiv] -> SetTitleFont(fTxt);
    hEffDiv[iDiv] -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
    hEffDiv[iDiv] -> GetXaxis() -> SetLabelFont(fTxt);
    hEffDiv[iDiv] -> GetXaxis() -> SetLabelSize(fLbl);
    hEffDiv[iDiv] -> GetXaxis() -> SetTitle(sTrg.Data());
    hEffDiv[iDiv] -> GetXaxis() -> SetTitleFont(fTxt);
    hEffDiv[iDiv] -> GetXaxis() -> SetTitleOffset(fOffX);
    hEffDiv[iDiv] -> GetXaxis() -> CenterTitle(fCnt);
    hEffDiv[iDiv] -> GetYaxis() -> SetLabelFont(fTxt);
    hEffDiv[iDiv] -> GetYaxis() -> SetLabelSize(fLbl);
    hEffDiv[iDiv] -> GetYaxis() -> SetTitle(sNum.Data());
    hEffDiv[iDiv] -> GetYaxis() -> SetTitleFont(fTxt);
    hEffDiv[iDiv] -> GetYaxis() -> SetTitleOffset(fOffY);
    hEffDiv[iDiv] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hSmeared[iTrg]      -> SetMarkerColor(fColSmear[iTrg]);
    hSmeared[iTrg]      -> SetMarkerStyle(fMarSmear[iTrg]);
    hSmeared[iTrg]      -> SetLineColor(fColSmear[iTrg]);
    hSmeared[iTrg]      -> SetLineStyle(fLin);
    hSmeared[iTrg]      -> SetLineWidth(fWid);
    hSmeared[iTrg]      -> SetFillColor(fColSmear[iTrg]);
    hSmeared[iTrg]      -> SetFillStyle(fFil);
    hSmeared[iTrg]      -> SetTitle(sTitleN.Data());
    hSmeared[iTrg]      -> SetTitleFont(fTxt);
    hSmeared[iTrg]      -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
    hSmeared[iTrg]      -> GetXaxis() -> SetLabelFont(fTxt);
    hSmeared[iTrg]      -> GetXaxis() -> SetLabelSize(fLbl);
    hSmeared[iTrg]      -> GetXaxis() -> SetTitle(sPar.Data());
    hSmeared[iTrg]      -> GetXaxis() -> SetTitleFont(fTxt);
    hSmeared[iTrg]      -> GetXaxis() -> SetTitleOffset(fOffX);
    hSmeared[iTrg]      -> GetXaxis() -> CenterTitle(fCnt);
    hSmeared[iTrg]      -> GetYaxis() -> SetLabelFont(fTxt);
    hSmeared[iTrg]      -> GetYaxis() -> SetLabelSize(fLbl);
    hSmeared[iTrg]      -> GetYaxis() -> SetTitle(sNum.Data());
    hSmeared[iTrg]      -> GetYaxis() -> SetTitleFont(fTxt);
    hSmeared[iTrg]      -> GetYaxis() -> SetTitleOffset(fOffY);
    hSmeared[iTrg]      -> GetYaxis() -> CenterTitle(fCnt);
    hDetector[iTrg]     -> SetMarkerColor(fColR[iTrg]);
    hDetector[iTrg]     -> SetMarkerStyle(fMarR[iTrg]);
    hDetector[iTrg]     -> SetLineColor(fColR[iTrg]);
    hDetector[iTrg]     -> SetLineStyle(fLin);
    hDetector[iTrg]     -> SetLineWidth(fWid);
    hDetector[iTrg]     -> SetFillColor(fColR[iTrg]);
    hDetector[iTrg]     -> SetFillStyle(fFil);
    hDetector[iTrg]     -> SetTitle(sTitleN.Data());
    hDetector[iTrg]     -> SetTitleFont(fTxt);
    hDetector[iTrg]     -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
    hDetector[iTrg]     -> GetXaxis() -> SetLabelFont(fTxt);
    hDetector[iTrg]     -> GetXaxis() -> SetLabelSize(fLbl);
    hDetector[iTrg]     -> GetXaxis() -> SetTitle(sTrg.Data());
    hDetector[iTrg]     -> GetXaxis() -> SetTitleFont(fTxt);
    hDetector[iTrg]     -> GetXaxis() -> SetTitleOffset(fOffX);
    hDetector[iTrg]     -> GetXaxis() -> CenterTitle(fCnt);
    hDetector[iTrg]     -> GetYaxis() -> SetLabelFont(fTxt);
    hDetector[iTrg]     -> GetYaxis() -> SetLabelSize(fLbl);
    hDetector[iTrg]     -> GetYaxis() -> SetTitle(sNum.Data());
    hDetector[iTrg]     -> GetYaxis() -> SetTitleFont(fTxt);
    hDetector[iTrg]     -> GetYaxis() -> SetTitleOffset(fOffY);
    hDetector[iTrg]     -> GetYaxis() -> CenterTitle(fCnt);
    hBias[iTrg]         -> SetMarkerColor(fColR[iTrg]);
    hBias[iTrg]         -> SetMarkerStyle(fMarW[iTrg]);
    hBias[iTrg]         -> SetLineColor(fColR[iTrg]);
    hBias[iTrg]         -> SetLineStyle(fLin);
    hBias[iTrg]         -> SetLineWidth(fWid);
    hBias[iTrg]         -> SetFillColor(fColR[iTrg]);
    hBias[iTrg]         -> SetFillStyle(fFil);
    hBias[iTrg]         -> SetTitle(sTitleB.Data());
    hBias[iTrg]         -> SetTitleFont(fTxt);
    hBias[iTrg]         -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
    hBias[iTrg]         -> GetXaxis() -> SetLabelFont(fTxt);
    hBias[iTrg]         -> GetXaxis() -> SetLabelSize(fLbl);
    hBias[iTrg]         -> GetXaxis() -> SetTitle(sTrg.Data());
    hBias[iTrg]         -> GetXaxis() -> SetTitleFont(fTxt);
    hBias[iTrg]         -> GetXaxis() -> SetTitleOffset(fOffX);
    hBias[iTrg]         -> GetXaxis() -> CenterTitle(fCnt);
    hBias[iTrg]         -> GetYaxis() -> SetLabelFont(fTxt);
    hBias[iTrg]         -> GetYaxis() -> SetLabelSize(fLbl);
    hBias[iTrg]         -> GetYaxis() -> SetTitle(sRat.Data());
    hBias[iTrg]         -> GetYaxis() -> SetTitleFont(fTxt);
    hBias[iTrg]         -> GetYaxis() -> SetTitleOffset(fOffY);
    hBias[iTrg]         -> GetYaxis() -> CenterTitle(fCnt);
    hSmearWeights[iTrg] -> SetMarkerColor(fColSmear[iTrg]);
    hSmearWeights[iTrg] -> SetMarkerStyle(fMarSmear[iTrg]);
    hSmearWeights[iTrg] -> SetLineColor(fColSmear[iTrg]);
    hSmearWeights[iTrg] -> SetLineStyle(fLin);
    hSmearWeights[iTrg] -> SetLineWidth(fWid);
    hSmearWeights[iTrg] -> SetFillColor(fColSmear[iTrg]);
    hSmearWeights[iTrg] -> SetFillStyle(fFil);
    hSmearWeights[iTrg] -> SetTitle(sTitleN.Data());
    hSmearWeights[iTrg] -> SetTitleFont(fTxt);
    hSmearWeights[iTrg] -> GetXaxis() -> SetRangeUser(yEtParRange[0], yEtParRange[1]);
    hSmearWeights[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hSmearWeights[iTrg] -> GetXaxis() -> SetLabelSize(fLbl);
    hSmearWeights[iTrg] -> GetXaxis() -> SetTitle(sPar.Data());
    hSmearWeights[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hSmearWeights[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hSmearWeights[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hSmearWeights[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hSmearWeights[iTrg] -> GetYaxis() -> SetLabelSize(fLbl);
    hSmearWeights[iTrg] -> GetYaxis() -> SetTitle(sWei.Data());
    hSmearWeights[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hSmearWeights[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hSmearWeights[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hEtProj[iTrg]       -> SetMarkerColor(fColR[iTrg]);
    hEtProj[iTrg]       -> SetMarkerStyle(fMarW[iTrg]);
    hEtProj[iTrg]       -> SetLineColor(fColR[iTrg]);
    hEtProj[iTrg]       -> SetLineStyle(fLin);
    hEtProj[iTrg]       -> SetLineWidth(fWid);
    hEtProj[iTrg]       -> SetFillColor(fColR[iTrg]);
    hEtProj[iTrg]       -> SetFillStyle(fFil);
    hEtProj[iTrg]       -> SetTitle(sTitleN.Data());
    hEtProj[iTrg]       -> SetTitleFont(fTxt);
    hEtProj[iTrg]       -> GetXaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
    hEtProj[iTrg]       -> GetXaxis() -> SetLabelFont(fTxt);
    hEtProj[iTrg]       -> GetXaxis() -> SetLabelSize(fLbl);
    hEtProj[iTrg]       -> GetXaxis() -> SetTitle(sDet.Data());
    hEtProj[iTrg]       -> GetXaxis() -> SetTitleFont(fTxt);
    hEtProj[iTrg]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hEtProj[iTrg]       -> GetXaxis() -> CenterTitle(fCnt);
    hEtProj[iTrg]       -> GetYaxis() -> SetLabelFont(fTxt);
    hEtProj[iTrg]       -> GetYaxis() -> SetLabelSize(fLbl);
    hEtProj[iTrg]       -> GetYaxis() -> SetTitle(sNum.Data());
    hEtProj[iTrg]       -> GetYaxis() -> SetTitleFont(fTxt);
    hEtProj[iTrg]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hEtProj[iTrg]       -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iProj = 0; iProj < NProj; iProj++) {
    hProj[iProj]   -> SetMarkerColor(fColProj[iProj]);
    hProj[iProj]   -> SetMarkerStyle(fMarProj[iProj]);
    hProj[iProj]   -> SetLineColor(fColProj[iProj]);
    hProj[iProj]   -> SetLineStyle(fLin);
    hProj[iProj]   -> SetLineWidth(fWid);
    hProj[iProj]   -> SetFillColor(fColProj[iProj]);
    hProj[iProj]   -> SetFillStyle(fFil);
    hProj[iProj]   -> SetTitle(sTitleN.Data());
    hProj[iProj]   -> SetTitleFont(fTxt);
    hProj[iProj]   -> GetXaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
    hProj[iProj]   -> GetXaxis() -> SetLabelFont(fTxt);
    hProj[iProj]   -> GetXaxis() -> SetLabelSize(fLbl);
    hProj[iProj]   -> GetXaxis() -> SetTitle(sPar.Data());
    hProj[iProj]   -> GetXaxis() -> SetTitleFont(fTxt);
    hProj[iProj]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hProj[iProj]   -> GetXaxis() -> CenterTitle(fCnt);
    hProj[iProj]   -> GetYaxis() -> SetLabelFont(fTxt);
    hProj[iProj]   -> GetYaxis() -> SetLabelSize(fLbl);
    hProj[iProj]   -> GetYaxis() -> SetTitle(sNum.Data());
    hProj[iProj]   -> GetYaxis() -> SetTitleFont(fTxt);
    hProj[iProj]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hProj[iProj]   -> GetYaxis() -> CenterTitle(fCnt);
    hProjW[iProj]  -> SetMarkerColor(fColProj[iProj]);
    hProjW[iProj]  -> SetMarkerStyle(fMarProj[iProj]);
    hProjW[iProj]  -> SetLineColor(fColProj[iProj]);
    hProjW[iProj]  -> SetLineStyle(fLin);
    hProjW[iProj]  -> SetLineWidth(fWid);
    hProjW[iProj]  -> SetFillColor(fColProj[iProj]);
    hProjW[iProj]  -> SetFillStyle(fFil);
    hProjW[iProj]  -> SetTitle(sTitleW.Data());
    hProjW[iProj]  -> SetTitleFont(fTxt);
    hProjW[iProj]  -> GetXaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
    hProjW[iProj]  -> GetXaxis() -> SetLabelFont(fTxt);
    hProjW[iProj]  -> GetXaxis() -> SetLabelSize(fLbl);
    hProjW[iProj]  -> GetXaxis() -> SetTitle(sPar.Data());
    hProjW[iProj]  -> GetXaxis() -> SetTitleFont(fTxt);
    hProjW[iProj]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hProjW[iProj]  -> GetXaxis() -> CenterTitle(fCnt);
    hProjW[iProj]  -> GetYaxis() -> SetLabelFont(fTxt);
    hProjW[iProj]  -> GetYaxis() -> SetLabelSize(fLbl);
    hProjW[iProj]  -> GetYaxis() -> SetTitle(sNum.Data());
    hProjW[iProj]  -> GetYaxis() -> SetTitleFont(fTxt);
    hProjW[iProj]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hProjW[iProj]  -> GetYaxis() -> CenterTitle(fCnt);
    hWeight[iProj] -> SetMarkerColor(fColProj[iProj]);
    hWeight[iProj] -> SetMarkerStyle(fMarProj[iProj]);
    hWeight[iProj] -> SetLineColor(fColProj[iProj]);
    hWeight[iProj] -> SetLineStyle(fLin);
    hWeight[iProj] -> SetLineWidth(fWid);
    hWeight[iProj] -> SetFillColor(fColProj[iProj]);
    hWeight[iProj] -> SetFillStyle(fFil);
    hWeight[iProj] -> SetTitle(sTitleW.Data());
    hWeight[iProj] -> SetTitleFont(fTxt);
    hWeight[iProj] -> GetXaxis() -> SetRangeUser(xyPlot[1], xyPlot[3]);
    hWeight[iProj] -> GetXaxis() -> SetLabelFont(fTxt);
    hWeight[iProj] -> GetXaxis() -> SetLabelSize(fLbl);
    hWeight[iProj] -> GetXaxis() -> SetTitle(sPar.Data());
    hWeight[iProj] -> GetXaxis() -> SetTitleFont(fTxt);
    hWeight[iProj] -> GetXaxis() -> SetTitleOffset(fOffX);
    hWeight[iProj] -> GetXaxis() -> CenterTitle(fCnt);
    hWeight[iProj] -> GetYaxis() -> SetLabelFont(fTxt);
    hWeight[iProj] -> GetYaxis() -> SetLabelSize(fLbl);
    hWeight[iProj] -> GetYaxis() -> SetTitle(sNum.Data());
    hWeight[iProj] -> GetYaxis() -> SetTitleFont(fTxt);
    hWeight[iProj] -> GetYaxis() -> SetTitleOffset(fOffY);
    hWeight[iProj] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legends
  const UInt_t  fColL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fAlnL(12);
  const UInt_t  nObjL(NTrgBins + 2);
  const UInt_t  nObjLT(NTrgBins + 1);
  const UInt_t  nObjLM(NTrgBins + 1);
  const UInt_t  nObjLP(NProj);
  const UInt_t  nObjLB(NTrgBins);
  const UInt_t  nObjLW(NProj + 1);
  const UInt_t  nObjLE(NEffHist);
  const UInt_t  nObjLC(NEtCuts + 1);
  const UInt_t  nObjLS(NTrgBins);
  const UInt_t  nObjLD(NTrgBins + 1);
  const UInt_t  nObjT(NTxt);
  const Float_t hObj(0.05);
  const Float_t hObjL(hObj * nObjL);
  const Float_t hObjLT(hObj * nObjLT);
  const Float_t hObjLM(hObj * nObjLM);
  const Float_t hObjLP(hObj * nObjLP);
  const Float_t hObjLB(hObj * nObjLB);
  const Float_t hObjLW(hObj * nObjLW);
  const Float_t hObjLE(hObj * nObjLE);
  const Float_t hObjLC(hObj * nObjLC);
  const Float_t hObjLS(hObj * nObjLS);
  const Float_t hObjLD(hObj * nObjLD);
  const Float_t hObjT(hObj * nObjT);
  const Float_t yObjL(0.1 + hObjL);
  const Float_t yObjLT(0.1 + hObjLT);
  const Float_t yObjLM(0.1 + hObjLM);
  const Float_t yObjLP(0.1 + hObjLP);
  const Float_t yObjLB(0.1 + hObjLB);
  const Float_t yObjLW(0.1 + hObjLW);
  const Float_t yObjLE(0.1 + hObjLE);
  const Float_t yObjLC(0.1 + hObjLC);
  const Float_t yObjLS(0.1 + hObjLS);
  const Float_t yObjLD(0.1 + hObjLD);
  const Float_t yObjT(0.1 + hObjT);
  const Float_t xyLeg[NVtx]  = {0.1, 0.1, 0.3, yObjL};
  const Float_t xyLegT[NVtx] = {0.1, 0.1, 0.3, yObjLT};
  const Float_t xyLegM[NVtx] = {0.1, 0.1, 0.3, yObjLM};
  const Float_t xyLegP[NVtx] = {0.1, 0.1, 0.3, yObjLP};
  const Float_t xyLegB[NVtx] = {0.1, 0.1, 0.3, yObjLB};
  const Float_t xyLegW[NVtx] = {0.1, 0.1, 0.3, yObjLW};
  const Float_t xyLegE[NVtx] = {0.1, 0.1, 0.3, yObjLE};
  const Float_t xyLegC[NVtx] = {0.3, 0.1, 0.5, yObjLC};
  const Float_t xyLegS[NVtx] = {0.1, 0.1, 0.3, yObjLS};
  const Float_t xyLegD[NVtx] = {0.1, 0.1, 0.3, yObjLD};
  const Float_t xyTxt[NVtx]  = {0.3, 0.1, 0.5, yObjT};

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3], sLabel.Data());
  leg -> SetLineColor(fColL);
  leg -> SetLineStyle(fLinL);
  leg -> SetFillColor(fColL);
  leg -> SetFillStyle(fFilL);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnL);
  leg -> AddEntry(hParticle, sLabelS.Data(), "pf");
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    leg -> AddEntry(hDetector[iTrg], sLabelR[iTrg].Data(), "pf");
  }

  TLegend *legT = new TLegend(xyLegT[0], xyLegT[1], xyLegT[2], xyLegT[3], sHeaderT.Data());
  legT -> SetLineColor(fColL);
  legT -> SetLineStyle(fLinL);
  legT -> SetFillColor(fColL);
  legT -> SetFillStyle(fFilL);
  legT -> SetTextFont(fTxt);
  legT -> SetTextAlign(fAlnL);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legT -> AddEntry(hBias[iTrg], sLabelT[iTrg].Data(), "pf");
  }

  // convert avg. integrals to strings
  TString sLabelI[NTrgBins] = {"", "", ""};
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {

    TString sRaw("");
    TString sTxt("");
    sRaw += avgNorm[iTrg];
    sTxt += sNorm;
    sTxt += " = ";

    const UInt_t nRaw = sRaw.First(".");
    const UInt_t nTxt = nRaw + (nDec + 1);

    sTxt.Append(sRaw, nTxt);
    sLabelI[iTrg] += sEnergy[iTrg];
    sLabelI[iTrg] += ", ";
    sLabelI[iTrg] += sTxt;
  }

  TLegend *legM = new TLegend(xyLegM[0], xyLegM[1], xyLegM[2], xyLegM[3], sHeaderI.Data());
  legM -> SetLineColor(fColL);
  legM -> SetLineStyle(fLinL);
  legM -> SetFillColor(fColL);
  legM -> SetFillStyle(fFilL);
  legM -> SetTextFont(fTxt);
  legM -> SetTextAlign(fAlnL);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legM -> AddEntry((TObject*)0, sLabelI[iTrg].Data(), "");
  }

  TLegend *legP = new TLegend(xyLegP[0], xyLegP[1], xyLegP[2], xyLegP[3]);
  TLegend *legW = new TLegend(xyLegW[0], xyLegW[1], xyLegW[2], xyLegW[3]);
  legP -> SetLineColor(fColL);
  legP -> SetLineStyle(fLinL);
  legP -> SetFillColor(fColL);
  legP -> SetFillStyle(fFilL);
  legP -> SetTextFont(fTxt);
  legP -> SetTextAlign(fAlnL);
  legW -> SetLineColor(fColL);
  legW -> SetLineStyle(fLinL);
  legW -> SetFillColor(fColL);
  legW -> SetFillStyle(fFilL);
  legW -> SetTextFont(fTxt);
  legW -> SetTextAlign(fAlnL);
  legW -> AddEntry(hWeightPar, sLabelS.Data(), "pf");
  for (UInt_t iProj = 0; iProj < NProj; iProj++) {

    // create label
    TString sLabelP(sLabPref.Data());
    sLabelP.Append(sProjSuff[iProj].Data());
    sLabelP.Append(sLabSuff.Data());

    // add to legend
    legP -> AddEntry(hProj[iProj], sLabelP.Data(), "pf");
    legW -> AddEntry(hProjW[iProj], sLabelP.Data(), "pf");
  }

  TLegend *legB = new TLegend(xyLegB[0], xyLegB[1], xyLegB[2], xyLegB[3]);
  legB -> SetLineColor(fColL);
  legB -> SetLineStyle(fLinL);
  legB -> SetFillColor(fColL);
  legB -> SetFillStyle(fFilL);
  legB -> SetTextFont(fTxt);
  legB -> SetTextAlign(fAlnL);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legB -> AddEntry(hEtProj[iTrg], sLabelT[iTrg].Data(), "pf");
  }

  TLegend *legE = new TLegend(xyLegE[0], xyLegE[1], xyLegE[2], xyLegE[3]);
  legE -> SetLineColor(fColL);
  legE -> SetLineStyle(fLinL);
  legE -> SetFillColor(fColL);
  legE -> SetFillStyle(fFilL);
  legE -> SetTextFont(fTxt);
  legE -> SetTextAlign(fAlnL);
  for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
    legE -> AddEntry(hEffDiv[iDiv], sLabelE[iDiv].Data(), "pf");
  }
  legE -> AddEntry(hEff, sLabelE[NEffDiv].Data(), "pf");

  TLegend *legS = new TLegend(xyLegS[0], xyLegS[1], xyLegS[2], xyLegS[3]);
  legS -> SetLineColor(fColL);
  legS -> SetLineStyle(fLinL);
  legS -> SetFillColor(fColL);
  legS -> SetFillStyle(fFilL);
  legS -> SetTextFont(fTxt);
  legS -> SetTextAlign(fAlnL);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legS -> AddEntry(hSmearWeights[iTrg], sLegSmear[iTrg].Data(), "pf");
  }

  TLegend *legD = new TLegend(xyLegD[0], xyLegD[1], xyLegD[2], xyLegD[3]);
  legD -> SetLineColor(fColL);
  legD -> SetLineStyle(fLinL);
  legD -> SetFillColor(fColL);
  legD -> SetFillStyle(fFilL);
  legD -> SetTextFont(fTxt);
  legD -> SetTextAlign(fAlnL);
  legD -> AddEntry(hParticle, sLegPar.Data(), "pf");
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legD -> AddEntry(hSmearWeights[iTrg], sLegSmear[iTrg].Data(), "pf");
  }

  TPaveText *ptTxt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  ptTxt -> SetLineColor(fColL);
  ptTxt -> SetLineStyle(fLinL);
  ptTxt -> SetFillColor(fColL);
  ptTxt -> SetFillStyle(fFilL);
  ptTxt -> SetTextFont(fTxt);
  ptTxt -> SetTextAlign(fAlnL);
  ptTxt -> AddText(sTxtSim.Data());
  ptTxt -> AddText(sTxtCut.Data());
  cout << "    Made legends and text." << endl;

  // make lines
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  const Float_t yLineW(1.);
  const Float_t yLineE(1.);
  const Float_t yLineS(1.);

  TLine *lMatrix = new TLine(xyPlot[0], xyPlot[1], xyPlot[2], xyPlot[3]);
  lMatrix -> SetLineColor(fColLi);
  lMatrix -> SetLineStyle(fLinLi);
  lMatrix -> SetLineWidth(fWidLi);

  TLine *lBias = new TLine(yEtParRange[0], yLineW, yEtParRange[1], yLineW);
  lBias -> SetLineColor(fColLi);
  lBias -> SetLineStyle(fLinLi);
  lBias -> SetLineWidth(fWidLi);

  TLine *lWeight = new TLine(yEtParRange[0], yLineW, yEtParRange[1], yLineW);
  lWeight -> SetLineColor(fColLi);
  lWeight -> SetLineStyle(fLinLi);
  lWeight -> SetLineWidth(fWidLi);

  TLine *lEff = new TLine(yEtParRange[0], yLineE, yEtParRange[1], yLineE);
  lEff -> SetLineColor(fColLi);
  lEff -> SetLineStyle(fLinLi);
  lEff -> SetLineWidth(fWidLi);

  TLine *lSmear = new TLine(yEtParRange[0], yLineS, yEtParRange[1], yLineS);
  lSmear -> SetLineColor(fColLi);
  lSmear -> SetLineStyle(fLinLi);
  lSmear -> SetLineWidth(fWidLi);
  cout << "    Made lines." << endl;

  // make cut indicators
  TBox *bEtCut = new TBox(xEtCutRange[0], xyPlot[1], xEtCutRange[1], xyPlot[3]);
  bEtCut -> SetFillColor(fEtCutCol);
  bEtCut -> SetFillStyle(fEtCutFil);
  bEtCut -> SetLineColor(fEtCutCol);
  bEtCut -> SetLineStyle(fEtCutLin);
  bEtCut -> SetLineWidth(fEtCutWid);
  cout << "    Made cut indicators." << endl;

  // make plots
  const UInt_t  widthM(950);
  const UInt_t  widthT(1500);
  const UInt_t  widthP(950);
  const UInt_t  widthB(950);
  const UInt_t  widthW(1500);
  const UInt_t  widthE(750);
  const UInt_t  widthS(950);
  const UInt_t  heightM(950);
  const UInt_t  heightT(750);
  const UInt_t  heightP(950);
  const UInt_t  heightB(950);
  const UInt_t  heightW(750);
  const UInt_t  heightE(950);
  const UInt_t  heightS(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogYN(1);
  const UInt_t  fLogYT(0);
  const UInt_t  fLogYM(0);
  const UInt_t  fLogYP(1);
  const UInt_t  fLogYB(1);
  const UInt_t  fLogYW(1);
  const UInt_t  fLogYE(1);
  const UInt_t  fLogYS(1);
  const UInt_t  fLogZ(1);
  const Float_t fMarginTN(0.02);
  const Float_t fMarginTT(0.02);
  const Float_t fMarginTP(0.02);
  const Float_t fMarginTB(0.02);
  const Float_t fMarginTW(0.02);
  const Float_t fMarginTD(0.02);
  const Float_t fMarginTE(0.005);
  const Float_t fMarginTS(0.02);
  const Float_t fMarginRN(0.02);
  const Float_t fMarginRT(0.02);
  const Float_t fMarginRP(0.02);
  const Float_t fMarginRB(0.02);
  const Float_t fMarginRW(0.02);
  const Float_t fMarginRD(0.02);
  const Float_t fMarginRE(0.02);
  const Float_t fMarginRS(0.02);
  const Float_t fMarginBN(0.15);
  const Float_t fMarginBT(0.15);
  const Float_t fMarginBP(0.15);
  const Float_t fMarginBB(0.15);
  const Float_t fMarginBW(0.15);
  const Float_t fMarginBD(0.005);
  const Float_t fMarginBE(0.25);
  const Float_t fMarginBS(0.15);
  const Float_t fMarginLN(0.15);
  const Float_t fMarginLT(0.15);
  const Float_t fMarginLP(0.15);
  const Float_t fMarginLB(0.15);
  const Float_t fMarginLW(0.15);
  const Float_t fMarginLD(0.15);
  const Float_t fMarginLE(0.15);
  const Float_t fMarginLS(0.15);
  const Float_t fPadXYN[NVtx] = {0., 0., 0.5, 1.};
  const Float_t fPadXYT[NVtx] = {0.5, 0., 1., 1.};
  const Float_t fPadXYP[NVtx] = {0., 0., 0.5, 1.};
  const Float_t fPadXYW[NVtx] = {0.5, 0., 1., 1.};
  const Float_t fPadXYD[NVtx] = {0., 0.3, 1., 1.};
  const Float_t fPadXYE[NVtx] = {0., 0., 1., 0.3};

  TCanvas *cMatrix = new TCanvas("cMatrix", "", widthM, heightM);
  cMatrix -> SetGrid(fGrid, fGrid);
  cMatrix -> SetTicks(fTick, fTick);
  cMatrix -> SetBorderMode(fMode);
  cMatrix -> SetBorderSize(fBord);
  cMatrix -> SetLogx(fLogX);
  cMatrix -> SetLogy(fLogYM);
  cMatrix -> SetLogz(fLogZ);
  cMatrix -> cd();
  hMatrix -> Draw("colz");
  lMatrix -> Draw();
  legM    -> Draw();
  bEtCut  -> Draw();
  ptTxt   -> Draw();
  fOut    -> cd();
  cMatrix -> Write();
  cMatrix -> Close();

  TCanvas *cMatrixW = new TCanvas("cWeightMatrix", "", widthM, heightM);
  cMatrixW -> SetGrid(fGrid, fGrid);
  cMatrixW -> SetTicks(fTick, fTick);
  cMatrixW -> SetBorderMode(fMode);
  cMatrixW -> SetBorderSize(fBord);
  cMatrixW -> SetLogx(fLogX);
  cMatrixW -> SetLogy(fLogYM);
  cMatrixW -> SetLogz(fLogZ);
  cMatrixW -> cd();
  hMatrixW -> Draw("colz");
  lMatrix  -> Draw();
  fOut     -> cd();
  cMatrixW -> Write();
  cMatrixW -> Close();

  TCanvas *cBias  = new TCanvas("cTspBias", "", widthT, heightT);
  TPad    *pDists = new TPad("pDistributions", "", fPadXYN[0], fPadXYN[1], fPadXYN[2], fPadXYN[3]);
  TPad    *pBias  = new TPad("pTspBias", "", fPadXYT[0], fPadXYT[1], fPadXYT[2], fPadXYT[3]);
  pDists    -> SetGrid(fGrid, fGrid);
  pDists    -> SetTicks(fTick, fTick);
  pDists    -> SetBorderMode(fMode);
  pDists    -> SetBorderSize(fBord);
  pDists    -> SetTopMargin(fMarginTN);
  pDists    -> SetRightMargin(fMarginRN);
  pDists    -> SetBottomMargin(fMarginBN);
  pDists    -> SetLeftMargin(fMarginLN);
  pDists    -> SetLogx(fLogX);
  pDists    -> SetLogy(fLogYN);
  pBias     -> SetGrid(fGrid, fGrid);
  pBias     -> SetTicks(fTick, fTick);
  pBias     -> SetBorderMode(fMode);
  pBias     -> SetBorderSize(fBord);
  pBias     -> SetTopMargin(fMarginTT);
  pBias     -> SetRightMargin(fMarginRT);
  pBias     -> SetBottomMargin(fMarginBT);
  pBias     -> SetLeftMargin(fMarginLT);
  pBias     -> SetLogx(fLogX);
  pBias     -> SetLogy(fLogYT);
  cBias     -> cd();
  pDists    -> Draw();
  pBias     -> Draw();
  pDists    -> cd();
  hParticle -> Draw();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hDetector[iTrg] -> Draw("same");
  }
  leg      -> Draw();
  pBias    -> cd();
  hBias[0] -> Draw();
  for (UInt_t iTrg = 1; iTrg < NTrgBins; iTrg++) {
    hBias[iTrg] -> Draw("same");
  }
  lBias -> Draw();
  legT  -> Draw();
  fOut  -> cd();
  cBias -> Write();
  cBias -> Close();

  TCanvas *cProj = new TCanvas("cFineProjections", "", widthP, heightP);
  cProj    -> SetGrid(fGrid, fGrid);
  cProj    -> SetTicks(fTick, fTick);
  cProj    -> SetBorderMode(fMode);
  cProj    -> SetBorderSize(fBord);
  cProj    -> SetTopMargin(fMarginTP);
  cProj    -> SetRightMargin(fMarginRP);
  cProj    -> SetBottomMargin(fMarginBP);
  cProj    -> SetLeftMargin(fMarginLP);
  cProj    -> SetLogx(fLogX);
  cProj    -> SetLogy(fLogYP);
  cProj    -> cd();
  hProj[0] -> Draw("hist p");
  hProj[0] -> Draw("hist l sames");
  for (UInt_t iProj = 1; iProj < NProj; iProj++) {
    hProj[iProj] -> Draw("hist p same");
    hProj[iProj] -> Draw("hist l sames");
  }
  legP  -> Draw();
  fOut  -> cd();
  cProj -> Write();
  cProj -> Close();

  TCanvas *cEtProj = new TCanvas("cEtBinProjections", "", widthB, heightB);
  cEtProj    -> SetGrid(fGrid, fGrid);
  cEtProj    -> SetTicks(fTick, fTick);
  cEtProj    -> SetBorderMode(fMode);
  cEtProj    -> SetBorderSize(fBord);
  cEtProj    -> SetTopMargin(fMarginTB);
  cEtProj    -> SetRightMargin(fMarginRB);
  cEtProj    -> SetBottomMargin(fMarginBB);
  cEtProj    -> SetLeftMargin(fMarginLB);
  cEtProj    -> SetLogx(fLogX);
  cEtProj    -> SetLogy(fLogYB);
  cEtProj    -> cd();
  hEtProj[0] -> Draw();
  for (UInt_t iTrg = 1; iTrg < NTrgBins; iTrg++) {
    hEtProj[iTrg] -> Draw("sames");
  }
  legB    -> Draw();
  fOut    -> cd();
  cEtProj -> Write();
  cEtProj -> Close();

  TCanvas *cWeight = new TCanvas("cWeights", "", widthW, heightW);
  TPad    *pProj   = new TPad("pProjections", "", fPadXYP[0], fPadXYP[1], fPadXYP[2], fPadXYP[3]);
  TPad    *pWeight = new TPad("pWeight", "", fPadXYW[0], fPadXYW[1], fPadXYW[2], fPadXYW[3]);
  pProj      -> SetGrid(fGrid, fGrid);
  pProj      -> SetTicks(fTick, fTick);
  pProj      -> SetBorderMode(fMode);
  pProj      -> SetBorderSize(fBord);
  pProj      -> SetTopMargin(fMarginTW);
  pProj      -> SetRightMargin(fMarginRW);
  pProj      -> SetBottomMargin(fMarginBW);
  pProj      -> SetLeftMargin(fMarginLW);
  pProj      -> SetLogx(fLogX);
  pProj      -> SetLogy(fLogYW);
  pWeight    -> SetGrid(fGrid, fGrid);
  pWeight    -> SetTicks(fTick, fTick);
  pWeight    -> SetBorderMode(fMode);
  pWeight    -> SetBorderSize(fBord);
  pWeight    -> SetTopMargin(fMarginTW);
  pWeight    -> SetRightMargin(fMarginRW);
  pWeight    -> SetBottomMargin(fMarginBW);
  pWeight    -> SetLeftMargin(fMarginLW);
  pWeight    -> SetLogx(fLogX);
  pWeight    -> SetLogy(fLogYW);
  cWeight    -> cd();
  pProj      -> Draw();
  pWeight    -> Draw();
  pProj      -> cd();
  hWeightPar -> Draw();
  for (UInt_t iProj = 0; iProj < NProj; iProj++) {
    hProjW[iProj] -> Draw("hist p same");
    hProjW[iProj] -> Draw("hist l same");
    hProjW[iProj] -> Draw("same");
  }
  legW       -> Draw();
  pWeight    -> cd();
  hWeight[0] -> Draw("hist p");
  hWeight[0] -> Draw("hist l same");
  hWeight[0] -> Draw("same");
  for (UInt_t iProj = 1; iProj < NProj; iProj++) {
    hWeight[iProj] -> Draw("hist p same");
    hWeight[iProj] -> Draw("hist l same");
    hWeight[iProj] -> Draw("same");
  }
  lWeight -> Draw();
  fOut    -> cd();
  cWeight -> Write();
  cWeight -> Close();

  TCanvas *cEff = new TCanvas("cEfficiency", "", widthE, heightE);
  TPad    *pDiv = new TPad("pDiv", "", fPadXYD[0], fPadXYD[1], fPadXYD[2], fPadXYD[3]);
  TPad    *pEff = new TPad("pEff", "", fPadXYE[0], fPadXYE[1], fPadXYE[2], fPadXYE[3]);
  pDiv -> SetGrid(fGrid, fGrid);
  pDiv -> SetTicks(fTick, fTick);
  pDiv -> SetBorderMode(fMode);
  pDiv -> SetBorderSize(fBord);
  pDiv -> SetTopMargin(fMarginTD);
  pDiv -> SetRightMargin(fMarginRD);
  pDiv -> SetBottomMargin(fMarginBD);
  pDiv -> SetLeftMargin(fMarginLD);
  pDiv -> SetLogx(fLogX);
  pDiv -> SetLogy(fLogYE);
  pEff -> SetGrid(fGrid, fGrid);
  pEff -> SetTicks(fTick, fTick);
  pEff -> SetBorderMode(fMode);
  pEff -> SetBorderSize(fBord);
  pEff -> SetTopMargin(fMarginTE);
  pEff -> SetRightMargin(fMarginRE);
  pEff -> SetBottomMargin(fMarginBE);
  pEff -> SetLeftMargin(fMarginLE);
  pEff -> SetLogx(fLogX);
  pEff -> SetLogy(fLogYE);
  cEff -> cd();
  pDiv -> Draw();
  pEff -> Draw();
  pDiv -> cd();
  for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
    if (iDiv == 0) {
      hEffDiv[iDiv] -> Draw();
    } else {
      hEffDiv[iDiv] -> Draw("same");
    }
  }
  legE -> Draw();
  pEff -> cd();
  hEff -> Draw();
  lEff -> Draw();
  fOut -> cd();
  cEff -> Write();
  cEff -> Close();

  TCanvas *cSmear = new TCanvas("cSmearWeights", "", widthS, heightS);
  cSmear           -> SetGrid(fGrid, fGrid);
  cSmear           -> SetTicks(fTick, fTick);
  cSmear           -> SetBorderMode(fMode);
  cSmear           -> SetBorderSize(fBord);
  cSmear           -> SetTopMargin(fMarginTS);
  cSmear           -> SetRightMargin(fMarginRS);
  cSmear           -> SetBottomMargin(fMarginBS);
  cSmear           -> SetLeftMargin(fMarginLS);
  cSmear           -> SetLogx(fLogX);
  cSmear           -> SetLogy(fLogYS);
  cSmear           -> cd();
  hSmearWeights[0] -> Draw();
  for (UInt_t iTrg = 1; iTrg < NTrgBins; iTrg++) {
    hSmearWeights[iTrg] -> Draw("same");
  }
  lSmear -> Draw();
  legS   -> Draw();
  ptTxt  -> Draw();
  fOut   -> cd();
  cSmear -> Write();
  cSmear -> Close();

  TCanvas *cSmeDists = new TCanvas("cSmearDistributions", "", widthS, heightS);
  cSmeDists -> SetGrid(fGrid, fGrid);
  cSmeDists -> SetTicks(fTick, fTick);
  cSmeDists -> SetBorderMode(fMode);
  cSmeDists -> SetBorderSize(fBord);
  cSmeDists -> SetTopMargin(fMarginTS);
  cSmeDists -> SetRightMargin(fMarginRS);
  cSmeDists -> SetBottomMargin(fMarginBS);
  cSmeDists -> SetLeftMargin(fMarginLS);
  cSmeDists -> SetLogx(fLogX);
  cSmeDists -> SetLogy(fLogYS);
  cSmeDists -> cd();
  hParticle -> Draw();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hSmeared[iTrg] -> Draw("same");
  }
  legD      -> Draw();
  ptTxt     -> Draw();
  fOut      -> cd();
  cSmeDists -> Write();
  cSmeDists -> Close();
  cout << "    Made plots." << endl;

  // save histograms and functions
  fOut       -> cd();
  hEff       -> Write();
  hMatrix    -> Write();
  hMatrixW   -> Write();
  hParticle  -> Write();
  hWeightPar -> Write();
  for (UInt_t iDiv = 0; iDiv < NEffDiv; iDiv++) {
    hEffDiv[iDiv] -> Write();
  }
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hSmeared[iTrg]      -> Write();
    hDetector[iTrg]     -> Write();
    hParForNorm[iTrg]   -> Write();
    hBias[iTrg]         -> Write();
    hSmearWeights[iTrg] -> Write();
    hEtProj[iTrg]       -> Write();
    if (smoothWeights) {
      for (UInt_t iSmooth = 0; iSmooth < nSmooth[iTrg]; iSmooth++) {
        fSmoothWeight[iTrg][iSmooth] -> Write();
      }
    }
  }
  for (UInt_t iProj = 0; iProj < NProj; iProj++) {
    hProj[iProj]   -> Write();
    hProjW[iProj]  -> Write();
    hWeight[iProj] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished TES/R matrix prepartion!\n" << endl;

}

// End ------------------------------------------------------------------------
