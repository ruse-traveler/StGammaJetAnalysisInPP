// 'DoMatrixProjections.C'
// Derek Anderson
// 11.29.2022
//
// Use this to project slices of the
// TES/R matrix, weight them, and
// calculate the corresponding
// smearing weights.
//
// NOTE: this assumes that the input
// matrix has already been normalized.
//
// NOTE: TrgType controls trigger
// species
//   TrgType = 0: pi0
//   TrgType = 1: gamma

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
const UInt_t NProjMax(11);
const UInt_t NTrgBins(3);
const UInt_t NRange(2);
const UInt_t NVtx(4);
const UInt_t NTxt(2);

// smoothing constants
static const UInt_t NParSmooth(3);
static const UInt_t NSmoothMax(4);
static const UInt_t NSmoothPi0[NTrgBins] = {4, 3, 2};
static const UInt_t NSmoothGam[NTrgBins] = {3, 2, 1};

// trigger species
static const UInt_t TrgType(0);



void DoMatrixProjections() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning matrix projection..." << endl;

  // i/o parameters
  const TString sOut("matrixProjections.forPaper_noEt1520.et650pi0.d30m11y2022.root");
  const TString sInMatrix("triggerMatrix.forWeightPlot_matrixErrorFixed.et650x650vz55tsp008pi0.d30m11y2022.root");
  const TString sInParticle("triggerMatrix.forWeightPlot_matrixErrorFixed.et650x650vz55tsp008pi0.d30m11y2022.root");
  const TString sMatrix("hTriggerMatrix");
  const TString sParticle("hEtPar");
  const TString sOutMat("hMatrix");
  const TString sOutPar("hEtParticle");

  // projection parameters
  const TString sAvgName[NTrgBins]   = {"hAverageEt911", "hAverageEt1115", "hAverageEt1520"};
  const TString sSmearName[NTrgBins] = {"hSmearEt911",   "hSmearEt1115",   "hSmearEt1520"};
  const TString sProjBase[NTrgBins]  = {"hProjEt911",    "hProjEt1115",    "hProjEt1520"};
  const Float_t eTdetMin[NTrgBins]   = {9.,  11., 15.};
  const Float_t eTdetMax[NTrgBins]   = {11., 15., 20.};

  // plot parameter
  const UInt_t  fColPar(923);
  const UInt_t  fMarPar(29);
  const UInt_t  fColPi0[NTrgBins] = {859, 899, 819};
  const UInt_t  fColGam[NTrgBins] = {799, 899, 859}; 
  const UInt_t  fMar[NTrgBins]    = {20,  22,  21};
  const Float_t xPlot[NRange]     = {6., 30.};
  const TString sParP("E_{T}^{#pi^{0},part-match} [GeV]");
  const TString sParG("E_{T}^{#gamma,part-match} [GeV]");

  // text parameters
  const TString sTxtSimPi0("STAR simulation #pi^{0}, 6 < E_{T}^{#pi^{0},part,det-clust} < 50 GeV"); 
  const TString sTxtSimGam("STAR simulation #gamma, 6 < E_{T}^{#gamma,part,det-clust} < 50 GeV");
  const TString sTxtCutPi0("|#eta^{#pi^{0}}| < 0.9, TSP #in (0, 0.08)");
  const TString sTxtCutGam("|#eta^{#gamma}| < 0.9, TSP #in (0.2, 0.6)");
  const TString sLegParSP("generated #pi^{0}");
  const TString sLegParSG("generated #gamma");
  const TString sLegP[NTrgBins]   = {"9 < E_{T}^{#pi^{0},det-clust} < 11 GeV", "11 < E_{T}^{#pi^{0},det-clust} < 15 GeV", "15 < E_{T}^{#pi^{0},det-clust} < 20 GeV"};
  const TString sLegG[NTrgBins]   = {"9 < E_{T}^{#gamma,det-clust} < 11 GeV",  "11 < E_{T}^{#gamma,det-clust} < 15 GeV",  "15 < E_{T}^{#gamma,det-clust} < 20 GeV"};

  // rebinning parameters
  const Bool_t rebinPar(true);
  const UInt_t nRebinPar(2);

  // generic smoothing parameters
  const Bool_t  smoothWeights(false);
  const UInt_t  fLinSmooth(1);
  const UInt_t  fWidSmooth(2);
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
  UInt_t  fCol[NTrgBins];
  UInt_t  nSmooth[NTrgBins];
  TString sPar("");
  TString sLegPar("");
  TString sTxtSim("");
  TString sTxtCut("");
  TString sLeg[NTrgBins];
  Float_t rSmoothRangesLo[NTrgBins][NSmoothMax];
  Float_t rSmoothRangesHi[NTrgBins][NSmoothMax];
  Float_t rCorrectRangesLo[NTrgBins][NSmoothMax];
  Float_t rCorrectRangesHi[NTrgBins][NSmoothMax];
  switch (TrgType) {
    case 0:
      sPar    = sParP;
      sLegPar = sLegParSP;
      sTxtSim = sTxtSimPi0;
      sTxtCut = sTxtCutPi0;
      for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
        fCol[iTrg] = fColPi0[iTrg];
        sLeg[iTrg] = sLegP[iTrg];
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
      sPar    = sParG;
      sLegPar = sLegParSG;
      sTxtSim = sTxtSimGam;
      sTxtCut = sTxtCutGam;
      for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
        fCol[iTrg] = fColGam[iTrg];
        sLeg[iTrg] = sLegG[iTrg];
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
      sPar    = sParP;
      sLegPar = sLegParSP;
      sTxtSim = sTxtSimPi0;
      sTxtCut = sTxtCutPi0;
      for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
        fCol[iTrg] = fColPi0[iTrg];
        sLeg[iTrg] = sLegP[iTrg];
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
  TFile *fOut      = new TFile(sOut.Data(), "recreate");
  TFile *fMatrix   = new TFile(sInMatrix.Data(), "read");
  TFile *fParticle = new TFile(sInParticle.Data(), "read");
  if (!fOut || !fMatrix || !fParticle) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fMatrix = " << fMatrix << ", fParticle = " << fParticle << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH2D *hMatrix   = (TH2D*) fMatrix   -> Get(sMatrix.Data());
  TH1D *hParticle = (TH1D*) fParticle -> Get(sParticle.Data());
  if (!hMatrix || !hParticle) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hMatrix = << " hMatrix << ", hParticle = " << hParticle << "\n"
         << endl;
    return;
  }
  hMatrix   -> SetName(sOutMat.Data());
  hParticle -> SetName(sOutPar.Data());
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  if (rebinPar) {
    hParticle -> Rebin(nRebinPar);
    cout << "    Rebinned particle histogram." << endl;
  }

  // determine number of projections per bin
  UInt_t nProj[NTrgBins];
  UInt_t iProj[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    iProj[iTrg] = 0;
    nProj[iTrg] = 0;
  }

  const UInt_t nBinX = hMatrix -> GetNbinsX();
  const UInt_t nBinY = hMatrix -> GetNbinsY();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    for (UInt_t iBinX = 1; iBinX < (nBinX + 1); iBinX++) {
      const Float_t xCenter = hMatrix -> GetXaxis() -> GetBinCenter(iBinX);
      const Bool_t  isInBin = ((xCenter > eTdetMin[iTrg]) && (xCenter < eTdetMax[iTrg]));
      if (isInBin) {
        ++nProj[iTrg];
      } else {
        continue;
      }
    }  // end eTdet loop
  }  // end trigger loop

  cout << "    Determined number of projections to do:" << endl;
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    cout << "      nProj[" << iTrg << "] = " << nProj[iTrg] << endl;
  }

  // project out matrix
  TH1D *hProj[NTrgBins][NProjMax];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    for (UInt_t iBinX = 1; iBinX < (nBinX + 1); iBinX++) {

      // check if in bin
      const Float_t xCenter = hMatrix -> GetXaxis() -> GetBinCenter(iBinX);
      const Bool_t  isInBin = ((xCenter > eTdetMin[iTrg]) && (xCenter < eTdetMax[iTrg]));
      if (!isInBin) continue;

      // create name
      TString sProjName(sProjBase[iTrg].Data());
      sProjName += "_";
      sProjName += iProj[iTrg];

      // do projection
      hProj[iTrg][iProj[iTrg]] = (TH1D*) hMatrix -> ProjectionY("", iBinX, iBinX, "e") -> Clone();
      hProj[iTrg][iProj[iTrg]] -> SetName(sProjName.Data());
      ++iProj[iTrg];

    }  // end eTdet loop
  }  // end trigger loop
  cout << "    Did projections." << endl;

  // do averages and normalize
  const Double_t intParAll = hParticle -> Integral();
  const Double_t binWidthP = hParticle -> GetBinWidth(1);
  const Double_t normPar   = 1. / intParAll;
  const Double_t xCheck    = binWidthP / 4.;
  hParticle -> Scale(normPar);

  TH1D *hAverage[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {

    // initialize histogram
    hAverage[iTrg] = (TH1D*) hProj[iTrg][0] -> Clone();
    hAverage[iTrg] -> SetName(sAvgName[iTrg].Data());
    hAverage[iTrg] -> Reset("ICES");

    // sum projections
    for (UInt_t iProjBin = 0; iProjBin < nProj[iTrg]; iProjBin++) {
      hAverage[iTrg] -> Add(hProj[iTrg][iProjBin]);
    }
    hAverage[iTrg] -> Scale(1. / (Double_t) nProj[iTrg]);

    // normalize
    const Float_t xStart = eTdetMin[iTrg] + xCheck;
    const Float_t xStop  = eTdetMax[iTrg] - xCheck;
    const UInt_t  iStart = hParticle -> FindBin(xStart);
    const UInt_t  iStop  = hParticle -> FindBin(xStop);

    const Double_t intPar = hParticle      -> Integral(iStart, iStop);
    const Double_t intAvg = hAverage[iTrg] -> Integral();
    const Double_t norm   = intPar / intAvg;
    hAverage[iTrg] -> Scale(norm);
  }
  cout << "    Calculated averages." << endl;

  // calculate smearing weights
  TH1D *hSmear[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {

    // initialize histogram
    hSmear[iTrg] = (TH1D*) hAverage[iTrg] -> Clone();
    hSmear[iTrg] -> SetName(sSmearName[iTrg].Data());

    // normalize
    const Double_t intSmear  = hSmear[iTrg] -> Integral();
    const Double_t normSmear = 1. / intSmear;
    hSmear[iTrg] -> Scale(normSmear);

  }  // end trigger loop
  cout << "    Calculated smearing weights." << endl;

  // smooth weights (if need be)
  TF1  *fSmooth[NTrgBins][NSmoothMax];
  if (smoothWeights) {
    for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
      for (UInt_t iSmooth = 0; iSmooth < nSmooth[iTrg]; iSmooth++) {

        // create name
        TString sSmoothFuncName(sSmoothFuncBase[iTrg].Data());
        sSmoothFuncName += "_";
        sSmoothFuncName += iSmooth;

        // initialize smoothing functions
        fSmooth[iTrg][iSmooth] = new TF1(sSmoothFuncName.Data(), sSmoothFunc.Data(), xPlot[0], xPlot[1]);
        fSmooth[iTrg][iSmooth] -> SetLineColor(fCol[iTrg]);
        fSmooth[iTrg][iSmooth] -> SetLineStyle(fLinSmooth);
        fSmooth[iTrg][iSmooth] -> SetLineWidth(fWidSmooth);
        for (UInt_t iSmoothPar = 0; iSmoothPar < NParSmooth; iSmoothPar++) {
          fSmooth[iTrg][iSmooth] -> SetParameter(iSmoothPar, pSmoothParGuess[iSmoothPar]);
        }

        // fit weights
        Int_t fitStatus(-1);
        if (iSmooth == 0) {
          fitStatus = (Int_t) hSmear[iTrg] -> Fit(sSmoothFuncName.Data(), "", "", rSmoothRangesLo[iTrg][iSmooth], rSmoothRangesHi[iTrg][iSmooth]);
        } else {
          fitStatus = (Int_t) hSmear[iTrg] -> Fit(sSmoothFuncName.Data(), "+", "", rSmoothRangesLo[iTrg][iSmooth], rSmoothRangesHi[iTrg][iSmooth]);
        }
        const Bool_t isGoodFit = (fitStatus == 0);

        // apply smoothing
        for (UInt_t iBinS = 1; iBinS < (nBinsP + 1); iBinS++) {

          // get bin info
          const Double_t xCenterS = hSmear[iTrg]           -> GetBinCenter(iBinS);
          const Double_t yRawS    = hSmear[iTrg]           -> GetBinContent(iBinS);
          const Double_t ySmooth  = fSmooth[iTrg][iSmooth] -> Eval(xCenterS);
          const Double_t eRawS    = hSmear[iTrg]           -> GetBinError(iBinS);
          const Double_t ePerRawS = eRawS / yRawS;
          const Double_t eSmooth  = ePerRawS * ySmooth;

          // if in range and fit is good, apply smoothing
          const Bool_t isBinNonzero     = ((yRawS > 0.) && (ySmooth > 0.));
          const Bool_t isInCorrectRange = ((xCenterS > rCorrectRangesLo[iTrg][iSmooth]) && (xCenterS < rCorrectRangesHi[iTrg][iSmooth]));
          if (isBinNonzero && isGoodFit && isInCorrectRange) {
            hSmear[iTrg] -> SetBinContent(iBinS, ySmooth);
            hSmear[iTrg] -> SetBinError(iBinS, eSmooth);
          }
        }  // end eTpar loop
      }  // end smooth loop
    }  // end trigger loop
    cout << "    Smoothed smear weights." << endl;
  }

  // set styles
  const TString sTitle("");
  const TString sNum("probability density");
  const TString sWei("arbitrary units");
  const Float_t fLbl(0.03);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.4);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fFil(0);
  hParticle -> SetMarkerColor(fColPar);
  hParticle -> SetMarkerStyle(fMarPar);
  hParticle -> SetLineColor(fColPar);
  hParticle -> SetLineStyle(fLin);
  hParticle -> SetLineWidth(fWid);
  hParticle -> SetFillColor(fColPar);
  hParticle -> SetFillStyle(fFil);
  hParticle -> SetTitle(sTitle.Data());
  hParticle -> SetTitleFont(fTxt);
  hParticle -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
  hParticle -> GetXaxis() -> SetLabelFont(fTxt);
  hParticle -> GetXaxis() -> SetLabelSize(fLbl);
  hParticle -> GetXaxis() -> SetTitle(sPar.Data());
  hParticle -> GetXaxis() -> SetTitleFont(fTxt);
  hParticle -> GetXaxis() -> SetTitleOffset(fOffX);
  hParticle -> GetXaxis() -> CenterTitle(fCnt);
  hParticle -> GetYaxis() -> SetLabelFont(fTxt);
  hParticle -> GetYaxis() -> SetLabelSize(fLbl);
  hParticle -> GetYaxis() -> SetTitle(sNum.Data());
  hParticle -> GetYaxis() -> SetTitleFont(fTxt);
  hParticle -> GetYaxis() -> SetTitleOffset(fOffY);
  hParticle -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hAverage[iTrg] -> SetMarkerColor(fCol[iTrg]);
    hAverage[iTrg] -> SetMarkerStyle(fMar[iTrg]);
    hAverage[iTrg] -> SetLineColor(fCol[iTrg]);
    hAverage[iTrg] -> SetLineStyle(fLin);
    hAverage[iTrg] -> SetLineWidth(fWid);
    hAverage[iTrg] -> SetFillColor(fCol[iTrg]);
    hAverage[iTrg] -> SetFillStyle(fFil);
    hAverage[iTrg] -> SetTitle(sTitle.Data());
    hAverage[iTrg] -> SetTitleFont(fTxt);
    hAverage[iTrg] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hAverage[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hAverage[iTrg] -> GetXaxis() -> SetLabelSize(fLbl);
    hAverage[iTrg] -> GetXaxis() -> SetTitle(sPar.Data());
    hAverage[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hAverage[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hAverage[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hAverage[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hAverage[iTrg] -> GetYaxis() -> SetLabelSize(fLbl);
    hAverage[iTrg] -> GetYaxis() -> SetTitle(sNum.Data());
    hAverage[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hAverage[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hAverage[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hSmear[iTrg]   -> SetMarkerColor(fCol[iTrg]);
    hSmear[iTrg]   -> SetMarkerStyle(fMar[iTrg]);
    hSmear[iTrg]   -> SetLineColor(fCol[iTrg]);
    hSmear[iTrg]   -> SetLineStyle(fLin);
    hSmear[iTrg]   -> SetLineWidth(fWid);
    hSmear[iTrg]   -> SetFillColor(fCol[iTrg]);
    hSmear[iTrg]   -> SetFillStyle(fFil);
    hSmear[iTrg]   -> SetTitle(sTitle.Data());
    hSmear[iTrg]   -> SetTitleFont(fTxt);
    hSmear[iTrg]   -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hSmear[iTrg]   -> GetXaxis() -> SetLabelFont(fTxt);
    hSmear[iTrg]   -> GetXaxis() -> SetLabelSize(fLbl);
    hSmear[iTrg]   -> GetXaxis() -> SetTitle(sPar.Data());
    hSmear[iTrg]   -> GetXaxis() -> SetTitleFont(fTxt);
    hSmear[iTrg]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hSmear[iTrg]   -> GetXaxis() -> CenterTitle(fCnt);
    hSmear[iTrg]   -> GetYaxis() -> SetLabelFont(fTxt);
    hSmear[iTrg]   -> GetYaxis() -> SetLabelSize(fLbl);
    hSmear[iTrg]   -> GetYaxis() -> SetTitle(sNum.Data());
    hSmear[iTrg]   -> GetYaxis() -> SetTitleFont(fTxt);
    hSmear[iTrg]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hSmear[iTrg]   -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make text boxes/legends
  const UInt_t  fColL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fAlnL(12);
  const UInt_t  nObjLA(NTrgBins + 1);
  const UInt_t  nObjLS(NTrgBins);
  const UInt_t  nObjT(NTxt);
  const Float_t hObj(0.05);
  const Float_t hObjLA(hObj * nObjLA);
  const Float_t hObjLS(hObj * nObjLS);
  const Float_t hObjT(hObj * nObjT);
  const Float_t yObjLA(0.1 + hObjLA);
  const Float_t yObjLS(0.1 + hObjLS);
  const Float_t yObjT(0.1 + hObjT);
  const Float_t xyLegA[NVtx] = {0.1, 0.1, 0.3, yObjLA};
  const Float_t xyLegS[NVtx] = {0.1, 0.1, 0.3, yObjLS};
  const Float_t xyTxt[NVtx]  = {0.3, 0.1, 0.5, yObjT};

  TLegend *legA = new TLegend(xyLegA[0], xyLegA[1], xyLegA[2], xyLegA[3]);
  legA -> SetLineColor(fColL);
  legA -> SetLineStyle(fLinL);
  legA -> SetFillColor(fColL);
  legA -> SetFillStyle(fFilL);
  legA -> SetTextFont(fTxt);
  legA -> SetTextAlign(fAlnL);
  legA -> AddEntry(hParticle, sLegPar.Data(), "pf");
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legA -> AddEntry(hAverage[iTrg], sLeg[iTrg].Data(), "pf");
  }

  TLegend *legS = new TLegend(xyLegS[0], xyLegS[1], xyLegS[2], xyLegS[3]);
  legS -> SetLineColor(fColL);
  legS -> SetLineStyle(fLinL);
  legS -> SetFillColor(fColL);
  legS -> SetFillStyle(fFilL);
  legS -> SetTextFont(fTxt);
  legS -> SetTextAlign(fAlnL);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    legS -> AddEntry(hAverage[iTrg], sLeg[iTrg].Data(), "pf");
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
  cout << "    Made text boxes." << endl;

  // make line
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  TLine *unity = new TLine(xPlot[0], 1., xPlot[1], 1.);
  unity -> SetLineColor(fColLi);
  unity -> SetLineStyle(fLinLi);
  unity -> SetLineWidth(fWidLi);

  // make plots
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const Float_t fMarginT(0.02);
  const Float_t fMarginR(0.02);
  const Float_t fMarginB(0.15);
  const Float_t fMarginL(0.15);

  TCanvas *cAverage = new TCanvas("cAverageVsParticle", "", width, height);
  cAverage  -> SetGrid(fGrid, fGrid);
  cAverage  -> SetTicks(fTick, fTick);
  cAverage  -> SetBorderMode(fMode);
  cAverage  -> SetBorderSize(fBord);
  cAverage  -> SetTopMargin(fMarginT);
  cAverage  -> SetRightMargin(fMarginR);
  cAverage  -> SetBottomMargin(fMarginB);
  cAverage  -> SetLeftMargin(fMarginL);
  cAverage  -> SetLogx(fLogX);
  cAverage  -> SetLogy(fLogY);
  cAverage  -> cd();
  hParticle -> Draw();
  //for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
  for (UInt_t iTrg = 0; iTrg < (NTrgBins - 1); iTrg++) {
    hAverage[iTrg] -> Draw("same");
  }
  legA     -> Draw();
  ptTxt    -> Draw();
  fOut     -> cd();
  cAverage -> Write();
  cAverage -> Close();

  TCanvas *cSmear = new TCanvas("cSmearingWeights", "", width, height);
  cSmear    -> SetGrid(fGrid, fGrid);
  cSmear    -> SetTicks(fTick, fTick);
  cSmear    -> SetBorderMode(fMode);
  cSmear    -> SetBorderSize(fBord);
  cSmear    -> SetTopMargin(fMarginT);
  cSmear    -> SetRightMargin(fMarginR);
  cSmear    -> SetBottomMargin(fMarginB);
  cSmear    -> SetLeftMargin(fMarginL);
  cSmear    -> SetLogx(fLogX);
  cSmear    -> SetLogy(fLogY);
  cSmear    -> cd();
  hSmear[0] -> Draw();
  //for (UInt_t iTrg = 1; iTrg < NTrgBins; iTrg++) {
  for (UInt_t iTrg = 1; iTrg < (NTrgBins - 1); iTrg++) {
    hSmear[iTrg] -> Draw("same");
  }
  unity  -> Draw();
  legS   -> Draw();
  ptTxt  -> Draw();
  fOut   -> cd();
  cSmear -> Write();
  cSmear -> Close();
  cout << "    Made plots." << endl;

  // save histograms
  fOut      -> cd();
  hMatrix   -> Write();
  hParticle -> Write();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hAverage[iTrg] -> Write();
    hSmear[iTrg]   -> Write();
    for (UInt_t iProjBin = 0; iProjBin < nProj[iTrg]; iProjBin++) {
      hProj[iTrg][iProjBin] -> Write();
    }
    if (smoothWeights) {
      for (UInt_t iSmooth = 0; iSmooth < nSmooth[iTrg]; iSmooth++) {
        fSmooth[iTrg][iSmooth] -> Write();
      }
    }
  }  // end trigger loop
  cout << "    Saved histograms." << endl;

  // close files
  fOut      -> cd();
  fOut      -> Close();
  fMatrix   -> cd();
  fMatrix   -> Close();
  fParticle -> cd();
  fParticle -> Close();
  cout << "  Finished projection script!\n" << endl;

}

// End ------------------------------------------------------------------------
