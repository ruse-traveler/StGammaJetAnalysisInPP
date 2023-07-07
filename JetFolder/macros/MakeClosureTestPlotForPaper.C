// MakeClosureTestPlotForPaper.C'
// Derek Anderson
// 03.03.2023
//
// Use this to plot the output of
// 'MakeClosureTestPlot.C' in
// the style of the long paper.
//
// Adapted from code by Nihar Sahoo.

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

using namespace std;

// global constants
static const UInt_t NVtx(4);
static const UInt_t NTxt(5);
static const UInt_t NVals(2);
static const UInt_t NHist(2);
static const UInt_t NRange(2);
static const UInt_t NBinMax(20);



void MakeClosureTestPlotForPaper() {

  // lower verbosity & set stat option
  gErrorIgnoreLevel = kError;
  gStyle            -> SetOptStat(0);
  cout << "\n  Plotting closure test for paper..." << endl;

  // io file parameters
  const TString sIn("closure/et911r05ff_rebinClosure/closureTestFF.forThesis_pTbinHuge.et911r05pi0.d9m2y2022.root");
  const TString sOut("closureTestForPaper_withDataErrorsAndLowPtCutoff_mayCollabComments.ffWithRff_pTbinHuge.et911r05pi0.d7m7y2023.root");

  // io histogram parameters
  const TString sParticle("hParticleWithStat");
  const TString sParUnity("hUnityWithStat");
  const TString sAverage[NHist] = {"hAverageLine",      "hVariationAverage"};
  const TString sRatio[NHist]   = {"hRatioAverageLine", "hRatioAverage"};

  // output histogram parameters
  const TString sNamePar("hParticleLevel");
  const TString sNameUnity("hUnityWithParErrors");
  const TString sNameParError("hParWithErrors");
  const TString sNameAvg[NHist] = {"hAvgValue",   "hAvgError"};
  const TString sNameRat[NHist] = {"hRatioValue", "hRatioError"};

  // style parameters
  const TString sTitleX("p_{T,jet}^{ch} [GeV/#it{c}]");
  const TString sTitleY("1/N_{trig} dN^{2}_{jets}/(dp_{T,jet}^{reco,ch} d#eta_{jet}) [GeV/#it{c}]^{-1}");
  const TString sTitleR("#frac{Corrected}{PYTHIA}");
  const UInt_t  fColPar(635);
  const UInt_t  fColErr(622);
  const UInt_t  fColUni(622);
  const UInt_t  fColOut(1);
  const UInt_t  fMarPar(1);
  const UInt_t  fMarErr(1);
  const UInt_t  fMarUni(1);
  const UInt_t  fMarAvg(1);
  const UInt_t  fMarRat(1);
  const UInt_t  fLinPar(1);
  const UInt_t  fLinErr(1);
  const UInt_t  fLinUni(1);
  const UInt_t  fLinOut(1);
  const UInt_t  fFilPar(0);
  const UInt_t  fFilErr(1001);
  const UInt_t  fFilUni(1001);
  const UInt_t  fFilOut(0);
  const UInt_t  fWidPar(2);
  const UInt_t  fWidErr(1);
  const UInt_t  fWidUni(1);
  const UInt_t  fWidOut(6);
  const UInt_t  fColAvg[NHist] = {1, 16};
  const UInt_t  fColRat[NHist] = {1, 16};
  const UInt_t  fLinAvg[NHist] = {2, 1};
  const UInt_t  fLinRat[NHist] = {2, 1};
  const UInt_t  fWidAvg[NHist] = {2, 1};
  const UInt_t  fWidRat[NHist] = {2, 1};
  const UInt_t  fFilAvg[NHist] = {0, 1001};
  const UInt_t  fFilRat[NHist] = {0, 1001};

  // frame parameters
  const Int_t    xNumBins(-10);
  const Int_t    yNumBins(1000);
  const Int_t    xBinPtRange[NRange] = {0,         29};
  const Double_t yRange[NRange]      = {0.0000004, 3.};
  const Double_t rRange[NRange]      = {0.02,      1.8};

  // plot parameters
  const UInt_t   width(700);
  const UInt_t   height(800);
  const Double_t xyTxt[NVtx]   = {0.52, 0.69,  0.79,  0.96};
  const Double_t xyRatio[NVtx] = {0.02, 0.007, 0.917, 0.377};
  const TString  sTxt[NTxt]    = {"#bf{STAR}", "Closure Test", "anti-k_{T}, R=0.5", "9 < E_{T}^{trig} < 11 GeV", "#pi^{0}+jet"};
  const TString  sJetTxt("Corrected");
  const TString  sParTxt("PYTHIA");

  // distribution parameters
  const Bool_t  doLowPtCutoff(true);
  const Bool_t  showTruthError(true);
  const Float_t minJetPt(3.);
  const Float_t maxJetPt(xBinPtRange[1]);

  // open files
  TFile *fIn  = new TFile(sIn.Data(),  "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut || !fIn) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fIn = " << fIn
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input
  TH1D *hInAverage[NHist];
  TH1D *hInRatio[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInAverage[iHist] = (TH1D*) fIn -> Get(sAverage[iHist].Data());
    hInRatio[iHist]   = (TH1D*) fIn -> Get(sRatio[iHist].Data());
    if (!hInAverage[iHist] || !hInRatio[iHist]) {
      cerr << "PANIC: couldn't grab an input histogram!\n"
           << "       hAverage[" << iHist << "] = " << hInAverage[iHist] << ", hRatio[" << iHist << "] = " << hInRatio[iHist] << "\n"
           << endl;
      return;
    }
  }

  TH1D *hInParticle = (TH1D*) fIn -> Get(sParticle.Data());
  if (!hInParticle) {
    cerr << "PANIC: couldn't grab input particle histogram!\n" << endl;
    return;
  }

  TH1D *hInUnity;
  if (showTruthError) {
    hInUnity = (TH1D*) fIn -> Get(sParUnity.Data());
    if (!hInUnity) {
      cerr << "PANIC: couldn't grab input unity histogram!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // copy input to output histograms
  TH1D *hUnity;
  TH1D *hParticle;
  TH1D *hParError;
  TH1D *hAverage[NHist];
  TH1D *hRatio[NHist];

  hParticle = (TH1D*) hInParticle -> Clone();
  if (showTruthError) {
    hParError = (TH1D*) hInParticle -> Clone();
    hUnity    = (TH1D*) hInUnity    -> Clone();
  }
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] = (TH1D*) hInAverage[iHist] -> Clone();
    hRatio[iHist]   = (TH1D*) hInRatio[iHist]   -> Clone();
  }

  UInt_t   nBinToDrawPar;
  UInt_t   nBinToDrawAvg[NHist];
  UInt_t   nBinToDrawRat[NHist];
  Double_t xyPar[NVals][NBinMax];
  Double_t xyErr[NVals][NBinMax];
  Double_t xyAvg[NHist][NVals][NBinMax];
  Double_t xyRat[NHist][NVals][NBinMax];

  // get no. of bins
  const UInt_t nBinsPar = hParticle   -> GetNbinsX();
  const UInt_t nBinsAvg = hAverage[0] -> GetNbinsX();
  const UInt_t nBinsRat = hRatio[0]   -> GetNbinsX();

  UInt_t nBinsUni = 0;
  if (showTruthError) {
    nBinsUni = hUnity -> GetNbinsX();
  }

  nBinToDrawPar = 0;
  for (UInt_t iBin = 0; iBin < NBinMax; iBin++) {
    xyPar[0][iBin] = 0.;
    xyPar[1][iBin] = 0.;
    xyErr[0][iBin] = 0.;
    xyErr[1][iBin] = 0.;
  }
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    nBinToDrawAvg[iHist] = 0;
    nBinToDrawRat[iHist] = 0;
    for (UInt_t iBin = 0; iBin < NBinMax; iBin++) {
      xyAvg[iHist][0][iBin] = 0.;
      xyAvg[iHist][1][iBin] = 0.;
      xyRat[iHist][0][iBin] = 0.;
      xyRat[iHist][1][iBin] = 0.;
    }
  }

  // apply low pt cutoff (if needed)
  if (doLowPtCutoff) {

    // apply cutoff to particle spectrum
    hParticle -> Reset("ICES");
    for (UInt_t iBinPar = 1; iBinPar < (nBinsPar + 1); iBinPar++) {
      const Double_t binWidthPar  = hInParticle -> GetBinWidth(iBinPar);
      const Double_t binCenterPar = hInParticle -> GetBinCenter(iBinPar);
      const Double_t binValuePar  = hInParticle -> GetBinContent(iBinPar);
      const Double_t binErrorPar  = hInParticle -> GetBinError(iBinPar);
      if ((binCenterPar >= minJetPt) && (binCenterPar < maxJetPt)) {

        // set histogram values
        hParticle -> SetBinContent(iBinPar, binValuePar);
        hParticle -> SetBinError(iBinPar,   binErrorPar);

        // store values
        xyPar[0][nBinToDrawPar] = binCenterPar;
        xyPar[1][nBinToDrawPar] = binValuePar;
        xyErr[0][nBinToDrawPar] = binWidthPar;
        xyErr[1][nBinToDrawPar] = binErrorPar;
        ++nBinToDrawPar;
      }
    }  // end particle bin loop

    // apply cutoff to unity histogram
    if (showTruthError) {

      hUnity -> Reset("ICES");
      for (UInt_t iBinUni = 1; iBinUni < (nBinsUni + 1); iBinUni++) {
        const Double_t binWidthUni  = hInUnity -> GetBinWidth(iBinUni);
        const Double_t binCenterUni = hInUnity -> GetBinCenter(iBinUni);
        const Double_t binValueUni  = hInUnity -> GetBinContent(iBinUni);
        const Double_t binErrorUni  = hInUnity -> GetBinError(iBinUni);
        if ((binCenterUni >= minJetPt) && (binCenterUni < maxJetPt)) {
          hUnity -> SetBinContent(iBinUni, binValueUni);
          hUnity -> SetBinError(iBinUni,   binErrorUni);
        }
      }  // end unity bin loop
    }

    // apply cutoff to averages and ratios
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      hAverage[iHist] -> Reset("ICES");
      hRatio[iHist]   -> Reset("ICES");
      for (UInt_t iBinAvg = 1; iBinAvg < (nBinsAvg + 1); iBinAvg++) {
        const Double_t binCenterAvg = hInAverage[iHist] -> GetBinCenter(iBinAvg);
        const Double_t binValueAvg  = hInAverage[iHist] -> GetBinContent(iBinAvg);
        const Double_t binErrorAvg  = hInAverage[iHist] -> GetBinError(iBinAvg);
        if ((binCenterAvg >= minJetPt) && (binCenterAvg < maxJetPt)) {

          // set histogram values
          hAverage[iHist] -> SetBinContent(iBinAvg, binValueAvg);
          hAverage[iHist] -> SetBinError(iBinAvg,   binErrorAvg);

          // store values
          xyAvg[iHist][0][nBinToDrawAvg[iHist]] = binCenterAvg;
          xyAvg[iHist][1][nBinToDrawAvg[iHist]] = binValueAvg;
          ++nBinToDrawAvg[iHist];
        }
      }  // end average bin loop
      for (UInt_t iBinRat = 1; iBinRat < (nBinsRat + 1); iBinRat++) {
        const Double_t binCenterRat = hInRatio[iHist] -> GetBinCenter(iBinRat);
        const Double_t binValueRat  = hInRatio[iHist] -> GetBinContent(iBinRat);
        const Double_t binErrorRat  = hInRatio[iHist] -> GetBinError(iBinRat);
        if ((binCenterRat >= minJetPt) && (binCenterRat < maxJetPt)) {

          // set histogram values
          hRatio[iHist] -> SetBinContent(iBinRat, binValueRat);
          hRatio[iHist] -> SetBinError(iBinRat,   binErrorRat);

          // store values
          xyRat[iHist][0][nBinToDrawRat[iHist]] = binCenterRat;
          xyRat[iHist][1][nBinToDrawRat[iHist]] = binValueRat;
          ++nBinToDrawRat[iHist];
        }
      }  // end average bin loop
    }  // end histogram loop
    cout << "    Cut off spectra at low pt." << endl;
  }  // end low pt cutoff

  // set histogram names
  hParticle -> SetName(sNamePar.Data());
  if (showTruthError) {
    hParError -> SetName(sNameParError.Data());
    hUnity    -> SetName(sNameUnity.Data());
  }
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> SetName(sNameAvg[iHist].Data());
    hRatio[iHist]   -> SetName(sNameRat[iHist].Data());
  }

  // set styles
  hParticle -> SetMarkerColor(fColPar);
  hParticle -> SetMarkerStyle(fMarPar);
  hParticle -> SetLineColor(fColPar);
  hParticle -> SetLineStyle(fLinPar);
  hParticle -> SetLineWidth(fWidPar);
  hParticle -> SetFillColor(fColPar);
  hParticle -> SetFillStyle(fFilPar);
  if (showTruthError) {
    hParError -> SetMarkerColor(fColErr);
    hParError -> SetMarkerStyle(fMarErr);
    hParError -> SetLineColor(fColErr);
    hParError -> SetLineStyle(fLinErr);
    hParError -> SetLineWidth(fWidErr);
    hParError -> SetFillColor(fColErr);
    hParError -> SetFillStyle(fFilErr);
    hUnity    -> SetMarkerColor(fColUni);
    hUnity    -> SetMarkerStyle(fMarUni);
    hUnity    -> SetLineColor(fColUni);
    hUnity    -> SetLineStyle(fLinUni);
    hUnity    -> SetLineWidth(fWidUni);
    hUnity    -> SetFillColor(fColUni);
    hUnity    -> SetFillStyle(fFilUni);
  }
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> SetMarkerColor(fColAvg[iHist]);
    hAverage[iHist] -> SetMarkerStyle(fMarAvg);
    hAverage[iHist] -> SetLineColor(fColAvg[iHist]);
    hAverage[iHist] -> SetLineStyle(fLinAvg[iHist]);
    hAverage[iHist] -> SetLineWidth(fWidAvg[iHist]);
    hAverage[iHist] -> SetFillColor(fColAvg[iHist]);
    hAverage[iHist] -> SetFillStyle(fFilAvg[iHist]);
    hRatio[iHist]   -> SetMarkerColor(fColRat[iHist]);
    hRatio[iHist]   -> SetMarkerStyle(fMarRat);
    hRatio[iHist]   -> SetLineColor(fColRat[iHist]);
    hRatio[iHist]   -> SetLineStyle(fLinRat[iHist]);
    hRatio[iHist]   -> SetLineWidth(fWidRat[iHist]);
    hRatio[iHist]   -> SetFillColor(fColRat[iHist]);
    hRatio[iHist]   -> SetFillStyle(fFilRat[iHist]);
  }
  cout << "    Set styles." << endl;

  // create outline histograms
  TH1D *hAvgOutline = (TH1D*) hAverage[1] -> Clone();
  hAvgOutline -> SetName("hAverageOutline");
  hAvgOutline -> SetLineColor(fColOut);
  hAvgOutline -> SetLineStyle(fLinOut);
  hAvgOutline -> SetLineWidth(fWidOut);
  hAvgOutline -> SetFillColor(fColOut);
  hAvgOutline -> SetFillStyle(fFilOut);

  TH1D *hRatOutline = (TH1D*) hRatio[1] -> Clone();
  hRatOutline -> SetName("hRatioOutline");
  hRatOutline -> SetLineColor(fColOut);
  hRatOutline -> SetLineStyle(fLinOut);
  hRatOutline -> SetLineWidth(fWidOut);
  hRatOutline -> SetFillColor(fColOut);
  hRatOutline -> SetFillStyle(fFilOut);

  TH1D *hParOutline = (TH1D*) hParticle -> Clone();
  hParOutline -> SetName("hParticleOutline");
  hParOutline -> SetLineColor(fColPar);
  hParOutline -> SetLineStyle(fLinOut);
  hParOutline -> SetLineWidth(fWidOut);
  hParOutline -> SetFillColor(fColPar);
  hParOutline -> SetFillStyle(fFilOut);
 
  TH1D *hUniOutline = (TH1D*) hUnity -> Clone();
  hUniOutline -> SetName("hUnityOutline");
  hUniOutline -> SetLineColor(fColPar);
  hUniOutline -> SetLineStyle(fLinOut);
  hUniOutline -> SetLineWidth(fWidOut);
  hUniOutline -> SetFillColor(fColPar);
  hUniOutline -> SetFillStyle(fFilOut);
  cout << "    Made outline histograms." << endl;

  // make legend histograms
  TH1D *hAvgLegend = (TH1D*) hAverage[1] -> Clone();
  hAvgLegend -> SetName("hAverageLegend");
  hAvgLegend -> SetLineColor(fColOut);

  TH1D *hParLegend = (TH1D*) hParticle -> Clone();
  hParLegend -> SetName("hParticleLegend");
  hParLegend -> SetLineColor(fColPar);
  hParLegend -> SetLineWidth(fWidPar);
  if (showTruthError) {
    hParLegend -> SetFillColor(fColErr);
    hParLegend -> SetFillStyle(fFilErr);
  } else {
    hParLegend -> SetFillColor(fColPar);
    hParLegend -> SetFillStyle(fFilPar);
  }
  cout << "    Made legend histograms." << endl;

  // make frame histograms
  TH2D *hFrameJets = new TH2D("hFrameJets", "" , xNumBins, xBinPtRange[0], xBinPtRange[1], yNumBins, yRange[0], yRange[1]);
  hFrameJets -> GetXaxis() -> SetLabelFont(42);
  hFrameJets -> GetXaxis() -> SetLabelSize(0.0);
  hFrameJets -> GetXaxis() -> SetTitleSize(0.033);
  hFrameJets -> GetXaxis() -> SetTitleFont(42);
  hFrameJets -> GetYaxis() -> SetTitle(sTitleY.Data());
  hFrameJets -> GetYaxis() -> CenterTitle(true);
  hFrameJets -> GetYaxis() -> SetLabelFont(42);
  hFrameJets -> GetYaxis() -> SetLabelSize(0.05);
  hFrameJets -> GetYaxis() -> SetTitleSize(0.04);
  hFrameJets -> GetYaxis() -> SetTitleOffset(1.93);
  hFrameJets -> GetYaxis() -> SetTitleFont(42);
  hFrameJets -> GetZaxis() -> SetLabelFont(42);
  hFrameJets -> GetZaxis() -> SetLabelSize(0.035);
  hFrameJets -> GetZaxis() -> SetTitleSize(0.035);
  hFrameJets -> GetZaxis() -> SetTitleFont(42);

  TH2F *hFrameRatio = new TH2F("hFrameRatio", "", xNumBins, xBinPtRange[0], xBinPtRange[1], yNumBins, rRange[0], rRange[1]);
  hFrameRatio -> GetXaxis() -> SetTitle(sTitleX.Data());
  hFrameRatio -> GetXaxis() -> CenterTitle(true);
  hFrameRatio -> GetXaxis() -> SetLabelFont(42);
  hFrameRatio -> GetXaxis() -> SetLabelOffset(0.001);
  hFrameRatio -> GetXaxis() -> SetLabelSize(0.14);
  hFrameRatio -> GetXaxis() -> SetTitleSize(0.14);
  hFrameRatio -> GetXaxis() -> SetTitleFont(42);
  hFrameRatio -> GetXaxis() -> SetTickLength(0.06);
  hFrameRatio -> GetYaxis() -> SetTickLength(0.04);
  hFrameRatio -> GetYaxis() -> SetTitle(sTitleR.Data());
  hFrameRatio -> GetYaxis() -> CenterTitle(true);
  hFrameRatio -> GetYaxis() -> SetLabelFont(42);
  hFrameRatio -> GetYaxis() -> SetLabelSize(0.12);
  hFrameRatio -> GetYaxis() -> SetTitleSize(0.1);
  hFrameRatio -> GetYaxis() -> SetNdivisions(505);
  hFrameRatio -> GetYaxis() -> SetTitleOffset(0.81);
  hFrameRatio -> GetYaxis() -> SetTitleFont(42);
  hFrameRatio -> GetZaxis() -> SetLabelFont(42);
  hFrameRatio -> GetZaxis() -> SetLabelSize(0.035);
  hFrameRatio -> GetZaxis() -> SetTitleSize(0.035);
  hFrameRatio -> GetZaxis() -> SetTitleFont(42);
  cout << "    Made frame histograms." << endl;

  // make text box
  TLegend *legTxt = new TLegend(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3]);
  legTxt -> SetTextFont(42);
  legTxt -> SetTextSize(0.036);  
  legTxt -> SetFillColor(0);   
  legTxt -> SetLineColor(0);
  legTxt -> AddEntry((TObject*) 0, sTxt[0].Data(), "");
  legTxt -> AddEntry((TObject*) 0, sTxt[1].Data(), "");
  legTxt -> AddEntry((TObject*) 0, sTxt[2].Data(), "");
  legTxt -> AddEntry((TObject*) 0, sTxt[3].Data(), "");
  legTxt -> AddEntry(hAvgLegend,   sJetTxt.Data(), "f");
  if (showTruthError) {
    legTxt -> AddEntry(hParLegend, sParTxt.Data(), "lf");
  } else {
    legTxt -> AddEntry(hParLegend, sParTxt.Data(), "l");
  }
  cout << "    Made text boxes." << endl;

  // make lines
  TLine *lUnity = new TLine(xBinPtRange[0], 1., xBinPtRange[1], 1.);
  lUnity -> SetLineColor(1);
  lUnity -> SetLineStyle(2);

  TGraph       *gParticle = new TGraph(nBinToDrawPar, xyPar[0], xyPar[1]);
  TGraphErrors *gParError = new TGraphErrors(nBinToDrawPar, xyPar[0], xyPar[1], xyErr[0], xyErr[1]);
  gParticle -> SetName("gParticle");
  gParticle -> SetMarkerColor(fColPar);
  gParticle -> SetMarkerStyle(fMarPar);
  gParticle -> SetLineColor(fColPar);
  gParticle -> SetLineStyle(fLinPar);
  gParticle -> SetLineWidth(fWidPar);
  gParticle -> SetFillColor(fColPar);
  gParticle -> SetFillStyle(fFilPar);
  gParError -> SetName("gParError");
  gParError -> SetMarkerColor(fColErr);
  gParError -> SetMarkerStyle(fMarErr);
  gParError -> SetLineColor(fColErr);
  gParError -> SetLineStyle(fLinErr);
  gParError -> SetLineWidth(fWidErr);
  gParError -> SetFillColor(fColErr);
  gParError -> SetFillStyle(fFilErr);

  TGraph *gAverage = new TGraph(nBinToDrawAvg[0], xyAvg[0][0], xyAvg[0][1]);
  TGraph *gRatio   = new TGraph(nBinToDrawRat[0], xyRat[0][0], xyRat[0][1]);
  gAverage -> SetName("gAverage");
  gAverage -> SetLineColor(fColAvg[0]);
  gAverage -> SetLineStyle(fLinAvg[0]);
  gAverage -> SetLineWidth(fWidAvg[0]);
  gRatio   -> SetName("gRatio");
  gRatio   -> SetLineColor(fColRat[0]);
  gRatio   -> SetLineStyle(fLinRat[0]);
  gRatio   -> SetLineWidth(fWidRat[0]);
  cout << "    Made lines." << endl;

  // make canvas
  TCanvas *cPlot = new TCanvas("cPlot", "canvas", width, height);
  cPlot       -> SetFillColor(0);
  cPlot       -> SetBorderMode(0);
  cPlot       -> SetBorderSize(2);
  cPlot       -> SetLogy();
  cPlot       -> SetTickx(1);
  cPlot       -> SetTicky(1);
  cPlot       -> SetLeftMargin(0.1794479);
  cPlot       -> SetRightMargin(0.0993865);
  cPlot       -> SetTopMargin(0.01536313);
  cPlot       -> SetBottomMargin(0.3743017);
  cPlot       -> SetFrameBorderMode(0);
  cPlot       -> SetFrameBorderMode(0);
  cPlot       -> cd();
  hFrameJets  -> DrawCopy("9");
  if (showTruthError) {
    if (doLowPtCutoff) {
      gParError   -> Draw("3");
      hParOutline -> Draw("E5 SAME");
      gParticle   -> Draw("L");
      hAverage[1] -> Draw("E5 SAME");
      hAvgOutline -> Draw("E5 SAME");
      gAverage    -> Draw("L");
    } else {
      hParError   -> Draw("E5 same");
      hParOutline -> Draw("E5 same");
      hParticle   -> Draw("][ L HIST SAME");
      hAverage[1] -> Draw("E5 SAME");
      hAvgOutline -> Draw("E5 SAME");
      hAverage[0] -> Draw("][ L HIST SAME");
    }
  } else {
    if (doLowPtCutoff) {
      hAverage[1]  -> Draw("E5 SAME");
      hAvgOutline  -> Draw("E5 SAME");
      gAverage     -> Draw("L");
      gParticle    -> Draw("L");
    } else {
      hAverage[1] -> Draw("E5 SAME");
      hAvgOutline -> Draw("E5 SAME");
      hAverage[0] -> Draw("][ L HIST SAME");
      hParticle   -> Draw("][ L HIST SAME");
    }
  }
  legTxt -> Draw("SAME");

  // create pad for ratio
  TPad *pRatio = new TPad("pRatio", "pRatio", xyRatio[0], xyRatio[1], xyRatio[2], xyRatio[3]);
  cPlot       -> cd();
  pRatio      -> Draw();
  pRatio      -> cd();  
  pRatio      -> SetFillColor(0);
  pRatio      -> SetBorderMode(0);
  pRatio      -> SetBorderSize(2);
  pRatio      -> SetTicky(1);
  pRatio      -> SetLeftMargin(0.1784615);
  pRatio      -> SetRightMargin(0.01692308);
  pRatio      -> SetTopMargin(0.007434944);
  pRatio      -> SetBottomMargin(0.3945725);
  pRatio      -> SetFrameBorderMode(0);
  pRatio      -> SetFrameBorderMode(0);
  hFrameRatio -> DrawCopy("9");
  if (showTruthError) {
    if (doLowPtCutoff) {
      hUnity      -> Draw("E5 SAME");
      hUniOutline -> Draw("E5 SAME");
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      gRatio      -> Draw("L");
    } else {
      hUnity      -> Draw("E5 SAME");
      hUniOutline -> Draw("E5 SAME");
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      hRatio[0]   -> Draw("][ L HIST SAME");
    }
  } else {
    if (doLowPtCutoff) {
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      gRatio      -> Draw("L");
    } else {
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      hRatio[0]   -> Draw("][ L HIST SAME");
    }
  }
  lUnity -> Draw();
  fOut   -> cd();
  cPlot  -> Write();
  cPlot  -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> Write();
    hRatio[iHist]   -> Write();
  }
  hParticle -> Write();
  if (showTruthError) {
    hParError   -> Write();
    hUnity      -> Write();
    hParOutline -> Write();
    hUniOutline -> Write();
  }
  if (doLowPtCutoff) {
    gParticle -> Write();
    gAverage  -> Write();
    gRatio    -> Write();
    if (showTruthError) {
      gParError -> Write();
    }
  }
  hAvgOutline -> Write();
  hRatOutline -> Write();
  hAvgLegend  -> Write();
  hParLegend  -> Write();
  hFrameJets  -> Write();
  hFrameRatio -> Write();

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished plotting closure test!\n" << endl;

}

// End ------------------------------------------------------------------------
