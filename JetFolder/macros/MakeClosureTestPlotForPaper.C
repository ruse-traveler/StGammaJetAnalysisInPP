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
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;

// global constants
const UInt_t NVtx(4);
const UInt_t NTxt(4);
const UInt_t NHist(2);
const UInt_t NRange(2);



void MakeClosureTestPlotForPaper() {

  // lower verbosity & set stat option
  gErrorIgnoreLevel = kError;
  gStyle            -> SetOptStat(0);
  cout << "\n  Plotting closure test for paper..." << endl;

  // io file parameters
  const TString sIn("closureTestFF.forThesis_pTbinHuge.et911r05pi0.d9m2y2022.root");
  const TString sOut("closureTestForPaper.ffWithRff_pTbinHuge.et911r05pi0.d5m3y2023.root");

  // io histogram parameters
  const TString sParticle("hParticleWithStat");
  const TString sNamePar("hParticleLevel");
  const TString sAverage[NHist] = {"hAverageLine",      "hVariationAverage"};
  const TString sRatio[NHist]   = {"hRatioAverageLine", "hRatioAverage"};
  const TString sNameAvg[NHist] = {"hAvgValue",         "hAvgError"};
  const TString sNameRat[NHist] = {"hRatioValue",       "hRatioError"};

  // style parameters
  const TString sTitleX("p_{T,jet}^{ch} [GeV/c]");
  const TString sTitleY("1/N_{trig} dN^{2}_{jets}/(dp_{T,jet}^{reco,ch} d#eta_{jet}) [GeV/c]^{-1}");
  const TString sTitleR("#frac{Corrected}{PYTHIA}");
  const UInt_t  fColPar(2);
  const UInt_t  fColOut(1);
  const UInt_t  fMarPar(1);
  const UInt_t  fMarAvg(1);
  const UInt_t  fMarRat(1);
  const UInt_t  fLinPar(1);
  const UInt_t  fLinOut(1);
  const UInt_t  fFilPar(0);
  const UInt_t  fFilOut(0);
  const UInt_t  fWidPar(2);
  const UInt_t  fWidOut(4);
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
  const Double_t rRange[NRange]      = {0.02,      2.8};

  // plot parameters
  const UInt_t   width(700);
  const UInt_t   height(800);
  const Double_t xyTxt[NVtx]   = {0.52, 0.69,  0.79,  0.96};
  const Double_t xyRatio[NVtx] = {0.02, 0.007, 0.917, 0.377};
  const TString  sTxt[NTxt]    = {"Closure Test", "anti-k_{T}, R=0.5", "9 < p_{T}^{trig} < 11 GeV/c", "#pi^{0}+jet"};
  const TString  sJetTxt("Corrected");
  const TString  sParTxt("PYTHIA");

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
  TH1D *hAverage[NHist];
  TH1D *hRatio[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] = (TH1D*) fIn -> Get(sAverage[iHist].Data());
    hRatio[iHist]   = (TH1D*) fIn -> Get(sRatio[iHist].Data());
    if (!hAverage[iHist] || !hRatio[iHist]) {
      cerr << "PANIC: couldn't grab an input histogram!\n"
           << "       hAverage[" << iHist << "] = " << hAverage[iHist] << ", hRatio[" << iHist << "] = " << hRatio[iHist] << "\n"
           << endl;
      return;
    }
    hAverage[iHist] -> SetName(sNameAvg[iHist].Data());
    hRatio[iHist]   -> SetName(sNameRat[iHist].Data());
  }

  TH1D *hParticle = (TH1D*) fIn -> Get(sParticle.Data());
  if (!hParticle) {
    cerr << "PANIC: couldn't gran input particle histogram!\n" << endl;
    return;
  }
  hParticle -> SetName(sNamePar.Data());
  cout << "    Grabbed histograms." << endl;

  // set styles
  hParticle   -> SetLineColor(fColPar);
  hParticle   -> SetLineStyle(fLinPar);
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

  TH1D *hRatOutline = (TH1D*) hRatio[1]   -> Clone();
  hRatOutline -> SetName("hRatioOutline");
  hRatOutline -> SetLineColor(fColOut);
  hRatOutline -> SetLineStyle(fLinOut);
  hRatOutline -> SetLineWidth(fWidOut);
  hRatOutline -> SetFillColor(fColOut);
  hRatOutline -> SetFillStyle(fFilOut);
  cout << "    Made outline histograms." << endl;

  // make legend histogram
  hAvgLegend = (TH1D*) hAverage[1] -> Clone();
  hAvgLegend -> SetName("hAverageLegend");
  hAvgLegend -> SetLineColor(fColOut);
  cout << "    Made legend histogram." << endl;

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

  TLegend *legTxt = new TLegend(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3]);
  legTxt -> SetTextSize(0.036);  
  legTxt -> SetFillColor(0);   
  legTxt -> SetLineColor(0);
  legTxt -> AddEntry((TObject*) 0, sTxt[0].Data(), "");
  legTxt -> AddEntry((TObject*) 0, sTxt[1].Data(), "");
  legTxt -> AddEntry((TObject*) 0, sTxt[2].Data(), "");
  legTxt -> AddEntry((TObject*) 0, sTxt[3].Data(), "");
  legTxt -> AddEntry(hAvgLegend,   sJetTxt.Data(), "f");
  legTxt -> AddEntry(hParticle,    sParTxt.Data(), "l");
  cout << "    Made text boxes." << endl;

  // make line
  TLine *lUnity = new TLine(xBinPtRange[0], 1., xBinPtRange[1], 1.);
  lUnity -> SetLineColor(1);
  lUnity -> SetLineStyle(2);
  cout << "    Made line." << endl;

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
  hAvgOutline -> Draw("E5 SAME");
  hAverage[1] -> Draw("E5 SAME");
  hAverage[0] -> Draw("][ L HIST SAME");
  hParticle   -> Draw("][ L HIST SAME");
  legTxt      -> Draw("SAME");

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
  pRatio      -> SetBottomMargin(0.3345725);
  pRatio      -> SetFrameBorderMode(0);
  pRatio      -> SetFrameBorderMode(0);
  hFrameRatio -> DrawCopy("9");
  hRatOutline -> Draw("E5 SAME");
  hRatio[1]   -> Draw("E5 SAME");
  hRatio[0]   -> Draw("][ L HIST SAME");
  lUnity      -> Draw();
  fOut        -> cd();
  cPlot       -> Write();
  cPlot       -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> Write();
    hRatio[iHist]   -> Write();
  }
  hParticle   -> Write();
  hAvgOutline -> Write();
  hRatOutline -> Write();
  hAvgLegend  -> Write();
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
