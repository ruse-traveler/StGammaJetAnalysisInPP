// 'QuickEfficiencyPlot.C'
// Derek Anderson
// 09.29.2020
//
// Grabs an efficiency histogram and plots it with the systematic
// uncertainty and a fit function.  Also produces a "smoothed"
// efficiency based on the fit function.

#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// io parameters
static const TString sIn("input/trackEffVsPseduoEff.onlyPtParGT0.et920vz55pi0.d2m3y2020.root");
static const TString sOut("smoothedEffWithPlot_pseudoEffForedToConvergeToEff_allOtherParametersFixed.et920vz55pi0.d6m2y2021.root");
static const TString sHist("hPseudoEff_AllVtx");

// plot parameters
static const Float_t eSys(0.04);
static const Float_t xPlotMin(0.);
static const Float_t xPlotMax(20.);
static const TString sName("hEff_forceEffAndPseudoEffToConverge");
static const TString sSmooth("hEffSmooth_forceEffAndPseudoEffToConverge");
static const TString sPlotName("cPseudoEff");
static const TString sPlotSmooth("cPseudoEffSmooth");

// fit parameters
static const UInt_t   nPar(5);
static const Float_t  fFitMin(0.2);
static const Float_t  fFitMax(20.);
static const TString  sFunc("[0]+[1]*TMath::Exp(-1.*[2]*x)+[3]*TMath::Exp(-1.*[4]*x*x)");
static const TString  sFuncName("fEff_forceEffAndPseudoEffToConverge");
static const TString  sPar[nPar]   = {"#epsilon_{0}", "#epsilon_{1}", "#sigma_{1}", "#epsilon_{2}", "#sigma_{2}"};
static const Double_t pGuess[nPar] = {0.82, -1.81, 13.0, 0.12, 0.67};



void QuickEfficiencyPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning quick efficiency plot..." << endl;


  // style parameters
  const UInt_t  fColE(1);
  const UInt_t  fColF(879);
  const UInt_t  fMar(8);
  const UInt_t  fLin(1);
  const UInt_t  fWid(2);
  const UInt_t  fTxt(42);
  const UInt_t  fFil(0);
  const UInt_t  fCnt(1);
  const Float_t fLbl(0.03);
  const Float_t fSiz(0.04);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.2);
  const TString sTitle("");
  const TString sTitleX("p_{T}^{trk} [GeV/c]");
  const TString sTitleY("#tilde{#epsilon}_{trk}");

  // canvas parameters
  const UInt_t  fColL(1);
  const UInt_t  fLinPlot(2);
  const UInt_t  fGridX(0);
  const UInt_t  fGridY(0);
  const UInt_t  width(750);
  const UInt_t  height(750);
  const Float_t lValueY(1.);
  const Float_t tMargin(0.02);
  const Float_t rMargin(0.02);
  const Float_t bMargin(1.28);
  const Float_t lMargin(1.28);


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn  = new TFile(sIn.Data(), "read");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // get histogram
  TH1D *hEff = (TH1D*) fIn -> Get(sHist.Data());
  if (!hEff) {
    cerr << "PANIC: couldn't grab input histogram!" << endl;
    return;
  }
  cout << "    Grabbed histogram." << endl;


  // set style
  hEff -> SetName(sName.Data());
  hEff -> SetTitle(sTitle.Data());
  hEff -> SetTitleFont(fTxt);
  hEff -> SetLineColor(fColE);
  hEff -> SetLineStyle(fLin);
  hEff -> SetFillColor(fColE);
  hEff -> SetFillStyle(fFil);
  hEff -> SetMarkerColor(fColE);
  hEff -> SetMarkerStyle(fMar);
  hEff -> GetXaxis() -> SetLabelFont(fTxt);
  hEff -> GetXaxis() -> SetLabelSize(fLbl);
  hEff -> GetXaxis() -> SetTitle(sTitleX.Data());
  hEff -> GetXaxis() -> SetTitleFont(fTxt);
  hEff -> GetXaxis() -> SetTitleSize(fSiz);
  hEff -> GetXaxis() -> SetTitleOffset(fOffX);
  hEff -> GetXaxis() -> CenterTitle(fCnt);
  hEff -> GetYaxis() -> SetLabelFont(fTxt);
  hEff -> GetYaxis() -> SetLabelSize(fLbl);
  hEff -> GetYaxis() -> SetTitle(sTitleY.Data());
  hEff -> GetYaxis() -> SetTitleFont(fTxt);
  hEff -> GetYaxis() -> SetTitleSize(fSiz);
  hEff -> GetYaxis() -> SetTitleOffset(fOffY);
  hEff -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set histogram styles." << endl;


  // fit efficiency
  TF1 *fEff = new TF1(sFuncName.Data(), sFunc.Data(), fFitMin, fFitMax);
  for (UInt_t iPar = 0; iPar < nPar; iPar++) {
    fEff -> SetParName(iPar, sPar[iPar]);
    fEff -> SetParameter(iPar, pGuess[iPar]);
  }
  // TEST [02.06.2021]
  fEff -> FixParameter(0, 0.82);
  //fEff -> FixParameter(1, -1.81);
  fEff -> FixParameter(2, 13.00);
  //fEff -> FixParameter(3, 0.12);
  fEff -> FixParameter(4, 0.67);
  fEff -> SetLineColor(fColF);
  fEff -> SetLineStyle(fLin);
  fEff -> SetLineWidth(fWid);
  hEff -> Fit(fEff, "R");
  hEff -> GetXaxis() -> SetRangeUser(xPlotMin, xPlotMax);
  cout << "    Fit efficiency." << endl;


  // set systematic uncertainties
  const UInt_t nBins = hEff -> GetNbinsX();
  for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
    hEff -> SetBinError(iBin, eSys);
  }
  cout << "    Set systematic uncertainties." << endl;

  // create smoothed histogram
  TH1D *hSmooth = (TH1D*) hEff -> Clone();
  hSmooth -> SetName(sSmooth.Data());
  hSmooth -> Reset("ICE");
  for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
    const Float_t xVal = hSmooth -> GetBinCenter(iBin);
    const Float_t yVal = fEff    -> Eval(xVal);
    hSmooth -> SetBinContent(iBin, yVal);
    hSmooth -> SetBinError(iBin, eSys);
  }
  cout << "    Created smoothed efficiency." << endl;


  // create line
  TLine *line = new TLine(xPlotMin, lValueY, xPlotMax, lValueY);
  line -> SetLineStyle(fLinPlot);
  line -> SetLineColor(fColL);

  // draw plots
  TCanvas *cEff = new TCanvas(sPlotName.Data(), "", width, height);
  cEff -> SetGrid(fGridX, fGridY);
  cEff -> SetTopMargin(tMargin);
  cEff -> SetRightMargin(rMargin);
  cEff -> SetBottomMargin(bMargin);
  cEff -> SetLeftMargin(lMargin);
  cEff -> cd();
  hEff -> Draw("E5");
  line -> Draw();
  fOut -> cd();
  cEff -> Write();
  cEff -> Close();

  TCanvas *cSmooth = new TCanvas(sPlotSmooth.Data(), "", width, height);
  cSmooth -> SetGrid(fGridX, fGridY);
  cSmooth -> SetTopMargin(tMargin);
  cSmooth -> SetRightMargin(rMargin);
  cSmooth -> SetBottomMargin(bMargin);
  cSmooth -> SetLeftMargin(lMargin);
  cSmooth -> cd();
  hSmooth -> Draw("E5");
  line    -> Draw();
  fOut    -> cd();
  cSmooth -> Write();
  cSmooth -> Close();
  cout << "    Made plots." << endl; 


  // save and close files
  fOut    -> cd();
  hEff    -> Write();
  fEff    -> Write();
  hSmooth -> Write();
  fOut    -> Close();
  fIn     -> Write();
  cout << "  Plotting script finished.\n" << endl;

}

// End ------------------------------------------------------------------------
