// 'QuickFitTest.C'
// Derek Anderson
// 04.26.2021
// 
// Fits specified histograms with a gaussian, and then uses the
// parameters from the guassian to inform a couple fits based
// on some representations of the crystal ball.
//
// Crystal Ball Function from:
//   https://root-forum.cern.ch/t/crystalball-fitting-function/26973


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const UInt_t NHist(4);
static const UInt_t NVtx(4);
static const UInt_t NPts(2);
static const UInt_t NLimit(2);
static const UInt_t NParGaus(3);
static const UInt_t NParGaus2(6);
static const UInt_t NParErf(6);
static const UInt_t NParXtal(5);
static const Bool_t DoXtalBallFit(false);



void QuickFitTest() {

  // lower verbosity
  gStyle -> SetOptStat("emrf");
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning quick fitting test..." << endl;

  // io parameters
  const TString sOut("testingLimits.root");
  //const TString sOut("fitToJetQt_forSmoothingCalc.et911pt0230x021kvz55pi0.r05a065rm1chrg.d28m5y2021.root");
  const TString sIn[NHist]   = {"jetQtVsPtPar_forSmoothingCalc_pTpar0206.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root", "jetQtVsPtPar_forSmoothingCalc_pTpar061.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root", "jetQtVsPtPar_forSmoothingCalc_pTpar12.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root", "jetQtVsPtPar_forSmoothingCalc_pTpar257.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root"};
  const TString sHist[NHist] = {"hNorm", "hNorm", "hNorm", "hNorm"};

  // histogram parameters
  const TString sTitle("");
  const TString sTitleX("q_{T}^{jet} = p_{T}^{det} / p_{T}^{par} [GeV/c]");
  const TString sTitleY("a. u.");
  const TString sPrefix("hQtEt911_pt");
  const TString sSuffix[NHist] = {"0p2", "0p6", "1", "2"};
  const UInt_t  fCol[NHist]    = {923, 819, 859, 899};
  const UInt_t  fColG0[NHist]  = {922, 816, 856, 896};
  const UInt_t  fColG1[NHist]  = {920, 814, 854, 894};
  const UInt_t  fMar[NHist]    = {24, 26, 32, 25};
  const Float_t xPlot[NPts]    = {0., 2.};

  // legend parameters
  const TString sLegPrefix("p_{T}^{par} #in ");
  const TString sLegSuffix(" GeV/c");
  const TString sLegBin[NHist] = {"(0.2, 0.6)", "(0.6, 1)", "(1, 2)", "(2, 57)"};

  // fit parameters
  const UInt_t  fWidFit(3);
  const UInt_t  fLinFit(9);
  const UInt_t  fLinSubGaus(1);
  const TString sGaus("gaus(0)");
  const TString sGaus2("gaus(0) + gaus(3)");
  const TString sModErf("(TMath::Erf((x - [0]) / [1]) + 1.) * ([2] * TMath::Power(x, [3]) / (TMath::Power(1 + (x / [4]), [5])))");
  const TString sXtalBall("[0] * (ROOT::Math::crystalball_function(x, [1], [2], [3], [4])}");
  const TString sGausPrefix("hGausPt");
  const TString sGaus0Prefix("hGausCorePt");
  const TString sGaus1Prefix("hGausTailPt");
  const TString sGaus2Prefix("hGausSumPt");
  const TString sErfPrefix("hErf");
  const TString sXtalPrefix("hXtalPt");
  const TString sGausFuncPrefix("fGausEt911_pt");
  const TString sGaus0FuncPrefix("fGausCoreEt911_pt");
  const TString sGaus1FuncPrefix("fGausTailEt911_pt");
  const TString sGaus2FuncPrefix("fGausSumEt911_pt");
  const TString sErfFuncPrefix("fErfEt911_pt");
  const TString sXtalFuncPrefix("fXtalEt911_pt");
  const TString sParGaus[NParGaus]    = {"Amp", "#mu", "#sigma"};
  const TString sParGaus2[NParGaus2]  = {"A_{0}", "#mu_{0}", "#sigma_{0}", "A_{1}", "#mu_{1}", "#sigma_{1}"};
  const TString sParErf[NParErf]      = {"#mu", "#sigma", "A", "B", "C", "D"};
  const TString sParXtal[NParXtal]    = {"Amp", "#alpha", "n", "#sigma", "#mu"};
  const Float_t pGuessErf[NParErf]    = {1., 1., 2.7, 5.5, 28.7, 34.3};
  const Float_t pGuessXtal[NParXtal]  = {1., 1.7, 0.23, 1., 1.};
  const Float_t pGuessGaus2[NParGaus] = {0.05, 0.7, 0.2};
  const Float_t pMeanLimits0[NLimit]  = {0.9, 1.1};
  const Float_t pMeanLimits1[NLimit]  = {0.6, 0.9};
  const Float_t pRmsLimits0[NLimit]   = {0., 0.2};
  const Float_t pRmsLimits1[NLimit]   = {0.1, 0.8};
  const Float_t xFuncDef[NPts]        = {0., 2.};

  // text parameters
  const TString sSystem("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrigger("#pi^{0} trigger, E_{T}^{trg} #in (9, 11) GeV");
  const TString sJets("anti-k_{T}, R = 0.5");
  const TString sType("#bf{charged jets}");


  // open output file
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  // open input files
  TFile *fIn[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] = new TFile(sIn[iHist].Data(), "read");
    if (!fIn[iHist]) {
      cerr << "PANIC: couldn't open input file #" << iHist << "!\n" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input histograms
  TH1D   *hInput[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInput[iHist] = (TH1D*) fIn[iHist] -> Get(sHist[iHist].Data());
    if (!hInput[iHist]) {
      cerr << "PANIC: couldn't grab input histogram #" << iHist << "!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed input histograms." << endl;


  // set names of and normalize input
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {

    // create histogram name and set
    TString sName(sPrefix.Data());
    sName.Append(sSuffix[iHist]);
    hInput[iHist] -> SetName(sName.Data());

    // check if integral is nonzero and normalize if so
    const Double_t intPtBin  = hInput[iHist] -> Integral();
    const Bool_t   intIsNot0 = (intPtBin > 0.);
    if (intIsNot0) hInput[iHist] -> Scale(1. / intPtBin);

  }  // end pTpar bin loop
  cout << "    Set input histogram names and normlized histograms." << endl;

  // set histogram styles
  const UInt_t  fLin(1);
  const UInt_t  fFil(0);
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInput[iHist] -> SetMarkerColor(fCol[iHist]);
    hInput[iHist] -> SetMarkerStyle(fMar[iHist]);
    hInput[iHist] -> SetLineColor(fCol[iHist]);
    hInput[iHist] -> SetLineStyle(fLin);
    hInput[iHist] -> SetFillColor(fCol[iHist]);
    hInput[iHist] -> SetFillStyle(fFil);
    hInput[iHist] -> SetTitle(sTitle.Data());
    hInput[iHist] -> SetTitleFont(fTxt);
    hInput[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hInput[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hInput[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hInput[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hInput[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hInput[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hInput[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hInput[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hInput[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hInput[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hInput[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hInput[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hInput[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end hist loop
  cout << "    Set histogram styles." << endl;


  // fit projections
  TF1  *fGaus[NHist];
  TF1  *fGaus0[NHist];
  TF1  *fGaus1[NHist];
  TF1  *fGaus2[NHist];
  TF1  *fErf[NHist];
  TF1  *fXtal[NHist];
  TH1D *hGaus[NHist];
  TH1D *hGaus2[NHist];
  TH1D *hErf[NHist];
  TH1D *hXtal[NHist];
  cout << "    Fitting histograms:" << endl;

  Double_t parGaus[NParGaus];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {

    cout << "      Fitting histogram " << (iHist + 1) << "/" << NHist << "..." << endl;

    // create histogram names
    TString sPtGaus(sGausPrefix.Data());
    TString sPtGaus0(sGaus0Prefix.Data());
    TString sPtGaus1(sGaus1Prefix.Data());
    TString sPtGaus2(sGaus2Prefix.Data());
    TString sPtErf(sErfPrefix.Data());
    TString sPtXtal(sXtalPrefix.Data());
    sPtGaus.Append(sSuffix[iHist]);
    sPtGaus0.Append(sSuffix[iHist]);
    sPtGaus1.Append(sSuffix[iHist]);
    sPtGaus2.Append(sSuffix[iHist]);
    sPtErf.Append(sSuffix[iHist]);
    sPtXtal.Append(sSuffix[iHist]);

    // create function names
    TString sFuncGaus(sGausFuncPrefix.Data());
    TString sFuncGaus0(sGaus0FuncPrefix.Data());
    TString sFuncGaus1(sGaus1FuncPrefix.Data());
    TString sFuncGaus2(sGaus2FuncPrefix.Data());
    TString sFuncErf(sErfFuncPrefix.Data());
    TString sFuncXtal(sXtalFuncPrefix.Data());
    sFuncGaus.Append(sSuffix[iHist]);
    sFuncGaus0.Append(sSuffix[iHist]);
    sFuncGaus1.Append(sSuffix[iHist]);
    sFuncGaus2.Append(sSuffix[iHist]);
    sFuncErf.Append(sSuffix[iHist]);
    sFuncXtal.Append(sSuffix[iHist]);

    // create histograms to fit
    hGaus[iHist]  = (TH1D*) hInput[iHist] -> Clone();
    hGaus2[iHist] = (TH1D*) hInput[iHist] -> Clone();
    hErf[iHist]   = (TH1D*) hInput[iHist] -> Clone();
    hXtal[iHist]  = (TH1D*) hInput[iHist] -> Clone();
    hGaus[iHist]  -> SetName(sPtGaus.Data());
    hGaus2[iHist] -> SetName(sPtGaus2.Data());
    hErf[iHist]   -> SetName(sPtErf.Data());
    hXtal[iHist]  -> SetName(sPtXtal.Data());

    // extract info for fits
    const UInt_t   iStart = hGaus[iHist] -> FindFirstBinAbove(0.);
    const UInt_t   iStop  = hGaus[iHist] -> FindLastBinAbove(0.);
    const Double_t xStart = hGaus[iHist] -> GetBinCenter(iStart);
    const Double_t xStop  = hGaus[iHist] -> GetBinCenter(iStop);
    parGaus[0] = hGaus[iHist] -> GetMaximum();
    parGaus[1] = hGaus[iHist] -> GetMean();
    parGaus[2] = hGaus[iHist] -> GetRMS();

    // create gaussian functions and fit
    fGaus[iHist] = new TF1(sFuncGaus.Data(), sGaus.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
      fGaus[iHist] -> SetParameter(iParGaus, parGaus[iParGaus]);
      fGaus[iHist] -> SetParName(iParGaus, sParGaus[iParGaus].Data());
    }  // end parameter loop
    fGaus[iHist] -> SetLineColor(fCol[iHist]);
    fGaus[iHist] -> SetLineStyle(fLinFit);
    fGaus[iHist] -> SetLineWidth(fWidFit);
    hGaus[iHist] -> Fit(sFuncGaus.Data(), "", "", xStart, xStop);

    // create sub gaussian functions
    fGaus0[iHist] = new TF1(sFuncGaus0.Data(), sGaus.Data(), xFuncDef[0], xFuncDef[1]);
    fGaus1[iHist] = new TF1(sFuncGaus1.Data(), sGaus.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
      fGaus0[iHist] -> SetParameter(iParGaus, fGaus[iHist] -> GetParameter(iParGaus));
      fGaus1[iHist] -> SetParameter(iParGaus, pGuessGaus2[iParGaus]);
      fGaus0[iHist] -> SetParName(iParGaus, sParGaus2[iParGaus].Data());
      fGaus1[iHist] -> SetParName(iParGaus, sParGaus2[iParGaus + 3].Data());
      if (iParGaus == 1) {
        fGaus0[iHist] -> SetParLimits(iParGaus, pMeanLimits0[0], pMeanLimits0[1]);
        fGaus1[iHist] -> SetParLimits(iParGaus, pMeanLimits1[0], pMeanLimits1[1]);
      }
    }
    fGaus0[iHist] -> SetLineColor(fColG0[iHist]);
    fGaus1[iHist] -> SetLineColor(fColG1[iHist]);
    fGaus0[iHist] -> SetLineStyle(fLinSubGaus);
    fGaus1[iHist] -> SetLineStyle(fLinSubGaus);
    fGaus0[iHist] -> SetLineWidth(fWidFit);
    fGaus1[iHist] -> SetLineWidth(fWidFit);
    hGaus2[iHist] -> Fit(sFuncGaus0.Data(), "", "", xStart, xStop);
    hGaus2[iHist] -> Fit(sFuncGaus1.Data(), "", "", xStart, xStop);

    // create double gaussian functions and fit
    fGaus2[iHist] = new TF1(sFuncGaus2.Data(), sGaus2.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParGaus2 = 0; iParGaus2 < NParGaus2; iParGaus2++) {
      fGaus2[iHist] -> SetParName(iParGaus2, sParGaus2[iParGaus2].Data());
      if (iParGaus2 < 3) {
        fGaus2[iHist] -> SetParameter(iParGaus2, fGaus[iHist] -> GetParameter(iParGaus2));
        if (iParGaus2 == 1) fGaus2[iHist] -> SetParLimits(iParGaus2, pMeanLimits0[0], pMeanLimits0[1]);
        if (iParGaus2 == 2) fGaus2[iHist] -> SetParLimits(iParGaus2, pRmsLimits0[0], pRmsLimits0[1]);
      }
      else {
        fGaus2[iHist] -> SetParameter(iParGaus2, pGuessGaus2[iParGaus2 - 3]);
        if (iParGaus2 == 4) fGaus2[iHist] -> SetParLimits(iParGaus2, pMeanLimits1[0], pMeanLimits1[1]);
        if (iParGaus2 == 5) fGaus2[iHist] -> SetParLimits(iParGaus2, pRmsLimits1[0], pRmsLimits1[1]);
      }
    }
    fGaus2[iHist] -> SetLineColor(fCol[iHist]);
    fGaus2[iHist] -> SetLineStyle(fLinFit);
    fGaus2[iHist] -> SetLineWidth(fWidFit);
    hGaus2[iHist] -> Fit(sFuncGaus2.Data(), "B", "", xStart, xStop);


    // create error functions and fit
    fErf[iHist] = new TF1(sFuncErf.Data(), sModErf.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParErf = 0; iParErf < NParErf; iParErf++) {
      fErf[iHist] -> SetParName(iParErf, sParErf[iParErf]);
      switch (iParErf) {
        case 0:
          fErf[iHist] -> SetParameter(iParErf, fGaus[iHist] -> GetParameter(1));
          break;
        case 1:
          fErf[iHist] -> SetParameter(iParErf, fGaus[iHist] -> GetParameter(2));
          break;
        default:
          fErf[iHist] -> SetParameter(iParErf, pGuessErf[iParErf]);
          break;
      }
    }  // end parameter loop
    fErf[iHist] -> SetLineColor(fCol[iHist]);
    fErf[iHist] -> SetLineStyle(fLinFit);
    fErf[iHist] -> SetLineWidth(fWidFit);
    hErf[iHist] -> Fit(sFuncErf.Data(), "", "", xStart, xStop);

    // create xtal ball functions and fit
    if (DoXtalBallFit) {
      fXtal[iHist] = new TF1(sFuncXtal.Data(), sXtalBall.Data(), xFuncDef[0], xFuncDef[1]);
      for (UInt_t iParXtal = 0; iParXtal < NParXtal; iParXtal++) {
        fXtal[iHist] -> SetLineColor(fCol[iHist]);
        fXtal[iHist] -> SetLineStyle(fLinFit);
        fXtal[iHist] -> SetLineWidth(fWidFit);
        fXtal[iHist] -> SetParName(iParXtal, sParXtal[iParXtal]);
        switch (iParXtal) {
          case 4:
            fXtal[iHist] -> SetParameter(iParXtal, fGaus[iHist] -> GetParameter(2));
            break;
          case 5:
            fXtal[iHist] -> SetParameter(iParXtal, fGaus[iHist] -> GetParameter(1));
            break;
          default:
            fXtal[iHist] -> SetParameter(iParXtal, pGuessXtal[iHist]);
            break;
        }
      }  // end parameter loop
      hXtal[iHist] -> Fit(sFuncXtal.Data(), "", "", xStart, xStop);
    }

  }  // end hist loop
  cout << "    Fit histograms." << endl;


  // make legends
  const UInt_t  fColL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fAlnL(12);
  const UInt_t  fAlnT(32);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  // create labels
  TString sLabel[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    sLabel[iHist] = sLegPrefix;
    sLabel[iHist].Append(sLegBin[iHist].Data());
    sLabel[iHist].Append(sLegSuffix.Data());
  }  // end hist loop

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  leg -> SetFillColor(fColL);
  leg -> SetFillStyle(fFilL);
  leg -> SetLineColor(fColL);
  leg -> SetLineStyle(fLinL);
  leg -> SetTextAlign(fAlnL);
  leg -> SetTextFont(fTxt);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    leg -> AddEntry(hInput[iHist], sLabel[iHist].Data());
  }  // end hist loop

  TLegend *legG = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legG -> SetFillColor(fColL);
  legG -> SetFillStyle(fFilL);
  legG -> SetLineColor(fColL);
  legG -> SetLineStyle(fLinL);
  legG -> SetTextAlign(fAlnL);
  legG -> SetTextFont(fTxt);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    legG -> AddEntry(hGaus[iHist], sLabel[iHist].Data());
  }  // end hist loop

  TLegend *legG2 = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legG2 -> SetFillColor(fColL);
  legG2 -> SetFillStyle(fFilL);
  legG2 -> SetLineColor(fColL);
  legG2 -> SetLineStyle(fLinL);
  legG2 -> SetTextAlign(fAlnL);
  legG2 -> SetTextFont(fTxt);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    legG2 -> AddEntry(hGaus2[iHist], sLabel[iHist].Data());
  }  // end hist loop

  TLegend *legE = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legE -> SetFillColor(fColL);
  legE -> SetFillStyle(fFilL);
  legE -> SetLineColor(fColL);
  legE -> SetLineStyle(fLinL);
  legE -> SetTextAlign(fAlnL);
  legE -> SetTextFont(fTxt);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    legE -> AddEntry(hErf[iHist], sLabel[iHist].Data());
  }  // end hist loop

  TLegend *legC = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  if (DoXtalBallFit) {
    legC -> SetFillColor(fColL);
    legC -> SetFillStyle(fFilL);
    legC -> SetLineColor(fColL);
    legC -> SetLineStyle(fLinL);
    legC -> SetTextAlign(fAlnL);
    legC -> SetTextFont(fTxt);
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      legC -> AddEntry(hXtal[iHist], sLabel[iHist].Data());
    }  // end hist loop
  }
  cout << "    Made legends." << endl;


  // make text box
  TPaveText *txt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  txt -> SetFillColor(fColL);
  txt -> SetFillStyle(fFilL);
  txt -> SetLineColor(fColL);
  txt -> SetLineStyle(fLinL);
  txt -> SetTextAlign(fAlnT);
  txt -> SetTextFont(fTxt);
  txt -> AddText(sSystem.Data());
  txt -> AddText(sTrigger.Data());
  txt -> AddText(sJets.Data());
  txt -> AddText(sType.Data());
  cout << "    Made text box." << endl;


  // make plots
  const UInt_t  fSizX(750);
  const UInt_t  fSizY(750);
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const Float_t fMargT(0.02);
  const Float_t fMargR(0.02);

  TCanvas *cInput = new TCanvas("cInput", "", fSizX, fSizY);
  cInput    -> SetTopMargin(fMargT);
  cInput    -> SetRightMargin(fMargR);
  cInput    -> SetLogy(fLog);
  cInput    -> SetGrid(fGrid, fGrid);
  cInput    -> cd();
  hInput[0] -> Draw("hist p l");
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hInput[iHist] -> Draw("hist p l same");
  }
  leg    -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cInput -> Write();
  cInput -> Close();

  TCanvas *cGaus = new TCanvas("cGaus", "", fSizX, fSizY);
  cGaus    -> SetTopMargin(fMargT);
  cGaus    -> SetRightMargin(fMargR);
  cGaus    -> SetLogy(fLog);
  cGaus    -> SetGrid(fGrid, fGrid);
  cGaus    -> cd();
  hGaus[0] -> Draw();
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hGaus[iHist] -> Draw("sames");
  }
  legG  -> Draw();
  txt   -> Draw();
  fOut  -> cd();
  cGaus -> Write();
  cGaus -> Close();

  TCanvas *cGaus2 = new TCanvas("cGausSum", "", fSizX, fSizY);
  cGaus2    -> SetTopMargin(fMargT);
  cGaus2    -> SetRightMargin(fMargR);
  cGaus2    -> SetLogy(fLog);
  cGaus2    -> SetGrid(fGrid, fGrid);
  cGaus2    -> cd();
  hGaus2[0] -> Draw();
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hGaus2[iHist] -> Draw("sames");
  }
  legG2  -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cGaus2 -> Write();
  cGaus2 -> Close();

  TCanvas *cErf = new TCanvas("cErf", "", fSizX, fSizY);
  cErf    -> SetTopMargin(fMargT);
  cErf    -> SetRightMargin(fMargR);
  cErf    -> SetLogy(fLog);
  cErf    -> SetGrid(fGrid, fGrid);
  cErf    -> cd();
  hErf[0] -> Draw("l");
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hErf[iHist] -> Draw("l sames");
  }
  legE -> Draw();
  txt  -> Draw();
  fOut -> cd();
  cErf -> Write();
  cErf -> Close();

  if (DoXtalBallFit) {
    TCanvas *cXtal = new TCanvas("cXtal", "", fSizX, fSizY);
    cXtal    -> SetTopMargin(fMargT);
    cXtal    -> SetRightMargin(fMargR);
    cXtal    -> SetLogy(fLog);
    cXtal    -> SetGrid(fGrid, fGrid);
    cXtal    -> cd();
    hXtal[0] -> Draw("l");
    for (UInt_t iHist = 1; iHist < NHist; iHist++) {
      hXtal[iHist] -> Draw("l sames");
    }
    legC  -> Draw();
    txt   -> Draw();
    fOut  -> cd();
    cXtal -> Write();
    cXtal -> Close();
  }
  cout << "    Made plots." << endl;


  // save histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInput[iHist] -> Write();
    hGaus[iHist]  -> Write();
    fGaus[iHist]  -> Write();
    hGaus2[iHist] -> Write();
    fGaus0[iHist] -> Write();
    fGaus1[iHist] -> Write();
    fGaus2[iHist] -> Write();
    hErf[iHist]   -> Write();
    fErf[iHist]   -> Write();
    if (DoXtalBallFit) {
      hXtal[iHist] -> Write();
      fXtal[iHist] -> Write();
    }
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Finished quick fitting test!\n" << endl;

}  // end 'PlotResponseMatrixProjections()'



Double_t EvalCrystalBallFunction(Double_t x, Double_t alpha, Double_t n, Double_t sigma, Double_t mean) {

  // check if width is negative
  if (sigma < 0.) return 0.;

  Double_t z = (x - mean) / sigma; 
  if (alpha < 0) z = -z; 

  // calculate value
  Double_t abs_alpha = TMath::Abs(alpha);
  if (z  > (-1. * abs_alpha)) {
    return TMath::Exp(-0.5 * z * z);
  } 
  else {
    Double_t nDivAlpha = n / abs_alpha;
    Double_t AA        = TMath::Exp(-0.5 * abs_alpha * abs_alpha);
    Double_t B         = nDivAlpha - abs_alpha;
    Double_t arg       = nDivAlpha /(B - z);
    Double_t returnVal = AA * TMath::Power(arg, n);
    return returnVal;
  }

}  // end 'EvalCrystalBallFunction(Double_t, Double_t, Double_t, Double_t, Double_t)'



Double_t CrystalBallFunction(const Double_t *x, const Double_t *par) {

  const Double_t xtal = par[0] * EvalCrystalBallFunction(x[0], par[1], par[2], par[3], par[4]);
  return xtal;

}  // end 'CrystalBallFunction(Double_t*, Double_t*)'

// End ------------------------------------------------------------------------
