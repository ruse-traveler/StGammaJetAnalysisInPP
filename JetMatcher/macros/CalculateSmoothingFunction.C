// 'CalculateSmoothingFunction.C'
// Derek Anderson
// 05.28.2021
//
// Fits specified histograms with either a single
// gaussian + a constant or a double gaussian to
// be used to smooth a response matrix.


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
static const UInt_t NLoPt(2);
static const UInt_t NHiPt(2);
static const UInt_t NHist(NLoPt + NHiPt);
static const UInt_t NVtx(4);
static const UInt_t NPts(2);
static const UInt_t NLimit(2);
static const UInt_t NParLo(4);
static const UInt_t NParHi(6);
static const UInt_t NSubFit(2);
static const UInt_t NSubPar(3);



void CalculateSmoothingFunction() {

  // lower verbosity
  gStyle -> SetOptStat("emr");
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning smoothing function calculation..." << endl;

  // io parameters
  //const TString sOut("fitToJetQt_forSmoothingCalc.et911pt0230x021kvz55pi0.r05a065rm1chrg.d28m5y2021.root");
  const TString sOut("test.root");
  const TString sInLo[NLoPt]   = {"jetQtVsPtPar_forSmoothingCalc_pTpar0206.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root", "jetQtVsPtPar_forSmoothingCalc_pTpar061.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root"};
  const TString sInHi[NHiPt]   = {"jetQtVsPtPar_forSmoothingCalc_pTpar12.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root", "jetQtVsPtPar_forSmoothingCalc_pTpar257.et911pt0230x021Kvz55pi0.r05a065rm1chrg.d28m5y2021.root"};
  const TString sHistLo[NLoPt] = {"hNorm", "hNorm"};
  const TString sHistHi[NHiPt] = {"hNorm", "hNorm"};

  // histogram parameters
  const TString sTitle("");
  const TString sTitleX("q_{T}^{jet} = p_{T}^{det} / p_{T}^{par} [GeV/c]");
  const TString sTitleY("a. u.");
  const TString sPrefix("hQtEt911_pt");
  const TString sSuffixLo[NLoPt]        = {"0p2", "0p6"};
  const TString sSuffixHi[NHiPt]        = {"1", "2"};
  const UInt_t  fMarLo[NLoPt]           = {24, 26};
  const UInt_t  fMarHi[NHiPt]           = {32, 25};
  const UInt_t  fColLo[NLoPt]           = {923, 819};
  const UInt_t  fColHi[NHiPt]           = {859, 899};
  const UInt_t  fColSub[NHiPt][NSubFit] = {{856, 896}, {854, 894}};
  const Float_t xPlot[NPts]             = {0., 2.};

  // legend parameters
  const TString sLegPrefix("p_{T}^{par} #in ");
  const TString sLegSuffix(" GeV/c");
  const TString sLegLo[NLoPt] = {"(0.2, 0.6)", "(0.6, 1)"};
  const TString sLegHi[NHiPt] = {"(1, 2)", "(2, 57)"};

  // fit parameters
  const UInt_t  fWidFit(3);
  const UInt_t  fLinFit(9);
  const UInt_t  fLinSubGaus(1);
  const TString sFitLo("gaus(0) + [3]");
  const TString sFitHi("gaus(0) + gaus(3)");
  const TString sLoPrefix("hLoEt911_pt");
  const TString sHiPrefix("hHiEt911_pt");
  const TString sLoFuncPrefix("fLoEt911_pt");
  const TString sHiFuncPrefix("fHiEt911_pt");
  const TString sSubPrefix[NSubFit]       = {"fTailHiEt911_pt", "fCoreHiEt911_pt"};
  const TString sSubFit[NSubFit]          = {"gaus(0)", "gaus(0)"};
  const TString sParLo[NParLo]            = {"Amp", "#mu", "#sigma", "Const"};
  const TString sParHi[NParHi]            = {"A_{0}", "#mu_{0}", "#sigma_{0}", "A_{1}", "#mu_{1}", "#sigma_{1}"};
  const TString sSubPar[NSubFit][NSubPar] = {{"A_{0}", "#mu_{0}", "#sigma_{0}"}, {"A_{1}", "#mu_{1}", "#sigma_{1}"}};
  const Float_t pSubGuess[NSubPar]        = {0.05, 0.7, 0.2};
  const Float_t pMuLim[NSubFit][NLimit]   = {{0.9, 1.1}, {0.6, 0.9}};
  const Float_t pRmsLimits1[NLimit]       = {0.1, 0.8};
  const Float_t xFuncDef[NPts]            = {0., 2.};

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

  // open low pT input files
  TFile *fInLo[NHist];
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    fInLo[iLoPt] = new TFile(sInLo[iLoPt].Data(), "read");
    if (!fInLo[iLoPt]) {
      cerr << "PANIC: couldn't open low pT input file #" << iLoPt << "!\n" << endl;
      return;
    }
  }

  // open high pT input files
  TFile *fInHi[NHist];
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    fInHi[iHiPt] = new TFile(sInHi[iHiPt].Data(), "read");
    if (!fInHi[iHiPt]) {
      cerr << "PANIC: couldn't open high pT input file #" << iHiPt << "!\n" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab low pT input histograms
  TH1D   *hInputLo[NLoPt];
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    hInputLo[iLoPt] = (TH1D*) fInLo[iLoPt] -> Get(sHistLo[iLoPt].Data());
    if (!hInputLo[iLoPt]) {
      cerr << "PANIC: couldn't grab low pT input histogram #" << iLoPt << "!\n" << endl;
      return;
    }
  }

  // grab high pT input histograms
  TH1D   *hInputHi[NHiPt];
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    hInputHi[iHiPt] = (TH1D*) fInHi[iHiPt] -> Get(sHistHi[iHiPt].Data());
    if (!hInputHi[iHiPt]) {
      cerr << "PANIC: couldn't grab high pT input histogram #" << iHiPt << "!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed input histograms." << endl;


  // set names of and normalize low pT input
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {

    // create histogram name and set
    TString sName(sPrefix.Data());
    sName.Append(sSuffixLo[iLoPt]);
    hInputLo[iLoPt] -> SetName(sName.Data());

    // check if integral is nonzero and normalize if so
    const Double_t intPtBin  = hInputLo[iLoPt] -> Integral();
    const Bool_t   intIsNot0 = (intPtBin > 0.);
    if (intIsNot0) hInputLo[iLoPt] -> Scale(1. / intPtBin);

  }  // end low pT bin loop

  // set names of and normalize high pT input
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {

    // create histogram name and set
    TString sName(sPrefix.Data());
    sName.Append(sSuffixHi[iHiPt]);
    hInputHi[iHiPt] -> SetName(sName.Data());

    // check if integral is nonzero and normalize if so
    const Double_t intPtBin  = hInputHi[iHiPt] -> Integral();
    const Bool_t   intIsNot0 = (intPtBin > 0.);
    if (intIsNot0) hInputHi[iHiPt] -> Scale(1. / intPtBin);

  }  // end low pT bin loop
  cout << "    Set input histogram names and normlized histograms." << endl;


  // set histogram styles
  const UInt_t  fLin(1);
  const UInt_t  fFil(0);
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    hInputLo[iLoPt] -> SetMarkerColor(fColLo[iLoPt]);
    hInputLo[iLoPt] -> SetMarkerStyle(fMarLo[iLoPt]);
    hInputLo[iLoPt] -> SetLineColor(fColLo[iLoPt]);
    hInputLo[iLoPt] -> SetLineStyle(fLin);
    hInputLo[iLoPt] -> SetFillColor(fColLo[iLoPt]);
    hInputLo[iLoPt] -> SetFillStyle(fFil);
    hInputLo[iLoPt] -> SetTitle(sTitle.Data());
    hInputLo[iLoPt] -> SetTitleFont(fTxt);
    hInputLo[iLoPt] -> GetXaxis() -> SetLabelSize(fLab);
    hInputLo[iLoPt] -> GetXaxis() -> SetLabelFont(fTxt);
    hInputLo[iLoPt] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hInputLo[iLoPt] -> GetXaxis() -> SetTitleFont(fTxt);
    hInputLo[iLoPt] -> GetXaxis() -> SetTitleOffset(fOffX);
    hInputLo[iLoPt] -> GetXaxis() -> CenterTitle(fCnt);
    hInputLo[iLoPt] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hInputLo[iLoPt] -> GetYaxis() -> SetLabelSize(fLab);
    hInputLo[iLoPt] -> GetYaxis() -> SetLabelFont(fTxt);
    hInputLo[iLoPt] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hInputLo[iLoPt] -> GetYaxis() -> SetTitleFont(fTxt);
    hInputLo[iLoPt] -> GetYaxis() -> SetTitleOffset(fOffY);
    hInputLo[iLoPt] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end low pT loop

  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    hInputHi[iHiPt] -> SetMarkerColor(fColHi[iHiPt]);
    hInputHi[iHiPt] -> SetMarkerStyle(fMarHi[iHiPt]);
    hInputHi[iHiPt] -> SetLineColor(fColHi[iHiPt]);
    hInputHi[iHiPt] -> SetLineStyle(fLin);
    hInputHi[iHiPt] -> SetFillColor(fColHi[iHiPt]);
    hInputHi[iHiPt] -> SetFillStyle(fFil);
    hInputHi[iHiPt] -> SetTitle(sTitle.Data());
    hInputHi[iHiPt] -> SetTitleFont(fTxt);
    hInputHi[iHiPt] -> GetXaxis() -> SetLabelSize(fLab);
    hInputHi[iHiPt] -> GetXaxis() -> SetLabelFont(fTxt);
    hInputHi[iHiPt] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hInputHi[iHiPt] -> GetXaxis() -> SetTitleFont(fTxt);
    hInputHi[iHiPt] -> GetXaxis() -> SetTitleOffset(fOffX);
    hInputHi[iHiPt] -> GetXaxis() -> CenterTitle(fCnt);
    hInputHi[iHiPt] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hInputHi[iHiPt] -> GetYaxis() -> SetLabelSize(fLab);
    hInputHi[iHiPt] -> GetYaxis() -> SetLabelFont(fTxt);
    hInputHi[iHiPt] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hInputHi[iHiPt] -> GetYaxis() -> SetTitleFont(fTxt);
    hInputHi[iHiPt] -> GetYaxis() -> SetTitleOffset(fOffY);
    hInputHi[iHiPt] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end high pT loop
  cout << "    Set histogram styles." << endl;


  // fit projections
  TF1  *fLoPt[NLoPt];
  TF1  *fHiPt[NHiPt];
  TF1  *fSub[NHiPt][NSubFit];
  TH1D *hLoPt[NLoPt];
  TH1D *hHiPt[NHiPt];
  cout << "    Fitting low pT histograms:" << endl;

  // fit low pT projections
  Double_t parLo[NParLo];
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {

    cout << "      Fitting low pT histogram " << (iLoPt + 1) << "/" << NLoPt << "..." << endl;

    // create histogram names
    TString sPtLo(sLoPrefix.Data());
    sPtLo.Append(sSuffixLo[iLoPt]);

    // create function names
    TString sFuncLo(sLoFuncPrefix.Data());
    sFuncLo.Append(sSuffixLo[iLoPt]);

    // create histograms to fit
    hLoPt[iLoPt] = (TH1D*) hInputLo[iLoPt] -> Clone();
    hLoPt[iLoPt] -> SetName(sPtLo.Data());

    // extract info for fits
    const UInt_t   iStart = hLoPt[iLoPt] -> FindFirstBinAbove(0.);
    const UInt_t   iStop  = hLoPt[iLoPt] -> FindLastBinAbove(0.);
    const Double_t xStart = hLoPt[iLoPt] -> GetBinCenter(iStart);
    const Double_t xStop  = hLoPt[iLoPt] -> GetBinCenter(iStop);
    parLo[0] = hLoPt[iLoPt] -> GetMaximum();
    parLo[1] = hLoPt[iLoPt] -> GetMean();
    parLo[2] = hLoPt[iLoPt] -> GetRMS();
    parLo[3] = hLoPt[iLoPt] -> GetMinimum(0.);
    cout << "CHECK: parLo[3] = " << parLo[3] << endl;

    // create gaussian functions and fit
    fLoPt[iLoPt] = new TF1(sFuncLo.Data(), sFitLo.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParLo = 0; iParLo < NParLo; iParLo++) {
      fGaus[iLoPt] -> SetParameter(iParLo, parLo[iParLo]);
      fGaus[iLoPt] -> SetParName(iParLo, sParLo[iParLo].Data());
    }  // end parameter loop
    fLoPt[iLoPt] -> SetLineColor(fColLo[iLoPt]);
    fLoPt[iLoPt] -> SetLineStyle(fLinFit);
    fLoPt[iLoPt] -> SetLineWidth(fWidFit);
    hLoPt[iLoPt] -> Fit(sFuncLoPt.Data(), "", "", xStart, xStop);

  }  // end low pT loop
  cout << "    Finished fitting low pT histograms.\n"
       << "    Fitting high pT histograms:"
       << endl;

  // fit high pT histograms
  Double_t parHi[NSubPar];
  TString  sFuncSub[NHiPt][NSubFit];
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {

    cout << "      Fitting high pT histogram " << (iHiPt + 1) << "/" << NHiPt << "..." << endl;

    // create histogram names
    TString sPtHi(sGaus2Prefix.Data());
    sPtHi.Append(sSuffixHi[iHiPt]);

    // create function names
    TString sFuncHi(sHiFuncPrefix.Data());
    sFuncHi.Append(sSuffixHi[iHiPt]);
    for (UInt_t iSubFit = 0; iSubFit < NSubFit; iSubFit++) {
      sFuncSub[iHiPt][iSubFit] = sSubPrefix[iSubFit];
      sFuncSub[iHiPt][iSubFit].Append(sSuffixHi[iHiPt]);
    }
    cout << "CHECK: sFuncSub = " << sFuncSub[iHiPt][iSubFit].Data() << endl;

    // create histograms to fit
    hHiPt[iHiPt] = (TH1D*) hInputHi[iHiPt] -> Clone();
    hHiPt[iHiPt] -> SetName(sPtHi.Data());

    // extract info for fits
    const UInt_t   iStart = hHiPt[iHiPt] -> FindFirstBinAbove(0.);
    const UInt_t   iStop  = hHiPt[iHiPt] -> FindLastBinAbove(0.);
    const Double_t xStart = hHiPt[iHiPt] -> GetBinCenter(iStart);
    const Double_t xStop  = hHiPt[iHiPt] -> GetBinCenter(iStop);
    parHi[0] = hHiPt[iHiPt] -> GetMaximum();
    parHi[1] = hHiPt[iHiPt] -> GetMean();
    parHi[2] = hHiPt[iHiPt] -> GetRMS();

    // create sub fits
    for (UInt_t iSubFit = 0; iSubFit < NSubFit; iSubFit++) {
      fSub[iHiPt][iSubFit] = new TF1(sFuncSub[iHiPt][iSubFit].Data(), sSubFit[iSubFit].Data(), xFuncDef[0], xFuncDef[1]);
      for (UInt_t iSubPar = 0; iSubPar < NSubPar; iSubPar++) {
        fSub[iHiPt][iSubFit] -> SetParName(iSubPar, sSubPar[iSubFit][iSubPar].Data());
        if (iSubFit == 0) {
          fSub[iHiPt][iSubFit] -> SetParameter(iSubPar, parHi[iSubPar]);
          if (iSubPar == 1) fSub[iHiPt][iSubFit] -> SetParLimits(iSubPar, pMuLim[iSubFit][0], pMuLim[iSubFit][1]);
        }
        else {
          fSub[iHiPt][iSubFit] -> SetParameter(iSubPar, pSubGuess[iSubPar]);
          if (iSubPar == 1) fSub[iHiPt][iSubFit] -> SetParLimits(iSubPar, pMuLim[iSubFit][0], pMuLim[iSubFit][1]);
        }
      }  // end sub fit parameter loop
      fSub[iHiPt][iSubFit] -> SetLineColor(fColSub[iHiPt][iSubFit]);
      fSub[iHiPt][iSubFit] -> SetLineStyle(fLinSubGaus);
      fSub[iHiPt][iSubFit] -> SetLineWidth(fWidFit);
      hHiPt[iHiPt]         -> Fit(sFunbcSub[iHiPt][iSubFit].Data(), "", "", xStart, xStop);
    }  // end sub fit loop

    // create double gaussian functions and fit
    fGaus2[iHiPt] = new TF1(sFuncGaus2.Data(), sGaus2.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParGaus2 = 0; iParGaus2 < NParGaus2; iParGaus2++) {
      fGaus2[iHiPt] -> SetParName(iParGaus2, sParGaus2[iParGaus2].Data());
      if (iParGaus2 < 3) {
        fGaus2[iHiPt] -> SetParameter(iParGaus2, fGaus[iHiPt] -> GetParameter(iParGaus2));
        if (iParGaus2 == 1) fGaus2[iHiPt] -> SetParLimits(iParGaus2, pMeanLimits0[0], pMeanLimits0[1]);
      }
      else {
        fGaus2[iHiPt] -> SetParameter(iParGaus2, pGuessGaus2[iParGaus2 - 3]);
        if (iParGaus2 == 4) fGaus2[iHiPt] -> SetParLimits(iParGaus2, pMeanLimits1[0], pMeanLimits1[1]);
        if (iParGaus2 == 5) fGaus2[iHiPt] -> SetParLimits(iParGaus2, pRmsLimits1[0], pRmsLimits1[1]);
      }
    }
    fGaus2[iHiPt] -> SetLineColor(fCol[iHiPt]);
    fGaus2[iHiPt] -> SetLineStyle(fLinFit);
    fGaus2[iHiPt] -> SetLineWidth(fWidFit);
    hGaus2[iHiPt] -> Fit(sFuncGaus2.Data(), "B", "", xStart, xStop);

  }  // end high pT loop
  cout << "    Finished fitting high pT histograms." << endl;


  // make legends
  const UInt_t  fColL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fAlnL(12);
  const UInt_t  fAlnT(32);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  // create labels
  TString sLoLabel[NLoPt];
  TString sHiLabel[NHiPt];
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    sLoLabel[iLoPt] = sLegPrefix;
    sLoLabel[iLoPt].Append(sLegLo[iLoPt].Data());
    sLoLabel[iLoPt].Append(sLegSuffix.Data());
  }
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    sHiLabel[iHiPt] = sLegPrefix;
    sHiLabel[iHiPt].Append(sLegHi[iHiPt].Data());
    sHiLabel[iHiPt].Append(sLegSuffix.Data());
  }

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  leg -> SetFillColor(fColL);
  leg -> SetFillStyle(fFilL);
  leg -> SetLineColor(fColL);
  leg -> SetLineStyle(fLinL);
  leg -> SetTextAlign(fAlnL);
  leg -> SetTextFont(fTxt);
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    leg -> AddEntry(hInputLo[iLoPt], sLabel[iLoPt].Data());
  }
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    leg -> AddEntry(hInputHi[iHiPt], sLabel[iHiPt].Data());
  }

  TLegend *legF = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legF -> SetFillColor(fColL);
  legF -> SetFillStyle(fFilL);
  legF -> SetLineColor(fColL);
  legF -> SetLineStyle(fLinL);
  legF -> SetTextAlign(fAlnL);
  legF -> SetTextFont(fTxt);
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    legF -> AddEntry(hLoPt[iLoPt], sLabel[iLoPt].Data());
  }
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    legF -> AddEntry(hHiPt[iHiPt], sLabel[iHiPt].Data());
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


  // save histograms
  fOut -> cd();
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    hInputLo[iLoPt] -> Write();
    fLoPt[iLoPt]    -> Write();
  }
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    hInputHi[iHiPt] -> Write();
    fHiPt[iHiPt]    -> Write();
    for (UInt_t iSubFit = 0; iSubFit < NSubFit; iSubFit++) {
      fSub[iHiPt][iSubFit] -> Write();
    }
  }
  cout << "    Saved histograms." << endl;


  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iLoPt = 0; iLoPt < NLoPt; iLoPt++) {
    fInLo[iLoPt] -> cd();
    fInLo[iLoPt] -> Close();
  }
  for (UInt_t iHiPt = 0; iHiPt < NHiPt; iHiPt++) {
    fInHi[iHiPt] -> cd();
    fInHi[iHiPt] -> Close();
  }
  cout << "  Finished smoothing function calculation!\n" << endl;

}

// End ------------------------------------------------------------------------
