// 'MakeEffVsPseudoEffPlot.C'
// Derek Anderson
// 10.17.2021
//
// Use this to plot the efficiency, the pseudo-
// efficiency, and their fit functions.

#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TPad.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NVtx(4);
static const UInt_t NHist(2);
static const UInt_t NFunc(3);
static const UInt_t NPlot(2);
static const UInt_t NParMax(5);



void MakeEffVsPseudoEffPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning eff vs. pseudoeff plot macro..." << endl;

  // io parameters
  const TString sOutput("effVsPseudoEff_fitVsHistPlot.et920vz55pi0.d19m10y2021.root");
  const TString sInput[NHist]     = {"trackEffVsPseduoEff.onlyPtParGT0.et920vz55pi0.d2m3y2020.root", "trackEffVsPseduoEff.onlyPtParGT0.et920vz55pi0.d2m3y2020.root"};
  const TString sHist[NHist]      = {"hEff_pTparOverPtPar", "hPseudoEff_AllVtx"};
  const TString sNameHist[NHist]  = {"hEff", "hPseudoEff"};
  const TString sLabelHist[NHist] = {"tracking efficiency, #epsilon_{trk}", "tracking pseudo-efficiency, #tilde{#epsilon}_{trk}"};

  // function parameters
  const UInt_t   nFuncPar[NFunc]       = {2, 5, 5};
  const UInt_t   fColFunc[NFunc]       = {854, 894, 874};
  const UInt_t   fColSys[NFunc]        = {856, 896, 876};
  const UInt_t   fLinFunc[NFunc]       = {9, 9, 9};
  const UInt_t   fLinSys[NFunc]        = {1, 1, 1};
  const UInt_t   fWidFunc[NFunc]       = {3, 3, 3};
  const UInt_t   fWidSys[NFunc]        = {2, 2, 2};
  const UInt_t   fFilFunc[NFunc]       = {0, 0, 0};
  const UInt_t   fFilSys[NFunc]        = {0, 0, 0};
  const UInt_t   fMarFunc[NFunc]       = {1, 1, 1};
  const UInt_t   fMarSys[NFunc]        = {1, 1, 1};
  const Float_t  xFuncDef[NPlot]       = {0.05, 19.95};
  const TString  sFuncName[NFunc]      = {"fEff", "fPseudoEff", "fInterpolate"};
  const TString  sLabelFunc[NFunc]     = {"fit to #epsilon_{trk}, E[#epsilon_{trk}]",
                                          "fit to #tilde{#epsilon}_{trk}, #tilde{E}[#tilde{#epsilon}_{trk}]",
                                          "interpolation of E and #tilde{E}, #tilde{E}^{*}[#tilde{#epsilon}_{trk}]"};
  const TString  sFunc[NFunc]          = {"[0] * (1 - TMath::Exp(-1. * [1] * x))",
                                          "[0] + ([1] * TMath::Exp(-1. * [2] * x)) + ([3] * TMath::Exp(-1. * [4] * x))",
                                          "[0] + ([1] * TMath::Exp(-1. * [2] * x)) + ([3] * TMath::Exp(-1. * [4] * x))"};
  const TString sParSmall[NParMax]     = {"#epsilon_{0}", "#sigma_{1}", "NA0", "NA1", "NA2"};
  const TString sParBig[NParMax]       = {"#epsilon_{0}", "#epsilon_{1}", "#sigma_{1}", "#epsilon_{2}", "#sigma_{2}"};
  const Double_t pFunc[NFunc][NParMax] = {{0.82, 7.59, 1., 1., 1.},
                                          {0.86, -1.81, 13., 0.12, 0.67},
                                          {0.82, -2.18, 13., 0.18, 0.67}};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("#color[603]{p_{T}^{MC}}, #color[635]{p_{T}} [GeV/c]");
  const TString sTitleY("#color[603]{#epsilon_{trk}}, #color[635]{#tilde{#epsilon}_{trk}}");
  const Float_t xPlotRange[NPlot] = {0., 20.};
  const Float_t fMarSize[NHist]   = {1.5, 1};
  const UInt_t  fMarHist[NHist]   = {20, 21};
  const UInt_t  fColHist[NHist]   = {603, 635};
  const UInt_t  fWidHist[NHist]   = {1, 1};

  // text parameters
  const TString sSys("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  const TString sErr("bands indicate #pm 4 sys. uncertainty");
  const TString sEtc("using only recoil tracks");

  // open output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n"
         << "       fOutput = " << fOutput
         << endl;
    return;
  }

  // open input files
  TFile *fInput[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInput[iHist] = new TFile(sInput[iHist].Data(), "read");
    if (!fInput[iHist]) {
      cerr << "PANIC: couldn't open input file #" << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input histograms
  TH1D *hEff[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hEff[iHist] = (TH1D*) fInput[iHist] -> Get(sHist[iHist]);
    if (!hEff[iHist]) {
      cerr << "PANIC: couldn't grab input histogram #" << iHist << "!" << endl;
      return;
    }
    hEff[iHist] -> SetName(sNameHist[iHist].Data());
  }
  cout << "    Grabbed input histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hEff[iHist] -> SetMarkerColor(fColHist[iHist]);
    hEff[iHist] -> SetMarkerStyle(fMarHist[iHist]);
    hEff[iHist] -> SetMarkerSize(fMarSize[iHist]);
    hEff[iHist] -> SetFillColor(fColHist[iHist]);
    hEff[iHist] -> SetFillStyle(fFil);
    hEff[iHist] -> SetLineColor(fColHist[iHist]);
    hEff[iHist] -> SetLineStyle(fLin);
    hEff[iHist] -> SetLineWidth(fWidHist[iHist]);
    hEff[iHist] -> SetTitle(sTitle.Data());
    hEff[iHist] -> SetTitleFont(fTxt);
    hEff[iHist] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hEff[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hEff[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hEff[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hEff[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hEff[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hEff[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hEff[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hEff[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hEff[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hEff[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hEff[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hEff[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hEff[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hEff[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make functions
  const TString sPlusSys(" + 0.04");
  const TString sMinusSys(" - 0.04");
  const TString sPlusName("PlusSys");
  const TString sMinusName("MinusSys");

  TF1 *fFunc[NFunc];
  TF1 *fSysPlus[NFunc];
  TF1 *fSysMinus[NFunc];
  for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {
    // create sys names/functions
    TString sFuncPlus(sFunc[iFunc].Data());
    TString sFuncMinus(sFunc[iFunc].Data());
    TString sNamePlus(sFuncName[iFunc].Data());
    TString sNameMinus(sFuncName[iFunc].Data());
    sFuncPlus.Append(sPlusSys.Data());
    sFuncMinus.Append(sMinusSys.Data());
    sNamePlus.Append(sPlusName.Data());
    sNameMinus.Append(sMinusName.Data());

    // create functions
    fFunc[iFunc]     = new TF1(sFuncName[iFunc].Data(), sFunc[iFunc].Data(), xFuncDef[0], xFuncDef[1]);
    fSysPlus[iFunc]  = new TF1(sNamePlus.Data(), sFuncPlus.Data(), xFuncDef[0], xFuncDef[1]);
    fSysMinus[iFunc] = new TF1(sNameMinus.Data(), sFuncMinus.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iPar = 0; iPar < nFuncPar[iFunc]; iPar++) {
      fFunc[iFunc]     -> SetParameter(iPar, pFunc[iFunc][iPar]);
      fSysPlus[iFunc]  -> SetParameter(iPar, pFunc[iFunc][iPar]);
      fSysMinus[iFunc] -> SetParameter(iPar, pFunc[iFunc][iPar]);
      if (iFunc == 0) {
        fFunc[iFunc]     -> SetParName(iPar, sParSmall[iPar]);
        fSysPlus[iFunc]  -> SetParName(iPar, sParSmall[iPar]);
        fSysMinus[iFunc] -> SetParName(iPar, sParSmall[iPar]);
      } else {
        fFunc[iFunc]     -> SetParName(iPar, sParBig[iPar]);
        fSysPlus[iFunc]  -> SetParName(iPar, sParBig[iPar]);
        fSysMinus[iFunc] -> SetParName(iPar, sParBig[iPar]);
      }
    }
    fFunc[iFunc]     -> SetLineColor(fColFunc[iFunc]);
    fFunc[iFunc]     -> SetLineStyle(fLinFunc[iFunc]);
    fFunc[iFunc]     -> SetLineWidth(fWidFunc[iFunc]);
    fFunc[iFunc]     -> SetFillColor(fColSys[iFunc]);
    fFunc[iFunc]     -> SetFillStyle(fFilFunc[iFunc]);
    fFunc[iFunc]     -> SetMarkerColor(fColFunc[iFunc]);
    fFunc[iFunc]     -> SetMarkerStyle(fMarFunc[iFunc]);
    fSysPlus[iFunc]  -> SetLineColor(fColSys[iFunc]);
    fSysPlus[iFunc]  -> SetLineStyle(fLinSys[iFunc]);
    fSysPlus[iFunc]  -> SetLineWidth(fWidSys[iFunc]);
    fSysPlus[iFunc]  -> SetFillColor(fColSys[iFunc]);
    fSysPlus[iFunc]  -> SetFillStyle(fFilSys[iFunc]);
    fSysPlus[iFunc]  -> SetMarkerColor(fColSys[iFunc]);
    fSysPlus[iFunc]  -> SetMarkerStyle(fMarSys[iFunc]);
    fSysMinus[iFunc] -> SetLineColor(fColSys[iFunc]);
    fSysMinus[iFunc] -> SetLineStyle(fLinSys[iFunc]);
    fSysMinus[iFunc] -> SetLineWidth(fWidSys[iFunc]);
    fSysMinus[iFunc] -> SetFillColor(fColSys[iFunc]);
    fSysMinus[iFunc] -> SetFillStyle(fFilSys[iFunc]);
    fSysMinus[iFunc] -> SetMarkerColor(fColSys[iFunc]);
    fSysMinus[iFunc] -> SetMarkerStyle(fMarSys[iFunc]);
  }
  cout << "    Made functions." << endl;

  // make legend
  const UInt_t  fColLe       = 0;
  const UInt_t  fFilLe       = 0;
  const UInt_t  fLinLe       = 0;
  const UInt_t  nObj         = NHist + NFunc;
  const UInt_t  nPlot        = TMath::Max(NHist, NFunc);
  const Float_t fHeightObj   = 0.05;
  const Float_t fHeightLe    = nObj * fHeightObj;
  const Float_t fTopLe       = 0.1 + fHeightLe;
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, fTopLe};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iPlot = 0; iPlot < nPlot; iPlot++) {
    if (iPlot < NFunc) leg -> AddEntry(fFunc[iPlot], sLabelFunc[iPlot], "lf");
    if (iPlot < NHist) leg -> AddEntry(hEff[iPlot], sLabelHist[iPlot], "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, 0.3};
  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sErr.Data());
  txt -> AddText(sEtc.Data());
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(3);
  const Float_t fLinXY[NVtx] = {xPlotRange[0], 1., xPlotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(0);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.15);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot   -> SetGrid(fGrid, fGrid);
  cPlot   -> SetTicks(fTick, fTick);
  cPlot   -> SetLogx(fLogX);
  cPlot   -> SetLogy(fLogY);
  cPlot   -> SetBorderMode(fMode);
  cPlot   -> SetBorderSize(fBord);
  cPlot   -> SetFrameBorderMode(fFrame);
  cPlot   -> SetLeftMargin(fMarginL);
  cPlot   -> SetRightMargin(fMarginR);
  cPlot   -> SetTopMargin(fMarginT);
  cPlot   -> SetBottomMargin(fMarginB);
  cPlot   -> cd();
  hEff[0] -> Draw();
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hEff[iHist] -> Draw("same");
  }
  for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {
    fFunc[iFunc]     -> Draw("same");
    fSysPlus[iFunc]  -> Draw("same");
    fSysMinus[iFunc] -> Draw("same");
  }
  line    -> Draw();
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;

  // save histograms and functions
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hEff[iHist] -> Write();
  }
  for (UInt_t iFunc = 0; iFunc < NFunc; iFunc++) {
    fFunc[iFunc]     -> Write();
    fSysPlus[iFunc]  -> Write();
    fSysMinus[iFunc] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInput[iHist] -> cd();
    fInput[iHist] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// End ------------------------------------------------------------------------
