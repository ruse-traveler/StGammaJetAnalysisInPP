// 'MakeSmearingWeightPlotForPaper.C'
// Derek Anderson
// 02.23.2023
//
// Use this to make the smearing
// weight figure for the long
// paper.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

// global constants
static const UInt_t NGam(3);
static const UInt_t NPi0(2);
static const UInt_t NVtx(4);
static const UInt_t NRange(2);



void MakeSmearingWeightPlotForPaper() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning smearing weight plot maker..." << endl;

  // io parameters
  const TString sOutput("smearingWeightsForPaper.et630pi0vsGam.d23m2y2023.root");
  const TString sInGam[NGam]   = {"matrixProjections.forPaper.et650gam.d30m11y2022.root",
                                  "matrixProjections.forPaper.et650gam.d30m11y2022.root",
                                  "matrixProjections.forPaper.et650gam.d30m11y2022.root"};
  const TString sInPi0[NPi0]   = {"FromSaskia/Weight_pi0_911.root",
                                  "FromSaskia/Weight_pi0_1115.root"};
  const TString sHistGam[NGam] = {"hSmearEt911",
                                  "hSmearEt1115",
                                  "hSmearEt1520"};
  const TString sHistPi0[NPi0] = {"hEtMatch911_NoWeight",
                                  "hEtMatch1115_NoWeight"};
  const TString sNameGam[NGam] = {"hGamWeightsEt911",
                                  "hGamWeightsEt1115",
                                  "hGamWeightsEt1520"};
  const TString sNamePi0[NPi0] = {"hPi0WeightsEt911",
                                  "hPi0WeightsEt1115"};

  // text parameters
  const TString sTitleY("Relative smearing weights (arbitrary normalization)");
  const TString sTitleGam("STAR simulation #gamma");
  const TString sTitlePi0("STAR simulation #pi^{0}");
  const TString sTitleGamX("E_{T}^{#gamma,part-match} [GeV]");
  const TString sTitlePi0X("E_{T}^{#gamma,part-match} [GeV]");
  const TString sLegGam[NGam] = {"9 < E_{T}^{#gamma,det-clust} < 11 GeV",
                                 "11 < E_{T}^{#gamma,det-clust} < 15 GeV",
                                 "15 < E_{T}^{#gamma,det-clust} < 20 GeV"};
  const TString sLegPi0[NPi0] = {"9 < E_{T}^{#pi^{0},det-clust} < 11 GeV",
                                 "11 < E_{T}^{#pi^{0},det-clust} < 15 GeV"};

  // plot parameters
  const Float_t yRange[NRange]    = {0.0000007, 77.};
  const Float_t xGamRange[NRange] = {6., 30.};
  const Float_t xPi0Range[NRange] = {6., 25.};
  const UInt_t  fColGam[NGam]     = {799, 899, 859};
  const UInt_t  fColPi0[NPi0]     = {859, 899};
  const UInt_t  fMarGam[NGam]     = {20,  22,  21};
  const UInt_t  fMarPi0[NPi0]     = {20,  22};

  // misc parameters
  const Bool_t doIntNormGam[NGam] = {false, false, false};
  const Bool_t doIntNormPi0[NPi0] = {true,  true};
  const Bool_t doRebinGam[NGam]   = {false, false, false};
  const Bool_t doRebinPi0[NPi0]   = {true,  true};
  const UInt_t nRebinGam[NGam]    = {2, 2, 2};
  const UInt_t nRebinPi0[NPi0]    = {2, 2};

  // open files
  TFile *fOut = new TFile(sOutput.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }

  TFile *fGam[NGam];
  TFile *fPi0[NPi0];
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    fGam[iGam] = new TFile(sInGam[iGam].Data(),  "read");
    if (!fGam[iGam]) {
      cerr << "PANIC: couldn't open input gamma file #" << iGam << "!" << endl;
      return;
    }
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    fPi0[iPi0] = new TFile(sInPi0[iPi0].Data(),  "read");
    if (!fPi0[iPi0]) {
      cerr << "PANIC: couldn't open input pi0 file #" << iPi0 << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hGam[NGam];
  TH1D *hPi0[NPi0];
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hGam[iGam] = (TH1D*) fGam[iGam] -> Get(sHistGam[iGam].Data());
    if (!hGam[iGam]) {
      cerr << "PANIC: couldn't grab input gamma histogram #" << iGam << "!" << endl;
      return;
    }
    hGam[iGam] -> SetName(sNameGam[iGam].Data());
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hPi0[iPi0] = (TH1D*) fPi0[iPi0] -> Get(sHistPi0[iPi0].Data());
    if (!hPi0[iPi0]) {
      cerr << "PANIC: couldn't grab an input pi0 histogram #" << iPi0 << "!" << endl;
      return;
    }
    hPi0[iPi0] -> SetName(sNamePi0[iPi0].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  Bool_t didRebin(false);
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    if (doRebinGam[iGam]) {
      hGam[iGam] -> Rebin(nRebinGam[iGam]);
      didRebin = true;
    }
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    if (doRebinPi0[iPi0]) {
      hPi0[iPi0] -> Rebin(nRebinPi0[iPi0]);
      didRebin = true;
    }
  }
  if (didRebin) cout << "    Rebinned histograms." << endl;

  // normalize histograms
  Bool_t   didNorm(false);
  Double_t integral(1.);
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    if (doIntNormGam[iGam]) {
      didNorm  = true;
      integral = hGam[iGam] -> Integral();
      if (integral > 0.) hGam[iGam] -> Scale(1. / integral);
    }
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    if (doIntNormPi0[iPi0]) {
      didNorm  = true;
      integral = hPi0[iPi0] -> Integral();
      if (integral > 0.) hPi0[iPi0] -> Scale(1. / integral);
    }
  }
  if (didNorm) cout << "    Normalized histograms." << endl;

  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fFil(0);
  const Float_t fLbl(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  const Float_t fOffL(0.07);
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hGam[iGam] -> SetMarkerColor(fColGam[iGam]);
    hGam[iGam] -> SetMarkerStyle(fMarGam[iGam]);
    hGam[iGam] -> SetLineColor(fColGam[iGam]);
    hGam[iGam] -> SetLineStyle(fLin);
    hGam[iGam] -> SetFillColor(fColGam[iGam]);
    hGam[iGam] -> SetFillStyle(fFil);
    hGam[iGam] -> SetTitle(sTitleGam.Data());
    hGam[iGam] -> SetTitleFont(fTxt);
    hGam[iGam] -> GetXaxis() -> SetRangeUser(xGamRange[0], xGamRange[1]);
    hGam[iGam] -> GetXaxis() -> SetLabelFont(fTxt);
    hGam[iGam] -> GetXaxis() -> SetLabelSize(fLbl);
    hGam[iGam] -> GetXaxis() -> SetTitle(sTitleGamX.Data());
    hGam[iGam] -> GetXaxis() -> SetTitleFont(fTxt);
    hGam[iGam] -> GetXaxis() -> SetTitleOffset(fOffX);
    hGam[iGam] -> GetXaxis() -> CenterTitle(fCnt);
    hGam[iGam] -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
    hGam[iGam] -> GetYaxis() -> SetLabelFont(fTxt);
    hGam[iGam] -> GetYaxis() -> SetLabelSize(fLbl);
    hGam[iGam] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hGam[iGam] -> GetYaxis() -> SetTitleFont(fTxt);
    hGam[iGam] -> GetYaxis() -> SetTitleOffset(fOffY);
    hGam[iGam] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hPi0[iPi0] -> SetMarkerColor(fColPi0[iPi0]);
    hPi0[iPi0] -> SetMarkerStyle(fMarPi0[iPi0]);
    hPi0[iPi0] -> SetLineColor(fColPi0[iPi0]);
    hPi0[iPi0] -> SetLineStyle(fLin);
    hPi0[iPi0] -> SetFillColor(fColPi0[iPi0]);
    hPi0[iPi0] -> SetFillStyle(fFil);
    hPi0[iPi0] -> SetTitle(sTitlePi0.Data());
    hPi0[iPi0] -> SetTitleFont(fTxt);
    hPi0[iPi0] -> GetXaxis() -> SetRangeUser(xPi0Range[0], xPi0Range[1]);
    hPi0[iPi0] -> GetXaxis() -> SetLabelFont(fTxt);
    hPi0[iPi0] -> GetXaxis() -> SetLabelSize(fLbl);
    hPi0[iPi0] -> GetXaxis() -> SetTitle(sTitlePi0X.Data());
    hPi0[iPi0] -> GetXaxis() -> SetTitleFont(fTxt);
    hPi0[iPi0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPi0[iPi0] -> GetXaxis() -> CenterTitle(fCnt);
    hPi0[iPi0] -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
    hPi0[iPi0] -> GetYaxis() -> SetLabelFont(fTxt);
    hPi0[iPi0] -> GetYaxis() -> SetLabelSize(fLbl);
    hPi0[iPi0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hPi0[iPi0] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPi0[iPi0] -> GetYaxis() -> SetTitleFont(fTxt);
    hPi0[iPi0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPi0[iPi0] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legends
  const UInt_t  fColLeg      = 0;
  const UInt_t  fLinLeg      = 0;
  const UInt_t  fAlnLeg      = 12;
  const Float_t hObj         = 0.05;
  const Float_t hObjG        = hObj * NGam;
  const Float_t hObjP        = hObj * NPi0;
  const Float_t yObjG        = 0.1 + hObjG;
  const Float_t yObjP        = 0.1 + hObjP;
  const Float_t xyLegG[NVtx] = {0.1, 0.1, 0.3, yObjG};
  const Float_t xyLegP[NVtx] = {0.1, 0.1, 0.3, yObjP};

  TLegend *legGam = new TLegend(xyLegG[0], xyLegG[1], xyLegG[2], xyLegG[3]);
  legGam -> SetLineColor(fColLeg);
  legGam -> SetLineStyle(fLinLeg);
  legGam -> SetFillColor(fColLeg);
  legGam -> SetLineColor(fLinLeg);
  legGam -> SetTextFont(fTxt);
  legGam -> SetTextAlign(fAlnLeg);
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    legGam -> AddEntry(hGam[iGam], sLegGam[iGam].Data(), "pf");
  }

  TLegend *legPi0 = new TLegend(xyLegP[0], xyLegP[1], xyLegP[2], xyLegP[3]);
  legPi0 -> SetLineColor(fColLeg);
  legPi0 -> SetLineStyle(fLinLeg);
  legPi0 -> SetFillColor(fColLeg);
  legPi0 -> SetLineColor(fLinLeg);
  legPi0 -> SetTextFont(fTxt);
  legPi0 -> SetTextAlign(fAlnLeg);
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    legPi0 -> AddEntry(hPi0[iPi0], sLegPi0[iPi0].Data(), "pf");
  }
  cout << "    Made legends." << endl;

  // make lines
  const UInt_t  fColOne(923);
  const UInt_t  fLinOne(9);
  const UInt_t  fWidOne(2);
  const Float_t yLine(1.);

  TLine *liGam = new TLine(xGamRange[0], yLine, xGamRange[1], yLine);
  TLine *liPi0 = new TLine(xPi0Range[0], yLine, xPi0Range[1], yLine);
  liGam -> SetLineColor(fColOne);
  liGam -> SetLineStyle(fLinOne);
  liGam -> SetLineWidth(fWidOne);
  liPi0 -> SetLineColor(fColOne);
  liPi0 -> SetLineStyle(fLinOne);
  liPi0 -> SetLineWidth(fWidOne);
  cout << "    Made lines." << endl;

  // make plot
  const UInt_t  width(1900);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const Float_t fMarginRG(0.015);
  const Float_t fMarginRP(0.1);
  const Float_t fMarginLG(0.1);
  const Float_t fMarginLP(0.015);
  const Float_t xyPadGam[NVtx] = {0.,  0., 0.5, 1.};
  const Float_t xyPadPi0[NVtx] = {0.5, 0., 1.,  1.};

  TCanvas *cWeight = new TCanvas("cWeights", "", width, height);
  TPad    *pGam    = new TPad("pGam", "", xyPadGam[0], xyPadGam[1], xyPadGam[2], xyPadGam[3]);
  TPad    *pPi0    = new TPad("pPi0", "", xyPadPi0[0], xyPadPi0[1], xyPadPi0[2], xyPadPi0[3]);
  cWeight -> SetGrid(fGrid, fGrid);
  cWeight -> SetTicks(fTick, fTick);
  cWeight -> SetBorderMode(fMode);
  cWeight -> SetBorderSize(fBord);
  pGam    -> SetGrid(fGrid, fGrid);
  pGam    -> SetTicks(fTick, fTick);
  pGam    -> SetBorderMode(fMode);
  pGam    -> SetBorderSize(fBord);
  pGam    -> SetLeftMargin(fMarginLG);
  pGam    -> SetRightMargin(fMarginRG);
  pGam    -> SetLogx(fLogX);
  pGam    -> SetLogy(fLogY);
  pPi0    -> SetGrid(fGrid, fGrid);
  pPi0    -> SetTicks(fTick, fTick);
  pPi0    -> SetBorderMode(fMode);
  pPi0    -> SetBorderSize(fBord);
  pPi0    -> SetLeftMargin(fMarginLP);
  pPi0    -> SetRightMargin(fMarginRP);
  pPi0    -> SetLogx(fLogX);
  pPi0    -> SetLogy(fLogY);
  cWeight -> cd();
  pGam    -> Draw();
  pPi0    -> Draw();
  pGam    -> cd();
  hGam[0] -> Draw();
  for (UInt_t iGam = 1; iGam < NGam; iGam++) {
    hGam[iGam] -> Draw("same");
  }
  liGam   -> Draw();
  legGam  -> Draw();
  pPi0    -> cd();
  hPi0[0] -> Draw();
  for (UInt_t iPi0 = 1; iPi0 < NPi0; iPi0++) {
    hPi0[iPi0] -> Draw("same");
  }
  liPi0   -> Draw();
  legPi0  -> Draw();
  fOut    -> cd();
  cWeight -> Write();
  cWeight -> Close();
  cout << "    Made plot." << endl;

  // save input
  fOut -> cd();
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hGam[iGam] -> Write();
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hPi0[iPi0] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    fGam[iGam] -> cd();
    fGam[iGam] -> Close();
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    fPi0[iPi0] -> cd();
    fPi0[iPi0] -> Close();
  }
  cout << "  Finished smearing weight plot maker!\n" << endl;

}

// End ------------------------------------------------------------------------
