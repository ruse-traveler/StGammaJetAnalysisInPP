// 'MakeTriggerEnergyPlotForPaper.C'
// Derek Anderson
// 02.24.2023
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



void MakeTriggerEnergyPlotForPaper() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning trigger energy plot maker..." << endl;

  // io parameters
  const TString sOutput("triggerEnergyPlotForPaper.et630pi0vsGam.d24m2y2023.root");
  const TString sInParG("trigEtComparison.forPaper_withShapeWeights.et650gam.d29m11y2022.root");
  const TString sInParP("FromSaskia/TriggerET_625_Pythia6_Combine_et6_and8.root");
  const TString sHistParG("hEtTrgParAllWeight");
  const TString sHistParP("hTrgSumProb");
  const TString sNameParG("hParGam");
  const TString sNameParP("hParPi0");
  const TString sInMatchG[NGam]   = {"trigEtComparison.forPaper_withShapeWeights.et650gam.d29m11y2022.root",
                                     "trigEtComparison.forPaper_withShapeWeights.et650gam.d29m11y2022.root",
                                     "trigEtComparison.forPaper_withShapeWeights.et650gam.d29m11y2022.root"};
  const TString sInMatchP[NPi0]   = {"FromSaskia/TriggerET_smearedPythia.root",
                                     "FromSaskia/TriggerET_smearedPythia.root"};
  const TString sHistMatchG[NGam] = {"hEtMatch911",
                                     "hEtMatch1115",
                                     "hEtMatch1520"};
  const TString sHistMatchP[NPi0] = {"hETMatch911",
                                     "hETMatch1115"};
  const TString sNameMatchG[NGam] = {"hMatchGamEt911",
                                     "hMatchGamEt1115",
                                     "hMatchGamEt1520"};
  const TString sNameMatchP[NPi0] = {"hMatchPi0Et911",
                                     "hMatchPi0Et1115"};

  // text parameters
  const TString sTitleY("Probability Density");
  const TString sTitleGam("");
  const TString sTitlePi0("");
  const TString sTitleGamX("E_{T}^{#gamma,part}, E_{T}^{#gamma,part-match} [GeV]");
  const TString sTitlePi0X("E_{T}^{#gamma,part}, E_{T}^{#gamma,part-match} [GeV]");
  const TString sLegPar("Generated");
  const TString sLegMatch("Smeared");
  const TString sLegGam[NGam] = {"9 < E_{T}^{#gamma,det-clust} < 11 GeV",
                                 "11 < E_{T}^{#gamma,det-clust} < 15 GeV",
                                 "15 < E_{T}^{#gamma,det-clust} < 20 GeV"};
  const TString sLegPi0[NPi0] = {"9 < E_{T}^{#pi^{0},det-clust} < 11 GeV",
                                 "11 < E_{T}^{#pi^{0},det-clust} < 15 GeV"};

  // plot parameters
  const UInt_t  fColHistL(923);
  const UInt_t  fColParG(923);
  const UInt_t  fColParP(923);
  const UInt_t  fMarParG(29);
  const UInt_t  fMarParP(29);
  const UInt_t  fColGam[NGam]     = {799, 899, 859};
  const UInt_t  fColPi0[NPi0]     = {859, 899};
  const UInt_t  fMarGam[NGam]     = {20,  22,  21};
  const UInt_t  fMarPi0[NPi0]     = {20,  22};
  const Float_t yRange[NRange]    = {0.000000007, 3333.};
  const Float_t xGamRange[NRange] = {6., 30.};
  const Float_t xPi0Range[NRange] = {6., 25.};

  // normalization parameters
  const Float_t xNormRangeG[NGam][NRange] = {{9., 11.}, {11., 15.}, {15., 20.}};
  const Float_t xNormRangeP[NPi0][NRange] = {{9., 11.}, {11., 15.}};

  // misc parameters
  const Bool_t doRebinParG(false);
  const Bool_t doRebinParP(false);
  const Bool_t doRebinMatchG[NGam] = {false, false, false};
  const Bool_t doRebinMatchP[NPi0] = {false, false};
  const UInt_t nRebinMatchG[NGam]  = {2, 2, 2};
  const UInt_t nRebinMatchP[NPi0]  = {2, 2};
  const UInt_t nRebinParG(2);
  const UInt_t nRebinParP(2);

  // open files
  TFile *fOut = new TFile(sOutput.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }

  TFile *fParG = new TFile(sInParG.Data(), "read");
  TFile *fParP = new TFile(sInParP.Data(), "read");
  if (!fParG || !fParP) {
    cerr << "PANIC: couldn't open input particle file!\n"
         << "       fParG = " << fParG << ", fParP = " << fParP
         << endl;
    return;
  }

  TFile *fMatchG[NGam];
  TFile *fMatchP[NPi0];
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    fMatchG[iGam] = new TFile(sInMatchG[iGam].Data(),  "read");
    if (!fMatchG[iGam]) {
      cerr << "PANIC: couldn't open input matched gamma file #" << iGam << "!" << endl;
      return;
    }
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    fMatchP[iPi0] = new TFile(sInMatchP[iPi0].Data(),  "read");
    if (!fMatchP[iPi0]) {
      cerr << "PANIC: couldn't open input matched pi0 file #" << iPi0 << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hParG = (TH1D*) fParG -> Get(sHistParG.Data());
  TH1D *hParP = (TH1D*) fParP -> Get(sHistParP.Data());
  if (!hParG || !hParP) {
    cerr << "PANIC: couldn't grab input particle histogram!\n"
         << "       hParG = " << hParG << ", hParP = " << hParP
         << endl;
    return;
  }
  hParG -> SetName(sNameParG.Data());
  hParP -> SetName(sNameParP.Data());


  TH1D *hMatchG[NGam];
  TH1D *hMatchP[NPi0];
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hMatchG[iGam] = (TH1D*) fMatchG[iGam] -> Get(sHistMatchG[iGam].Data());
    if (!hMatchG[iGam]) {
      cerr << "PANIC: couldn't grab input matched gamma histogram #" << iGam << "!" << endl;
      return;
    }
    hMatchG[iGam] -> SetName(sNameMatchG[iGam].Data());
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hMatchP[iPi0] = (TH1D*) fMatchP[iPi0] -> Get(sHistMatchP[iPi0].Data());
    if (!hMatchP[iPi0]) {
      cerr << "PANIC: couldn't grab an input matched pi0 histogram #" << iPi0 << "!" << endl;
      return;
    }
    hMatchP[iPi0] -> SetName(sNameMatchP[iPi0].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  Bool_t didRebin(false);
  if (doRebinParG) {
    hParG -> Rebin(nRebinParG);
    didRebin = true;
  }
  if (doRebinParP) {
    hParP -> Rebin(nRebinParP);
    didRebin = true;
  }

  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    if (doRebinMatchG[iGam]) {
      hMatchG[iGam] -> Rebin(nRebinMatchG[iGam]);
      didRebin = true;
    }
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    if (doRebinMatchP[iPi0]) {
      hMatchP[iPi0] -> Rebin(nRebinMatchP[iPi0]);
      didRebin = true;
    }
  }
  if (didRebin) cout << "    Rebinned histograms." << endl;

  // normalize histograms
  const Double_t intParAllG = hParG -> Integral();
  const Double_t intParAllP = hParP -> Integral();
  cout << "    Grabbing integrals:" << endl;

  Double_t intMatchG[NGam];
  Double_t intMatchP[NPi0];
  Double_t intParSelG[NGam];
  Double_t intParSelP[NPi0];
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {

    // locate bins
    const UInt_t iNormStart = hParG -> FindBin(xNormRangeG[iGam][0]);
    const UInt_t iNormStop  = hParG -> FindBin(xNormRangeG[iGam][1]);

    // grab integrals
    intMatchG[iGam]  = hMatchG[iGam] -> Integral();
    intParSelG[iGam] = hParG         -> Integral(iNormStart, iNormStop);
    cout << "      Gamma Bin [" << iGam << "]: match integral = " << intMatchG[iGam] << ", particle integral = " << intParSelG[iGam] << endl;
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {

    // locate bins
    const UInt_t iNormStart = hParP -> FindBin(xNormRangeP[iPi0][0]);
    const UInt_t iNormStop  = hParP -> FindBin(xNormRangeP[iPi0][1]);

    // grab integrals
    intMatchP[iPi0]  = hMatchP[iPi0] -> Integral();
    intParSelP[iPi0] = hParP         -> Integral(iNormStart, iNormStop);
    cout << "      Pi0 Bin [" << iPi0 << "]: match integral = " << intMatchP[iPi0] << ", particle integral = " << intParSelP[iPi0] << endl;
  }

  // do normalization
  if (intParAllG > 0.) {
    hParG -> Scale(1. / intParAllG);
    for (UInt_t iGam = 0; iGam < NGam; iGam++) {
      const Double_t intSelect = (intParSelG[iGam] / intParAllG);
      const Double_t normScale = (intSelect / intMatchG[iGam]);
      hMatchG[iGam] -> Scale(normScale);
    }
  }
  if (intParAllP > 0.) {
    hParP -> Scale(1. / intParAllP);
    for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
      const Double_t intSelect = (intParSelP[iPi0] / intParAllP);
      const Double_t normScale = (intSelect / intMatchP[iPi0]);
      hMatchP[iPi0] -> Scale(normScale);
    }
  }
  cout << "    Normalized histograms." << endl;

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
  hParG -> SetMarkerColor(fColParG);
  hParG -> SetMarkerStyle(fMarParG);
  hParG -> SetLineColor(fColParG);
  hParG -> SetLineStyle(fLin);
  hParG -> SetFillColor(fColParG);
  hParG -> SetFillStyle(fFil);
  hParG -> SetTitle(sTitleGam.Data());
  hParG -> SetTitleFont(fTxt);
  hParG -> GetXaxis() -> SetRangeUser(xGamRange[0], xGamRange[1]);
  hParG -> GetXaxis() -> SetLabelFont(fTxt);
  hParG -> GetXaxis() -> SetLabelSize(fLbl);
  hParG -> GetXaxis() -> SetTitle(sTitleGamX.Data());
  hParG -> GetXaxis() -> SetTitleFont(fTxt);
  hParG -> GetXaxis() -> SetTitleOffset(fOffX);
  hParG -> GetXaxis() -> CenterTitle(fCnt);
  hParG -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
  hParG -> GetYaxis() -> SetLabelFont(fTxt);
  hParG -> GetYaxis() -> SetLabelSize(fLbl);
  hParG -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParG -> GetYaxis() -> SetTitleFont(fTxt);
  hParG -> GetYaxis() -> SetTitleOffset(fOffY);
  hParG -> GetYaxis() -> CenterTitle(fCnt);
  hParP -> SetMarkerColor(fColParP);
  hParP -> SetMarkerStyle(fMarParP);
  hParP -> SetLineColor(fColParP);
  hParP -> SetLineStyle(fLin);
  hParP -> SetFillColor(fColParP);
  hParP -> SetFillStyle(fFil);
  hParP -> SetTitle(sTitlePi0.Data());
  hParP -> SetTitleFont(fTxt);
  hParP -> GetXaxis() -> SetRangeUser(xPi0Range[0], xPi0Range[1]);
  hParP -> GetXaxis() -> SetLabelFont(fTxt);
  hParP -> GetXaxis() -> SetLabelSize(fLbl);
  hParP -> GetXaxis() -> SetTitle(sTitlePi0X.Data());
  hParP -> GetXaxis() -> SetTitleFont(fTxt);
  hParP -> GetXaxis() -> SetTitleOffset(fOffX);
  hParP -> GetXaxis() -> CenterTitle(fCnt);
  hParP -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
  hParP -> GetYaxis() -> SetLabelFont(fTxt);
  hParP -> GetYaxis() -> SetLabelSize(fLbl);
  hParP -> GetYaxis() -> SetLabelOffset(fOffL);
  hParP -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParP -> GetYaxis() -> SetTitleFont(fTxt);
  hParP -> GetYaxis() -> SetTitleOffset(fOffY);
  hParP -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hMatchG[iGam] -> SetMarkerColor(fColGam[iGam]);
    hMatchG[iGam] -> SetMarkerStyle(fMarGam[iGam]);
    hMatchG[iGam] -> SetLineColor(fColGam[iGam]);
    hMatchG[iGam] -> SetLineStyle(fLin);
    hMatchG[iGam] -> SetFillColor(fColGam[iGam]);
    hMatchG[iGam] -> SetFillStyle(fFil);
    hMatchG[iGam] -> SetTitle(sTitleGam.Data());
    hMatchG[iGam] -> SetTitleFont(fTxt);
    hMatchG[iGam] -> GetXaxis() -> SetRangeUser(xGamRange[0], xGamRange[1]);
    hMatchG[iGam] -> GetXaxis() -> SetLabelFont(fTxt);
    hMatchG[iGam] -> GetXaxis() -> SetLabelSize(fLbl);
    hMatchG[iGam] -> GetXaxis() -> SetTitle(sTitleGamX.Data());
    hMatchG[iGam] -> GetXaxis() -> SetTitleFont(fTxt);
    hMatchG[iGam] -> GetXaxis() -> SetTitleOffset(fOffX);
    hMatchG[iGam] -> GetXaxis() -> CenterTitle(fCnt);
    hMatchG[iGam] -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
    hMatchG[iGam] -> GetYaxis() -> SetLabelFont(fTxt);
    hMatchG[iGam] -> GetYaxis() -> SetLabelSize(fLbl);
    hMatchG[iGam] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hMatchG[iGam] -> GetYaxis() -> SetTitleFont(fTxt);
    hMatchG[iGam] -> GetYaxis() -> SetTitleOffset(fOffY);
    hMatchG[iGam] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hMatchP[iPi0] -> SetMarkerColor(fColPi0[iPi0]);
    hMatchP[iPi0] -> SetMarkerStyle(fMarPi0[iPi0]);
    hMatchP[iPi0] -> SetLineColor(fColPi0[iPi0]);
    hMatchP[iPi0] -> SetLineStyle(fLin);
    hMatchP[iPi0] -> SetFillColor(fColPi0[iPi0]);
    hMatchP[iPi0] -> SetFillStyle(fFil);
    hMatchP[iPi0] -> SetTitle(sTitlePi0.Data());
    hMatchP[iPi0] -> SetTitleFont(fTxt);
    hMatchP[iPi0] -> GetXaxis() -> SetRangeUser(xPi0Range[0], xPi0Range[1]);
    hMatchP[iPi0] -> GetXaxis() -> SetLabelFont(fTxt);
    hMatchP[iPi0] -> GetXaxis() -> SetLabelSize(fLbl);
    hMatchP[iPi0] -> GetXaxis() -> SetTitle(sTitlePi0X.Data());
    hMatchP[iPi0] -> GetXaxis() -> SetTitleFont(fTxt);
    hMatchP[iPi0] -> GetXaxis() -> SetTitleOffset(fOffX);
    hMatchP[iPi0] -> GetXaxis() -> CenterTitle(fCnt);
    hMatchP[iPi0] -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
    hMatchP[iPi0] -> GetYaxis() -> SetLabelFont(fTxt);
    hMatchP[iPi0] -> GetYaxis() -> SetLabelSize(fLbl);
    hMatchP[iPi0] -> GetYaxis() -> SetLabelOffset(fOffL);
    hMatchP[iPi0] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hMatchP[iPi0] -> GetYaxis() -> SetTitleFont(fTxt);
    hMatchP[iPi0] -> GetYaxis() -> SetTitleOffset(fOffY);
    hMatchP[iPi0] -> GetYaxis() -> CenterTitle(fCnt);
  }

  // make histogram for legend
  hParL   = (TH1D*) hParG      -> Clone();
  hMatchL = (TH1D*) hMatchG[0] -> Clone();
  hParL   -> SetName("hParLegend");
  hMatchL -> SetName("hMatchLegend");
  hParL   -> SetMarkerColor(fColHistL);
  hMatchL -> SetMarkerColor(fColHistL);
  hParL   -> SetFillColor(fColHistL);
  hMatchL -> SetFillColor(fColHistL);
  hParL   -> SetLineColor(fColHistL);
  hMatchL -> SetLineColor(fColHistL);
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLeg(0);
  const UInt_t  fLinLeg(0);
  const UInt_t  fColLin(1);
  const UInt_t  fAlnLeg(12);
  const UInt_t  nColumn(3);
  const Float_t xyLeg[NVtx] = {0.05, 0.015, 0.95, 0.9};

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]); 
  leg -> SetFillColor(fColLeg);
  leg -> SetLineColor(fColLin);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnLeg);
  leg -> SetNColumns(nColumn);
  leg -> AddEntry(hMatchG[0], sLegGam[0].Data(), "pf");
  leg -> AddEntry(hMatchP[0], sLegPi0[0].Data(), "pf");
  leg -> AddEntry(hParL,      sLegPar.Data(),    "pf");
  leg -> AddEntry(hMatchG[1], sLegGam[1].Data(), "pf");
  leg -> AddEntry(hMatchP[1], sLegPi0[1].Data(), "pf");
  leg -> AddEntry(hMatchL,    sLegMatch.Data(),  "pf");
  leg -> AddEntry(hMatchG[2], sLegGam[2].Data(), "pf");
  cout << "    Made legend." << endl;

  // make plot
  const UInt_t  width(1900);
  const UInt_t  height(1056);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fGridL(0);
  const UInt_t  fTick(1);
  const UInt_t  fTickL(0);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fLogYL(0);
  const Float_t fMarginTG(0.015);
  const Float_t fMarginTP(0.015);
  const Float_t fMarginTL(0.015);
  const Float_t fMarginBG(0.1);
  const Float_t fMarginBP(0.1);
  const Float_t fMarginBL(0.015);
  const Float_t fMarginRG(0.015);
  const Float_t fMarginRP(0.1);
  const Float_t fMarginRL(0.1);
  const Float_t fMarginLG(0.1);
  const Float_t fMarginLP(0.015);
  const Float_t fMarginLL(0.1);
  const Float_t xyPadGam[NVtx] = {0.,  0.,  0.5, 0.9};
  const Float_t xyPadPi0[NVtx] = {0.5, 0.,  1.,  0.9};
  const Float_t xyPadLeg[NVtx] = {0.,  0.9, 1.,  1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pGam    = new TPad("pGam", "", xyPadGam[0], xyPadGam[1], xyPadGam[2], xyPadGam[3]);
  TPad    *pPi0    = new TPad("pPi0", "", xyPadPi0[0], xyPadPi0[1], xyPadPi0[2], xyPadPi0[3]);
  TPad    *pLeg    = new TPad("pLeg", "", xyPadLeg[0], xyPadLeg[1], xyPadLeg[2], xyPadLeg[3]);
  cPlot -> SetGrid(fGrid, fGrid);
  cPlot -> SetTicks(fTick, fTick);
  cPlot -> SetBorderMode(fMode);
  cPlot -> SetBorderSize(fBord);
  pGam  -> SetGrid(fGrid, fGrid);
  pGam  -> SetTicks(fTick, fTick);
  pGam  -> SetBorderMode(fMode);
  pGam  -> SetBorderSize(fBord);
  pGam  -> SetTopMargin(fMarginTG);
  pGam  -> SetBottomMargin(fMarginBG);
  pGam  -> SetLeftMargin(fMarginLG);
  pGam  -> SetRightMargin(fMarginRG);
  pGam  -> SetLogx(fLogX);
  pGam  -> SetLogy(fLogY);
  pPi0  -> SetGrid(fGrid, fGrid);
  pPi0  -> SetTicks(fTick, fTick);
  pPi0  -> SetBorderMode(fMode);
  pPi0  -> SetBorderSize(fBord);
  pPi0  -> SetTopMargin(fMarginTP);
  pPi0  -> SetBottomMargin(fMarginBP);
  pPi0  -> SetLeftMargin(fMarginLP);
  pPi0  -> SetRightMargin(fMarginRP);
  pPi0  -> SetLogx(fLogX);
  pPi0  -> SetLogy(fLogY);
  pLeg  -> SetGrid(fGridL, fGridL);
  pLeg  -> SetTicks(fTickL, fTickL);
  pLeg  -> SetBorderMode(fMode);
  pLeg  -> SetBorderSize(fBord);
  pLeg  -> SetTopMargin(fMarginTL);
  pLeg  -> SetBottomMargin(fMarginBL);
  pLeg  -> SetLeftMargin(fMarginLL);
  pLeg  -> SetRightMargin(fMarginRL);
  pLeg  -> SetLogx(fLogX);
  pLeg  -> SetLogy(fLogYL);
  cPlot -> cd();
  pGam  -> Draw();
  pPi0  -> Draw();
  pLeg  -> Draw();
  pGam  -> cd();
  hParG -> Draw();
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hMatchG[iGam] -> Draw("same");
  }
  pPi0  -> cd();
  hParP -> Draw();
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hMatchP[iPi0] -> Draw("same");
  }
  pLeg  -> cd();
  leg   -> Draw();
  fOut  -> cd();
  cPlot -> Write();
  cPlot -> Close();
  cout << "    Made plot." << endl;

  // save input
  fOut    -> cd();
  hParG   -> Write();
  hParP   -> Write();
  hParL   -> Write();
  hMatchL -> Write();
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    hMatchG[iGam] -> Write();
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    hMatchP[iPi0] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut  -> cd();
  fOut  -> Close();
  fParG -> cd();
  fParG -> Close();
  fParP -> cd();
  fParP -> Close();
  for (UInt_t iGam = 0; iGam < NGam; iGam++) {
    fMatchG[iGam] -> cd();
    fMatchG[iGam] -> Close();
  }
  for (UInt_t iPi0 = 0; iPi0 < NPi0; iPi0++) {
    fMatchP[iPi0] -> cd();
    fMatchP[iPi0] -> Close();
  }
  cout << "  Finished trigger energy plot maker!\n" << endl;

}

// End ------------------------------------------------------------------------
