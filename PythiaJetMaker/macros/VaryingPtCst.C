// 'VaryingPtCst.C'
// Derek Anderson
//
// Use this to compare jet pT spectra after varying pTcst

#include <stdio>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// constants
static const TString sFileO("pythia23.varyingPtCst.d19m1y2017.root");
static const TString sFile1("Pythia23p.r03a02rm1.pChrg.root");
static const TString sFile2("Pythia23p.r03a02rm1p07.pChrg.root");
static const TString sFile3("Pythia23p.r03a02rm1p12.pChrg.root");
static const TString sFileD("../QuarkMatterPlots/data/pp200corr.r03rm1chrg.d16m1y2017.differentPlots.root");
static const TString sNameP("QA/hPtRE");
static const TString sNameD("Pi0/hJetPtCorrP");
static const TString sNameX("p_{T,jet}^{reco,ch}");
static const TString sNameY("(1/N_{trg}) dN_{jet}/(dp_{T,jet}^{reco,ch} d#eta_{jet})");
static const TString sTitle("");
static const TString sPtCst1("0.2 GeV/c");
static const TString sPtCst2("0.7 GeV/c");
static const TString sPtCst3("1.2 GeV/c");
static const TString sPtCstD("0.2 GeV/c");
static const TString sJetType("Charged Jets");



void VaryingPtCst() {

  gErrorIgnoreLevel = kFatal;
  gStyle -> SetOptStat(0);


  TFile *fO = new TFile(sFileO.Data(), "recreate");
  TFile *f1 = new TFile(sFile1.Data(), "read");
  TFile *f2 = new TFile(sFile2.Data(), "read");
  TFile *f3 = new TFile(sFile3.Data(), "read");
  TFile *fD = new TFile(sFileD.Data(), "read");

  TH1D *h1, *h2, *h3, *hD;
  if (f1) h1 = (TH1D*) f1 -> Get(sNameP.Data());
  if (f2) h2 = (TH1D*) f2 -> Get(sNameP.Data());
  if (f3) h3 = (TH1D*) f3 -> Get(sNameP.Data());
  if (fD) hD = (TH1D*) fD -> Get(sNameD.Data());


  const Double_t xMin = -2.;
  const Double_t xMax = 20.;
  const Double_t yMin = 0.00005;
  const Double_t yMax = 5.;
  h1 -> GetXaxis() -> SetRangeUser(xMin, xMax);
  h1 -> GetXaxis() -> SetTitle(sNameX);
  h1 -> GetXaxis() -> CenterTitle(true);
  h1 -> GetYaxis() -> SetRangeUser(yMin, yMax);
  h1 -> GetYaxis() -> SetTitle(sNameY);
  h1 -> GetYaxis() -> CenterTitle(true);
  h1 -> SetTitle(sTitle);
  h1 -> SetLineColor(kRed);
  h1 -> SetMarkerColor(kRed);
  h2 -> GetXaxis() -> SetRangeUser(xMin, xMax);
  h2 -> GetXaxis() -> SetTitle(sNameX);
  h2 -> GetXaxis() -> CenterTitle(true);
  h2 -> GetYaxis() -> SetRangeUser(yMin, yMax);
  h2 -> GetYaxis() -> SetRangeUser(yMin, yMax);
  h2 -> GetYaxis() -> SetTitle(sNameY);
  h2 -> GetYaxis() -> CenterTitle(true);
  h2 -> SetTitle(sTitle);
  h2 -> SetLineColor(kGreen);
  h2 -> SetMarkerColor(kGreen);
  h3 -> GetXaxis() -> SetRangeUser(xMin, xMax);
  h3 -> GetXaxis() -> SetTitle(sNameX);
  h3 -> GetXaxis() -> CenterTitle(true);
  h3 -> GetYaxis() -> SetRangeUser(yMin, yMax);
  h3 -> GetYaxis() -> SetRangeUser(yMin, yMax);
  h3 -> GetYaxis() -> SetTitle(sNameY);
  h3 -> GetYaxis() -> CenterTitle(true);
  h3 -> SetTitle(sTitle);
  h3 -> SetLineColor(kBlue);
  h3 -> SetMarkerColor(kBlue);
  hD -> GetXaxis() -> SetRangeUser(xMin, xMax);
  hD -> GetXaxis() -> SetTitle(sNameX);
  hD -> GetXaxis() -> CenterTitle(true);
  hD -> GetYaxis() -> SetRangeUser(yMin, yMax);
  hD -> GetYaxis() -> SetRangeUser(yMin, yMax);
  hD -> GetYaxis() -> SetTitle(sNameY);
  hD -> GetYaxis() -> CenterTitle(true);
  hD -> SetTitle(sTitle);
  hD -> SetLineColor(kBlack);
  hD -> SetMarkerColor(kBlack);


  TString entry1("Pythia, p_{T}^{cst} > ");
  TString entry2("Pythia, p_{T}^{cst} > ");
  TString entry3("Pythia, p_{T}^{cst} > ");
  TString entryD("Data, p_{T}^{cst} > ");
  entry1 += sPtCst1;
  entry2 += sPtCst2;
  entry3 += sPtCst3;
  entryD += sPtCstD;
  TLegend *lg = new TLegend(0.1, 0.1, 0.3, 0.3);
  lg -> SetFillColor(kWhite);
  lg -> SetLineColor(kWhite);
  lg -> AddEntry(h1, entry1.Data());
  lg -> AddEntry(h2, entry2.Data());
  lg -> AddEntry(h3, entry3.Data());
  lg -> AddEntry(hD, entryD.Data());

  TPaveText *pt = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pt -> SetFillColor(kWhite);
  pt -> SetLineColor(kWhite);
  pt -> SetTextColor(kRed);
  pt -> AddText(sJetType);


  fO -> cd();
  TCanvas *cVaryPtCst = new TCanvas("cVaryPtCst", "Recoil jet pTreco after varying pTcst", 800, 800);
  cVaryPtCst -> cd();
  cVaryPtCst -> SetGrid(0, 0);
  cVaryPtCst -> SetLogy(1);
  cVaryPtCst -> SetTickx(1);
  cVaryPtCst -> SetTicky(1);
  h1         -> Draw();
  h2         -> Draw("same");
  h3         -> Draw("same");
  hD         -> Draw("same");
  lg         -> Draw();
  pt         -> Draw();
  cVaryPtCst -> Write();
  cVaryPtCst -> Close();


  h1 -> Write();
  h2 -> Write();
  h3 -> Write();
  hD -> Write();
  fO -> Close();

  f1 -> cd();
  f1 -> Close();
  f2 -> cd();
  f2 -> Close();
  f3 -> cd();
  f3 -> Close();
  fD -> cd();
  fD -> Close();

}

// End ------------------------------------------------------------------------
