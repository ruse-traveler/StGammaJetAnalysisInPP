// 'MakeJetEfficiencyPlot.C'
// Derek Anderson
// 10.24.2021
//
// Makes a plot of the jet-matching
// efficiency with a systematic
// uncertainty


#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NRes(2);
static const UInt_t NSys(3);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NPlot(2);



void MakeJetEfficiencyPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning jet efficiency plot macro..." << endl;

  // output and denominator parameters
  const TString sOutput("jetEfficiency.pTbinUltra.et920r02x05pi0.d24m10y2021.root");
  const TString sHistDef("hEfficiencyAll");
  const TString sHistSysM("hEfficiency");
  const TString sHistSysP("hEfficiency");
  const TString sInDef[NRes]  = {"pp200r9embed.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r02a005rm1chrg.dr02qt05130.d14m10y2021.root", "pp200r9embed.forJetEffPlot_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root"};
  const TString sInSysM[NRes] = {"pp200py8.forJetEffPlotM04_pTbinUltra.et920pt0230x021Kvz55pi0.r02a005rm1chrg.dr02qt05130.d14m10y2021.root", "pp200py8.forJetEffPlotM04_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root"};
  const TString sInSysP[NRes] = {"pp200py8.forJetEffPlotP04_pTbinUltra.et920pt0230x021Kvz55pi0.r02a005rm1chrg.dr02qt05130.d14m10y2021.root", "pp200py8.forJetEffPlotP04_pTbinUltra.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root"};

  // histogram parameters
  const TString sNameDef[NRes]  = {"hEffDefR02", "hEffDefR05"};
  const TString sNameSysM[NRes] = {"hEffSysM02", "hEffSysM05"};
  const TString sNameSysP[NRes] = {"hEffSysP02", "hEffSysP05"};
  const TString sNameSysA[NRes] = {"hEffSysA02", "hEffSysA05"};
  const TString sNameErr[NRes]  = {"hSysErrR02", "hErrSysR05"};
  const TString sLabelDef[NRes] = {"R = 0.2", "R = 0.5"};
  const TString sLabelSys[NSys] = {"#epsilon_{jet}(Py8#oplusParam) [#epsilon_{trk} - 4%]", "#epsilon_{jet}(Py8#oplusParam) [#epsilon_{trk} + 4%]", "average of #epsilon_{trk} #pm 4%"};

  // canvas parameters
  const TString sCanDef("cJetEfficiency");
  const TString sCanSys[NRes] = {"cJetEffSystematicR02", "cJetEffSystematicR05"};
  const TString sPadSys[NPad] = {"pEffVar", "pEffSys"};

  // text parameters
  const TString sSys("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  const TString sTyp("#bf{charged jets}");
  const TString sJetD("anti-k_{T} algo., #bf{charged jets}");
  const TString sJetS[NRes] = {"anti-k_{T} algo., R = 0.2", "anti-k_{T} algo., R = 0.5"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{jet} [GeV/c]");
  const TString sTitleE("jet matching efficiency");
  const TString sTitleV("relative #epsilon_{jet} uncertainty");
  const Float_t xCalcRange[NPlot] = {0.2, 30.};
  const Float_t xPlotRange[NPlot] = {-1., 20.};
  const UInt_t  fColLinDef[NRes]  = {923, 893};
  const UInt_t  fColFilDef[NRes]  = {920, 623};
  const UInt_t  fColLinSys[NSys]  = {859, 899, 923};
  const UInt_t  fColFilSys[NSys]  = {859, 899, 923};
  const UInt_t  fLinDef[NRes]     = {9, 9};
  const UInt_t  fLinSys[NSys]     = {1, 1, 1};
  const UInt_t  fFilDef[NRes]     = {1001, 1001};
  const UInt_t  fFilSys[NSys]     = {0, 0, 0};
  const UInt_t  fMarDef[NRes]     = {1, 1};
  const UInt_t  fMarSys[NSys]     = {22, 23, 20};
  const UInt_t  fColE(920);
  const UInt_t  fMarE(1);
  const UInt_t  fFilE(1001);

  // open output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  // open input files
  TFile *fInDef[NRes];
  TFile *fInSysM[NRes];
  TFile *fInSysP[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    fInDef[iRes]  = new TFile(sInDef[iRes].Data(), "read");
    fInSysM[iRes] = new TFile(sInSysM[iRes].Data(), "read");
    fInSysP[iRes] = new TFile(sInSysP[iRes].Data(), "read");
    if (!fInDef[iRes] || !fInSysM[iRes] || !fInSysP[iRes]) {
      cerr << "PANIC: couldn't open an input file!\n"
           << "       fInDef[" << iRes << "] = " << fInDef[iRes] << ", fInSysM[" << iRes << "] = " << fInSysM[iRes] << ", fInSysP[" << iRes << "] = " << fInSysP[iRes] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input histograms
  TH1D *hEffDef[NRes];
  TH1D *hEffSysM[NRes];
  TH1D *hEffSysP[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    hEffDef[iRes]  = (TH1D*) fInDef[iRes]  -> Get(sHistDef.Data());
    hEffSysM[iRes] = (TH1D*) fInSysM[iRes] -> Get(sHistSysM.Data());
    hEffSysP[iRes] = (TH1D*) fInSysP[iRes] -> Get(sHistSysP.Data());
    if (!hEffDef[iRes] || !hEffSysM[iRes] || !hEffSysP[iRes]) {
      cerr << "PANIC: couldn't grab input histogram!\n"
           << "       hEffDef[" << iRes << "] = " << hEffDef[iRes] << ", hEffSysM[" << iRes << "] = " << hEffSysM[iRes] << ", hEffSysP[" << iRes << "] = " << hEffSysP[iRes] << "\n"
           << endl;
      return;
    }
    hEffDef[iRes]  -> SetName(sNameDef[iRes].Data());
    hEffSysM[iRes] -> SetName(sNameSysM[iRes].Data());
    hEffSysP[iRes] -> SetName(sNameSysP[iRes].Data());
  }
  cout << "    Grabbed input histograms." << endl;

  // calculate averages
  TH1D *hEffAvg[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    hEffAvg[iRes] = (TH1D*) hEffSysM[iRes] -> Clone();
    hEffAvg[iRes] -> SetName(sNameSysA[iRes].Data());
    hEffAvg[iRes] -> Reset("ICES");
    hEffAvg[iRes] -> Add(hEffSysM[iRes]);
    hEffAvg[iRes] -> Add(hEffSysP[iRes]);
    hEffAvg[iRes] -> Scale(0.5);
  }
  cout << "    Calculated ratios." << endl;

  // calculate uncertainties
  TH1D *hEffErr[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {

    hEffErr[iRes] = (TH1D*) hEffAvg[iRes] -> Clone();
    hEffErr[iRes] -> SetName(sNameErr[iRes].Data());
    hEffErr[iRes] -> Reset("ICES");

    const UInt_t nBins = hEffAvg[iRes] -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {

      const Double_t binCent = hEffAvg[iRes]  -> GetBinCenter(iBin);
      const Double_t binValA = hEffAvg[iRes]  -> GetBinContent(iBin);
      const Double_t binValM = hEffSysM[iRes] -> GetBinContent(iBin);
      const Double_t binValP = hEffSysP[iRes] -> GetBinContent(iBin);
      const Double_t binValD = hEffDef[iRes]  -> GetBinContent(iBin);
      const Double_t diffM   = TMath::Abs(binValM - binValA);
      const Double_t diffP   = TMath::Abs(binValP - binValA);
      const Double_t maxDiff = TMath::Max(diffM, diffP);
      const Double_t relDiff = maxDiff / binValA;
      const Double_t defDiff = binValD * relDiff;

      const Bool_t isNonzeroM = (binValM > 0.);
      const Bool_t isNonzeroP = (binValP > 0.);
      const Bool_t isNonzeroA = (binValA > 0.);
      const Bool_t isNonzeroD = (binValD > 0.);
      const Bool_t isInRange  = ((binCent > xCalcRange[0]) && (binCent < xCalcRange[1]));
      if (isNonzeroM && isNonzeroP && isNonzeroA && isNonzeroD && isInRange) {
        hEffErr[iRes] -> SetBinContent(iBin, 1.);
        hEffErr[iRes] -> SetBinError(iBin, relDiff);
        hEffDef[iRes] -> SetBinError(iBin, defDiff);
      } else {
        hEffErr[iRes]  -> SetBinContent(iBin, 1.);
        hEffDef[iRes]  -> SetBinContent(iBin, 0.);
        hEffAvg[iRes]  -> SetBinContent(iBin, 0.);
        hEffSysM[iRes] -> SetBinContent(iBin, 0.);
        hEffSysP[iRes] -> SetBinContent(iBin, 0.);
        hEffErr[iRes]  -> SetBinError(iBin, 0.);
        hEffDef[iRes]  -> SetBinError(iBin, 0.);
        hEffAvg[iRes]  -> SetBinError(iBin, 0.);
        hEffSysM[iRes] -> SetBinError(iBin, 0.);
        hEffSysP[iRes] -> SetBinError(iBin, 0.);
      }
    }  // end bin loop
  }  // end resolution loop
  cout << "    Calculated and set uncertainties." << endl;

  // set styles
  const UInt_t  fWidD(3);
  const UInt_t  fWidS(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLabD(0.04);
  const Float_t fTtlD(0.04);
  const Float_t fOffXD(1.1);
  const Float_t fOffYD(1.5);
  const Float_t fLabS[NPad]  = {0.04, 0.04};
  const Float_t fTtlS[NPad]  = {0.04, 0.04};
  const Float_t fOffXS[NPad] = {1.1, 1.1};
  const Float_t fOffYS[NPad] = {1.5, 1.5};
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    hEffDef[iRes]  -> SetMarkerColor(fColLinDef[iRes]);
    hEffDef[iRes]  -> SetMarkerStyle(fMarDef[iRes]);
    hEffDef[iRes]  -> SetFillColor(fColFilDef[iRes]);
    hEffDef[iRes]  -> SetFillStyle(fFilDef[iRes]);
    hEffDef[iRes]  -> SetLineColor(fColLinDef[iRes]);
    hEffDef[iRes]  -> SetLineStyle(fLinDef[iRes]);
    hEffDef[iRes]  -> SetLineWidth(fWidD);
    hEffDef[iRes]  -> SetTitle(sTitle.Data());
    hEffDef[iRes]  -> SetTitleFont(fTxt);
    hEffDef[iRes]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hEffDef[iRes]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hEffDef[iRes]  -> GetXaxis() -> SetTitleFont(fTxt);
    hEffDef[iRes]  -> GetXaxis() -> SetTitleSize(fTtlD);
    hEffDef[iRes]  -> GetXaxis() -> SetTitleOffset(fOffXD);
    hEffDef[iRes]  -> GetXaxis() -> SetLabelFont(fTxt);
    hEffDef[iRes]  -> GetXaxis() -> SetLabelSize(fLabD);
    hEffDef[iRes]  -> GetXaxis() -> CenterTitle(fCnt);
    hEffDef[iRes]  -> GetYaxis() -> SetTitle(sTitleE.Data());
    hEffDef[iRes]  -> GetYaxis() -> SetTitleFont(fTxt);
    hEffDef[iRes]  -> GetYaxis() -> SetTitleSize(fTtlD);
    hEffDef[iRes]  -> GetYaxis() -> SetTitleOffset(fOffYD);
    hEffDef[iRes]  -> GetYaxis() -> SetLabelFont(fTxt);
    hEffDef[iRes]  -> GetYaxis() -> SetLabelSize(fLabD);
    hEffDef[iRes]  -> GetYaxis() -> CenterTitle(fCnt);
    hEffSysM[iRes] -> SetMarkerColor(fColLinSys[0]);
    hEffSysM[iRes] -> SetMarkerStyle(fMarSys[0]);
    hEffSysM[iRes] -> SetFillColor(fColFilSys[0]);
    hEffSysM[iRes] -> SetFillStyle(fFilSys[0]);
    hEffSysM[iRes] -> SetLineColor(fColLinSys[0]);
    hEffSysM[iRes] -> SetLineStyle(fLinSys[0]);
    hEffSysM[iRes] -> SetLineWidth(fWidS);
    hEffSysM[iRes] -> SetTitle(sTitle.Data());
    hEffSysM[iRes] -> SetTitleFont(fTxt);
    hEffSysM[iRes] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hEffSysM[iRes] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hEffSysM[iRes] -> GetXaxis() -> SetTitleFont(fTxt);
    hEffSysM[iRes] -> GetXaxis() -> SetTitleSize(fTtlS[0]);
    hEffSysM[iRes] -> GetXaxis() -> SetTitleOffset(fOffXS[0]);
    hEffSysM[iRes] -> GetXaxis() -> SetLabelFont(fTxt);
    hEffSysM[iRes] -> GetXaxis() -> SetLabelSize(fLabS[0]);
    hEffSysM[iRes] -> GetXaxis() -> CenterTitle(fCnt);
    hEffSysM[iRes] -> GetYaxis() -> SetTitle(sTitleE.Data());
    hEffSysM[iRes] -> GetYaxis() -> SetTitleFont(fTxt);
    hEffSysM[iRes] -> GetYaxis() -> SetTitleSize(fTtlS[0]);
    hEffSysM[iRes] -> GetYaxis() -> SetTitleOffset(fOffYS[0]);
    hEffSysM[iRes] -> GetYaxis() -> SetLabelFont(fTxt);
    hEffSysM[iRes] -> GetYaxis() -> SetLabelSize(fLabS[0]);
    hEffSysM[iRes] -> GetYaxis() -> CenterTitle(fCnt);
    hEffSysP[iRes] -> SetMarkerColor(fColLinSys[1]);
    hEffSysP[iRes] -> SetMarkerStyle(fMarSys[1]);
    hEffSysP[iRes] -> SetFillColor(fColFilSys[1]);
    hEffSysP[iRes] -> SetFillStyle(fFilSys[1]);
    hEffSysP[iRes] -> SetLineColor(fColLinSys[1]);
    hEffSysP[iRes] -> SetLineStyle(fLinSys[1]);
    hEffSysP[iRes] -> SetLineWidth(fWidS);
    hEffSysP[iRes] -> SetTitle(sTitle.Data());
    hEffSysP[iRes] -> SetTitleFont(fTxt);
    hEffSysP[iRes] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hEffSysP[iRes] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hEffSysP[iRes] -> GetXaxis() -> SetTitleFont(fTxt);
    hEffSysP[iRes] -> GetXaxis() -> SetTitleSize(fTtlS[0]);
    hEffSysP[iRes] -> GetXaxis() -> SetTitleOffset(fOffXS[0]);
    hEffSysP[iRes] -> GetXaxis() -> SetLabelFont(fTxt);
    hEffSysP[iRes] -> GetXaxis() -> SetLabelSize(fLabS[0]);
    hEffSysP[iRes] -> GetXaxis() -> CenterTitle(fCnt);
    hEffSysP[iRes] -> GetYaxis() -> SetTitle(sTitleE.Data());
    hEffSysP[iRes] -> GetYaxis() -> SetTitleFont(fTxt);
    hEffSysP[iRes] -> GetYaxis() -> SetTitleSize(fTtlS[0]);
    hEffSysP[iRes] -> GetYaxis() -> SetTitleOffset(fOffYS[0]);
    hEffSysP[iRes] -> GetYaxis() -> SetLabelFont(fTxt);
    hEffSysP[iRes] -> GetYaxis() -> SetLabelSize(fLabS[0]);
    hEffSysP[iRes] -> GetYaxis() -> CenterTitle(fCnt);
    hEffAvg[iRes]  -> SetMarkerColor(fColLinSys[2]);
    hEffAvg[iRes]  -> SetMarkerStyle(fMarSys[2]);
    hEffAvg[iRes]  -> SetFillColor(fColFilSys[2]);
    hEffAvg[iRes]  -> SetFillStyle(fFilSys[2]);
    hEffAvg[iRes]  -> SetLineColor(fColLinSys[2]);
    hEffAvg[iRes]  -> SetLineStyle(fLinSys[2]);
    hEffAvg[iRes]  -> SetLineWidth(fWidS);
    hEffAvg[iRes]  -> SetTitle(sTitle.Data());
    hEffAvg[iRes]  -> SetTitleFont(fTxt);
    hEffAvg[iRes]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hEffAvg[iRes]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hEffAvg[iRes]  -> GetXaxis() -> SetTitleFont(fTxt);
    hEffAvg[iRes]  -> GetXaxis() -> SetTitleSize(fTtlS[0]);
    hEffAvg[iRes]  -> GetXaxis() -> SetTitleOffset(fOffXS[0]);
    hEffAvg[iRes]  -> GetXaxis() -> SetLabelFont(fTxt);
    hEffAvg[iRes]  -> GetXaxis() -> SetLabelSize(fLabS[0]);
    hEffAvg[iRes]  -> GetXaxis() -> CenterTitle(fCnt);
    hEffAvg[iRes]  -> GetYaxis() -> SetTitle(sTitleE.Data());
    hEffAvg[iRes]  -> GetYaxis() -> SetTitleFont(fTxt);
    hEffAvg[iRes]  -> GetYaxis() -> SetTitleSize(fTtlS[0]);
    hEffAvg[iRes]  -> GetYaxis() -> SetTitleOffset(fOffYS[0]);
    hEffAvg[iRes]  -> GetYaxis() -> SetLabelFont(fTxt);
    hEffAvg[iRes]  -> GetYaxis() -> SetLabelSize(fLabS[0]);
    hEffAvg[iRes]  -> GetYaxis() -> CenterTitle(fCnt);
    hEffErr[iRes]  -> SetMarkerColor(fColE);
    hEffErr[iRes]  -> SetMarkerStyle(fMarE);
    hEffErr[iRes]  -> SetFillColor(fColE);
    hEffErr[iRes]  -> SetFillStyle(fFilE);
    hEffErr[iRes]  -> SetLineColor(fColLinSys[2]);
    hEffErr[iRes]  -> SetLineStyle(fLinSys[2]);
    hEffErr[iRes]  -> SetLineWidth(fWidD);
    hEffErr[iRes]  -> SetTitle(sTitle.Data());
    hEffErr[iRes]  -> SetTitleFont(fTxt);
    hEffErr[iRes]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hEffErr[iRes]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hEffErr[iRes]  -> GetXaxis() -> SetTitleFont(fTxt);
    hEffErr[iRes]  -> GetXaxis() -> SetTitleSize(fTtlS[1]);
    hEffErr[iRes]  -> GetXaxis() -> SetTitleOffset(fOffXS[1]);
    hEffErr[iRes]  -> GetXaxis() -> SetLabelFont(fTxt);
    hEffErr[iRes]  -> GetXaxis() -> SetLabelSize(fLabS[1]);
    hEffErr[iRes]  -> GetXaxis() -> CenterTitle(fCnt);
    hEffErr[iRes]  -> GetYaxis() -> SetTitle(sTitleV.Data());
    hEffErr[iRes]  -> GetYaxis() -> SetTitleFont(fTxt);
    hEffErr[iRes]  -> GetYaxis() -> SetTitleSize(fTtlS[1]);
    hEffErr[iRes]  -> GetYaxis() -> SetTitleOffset(fOffYS[1]);
    hEffErr[iRes]  -> GetYaxis() -> SetLabelFont(fTxt);
    hEffErr[iRes]  -> GetYaxis() -> SetLabelSize(fLabS[1]);
    hEffErr[iRes]  -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make line and outline histograms
  const UInt_t  fFilOut(0);
  const UInt_t  fLinOut(1);
  const UInt_t  fWidLin(2);
  const UInt_t  fWidOut(3);
  const TString sLinSuff("_line");
  const TString sOutSuff("_outline");

  TH1D *hLinDef[NRes];
  TH1D *hLinErr[NRes];
  TH1D *hOutDef[NRes];
  TH1D *hOutErr[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {

    // create names
    TString sLinNameD = hEffDef[iRes] -> GetName();
    TString sOutNameD = hEffDef[iRes] -> GetName();
    TString sLinNameE = hEffErr[iRes] -> GetName();
    TString sOutNameE = hEffErr[iRes] -> GetName();
    sLinNameD.Append(sLinSuff.Data());
    sOutNameD.Append(sOutSuff.Data());
    sLinNameE.Append(sLinSuff.Data());
    sOutNameE.Append(sOutSuff.Data());

    // create histograms
    hLinDef[iRes] = (TH1D*) hEffDef[iRes] -> Clone();
    hOutDef[iRes] = (TH1D*) hEffDef[iRes] -> Clone();
    hLinErr[iRes] = (TH1D*) hEffErr[iRes] -> Clone();
    hOutErr[iRes] = (TH1D*) hEffErr[iRes] -> Clone();
    hLinDef[iRes] -> SetName(sLinNameD.Data());
    hOutDef[iRes] -> SetName(sOutNameD.Data());
    hLinErr[iRes] -> SetName(sLinNameE.Data());
    hOutErr[iRes] -> SetName(sOutNameE.Data());
    hLinDef[iRes] -> SetFillColor(fColLinDef[iRes]);
    hLinDef[iRes] -> SetFillStyle(fFilOut);
    hLinDef[iRes] -> SetLineWidth(fWidLin);
    hOutDef[iRes] -> SetFillColor(fColLinDef[iRes]);
    hOutDef[iRes] -> SetFillStyle(fFilOut);
    hOutDef[iRes] -> SetLineStyle(fLinOut);
    hOutDef[iRes] -> SetLineWidth(fWidOut);
    hLinErr[iRes] -> SetFillColor(fColLinSys[2]);
    hLinErr[iRes] -> SetFillStyle(fFilOut);
    hLinErr[iRes] -> SetLineWidth(fWidLin);
    hOutErr[iRes] -> SetFillColor(fColLinSys[2]);
    hOutErr[iRes] -> SetFillStyle(fFilOut);
    hOutErr[iRes] -> SetLineStyle(fLinOut);
    hOutErr[iRes] -> SetLineWidth(fWidOut);
  }
  cout << "    Made outline histograms." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXyD[NVtx] = {0.1, 0.1, 0.3, 0.2};
  const Float_t fLegXyS[NVtx] = {0.1, 0.1, 0.3, 0.25};

  TLegend *legD = new TLegend(fLegXyD[0], fLegXyD[1], fLegXyD[2], fLegXyD[3]);
  legD -> SetFillColor(fColLe);
  legD -> SetFillStyle(fFilLe);
  legD -> SetLineColor(fColLe);
  legD -> SetLineStyle(fLinLe);
  legD -> SetTextFont(fTxt);
  legD -> SetTextAlign(fAln);
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    legD -> AddEntry(hEffDef[iRes], sLabelDef[iRes], "lf");
  }

  TLegend *legS[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    legS[iRes] = new TLegend(fLegXyS[0], fLegXyS[1], fLegXyS[2], fLegXyS[3]);
    legS[iRes] -> SetFillColor(fColLe);
    legS[iRes] -> SetFillStyle(fFilLe);
    legS[iRes] -> SetLineColor(fColLe);
    legS[iRes] -> SetLineStyle(fLinLe);
    legS[iRes] -> SetTextFont(fTxt);
    legS[iRes] -> SetTextAlign(fAln);
    legS[iRes] -> AddEntry(hEffSysM[iRes], sLabelSys[0].Data(), "pf");
    legS[iRes] -> AddEntry(hEffSysP[iRes], sLabelSys[1].Data(), "pf");
    legS[iRes] -> AddEntry(hEffAvg[iRes], sLabelSys[2].Data(), "pf");
  }
  cout << "    Made legends." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXyD[NVtx] = {0.3, 0.1, 0.5, 0.25};
  const Float_t fTxtXyS[NVtx] = {0.3, 0.1, 0.5, 0.3};

  TPaveText *txtD = new TPaveText(fTxtXyD[0], fTxtXyD[1], fTxtXyD[2], fTxtXyD[3], "NDC NB");
  txtD -> SetFillColor(fColTx);
  txtD -> SetFillStyle(fFilTx);
  txtD -> SetLineColor(fColTx);
  txtD -> SetLineStyle(fLinTx);
  txtD -> SetTextFont(fTxt);
  txtD -> SetTextAlign(fAln);
  txtD -> AddText(sSys.Data());
  txtD -> AddText(sTrg.Data());
  txtD -> AddText(sJetD.Data());

  TPaveText *txtS[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    txtS[iRes] = new TPaveText(fTxtXyS[0], fTxtXyS[1], fTxtXyS[2], fTxtXyS[3], "NDC NB");
    txtS[iRes] -> SetFillColor(fColTx);
    txtS[iRes] -> SetFillStyle(fFilTx);
    txtS[iRes] -> SetLineColor(fColTx);
    txtS[iRes] -> SetLineStyle(fLinTx);
    txtS[iRes] -> SetTextFont(fTxt);
    txtS[iRes] -> SetTextAlign(fAln);
    txtS[iRes] -> AddText(sSys.Data());
    txtS[iRes] -> AddText(sTrg.Data());
    txtS[iRes] -> AddText(sJetS[iRes].Data());
    txtS[iRes] -> AddText(sTyp.Data());
  }
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

  // make plots
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  bigWidth(1500);
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
  const Float_t fPadXY1[NVtx] = {0., 0., 0.5, 1.};
  const Float_t fPadXY2[NVtx] = {0.5, 0., 1., 1.};

  // jet efficiencies
  TCanvas *cEff = new TCanvas(sCanDef.Data(), "", width, height);
  cEff       -> SetGrid(fGrid, fGrid);
  cEff       -> SetTicks(fTick, fTick);
  cEff       -> SetLogx(fLogX);
  cEff       -> SetLogy(fLogY);
  cEff       -> SetBorderMode(fMode);
  cEff       -> SetBorderSize(fBord);
  cEff       -> SetFrameBorderMode(fFrame);
  cEff       -> SetLeftMargin(fMarginL);
  cEff       -> SetRightMargin(fMarginR);
  cEff       -> SetTopMargin(fMarginT);
  cEff       -> SetBottomMargin(fMarginB);
  cEff       -> cd();
  hEffDef[0] -> Draw("E5");
  hOutDef[0] -> Draw("E5 same");
  hLinDef[0] -> Draw("hist ][ L same");
  for (UInt_t iRes = 1; iRes < NRes; iRes++) {
    hEffDef[iRes] -> Draw("E5 same");
    hOutDef[iRes] -> Draw("E5 same");
    hLinDef[iRes] -> Draw("hist ][ L same");
  }
  line    -> Draw();
  legD    -> Draw();
  txtD    -> Draw();
  fOutput -> cd();
  cEff    -> Write();
  cEff    -> Close();

  TPad    *pEff[NRes];
  TPad    *pSys[NRes];
  TCanvas *cSys[NRes];
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    cSys[iRes] = new TCanvas(sCanSys[iRes].Data(), "", bigWidth, height);
    pEff[iRes] = new TPad(sPadSys[0].Data(), "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
    pSys[iRes] = new TPad(sPadSys[1].Data(), "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
    pEff[iRes]     -> SetGrid(fGrid, fGrid);
    pEff[iRes]     -> SetTicks(fTick, fTick);
    pEff[iRes]     -> SetLogx(fLogX);
    pEff[iRes]     -> SetLogy(fLogY);
    pEff[iRes]     -> SetBorderMode(fMode);
    pEff[iRes]     -> SetBorderSize(fBord);
    pEff[iRes]     -> SetFrameBorderMode(fFrame);
    pEff[iRes]     -> SetLeftMargin(fMarginL);
    pEff[iRes]     -> SetRightMargin(fMarginR);
    pEff[iRes]     -> SetTopMargin(fMarginT);
    pEff[iRes]     -> SetBottomMargin(fMarginB);
    pSys[iRes]     -> SetGrid(fGrid, fGrid);
    pSys[iRes]     -> SetTicks(fTick, fTick);
    pSys[iRes]     -> SetLogx(fLogX);
    pSys[iRes]     -> SetLogy(fLogY);
    pSys[iRes]     -> SetBorderMode(fMode);
    pSys[iRes]     -> SetBorderSize(fBord);
    pSys[iRes]     -> SetFrameBorderMode(fFrame);
    pSys[iRes]     -> SetLeftMargin(fMarginL);
    pSys[iRes]     -> SetRightMargin(fMarginR);
    pSys[iRes]     -> SetTopMargin(fMarginT);
    pSys[iRes]     -> SetBottomMargin(fMarginB);
    cSys[iRes]     -> cd();
    pEff[iRes]     -> Draw();
    pSys[iRes]     -> Draw();
    pEff[iRes]     -> cd();
    hEffSysM[iRes] -> Draw();
    hEffSysP[iRes] -> Draw("same");
    hEffAvg[iRes]  -> Draw("same");
    line           -> Draw();
    legS[iRes]     -> Draw();
    pSys[iRes]     -> cd();
    hOutErr[iRes]  -> Draw("E2");
    hEffErr[iRes]  -> Draw("E2 same");
    line           -> Draw();
    txtS[iRes]     -> Draw();
    fOutput        -> cd();
    cSys[iRes]     -> Write();
    cSys[iRes]     -> Close();
  }
  cout << "    Made plot." << endl;

  // save histograms
  fOutput -> cd();
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    hEffDef[iRes]  -> Write();
    hEffSysM[iRes] -> Write();
    hEffSysP[iRes] -> Write();
    hEffAvg[iRes]  -> Write();
    hEffErr[iRes]  -> Write();
    hLinDef[iRes]  -> Write();
    hOutDef[iRes]  -> Write();
    hLinErr[iRes]  -> Write();
    hOutErr[iRes]  -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iRes = 0; iRes < NRes; iRes++) {
    fInDef[iRes]  -> cd();
    fInDef[iRes]  -> Close();
    fInSysM[iRes] -> cd();
    fInSysP[iRes] -> Close();
  }
  cout << "  Finished jet efficiency plot!\n" << endl;

}

// End ------------------------------------------------------------------------
