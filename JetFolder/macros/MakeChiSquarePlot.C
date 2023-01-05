// 'MakeChiSquarePlot.C'
// Derek Anderson
// 10.20.2021
//
// Use this to plot a set of chi2/ndf
// as a function of regularization.

#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NChi(3);
static const UInt_t NPlot(2);
static const UInt_t NVtx(4);



void MakeChiSquarePlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning chi2 plotting macro..." << endl;

  // io/hist parameters
  const TString sOutput("backfoldChiSquareFromData.et911r02pi0.d27m10y2021.root");
  const TString sInput[NChi]    = {"et911r02pi0_rebinUnfold/pp200r9pi0.forChi2doubleCheck_pTbinHuge.et911r02qt05130.p0n38t5.performance.root",
                                   "et1115r02pi0_rebinUnfold/pp200r9pi0.forChi2doubleCheck_pTbinHuge.et1115r02qt05130.p0n38t5.performance.root",
                                   "et1520r02pi0_rebinUnfold/pp200r9pi0.forChi2doubleCheck_pTbinHuge.et1520r02qt05130.p0n38t5.performance.root"};
  const TString sHistChi[NChi]  = {"hBayBackfold", "hBayBackfold", "hBayBackfold"};
  const TString sNameChi[NChi]  = {"hBackfoldChi911", "hBackfoldChi1115", "hBackfoldChi1520"};
  const TString sLabelChi[NChi] = {"9 - 11 GeV", "11 - 15 GeV", "15 - 20 GeV"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("n_{iter}");
  const TString sTitleY("#chi^{2}/ndf(backfold, data)");
  const Float_t xPlotRange[NPlot] = {1., 7.};
  const UInt_t  fColChi[NChi]     = {602, 434, 418};
  const UInt_t  fMarChi[NChi]     = {20, 21, 22};
  const UInt_t  fWidChi[NChi]     = {1, 1, 1};

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.2");
  const TString sTyp("#bf{charged jets}");
  const TString sTest("");

  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  TFile *fInput[NChi];
  for (UInt_t iChi = 0; iChi < NChi; iChi++) {
    fInput[iChi] = new TFile(sInput[iChi].Data(), "read");
    if (!fInput[iChi]) {
      cerr << "PANIC: couldn't open input file #" << iChi << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hChi[NChi];
  for (UInt_t iChi = 0; iChi < NChi; iChi++) {
    hChi[iChi] = (TH1D*) fInput[iChi] -> Get(sHistChi[iChi]);
    if (!hChi[iChi]) {
      cerr << "PANIC: couldn't grab histogram #" << iChi << "!" << endl;
      return;
    }
    hChi[iChi] -> SetName(sNameChi[iChi].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fSiz(1.15);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  for (UInt_t iChi = 0; iChi < NChi; iChi++) {
    hChi[iChi] -> SetMarkerColor(fColChi[iChi]);
    hChi[iChi] -> SetMarkerStyle(fMarChi[iChi]);
    hChi[iChi] -> SetMarkerSize(fSiz);
    hChi[iChi] -> SetFillColor(fColChi[iChi]);
    hChi[iChi] -> SetFillStyle(fFil);
    hChi[iChi] -> SetLineColor(fColChi[iChi]);
    hChi[iChi] -> SetLineStyle(fLin);
    hChi[iChi] -> SetLineWidth(fWidChi[iChi]);
    hChi[iChi] -> SetTitle(sTitle.Data());
    hChi[iChi] -> SetTitleFont(fTxt);
    hChi[iChi] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hChi[iChi] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hChi[iChi] -> GetXaxis() -> SetTitleFont(fTxt);
    hChi[iChi] -> GetXaxis() -> SetTitleSize(fTit);
    hChi[iChi] -> GetXaxis() -> SetTitleOffset(fOffX);
    hChi[iChi] -> GetXaxis() -> SetLabelFont(fTxt);
    hChi[iChi] -> GetXaxis() -> SetLabelSize(fLab);
    hChi[iChi] -> GetXaxis() -> CenterTitle(fCnt);
    hChi[iChi] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hChi[iChi] -> GetYaxis() -> SetTitleFont(fTxt);
    hChi[iChi] -> GetYaxis() -> SetTitleSize(fTit);
    hChi[iChi] -> GetYaxis() -> SetTitleOffset(fOffY);
    hChi[iChi] -> GetYaxis() -> SetLabelFont(fTxt);
    hChi[iChi] -> GetYaxis() -> SetLabelSize(fLab);
    hChi[iChi] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t hObj(0.05);
  const Float_t hLeg((NChi + 1) * hObj);
  const Float_t yLeg(0.1 + hLeg);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yLeg};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sTest.Data());
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iChi = 0; iChi < NChi; iChi++) {
    leg -> AddEntry(hChi[iChi], sLabelChi[iChi], "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, 0.25};
  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(2);
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
  hChi[0] -> Draw("hist pl ][");
  for(UInt_t iChi = 1; iChi < NChi; iChi++) {
    hChi[iChi] -> Draw("hist pl ][ same");
  }
  line    -> Draw();
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOutput -> cd();
  for (UInt_t iChi = 0; iChi < NChi; iChi++) {
    hChi[iChi] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iChi = 0; iChi < NChi; iChi++) {
    fInput[iChi] -> cd();
    fInput[iChi] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// End ------------------------------------------------------------------------
