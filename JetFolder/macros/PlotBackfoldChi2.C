// 'PlotBackfoldChi2.C'
// Derek Anderson
// 08.02.2021
//
// Use this to plot a set of chi2/ndf's
// between backfolded and measured/
// detector-level distributions.

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
static const UInt_t NRange(2);
static const UInt_t NHist(3);
static const UInt_t NVtx(4);
 


void PlotBackfoldChi2() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning chi2/ndf plot routine..." << endl;

  // io parameters
  const TString sOutput("backfoldChi2RFF.nIterCheck_unfoldedWithFF.et920r02qt05130.d2m8y2021.root");
  const TString sChi2[NHist]      = {"output/August2021/pp200r9rff.kRegCheck_unfoldWithFF_pTbinFine.et911r02qt05130.p0n38t5.performance.root", "output/August2021/pp200r9rff.kRegCheck_unfoldWithFF_pTbinFine.et1115r02qt05130.p0n38t5.performance.root", "output/August2021/pp200r9rff.kRegCheck_unfoldWithFF_pTbinBig.et1520r02qt05130.p0n38t5.performance.root"};
  const TString sHistChi2[NHist]  = {"hBayBackfold", "hBayBackfold", "hBayBackfold"};
  const TString sNameChi2[NHist]  = {"hBackfoldChi2_et911", "hBackfoldChi2_et1115", "hBackfoldChi2_et1520"};
  const TString sLabelChi2[NHist] = {"backfolded Py6#oplusGeant(RFF) #chi^{2}/ndf vs. n_{iter} #color[899]{[E_{T}^{trg} #in (9, 11) GeV]}", "backfolded Py6#oplusGeant(RFF) #chi^{2}/ndf vs. n_{iter} #color[859]{[E_{T}^{trg} #in (11, 15) GeV]}", "backfolded Py6#oplusGeant(RFF) #chi^{2}/ndf vs. n_{iter} #color[819]{[E_{T}^{trg} #in (15, 20) GeV]}"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("n_{iter}");
  const TString sTitleY("#chi^{2}/ndf(backfold, detector level)");
  const Float_t xPlotRange[NRange] = {1., 6.};
  const UInt_t  fColChi2[NHist]    = {899, 859, 819};
  const UInt_t  fMarChi2[NHist]    = {24, 26, 25};
  const UInt_t  fLinChi2[NHist]    = {1, 1, 1};

  // text parameters
  const TString sSys("Py6#oplusGeant(RFF), #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.2");
  const TString sTyp("#bf{charged jets}");
  const TString sTest("Py6#oplusGeant(RFF) unfolded w/ Py6#oplusGeant(FF)");

  // open output file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }
  cout << "    Opened output and target files." << endl;

  // open chi2/ndf files
  TFile *fChi2[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fChi2[iHist] = new TFile(sChi2[iHist].Data(), "read");
    if (!fChi2[iHist]) {
      cerr << "PANIC: couldn't open chi2/ndf file #" << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened unfolding files." << endl;

  // grab chi2/ndf histograms
  TH1D *hChi2[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hChi2[iHist] = (TH1D*) fChi2[iHist] -> Get(sHistChi2[iHist]);
    if (!hChi2[iHist]) {
      cerr << "PANIC: couldn't grab chi2/ndf histogram #" << iHist << "!" << endl;
      return;
    }
    hChi2[iHist] -> SetName(sNameChi2[iHist].Data());
  }
  cout << "    Grabbed unfolding histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.3);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hChi2[iHist] -> SetMarkerColor(fColChi2[iHist]);
    hChi2[iHist] -> SetMarkerStyle(fMarChi2[iHist]);
    hChi2[iHist] -> SetLineColor(fColChi2[iHist]);
    hChi2[iHist] -> SetLineStyle(fLinChi2[iHist]);
    hChi2[iHist] -> SetLineWidth(fWid);
    hChi2[iHist] -> SetFillColor(fColChi2[iHist]);
    hChi2[iHist] -> SetFillStyle(fFil);
    hChi2[iHist] -> SetTitle(sTitle.Data());
    hChi2[iHist] -> SetTitleFont(fTxt);
    hChi2[iHist] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hChi2[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hChi2[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hChi2[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hChi2[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hChi2[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hChi2[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hChi2[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hChi2[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hChi2[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hChi2[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hChi2[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hChi2[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hChi2[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hChi2[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, 0.3};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sTest.Data());
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    leg -> AddEntry(hChi2[iHist], sLabelChi2[iHist]);
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
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(1);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NVtx] = {xPlotRange[0], 1., xPlotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

  // make plot
  const UInt_t  width(950);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.1);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.1);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot    -> SetGrid(fGrid, fGrid);
  cPlot    -> SetTicks(fTick, fTick);
  cPlot    -> SetBorderMode(fMode);
  cPlot    -> SetBorderSize(fBord);
  cPlot    -> SetFrameBorderMode(fFrame);
  cPlot    -> SetLeftMargin(fMarginL);
  cPlot    -> SetRightMargin(fMarginR);
  cPlot    -> SetTopMargin(fMarginT);
  cPlot    -> SetBottomMargin(fMarginB);
  cPlot    -> cd();
  hChi2[0] -> Draw("hist p l");
  for(UInt_t iHist = 1; iHist < NHist; iHist++) {
    hChi2[iHist] -> Draw("hist p l same");
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
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hChi2[iHist] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fChi2[iHist] -> cd();
    fChi2[iHist] -> Close();
  }
  cout << "  Finished chi2/ndf plot!\n" << endl;

}

// End ------------------------------------------------------------------------
