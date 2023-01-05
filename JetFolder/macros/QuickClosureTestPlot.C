// 'QuickClosureTestPlot.C'
// Derek Anderson
// 04.17.2021
//
// Use this quickly plot the results of a
// closure test (i.e. the unfolding results
// vs. their target particle-level
// distribution).


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
static const UInt_t NIter(6);
static const UInt_t NPlot(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);



void QuickClosureTestPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning closure test plot..." << endl;

  // output and target parameters
  const TString sOutput("numIterVsParCheckRFF_pTbinHuge.et911r02pi0.d29m10y2021.root");
  const TString sTarget("et911r02rff_rebinClosure/pp200r9rff.modifiedStats_parLvlRFF_pTbinHuge.et911r02pi0.d4m10y2021.root");
  const TString sHistTarget("hSimJetOutput");
  const TString sNameTarget("hParticle");
  const TString sLabelTarget("par.-lvl. [with mod. #sigma_{stat}]");

  // unfolding parameters
  const TString sUnfold[NIter]      = {"et911r02rff_rebinClosure/pp200r9rff.forNumIterCheck_UnfoldWithFF_pTbinHuge.et911r02qt05130.p0m1k1n38t5.root",
                                       "et911r02rff_rebinClosure/pp200r9rff.forNumIterCheck_UnfoldWithFF_pTbinHuge.et911r02qt05130.p0m1k2n38t5.root",
                                       "et911r02rff_rebinClosure/pp200r9rff.forNumIterCheck_UnfoldWithFF_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root",
                                       "et911r02rff_rebinClosure/pp200r9rff.forNumIterCheck_UnfoldWithFF_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root",
                                       "et911r02rff_rebinClosure/pp200r9rff.forNumIterCheck_UnfoldWithFF_pTbinHuge.et911r02qt05130.p0m1k5n38t5.root",
                                       "et911r02rff_rebinClosure/pp200r9rff.forNumIterCheck_UnfoldWithFF_pTbinHuge.et911r02qt05130.p0m1k6n38t5.root"};
  const TString sHistUnfold[NIter]  = {"hUnfolded", "hUnfolded", "hUnfolded", "hUnfolded", "hUnfolded", "hUnfolded"};
  const TString sNameUnfold[NIter]  = {"hUnfoldN1", "hUnfoldN2", "hUnfoldN3", "hUnfoldN4", "hUnfoldN5", "hUnfoldedN6"};
  const TString sNameRatio[NIter]   = {"hRatioN1", "hRatioN2", "hRatioN3", "hRatioN4", "hRatioN5", "hRatioN6"};
  const TString sLabelUnfold[NIter] = {"unfolded det.-lvl. #bf{#color[799]{[n_{iter} = 1]}}",
                                       "unfolded det.-lvl. #bf{#color[899]{[n_{iter} = 2]}}",
                                       "unfolded det.-lvl. #bf{#color[879]{[n_{iter} = 3]}}",
                                       "unfolded det.-lvl. #bf{#color[859]{[n_{iter} = 4]}}",
                                       "unfolded det.-lvl. #bf{#color[839]{[n_{iter} = 5]}}",
                                       "unfolded det.-lvl. #bf{#color[819]{[n_{iter} = 6]}}"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{par} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{3}N^{jet}/d(p_{T}^{par} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("unfolded / particle");
  const Float_t xPlotRange[NPlot] = {-1., 30.};
  const UInt_t  fColTar(923);
  const UInt_t  fMarTar(8);
  const UInt_t  fColUnf[NIter] = {799, 899, 879, 859, 839, 819};
  const UInt_t  fMarUnf[NIter] = {24, 26, 32, 25, 28, 30};

  // text parameters
  const TString sSys("Py6#oplusGeant(RFF), #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.2");
  const TString sTyp("#bf{charged jets}");
  const TString sTest("Py6#oplusGeant(RFF) unfolded w/ Py6#oplusGeant(FF)");


  // open output and target files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fTarget = new TFile(sTarget.Data(), "read");
  if (!fOutput || !fTarget) {
    cerr << "PANIC: couldn't open output or target file!\n"
         << "       fOutput = " << fOutput << ", fTarget = " << fTarget
         << endl;
    return;
  }
  cout << "    Opened output and target files." << endl;

  // open unfolding files
  TFile *fUnfold[NIter];
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    fUnfold[iIter] = new TFile(sUnfold[iIter].Data(), "read");
    if (!fUnfold[iIter]) {
      cerr << "PANIC: couldn't open unfolding file #" << iIter << "!" << endl;
      return;
    }
  }
  cout << "    Opened unfolding files." << endl;


  // grab target histogram
  TH1D *hTarget = (TH1D*) fTarget -> Get(sHistTarget.Data());
  if (!hTarget) {
    cerr << "PANIC: couldn't grab target histogram!" << endl;
    return;
  }
  hTarget -> SetName(sNameTarget.Data());
  cout << "    Grabbed target histogram." << endl;

  // grab unfolding histograms
  TH1D *hUnfold[NIter];
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    hUnfold[iIter] = (TH1D*) fUnfold[iIter] -> Get(sHistUnfold[iIter]);
    if (!hUnfold[iIter]) {
      cerr << "PANIC: couldn't grab unfolding histogram #" << iIter << "!" << endl;
      return;
    }
    hUnfold[iIter] -> SetName(sNameUnfold[iIter].Data());
  }
  cout << "    Grabbed unfolding histograms." << endl;


  // calculate ratios
  TH1D *hRatio[NIter];
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    hRatio[iIter] = (TH1D*) hTarget -> Clone();
    hRatio[iIter] -> Reset("ICE");
    hRatio[iIter] -> Divide(hUnfold[iIter], hTarget, 1., 1.);
    hRatio[iIter] -> SetName(sNameRatio[iIter]);
  }
  cout << "    Calculated ratios." << endl;


  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fSiz(1.25);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hTarget -> SetMarkerColor(fColTar);
  hTarget -> SetMarkerStyle(fMarTar);
  hTarget -> SetMarkerSize(fSiz);
  hTarget -> SetFillColor(fColTar);
  hTarget -> SetFillStyle(fFil);
  hTarget -> SetLineColor(fColTar);
  hTarget -> SetLineStyle(fLin);
  hTarget -> SetLineWidth(fWid);
  hTarget -> SetTitle(sTitle.Data());
  hTarget -> SetTitleFont(fTxt);
  hTarget -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hTarget -> GetXaxis() -> SetTitle(sTitleX.Data());
  hTarget -> GetXaxis() -> SetTitleFont(fTxt);
  hTarget -> GetXaxis() -> SetTitleSize(fTit[1]);
  hTarget -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hTarget -> GetXaxis() -> SetLabelFont(fTxt);
  hTarget -> GetXaxis() -> SetLabelSize(fLab[1]);
  hTarget -> GetXaxis() -> CenterTitle(fCnt);
  hTarget -> GetYaxis() -> SetTitle(sTitleY.Data());
  hTarget -> GetYaxis() -> SetTitleFont(fTxt);
  hTarget -> GetYaxis() -> SetTitleSize(fTit[1]);
  hTarget -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hTarget -> GetYaxis() -> SetLabelFont(fTxt);
  hTarget -> GetYaxis() -> SetLabelSize(fLab[1]);
  hTarget -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    hUnfold[iIter] -> SetMarkerColor(fColUnf[iIter]);
    hUnfold[iIter] -> SetMarkerStyle(fMarUnf[iIter]);
    hUnfold[iIter] -> SetMarkerSize(fSiz);
    hUnfold[iIter] -> SetFillColor(fColUnf[iIter]);
    hUnfold[iIter] -> SetFillStyle(fFil);
    hUnfold[iIter] -> SetLineColor(fColUnf[iIter]);
    hUnfold[iIter] -> SetLineStyle(fLin);
    hUnfold[iIter] -> SetLineWidth(fWid);
    hUnfold[iIter] -> SetTitle(sTitle.Data());
    hUnfold[iIter] -> SetTitleFont(fTxt);
    hUnfold[iIter] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hUnfold[iIter] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hUnfold[iIter] -> GetXaxis() -> SetTitleFont(fTxt);
    hUnfold[iIter] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hUnfold[iIter] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hUnfold[iIter] -> GetXaxis() -> SetLabelFont(fTxt);
    hUnfold[iIter] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hUnfold[iIter] -> GetXaxis() -> CenterTitle(fCnt);
    hUnfold[iIter] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hUnfold[iIter] -> GetYaxis() -> SetTitleFont(fTxt);
    hUnfold[iIter] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hUnfold[iIter] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hUnfold[iIter] -> GetYaxis() -> SetLabelFont(fTxt);
    hUnfold[iIter] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hUnfold[iIter] -> GetYaxis() -> CenterTitle(fCnt);
    hRatio[iIter]  -> SetMarkerColor(fColUnf[iIter]);
    hRatio[iIter]  -> SetMarkerStyle(fMarUnf[iIter]);
    hRatio[iIter]  -> SetMarkerSize(fSiz);
    hRatio[iIter]  -> SetFillColor(fColUnf[iIter]);
    hRatio[iIter]  -> SetFillStyle(fFil);
    hRatio[iIter]  -> SetLineColor(fColUnf[iIter]);
    hRatio[iIter]  -> SetLineStyle(fLin);
    hRatio[iIter]  -> SetLineWidth(fWid);
    hRatio[iIter]  -> SetTitle(sTitle.Data());
    hRatio[iIter]  -> SetTitleFont(fTxt);
    hRatio[iIter]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatio[iIter]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatio[iIter]  -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio[iIter]  -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatio[iIter]  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatio[iIter]  -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio[iIter]  -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatio[iIter]  -> GetXaxis() -> CenterTitle(fCnt);
    hRatio[iIter]  -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatio[iIter]  -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio[iIter]  -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatio[iIter]  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatio[iIter]  -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio[iIter]  -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatio[iIter]  -> GetYaxis() -> CenterTitle(fCnt);
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
  leg -> AddEntry(hTarget, sLabelTarget.Data(), "pf");
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    leg -> AddEntry(hUnfold[iIter], sLabelUnfold[iIter], "pf");
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
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY1(0);
  const UInt_t  fLogY2(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT1(0.005);
  const Float_t fMarginT2(0.02);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.005);
  const Float_t fPadXY1[NVtx] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NVtx] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot     -> SetGrid(fGrid, fGrid);
  cPlot     -> SetTicks(fTick, fTick);
  cPlot     -> SetBorderMode(fMode);
  cPlot     -> SetBorderSize(fBord);
  pPad1     -> SetGrid(fGrid, fGrid);
  pPad1     -> SetTicks(fTick, fTick);
  pPad1     -> SetLogx(fLogX);
  pPad1     -> SetLogy(fLogY1);
  pPad1     -> SetBorderMode(fMode);
  pPad1     -> SetBorderSize(fBord);
  pPad1     -> SetFrameBorderMode(fFrame);
  pPad1     -> SetLeftMargin(fMarginL);
  pPad1     -> SetRightMargin(fMarginR);
  pPad1     -> SetTopMargin(fMarginT1);
  pPad1     -> SetBottomMargin(fMarginB1);
  pPad2     -> SetGrid(fGrid, fGrid);
  pPad2     -> SetTicks(fTick, fTick);
  pPad2     -> SetLogx(fLogX);
  pPad2     -> SetLogy(fLogY2);
  pPad2     -> SetBorderMode(fMode);
  pPad2     -> SetBorderSize(fBord);
  pPad2     -> SetFrameBorderMode(fFrame);
  pPad2     -> SetLeftMargin(fMarginL);
  pPad2     -> SetRightMargin(fMarginR);
  pPad2     -> SetTopMargin(fMarginT2);
  pPad2     -> SetBottomMargin(fMarginB2);
  cPlot     -> cd();
  pPad1     -> Draw();
  pPad2     -> Draw();
  pPad1     -> cd();
  hRatio[0] -> Draw("E5");
  for (UInt_t iIter = 1; iIter < NIter; iIter++) {
    hRatio[iIter] -> Draw("E5 SAME");
  }
  line    -> Draw();
  pPad2   -> cd();
  hTarget -> Draw();
  for(UInt_t iIter = 0; iIter < NIter; iIter++) {
    hUnfold[iIter] -> Draw("E5 SAME");
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;


  // save histograms
  fOutput -> cd();
  hTarget -> Write();
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    hUnfold[iIter] -> Write();
    hRatio[iIter]  -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fTarget -> cd();
  fTarget -> Close();
  for (UInt_t iIter = 0; iIter < NIter; iIter++) {
    fUnfold[iIter] -> cd();
    fUnfold[iIter] -> Close();
  }
  cout << "  Finished closure test plot!\n" << endl;

}

// End ------------------------------------------------------------------------
