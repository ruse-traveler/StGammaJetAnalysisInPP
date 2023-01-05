// 'MakePrettyPlot.C'
// Derek Anderson
// 04.18.2021
//
// Use this quickly plot a set of numerator
// distributions against a denominator
// distribution

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
static const UInt_t NNumer(2);
static const UInt_t NPlot(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const Bool_t DoOutline(false);



void MakePrettyPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning plot macro..." << endl;

  // output and denominator parameters
  const TString sOutput("weightVsUnweightPythia_pTbinHuge.et911r02gam.d7m11y2021.root");
  const TString sDenom("input/pythia/pp200py8par.forComparison_pTbinHuge.et911nTrg100Kpt021Kgam.r02a005rm1chrg.root");
  const TString sHistDenom("Gam/hJetPtCorrG");
  const TString sNameDenom("hUnweighted");
  const TString sLabelDenom("Unweighted Pythia8");

  // numerator parameters
  const TString sNumer[NNumer]      = {"input/weightedPythia/pp200py8par.forComparison_fitWeighted911_pTbinHuge.et0100r02a005gam.root", "input/weightedPythia/pp200py8par.forComparison_histWeighted911_pTbinHuge.et0100r02a005gam.root"};
  const TString sHistNumer[NNumer]  = {"Gam/hJetPtCorrG", "Gam/hJetPtCorrG"};
  const TString sNameNumer[NNumer]  = {"hFitWeighted", "hHistWeighted"};
  const TString sNameRatio[NNumer]  = {"hRatioFit", "hRatioHist"};
  const TString sLabelNumer[NNumer] = {"Weighted with #Delta^{Py8}_{F}", "Weighted with #Delta^{Py8}_{D}"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{3}N^{jet}/d(p_{T}^{jet} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("weighted / unweighted");
  const Float_t xPlotRange[NPlot] = {-1., 30.};
  const UInt_t  fColMarD(923);
  const UInt_t  fColFilD(923);
  const UInt_t  fMarDen(29);
  const UInt_t  fSizDen(1.);
  const UInt_t  fFilDen(0);
  const UInt_t  fLinDen(1);
  const UInt_t  fWidDen(1);
  const UInt_t  fColMarN[NNumer] = {899, 896};
  const UInt_t  fColFilN[NNumer] = {899, 896};
  const UInt_t  fMarNum[NNumer]  = {20, 24};
  const UInt_t  fSizNum[NNumer]  = {1., 1.};
  const UInt_t  fFilNum[NNumer]  = {0, 0};
  const UInt_t  fLinNum[NNumer]  = {1, 1};
  const UInt_t  fWidNum[NNumer]  = {1, 1};

  // text parameters
  const TString sSys("Pythia8, #sqrt{s} = 200 GeV");
  const TString sTrg("#gamma_{dir} trig., E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.2");
  const TString sTyp("#bf{charged jets}");

  // open output and target files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fDenom  = new TFile(sDenom.Data(), "read");
  if (!fOutput || !fDenom) {
    cerr << "PANIC: couldn't open output or denominator file!\n"
         << "       fOutput = " << fOutput << ", fDenom = " << fDenom
         << endl;
    return;
  }
  cout << "    Opened output and denominator files." << endl;

  // open unfolding files
  TFile *fNumer[NNumer];
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    fNumer[iNumer] = new TFile(sNumer[iNumer].Data(), "read");
    if (!fNumer[iNumer]) {
      cerr << "PANIC: couldn't open numerator file #" << iNumer << "!" << endl;
      return;
    }
  }
  cout << "    Opened numerator files." << endl;

  // grab target histogram
  TH1D *hDenom = (TH1D*) fDenom -> Get(sHistDenom.Data());
  if (!hDenom) {
    cerr << "PANIC: couldn't grab denominator histogram!" << endl;
    return;
  }
  hDenom -> SetName(sNameDenom.Data());
  cout << "    Grabbed denominator histogram." << endl;

  // grab unfolding histograms
  TH1D *hNumer[NNumer];
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    hNumer[iNumer] = (TH1D*) fNumer[iNumer] -> Get(sHistNumer[iNumer]);
    if (!hNumer[iNumer]) {
      cerr << "PANIC: couldn't grab numerator histogram #" << iNumer << "!" << endl;
      return;
    }
    hNumer[iNumer] -> SetName(sNameNumer[iNumer].Data());
  }
  cout << "    Grabbed numerator histograms." << endl;

  // calculate ratios
  TH1D *hRatio[NNumer];
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    hRatio[iNumer] = (TH1D*) hDenom -> Clone();
    hRatio[iNumer] -> Reset("ICE");
    hRatio[iNumer] -> Divide(hNumer[iNumer], hDenom, 1., 1.);
    hRatio[iNumer] -> SetName(sNameRatio[iNumer]);
  }
  cout << "    Calculated ratios." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hDenom -> SetMarkerColor(fColMarD);
  hDenom -> SetMarkerStyle(fMarDen);
  hDenom -> SetMarkerSize(fSizDen);
  hDenom -> SetFillColor(fColFilD);
  hDenom -> SetFillStyle(fFilDen);
  hDenom -> SetLineColor(fColFilD);
  hDenom -> SetLineStyle(fLinDen);
  hDenom -> SetLineWidth(fWidDen);
  hDenom -> SetTitle(sTitle.Data());
  hDenom -> SetTitleFont(fTxt);
  hDenom -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hDenom -> GetXaxis() -> SetTitle(sTitleX.Data());
  hDenom -> GetXaxis() -> SetTitleFont(fTxt);
  hDenom -> GetXaxis() -> SetTitleSize(fTit[1]);
  hDenom -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hDenom -> GetXaxis() -> SetLabelFont(fTxt);
  hDenom -> GetXaxis() -> SetLabelSize(fLab[1]);
  hDenom -> GetXaxis() -> CenterTitle(fCnt);
  hDenom -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDenom -> GetYaxis() -> SetTitleFont(fTxt);
  hDenom -> GetYaxis() -> SetTitleSize(fTit[1]);
  hDenom -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hDenom -> GetYaxis() -> SetLabelFont(fTxt);
  hDenom -> GetYaxis() -> SetLabelSize(fLab[1]);
  hDenom -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    hNumer[iNumer] -> SetMarkerColor(fColMarN[iNumer]);
    hNumer[iNumer] -> SetMarkerStyle(fMarNum[iNumer]);
    hNumer[iNumer] -> SetMarkerSize(fSizNum[iNumer]);
    hNumer[iNumer] -> SetFillColor(fColFilN[iNumer]);
    hNumer[iNumer] -> SetFillStyle(fFilNum[iNumer]);
    hNumer[iNumer] -> SetLineColor(fColFilN[iNumer]);
    hNumer[iNumer] -> SetLineStyle(fLinNum[iNumer]);
    hNumer[iNumer] -> SetLineWidth(fWidNum[iNumer]);
    hNumer[iNumer] -> SetTitle(sTitle.Data());
    hNumer[iNumer] -> SetTitleFont(fTxt);
    hNumer[iNumer] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hNumer[iNumer] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hNumer[iNumer] -> GetXaxis() -> SetTitleFont(fTxt);
    hNumer[iNumer] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hNumer[iNumer] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hNumer[iNumer] -> GetXaxis() -> SetLabelFont(fTxt);
    hNumer[iNumer] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hNumer[iNumer] -> GetXaxis() -> CenterTitle(fCnt);
    hNumer[iNumer] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hNumer[iNumer] -> GetYaxis() -> SetTitleFont(fTxt);
    hNumer[iNumer] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hNumer[iNumer] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hNumer[iNumer] -> GetYaxis() -> SetLabelFont(fTxt);
    hNumer[iNumer] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hNumer[iNumer] -> GetYaxis() -> CenterTitle(fCnt);
    hRatio[iNumer] -> SetMarkerColor(fColMarN[iNumer]);
    hRatio[iNumer] -> SetMarkerStyle(fMarNum[iNumer]);
    hRatio[iNumer] -> SetMarkerSize(fSizNum[iNumer]);
    hRatio[iNumer] -> SetFillColor(fColFilN[iNumer]);
    hRatio[iNumer] -> SetFillStyle(fFilNum[iNumer]);
    hRatio[iNumer] -> SetLineColor(fColFilN[iNumer]);
    hRatio[iNumer] -> SetLineStyle(fLinNum[iNumer]);
    hRatio[iNumer] -> SetLineWidth(fWidNum[iNumer]);
    hRatio[iNumer] -> SetTitle(sTitle.Data());
    hRatio[iNumer] -> SetTitleFont(fTxt);
    hRatio[iNumer] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatio[iNumer] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatio[iNumer] -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio[iNumer] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatio[iNumer] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatio[iNumer] -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio[iNumer] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatio[iNumer] -> GetXaxis() -> CenterTitle(fCnt);
    hRatio[iNumer] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatio[iNumer] -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio[iNumer] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatio[iNumer] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatio[iNumer] -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio[iNumer] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatio[iNumer] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // create outlines
  const UInt_t  fFilOut(0);
  const UInt_t  fWidOut(3);
  const TString sOutD("hOutlineDenom");
  const TString sOutBaseN("hOutlineNum");
  const TString sOutBaseR("hOutlineRatio");

  TH1D *hOutD;
  TH1D *hOutN[NNumer];
  TH1D *hOutR[NNumer];
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {

    // create name
    TString sOutN(sOutBaseN.Data());
    TString sOutR(sOutBaseR.Data());
    sOutN += iNumer;
    sOutR += iNumer;

    // create outline
    hOutN[iNumer] = (TH1D*) hNumer[iNumer] -> Clone();
    hOutR[iNumer] = (TH1D*) hRatio[iNumer] -> Clone();
    hOutN[iNumer] -> SetName(sOutN.Data());
    hOutR[iNumer] -> SetName(sOutR.Data());
    hOutN[iNumer] -> SetFillStyle(fFilOut);
    hOutR[iNumer] -> SetFillStyle(fFilOut);
    hOutN[iNumer] -> SetFillColor(fColMarN[iNumer]);
    hOutR[iNumer] -> SetFillColor(fColMarN[iNumer]);
    hOutN[iNumer] -> SetLineWidth(fWidOut);
    hOutR[iNumer] -> SetLineWidth(fWidOut);
  }
  hOutD = (TH1D*) hDenom -> Clone();
  hOutD -> SetName(sOutD.Data());
  hOutD -> SetFillStyle(fFilOut);
  hOutD -> SetFillColor(fColMarD);
  hOutD -> SetLineWidth(fWidOut);
  cout << "    Created outlines." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, 0.3};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  if (DoOutline) {
    leg -> AddEntry(hDenom, sLabelDenom.Data(), "f");
    for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
      leg -> AddEntry(hNumer[iNumer], sLabelNumer[iNumer], "f");
    }
  } else {
    leg -> AddEntry(hDenom, sLabelDenom.Data(), "pf");
    for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
      leg -> AddEntry(hNumer[iNumer], sLabelNumer[iNumer], "pf");
    }
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
  const UInt_t  fColLi(923);
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
  hRatio[0] -> Draw("hist p");
  if (DoOutline) hOutR[0] -> Draw("hist p same");
  for (UInt_t iNumer = 1; iNumer < NNumer; iNumer++) {
    hRatio[iNumer] -> Draw("hist p same");
    if (DoOutline) hOutR[iNumer] -> Draw("hist p same");
  }
  line      -> Draw();
  pPad2     -> cd();
  hDenom    -> Draw();
  if (DoOutline) hOutD -> Draw("hist p same");
  hNumer[0] -> Draw("hist p same");
  if (DoOutline) hOutN[0] -> Draw("E5 same");
  for (UInt_t iNumer = 1; iNumer < NNumer; iNumer++) {
    hNumer[iNumer] -> Draw("hist p same");
    if (DoOutline) hOutN[iNumer] -> Draw("hist p same");
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOutput -> cd();
  hDenom  -> Write();
  if (DoOutline) hOutD -> Write();
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    hNumer[iNumer] -> Write();
    hRatio[iNumer] -> Write();
    if (DoOutline) {
      hOutN[iNumer] -> Write();
      hOutR[iNumer] -> Write();
    }
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fDenom  -> cd();
  fDenom  -> Close();
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    fNumer[iNumer] -> cd();
    fNumer[iNumer] -> Close();
  }
  cout << "  Finished plot!\n" << endl;

}

// End ------------------------------------------------------------------------
