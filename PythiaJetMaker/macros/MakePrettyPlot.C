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



void MakePrettyPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning plot macro..." << endl;

  // output and denominator parameters
  const TString sOutput("pp200r9vsPy8and6.et911r05pi0.d5m3y2022.root");
  const TString sDenom("../DataJetMaker/pp200r9data.minBiasPi0check_pTbinOne.et911r05.d5m3y2022.root");
  const TString sHistDenom("Pi0/hJetPtCorrP");
  const TString sNameDenom("hPtJetData");
  const TString sLabelDenom("Raw data");

  // numerator parameters
  const TString sNumer[NNumer]      = {"./pp200py8par.minBiasPi0check_pTbinOne.et911r05pi0.d4m3y2022.root", "../StarPythiaJetMaker/pp200py6par.minBiasPi0check_pTbinOne.et911r05pi0.d4m3y2022.root"};
  const TString sHistNumer[NNumer]  = {"Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP"};
  const TString sNameNumer[NNumer]  = {"hPtJetPy8", "hPtJetPy6"};
  const TString sNameRatio[NNumer]  = {"hRatioPy8", "hRatioPy6"};
  const TString sLabelNumer[NNumer] = {"Pythia8 [#hat{p}_{T} > 4 GeV/c]", "Pythia6 [#hat{p}_{T} > 4 GeV/c]"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - (#rho #upoint A_{jet}) [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{2}N^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("pythia / data");
  const Float_t xPlotRange[NPlot] = {-1., 30.};
  const UInt_t  fColDen(923);
  const UInt_t  fMarDen(29);
  const UInt_t  fColNum[NNumer] = {859, 899};
  const UInt_t  fMarNum[NNumer] = {24, 25};

  // text parameters
  const TString sSys("data vs. pythia8/6, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
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

  // QUICK FIX [02.28.2022]
  const Double_t loNormBin(9.5);
  const Double_t hiNormBin(29.5);
  const UInt_t   iLoNormD = hDenom    -> FindBin(loNormBin);
  const UInt_t   iLoNormN = hNumer[0] -> FindBin(loNormBin);
  const UInt_t   iHiNormD = hDenom    -> FindBin(hiNormBin);
  const UInt_t   iHiNormN = hNumer[0] -> FindBin(hiNormBin);

  const Double_t iNormD = hDenom    -> Integral(iLoNormD, iHiNormD);
  const Double_t iNormN = hNumer[0] -> Integral(iLoNormN, iHiNormN);
  const Double_t norm   = iNormN / iNormD;
  hDenom -> Scale(norm);


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
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hDenom -> SetMarkerColor(fColDen);
  hDenom -> SetMarkerStyle(fMarDen);
  hDenom -> SetFillColor(fColDen);
  hDenom -> SetFillStyle(fFil);
  hDenom -> SetLineColor(fColDen);
  hDenom -> SetLineStyle(fLin);
  hDenom -> SetLineWidth(fWid);
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
    hNumer[iNumer] -> SetMarkerColor(fColNum[iNumer]);
    hNumer[iNumer] -> SetMarkerStyle(fMarNum[iNumer]);
    hNumer[iNumer] -> SetFillColor(fColNum[iNumer]);
    hNumer[iNumer] -> SetFillStyle(fFil);
    hNumer[iNumer] -> SetLineColor(fColNum[iNumer]);
    hNumer[iNumer] -> SetLineStyle(fLin);
    hNumer[iNumer] -> SetLineWidth(fWid);
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
    hRatio[iNumer]  -> SetMarkerColor(fColNum[iNumer]);
    hRatio[iNumer]  -> SetMarkerStyle(fMarNum[iNumer]);
    hRatio[iNumer]  -> SetFillColor(fColNum[iNumer]);
    hRatio[iNumer]  -> SetFillStyle(fFil);
    hRatio[iNumer]  -> SetLineColor(fColNum[iNumer]);
    hRatio[iNumer]  -> SetLineStyle(fLin);
    hRatio[iNumer]  -> SetLineWidth(fWid);
    hRatio[iNumer]  -> SetTitle(sTitle.Data());
    hRatio[iNumer]  -> SetTitleFont(fTxt);
    hRatio[iNumer]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatio[iNumer]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatio[iNumer]  -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio[iNumer]  -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatio[iNumer]  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatio[iNumer]  -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio[iNumer]  -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatio[iNumer]  -> GetXaxis() -> CenterTitle(fCnt);
    hRatio[iNumer]  -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatio[iNumer]  -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio[iNumer]  -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatio[iNumer]  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatio[iNumer]  -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio[iNumer]  -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatio[iNumer]  -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

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
  leg -> AddEntry(hDenom, sLabelDenom.Data(), "pf");
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    leg -> AddEntry(hNumer[iNumer], sLabelNumer[iNumer], "pf");
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
  const UInt_t  fLinLi(2);
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
  hRatio[0] -> Draw();
  for (UInt_t iNumer = 1; iNumer < NNumer; iNumer++) {
    hRatio[iNumer] -> Draw("SAME");
  }
  line   -> Draw();
  pPad2  -> cd();
  hDenom -> Draw();
  for(UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    hNumer[iNumer] -> Draw("SAME");
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
  for (UInt_t iNumer = 0; iNumer < NNumer; iNumer++) {
    hNumer[iNumer] -> Write();
    hRatio[iNumer] -> Write();
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
