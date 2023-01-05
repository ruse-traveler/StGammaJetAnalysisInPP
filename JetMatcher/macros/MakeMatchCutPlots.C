// 'MakeMatchCutPlots.C'
// Derek Anderson
// 10.08.2021
//
// Use this to make two plots showing
// the qTjet and dRjet of the match vs.
// all match candidates

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TError.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

// global constants
static const UInt_t NHist(2);
static const UInt_t NVtx(4);
static const UInt_t NQt(2);


void MakeMatchCutPlots() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning match cut plot making script..." << endl;

  // io parameters
  const TString sOut("unfoldMatchCuts.et920r02qt05130.d9m10y2021.root");
  const TString sInDr[NHist]   = {"dRall.et920r02qt05130.d9m10y2021.root", "dRsel.et920r02qt05130.d9m10y2021.root"};
  const TString sInQt[NHist]   = {"qTall.et920r02qt05130.d9m10y2021.root", "qTsel.et920r02qt05130.d9m10y2021.root"};
  const TString sHistDr[NHist] = {"hNorm", "hNorm"};
  const TString sHistQt[NHist] = {"hNorm", "hNorm"};

  // match parameters
  const Float_t rJet(0.2);
  const Float_t qTmatch[NQt] = {0.5, 1.3};

  // text parameters
  const TString sSys("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.2");
  const TString sTyp("#bf{charged jets}");

  // histogram parameters
  const TString sNameDr[NHist]  = {"hJetDrAll", "hJetDrMatch"};
  const TString sNameQt[NHist]  = {"hJetQtAll", "hJetQtMatch"};
  const TString sLabelDr[NHist] = {"all accepted jets", "match candidates"};
  const TString sLabelQt[NHist] = {"all accepted jets", "match candidates"};

  // open histograms
  TFile *fInDr[NHist];
  TFile *fInQt[NHist];
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInDr[iHist] = new TFile(sInDr[iHist].Data(), "read");
    fInQt[iHist] = new TFile(sInQt[iHist].Data(), "read");
    if (!fInDr[iHist] || !fInQt[iHist]) {
      cerr << "PANIC: couldn't open an input file!\n"
           << "       fInDr[" << iHist << "] = " << fInDr[iHist] << ", fInQt[" << iHist << "] = " << fInQt[iHist]
           << endl;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hJetDr[NHist];
  TH1D *hJetQt[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hJetDr[iHist] = (TH1D*) fInDr[iHist] -> Get(sHistDr[iHist].Data());
    hJetQt[iHist] = (TH1D*) fInQt[iHist] -> Get(sHistQt[iHist].Data());
    if (!hJetDr[iHist] || !hJetQt[iHist]) {
      cerr << "PANIC: couldn't open an input histogram!\n"
           << "       hJetDr[" << iHist << "] = " << hJetDr[iHist] << ", hJetQt[" << iHist << "] = " << hJetQt[iHist]
           << endl;
    }
    hJetDr[iHist] -> SetName(sNameDr[iHist].Data());
    hJetQt[iHist] -> SetName(sNameQt[iHist].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // normalize histograms
  Double_t intDr[NHist];
  Double_t intQt[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    intDr[iHist] = hJetDr[iHist] -> Integral();
    intQt[iHist] = hJetQt[iHist] -> Integral();
  }
  const Double_t normDrCand  = 1. / intDr[0];
  const Double_t normQtCand  = 1. / intQt[0];
  hJetDr[0] -> Scale(normDrCand);
  hJetQt[0] -> Scale(normQtCand);

  UInt_t iQtMatch[NQt];
  for (UInt_t iQt = 0; iQt < NQt; iQt++) {
    iQtMatch[iQt] = hJetQt[iQt] -> FindBin(qTmatch[iQt]);
  }
  const UInt_t   iDrMatch    = hJetDr[0] -> FindBin(rJet);
  const Double_t intDrMatchC = hJetDr[0] -> Integral(1, iDrMatch - 1);
  const Double_t intQtMatchC = hJetQt[0] -> Integral(iQtMatch[0], iQtMatch[1]);
  const Double_t normDrMatch = intDrMatchC / intDr[1];
  const Double_t normQtMatch = intQtMatchC / intQt[1];
  hJetDr[1] -> Scale(normDrMatch);
  hJetQt[1] -> Scale(normQtMatch);
  cout << "    Normalized histograms:\n"
       << "      normDr(cand., match) = (" << normDrCand << ", " << normDrMatch << ")\n"
       << "      normQt(cand., match) = (" << normQtCand << ", " << normQtMatch << ")"
       << endl;

  // set styles
  const UInt_t  fCol[NHist] = {923, 899};
  const UInt_t  fLin[NHist] = {1, 1};
  const UInt_t  fFil[NHist] = {0, 3345};
  const UInt_t  fMar[NHist] = {8, 8};
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const TString sTitle("");
  const TString sJetDr("#Deltar = #sqrt{#Delta#eta^{2} + #Delta#varphi^{2}}");
  const TString sJetQt("q_{T}^{jet} = p_{T}^{det} / p_{T}^{par}");
  const TString sJetY("arbitrary units");
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.3);
  const Float_t fOffL(0.005);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hJetDr[iHist] -> SetFillColor(fCol[iHist]);
    hJetDr[iHist] -> SetFillStyle(fFil[iHist]);
    hJetDr[iHist] -> SetLineColor(fCol[iHist]);
    hJetDr[iHist] -> SetLineStyle(fLin[iHist]);
    hJetDr[iHist] -> SetMarkerColor(fCol[iHist]);
    hJetDr[iHist] -> SetMarkerStyle(fMar[iHist]);
    hJetDr[iHist] -> SetTitle(sTitle.Data());
    hJetDr[iHist] -> SetTitleFont(fTxt);
    hJetDr[iHist] -> GetXaxis() -> SetTitle(sJetDr.Data());
    hJetDr[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetDr[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hJetDr[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetDr[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetDr[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hJetDr[iHist] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetDr[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hJetDr[iHist] -> GetYaxis() -> SetTitle(sJetY.Data());
    hJetDr[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetDr[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hJetDr[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetDr[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetDr[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hJetDr[iHist] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetDr[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hJetQt[iHist] -> SetFillColor(fCol[iHist]);
    hJetQt[iHist] -> SetFillStyle(fFil[iHist]);
    hJetQt[iHist] -> SetLineColor(fCol[iHist]);
    hJetQt[iHist] -> SetLineStyle(fLin[iHist]);
    hJetQt[iHist] -> SetMarkerColor(fCol[iHist]);
    hJetQt[iHist] -> SetMarkerStyle(fMar[iHist]);
    hJetQt[iHist] -> SetTitle(sTitle.Data());
    hJetQt[iHist] -> SetTitleFont(fTxt);
    hJetQt[iHist] -> GetXaxis() -> SetTitle(sJetQt.Data());
    hJetQt[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hJetQt[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hJetQt[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hJetQt[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hJetQt[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hJetQt[iHist] -> GetXaxis() -> SetLabelOffset(fOffL);
    hJetQt[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hJetQt[iHist] -> GetYaxis() -> SetTitle(sJetY.Data());
    hJetQt[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hJetQt[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hJetQt[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hJetQt[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hJetQt[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hJetQt[iHist] -> GetYaxis() -> SetLabelOffset(fOffL);
    hJetQt[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // create legends
  const UInt_t  fColLeg(0);
  const UInt_t  fFilLeg(0);
  const UInt_t  fLinLeg(0);
  const UInt_t  fAlnLeg(12);
  const UInt_t  fAlnTxt(32);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.2};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  TLegend *lJetDr = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  TLegend *lJetQt = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  lJetDr -> SetFillColor(fColLeg);
  lJetDr -> SetFillStyle(fFilLeg);
  lJetDr -> SetLineColor(fColLeg);
  lJetDr -> SetLineStyle(fLinLeg);
  lJetDr -> SetTextFont(fTxt);
  lJetDr -> SetTextAlign(fAlnLeg);
  lJetQt -> SetFillColor(fColLeg);
  lJetQt -> SetFillStyle(fFilLeg);
  lJetQt -> SetLineColor(fColLeg);
  lJetQt -> SetLineStyle(fLinLeg);
  lJetQt -> SetTextFont(fTxt);
  lJetQt -> SetTextAlign(fAlnLeg);
  lJetDr -> AddEntry(hJetDr[0], sLabelDr[0].Data(), "pf");
  lJetDr -> AddEntry(hJetDr[1], sLabelDr[1].Data(), "f");
  lJetQt -> AddEntry(hJetQt[0], sLabelQt[0].Data(), "pf");
  lJetQt -> AddEntry(hJetQt[1], sLabelQt[1].Data(), "f");

  TPaveText *pt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  pt -> SetFillColor(fColLeg);
  pt -> SetFillStyle(fFilLeg);
  pt -> SetLineColor(fColLeg);
  pt -> SetLineStyle(fLinLeg);
  pt -> SetTextFont(fTxt);
  pt -> SetTextAlign(fAlnTxt);
  pt -> AddText(sSys.Data());
  pt -> AddText(sTrg.Data());
  pt -> AddText(sJet.Data());
  pt -> AddText(sTyp.Data());
  cout << "    Made text boxes." << endl;

  // draw plots
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fFrame(0);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fTick(1);
  const Float_t fMarginBig(0.15);
  const Float_t fMarginSmall(0.02);

  TCanvas *cJetDr = new TCanvas("cJetDr", "", width, height);
  cJetDr    -> SetLogx(fLogX);
  cJetDr    -> SetLogy(fLogY);
  cJetDr    -> SetGrid(fGrid, fGrid);
  cJetDr    -> SetTicks(fTick, fTick);
  cJetDr    -> SetBorderMode(fMode);
  cJetDr    -> SetBorderSize(fBord);
  cJetDr    -> SetFrameBorderMode(fFrame);
  cJetDr    -> SetLeftMargin(fMarginBig);
  cJetDr    -> SetTopMargin(fMarginSmall);
  cJetDr    -> SetRightMargin(fMarginSmall);
  cJetDr    -> SetBottomMargin(fMarginBig);
  cJetDr    -> cd();
  hJetDr[0] -> Draw();
  hJetDr[1] -> Draw("hist same");
  hJetDr[1] -> Draw("same");
  lJetDr    -> Draw();
  pt        -> Draw();
  fOut      -> cd();
  cJetDr    -> Write();
  cJetDr    -> Close();

  TCanvas *cJetQt = new TCanvas("cJetQt", "", width, height);
  cJetQt    -> SetLogx(fLogX);
  cJetQt    -> SetLogy(fLogY);
  cJetQt    -> SetGrid(fGrid, fGrid);
  cJetQt    -> SetTicks(fTick, fTick);
  cJetQt    -> SetBorderMode(fMode);
  cJetQt    -> SetBorderSize(fBord);
  cJetQt    -> SetFrameBorderMode(fFrame);
  cJetQt    -> SetLeftMargin(fMarginBig);
  cJetQt    -> SetTopMargin(fMarginSmall);
  cJetQt    -> SetRightMargin(fMarginSmall);
  cJetQt    -> SetBottomMargin(fMarginBig);
  cJetQt    -> cd();
  hJetQt[0] -> Draw();
  hJetQt[1] -> Draw("hist same");
  hJetQt[1] -> Draw("same");
  lJetQt    -> Draw();
  pt        -> Draw();
  fOut      -> cd();
  cJetQt    -> Write();
  cJetQt    -> Close();
  cout << "    Made plots." << endl;

  // save histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hJetDr[iHist] -> Write();
    hJetQt[iHist] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInDr[iHist] -> cd();
    fInDr[iHist] -> Close();
    fInQt[iHist] -> cd();
    fInQt[iHist] -> Close();
  }
  cout << "  Finished match cut plot making script!\n" << endl;

}

// End ------------------------------------------------------------------------
