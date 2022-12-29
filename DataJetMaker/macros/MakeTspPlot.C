// 'MakeTspPlot.C'
// Derek Anderson
// 09.25.2021
//
// Plots two simulated TSP
// distributions on top of
// a measured TSP

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;



void MakeTspPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning TSP plot script..." << endl;

  // input/output files
  const TString sOut("measVsSimTsp.et920vz55.d25m9y2021.root");
  const TString sInMeas("eventQaPlots.et920vz55pt0230.d8m9y2021.root");
  const TString sInSimPi0("particleGun.noEtCutCheck_withNoTsp.et0100x0100x0100vz55tsp0100pi0.d24m11y2020.root");
  const TString sInSimGam("particleGun.noEtCutCheck_withNoTsp.et0100x0100x0100vz55tsp0206gam.d24m11y2020.root");

  // input histograms
  const TString sHistMeasA("TriggerQA/hTrgTspAll");
  const TString sHistMeasP("TriggerQA/hTrgTspPi0");
  const TString sHistMeasG("TriggerQA/hTrgTspGam");
  const TString sHistSimPi0("match/hTspTrgDetMatchWeight");
  const TString sHistSimGam("match/hTspTrgDetMatchWeight");

  // histogram parameters
  const UInt_t  fColMeasA(923);
  const UInt_t  fColMeasP(859);
  const UInt_t  fColMeasG(899);
  const UInt_t  fColSimP(603);
  const UInt_t  fColSimG(635);
  const UInt_t  fFilMeasA(0);
  const UInt_t  fFilMeasP(3345);
  const UInt_t  fFilMeasG(3354);
  const UInt_t  fFilSimP(0);
  const UInt_t  fFilSimG(0);
  const UInt_t  fMarMeas(20);
  const UInt_t  fMarSim(1);
  const UInt_t  fSizMeas(1.7);
  const UInt_t  fSizSim(1.);
  const UInt_t  fWidMeas(1);
  const UInt_t  fWidSim(2);
  const TString sTitle("");
  const TString sTitleX("TSP");
  const TString sTitleY("arbitrary units");
  const Float_t xPlot[2] = {0., 1.};

  // text parameters
  const TString sLabelMA("all measured triggers");
  const TString sLabelMP("measured #pi^{0}");
  const TString sLabelMG("measured #gamma_{rich}");
  const TString sLabelSP("simulated #pi^{0}");
  const TString sLabelSG("simulated #gamma_{dir}");
  const TString sCollision("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrigger("#pi^{0}/#gamma_{rich} trig., E_{T}^{trg} #in (9, 20) GeV");

  // open files
  TFile *fOut  = new TFile(sOut.Data(), "recreate");
  TFile *fMeas = new TFile(sInMeas.Data(), "read");
  TFile *fSimP = new TFile(sInSimPi0.Data(), "read");
  TFile *fSimG = new TFile(sInSimGam.Data(), "read");
  if (!fOut || !fMeas || !fSimP || !fSimG) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut  = " << fOut  << ", fMeas = " << fMeas << "\n"
         << "       fSimP = " << fSimP << ", fSimG = " << fSimG
         << endl;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hMeasA = (TH1D*) fMeas -> Get(sHistMeasA.Data());
  TH1D *hMeasP = (TH1D*) fMeas -> Get(sHistMeasP.Data());
  TH1D *hMeasG = (TH1D*) fMeas -> Get(sHistMeasG.Data());
  TH1D *hSimP  = (TH1D*) fSimP -> Get(sHistSimPi0.Data());
  TH1D *hSimG  = (TH1D*) fSimG -> Get(sHistSimGam.Data());
  if (!hMeasA || !hMeasP || !hMeasG || !hSimP || !hSimG) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hMeasA = " << hMeasA << ", hMeasP = " << hMeasP << ", hMeasG = " << hMeasG << "\n"
         << "       hSimP  = " << hSimP  << ", hSimG  = " << hSimG
         << endl;
  }
  hMeasA -> SetName("hMeasAll");
  hMeasP -> SetName("hMeasPi0");
  hMeasG -> SetName("hMeasGam");
  hSimP  -> SetName("hSimPi0");
  hSimG  -> SetName("hSimGam");
  cout << "    Grabbed histograms." << endl;

  // normalize histograms
  const Double_t intMA = hMeasA -> Integral();
  const Double_t intSP = hSimP  -> Integral();
  const Double_t intSG = hSimG  -> Integral();
  hMeasA -> Scale(1. / intMA);
  hMeasP -> Scale(1. / intMA);
  hMeasG -> Scale(1. / intMA);
  hSimP  -> Scale(1. / intSP);
  hSimG  -> Scale(1. / intSG); 
  cout << "    Scaled histograms." << endl;

  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fLin(1);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.35);
  const Float_t fOffL(0.);
  const Float_t fSizX(0.04);
  const Float_t fSizY(0.04);
  const Float_t fSizL(0.04);
  hMeasA -> SetLineColor(fColMeasA);
  hMeasA -> SetLineStyle(fLin);
  hMeasA -> SetLineWidth(fWidMeas);
  hMeasA -> SetFillColor(fColMeasA);
  hMeasA -> SetFillStyle(fFilMeasA);
  hMeasA -> SetMarkerColor(fColMeasA);
  hMeasA -> SetMarkerStyle(fMarMeas);
  hMeasA -> SetMarkerSize(fSizMeas);
  hMeasA -> SetTitle(sTitle.Data());
  hMeasA -> SetTitleFont(fTxt);
  hMeasA -> GetXaxis() -> SetTitle(sTitleX.Data());
  hMeasA -> GetXaxis() -> SetTitleFont(fTxt);
  hMeasA -> GetXaxis() -> SetTitleSize(fSizX);
  hMeasA -> GetXaxis() -> SetTitleOffset(fOffX);
  hMeasA -> GetXaxis() -> CenterTitle(fCnt);
  hMeasA -> GetXaxis() -> SetLabelFont(fTxt);
  hMeasA -> GetXaxis() -> SetLabelSize(fSizL);
  hMeasA -> GetYaxis() -> SetTitle(sTitleY.Data());
  hMeasA -> GetYaxis() -> SetTitleFont(fTxt);
  hMeasA -> GetYaxis() -> SetTitleSize(fSizY);
  hMeasA -> GetYaxis() -> SetTitleOffset(fOffY);
  hMeasA -> GetYaxis() -> CenterTitle(fCnt);
  hMeasA -> GetYaxis() -> SetLabelFont(fTxt);
  hMeasA -> GetYaxis() -> SetLabelSize(fSizL);
  hMeasA -> GetYaxis() -> SetLabelOffset(fOffL);
  hMeasP -> SetLineColor(fColMeasP);
  hMeasP -> SetLineStyle(fLin);
  hMeasP -> SetLineWidth(fWidMeas);
  hMeasP -> SetFillColor(fColMeasP);
  hMeasP -> SetFillStyle(fFilMeasP);
  hMeasP -> SetMarkerColor(fColMeasP);
  hMeasP -> SetMarkerStyle(fMarMeas);
  hMeasP -> SetMarkerSize(fSizMeas);
  hMeasP -> SetTitle(sTitle.Data());
  hMeasP -> SetTitleFont(fTxt);
  hMeasP -> GetXaxis() -> SetTitle(sTitleX.Data());
  hMeasP -> GetXaxis() -> SetTitleFont(fTxt);
  hMeasP -> GetXaxis() -> SetTitleSize(fSizX);
  hMeasP -> GetXaxis() -> SetTitleOffset(fOffX);
  hMeasP -> GetXaxis() -> CenterTitle(fCnt);
  hMeasP -> GetXaxis() -> SetLabelFont(fTxt);
  hMeasP -> GetXaxis() -> SetLabelSize(fSizL);
  hMeasP -> GetYaxis() -> SetTitle(sTitleY.Data());
  hMeasP -> GetYaxis() -> SetTitleFont(fTxt);
  hMeasP -> GetYaxis() -> SetTitleSize(fSizY);
  hMeasP -> GetYaxis() -> SetTitleOffset(fOffY);
  hMeasP -> GetYaxis() -> CenterTitle(fCnt);
  hMeasP -> GetYaxis() -> SetLabelFont(fTxt);
  hMeasP -> GetYaxis() -> SetLabelSize(fSizL);
  hMeasP -> GetYaxis() -> SetLabelOffset(fOffL);
  hMeasG -> SetLineColor(fColMeasG);
  hMeasG -> SetLineStyle(fLin);
  hMeasG -> SetLineWidth(fWidMeas);
  hMeasG -> SetFillColor(fColMeasG);
  hMeasG -> SetFillStyle(fFilMeasG);
  hMeasG -> SetMarkerColor(fColMeasG);
  hMeasG -> SetMarkerStyle(fMarMeas);
  hMeasG -> SetMarkerSize(fSizMeas);
  hMeasG -> SetTitle(sTitle.Data());
  hMeasG -> SetTitleFont(fTxt);
  hMeasG -> GetXaxis() -> SetTitle(sTitleX.Data());
  hMeasG -> GetXaxis() -> SetTitleFont(fTxt);
  hMeasG -> GetXaxis() -> SetTitleSize(fSizX);
  hMeasG -> GetXaxis() -> SetTitleOffset(fOffX);
  hMeasG -> GetXaxis() -> CenterTitle(fCnt);
  hMeasG -> GetXaxis() -> SetLabelFont(fTxt);
  hMeasG -> GetXaxis() -> SetLabelSize(fSizL);
  hMeasG -> GetYaxis() -> SetTitle(sTitleY.Data());
  hMeasG -> GetYaxis() -> SetTitleFont(fTxt);
  hMeasG -> GetYaxis() -> SetTitleSize(fSizY);
  hMeasG -> GetYaxis() -> SetTitleOffset(fOffY);
  hMeasG -> GetYaxis() -> CenterTitle(fCnt);
  hMeasG -> GetYaxis() -> SetLabelFont(fTxt);
  hMeasG -> GetYaxis() -> SetLabelSize(fSizL);
  hMeasG -> GetYaxis() -> SetLabelOffset(fOffL);
  hSimP  -> SetLineColor(fColSimP);
  hSimP  -> SetLineStyle(fLin);
  hSimP  -> SetLineWidth(fWidSim);
  hSimP  -> SetFillColor(fColSimP);
  hSimP  -> SetFillStyle(fFilSimP);
  hSimP  -> SetMarkerColor(fColSimP);
  hSimP  -> SetMarkerStyle(fMarSim);
  hSimP  -> SetMarkerSize(fSizSim);
  hSimP  -> SetTitle(sTitle.Data());
  hSimP  -> SetTitleFont(fTxt);
  hSimP  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSimP  -> GetXaxis() -> SetTitleFont(fTxt);
  hSimP  -> GetXaxis() -> SetTitleSize(fSizX);
  hSimP  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSimP  -> GetXaxis() -> CenterTitle(fCnt);
  hSimP  -> GetXaxis() -> SetLabelFont(fTxt);
  hSimP  -> GetXaxis() -> SetLabelSize(fSizL);
  hSimP  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSimP  -> GetYaxis() -> SetTitleFont(fTxt);
  hSimP  -> GetYaxis() -> SetTitleSize(fSizY);
  hSimP  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSimP  -> GetYaxis() -> CenterTitle(fCnt);
  hSimP  -> GetYaxis() -> SetLabelFont(fTxt);
  hSimP  -> GetYaxis() -> SetLabelSize(fSizL);
  hSimP  -> GetYaxis() -> SetLabelOffset(fOffL);
  hSimG  -> SetLineColor(fColSimG);
  hSimG  -> SetLineStyle(fLin);
  hSimG  -> SetLineWidth(fWidSim);
  hSimG  -> SetFillColor(fColSimG);
  hSimG  -> SetFillStyle(fFilSimG);
  hSimG  -> SetMarkerColor(fColSimG);
  hSimG  -> SetMarkerStyle(fMarSim);
  hSimG  -> SetMarkerSize(fSizSim);
  hSimG  -> SetTitle(sTitle.Data());
  hSimG  -> SetTitleFont(fTxt);
  hSimG  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSimG  -> GetXaxis() -> SetTitleFont(fTxt);
  hSimG  -> GetXaxis() -> SetTitleSize(fSizX);
  hSimG  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSimG  -> GetXaxis() -> CenterTitle(fCnt);
  hSimG  -> GetXaxis() -> SetLabelFont(fTxt);
  hSimG  -> GetXaxis() -> SetLabelSize(fSizL);
  hSimG  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSimG  -> GetYaxis() -> SetTitleFont(fTxt);
  hSimG  -> GetYaxis() -> SetTitleSize(fSizY);
  hSimG  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSimG  -> GetYaxis() -> CenterTitle(fCnt);
  hSimG  -> GetYaxis() -> SetLabelFont(fTxt);
  hSimG  -> GetYaxis() -> SetLabelSize(fSizL);
  hSimG  -> GetYaxis() -> SetLabelOffset(fOffL);
  cout << "    Set styles." << endl;

  // create legend
  const UInt_t  fFilLe(0);
  const UInt_t  fColLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  fAlnLe(12);
  const Float_t vtxLe[4] = {0.1, 0.1, 0.3, 0.35};

  TLegend *leg = new TLegend(vtxLe[0], vtxLe[1], vtxLe[2], vtxLe[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnLe);
  leg -> AddEntry(hMeasA, sLabelMA.Data(), "fp");
  leg -> AddEntry(hMeasP, sLabelMP.Data(), "fp");
  leg -> AddEntry(hSimP,  sLabelSP.Data(), "l");
  leg -> AddEntry(hMeasG, sLabelMG.Data(), "fp");
  leg -> AddEntry(hSimG,  sLabelSG.Data(), "l");

  // create text box
  const UInt_t  fFilTx(0);
  const UInt_t  fColTx(0);
  const UInt_t  fLinTx(0);
  const UInt_t  fAlnTx(32);
  const Float_t vtxTx[4] = {0.3, 0.1, 0.5, 0.2};

  TPaveText *txt = new TPaveText(vtxTx[0], vtxTx[1], vtxTx[2], vtxTx[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlnTx);
  txt -> AddText(sCollision.Data());
  txt -> AddText(sTrigger.Data());
  cout << "    Created legend and text box." << endl;

  // create plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fLogY(0);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fBord(2);
  const UInt_t  fMode(0);
  const UInt_t  fFrame(0);
  const Float_t fMarginBig(1.15);
  const Float_t fMarginSmall(0.02);

  TCanvas *cPlot = new TCanvas("cMeasVsSimTSP", "", width, height);
  cPlot  -> SetGrid(fGrid, fGrid);
  cPlot  -> SetTicks(fTick, fTick);
  cPlot  -> SetBorderMode(fMode);
  cPlot  -> SetBorderSize(fBord);
  cPlot  -> SetFrameBorderMode(fFrame);
  cPlot  -> SetLogy(fLogY);
  cPlot  -> SetLeftMargin(fMarginBig);
  cPlot  -> SetBottomMargin(fMarginBig);
  cPlot  -> SetRightMargin(fMarginSmall);
  cPlot  -> SetTopMargin(fMarginSmall);
  cPlot  -> cd();
  hMeasA -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
  hMeasA -> Draw();
  hMeasP -> Draw("hist same");
  hMeasP -> Draw("same");
  hMeasG -> Draw("hist same");
  hMeasG -> Draw("same");
  hSimP  -> Draw("hist ][ same");
  hSimG  -> Draw("hist ][ same");
  leg    -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cPlot  -> Write();
  cPlot  -> Close();

  // save histograms
  fOut   -> cd();
  hMeasA -> Write();
  hMeasP -> Write();
  hMeasG -> Write();
  hSimP  -> Write();
  hSimG  -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut  -> cd();
  fOut  -> Close();
  fMeas -> cd();
  fMeas -> Close();
  fSimP -> cd();
  fSimP -> Close();
  fSimG -> cd();
  fSimG -> Close();
  cout << "  Script finished.\n" << endl;

}

// End ------------------------------------------------------------------------
