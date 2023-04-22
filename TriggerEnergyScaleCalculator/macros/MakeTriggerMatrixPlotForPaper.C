// 'MakeTriggerMatrixPlotForPaper.C'
// Derek Anderson
// 02.21.2023
//
// Use this to make the trigger eT
// matrix for the long paper.

#include <fstream>
#include <iostream>
#include "TH2.h"
#include "TPad.h"
#include "TBox.h"
#include "TMath.h"
#include "TFile.h"
#include "TError.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

// global constants
static const UInt_t NVtx(4);
static const UInt_t NRange(2);



void MakeTriggerMatrixPlotForPaper() {

  // lower verbosity and set stat option
  gErrorIgnoreLevel = kFatal;
  gStyle            -> SetOptStat(0);
  cout << "\n  Beginning matrix plot maker..." << endl;

  // io parameters
  const TString sOutput("trigMatrixForLongPaper_withAprGpcComments.et630pi0vsGam.d20m4y2023.root");
  const TString sInGam("output/2022/November2022/triggerMatrix.forPaper_noShapeWeights_withWeightSmooth.et650x650vz55tsp0206gam.d22m11y2022.root");
  const TString sInPi0("output/2022/November2022/triggerMatrix.forPaper_noShapeWeights_withWeightSmooth.et650x650vz55tsp008pi0.d22m11y2022.root");
  const TString sHistGam("hTriggerMatrix");
  const TString sHistPi0("hTriggerMatrix");
  const TString sNameGam("hMatrixGam");
  const TString sNamePi0("hMatrixPi0");

  // plot parameters
  const TString sTitleGam("#bf{STAR} simulation #gamma");
  const TString sTitlePi0("#bf{STAR} simulation #pi^{0}");
  const TString sTitleGamX("E_{T}^{det-clust} [GeV]");
  const TString sTitlePi0X("E_{T}^{det-clust} [GeV]");
  const TString sTitleGamY("E_{T}^{part-match} [GeV]");
  const TString sTitlePi0Y("E_{T}^{part-match} [GeV]");
  const Float_t xCutRange[NRange]  = {9., 20.};
  const Float_t xyGamRange[NRange] = {6., 30.};
  const Float_t xyPi0Range[NRange] = {6., 30.};

  // open files
  TFile *fOut = new TFile(sOutput.Data(), "recreate");
  TFile *fGam = new TFile(sInGam.Data(),  "read");
  TFile *fPi0 = new TFile(sInPi0.Data(),  "read");
  if (!fPi0 || !fGam || !fOut) {
    cerr << "PANIC: couldn't open input or output file!\n"
         << "       fGam = " << fGam << ", fPi0 = " << fPi0 << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH2D *hGam = (TH2D*) fGam -> Get(sHistGam.Data());
  TH2D *hPi0 = (TH2D*) fPi0 -> Get(sHistPi0.Data());
  if (!hGam || !hPi0) {
    cerr << "PANIC: couldn't grab an input histogram!\n"
         << "       hGam = " << hGam << ", hPi0 = " << hPi0 << "\n"
         << endl;
    return;
  }
  hGam -> SetName(sNameGam.Data());
  hPi0 -> SetName(sNamePi0.Data());
  cout << "    Grabbed matrices." << endl;

  // get matrix minima and maxima
  const Double_t minGam    = hGam -> GetMinimum(0.);
  const Double_t minPi0    = hPi0 -> GetMinimum(0.);
  const Double_t maxGam    = hGam -> GetMaximum();
  const Double_t maxPi0    = hPi0 -> GetMaximum();
  const Double_t minMatrix = TMath::Min(minGam, minPi0);
  const Double_t maxMatrix = TMath::Max(maxGam, maxPi0);
  cout << "    Got minimum & maximum:\n"
       << "      min = " << minMatrix << ", max = " << maxMatrix
       << endl;

  // set styles
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fFil(0);
  const UInt_t  nDivX(509);
  const Float_t fTickL(0.06);
  const Float_t fTtlX(0.05);
  const Float_t fTtlY(0.05);
  const Float_t fTtlZ(0.035);
  const Float_t fLblX(0.06);
  const Float_t fLblY(0.05);
  const Float_t fLblZ(0.035);
  const Float_t fOffX(1.2);
  const Float_t fOffY(1.13);
  const Float_t fOffZ(1.2);
  const Float_t fOffLX(0.001);
  const Float_t fOffLY(0.005);
  const Float_t fOffLZ(0.001);
  const Float_t fOffLblNudge(0.001);
  const TString sTitleZ("Probability Density");

  // set gamma options
  hGam -> SetTitle(sTitleGam.Data());
  hGam -> SetTitleFont(fTxt);
  hGam -> GetXaxis() -> SetRangeUser(xyGamRange[0], xyGamRange[1]);
  hGam -> GetXaxis() -> SetNdivisions(nDivX);
  hGam -> GetXaxis() -> SetTickLength(fTickL);
  hGam -> GetXaxis() -> SetLabelFont(fTxt);
  hGam -> GetXaxis() -> SetLabelSize(fLblX);
  hGam -> GetXaxis() -> SetLabelOffset(fOffLX);
  hGam -> GetXaxis() -> SetTitle(sTitleGamX.Data());
  hGam -> GetXaxis() -> SetTitleFont(fTxt);
  hGam -> GetXaxis() -> SetTitleSize(fTtlX);
  hGam -> GetXaxis() -> SetTitleOffset(fOffX);
  hGam -> GetXaxis() -> CenterTitle(fCnt);
  hGam -> GetYaxis() -> SetRangeUser(xyGamRange[0], xyGamRange[1]);
  hGam -> GetYaxis() -> SetLabelFont(fTxt);
  hGam -> GetYaxis() -> SetLabelSize(fLblY);
  hGam -> GetYaxis() -> SetLabelOffset(fOffLY);
  hGam -> GetYaxis() -> SetTitle(sTitleGamY.Data());
  hGam -> GetYaxis() -> SetTitleFont(fTxt);
  hGam -> GetYaxis() -> SetTitleSize(fTtlY);
  hGam -> GetYaxis() -> SetTitleOffset(fOffY);
  hGam -> GetYaxis() -> CenterTitle(fCnt);
  hGam -> GetZaxis() -> SetRangeUser(minMatrix, maxMatrix);
  hGam -> GetZaxis() -> SetLabelFont(fTxt);
  hGam -> GetZaxis() -> SetLabelSize(fLblZ);
  hGam -> GetZaxis() -> SetLabelOffset(fOffLZ);
  hGam -> GetZaxis() -> SetTitle(sTitleZ.Data());
  hGam -> GetZaxis() -> SetTitleFont(fTxt);
  hGam -> GetZaxis() -> SetTitleSize(fTtlZ);
  hGam -> GetZaxis() -> SetTitleOffset(fOffZ);
  hGam -> GetZaxis() -> CenterTitle(fCnt);

  // set pi0 options
  hPi0 -> SetTitle(sTitlePi0.Data());
  hPi0 -> SetTitleFont(fTxt);
  hPi0 -> GetXaxis() -> SetRangeUser(xyPi0Range[0], xyPi0Range[1]);
  hPi0 -> GetXaxis() -> SetNdivisions(nDivX);
  hPi0 -> GetXaxis() -> SetTickLength(fTickL);
  hPi0 -> GetXaxis() -> SetLabelFont(fTxt);
  hPi0 -> GetXaxis() -> SetLabelSize(fLblX);
  hPi0 -> GetXaxis() -> SetLabelOffset(fOffLX);
  hPi0 -> GetXaxis() -> SetTitle(sTitlePi0X.Data());
  hPi0 -> GetXaxis() -> SetTitleFont(fTxt);
  hPi0 -> GetXaxis() -> SetTitleSize(fTtlX);
  hPi0 -> GetXaxis() -> SetTitleOffset(fOffX);
  hPi0 -> GetXaxis() -> CenterTitle(fCnt);
  hPi0 -> GetYaxis() -> SetRangeUser(xyPi0Range[0], xyPi0Range[1]);
  hPi0 -> GetYaxis() -> SetLabelFont(fTxt);
  hPi0 -> GetYaxis() -> SetLabelSize(fLblY);
  hPi0 -> GetYaxis() -> SetLabelOffset(fOffLY + fOffLblNudge);
  hPi0 -> GetYaxis() -> SetTitle(sTitlePi0Y.Data());
  hPi0 -> GetYaxis() -> SetTitleFont(fTxt);
  hPi0 -> GetYaxis() -> SetTitleSize(fTtlY);
  hPi0 -> GetYaxis() -> SetTitleOffset(fOffY);
  hPi0 -> GetYaxis() -> CenterTitle(fCnt);
  hPi0 -> GetZaxis() -> SetRangeUser(minMatrix, maxMatrix);
  hPi0 -> GetZaxis() -> SetLabelFont(fTxt);
  hPi0 -> GetZaxis() -> SetLabelSize(fLblZ);
  hPi0 -> GetZaxis() -> SetLabelOffset(fOffLZ);
  hPi0 -> GetZaxis() -> SetTitle(sTitleZ.Data());
  hPi0 -> GetZaxis() -> SetTitleFont(fTxt);
  hPi0 -> GetZaxis() -> SetTitleSize(fTtlZ);
  hPi0 -> GetZaxis() -> SetTitleOffset(fOffZ);
  hPi0 -> GetZaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // make cut indicators
  const UInt_t fEtCutCol(1);
  const UInt_t fEtCutFil(3345);
  const UInt_t fEtCutLin(1);
  const UInt_t fEtCutWid(3);

  TBox *bCutGam = new TBox(xCutRange[0], xyGamRange[0], xCutRange[1], xyGamRange[1]);
  TBox *bCutPi0 = new TBox(xCutRange[0], xyPi0Range[0], xCutRange[1], xyPi0Range[1]);
  bCutGam -> SetFillColor(fEtCutCol);
  bCutGam -> SetFillStyle(fEtCutFil);
  bCutGam -> SetLineColor(fEtCutCol);
  bCutGam -> SetLineStyle(fEtCutLin);
  bCutGam -> SetLineWidth(fEtCutWid);
  bCutPi0 -> SetFillColor(fEtCutCol);
  bCutPi0 -> SetFillStyle(fEtCutFil);
  bCutPi0 -> SetLineColor(fEtCutCol);
  bCutPi0 -> SetLineStyle(fEtCutLin);
  bCutPi0 -> SetLineWidth(fEtCutWid);
  cout << "    Made cut indicators." << endl;

  // make plot
  const UInt_t  padWidth       = 850;
  const UInt_t  padHeight      = 700;
  const UInt_t  width          = 2 * padWidth;
  const UInt_t  height         = padHeight;
  const UInt_t  fMode          = 0;
  const UInt_t  fBord          = 2;
  const UInt_t  fGrid          = 0;
  const UInt_t  fTick          = 1;
  const UInt_t  fLogX          = 0;
  const UInt_t  fLogY          = 0;
  const UInt_t  fLogZ          = 1;
  const Float_t fMarginT       = 0.0536313;
  const Float_t fMarginB       = 0.15043017;
  const Float_t fMarginRG      = 0.025;
  const Float_t fMarginRP      = 0.170993865;
  const Float_t fMarginLG      = 0.1504479;
  const Float_t fMarginLP      = 0.005;
  const Float_t xyPadGam[NVtx] = {0.,  0., 0.5, 1.};
  const Float_t xyPadPi0[NVtx] = {0.5, 0., 1.,  1.};

  TCanvas *cMatrix = new TCanvas("cMatrix", "", width, height);
  TPad    *pGam    = new TPad("pGam", "", xyPadGam[0], xyPadGam[1], xyPadGam[2], xyPadGam[3]);
  TPad    *pPi0    = new TPad("pPi0", "", xyPadPi0[0], xyPadPi0[1], xyPadPi0[2], xyPadPi0[3]);

  // set canvas options
  cMatrix -> SetGrid(fGrid, fGrid);
  cMatrix -> SetTicks(fTick, fTick);
  cMatrix -> SetBorderMode(fMode);
  cMatrix -> SetBorderSize(fBord);
  cMatrix -> SetLogx(fLogX);
  cMatrix -> SetLogy(fLogY);
  cMatrix -> SetLogz(fLogZ);
  cMatrix -> SetTopMargin(fMarginT);
  cMatrix -> SetBottomMargin(fMarginB);
  cMatrix -> SetLeftMargin(fMarginLG);
  cMatrix -> SetRightMargin(fMarginRG);

  // set pad options
  pGam    -> SetGrid(fGrid, fGrid);
  pGam    -> SetTicks(fTick, fTick);
  pGam    -> SetBorderMode(fMode);
  pGam    -> SetBorderSize(fBord);
  pGam    -> SetTopMargin(fMarginT);
  pGam    -> SetBottomMargin(fMarginB);
  pGam    -> SetLeftMargin(fMarginLG);
  pGam    -> SetRightMargin(fMarginRG);
  pGam    -> SetLogx(fLogX);
  pGam    -> SetLogy(fLogY);
  pGam    -> SetLogz(fLogZ);
  pPi0    -> SetGrid(fGrid, fGrid);
  pPi0    -> SetTicks(fTick, fTick);
  pPi0    -> SetBorderMode(fMode);
  pPi0    -> SetBorderSize(fBord);
  pPi0    -> SetTopMargin(fMarginT);
  pPi0    -> SetBottomMargin(fMarginB);
  pPi0    -> SetLeftMargin(fMarginLP);
  pPi0    -> SetRightMargin(fMarginRP);
  pPi0    -> SetLogx(fLogX);
  pPi0    -> SetLogy(fLogY);
  pPi0    -> SetLogz(fLogZ);

  // draw plot
  cMatrix -> cd();
  pGam    -> Draw();
  pPi0    -> Draw();
  pGam    -> cd();
  hGam    -> Draw("col");
  bCutGam -> Draw();
  pPi0    -> cd();
  hPi0    -> Draw("colz");
  bCutPi0 -> Draw();
  fOut    -> cd();
  cMatrix -> Write();
  cMatrix -> Close();
  cout << "    Made plot." << endl;

  // save input
  fOut -> cd();
  hGam -> Write();
  hPi0 -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fGam -> cd();
  fGam -> Close();
  fPi0 -> cd();
  fPi0 -> Close();
  cout << "  Finished matrix plot maker!\n" << endl;

}

// End ------------------------------------------------------------------------
