// 'MakeTriggerEnergyPlot.C'
// Derek Anderson
// 11.23.2022
//
// Plots the eTpart, eTdet,
// and eTmatch distributions
// from the TES/R calculation.

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NTrgBins(3);
static const UInt_t NRange(2);
static const UInt_t NTxt(2);
static const UInt_t NLvl(3);
static const UInt_t NVtx(4);



void MakeTriggerEnergyPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning TES/R energy comparison plotting..." << endl;

  // file parameters
  const TString sIn("./output/particleGun.forPaper_withNoWeightMatrix_withTspNoEt.et650x650x0100vz55tsp008pi0.d14m11y2022.root");
  //const TString sIn("./output/particleGun.forPaper_withNoWeightMatrix_withTspNoEt.et650x650x0100vz55tsp0206gam.d14m11y2022.root");
  const TString sOut("trigEtComparison.forPaper_withShapeWeights.et650pi0.d29m11y2022.root");
  //const TString sOut("trigEtComparison.forPaper_withShapeWeights.et650gam.d29m11y2022.root");

  // histogram parameters
  const TString sInPar("particle/hEtTrgParAllWeight");
  const TString sInDet[NTrgBins] = {"match/hEtDet911",   "match/hEtDet1115",   "match/hEtDet1520"};
  const TString sInMat[NTrgBins] = {"match/hEtMatch911", "match/hEtMatch1115", "match/hEtMatch1520"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("E_{T}^{#pi^{0},part}, E_{T}^{#pi^{0},det-clust}, E_{T}^{#pi^{0},part-match} [GeV]");
  //const TString sTitleX("E_{T}^{#gamma,part}, E_{T}^{#gamma,det-clust}, E_{T}^{#gamma,part-match} [GeV]");
  const TString sTitleY("probability density");
  const UInt_t  fColPar(923);
  const UInt_t  fMarPar(29);
  const UInt_t fColDet[NTrgBins]  = {859, 899, 819};
  const UInt_t fColMat[NTrgBins]  = {859, 899, 819};
  //const UInt_t  fColDet[NTrgBins]  = {799, 899, 859};
  //const UInt_t  fColMat[NTrgBins]  = {799, 899, 859};
  const UInt_t  fMarDet[NTrgBins]  = {20, 22, 21};
  const UInt_t  fMarMat[NTrgBins]  = {24, 26, 25};
  const Float_t xPlotRange[NRange] = {6, 30.};

  // norm parameters
  const Double_t norm(1.);
  const Double_t fCheck(0.25);
  const Double_t xNormRange[NTrgBins][NRange] = {{9., 11.}, {11., 15.}, {15., 20.}};

  // text parameters
  const UInt_t fColLeg(923);
  const TString sTxt[NTxt]     = {"STAR simulation #pi^{0}, 6 < E_{T}^{#pi^{0},part,det-clust} < 50 GeV", "|#eta^{#pi^{0}}| < 0.9, TSP #in (0., 0.08)"};
  //const TString sTxt[NTxt]     = {"STAR simulation #gamma, 6 < E_{T}^{#gamma,part,det-clust} < 50 GeV", "|#eta^{#gamma}| < 0.9, TSP #in (0.2, 0.6)"};
  const TString sLvl[NLvl]     = {"Generated", "Reco.", "Matched"};
  const TString sTrg[NTrgBins] = {"E_{T}^{#pi^{0},det-clust} #in (9, 11) GeV", "E_{T}^{#pi^{0},det-clust} #in (11, 15) GeV", "E_{T}^{#pi^{0},det-clust} #in (15, 20) GeV"};
  //const TString sTrg[NTrgBins] = {"E_{T}^{#gamma,det-clust} #in (9, 11) GeV", "E_{T}^{#gamma,det-clust} #in (11, 15) GeV", "E_{T}^{#gamma,det-clust} #in (15, 20) GeV"};

  // open files
  TFile *fIn  = new TFile(sIn.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn || !fOut) {
    cerr << "PANIC: couldn't open input or output file!\n"
         << "       fIn = " << fIn << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hPar = (TH1D*) fIn -> Get(sInPar.Data());
  if (!hPar) {
    cerr << "PANIC: couldn't grab input particle histogram!\n" << endl;
    return;
  }

  TH1D *hDet[NTrgBins];
  TH1D *hMat[NTrgBins];
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hDet[iTrg] = (TH1D*) fIn -> Get(sInDet[iTrg].Data());
    hMat[iTrg] = (TH1D*) fIn -> Get(sInMat[iTrg].Data());
    if (!hDet[iTrg] || !hMat[iTrg]) {
      cerr << "PANIC: couldn't grab input detector or matched histogram!\n"
           << "       hDet[" << iTrg << "] = " << hDet[iTrg] << ", hMat[" << iTrg << "] = " << hMat[iTrg] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // normalize particle histogram
  const Double_t intPar  = hPar -> Integral();
  const Double_t normPar = norm / intPar;
  hPar -> Scale(normPar);

  // normalize detector/matched histograms
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {

    // determine norm range
    const Double_t xWidth = hPar -> GetBinWidth(1);
    const Double_t xCheck = xWidth * fCheck;
    const UInt_t   iStart = hPar -> FindBin(xNormRange[iTrg][0] + xCheck);
    const UInt_t   iStop  = hPar -> FindBin(xNormRange[iTrg][1] - xCheck);

    // get norms
    const Double_t intParBin = hPar       -> Integral(iStart, iStop);
    const Double_t intDet    = hDet[iTrg] -> Integral();
    const Double_t intMat    = hMat[iTrg] -> Integral();
    const Double_t normDet   = intParBin / intDet;
    const Double_t normMat   = intParBin / intMat;
    hDet[iTrg] -> Scale(normDet);
    hMat[iTrg] -> Scale(normMat);
  }
  cout << "    Normalized histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLbl(0.04);
  const Float_t fTtl(0.04);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.2);
  hPar -> SetMarkerColor(fColPar);
  hPar -> SetMarkerStyle(fMarPar);
  hPar -> SetLineColor(fColPar);
  hPar -> SetLineStyle(fLin);
  hPar -> SetLineWidth(fWid);
  hPar -> SetFillColor(fColPar);
  hPar -> SetFillStyle(fFil);
  hPar -> SetTitle(sTitle.Data());
  hPar -> SetTitleFont(fTxt);
  hPar -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hPar -> GetXaxis() -> SetLabelFont(fTxt);
  hPar -> GetXaxis() -> SetLabelSize(fLbl);
  hPar -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPar -> GetXaxis() -> SetTitleSize(fTtl);
  hPar -> GetXaxis() -> SetTitleFont(fTxt);
  hPar -> GetXaxis() -> SetTitleOffset(fOffX);
  hPar -> GetXaxis() -> CenterTitle(fCnt);
  hPar -> GetYaxis() -> SetLabelFont(fTxt);
  hPar -> GetYaxis() -> SetLabelSize(fLbl);
  hPar -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPar -> GetYaxis() -> SetTitleSize(fTtl);
  hPar -> GetYaxis() -> SetTitleFont(fTxt);
  hPar -> GetYaxis() -> SetTitleOffset(fOffY);
  hPar -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hDet[iTrg] -> SetMarkerColor(fColDet[iTrg]);
    hDet[iTrg] -> SetMarkerStyle(fMarDet[iTrg]);
    hDet[iTrg] -> SetLineColor(fColDet[iTrg]);
    hDet[iTrg] -> SetLineStyle(fLin);
    hDet[iTrg] -> SetLineWidth(fWid);
    hDet[iTrg] -> SetFillColor(fColDet[iTrg]);
    hDet[iTrg] -> SetFillStyle(fFil);
    hDet[iTrg] -> SetTitle(sTitle.Data());
    hDet[iTrg] -> SetTitleFont(fTxt);
    hDet[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hDet[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hDet[iTrg] -> GetXaxis() -> SetLabelSize(fLbl);
    hDet[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDet[iTrg] -> GetXaxis() -> SetTitleSize(fTtl);
    hDet[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hDet[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hDet[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hDet[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hDet[iTrg] -> GetYaxis() -> SetLabelSize(fLbl);
    hDet[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hDet[iTrg] -> GetYaxis() -> SetTitleSize(fTtl);
    hDet[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hDet[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hDet[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hMat[iTrg] -> SetMarkerColor(fColMat[iTrg]);
    hMat[iTrg] -> SetMarkerStyle(fMarMat[iTrg]);
    hMat[iTrg] -> SetLineColor(fColMat[iTrg]);
    hMat[iTrg] -> SetLineStyle(fLin);
    hMat[iTrg] -> SetLineWidth(fWid);
    hMat[iTrg] -> SetFillColor(fColMat[iTrg]);
    hMat[iTrg] -> SetFillStyle(fFil);
    hMat[iTrg] -> SetTitle(sTitle.Data());
    hMat[iTrg] -> SetTitleFont(fTxt);
    hMat[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hMat[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hMat[iTrg] -> GetXaxis() -> SetLabelSize(fLbl);
    hMat[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hMat[iTrg] -> GetXaxis() -> SetTitleSize(fTtl);
    hMat[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hMat[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hMat[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hMat[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hMat[iTrg] -> GetYaxis() -> SetLabelSize(fLbl);
    hMat[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hMat[iTrg] -> GetYaxis() -> SetTitleSize(fTtl);
    hMat[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hMat[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hMat[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }

  // create histograms for legend
  TH1D *hParL = (TH1D*) hPar    -> Clone();
  TH1D *hDetL = (TH1D*) hDet[0] -> Clone();
  TH1D *hMatL = (TH1D*) hMat[0] -> Clone();
  hParL -> SetName("hParticleLegend");
  hDetL -> SetName("hDetectorLegend");
  hMatL -> SetName("hMatchedLegend");
  hParL -> SetMarkerColor(fColLeg);
  hParL -> SetFillColor(fColLeg);
  hParL -> SetLineColor(fColLeg);
  hDetL -> SetMarkerColor(fColLeg);
  hDetL -> SetFillColor(fColLeg);
  hDetL -> SetLineColor(fColLeg);
  hMatL -> SetMarkerColor(fColLeg);
  hMatL -> SetFillColor(fColLeg);
  hMatL -> SetLineColor(fColLeg);
  cout << "    Set styles." << endl;

  // create legend/text box
  const UInt_t  fAln        = 12;
  const UInt_t  fColBox     = 0;
  const UInt_t  fFilBox     = 0;
  const UInt_t  fLinBox     = 1;
  const UInt_t  nColumn     = 2;
  const UInt_t  nObjL       = NTrgBins;
  const UInt_t  nObjT       = NTxt;
  const Float_t fObj        = 0.05;
  const Float_t hObjL       = fObj * nObjL;
  const Float_t hObjT       = fObj * nObjT;
  const Float_t yLeg        = 0.1 + hObjL;
  const Float_t yTxt        = 0.1 + hObjT;
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, yLeg};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, yTxt};

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  leg -> SetFillColor(fColBox);
  leg -> SetFillStyle(fFilBox);
  leg -> SetLineColor(fColBox);
  leg -> SetLineStyle(fLinBox);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> SetNColumns(nColumn);
  leg -> AddEntry(hParL,   sLvl[0].Data(), "pf");
  leg -> AddEntry(hDet[0], sTrg[0].Data(), "pf");
  leg -> AddEntry(hDetL,   sLvl[1].Data(), "pf");
  leg -> AddEntry(hDet[1], sTrg[1].Data(), "pf");
  leg -> AddEntry(hMatL,   sLvl[2].Data(), "pf");
  leg -> AddEntry(hDet[2], sTrg[2].Data(), "pf");

  TPaveText *txt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  txt -> SetFillColor(fColBox);
  txt -> SetFillStyle(fFilBox);
  txt -> SetLineColor(fColBox);
  txt -> SetLineStyle(fLinBox);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  for (UInt_t iTxt = 0; iTxt < NTxt; iTxt++) {
    txt -> AddText(sTxt[iTxt].Data());
  }
  cout << "    Made legend and text box." << endl;

  // save histograms
  fOut -> cd();
  hPar -> Write();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hDet[iTrg] -> Write();
    hMat[iTrg] -> Write();
  }
  hParL -> Write();
  hDetL -> Write();
  hMatL -> Write();
  cout << "    Saved histograms." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogY(1);
  const Float_t fMarginT(0.02);
  const Float_t fMarginR(0.02);
  const Float_t fMarginB(0.10);
  const Float_t fMarginL(0.12);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot -> SetGrid(fGrid, fGrid);
  cPlot -> SetTicks(fTick, fTick);
  cPlot -> SetBorderMode(fMode);
  cPlot -> SetBorderSize(fBord);
  cPlot -> SetTopMargin(fMarginT);
  cPlot -> SetRightMargin(fMarginR);
  cPlot -> SetBottomMargin(fMarginB);
  cPlot -> SetLeftMargin(fMarginL);
  cPlot -> SetLogy(fLogY);
  cPlot -> cd();
  hPar  -> Draw();
  for (UInt_t iTrg = 0; iTrg < NTrgBins; iTrg++) {
    hDet[iTrg] -> Draw("same");
    hMat[iTrg] -> Draw("same");
  }
  leg   -> Draw();
  txt   -> Draw();
  fOut  -> cd();
  cPlot -> Write();
  cPlot -> Close();
  cout << "    Made plot." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished making TES/R energy comparison plot!\n" << endl;

}

// end ------------------------------------------------------------------------
