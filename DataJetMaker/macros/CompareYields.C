// 'CompareYields.C'
// Derek Anderson
// 05.10.2020
//
// Compares yields extracted from Run 9
// pp data for gamma-jet analysis against
// yields measured for gamma-hadron analysis.


#include <iostream>
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"

using namespace std;


// global constants
static const UInt_t  NPoints(8);
static const Bool_t  DoScaleAS(true);
static const Bool_t  DoScaleNS(true);
static const Float_t ScaleAS(1. / 0.87);
static const Float_t ScaleNS(1. / 0.87);



void CompareYields() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning yield comparison..." << endl;


  // io parameters
  const TString sOutput("pp200r9.comparingYields_effCorrectedNewYields087_withRatios.et1220pt1220vz55pi0.d12m5y2020.root");
  const TString sNewYieldAS("gNewYieldAS");
  const TString sNewYieldNS("gNewYieldNS");
  const TString sOldStatAS("gOldStatAS");
  const TString sOldStatNS("gOldStatNS");
  const TString sOldSysAS("gOldSysAS");
  const TString sOldSysNS("gOldSysNS");
  const TString sRatioStatAS("gRatioStatAS");
  const TString sRatioStatNS("gRatioStatNS");
  const TString sRatioSysAS("gRatioSysAS");
  const TString sRatioSysNS("gRatioSysNS");

  // AS graph parameters
  const UInt_t  fColNewAS(623);
  const UInt_t  fColOldAS(899);
  const UInt_t  fMarNewAS(4);
  const UInt_t  fMarOldAS(8);
  const UInt_t  fFilNewAS(0);
  const UInt_t  fFilStatAS(0);
  const UInt_t  fFilSysAS(3001);
  const UInt_t  fLinNewAS(1);
  const UInt_t  fLinOldAS(1);
  const TString sTitleAS("Away Side, |#Delta#varphi - #pi| < 1.3");
  const TString sTitleXas("z_{T}^{trk}");
  const TString sTitleYas("D^{AS}_{pp}");
  const TString sTitleRas("new / old");

  // NS graph parameters
  const UInt_t  fColNewNS(591);
  const UInt_t  fColOldNS(859);
  const UInt_t  fMarNewNS(4);
  const UInt_t  fMarOldNS(8);
  const UInt_t  fFilNewNS(0);
  const UInt_t  fFilStatNS(0);
  const UInt_t  fFilSysNS(3001);
  const UInt_t  fLinNewNS(1);
  const UInt_t  fLinOldNS(1);
  const TString sTitleNS("Near Side, |#Delta#varphi| < 1.3");
  const TString sTitleXns("z_{T}^{trk}");
  const TString sTitleYns("D^{NS}_{pp}");
  const TString sTitleRns("new / old");

  // text parameters
  const UInt_t  fLeg(0);
  const UInt_t  fTxt(42);
  const UInt_t  fColTxt(1);
  const UInt_t  fAlign(12);
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrigger("#pi^{0} trig., E_{T}^{trg} #in (12, 20) GeV");
  const TString sTrack("dca < 3 cm, p_{T}^{trk} #in (1.2, 20) GeV/c");
  const TString sNewYields("new data, #gamma-jet [stat.] (corrected for #epsilon_{trk} ~ 0.87)");
  const TString sOldStat("old data, #gamma-h^{#pm} [stat.] (efficiency corrected)");
  const TString sOldSys("old data, #gamma-h^{#pm} [sys.] (efficiency corrected)");
  const TString sRatStatAS("Away Side ratio [stat.]");
  const TString sRatStatNS("Near Side ratio [stat.]");
  const TString sRatSysAS("Away Side ratio [sys.]");
  const TString sRatSysNS("Near Side ratio [sys.]");
  const Float_t xyText[4]  = {0.3, 0.1, 0.5, 0.3};
  const Float_t xyLegAS[4] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyLegNS[4] = {0.1, 0.1, 0.3, 0.3};

  // frame parameters
  const UInt_t  fCnt(1);
  const UInt_t  nXplot(20);
  const UInt_t  nYplot(300);
  const UInt_t  nRplot(20);
  const TString sYieldAS("hYieldAS");
  const TString sYieldNS("hYieldNS");
  const TString sRatioAS("hRatioAS");
  const TString sRatioNS("hRatioNS");
  const Float_t fTitleSize(0.04);
  const Float_t fTitleOffX(1.);
  const Float_t fTitleOffY(1.1);
  const Float_t fLabSize(0.03);
  const Float_t xPlot[2]  = {0., 1.};
  const Float_t yPlot[2]  = {0.003, 33.};
  const Float_t rPlot[2]  = {0.47, 1.53};
  const Float_t xRange[2] = {0., 1.};
  const Float_t yRange[2] = {0.003, 33.};
  const Float_t rRange[2] = {0.47, 1.53};

  // plot parameters
  const UInt_t  fTick(1);
  const UInt_t  fGrid(0);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const TString sYields("cYields");
  const TString sRatios("cRatios");
  const TString sPadYieldAS("pYieldsAS");
  const TString sPadYieldNS("pYieldsNS");
  const TString sPadRatioAS("pRatiosAS");
  const TString sPadRatioNS("pRatiosNS");
  const Float_t width(1500.);
  const Float_t height(750.);
  const Float_t fMargin(0.02);
  const Float_t xyPadAS[4] = {0.5, 0., 1., 1.};
  const Float_t xyPadNS[4] = {0., 0., 0.5, 1.};

  // line parameters
  const UInt_t  fLinCol(1);
  const UInt_t  fLinSty(9);
  const UInt_t  fLinSiz(2);
  const Float_t lineVal(1.);


  // new AS yields
  Double_t xNewAS[NPoints]      = {0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85};
  Double_t yNewAS[NPoints]      = {7.854, 3.130, 1.496, 0.683, 0.521, 0.208, 0.142, 0.099};
  Double_t xErrLoNewAS[NPoints] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  Double_t xErrHiNewAS[NPoints] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  Double_t yErrLoNewAS[NPoints] = {0.166, 0.096, 0.067, 0.045, 0.039, 0.026, 0.021, 0.019};
  Double_t yErrHiNewAS[NPoints] = {0.166, 0.096, 0.067, 0.045, 0.039, 0.026, 0.021, 0.019};

  // new NS yields
  Double_t xNewNS[NPoints]      = {0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85};
  Double_t yNewNS[NPoints]      = {3.234, 1.338, 0.435, 0.279, 0.096, 0.054, 0.068, 0.010};
  Double_t xErrLoNewNS[NPoints] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  Double_t xErrHiNewNS[NPoints] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  Double_t yErrLoNewNS[NPoints] = {0.120, 0.064, 0.037, 0.029, 0.017, 0.013, 0.015, 0.008};
  Double_t yErrHiNewNS[NPoints] = {0.120, 0.064, 0.037, 0.029, 0.017, 0.013, 0.015, 0.008};

  // old AS yields
  const Double_t xOldAS[NPoints]       = {0.126, 0.23, 0.33, 0.44, 0.54, 0.65, 0.75, 0.85};
  const Double_t yOldAS[NPoints]       = {9.45625, 3.81209, 1.82786, 0.865953, 0.646539, 0.315956, 0.204786, 0.118483};
  const Double_t xStatLoOldAS[NPoints] = {0.026, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.05};
  const Double_t xStatHiOldAS[NPoints] = {0.074, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.05};
  const Double_t yStatLoOldAS[NPoints] = {0.061169, 0.0348081, 0.0233488, 0.0511065, 0.0198469, 0.010498, 0.00783924, 0.00589885};
  const Double_t yStatHiOldAS[NPoints] = {0.061169, 0.0348081, 0.0233488, 0.0511065, 0.0198469, 0.010498, 0.00783924, 0.00589885};
  const Double_t xSysLoOldAS[NPoints]  = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  const Double_t xSysHiOldAS[NPoints]  = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  const Double_t ySysLoOldAS[NPoints]  = {1.13004, 0.466637, 0.136008, 0.0804235, 0.0696328, 0.0229428, 0.0523296, 0.};
  const Double_t ySysHiOldAS[NPoints]  = {0.302744, 0.11572, 0.0694605, 0.0242231, 0.0095178, 0.0219204, 0.00819245, 0.0551975};

  // old NS yields
  const Double_t xOldNS[NPoints]       = {0.126, 0.23, 0.33, 0.44, 0.54, 0.65, 0.75, 0.85};
  const Double_t yOldNS[NPoints]       = {4.25761, 1.55359, 0.678129, 0.351062, 0.134574, 0.0672869, 0.0643614};
  const Double_t xStatLoOldNS[NPoints] = {0.026, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.05};
  const Double_t xStatHiOldNS[NPoints] = {0.074, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.05};
  const Double_t yStatLoOldNS[NPoints] = {0.0472879, 0.0234861, 0.0145112, 0.0496109, 0.0156244, 0.00612046, 0.00451355, 0.00247488};
  const Double_t yStatHiOldNS[NPoints] = {0.0472879, 0.0234861, 0.0145112, 0.0496109, 0.0156244, 0.00612046, 0.00451355, 0.00247488};
  const Double_t xSysLoOldNS[NPoints]  = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  const Double_t xSysHiOldNS[NPoints]  = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  const Double_t ySysLoOldNS[NPoints]  = {0.339517, 0.0964103, 0.0403867, 0.0258295, 0.0172783, 0.00348362, 0.00879478, 0.};
  const Double_t ySysHiOldNS[NPoints]  = {0.166878, 0.0521669, 0.0659449, 0.0542094, 0.010947, 0.00705247, 0.000281727, 0.00949034};


  // scale new yields if necessary
  if (DoScaleAS) {
    for (UInt_t iPointAS = 0; iPointAS < NPoints; iPointAS++) {
      yNewAS[iPointAS]      = yNewAS[iPointAS] * ScaleAS;
      yErrLoNewAS[iPointAS] = yErrLoNewAS[iPointAS] * ScaleAS;
      yErrHiNewAS[iPointAS] = yErrHiNewAS[iPointAS] * ScaleAS;
    }
    cout << "    Scaled new AS yields by " << ScaleAS << endl;
  }
  if (DoScaleNS) {
    for (UInt_t iPointNS = 0; iPointNS < NPoints; iPointNS++) {
      yNewNS[iPointNS]      = yNewNS[iPointNS] * ScaleNS;
      yErrLoNewNS[iPointNS] = yErrLoNewNS[iPointNS] * ScaleNS;
      yErrHiNewNS[iPointNS] = yErrHiNewNS[iPointNS] * ScaleNS;
    }
    cout << "    Scaled new NS yields by " << ScaleNS << endl;
  }

  // calculate ratios
  Double_t yRatioAS[NPoints];
  Double_t yRatioNS[NPoints];
  Double_t yRatioStatLoAS[NPoints];
  Double_t yRatioStatHiAS[NPoints];
  Double_t yRatioStatLoNS[NPoints];
  Double_t yRatioStatHiNS[NPoints];
  Double_t yRatioSysLoAS[NPoints];
  Double_t yRatioSysHiAS[NPoints];
  Double_t yRatioSysLoNS[NPoints];
  Double_t yRatioSysHiNS[NPoints];
  for (UInt_t iPointAS = 0; iPointAS < NPoints; iPointAS++) {
    const Double_t valAS       = yNewAS[iPointAS] / yOldAS[iPointAS];
    const Double_t errStatLoAS = valAS * TMath::Sqrt(TMath::Power((yErrLoNewAS[iPointAS] / yNewAS[iPointAS]), 2.) + TMath::Power((yStatLoOldAS[iPointAS] / yOldAS[iPointAS]), 2.));
    const Double_t errStatHiAS = valAS * TMath::Sqrt(TMath::Power((yErrHiNewAS[iPointAS] / yNewAS[iPointAS]), 2.) + TMath::Power((yStatHiOldAS[iPointAS] / yOldAS[iPointAS]), 2.));
    const Double_t errSysLoAS  = valAS * TMath::Sqrt(TMath::Power((yErrLoNewAS[iPointAS] / yNewAS[iPointAS]), 2.) + TMath::Power((ySysLoOldAS[iPointAS] / yOldAS[iPointAS]), 2.));
    const Double_t errSysHiAS  = valAS * TMath::Sqrt(TMath::Power((yErrHiNewAS[iPointAS] / yNewAS[iPointAS]), 2.) + TMath::Power((ySysHiOldAS[iPointAS] / yOldAS[iPointAS]), 2.));
    yRatioAS[iPointAS] = valAS;
    yRatioStatLoAS[iPointAS] = errStatLoAS;
    yRatioStatHiAS[iPointAS] = errStatHiAS;
    yRatioSysLoAS[iPointAS]  = errSysLoAS;
    yRatioSysHiAS[iPointAS]  = errSysHiAS;
  }
  for (UInt_t iPointNS = 0; iPointNS < NPoints; iPointNS++) {
    const Double_t valNS       = yNewNS[iPointNS] / yOldNS[iPointNS];
    const Double_t errStatLoNS = valNS * TMath::Sqrt(TMath::Power((yErrLoNewNS[iPointNS] / yNewNS[iPointNS]), 2.) + TMath::Power((yStatLoOldNS[iPointNS] / yOldNS[iPointNS]), 2.));
    const Double_t errStatHiNS = valNS * TMath::Sqrt(TMath::Power((yErrHiNewNS[iPointNS] / yNewNS[iPointNS]), 2.) + TMath::Power((yStatHiOldNS[iPointNS] / yOldNS[iPointNS]), 2.));
    const Double_t errSysLoNS  = valNS * TMath::Sqrt(TMath::Power((yErrLoNewNS[iPointNS] / yNewNS[iPointNS]), 2.) + TMath::Power((ySysLoOldNS[iPointNS] / yOldNS[iPointNS]), 2.));
    const Double_t errSysHiNS  = valNS * TMath::Sqrt(TMath::Power((yErrHiNewNS[iPointNS] / yNewNS[iPointNS]), 2.) + TMath::Power((ySysHiOldNS[iPointNS] / yOldNS[iPointNS]), 2.));
    yRatioNS[iPointNS] = valNS;
    yRatioStatLoNS[iPointNS] = errStatLoNS;
    yRatioStatHiNS[iPointNS] = errStatHiNS;
    yRatioSysLoNS[iPointNS]  = errSysLoNS;
    yRatioSysHiNS[iPointNS]  = errSysHiNS;
  }
  cout << "    Calculated ratios." << endl;


  // declare yield graphs
  TGraphAsymmErrors *gNewYieldAS = new TGraphAsymmErrors(NPoints, xNewAS, yNewAS, xErrLoNewAS, xErrHiNewAS, yErrLoNewAS, yErrHiNewAS);
  TGraphAsymmErrors *gNewYieldNS = new TGraphAsymmErrors(NPoints, xNewNS, yNewNS, xErrLoNewNS, xErrHiNewNS, yErrLoNewNS, yErrHiNewNS);
  TGraphAsymmErrors *gOldStatAS  = new TGraphAsymmErrors(NPoints, xOldAS, yOldAS, xStatLoOldAS, xStatHiOldAS, yStatLoOldAS, yStatHiOldAS);
  TGraphAsymmErrors *gOldStatNS  = new TGraphAsymmErrors(NPoints, xOldNS, yOldNS, xStatLoOldNS, xStatHiOldNS, yStatLoOldNS, yStatHiOldNS);
  TGraphAsymmErrors *gOldSysAS   = new TGraphAsymmErrors(NPoints, xOldAS, yOldAS, xSysLoOldAS, xSysHiOldAS, ySysLoOldAS, ySysHiOldAS);
  TGraphAsymmErrors *gOldSysNS   = new TGraphAsymmErrors(NPoints, xOldNS, yOldNS, xSysLoOldNS, xSysHiOldNS, ySysLoOldNS, ySysLoOldNS);

  // declare ratio graphs
  TGraphAsymmErrors *gRatioStatAS = new TGraphAsymmErrors(NPoints, xOldAS, yRatioAS, xStatLoOldAS, xStatHiOldAS, yRatioStatLoAS, yRatioStatHiAS);
  TGraphAsymmErrors *gRatioStatNS = new TGraphAsymmErrors(NPoints, xOldNS, yRatioNS, xStatLoOldNS, xStatHiOldNS, yRatioStatLoNS, yRatioStatHiNS);
  TGraphAsymmErrors *gRatioSysAS  = new TGraphAsymmErrors(NPoints, xOldAS, yRatioAS, xSysLoOldAS, xSysHiOldAS, yRatioSysLoAS, yRatioSysHiAS);
  TGraphAsymmErrors *gRatioSysNS  = new TGraphAsymmErrors(NPoints, xOldNS, yRatioNS, xSysLoOldNS, xSysHiOldNS, yRatioSysLoNS, yRatioSysHiNS);
  cout << "    Declared graphs." << endl;


  // set new styles
  gNewYieldAS -> SetName(sNewYieldAS.Data());
  gNewYieldAS -> SetTitle(sTitleAS.Data());
  gNewYieldAS -> SetMarkerColor(fColNewAS);
  gNewYieldAS -> SetMarkerStyle(fMarNewAS);
  gNewYieldAS -> SetLineColor(fColNewAS);
  gNewYieldAS -> SetLineStyle(fLinNewAS);
  gNewYieldAS -> SetFillColor(fColNewAS);
  gNewYieldAS -> SetFillStyle(fFilNewAS);
  gNewYieldNS -> SetName(sNewYieldNS.Data());
  gNewYieldNS -> SetTitle(sTitleNS.Data());
  gNewYieldNS -> SetMarkerColor(fColNewNS);
  gNewYieldNS -> SetMarkerStyle(fMarNewNS);
  gNewYieldNS -> SetLineColor(fColNewNS);
  gNewYieldNS -> SetLineStyle(fLinNewNS);
  gNewYieldNS -> SetFillColor(fColNewNS);
  gNewYieldNS -> SetFillStyle(fFilNewNS);

  // set old styles
  gOldStatAS -> SetName(sOldStatAS.Data());
  gOldStatAS -> SetTitle(sTitleAS.Data());
  gOldStatAS -> SetMarkerColor(fColOldAS);
  gOldStatAS -> SetMarkerStyle(fMarOldAS);
  gOldStatAS -> SetLineColor(fColOldAS);
  gOldStatAS -> SetLineStyle(fLinOldAS);
  gOldStatAS -> SetFillColor(fColOldAS);
  gOldStatAS -> SetFillStyle(fFilStatAS);
  gOldStatNS -> SetName(sOldStatNS.Data());
  gOldStatNS -> SetTitle(sTitleNS.Data());
  gOldStatNS -> SetMarkerColor(fColOldNS);
  gOldStatNS -> SetMarkerStyle(fMarOldNS);
  gOldStatNS -> SetLineColor(fColOldNS);
  gOldStatNS -> SetLineStyle(fLinOldNS);
  gOldStatNS -> SetFillColor(fColOldNS);
  gOldStatNS -> SetFillStyle(fFilStatNS);
  gOldSysAS  -> SetName(sOldSysAS.Data());
  gOldSysAS  -> SetTitle(sTitleAS.Data());
  gOldSysAS  -> SetMarkerColor(fColOldAS);
  gOldSysAS  -> SetMarkerStyle(fMarOldAS);
  gOldSysAS  -> SetLineColor(fColOldAS);
  gOldSysAS  -> SetLineStyle(fLinOldAS);
  gOldSysAS  -> SetFillColor(fColOldAS);
  gOldSysAS  -> SetFillStyle(fFilSysAS);
  gOldSysNS  -> SetName(sOldSysNS.Data());
  gOldSysNS  -> SetTitle(sTitleNS.Data());
  gOldSysNS  -> SetMarkerColor(fColOldNS);
  gOldSysNS  -> SetMarkerStyle(fMarOldNS);
  gOldSysNS  -> SetLineColor(fColOldNS);
  gOldSysNS  -> SetLineStyle(fLinOldNS);
  gOldSysNS  -> SetFillColor(fColOldNS);
  gOldSysNS  -> SetFillStyle(fFilSysNS);

  // set ratio styles
  gRatioStatAS -> SetName(sRatioStatAS.Data());
  gRatioStatAS -> SetTitle(sTitleAS.Data());
  gRatioStatAS -> SetMarkerColor(fColOldAS);
  gRatioStatAS -> SetMarkerStyle(fMarOldAS);
  gRatioStatAS -> SetLineColor(fColOldAS);
  gRatioStatAS -> SetLineStyle(fLinOldAS);
  gRatioStatAS -> SetFillColor(fColOldAS);
  gRatioStatAS -> SetFillStyle(fFilStatAS);
  gRatioStatNS -> SetName(sRatioStatNS.Data());
  gRatioStatNS -> SetTitle(sTitleNS.Data());
  gRatioStatNS -> SetMarkerColor(fColOldNS);
  gRatioStatNS -> SetMarkerStyle(fMarOldNS);
  gRatioStatNS -> SetLineColor(fColOldNS);
  gRatioStatNS -> SetLineStyle(fLinOldNS);
  gRatioStatNS -> SetFillColor(fColOldNS);
  gRatioStatNS -> SetFillStyle(fFilStatNS);
  gRatioSysAS  -> SetName(sRatioSysAS.Data());
  gRatioSysAS  -> SetTitle(sTitleAS.Data());
  gRatioSysAS  -> SetMarkerColor(fColOldAS);
  gRatioSysAS  -> SetMarkerStyle(fMarOldAS);
  gRatioSysAS  -> SetLineColor(fColOldAS);
  gRatioSysAS  -> SetLineStyle(fLinOldAS);
  gRatioSysAS  -> SetFillColor(fColOldAS);
  gRatioSysAS  -> SetFillStyle(fFilSysAS);
  gRatioSysNS  -> SetName(sRatioSysNS.Data());
  gRatioSysNS  -> SetTitle(sTitleNS.Data());
  gRatioSysNS  -> SetMarkerColor(fColOldNS);
  gRatioSysNS  -> SetMarkerStyle(fMarOldNS);
  gRatioSysNS  -> SetLineColor(fColOldNS);
  gRatioSysNS  -> SetLineStyle(fLinOldNS);
  gRatioSysNS  -> SetFillColor(fColOldNS);
  gRatioSysNS  -> SetFillStyle(fFilSysNS);
  cout << "    Set styles." << endl;


  // make legend and text
  TPaveText *pt = new TPaveText(xyText[0], xyText[1], xyText[2], xyText[3], "NDC NB");
  pt -> SetFillColor(fLeg);
  pt -> SetFillStyle(fLeg);
  pt -> SetLineColor(fLeg);
  pt -> SetLineStyle(fLeg);
  pt -> SetTextFont(fTxt);
  pt -> SetTextColor(fColTxt);
  pt -> SetTextAlign(fAlign);
  pt -> AddText(sSystem.Data());
  pt -> AddText(sTrigger.Data());
  pt -> AddText(sTrack.Data());

  TLegend *legAS = new TLegend(xyLegAS[0], xyLegAS[1], xyLegAS[2], xyLegAS[3]);
  legAS -> SetFillColor(fLeg);
  legAS -> SetFillStyle(fLeg);
  legAS -> SetLineColor(fLeg);
  legAS -> SetLineStyle(fLeg);
  legAS -> SetTextFont(fTxt);
  legAS -> SetTextColor(fColTxt);
  legAS -> SetTextAlign(fAlign);
  legAS -> AddEntry(gNewYieldAS, sNewYields.Data());
  legAS -> AddEntry(gOldStatAS, sOldStat.Data());
  legAS -> AddEntry(gOldSysAS, sOldSys.Data());

  TLegend *legNS = new TLegend(xyLegNS[0], xyLegNS[1], xyLegNS[2], xyLegNS[3]);
  legNS -> SetFillColor(fLeg);
  legNS -> SetFillStyle(fLeg);
  legNS -> SetLineColor(fLeg);
  legNS -> SetLineStyle(fLeg);
  legNS -> SetTextFont(fTxt);
  legNS -> SetTextColor(fColTxt);
  legNS -> SetTextAlign(fAlign);
  legNS -> AddEntry(gNewYieldNS, sNewYields.Data());
  legNS -> AddEntry(gOldStatNS, sOldStat.Data());
  legNS -> AddEntry(gOldSysNS, sOldSys.Data());

  TLegend *legRatio = new TLegend(xyText[0], xyText[1], xyText[2], xyText[3]);
  legRatio -> SetFillColor(fLeg);
  legRatio -> SetFillStyle(fLeg);
  legRatio -> SetLineColor(fLeg);
  legRatio -> SetLineStyle(fLeg);
  legRatio -> SetTextFont(fTxt);
  legRatio -> SetTextColor(fColTxt);
  legRatio -> SetTextAlign(fAlign);
  legRatio -> AddEntry(gRatioStatNS, sRatStatNS.Data());
  legRatio -> AddEntry(gRatioSysNS, sRatSysNS.Data());
  legRatio -> AddEntry(gRatioStatAS, sRatStatAS.Data());
  legRatio -> AddEntry(gRatioSysAS, sRatSysAS.Data());
  cout << "    Made labels." << endl;


  // create line
  TLine *lRatio = new TLine(xRange[0], lineVal, xRange[1], lineVal);
  lRatio -> SetLineColor(fLinCol);
  lRatio -> SetLineStyle(fLinSty);
  lRatio -> SetLineWidth(fLinSiz);
  cout << "    Made line." << endl;


  // create yield frames
  TH2D *hYieldsAS = new TH2D(sYieldAS.Data(), sTitleAS.Data(), nXplot, xPlot[0], xPlot[1], nYplot, yPlot[0], yPlot[1]);
  TH2D *hYieldsNS = new TH2D(sYieldNS.Data(), sTitleNS.Data(), nXplot, xPlot[0], xPlot[1], nYplot, yPlot[0], yPlot[1]);
  hYieldsAS -> SetTitleFont(fTxt);
  hYieldsAS -> GetXaxis() -> SetTitle(sTitleXas.Data());
  hYieldsAS -> GetXaxis() -> SetTitleFont(fTxt);
  hYieldsAS -> GetXaxis() -> SetTitleSize(fTitleSize);
  hYieldsAS -> GetXaxis() -> SetTitleOffset(fTitleOffX);
  hYieldsAS -> GetXaxis() -> SetLabelSize(fLabSize);
  hYieldsAS -> GetXaxis() -> SetLabelFont(fTxt);
  hYieldsAS -> GetXaxis() -> CenterTitle(fCnt);
  hYieldsAS -> GetYaxis() -> SetTitle(sTitleYas.Data());
  hYieldsAS -> GetYaxis() -> SetTitleFont(fTxt);
  hYieldsAS -> GetYaxis() -> SetTitleSize(fTitleSize);
  hYieldsAS -> GetYaxis() -> SetTitleOffset(fTitleOffY);
  hYieldsAS -> GetYaxis() -> SetLabelSize(fLabSize);
  hYieldsAS -> GetYaxis() -> SetLabelFont(fTxt);
  hYieldsAS -> GetYaxis() -> CenterTitle(fCnt);
  hYieldsNS -> SetTitleFont(fTxt);
  hYieldsNS -> GetXaxis() -> SetTitle(sTitleXns.Data());
  hYieldsNS -> GetXaxis() -> SetTitleFont(fTxt);
  hYieldsNS -> GetXaxis() -> SetTitleSize(fTitleSize);
  hYieldsNS -> GetXaxis() -> SetTitleOffset(fTitleOffX);
  hYieldsNS -> GetXaxis() -> SetLabelSize(fLabSize);
  hYieldsNS -> GetXaxis() -> SetLabelFont(fTxt);
  hYieldsNS -> GetXaxis() -> CenterTitle(fCnt);
  hYieldsNS -> GetYaxis() -> SetTitle(sTitleYns.Data());
  hYieldsNS -> GetYaxis() -> SetTitleFont(fTxt);
  hYieldsNS -> GetYaxis() -> SetTitleSize(fTitleSize);
  hYieldsNS -> GetYaxis() -> SetTitleOffset(fTitleOffY);
  hYieldsNS -> GetYaxis() -> SetLabelSize(fLabSize);
  hYieldsNS -> GetYaxis() -> SetLabelFont(fTxt);
  hYieldsNS -> GetYaxis() -> CenterTitle(fCnt);

  // create ratio frames
  TH2D *hRatiosAS = new TH2D(sRatioAS.Data(), sTitleAS.Data(), nXplot, xPlot[0], xPlot[1], nRplot, rPlot[0], rPlot[1]);
  TH2D *hRatiosNS = new TH2D(sRatioNS.Data(), sTitleNS.Data(), nXplot, xPlot[0], xPlot[1], nRplot, rPlot[0], rPlot[1]);
  hRatiosAS -> SetTitleFont(fTxt);
  hRatiosAS -> GetXaxis() -> SetTitle(sTitleXas.Data());
  hRatiosAS -> GetXaxis() -> SetTitleFont(fTxt);
  hRatiosAS -> GetXaxis() -> SetTitleSize(fTitleSize);
  hRatiosAS -> GetXaxis() -> SetTitleOffset(fTitleOffX);
  hRatiosAS -> GetXaxis() -> SetLabelSize(fLabSize);
  hRatiosAS -> GetXaxis() -> SetLabelFont(fTxt);
  hRatiosAS -> GetXaxis() -> CenterTitle(fCnt);
  hRatiosAS -> GetYaxis() -> SetTitle(sTitleRas.Data());
  hRatiosAS -> GetYaxis() -> SetTitleFont(fTxt);
  hRatiosAS -> GetYaxis() -> SetTitleSize(fTitleSize);
  hRatiosAS -> GetYaxis() -> SetTitleOffset(fTitleOffY);
  hRatiosAS -> GetYaxis() -> SetLabelSize(fLabSize);
  hRatiosAS -> GetYaxis() -> SetLabelFont(fTxt);
  hRatiosAS -> GetYaxis() -> CenterTitle(fCnt);
  hRatiosNS -> SetTitleFont(fTxt);
  hRatiosNS -> GetXaxis() -> SetTitle(sTitleXns.Data());
  hRatiosNS -> GetXaxis() -> SetTitleFont(fTxt);
  hRatiosNS -> GetXaxis() -> SetTitleSize(fTitleSize);
  hRatiosNS -> GetXaxis() -> SetTitleOffset(fTitleOffX);
  hRatiosNS -> GetXaxis() -> SetLabelSize(fLabSize);
  hRatiosNS -> GetXaxis() -> SetLabelFont(fTxt);
  hRatiosNS -> GetXaxis() -> CenterTitle(fCnt);
  hRatiosNS -> GetYaxis() -> SetTitle(sTitleRns.Data());
  hRatiosNS -> GetYaxis() -> SetTitleFont(fTxt);
  hRatiosNS -> GetYaxis() -> SetTitleSize(fTitleSize);
  hRatiosNS -> GetYaxis() -> SetTitleOffset(fTitleOffY);
  hRatiosNS -> GetYaxis() -> SetLabelSize(fLabSize);
  hRatiosNS -> GetYaxis() -> SetLabelFont(fTxt);
  hRatiosNS -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Created frames." << endl;


  // make yield plot
  TCanvas *cYields  = new TCanvas(sYields.Data(), "", width, height);
  TPad    *pYieldNS = new TPad(sPadYieldNS.Data(), "", xyPadNS[0], xyPadNS[1], xyPadNS[2], xyPadNS[3]);
  TPad    *pYieldAS = new TPad(sPadYieldAS.Data(), "", xyPadAS[0], xyPadAS[1], xyPadAS[2], xyPadAS[3]);
  pYieldNS    -> SetTicks(fTick, fTick);
  pYieldNS    -> SetGrid(fGrid, fGrid);
  pYieldNS    -> SetLogx(fLogX);
  pYieldNS    -> SetLogy(fLogY);
  pYieldNS    -> SetRightMargin(fMargin);
  pYieldAS    -> SetTicks(fTick, fTick);
  pYieldAS    -> SetGrid(fGrid, fGrid);
  pYieldAS    -> SetLogx(fLogX);
  pYieldAS    -> SetLogy(fLogY);
  pYieldAS    -> SetRightMargin(fMargin);
  cYields     -> cd();
  pYieldNS    -> Draw();
  pYieldAS    -> Draw();
  pYieldNS    -> cd();
  hYieldsNS   -> Draw();
  gOldSysNS   -> Draw("P2 same");
  gOldStatNS  -> Draw("P2 same");
  gNewYieldNS -> Draw("P2 same");
  hYieldsNS   -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hYieldsNS   -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
  legNS       -> Draw();
  pt          -> Draw();
  pYieldAS    -> cd();
  hYieldsAS   -> Draw();
  gOldSysAS   -> Draw("P2 same");
  gOldStatAS  -> Draw("P2 same");
  gNewYieldAS -> Draw("P2 same");
  hYieldsAS   -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hYieldsAS   -> GetYaxis() -> SetRangeUser(yRange[0], yRange[1]);
  legAS       -> Draw();

  // make ratio plot
  TCanvas *cRatios  = new TCanvas(sRatios.Data(), "", width, height);
  TPad    *pRatioNS = new TPad(sPadRatioNS.Data(), "", xyPadNS[0], xyPadNS[1], xyPadNS[2], xyPadNS[3]);
  TPad    *pRatioAS = new TPad(sPadRatioAS.Data(), "", xyPadAS[0], xyPadAS[1], xyPadAS[2], xyPadAS[3]);
  pRatioNS     -> SetTicks(fTick, fTick);
  pRatioNS     -> SetGrid(fGrid, fGrid);
  pRatioNS     -> SetRightMargin(fMargin);
  pRatioAS     -> SetTicks(fTick, fTick);
  pRatioAS     -> SetGrid(fGrid, fGrid);
  pRatioAS     -> SetRightMargin(fMargin);
  cRatios      -> cd();
  pRatioNS     -> Draw();
  pRatioAS     -> Draw();
  pRatioNS     -> cd();
  hRatiosNS    -> Draw();
  gRatioSysNS  -> Draw("P2 same");
  gRatioStatNS -> Draw("P2 same");
  hRatiosNS    -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hRatiosNS    -> GetYaxis() -> SetRangeUser(rRange[0], rRange[1]);
  legRatio     -> Draw();
  lRatio       -> Draw();
  pRatioAS     -> cd();
  hRatiosAS    -> Draw();
  gRatioSysAS  -> Draw("P2 same");
  gRatioStatAS -> Draw("P2 same");
  hRatiosAS    -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hRatiosAS    -> GetYaxis() -> SetRangeUser(rRange[0], rRange[1]);
  pt           -> Draw();
  lRatio       -> Draw();
  cout << "    Made plots." << endl;


  // save graphs and close file
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  fOutput      -> cd();
  cYields      -> Write();
  cYields      -> Close();
  cRatios      -> Write();
  cRatios      -> Close();
  gNewYieldAS  -> Write();
  gNewYieldNS  -> Write();
  gOldStatAS   -> Write();
  gOldStatNS   -> Write();
  gOldSysAS    -> Write();
  gOldSysNS    -> Write();
  gRatioStatAS -> Write();
  gRatioStatNS -> Write();
  gRatioSysAS  -> Write();
  gRatioSysNS  -> Write();
  hYieldsAS    -> Write();
  hYieldsNS    -> Write();
  hRatiosAS    -> Write();
  hRatiosNS    -> Write();
  fOutput      -> Close();
  cout << "  Finished comparing yields!\n" << endl;

}

// End ------------------------------------------------------------------------
