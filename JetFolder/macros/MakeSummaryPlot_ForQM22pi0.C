// 'MakeSummaryPlot_ForQM22pi0.C'
// Derek Anderson
// 03.31.2022
//
// Creates a plot summarizing
// the fully corrected recoil
// jet distributions.
//
// NOTE: use the 'JetRes' variable
//       to change how certain spectra
//       are handled.
//         JetRes = 0 -- R = 0.2
//         JetRes = 1 -- R = 0.5

#include <iostream>
#include "TLine.h"
#include "TError.h"
#include "TGraph.h"
#include "TArrow.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"

// global constants
static const UInt_t NVtx(4);
static const UInt_t NErr(4);
static const UInt_t NPlot(2);
static const UInt_t NPtMax(40);
static const UInt_t NJetRes(2);
static const UInt_t NTrgTot(2);

// jet parameters
static const UInt_t JetRes(1);



void MakeSummaryPlot_ForQM22pi0() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distributions..." << endl;

  // io parameters
  const TString sOut("pionSummaryPP.onlyPy8_forQM22_betterRes_pTbinHuge.et920r05pi0.d2m4y2022.root");
  const TString sGraphStat("gIter1");
  const TString sGraphUnf("gIter1");
  const TString sGraphPyF("gIter1");
  const TString sGraphPyH("gIter1");
  const TString sInStat[NTrgTot] = {"et911r05pi0_forQM22/statErr.shiftedDataPoints_forQM22_pTbinHuge.et911r05pi0.d27m3y2022.root",
                                    "et1115r05pi0_forQM22/statErr.shiftedDataPoints_forQM22_pTbinHuge.et1115r05pi0.d27m3y2022.root"};
  const TString sInUnf[NTrgTot]  = {"et911r05pi0_forQM22/sysErr.shiftedDataPoints_forQM22_pTbinHuge.et911r05pi0.d27m3y2022.root",
                                    "et1115r05pi0_forQM22/sysErr.shiftedDataPoints_forQM22_pTbinHuge.et1115r05pi0.d27m3y2022.root"};
  const TString sInPyF[NTrgTot]  = {"et911r05pi0_forQM22/pythia.shiftedDataPoints_forQM22_pTbinHuge.et6100x911r05pi0.d28m3y2022.root",
                                    "et1115r05pi0_forQM22/pythia.shiftedDataPoints_forQM22_pTbinHuge.et6100x1115r05pi0.d28m3y2022.root"};
  const TString sInPyH[NTrgTot]  = {"et911r05pi0_forQM22/pythia.shiftedDataPoints_forQM22_pTbinHuge.et6100x911r05pi0.d28m3y2022.root",
                                    "et1115r05pi0_forQM22/pythia.shiftedDataPoints_forQM22_pTbinHuge.et6100x1115r05pi0.d28m3y2022.root"};

  // text parameters
  const TString sSys("p+p #sqrt{s} = 200 GeV");
  const TString sTyp("#pi^{0}+jet, anti-k_{T}");
  const TString sEtc("R=0.5");

  // plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T,jet}^{reco,ch} [GeV/c]");
  const TString  sTitleY("(1/N^{trg}) d^{2}N_{jets}/(dp_{T,jet}^{reco,ch} d#eta_{jet}) [GeV/c]^{-1}");
  const Double_t xPlotRange[NPlot] = {-1., 25.};
  const Double_t yPlotRange[NPlot] = {0.0000001, 3.};

  // graph parameters
  const Float_t xMinJet(4.);
  const TString sNameStat[NTrgTot] = {"gStat911pi",     "gStat1115pi"};
  const TString sNameUnf[NTrgTot]  = {"gSys911pi",      "gSys1115pi"};
  const TString sPlotStat[NTrgTot] = {"gPlotStat911pi", "gPlotStat1115pi"};
  const TString sPlotUnf[NTrgTot]  = {"gPlotSys911pi",  "gPlotSys1115pi"};
  const TString sLabels[NTrgTot]   = {"9 < E_{T}^{trig} < 11 GeV",
                                      "11 < E_{T}^{trig} < 15 GeV"};

  // trigger and jet parameters
  const UInt_t  fColSt[NJetRes][NTrgTot]      = {{923, 634}, {923, 634}};
  const UInt_t  fColUn[NJetRes][NTrgTot]      = {{920, 623}, {920, 623}};
  const UInt_t  fColOu[NJetRes][NTrgTot]      = {{923, 634}, {923, 634}};
  const UInt_t  fMarSt[NJetRes][NTrgTot]      = {{1, 1}, {1, 1}};
  const UInt_t  fMarUn[NJetRes][NTrgTot]      = {{1, 1}, {1, 1}};
  const UInt_t  fMarOu[NJetRes][NTrgTot]      = {{1, 1}, {1, 1}};
  const UInt_t  nJetPtBins[NJetRes][NTrgTot]  = {{7, 7}, {7, 7}};
  const Bool_t  goesToZero[NJetRes][NTrgTot]  = {{false, false}, {false, false}};
  const Float_t xMaxJetBins[NJetRes][NTrgTot] = {{30., 30.}, {30., 30.}};
  const Float_t fTrgScale[NJetRes][NTrgTot]   = {{1., 1.},   {1., 1.}};

  // pythia parameters
  const TString sNamePyF[NTrgTot]         = {"gPyFit911pi",      "gPyFit1115pi"};
  const TString sNamePyH[NTrgTot]         = {"gPyHist911pi",     "gPyHist1115pi"};
  const TString sPlotPyF[NTrgTot]         = {"gPlotPyFit911pi",  "gPlotPyFit1115pi"};
  const TString sPlotPyH[NTrgTot]         = {"gPlotPyHist911pi", "gPlotPyHist1115pi"};
  const UInt_t  fColPyF[NJetRes][NTrgTot] = {{1, 899}, {1, 899}};
  const UInt_t  fColPyH[NJetRes][NTrgTot] = {{1, 899}, {1, 899}};
  const UInt_t  fMarPyF[NJetRes][NTrgTot] = {{1, 1}, {1, 1}};
  const UInt_t  fMarPyH[NJetRes][NTrgTot] = {{1, 1}, {1, 1}};

  // style parameters
  const UInt_t  fLinStat(1);
  const UInt_t  fLinUnf(1);
  const UInt_t  fLinDet(1);
  const UInt_t  fLinOut(1);
  const UInt_t  fLinArr(1);
  const UInt_t  fLinPyF(9);
  const UInt_t  fLinPyH(9);
  const UInt_t  fWidStat(2);
  const UInt_t  fWidUnf(2);
  const UInt_t  fWidDet(2);
  const UInt_t  fWidOut(2);
  const UInt_t  fWidArr(2);
  const UInt_t  fWidPyF(2);
  const UInt_t  fWidPyH(2);
  const UInt_t  fFilStat(1001);
  const UInt_t  fFilUnf(1001);
  const UInt_t  fFilDet(1001);
  const UInt_t  fFilOut(0);
  const UInt_t  fFilPy(0);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  const Float_t fSizArr(0.05);
  const Float_t fSizMar(1.25);
  const Float_t fSizHat(1.5);

  // select trigger and jets
  UInt_t  fColStat[NTrgTot];
  UInt_t  fColUnf[NTrgTot];
  UInt_t  fColPyWF[NTrgTot];
  UInt_t  fColPyWH[NTrgTot];
  UInt_t  fColOut[NTrgTot];
  UInt_t  fMarStat[NTrgTot];
  UInt_t  fMarUnf[NTrgTot];
  UInt_t  fMarPyWF[NTrgTot];
  UInt_t  fMarPyWH[NTrgTot];
  UInt_t  fMarOut[NTrgTot];
  UInt_t  nJetPoints[NTrgTot];
  Bool_t  endIsArrow[NTrgTot];
  Float_t xMaxJet[NTrgTot];
  Float_t fScalers[NTrgTot];
  TString sTrg("");
  TString sJet("");
  switch (JetRes) {
    case 0:
      sTrg = "#pi^{0}/#gamma_{dir} trig.";
      sJet = "anti-k_{T}, R = 0.2";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColStat[iTrg]   = fColSt[0][iTrg];
        fColUnf[iTrg]    = fColUn[0][iTrg];
        fColPyWF[iTrg]   = fColPyF[0][iTrg];
        fColPyWH[iTrg]   = fColPyH[0][iTrg];
        fColOut[iTrg]    = fColOu[0][iTrg];
        fMarStat[iTrg]   = fMarSt[0][iTrg];
        fMarUnf[iTrg]    = fMarUn[0][iTrg];
        fMarPyWF[iTrg]   = fMarPyF[0][iTrg];
        fMarPyWH[iTrg]   = fMarPyH[0][iTrg];
        fMarOut[iTrg]    = fMarOu[0][iTrg];
        nJetPoints[iTrg] = nJetPtBins[0][iTrg];
        endIsArrow[iTrg] = goesToZero[0][iTrg];
        xMaxJet[iTrg]    = xMaxJetBins[0][iTrg];
        fScalers[iTrg]   = fTrgScale[0][iTrg];
      }
      break;
    case 1:
      sTrg = "#pi^{0}/#gamma_{dir} trig.";
      sJet = "anti-k_{T}, R = 0.5";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColStat[iTrg]   = fColSt[1][iTrg];
        fColUnf[iTrg]    = fColUn[1][iTrg];
        fColPyWF[iTrg]   = fColPyF[1][iTrg];
        fColPyWH[iTrg]   = fColPyH[1][iTrg];
        fColOut[iTrg]    = fColOu[1][iTrg];
        fMarStat[iTrg]   = fMarSt[1][iTrg];
        fMarUnf[iTrg]    = fMarUn[1][iTrg];
        fMarPyWF[iTrg]   = fMarPyF[1][iTrg];
        fMarPyWH[iTrg]   = fMarPyH[1][iTrg];
        fMarOut[iTrg]    = fMarOu[1][iTrg];
        nJetPoints[iTrg] = nJetPtBins[1][iTrg];
        endIsArrow[iTrg] = goesToZero[1][iTrg];
        xMaxJet[iTrg]    = xMaxJetBins[1][iTrg];
        fScalers[iTrg]   = fTrgScale[1][iTrg];
      }
      break;
    default:
      sTrg = "#pi^{0}/#gamma_{dir} trig.";
      sJet = "anti-k_{T}, R = 0.2";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColStat[iTrg]   = fColSt[0][iTrg];
        fColUnf[iTrg]    = fColUn[0][iTrg];
        fColPyWF[iTrg]   = fColPyF[0][iTrg];
        fColPyWH[iTrg]   = fColPyH[0][iTrg];
        fColOut[iTrg]    = fColOu[0][iTrg];
        fMarStat[iTrg]   = fMarSt[0][iTrg];
        fMarUnf[iTrg]    = fMarUn[0][iTrg];
        fMarPyWF[iTrg]   = fMarPyF[0][iTrg];
        fMarPyWH[iTrg]   = fMarPyH[0][iTrg];
        fMarOut[iTrg]    = fMarOu[0][iTrg];
        nJetPoints[iTrg] = nJetPtBins[0][iTrg];
        endIsArrow[iTrg] = goesToZero[0][iTrg];
        xMaxJet[iTrg]    = xMaxJetBins[0][iTrg];
        fScalers[iTrg]   = fTrgScale[0][iTrg];
      }
      break;
  }  // end switch-case

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }

  TFile *fInStat[NTrgTot];
  TFile *fInUnf[NTrgTot];
  TFile *fInPyF[NTrgTot];
  TFile *fInPyH[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    fInStat[iTrg] = new TFile(sInStat[iTrg].Data(), "read");
    fInUnf[iTrg]  = new TFile(sInUnf[iTrg].Data(), "read");
    fInPyF[iTrg]  = new TFile(sInPyF[iTrg].Data(), "read");
    fInPyH[iTrg]  = new TFile(sInPyH[iTrg].Data(), "read");
    if (!fInStat[iTrg] || !fInUnf[iTrg] || !fInPyF[iTrg] || !fInPyH[iTrg]) {
      cerr << "PANIC: couldn't open an input file!\n"
           << "       fInStat[" << iTrg << "] = " << fInStat[iTrg] << ", fInUnf[" << iTrg << "] = " << fInUnf[iTrg] << ", fInPyF[" << iTrg << "] = " << fInPyF[iTrg] << ", fInPyH[" << iTrg << "] = " << fInPyH[iTrg] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input graphs
  TGraphAsymmErrors *gStat[NTrgTot];
  TGraphAsymmErrors *gUnf[NTrgTot];
  TGraphAsymmErrors *gPyF[NTrgTot];
  TGraphAsymmErrors *gPyH[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gStat[iTrg] = (TGraphAsymmErrors*) fInStat[iTrg] -> Get(sGraphStat.Data());
    gUnf[iTrg]  = (TGraphAsymmErrors*) fInUnf[iTrg]  -> Get(sGraphUnf.Data());
    gPyF[iTrg]  = (TGraphAsymmErrors*) fInPyF[iTrg]  -> Get(sGraphPyF.Data());
    gPyH[iTrg]  = (TGraphAsymmErrors*) fInPyH[iTrg]  -> Get(sGraphPyH.Data());
    if (!gStat[iTrg] || !gUnf[iTrg] || !gPyF[iTrg] || !gPyH[iTrg]) {
      cerr << "PANIC: couldn't open an input graph!\n"
           << "       gStat[" << iTrg << "] = " << gStat[iTrg] << ", gUnf[" << iTrg << "] = " << gUnf[iTrg] << ", gPyF[" << iTrg << "] = " << gPyF[iTrg] << ", gPyH[" << iTrg << "] = " << gPyH[iTrg] << "\n"
           << endl;
      return;
    }
    gStat[iTrg] -> SetName(sNameStat[iTrg].Data());
    gUnf[iTrg]  -> SetName(sNameUnf[iTrg].Data());
    gPyF[iTrg]  -> SetName(sNamePyF[iTrg].Data());
    gPyH[iTrg]  -> SetName(sNamePyH[iTrg].Data());
  }
  cout << "    Grabbed graphs." << endl;


  // create graphs for plotting
  Double_t xPlot[NTrgTot][NPtMax];
  Double_t yPlot[NTrgTot][NPtMax];
  Double_t yPlotLoU[NTrgTot][NPtMax];
  Double_t yPlotHiU[NTrgTot][NPtMax];
  Double_t xPlotPyF[NTrgTot][NPtMax];
  Double_t xPlotPyH[NTrgTot][NPtMax];
  Double_t yPlotPyF[NTrgTot][NPtMax];
  Double_t yPlotPyH[NTrgTot][NPtMax];
  Double_t ePlotZe[NTrgTot][NErr][NPtMax];
  Double_t ePlotSt[NTrgTot][NErr][NPtMax];
  Double_t ePlotUn[NTrgTot][NErr][NPtMax];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    for (UInt_t iPt = 0; iPt < NPtMax; iPt++) {
      xPlot[iTrg][iPt]    = 0.;
      yPlot[iTrg][iPt]    = 0.;
      yPlotLoU[iTrg][iPt] = 0.;
      yPlotHiU[iTrg][iPt] = 0.;
      xPlotPyF[iTrg][iPt] = 0.;
      yPlotPyH[iTrg][iPt] = 0.;
      for (UInt_t iErr = 0; iErr < NErr; iErr++) {
        ePlotZe[iTrg][iErr][iPt] = 0.;
        ePlotSt[iTrg][iErr][iPt] = 0.;
        ePlotUn[iTrg][iErr][iPt] = 0.;
      }
    }  // end bin loop
  }  // end trigger loop

  TLine             *lEndPoint[NTrgTot];
  TArrow            *aEndArrow[NTrgTot];
  TGraph            *gPlotLoU[NTrgTot];
  TGraph            *gPlotHiU[NTrgTot];
  TGraph            *gPlotPyF[NTrgTot];
  TGraph            *gPlotPyH[NTrgTot];
  TGraphAsymmErrors *gPlotSt[NTrgTot];
  TGraphAsymmErrors *gPlotUn[NTrgTot];
  cout << "    Reading in data points..." << endl;

  UInt_t iJetPoint(0);
  UInt_t nJetCheck(0);
  UInt_t nPoints(0);
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {

    // loop over points
    iJetPoint = 0;
    nJetCheck = 0;
    nPoints   = gStat[iTrg] -> GetN();
    for (UInt_t iPoint = 0; iPoint < nPoints; iPoint++) {

      // get input points
      Double_t xStat(0.);
      Double_t xUnf(0.);
      Double_t xDet(0.);
      Double_t yStat(0.);
      Double_t yUnf(0.);
      Double_t yDet(0.);
      gStat[iTrg] -> GetPoint(iPoint, xStat, yStat);
      gUnf[iTrg]  -> GetPoint(iPoint, xUnf, yUnf);

      // check if point is in range
      const Bool_t isInRange = ((xStat > xMinJet) && (xStat < xMaxJet[iTrg]));
      if (!isInRange) {
        continue;
      } else {
        nJetCheck++;
      }

      // check if stat and sys yields match
      const Bool_t unfYieldsMatch = (yUnf == yStat);
      if (!unfYieldsMatch) {
        cerr << "      WARNING: stat. and sys. yields don't match at bin #" << iPoint << " for trigger bin #" << iTrg << "!\n"
             << "               yStat = " << yStat << ", yUnf = " << yUnf
             << endl;
      }
      xPlot[iTrg][iJetPoint] = xStat;
      yPlot[iTrg][iJetPoint] = yStat * fScalers[iTrg];

      // get errors
      ePlotSt[iTrg][0][iJetPoint] = gStat[iTrg] -> GetErrorXlow(iPoint);
      ePlotUn[iTrg][0][iJetPoint] = gUnf[iTrg]  -> GetErrorXlow(iPoint);
      ePlotSt[iTrg][1][iJetPoint] = gStat[iTrg] -> GetErrorXhigh(iPoint);
      ePlotUn[iTrg][1][iJetPoint] = gUnf[iTrg]  -> GetErrorXhigh(iPoint);
      ePlotSt[iTrg][2][iJetPoint] = fScalers[iTrg] * (gStat[iTrg] -> GetErrorYlow(iPoint));
      ePlotUn[iTrg][2][iJetPoint] = fScalers[iTrg] * (gUnf[iTrg]  -> GetErrorYlow(iPoint));
      ePlotSt[iTrg][3][iJetPoint] = fScalers[iTrg] * (gStat[iTrg] -> GetErrorYhigh(iPoint));
      ePlotUn[iTrg][3][iJetPoint] = fScalers[iTrg] * (gUnf[iTrg]  -> GetErrorYhigh(iPoint));

      // get upper and lower bands
      yPlotLoU[iTrg][iJetPoint] = yPlot[iTrg][iJetPoint] - ePlotUn[iTrg][2][iJetPoint];
      yPlotHiU[iTrg][iJetPoint] = yPlot[iTrg][iJetPoint] + ePlotUn[iTrg][3][iJetPoint];

      // get pythia points
      gPyF[iTrg] -> GetPoint(iPoint, xPlotPyF[iTrg][iJetPoint], yPlotPyF[iTrg][iJetPoint]);
      gPyH[iTrg] -> GetPoint(iPoint, xPlotPyH[iTrg][iJetPoint], yPlotPyH[iTrg][iJetPoint]);
      yPlotPyF[iTrg][iJetPoint] = yPlotPyF[iTrg][iJetPoint] * fScalers[iTrg];
      yPlotPyH[iTrg][iJetPoint] = yPlotPyH[iTrg][iJetPoint] * fScalers[iTrg];

      // check if arrow needs to be created
      const Bool_t isLastPoint = ((iJetPoint + 1) == nJetPoints[iTrg]);
      if (endIsArrow[iTrg] && isLastPoint) {

        // get start and stop points
        const Double_t xStart = xPlot[iTrg][iJetPoint] - ePlotSt[iTrg][0][iJetPoint];
        const Double_t xStop  = xPlot[iTrg][iJetPoint] + ePlotSt[iTrg][1][iJetPoint];
        const Double_t yLoSt  = yPlot[iTrg][iJetPoint] - ePlotSt[iTrg][2][iJetPoint];
        const Double_t yLoUn  = yPlot[iTrg][iJetPoint] - ePlotUn[iTrg][2][iJetPoint];
        const Double_t yHiSt  = yPlot[iTrg][iJetPoint] + ePlotSt[iTrg][3][iJetPoint];
        const Double_t yHiUn  = yPlot[iTrg][iJetPoint] + ePlotUn[iTrg][3][iJetPoint];
        const Double_t yLoSy  = yLoUn;
        const Double_t yHiSy  = yHiUn;
        const Double_t yStart = TMath::Max(yHiSt, yHiSy);
        const Double_t yStop  = TMath::Max(yLoSt, yLoSy);

        // make line
        lEndPoint[iTrg] = new TLine(xStart, yStart, xStop, yStart);
        lEndPoint[iTrg] -> SetLineColor(fColStat[iTrg]);
        lEndPoint[iTrg] -> SetLineWidth(fWidArr);
        lEndPoint[iTrg] -> SetLineStyle(fLinArr);

        // make arrow
        aEndArrow[iTrg] = new TArrow(xPlot[iTrg][iJetPoint], yStart, xPlot[iTrg][iJetPoint], yStop, fSizArr, ">");
        aEndArrow[iTrg] -> SetLineColor(fColStat[iTrg]);
        aEndArrow[iTrg] -> SetLineWidth(fWidArr);
        aEndArrow[iTrg] -> SetLineStyle(fLinArr);
      }
      iJetPoint++;
    }  // end point loop

    // create graphs
    const Int_t nGraph = (Int_t) nJetPoints[iTrg];
    gPlotSt[iTrg]  = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotZe[iTrg][0], ePlotZe[iTrg][1], ePlotSt[iTrg][2], ePlotSt[iTrg][3]);
    gPlotUn[iTrg]  = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotZe[iTrg][0], ePlotZe[iTrg][1], ePlotUn[iTrg][2], ePlotUn[iTrg][3]);
    gPlotLoU[iTrg] = new TGraph(nGraph, xPlot[iTrg], yPlotLoU[iTrg]);
    gPlotHiU[iTrg] = new TGraph(nGraph, xPlot[iTrg], yPlotHiU[iTrg]);
    gPlotPyF[iTrg] = new TGraph(nGraph, xPlotPyF[iTrg], yPlotPyF[iTrg]);
    gPlotPyH[iTrg] = new TGraph(nGraph, xPlotPyH[iTrg], yPlotPyH[iTrg]);

/* will fix later
    if (!endIsArrow[iTrg]) {
      gPlotOu[iTrg] = new TGraph(nGraph, xPlot[iTrg], yPlot[iTrg]);
      gPlotSt[iTrg] = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotSt[iTrg][2], ePlotSt[iTrg][3]);
      gPlotUn[iTrg] = new TGraphAsymmErrors(nGraph, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotUn[iTrg][2], ePlotUn[iTrg][3]);
    } else {
      gPlotOu[iTrg] = new TGraphAsymmErrors(nGraph - 1, xPlot[iTrg], yPlot[iTrg]);
      gPlotSt[iTrg] = new TGraphAsymmErrors(nGraph - 1, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotSt[iTrg][2], ePlotSt[iTrg][3]);
      gPlotUn[iTrg] = new TGraphAsymmErrors(nGraph - 1, xPlot[iTrg], yPlot[iTrg], ePlotSt[iTrg][0], ePlotSt[iTrg][1], ePlotUn[iTrg][2], ePlotUn[iTrg][3]);
    }
*/

    // set data and pythia names
    gPlotSt[iTrg]  -> SetName(sPlotStat[iTrg].Data());
    gPlotUn[iTrg]  -> SetName(sPlotUnf[iTrg].Data());
    gPlotPyF[iTrg] -> SetName(sPlotPyF[iTrg].Data());
    gPlotPyH[iTrg] -> SetName(sPlotPyH[iTrg].Data());

    // create bound names
    TString sNameLoU = gPlotUn[iTrg] -> GetName();
    TString sNameHiU = gPlotUn[iTrg] -> GetName();
    sNameLoU.Append("_LoErr");
    sNameHiU.Append("_HiErr");

    // set bound names
    gPlotLoU[iTrg] -> SetName(sNameLoU.Data());
    gPlotHiU[iTrg] -> SetName(sNameHiU.Data());

    // check if no. of points is reasonable
    if (nJetCheck != nJetPoints[iTrg]) {
      cerr << "      WARNING: no. of accepted points is off in trigger bin #" << iTrg << "!" << endl;
    }
    cout << "      Created plotting graphs for bin #" << iTrg << "..." << endl;

  }  // end trigger loop
  cout << "    Made graphs for plotting." << endl;

  // set styles
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gPlotSt[iTrg]  -> SetLineColor(fColStat[iTrg]);
    gPlotSt[iTrg]  -> SetLineStyle(fLinStat);
    gPlotSt[iTrg]  -> SetLineWidth(fWidStat);
    gPlotSt[iTrg]  -> SetFillColor(fColStat[iTrg]);
    gPlotSt[iTrg]  -> SetFillStyle(fFilStat);
    gPlotSt[iTrg]  -> SetMarkerColor(fColStat[iTrg]);
    gPlotSt[iTrg]  -> SetMarkerStyle(fMarStat[iTrg]);
    gPlotSt[iTrg]  -> SetMarkerSize(fSizMar);
    gPlotSt[iTrg]  -> SetTitle(sTitle.Data());
    gPlotSt[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotSt[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotSt[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotSt[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotSt[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotSt[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    gPlotSt[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    gPlotSt[iTrg]  -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    gPlotSt[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotSt[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotSt[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotSt[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotSt[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    gPlotSt[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    gPlotUn[iTrg]  -> SetLineColor(fColUnf[iTrg]);
    gPlotUn[iTrg]  -> SetLineStyle(fLinUnf);
    gPlotUn[iTrg]  -> SetLineWidth(fWidUnf);
    gPlotUn[iTrg]  -> SetFillColor(fColUnf[iTrg]);
    gPlotUn[iTrg]  -> SetFillStyle(fFilUnf);
    gPlotUn[iTrg]  -> SetMarkerColor(fColUnf[iTrg]);
    gPlotUn[iTrg]  -> SetMarkerStyle(fMarUnf[iTrg]);
    gPlotUn[iTrg]  -> SetMarkerSize(fSizMar);
    gPlotUn[iTrg]  -> SetTitle(sTitle.Data());
    gPlotUn[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotUn[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotUn[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotUn[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotUn[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotUn[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    gPlotUn[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    gPlotUn[iTrg]  -> GetYaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotUn[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotUn[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotUn[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotUn[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotUn[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    gPlotUn[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    gPlotPyF[iTrg]  -> SetLineColor(fColPyWF[iTrg]);
    gPlotPyF[iTrg]  -> SetLineStyle(fLinPyF);
    gPlotPyF[iTrg]  -> SetLineWidth(fWidPyF);
    gPlotPyF[iTrg]  -> SetFillColor(fColPyWF[iTrg]);
    gPlotPyF[iTrg]  -> SetFillStyle(fFilPy);
    gPlotPyF[iTrg]  -> SetMarkerColor(fColPyWF[iTrg]);
    gPlotPyF[iTrg]  -> SetMarkerStyle(fMarPyWF[iTrg]);
    gPlotPyF[iTrg]  -> SetMarkerSize(fSizMar);
    gPlotPyF[iTrg]  -> SetTitle(sTitle.Data());
    gPlotPyF[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotPyF[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotPyF[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotPyF[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotPyF[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotPyF[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    gPlotPyF[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    gPlotPyF[iTrg]  -> GetYaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotPyF[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotPyF[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotPyF[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotPyF[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotPyF[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    gPlotPyF[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    gPlotPyH[iTrg]  -> SetLineColor(fColPyWH[iTrg]);
    gPlotPyH[iTrg]  -> SetLineStyle(fLinPyH);
    gPlotPyH[iTrg]  -> SetLineWidth(fWidPyH);
    gPlotPyH[iTrg]  -> SetFillColor(fColPyWH[iTrg]);
    gPlotPyH[iTrg]  -> SetFillStyle(fFilPy);
    gPlotPyH[iTrg]  -> SetMarkerColor(fColPyWH[iTrg]);
    gPlotPyH[iTrg]  -> SetMarkerStyle(fMarPyWH[iTrg]);
    gPlotPyH[iTrg]  -> SetMarkerSize(fSizMar);
    gPlotPyH[iTrg]  -> SetTitle(sTitle.Data());
    gPlotPyH[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotPyH[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotPyH[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotPyH[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotPyH[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotPyH[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    gPlotPyH[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    gPlotPyH[iTrg]  -> GetYaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotPyH[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotPyH[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotPyH[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotPyH[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotPyH[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    gPlotPyH[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    gPlotLoU[iTrg] -> SetLineColor(fColOut[iTrg]);
    gPlotLoU[iTrg] -> SetLineStyle(fLinOut);
    gPlotLoU[iTrg] -> SetLineWidth(fWidOut);
    gPlotLoU[iTrg] -> SetFillColor(fColOut[iTrg]);
    gPlotLoU[iTrg] -> SetFillStyle(fFilOut);
    gPlotLoU[iTrg] -> SetMarkerColor(fColOut[iTrg]);
    gPlotLoU[iTrg] -> SetMarkerStyle(fMarOut[iTrg]);
    gPlotLoU[iTrg] -> SetMarkerSize(fSizMar);
    gPlotLoU[iTrg] -> SetTitle(sTitle.Data());
    gPlotLoU[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotLoU[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotLoU[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotLoU[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotLoU[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotLoU[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    gPlotLoU[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    gPlotLoU[iTrg] -> GetYaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotLoU[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotLoU[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotLoU[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotLoU[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotLoU[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    gPlotLoU[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    gPlotHiU[iTrg] -> SetLineColor(fColOut[iTrg]);
    gPlotHiU[iTrg] -> SetLineStyle(fLinOut);
    gPlotHiU[iTrg] -> SetLineWidth(fWidOut);
    gPlotHiU[iTrg] -> SetFillColor(fColOut[iTrg]);
    gPlotHiU[iTrg] -> SetFillStyle(fFilOut);
    gPlotHiU[iTrg] -> SetMarkerColor(fColOut[iTrg]);
    gPlotHiU[iTrg] -> SetMarkerStyle(fMarOut[iTrg]);
    gPlotHiU[iTrg] -> SetMarkerSize(fSizMar);
    gPlotHiU[iTrg] -> SetTitle(sTitle.Data());
    gPlotHiU[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotHiU[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    gPlotHiU[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    gPlotHiU[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    gPlotHiU[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    gPlotHiU[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    gPlotHiU[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    gPlotHiU[iTrg] -> GetYaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    gPlotHiU[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    gPlotHiU[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    gPlotHiU[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    gPlotHiU[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    gPlotHiU[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    gPlotHiU[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end trigger loop

  // for legend
  const UInt_t fColPyWfLeg(923);
  const UInt_t fColPyWhLeg(923);

  TGraph *gPyWfLeg = (TGraph*) gPlotPyF[0] -> Clone();
  TGraph *gPyWhLeg = (TGraph*) gPlotPyH[0] -> Clone();
  gPyWfLeg -> SetFillColor(fColPyWfLeg);
  gPyWfLeg -> SetLineColor(fColPyWfLeg);
  gPyWhLeg -> SetFillColor(fColPyWhLeg);
  gPyWhLeg -> SetLineColor(fColPyWhLeg);

  TGraphAsymmErrors *gStatLeg[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gStatLeg[iTrg] = (TGraphAsymmErrors*) gPlotSt[iTrg] -> Clone();
    gStatLeg[iTrg] -> SetLineColor(fColUnf[iTrg]);
  }

  // for QM22
  TGraph *gDataLeg[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gDataLeg[iTrg] = (TGraphAsymmErrors*) gPlotSt[iTrg] -> Clone();
    gDataLeg[iTrg] -> SetFillColor(fColUnf[iTrg]);
    gDataLeg[iTrg] -> SetFillStyle(1001);
    gDataLeg[iTrg] -> SetLineColor(fColStat[iTrg]);
    gDataLeg[iTrg] -> SetLineWidth(5);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  fAlnLe(12);
  const UInt_t  nColLe(1);
  const Float_t hObj(0.05);
  //const Float_t hLeg(((NTrgTot / 2) + 2) * hObj);
  const Float_t hLeg((3 + NTrgTot) * hObj);
  const Float_t yLeg(0.1 + hLeg);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.7, yLeg};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnLe);
  leg -> SetNColumns(nColLe);
  leg -> AddEntry((TObject*)0, sSys.Data(), "");
  leg -> AddEntry((TObject*)0, sTyp.Data(), "");
  leg -> AddEntry((TObject*)0, sEtc.Data(), "");
  leg -> AddEntry(gPyWfLeg, "PYTHIA-8", "l");
  leg -> AddEntry(gDataLeg[0], sLabels[0].Data(), "lf");
  leg -> AddEntry(gDataLeg[1], sLabels[1].Data(), "lf");
  //leg -> AddEntry(gPyWhLeg, "Pythia8 [weighted w/ #Delta^{Py8}_{D}]", "l");
  //for (UInt_t iPlot = 0; iPlot < (NTrgTot / 2); iPlot++) {
  //  const UInt_t iTrgPi = iPlot;
  //  const UInt_t iTrgGa = iPlot + (NTrgTot / 2);
  //  leg -> AddEntry(gStatLeg[iTrgPi], sLabels[iTrgPi].Data(), "lf");
  //  leg -> AddEntry(gStatLeg[iTrgGa], sLabels[iTrgGa].Data(), "lf");
  //}
/* will fix later
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    leg -> AddEntry(gPlotUn[iTrg], sLabels[iTrg].Data(), "fl");
  }
*/
  cout << "    Made legend." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const UInt_t  fAlnTx(32);
  const UInt_t  fAlnLb(22);
  const UInt_t  fFonLb(62);
  const Float_t fTxtTx[NVtx] = {0.5, 0.5, 0.7, 0.7};
  TPaveText *txt = new TPaveText(fTxtTx[0], fTxtTx[1], fTxtTx[2], fTxtTx[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlnTx);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sTyp.Data());
  txt -> AddText(sJet.Data());

  // for QM22
  TPaveText *prelim = new TPaveText(0.7, 0.7, 0.9, 0.75, "NDC NB");
  prelim -> SetFillColor(0);
  prelim -> SetFillStyle(0);
  prelim -> SetLineColor(0);
  prelim -> SetLineStyle(0);
  prelim -> SetTextFont(42);
  prelim -> SetTextAlign(22);
  prelim -> SetTextColor(2);
  prelim -> AddText("#bf{STAR} preliminary");
  cout << "    Made text box." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.15);

  TCanvas *cSummary = new TCanvas("cCorrectedData", "", width, height);
  cSummary    -> SetGrid(fGrid, fGrid);
  cSummary    -> SetTicks(fTick, fTick);
  cSummary    -> SetLogx(fLogX);
  cSummary    -> SetLogy(fLogY);
  cSummary    -> SetBorderMode(fMode);
  cSummary    -> SetBorderSize(fBord);
  cSummary    -> SetFrameBorderMode(fFrame);
  cSummary    -> SetLeftMargin(fMarginL);
  cSummary    -> SetRightMargin(fMarginR);
  cSummary    -> SetTopMargin(fMarginT);
  cSummary    -> SetBottomMargin(fMarginB);
  cSummary    -> cd();
  gPlotUn[1]  -> Draw("AZ3");
  gPlotLoU[1] -> Draw("L");
  gPlotHiU[1] -> Draw("L");
  gPlotSt[1]  -> Draw("Z3");
  gPlotUn[0]  -> Draw("Z3");
  gPlotLoU[0] -> Draw("L");
  gPlotHiU[0] -> Draw("L");
  gPlotSt[0]  -> Draw("Z3");
  gPlotPyH[1] -> Draw("L");
  gPlotPyH[0] -> Draw("L");
/* will fix later
  for (Int_t iTrgPlot = (NTrgTot - 1); iTrgPlot > -1; iTrgPlot--) {
    if (iTrgPlot == (NTrgTot - 1)) {
      gPlotUn[iTrgPlot]  -> Draw("AZ3");
      gPlotLoU[iTrgPlot] -> Draw("L");
      gPlotHiU[iTrgPlot] -> Draw("L");
      gPlotPyF[iTrgPlot] -> Draw("L");
      gPlotPyH[iTrgPlot] -> Draw("L");
    } else {
      gPlotUn[iTrgPlot]  -> Draw("Z3");
      gPlotLoU[iTrgPlot] -> Draw("L");
      gPlotHiU[iTrgPlot] -> Draw("L");
      gPlotPyF[iTrgPlot] -> Draw("L");
      gPlotPyH[iTrgPlot] -> Draw("L");
    }
    gPlotSt[iTrgPlot] -> Draw("P");
    if (endIsArrow[iTrgPlot]) {
      lEndPoint[iTrgPlot] -> Draw();
      aEndArrow[iTrgPlot] -> Draw();
    }
  }
*/
  //txt      -> Draw();
  prelim   -> Draw();
  leg      -> Draw();
  fOut     -> cd();
  cSummary -> Write();
  cSummary -> Close();
  cout << "    Made summary plot.\n"
       << "    Preparing Quark Matter plot..."
       << endl;

  // for QM22
  const Double_t start_bin            = 1.0;
  const Double_t stop_bin             = 30.5;
  const Double_t start_bin_r5         = 2.5;
  const Double_t stop_bin_r5          = 30.5;
  const Double_t start_bin_pi0        = 1.0;
  const Double_t stop_bin_pi0         = 25.5;
  const Double_t start_bin_pi0_r5     = 3.0;
  const Double_t stop_bin_pi0_r5      = 25.5;
  const Double_t unfold_color         = kBlue + 3;
  const Double_t backfold_color       = kGreen - 3;
  const Double_t YBin                 = 1000;
  const Double_t Ymin                 = 0.0001;
  const Double_t Ymax                 = 0.2;
  const Double_t XBin                 = 100;
  const Double_t Xmin                 = 0;
  const Double_t Xmax                 = 22;
  const Double_t YBin_pi0             = 1000;
  const Double_t Ymin_pi0             = 0.00006;
  const Double_t Ymax_pi0             = 0.8;
  const Double_t XBin_pi0             = 100;
  const Double_t Xmin_pi0             = -2.0;
  const Double_t Xmax_pi0             = 26.6;
  const Int_t    PYTHIALINESTYLE      = 9;
  const Int_t    SysColor_15          = kRed - 6;  //kRed-3;
  const Int_t    StatColor_15         = kRed - 2;  //kRed+2;
  const Int_t    SysLineColor_15      = 1;  //kRed+4;
  const Int_t    SysLineColor_pyth_15 = kRed + 2;
  const Int_t    SysColor_11          = kSpring + 3;  //kRed+3;
  const Int_t    StatColor_11         = kSpring - 7;  //kRed+4;
  const Int_t    SysLineColor_11      = 1;  //kRed+4;
  const Int_t    SysLineColor_pyth_11 = kSpring-7;   //kRed+4;
  const Int_t    SysColor_9           = kGray;  //Orange-4;  //kOrange+5;
  const Int_t    StatColor_9          = 14;  //kOrange-7;
  const Int_t    SysLineColor_9       = 1;
  const Int_t    SysLineColor_pyth_9  = 1;  //kOrange-1;  //kRed+4;
  const Int_t    PythColor_pi0_9      = kBlue + 3;
  const Int_t    SysColor_pi0_9       = kBlue - 9;
  const Int_t    StatColor_pi0_9      = kBlue - 6;
  const Int_t    LineColor_pi0_9      = kBlue + 3; 
  const Int_t    PythColor_pi0_11     = kMagenta + 1;
  const Int_t    SysColor_pi0_11      = kViolet - 9;
  const Int_t    StatColor_pi0_11     = kViolet - 4; 
  const Int_t    LineColor_pi0_11     = kViolet - 5;
  const Int_t    MainColor            = kRed - 3;
  const Double_t plot_NStart          = 0;
  const Double_t plot_NEnd            = 16.7; 
  const Double_t Nstart_bin           = 1.0;
  const Double_t Nstop_bin            = 16.1;
  const Double_t upperY               = 0.025;  //8.8;  //3.7
  const Double_t lowerY               = -0.01;  //0599;
  const Double_t MarkClr510           = kYellow + 4;
  const Double_t MarkClr1015          = kOrange + 9;
  const Double_t MarkClr1520          = kAzure - 6;
  cout << "      Declared variables." << endl;

  TCanvas *c3 = new TCanvas("cQuarkMatterPlot", "", 870, 870);
  gStyle -> SetOptStat(0);
  c3     -> Range(0, 0, 1, 1);
  c3     -> SetFillColor(0);
  c3     -> SetBorderMode(0);
  c3     -> SetBorderSize(2);
  c3     -> SetLogy();
  c3     -> SetLogz();
  c3     -> SetTickx(1);
  c3     -> SetTicky(1);
  c3     -> SetGridx(0);
  c3     -> SetGridy(0);
  c3     -> SetLeftMargin(0.196692);
  c3     -> SetRightMargin(0.03933966);
  c3     -> SetTopMargin(0.03678738);
  c3     -> SetBottomMargin(0.3);
  c3     -> SetFrameLineWidth(2);
  c3     -> SetFrameBorderMode(0);
  cout << "      Created canvas." << endl;

  TPad *c3_1 = new TPad("cQuarkMatterPlot_1", "_1", 0, 0.08174387, 0.9153605, 0.9659401);
  c3_1 -> Draw();
  c3_1 -> cd();
  c3_1 -> Range(1.200484, -6.512013, 3.205194, 0.7091003);
  c3_1 -> SetFillColor(0);
  c3_1 -> SetBorderMode(0);
  c3_1 -> SetBorderSize(2);
  c3_1 -> SetLogy();
  c3_1 -> SetLogz();
  c3_1 -> SetTickx(1);
  c3_1 -> SetTicky(1);
  c3_1 -> SetGridx(0);
  c3_1 -> SetGridy(0);
  c3_1 -> SetLeftMargin(0.1943005);
  c3_1 -> SetRightMargin(0.002590665);
  c3_1 -> SetTopMargin(0);
  c3_1 -> SetBottomMargin(0.1125925);
  c3_1 -> SetFrameLineWidth(2);
  c3_1 -> SetFrameBorderMode(0);
  c3_1 -> SetFrameLineWidth(2);
  c3_1 -> SetFrameBorderMode(0);
  cout << "      Created pad." << endl;
  
  TH2D *h1 = new TH2D("h1","", XBin, Xmin, Xmax, YBin, Ymin, Ymax);
  h1 -> SetDirectory(0);
  h1 -> SetStats(0);
  h1 -> GetXaxis() -> SetTitle("p_{T,jet}^{ch} [GeV/c]");
  h1 -> GetXaxis() -> CenterTitle(true);
  h1 -> GetXaxis() -> SetNdivisions(508);
  h1 -> GetXaxis() -> SetLabelFont(42);
  h1 -> GetXaxis() -> SetLabelSize(0.05);
  h1 -> GetXaxis() -> SetTitleOffset(1.04);
  h1 -> GetXaxis() -> SetTitleSize(0.05);
  h1 -> GetXaxis() -> SetTitleFont(42);
  h1 -> GetYaxis() -> SetTitle("#frac{1}{N_{trig}} #frac{dN^{2}_{jet}}{dp_{T,jet}^{ch} d#eta_{jet}} [GeV/c]^{-1}");
  h1 -> GetYaxis() -> CenterTitle(true);
  h1 -> GetYaxis() -> SetNdivisions(504);
  h1 -> GetYaxis() -> SetLabelFont(42);
  h1 -> GetYaxis() -> SetLabelSize(0.043);
  h1 -> GetYaxis() -> SetTitleSize(0.035);
  h1 -> GetYaxis() -> SetTitleOffset(2.29);
  h1 -> GetYaxis() -> SetTitleFont(42);
  h1 -> GetZaxis() -> SetLabelFont(42);
  h1 -> GetZaxis() -> SetLabelSize(0.03);
  h1 -> GetZaxis() -> SetTitleSize(0.35);
  h1 -> GetZaxis() -> SetTitleOffset(1);
  h1 -> GetZaxis() -> SetTitleFont(42);
  h1 -> Draw("9");
  cout << "      Created histogram frame." << endl;

  gPlotUn[1]  -> Draw("Z3");
  gPlotLoU[1] -> Draw("L");
  gPlotHiU[1] -> Draw("L");
  gPlotSt[1]  -> Draw("Z3");
  gPlotUn[0]  -> Draw("Z3");
  gPlotLoU[0] -> Draw("L");
  gPlotHiU[0] -> Draw("L");
  gPlotSt[0]  -> Draw("Z3");
  gPlotPyH[1] -> Draw("L");
  gPlotPyH[0] -> Draw("L");
  cout << "      Drew graphs." << endl;

  TLegend *legend3 = new TLegend(0.23, 0.14, 0.53, 0.45);
  legend3 -> SetTextSize(0.035);  
  legend3 -> SetFillColor(0);   
  legend3 -> SetLineColor(0);
  legend3 -> SetTextFont(42);
  legend3 -> AddEntry((TObject*)0, "p+p, #surds = 200 GeV", "");
  legend3 -> AddEntry((TObject*)0, "#pi^{0}+jet, anti-k_{T}, R = 0.5", "");   
  legend3 -> AddEntry(gDataLeg[0], "9 < E_{T}^{trig} < 11 GeV", "lf");    
  legend3 -> AddEntry(gDataLeg[1], "11 < E_{T}^{trig} < 15 GeV", "lf");
  legend3 -> AddEntry(gPyWfLeg, "PYTHIA-8 Monash tune", "l");
  legend3 -> Draw("SAME");
  cout << "      Created legend." << endl;

  TLegend *legend4 = new TLegend(0.53, 0.14, 0.83, 0.24);
  legend4 -> SetTextSize(0.035);  
  legend4 -> SetFillColor(0);   
  legend4 -> SetLineColor(0);
  legend4 -> SetTextFont(42);
  legend4 -> AddEntry((TObject*)0, "p+p, #surds = 200 GeV", "");
  legend4 -> AddEntry((TObject*)0, "#pi^{0}+jet, anti-k_{T}, R = 0.5", "");   
  //legend4 -> Draw("SAME");

  TLatex T1;
  T1.SetTextSize(0.05);
  T1.SetTextFont(42);
  T1.DrawLatex(10, 0.08, "STAR Preliminary");
  cout << "      Created preliminary label." << endl;

  fOut -> cd();
  c3   -> Write();
  c3   -> Close();
  cout << "    Made Quark Matter plot." << endl;

  // save graphs
  fOut -> cd();
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    gStat[iTrg]    -> Write();
    gUnf[iTrg]     -> Write();
    gPyF[iTrg]     -> Write();
    gPyH[iTrg]     -> Write();
    gPlotSt[iTrg]  -> Write();
    gPlotUn[iTrg]  -> Write();
    gPlotLoU[iTrg] -> Write();
    gPlotHiU[iTrg] -> Write();
    gPlotPyF[iTrg] -> Write();
    gPlotPyH[iTrg] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    fInStat[iTrg] -> cd();
    fInStat[iTrg] -> Close();
    fInUnf[iTrg]  -> cd();
    fInUnf[iTrg]  -> Close();
    fInPyF[iTrg]  -> cd();
    fInPyF[iTrg]  -> Close();
    fInPyH[iTrg]  -> cd();
    fInPyH[iTrg] -> Close();
  }
  cout << "  Finished making summary plot!" << endl;

}

// End ------------------------------------------------------------------------
