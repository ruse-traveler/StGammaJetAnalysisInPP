// MakeClosureTestPlotForPaper.C'
// Derek Anderson
// 03.03.2023
//
// Use this to plot the output of
// 'MakeClosureTestPlot.C' in
// the style of the long paper.
//
// Adapted from code by Nihar Sahoo.

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

using namespace std;

// global constants
static const UInt_t NVtx(4);
static const UInt_t NTxt(5);
static const UInt_t NVals(2);
static const UInt_t NHist(2);
static const UInt_t NRange(2);
static const UInt_t NPoints(10);
static const UInt_t NFields(6);
static const UInt_t NBinMax(20);



void MakeClosureTestPlotForPaper() {

  // lower verbosity & set stat option
  gErrorIgnoreLevel = kError;
  gStyle            -> SetOptStat(0);
  cout << "\n  Plotting closure test for paper..." << endl;

  // io file parameters
  const TString sIn("closure/et911r05ff_rebinClosure/closureTestFF.forThesis_pTbinHuge.et911r05pi0.d9m2y2022.root");
  const TString sOut("closureTestForPaper_totErrorOnlyOnUnity.ffWithRff_pTbinHuge.et911r05pi0.d10m7y2023.root");

  // io histogram parameters
  const TString sParticle("hParticleWithStat");
  const TString sParUnity("hUnityWithStat");
  const TString sAverage[NHist] = {"hAverageLine",      "hVariationAverage"};
  const TString sRatio[NHist]   = {"hRatioAverageLine", "hRatioAverage"};

  // output histogram parameters
  const TString sNamePar("hParticleLevel");
  const TString sNameUnity("hUnityWithParErrors");
  const TString sNameParError("hParWithErrors");
  const TString sNameAvg[NHist] = {"hAvgValue",   "hAvgError"};
  const TString sNameRat[NHist] = {"hRatioValue", "hRatioError"};

  // style parameters
  const TString sTitleX("p_{T,jet}^{ch} [GeV/#it{c}]");
  const TString sTitleY("1/N_{trig} dN^{2}_{jets}/(dp_{T,jet}^{reco,ch} d#eta_{jet}) [GeV/#it{c}]^{-1}");
  const TString sTitleR("Ratio to PYTHIA");
  const UInt_t  fColPar(635);
  const UInt_t  fColErr(622);
  const UInt_t  fColUni(622);
  const UInt_t  fColOut(1);
  const UInt_t  fMarPar(1);
  const UInt_t  fMarErr(1);
  const UInt_t  fMarUni(1);
  const UInt_t  fMarAvg(1);
  const UInt_t  fMarRat(1);
  const UInt_t  fLinPar(1);
  const UInt_t  fLinErr(1);
  const UInt_t  fLinUni(1);
  const UInt_t  fLinOut(1);
  const UInt_t  fFilPar(0);
  const UInt_t  fFilErr(1001);
  const UInt_t  fFilUni(1001);
  const UInt_t  fFilOut(0);
  const UInt_t  fWidPar(2);
  const UInt_t  fWidErr(1);
  const UInt_t  fWidUni(1);
  const UInt_t  fWidOut(6);
  const UInt_t  fColAvg[NHist] = {1, 16};
  const UInt_t  fColRat[NHist] = {1, 16};
  const UInt_t  fLinAvg[NHist] = {2, 1};
  const UInt_t  fLinRat[NHist] = {2, 1};
  const UInt_t  fWidAvg[NHist] = {2, 1};
  const UInt_t  fWidRat[NHist] = {2, 1};
  const UInt_t  fFilAvg[NHist] = {0, 1001};
  const UInt_t  fFilRat[NHist] = {0, 1001};

  // frame parameters
  const Int_t    xNumBins(-10);
  const Int_t    yNumBins(1000);
  const Int_t    xBinPtRange[NRange] = {0,         29};
  const Double_t yRange[NRange]      = {0.0000004, 3.};
  //const Double_t rRange[NRange]      = {0.02,      1.8};
  const Double_t rRange[NRange]      = {0.37,      2.37};

  // plot parameters
  const UInt_t   width(700);
  const UInt_t   height(800);
  const Double_t xyLegUp[NVtx]    = {0.52, 0.69,  0.79,  0.96};  // 0.27 height, 0.27 width
  const Double_t xyLegDown[NVtx]  = {0.52, 0.69,  0.82,  0.93};  // (1/3) of the up lines => 0.09 =(scale by 37%)=> 0.24 height, =(scale by 89.7%)=> 0.30
  const Double_t xyPadRatio[NVtx] = {0.02, 0.007, 0.917, 0.377};  // 37% the height of the canvas, 89.7% the width of the canvas
  const TString  sTxt[NTxt]       = {"#bf{STAR}", "Closure Test", "anti-k_{T}, R=0.5", "9 < E_{T}^{trig} < 11 GeV", "#pi^{0}+jet"};
  const TString  sJetTxt("Unfolded");
  const TString  sParTxt("PYTHIA");
  const TString  sUniTxt("Data uncertainty");

  // distribution parameters
  const Bool_t  doLowPtCutoff(true);
  const Bool_t  showUnityError(true);
  const Bool_t  showPythiaError(false);
  const Bool_t  putDataErrorsOnUnity(true);
  const Bool_t  putDataErrorsOnPythia(false);
  const UInt_t  fDataset(1);
  const Float_t minJetPt(3.);
  const Float_t maxJetPt(xBinPtRange[1]);

  // data points
  const Double_t sysData_et911r02[NPoints][NFields] = {
    {0.41991, 1.67912,     0.41991,  0.58009, 0.0727417,   0.0727417},
    {1.41991, 0.565804,    0.41991,  0.58009, 0.0227303,   0.0227303},
    {2.69355, 0.157091,    0.693549, 1.30645, 0.00952979,  0.00952979},
    {4.92234, 0.0429032,   0.922344, 1.07766, 0.00197916,  0.00197916},
    {6.92234, 0.0251721,   0.922344, 1.07766, 0.00126535,  0.00126535},
    {9.32606, 0.0148911,   1.32606,  1.67394, 0.00101062,  0.00101062},
    {12.6927, 0.00661861,  1.6927,   2.3073,  0.000567637, 0.000567637},
    {16.9582, 0.00230692,  1.95822,  3.04178, 0.000276657, 0.000276657},
    {21.9582, 0.000774539, 1.95822,  3.04178, 0.000182686, 0.000182686},
    {26.9582, 0.000162226, 1.95822,  3.04178, 4.17115e-05, 4.17115e-05}
  };
  const Double_t statData_et911r02[NPoints][NFields] = {
    {0.419903, 1.67912,     0.419903, 0.580097, 0.0103549,   0.0103549},
    {1.4199,   0.565804,    0.419903, 0.580097, 0.00601886,  0.00601886},
    {2.69352,  0.157091,    0.693524, 1.30648,  0.00232342,  0.00232342},
    {4.92234,  0.0429032,   0.922337, 1.07766,  0.00129114,  0.00129114},
    {6.92234,  0.0251721,   0.922337, 1.07766,  0.00103774,  0.00103774},
    {9.32604,  0.0148911,   1.32604,  1.67396,  0.000656481, 0.000656481},
    {12.6927,  0.00661861,  1.69267,  2.30733,  0.000376287, 0.000376287},
    {16.9582,  0.00230692,  1.95824,  3.04176,  0.00020196,  0.00020196},
    {21.9582,  0.000774539, 1.95824,  3.04176,  0.000112017, 0.000112017},
    {26.9582,  0.000162226, 1.95824,  3.04176,  4.83076e-05, 4.83076e-05}
  };
  const Double_t sysData_et911r05[NPoints][NFields] = {
    {0.415792, 1.07316,     0.415792, 0.584208, 0.0423119,   0.0423119},
    {1.41579,  0.289051,    0.415792, 0.584208, 0.00695692,  0.00695692},
    {2.67923,  0.0981596,   0.679232, 1.32077,  0.00635806,  0.00635806},
    {4.94215,  0.0543701,   0.942145, 1.05785,  0.00259349,  0.00259349},
    {6.94215,  0.038216,    0.942145, 1.05785,  0.00218521,  0.00218521},
    {9.37015,  0.0251227,   1.37015,  1.62985,  0.00172022,  0.00172022},
    {12.77,    0.0134148,   1.76996,  2.23004,  0.00136424,  0.00136424},
    {17.0654,  0.00443086,  2.06537,  2.93463,  0.000670039, 0.000670039},
    {22.0654,  0.0014075,   2.06537,  2.93463,  0.000339134, 0.000339134},
    {27.0654,  0.000544651, 2.06537,  2.93463,  0.00029686,  0.00029686}
  };
  const Double_t statData_et911r05[NPoints][NFields] = {
    {0.415807, 1.07316,     0.415807, 0.584193, 0.0169683,   0.0169683},
    {1.41581,  0.289051,    0.415807, 0.584193, 0.00670517,  0.00670517},
    {2.67928,  0.0981596,   0.679283, 1.32072,  0.00253153,  0.00253153},
    {4.94214,  0.0543701,   0.942141, 1.05786,  0.00180289,  0.00180289},
    {6.94214,  0.038216,    0.942141, 1.05786,  0.00148595,  0.00148595},
    {9.37014,  0.0251227,   1.37014,  1.62986,  0.000988172, 0.000988172},
    {12.7699,  0.0134148,   1.76995,  2.23005,  0.000616403, 0.000616403},
    {17.0693,  0.00443086,  2.06935,  2.93065,  0.000302308, 0.000302308},
    {22.0693,  0.0014075,   2.06935,  2.93065,  0.000162512, 0.000162512},
    {27.0693,  0.000544651, 2.06935,  2.93065,  0.000105142, 0.000105142}
  };

  // open files
  TFile *fIn  = new TFile(sIn.Data(),  "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut || !fIn) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fIn = " << fIn
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input
  TH1D *hInAverage[NHist];
  TH1D *hInRatio[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInAverage[iHist] = (TH1D*) fIn -> Get(sAverage[iHist].Data());
    hInRatio[iHist]   = (TH1D*) fIn -> Get(sRatio[iHist].Data());
    if (!hInAverage[iHist] || !hInRatio[iHist]) {
      cerr << "PANIC: couldn't grab an input histogram!\n"
           << "       hAverage[" << iHist << "] = " << hInAverage[iHist] << ", hRatio[" << iHist << "] = " << hInRatio[iHist] << "\n"
           << endl;
      return;
    }
  }

  TH1D *hInParticle = (TH1D*) fIn -> Get(sParticle.Data());
  if (!hInParticle) {
    cerr << "PANIC: couldn't grab input particle histogram!\n" << endl;
    return;
  }

  TH1D *hInUnity;
  if (showUnityError) {
    hInUnity = (TH1D*) fIn -> Get(sParUnity.Data());
    if (!hInUnity) {
      cerr << "PANIC: couldn't grab input unity histogram!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // copy input to output histograms
  TH1D *hUnity;
  TH1D *hParticle;
  TH1D *hParError;
  TH1D *hAverage[NHist];
  TH1D *hRatio[NHist];

  hParticle = (TH1D*) hInParticle -> Clone();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] = (TH1D*) hInAverage[iHist] -> Clone();
    hRatio[iHist]   = (TH1D*) hInRatio[iHist]   -> Clone();
  }
  if (showUnityError)  hUnity    = (TH1D*) hInUnity    -> Clone();
  if (showPythiaError) hParError = (TH1D*) hInParticle -> Clone();

  // parse dataset selection
  Double_t sysDataPoints[NPoints][NFields];
  Double_t statDataPoints[NPoints][NFields];
  for (UInt_t iPoint = 0; iPoint < NPoints; iPoint++) {
    for (UInt_t iField = 0; iField < NFields; iField++) {
      switch (fDataset) {
        case 0:
          sysDataPoints[iPoint][iField]  = sysData_et911r02[iPoint][iField];
          statDataPoints[iPoint][iField] = statData_et911r02[iPoint][iField];
          break;
        case 1:
          sysDataPoints[iPoint][iField]  = sysData_et911r05[iPoint][iField];
          statDataPoints[iPoint][iField] = statData_et911r05[iPoint][iField];
          break;
        default:
          sysDataPoints[iPoint][iField]  = sysData_et911r02[iPoint][iField];
          statDataPoints[iPoint][iField] = statData_et911r02[iPoint][iField];
          break;
      }
    }  // end field loop
  }  // end point loop

  // calculate total data error
  Double_t totDataPoints[NPoints][NFields];
  for (UInt_t iPoint = 0; iPoint < NPoints; iPoint++) {

    // add errors in quadrature
    const Double_t sysFracErr  = sysDataPoints[iPoint][4]  / sysDataPoints[iPoint][1];
    const Double_t statFracErr = statDataPoints[iPoint][4] / statDataPoints[iPoint][1];
    const Double_t totFracErr  = TMath::Sqrt((sysFracErr * sysFracErr) + (statFracErr * statFracErr));
    const Double_t totErr      = totFracErr * sysDataPoints[iPoint][1];

    // set total data points
    totDataPoints[iPoint][0] = sysDataPoints[iPoint][0];
    totDataPoints[iPoint][1] = sysDataPoints[iPoint][1];
    totDataPoints[iPoint][2] = sysDataPoints[iPoint][2];
    totDataPoints[iPoint][3] = sysDataPoints[iPoint][3];
    totDataPoints[iPoint][4] = totErr;
    totDataPoints[iPoint][5] = totErr;
  }  // end point loop

  // get no. of bins
  const UInt_t nBinsPar = hInParticle   -> GetNbinsX();
  const UInt_t nBinsAvg = hInAverage[0] -> GetNbinsX();
  const UInt_t nBinsRat = hInRatio[0]   -> GetNbinsX();

  UInt_t nBinsUni = 0;
  if (showUnityError) {
    nBinsUni = hInUnity -> GetNbinsX();
  }

  // adjust particle-level errors if needed
  const Bool_t useDataErrors = (putDataErrorsOnUnity || putDataErrorsOnPythia);
  if (useDataErrors) {
    for (UInt_t iBinPar = 0; iBinPar < nBinsPar; iBinPar++) {

      // get particle-level info
      const Double_t binCenterPar = hInParticle -> GetBinCenter(iBinPar);
      const Double_t binValuePar  = hInParticle -> GetBinContent(iBinPar);

      // get relevant datapoint
      for (UInt_t iPoint = 0; iPoint < NPoints; iPoint++) {
        const Double_t binLowEdgeDat  = totDataPoints[iPoint][0] - totDataPoints[iPoint][2];
        const Double_t binHighEdgeDat = totDataPoints[iPoint][0] + totDataPoints[iPoint][3];
        const Double_t binFracDat     = totDataPoints[iPoint][4] / totDataPoints[iPoint][1];
        const Bool_t   isInDataBin    = ((binCenterPar >= binLowEdgeDat) && (binCenterPar < binHighEdgeDat));
        if (isInDataBin && (totDataPoints[iPoint][1] > 0.)) {
          if (putDataErrorsOnUnity)  hInUnity    -> SetBinError(iBinPar, binFracDat);
          if (putDataErrorsOnPythia) hInParticle -> SetBinError(iBinPar, binFracDat * binValuePar);
          break;
        }
      }  // end point loop
    }  // end particle bin loop
    cout << "    Adjusted particle-level errors." << endl;
  }  // end error adjustment

  UInt_t   nBinToDrawPar;
  UInt_t   nBinToDrawAvg[NHist];
  UInt_t   nBinToDrawRat[NHist];
  Double_t xyPar[NVals][NBinMax];
  Double_t xyErr[NVals][NBinMax];
  Double_t xyAvg[NHist][NVals][NBinMax];
  Double_t xyRat[NHist][NVals][NBinMax];

  nBinToDrawPar = 0;
  for (UInt_t iBin = 0; iBin < NBinMax; iBin++) {
    xyPar[0][iBin] = 0.;
    xyPar[1][iBin] = 0.;
    xyErr[0][iBin] = 0.;
    xyErr[1][iBin] = 0.;
  }
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    nBinToDrawAvg[iHist] = 0;
    nBinToDrawRat[iHist] = 0;
    for (UInt_t iBin = 0; iBin < NBinMax; iBin++) {
      xyAvg[iHist][0][iBin] = 0.;
      xyAvg[iHist][1][iBin] = 0.;
      xyRat[iHist][0][iBin] = 0.;
      xyRat[iHist][1][iBin] = 0.;
    }
  }

  // apply low pt cutoff (if needed)
  if (doLowPtCutoff) {

    // apply cutoff to particle spectrum
    hParticle -> Reset("ICES");
    for (UInt_t iBinPar = 1; iBinPar < (nBinsPar + 1); iBinPar++) {
      const Double_t binWidthPar  = hInParticle -> GetBinWidth(iBinPar);
      const Double_t binCenterPar = hInParticle -> GetBinCenter(iBinPar);
      const Double_t binValuePar  = hInParticle -> GetBinContent(iBinPar);
      const Double_t binErrorPar  = hInParticle -> GetBinError(iBinPar);
      if ((binCenterPar >= minJetPt) && (binCenterPar < maxJetPt)) {

        // set histogram values
        hParticle -> SetBinContent(iBinPar, binValuePar);
        hParticle -> SetBinError(iBinPar,   binErrorPar);

        // store values
        xyPar[0][nBinToDrawPar] = binCenterPar;
        xyPar[1][nBinToDrawPar] = binValuePar;
        xyErr[0][nBinToDrawPar] = binWidthPar;
        xyErr[1][nBinToDrawPar] = binErrorPar;
        ++nBinToDrawPar;
      }
    }  // end particle bin loop

    // apply cutoff to unity histogram
    if (showUnityError) {

      hUnity -> Reset("ICES");
      for (UInt_t iBinUni = 1; iBinUni < (nBinsUni + 1); iBinUni++) {
        const Double_t binWidthUni  = hInUnity -> GetBinWidth(iBinUni);
        const Double_t binCenterUni = hInUnity -> GetBinCenter(iBinUni);
        const Double_t binValueUni  = hInUnity -> GetBinContent(iBinUni);
        const Double_t binErrorUni  = hInUnity -> GetBinError(iBinUni);
        if ((binCenterUni >= minJetPt) && (binCenterUni < maxJetPt)) {
          hUnity -> SetBinContent(iBinUni, binValueUni);
          hUnity -> SetBinError(iBinUni,   binErrorUni);
        }
      }  // end unity bin loop
    }

    // apply cutoff to averages and ratios
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      hAverage[iHist] -> Reset("ICES");
      hRatio[iHist]   -> Reset("ICES");
      for (UInt_t iBinAvg = 1; iBinAvg < (nBinsAvg + 1); iBinAvg++) {
        const Double_t binCenterAvg = hInAverage[iHist] -> GetBinCenter(iBinAvg);
        const Double_t binValueAvg  = hInAverage[iHist] -> GetBinContent(iBinAvg);
        const Double_t binErrorAvg  = hInAverage[iHist] -> GetBinError(iBinAvg);
        if ((binCenterAvg >= minJetPt) && (binCenterAvg < maxJetPt)) {

          // set histogram values
          hAverage[iHist] -> SetBinContent(iBinAvg, binValueAvg);
          hAverage[iHist] -> SetBinError(iBinAvg,   binErrorAvg);

          // store values
          xyAvg[iHist][0][nBinToDrawAvg[iHist]] = binCenterAvg;
          xyAvg[iHist][1][nBinToDrawAvg[iHist]] = binValueAvg;
          ++nBinToDrawAvg[iHist];
        }
      }  // end average bin loop
      for (UInt_t iBinRat = 1; iBinRat < (nBinsRat + 1); iBinRat++) {
        const Double_t binCenterRat = hInRatio[iHist] -> GetBinCenter(iBinRat);
        const Double_t binValueRat  = hInRatio[iHist] -> GetBinContent(iBinRat);
        const Double_t binErrorRat  = hInRatio[iHist] -> GetBinError(iBinRat);
        if ((binCenterRat >= minJetPt) && (binCenterRat < maxJetPt)) {

          // set histogram values
          hRatio[iHist] -> SetBinContent(iBinRat, binValueRat);
          hRatio[iHist] -> SetBinError(iBinRat,   binErrorRat);

          // store values
          xyRat[iHist][0][nBinToDrawRat[iHist]] = binCenterRat;
          xyRat[iHist][1][nBinToDrawRat[iHist]] = binValueRat;
          ++nBinToDrawRat[iHist];
        }
      }  // end average bin loop
    }  // end histogram loop
    cout << "    Cut off spectra at low pt." << endl;
  }  // end low pt cutoff

  // set histogram names
  hParticle -> SetName(sNamePar.Data());
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> SetName(sNameAvg[iHist].Data());
    hRatio[iHist]   -> SetName(sNameRat[iHist].Data());
  }
  if (showUnityError)  hUnity    -> SetName(sNameUnity.Data());
  if (showPythiaError) hParError -> SetName(sNameParError.Data());

  // set styles
  if (showPythiaError) {
    hParticle -> SetMarkerColor(fColPar);
    hParticle -> SetMarkerStyle(fMarPar);
    hParticle -> SetLineColor(fColPar);
    hParticle -> SetLineStyle(fLinPar);
    hParticle -> SetLineWidth(fWidPar);
    hParticle -> SetFillColor(fColPar);
    hParticle -> SetFillStyle(fFilPar);
    hParError -> SetMarkerColor(fColErr);
    hParError -> SetMarkerStyle(fMarErr);
    hParError -> SetLineColor(fColErr);
    hParError -> SetLineStyle(fLinErr);
    hParError -> SetLineWidth(fWidErr);
    hParError -> SetFillColor(fColErr);
    hParError -> SetFillStyle(fFilErr);
  } else {
    hParticle -> SetMarkerColor(fColPar);
    hParticle -> SetMarkerStyle(fMarPar);
    hParticle -> SetLineColor(fColPar);
    hParticle -> SetLineStyle(fLinPar);
    hParticle -> SetLineWidth(fWidPar);
    hParticle -> SetFillColor(fColPar);
    hParticle -> SetFillStyle(fFilPar);
  }
  if (showUnityError) {
    hUnity -> SetMarkerColor(fColUni);
    hUnity -> SetMarkerStyle(fMarUni);
    hUnity -> SetLineColor(fColUni);
    hUnity -> SetLineStyle(fLinUni);
    hUnity -> SetLineWidth(fWidUni);
    hUnity -> SetFillColor(fColUni);
    hUnity -> SetFillStyle(fFilUni);
  }
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> SetMarkerColor(fColAvg[iHist]);
    hAverage[iHist] -> SetMarkerStyle(fMarAvg);
    hAverage[iHist] -> SetLineColor(fColAvg[iHist]);
    hAverage[iHist] -> SetLineStyle(fLinAvg[iHist]);
    hAverage[iHist] -> SetLineWidth(fWidAvg[iHist]);
    hAverage[iHist] -> SetFillColor(fColAvg[iHist]);
    hAverage[iHist] -> SetFillStyle(fFilAvg[iHist]);
    hRatio[iHist]   -> SetMarkerColor(fColRat[iHist]);
    hRatio[iHist]   -> SetMarkerStyle(fMarRat);
    hRatio[iHist]   -> SetLineColor(fColRat[iHist]);
    hRatio[iHist]   -> SetLineStyle(fLinRat[iHist]);
    hRatio[iHist]   -> SetLineWidth(fWidRat[iHist]);
    hRatio[iHist]   -> SetFillColor(fColRat[iHist]);
    hRatio[iHist]   -> SetFillStyle(fFilRat[iHist]);
  }
  cout << "    Set styles." << endl;

  // create outline histograms
  TH1D *hAvgOutline = (TH1D*) hAverage[1] -> Clone();
  hAvgOutline -> SetName("hAverageOutline");
  hAvgOutline -> SetLineColor(fColOut);
  hAvgOutline -> SetLineStyle(fLinOut);
  hAvgOutline -> SetLineWidth(fWidOut);
  hAvgOutline -> SetFillColor(fColOut);
  hAvgOutline -> SetFillStyle(fFilOut);

  TH1D *hRatOutline = (TH1D*) hRatio[1] -> Clone();
  hRatOutline -> SetName("hRatioOutline");
  hRatOutline -> SetLineColor(fColOut);
  hRatOutline -> SetLineStyle(fLinOut);
  hRatOutline -> SetLineWidth(fWidOut);
  hRatOutline -> SetFillColor(fColOut);
  hRatOutline -> SetFillStyle(fFilOut);

  TH1D *hParOutline = (TH1D*) hParticle -> Clone();
  hParOutline -> SetName("hParticleOutline");
  hParOutline -> SetLineColor(fColPar);
  hParOutline -> SetLineStyle(fLinOut);
  hParOutline -> SetLineWidth(fWidOut);
  hParOutline -> SetFillColor(fColPar);
  hParOutline -> SetFillStyle(fFilOut);
 
  TH1D *hUniOutline = (TH1D*) hUnity -> Clone();
  hUniOutline -> SetName("hUnityOutline");
  hUniOutline -> SetLineColor(fColPar);
  hUniOutline -> SetLineStyle(fLinOut);
  hUniOutline -> SetLineWidth(fWidOut);
  hUniOutline -> SetFillColor(fColPar);
  hUniOutline -> SetFillStyle(fFilOut);
  cout << "    Made outline histograms." << endl;

  // make legend histograms
  TH1D *hAvgLegend = (TH1D*) hAverage[1] -> Clone();
  hAvgLegend -> SetName("hAverageLegend");
  hAvgLegend -> SetLineColor(fColOut);

  TH1D *hParLegend = (TH1D*) hParticle -> Clone();
  hParLegend -> SetName("hParticleLegend");
  if (showPythiaError) {
    hParLegend -> SetLineColor(fColPar);
    hParLegend -> SetLineWidth(fWidPar);
    hParLegend -> SetFillColor(fColErr);
    hParLegend -> SetFillStyle(fFilErr);
  } else {
    hParLegend -> SetLineColor(fColPar);
    hParLegend -> SetLineWidth(fWidPar);
    hParLegend -> SetFillColor(fColPar);
    hParLegend -> SetFillStyle(fFilPar);
  }

  TH1D *hUniLegend = (TH1D*) hUnity -> Clone();
  hUniLegend -> SetName("hUnityLegend");
  if (showUnityError) {
    hUniLegend -> SetLineColor(fColPar);
    hUniLegend -> SetLineWidth(fWidPar);
    hUniLegend -> SetFillColor(fColErr);
    hUniLegend -> SetFillStyle(fFilErr);
  } else {
    hUniLegend -> SetLineColor(fColPar);
    hUniLegend -> SetLineWidth(fWidPar);
    hUniLegend -> SetFillColor(fColPar);
    hUniLegend -> SetFillStyle(fFilPar);
  }
  cout << "    Made legend histograms." << endl;

  // make frame histograms
  TH2D *hFrameJets = new TH2D("hFrameJets", "" , xNumBins, xBinPtRange[0], xBinPtRange[1], yNumBins, yRange[0], yRange[1]);
  hFrameJets -> GetXaxis() -> SetLabelFont(42);
  hFrameJets -> GetXaxis() -> SetLabelSize(0.0);
  hFrameJets -> GetXaxis() -> SetTitleSize(0.033);
  hFrameJets -> GetXaxis() -> SetTitleFont(42);
  hFrameJets -> GetYaxis() -> SetTitle(sTitleY.Data());
  hFrameJets -> GetYaxis() -> CenterTitle(true);
  hFrameJets -> GetYaxis() -> SetLabelFont(42);
  hFrameJets -> GetYaxis() -> SetLabelSize(0.05);
  hFrameJets -> GetYaxis() -> SetTitleSize(0.035);
  hFrameJets -> GetYaxis() -> SetTitleOffset(2.2);
  hFrameJets -> GetYaxis() -> SetTitleFont(42);
  hFrameJets -> GetZaxis() -> SetLabelFont(42);
  hFrameJets -> GetZaxis() -> SetLabelSize(0.035);
  hFrameJets -> GetZaxis() -> SetTitleSize(0.035);
  hFrameJets -> GetZaxis() -> SetTitleFont(42);

  TH2F *hFrameRatio = new TH2F("hFrameRatio", "", xNumBins, xBinPtRange[0], xBinPtRange[1], yNumBins, rRange[0], rRange[1]);
  hFrameRatio -> GetXaxis() -> SetTitle(sTitleX.Data());
  hFrameRatio -> GetXaxis() -> CenterTitle(true);
  hFrameRatio -> GetXaxis() -> SetLabelFont(42);
  hFrameRatio -> GetXaxis() -> SetLabelOffset(0.001);
  hFrameRatio -> GetXaxis() -> SetLabelSize(0.14);
  hFrameRatio -> GetXaxis() -> SetTitleSize(0.1);
  hFrameRatio -> GetXaxis() -> SetTitleOffset(1.4);
  hFrameRatio -> GetXaxis() -> SetTitleFont(42);
  hFrameRatio -> GetXaxis() -> SetTickLength(0.06);
  hFrameRatio -> GetYaxis() -> SetTickLength(0.04);
  hFrameRatio -> GetYaxis() -> SetTitle(sTitleR.Data());
  hFrameRatio -> GetYaxis() -> CenterTitle(true);
  hFrameRatio -> GetYaxis() -> SetLabelFont(42);
  hFrameRatio -> GetYaxis() -> SetLabelSize(0.12);
  hFrameRatio -> GetYaxis() -> SetTitleSize(0.09);
  hFrameRatio -> GetYaxis() -> SetNdivisions(505);
  hFrameRatio -> GetYaxis() -> SetTitleOffset(0.9);
  hFrameRatio -> GetYaxis() -> SetTitleFont(42);
  hFrameRatio -> GetZaxis() -> SetLabelFont(42);
  hFrameRatio -> GetZaxis() -> SetLabelSize(0.035);
  hFrameRatio -> GetZaxis() -> SetTitleSize(0.035);
  hFrameRatio -> GetZaxis() -> SetTitleFont(42);
  cout << "    Made frame histograms." << endl;

  // make text boxes
  TLegend *legTxtUp = new TLegend(xyLegUp[0], xyLegUp[1], xyLegUp[2], xyLegUp[3]);
  legTxtUp -> SetTextFont(42);
  legTxtUp -> SetTextSize(0.036);  
  legTxtUp -> SetFillColor(0);   
  legTxtUp -> SetLineColor(0);
  legTxtUp -> AddEntry((TObject*) 0, sTxt[0].Data(), "");
  legTxtUp -> AddEntry((TObject*) 0, sTxt[1].Data(), "");
  legTxtUp -> AddEntry((TObject*) 0, sTxt[2].Data(), "");
  legTxtUp -> AddEntry((TObject*) 0, sTxt[3].Data(), "");
  legTxtUp -> AddEntry(hAvgLegend,   sJetTxt.Data(), "f");
  if (showPythiaError) {
    legTxtUp -> AddEntry(hParLegend, sParTxt.Data(), "lf");
  } else {
    legTxtUp -> AddEntry(hParLegend, sParTxt.Data(), "l");
  }

  TLegend *legTxtDown = new TLegend(xyLegDown[0], xyLegDown[1], xyLegDown[2], xyLegDown[3]);
  legTxtDown -> SetTextFont(42);
  legTxtDown -> SetTextSize(0.085);  
  legTxtDown -> SetFillColor(0);   
  legTxtDown -> SetLineColor(0);
  legTxtDown -> AddEntry(hAvgLegend, sJetTxt.Data(), "f");
  if (showUnityError) {
    legTxtDown -> AddEntry(hUniLegend, sUniTxt.Data(), "lf");
  }
  cout << "    Made text boxes." << endl;

  // make lines
  TLine *lUnity = new TLine(xBinPtRange[0], 1., xBinPtRange[1], 1.);
  lUnity -> SetLineColor(1);
  lUnity -> SetLineStyle(2);

  TGraph       *gParticle = new TGraph(nBinToDrawPar, xyPar[0], xyPar[1]);
  TGraphErrors *gParError = new TGraphErrors(nBinToDrawPar, xyPar[0], xyPar[1], xyErr[0], xyErr[1]);
  gParticle -> SetName("gParticle");
  if (showPythiaError) {
    gParticle -> SetMarkerColor(fColPar);
    gParticle -> SetMarkerStyle(fMarPar);
    gParticle -> SetLineColor(fColPar);
    gParticle -> SetLineStyle(fLinPar);
    gParticle -> SetLineWidth(fWidPar);
    gParticle -> SetFillColor(fColPar);
    gParticle -> SetFillStyle(fFilPar);
  } else {
    gParticle -> SetMarkerColor(fColPar);
    gParticle -> SetMarkerStyle(fMarPar);
    gParticle -> SetLineColor(fColPar);
    gParticle -> SetLineStyle(fLinPar);
    gParticle -> SetLineWidth(fWidPar);
    gParticle -> SetFillColor(fColPar);
    gParticle -> SetFillStyle(fFilPar);
  }
  gParError -> SetName("gParError");
  gParError -> SetMarkerColor(fColErr);
  gParError -> SetMarkerStyle(fMarErr);
  gParError -> SetLineColor(fColErr);
  gParError -> SetLineStyle(fLinErr);
  gParError -> SetLineWidth(fWidErr);
  gParError -> SetFillColor(fColErr);
  gParError -> SetFillStyle(fFilErr);

  TGraph *gAverage = new TGraph(nBinToDrawAvg[0], xyAvg[0][0], xyAvg[0][1]);
  TGraph *gRatio   = new TGraph(nBinToDrawRat[0], xyRat[0][0], xyRat[0][1]);
  gAverage -> SetName("gAverage");
  gAverage -> SetLineColor(fColAvg[0]);
  gAverage -> SetLineStyle(fLinAvg[0]);
  gAverage -> SetLineWidth(fWidAvg[0]);
  gRatio   -> SetName("gRatio");
  gRatio   -> SetLineColor(fColRat[0]);
  gRatio   -> SetLineStyle(fLinRat[0]);
  gRatio   -> SetLineWidth(fWidRat[0]);
  cout << "    Made lines." << endl;

  // make canvas
  TCanvas *cPlot = new TCanvas("cPlot", "canvas", width, height);
  cPlot       -> SetFillColor(0);
  cPlot       -> SetBorderMode(0);
  cPlot       -> SetBorderSize(2);
  cPlot       -> SetLogy(1);
  cPlot       -> SetTickx(1);
  cPlot       -> SetTicky(1);
  cPlot       -> SetLeftMargin(0.1794479);
  cPlot       -> SetRightMargin(0.0993865);
  cPlot       -> SetTopMargin(0.01536313);
  cPlot       -> SetBottomMargin(0.3743017);
  cPlot       -> SetFrameBorderMode(0);
  cPlot       -> SetFrameBorderMode(0);
  cPlot       -> cd();
  hFrameJets  -> DrawCopy("9");
  if (showPythiaError) {
    if (doLowPtCutoff) {
      gParError   -> Draw("3");
      hParOutline -> Draw("E5 SAME");
      gParticle   -> Draw("L");
      hAverage[1] -> Draw("E5 SAME");
      hAvgOutline -> Draw("E5 SAME");
      gAverage    -> Draw("L");
    } else {
      hParError   -> Draw("E5 same");
      hParOutline -> Draw("E5 same");
      hParticle   -> Draw("][ L HIST SAME");
      hAverage[1] -> Draw("E5 SAME");
      hAvgOutline -> Draw("E5 SAME");
      hAverage[0] -> Draw("][ L HIST SAME");
    }
  } else {
    if (doLowPtCutoff) {
      hAverage[1]  -> Draw("E5 SAME");
      hAvgOutline  -> Draw("E5 SAME");
      gAverage     -> Draw("L");
      gParticle    -> Draw("L");
    } else {
      hAverage[1] -> Draw("E5 SAME");
      hAvgOutline -> Draw("E5 SAME");
      hAverage[0] -> Draw("][ L HIST SAME");
      hParticle   -> Draw("][ L HIST SAME");
    }
  }
  legTxtUp -> Draw("SAME");

  // create pad for ratio
  TPad *pRatio = new TPad("pRatio", "pRatio", xyPadRatio[0], xyPadRatio[1], xyPadRatio[2], xyPadRatio[3]);
  cPlot       -> cd();
  pRatio      -> Draw();
  pRatio      -> cd();  
  pRatio      -> SetFillColor(0);
  pRatio      -> SetBorderMode(0);
  pRatio      -> SetBorderSize(2);
  pRatio      -> SetLogy(1);
  pRatio      -> SetTicky(1);
  pRatio      -> SetLeftMargin(0.1784615);
  pRatio      -> SetRightMargin(0.01692308);
  pRatio      -> SetTopMargin(0.007434944);
  pRatio      -> SetBottomMargin(0.3945725);
  pRatio      -> SetFrameBorderMode(0);
  pRatio      -> SetFrameBorderMode(0);
  hFrameRatio -> DrawCopy("9");
  if (showUnityError) {
    if (doLowPtCutoff) {
      hUnity      -> Draw("E5 SAME");
      hUniOutline -> Draw("E5 SAME");
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      gRatio      -> Draw("L");
    } else {
      hUnity      -> Draw("E5 SAME");
      hUniOutline -> Draw("E5 SAME");
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      hRatio[0]   -> Draw("][ L HIST SAME");
    }
  } else {
    if (doLowPtCutoff) {
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      gRatio      -> Draw("L");
    } else {
      hRatio[1]   -> Draw("E5 SAME");
      hRatOutline -> Draw("E5 SAME");
      hRatio[0]   -> Draw("][ L HIST SAME");
    }
  }
  lUnity     -> Draw();
  legTxtDown -> Draw();
  fOut       -> cd();
  cPlot      -> Write();
  cPlot      -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hAverage[iHist] -> Write();
    hRatio[iHist]   -> Write();
  }
  hParticle -> Write();
  if (showUnityError) {
    hUnity      -> Write();
    hUniOutline -> Write();
  }
  if (showPythiaError) {
    hParError   -> Write();
    hParOutline -> Write();
  }
  if (doLowPtCutoff) {
    gParticle -> Write();
    gAverage  -> Write();
    gRatio    -> Write();
    if (showPythiaError) {
      gParError -> Write();
    }
  }
  hAvgOutline -> Write();
  hRatOutline -> Write();
  hAvgLegend  -> Write();
  hParLegend  -> Write();
  if (showUnityError) {
    hUniLegend  -> Write();
  }
  hFrameJets  -> Write();
  hFrameRatio -> Write();

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished plotting closure test!\n" << endl;

}

// End ------------------------------------------------------------------------
