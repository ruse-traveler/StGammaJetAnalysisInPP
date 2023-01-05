// 'MakeUncertaintySummaryPlot.C'
// Derek Anderson
// 12.16.2021
//
// Use this produce a plot summarizing
// the systematic uncertainties on
// the data.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NTrg(3);
static const UInt_t NSys(2);
static const UInt_t NUnc(NSys + 1);
static const UInt_t NAvg(6);
static const UInt_t NVtx(4);
static const UInt_t NPlot(2);
static const UInt_t NMaxBin(20);
static const UInt_t NAvgBin(NAvg - 1);



void MakeUncertaintySummary() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning uncertainty summary plot making..." << endl;

  // io parameters
  const TString sOut("test.root");
  const TString sInP[NTrg]    = {"errorSumWithHists.et911r02pi0.d16m12y2021.root", "errorSumWithHists.et1115r02pi0.d16m12y2021.root", "errorSumWithHists.et1520r02pi0.d16m12y2021.root"};
  const TString sInG[NTrg]    = {"errorSumWithHists.et911r02gam.d16m12y2021.root", "errorSumWithHists.et1115r02gam.d16m12y2021.root", "errorSumWithHists.et1520r02gam.d16m12y2021.root"};
  const TString sInHist[NUnc] = {"hStatistics", "hSysUnfoldSmooth", "hSysDet"};

  // name parameters
  const TString sSuffP[NTrg] = {"Et911pi", "Et1115pi", "Et1520pi"};
  const TString sSuffG[NTrg] = {"Et911ga", "Et1115ga", "Et1520ga"};
  const TString sHist[NUnc]  = {"hStats", "hUnfold", "hDet"};
  const TString sBase[NSys]  = {"hPlotSysUnf", "hPlotSysDet"};
  const TString sStat("hPlotStat");

  // calculation parameters
  const Float_t yValue(1.);
  const Float_t yScales[NTrg]        = {1., 10., 100.};
  const Float_t xAvg[NTrg][NAvg]     = {0., 6., 11., 15., 20., 25., 30.};
  const Float_t xRangeP[NTrg][NPlot] = {{0.2, 30.}, {0.2, 30.}, {0.2, 30.}};
  const Float_t xRangeG[NTrg][NPlot] = {{0.2, 11.}, {0.2, 15.}, {0.2, 20.}};

  // histogram parameters
  const TString sTitle("");
  const TString sXtitle("p_{T}^{unfold} [GeV/c]");
  const TString sYtitle("percent uncertainty");

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  TFile *fInP[NTrg];
  TFile *fInG[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fInP[iTrg] = new TFile(sInP[iTrg].Data(), "read");
    fInG[iTrg] = new TFile(sInG[iTrg].Data(), "read");
    if (!fInP[iTrg] || !fInG[iTrg]) {
      cerr << "PANIC: couldn't open pi0 or gamma input file #" << iTrg << "!\n"
           << "       fInP = " << fInP[iTrg] << ", fInG = " << fInG[iTrg] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hInP[NTrg][NUnc];
  TH1D *hInG[NTrg][NUnc];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    for (UInt_t iUnc = 0; iUnc < NUnc; iUnc++) {
      hInP[iTrg][iUnc] = (TH1D*) fInP[iTrg] -> Get(sInHist[iUnc]);
      hInG[iTrg][iUnc] = (TH1D*) fInG[iTrg] -> Get(sInHist[iUnc]);
      if (!hInP[iTrg][iUnc] || !hInG[iTrg][iUnc]) {
        cerr << "PANIC: couldn't grab pi0 or gamma input histogram #" << iTrg << "!\n"
             << "       hInP = " << hInP[iTrg][iUnc] << ", hInG = " << hInG[iTrg][iUnc] << "\n"
             << endl;
        return;
      }

      // create names
      TString sNameP(sHist[iUnc].Data());
      TString sNameG(sHist[iUnc].Data());
      sNameP.Append(sSuffP[iTrg].Data());
      sNameG.Append(sSuffG[iTrg].Data());
      hInP[iTrg][iUnc] -> SetName(sNameP.Data());
      hInG[iTrg][iUnc] -> SetName(sNameG.Data());
    }
  }
  cout << "    Grabbed histograms." << endl;

  // create plot histograms
  TH1D *hStatP[NTrg];
  TH1D *hStatG[NTrg];
  TH1D *hSysP[NTrg][NSys];
  TH1D *hSysG[NTrg][NSys];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      hSysP[iTrg][iSys] = (TH1D*) hInP[iTrg][iSys + 1] -> Clone();
      hSysG[iTrg][iSys] = (TH1D*) hInG[iTrg][iSys + 1] -> Clone();
      hSysP[iTrg][iSys] -> Reset("ICES");
      hSysG[iTrg][iSys] -> Reset("ICES");

      // create systematic names
      TString sSysP(sBase[iSys].Data());
      TString sSysG(sBase[iSys].Data());
      sSysP.Append(sSuffP[iTrg].Data());
      sSysG.Append(sSuffG[iTrg].Data());
      hSysP[iTrg][iSys] -> SetName(sSysP.Data());
      hSysG[iTrg][iSys] -> SetName(sSysG.Data());
    }
    hStatP[iTrg] = (TH1D*) hInP[iTrg][0] -> Clone();
    hStatG[iTrg] = (TH1D*) hInG[iTrg][0] -> Clone();
    hStatP[iTrg] -> Reset("ICES");
    hStatG[iTrg] -> Reset("ICES");

    // create stat names
    TString sStatP(sStat.Data());
    TString sStatG(sStat.Data());
    sStatP.Append(sSuffP[iTrg].Data());
    sStatG.Append(sSuffG[iTrg].Data());
    hStatP[iTrg] -> SetName(sStatP.Data());
    hStatG[iTrg] -> SetName(sStatG.Data());
  }
  cout << "    Created plotting histograms." << endl;

  // sum uncertainties
  UInt_t   nStaBinP[NTrg];
  UInt_t   nStaBinG[NTrg];
  UInt_t   nSysBinP[NTrg][NSys];
  UInt_t   nSysBinG[NTrg][NSys];
  UInt_t   nAvgBinP[NTrg][NUnc];
  UInt_t   nAvgBinG[NTrg][NUnc];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    nStaBinP[iTrg] = 0;
    nStaBinG[iTrg] = 0;
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      nSysBinP[iTrg][iSys] = 0;
      nSysBinG[iTrg][iSys] = 0;
    }
    for (UInt_t iUnc = 0; iUnc < NUnc; iUnc++) {
      nAvgBinP[iTrg][iUnc] = 0;
      nAvgBinG[iTrg][iUnc] = 0;
    }
  }

  Double_t staErrP[NTrg][NMaxBin];
  Double_t staErrG[NTrg][NMaxBin];
  Double_t sysErrP[NTrg][NSys][NMaxBin];
  Double_t sysErrG[NTrg][NSys][NMaxBin];
  Double_t avgErrP[NTrg][NUnc][NAvgBin];
  Double_t avgErrG[NTrg][NUnc][NAvgBin];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    for (UInt_t iBinArr = 0; iBinArr < NMaxBin; iBinArr++) {
      staErrP[iTrg][iBinArr] = 0.;
      staErrG[iTrg][iBinArr] = 0.;
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        sysErrP[iTrg][iSys][iBinArr] = 0.;
        sysErrG[iTrg][iSys][iBinArr] = 0.;
      }
    }
    for (UInt_t iAvgArr = 0; iAvgArr < NAvgBin; iAvgArr++) {
      for (UInt_t iUnc = 0; iUnc < NUnc; iUnc++) {
        avgErrP[iTrg][iUnc][iAvgArr] = 0.;
        avgErrG[iTrg][iUnc][iAvgArr] = 0.;
      }
    }
  }

  // loop over triggers
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {

    // loop over bins
    const UInt_t nBins = (UInt_t) hStatP[iTrg] -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {

      // get stat uncertainties
      const Double_t xStat  = hInP[iTrg][0] -> GetBinCenter(iBin);
      const Double_t yStatP = hInP[iTrg][0] -> GetBinContent(iBin);
      const Double_t yStatG = hInG[iTrg][0] -> GetBinContent(iBin);
      const Double_t eStatP = hInP[iTrg][0] -> GetBinError(iBin);
      const Double_t eStatG = hInG[iTrg][0] -> GetBinError(iBin);

      // calculate stat percent
      Double_t pStatP(0.);
      Double_t pStatG(0.);
      if (yStatP > 0.) pStatP = eStatP / yStatP;
      if (yStatG > 0.) pStatG = eStatG / yStatG;

      // set stat bins
      const Bool_t isInStatRangeP = ((xStat > xRangeP[iTrg][0]) && (xStat < xRangeP[iTrg][1]));
      const Bool_t isInStatRangeG = ((xStat > xRangeG[iTrg][0]) && (xStat < xRangeG[iTrg][1]));
      if (isInStatRangeP) {
        staErrP[iTrg][nStaBinP[iTrg]] = pStatP;
        nStaBinP[iTrg]++;
      }
      if (isInStatRangeG) {
        staErrG[iTrg][nStaBinG[iTrg]] = pStatG;
        nStaBinG[iTrg]++;
      }

      // loop over systematics
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {

        // get sys uncertainties
        const Double_t xSys  = hInP[iTrg][iSys + 1] -> GetBinCenter(iBin);
        const Double_t ySysP = hInP[iTrg][iSys + 1] -> GetBinContent(iBin);
        const Double_t ySysG = hInG[iTrg][iSys + 1] -> GetBinContent(iBin);
        const Double_t eSysP = hInP[iTrg][iSys + 1] -> GetBinError(iBin);
        const Double_t eSysG = hInG[iTrg][iSys + 1] -> GetBinError(iBin);

        // calculate sys percent
        Double_t pSysP(0.);
        Double_t pSysG(0.);
        if (ySysP > 0.) pSysP = eSysP / ySysP;
        if (ySysG > 0.) pSysG = eSysG / ySysG;

        // set sys bins
        const Bool_t isInSysRangeP = ((xSys > xRangeP[iTrg][0]) && (xSys < xRangeP[iTrg][1]));
        const Bool_t isInSysRangeG = ((xSys > xRangeG[iTrg][0]) && (xSys < xRangeG[iTrg][1]));
        if (isInSysRangeP) {
          sysErrP[iTrg][nSysBinP[iTrg][iSys]] = pSysP;
          nSysBinP[iTrg][iSys]++;
        }
        if (isInSysRangeG) {
          sysErrG[iTrg][nSysBinG[iTrg][iSys]] = pSysG;
          nSysBinG[iTrg][iSys]++;
        }

      }  // end sys loop
    }  // end bin loop
  }  // end trigger loop
  cout << "    Calculated uncertainties." << endl;

  // check
  cout << "CHECK:" << endl;
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    cout << "  nStatBinP[" << iTrg << "] = " << nStaBinP[iTrg] << endl;
    cout << "  nStatBinG[" << iTrg << "] = " << nStaBinG[iTrg] << endl;
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      cout << "  nSysBinP[" << iTrg << " " << iSys << "] = " << nSysBinP[iTrg][iSys] << endl;
      cout << "  nSysBinG[" << iTrg << " " << iSys << "] = " << nSysBinG[iTrg][iSys] << endl;
    }
  }

  // save histograms
  fOut -> cd();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    for (UInt_t iUnc = 0; iUnc < NUnc; iUnc++) {
      hInP[iTrg][iUnc] -> Write();
      hInG[iTrg][iUnc] -> Write();
    }
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      hSysP[iTrg][iSys] -> Write();
      hSysG[iTrg][iSys] -> Write();
    }
    hStatP[iTrg] -> Write();
    hStatG[iTrg] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fInP[iTrg] -> cd();
    fInP[iTrg] -> Close();
    fInG[iTrg] -> cd();
    fInG[iTrg] -> Close();
  }
  cout << "  Finished making uncertainty summary plot!\n" << endl;

}

// End ------------------------------------------------------------------------
