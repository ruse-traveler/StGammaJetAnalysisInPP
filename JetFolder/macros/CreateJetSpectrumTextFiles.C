// 'CreateJetSpectrumTextFiles.C'
// Derek Anderson
// 07.08.2023
// 
// Reads in rootfiles for fig.s 20 in the gamma+jet
// PRC (or something equivalent) and sgats out text
// files with their data in a table.

#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TGraphErrors.h"

using namespace std;

static const Ssiz_t NStat(3);
static const Ssiz_t NSys(3);



void CreateJetSpectrumTextFiles() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Begining data tables creation script..." << endl;

  // input parameters
  const TString sInStat[NStat] = {
    "../../misc/ThesisAndPaperFiles/RootFilesForPlotsPP/Results/gammaVsPionSummary.onlyPy8_pTbinHuge.et920r05.d8m11y2021.root",
    "../../misc/ThesisAndPaperFiles/RootFilesForPlotsPP/Results/gammaVsPionSummary.onlyPy8_pTbinHuge.et920r05.d8m11y2021.root",
    "../../misc/ThesisAndPaperFiles/RootFilesForPlotsPP/Results/gammaVsPionSummary.onlyPy8_pTbinHuge.et920r05.d8m11y2021.root"
  };
  const TString sInSys[NSys] = {
    "../../misc/ThesisAndPaperFiles/RootFilesForPlotsPP/Results/gammaVsPionSummary.onlyPy8_pTbinHuge.et920r05.d8m11y2021.root",
    "../../misc/ThesisAndPaperFiles/RootFilesForPlotsPP/Results/gammaVsPionSummary.onlyPy8_pTbinHuge.et920r05.d8m11y2021.root",
    "../../misc/ThesisAndPaperFiles/RootFilesForPlotsPP/Results/gammaVsPionSummary.onlyPy8_pTbinHuge.et920r05.d8m11y2021.root"
  };
  const TString sGraphStat[NStat] = {
    "gStat911ga",
    "gStat1115ga",
    "gStat1520ga"
  };
  const TString sGraphSys[NSys] = {
    "gSys911ga",
    "gSys1115ga",
    "gSys1520ga"
  };

  // output parameters
  const TString sTxtStat[NStat] = {
    "correctedJets_statErr.et911r05gam.d8m7y2023.txt",
    "correctedJets_statErr.et1115r05gam.d8m7y2023.txt",
    "correctedJets_statErr.et1520r05gam.d8m7y2023.txt"
  };
  const TString sTxtSys[NSys] = {
    "correctedJets_sysErr.et911r05gam.d8m7y2023.txt",
    "correctedJets_sysErr.et1115r05gam.d8m7y2023.txt",
    "correctedJets_sysErr.et1520r05gam.d8m7y2023.txt"
  };

  // open input files
  TFile *fInStat[NStat];
  for (Ssiz_t iStat = 0 ; iStat < NStat; iStat++) {
    fInStat[iStat] = new TFile(sInStat[iStat].Data(),   "read");
    if (!fInStat[iStat]) {
      cerr << "PANIC: couldn't open input statistics file #" << iStat << "!\n" << endl;
      return;
    }
  }

  TFile *fInSys[NSys];
  for (Ssiz_t iSys = 0 ; iSys < NSys; iSys++) {
    fInSys[iSys] = new TFile(sInSys[iSys].Data(),  "read");
    if (!fInSys[iSys]) {
      cerr << "PANIC: couldn't open input systematics file #" << iSys << "!\n" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab graphs
  TGraphAsymmErrors *gStat[NStat];
  for (Ssiz_t iStat = 0; iStat < NStat; iStat++) {
    gStat[iStat] = (TGraphAsymmErrors*) fInStat[iStat] -> Get(sGraphStat[iStat].Data());
    if (!gStat[iStat]) {
      cerr << "PANIC: couldn't grab input statistic graph #" << iStat << "!\n" << endl;
      return;
    }
  }

  TGraphAsymmErrors *gSys[NSys];
  for (Ssiz_t iSys = 0; iSys < NSys; iSys++) {
    gSys[iSys] = (TGraphAsymmErrors*) fInSys[iSys] -> Get(sGraphSys[iSys].Data());
    if (!gSys[iSys]) {
      cerr << "PANIC: couldn't grab input systematics graph #" << iSys << "!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed graphs." << endl;

  // create statistics data tableis
  const UInt_t nStatBins = gStat[0] -> GetN();
  for (Ssiz_t iStat = 0; iStat < NStat; iStat++) {

    // open stream
    ofstream ofStat(sTxtStat[iStat].Data());
    ofStat << "x_value";
    ofStat << " ";
    ofStat << "y_value";
    ofStat << " ";
    ofStat << "x_error_down";
    ofStat << " ";
    ofStat << "x_error_up";
    ofStat << " ";
    ofStat << "y_error_down";
    ofStat << " ";
    ofStat << "y_error_up";
    ofStat << endl;

    // loop over bins
    for (UInt_t iBin = 1; iBin < (nStatBins + 1); iBin++) {

      // get bin values
      Double_t valX(0.);
      Double_t valY(0.);
      gStat[iStat] -> GetPoint(iBin, valX, valY);

      // get bin errors
      const Double_t errLoX = gStat[iStat] -> GetErrorXlow(iBin);
      const Double_t errHiX = gStat[iStat] -> GetErrorXhigh(iBin);
      const Double_t errLoY = gStat[iStat] -> GetErrorYlow(iBin);
      const Double_t errHiY = gStat[iStat] -> GetErrorYhigh(iBin);

      // print contents
      ofStat << valX;
      ofStat << " ";
      ofStat << valY;
      ofStat << " ";
      ofStat << errLoX;
      ofStat << " ";
      ofStat << errHiX;
      ofStat << " ";
      ofStat << errLoY;
      ofStat << " ";
      ofStat << errHiY;
      ofStat << endl;

    }  // end bin loop
  }  // end weight loop

  // create systematics data tableis
  const UInt_t nSysBins = gSys[0] -> GetN();
  for (Ssiz_t iSys = 0; iSys < NSys; iSys++) {

    // open stream
    ofstream ofSys(sTxtSys[iSys].Data());
    ofSys << "x_value";
    ofSys << " ";
    ofSys << "y_value";
    ofSys << " ";
    ofSys << "x_error_down";
    ofSys << " ";
    ofSys << "x_error_up";
    ofSys << " ";
    ofSys << "y_error_down";
    ofSys << " ";
    ofSys << "y_error_up";
    ofSys << endl;

    // loop over bins
    for (UInt_t iBin = 1; iBin < (nSysBins + 1); iBin++) {

      // get bin values
      Double_t valX(0.);
      Double_t valY(0.);
      gSys[iSys] -> GetPoint(iBin, valX, valY);

      // get bin errors
      const Double_t errLoX = gSys[iSys] -> GetErrorXlow(iBin);
      const Double_t errHiX = gSys[iSys] -> GetErrorXhigh(iBin);
      const Double_t errLoY = gSys[iSys] -> GetErrorYlow(iBin);
      const Double_t errHiY = gSys[iSys] -> GetErrorYhigh(iBin);

      // print contents
      ofSys << valX;
      ofSys << " ";
      ofSys << valY;
      ofSys << " ";
      ofSys << errLoX;
      ofSys << " ";
      ofSys << errHiX;
      ofSys << " ";
      ofSys << errLoY;
      ofSys << " ";
      ofSys << errHiY;
      ofSys << endl;

    }  // end bin loop
  }  // end weight loop
  cout << "    Created data tables." << endl;

  // close files
  for (Ssiz_t iStat = 0; iStat < NStat; iStat++) {
    fInStat[iStat] -> cd();
    fInStat[iStat] -> Close();
  }
  for (Ssiz_t iSys = 0; iSys < NSys; iSys++) {
    fInSys[iSys]  -> cd();
    fInSys[iSys]  -> Close();
  }
  cout << "  Data tables creation finished.\n" << endl;

}

// End ------------------------------------------------------------------------
