// 'CreateDataTablesForPaper.C'
// Derek Anderson
// 05.28.2023
// 
// Reads in rootfiles for fig.s 10 and 11 in the
// gamma+jet PRC and spits out text files with
// their data in table and saves the canvases
// to their own files.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

const Ssiz_t NMatrix(2);
const Ssiz_t NWeights(5);


void CreateDataTablesForPaper() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Begining data tables creation script..." << endl;

  // io parameters
  const TString sInMatrix("trigMatrixForLongPaper_betterMargins_withAprGpcComments.et630pi0vsGam.d20m4y2023.root");
  const TString sInWeights("smearingWeightsForPaper_betterTitles_withAprGpcComments.et630pi0vsGam.d20m4y2023.root");
  const TString sOutMatrix("trigMatrixForLongPaper_onlyCanvas.d28m5y2023.root");
  const TString sOutWeights("smearingWeightsForLongPaper_onlyCanvas.d28m5y2023.root");

  // output text parameters
  const TString sTxtMatrixValues[NMatrix] = {
    "trigMatrixForLongPaper_values_dataTable.et630gam.d28m5y2023.txt",
    "trigMatrixForLongPaper_values_dataTable.et630pi0.d28m5y2023.txt"
  };
  const TString sTxtMatrixErrors[NMatrix] = {
    "trigMatrixForLongPaper_errors_dataTable.et630gam.d28m5y2023.txt",
    "trigMatrixForLongPaper_errors_dataTable.et630pi0.d28m5y2023.txt"
  };
  const TString sTxtWeights[NWeights] = {
    "smearingWeightsForLongPaper_dataTable.et911gam.d28m5y2023.txt",
    "smearingWeightsForLongPaper_dataTable.et1115gam.d28m5y2023.txt",
    "smearingWeightsForLongPaper_dataTable.et1520gam.d28m5y2023.txt",
    "smearingWeightsForLongPaper_dataTable.et911pi0.d28m5y2023.txt",
    "smearingWeightsForLongPaper_dataTable.et1115pi0.d28m5y2023.txt"
  };

  // input canvas/hist parameters
  const TString sCanvasMatrix("cMatrix");
  const TString sCanvasWeights("cWeights");
  const TString sHistMatrix[NMatrix] = {
    "hMatrixGam",
    "hMatrixPi0"
  };
  const TString sHistWeights[NWeights] = {
    "hGamWeightsEt911",
    "hGamWeightsEt1115",
    "hGamWeightsEt1520",
    "hPi0WeightsEt911",
    "hPi0WeightsEt1115"
  };

  // open files
  TFile *fInMatrix   = new TFile(sInMatrix.Data(),   "read");
  TFile *fInWeights  = new TFile(sInWeights.Data(),  "read");
  TFile *fOutMatrix  = new TFile(sOutMatrix.Data(),  "recreate");
  TFile *fOutWeights = new TFile(sOutWeights.Data(), "recreate");
  if (!fInMatrix || !fInWeights || !fOutMatrix || !fOutWeights) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fInMatrix  = " << fInMatrix  << ", fInWeights = "  << fInWeights  << "\n"
         << "       fOutMatrix = " << fOutMatrix << ", fOutWeights = " << fOutWeights << "\n"
         << endl;
    return;
  }
  cout << "    Opened histograms." << endl;

  // grab canvases
  TCanvas *cMatrix  = (TCanvas*) fInMatrix  -> Get(sCanvasMatrix.Data())  -> Clone();
  TCanvas *cWeights = (TCanvas*) fInWeights -> Get(sCanvasWeights.Data()) -> Clone();
  if (!cMatrix || !cWeights) {
    cerr << "PANIC: couldn't grab a canvas!\n"
         << "       cMatrix = " << cMatrix << ", cWeights = " << cWeights << "\n"
         << endl;
    return;
  }
  cout << "    Grabbed canvases." << endl;

  // grab histograms
  TH2D *hMatrix[NMatrix];
  for (Ssiz_t iMatrix = 0; iMatrix < NMatrix; iMatrix++) {
    hMatrix[iMatrix] = (TH2D*) fInMatrix -> Get(sHistMatrix[iMatrix].Data());
    if (!hMatrix[iMatrix]) {
      cerr << "PANIC: couldn't grab input matrix #" << iMatrix << "!\n" << endl;
      return;
    }
  }

  TH1D *hWeights[NWeights];
  for (Ssiz_t iWeights = 0; iWeights < NWeights; iWeights++) {
    hWeights[iWeights] = (TH1D*) fInWeights -> Get(sHistWeights[iWeights].Data());
    if (!hWeights[iWeights]) {
      cerr << "PANIC: couldn't grab input matrix #" << iWeights << "!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // create matrix data tables
  const UInt_t nMatrixBinsX = hMatrix[0] -> GetNbinsX();
  const UInt_t nMatrixBinsY = hMatrix[0] -> GetNbinsY();
  for (Ssiz_t iMatrix = 0; iMatrix < NMatrix; iMatrix++) {

    // open streams
    ofstream ofMatrixValues(sTxtMatrixValues[iMatrix].Data());
    ofstream ofMatrixErrors(sTxtMatrixErrors[iMatrix].Data());

    // loop over bins
    for (UInt_t iBinY = 0; iBinY < (nMatrixBinsY + 1); iBinY++) {
      for (UInt_t iBinX = 0; iBinX < (nMatrixBinsX + 1); iBinX++) {

        // determine location
        const Bool_t isBeforeBinsY = (iBinY == 0);
        const Bool_t isBeforeBinsX = (iBinX == 0);
        const Bool_t isAtLastBinX  = (iBinX == nMatrixBinsX);
        const Bool_t isInBinsY     = (iBinY > 0);
        const Bool_t isInBinsX     = (iBinX > 0);

        // put nothing in top left corner
        if (isBeforeBinsY && isBeforeBinsX) {
          ofMatrixValues << " ";
          ofMatrixErrors << " ";
        }

        // get bin center and content
        const Double_t xCenter = hMatrix[iMatrix] -> GetXaxis() -> GetBinCenter(iBinX);
        const Double_t yCenter = hMatrix[iMatrix] -> GetYaxis() -> GetBinCenter(iBinY);
        const Double_t value   = hMatrix[iMatrix] -> GetBinContent(iBinX, iBinY);
        const Double_t error   = hMatrix[iMatrix] -> GetBinError(iBinX, iBinY);

        // print axes
        if (isBeforeBinsY && isInBinsX) {
          ofMatrixValues << xCenter;
          ofMatrixValues << " ";
          ofMatrixErrors << xCenter;
          ofMatrixErrors << " ";
        }
        if (isBeforeBinsX && isInBinsY) {
          ofMatrixValues << yCenter;
          ofMatrixValues << " ";
          ofMatrixErrors << yCenter;
          ofMatrixErrors << " ";
        }

        // print value and error
        if (isInBinsY && isInBinsX) {
          ofMatrixValues << value;
          ofMatrixValues << " ";
          ofMatrixErrors << error;
          ofMatrixErrors << " ";
        }

        // print out endlines
        if (isAtLastBinX) {
          ofMatrixValues << endl;
          ofMatrixErrors << endl;
        }
      }  // end x bin loop
    }  // end y bin loop
  }  // end matrix loop
  cout << "    Created matrix data tables." << endl;

  // create weight data tables
  const UInt_t nWeightsBins = hWeights[0] -> GetNbinsX();
  for (Ssiz_t iWeights = 0; iWeights < NWeights; iWeights++) {

    // open stream
    ofstream ofWeights(sTxtWeights[iWeights].Data());
    ofWeights << "bin_center";
    ofWeights << " ";
    ofWeights << "bin_value";
    ofWeights << " ";
    ofWeights << "bin_error";
    ofWeights << endl;

    // loop over bins
    for (UInt_t iBin = 1; iBin < (nWeightsBins + 1); iBin++) {

      // get bin center and content
      const Double_t center = hWeights[iWeights] -> GetBinCenter(iBin);
      const Double_t value  = hWeights[iWeights] -> GetBinContent(iBin);
      const Double_t error  = hWeights[iWeights] -> GetBinError(iBin);

      // print contents
      ofWeights << center;
      ofWeights << " ";
      ofWeights << value;
      ofWeights << " ";
      ofWeights << error;
      ofWeights << endl;

    }  // end bin loop
  }  // end weight loop
  cout << "    Created matrix data tables." << endl;

  // save canvases to file
  fOutMatrix  -> cd();
  cMatrix     -> Write();
  fOutWeights -> cd();
  cWeights    -> Write();
  cout << "    Saved canvases." << endl;

  // close files
  fOutMatrix  -> cd();
  fOutMatrix  -> Close();
  fOutWeights -> cd();
  fOutWeights -> Close();
  fInMatrix   -> cd();
  fInMatrix   -> Close();
  fInWeights  -> cd();
  fInWeights  -> Close();
  cout << "  Data tables creation finished.\n" << endl;

}

// End ------------------------------------------------------------------------
