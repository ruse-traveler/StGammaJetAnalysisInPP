// 'MergeHistograms.C'
// Derek Anderson
//
// Use this to merge and plot the debugging
// histograms from the ThirdJetMaker.


#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// constants
static const Int_t nHist  = 19;
static const Int_t nClass = 2;
// i/o parameters
static const TString sFiles("hist2.1fpp.r10120014-10124011.d20m4y2017.list");
static const TString sOut("hist.1fpp.r10120014-10124011.d20m4y2017.root");
static const TString sEvt[nHist] = {"hEvtID_s", "hRunID_s", "hVx_s", "hVy_s", "hVz_s", "hMag_s", "hZdcCo_s", "hBbcCo_s", "hBkgdR_s", "hBlueBR_s", "hYellBR_s", "hBbcZ_s", "hVpdZ_s", "hPMult_s", "hNMult_s", "hVtxR_s", "hNvtx_s", "hNprim_s", "hNglob_s"};
static const TString sGen[nHist] = {"hEvtID_g", "hRunID_g", "hVx_g", "hVy_g", "hVz_g", "hMag_g", "hZdcCo_g", "hBbcCo_g", "hBkgdR_g", "hBlueBR_g", "hYellBR_g", "hBbcZ_g", "hVpdZ_g", "hPMult_g", "hNMult_g", "hVtxR_g", "hNvtx_g", "hNprim_g", "hNglob_g"};
static const TString sCan[nHist] = {"cEvtID", "cRunID", "cVx", "cVy", "cVz", "cMag", "cZdcCo", "cBbcCo", "cBkgdR", "cBlueBR", "cYellBR", "cBbcZ", "cVpdZ", "cPMult", "cNMult", "cVtxR", "cNvtx", "cNprim", "cNglob"};
// histogram parameters
static const Int_t    iColor[nClass]  = {890, 810};
static const Int_t    iFill[nClass]   = {3004, 3005};
static const Int_t    iLine[nClass]   = {1, 1};
static const Int_t    iWidth[nClass]  = {1, 1};
static const Int_t    iText[nClass]   = {42, 42};
static const Int_t    iCent[nClass]   = {1, 1};
static const TString  sTitle[nHist]   = {"Event ID", "Run ID", "Primary V_{x}", "Primary V_{y}", "Primary V_{z}", "Magnetic Field", "ZDC Coincidence", "BBC Coincidence", "Background Rate", "Blue Background Rate", "Yellow Background Rate", "BBC Z", "VPD Z", "Positive Multiplicity", "Negative Multiplicity", "Vertex Rank", "Number of Vertices", "Number of Primary Tracks", "Number of Global Tracks"};
static const TString  sTitleX[nHist]  = {"#DeltaEvtID", "#DeltaRunID", "#DeltaV_{x}", "#DeltaV_{y}", "#DeltaV_{z}", "#DeltaB", "#DeltaZDC", "#DeltaBBC", "#DeltaR_{bkgd}", "#DeltaR_{blue}", "#DeltaR_{yellow}", "#DeltaZ_{bbc}", "#DeltaZ_{vpd}", "#DeltaM_{pos}", "#DeltaM_{neg}", "#DeltaK_{vtx}", "#DeltaN_{vtx}", "#DeltaN_{prim}", "#DeltaN_{glob}"};
static const TString  sTitleY[nHist]  = {"counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts", "counts"};
static const TString  sLegend[nClass] = {"StMuEvent - StEvent", "StMuEvent - Generated StEvent"};
static const Double_t dSize[nClass]   = {0.03, 0.03};



void MergeHistograms(const Bool_t inBatchMode=false) {

  cout << "\n  Beginning merge script..." << endl;
  gErrorIgnoreLevel = kFatal;


  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TH1D  *hist[nHist][nClass];
  // open stream
  ifstream files(sFiles.Data());
  if (!files) {
    cerr << "PANIC: couldn't open input stream!" << endl;
    assert(files);
  }


  // stream files, merge histograms
  cout << "    Streaming files..." << endl;

  Int_t  nF(0);
  string in("");
  TFile  *fIn;
  while (getline(files, in)) {
    // open files
    TString sIn(in);
    fIn = new TFile(sIn.Data(), "read");
    if (!fIn) {
      cerr << "PANIC: couldn't open input file\n"
           << "       '" << sIn.Data() << "'"
           << endl;
      assert(fIn);
    }
    // merge histograms
    if (!inBatchMode) {
      cout << "      Processing file " << nF << "...\r" << flush;
    }
    else {
      cout << "      Processing file " << nF << "..." << endl;
    }
    for (Int_t h = 0; h < nHist; h++) {
      if (nF == 0) {
        hist[h][0] = (TH1D*) fIn -> Get(sEvt[h].Data()) -> Clone();
        hist[h][1] = (TH1D*) fIn -> Get(sGen[h].Data()) -> Clone();
        hist[h][0] -> SetDirectory(0);
        hist[h][1] -> SetDirectory(0);
      }
      else {
        hist[h][0] -> Add((TH1D*) fIn -> Get(sEvt[h].Data()));
        hist[h][1] -> Add((TH1D*) fIn -> Get(sGen[h].Data()));
      }
      // check histogram
      if (!hist[h][0]) {
        cerr << "PANIC: couldn't grab StEvent histogram " << h
             << endl;
        assert(hist[h][0]);
      }
      if (!hist[h][1]) {
        cerr << "PANIC: couldn't grab generated histogram " << h
             << endl;
        assert(hist[h][1]);
      }
    }
    fIn -> Close();
    nF++;
  }  // end stream

  if (!inBatchMode) {
    cout << endl;
  }
  files.close();


  // set styles
  cout << "    Setting styles..." << endl;

  for (Int_t h = 0; h < nHist; h++) {
    for (Int_t c = 0; c < nClass; c++) {
      hist[h][c] -> SetTitle(sTitle[h].Data());
      hist[h][c] -> SetTitleFont(iText[c]);
      hist[h][c] -> SetLineColor(iColor[c]);
      hist[h][c] -> SetLineStyle(iLine[c]);
      hist[h][c] -> SetLineWidth(iWidth[c]);
      hist[h][c] -> SetFillColor(iColor[c]);
      hist[h][c] -> SetFillStyle(iFill[c]);
      hist[h][c] -> GetXaxis() -> SetTitle(sTitleX[h].Data());
      hist[h][c] -> GetXaxis() -> SetTitleFont(iText[c]);
      hist[h][c] -> GetXaxis() -> CenterTitle(iCent[c]);
      hist[h][c] -> GetXaxis() -> SetLabelSize(dSize[c]);
      hist[h][c] -> GetYaxis() -> SetTitle(sTitleY[h].Data());
      hist[h][c] -> GetYaxis() -> SetTitleFont(iText[c]);
      hist[h][c] -> GetYaxis() -> CenterTitle(iCent[c]);
      hist[h][c] -> GetYaxis() -> SetLabelSize(dSize[c]);
    }
  }


  // create text boxes
  cout << "    Creating text boxes..." << endl;

  TPaveText *pHist[nHist][nClass];
  for (Int_t h = 0; h < nHist; h++) {
    for (Int_t c = 0; c < nClass; c++) {
      pHist[h][c] = new TPaveText(0.1, 0.1, 0.3, 0.3, "NDC NB");
      pHist[h][c] -> SetLineColor(0);
      pHist[h][c] -> SetLineWidth(1);
      pHist[h][c] -> SetFillColor(0);
      pHist[h][c] -> SetFillStyle(0);
      pHist[h][c] -> SetTextColor(1);
      pHist[h][c] -> SetTextFont(42);
      pHist[h][c] -> SetTextAlign(22);

      Double_t mean = hist[h][c] -> GetMean();
      Double_t rms  = hist[h][c] -> GetRMS();
      TString  sMean("");
      TString  sRMS("");
      sMean += "#color[";
      sMean += iColor[c];
      sMean += "]{mean = ";
      sMean += mean;
      sMean += "}";
      sRMS  += "#color[";
      sRMS  += iColor[c];
      sRMS  += "]{RMS = ";
      sRMS  += rms;
      sRMS  += "}";
      pHist[h][c] -> AddText(sMean.Data());
      pHist[h][c] -> AddText(sRMS.Data());
    }
  }

  TLegend *lHist[nHist];
  for (Int_t h = 0; h < nHist; h++) {
    lHist[h] = new TLegend(0.5, 0.1, 0.7, 0.3);
    lHist[h] -> SetLineColor(0);
    lHist[h] -> SetLineWidth(1);
    lHist[h] -> SetFillColor(0);
    lHist[h] -> SetFillStyle(0);
    lHist[h] -> SetTextColor(1);
    lHist[h] -> SetTextFont(42);
    lHist[h] -> SetTextAlign(22);
    for (Int_t c = 0; c < nClass; c++) {
      lHist[h] -> AddEntry(hist[h][c], sLegend[c].Data());
    }
  }


  // create canvases
  cout << "    Drawing plots..." << endl;
  fOut -> cd();

  TCanvas *cHist[nHist];
  for (Int_t h = 0; h < nHist; h++) {
    cHist[h] = new TCanvas(sCan[h].Data(), "", 800, 800);
    cHist[h] -> SetGrid(0, 0);
    cHist[h] -> SetLogy(1);
    for (Int_t c = 0; c < nClass; c++) {
      if (c == 0) {
        hist[h][c] -> Draw("");
      }
      else {
        hist[h][c] -> Draw("same");
      }
      pHist[h][c] -> Draw();
    }
    lHist[h] -> Draw();
    cHist[h] -> Write();
    cHist[h] -> Close();
  }

  // save and close
  cout << "    Closing file..." << endl;

  fOut -> cd();
  for (Int_t h = 0; h < nHist; h++) {
    hist[h][0] -> Write();
    hist[h][1] -> Write();
  }
  fOut -> Close();


  cout << "  Merging script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
