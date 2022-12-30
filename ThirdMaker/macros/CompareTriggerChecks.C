// 'CompareTriggerChecks.C'
// Derek Anderson
// 03.01.2017
//
// Short plotting macro to compare ouput of
// 'CheckTreeTriggers.C' against output of
// 'CheckTupleTriggers.C'.

#include <assert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const Int_t    nDec   = 2;
static const Int_t    nHist  = 3;
static const Int_t    nMaker = 2;
static const Int_t    font   = 42;
static const Int_t    center = 1;
static const Double_t label  = 0.03;
// file-paths
static const TString sMaker[nMaker] = {"treeTrigs.beforeFix.r10145070.d12m4y2017.root", "treeTrigs.afterFix.r10145070.d12m4y2017.root"};
static const TString sOut("afterXbefore.trigs.r10145070.d12m4y2017.root");
// hist. labels
static const TString sLeg[nMaker]  = {"Before fix", "After fix"};
static const TString sHist[nHist]  = {"hNumTrg", "hNumPass", "hNumRatio"};
static const TString sCanv[nHist]  = {"cNumTrg", "cNumPass", "cNumRatio"};
static const TString sHistT[nHist] = {"No of triggers per event", "No of triggers passing acceptance cuts", "Ratio of accepted triggers to total triggers"};
static const TString sHistX[nHist] = {"N_{trg}", "N_{pass}", "N_{pass} / N_{ratio}"};
static const TString sHistY[nHist] = {"counts", "counts", "counts"};
// hist. style
static const Int_t color[nMaker]  = {890, 810};
static const Int_t marker[nMaker] = {7, 4};
static const Int_t line[nMaker]   = {1, 2};
static const Int_t width[nMaker]  = {2, 1};



void CompareTriggerChecks() {

  cout << "\n  Beginning comparison script..." << endl;
  gErrorIgnoreLevel = kError;


  // open files and grab histograms
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fMaker[nMaker];
  for (Int_t m = 0; m < nMaker; m++) {
    fMaker[m] = new TFile(sMaker[m].Data(), "read");
    if (!fMaker[m]) {
      cerr << "PANIC: couldn't open input file " << m << endl;
      assert(fMaker[m]);
    }
  }

  TH1D *hHist[nHist][nMaker];
  for (Int_t h = 0; h < nHist; h++) {
    for (Int_t m = 0; m < nMaker; m++) {
      hHist[h][m] = (TH1D*) fMaker[m] -> Get(sHist[h].Data());
      if (!hHist[h][m]) {
        cerr << "PANIC: couldn't grab histogram " << h << " in file " << m << endl;
        assert(hHist[h][m]);
      }
    }
  }
  cout << "    Histograms grabbed..." << endl;


  // set histogram styles
  for (Int_t h = 0; h < nHist; h++) {
    for (Int_t m = 0; m < nMaker; m++) {
      hHist[h][m] -> SetTitle(sHistT[h].Data());
      hHist[h][m] -> SetTitleFont(font);
      hHist[h][m] -> SetMarkerColor(color[m]);
      hHist[h][m] -> SetMarkerStyle(marker[m]);
      hHist[h][m] -> SetLineColor(color[m]);
      hHist[h][m] -> SetLineStyle(line[m]);
      hHist[h][m] -> SetLineWidth(width[m]);
      hHist[h][m] -> GetXaxis() -> SetTitle(sHistX[h].Data());
      hHist[h][m] -> GetXaxis() -> SetTitleFont(font);
      hHist[h][m] -> GetXaxis() -> CenterTitle(center);
      hHist[h][m] -> GetXaxis() -> SetLabelSize(label);
      hHist[h][m] -> GetYaxis() -> SetTitle(sHistY[h].Data());
      hHist[h][m] -> GetYaxis() -> SetTitleFont(font);
      hHist[h][m] -> GetYaxis() -> CenterTitle(center);
      hHist[h][m] -> GetYaxis() -> SetLabelSize(label);
    }
  }
  cout << "    Styles set..." << endl;


  // text box styles
  const Int_t colorB = 0;
  const Int_t colorT = 1;
  const Int_t widthB = 1;
  const Int_t fillB  = 0;

  // calculate integrals
  Int_t     nDig(0);
  Int_t     size(0);
  TString   raw("");
  TString   sInt[nHist][nMaker];
  TString   sTxt[nHist][nMaker];
  Double_t  ints[nHist][nMaker];
  TPaveText *pHist[nHist];
  for (Int_t h = 0; h < nHist; h++) {
    for (Int_t m = 0; m < nMaker; m++) {
      ints[h][m] = hHist[h][m] -> Integral();
      raw  = "";
      raw += ints[h][m];
      nDig = raw.First(".");
      if (nDig > 0) {
        size = (nDig + 1) + nDec;
      }
      else {
        size = raw.Length();
      }
      sInt[h][m]  = "";
      sInt[h][m].Append(raw, size);
      sTxt[h][m]  = "#color[";
      sTxt[h][m] += color[m];
      sTxt[h][m] += "]{I = ";
      sTxt[h][m] += sInt[h][m];
      sTxt[h][m] += "}";
    }
    pHist[h] = new TPaveText(0.1, 0.1, 0.3, 0.3, "NDC NB");
    pHist[h] -> SetLineColor(colorB);
    pHist[h] -> SetLineWidth(widthB);
    pHist[h] -> SetFillColor(colorB);
    pHist[h] -> SetFillStyle(fillB);
    pHist[h] -> SetTextColor(colorT);
    pHist[h] -> SetTextFont(font);
    for (Int_t m = 0; m < nMaker; m++) {
      pHist[h] -> AddText(sTxt[h][m]);
    }
  }
  cout << "    Integrals calculated..." << endl;

  // create legends
  TLegend *lHist[nHist];
  for (Int_t h = 0; h < nHist; h++) {
    lHist[h] = new TLegend(0.3, 0.1, 0.5, 0.3);
    lHist[h] -> SetLineColor(colorB);
    lHist[h] -> SetLineWidth(widthB);
    lHist[h] -> SetFillColor(colorB);
    lHist[h] -> SetFillStyle(fillB);
    lHist[h] -> SetTextColor(colorT);
    lHist[h] -> SetTextFont(font);
    for (Int_t m = 0; m < nMaker; m++) {
      lHist[h] -> AddEntry(hHist[h][m], sLeg[m].Data());
    }
  }
  cout << "    Legends created..." << endl;


  // canvas size
  const Int_t grid    = 0;
  const Int_t logY    = 1;
  const Int_t widthC  = 800;
  const Int_t heightC = 800;

  // draw histograms
  fOut -> cd();
  TCanvas *cHist[nHist];
  for (Int_t h = 0; h < nHist; h++) {
    cHist[h] = new TCanvas(sCanv[h].Data(), "", widthC, heightC);
    cHist[h] -> SetGrid(grid, grid);
    cHist[h] -> SetLogy(logY);
    for (Int_t m = 0; m < nMaker; m++) {
      if (m == 0)
        hHist[h][m] -> Draw();
      else
        hHist[h][m] -> Draw("same");
    }
    pHist[h] -> Draw();
    lHist[h] -> Draw();
    cHist[h] -> Write();
    cHist[h] -> Close(); 
  }
  cout << "    Canvases painted..." << endl;


  // close files
  fOut -> cd();
  fOut -> Close();
  for (Int_t m = 0; m < nMaker; m++) {
    fMaker[m] -> cd();
    fMaker[m] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
