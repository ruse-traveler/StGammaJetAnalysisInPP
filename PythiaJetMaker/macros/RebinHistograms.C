// 'RebinHistograms.C'
// Derek Anderson
// 05.15.2017
//
// Rebin your histograms.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace std;


// input
static const TString sOut("rebin.pTtrk.n2d2.d15m5y2017.root");
static const TString sInA("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Embedding/Run12pp/PythiaMatching/output/merged/pp200r12.pythia.pt2_-1.root");
static const TString sInB("/global/project/projectdirs/star/pwg/starjetc/dmawxc/Embedding/Run12pp/JetMaker/mudst/cut.gam.d8m5y2017.root");
static const TString sHistA("hMcFinalPtS1_Inside");
static const TString sHistB("FitCut2/hPtTrkG_n2d2");
// output
static const TString sNewA("hMcParPt");
static const TString sNewB("hMcDetPt");



void RebinHistograms() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning rebinning script..." << endl;


  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInA = new TFile(sInA.Data(), "read");
  TFile *fInB = new TFile(sInB.Data(), "read");
  if (!fOut || !fInA || !fInB) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOut);
    assert(fInA);
    assert(fInB);
  }
  cout << "    Files opened." << endl;

  TH1D *hA = (TH1D*) fInA -> Get(sHistA.Data());
  TH1D *hB = (TH1D*) fInB -> Get(sHistB.Data());
  if (!hA || !hB) {
    cerr << "PANIC: couldn't grab a histogram!" << endl;
    assert(hA);
    assert(hB);
  }
  cout << "    Histograms grabbed." << endl;


  // determine binning
  Int_t    nXA = hA -> GetNbinsX();
  Int_t    nXB = hB -> GetNbinsX();
  Double_t xA1 = hA -> GetBinLowEdge(1);
  Double_t xB1 = hB -> GetBinLowEdge(1);
  Double_t xA2 = hA -> GetBinLowEdge(nXA + 1);
  Double_t xB2 = hB -> GetBinLowEdge(nXB + 1);
  Double_t dXA = (xA2 - xA1) / nXA;
  Double_t dXB = (xB2 - xB1) / nXB;
  cout << "    Input histogram parameters:\n"
       << "      A) nX = " << nXA << "; (x1, x2) = (" << xA1 << ", " << xA2 << "); dX = " << dXA << "\n"
       << "      B) nX = " << nXB << "; (x1, x2) = (" << xB1 << ", " << xB2 << "); dX = " << dXB
       << endl;

  Int_t   iRebin = 0;
  TString sRebin("");
  if (dXA > dXB) {
    iRebin = 1;
    sRebin = "B";
  }
  else if (dXA < dXB) {
    iRebin = 2;
    sRebin = "A";
  }
  else {
    cout << "    Neither histogram needs to be rebinned!" << endl;
  }

  const Double_t dXold  = TMath::Min(dXA, dXB);
  const Double_t dXnew  = TMath::Max(dXA, dXB);
  const Int_t    nGroup = (Int_t) (dXnew / dXold);
  cout << "    Rebinning histogram " << sRebin.Data() << ":\n"
       << "      dXnew = " << dXnew << ", nGroup = " << nGroup
       << endl;

  TH1D *hNew;
  if (iRebin == 1) {
    hNew = (TH1D*) hB -> Rebin(nGroup, sNewB.Data());
    hA -> SetName(sNewA.Data());
  }
  else if (iRebin == 2) {
    hNew = (TH1D*) hA -> Rebin(nGroup, sNewA.Data());
    hB -> SetName(sNewB.Data());
  }
  else {
    hA -> SetName(sNewA.Data());
    hB -> SetName(sNewB.Data());
  }


  // display new dimensions
  if (iRebin == 1) {
    nXA = hA   -> GetNbinsX();
    nXB = hNew -> GetNbinsX();
    xA1 = hA   -> GetBinLowEdge(1);
    xB1 = hNew -> GetBinLowEdge(1);
    xA2 = hA   -> GetBinLowEdge(nXA + 1);
    xB2 = hNew -> GetBinLowEdge(nXB + 1);
    dXA = (xA2 - xA1) / nXA;
    dXB = (xB2 - xB1) / nXB;
    cout << "    New histogram parameters:\n"
         << "      A) nX = " << nXA << "; (x1, x2) = (" << xA1 << ", " << xA2 << "); dX = " << dXA << "\n"
         << "      B) nX = " << nXB << "; (x1, x2) = (" << xB1 << ", " << xB2 << "); dX = " << dXB
         << endl;
  }
  else if (iRebin == 2) {
    nXA = hNew -> GetNbinsX();
    nXB = hB   -> GetNbinsX();
    xA1 = hNew -> GetBinLowEdge(1);
    xB1 = hB   -> GetBinLowEdge(1);
    xA2 = hNew -> GetBinLowEdge(nXA + 1);
    xB2 = hB   -> GetBinLowEdge(nXB + 1);
    dXA = (xA2 - xA1) / nXA;
    dXB = (xB2 - xB1) / nXB;
    cout << "    New histogram parameters:\n"
         << "      A) nX = " << nXA << "; (x1, x2) = (" << xA1 << ", " << xA2 << "); dX = " << dXA << "\n"
         << "      B) nX = " << nXB << "; (x1, x2) = (" << xB1 << ", " << xB2 << "); dX = " << dXB
         << endl;
  }


  // close files
  fOut -> cd();
  if (iRebin == 1) {
    hA   -> Write();
    hNew -> Write();
  }
  else if (iRebin == 2) {
    hNew -> Write();
    hB   -> Write();
  }
  else {
    hA -> Write();
    hB -> Write();
  }
  fOut -> Close();
  fInA -> cd();
  fInA -> Close();
  fInB -> cd();
  fInB -> Close();


  cout << "  Rebinning script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
