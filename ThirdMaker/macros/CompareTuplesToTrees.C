// 'CompareTuplesToTrees.C'
// Derek Anderson
// 02.28.2017
//
// Use this to compare the output of 'ReadTuples.C'
// to the output of 'ReadFemtoDst.C'.

#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const TString sPythia("/star/u/dmawxc/JetRecon_pp/FullJetTreeTest/input/QuickCheck.r03a02rm1chrg.root");
static const TString sTuple("tuple.d4m5y2017.root");
static const TString sTree("tree.d4m5y2017.root");
static const TString sOut("treeXtuple.d4m5y2017.root");
// global constants
static const Int_t  nDec   = 2;
static const Int_t  nMaker = 2;
static const Int_t  nHist  = 7;
static const Int_t  nProf  = 5;
static const Int_t  nPyth  = 3;
static const Int_t  nEvt   = 3;
static const Int_t  iPyth  = nEvt + 1;
static const Bool_t pythia = true;
// hist / prof names
static const TString sHist[nHist] = {"hPrimeVz", "hEtTrg", "hFtrg", "hHtrg", "hPtTrk0", "hDfTrk0", "hHtrk0"};
static const TString sProf[nProf] = {"hVzVsRun", "hVzVsEvt", "hVzVsInd", "hRunVsInd", "hEvtVsInd"};
static const TString sPyth[nPyth] = {"QA/hPtTrk", "QA/hDfTrk", "QA/hHtrk"};



void CompareTuplesToTrees() {

  cout << "\n  Beginning comparison script..." << endl;
  gErrorIgnoreLevel = kError;


  TFile *fOut   = new TFile(sOut.Data(), "recreate");
  TFile *fTree  = new TFile(sTree.Data(), "read");
  TFile *fTuple = new TFile(sTuple.Data(), "read");
  if (!fTree || !fTuple) {
    cerr << "PANIC: couldn't open an input file!" << endl;
    assert(fTree && fTuple);
  }

  TFile *fPyth;
  if (pythia) {
    fPyth = new TFile(sPythia.Data(), "read");
    if (!fPyth) {
      cerr << "PANIC: couldn't open pythia file!" << endl;
      assert(fPyth);
    }
  } 

  // grab histograms / profiles
  TH1D *hist[nHist][nMaker];
  TH1D *prof[nProf][nMaker];
  for (Int_t i = 0; i < nHist; i++) {
    hist[i][0] = (TH1D*) fTree  -> Get(sHist[i].Data());
    hist[i][1] = (TH1D*) fTuple -> Get(sHist[i].Data());
    if (!hist[i][0]) {
      cerr << "PANIC: couldn't grab tree histogram no. " << i << "!" << endl;
      assert(hist[i][0]);
    }
    if (!hist[i][1]) {
      cerr << "PANIC: couldn't grab tuple histogram no. " << i << "!" << endl;
      assert(hist[i][1]);
    }
  }
  for (Int_t i = 0; i < nProf; i++) {
    prof[i][0] = (TH1D*) fTree  -> Get(sProf[i].Data());
    prof[i][1] = (TH1D*) fTuple -> Get(sProf[i].Data());
    if (!prof[i][0]) {
      cerr << "PANIC: couldn't grab tree profile no. " << i << "!" << endl;
      assert(prof[i][0]);
    }
    if (!prof[i][1]) {
      cerr << "PANIC: couldn't grab tuple profile no. " << i << "!" << endl;
      assert(prof[i][1]);
    }
  }

  TH1D *pyth[nPyth];
  if (pythia) {
    for (Int_t i = 0; i < nPyth; i++) {
      pyth[i] = (TH1D*) fPyth -> Get(sPyth[i].Data());
      if (!pyth[i]) {
        cerr << "PANIC: couldn't grab pythia histogram no. " << i << "!" << endl;
        assert(pyth[i]);
      }
    }
  }

  cout << "    Histograms and profiles grabbed..." << endl;


  // compute differences (from profiles) and create lines
  TLine *zero[nProf];
  TH1D  *diff[nProf];
  for (Int_t i = 0; i < nProf; i++) {
    TString diffName("hDiff");
    diffName += i;

    const Int_t    nX = prof[i][0] -> GetNbinsX();
    const Double_t x1 = prof[i][0] -> GetBinLowEdge(1);
    const Double_t x2 = prof[i][0] -> GetBinLowEdge(nX + 1);
    zero[i] = new TLine(x1, 0., x2, 0.);
    diff[i] = new TH1D(diffName.Data(), "", nX, x1, x2);
    diff[i] -> Sumw2();

    for (Int_t j = 1; j < nX + 1; j++) {
      const Double_t yTree  = prof[i][0] -> GetBinContent(j);
      const Double_t yTuple = prof[i][1] -> GetBinContent(j);
      const Double_t yDiff  = yTree - yTuple;
      diff[i] -> SetBinContent(j, yDiff);
      diff[i] -> SetBinError(j, 0);
    } 
  }

  cout << "    Differences computed..." << endl;


  // axis labels
  const TString sHistT[nHist] = {"Primary Vertex z-Coordinate", "Trigger E_{T}", "Trigger #varphi", "Trigger #eta", "Track p_{T}", "Track #Delta#varphi", "Track #eta"};
  const TString sHistX[nHist] = {"v_{z}", "E_{T}^{trg}", "#varphi_{trg}", "#eta_{trg}", "p_{T}^{trk}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#eta_{trk}"};
  const TString sHistY[nHist] = {"counts", "counts", "counts", "counts", "(1/N_{trg})dN_{trk}/dp_{T}^{trk}", "(1/N_{trg})dN_{trk}/d#Delta#varphi", "(1/N_{trg})dN_{trk}/d#eta_{trk}"};
  const TString sProfT[nProf] = {"v_{z} vs. run ID", "v_{z} vs. event ID", "v_{z} vs. tree-index", "run ID vs. tree-index", "event ID vs. tree-index"};
  const TString sProfX[nProf] = {"run", "event", "index", "index", "index"};
  const TString sProfY[nProf] = {"v_{z}", "v_{z}", "v_{z}", "run", "event"};
  const TString sPythX[nPyth] = {"p_{T}^{trk}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#eta_{trk}"};
  const TString sPythY[nPyth] = {"(1/N_{trg})dN_{trk}/dp_{T}^{trk}", "(1/N_{trg})dN_{trk}/d#Delta#varphi", "(1/N_{trk})dN_{trk}/d#eta_{trk}"};
  const TString sDiff("Tuple - Tree");

  // set styles
  const Int_t font               = 42;
  const Int_t center             = 1;
  const Int_t color[nMaker + 1]  = {860, 810, 850};
  const Int_t marker[nMaker + 1] = {7, 4, 1};
  const Int_t line[nMaker + 1]   = {1, 2, 1};
  const Int_t width[nMaker + 1]  = {2, 1, 2};
  for (Int_t i = 0; i < nHist; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      hist[i][j] -> SetTitle(sHistT[i].Data());
      hist[i][j] -> SetMarkerColor(color[j]);
      hist[i][j] -> SetMarkerStyle(marker[j]);
      hist[i][j] -> SetLineColor(color[j]);
      hist[i][j] -> SetLineStyle(line[j]);
      hist[i][j] -> SetLineWidth(width[j]);
      hist[i][j] -> GetXaxis() -> SetTitle(sHistX[i].Data());
      hist[i][j] -> GetXaxis() -> SetTitleFont(font);
      hist[i][j] -> GetXaxis() -> SetLabelFont(font);
      hist[i][j] -> GetXaxis() -> CenterTitle(center);
      hist[i][j] -> GetYaxis() -> SetTitle(sHistY[i].Data());
      hist[i][j] -> GetYaxis() -> SetTitleFont(font);
      hist[i][j] -> GetYaxis() -> SetLabelFont(font);
      hist[i][j] -> GetYaxis() -> CenterTitle(center);
    }
  }
  for (Int_t i = 0; i < nProf; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      prof[i][j] -> SetTitle(sProfT[i].Data());
      prof[i][j] -> SetMarkerColor(color[j]);
      prof[i][j] -> SetMarkerStyle(marker[j]);
      prof[i][j] -> SetLineColor(color[j]);
      prof[i][j] -> SetLineStyle(line[j]);
      prof[i][j] -> SetLineWidth(width[j]);
      prof[i][j] -> GetXaxis() -> SetTitle(sProfX[i].Data());
      prof[i][j] -> GetXaxis() -> SetTitleFont(font);
      prof[i][j] -> GetXaxis() -> SetLabelFont(font);
      prof[i][j] -> GetXaxis() -> CenterTitle(center);
      prof[i][j] -> GetYaxis() -> SetTitle(sProfY[i].Data());
      prof[i][j] -> GetYaxis() -> SetTitleFont(font);
      prof[i][j] -> GetYaxis() -> SetLabelFont(font);
      prof[i][j] -> GetYaxis() -> CenterTitle(center);
    }
  }
  for (Int_t i = 0; i < nProf; i++) {
    diff[i] -> SetTitle("");
    diff[i] -> SetMarkerColor(color[2]);
    diff[i] -> SetMarkerStyle(marker[2]);
    diff[i] -> SetLineColor(color[2]);
    diff[i] -> SetLineStyle(line[2]);
    diff[i] -> SetLineWidth(width[2]);
    diff[i] -> GetXaxis() -> SetTitle(sProfX[i]);
    diff[i] -> GetXaxis() -> SetTitleFont(font);
    diff[i] -> GetXaxis() -> SetLabelFont(font);
    diff[i] -> GetXaxis() -> CenterTitle(center);
    diff[i] -> GetYaxis() -> SetTitle(sDiff);
    diff[i] -> GetYaxis() -> SetTitleFont(font);
    diff[i] -> GetYaxis() -> SetLabelFont(font);
    diff[i] -> GetYaxis() -> CenterTitle(center);
  }

  // set pythia styles
  const Int_t colorP  = 850;
  const Int_t markerP = 1;
  const Int_t lineP   = 2;
  const Int_t widthP  = 1;
  if (pythia) {
    for (Int_t i = 0; i < nPyth; i++) {
      pyth[i] -> SetTitle("");
      pyth[i] -> SetMarkerColor(color[j]);
      pyth[i] -> SetMarkerStyle(marker[j]);
      pyth[i] -> SetLineColor(color[j]);
      pyth[i] -> SetLineStyle(line[j]);
      pyth[i] -> SetLineWidth(width[j]);
      pyth[i] -> GetXaxis() -> SetTitle(sPythX[i]);
      pyth[i] -> GetXaxis() -> SetTitleFont(font);
      pyth[i] -> GetXaxis() -> SetLabelFont(font);
      pyth[i] -> GetXaxis() -> CenterTitle(center);
      pyth[i] -> GetYaxis() -> SetTitle(sPythY[i]);
      pyth[i] -> GetYaxis() -> SetTitleFont(font);
      pyth[i] -> GetYaxis() -> SetLabelFont(font);
      pyth[i] -> GetYaxis() -> CenterTitle(center);
    }
  }

  // set line styles
  const Int_t colorL = 910;
  const Int_t widthL = 1;
  const Int_t styleL = 2;
  for (Int_t i = 0; i < nProf; i++) {
    zero[i] -> SetLineColor(colorL);
    zero[i] -> SetLineStyle(styleL);
    zero[i] -> SetLineWidth(widthL);
  }

  cout << "    Styles set..." << endl;


  // compute integrals
  Double_t hIntegrals[nHist][nMaker];
  Double_t pIntegrals[nProf][nMaker];
  for (Int_t i = 0; i < nHist; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      hIntegrals[i][j] = hist[i][j] -> Integral();
    }
  }
  for (Int_t i = 0; i < nProf; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      pIntegrals[i][j] = prof[i][j] -> Integral();
    }
  }

  // display only nDec decimals
  Int_t   nDig(0);
  Int_t   size(0);
  TString raw("");
  TString sIntH[nHist][nMaker];
  TString sIntP[nProf][nMaker];
  for (Int_t i = 0; i < nHist; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      raw   = "";
      raw  += hIntegrals[i][j];
      nDig  = raw.First(".");
      if (nDig > 0)
        size = (nDig + 1) + nDec;
      else
        size = raw.Length();
      sIntH[i][j] = "";
      sIntH[i][j].Append(raw, size);
    }
  }
  for (Int_t i = 0; i < nProf; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      raw   = "";
      raw  += pIntegrals[i][j];
      nDig  = raw.First(".");
      if (nDig > 0)
        size = (nDig + 1) + nDec;
      else
        size = raw.Length();
      sIntP[i][j] = "";
      sIntP[i][j].Append(raw, size);
    }
  }

  // set up text
  TString sTxtH[nHist][nMaker];
  TString sTxtP[nProf][nMaker];
  for (Int_t i = 0; i < nHist; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      sTxtH[i][j]  = "#color[";
      sTxtH[i][j] += color[j];
      sTxtH[i][j] += "]{I = ";
      sTxtH[i][j] += sIntH[i][j];
      sTxtH[i][j] += "}"; 
    }
  }
  for (Int_t i = 0; i < nProf; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      sTxtP[i][j]  = "#color[";
      sTxtP[i][j] += color[j];
      sTxtP[i][j] += "]{I = ";
      sTxtP[i][j] += sIntP[i][j];
      sTxtP[i][j] += "}"; 
    }
  }

  // create text boxes
  TPaveText *pHist[nHist];
  TPaveText *pProf[nProf];
  for (Int_t i = 0; i < nHist; i++) {
    pHist[i] = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
    pHist[i] -> SetLineColor(0);
    pHist[i] -> SetFillColor(0);
    pHist[i] -> SetFillStyle(0);
    pHist[i] -> SetTextColor(1);
    pHist[i] -> SetTextFont(42);
    pHist[i] -> AddText(sTxtH[i][0]);
    pHist[i] -> AddText(sTxtH[i][1]);
  }
  for (Int_t i = 0; i < nProf; i++) {
    pProf[i] = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
    pProf[i] -> SetLineColor(0);
    pProf[i] -> SetFillColor(0);
    pProf[i] -> SetFillStyle(0);
    pProf[i] -> SetTextColor(1);
    pProf[i] -> SetTextFont(42);
    pProf[i] -> AddText(sTxtP[i][0]);
    pProf[i] -> AddText(sTxtP[i][1]);
  }

  cout << "    Integrals computed..." << endl;


  TLegend *lHist[nHist];
  TLegend *lProf[nProf];
  for (Int_t i = 0; i < nHist; i++) {
    lHist[i] = new TLegend(0.1, 0.1, 0.3, 0.3);
    lHist[i] -> SetLineColor(0);
    lHist[i] -> SetFillColor(0);
    lHist[i] -> SetFillStyle(0);
    lHist[i] -> SetTextColor(1);
    lHist[i] -> SetTextFont(42);
    lHist[i] -> AddEntry(hist[i][0], "StJetTreeThirdMaker");
    lHist[i] -> AddEntry(hist[i][1], "StThirdMaker");
    if (pythia && (i > nEvt))
      lHist[i] -> AddEntry(pyth[i - iPyth], "Pythia, #pi^{0} trigger");
  }
  for (Int_t i = 0; i < nProf; i++) {
    lProf[i] = new TLegend(0.1, 0.1, 0.3, 0.3);
    lProf[i] -> SetLineColor(0);
    lProf[i] -> SetFillColor(0);
    lProf[i] -> SetFillStyle(0);
    lProf[i] -> SetTextColor(1);
    lProf[i] -> SetTextFont(42);
    lProf[i] -> AddEntry(prof[i][0], "StJetTreeThirdMaker");
    lProf[i] -> AddEntry(prof[i][1], "StThirdMaker");
  }

  cout << "    Legends created..." << endl;


  // canvas names
  const TString sCanNamesH[nHist] = {"cPrimeVz", "cEtTrg", "cFtrg", "cHtrg", "cPtTrk", "cDfTrk", "cHtrk"};
  const TString sCanNamesP[nProf] = {"cVzVsRun", "cVzVsEvt", "cVzVsInd", "cRunVsInd", "cEvtVsInd"};
  const TString sPadNamesP[nProf] = {"pVzRunPlot", "pVzEvtPlot", "pVzIndPlot", "pRunIndPlot", "pEvtIndPlot"};
  const TString sPadNamesD[nProf] = {"pVzRunDiff", "pVzEvtDiff", "pVzIndDiff", "pRunIndDiff", "pEvtIndDiff"};

  // create canvases
  fOut -> cd();
  TPad    *pPlot[nProf];
  TPad    *pDiff[nProf];
  TCanvas *cHist[nHist];
  TCanvas *cProf[nProf];
  for (Int_t i = 0; i < nHist; i++) {
    cHist[i] = new TCanvas(sCanNamesH[i].Data(), "", 800, 800);
    cHist[i]   -> SetGrid(0, 0);
    cHist[i]   -> SetLogy(1);
    hist[i][0] -> Draw("");
    hist[i][1] -> Draw("same");
    if (pythia && (i > nEvt))
      pyth[i - iPyth]  -> Draw("same hist");
    lHist[i]   -> Draw();
    pHist[i]   -> Draw();
    cHist[i]   -> Write();
    cHist[i]   -> Close();
  }
  for (Int_t i = 0; i < nProf; i++) {
    cProf[i] = new TCanvas(sCanNamesP[i].Data(), "", 800, 800);
    pDiff[i] = new TPad(sPadNamesD[i].Data(), "", 0, 0, 1, 0.33);
    pPlot[i] = new TPad(sPadNamesP[i].Data(), "", 0, 0.33, 1, 1);
    pDiff[i]   -> SetBorderMode(0);
    pDiff[i]   -> SetBorderSize(2);
    pDiff[i]   -> SetFrameBorderMode(0);
    pDiff[i]   -> SetTopMargin(0);
    pDiff[i]   -> SetGrid(0, 0);
    pPlot[i]   -> SetBorderMode(0);
    pPlot[i]   -> SetBorderSize(2);
    pPlot[i]   -> SetFrameBorderMode(0);
    pPlot[i]   -> SetBottomMargin(0);
    pPlot[i]   -> SetGrid(0, 0);
    pDiff[i]   -> Draw();
    pPlot[i]   -> Draw();
    pDiff[i]   -> cd();
    diff[i]    -> Draw("hist");
    zero[i]    -> Draw("");
    pPlot[i]   -> cd();
    prof[i][0] -> Draw("");
    prof[i][1] -> Draw("same");
    lProf[i]   -> Draw();
    pProf[i]   -> Draw();
    cProf[i]   -> Write();
    cProf[i]   -> Close();
  }

  cout << "    Canvases created..." << endl;


  fOut   -> cd();
  fOut   -> Close();
  fTree  -> cd();
  fTree  -> Close();
  fTuple -> cd();
  fTuple -> Close();

  cout << "  Comparison script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
