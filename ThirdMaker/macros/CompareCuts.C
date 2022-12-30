// 'CompareCuts.C'
// Derek Anderson
// 03.23.2017
//
// Use this to compare the effects of QA cuts
// on the output of the ThirdMaker to the
// effects on the output of the ThirdJetMaker.
//
// NOTE: Please make sure the profiles in
// 'sProf2[nProf2]' below have the *same*
// no. of bins.  Otherwise, you'll need
// to modify how the differences are
// calculated...

#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const TString sTuple("tuple.r10120014-10124011.d18m4y2017.root");
static const TString sTree("tree.after.r10120014-10124011.d17m4y2017.root");
static const TString sOut("treeXtuple.after.r10120014-10124011.d18m4y2017.root");
// global constants
static const Int_t  nDec   = 2;
static const Int_t  nMaker = 2;
static const Int_t  nHist  = 26;
static const Int_t  nProf1 = 5;
static const Int_t  nProf2 = 2;
static const Bool_t doDiff = true;
// hist / prof names
static const TString sHist[nHist]   = {"hPrimeVz", "hEtTrg", "hFtrg", "hHtrg", "hPtTrk0", "hPtTrk1", "hPtTrk2", "hPtTrk3", "hPtTrk4", "hPtTrk5", "hDfTrk0", "hDfTrk1", "hDfTrk2", "hDfTrk3", "hDfTrk4", "hDfTrk5", "hHtrk0", "hHtrk1", "hHtrk2", "hHtrk3", "hHtrk4", "hHtrk5", "hFpTrk", "hPpTrk", "hRpTrk", "hDcaTrk"};
static const TString sProf1[nProf1] = {"hVzVsRun", "hVzVsEvt", "hVzVsInd", "hRunVsInd", "hEvtVsInd"};
static const TString sProf2[nProf2] = {"pNtrkVsCut", "pRatioVsCut"};



void CompareCuts() {

  cout << "\n  Beginning comparison script..." << endl;
  gErrorIgnoreLevel = kError;


  TFile *fOut   = new TFile(sOut.Data(), "recreate");
  TFile *fTree  = new TFile(sTree.Data(), "read");
  TFile *fTuple = new TFile(sTuple.Data(), "read");
  if (!fTree || !fTuple) {
    cerr << "PANIC: couldn't open an input file!" << endl;
    assert(fTree && fTuple);
  }

  // grab histograms / profiles
  TH1D     *hist[nHist][nMaker];
  TH1D     *prof1[nProf1][nMaker];
  TProfile *prof2[nProf2][nMaker];
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
  for (Int_t i = 0; i < nProf1; i++) {
    prof1[i][0] = (TH1D*) fTree  -> Get(sProf1[i].Data());
    prof1[i][1] = (TH1D*) fTuple -> Get(sProf1[i].Data());
    if (!prof1[i][0]) {
      cerr << "PANIC: couldn't grab tree profile(1) no. " << i << "!" << endl;
      assert(prof1[i][0]);
    }
    if (!prof1[i][1]) {
      cerr << "PANIC: couldn't grab tuple profile(1) no. " << i << "!" << endl;
      assert(prof1[i][1]);
    }
  }
  for (Int_t i = 0; i < nProf2; i++) {
    prof2[i][0] = (TProfile*) fTree  -> Get(sProf2[i].Data());
    prof2[i][1] = (TProfile*) fTuple -> Get(sProf2[i].Data());
    if (!prof2[i][0]) {
      cerr << "PANIC: couldn't grab tree profile(2) no. " << i << "!" << endl;
      assert(prof2[i][0]);
    }
    if (!prof2[i][1]) {
      cerr << "PANIC: couldn't grab tuple profile(2) no. " << i << "!" << endl;
      assert(prof2[i][1]);
    }
  }

  cout << "    Histograms and profiles grabbed..." << endl;


  // compute differences and create lines
  TLine *zero0[nHist];
  TH1D  *diff0[nHist];
  for (Int_t i = 0; i < nHist; i++) {
    TString diff0Name("hDiff0_");
    diff0Name += i;

    const Int_t    nX0 = hist[i][0] -> GetNbinsX();
    const Double_t x10 = hist[i][0] -> GetBinLowEdge(1);
    const Double_t x20 = hist[i][0] -> GetBinLowEdge(nX0 + 1);
    zero0[i] = new TLine(x10, 0., x20, 0.);
    diff0[i] = new TH1D(diff0Name.Data(), "", nX0, x10, x20);
    diff0[i] -> Sumw2();

    for (Int_t j = 1; j < nX0 + 1; j++) {
      const Double_t yTree  = hist[i][0] -> GetBinContent(j);
      const Double_t yTuple = hist[i][1] -> GetBinContent(j);
      const Double_t yDiff  = yTree - yTuple;
      diff0[i] -> SetBinContent(j, yDiff);
      diff0[i] -> SetBinError(j, 0);
    } 
  }

  TLine *zero1[nProf1];
  TH1D  *diff1[nProf1];
  for (Int_t i = 0; i < nProf1; i++) {
    TString diff1Name("hDiff1_");
    diff1Name += i;

    const Int_t    nX1 = prof1[i][0] -> GetNbinsX();
    const Double_t x11 = prof1[i][0] -> GetBinLowEdge(1);
    const Double_t x21 = prof1[i][0] -> GetBinLowEdge(nX1 + 1);
    zero1[i] = new TLine(x11, 0., x21, 0.);
    diff1[i] = new TH1D(diff1Name.Data(), "", nX1, x11, x21);
    diff1[i] -> Sumw2();

    for (Int_t j = 1; j < nX1 + 1; j++) {
      const Double_t yTree  = prof1[i][0] -> GetBinContent(j);
      const Double_t yTuple = prof1[i][1] -> GetBinContent(j);
      const Double_t yDiff  = yTree - yTuple;
      diff1[i] -> SetBinContent(j, yDiff);
      diff1[i] -> SetBinError(j, 0);
    } 
  }

  TLine *zero2[nProf2];
  TH1D  *diff2[nProf2];
  for (Int_t i = 0; i < nProf2; i++) {
    TString diff2Name("hDiff2_");
    diff2Name += i;

    const Int_t    nX2 = prof2[i][0] -> GetNbinsX();
    const Double_t x12 = prof2[i][0] -> GetBinLowEdge(1);
    const Double_t x22 = prof2[i][0] -> GetBinLowEdge(nX2 + 1);
    zero2[i] = new TLine(x12, 0., x22, 0.);
    diff2[i] = new TH1D(diff2Name.Data(), "", nX2, x12, x22);
    diff2[i] -> Sumw2();

    for (Int_t j = 1; j < nX2 + 1; j++) {
      const Double_t yTree  = prof2[i][0] -> GetBinContent(j);
      const Double_t yTuple = prof2[i][1] -> GetBinContent(j);
      const Double_t yDiff  = yTree - yTuple;
      diff2[i] -> SetBinContent(j, yDiff);
      diff2[i] -> SetBinError(j, 0);
    }
  }

  cout << "    Differences computed..." << endl;


  // axis labels
  const TString sHistT[nHist]  = {"Primary Vertex z-Coordinate", "Trigger E_{T}", "Trigger #varphi", "Trigger #eta", "Track p_{T}: no cuts", "Track p_{T}: fit cut applied", "Track p_{T}: fit, ratio cuts applied", "Track p_{T}: fit, ratio, dca cuts applied", "Track p_{T}: fit, ratio, dca, p_{T} cuts applied", "Track p_{T}: fit, ratio, dca, p_{T}, #eta cuts applied", "Track #Delta#varphi: no cuts", "Track #Delta#varphi: fit cut applied", "Track #Delta#varphi: fit, ratio cuts applied", "Track #Delta#varphi: fit, ratio, dca cuts applied", "Track #Delta#varphi: fit, ratio, dca, p_{T} cuts applied", "Track #Delta#varphi: fit, ratio, dca, p_{T}, #eta cuts applied", "Track #eta: no cuts", "Track #eta: fit cut applied", "Track #eta: fit, ratio cuts applied", "Track #eta: fit, ratio, dca cuts applied", "Track #eta: fit, ratio, dca, p_{T} cuts applied", "Track #eta: fit, ratio, dca, p_{T}, #eta cuts applied", "Track fit points", "Track possible points", "Track fit points / possible points", "Track DCA"};
  const TString sHistX[nHist]  = {"v_{z}", "E_{T}^{trg}", "#varphi_{trg}", "#eta_{trg}", "p_{T}^{trk}", "p_{T}^{trk}", "p_{T}^{trk}", "p_{T}^{trk}", "p_{T}^{trk}", "p_{T}^{trk}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#Delta#varphi = #varphi_{trk} - #varphi_{trg}", "#eta_{trk}", "#eta_{trk}", "eta_{trk}", "eta_{trk}", "eta_{trk}", "eta_{trk}", "nFit", "nPoss", "nFit / nPoss", "dca"};
  const TString sHistY[nHist]  = {"counts", "counts", "counts", "counts", "dN_{trk}/dp_{T}^{trk}", "dN_{trk}/dp_{T}^{trk}", "dN_{trk}/dp_{T}^{trk}", "dN_{trk}/dp_{T}^{trk}", "dN_{trk}/dp_{T}^{trk}", "dN_{trk}/dp_{T}^{trk}", "dN_{trk}/d#Delta#varphi", "dN_{trk}/d#Delta#varphi", "dN_{trk}/d#Delta#varphi", "dN_{trk}/d#Delta#varphi", "dN_{trk}/d#Delta#varphi", "dN_{trk}/d#Delta#varphi", "dN_{trk}/d#eta_{trk}", "dN_{trk}/d#eta_{trk}", "dN_{trk}/d#eta_{trk}", "dN_{trk}/d#eta_{trk}", "dN_{trk}/d#eta_{trk}", "dN_{trk}/d#eta_{trk}", "counts", "counts", "counts", "counts"};
  const TString sProf1T[nProf1] = {"v_{z} vs. run ID", "v_{z} vs. event ID", "v_{z} vs. tree-index", "run ID vs. tree-index", "event ID vs. tree-index"};
  const TString sProf1X[nProf1] = {"run", "event", "index", "index", "index"};
  const TString sProf1Y[nProf1] = {"v_{z}", "v_{z}", "v_{z}", "run", "event"};
  const TString sProf2T[nProf2] = {"Number of tracks vs. cut number", "R(cut) = N_{trk}(cut) / N_{trk}^{total}"};
  const TString sProf2X[nProf2] = {"cut", "cut"};
  const TString sProf2Y[nProf2] = {"N_{trk}(cut)", "R(cut)"};
  const TString sDiff("Tree - Tuple");

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
  for (Int_t i = 0; i < nProf1; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      prof1[i][j] -> SetTitle(sProf1T[i].Data());
      prof1[i][j] -> SetMarkerColor(color[j]);
      prof1[i][j] -> SetMarkerStyle(marker[j]);
      prof1[i][j] -> SetLineColor(color[j]);
      prof1[i][j] -> SetLineStyle(line[j]);
      prof1[i][j] -> SetLineWidth(width[j]);
      prof1[i][j] -> GetXaxis() -> SetTitle(sProf1X[i].Data());
      prof1[i][j] -> GetXaxis() -> SetTitleFont(font);
      prof1[i][j] -> GetXaxis() -> SetLabelFont(font);
      prof1[i][j] -> GetXaxis() -> CenterTitle(center);
      prof1[i][j] -> GetYaxis() -> SetTitle(sProf1Y[i].Data());
      prof1[i][j] -> GetYaxis() -> SetTitleFont(font);
      prof1[i][j] -> GetYaxis() -> SetLabelFont(font);
      prof1[i][j] -> GetYaxis() -> CenterTitle(center);
    }
  }
  for (Int_t i = 0; i < nProf2; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      prof2[i][j] -> SetTitle(sProf2T[i].Data());
      prof2[i][j] -> SetMarkerColor(color[j]);
      prof2[i][j] -> SetMarkerStyle(marker[j]);
      prof2[i][j] -> SetLineColor(color[j]);
      prof2[i][j] -> SetLineStyle(line[j]);
      prof2[i][j] -> SetLineWidth(width[j]);
      prof2[i][j] -> GetXaxis() -> SetTitle(sProf2X[i].Data());
      prof2[i][j] -> GetXaxis() -> SetTitleFont(font);
      prof2[i][j] -> GetXaxis() -> SetLabelFont(font);
      prof2[i][j] -> GetXaxis() -> CenterTitle(center);
      prof2[i][j] -> GetYaxis() -> SetTitle(sProf2Y[i].Data());
      prof2[i][j] -> GetYaxis() -> SetTitleFont(font);
      prof2[i][j] -> GetYaxis() -> SetLabelFont(font);
      prof2[i][j] -> GetYaxis() -> CenterTitle(center);
    }
  }
  for (Int_t i = 0; i < nHist; i++) {
    diff0[i] -> SetTitle("");
    diff0[i] -> SetMarkerColor(color[2]);
    diff0[i] -> SetMarkerStyle(marker[2]);
    diff0[i] -> SetLineColor(color[2]);
    diff0[i] -> SetLineStyle(line[2]);
    diff0[i] -> SetLineWidth(width[2]);
    diff0[i] -> GetXaxis() -> SetTitle(sHistX[i].Data());
    diff0[i] -> GetXaxis() -> SetTitleFont(font);
    diff0[i] -> GetXaxis() -> SetLabelFont(font);
    diff0[i] -> GetXaxis() -> CenterTitle(center);
    diff0[i] -> GetYaxis() -> SetTitle(sDiff.Data());
    diff0[i] -> GetYaxis() -> SetTitleFont(font);
    diff0[i] -> GetYaxis() -> SetLabelFont(font);
    diff0[i] -> GetYaxis() -> CenterTitle(center);
  }
  for (Int_t i = 0; i < nProf1; i++) {
    diff1[i] -> SetTitle("");
    diff1[i] -> SetMarkerColor(color[2]);
    diff1[i] -> SetMarkerStyle(marker[2]);
    diff1[i] -> SetLineColor(color[2]);
    diff1[i] -> SetLineStyle(line[2]);
    diff1[i] -> SetLineWidth(width[2]);
    diff1[i] -> GetXaxis() -> SetTitle(sProf1X[i].Data());
    diff1[i] -> GetXaxis() -> SetTitleFont(font);
    diff1[i] -> GetXaxis() -> SetLabelFont(font);
    diff1[i] -> GetXaxis() -> CenterTitle(center);
    diff1[i] -> GetYaxis() -> SetTitle(sDiff.Data());
    diff1[i] -> GetYaxis() -> SetTitleFont(font);
    diff1[i] -> GetYaxis() -> SetLabelFont(font);
    diff1[i] -> GetYaxis() -> CenterTitle(center);
  }
  for (Int_t i = 0; i < nProf2; i++) {
    diff2[i] -> SetTitle("");
    diff2[i] -> SetMarkerColor(color[2]);
    diff2[i] -> SetMarkerStyle(marker[2]);
    diff2[i] -> SetLineColor(color[2]);
    diff2[i] -> SetLineStyle(line[2]);
    diff2[i] -> SetLineWidth(width[2]);
    diff2[i] -> GetXaxis() -> SetTitle(sProf2X[i].Data());
    diff2[i] -> GetXaxis() -> SetTitleFont(font);
    diff2[i] -> GetXaxis() -> SetLabelFont(font);
    diff2[i] -> GetXaxis() -> CenterTitle(center);
    diff2[i] -> GetYaxis() -> SetTitle(sDiff.Data());
    diff2[i] -> GetYaxis() -> SetTitleFont(font);
    diff2[i] -> GetYaxis() -> SetLabelFont(font);
    diff2[i] -> GetYaxis() -> CenterTitle(center);
  }

  // set line styles
  const Int_t colorL = 910;
  const Int_t widthL = 1;
  const Int_t styleL = 2;
  for (Int_t i = 0; i < nHist; i++) {
    zero0[i] -> SetLineColor(colorL);
    zero0[i] -> SetLineStyle(styleL);
    zero0[i] -> SetLineWidth(widthL);
  }
  for (Int_t i = 0; i < nProf1; i++) {
    zero1[i] -> SetLineColor(colorL);
    zero1[i] -> SetLineStyle(styleL);
    zero1[i] -> SetLineWidth(widthL);
  }
  for (Int_t i = 0; i < nProf2; i++) {
    zero2[i] -> SetLineColor(colorL);
    zero2[i] -> SetLineStyle(styleL);
    zero2[i] -> SetLineWidth(widthL);
  }

  cout << "    Styles set..." << endl;


  // compute integrals
  Double_t hIntegrals[nHist][nMaker];
  Double_t pIntegrals[nProf1][nMaker];
  for (Int_t i = 0; i < nHist; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      hIntegrals[i][j] = hist[i][j] -> Integral();
    }
  }
  for (Int_t i = 0; i < nProf1; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      pIntegrals[i][j] = prof1[i][j] -> Integral();
    }
  }

  // display only nDec decimals
  Int_t   nDig(0);
  Int_t   size(0);
  TString raw("");
  TString sIntH[nHist][nMaker];
  TString sIntP[nProf1][nMaker];
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
  for (Int_t i = 0; i < nProf1; i++) {
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
  TString sTxtP[nProf1][nMaker];
  for (Int_t i = 0; i < nHist; i++) {
    for (Int_t j = 0; j < nMaker; j++) {
      sTxtH[i][j]  = "#color[";
      sTxtH[i][j] += color[j];
      sTxtH[i][j] += "]{I = ";
      sTxtH[i][j] += sIntH[i][j];
      sTxtH[i][j] += "}"; 
    }
  }
  for (Int_t i = 0; i < nProf1; i++) {
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
  TPaveText *pProf[nProf1];
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
  for (Int_t i = 0; i < nProf1; i++) {
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
  TLegend *lProf1[nProf1];
  TLegend *lProf2[nProf2];
  for (Int_t i = 0; i < nHist; i++) {
    lHist[i] = new TLegend(0.1, 0.1, 0.3, 0.3);
    lHist[i] -> SetLineColor(0);
    lHist[i] -> SetFillColor(0);
    lHist[i] -> SetFillStyle(0);
    lHist[i] -> SetTextColor(1);
    lHist[i] -> SetTextFont(42);
    lHist[i] -> AddEntry(hist[i][0], "StThirdJetMaker");
    lHist[i] -> AddEntry(hist[i][1], "StThirdMaker");
  }
  for (Int_t i = 0; i < nProf1; i++) {
    lProf1[i] = new TLegend(0.1, 0.1, 0.3, 0.3);
    lProf1[i] -> SetLineColor(0);
    lProf1[i] -> SetFillColor(0);
    lProf1[i] -> SetFillStyle(0);
    lProf1[i] -> SetTextColor(1);
    lProf1[i] -> SetTextFont(42);
    lProf1[i] -> AddEntry(prof1[i][0], "StThirdJetMaker");
    lProf1[i] -> AddEntry(prof1[i][1], "StThirdMaker");
  }
  for (Int_t i = 0; i < nProf2; i++) {
    lProf2[i] = new TLegend(0.1, 0.1, 0.3, 0.3);
    lProf2[i] -> SetLineColor(0);
    lProf2[i] -> SetFillColor(0);
    lProf2[i] -> SetFillStyle(0);
    lProf2[i] -> SetTextColor(1);
    lProf2[i] -> SetTextFont(42);
    lProf2[i] -> AddEntry(prof2[i][0], "StThirdJetMaker");
    lProf2[i] -> AddEntry(prof2[i][1], "StThirdMaker");
  }

  cout << "    Legends created..." << endl;


  // canvas names
  const TString sCanNamesH[nHist]   = {"cPrimeVz", "cEtTrg", "cFtrg", "cHtrg", "cPtTrk0", "cPtTrk1", "cPtTrk2", "cPtTrk3", "cPtTrk4", "cPtTrk5", "cDfTrk0", "cDfTrk1", "cDfTrk2", "cDfTrk3", "cDfTrk4", "cDfTrk5", "cHtrk0", "cHtrk1", "cHtrk2", "cHtrk3", "cHtrk4", "cHtrk5", "cFpTrk", "cPpTrk", "cRpTrk", "cDcaTrk"};
  const TString sCanNamesP1[nProf1] = {"cVzVsRun", "cVzVsEvt", "cVzVsInd", "cRunVsInd", "cEvtVsInd"};
  const TString sCanNamesP2[nProf2] = {"cNtrkVsCut", "cRatioVsCut"};
  const TString sPadNamesP0[nHist]  = {"pVzPlot", "pEtPlot", "pFtrgPlot", "pHtrgPlot", "pPtPlot0", "pPtPlot1", "pPtPlot2", "pPtPlot3", "pPtPlot4", "pPtPlot5", "pDfPlot0", "pDfPlot1","pDfPlot2", "pDfPlot3", "pDfPlot4", "pDfPlot5", "pHtrkPlot0", "pHtrkPlot1", "pHtrkPlot2", "pHtrkPlot3", "pHtrkPlot4", "pHtrkPlot5", "pFpPlot", "pPpPlot", "pRpPlot", "pDcaPlot"};
  const TString sPadNamesD0[nHist]  = {"pVzDiff", "pEtDiff", "pFtrgDiff", "pHtrgDiff", "pPtDiff0", "pPtDiff1", "pPtDiff2", "pPtDiff3", "pPtDiff4", "pPtDiff5", "pDfDiff0", "pDfDiff1","pDfDiff2", "pDfDiff3", "pDfDiff4", "pDfDiff5", "pHtrkDiff0", "pHtrkDiff1", "pHtrkDiff2", "pHtrkDiff3", "pHtrkDiff4", "pHtrkDiff5", "pFpDiff", "pPpDiff", "pRpDiff", "pDcaDiff"};
  const TString sPadNamesP1[nProf1] = {"pVzRunPlot", "pVzEvtPlot", "pVzIndPlot", "pRunIndPlot", "pEvtIndPlot"};
  const TString sPadNamesD1[nProf1] = {"pVzRunDiff", "pVzEvtDiff", "pVzIndDiff", "pRunIndDiff", "pEvtIndDiff"};
  const TString sPadNamesP2[nProf2] = {"pNtrkPlot", "pRatioPlot"};
  const TString sPadNamesD2[nProf2] = {"pNtrkDiff", "pRatioDiff"};

  // create canvases
  fOut -> cd();
  TPad    *pPlot0[nHist];
  TPad    *pDiff0[nHist];
  TPad    *pPlot1[nProf1];
  TPad    *pDiff1[nProf1];
  TPad    *pPlot2[nProf2];
  TPad    *pDiff2[nProf2];
  TCanvas *cHist[nHist];
  TCanvas *cProf1[nProf1];
  TCanvas *cProf2[nProf2];
  for (Int_t i = 0; i < nHist; i++) {
    cHist[i] = new TCanvas(sCanNamesH[i].Data(), "", 800, 800);
    if (doDiff) {
      pDiff0[i] = new TPad(sPadNamesD0[i].Data(), "", 0, 0, 1, 0.33);
      pPlot0[i] = new TPad(sPadNamesP0[i].Data(), "", 0, 0.33, 1, 1);
      pDiff0[i]   -> SetBorderMode(0);
      pDiff0[i]   -> SetBorderSize(2);
      pDiff0[i]   -> SetFrameBorderMode(0);
      pDiff0[i]   -> SetTopMargin(0);
      pDiff0[i]   -> SetGrid(0, 0);
      pPlot0[i]   -> SetBorderMode(0);
      pPlot0[i]   -> SetBorderSize(2);
      pPlot0[i]   -> SetFrameBorderMode(0);
      pPlot0[i]   -> SetBottomMargin(0);
      pPlot0[i]   -> SetGrid(0, 0);
      pDiff0[i]   -> Draw();
      pPlot0[i]   -> Draw();
      pDiff0[i]   -> cd();
      diff0[i]    -> Draw("");
      zero0[i]    -> Draw("");
      pPlot0[i]   -> cd();
      hist[i][0]  -> Draw("");
      hist[i][1]  -> Draw("same");
      lHist[i]    -> Draw();
      pHist[i]    -> Draw();
      cHist[i]    -> Write();
      cHist[i]    -> Close();
    }
    else {
      cHist[i]   -> SetGrid(0, 0);
      cHist[i]   -> SetLogy(1);
      hist[i][0] -> Draw("");
      hist[i][1] -> Draw("same");
      lHist[i]   -> Draw();
      pHist[i]   -> Draw();
      cHist[i]   -> Write();
      cHist[i]   -> Close();
    }
  }
  for (Int_t i = 0; i < nProf1; i++) {
    cProf1[i] = new TCanvas(sCanNamesP1[i].Data(), "", 800, 800);
    pDiff1[i] = new TPad(sPadNamesD1[i].Data(), "", 0, 0, 1, 0.33);
    pPlot1[i] = new TPad(sPadNamesP1[i].Data(), "", 0, 0.33, 1, 1);
    pDiff1[i]   -> SetBorderMode(0);
    pDiff1[i]   -> SetBorderSize(2);
    pDiff1[i]   -> SetFrameBorderMode(0);
    pDiff1[i]   -> SetTopMargin(0);
    pDiff1[i]   -> SetGrid(0, 0);
    pPlot1[i]   -> SetBorderMode(0);
    pPlot1[i]   -> SetBorderSize(2);
    pPlot1[i]   -> SetFrameBorderMode(0);
    pPlot1[i]   -> SetBottomMargin(0);
    pPlot1[i]   -> SetGrid(0, 0);
    pDiff1[i]   -> Draw();
    pPlot1[i]   -> Draw();
    pDiff1[i]   -> cd();
    diff1[i]    -> Draw("hist");
    zero1[i]    -> Draw("");
    pPlot1[i]   -> cd();
    prof1[i][0] -> Draw("");
    prof1[i][1] -> Draw("same");
    lProf1[i]   -> Draw();
    pProf[i]    -> Draw();
    cProf1[i]   -> Write();
    cProf1[i]   -> Close();
  }
  for (Int_t i = 0; i < nProf2; i++) {
    cProf2[i] = new TCanvas(sCanNamesP2[i].Data(), "", 800, 800);
    if (doDiff) {
      pDiff2[i] = new TPad(sPadNamesD2[i].Data(), "", 0, 0, 1, 0.33);
      pPlot2[i] = new TPad(sPadNamesP2[i].Data(), "", 0, 0.33, 1, 1);
      pDiff2[i]   -> SetBorderMode(0);
      pDiff2[i]   -> SetBorderSize(2);
      pDiff2[i]   -> SetFrameBorderMode(0);
      pDiff2[i]   -> SetTopMargin(0);
      pDiff2[i]   -> SetGrid(0, 0);
      pPlot2[i]   -> SetBorderMode(0);
      pPlot2[i]   -> SetBorderSize(2);
      pPlot2[i]   -> SetFrameBorderMode(0);
      pPlot2[i]   -> SetBottomMargin(0);
      pPlot2[i]   -> SetGrid(0, 0);
      pDiff2[i]   -> Draw();
      pPlot2[i]   -> Draw();
      pDiff2[i]   -> cd();
      diff2[i]    -> Draw("");
      zero2[i]    -> Draw("");
      pPlot2[i]   -> cd();
      prof2[i][0] -> Draw("");
      prof2[i][1] -> Draw("same");
      lProf2[i]   -> Draw();
      cProf2[i]   -> Write();
      cProf2[i]   -> Close();
    }
    else {
      cProf2[i]   -> SetGrid(0, 0);
      cProf2[i]   -> SetLogy(1);
      prof2[i][0] -> Draw("");
      prof2[i][1] -> Draw("same");
      lProf2[i]   -> Draw();
      cProf2[i]   -> Write();
      cProf2[i]   -> Close();
    }
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
