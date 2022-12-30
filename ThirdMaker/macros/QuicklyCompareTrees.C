// 'QuicklyCompareTrees.C'
// Derek Anderseon
// 02.22.2017
//
// Use this to (quickly) compare the output of 'StThirdMaker' to
// the output of 'StThirdJetMaker'.  Run 'MergeTrees.C' first.
//
// NOTE: that this assumes both makers were run on the *exact*
// same set of events.
//
// NOTE: also that if any of the leaves drawn are empty (for
// one reason or another) this will break.

#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const Int_t   nDec(2);
static const Int_t   nTuples(2);
static const Int_t   nMakers(2);
static const Int_t   nVariables(6);
static const TString sIn("merge.d13m3y2017.root");
static const TString sOut("comparison.d13m3y2017.root");



void QuicklyCompareTrees() {

  cout << "\n  Beginning comparison script..." << endl;
  gErrorIgnoreLevel = kError;


  // open files and load chains
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn  = new TFile(sIn.Data(), "read");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    assert(fIn);
  }

  cout << "    Loading chains..." << endl;
  TChain *tTree;
  TChain *tTuple[nTuples];
  if (fIn) {
    fIn -> GetObject("Gfmtodst", tTree);
    fIn -> GetObject("QAAB", tTuple[0]);
    fIn -> GetObject("Track", tTuple[1]);
  }

  // check if chains were loaded
  Bool_t tupleOpen  = true;
  Bool_t tuplesGood = true;
  for (Int_t i = 0; i < nTuples; i++) {
    tupleOpen   = (Bool_t) tTuple[i];
    tuplesGood *= tupleOpen;
  }
  if (!tTree || !tuplesGood) {
    cerr << "PANIC: couldn't grab a tree or tuple!" << endl;
    assert(tTree && tuplesGood);
  }


  cout << "    Drawing tree histograms..." << endl;
  TCanvas *cTree = new TCanvas("cTree", "", 200, 200);

  TH1D *hEtrgTree;
  TH1D *hFtrgTree;
  TH1D *hHtrgTree;
  TH1D *hPtTrkTree;
  TH1D *hFtrkTree;
  TH1D *hHtrkTree;
  tTree -> Draw("EClustEneT0>>hEtrg(500,0.,50.)");
  hEtrgTree  = (TH1D*) gDirectory -> Get("hEtrg");
  tTree -> Draw("EClustphiv1>>hFtrg(80,-4.,4.)");
  hFtrgTree  = (TH1D*) gDirectory -> Get("hFtrg");
  tTree -> Draw("EClustetav1>>hHtrg(100,-5.,5.)");
  hHtrgTree  = (TH1D*) gDirectory -> Get("hHtrg");
  tTree -> Draw("PrimaryTrackArray.pT>>hPtTrk(500,0.,50.)");
  hPtTrkTree = (TH1D*) gDirectory -> Get("hPtTrk");
  tTree -> Draw("PrimaryTrackArray.phi>>hFtrk(80,-4.,4.)");
  hFtrkTree  = (TH1D*) gDirectory -> Get("hFtrk");
  tTree -> Draw("PrimaryTrackArray.eta>>hHtrk(100,-5.,5.)");
  hHtrkTree  = (TH1D*) gDirectory -> Get("hHtrk");

  cout << "    Tree histograms drawn." << endl;
  cTree -> Close();

  cout << "    Setting tree names and styles..." << endl;
  const Int_t oTree  = 860;
  const Int_t mTree  = 7;
  const Int_t lTree  = 1;
  const Int_t wTree  = 2;
  const Int_t fTree  = 3004;
  const Int_t font   = 42;
  const Int_t center = 1;
  hEtrgTree  -> SetNameTitle("hEtrgTree", "E^{trg}, StThirdJetMaker");
  hEtrgTree  -> SetMarkerColor(oTree);
  hEtrgTree  -> SetMarkerStyle(mTree);
  hEtrgTree  -> SetLineColor(oTree);
  hEtrgTree  -> SetLineStyle(lTree);
  hEtrgTree  -> SetLineWidth(wTree);
  hEtrgTree  -> SetFillColor(oTree);
  hEtrgTree  -> SetFillStyle(fTree);
  hEtrgTree  -> GetXaxis() -> SetTitle("E^{trg}");
  hEtrgTree  -> GetXaxis() -> SetTitleFont(font);
  hEtrgTree  -> GetXaxis() -> CenterTitle(center);
  hEtrgTree  -> GetYaxis() -> SetTitle("counts");
  hEtrgTree  -> GetYaxis() -> SetTitleFont(font);
  hEtrgTree  -> GetYaxis() -> CenterTitle(center);
  hFtrgTree  -> SetNameTitle("hFtrgTree", "#varphi^{trg}, StThirdJetMaker");
  hFtrgTree  -> SetMarkerColor(oTree);
  hFtrgTree  -> SetMarkerStyle(mTree);
  hFtrgTree  -> SetLineColor(oTree);
  hFtrgTree  -> SetLineStyle(lTree);
  hFtrgTree  -> SetLineWidth(wTree);
  hFtrgTree  -> SetFillColor(oTree);
  hFtrgTree  -> SetFillStyle(fTree);
  hFtrgTree  -> GetXaxis() -> SetTitle("#varphi^{trg}");
  hFtrgTree  -> GetXaxis() -> SetTitleFont(font);
  hFtrgTree  -> GetXaxis() -> CenterTitle(center);
  hFtrgTree  -> GetYaxis() -> SetTitle("counts");
  hFtrgTree  -> GetYaxis() -> SetTitleFont(font);
  hFtrgTree  -> GetYaxis() -> CenterTitle(center);
  hHtrgTree  -> SetNameTitle("hHtrgTree", "#eta^{trg}, StThirdJetMaker");
  hHtrgTree  -> SetMarkerColor(oTree);
  hHtrgTree  -> SetMarkerStyle(mTree);
  hHtrgTree  -> SetLineColor(oTree);
  hHtrgTree  -> SetLineStyle(lTree);
  hHtrgTree  -> SetLineWidth(wTree);
  hHtrgTree  -> SetFillColor(oTree);
  hHtrgTree  -> SetFillStyle(fTree);
  hHtrgTree  -> GetXaxis() -> SetTitle("#eta^{trg}");
  hHtrgTree  -> GetXaxis() -> SetTitleFont(font);
  hHtrgTree  -> GetXaxis() -> CenterTitle(center);
  hHtrgTree  -> GetYaxis() -> SetTitle("counts");
  hHtrgTree  -> GetYaxis() -> SetTitleFont(font);
  hHtrgTree  -> GetYaxis() -> CenterTitle(center);
  hPtTrkTree -> SetNameTitle("hPtTrkTree", "p_{T}^{trk}, StThirdJetMaker");
  hPtTrkTree -> SetMarkerColor(oTree);
  hPtTrkTree -> SetMarkerStyle(mTree);
  hPtTrkTree -> SetLineColor(oTree);
  hPtTrkTree -> SetLineStyle(lTree);
  hPtTrkTree -> SetLineWidth(wTree);
  hPtTrkTree -> SetFillColor(oTree);
  hPtTrkTree -> SetFillStyle(fTree);
  hPtTrkTree -> GetXaxis() -> SetTitle("p_{T}^{trk}");
  hPtTrkTree -> GetXaxis() -> SetTitleFont(font);
  hPtTrkTree -> GetXaxis() -> CenterTitle(center);
  hPtTrkTree -> GetYaxis() -> SetTitle("counts");
  hPtTrkTree -> GetYaxis() -> SetTitleFont(font);
  hPtTrkTree -> GetYaxis() -> CenterTitle(center);
  hFtrkTree  -> SetNameTitle("hFtrkTree", "#varphi^{trk}, StThirdJetMaker");
  hFtrkTree  -> SetMarkerColor(oTree);
  hFtrkTree  -> SetMarkerStyle(mTree);
  hFtrkTree  -> SetLineColor(oTree);
  hFtrkTree  -> SetLineStyle(lTree);
  hFtrkTree  -> SetLineWidth(wTree);
  hFtrkTree  -> SetFillColor(oTree);
  hFtrkTree  -> SetFillStyle(fTree);
  hFtrkTree  -> GetXaxis() -> SetTitle("#varphi^{trk}");
  hFtrkTree  -> GetXaxis() -> SetTitleFont(font);
  hFtrkTree  -> GetXaxis() -> CenterTitle(center);
  hFtrkTree  -> GetYaxis() -> SetTitle("counts");
  hFtrkTree  -> GetYaxis() -> SetTitleFont(font);
  hFtrkTree  -> GetYaxis() -> CenterTitle(center);
  hHtrkTree  -> SetNameTitle("hHtrkTree", "#eta^{trk}, StThirdJetMaker");
  hHtrkTree  -> SetMarkerColor(oTree);
  hHtrkTree  -> SetMarkerStyle(mTree);
  hHtrkTree  -> SetLineColor(oTree);
  hHtrkTree  -> SetLineStyle(lTree);
  hHtrkTree  -> SetLineWidth(wTree);
  hHtrkTree  -> SetFillColor(oTree);
  hHtrkTree  -> SetFillStyle(fTree);
  hHtrkTree  -> GetXaxis() -> SetTitle("#eta^{trk}");
  hHtrkTree  -> GetXaxis() -> SetTitleFont(font);
  hHtrkTree  -> GetXaxis() -> CenterTitle(center);
  hHtrkTree  -> GetYaxis() -> SetTitle("counts");
  hHtrkTree  -> GetYaxis() -> SetTitleFont(font);
  hHtrkTree  -> GetYaxis() -> CenterTitle(center);


  cout  << "    Drawing tuple histograms..." << endl;
  TCanvas *cTuple = new TCanvas("cTuple", "", 200, 200);
        
  TH1D *hEtrgTuple;
  TH1D *hFtrgTuple;
  TH1D *hHtrgTuple;
  TH1D *hPtTrkTuple;
  TH1D *hFtrkTuple;
  TH1D *hHtrkTuple;
  tTuple[0] -> Draw("EneT0>>hEtrg(500,0.,50.)");
  hEtrgTuple  = (TH1D*) gDirectory -> Get("hEtrg");
  tTuple[0] -> Draw("phiv1>>hFtrg(80,-4.,4.)");
  hFtrgTuple  = (TH1D*) gDirectory -> Get("hFtrg");
  tTuple[0] -> Draw("etav1>>hHtrg(100,-5.,5.)");
  hHtrgTuple  = (TH1D*) gDirectory -> Get("hHtrg");
  tTuple[1] -> Draw("ptTPC>>hPtTrk(500,0.,50.)");
  hPtTrkTuple = (TH1D*) gDirectory -> Get("hPtTrk");
  tTuple[1] -> Draw("epsi>>hFtrk(80,-4.,4.)");
  hFtrkTuple  = (TH1D*) gDirectory -> Get("hFtrk");
  tTuple[1] -> Draw("EtaTr>>hHtrk(100,-5.,5.)");
  hHtrkTuple  = (TH1D*) gDirectory -> Get("hHtrk");

  cout << "    Tuple histograms drawn." << endl;
  cTuple -> Close();

  cout << "    Setting tuple names and style..." << endl;
  const Int_t oTuple = 810;
  const Int_t mTuple = 4;
  const Int_t lTuple = 2;
  const Int_t wTuple = 1;
  const Int_t fTuple = 3005;
  hEtrgTuple  -> SetNameTitle("hEtrgTuple", "E^{trg}, StThirdMaker");
  hEtrgTuple  -> SetMarkerColor(oTuple);
  hEtrgTuple  -> SetMarkerStyle(mTuple);
  hEtrgTuple  -> SetLineColor(oTuple);
  hEtrgTuple  -> SetLineStyle(lTuple);
  hEtrgTuple  -> SetLineWidth(wTuple);
  hEtrgTuple  -> SetFillColor(oTuple);
  hEtrgTuple  -> SetFillStyle(fTuple);
  hEtrgTuple  -> GetXaxis() -> SetTitle("E^{trg}");
  hEtrgTuple  -> GetXaxis() -> SetTitleFont(font);
  hEtrgTuple  -> GetXaxis() -> CenterTitle(center);
  hEtrgTuple  -> GetYaxis() -> SetTitle("counts");
  hEtrgTuple  -> GetYaxis() -> SetTitleFont(font);
  hEtrgTuple  -> GetYaxis() -> CenterTitle(center);
  hFtrgTuple  -> SetNameTitle("hFtrgTuple", "#varphi^{trg}, StThirdMaker");
  hFtrgTuple  -> SetMarkerColor(oTuple);
  hFtrgTuple  -> SetMarkerStyle(mTuple);
  hFtrgTuple  -> SetLineColor(oTuple);
  hFtrgTuple  -> SetLineStyle(lTuple);
  hFtrgTuple  -> SetLineWidth(wTuple);
  hFtrgTuple  -> SetFillColor(oTuple);
  hFtrgTuple  -> SetFillStyle(fTuple);
  hFtrgTuple  -> GetXaxis() -> SetTitle("#varphi^{trg}");
  hFtrgTuple  -> GetXaxis() -> SetTitleFont(font);
  hFtrgTuple  -> GetXaxis() -> CenterTitle(center);
  hFtrgTuple  -> GetYaxis() -> SetTitle("counts");
  hFtrgTuple  -> GetYaxis() -> SetTitleFont(font);
  hFtrgTuple  -> GetYaxis() -> CenterTitle(center);
  hHtrgTuple  -> SetNameTitle("hHtrgTuple", "#eta^{trg}, StThirdMaker");
  hHtrgTuple  -> SetMarkerColor(oTuple);
  hHtrgTuple  -> SetMarkerStyle(mTuple);
  hHtrgTuple  -> SetLineColor(oTuple);
  hHtrgTuple  -> SetLineStyle(lTuple);
  hHtrgTuple  -> SetLineWidth(wTuple);
  hHtrgTuple  -> SetFillColor(oTuple);
  hHtrgTuple  -> SetFillStyle(fTuple);
  hHtrgTuple  -> GetXaxis() -> SetTitle("#eta^{trg}");
  hHtrgTuple  -> GetXaxis() -> SetTitleFont(font);
  hHtrgTuple  -> GetXaxis() -> CenterTitle(center);
  hHtrgTuple  -> GetYaxis() -> SetTitle("counts");
  hHtrgTuple  -> GetYaxis() -> SetTitleFont(font);
  hHtrgTuple  -> GetYaxis() -> CenterTitle(center);
  hPtTrkTuple -> SetNameTitle("hPtTrkTuple", "p_{T}^{trk}, StThirdMaker");
  hPtTrkTuple -> SetMarkerColor(oTuple);
  hPtTrkTuple -> SetMarkerStyle(mTuple);
  hPtTrkTuple -> SetLineColor(oTuple);
  hPtTrkTuple -> SetLineStyle(lTuple);
  hPtTrkTuple -> SetLineWidth(wTuple);
  hPtTrkTuple -> SetFillColor(oTuple);
  hPtTrkTuple -> SetFillStyle(fTuple);
  hPtTrkTuple -> GetXaxis() -> SetTitle("p_{T}^{trk}");
  hPtTrkTuple -> GetXaxis() -> SetTitleFont(font);
  hPtTrkTuple -> GetXaxis() -> CenterTitle(center);
  hPtTrkTuple -> GetYaxis() -> SetTitle("counts");
  hPtTrkTuple -> GetYaxis() -> SetTitleFont(font);
  hPtTrkTuple -> GetYaxis() -> CenterTitle(center);
  hFtrkTuple  -> SetNameTitle("hFtrkTuple", "#varphi^{trk}, StThirdMaker");
  hFtrkTuple  -> SetMarkerColor(oTuple);
  hFtrkTuple  -> SetMarkerStyle(mTuple);
  hFtrkTuple  -> SetLineColor(oTuple);
  hFtrkTuple  -> SetLineStyle(lTuple);
  hFtrkTuple  -> SetLineWidth(wTuple);
  hFtrkTuple  -> SetFillColor(oTuple);
  hFtrkTuple  -> SetFillStyle(fTuple);
  hFtrkTuple  -> GetXaxis() -> SetTitle("#varphi^{trk}");
  hFtrkTuple  -> GetXaxis() -> SetTitleFont(font);
  hFtrkTuple  -> GetXaxis() -> CenterTitle(center);
  hFtrkTuple  -> GetYaxis() -> SetTitle("counts");
  hFtrkTuple  -> GetYaxis() -> SetTitleFont(font);
  hFtrkTuple  -> GetYaxis() -> CenterTitle(center);
  hHtrkTuple  -> SetNameTitle("hHtrkTuple", "#eta^{trk}, StThirdMaker");
  hHtrkTuple  -> SetMarkerColor(oTuple);
  hHtrkTuple  -> SetMarkerStyle(mTuple);
  hHtrkTuple  -> SetLineColor(oTuple);
  hHtrkTuple  -> SetLineStyle(lTuple);
  hHtrkTuple  -> SetLineWidth(wTuple);
  hHtrkTuple  -> SetFillColor(oTuple);
  hHtrkTuple  -> SetFillStyle(fTuple);
  hHtrkTuple  -> GetXaxis() -> SetTitle("#eta^{trk}");
  hHtrkTuple  -> GetXaxis() -> SetTitleFont(font);
  hHtrkTuple  -> GetXaxis() -> CenterTitle(center);
  hHtrkTuple  -> GetYaxis() -> SetTitle("counts");
  hHtrkTuple  -> GetYaxis() -> SetTitleFont(font);
  hHtrkTuple  -> GetYaxis() -> CenterTitle(center);


  cout << "    Computing integrals..." << endl;
  Double_t integrals[nVariables][nMakers];
  integrals[0][0] = hEtrgTree   -> Integral();
  integrals[0][1] = hEtrgTuple  -> Integral();
  integrals[1][0] = hFtrgTree   -> Integral();
  integrals[1][1] = hFtrgTuple  -> Integral();
  integrals[2][0] = hHtrgTree   -> Integral();
  integrals[2][1] = hHtrgTuple  -> Integral();
  integrals[3][0] = hPtTrkTree  -> Integral();
  integrals[3][1] = hPtTrkTuple -> Integral();
  integrals[4][0] = hFtrkTree   -> Integral();
  integrals[4][1] = hFtrkTuple  -> Integral();
  integrals[5][0] = hHtrkTree   -> Integral();
  integrals[5][1] = hHtrkTuple  -> Integral();

  // display only nDec decimals
  Int_t   nDig(0);
  Int_t   size(0);
  TString raw("");
  TString sInt[nVariables][nMakers];
  for (Int_t i = 0; i < nVariables; i++) {
    for (Int_t j = 0; j < nMakers; j++) {
      raw        = "";
      raw       += integrals[i][j];
      nDig       = raw.First(".");
      if (nDig > 0)
        size     = (nDig + 1) + nDec;
      else
        size     = raw.Length();
      sInt[i][j] = "";
      sInt[i][j].Append(raw, size);
    }
  }

  // set up text
  TString sEnd("}");
  TString sStart[nVariables];
  TString sColTree("#color[");
  TString sColTuple("#color[");
  sStart[0]  = "I(E^{trg}) = ";
  sStart[1]  = "I(#varphi^{trg}) = ";
  sStart[2]  = "I(#eta^{trg}) = ";
  sStart[3]  = "I(p_{T}^{trk}) = ";
  sStart[4]  = "I(#varphi^{trk}) = ";
  sStart[5]  = "I(#eta^{trk}) = ";
  sColTree  += oTree;
  sColTree  += "]{";
  sColTuple += oTuple;
  sColTuple += "]{";

  // create text boxes
  const Int_t oLine = 0;
  const Int_t oFill = 0;
  const Int_t fTxt  = 0;
  const Int_t oTxt  = 1;
  TString   sEtrgTree(sColTree);
  TString   sEtrgTuple(sColTuple);
  TPaveText *pEtrg = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pEtrg -> SetLineColor(oLine);
  pEtrg -> SetFillColor(oFill);
  pEtrg -> SetFillStyle(fTxt);
  pEtrg -> SetTextColor(oTxt);
  pEtrg -> SetTextFont(font);
  sEtrgTree  += sStart[0];
  sEtrgTree  += sInt[0][0];
  sEtrgTree  += sEnd;
  sEtrgTuple += sStart[0];
  sEtrgTuple += sInt[0][1];
  sEtrgTuple += sEnd;
  pEtrg -> AddText(sEtrgTree.Data());
  pEtrg -> AddText(sEtrgTuple.Data());

  TString   sFtrgTree(sColTree);
  TString   sFtrgTuple(sColTuple);
  TPaveText *pFtrg = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pFtrg -> SetLineColor(oLine);
  pFtrg -> SetFillColor(oFill);
  pFtrg -> SetFillStyle(fTxt);
  pFtrg -> SetTextColor(oTxt);
  pFtrg -> SetTextFont(font);
  sFtrgTree  += sStart[1];
  sFtrgTree  += sInt[1][0];
  sFtrgTree  += sEnd;
  sFtrgTuple += sStart[1];
  sFtrgTuple += sInt[1][1];
  sFtrgTuple += sEnd;
  pFtrg -> AddText(sFtrgTree.Data());
  pFtrg -> AddText(sFtrgTuple.Data());

  TString   sHtrgTree(sColTree);
  TString   sHtrgTuple(sColTuple);
  TPaveText *pHtrg = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pHtrg -> SetLineColor(oLine);
  pHtrg -> SetFillColor(oFill);
  pHtrg -> SetFillStyle(fTxt);
  pHtrg -> SetTextColor(oTxt);
  pHtrg -> SetTextFont(font);
  sHtrgTree  += sStart[2];
  sHtrgTree  += sInt[2][0];
  sHtrgTree  += sEnd;
  sHtrgTuple += sStart[2];
  sHtrgTuple += sInt[2][1];
  sHtrgTuple += sEnd;
  pHtrg -> AddText(sHtrgTree.Data());
  pHtrg -> AddText(sHtrgTuple.Data());

  TString   sPtTrkTree(sColTree);
  TString   sPtTrkTuple(sColTuple);
  TPaveText *pPtTrk = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pPtTrk -> SetLineColor(oLine);
  pPtTrk -> SetFillColor(oFill);
  pPtTrk -> SetFillStyle(fTxt);
  pPtTrk -> SetTextColor(oTxt);
  pPtTrk -> SetTextFont(font);
  sPtTrkTree  += sStart[3];
  sPtTrkTree  += sInt[3][0];
  sPtTrkTree  += sEnd;
  sPtTrkTuple += sStart[3];
  sPtTrkTuple += sInt[3][1];
  sPtTrkTuple += sEnd;
  pPtTrk -> AddText(sPtTrkTree.Data());
  pPtTrk -> AddText(sPtTrkTuple.Data());

  TString   sFtrkTree(sColTree);
  TString   sFtrkTuple(sColTuple);
  TPaveText *pFtrk = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pFtrk -> SetLineColor(oLine);
  pFtrk -> SetFillColor(oFill);
  pFtrk -> SetFillStyle(fTxt);
  pFtrk -> SetTextColor(oTxt);
  pFtrk -> SetTextFont(font);
  sFtrkTree  += sStart[4];
  sFtrkTree  += sInt[4][0];
  sFtrkTree  += sEnd;
  sFtrkTuple += sStart[4];
  sFtrkTuple += sInt[4][1];
  sFtrkTuple += sEnd;
  pFtrk -> AddText(sFtrkTree.Data());
  pFtrk -> AddText(sFtrkTuple.Data());

  TString   sHtrkTree(sColTree);
  TString   sHtrkTuple(sColTuple);
  TPaveText *pHtrk = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pHtrk -> SetLineColor(oLine);
  pHtrk -> SetFillColor(oFill);
  pHtrk -> SetFillStyle(fTxt);
  pHtrk -> SetTextColor(oTxt);
  pHtrk -> SetTextFont(font);
  sHtrkTree  += sStart[4];
  sHtrkTree  += sInt[4][0];
  sHtrkTree  += sEnd;
  sHtrkTuple += sStart[4];
  sHtrkTuple += sInt[4][1];
  sHtrkTuple += sEnd;
  pHtrk -> AddText(sHtrkTree.Data());
  pHtrk -> AddText(sHtrkTuple.Data());


  cout << "    Creating legend..." << endl;
  TLegend *lg = new TLegend(0.1, 0.1, 0.3, 0.3);
  lg -> SetLineColor(oLine);
  lg -> SetFillColor(oFill);
  lg -> SetFillStyle(fTxt);
  lg -> SetTextColor(oTxt);
  lg -> SetTextFont(font);
  lg -> AddEntry(hEtrgTree, "StThirdJetMaker");
  lg -> AddEntry(hEtrgTuple, "StThirdMaker");


  cout << "    Drawing plots..." << endl;
  fOut -> cd();

  TCanvas *cEtrg = new TCanvas("cEtrg", "E^{trg}, StThirdMaker vs. StThirdJetMaker", 800, 800);
  cEtrg      -> SetGrid(0, 0);
  cEtrg      -> SetLogy(1);
  hEtrgTree  -> SetTitle("Trigger energy, E^{trg}");
  hEtrgTree  -> Draw("");
  hEtrgTuple -> Draw("same");
  lg         -> Draw();
  pEtrg      -> Draw();
  cEtrg      -> Write();
  cEtrg      -> Close();

  TCanvas *cFtrg = new TCanvas("cFtrg", "#varphi^{trg}, StThirdMaker vs. StThirdJetMaker", 800, 800);
  cFtrg      -> SetGrid(0, 0);
  cFtrg      -> SetLogy(1);
  hFtrgTree  -> SetTitle("Trigger azimuth, #varphi^{trg}");
  hFtrgTree  -> Draw("");
  hFtrgTuple -> Draw("same");
  lg         -> Draw();
  pFtrg      -> Draw();
  cFtrg      -> Write();
  cFtrg      -> Close();

  TCanvas *cHtrg = new TCanvas("cHtrg", "#eta^{trg}, StThirdMaker vs. StThirdJetMaker", 800, 800);
  cHtrg      -> SetGrid(0, 0);
  cHtrg      -> SetLogy(1);
  hHtrgTree  -> SetTitle("Trigger pseudorapdity, #eta^{trg}");
  hHtrgTree  -> Draw("");
  hHtrgTuple -> Draw("same");
  lg         -> Draw();
  pHtrg      -> Draw();
  cHtrg      -> Write();
  cHtrg      -> Close();

  TCanvas *cPtTrk = new TCanvas("cPtTrk", "p_{T}^{trk}, StThirdMaker vs. StThirdJetMaker", 800, 800);
  cPtTrk      -> SetGrid(0, 0);
  cPtTrk      -> SetLogy(1);
  hPtTrkTree  -> SetTitle("Primary track p_{T}, p_{T}^{trk}");
  hPtTrkTree  -> Draw("");
  hPtTrkTuple -> Draw("same");
  lg          -> Draw();
  pPtTrk      -> Draw();
  cPtTrk      -> Write();
  cPtTrk      -> Close();

  TCanvas *cFtrk = new TCanvas("cFtrk", "#varphi^{trk}, StThirdMaker vs. StThirdJetMaker", 800, 800);
  cFtrk      -> SetGrid(0, 0);
  cFtrk      -> SetLogy(1);
  hFtrkTree  -> SetTitle("Primary track azimuth, #varphi^{trk}");
  hFtrkTree  -> Draw("");
  hFtrkTuple -> Draw("same");
  lg         -> Draw();
  pFtrk      -> Draw();
  cFtrk      -> Write();
  cFtrk      -> Close();

  TCanvas *cHtrk = new TCanvas("cHtrk", "#eta^{trk}, StThirdMaker vs. StThirdJetMaker", 800, 800);
  cHtrk      -> SetGrid(0, 0);
  cHtrk      -> SetLogy(1);
  hHtrkTree  -> SetTitle("Primary track pseudorapidity, #eta^{trk}");
  hHtrkTree  -> Draw("");
  hHtrkTuple -> Draw("same");
  lg         -> Draw();
  pHtrk      -> Draw();
  cHtrk      -> Write();
  cHtrk      -> Close();


  cout  << "    Saving and closing files..." << endl;

  fOut        -> cd();
  hEtrgTree   -> Write();
  hFtrgTree   -> Write();
  hHtrgTree   -> Write();
  hPtTrkTree  -> Write();
  hFtrkTree   -> Write();
  hHtrkTree   -> Write();
  hEtrgTuple  -> Write();
  hFtrgTuple  -> Write();
  hHtrgTuple  -> Write();
  hPtTrkTuple -> Write();
  hFtrkTuple  -> Write();
  hHtrkTuple  -> Write();
  fOut        -> Close();
  fIn         -> cd();
  fIn         -> Close();

  cout << "  Comparison script finished.\n" << endl;

}

// End ------------------------------------------------------------------------
