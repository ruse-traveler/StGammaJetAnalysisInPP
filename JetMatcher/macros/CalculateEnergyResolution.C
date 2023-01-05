// 'CalculateEnergyResolution.C'
// Derek Anderson
// 07.11.2018
//
// Use this to calculate the JER as a function
// of particle-lvl. jet pT from a response
// matrix.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;


// global constants
static const UInt_t  NDec(2);
static const UInt_t  NPtBins(3);
static const UInt_t  PtCols[NPtBins] = {899, 859, 819};
static const Float_t PtBins[NPtBins] = {5., 10., 20.};
static const Float_t FitLimits(3.);



void CalculateEnergyResolution() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning JER calculation..." << endl;

  // input parameters
  const TString sIn("output/October2021/pp200r9embed.forPlots_pTbinOne.et920pt0230x021Kvz55pi0.r02a005rm1chrg.dr02qt05130.d7m10y2021.root");
  const TString sOut("jetEnergyResolution.forLongPaper_pTbinOne.et911r02pi0.dr02qt05130.d8m1y2022.root");
  const TString sMatrix("hResponseAll");

  // style parameters
  const UInt_t  fMar(4);
  const UInt_t  fLin(1);
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const UInt_t  fColL(0);
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const UInt_t  fSizX(750);
  const UInt_t  fSizY(750);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  const Float_t fMargin(0.02);
  const TString sTitle("");
  const TString sTitleX("p_{T}^{det} [GeV/c]");
  const TString sTitleY("a. u.");

  // misc parameters
  const UInt_t  fAln(12);
  const Float_t fLeg1(0.1);
  const Float_t fLeg2(0.3);
  const Float_t fTxt1(0.3);
  const Float_t fTxt2(0.5);
  const TString sSystem("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrigger("#pi^{0} trigger, E_{T}^{trg} #in (9, 11) GeV");
  const TString sJets("anti-k_{T}, R = 0.2");
  const TString sType("#bf{charged jets}");


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn  = new TFile(sIn.Data(), "read");
  if (!fOut || !fIn) {
    cerr << "PANIC: couldn't open file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input matrix
  TH2D *hMatrix = (TH2D*) fIn -> Get(sMatrix.Data());
  if (!hMatrix) {
    cerr << "PANIC: couldn't grab input matrix!" << endl;
    return;
  }
  cout << "    Grabbed input matrix." << endl;


  // project and fit bins
  TF1     *fPtBins[NPtBins];
  TH1D    *hPtBins[NPtBins];
  for (UInt_t iPtBin = 0; iPtBin < NPtBins; iPtBin++) {

    // create name
    TString sPtHist("hPt");
    TString sPtFunc("fPt");
    sPtHist += PtBins[iPtBin];
    sPtFunc += PtBins[iPtBin];

    // create functions
    const Float_t guess0 = 1.;
    const Float_t guess1 = PtBins[iPtBin];
    const Float_t guess2 = 2. * FitLimits;
    const Float_t fitMin = PtBins[iPtBin] - FitLimits;
    const Float_t fitMax = PtBins[iPtBin] + FitLimits;
    fPtBins[iPtBin] = new TF1(sPtFunc.Data(), "gaus", fitMin, fitMax);
    fPtBins[iPtBin] -> SetParameter(0, guess0);
    fPtBins[iPtBin] -> SetParameter(1, guess1);
    fPtBins[iPtBin] -> SetParameter(2, guess2);
    fPtBins[iPtBin] -> SetLineStyle(fLin);
    fPtBins[iPtBin] -> SetLineColor(PtCols[iPtBin]);

     // project bin
    const UInt_t iMatrixBin = hMatrix -> GetYaxis() -> FindBin(PtBins[iPtBin]);
    hPtBins[iPtBin] = (TH1D*) hMatrix -> ProjectionX("", iMatrixBin, iMatrixBin, "") -> Clone();
    hPtBins[iPtBin] -> SetName(sPtHist.Data());
    hPtBins[iPtBin] -> Fit(sPtFunc.Data(), "QR0");
    hPtBins[iPtBin] -> GetFunction(sPtFunc.Data()) -> ResetBit(1<<9);

    // set styles
    hPtBins[iPtBin] -> SetLineStyle(fLin);
    hPtBins[iPtBin] -> SetLineColor(PtCols[iPtBin]);
    hPtBins[iPtBin] -> SetMarkerStyle(fMar);
    hPtBins[iPtBin] -> SetMarkerColor(PtCols[iPtBin]);
    hPtBins[iPtBin] -> SetTitle(sTitle.Data());
    hPtBins[iPtBin] -> SetTitleFont(fTxt);
    hPtBins[iPtBin] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPtBins[iPtBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtBins[iPtBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPtBins[iPtBin] -> GetXaxis() -> SetLabelFont(fTxt);
    hPtBins[iPtBin] -> GetXaxis() -> SetLabelSize(fLab);
    hPtBins[iPtBin] -> GetXaxis() -> CenterTitle(fCnt);
    hPtBins[iPtBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPtBins[iPtBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtBins[iPtBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPtBins[iPtBin] -> GetYaxis() -> SetLabelFont(fTxt);
    hPtBins[iPtBin] -> GetYaxis() -> SetLabelSize(fLab);
    hPtBins[iPtBin] -> GetYaxis() -> CenterTitle(fCnt);

  }
  cout << "    Projected pT bins." << endl;


  // make legend
  TLegend *lPtBins = new TLegend(fLeg1, fLeg1, fLeg2, fLeg2);
  lPtBins -> SetFillColor(fColL);
  lPtBins -> SetLineColor(fColL);
  lPtBins -> SetTextFont(fTxt);
  lPtBins -> SetTextAlign(fAln);
  for (UInt_t iPtBin = 0; iPtBin < NPtBins; iPtBin++) {

    const Double_t mu  = fPtBins[iPtBin] -> GetParameter(1);
    const Double_t sig = fPtBins[iPtBin] -> GetParameter(2);
    const Double_t rms = hPtBins[iPtBin] -> GetRMS();

    // display only a few decimals
    TString sMuRaw("");
    TString sSigRaw("");
    TString sRmsRaw("");
    sMuRaw  += mu;
    sSigRaw += sig;
    sRmsRaw += rms;

    const UInt_t nMuRaw  = sMuRaw.First(".");
    const UInt_t nSigRaw = sSigRaw.First(".");
    const UInt_t nRmsRaw = sRmsRaw.First(".");
    const UInt_t nMuTxt  = (nMuRaw + NDec) + 1;
    const UInt_t nSigTxt = (nSigRaw + NDec) + 1;
    const UInt_t nRmsTxt = (nRmsRaw + NDec) + 1;

    // create label
    TString sMuTxt("#mu = ");
    TString sSigTxt("#sigma = ");
    TString sRmsTxt("RMS = ");
    sMuTxt.Append(sMuRaw.Data(), nMuTxt);
    sSigTxt.Append(sSigRaw.Data(), nSigTxt);
    sRmsTxt.Append(sRmsRaw.Data(), nRmsTxt);

    TString sLabel("p_{T}^{par.} = ");
    sLabel += PtBins[iPtBin];
    sLabel += ": ";
    sLabel += sRmsTxt.Data();
    sLabel += ", ";
    sLabel += sMuTxt.Data();
    sLabel += ", ";
    sLabel += sSigTxt.Data();

    // add entry
    lPtBins -> AddEntry(hPtBins[iPtBin], sLabel.Data());

  }
  cout << "    Made legend." << endl;

  // make text
  TPaveText *pInfo = new TPaveText(fTxt1, fTxt1, fTxt2, fTxt2, "NDC NB");
  pInfo -> SetFillColor(fColL);
  pInfo -> SetLineColor(fColL);
  pInfo -> SetTextFont(fTxt);
  pInfo -> SetTextAlign(fAln);
  pInfo -> AddText(sSystem.Data());
  pInfo -> AddText(sTrigger.Data());
  pInfo -> AddText(sJets.Data());
  pInfo -> AddText(sType.Data());
  cout << "    Made text." << endl;


  // plot results
  TCanvas *cResolution = new TCanvas("cJetEnergyResolution", "", fSizX, fSizY);
  cResolution -> SetGrid(fGrid, fGrid);
  cResolution -> SetLogy(fLog);
  cResolution -> SetTopMargin(fMargin);
  cResolution -> SetRightMargin(fMargin);
  hPtBins[0]  -> Draw();
  for (UInt_t iPtBin = 1; iPtBin < NPtBins; iPtBin++) {
    hPtBins[iPtBin] -> Draw("sames");
  }
  lPtBins     -> Draw();
  pInfo       -> Draw();
  fOut        -> cd();
  cResolution -> Write();
  cResolution -> Close();
  cout << "    Made plot." << endl;


  // save and close histograms
  fOut -> cd();
  for (UInt_t iPtBin = 0; iPtBin < NPtBins; iPtBin++) {
    hPtBins[iPtBin] -> Write();
  }
  cout << "    Saved histograms." << endl;

  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished calculation!\n" << endl;

}

// End ------------------------------------------------------------------------
