// 'PriorExponentFitTest.C'
// Derek Anderson
// 04.21.2021
//
// Fits some exponentials merged
// with a hypertangent to a
// provided prior.

#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

// global constants
//static const UInt_t NExp(2);  // 9 - 11 GeV
static const UInt_t NExp(3);  // 11 - 15 GeV
//static const UInt_t NExp(4);  // 15 - 20 GeV
//static const UInt_t NTan(1);  // 9 - 11 GeV
static const UInt_t NTan(2);  // 11 - 15, 15 - 20 GeV
static const UInt_t NPts(2);
static const UInt_t NPar(2);



void PriorExponentFitTest() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning exponential fit test..." << endl;

  // io parameters
  const TString sIn("input/embed/pp200r9embed.forUnfolding_pTbinGam11.et1115pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.root");
  const TString sOut("priorExponentFitTest_twoTanhTest_twoFitCallAndNewBin_pTbinFine.et1115r05pi0.d30m8y2021.root");
  const TString sHistIn("hSumParAll");

  // fit parameters
  const TString sTotName("fTotalFit");
  const UInt_t  fColTot(6);

/*
  // fit parameters [9 - 11 GeV]
  const TString sTotFunc("(expo(0) + expo(2)) * tanh((x-[4])/[5])");
  const TString sExpNames[NExp]        = {"fExpSoft", "fExpHard"};
  const TString sExpFuncs[NExp]        = {"expo(0)", "expo(0)"};
  const Float_t xExpRange[NExp][NPts]  = {{0.6, 2.}, {2., 40.}};
  const Float_t pExpGuess[NExp][NPar]  = {{0.25, -1.03}, {-1.76, -0.26}};
  const Float_t pTanhGuess[NTan][NPar] = {2., 1.};
  const UInt_t  fColExp[NExp]          = {2, 4};
*/

  // fit parameters [11 - 15 GeV]
  const TString sTotFunc("(expo(0) + expo(2) + expo(4)) * (tanh((x-[6])/[7]) + tanh((x-[8])/[9]))");
  const TString sExpNames[NExp]        = {"fExpSoft", "fExpMed", "fExpHard"};
  const TString sExpFuncs[NExp]        = {"expo(0)", "expo(0)", "expo(0)"};
  const Float_t xExpRange[NExp][NPts]  = {{0.6, 3.}, {3., 12.}, {12., 40.}};
  const Float_t pExpGuess[NExp][NPar]  = {{0.19, -1.01}, {-2.44, -0.14}, {-1.06, -0.25}};
  const Float_t pTanhGuess[NTan][NPar] = {{3., 1.}, {12., 1.}};
  const UInt_t  fColExp[NExp]          = {2, 3, 4};

/*
  // fit parameters [15 - 20 GeV]
  const TString sTotFunc("(expo(0) + expo(2) + expo(4) + expo(6)) * (tanh((x-[8])/[9]) + tanh((x-[10])/[11]))");
  const TString sExpNames[NExp]        = {"fExpSoft", "fExpMedSoft", "fExpMedHard", "fExpHard"};
  const TString sExpFuncs[NExp]        = {"expo(0)", "expo(0)", "expo(0)", "expo(0)"};
  const Float_t xExpRange[NExp][NPts]  = {{0.2, 2.}, {2., 5.}, {5., 14.5}, {14.5, 40.}};
  const Float_t pExpGuess[NExp][NPar]  = {{0.31, -1.18}, {-1.36, -0.39}, {-3.81, -0.01}, {-0.49, -0.22}};
  const Float_t pTanhGuess[NTan][NPar] = {{7.5, 13.}, {7.5, -13.}};
  const UInt_t  fColExp[NExp]          = {2, 3, 4, 6};
*/

  // histogram parameters
  const TString sFitHistExp("hFitPrior_onlyExpos");
  const TString sFitHistTot("hFitPrior_totalFit");
  const TString sTitle("");
  const TString sTitleX("p_{T}^{par} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{2}N^{jet}/d(p_{T}^{par} #eta^{jet}) [GeV/c]^{-1}");
  const Float_t fOffHistX(1.0);
  const Float_t fOffHistY(1.05);
  const Float_t fLabHist(0.03);
  const UInt_t  fColHist(1);
  const UInt_t  fMarHist(8);
  const UInt_t  fLinHist(1);
  const UInt_t  fFilHist(0);
  const UInt_t  fTxtHist(42);
  const UInt_t  fCntHist(1);

  // canvas parameters
  const TString sCanExp("cOnlyExpos");
  const TString sCanTot("cTotalFit");
  const Float_t fWidth(750);
  const Float_t fHeight(750);
  const Float_t fTopMarg(0.02);
  const Float_t fRightMarg(0.02);
  const UInt_t  fLogY(1);
  const UInt_t  fGrid(0);

  // open files
  TFile *fIn  = new TFile(sIn.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn || !fOut) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fIn = " << fIn << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histogram
  TH1D *hPrior = (TH1D*) fIn -> Get(sHistIn.Data());
  if (!hPrior) {
    cerr << "PANIC: couldn't grab input histogram!\n" << endl;
    return;
  }
  cout << "    Grabbed input histogram." << endl;

  // set histogram style
  hPrior -> SetMarkerColor(fColHist);
  hPrior -> SetMarkerStyle(fMarHist);
  hPrior -> SetLineColor(fColHist);
  hPrior -> SetLineStyle(fLinHist);
  hPrior -> SetFillColor(fColHist);
  hPrior -> SetFillStyle(fFilHist);
  hPrior -> SetTitle(sTitle.Data());
  hPrior -> SetTitleFont(fTxtHist);
  hPrior -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPrior -> GetXaxis() -> SetTitleFont(fTxtHist);
  hPrior -> GetXaxis() -> SetTitleOffset(fOffHistX);
  hPrior -> GetXaxis() -> SetLabelSize(fLabHist);
  hPrior -> GetXaxis() -> SetLabelFont(fTxtHist);
  hPrior -> GetXaxis() -> CenterTitle(fCntHist);
  hPrior -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPrior -> GetYaxis() -> SetTitleOffset(fOffHistY);
  hPrior -> GetYaxis() -> SetTitleFont(fTxtHist);
  hPrior -> GetYaxis() -> SetLabelSize(fLabHist);
  hPrior -> GetYaxis() -> SetLabelFont(fTxtHist);
  hPrior -> GetYaxis() -> CenterTitle(fCntHist);
  cout << "    Set histogram style." << endl;

  // declare histograms for fitting
  TH1D *hFitExp = (TH1D*) hPrior -> Clone();
  TH1D *hFitTot = (TH1D*) hPrior -> Clone();
  hFitExp -> SetName(sFitHistExp.Data());
  hFitTot -> SetName(sFitHistTot.Data());
  cout << "    Declared histograms for fitting." << endl;

  // initialize 1st set of exponentials
  TF1 *fExp[NExp];
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    fExp[iExp] = new TF1(sExpNames[iExp].Data(), sExpFuncs[iExp].Data(), xExpRange[iExp][0], xExpRange[iExp][1]);
    fExp[iExp] -> SetParameter(0, pExpGuess[iExp][0]);
    fExp[iExp] -> SetParameter(1, pExpGuess[iExp][1]);
    fExp[iExp] -> SetLineColor(fColExp[iExp]);
  }
  cout << "    Initialized individual exponentials." << endl;

  // determine exponential parameters
  Double_t pExpFit[NExp][NPar];
  Double_t eExpFit[NExp][NPar];
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    if (iExp == 0) {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R");
    } else {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R+");
    }
    pExpFit[iExp][0] = fExp[iExp] -> GetParameter(0);
    eExpFit[iExp][0] = fExp[iExp] -> GetParError(0);
    pExpFit[iExp][1] = fExp[iExp] -> GetParameter(1);
    eExpFit[iExp][1] = fExp[iExp] -> GetParError(1);
  }

  // announce parameters
  cout << "    Extracted exponential parameters:" << endl;
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    cout << "      " << fExp[iExp] -> GetName() << ": const. = " << pExpFit[iExp][0] << " +- " << eExpFit[iExp][0]
         << ", slope = " << pExpFit[iExp][1] << " +- " << eExpFit[iExp][1]
         << endl;
  }

  // initialize total fit function
  TF1 *fTot = new TF1(sTotName.Data(), sTotFunc.Data(), xExpRange[0][0], xExpRange[NExp - 1][1]);
  fTot -> SetLineColor(fColTot);

  const UInt_t iTanMin = 2 * NExp;
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    const UInt_t  iParExp = 2 * iExp;
    //fTot -> FixParameter(iParExp + 0, pExpFit[iExp][0]);  // 9 - 11 GeV
    fTot -> SetParameter(iParExp + 0, pExpFit[iExp][0]);  // 11 - 15, 15 - 20 GeV
    fTot -> FixParameter(iParExp + 1, pExpFit[iExp][1]);
  }
  for (UInt_t iTan = 0; iTan < NTan; iTan++) {
    const UInt_t iParTan = iTanMin + (2 * iTan);
    fTot -> SetParameter(iParTan + 0, pTanhGuess[iTan][0]);
    fTot -> SetParameter(iParTan + 1, pTanhGuess[iTan][1]);
  }
  cout << "    Initialized total fit function." << endl;

/*
  // set tanh names [9 - 11 GeV]
  fTot -> SetParName(4, "p_{T}^{0}");
  fTot -> SetParName(5, "a_{0}");
*/

  // set tanh names [11 - 15 GeV]
  fTot -> SetParName(6, "p_{T}^{0}");
  fTot -> SetParName(7, "a_{0}");
  fTot -> SetParName(8, "p_{T}^{1}");
  fTot -> SetParName(9, "a_{1}");

/*
  // set tanh names [15 - 20 GeV]
  fTot -> SetParName(8, "p_{T}^{0}");
  fTot -> SetParName(9, "a_{0}");
  fTot -> SetParName(10, "p_{T}^{1}");
  fTot -> SetParName(11, "a_{1}");
*/
  cout << "    Set parameter names." << endl;

  // fit function
  Int_t fitStat(-1);
  fitStat = (Int_t) hFitTot -> Fit(sTotName.Data(), "RB");
  fitStat = (Int_t) hFitTot -> Fit(sTotName.Data(), "RB");

  // 11 - 15 GeV
  if (fitStat != 0) {
    cout << "    11 - 15 GeV fit failed, using default fit values." << endl;
    fTot    -> SetParameter(0, 3.81);
    fTot    -> SetParameter(2, 1.03);
    fTot    -> SetParameter(4, -14.53);
    fTot    -> FixParameter(6, 7.14);
    fTot    -> FixParameter(7, 13.29);
    fTot    -> FixParameter(8, 7.33);
    fTot    -> FixParameter(9, -13.34);
    hFitTot -> Fit(sTotName.Data(), "RB");
  }

  cout << "    Fit prior with total function:" << endl;
  for (UInt_t iTan = 0; iTan < NTan; iTan++) {
    const UInt_t iParTan = (2 * NExp) + (2 * iTan);
    cout << "      " << fTot -> GetParName(iParTan) << " = " << fTot -> GetParameter(iParTan) << ", " << fTot -> GetParName(iParTan + 1) << " = " << fTot -> GetParameter(iParTan + 1)
         << endl;
  }

  // make canvases
  TCanvas *cExp = new TCanvas(sCanExp.Data(), "", fWidth, fHeight);
  cExp    -> SetTopMargin(fTopMarg);
  cExp    -> SetRightMargin(fRightMarg);
  cExp    -> SetGrid(fGrid, fGrid);
  cExp    -> SetLogy(fLogY);
  cExp    -> cd();
  hFitExp -> Draw();
  fOut    -> cd();
  cExp    -> Write();
  cExp    -> Close();

  TCanvas *cTot = new TCanvas(sCanTot.Data(), "", fWidth, fHeight);
  cTot    -> SetTopMargin(fTopMarg);
  cTot    -> SetRightMargin(fRightMarg);
  cTot    -> SetGrid(fGrid, fGrid);
  cTot    -> SetLogy(fLogY);
  cTot    -> cd();
  hFitTot -> Draw();
  fOut    -> cd();
  cTot    -> Write();
  cTot    -> Close();
  cout << "    Made canvases." << endl;

  // save histograms and functions
  fOut    -> cd();
  hPrior  -> Write();
  hFitExp -> Write();
  hFitTot -> Write();
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    fExp[iExp] -> Write();
  }
  fTot -> Write();
  cout << "    Saved histograms and functions." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished exponential fit test!\n" << endl;

}

// End ------------------------------------------------------------------------
