// 'PlotResponseMatrixProjections.C'
// Derek Anderson
// 04.26.2021
// 
// Plots projections of x-axis (i.e. detector pT)
// from a specified response matrix, and then
// fits projections with a variety of functions
// (for smoothing).
//
// Crystal Ball Function from:
//   https://root-forum.cern.ch/t/crystalball-fitting-function/26973


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
//static const UInt_t  NPtBin(18);
static const UInt_t  NPtBin(12);
static const UInt_t  NVtx(4);
static const UInt_t  NPts(2);
static const UInt_t  NParGaus(3);
static const UInt_t  NParErf(6);
static const UInt_t  NParXtal(5);
static const Bool_t  DoXtalBallFit(false);
static const Float_t PtStart(0.2);
//static const Float_t PtStop(34.);
static const Float_t PtStop(38.);



void PlotResponseMatrixProjections() {

  // lower verbosity
  gStyle -> SetOptStat("emrf");
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning response matrix projections plot..." << endl;

  // io parameters
  const TString sIn("smearedPriorWithSmoothRes_forGammaUnfoldTest_pTbinBig.et1520gam.r05a065rm1chrg.d21m7y2021.root");
  const TString sOut("pi0TrgResponseProjections_smearWithSmoothRes_pTbinBig.et1520r05pi0.d21m7y2021.root");
  const TString sMatrix("hResponseInput");

  // histogram parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{det} [GeV/c]");
  const TString sTitleY("a. u.");
  const TString sPtPrefix("hProjPt");
  //const TString sPtSuffix[NPtBin] = {"0p2", "0p6", "1", "1p5", "2", "3", "4", "5", "6p5", "8", "10", "12", "14p5", "17", "20", "23", "26p5", "30"};
  const TString sPtSuffix[NPtBin] = {"0p2", "0p6", "1", "1p5", "2", "3", "5", "8", "12", "17", "23", "30"};
  //const UInt_t  fPtCol[NPtBin]    = {633, 809, 799, 401, 829, 819, 417, 849, 839, 433, 869, 859, 601, 889, 879, 617, 909, 899};
  const UInt_t  fPtCol[NPtBin]    = {809, 799, 829, 819, 849, 839, 869, 859, 889, 879, 909, 899};
  //const UInt_t  fPtMar[NPtBin]    = {24, 26, 32, 27, 25, 42, 28, 46, 30, 24, 26, 32, 27, 25, 42, 28, 46, 30};
  const UInt_t  fPtMar[NPtBin]    = {24, 26, 32, 27, 25, 42, 28, 46, 30, 24, 26, 32};
  const Float_t xPlot[NPts]       = {0., 38.};

  // legend parameters
  const TString sLegPrefix("p_{T}^{par} #in ");
  const TString sLegSuffix(" GeV/c");
  //const TString sLegBin[NPtBin] = {"(0.2, 0.6)", "(0.6, 1)", "(1, 1.5)", "(1.5, 2)", "(2, 3)", "(3, 4)", "(4, 5)", "(5, 6.5)", "(6.5, 8)", "(8, 10)", "(10, 12)", "(12, 14.5)", "(14.5, 17)", "(17, 20)", "(20, 23)", "(23, 26.5)", "(26.5, 30)", "(30, 34)"};
  const TString sLegBin[NPtBin] = {"(0.2, 0.6)", "(0.6, 1)", "(1, 1.5)", "(1.5, 2)", "(2, 3)", "(3, 5)", "(5, 8)", "(8, 12)", "(12, 17)", "(17, 23)", "(23, 30)", "(30, 38)"};

  // fit parameters
  const UInt_t  fWidFit(3);
  const UInt_t  fLinFit(9);
  const TString sGaus("gaus(0)");
  const TString sModErf("(TMath::Erf((x - [0]) / [1]) + 1.) * ([2] * TMath::Power(x, [3]) / (TMath::Power(1 + (x / [4]), [5])))");
  const TString sXtalBall("ROOT::Math::crystalball_function(x, [0], [1], [2], [3])");
  const TString sGausPrefix("hGausPt");
  const TString sErfPrefix("hErfPt");
  const TString sXtalPrefix("hXtalPt");
  const TString sGausFuncPrefix("fGausPt");
  const TString sErfFuncPrefix("fErfPt");
  const TString sXtalFuncPrefix("fXtalPt");
  const TString sParGaus[NParGaus]   = {"Amp", "#mu", "#sigma"};
  const TString sParErf[NParErf]     = {"#mu", "#sigma", "A", "B", "C", "D"};
  const TString sParXtal[NParXtal]   = {"Amp", "#alpha", "n", "#sigma", "#mu"};
  const Float_t pGuessErf[NParErf]   = {1., 1., 2.7, 5.5, 28.7, 34.3};
  const Float_t pGuessXtal[NParXtal] = {1., 2.7, 5.5, 1., 1.};
  const Float_t xFuncDef[NPts]       = {0., 57.};

  // text parameters
  const TString sSystem("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrigger("#pi^{0} trigger, E_{T}^{trg} #in (15, 20) GeV");
  const TString sJets("anti-k_{T}, R = 0.5");
  const TString sType("#bf{charged jets}");


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

  // grab response matrix
  TH2D *hMatrix = (TH2D*) fIn -> Get(sMatrix.Data());
  if (!hMatrix) {
    cerr << "PANIC: couldn't grab response matrix!\n" << endl;
    return;
  }
  cout << "    Grabbed response matrix." << endl;


  // make projections
  TH1D   *hPtBins[NPtBin];
  Bool_t hasEntries[NPtBin];
  UInt_t iPtHist(0);

  const UInt_t nPtParBin = hMatrix -> GetNbinsY();
  for (UInt_t iPtParBin = 1; iPtParBin < (nPtParBin + 1); iPtParBin++) {

    // check if bin is in specified range
    const Float_t pTparCnt    = hMatrix -> GetYaxis() -> GetBinCenter(iPtParBin);
    const Bool_t  isInPtRange = ((pTparCnt > PtStart) && (pTparCnt < PtStop));
    if (!isInPtRange) continue;

    // create histogram name
    TString sPtName(sPtPrefix.Data());
    sPtName.Append(sPtSuffix[iPtHist]);

    // project pTdet distribution
    hPtBins[iPtHist] = (TH1D*) hMatrix -> ProjectionX("", iPtParBin, iPtParBin, "") -> Clone();
    hPtBins[iPtHist] -> SetName(sPtName.Data());

    // check if projection has entries
    const Double_t intPtBin  = hPtBins[iPtHist] -> Integral();
    const Bool_t   intIsNot0 = (intPtBin > 0.);

    // if not empty, normalize projection
    hasEntries[iPtHist] = intIsNot0;
    if (hasEntries[iPtHist]) {
      hPtBins[iPtHist] -> Scale(1. / intPtBin);
    }
    iPtHist++;

  }  // end pTpar bin loop
  cout << "    Projected pTdet distributions." << endl;

  // set histogram styles
  const UInt_t  fLin(1);
  const UInt_t  fFil(0);
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.1);
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
    hPtBins[iPtBin] -> SetMarkerColor(fPtCol[iPtBin]);
    hPtBins[iPtBin] -> SetMarkerStyle(fPtMar[iPtBin]);
    hPtBins[iPtBin] -> SetLineColor(fPtCol[iPtBin]);
    hPtBins[iPtBin] -> SetLineStyle(fLin);
    hPtBins[iPtBin] -> SetFillColor(fPtCol[iPtBin]);
    hPtBins[iPtBin] -> SetFillStyle(fFil);
    hPtBins[iPtBin] -> SetTitle(sTitle.Data());
    hPtBins[iPtBin] -> SetTitleFont(fTxt);
    hPtBins[iPtBin] -> GetXaxis() -> SetLabelSize(fLab);
    hPtBins[iPtBin] -> GetXaxis() -> SetLabelFont(fTxt);
    hPtBins[iPtBin] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPtBins[iPtBin] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtBins[iPtBin] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPtBins[iPtBin] -> GetXaxis() -> CenterTitle(fCnt);
    hPtBins[iPtBin] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hPtBins[iPtBin] -> GetYaxis() -> SetLabelSize(fLab);
    hPtBins[iPtBin] -> GetYaxis() -> SetLabelFont(fTxt);
    hPtBins[iPtBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPtBins[iPtBin] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtBins[iPtBin] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPtBins[iPtBin] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end pTbin loop
  cout << "    Set histogram styles." << endl;


  // fit projections
  TF1      *fPtGaus[NPtBin];
  TF1      *fPtErf[NPtBin];
  TF1      *fPtXtal[NPtBin];
  TH1D     *hPtGaus[NPtBin];
  TH1D     *hPtErf[NPtBin];
  TH1D     *hPtXtal[NPtBin];
  cout << "    Fitting histograms:" << endl;

  Double_t parGausPt[NParGaus];
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {

    cout << "      Fitting histogram " << (iPtBin + 1) << "/" << NPtBin << "..." << endl;

    // create histogram names
    TString sPtGaus(sGausPrefix.Data());
    TString sPtErf(sErfPrefix.Data());
    TString sPtXtal(sXtalPrefix.Data());
    sPtGaus.Append(sPtSuffix[iPtBin]);
    sPtErf.Append(sPtSuffix[iPtBin]);
    sPtXtal.Append(sPtSuffix[iPtBin]);

    // create function names
    TString sFuncGaus(sGausFuncPrefix.Data());
    TString sFuncErf(sErfFuncPrefix.Data());
    TString sFuncXtal(sXtalFuncPrefix.Data());
    sFuncGaus.Append(sPtSuffix[iPtBin]);
    sFuncErf.Append(sPtSuffix[iPtBin]);
    sFuncXtal.Append(sPtSuffix[iPtBin]);

    // create histograms to fit
    hPtGaus[iPtBin] = (TH1D*) hPtBins[iPtBin] -> Clone();
    hPtErf[iPtBin]  = (TH1D*) hPtBins[iPtBin] -> Clone();
    hPtXtal[iPtBin] = (TH1D*) hPtBins[iPtBin] -> Clone();
    hPtGaus[iPtBin] -> SetName(sPtGaus.Data());
    hPtErf[iPtBin]  -> SetName(sPtErf.Data());
    hPtXtal[iPtBin] -> SetName(sPtXtal.Data());

    // extract info for fits
    const UInt_t   iStart = hPtGaus[iPtBin] -> FindFirstBinAbove(0.);
    const UInt_t   iStop  = hPtGaus[iPtBin] -> FindLastBinAbove(0.);
    const Double_t xStart = hPtGaus[iPtBin] -> GetBinCenter(iStart);
    const Double_t xStop  = hPtGaus[iPtBin] -> GetBinCenter(iStop);
    parGausPt[0] = hPtGaus[iPtBin] -> GetMaximum();
    parGausPt[1] = hPtGaus[iPtBin] -> GetMean();
    parGausPt[2] = hPtGaus[iPtBin] -> GetRMS();

    // create gaussian functions and fit
    fPtGaus[iPtBin] = new TF1(sFuncGaus.Data(), sGaus.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
      fPtGaus[iPtBin] -> SetParameter(iParGaus, parGausPt[iParGaus]);
      fPtGaus[iPtBin] -> SetParName(iParGaus, sParGaus[iParGaus]);
      fPtGaus[iPtBin] -> SetLineColor(fPtCol[iPtBin]);
      fPtGaus[iPtBin] -> SetLineStyle(fLinFit);
      fPtGaus[iPtBin] -> SetLineWidth(fWidFit);
    }  // end parameter loop
    hPtGaus[iPtBin] -> Fit(sFuncGaus.Data(), "", "", xStart, xStop);

    // create error functions and fit
    fPtErf[iPtBin] = new TF1(sFuncErf.Data(), sModErf.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParErf = 0; iParErf < NParErf; iParErf++) {
      fPtErf[iPtBin] -> SetLineColor(fPtCol[iPtBin]);
      fPtErf[iPtBin] -> SetLineStyle(fLinFit);
      fPtErf[iPtBin] -> SetLineWidth(fWidFit);
      fPtErf[iPtBin] -> SetParName(iParErf, sParErf[iParErf]);
      switch (iParErf) {
        case 0:
          fPtErf[iPtBin] -> FixParameter(iParErf, fPtGaus[iPtBin] -> GetParameter(1));
          break;
        case 1:
          fPtErf[iPtBin] -> FixParameter(iParErf, fPtGaus[iPtBin] -> GetParameter(2));
          break;
        default:
          fPtErf[iPtBin] -> SetParameter(iParErf, pGuessErf[iParErf]);
          break;
      }
    }  // end parameter loop
    hPtErf[iPtBin] -> Fit(sFuncErf.Data(), "B", "", xStart, xStop);

    // create xtal ball functions and fit
    if (DoXtalBallFit) {
      fPtXtal[iPtBin] = new TF1(sFuncXtal.Data(), sXtalBall.Data(), xFuncDef[0], xFuncDef[1]);
      for (UInt_t iParXtal = 0; iParXtal < NParXtal; iParXtal++) {
        fPtXtal[iPtBin] -> SetLineColor(fPtCol[iPtBin]);
        fPtXtal[iPtBin] -> SetLineStyle(fLinFit);
        fPtXtal[iPtBin] -> SetLineWidth(fWidFit);
        fPtXtal[iPtBin] -> SetParName(iParXtal, sParXtal[iParXtal]);
        switch (iParXtal) {
          case 4:
            fPtXtal[iPtBin] -> FixParameter(iParXtal, fPtGaus[iPtBin] -> GetParameter(2));
            break;
          case 5:
            fPtXtal[iPtBin] -> FixParameter(iParXtal, fPtGaus[iPtBin] -> GetParameter(1));
            break;
          default:
            fPtXtal[iPtBin] -> SetParameter(iParXtal, pGuessXtal[iPtBin]);
            break;
        }
      }  // end parameter loop
      hPtXtal[iPtBin] -> Fit(sFuncXtal.Data(), "B", "", xStart, xStop);
    }

  }  // end pTbin loop
  cout << "    Fit histograms." << endl;


  // make legends
  const UInt_t  fColL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fAlnL(12);
  const UInt_t  fAlnT(32);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  // create labels
  TString sPtLabel[NPtBin];
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
    sPtLabel[iPtBin] = sLegPrefix;
    sPtLabel[iPtBin].Append(sLegBin[iPtBin].Data());
    sPtLabel[iPtBin].Append(sLegSuffix.Data());
  }  // end pTbin loop

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  leg -> SetFillColor(fColL);
  leg -> SetFillStyle(fFilL);
  leg -> SetLineColor(fColL);
  leg -> SetLineStyle(fLinL);
  leg -> SetTextAlign(fAlnL);
  leg -> SetTextFont(fTxt);
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
    if (hasEntries[iPtBin]) leg -> AddEntry(hPtBins[iPtBin], sPtLabel[iPtBin].Data());
  }  // end pTbin loop

  TLegend *legG = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legG -> SetFillColor(fColL);
  legG -> SetFillStyle(fFilL);
  legG -> SetLineColor(fColL);
  legG -> SetLineStyle(fLinL);
  legG -> SetTextAlign(fAlnL);
  legG -> SetTextFont(fTxt);
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
    if (hasEntries[iPtBin]) legG -> AddEntry(hPtGaus[iPtBin], sPtLabel[iPtBin].Data());
  }  // end pTbin loop

  TLegend *legE = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legE -> SetFillColor(fColL);
  legE -> SetFillStyle(fFilL);
  legE -> SetLineColor(fColL);
  legE -> SetLineStyle(fLinL);
  legE -> SetTextAlign(fAlnL);
  legE -> SetTextFont(fTxt);
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
    if (hasEntries[iPtBin]) legE -> AddEntry(hPtErf[iPtBin], sPtLabel[iPtBin].Data());
  }  // end pTbin loop

  TLegend *legC = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  if (DoXtalBallFit) {
    legC -> SetFillColor(fColL);
    legC -> SetFillStyle(fFilL);
    legC -> SetLineColor(fColL);
    legC -> SetLineStyle(fLinL);
    legC -> SetTextAlign(fAlnL);
    legC -> SetTextFont(fTxt);
    for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
      if (hasEntries[iPtBin]) legC -> AddEntry(hPtXtal[iPtBin], sPtLabel[iPtBin].Data());
    }  // end pTbin loop
  }
  cout << "    Made legends." << endl;


  // make text box
  TPaveText *txt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  txt -> SetFillColor(fColL);
  txt -> SetFillStyle(fFilL);
  txt -> SetLineColor(fColL);
  txt -> SetLineStyle(fLinL);
  txt -> SetTextAlign(fAlnT);
  txt -> SetTextFont(fTxt);
  txt -> AddText(sSystem.Data());
  txt -> AddText(sTrigger.Data());
  txt -> AddText(sJets.Data());
  txt -> AddText(sType.Data());
  cout << "    Made text box." << endl;


  // make plot
  const UInt_t  fSizX(750);
  const UInt_t  fSizY(750);
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const Float_t fMargT(0.02);
  const Float_t fMargR(0.02);

  TCanvas *cPtProj = new TCanvas("cPtProj", "", fSizX, fSizY);
  cPtProj -> SetTopMargin(fMargT);
  cPtProj -> SetRightMargin(fMargR);
  cPtProj -> SetLogy(fLog);
  cPtProj -> SetGrid(fGrid, fGrid);
  cPtProj -> cd();
  if (hasEntries[0])
    hPtBins[0] -> Draw("hist p l");
  for (UInt_t iPtBin = 1; iPtBin < NPtBin; iPtBin++) {
    if (hasEntries[iPtBin])
      hPtBins[iPtBin] -> Draw("hist p l same");
  }
  leg     -> Draw();
  txt     -> Draw();
  fOut    -> cd();
  cPtProj -> Write();
  cPtProj -> Close();

  TCanvas *cPtGaus = new TCanvas("cPtGaus", "", fSizX, fSizY);
  cPtGaus -> SetTopMargin(fMargT);
  cPtGaus -> SetRightMargin(fMargR);
  cPtGaus -> SetLogy(fLog);
  cPtGaus -> SetGrid(fGrid, fGrid);
  cPtGaus -> cd();
  if (hasEntries[0])
    hPtGaus[0] -> Draw("l");
  for (UInt_t iPtBin = 1; iPtBin < NPtBin; iPtBin++) {
    if (hasEntries[iPtBin])
      hPtGaus[iPtBin] -> Draw("l sames");
  }
  legG    -> Draw();
  txt     -> Draw();
  fOut    -> cd();
  cPtGaus -> Write();
  cPtGaus -> Close();

  TCanvas *cPtErf = new TCanvas("cPtErf", "", fSizX, fSizY);
  cPtErf -> SetTopMargin(fMargT);
  cPtErf -> SetRightMargin(fMargR);
  cPtErf -> SetLogy(fLog);
  cPtErf -> SetGrid(fGrid, fGrid);
  cPtErf -> cd();
  if (hasEntries[0])
    hPtErf[0] -> Draw("l");
  for (UInt_t iPtBin = 1; iPtBin < NPtBin; iPtBin++) {
    if (hasEntries[iPtBin])
      hPtErf[iPtBin] -> Draw("l sames");
  }
  legE   -> Draw();
  txt    -> Draw();
  fOut   -> cd();
  cPtErf -> Write();
  cPtErf -> Close();

  if (DoXtalBallFit) {
    TCanvas *cPtXtal = new TCanvas("cPtXtal", "", fSizX, fSizY);
    cPtXtal -> SetTopMargin(fMargT);
    cPtXtal -> SetRightMargin(fMargR);
    cPtXtal -> SetLogy(fLog);
    cPtXtal -> SetGrid(fGrid, fGrid);
    cPtXtal -> cd();
    if (hasEntries[0])
      hPtXtal[0] -> Draw("l");
    for (UInt_t iPtBin = 1; iPtBin < NPtBin; iPtBin++) {
      if (hasEntries[iPtBin])
        hPtXtal[iPtBin] -> Draw("l sames");
    }
    legC    -> Draw();
    txt     -> Draw();
    fOut    -> cd();
    cPtXtal -> Write();
    cPtXtal -> Close();
  }
  cout << "    Made plots." << endl;


  // save histograms
  fOut -> cd();
  for (UInt_t iPtBin = 0; iPtBin < NPtBin; iPtBin++) {
    hPtBins[iPtBin] -> Write();
    hPtGaus[iPtBin] -> Write();
    hPtErf[iPtBin]  -> Write();
    if (DoXtalBallFit) hPtXtal[iPtBin] -> Write();
  }
  hMatrix -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished plotting response matrix projection!\n" << endl;

}  // end 'PlotResponseMatrixProjections()'



Double_t EvalCrystalBallFunction(Double_t x, Double_t alpha, Double_t n, Double_t sigma, Double_t mean) {

  // check if width is negative
  if (sigma < 0.) return 0.;

  Double_t z = (x - mean) / sigma; 
  if (alpha < 0) z = -z; 

  // calculate value
  Double_t abs_alpha = TMath::Abs(alpha);
  if (z  > (-1. * abs_alpha)) {
    return TMath::Exp(-0.5 * z * z);
  } 
  else {
    Double_t nDivAlpha = n / abs_alpha;
    Double_t AA        =  TMath::Exp(-0.5 * abs_alpha * abs_alpha);
    Double_t B         = nDivAlpha - abs_alpha;
    Double_t arg       = nDivAlpha /(B - z);
    Double_t returnVal = AA * TMath::Power(arg, n);
    return returnVal;
  }

}  // end 'EvalCrystalBallFunction(Double_t, Double_t, Double_t, Double_t, Double_t)'



Double_t CrystalBallFunction(const Double_t *x, const Double_t *par) {

  const Double_t xtal = p[0] * EvalCrystalBallFunction(x[0], p[1], p[2], p[3], p[4]);
  return xtal;

}  // end 'CrystalBallFunction(Double_t*, Double_t*)'

// End ------------------------------------------------------------------------
