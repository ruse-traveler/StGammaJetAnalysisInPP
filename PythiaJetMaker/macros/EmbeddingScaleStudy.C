// 'EmbeddingScaleStudy.C'
// Derek Anderson
// 09.07.2017 (v3)
// 
// This takes a set histograms, and
// using a set of known fudge factors,
// x-sections, and no. of events scale
// each histogram accordingly and sum
// them.
//
// NOTE: don't forget to change the io-
// parameters!

#include <vector>
#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const UInt_t nBins(8);
static const UInt_t nHist(nBins + 1);
static const UInt_t nTotal(11);
// scales
static const Double_t fudge[nTotal] = {1.228, 1.051, 1.014, 1., 1., 1., 1., 1., 1., 1., 1.};
static const Double_t xSctn[nTotal] = {9.005805, 1.461907, 0.3544354, 0.1513760, 2.488644e-02, 5.845846e-03, 2.304880e-03, 3.426618e-04, 4.562988e-05, 9.738044e-06, 5.019977e-07};
static const Double_t nEvts[nTotal] = {2100295., 600300., 600300., 300289., 300289., 300289., 160295., 100302., 80293., 76303., 23307.};



void EmbeddingScaleStudy() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Starting scale study..." << endl;


  // io-parameters
  const TString sInput("dataXembedding.eTtrg.pi0vsHighTwrs.d3m9y2017.root");
  const TString sOutput("dataXembedding.eTtrg.pi0vsHighTwrs.scaledDistributions.withHistograms.d7m9y2017.root");
  const TString sHist[nHist]     = {"hTrgEt_pt57", "hTrgEt_pt79", "hTrgEt_pt911", "hTrgEt_pt1115", "hTrgEt_pt1520", "hTrgEt_pt2025", "hTrgEt_pt2535", "hTrgEt_pt35-1", "hTrgEt_data"};
  const TString sBins[nBins]     = {"(5,7)", "(7,9)", "(9,11)", "(11,15)", "(15,20)", "(20,25)", "(25,35)", "(35,-1)"};
  const TString sLegs[nHist + 1] = {"p_{T}^{part}#in(5,7)", "p_{T}^{part}#in(7,9)", "p_{T}^{part}#in(9,11)", "p_{T}^{part}#in(11,15)", "p_{T}^{part}#in(15,20)", "p_{T}^{part}#in(20,25)", "p_{T}^{part}#in(25,35)", "p_{T}^{part}#in(35,-1)", "sum", "data"};
  const TString sTitle("Trigger E_{T}");
  const TString sXtitle("E_{T}^{trg}");
  const TString sYtitle("dN^{trg}/dE_{T}^{trg}");
  const TString sYtitleR("embedding / data");

  // constants
  const UInt_t  nDec(2);
  const UInt_t  nRebin(5);
  const UInt_t  wSize(775);
  const UInt_t  hSize(775);
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const UInt_t  fColLbl(0);
  const UInt_t  fColLne(1);
  const UInt_t  fStyLbl(0);
  const UInt_t  fStyLne(2);
  const UInt_t  fGrid(0);
  const UInt_t  fLog(0);
  const UInt_t  fCol[nHist + 1] = {810, 800, 830, 850, 870, 860, 890, 880, 910, 1};
  const UInt_t  fMar[nHist + 1] = {7, 24, 26, 32, 27, 25, 28, 30, 24, 7};
  const Bool_t  doRebin(false);
  const Float_t fLbl(0.02);
  const Float_t fMargin(0.);
  const Float_t weight(1.);
  const Float_t xRange[2] = {9., 20.};
  const Float_t pRange[2] = {5., 25.};
  const Float_t xTxt[2]   = {0.1, 0.3};
  const Float_t xLeg[2]   = {0.3, 0.5};
  const Float_t yTxt[2]   = {0.1, 0.3};
  const Float_t yLeg[2]   = {0.1, 0.3};
  const Float_t xPadR[2]  = {0., 1.};
  const Float_t xPadD[2]  = {0., 1.};
  const Float_t yPadR[2]  = {0., 0.35};
  const Float_t yPadD[2]  = {0.35, 1.};
  const TString sCan("cScaledDistributions");
  const TString sPadD("pDistributions");
  const TString sPadR("pRatios");


  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput  = new TFile(sInput.Data(), "read");
  if (!fInput) {
    cerr << "PANIC: couldn't open input file!\n" << endl;
    return;
  }
  cout << "    Files opened." << endl;

  // grab histograms
  TH1D *hHist[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] = (TH1D*) fInput -> Get(sHist[iHist].Data());
    if (!hHist[iHist]) {
      cerr << "PANIC: couldn't grab histogram no. " << iHist << "!\n" << endl;
      return;
    }
  }
  fOutput -> cd();
  cout << "    Histograms grabbed." << endl;


  // rebin histograms (if necessary)
  if (doRebin) {
    for (UInt_t iHist = 0; iHist < nHist; iHist++) {
      hHist[iHist] -> Rebin(nRebin);
    }
    cout << "    Histograms rebinned." << endl;
  }


  // calculate scales
  Double_t weights[nBins];

  const UInt_t iStart = nTotal - nBins;
  for (UInt_t iBins = iStart; iBins < nTotal; iBins++) {
    const Double_t binLumi   = nEvts[iBins] / xSctn[iBins];
    const Double_t binWeight = 1. / (fudge[iBins] * binLumi);
    weights[iBins - iStart]  = binWeight;
  }
  cout << "    Weights calculated." << endl;


  // sum histograms
  const UInt_t  nBinsX = hHist[nBins] -> GetNbinsX();
  const Float_t xBin0  = hHist[nBins] -> GetBinLowEdge(1);
  const Float_t xBin1  = hHist[nBins] -> GetBinLowEdge(nBinsX + 1);

  TH1D *hSumU = new TH1D("hSumUnormalized", "", nBinsX, xBin0, xBin1);
  TH1D *hSumN = new TH1D("hSumNormalized", "", nBinsX, xBin0, xBin1);
  hSumU -> Sumw2();
  hSumN -> Sumw2();
  for (UInt_t iBins = 0; iBins < nBins; iBins++) {
    hSumU -> Add(hHist[iBins], weights[iBins]);
    hSumN -> Add(hHist[iBins], weights[iBins]);
  }
  cout << "    Histograms summed." << endl;

  // calculate normalization
  const Double_t intSum   = hSumN        -> Integral();
  const Double_t intData  = hHist[nBins] -> Integral();
  const Double_t normalization = intData / intSum;
  hSumN -> Scale(normalization);
  cout << "    Histograms normalized." << endl;


  // calculate chi2
  UInt_t   nChi  = 0;
  UInt_t   iChi0 = hHist[nBins] -> FindBin(xRange[0]);
  UInt_t   iChi1 = hHist[nBins] -> FindBin(xRange[1]);
  Double_t chi2  = 0.;
  for (UInt_t iBinX = iChi0; iBinX <= iChi1; iBinX++) {
    const Double_t xData = hHist[nBins] -> GetBinCenter(iBinX);
    const Double_t yData = hHist[nBins] -> GetBinContent(iBinX);
    const Double_t eData = hHist[nBins] -> GetBinError(iBinX);

    const UInt_t   iSumX = hSumN -> FindBin(xData);
    const Double_t ySum  = hSumN -> GetBinContent(iSumX);
    const Double_t eSum  = hSumN -> GetBinError(iSumX);
    if ((yData < 0.) || (ySum < 0.))
      continue;

    if ((eData > 0.) && (eSum > 0.)) {
      const Double_t dErr = TMath::Sqrt((eSum * eSum) + (eData * eData));
      const Double_t num  = TMath::Power(yData - ySum, 2.);
      const Double_t den  = TMath::Power(dErr, 2.);
      const Double_t c2   = num / den;
      chi2 += c2;
      nChi++;
    }
  }  // end bin loop
  if (nChi > 0) chi2 /= (Double_t) nChi;
  cout << "    Chi2 calculated." << endl;


  // scale distributions and calculate ratios
  TH1D *hRatios[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    TString sRatioName("hRatio");
    sRatioName += iHist;

    hRatios[iHist] = new TH1D(sRatioName.Data(), "", nBinsX, xBin0, xBin1);
    hRatios[iHist] -> Sumw2();
    if (iHist < nBins) {
      hRatios[iHist] -> Divide(hHist[iHist], hHist[nBins], weight, weight);
    }
    else if (iHist == nBins)
      hRatios[iHist] -> Divide(hSumN, hHist[nBins], weight, weight);
  }
  cout << "    Ratios calculated." << endl;


  // set styles
  for (UInt_t iHist = 0; iHist < nHist + 1; iHist++) {
    if (iHist < nBins) {
      // scaled distributions
      hHist[iHist] -> SetLineColor(fCol[iHist]);
      hHist[iHist] -> SetMarkerColor(fCol[iHist]);
      hHist[iHist] -> SetMarkerStyle(fMar[iHist]);
      hHist[iHist] -> SetTitle(sTitle.Data());
      hHist[iHist] -> SetTitleFont(fTxt);
      hHist[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
      hHist[iHist] -> GetXaxis() -> SetTitle(sXtitle.Data());
      hHist[iHist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
      hHist[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
      hHist[iHist] -> GetXaxis() -> CenterTitle(fCnt);
      hHist[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
      hHist[iHist] -> GetYaxis() -> SetTitle(sYtitle.Data());
      hHist[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
      hHist[iHist] -> GetYaxis() -> CenterTitle(fCnt);
      // ratios
      hRatios[iHist] -> SetLineColor(fCol[iHist]);
      hRatios[iHist] -> SetMarkerColor(fCol[iHist]);
      hRatios[iHist] -> SetMarkerStyle(fMar[iHist]);
      hRatios[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
      hRatios[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetXaxis() -> SetTitle(sXtitle.Data());
      hRatios[iHist] -> GetXaxis() -> CenterTitle(fCnt);
      hRatios[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetYaxis() -> SetTitle(sYtitleR.Data());
      hRatios[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    }
    else if (iHist == nBins) {
      // summed distributions
      hSumN -> SetLineColor(fCol[iHist]);
      hSumN -> SetMarkerColor(fCol[iHist]);
      hSumN -> SetMarkerStyle(fMar[iHist]);
      hSumN -> SetTitle(sTitle.Data());
      hSumN -> SetTitleFont(fTxt);
      hSumN -> GetXaxis() -> SetLabelSize(fLbl);
      hSumN -> GetXaxis() -> SetTitle(sXtitle.Data());
      hSumN -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
      hSumN -> GetXaxis() -> SetLabelSize(fLbl);
      hSumN -> GetXaxis() -> CenterTitle(fCnt);
      hSumN -> GetYaxis() -> SetLabelSize(fLbl);
      hSumN -> GetYaxis() -> SetTitle(sYtitle.Data());
      hSumN -> GetYaxis() -> SetLabelSize(fLbl);
      hSumN -> GetYaxis() -> CenterTitle(fCnt);
      // ratios
      hRatios[iHist] -> SetLineColor(fCol[iHist]);
      hRatios[iHist] -> SetMarkerColor(fCol[iHist]);
      hRatios[iHist] -> SetMarkerStyle(fMar[iHist]);
      hRatios[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetXaxis() -> SetTitle(sXtitle.Data());
      hRatios[iHist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
      hRatios[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetXaxis() -> CenterTitle(fCnt);
      hRatios[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetYaxis() -> SetTitle(sYtitleR.Data());
      hRatios[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
      hRatios[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    }
    else {
      hHist[iHist - 1] -> SetLineColor(fCol[iHist]);
      hHist[iHist - 1] -> SetMarkerColor(fCol[iHist]);
      hHist[iHist - 1] -> SetMarkerStyle(fMar[iHist]);
      hHist[iHist - 1] -> SetTitle(sTitle.Data());
      hHist[iHist - 1] -> SetTitleFont(fTxt);
      hHist[iHist - 1] -> GetXaxis() -> SetLabelSize(fLbl);
      hHist[iHist - 1] -> GetXaxis() -> SetTitle(sXtitle.Data());
      hHist[iHist - 1] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
      hHist[iHist - 1] -> GetXaxis() -> SetLabelSize(fLbl);
      hHist[iHist - 1] -> GetXaxis() -> CenterTitle(fCnt);
      hHist[iHist - 1] -> GetYaxis() -> SetLabelSize(fLbl);
      hHist[iHist - 1] -> GetYaxis() -> SetTitle(sYtitle.Data());
      hHist[iHist - 1] -> GetYaxis() -> SetLabelSize(fLbl);
      hHist[iHist - 1] -> GetYaxis() -> CenterTitle(fCnt);
    }
  }  // end bin loop
  cout << "    Styles set." << endl;


  // create labels
  TPaveText *ptScales = new TPaveText(xTxt[0], yTxt[0], xTxt[1], yTxt[1], "NDC NB");
  TPaveText *ptChi2   = new TPaveText(xTxt[0], yTxt[0], xTxt[1], yTxt[1], "NDC NB");
  ptScales -> SetFillColor(fColLbl);
  ptChi2   -> SetFillColor(fColLbl);
  ptScales -> SetFillStyle(fStyLbl);
  ptChi2   -> SetFillStyle(fStyLbl);
  ptScales -> SetLineColor(fColLbl);
  ptChi2   -> SetLineColor(fColLbl);
  ptScales -> SetLineStyle(fStyLbl);
  ptChi2   -> SetLineStyle(fStyLbl);

  // scale factors
  for (UInt_t iBin = 0; iBin < nHist; iBin++) {
    TString sRawScale("");
    TString sTxtScale("#color[");
    if (iBin < nBins)
      sRawScale += weights[iBin];
    else
      sRawScale += normalization;
    sTxtScale += fCol[iBin];
    if (iBin < nBins) {
      sTxtScale += "]{#omega";
      sTxtScale += sBins[iBin].Data();
      sTxtScale += " = ";
    }
    else
      sTxtScale += "]{Norm = ";

    const UInt_t nRawScale = sRawScale.First(".");
    const UInt_t nTxtScale = (nRawScale + nDec) + 1;
    sTxtScale.Append(sRawScale, nTxtScale);
    sTxtScale += "}";
    ptScales -> AddText(sTxtScale.Data());
  }

  // chi2 value
  {
    TString sRawChi2("");
    TString sTxtChi2("#chi^{2}(data, sum) = ");
    sRawChi2 += chi2;

    const UInt_t nRawChi2 = sRawChi2.First(".");
    const UInt_t nTxtChi2 = (nRawChi2 + nDec) + 1;
    sTxtChi2.Append(sRawChi2, nTxtChi2);
    ptChi2 -> AddText(sTxtChi2.Data());
  }
  cout << "    Text created." << endl;


  // create legend
  TLegend *leHist = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
  leHist -> SetFillColor(fColLbl);
  leHist -> SetFillStyle(fStyLbl);
  leHist -> SetLineColor(fColLbl);
  leHist -> SetLineStyle(fStyLbl);

  for (UInt_t iHist = 0; iHist < nHist + 1; iHist++) {
    if (iHist < nBins)
      leHist -> AddEntry(hHist[iHist], sLegs[iHist].Data());
    else if (iHist == nBins)
      leHist -> AddEntry(hSumN, sLegs[iHist].Data());
    else
      leHist -> AddEntry(hHist[iHist -1], sLegs[iHist].Data());
  }
  cout << "    Legend created." << endl;


  // create line
  TLine *liOne = new TLine(pRange[0], 1., pRange[1], 1.);
  liOne -> SetLineColor(fColLne);
  liOne -> SetLineStyle(fStyLne);


  // draw plot
  TCanvas *cHist  = new TCanvas(sCan.Data(), "", wSize, hSize);
  TPad    *pRatio = new TPad(sPadR.Data(), "", xPadR[0], yPadR[0], xPadR[1], yPadR[1]);
  TPad    *pDists = new TPad(sPadD.Data(), "", xPadD[0], yPadD[0], xPadD[1], yPadD[1]);
  cHist  -> cd();
  pRatio -> SetGrid(fGrid, fGrid);
  pRatio -> SetTopMargin(fMargin);
  pDists -> SetGrid(fGrid, fGrid);
  pDists -> SetLogy(fLog);
  pDists -> SetBottomMargin(fMargin);
  cHist  -> cd();
  pRatio -> Draw();
  pDists -> Draw();
  pRatio -> cd();
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    if (iHist == 0)
      hRatios[iHist] -> Draw();
    else
      hRatios[iHist] -> Draw("same");
  }
  liOne  -> Draw();
  ptChi2 -> Draw();
  pDists -> cd();
  for (UInt_t iHist = 0; iHist < nHist + 1; iHist++) {
    if (iHist == 0)
      hHist[iHist] -> Draw();
    else if (iHist < nBins)
      hHist[iHist] -> Draw("same");
    else if (iHist == nBins)
      hSumN -> Draw("same");
    else
      hHist[iHist - 1] -> Draw("same");
  }
  ptScales -> Draw();
  leHist   -> Draw();
  cHist    -> Write();
  cHist    -> Close();
  cout << "    Plots drawn." << endl;


  // set histogram names
  const TString sBin("hBin");
  const TString sData("hData");
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    if (iHist < nBins) {
      TString sBinName(sBin.Data());
      sBinName += iHist;
      hHist[iHist] -> SetName(sBinName.Data());
    }
    else
      hHist[iHist] -> SetName(sData.Data());
  }


  // close files
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] -> Write();
  }
  hSumU   -> Write();
  hSumN   -> Write();
  fOutput -> Close();
  fInput  -> cd();
  fInput  -> Close();
  cout << "  Scale study finished!\n" << endl;

}

// End ------------------------------------------------------------------------
