// 'PlotScaledDistributions.C'
// Derek Anderson
// 08.25.2017
//
// Use this to scale some distributions, add em' up,
// and compare them to another distribution.
//
// NOTE: if 'embedNorm' is set to -1., then it will
//       scale each element of 'nEmbedTrg' by the
//       corresponding scale and add them up to give
//       the embedding normalization.
//
// NOTE: if 'fScale' is set to 5, it will calculate
//       the weights corresponding to each pTparton
//       bin and weight the distributions accordingly.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

static const UInt_t fScale(5);
static const UInt_t nDist(8);
static const UInt_t nHist(nDist + 1);
static const UInt_t nTotal(11);



void PlotScaledDistributions() {

  //gErrorIgnoreLevel = kError;
  cout << "\n  Beginning plot script..." << endl;


  // i/o parameters
  const TString sOutput("pythiaXembedding.pTcorr.geant.r03a02rm1chrg.d27m11y2017.root");
  const TString sInput[nHist] = {"pp200r12pt5.test.r03a02rm1chrg.plots.root", "pp200r12pt7.test.r03a02rm1chrg.plots.root", "pp200r12pt9.test.r03a02rm1chrg.plots.root", "pp200r12pt11.test.r03a02rm1chrg.plots.root", "pp200r12pt15.test.r03a02rm1chrg.plots.root", "pp200r12pt20.test.r03a02rm1chrg.plots.root", "pp200r12pt25.test.r03a02rm1chrg.plots.root", "pp200r12pt35.test.r03a02rm1chrg.plots.root", "pythia.pi0eTtrg9.r03a02rm1chrg.plots.d27m11y2017.root"};

  // histogram parameters
  const TString sHist[nHist] = {"Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP", "Pi0/hJetPtCorrP"};
  const TString sTitle("Recoil jet p_{T}^{corr}: |#Delta#varphi - #pi| < #pi/4");
  const TString sTitleX("p_{T}^{corr} = p_{T}^{jet} - #rhoA^{jet}");
  const TString sTitleYR("embedding / data");
  const TString sTitleYD("dN^{jet}/dp_{T}^{corr}");
  const TString sTitleYDN("(1/N^{trg}_{eff}) dN^{jet}/dp_{T}^{corr}");

  // scales and misc. parameters
  const Double_t fudge[nTotal]    = {1.228, 1.051, 1.014, 1., 1., 1., 1., 1., 1., 1., 1.};
  const Double_t xSctn[nTotal]    = {9.005805, 1.461907, 0.3544354, 0.1513760, 2.488644e-02, 5.845846e-03, 2.304880e-03, 3.426618e-04, 4.562988e-05, 9.738044e-06, 5.019977e-07};
  const Double_t nEvts[nTotal]    = {2100295., 600300., 600300., 300289., 300289., 300289., 160295., 100302., 80293., 76303., 23307.};
  const Double_t nTrgEmbed[nDist] = {2., 34., 245., 659., 1982., 4483., 9044., 5600.};
  const Double_t dataNorm(250000.);
  const Double_t embedNorm(-1.);
  const Bool_t   doNorm(false);
  const Bool_t   doIntNorm(false);
  const Bool_t   doBinNorm(false);
  const Bool_t   doRebin(false);
  const UInt_t   nRebin(3);
  const UInt_t   nRebinS(3);

  // constants
  const UInt_t  fCol[3] = {810, 890, 850};
  const UInt_t  fMar[3] = {7, 4, 4};
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fColLeg(0);
  const UInt_t  fStyLeg(0);
  const UInt_t  fColLin(1);
  const UInt_t  fStyLin(2);
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const UInt_t  nDec(2);
  const Float_t weight(1.);
  const Float_t fLblSizeD(0.02);
  const Float_t fLblSizeR(0.035);
  const Float_t fTtlSizeR(0.06);
  const Float_t fLblOffR(0.0035);
  const Float_t fTtlOffR(0.6);
  const Float_t wMar(0.);
  const Float_t pRange[2] = {-3., 43.};
  const Float_t xLeg[2]   = {0.1, 0.3};
  const Float_t yLeg[2]   = {0.1, 0.3};
  const Float_t xPav[2]   = {0.3, 0.5};
  const Float_t yPav[2]   = {0.1, 0.3};
  const Float_t xPavN[2]  = {0.5, 0.7};
  const Float_t yPavN[2]  = {0.1, 0.3};
  const Float_t xPad[2]   = {0., 1.};
  const Float_t yPadR[2]  = {0., 0.35};
  const Float_t yPadD[2]  = {0.35, 1.};
  const Float_t width(700);
  const Float_t height(950);
  const TString sSum("hSum");
  const TString sSumN("hSumNorm");
  const TString sCan("cHist");
  const TString sCanN("cNorm");
  const TString sPadR("pRatio");
  const TString sPadRN("pRatioNorm");
  const TString sPadD("pDist");
  const TString sPadDN("pDistNorm");
  const TString sRatio[nHist]  = {"hRatio0", "hRatio1", "hRatio2", "hRatio3", "hRatio4", "hRatio5", "hRatio6", "hRatio7", "hRatioSum"};
  const TString sRatioN[nHist] = {"hRatio0N", "hRatio1N", "hRatio2N", "hRatio3N", "hRatio4N", "hRatio5N", "hRatio6N", "hRatio7N", "hRatioSumN"};


  // open files and grab histograms
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    fInput[iHist]  = new TFile(sInput[iHist].Data(), "read");
    if (!fInput[iHist]) {
      cerr << "PANIC: couldn't open input file " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Files opened." << endl;

  TH1D *hHist[nHist];
  TH1D *hNorm[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hHist[iHist] = (TH1D*) fInput[iHist] -> Get(sHist[iHist].Data()) -> Clone();
    hNorm[iHist] = (TH1D*) fInput[iHist] -> Get(sHist[iHist].Data()) -> Clone();
    if (!hHist[iHist]) {
      cerr << "PANIC: couldn't grab histogram " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Histograms grabbed." << endl;


  // calculate scale set 5
  Double_t scales5[nDist];
  UInt_t   iStart = nTotal - nDist;
  for (UInt_t iDist = iStart; iDist < nTotal; iDist++) {
    const Double_t binLumi   = nEvts[iDist] / xSctn[iDist];
    const Double_t binWeight = 1. / (fudge[iDist] * binLumi);
    scales5[iDist - iStart]  = binWeight;
  }

  // select scales
  Double_t *scales;
  switch (fScale) {
    case 1:
      scales = scales1;
      break;
    case 2:
      scales = scales2;
      break;
    case 3:
      scales = scales3;
      break;
    case 4:
      scales = scales4;
      break;
    case 5:
      scales = scales5;
      break;
  }

  // norm check
  Double_t normer = 1.;
  if (doNorm) {
    TFile *fNormer = new TFile("scaleFactors.eTtrg.thirdJetMaker.allFiles.d2m11y2017.root", "read");
    TH1D  *hEmbedN = (TH1D*) fNormer -> Get("hSumUnormalized");
    TH1D  *hDataN  = (TH1D*) fNormer -> Get("hData");

    const Double_t embInt = hEmbedN -> Integral();
    const Double_t datInt = hDataN  -> Integral();
    normer = datInt / embInt;
    fNormer -> Close();
    fOutput -> cd();
  }

  // scale histograms
  Double_t nTrgTotal(0.);
  Double_t nTrgScale[nDist];
  for (UInt_t iHist = 0; iHist < nDist; iHist++) {
    hHist[iHist] -> Scale(scales[iHist]);
    hNorm[iHist] -> Scale(scales[iHist]);
    hHist[iHist] -> Scale(normer);
    hNorm[iHist] -> Scale(normer);
    nTrgScale[iHist]  = nTrgEmbed[iHist] * scales[iHist] * normer;
    nTrgTotal        += nTrgScale[iHist];
  }
  if (embedNorm != -1.) nTrgTotal = embedNorm;

  // normalize histograms
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    if (iHist == nDist)
      hNorm[iHist] -> Scale(1. / dataNorm);
    else
      hNorm[iHist] -> Scale(1. / nTrgTotal);
  }
  cout << "    Normalization:\n"
       << "      data  -- " << dataNorm <<"\n"
       << "      embed -- " << nTrgTotal
       << endl;

  fOutput -> cd();
  cout << "    Histograms scaled." << endl;


  // sum histograms
  const UInt_t  nBins = hHist[0] -> GetNbinsX();
  const Float_t xBin1 = hHist[0] -> GetBinLowEdge(1);
  const Float_t xBin2 = hHist[0] -> GetBinLowEdge(nBins + 1);

  TH1D *hSum  = new TH1D(sSum.Data(), "", nBins, xBin1, xBin2);
  TH1D *hSumN = new TH1D(sSumN.Data(), "", nBins, xBin1, xBin2);
  hSum  -> Sumw2();
  hSumN -> Sumw2();
  for (UInt_t iDist = 0; iDist < nDist; iDist++) {
    hSum  -> Add(hHist[iDist]);
    hSumN -> Add(hHist[iDist]);
  }
  hSumN -> Scale(1. / nTrgTotal);
  cout << "    Histograms summed." << endl;

  // rebin histograms (if necessary)
  if (doRebin) {
    for (UInt_t iHist = 0; iHist < nHist; iHist++) {
      hHist[iHist] -> Rebin(nRebin);
      hNorm[iHist] -> Rebin(nRebin);
    }
    hSum  -> Rebin(nRebinS);
    hSumN -> Rebin(nRebinS);
    cout << "    Histograms rebinned." << endl;
  }

  // quick fix [10.09.2017]
  const Double_t intDataN = hNorm[nDist] -> Integral();
  const Double_t intEmbdN = hSumN        -> Integral();
  if (doIntNorm) {
    hNorm[nDist] -> Scale(1. / intDataN);
    hSumN        -> Scale(1. / intEmbdN);
  }

  // quick fix [11.07.2017]
  if (doBinNorm) {
    const Double_t binData  = hHist[nDist] -> GetBinWidth(17);
    const Double_t binDataN = hNorm[nDist] -> GetBinWidth(17);
    const Double_t binEmbd  = hSum         -> GetBinWidth(17);
    const Double_t binEmbdN = hSumN        -> GetBinWidth(17);
    hHist[nDist] -> Scale(1. / binData);
    hNorm[nDist] -> Scale(1. / binDataN);
    hSum         -> Scale(1. / binEmbd);
    hSumN        -> Scale(1. / binEmbdN);
  }


  // calculate ratios
  const UInt_t  nBinsRatio = hHist[nDist] -> GetNbinsX();
  const Float_t xBin1Ratio = hHist[nDist] -> GetBinLowEdge(1);
  const Float_t xBin2Ratio = hHist[nDist] -> GetBinLowEdge(nBinsRatio + 1);

  TH1D *hRatio[nHist];
  TH1D *hRatioN[nHist];
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    hRatio[iHist]  = new TH1D(sRatio[iHist].Data(), "", nBinsRatio, xBin1Ratio, xBin2Ratio);
    hRatioN[iHist] = new TH1D(sRatioN[iHist].Data(), "", nBinsRatio, xBin1Ratio,xBin2Ratio);
    hRatio[iHist]  -> Sumw2();
    hRatioN[iHist] -> Sumw2();
    if (iHist < nDist) {
      hRatio[iHist]  -> Divide(hHist[iHist], hHist[nDist], weight, weight);
      hRatioN[iHist] -> Divide(hNorm[iHist], hNorm[nDist], weight, weight); 
    }
    else {
      hRatio[iHist]  -> Divide(hSum, hHist[nDist], weight, weight);
      hRatioN[iHist] -> Divide(hSumN, hNorm[nDist], weight, weight);
    }
  }
  cout << "    Ratios calculated." << endl;

  // set styles
  hHist[nDist] -> SetLineColor(fCol[0]);
  hHist[nDist] -> SetMarkerColor(fCol[0]);
  hHist[nDist] -> SetMarkerStyle(fMar[0]);
  hHist[nDist] -> SetTitle(sTitle.Data());
  hHist[nDist] -> SetTitleFont(fTxt);
  hHist[nDist] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHist[nDist] -> GetXaxis() -> SetTitleFont(fTxt);
  hHist[nDist] -> GetXaxis() -> CenterTitle(fCnt);
  hHist[nDist] -> GetXaxis() -> SetLabelSize(fLblSizeD);
  hHist[nDist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hHist[nDist] -> GetYaxis() -> SetTitle(sTitleYD.Data());
  hHist[nDist] -> GetYaxis() -> SetTitleFont(fTxt);
  hHist[nDist] -> GetYaxis() -> CenterTitle(fCnt);
  hHist[nDist] -> GetYaxis() -> SetLabelSize(fLblSizeD);
  // normalized distributions
  hNorm[nDist] -> SetLineColor(fCol[0]);
  hNorm[nDist] -> SetMarkerColor(fCol[0]);
  hNorm[nDist] -> SetMarkerStyle(fMar[0]);
  hNorm[nDist] -> SetTitle(sTitle.Data());
  hNorm[nDist] -> SetTitleFont(fTxt);
  hNorm[nDist] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hNorm[nDist] -> GetXaxis() -> SetTitleFont(fTxt);
  hNorm[nDist] -> GetXaxis() -> CenterTitle(fCnt);
  hNorm[nDist] -> GetXaxis() -> SetLabelSize(fLblSizeD);
  hNorm[nDist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hNorm[nDist] -> GetYaxis() -> SetTitle(sTitleYDN.Data());
  hNorm[nDist] -> GetYaxis() -> SetTitleFont(fTxt);
  hNorm[nDist] -> GetYaxis() -> CenterTitle(fCnt);
  hNorm[nDist] -> GetYaxis() -> SetLabelSize(fLblSizeD);
  // summed distribution
  hSum -> SetLineColor(fCol[1]);
  hSum -> SetMarkerColor(fCol[1]);
  hSum -> SetMarkerStyle(fMar[1]);
  hSum -> SetTitle(sTitle.Data());
  hSum -> SetTitleFont(fTxt);
  hSum -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSum -> GetXaxis() -> SetTitleFont(fTxt);
  hSum -> GetXaxis() -> CenterTitle(fCnt);
  hSum -> GetXaxis() -> SetLabelSize(fLblSizeD);
  hSum -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hSum -> GetYaxis() -> SetTitle(sTitleYD.Data());
  hSum -> GetYaxis() -> SetTitleFont(fTxt);
  hSum -> GetYaxis() -> CenterTitle(fCnt);
  hSum -> GetYaxis() -> SetLabelSize(fLblSizeD);
  // summed (normalized) distribution
  hSumN -> SetLineColor(fCol[1]);
  hSumN -> SetMarkerColor(fCol[1]);
  hSumN -> SetMarkerStyle(fMar[1]);
  hSumN -> SetTitle(sTitle.Data());
  hSumN -> SetTitleFont(fTxt);
  hSumN -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSumN -> GetXaxis() -> SetTitleFont(fTxt);
  hSumN -> GetXaxis() -> CenterTitle(fCnt);
  hSumN -> GetXaxis() -> SetLabelSize(fLblSizeD);
  hSumN -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hSumN -> GetYaxis() -> SetTitle(sTitleYDN.Data());
  hSumN -> GetYaxis() -> SetTitleFont(fTxt);
  hSumN -> GetYaxis() -> CenterTitle(fCnt);
  hSumN -> GetYaxis() -> SetLabelSize(fLblSizeD);
  // ratios
  hRatio[nDist] -> SetLineColor(fCol[2]);
  hRatio[nDist] -> SetMarkerColor(fCol[2]);
  hRatio[nDist] -> SetMarkerStyle(fMar[2]);
  hRatio[nDist] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatio[nDist] -> GetXaxis() -> SetTitleFont(fTxt);
  hRatio[nDist] -> GetXaxis() -> SetTitleSize(fTtlSizeR);
  hRatio[nDist] -> GetXaxis() -> SetTitleOffset(fTtlOffR);
  hRatio[nDist] -> GetXaxis() -> CenterTitle(fCnt);
  hRatio[nDist] -> GetXaxis() -> SetLabelSize(fLblSizeR);
  hRatio[nDist] -> GetXaxis() -> SetLabelOffset(fLblOffR);
  hRatio[nDist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hRatio[nDist] -> GetYaxis() -> SetTitle(sTitleYR.Data());
  hRatio[nDist] -> GetYaxis() -> SetTitleFont(fTxt);
  hRatio[nDist] -> GetYaxis() -> SetTitleSize(fTtlSizeR);
  hRatio[nDist] -> GetYaxis() -> SetTitleOffset(fTtlOffR);
  hRatio[nDist] -> GetYaxis() -> CenterTitle(fCnt);
  hRatio[nDist] -> GetYaxis() -> SetLabelSize(fLblSizeR);
  hRatio[nDist] -> GetYaxis() -> SetLabelOffset(fLblOffR);
  // normalized ratios
  hRatioN[nDist] -> SetLineColor(fCol[2]);
  hRatioN[nDist] -> SetMarkerColor(fCol[2]);
  hRatioN[nDist] -> SetMarkerStyle(fMar[2]);
  hRatioN[nDist] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioN[nDist] -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioN[nDist] -> GetXaxis() -> SetTitleSize(fTtlSizeR);
  hRatioN[nDist] -> GetXaxis() -> SetTitleOffset(fTtlOffR);
  hRatioN[nDist] -> GetXaxis() -> CenterTitle(fCnt);
  hRatioN[nDist] -> GetXaxis() -> SetLabelSize(fLblSizeR);
  hRatioN[nDist] -> GetXaxis() -> SetLabelOffset(fLblOffR);
  hRatioN[nDist] -> GetXaxis() -> SetRangeUser(pRange[0], pRange[1]);
  hRatioN[nDist] -> GetYaxis() -> SetTitle(sTitleYR.Data());
  hRatioN[nDist] -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioN[nDist] -> GetYaxis() -> SetTitleSize(fTtlSizeR);
  hRatioN[nDist] -> GetYaxis() -> SetTitleOffset(fTtlOffR);
  hRatioN[nDist] -> GetYaxis() -> CenterTitle(fCnt);
  hRatioN[nDist] -> GetYaxis() -> SetLabelSize(fLblSizeR);
  hRatioN[nDist] -> GetYaxis() -> SetLabelOffset(fLblOffR);
  cout << "    Styles set." << endl;


  // make legend
  TLegend *lHist = new TLegend(xLeg[0], yLeg[0], xLeg[1], yLeg[1]);
  lHist -> SetFillColor(fColLeg);
  lHist -> SetFillStyle(fStyLeg);
  lHist -> SetLineColor(fColLeg);
  lHist -> SetLineStyle(fStyLeg);
  lHist -> SetTextFont(fTxt);
  lHist -> AddEntry(hHist[nDist], "Pythia 8 (standalone)");
  lHist -> AddEntry(hSum, "Pythia 6 (embedding)");
  cout << "    Legend made." << endl;


  // make labels
  cout << "    Labels made." << endl;


  // make line
  TLine *lOne = new TLine(pRange[0], 1., pRange[1], 1.);
  lOne -> SetLineColor(fColLin);
  lOne -> SetLineStyle(fStyLin);


  // make plots
  TCanvas *cHist  = new TCanvas(sCan.Data(), "", width, height);
  TPad    *pRatio = new TPad(sPadR.Data(), "", xPad[0], yPadR[0], xPad[1], yPadR[1]);
  TPad    *pHist  = new TPad(sPadD.Data(), "", xPad[0], yPadD[0], xPad[1], yPadD[1]);
  pRatio        -> SetLogy(fLog);
  pRatio        -> SetGrid(fGrid, fGrid);
  pRatio        -> SetTopMargin(wMar);
  pHist         -> SetLogy(fLog);
  pHist         -> SetGrid(fGrid, fGrid);
  pHist         -> SetBottomMargin(wMar);
  pRatio        -> Draw();
  pHist         -> Draw();
  pRatio        -> cd();
  hRatio[nDist] -> Draw();
  lOne          -> Draw();
  pHist         -> cd();
  hHist[nDist]  -> Draw();
  hSum          -> Draw("same");
  lHist         -> Draw();
  cHist         -> Write();
  cHist         -> Close();

  TCanvas *cNorm   = new TCanvas(sCanN.Data(), "", width, height);
  TPad    *pRatioN = new TPad(sPadRN.Data(), "", xPad[0], yPadR[0], xPad[1], yPadR[1]);
  TPad    *pNorm   = new TPad(sPadDN.Data(), "", xPad[0], yPadD[0], xPad[1], yPadD[1]);
  pRatioN        -> SetLogy(fLog);
  pRatioN        -> SetGrid(fGrid, fGrid);
  pRatioN        -> SetTopMargin(wMar);
  pNorm          -> SetLogy(fLog);
  pNorm          -> SetGrid(fGrid, fGrid);
  pNorm          -> SetBottomMargin(wMar);
  pRatioN        -> Draw();
  pNorm          -> Draw();
  pRatioN        -> cd();
  hRatioN[nDist] -> Draw();
  lOne           -> Draw();
  pNorm          -> cd();
  hNorm[nDist]   -> Draw();
  hSumN          -> Draw("same");
  lHist          -> Draw();
  cNorm          -> Write();
  cNorm          -> Close();


  // set names of data and sum
  hSum         -> SetName("hSum");
  hSumN        -> SetName("hSumPerTrigger");
  hHist[nDist] -> SetName("hData");
  hNorm[nDist] -> SetName("hDataPerTrigger");


  // close files
  fOutput      -> cd();
  hSum         -> Write();
  hSumN        -> Write();
  hHist[nDist] -> Write();
  hNorm[nDist] -> Write();
  fOutput      -> Close();
  for (UInt_t iHist = 0; iHist < nHist; iHist++) {
    fInput[iHist] -> cd();
    fInput[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
