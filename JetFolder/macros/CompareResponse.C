// 'CompareResponse.C'
// Derek Anderson
// 05.14.2018
//
// Use this to compare two response matrices.

#include <cassert>
#include <iostream>
#include "TH2.h"
#include "TPad.h"
#include "TMath.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TPaveText.h"

using namespace std;


// global constants
static const UInt_t NHist(2);
static const UInt_t NAxes(2);
static const UInt_t NPts(4);



void CompareResponse() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning comparison script..." << endl;


  // i/o parameters
  const TString sOutput("folderCrossCheck_responseWithMacroVsFolder.et911r05qt05130.d19m6y2021.root");
  const TString sInput[NHist] = {"output/June2021/unfoldWithSmoothResponse_usingMacro_forFolderCrossCheck.et911r05qt05130.p0m1k4n38t5.root", "output/June2021/unfoldWithSmoothResponse_usingFolder_forFolderCrossCheck.et911r05qt05130.p0m1k4n38t5.root"};
  const TString sHist[NHist]  = {"hResponse", "hResponseSmooth"};

  // hist parameters
  const TString sRatio("|macro - folder| / macro");
  const TString sRatioName("hRatio");
  const TString sTitle[NHist] = {"smoothed R_{ij} from macro", "smoothed R_{ij} from StJetFolder"};
  const TString sNameH[NHist] = {"hResponseFromMacro", "hResponseFromFolder"};
  const TString sNameP[NHist] = {"pResponseFromMacro", "pResponseFromFolder"};
  const TString sAxis[NAxes]  = {"p_{T}^{det} [GeV/c]", "p_{T}^{par} [GeV/c]"};
  const TString sLabel[NHist] = {"smoothed R_{ij} from macro", "smoothed R_{ij} from StJetFolder"};
  const Float_t xyRange[NPts] = {-1., -1., 57., 57.};
  const Float_t fLabelXY(0.03);
  const Float_t fLabelXYR(0.03);
  const Float_t fLabelZ(0.03);
  const Float_t fLabelZR(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.15);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);

  // profile parameters
  const TString sProfTitle("smoothed R_{ij} profiles");
  const TString sProfBase("pResponse");
  const UInt_t  fAlign(12);
  const UInt_t  cProf[NHist]     = {923, 923};
  const UInt_t  mProf[NHist]     = {20, 21};
  const UInt_t  cProfPlot[NHist] = {923, 899};
  const UInt_t  mProfPlot[NHist] = {20, 25};
  const UInt_t  lProfPlot[NHist] = {1, 1};
  const UInt_t  fProfPlot[NHist] = {0, 0};

  // text parameters
  const TString sSys("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // other parameters
  const UInt_t nRebinX[NHist] = {5, 5};
  const UInt_t nRebinY[NHist] = {5, 5};
  const Bool_t doNorm[NHist]  = {false, false};
  const Bool_t doRebin[NHist] = {false, false};


  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInput[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInput[iHist] = new TFile(sInput[iHist].Data(), "read");
    if (!fInput[iHist]) {
      cerr << "PANIC: couldn't open file " << iHist << "!" << endl;
      assert(0);
    }
  }  // end hist loop
  cout << "    Opened files." << endl;

  // grab histograms
  TH2D *hInput[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInput[iHist] = (TH2D*) fInput[iHist] -> Get(sHist[iHist].Data());
    if (!hInput[iHist]) {
      cerr << "PANIC: couldn't grab histogram " << iHist << "!" << endl;
      assert(0);
    }
    hInput[iHist] -> SetName(sNameH[iHist].Data());
  }  // end hist loop
  cout << "    Grabbed histograms." << endl;


  // do rebin
  cout << "    Rebinning..." << endl;
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (doRebin[iHist]) {
      hInput[iHist] -> Rebin2D(nRebinX[iHist], nRebinY[iHist]);
      cout << "      histogram " << iHist << endl;
    }
  }  // end hist loop


  // do normalization
  cout << "    Normalizing..." << endl;
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (doNorm[iHist]) {
      const UInt_t nBinsX = hInput[iHist] -> GetNbinsX();
      const UInt_t nBinsY = hInput[iHist] -> GetNbinsY();
      for (UInt_t iBinY = 1; iBinY < (nBinsY + 1); iBinY++) {
        const Double_t norm = hInput[iHist] -> Integral(1, nBinsX + 1, iBinY, iBinY);
        if (norm > 0.) {
          for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {
            const Double_t oldVal = hInput[iHist] -> GetBinContent(iBinX, iBinY);
            const Double_t oldErr = hInput[iHist] -> GetBinError(iBinX, iBinY);
            const Double_t newVal = oldVal / norm;
            const Double_t newErr = oldErr / norm;
            hInput[iHist] -> SetBinContent(iBinX, iBinY, newVal);
            hInput[iHist] -> SetBinError(iBinX, iBinY, newErr);
          }  // end x bin loop
        }
      }  // end y bin loop
    }
    cout << "      histogram " << iHist << endl;
  }  // end hist loop


  // calculate ratio
  TH2D *hRatio = (TH2D*) hInput[0] -> Clone();
  hRatio -> SetName(sRatioName.Data());

  const UInt_t  nTotX = hRatio -> GetNbinsX();
  const UInt_t  nTotY = hRatio -> GetNbinsY();
  const UInt_t  nCell = hRatio -> GetBin(nTotX, nTotY);
  for (UInt_t iCell = 0; iCell < nCell; iCell++) {
    const Double_t numVal = hInput[0] -> GetBinContent(iCell);
    const Double_t numErr = hInput[0] -> GetBinError(iCell);
    const Double_t denVal = hInput[1] -> GetBinContent(iCell);
    const Double_t denErr = hInput[1] -> GetBinError(iCell);
    if (numVal > 0) {
      const Double_t ratioV = TMath::Abs(denVal - numVal) / numVal;
      const Double_t ratioE = 0.;
      hRatio -> SetBinContent(iCell, ratioV);
      hRatio -> SetBinError(iCell, ratioE);
    }
    else {
      hRatio -> SetBinContent(iCell, -1.);
      hRatio -> SetBinError(iCell, 0.);
    }
  }  // end cell loop
  const Float_t zMin = hRatio -> GetMinimum(0.);
  const Float_t zMax = hRatio -> GetMaximum();
  hRatio -> GetZaxis() -> SetRangeUser(zMin, zMax);
  cout << "    Calculated ratio...\n"
       << "      max value = " << hRatio -> GetMaximum()
       << endl;


  // find minimum / maximum
  Double_t minimum(0.);
  Double_t maximum(0.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (iHist == 0) {
      minimum = hInput[iHist] -> GetMinimum(0.);
      maximum = hInput[iHist] -> GetMaximum();
    }
    else {
      const Double_t minCheck = hInput[iHist] -> GetMinimum(0.);
      const Double_t maxCheck = hInput[iHist] -> GetMaximum();
      if (minCheck < minimum) minimum = minCheck;
      if (maxCheck > maximum) maximum = maxCheck;
    }
  }
  cout << "    Determined z-axis range:\n"
       << "      min = " << minimum << ", max = " << maximum
       << endl;


  // set styles
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInput[iHist] -> SetTitle(sTitle[iHist].Data());
    hInput[iHist] -> SetTitleFont(fTxt);
    hInput[iHist] -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
    hInput[iHist] -> GetXaxis() -> SetTitle(sAxis[0].Data());
    hInput[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hInput[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hInput[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hInput[iHist] -> GetXaxis() -> SetLabelSize(fLabelXY);
    hInput[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hInput[iHist] -> GetYaxis() -> SetRangeUser(xyRange[1], xyRange[3]);
    hInput[iHist] -> GetYaxis() -> SetTitle(sAxis[1].Data());
    hInput[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hInput[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hInput[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hInput[iHist] -> GetYaxis() -> SetLabelSize(fLabelXY);
    hInput[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hInput[iHist] -> GetZaxis() -> SetRangeUser(minimum, maximum);
    hInput[iHist] -> GetZaxis() -> SetLabelFont(fTxt);
    hInput[iHist] -> GetZaxis() -> SetLabelSize(fLabelZ);
  }
  hRatio -> SetTitle(sRatio.Data());
  hRatio -> SetTitleFont(fTxt);
  hRatio -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
  hRatio -> GetXaxis() -> SetTitle(sAxis[0].Data());
  hRatio -> GetXaxis() -> SetTitleFont(fTxt);
  hRatio -> GetXaxis() -> SetTitleOffset(fOffX);
  hRatio -> GetXaxis() -> SetLabelFont(fTxt);
  hRatio -> GetXaxis() -> SetLabelSize(fLabelXYR);
  hRatio -> GetXaxis() -> CenterTitle(fCnt);
  hRatio -> GetYaxis() -> SetRangeUser(xyRange[1], xyRange[3]);
  hRatio -> GetYaxis() -> SetTitle(sAxis[1].Data());
  hRatio -> GetYaxis() -> SetTitleFont(fTxt);
  hRatio -> GetYaxis() -> SetTitleOffset(fOffY);
  hRatio -> GetYaxis() -> SetLabelFont(fTxt);
  hRatio -> GetYaxis() -> SetLabelSize(fLabelXYR);
  hRatio -> GetYaxis() -> CenterTitle(fCnt);
  hRatio -> GetZaxis() -> SetLabelFont(fTxt);
  hRatio -> GetZaxis() -> SetLabelSize(fLabelZR);
  cout << "    Set styles." << endl;


  // make profiles
  TProfile *pInput[NHist];
  TProfile *pPlot[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    pInput[iHist] = hInput[iHist] -> ProfileX(sNameP[iHist].Data(), 1, -1, "S");
    pInput[iHist] -> SetLineColor(cProf[iHist]);
    pInput[iHist] -> SetMarkerColor(cProf[iHist]);
    pInput[iHist] -> SetMarkerStyle(mProf[iHist]);
    pInput[iHist] -> SetTitle(sTitle[iHist].Data());
    pInput[iHist] -> SetTitleFont(fTxt);
    pInput[iHist] -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
    pInput[iHist] -> GetXaxis() -> SetTitle(sAxis[0].Data());
    pInput[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    pInput[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    pInput[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    pInput[iHist] -> GetXaxis() -> SetLabelSize(fLabelXYR);
    pInput[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    pInput[iHist] -> GetYaxis() -> SetRangeUser(xyRange[2], xyRange[3]);
    pInput[iHist] -> GetYaxis() -> SetTitle(sAxis[1].Data());
    pInput[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    pInput[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    pInput[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    pInput[iHist] -> GetYaxis() -> SetLabelSize(fLabelXYR);
    pInput[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end hist loop
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    // make name
    TString sPlot("pPlot");
    sPlot += iHist;

    // make profile for plot
    pPlot[iHist] = (TProfile*) pInput[iHist] -> Clone();
    pPlot[iHist] -> SetName(sPlot.Data());
    pPlot[iHist] -> SetTitle(sProfTitle.Data());
    pPlot[iHist] -> SetLineColor(cProfPlot[iHist]);
    pPlot[iHist] -> SetLineStyle(lProfPlot[iHist]);
    pPlot[iHist] -> SetFillColor(cProfPlot[iHist]);
    pPlot[iHist] -> SetFillStyle(fProfPlot[iHist]);
    pPlot[iHist] -> SetMarkerColor(cProfPlot[iHist]);
    pPlot[iHist] -> SetMarkerStyle(mProfPlot[iHist]);
    pPlot[iHist] -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[1]);
    pPlot[iHist] -> GetYaxis() -> SetRangeUser(xyRange[2], xyRange[3]);
  }
  cout << "    Made profiles." << endl;


  // make legend and text box
  const UInt_t  cProfFill(0);
  const UInt_t  cProfLine(0);
  const Float_t xyProfLeg[NPts] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyProfTxt[NPts] = {0.3, 0.1, 0.5, 0.3};

  TLegend *lProf = new TLegend(xyProfLeg[0], xyProfLeg[1], xyProfLeg[2], xyProfLeg[3]);
  lProf -> SetFillColor(cProfFill);
  lProf -> SetLineColor(cProfLine);
  lProf -> SetTextFont(fTxt);
  lProf -> SetTextAlign(fAlign);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    lProf -> AddEntry(pPlot[iHist], sLabel[iHist].Data());
  }

  TPaveText *ptProf = new TPaveText(xyProfTxt[0], xyProfTxt[1], xyProfTxt[2], xyProfTxt[3], "NDC NB");
  ptProf -> SetFillColor(cProfFill);
  ptProf -> SetLineColor(cProfLine);
  ptProf -> SetTextFont(fTxt);
  ptProf -> SetTextAlign(fAlign);
  ptProf -> AddText(sSys.Data());
  ptProf -> AddText(sTrg.Data());
  ptProf -> AddText(sJet.Data());
  ptProf -> AddText(sTyp.Data());
  cout << "    Made legend and text boxes." << endl;


  // make response plot
  const UInt_t  widthRes(1500);
  const UInt_t  heightRes(1500);
  const UInt_t  widthProf(750);
  const UInt_t  heightProf(750);
  const UInt_t  log(1);
  const UInt_t  grid(0);
  const UInt_t  ticks(1);
  const Float_t xyNum[NPts] = {0., 0.5, 0.5, 1.};
  const Float_t xyDen[NPts] = {0.5, 0.5, 1., 1.};
  const Float_t xyDiv[NPts] = {0., 0., 0.5, 0.5};
  const Float_t xyPro[NPts] = {0.5, 0., 1., 0.5};

  TCanvas *cResponse = new TCanvas("cResponse", "", widthRes, heightRes);
  TPad    *pNum      = new TPad("pNum", "", xyNum[0], xyNum[1], xyNum[2], xyNum[3]);
  TPad    *pDen      = new TPad("pDen", "", xyDen[0], xyDen[1], xyDen[2], xyDen[3]);
  TPad    *pDiv      = new TPad("pDiv", "", xyDiv[0], xyDiv[1], xyDiv[2], xyDiv[3]);
  TPad    *pPro      = new TPad("pPro", "", xyPro[0], xyPro[1], xyPro[2], xyPro[3]);
  pNum      -> SetGrid(grid, grid);
  pNum      -> SetTicks(ticks, ticks);
  pNum      -> SetLogz(log);
  pDen      -> SetGrid(grid, grid);
  pDen      -> SetTicks(ticks, ticks);
  pDen      -> SetLogz(log);
  pDiv      -> SetGrid(grid, grid);
  pPro      -> SetGrid(grid, grid);
  pPro      -> SetTicks(ticks, ticks);
  pDiv      -> SetTicks(ticks, ticks);
  cResponse -> cd();
  pNum      -> Draw();
  pDen      -> Draw();
  pDiv      -> Draw();
  pPro      -> Draw();
  pNum      -> cd();
  hInput[0] -> Draw("colz");
  pInput[0] -> Draw("same");
  pDen      -> cd();
  hInput[1] -> Draw("colz");
  pInput[1] -> Draw("same");
  pDiv      -> cd();
  hRatio    -> Draw("colz");
  pPro      -> cd();
  pPlot[0]  -> Draw();
  pPlot[1]  -> Draw("sames");
  lProf     -> Draw();
  ptProf    -> Draw();
  fOutput   -> cd();
  cResponse -> Write();
  cResponse -> Close();
  cout << "    Made response plot." << endl;


  // make profile plot
  TCanvas *cProfile = new TCanvas("cProfile", "", widthProf, heightProf);
  cProfile -> SetGrid(grid, grid);
  cProfile -> SetTicks(ticks, ticks);
  pPlot[0] -> Draw();
  pPlot[1] -> Draw("sames");
  lProf    -> Draw();
  ptProf   -> Draw();
  cProfile -> Write();
  cProfile -> Close();
  cout << "    Made profile plot." << endl;


  // save histograms and close files
  fOutput -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hInput[iHist] -> Write();
    pInput[iHist] -> Write();
  }
  hRatio  -> Write();
  fOutput -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInput[iHist] -> cd();
    fInput[iHist] -> Close();
  }
  cout << "  Finished comparison script.\n" << endl;

}

// End ------------------------------------------------------------------------
