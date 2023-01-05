// 'QuickSmear.C'
// Derek Anderson
// 07.08.2021
//
// Quick macro that grabs a prior,
// efficiency, and response matrix
// and generates a corresponding
// smeared prior and re-trained
// response matrix.


#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

// global constants
static const UInt_t NPts(2);
static const UInt_t NPar(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NHist(2);
static const UInt_t NAxes(2);
static const UInt_t NMc(50000000);


void QuickSmear() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning quick smearing macro..." << endl;

  // io parameters
  const TString sOut("smearedPriorWithSmoothResRFF_forClosureTest_pTbinHuge.et911gam.r05a065rm1chrg.d25m10y2021.root");
  const TString sParIn("input/pythia/pp200py8par.forComparison_pTbinHuge.et911nTrg100Kpt021Kgam.r05a065rm1chrg.root");
  const TString sEffIn("et911r05ff_rebinClosure/pp200r9ff.default_modStatsUnfoldWithRFF_pTbinHuge.et911r05qt05130.p0m1k3n38t5.root");
  const TString sResIn("et911r05ff_rebinClosure/pp200r9ff.default_modStatsUnfoldWithRFF_pTbinHuge.et911r05qt05130.p0m1k3n38t5.root");
  const TString sPar("Gam/hJetPtCorrG");
  const TString sEff("hEfficiency");
  const TString sRes("hResponseSmooth");

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg("#gamma_{dir} trigger, E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // efficiency fit parameters
  const TString sEffName("fEfficiency");
  const TString sEffTotName("fEffTotalFit");
  const TString sEffFunc("[0] * (1 - TMath::Exp(-1.*[1]*x))");
  const TString sEffParNames[NPar] = {"#epsilon_{0}", "#sigma"};
  const Float_t xEffRange[NPts]    = {0., 100.};
  const Float_t xFitRange[NPts]    = {0.2, 15.};
  const Float_t pEffGuess[NPar]    = {0.87, 4.};
  const UInt_t  fColEff(2);

  // distribution parameters
  const TString sTitle("");
  const TString sTitleY("a. u.");
  const TString sTitleRY("out / in");
  const TString sTitleP("p_{T}^{par}");
  const TString sTitleD("p_{T}^{det}");
  const TString sLegPar("input particle-level distribution");
  const TString sLegDet("output detector-level distribution");
  const Float_t xRange[NPts] = {-1., 57.};
  const Float_t fLbl(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.15);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fCol[NHist] = {923, 899};
  const UInt_t  fMar[NHist] = {20, 25};
  const UInt_t  fLin[NHist] = {1, 1};
  const UInt_t  fFil[NHist] = {0, 0};

  // martrix parameters
  const TString sRatio("|input - output| / input");
  const TString sProfTitle("R_{ij} profiles");
  const TString sTitleM[NHist] = {"input R_{ij}", "output R_{ij}"};
  const TString sAxis[NAxes]   = {"p_{T}^{det} [GeV/c]", "p_{T}^{par} [GeV/c]"};
  const TString sLegM[NHist]   = {"input R_{ij}", "output R_{ij}"};
  const Float_t xyRange[NVtx]  = {-1., -1., 57., 57.};
  const Float_t fLabelXY(0.03);
  const Float_t fLabelXYR(0.03);
  const Float_t fLabelZ(0.03);
  const Float_t fLabelZR(0.03);
  const UInt_t  cProf[NHist] = {923, 923};
  const UInt_t  cPlot[NHist] = {923, 899};
  const UInt_t  mProf[NHist] = {20, 21};
  const UInt_t  mPlot[NHist] = {20, 25};
  const UInt_t  fLinP[NHist] = {1, 1};
  const UInt_t  fFilP[NHist] = {0, 0};


  // open files
  TFile *fOut   = new TFile(sOut.Data(), "recreate");
  TFile *fParIn = new TFile(sParIn.Data(), "read");
  TFile *fEffIn = new TFile(sEffIn.Data(), "read");
  TFile *fResIn = new TFile(sResIn.Data(), "read");
  if (!fOut || !fParIn || !fEffIn || !fResIn) {
    cerr << "PANIC: couldn't open a file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hParIn = (TH1D*) fParIn -> Get(sPar.Data());
  TH1D *hEffIn = (TH1D*) fEffIn -> Get(sEff.Data());
  TH2D *hResIn = (TH2D*) fResIn -> Get(sRes.Data());
  if (!hParIn || !hEffIn || !hResIn) {
    cerr << "PANIC: couldn't grab a histogram!" << endl;
    return;
  }
  hParIn -> SetNameTitle("hParticleInput", "");
  hEffIn -> SetNameTitle("hEfficiencyInput", "");
  hResIn -> SetNameTitle("hResponseInput", "");
  cout << "    Grabbed histograms." << endl;


  // get efficiency function
  TF1 *fEff = new TF1(sEffName.Data(), sEffFunc.Data(), xEffRange[0], xEffRange[1]);
  fEff   -> SetParameter(0, pEffGuess[0]);
  fEff   -> SetParameter(1, pEffGuess[1]);
  fEff   -> SetParName(0, sEffParNames[0].Data());
  fEff   -> SetParName(1, sEffParNames[1].Data());
  fEff   -> SetLineColor(fColEff);
  hEffIn -> Fit(sEffName.Data(), "", "", xFitRange[0], xFitRange[1]);
  cout << "    Fit efficiency." << endl;


  // create output histograms
  hNorm   = (TH1D*) hParIn -> Clone();
  hDetOut = (TH1D*) hParIn -> Clone();
  hResOut = (TH2D*) hResIn -> Clone();
  hNorm   -> SetNameTitle("hNormalization", "");
  hDetOut -> SetNameTitle("hDetectorOutput", "");
  hResOut -> SetNameTitle("hResponseOutput", "");
  hNorm   -> Reset("ICES");
  hDetOut -> Reset("ICES");
  hResOut -> Reset("ICES");
  cout << "    Created output histograms." << endl;

  // create histograms for sampling
  TH1D *hProj;
  TH1D *hSample = (TH1D*) hParIn -> Clone();
  hSample -> SetName("hForSampling");
  hSample -> Multiply(fEff);
  cout << "    Created sampling histograms.\n"
       << "    Beginning MC loop..."
       << endl;


  // MC loop
  UInt_t   iPar(0);
  UInt_t   nProj(0);
  Double_t par(0.);
  Double_t det(0.);
  for (Int_t iMC = 0; iMC < NMc; iMC++) {

    // announce progress
    if ((iMC + 1) < NMc) {
      cout << "      Processing iteration " << (iMC + 1) << "/" << NMc << "\r" << flush;
    } else {
      cout << endl;
    }

    // apply smearing
    nProj = 0;
    do {
      par   = hSample -> GetRandom();
      iPar  = hResIn  -> GetYaxis() -> FindBin(par);
      hProj = hResIn  -> ProjectionX("", iPar, iPar);
      nProj = hProj   -> GetEntries();
      det   = hProj   -> GetRandom();
    }  while (nProj < 1);
    hNorm   -> Fill(par);
    hDetOut -> Fill(det);
    hResOut -> Fill(det, par);

  }  // end MC loop
  cout << "    Finished MC loop." << endl;


  // set negative bins to 0
  const UInt_t nNorm = hNorm   -> GetNbinsX();
  const UInt_t nDet  = hDetOut -> GetNbinsX();
  for (UInt_t iNorm = 1; iNorm < nNorm + 1; iNorm++) {
    const Double_t nVal = hNorm -> GetBinContent(iNorm);
    if (nVal < 0.) {
      hNorm -> SetBinContent(iNorm, 0.);
      hNorm -> SetBinError(iNorm, 0.);
    }
  }  // end norm bin loop
  for (UInt_t iDet = 1; iDet < nDet + 1; iDet++) {
    const Double_t bVal = hDetOut -> GetBinContent(iDet);
    if (bVal < 0.) {
      hDetOut -> SetBinContent(iDet, 0.);
      hDetOut -> SetBinError(iDet, 0.);
    }
  }  // end detector bin loop
  cout << "    Removed any negative bins." << endl;


  // normalize backfolded spectrum
  const Double_t intNorm   = hNorm   -> Integral();
  const Double_t intSample = hSample -> Integral();
  if (intSample > 0.) {
    const Double_t scale = intSample / intNorm;
    hDetOut -> Scale(scale);
  }

  // normalize response
  const UInt_t nResX = hResOut -> GetNbinsX();
  const UInt_t nResY = hResOut -> GetNbinsY();
  for (UInt_t iResY = 0; iResY < nResY; iResY++) {
    const Double_t binNorm = hResOut -> Integral(1, nResX, iResY, iResY);
    if (binNorm != 0.) {
      for (UInt_t iResX = 0; iResX < nResX; iResX++) {
        const Double_t binVal = hResOut -> GetBinContent(iResX, iResY);
        const Double_t binErr = hResOut -> GetBinError(iResX, iResY);
        const Double_t newVal = (binVal / binNorm);
        const Double_t newErr = (binErr / binNorm);
        hResOut -> SetBinContent(iResX, iResY, newVal);
        hResOut -> SetBinError(iResX, iResY, newErr);
      }  // end x loop
    }
  }  // end y loop
  cout << "    Normalized output." << endl;


  // create X-axis title for distribution plot
  TString sTitleX("#color[");
  sTitleX += fCol[0];
  sTitleX.Append("]{");
  sTitleX.Append(sTitleP.Data());
  sTitleX.Append("}, #color[");
  sTitleX += fCol[1];
  sTitleX.Append("]{");
  sTitleX.Append(sTitleD.Data());
  sTitleX.Append("} [GeV/c]");

  // set histogram styles
  hParIn  -> SetFillColor(fCol[0]);
  hParIn  -> SetFillStyle(fFil[0]);
  hParIn  -> SetLineColor(fCol[0]);
  hParIn  -> SetLineStyle(fLin[0]);
  hParIn  -> SetMarkerColor(fCol[0]);
  hParIn  -> SetMarkerStyle(fMar[0]);
  hParIn  -> SetTitle(sTitle.Data());
  hParIn  -> SetTitleFont(fTxt);
  hParIn  -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hParIn  -> GetXaxis() -> SetLabelSize(fLbl);
  hParIn  -> GetXaxis() -> SetLabelFont(fTxt);
  hParIn  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hParIn  -> GetXaxis() -> SetTitleFont(fTxt);
  hParIn  -> GetXaxis() -> SetTitleOffset(fOffX);
  hParIn  -> GetXaxis() -> CenterTitle(fCnt);
  hParIn  -> GetYaxis() -> SetLabelSize(fLbl);
  hParIn  -> GetYaxis() -> SetLabelFont(fTxt);
  hParIn  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParIn  -> GetYaxis() -> SetTitleFont(fTxt);
  hParIn  -> GetYaxis() -> SetTitleOffset(fOffY);
  hParIn  -> GetYaxis() -> CenterTitle(fCnt);
  hDetOut -> SetFillColor(fCol[1]);
  hDetOut -> SetFillStyle(fFil[1]);
  hDetOut -> SetLineColor(fCol[1]);
  hDetOut -> SetLineStyle(fLin[1]);
  hDetOut -> SetMarkerColor(fCol[1]);
  hDetOut -> SetMarkerStyle(fMar[1]);
  hDetOut -> SetTitle(sTitle.Data());
  hDetOut -> SetTitleFont(fTxt);
  hDetOut -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hDetOut -> GetXaxis() -> SetLabelSize(fLbl);
  hDetOut -> GetXaxis() -> SetLabelFont(fTxt);
  hDetOut -> GetXaxis() -> SetTitle(sTitleX.Data());
  hDetOut -> GetXaxis() -> SetTitleFont(fTxt);
  hDetOut -> GetXaxis() -> SetTitleOffset(fOffX);
  hDetOut -> GetXaxis() -> CenterTitle(fCnt);
  hDetOut -> GetYaxis() -> SetLabelSize(fLbl);
  hDetOut -> GetYaxis() -> SetLabelFont(fTxt);
  hDetOut -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDetOut -> GetYaxis() -> SetTitleFont(fTxt);
  hDetOut -> GetYaxis() -> SetTitleOffset(fOffY);
  hDetOut -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set distribution styles." << endl;


  // calculate ratios
  TH1D *hRatio  = (TH1D*) hParIn -> Clone();
  TH2D *hRatioM = (TH2D*) hResIn -> Clone();
  hRatio  -> SetName("hRatioParVsDet");
  hRatioM -> SetName("hMatrixRatio");
  hRatio  -> Reset("ICES");
  hRatioM -> Reset("ICES");
  hRatio  -> Divide(hDetOut, hParIn, 1., 1.);

  // take ratio of matrices
  const UInt_t nTotX = hRatioM -> GetNbinsX();
  const UInt_t nTotY = hRatioM -> GetNbinsY();
  const UInt_t nCell = hRatioM -> GetBin(nTotX, nTotY);
  for (UInt_t iCell = 0; iCell < nCell; iCell++) {
    const Double_t numVal = hResIn  -> GetBinContent(iCell);
    const Double_t numErr = hResIn  -> GetBinError(iCell);
    const Double_t denVal = hResOut -> GetBinContent(iCell);
    const Double_t denErr = hResOut -> GetBinError(iCell);
    if (numVal > 0) {
      const Double_t ratioV = TMath::Abs(denVal - numVal) / numVal;
      const Double_t ratioE = 0.;
      hRatioM -> SetBinContent(iCell, ratioV);
      hRatioM -> SetBinError(iCell, ratioE);
    } else {
      hRatioM -> SetBinContent(iCell, -1.);
      hRatioM -> SetBinError(iCell, 0.);
    }
  }  // end cell loop
  const Float_t zMin = hRatioM -> GetMinimum(0.);
  const Float_t zMax = hRatioM -> GetMaximum();
  hRatioM -> GetZaxis() -> SetRangeUser(zMin, zMax);
  cout << "    Calculate ratios." << endl;


  // find minimum / maximum
  Double_t minimum(0.);
  Double_t maximum(0.);
  minimum = hResIn -> GetMinimum(0.);
  maximum = hResIn -> GetMaximum();

  const Double_t minCheck = hResOut -> GetMinimum(0.);
  const Double_t maxCheck = hResOut -> GetMaximum();
  if (minCheck < minimum) minimum = minCheck;
  if (maxCheck > maximum) maximum = maxCheck;

  // set matrix styles
  hResIn  -> SetTitle(sTitleM[0].Data());
  hResIn  -> SetTitleFont(fTxt);
  hResIn  -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
  hResIn  -> GetXaxis() -> SetTitle(sAxis[0].Data());
  hResIn  -> GetXaxis() -> SetTitleFont(fTxt);
  hResIn  -> GetXaxis() -> SetTitleOffset(fOffX);
  hResIn  -> GetXaxis() -> SetLabelFont(fTxt);
  hResIn  -> GetXaxis() -> SetLabelSize(fLabelXY);
  hResIn  -> GetXaxis() -> CenterTitle(fCnt);
  hResIn  -> GetYaxis() -> SetRangeUser(xyRange[1], xyRange[3]);
  hResIn  -> GetYaxis() -> SetTitle(sAxis[1].Data());
  hResIn  -> GetYaxis() -> SetTitleFont(fTxt);
  hResIn  -> GetYaxis() -> SetTitleOffset(fOffY);
  hResIn  -> GetYaxis() -> SetLabelFont(fTxt);
  hResIn  -> GetYaxis() -> SetLabelSize(fLabelXY);
  hResIn  -> GetYaxis() -> CenterTitle(fCnt);
  hResIn  -> GetZaxis() -> SetRangeUser(minimum, maximum);
  hResIn  -> GetZaxis() -> SetLabelFont(fTxt);
  hResIn  -> GetZaxis() -> SetLabelSize(fLabelZ);
  hResOut -> SetTitle(sTitleM[1].Data());
  hResOut -> SetTitleFont(fTxt);
  hResOut -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
  hResOut -> GetXaxis() -> SetTitle(sAxis[0].Data());
  hResOut -> GetXaxis() -> SetTitleFont(fTxt);
  hResOut -> GetXaxis() -> SetTitleOffset(fOffX);
  hResOut -> GetXaxis() -> SetLabelFont(fTxt);
  hResOut -> GetXaxis() -> SetLabelSize(fLabelXY);
  hResOut -> GetXaxis() -> CenterTitle(fCnt);
  hResOut -> GetYaxis() -> SetRangeUser(xyRange[1], xyRange[3]);
  hResOut -> GetYaxis() -> SetTitle(sAxis[1].Data());
  hResOut -> GetYaxis() -> SetTitleFont(fTxt);
  hResOut -> GetYaxis() -> SetTitleOffset(fOffY);
  hResOut -> GetYaxis() -> SetLabelFont(fTxt);
  hResOut -> GetYaxis() -> SetLabelSize(fLabelXY);
  hResOut -> GetYaxis() -> CenterTitle(fCnt);
  hResOut -> GetZaxis() -> SetRangeUser(minimum, maximum);
  hResOut -> GetZaxis() -> SetLabelFont(fTxt);
  hResOut -> GetZaxis() -> SetLabelSize(fLabelZ);
  cout << "    Set response matrix styles." << endl;

  // set ratio styles
  const Float_t fTitlR(0.074);
  const Float_t fLabR(0.056);
  const Float_t fOffXR(1.1);
  const Float_t fOffYR(0.7);
  hRatio  -> GetXaxis() -> SetTitleSize(fTitlR);
  hRatio  -> GetXaxis() -> SetTitleOffset(fOffXR);
  hRatio  -> GetXaxis() -> SetLabelSize(fLabR);
  hRatio  -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatio  -> GetYaxis() -> SetTitleSize(fTitlR);
  hRatio  -> GetYaxis() -> SetTitleOffset(fOffYR);
  hRatio  -> GetYaxis() -> SetLabelSize(fLabR);
  hRatioM -> SetTitle(sRatio.Data());
  hRatioM -> SetTitleFont(fTxt);
  hRatioM -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
  hRatioM -> GetXaxis() -> SetTitle(sAxis[0].Data());
  hRatioM -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioM -> GetXaxis() -> SetTitleOffset(fOffX);
  hRatioM -> GetXaxis() -> SetLabelFont(fTxt);
  hRatioM -> GetXaxis() -> SetLabelSize(fLabelXYR);
  hRatioM -> GetXaxis() -> CenterTitle(fCnt);
  hRatioM -> GetYaxis() -> SetRangeUser(xyRange[1], xyRange[3]);
  hRatioM -> GetYaxis() -> SetTitle(sAxis[1].Data());
  hRatioM -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioM -> GetYaxis() -> SetTitleOffset(fOffY);
  hRatioM -> GetYaxis() -> SetLabelFont(fTxt);
  hRatioM -> GetYaxis() -> SetLabelSize(fLabelXYR);
  hRatioM -> GetYaxis() -> CenterTitle(fCnt);
  hRatioM -> GetZaxis() -> SetLabelFont(fTxt);
  hRatioM -> GetZaxis() -> SetLabelSize(fLabelZR);
  cout << "    Set ratio styles." << endl;


  // make profiles
  TProfile *pResIn  = hResIn - > ProfileX("pResIn", 1, -1, "S");
  TProfile *pResOut = hResOut -> ProfileX("pResOut", 1, -1, "S");

  // set profile styles
  pResIn  -> SetFillColor(cProf[0]);
  pResIn  -> SetFillStyle(fFilP[0]);
  pResIn  -> SetLineColor(cProf[0]);
  pResIn  -> SetLineStyle(fLinP[0]);
  pResIn  -> SetMarkerColor(cProf[0]);
  pResIn  -> SetMarkerStyle(mProf[0]);
  pResIn  -> SetTitle(sProfTitle.Data());
  pResIn  -> SetTitleFont(fTxt);
  pResIn  -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
  pResIn  -> GetXaxis() -> SetTitle(sAxis[0].Data());
  pResIn  -> GetXaxis() -> SetTitleFont(fTxt);
  pResIn  -> GetXaxis() -> SetTitleOffset(fOffX);
  pResIn  -> GetXaxis() -> SetLabelFont(fTxt);
  pResIn  -> GetXaxis() -> SetLabelSize(fLabelXYR);
  pResIn  -> GetXaxis() -> CenterTitle(fCnt);
  pResIn  -> GetYaxis() -> SetRangeUser(xyRange[2], xyRange[3]);
  pResIn  -> GetYaxis() -> SetTitle(sAxis[1].Data());
  pResIn  -> GetYaxis() -> SetTitleFont(fTxt);
  pResIn  -> GetYaxis() -> SetTitleOffset(fOffY);
  pResIn  -> GetYaxis() -> SetLabelFont(fTxt);
  pResIn  -> GetYaxis() -> SetLabelSize(fLabelXYR);
  pResIn  -> GetYaxis() -> CenterTitle(fCnt);
  pResIn  -> SetFillColor(cProf[1]);
  pResIn  -> SetFillStyle(fFilP[1]);
  pResIn  -> SetLineColor(cProf[1]);
  pResIn  -> SetLineStyle(fLinP[1]);
  pResOut -> SetMarkerColor(cProf[1]);
  pResOut -> SetMarkerStyle(mProf[1]);
  pResOut -> SetTitle(sProfTitle.Data());
  pResOut -> SetTitleFont(fTxt);
  pResOut -> GetXaxis() -> SetRangeUser(xyRange[0], xyRange[2]);
  pResOut -> GetXaxis() -> SetTitle(sAxis[0].Data());
  pResOut -> GetXaxis() -> SetTitleFont(fTxt);
  pResOut -> GetXaxis() -> SetTitleOffset(fOffX);
  pResOut -> GetXaxis() -> SetLabelFont(fTxt);
  pResOut -> GetXaxis() -> SetLabelSize(fLabelXYR);
  pResOut -> GetXaxis() -> CenterTitle(fCnt);
  pResOut -> GetYaxis() -> SetRangeUser(xyRange[2], xyRange[3]);
  pResOut -> GetYaxis() -> SetTitle(sAxis[1].Data());
  pResOut -> GetYaxis() -> SetTitleFont(fTxt);
  pResOut -> GetYaxis() -> SetTitleOffset(fOffY);
  pResOut -> GetYaxis() -> SetLabelFont(fTxt);
  pResOut -> GetYaxis() -> SetLabelSize(fLabelXYR);
  pResOut -> GetYaxis() -> CenterTitle(fCnt);

  // profiles for plot
  TProfile *pPlotIn  = (TProfile*) pResIn -> Clone();
  TProfile *pPlotOut = (TProfile*) pResOut -> Clone();
  pPlotIn  -> SetName("pResForPlotIn");
  pPlotIn  -> SetFillColor(cPlot[0]);
  pPlotIn  -> SetLineColor(cPlot[0]);
  pPlotIn  -> SetMarkerColor(cPlot[0]);
  pPlotOut -> SetName("pResForPlotOut");
  pPlotOut -> SetFillColor(cPlot[1]);
  pPlotOut -> SetLineColor(cPlot[1]);
  pPlotOut -> SetMarkerColor(cPlot[1]);
  cout << "    Made profiles and set styles." << endl;


  // make legends and text box
  const UInt_t  fColL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fAlnL(12);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  leg -> SetFillColor(fColL);
  leg -> SetFillStyle(fFilL);
  leg -> SetLineColor(fColL);
  leg -> SetLineStyle(fLinL);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnL);
  leg -> AddEntry(hParIn, sLegPar.Data());
  leg -> AddEntry(hDetOut, sLegDet.Data());

  TLegend *legP = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legP -> SetFillColor(fColL);
  legP -> SetFillStyle(fFilL);
  legP -> SetLineColor(fColL);
  legP -> SetLineStyle(fLinL);
  legP -> SetTextFont(fTxt);
  legP -> SetTextAlign(fAlnL);
  legP -> AddEntry(pPlotIn, sLegM[0].Data());
  legP -> AddEntry(pPlotOut, sLegM[1].Data());

  TPaveText *ptTxt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  ptTxt -> SetFillColor(fColL);
  ptTxt -> SetFillStyle(fFilL);
  ptTxt -> SetLineColor(fColL);
  ptTxt -> SetLineStyle(fLinL);
  ptTxt -> SetTextFont(fTxt);
  ptTxt -> SetTextAlign(fAlnL);
  ptTxt -> AddText(sSys.Data());
  ptTxt -> AddText(sTrg.Data());
  ptTxt -> AddText(sJet.Data());
  ptTxt -> AddText(sTyp.Data());
  cout << "    Made legends and text box." << endl;


  // make line
  const UInt_t  fColLi(1);
  const UInt_t  fStyLi(9);
  const Float_t yLine(1.);

  TLine *line = new TLine(xRange[0], yLine, xRange[1], yLine);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fStyLi);
  cout << "    Made line." << endl;


  // canvas parameters
  const UInt_t  fWidth(750);
  const UInt_t  fHeight(750);
  const UInt_t  fWidthRes(1500);
  const UInt_t  fHeightRes(1500);
  const UInt_t  fBigHeight(950);
  const UInt_t  fLogY(1);
  const UInt_t  fLogZ(1);
  const UInt_t  fGrid(0);
  const Float_t fMargin(0.02);
  const Float_t fSmallMargin(0.005);
  const Float_t fBigMargin(0.2);
  const Float_t xyNum[NVtx] = {0., 0.5, 0.5, 1.};
  const Float_t xyDen[NVtx] = {0.5, 0.5, 1., 1.};
  const Float_t xyDiv[NVtx] = {0., 0., 0.5, 0.5};
  const Float_t xyPro[NVtx] = {0.5, 0., 1., 0.5};
  const Float_t xyPad[NPad][NVtx] = {{0., 0., 1., 0.35}, {0., 0.35, 1., 1.}};
  const TString sPadNames[NPad]   = {"pRatio", "pDists"};

  // make plots
  TCanvas *cDist  = new TCanvas("cDistributions", "", fWidth, fBigHeight);
  TPad    *pDistR = new TPad(sPadNames[0].Data(), "", xyPad[0][0], xyPad[0][1], xyPad[0][2], xyPad[0][3]);
  TPad    *pDistD = new TPad(sPadNames[1].Data(), "", xyPad[1][0], xyPad[1][1], xyPad[1][2], xyPad[1][3]);
  pDistR  -> SetGrid(fGrid, fGrid);
  pDistR  -> SetTopMargin(fSmallMargin);
  pDistR  -> SetRightMargin(fMargin);
  pDistR  -> SetBottomMargin(fBigMargin);
  pDistD  -> SetGrid(fGrid, fGrid);
  pDistD  -> SetLogy(fLogY);
  pDistD  -> SetTopMargin(fMargin);
  pDistD  -> SetRightMargin(fMargin);
  pDistD  -> SetBottomMargin(fSmallMargin);
  cDist   -> cd();
  pDistR  -> Draw();
  pDistD  -> Draw();
  pDistR  -> cd();
  hRatio  -> Draw();
  line    -> Draw();
  pDistD  -> cd();
  hParIn  -> Draw();
  hDetOut -> Draw("sames");
  leg     -> Draw();
  ptTxt   -> Draw();
  fOut    -> cd();
  cDist   -> Write();
  cDist   -> Close();

  TCanvas *cMatrix = new TCanvas("cMatrices", "", fWidthRes, fHeightRes);
  TPad    *pNum      = new TPad("pNum", "", xyNum[0], xyNum[1], xyNum[2], xyNum[3]);
  TPad    *pDen      = new TPad("pDen", "", xyDen[0], xyDen[1], xyDen[2], xyDen[3]);
  TPad    *pDiv      = new TPad("pDiv", "", xyDiv[0], xyDiv[1], xyDiv[2], xyDiv[3]);
  TPad    *pPro      = new TPad("pPro", "", xyPro[0], xyPro[1], xyPro[2], xyPro[3]);
  pNum     -> SetGrid(fGrid, fGrid);
  pNum     -> SetLogz(fLogZ);
  pDen     -> SetGrid(fGrid, fGrid);
  pDen     -> SetLogz(fLogZ);
  pDiv     -> SetGrid(fGrid, fGrid);
  pPro     -> SetGrid(fGrid, fGrid);
  cMatrix  -> cd();
  pNum     -> Draw();
  pDen     -> Draw();
  pDiv     -> Draw();
  pPro     -> Draw();
  pNum     -> cd();
  hResIn   -> Draw("colz");
  pResIn   -> Draw("same");
  pDen     -> cd();
  hResOut  -> Draw("colz");
  pResOut  -> Draw("same");
  pDiv     -> cd();
  hRatioM  -> Draw("colz");
  pPro     -> cd();
  pPlotIn  -> Draw();
  pPlotOut -> Draw("sames");
  legP     -> Draw();
  ptTxt    -> Draw();
  fOut     -> cd();
  cMatrix  -> Write();
  cMatrix  -> Close();
  cout << "    Made plots." << endl;


  // save histograms
  fOut     -> cd();
  hParIn   -> Write();
  hEffIn   -> Write();
  hResIn   -> Write();
  pResIn   -> Write();
  pPlotIn  -> Write();
  hSample  -> Write();
  hDetOut  -> Write();
  hResOut  -> Write();
  pResOut  -> Write();
  pPlotOut -> Write();
  hRatioM  -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut   -> cd();
  fOut   -> Close();
  fParIn -> cd();
  fParIn -> Close();
  fEffIn -> cd();
  fEffIn -> Close();
  fResIn -> cd();
  fResIn -> Close();
  cout << "  Finished quick smearing macro!\n" << endl;

}

// End ------------------------------------------------------------------------
