// 'MakeClosureTestPlotWithAverage.C'
// Derek Anderson
// 01.10.2022
//
// Use this to plot a particle-level (with stat
// and sys errors) vs. unfolded variations
// for closure test.
//
// NOTE: use the 'TrigId' and 'TrigBin'
//       variables change the color palette
//       and labels to fit the trigger used.
//       'TrigId' sets the trigger species:
//         TrigId  = 0 -- pi0
//         TrigId  = 1 -- gamma
//       And 'TrgBin' sets the eTtrg bin:
//         TrigBin = 0 -- (9, 11) GeV
//         TrigBin = 1 -- (11, 15) GeV
//         TrigBin = 2 -- (15, 20) GeV

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraphErrors.h"

using namespace std;

// global constants
static const UInt_t NPad(2);
static const UInt_t NPlot(2);
static const UInt_t NData(2);
static const UInt_t NVars(5);
static const UInt_t NTrgIds(2);
static const UInt_t NTrgBins(3);
static const UInt_t NTrgs(NTrgIds * NTrgBins);

// trigger parameters
static const UInt_t TrigId(0);
static const UInt_t TrigBin(2);



void MakeClosureTestPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting closure test..." << endl;

  // io particle parameters
  const TString sOut("closureTestRFF.forDetCutoffRedo_pTbinHuge.et1520r05pi0.d17m5y2022.root");
  const TString sInSys("et1520r05rff_cutoffClosure/summedErrorsRFF.modStats_forClosureTest_pTbinHuge.et1520r05pi0.d26m10y2021.root");
  const TString sInStat("et1520r05rff_cutoffClosure/summedErrorsRFF.modStats_forClosureTest_pTbinHuge.et1520r05pi0.d26m10y2021.root");
  const TString sHistSys("hPlotSysUnfoldSmooth");
  const TString sHistStat("hStatistics");

  // io variation parameters
  const TString sInVar[NVars]   = {"et1520r05rff_cutoffClosure/pp200r9rff.default_detCutButNormalRes_pTbinHuge.et1520r05qt05130.p0m1k3n38t5.root",
                                   "et1520r05rff_cutoffClosure/pp200r9rff.forRegSysM1_detCutButNormalRes_pTbinHuge.et1520r05qt05130.p0m1k2n38t5.root",
                                   "et1520r05rff_cutoffClosure/pp200r9rff.forRegSysP1_detCutButNormalRes_pTbinHuge.et1520r05qt05130.p0m1k4n38t5.root",
                                   "et1520r05rff_cutoffClosure/pp200r9rff.forDefLevySys_detCutButNormalRes_pTbinHuge.et1520r05qt05130.p1m1k3n38t5.root",
                                   "et1520r05rff_cutoffClosure/pp200r9rff.forAltLevySys_detCutButNormalRes_pTbinHuge.et1520r05qt05130.p1m1k3n64t8.root"};
  const TString sHistVar[NVars] = {"hUnfolded", "hUnfolded", "hUnfolded", "hUnfolded", "hUnfolded"};

  // particle plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T}^{par} [GeV/c]");
  const TString  sTitleY("(1/N^{trg}) d^{2}N^{jet}/d(p_{T}^{par} #eta^{jet}) [GeV/c]^{-1}");
  const TString  sNameSys("hParticleWithDataSys");
  const TString  sNameStat("hParticleWithStat");
  const TString  sLabelSys("par.-lvl. [with data #sigma_{sys}^{unfold}]");
  const TString  sLabelStat("par.-lvl. [with modified #sigma_{stat}]");
  const Float_t  xPlotRange[NPlot] = {0., 30.};

  // variation and ratio plot parameters
  const TString sTitleRY("unfolded / par.-lvl.");
  const TString sNameVar[NVars]    = {"hVarDefault", "hVarRegM1", "hVarRegP1", "hVarDefLevy", "hVarAltLevy"};
  const TString sGraphVar[NVars]   = {"gVarDefault", "gVarRegM1", "gVarRegP1", "gVarDefLevy", "gVarAltLevy"};
  const TString sNameRatio[NVars]  = {"hRatDefault", "hRatRegM1", "hRatRegP1", "hRatDefLevy", "hRatAltLevy"};
  const TString sGraphRatio[NVars] = {"gRatDefault", "gRatRegM1", "gRatRegP1", "gRatDefLevy", "gRatAltLevy"};
  const TString sLabelVar[NVars]   = {"unfolded det.-lvl. [default, n^{*}_{iter} = 3]",
                                      "unfolded det.-lvl. [n^{*}_{iter} - 1]",
                                      "unfolded det.-lvl. [n^{*}_{iter} + 1]",
                                      "unfolded det.-lvl. [levy prior: n = 3.8, t = 0.5]",
                                      "unfolded det.-lvl. [levy prior: n = 6.4, t = 0.8]"};
  const UInt_t  fMarSty[NVars]     = {21, 20, 22, 23, 34};
  const UInt_t  fLinSty[NVars]     = {1, 1, 1, 1, 1};
  const UInt_t  fFilStyR[NVars]    = {0, 0, 0, 0, 0};
  const Float_t fPerShiftX[NVars]  = {0.08, 0.18, 0.28, 0.38, 0.48};

  // average parameters
  const UInt_t  fMarAvg(1);
  const UInt_t  fLinAvg(9);
  const UInt_t  fWidAvg(2);
  const UInt_t  fFilAvg(0);
  const TString sAvgLabel("corrected");
  const TString sAvgStat("Unfolding stat. error");
  const TString sAvgSys("Unfolding sys. error");
  const Float_t xCalcRange[NPlot] = {0., 30.};

  // text parameters
  const TString sSys("Py6#oplusGeant(RFF), #sqrt{s} = 200 GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");
  const TString sTest("Py6#oplusGeant(RFF) unfolded w/ Py6#oplusGeant(FF)");
  const TString sTestAvg("p+p closure test");
  const TString sTrgAvg("15 < p_{T}^{trg} < 20 GeV/c");
  const TString sSysAvg("#pi^{0}+jet");

  // define color schemes
  const UInt_t fColSy[NTrgs]       = {593, 425, 409, 625, 609, 593};
  const UInt_t fColOu[NTrgs]       = {604, 436, 420, 636, 620, 604};
  const UInt_t fColSt[NTrgs]       = {923, 923, 923, 923, 923, 923};
  const UInt_t fColV[NTrgs][NVars] = {{602, 883, 618, 893, 634}, {434, 863, 602, 883, 618}, {418, 843, 434, 863, 602}, {634, 803, 402, 813, 418}, {618, 893, 634, 803, 402}, {602, 883, 618, 893, 634}};
  const UInt_t fColR[NTrgs][NVars] = {{602, 883, 618, 893, 634}, {434, 863, 602, 883, 618}, {418, 843, 434, 863, 602}, {634, 803, 402, 813, 418}, {618, 893, 634, 803, 402}, {602, 883, 618, 893, 634}};

  // parse trigger selection
  UInt_t fTrg(0);
  if (TrigId == 0) {
    fTrg = TrigBin;
  } else {
    fTrg = TrigBin + 3;
  }

  // select trigger
  UInt_t  fColSys(0);
  UInt_t  fColOut(0);
  UInt_t  fColStat(0);
  UInt_t  fColVar[NVars];
  UInt_t  fColRatio[NVars];
  TString sTrg("");
  switch(fTrg) {
    case 0:
      sTrg     = "#pi^{0} trig., E_{T}^{trg} #in (9, 11) GeV";
      fColSys  = fColSy[0];
      fColOut  = fColOu[0];
      fColStat = fColSt[0];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[0][iVar];
        fColRatio[iVar] = fColR[0][iVar];
      }
      break;
    case 1:
      sTrg     = "#pi^{0} trig., E_{T}^{trg} #in (11, 15) GeV";
      fColSys  = fColSy[1];
      fColOut  = fColOu[1];
      fColStat = fColSt[1];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[1][iVar];
        fColRatio[iVar] = fColR[1][iVar];
      }
      break;
    case 2:
      sTrg     = "#pi^{0} trig., E_{T}^{trg} #in (15, 20) GeV";
      fColSys  = fColSy[2];
      fColOut  = fColOu[2];
      fColStat = fColSt[2];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[2][iVar];
        fColRatio[iVar] = fColR[2][iVar];
      }
      break;
    case 3:
      sTrg     = "#gamma^{dir} trig., E_{T}^{trg} #in (9, 11) GeV";
      fColSys  = fColSy[3];
      fColOut  = fColOu[3];
      fColStat = fColSt[3];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[3][iVar];
        fColRatio[iVar] = fColR[3][iVar];
      }
      break;
    case 4:
      sTrg     = "#gamma^{dir} trig., E_{T}^{trg} #in (11, 15) GeV";
      fColSys  = fColSy[4];
      fColOut  = fColOu[4];
      fColStat = fColSt[4];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[4][iVar];
        fColRatio[iVar] = fColR[4][iVar];
      }
      break;
    case 5:
      sTrg     = "#gamma^{dir} trig., E_{T}^{trg} #in (15, 20) GeV";
      fColSys  = fColSy[5];
      fColOut  = fColOu[5];
      fColStat = fColSt[5];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[5][iVar];
        fColRatio[iVar] = fColR[5][iVar];
      }
      break;
    default:
      sTrg     = "#pi^{0} trig., E_{T}^{trg} #in (9, 11) GeV";
      fColSys  = fColSy[0];
      fColOut  = fColOu[0];
      fColStat = fColSt[0];
      for (UInt_t iVar = 0; iVar < NVars; iVar++) {
        fColVar[iVar]   = fColV[0][iVar];
        fColRatio[iVar] = fColR[0][iVar];
      }
      break;
  }  // end switch-case

  // open particle files
  TFile *fOut    = new TFile(sOut.Data(), "recreate");
  TFile *fInSys  = new TFile(sInSys.Data(), "read");
  TFile *fInStat = new TFile(sInStat.Data(), "read");
  if (!fOut || !fInSys || !fInStat) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInSys = " << fInSys << ", fInStat" << fInStat
         << endl;
    return;
  }
  cout << "    Opened particle files." << endl;

  // open variation files
  TFile *fInVar[NVars];
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {

    fInVar[iVar] = new TFile(sInVar[iVar], "read");
    if (!fInVar[iVar]) {
      cerr << "PANIC: couldn't open a variation file!\n"
           << "       fInVar[" << iVar << "] = " << fInVar[iVar]
           << endl;
      return;
    }

  }  // end variation loop
  cout << "    Opened variation files." << endl;

  // grab particle histograms
  TH1D *hParticleSys  = (TH1D*) fInSys  -> Get(sHistSys.Data());
  TH1D *hParticleStat = (TH1D*) fInStat -> Get(sHistStat.Data());
  if (!hParticleSys || !hParticleStat) {
    cerr << "PANIC: couldn't grab a particle histogram!\n"
         << "       hParticleSys = " << hParticleSys << ", hParticleStat = " << hParticleStat
         << endl;
    return;
  }
  hParticleSys  -> SetName(sNameSys.Data());
  hParticleStat -> SetName(sNameStat.Data());
  cout << "    Grabbed particle histograms." << endl;

  // grab variation histograms
  TH1D *hVariation[NVars];
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {

    hVariation[iVar] = (TH1D*) fInVar[iVar] -> Get(sHistVar[iVar].Data());
    if (!hVariation[iVar]) {
      cerr << "PANIC: couldn't grab a variation histogram!\n"
           << "        hVariation[" << iVar << "] = " << hVariation[iVar]
           << endl;
      return;
    }
    hVariation[iVar] -> SetName(sNameVar[iVar].Data());

  }  // end variation loop
  cout << "    Grabbed variation histograms." << endl;

  // calculate ratio(s)
  TH1D *hRatioSys[NVars];
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {
    hRatioSys[iVar] = (TH1D*) hParticleSys  -> Clone();
    hRatioSys[iVar] -> Divide(hVariation[iVar], hParticleSys, 1., 1.);
    hRatioSys[iVar] -> SetName(sNameRatio[iVar].Data());
  }
  cout << "    Calculated ratio(s)." << endl;

  // calculate average histograms
  TH1D *hAverage  = (TH1D*) hVariation[0] -> Clone();
  TH1D *hRatioAvg = (TH1D*) hRatioSys[0]  -> Clone();
  hAverage  -> SetName("hVariationAverage");
  hRatioAvg -> SetName("hRatioAverage");
  hAverage  -> Reset("ICES");
  hRatioAvg -> Reset("ICES");

  const UInt_t nBins = hParticleSys -> GetNbinsX();
  for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {

    // check if in range
    const Double_t binCenter = hAverage -> GetBinCenter(iBin);
    const Bool_t   isInRange = ((binCenter > xCalcRange[0]) && (binCenter < xCalcRange[1]));
    if (!isInRange) continue;

    // get averages
    Double_t varAvg(0.);
    Double_t ratAvg(0.);
    for (UInt_t iVar = 0; iVar < NVars; iVar++) {
      const Double_t varVal = hVariation[iVar] -> GetBinContent(iBin);
      const Double_t ratVal = hRatioSys[iVar]  -> GetBinContent(iBin);
      varAvg += varVal;
      ratAvg += ratVal;
    }
    varAvg = varAvg / NVars;
    ratAvg = ratAvg / NVars;

    // get max differences
    Double_t varMaxD(0.);
    Double_t ratMaxD(0.);
    for (UInt_t iVar = 0; iVar < NVars; iVar++) {
      const Double_t varVal  = hVariation[iVar] -> GetBinContent(iBin);
      const Double_t ratVal  = hRatioSys[iVar]  -> GetBinContent(iBin);
      const Double_t varDiff = TMath::Abs(varVal - varAvg);
      const Double_t ratDiff = TMath::Abs(ratVal - ratAvg);
      if (varDiff > varMaxD) varMaxD = varDiff;
      if (ratDiff > ratMaxD) ratMaxD = ratDiff;
    }

    // set bins
    hAverage  -> SetBinContent(iBin, varAvg);
    hRatioAvg -> SetBinContent(iBin, ratAvg);
    hAverage  -> SetBinError(iBin, varMaxD);
    hRatioAvg -> SetBinError(iBin, ratMaxD);
  }
  cout << "    Calculated averages." << endl;

  // create uncertainty bands on unity
  TH1D *hUnitySys    = (TH1D*) hParticleSys  -> Clone();
  TH1D *hUnityStat   = (TH1D*) hParticleStat -> Clone();
  for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
    const Double_t sysVal  = hParticleSys  -> GetBinContent(iBin);
    const Double_t statVal = hParticleStat -> GetBinContent(iBin);
    const Double_t sysErr  = hParticleSys  -> GetBinError(iBin);
    const Double_t statErr = hParticleStat -> GetBinError(iBin);
    if (sysVal > 0.) {
      const Double_t sysPer = sysErr / sysVal;
      hUnitySys -> SetBinContent(iBin, 1.);
      hUnitySys -> SetBinError(iBin, sysPer);
    } else {
      hUnitySys -> SetBinContent(iBin, 1.);
      hUnitySys -> SetBinError(iBin, 0.);
    }

    if (statVal > 0.) {
      const Double_t statPer = statErr / statVal;
      hUnityStat -> SetBinContent(iBin, 1.);
      hUnityStat -> SetBinError(iBin, statPer);
    } else {
      hUnityStat -> SetBinContent(iBin, 1.);
      hUnityStat -> SetBinError(iBin, 0.);
    }
  }  // end bin loop
  hUnitySys  -> SetName("hUnityWithDataSys");
  hUnityStat -> SetName("hUnityWithStat"); 
  cout << "    Calculated uncertainty on unity." << endl;

  // set styles
  const UInt_t  fMarPar(29);
  const UInt_t  fMarRat(1);
  const UInt_t  fFilSys(1001);
  const UInt_t  fFilStat(0);
  const UInt_t  fFilUnityS(1001);
  const UInt_t  fFilVar(0);
  const UInt_t  fLinPar(1);
  const UInt_t  fWidPar(2);
  const UInt_t  fWidVar(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fSiz(1.15);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hParticleSys  -> SetMarkerColor(fColSys);
  hParticleSys  -> SetMarkerStyle(fMarPar);
  hParticleSys  -> SetMarkerSize(fSiz);
  hParticleSys  -> SetFillColor(fColSys);
  hParticleSys  -> SetFillStyle(fFilSys);
  hParticleSys  -> SetLineColor(fColSys);
  hParticleSys  -> SetLineStyle(fLinPar);
  hParticleSys  -> SetLineWidth(fWidPar);
  hParticleSys  -> SetTitle(sTitle.Data());
  hParticleSys  -> SetTitleFont(fTxt);
  hParticleSys  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hParticleSys  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hParticleSys  -> GetXaxis() -> SetTitleFont(fTxt);
  hParticleSys  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hParticleSys  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hParticleSys  -> GetXaxis() -> SetLabelFont(fTxt);
  hParticleSys  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hParticleSys  -> GetXaxis() -> CenterTitle(fCnt);
  hParticleSys  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParticleSys  -> GetYaxis() -> SetTitleFont(fTxt);
  hParticleSys  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hParticleSys  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hParticleSys  -> GetYaxis() -> SetLabelFont(fTxt);
  hParticleSys  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hParticleSys  -> GetYaxis() -> CenterTitle(fCnt);
  hParticleStat -> SetMarkerColor(fColStat);
  hParticleStat -> SetMarkerStyle(fMarPar);
  hParticleStat -> SetMarkerSize(fSiz);
  hParticleStat -> SetFillColor(fColStat);
  hParticleStat -> SetFillStyle(fFilStat);
  hParticleStat -> SetLineColor(fColStat);
  hParticleStat -> SetLineStyle(fLinPar);
  hParticleStat -> SetLineWidth(fWidPar);
  hParticleStat -> SetTitle(sTitle.Data());
  hParticleStat -> SetTitleFont(fTxt);
  hParticleStat -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hParticleStat -> GetXaxis() -> SetTitle(sTitleX.Data());
  hParticleStat -> GetXaxis() -> SetTitleFont(fTxt);
  hParticleStat -> GetXaxis() -> SetTitleSize(fTit[1]);
  hParticleStat -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hParticleStat -> GetXaxis() -> SetLabelFont(fTxt);
  hParticleStat -> GetXaxis() -> SetLabelSize(fLab[1]);
  hParticleStat -> GetXaxis() -> CenterTitle(fCnt);
  hParticleStat -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParticleStat -> GetYaxis() -> SetTitleFont(fTxt);
  hParticleStat -> GetYaxis() -> SetTitleSize(fTit[1]);
  hParticleStat -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hParticleStat -> GetYaxis() -> SetLabelFont(fTxt);
  hParticleStat -> GetYaxis() -> SetLabelSize(fLab[1]);
  hParticleStat -> GetYaxis() -> CenterTitle(fCnt);
  hUnitySys     -> SetMarkerColor(fColSys);
  hUnitySys     -> SetMarkerStyle(fMarRat);
  hUnitySys     -> SetMarkerSize(fSiz);
  hUnitySys     -> SetFillColor(fColSys);
  hUnitySys     -> SetFillStyle(fFilUnityS);
  hUnitySys     -> SetLineColor(fColSys);
  hUnitySys     -> SetLineStyle(fLinPar);
  hUnitySys     -> SetLineWidth(fWidPar);
  hUnitySys     -> SetTitle(sTitle.Data());
  hUnitySys     -> SetTitleFont(fTxt);
  hUnitySys     -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hUnitySys     -> GetXaxis() -> SetTitle(sTitleX.Data());
  hUnitySys     -> GetXaxis() -> SetTitleFont(fTxt);
  hUnitySys     -> GetXaxis() -> SetTitleSize(fTit[0]);
  hUnitySys     -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hUnitySys     -> GetXaxis() -> SetLabelFont(fTxt);
  hUnitySys     -> GetXaxis() -> SetLabelSize(fLab[0]);
  hUnitySys     -> GetXaxis() -> CenterTitle(fCnt);
  hUnitySys     -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hUnitySys     -> GetYaxis() -> SetTitleFont(fTxt);
  hUnitySys     -> GetYaxis() -> SetTitleSize(fTit[0]);
  hUnitySys     -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hUnitySys     -> GetYaxis() -> SetLabelFont(fTxt);
  hUnitySys     -> GetYaxis() -> SetLabelSize(fLab[0]);
  hUnitySys     -> GetYaxis() -> CenterTitle(fCnt);
  hUnityStat    -> SetMarkerColor(fColStat);
  hUnityStat    -> SetMarkerStyle(fMarRat);
  hUnityStat    -> SetMarkerSize(fSiz);
  hUnityStat    -> SetFillColor(fColStat);
  hUnityStat    -> SetFillStyle(fFilStat);
  hUnityStat    -> SetLineColor(fColStat);
  hUnityStat    -> SetLineStyle(fLinPar);
  hUnityStat    -> SetLineWidth(fWidPar);
  hUnityStat    -> SetTitle(sTitle.Data());
  hUnityStat    -> SetTitleFont(fTxt);
  hUnityStat    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hUnityStat    -> GetXaxis() -> SetTitle(sTitleX.Data());
  hUnityStat    -> GetXaxis() -> SetTitleFont(fTxt);
  hUnityStat    -> GetXaxis() -> SetTitleSize(fTit[0]);
  hUnityStat    -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hUnityStat    -> GetXaxis() -> SetLabelFont(fTxt);
  hUnityStat    -> GetXaxis() -> SetLabelSize(fLab[0]);
  hUnityStat    -> GetXaxis() -> CenterTitle(fCnt);
  hUnityStat    -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hUnityStat    -> GetYaxis() -> SetTitleFont(fTxt);
  hUnityStat    -> GetYaxis() -> SetTitleSize(fTit[0]);
  hUnityStat    -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hUnityStat    -> GetYaxis() -> SetLabelFont(fTxt);
  hUnityStat    -> GetYaxis() -> SetLabelSize(fLab[0]);
  hUnityStat    -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {
    hVariation[iVar] -> SetMarkerColor(fColVar[iVar]);
    hVariation[iVar] -> SetMarkerStyle(fMarSty[iVar]);
    hVariation[iVar] -> SetMarkerSize(fSiz);
    hVariation[iVar] -> SetFillColor(fColVar[iVar]);
    hVariation[iVar] -> SetFillStyle(fFilVar);
    hVariation[iVar] -> SetLineColor(fColVar[iVar]);
    hVariation[iVar] -> SetLineStyle(fLinSty[iVar]);
    hVariation[iVar] -> SetLineWidth(fWidVar);
    hVariation[iVar] -> SetTitle(sTitle.Data());
    hVariation[iVar] -> SetTitleFont(fTxt);
    hVariation[iVar] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hVariation[iVar] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hVariation[iVar] -> GetXaxis() -> SetTitleFont(fTxt);
    hVariation[iVar] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hVariation[iVar] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hVariation[iVar] -> GetXaxis() -> SetLabelFont(fTxt);
    hVariation[iVar] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hVariation[iVar] -> GetXaxis() -> CenterTitle(fCnt);
    hVariation[iVar] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hVariation[iVar] -> GetYaxis() -> SetTitleFont(fTxt);
    hVariation[iVar] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hVariation[iVar] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hVariation[iVar] -> GetYaxis() -> SetLabelFont(fTxt);
    hVariation[iVar] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hVariation[iVar] -> GetYaxis() -> CenterTitle(fCnt);
    hRatioSys[iVar]  -> SetMarkerColor(fColRatio[iVar]);
    hRatioSys[iVar]  -> SetMarkerStyle(fMarSty[iVar]);
    hRatioSys[iVar]  -> SetMarkerSize(fSiz);
    hRatioSys[iVar]  -> SetFillColor(fColRatio[iVar]);
    hRatioSys[iVar]  -> SetFillStyle(fFilStyR[iVar]);
    hRatioSys[iVar]  -> SetLineColor(fColRatio[iVar]);
    hRatioSys[iVar]  -> SetLineStyle(fLinSty[iVar]);
    hRatioSys[iVar]  -> SetLineWidth(fWidVar);
    hRatioSys[iVar]  -> SetTitle(sTitle.Data());
    hRatioSys[iVar]  -> SetTitleFont(fTxt);
    hRatioSys[iVar]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatioSys[iVar]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatioSys[iVar]  -> GetXaxis() -> SetTitleFont(fTxt);
    hRatioSys[iVar]  -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatioSys[iVar]  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatioSys[iVar]  -> GetXaxis() -> SetLabelFont(fTxt);
    hRatioSys[iVar]  -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatioSys[iVar]  -> GetXaxis() -> CenterTitle(fCnt);
    hRatioSys[iVar]  -> GetYaxis() -> SetTitle(sTitleRY.Data());
    hRatioSys[iVar]  -> GetYaxis() -> SetTitleFont(fTxt);
    hRatioSys[iVar]  -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatioSys[iVar]  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatioSys[iVar]  -> GetYaxis() -> SetLabelFont(fTxt);
    hRatioSys[iVar]  -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatioSys[iVar]  -> GetYaxis() -> CenterTitle(fCnt);
  }  // end variation loop
  hAverage  -> SetMarkerColor(fColVar[NVars - 1]);
  hAverage  -> SetMarkerStyle(fMarAvg);
  hAverage  -> SetMarkerSize(fSiz);
  hAverage  -> SetFillColor(fColVar[NVars - 1]);
  hAverage  -> SetFillStyle(fFilAvg);
  hAverage  -> SetLineColor(fColVar[NVars - 1]);
  hAverage  -> SetLineStyle(fLinAvg);
  hAverage  -> SetLineWidth(fWidAvg);
  hAverage  -> SetTitle(sTitle.Data());
  hAverage  -> SetTitleFont(fTxt);
  hAverage  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hAverage  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hAverage  -> GetXaxis() -> SetTitleFont(fTxt);
  hAverage  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hAverage  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hAverage  -> GetXaxis() -> SetLabelFont(fTxt);
  hAverage  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hAverage  -> GetXaxis() -> CenterTitle(fCnt);
  hAverage  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hAverage  -> GetYaxis() -> SetTitleFont(fTxt);
  hAverage  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hAverage  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hAverage  -> GetYaxis() -> SetLabelFont(fTxt);
  hAverage  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hAverage  -> GetYaxis() -> CenterTitle(fCnt);
  hRatioAvg -> SetMarkerColor(fColRatio[NVars - 1]);
  hRatioAvg -> SetMarkerStyle(fMarAvg);
  hRatioAvg -> SetMarkerSize(fSiz);
  hRatioAvg -> SetFillColor(fColRatio[NVars - 1]);
  hRatioAvg -> SetFillStyle(fFilAvg);
  hRatioAvg -> SetLineColor(fColRatio[NVars - 1]);
  hRatioAvg -> SetLineStyle(fLinAvg);
  hRatioAvg -> SetLineWidth(fWidAvg);
  hRatioAvg -> SetTitle(sTitle.Data());
  hRatioAvg -> SetTitleFont(fTxt);
  hRatioAvg -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hRatioAvg -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioAvg -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioAvg -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatioAvg -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatioAvg -> GetXaxis() -> SetLabelFont(fTxt);
  hRatioAvg -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatioAvg -> GetXaxis() -> CenterTitle(fCnt);
  hRatioAvg -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioAvg -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioAvg -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatioAvg -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatioAvg -> GetYaxis() -> SetLabelFont(fTxt);
  hRatioAvg -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatioAvg -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // create outlines
  const UInt_t fMarOut(1);
  const UInt_t fFilOut(0);
  const UInt_t fLinOut(1);
  const UInt_t fWidOut(2);

  TH1D *hOutSys = (TH1D*) hParticleSys -> Clone();
  TH1D *hOutRat = (TH1D*) hUnitySys    -> Clone();
  hOutSys -> SetName("hSysOutline");
  hOutRat -> SetName("hRatioOutline");
  hOutSys -> SetFillColor(fColOut);
  hOutRat -> SetFillColor(fColOut);
  hOutSys -> SetFillStyle(fFilOut);
  hOutRat -> SetFillStyle(fFilOut);
  hOutSys -> SetLineColor(fColOut);
  hOutRat -> SetLineColor(fColOut);
  hOutSys -> SetLineStyle(fLinOut);
  hOutRat -> SetLineStyle(fLinOut);
  hOutSys -> SetLineWidth(fWidOut);
  hOutRat -> SetLineWidth(fWidOut);
  hOutSys -> SetMarkerColor(fColOut);
  hOutRat -> SetMarkerColor(fColOut);
  hOutSys -> SetMarkerStyle(fMarOut);
  hOutRat -> SetMarkerStyle(fMarOut);
  cout << "    Created outlines." << endl;

  // create graphs
  TGraphErrors *gVariation[NVars];
  TGraphErrors *gRatioSys[NVars];
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {

    // initialize graphs
    gVariation[iVar] = new TGraphErrors(hVariation[iVar]);
    gRatioSys[iVar]  = new TGraphErrors(hRatioSys[iVar]);
    gVariation[iVar] -> SetName(sGraphVar[iVar].Data());
    gRatioSys[iVar]  -> SetName(sGraphRatio[iVar].Data());

    // shift points
    const UInt_t nPoints = gVariation[iVar] -> GetN();
    for (UInt_t iPoint = 0; iPoint < nPoints; iPoint++) {

      // get points
      Double_t xValueV(0.);
      Double_t xValueR(0.);
      Double_t yValueV(0.);
      Double_t yValueR(0.);
      gVariation[iVar] -> GetPoint(iPoint, xValueV, yValueV);
      gRatioSys[iVar]  -> GetPoint(iPoint, xValueR, yValueR);

      // set new points
      const UInt_t  iBin   = hVariation[iVar] -> FindBin(xValueV);
      const Float_t wBin   = hVariation[iVar] -> GetBinWidth(iBin);
      const Float_t xShift = fPerShiftX[iVar] * wBin;
      const Float_t xNew   = xValueV - xShift;
      gVariation[iVar] -> SetPoint(iPoint, xNew, yValueV);
      gRatioSys[iVar]  -> SetPoint(iPoint, xNew, yValueR);
    }
  }
  cout << "    Created graphs." << endl;

  // make legends
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.4};
  const Float_t fAvgXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.3};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sTest.Data());
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hParticleStat, sLabelStat.Data(), "pf");
  leg -> AddEntry(hParticleSys, sLabelSys.Data(), "f");
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {
    leg -> AddEntry(hVariation[iVar], sLabelVar[iVar].Data(), "p");
  }  // end variation loop

  TLegend *avg = new TLegend(fAvgXY[0], fAvgXY[1], fAvgXY[2], fAvgXY[3], sTest.Data());
  avg -> SetFillColor(fColLe);
  avg -> SetFillStyle(fFilLe);
  avg -> SetLineColor(fColLe);
  avg -> SetLineStyle(fLinLe);
  avg -> SetTextFont(fTxt);
  avg -> SetTextAlign(fAln);
  avg -> AddEntry(hParticleStat, sAvgStat.Data(), "pf");
  avg -> AddEntry(hParticleSys, sAvgSys.Data(), "f");
  avg -> AddEntry(hAverage, sAvgLabel.Data(), "l");
  cout << "    Made legends." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const Float_t fTxtXY[NPlot * NPlot]    = {0.3, 0.1, 0.5, 0.3};
  const Float_t fTxtAvgXY[NPlot * NPlot] = {0.3, 0.1, 0.5, 0.3};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());

  TPaveText *txtAvg = new TPaveText(fTxtAvgXY[0], fTxtAvgXY[1], fTxtAvgXY[2], fTxtAvgXY[3], "NDC NB");
  txtAvg -> SetFillColor(fColTx);
  txtAvg -> SetFillStyle(fFilTx);
  txtAvg -> SetLineColor(fColTx);
  txtAvg -> SetLineStyle(fLinTx);
  txtAvg -> SetTextFont(fTxt);
  txtAvg -> SetTextAlign(fAln);
  txtAvg -> AddText(sTestAvg.Data());
  txtAvg -> AddText(sJet.Data());
  txtAvg -> AddText(sTrgAvg.Data());
  txtAvg -> AddText(sSysAvg.Data());
  cout << "    Made text boxes." << endl;

  // make lines
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(2);
  const Float_t fLinXY[NPlot * NPlot] = {xPlotRange[0], 1., xPlotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);

  TH1D *hAvgLine   = (TH1D*) hAverage  -> Clone();
  TH1D *hRatioLine = (TH1D*) hRatioAvg -> Clone();
  hAvgLine   -> SetName("hAverageLine");
  hRatioLine -> SetName("hRatioAverageLine");
  cout << "    Made lines." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT1(0.005);
  const Float_t fMarginT2(0.02);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.005);
  const Float_t fPadXY1[NPlot * NPlot] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NPlot * NPlot] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot         -> SetGrid(fGrid, fGrid);
  cPlot         -> SetTicks(fTick, fTick);
  cPlot         -> SetBorderMode(fMode);
  cPlot         -> SetBorderSize(fBord);
  pPad1         -> SetGrid(fGrid, fGrid);
  pPad1         -> SetTicks(fTick, fTick);
  pPad1         -> SetLogx(fLogX);
  pPad1         -> SetLogy(fLogY);
  pPad1         -> SetBorderMode(fMode);
  pPad1         -> SetBorderSize(fBord);
  pPad1         -> SetFrameBorderMode(fFrame);
  pPad1         -> SetLeftMargin(fMarginL);
  pPad1         -> SetRightMargin(fMarginR);
  pPad1         -> SetTopMargin(fMarginT1);
  pPad1         -> SetBottomMargin(fMarginB1);
  pPad2         -> SetGrid(fGrid, fGrid);
  pPad2         -> SetTicks(fTick, fTick);
  pPad2         -> SetLogx(fLogX);
  pPad2         -> SetLogy(fLogY);
  pPad2         -> SetBorderMode(fMode);
  pPad2         -> SetBorderSize(fBord);
  pPad2         -> SetFrameBorderMode(fFrame);
  pPad2         -> SetLeftMargin(fMarginL);
  pPad2         -> SetRightMargin(fMarginR);
  pPad2         -> SetTopMargin(fMarginT2);
  pPad2         -> SetBottomMargin(fMarginB2);
  cPlot         -> cd();
  pPad1         -> Draw();
  pPad2         -> Draw();
  pPad1         -> cd();
  hUnitySys     -> Draw("E2");
  hOutRat       -> Draw("E2 SAME");
  hUnityStat    -> Draw("SAME");
  gRatioSys[0]  -> Draw("PX");
  for (UInt_t iVar = 1; iVar < NVars; iVar++) {
    gRatioSys[iVar] -> Draw("PX");
  }  // end variation loop
  line          -> Draw();
  pPad2         -> cd();
  hParticleSys  -> Draw("E2");
  hOutSys       -> Draw("E2 SAME");
  hParticleStat -> Draw("SAME");
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {
    gVariation[iVar] -> Draw("PX");
  }  // end variation loop
  leg           -> Draw();
  txt           -> Draw();
  fOut          -> cd();
  cPlot         -> Write();
  cPlot         -> Close();

  TCanvas *cAverage = new TCanvas("cAverage", "", width, height);
  TPad    *pPad3 = new TPad("pPad3", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad4 = new TPad("pPad4", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cAverage      -> SetGrid(fGrid, fGrid);
  cAverage      -> SetTicks(fTick, fTick);
  cAverage      -> SetBorderMode(fMode);
  cAverage      -> SetBorderSize(fBord);
  pPad3         -> SetGrid(fGrid, fGrid);
  pPad3         -> SetTicks(fTick, fTick);
  pPad3         -> SetLogx(fLogX);
  pPad3         -> SetLogy(fLogY);
  pPad3         -> SetBorderMode(fMode);
  pPad3         -> SetBorderSize(fBord);
  pPad3         -> SetFrameBorderMode(fFrame);
  pPad3         -> SetLeftMargin(fMarginL);
  pPad3         -> SetRightMargin(fMarginR);
  pPad3         -> SetTopMargin(fMarginT1);
  pPad3         -> SetBottomMargin(fMarginB1);
  pPad4         -> SetGrid(fGrid, fGrid);
  pPad4         -> SetTicks(fTick, fTick);
  pPad4         -> SetLogx(fLogX);
  pPad4         -> SetLogy(fLogY);
  pPad4         -> SetBorderMode(fMode);
  pPad4         -> SetBorderSize(fBord);
  pPad4         -> SetFrameBorderMode(fFrame);
  pPad4         -> SetLeftMargin(fMarginL);
  pPad4         -> SetRightMargin(fMarginR);
  pPad4         -> SetTopMargin(fMarginT2);
  pPad4         -> SetBottomMargin(fMarginB2);
  cAverage      -> cd();
  pPad3         -> Draw();
  pPad4         -> Draw();
  pPad3         -> cd();
  hOutRat       -> Draw("E2");
  hUnitySys     -> Draw("E2 SAME");
  hUnityStat    -> Draw("E0 X0 SAME");
  hRatioAvg     -> Draw("E5 SAME");
  hRatioLine    -> Draw("HIST ][ L SAME");
  line          -> Draw();
  pPad4         -> cd();
  hOutSys       -> Draw("E2");
  hParticleSys  -> Draw("E2 SAME");
  hParticleStat -> Draw("E0 X0 SAME");
  hAverage      -> Draw("E5 SAME");
  hAvgLine      -> Draw("HIST ][ L SAME");
  avg           -> Draw();
  txtAvg        -> Draw();
  fOut          -> cd();
  cAverage      -> Write();
  cAverage      -> Close();
  cout << "    Made plot." << endl;

  // close files
  fOut          -> cd();
  hParticleSys  -> Write();
  hParticleStat -> Write();
  hUnitySys     -> Write();
  hUnityStat    -> Write();
  hOutSys       -> Write();
  hOutRat       -> Write();
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {
    hVariation[iVar] -> Write();
    gVariation[iVar] -> Write();
    hRatioSys[iVar]  -> Write();
    gRatioSys[iVar]  -> Write();
  }  // end variation loop
  hAverage   -> Write();
  hAvgLine   -> Write();
  hRatioAvg  -> Write();
  hRatioLine -> Write();
  fOut       -> Close();
  fInSys     -> cd();
  fInSys     -> Close();
  fInStat    -> cd();
  fInStat    -> Close();
  for (UInt_t iVar = 0; iVar < NVars; iVar++) {
    fInVar[iVar] -> cd();
    fInVar[iVar] -> Close();
  }  // end variation loop
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
