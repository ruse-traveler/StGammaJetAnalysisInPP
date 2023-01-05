// 'QuickResponseSmoothTest.C'
// Derek Anderson
// 06.03.2021
//
// Quick test macro to try out
// response matrix smoothing
// scheme.


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
static const UInt_t NExp(2);
static const UInt_t NPts(2);
static const UInt_t NPar(2);
static const UInt_t NPad(2);
static const UInt_t NVtx(4);
static const UInt_t NHist(2);
static const UInt_t NRange(2);
static const UInt_t NPtPar(5);
static const UInt_t NQtJet(6);
static const UInt_t NParam(6);

// calculation parameters
static const UInt_t NMcIter(1000);
//static const UInt_t NMcIter(5000000);
static const Bool_t ApplyEff(false);



void QuickResponseSmoothTest() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning quick response matrix smoothing test..." << endl;

  // io parameters
  const TString sOut("quickResponseSmoothTestRFF_forClosureTestWithSmoothHists.et911r05qt05130pi0.d29m6y2021.root");
  const TString sResIn("input/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d15m6y2021.root");
  const TString sEffIn("input/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d15m6y2021.root");
  const TString sParIn("input/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d15m6y2021.root");
  const TString sDetIn("input/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d15m6y2021.root");
  const TString sRes("hResponseRFF");
  const TString sEff("hEfficiencyAll");
  const TString sPar("hSumParRFF");
  const TString sDet("hSumDetRFF");

  // histogram parameters
  const TString sTitle("");
  const TString sTitleY("a. u.");
  const TString sTitleRY("ratio");
  const TString sTitlePX("p_{T}^{par} [GeV/c]");
  const TString sTitleDX("p_{T}^{det} [GeV/c]");
  const TString sLegCheck("p_{T}^{det} sampled from embed. matrix");
  const TString sLegParIn("p_{T}^{par} from embedding");
  const TString sLegDetIn("p_{T}^{det} from embedding");
  const TString sLegParOut("p_{T}^{par} sampled from prior");
  const TString sLegDetOut("p_{T}^{det} from sampling p_{T}^{par} and q_{T}^{jet}");
  const TString sParName[NPtPar] = {"hPtPar_pt0206", "hPtPar_pt061", "hPtPar_pt12", "hPtPar_pt210", "hPtPar_pt1057"};
  const TString sDetName[NPtPar] = {"hPtDet_pt0206", "hPtDet_pt061", "hPtDet_pt12", "hPtDet_pt210", "hPtDet_pt1057"};
  const Float_t xRange[NRange]   = {-1., 57.};
  const Float_t fLbl(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.15);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const UInt_t  fCol[NHist]   = {923, 899};
  const UInt_t  fColP[NPtPar] = {923, 819, 859, 899, 879};
  const UInt_t  fColD[NPtPar] = {923, 819, 859, 899, 879};
  const UInt_t  fMar[NHist]   = {20, 25};
  const UInt_t  fMarP[NPtPar] = {24, 26, 32, 25, 28};
  const UInt_t  fMarD[NPtPar] = {24, 26, 32, 25, 28};
  const UInt_t  fLin[NHist]   = {1, 1};
  const UInt_t  fFil[NHist]   = {0, 0};

  // resposne smoothing parameters
  const TString sFunc("gaus(0) + gaus(3)");
  const TString sQtNames[NPtPar]         = {"fQt_pt0206", "fQt_pt061", "fQt_pt12", "fQt_pt210", "fQt_pt1057"};
  const Float_t xFuncDef[NRange]         = {0., 2.};
  const Float_t pTdetRange[NRange]       = {0.2, 30.};
  const Float_t xPtParLo[NPtPar]         = {0.2, 0.6, 1., 2., 10.};
  const Float_t xPtParHi[NPtPar]         = {0.6, 1., 2., 10., 57.};
  const Float_t xQtJetLo[NQtJet]         = {0.5, 0.7, 0.8, 0.9, 1., 1.1};
  const Float_t xQtJetHi[NQtJet]         = {0.7, 0.8, 0.9, 1., 1.1, 1.3};
  const Float_t fFuncPar[NPtPar][NParam] = {{0.59, 1.01, 0.06, 0.007, 0.62, 0.06}, {0.68, 1.02, 0.06, 0.04, 0.67, 0.06}, {0.57, 1.02, 0.07, 0.06, 0.68, 0.06}, {0.47, 1.01, 0.07, 0.15, 0.73, 0.08}, {0.48, 0.99, 0.07, 0.12, 0.71, 0.08}};  // 9 - 11 GeV parameters

  // prior smoothing parameters
  const TString sFitHistExp("hFitPrior_onlyExpos");
  const TString sFitHistTot("hFitPrior_totalFit");
  const TString sTotName("fTotalFit");
  const TString sTotFunc("(expo(0) + expo(2)) * TMath::TanH((x-[4])/[5])");
  const TString sExpNames[NExp]       = {"fExpSoft", "fExpHard"};
  const TString sExpFuncs[NExp]       = {"expo(0)", "expo(0)"};
  const TString sParTanhNames[NPar]   = {"p_{T}^{0}", "a"};
  const Float_t xFuncTotRange[NRange] = {0.2, 57.};
  const Float_t xExpRange[NExp][NPts] = {{0.6, 2.}, {2., 40.}};
  const Float_t pExpGuess[NExp][NPar] = {{0.25, -1.03}, {-1.76, -0.26}};
  const Float_t pTanhGuess[NPar]      = {2., 1.};
  const UInt_t  fColExp[NExp]         = {2, 4};
  const UInt_t  fColTot(6);

  // efficiency fit parameters
  const TString sEffName("fEfficiency");
  const TString sEffTotName("fEffTotalFit");
  const TString sEffFunc("[0] * (1 - TMath::Exp(-1.*[1]*x))");
  const TString sEffTot("(expo(0) + expo(2)) * TMath::TanH((x-[4])/[5]) * [6] * (1 - TMath::Exp(-1.*[7]*x))");
  const TString sEffParNames[NPar] = {"#epsilon_{0}", "#sigma"};
  const Float_t xEffRange[NPts]    = {0., 100.};
  const Float_t xFitRange[NPts]    = {0.2, 15.};
  const Float_t pEffGuess[NPar]    = {0.87, 4.};
  const UInt_t  fColEff(2);

  // text parameters
  const TString sSys("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trigger, E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.5");
  const TString sTyp("#bf{charged jets}");
  const TString sLegPt[NPtPar] = {"p_{T}^{par} #in (0.2, 0.6) GeV/c", "p_{T}^{par} #in (0.6, 1) GeV/c", "p_{T}^{par} #in (1, 2) GeV/c", "p_{T}^{par} #in (2, 10) GeV/c", "p_{T}^{par} #in (10, 57) GeV/c"};


  // open files
  TFile *fOut   = new TFile(sOut.Data(), "recreate");
  TFile *fResIn = new TFile(sResIn.Data(), "read");
  TFile *fParIn = new TFile(sParIn.Data(), "read");
  TFile *fDetIn = new TFile(sDetIn.Data(), "read");
  if (!fOut || !fResIn || !fParIn || !fDetIn) {
    cerr << "PANIC: couldn't open a file!" << endl;
    return;
  }

  TFile *fEffIn;
  if (ApplyEff) {
    fEffIn = new TFile(sEffIn.Data(), "read");
    if (!fEffIn) {
      cerr << "PANIC: couldn't open efficiency file!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH2D *hResIn = (TH2D*) fResIn -> Get(sRes.Data());
  TH1D *hParIn = (TH1D*) fParIn -> Get(sPar.Data());
  TH1D *hDetIn = (TH1D*) fDetIn -> Get(sDet.Data());
  if (!hResIn || !hParIn || !hDetIn) {
    cerr << "PANIC: couldn't grab a histogram!" << endl;
    return;
  }
  hResIn -> SetName("hResponseInput");
  hParIn -> SetName("hParticleInput");
  hDetIn -> SetName("hDetectorInput");

  TH1D *hEffIn;
  if (ApplyEff) {
    hEffIn = (TH1D*) fEffIn -> Get(sEff.Data());
    if (!hEffIn) {
      cerr << "PANIC: couldn't grab efficiency histogram!" << endl;
      return;
    }
    hEffIn -> SetName("hEfficiencyInput");
  }
  cout << "    Grabbed histograms." << endl;


  // create smoothing functions
  TF1 *fSmooth[NPtPar];
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    fSmooth[iPtPar] = new TF1(sQtNames[iPtPar].Data(), sFunc.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParam = 0; iParam < NParam; iParam++) {
      fSmooth[iPtPar] -> SetParameter(iParam, fFuncPar[iPtPar][iParam]);
    }
  }  // end pTPar bin loop

  // create histograms
  TH1D *hPtPar[NPtPar];
  TH1D *hPtDet[NPtPar];
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    hPtPar[iPtPar] = (TH1D*) hParIn -> Clone();
    hPtDet[iPtPar] = (TH1D*) hDetIn -> Clone();
    hPtPar[iPtPar] -> SetName(sParName[iPtPar].Data());
    hPtDet[iPtPar] -> SetName(sDetName[iPtPar].Data());
    hPtPar[iPtPar] -> Reset("ICES");
    hPtDet[iPtPar] -> Reset("ICES");
  }
  cout << "    Created smoothing functions and histograms." << endl;


  // declare histograms for fitting prior
  TH1D *hFitExp = (TH1D*) hParIn -> Clone();
  TH1D *hFitTot = (TH1D*) hParIn -> Clone();
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
    }
    else {
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
  TF1 *fTot = new TF1(sTotName.Data(), sTotFunc.Data(), xFuncTotRange[0], xFuncTotRange[1]);
  fTot -> SetLineColor(fColTot);

  const UInt_t iParTanh = 2 * NExp;
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    const Float_t par0    = pExpFit[iExp][0];
    const Float_t par1    = pExpFit[iExp][1];
    const Float_t par0m   = par0 - eExpFit[iExp][0];
    const Float_t par0p   = par0 + eExpFit[iExp][0];
    const Float_t par1m   = par1 - eExpFit[iExp][1];
    const Float_t par1p   = par1 + eExpFit[iExp][1];
    const UInt_t  iParExp = 2 * iExp;
    fTot -> SetParameter(iParExp + 0, pExpFit[iExp][0]);
    fTot -> SetParameter(iParExp + 1, pExpFit[iExp][1]);
    fTot -> SetParLimits(iParExp + 0, TMath::Min(par0m, par0p), TMath::Max(par0m, par0p));
    fTot -> SetParLimits(iParExp + 1, TMath::Min(par1m, par1p), TMath::Max(par1m, par1p));
  }
  fTot -> SetParameter(iParTanh + 0, pTanhGuess[0]);
  fTot -> SetParameter(iParTanh + 1, pTanhGuess[1]);
  fTot -> SetParName(iParTanh + 0, sParTanhNames[0].Data());
  fTot -> SetParName(iParTanh + 1, sParTanhNames[1].Data());
  cout << "    Initialized total fit function." << endl;

  // fit prior
  hFitTot -> Fit(sTotName.Data(), "RB", "", xExpRange[0][0], xExpRange[NExp - 1][1]);
  cout << "    Fit prior with total function:\n"
       << "      " << sParTanhNames[0].Data() << " = " << fTot -> GetParameter(iParTanh) << ", " << sParTanhNames[1].Data() << " = " << fTot -> GetParameter(iParTanh + 1)
       << endl;


  // get efficiency function
  TF1 *fEff;
  TF1 *fEffTot;
  if (ApplyEff) {
    fEff = new TF1(sEffName.Data(), sEffFunc.Data(), xEffRange[0], xEffRange[1]);
    fEff   -> SetParameter(0, pEffGuess[0]);
    fEff   -> SetParameter(1, pEffGuess[1]);
    fEff   -> SetParName(0, sEffParNames[0].Data());
    fEff   -> SetParName(1, sEffParNames[1].Data());
    fEff   -> SetLineColor(fColEff);
    hEffIn -> Fit(sEffName.Data(), "", "", xFitRange[0], xFitRange[1]);
    cout << "    Fit efficiency." << endl;

    hParIn -> Multiply(fEff, 1.);
    cout << "    Applied efficiency to prior." << endl;

    fEffTot = new TF1(sEffTotName.Data(), sEffTot.Data(), xFuncTotRange[0], xFuncTotRange[1]);
    fEffTot -> SetParameter(0, fTot -> GetParameter(0));
    fEffTot -> SetParameter(1, fTot -> GetParameter(1));
    fEffTot -> SetParameter(2, fTot -> GetParameter(2));
    fEffTot -> SetParameter(3, fTot -> GetParameter(3));
    fEffTot -> SetParameter(4, fTot -> GetParameter(4));
    fEffTot -> SetParameter(5, fTot -> GetParameter(5));
    fEffTot -> SetParameter(6, fEff -> GetParameter(0));
    fEffTot -> SetParameter(7, fEff -> GetParameter(1));
    fEffTot -> SetLineColor(fColTot);
    cout << "    Made prior function with efficiency applied." << endl;
  }

  // determine which function to sample
  TF1 *fParSample;
  if (ApplyEff) {
    fParSample = (TF1*) fEffTot -> Clone();
    fParSample -> SetName("fParSample");
  }
  else {
    fParSample = (TF1*) fTot -> Clone();
    fParSample -> SetName("fParSample");
  }
  cout << "    Selected function to sample." << endl;


  // create smoothed response
  TH2D *hResOut = (TH2D*) hResIn -> Clone();
  TH1D *hParOut = (TH1D*) hParIn -> Clone();
  TH1D *hDetOut = (TH1D*) hDetIn -> Clone();
  hResOut -> SetName("hResponseOutput");
  hParOut -> SetName("hParticleOutput");
  hDetOut -> SetName("hDetectorOutput");
  hResOut -> Reset("ICES");
  hParOut -> Reset("ICES");
  hDetOut -> Reset("ICES");

  // for detector-level check
  TH1D *hDetCheck = (TH1D*) hDetIn -> Clone();
  hDetCheck -> SetName("hDetectorCheck");
  hDetCheck -> Reset("ICES");

  // mc loop
  cout << "    Beginning MC loop:" << endl;
  for (UInt_t iMC = 0; iMC < NMcIter; iMC++) {

    cout << "      Processing iteration " << (iMC + 1) << "/" << NMcIter << "...\r" << flush;
    if ((iMC + 1) == NMcIter) cout << endl;

    // sample from particle distribution
    Bool_t   isAboveMin(false);
    Double_t pTpar(0.);
    do {
      pTpar      = fParSample -> GetRandom();
      isAboveMin = (pTpar > 0.2);
    } while (!isAboveMin);
    hParOut -> Fill(pTpar);

    // determine pTpar bin and fill histograms
    const UInt_t iPtBig = NPtPar - 1;
    for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
      if ((pTpar >= xPtParLo[iPtPar]) && (pTpar < xPtParHi[iPtPar])) {
        const UInt_t   iPar  = hResIn -> GetYaxis() -> FindBin(pTpar);
        const Double_t pTchk = hResIn -> ProjectionX("", iPar, iPar) -> GetRandom();
        const Double_t qTjet = fSmooth[iPtPar] -> GetRandom(xQtJetLo[0], xQtJetHi[NQtJet -1]);
        const Double_t pTdet = pTpar * qTjet;
        hPtPar[iPtPar] -> Fill(pTpar);
        hPtDet[iPtPar] -> Fill(pTdet);
        hResOut        -> Fill(pTdet, pTpar);
        hDetOut        -> Fill(pTdet);
        hDetCheck      -> Fill(pTchk);
        break;
      }
      else if (pTpar > xPtParHi[iPtBig]) {
        const UInt_t   iPar  = hResIn -> GetYaxis() -> FindBin(pTpar);
        const Double_t pTchk = hResIn -> ProjectionX("", iPar, iPar) -> GetRandom();
        const Double_t qTjet = fSmooth[iPtBig] -> GetRandom(xQtJetLo[0], xQtJetHi[NQtJet -1]);
        const Double_t pTdet = pTpar * qTjet;
        hPtPar[iPtBig] -> Fill(pTpar);
        hPtDet[iPtBig] -> Fill(pTdet);
        hResOut        -> Fill(pTdet, pTpar);
        hDetOut        -> Fill(pTdet);
        hDetCheck      -> Fill(pTchk);
        break;
      } 
    }  // end pTpar bin loop

  }  // end mc loop
  cout << "    Finished MC loop." << endl;


  // normalize response
  const UInt_t nResX = hResOut -> GetNbinsX();
  const UInt_t nResY = hResOut -> GetNbinsY();
  for (UInt_t iResY = 1; iResY < (nResY + 1); iResY++) {
    const Double_t intResY = hResOut -> Integral(1, nResX, iResY, iResY);
    for (UInt_t iResX = 1; iResX < (nResX + 1); iResX++) {
      const Double_t resBin  = hResOut -> GetBinContent(iResX, iResY);
      const Double_t resErr  = hResOut -> GetBinError(iResX, iResY);
      if (intResY > 0.) {
        hResOut -> SetBinContent(iResX, iResY, resBin / intResY);
        hResOut -> SetBinError(iResX, iResY, resErr / intResY);
      }
    }
  }
  cout << "    Normalized smoothed response matrix." << endl;

  // normalize particle and detector histograms
  const UInt_t nPtX = hParIn -> GetNbinsX();
  for (UInt_t iPt = 1; iPt < (nPtX + 1); iPt++) {
    const Double_t binParOut  = hParOut   -> GetBinContent(iPt);
    const Double_t binDetOut  = hDetOut   -> GetBinContent(iPt);
    const Double_t binDetChk  = hDetCheck -> GetBinContent(iPt);
    const Double_t errParOut  = hParOut   -> GetBinError(iPt);
    const Double_t errDetOut  = hDetOut   -> GetBinError(iPt);
    const Double_t errDetChk  = hDetCheck -> GetBinError(iPt);
    const Double_t binWidth = hParOut   -> GetBinWidth(iPt);
    hParOut   -> SetBinContent(iPt, binParOut / binWidth);
    hDetOut   -> SetBinContent(iPt, binDetOut / binWidth);
    hDetCheck -> SetBinContent(iPt, binDetChk / binWidth);
    hParOut   -> SetBinError(iPt, errParOut / binWidth);
    hDetOut   -> SetBinError(iPt, errDetOut / binWidth);
    hDetCheck -> SetBinError(iPt, errDetChk / binWidth);
  }

  const UInt_t   iNormStart = hParOut   -> FindFirstBinAbove(0.);
  const UInt_t   iNormStop  = hParOut   -> FindLastBinAbove(0.);
  const Double_t intParIn   = hParIn    -> Integral(iNormStart, iNormStop);
  const Double_t intParOut  = hParOut   -> Integral();
  const Double_t intNorm    = intParIn / intParOut;
  hParOut   -> Scale(intNorm);
  hDetOut   -> Scale(intNorm);
  hDetCheck -> Scale(intNorm);
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    const UInt_t   iPtNormStart = hPtPar[iPtPar] -> FindFirstBinAbove(0.);
    const UInt_t   iPtNormStop  = hPtPar[iPtPar] -> FindLastBinAbove(0.);
    const Double_t intPtParIn   = hParIn         -> Integral(iPtNormStart, iPtNormStop);
    const Double_t intPtParOut  = hPtPar[iPtPar] -> Integral();
    const Double_t intPtNorm    = intPtParIn / intPtParOut;
    if (intPtParOut > 0.) {
      hPtPar[iPtPar] -> Scale(intPtNorm);
      hPtDet[iPtPar] -> Scale(intPtNorm);
    }
  }
  cout << "    Normalized particle and detector distributions." << endl;


  // set particle histogram styles
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
  hParIn  -> GetXaxis() -> SetTitle(sTitlePX.Data());
  hParIn  -> GetXaxis() -> SetTitleFont(fTxt);
  hParIn  -> GetXaxis() -> SetTitleOffset(fOffX);
  hParIn  -> GetXaxis() -> CenterTitle(fCnt);
  hParIn  -> GetYaxis() -> SetLabelSize(fLbl);
  hParIn  -> GetYaxis() -> SetLabelFont(fTxt);
  hParIn  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParIn  -> GetYaxis() -> SetTitleFont(fTxt);
  hParIn  -> GetYaxis() -> SetTitleOffset(fOffY);
  hParIn  -> GetYaxis() -> CenterTitle(fCnt);
  hParOut -> SetFillColor(fCol[1]);
  hParOut -> SetFillStyle(fFil[1]);
  hParOut -> SetLineColor(fCol[1]);
  hParOut -> SetLineStyle(fLin[1]);
  hParOut -> SetMarkerColor(fCol[1]);
  hParOut -> SetMarkerStyle(fMar[1]);
  hParOut -> SetTitle(sTitle.Data());
  hParOut -> SetTitleFont(fTxt);
  hParOut -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hParOut -> GetXaxis() -> SetLabelSize(fLbl);
  hParOut -> GetXaxis() -> SetLabelFont(fTxt);
  hParOut -> GetXaxis() -> SetTitle(sTitlePX.Data());
  hParOut -> GetXaxis() -> SetTitleFont(fTxt);
  hParOut -> GetXaxis() -> SetTitleOffset(fOffX);
  hParOut -> GetXaxis() -> CenterTitle(fCnt);
  hParOut -> GetYaxis() -> SetLabelSize(fLbl);
  hParOut -> GetYaxis() -> SetLabelFont(fTxt);
  hParOut -> GetYaxis() -> SetTitle(sTitleY.Data());
  hParOut -> GetYaxis() -> SetTitleFont(fTxt);
  hParOut -> GetYaxis() -> SetTitleOffset(fOffY);
  hParOut -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    hPtPar[iPtPar] -> SetFillColor(fColP[iPtPar]);
    hPtPar[iPtPar] -> SetFillStyle(fFil[1]);
    hPtPar[iPtPar] -> SetLineColor(fColP[iPtPar]);
    hPtPar[iPtPar] -> SetLineStyle(fLin[1]);
    hPtPar[iPtPar] -> SetMarkerColor(fColP[iPtPar]);
    hPtPar[iPtPar] -> SetMarkerStyle(fMarP[iPtPar]);
    hPtPar[iPtPar] -> SetTitle(sTitle.Data());
    hPtPar[iPtPar] -> SetTitleFont(fTxt);
    hPtPar[iPtPar] -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
    hPtPar[iPtPar] -> GetXaxis() -> SetLabelSize(fLbl);
    hPtPar[iPtPar] -> GetXaxis() -> SetLabelFont(fTxt);
    hPtPar[iPtPar] -> GetXaxis() -> SetTitle(sTitlePX.Data());
    hPtPar[iPtPar] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtPar[iPtPar] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPtPar[iPtPar] -> GetXaxis() -> CenterTitle(fCnt);
    hPtPar[iPtPar] -> GetYaxis() -> SetLabelSize(fLbl);
    hPtPar[iPtPar] -> GetYaxis() -> SetLabelFont(fTxt);
    hPtPar[iPtPar] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPtPar[iPtPar] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtPar[iPtPar] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPtPar[iPtPar] -> GetYaxis() -> CenterTitle(fCnt);
  }
  hFitExp  -> SetFillColor(fCol[0]);
  hFitExp  -> SetFillStyle(fFil[0]);
  hFitExp  -> SetLineColor(fCol[0]);
  hFitExp  -> SetLineStyle(fLin[0]);
  hFitExp  -> SetMarkerColor(fCol[0]);
  hFitExp  -> SetMarkerStyle(fMar[0]);
  hFitExp  -> SetTitle(sTitle.Data());
  hFitExp  -> SetTitleFont(fTxt);
  hFitExp  -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hFitExp  -> GetXaxis() -> SetLabelSize(fLbl);
  hFitExp  -> GetXaxis() -> SetLabelFont(fTxt);
  hFitExp  -> GetXaxis() -> SetTitle(sTitlePX.Data());
  hFitExp  -> GetXaxis() -> SetTitleFont(fTxt);
  hFitExp  -> GetXaxis() -> SetTitleOffset(fOffX);
  hFitExp  -> GetXaxis() -> CenterTitle(fCnt);
  hFitExp  -> GetYaxis() -> SetLabelSize(fLbl);
  hFitExp  -> GetYaxis() -> SetLabelFont(fTxt);
  hFitExp  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hFitExp  -> GetYaxis() -> SetTitleFont(fTxt);
  hFitExp  -> GetYaxis() -> SetTitleOffset(fOffY);
  hFitExp  -> GetYaxis() -> CenterTitle(fCnt);
  hFitTot  -> SetFillColor(fCol[0]);
  hFitTot  -> SetFillStyle(fFil[0]);
  hFitTot  -> SetLineColor(fCol[0]);
  hFitTot  -> SetLineStyle(fLin[0]);
  hFitTot  -> SetMarkerColor(fCol[0]);
  hFitTot  -> SetMarkerStyle(fMar[0]);
  hFitTot  -> SetTitle(sTitle.Data());
  hFitTot  -> SetTitleFont(fTxt);
  hFitTot  -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hFitTot  -> GetXaxis() -> SetLabelSize(fLbl);
  hFitTot  -> GetXaxis() -> SetLabelFont(fTxt);
  hFitTot  -> GetXaxis() -> SetTitle(sTitlePX.Data());
  hFitTot  -> GetXaxis() -> SetTitleFont(fTxt);
  hFitTot  -> GetXaxis() -> SetTitleOffset(fOffX);
  hFitTot  -> GetXaxis() -> CenterTitle(fCnt);
  hFitTot  -> GetYaxis() -> SetLabelSize(fLbl);
  hFitTot  -> GetYaxis() -> SetLabelFont(fTxt);
  hFitTot  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hFitTot  -> GetYaxis() -> SetTitleFont(fTxt);
  hFitTot  -> GetYaxis() -> SetTitleOffset(fOffY);
  hFitTot  -> GetYaxis() -> CenterTitle(fCnt);

  // set detector histogram styles
  hDetIn  -> SetFillColor(fCol[0]);
  hDetIn  -> SetFillStyle(fFil[0]);
  hDetIn  -> SetLineColor(fCol[0]);
  hDetIn  -> SetLineStyle(fLin[0]);
  hDetIn  -> SetMarkerColor(fCol[0]);
  hDetIn  -> SetMarkerStyle(fMar[0]);
  hDetIn  -> SetTitle(sTitle.Data());
  hDetIn  -> SetTitleFont(fTxt);
  hDetIn  -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hDetIn  -> GetXaxis() -> SetLabelSize(fLbl);
  hDetIn  -> GetXaxis() -> SetLabelFont(fTxt);
  hDetIn  -> GetXaxis() -> SetTitle(sTitleDX.Data());
  hDetIn  -> GetXaxis() -> SetTitleFont(fTxt);
  hDetIn  -> GetXaxis() -> SetTitleOffset(fOffX);
  hDetIn  -> GetXaxis() -> CenterTitle(fCnt);
  hDetIn  -> GetYaxis() -> SetLabelSize(fLbl);
  hDetIn  -> GetYaxis() -> SetLabelFont(fTxt);
  hDetIn  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDetIn  -> GetYaxis() -> SetTitleFont(fTxt);
  hDetIn  -> GetYaxis() -> SetTitleOffset(fOffY);
  hDetIn  -> GetYaxis() -> CenterTitle(fCnt);
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
  hDetOut -> GetXaxis() -> SetTitle(sTitleDX.Data());
  hDetOut -> GetXaxis() -> SetTitleFont(fTxt);
  hDetOut -> GetXaxis() -> SetTitleOffset(fOffX);
  hDetOut -> GetXaxis() -> CenterTitle(fCnt);
  hDetOut -> GetYaxis() -> SetLabelSize(fLbl);
  hDetOut -> GetYaxis() -> SetLabelFont(fTxt);
  hDetOut -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDetOut -> GetYaxis() -> SetTitleFont(fTxt);
  hDetOut -> GetYaxis() -> SetTitleOffset(fOffY);
  hDetOut -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    hPtDet[iPtPar] -> SetFillColor(fColD[iPtPar]);
    hPtDet[iPtPar] -> SetFillStyle(fFil[1]);
    hPtDet[iPtPar] -> SetLineColor(fColD[iPtPar]);
    hPtDet[iPtPar] -> SetLineStyle(fLin[1]);
    hPtDet[iPtPar] -> SetMarkerColor(fColD[iPtPar]);
    hPtDet[iPtPar] -> SetMarkerStyle(fMarD[iPtPar]);
    hPtDet[iPtPar] -> SetTitle(sTitle.Data());
    hPtDet[iPtPar] -> SetTitleFont(fTxt);
    hPtDet[iPtPar] -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
    hPtDet[iPtPar] -> GetXaxis() -> SetLabelSize(fLbl);
    hPtDet[iPtPar] -> GetXaxis() -> SetLabelFont(fTxt);
    hPtDet[iPtPar] -> GetXaxis() -> SetTitle(sTitleDX.Data());
    hPtDet[iPtPar] -> GetXaxis() -> SetTitleFont(fTxt);
    hPtDet[iPtPar] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPtDet[iPtPar] -> GetXaxis() -> CenterTitle(fCnt);
    hPtDet[iPtPar] -> GetYaxis() -> SetLabelSize(fLbl);
    hPtDet[iPtPar] -> GetYaxis() -> SetLabelFont(fTxt);
    hPtDet[iPtPar] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPtDet[iPtPar] -> GetYaxis() -> SetTitleFont(fTxt);
    hPtDet[iPtPar] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPtDet[iPtPar] -> GetYaxis() -> CenterTitle(fCnt);
  }
  hDetCheck -> SetFillColor(fCol[0]);
  hDetCheck -> SetFillStyle(fFil[0]);
  hDetCheck -> SetLineColor(fCol[0]);
  hDetCheck -> SetLineStyle(fLin[0]);
  hDetCheck -> SetMarkerColor(fCol[0]);
  hDetCheck -> SetMarkerStyle(fMar[0]);
  hDetCheck -> SetTitle(sTitle.Data());
  hDetCheck -> SetTitleFont(fTxt);
  hDetCheck -> GetXaxis() -> SetRangeUser(xRange[0], xRange[1]);
  hDetCheck -> GetXaxis() -> SetLabelSize(fLbl);
  hDetCheck -> GetXaxis() -> SetLabelFont(fTxt);
  hDetCheck -> GetXaxis() -> SetTitle(sTitleDX.Data());
  hDetCheck -> GetXaxis() -> SetTitleFont(fTxt);
  hDetCheck -> GetXaxis() -> SetTitleOffset(fOffX);
  hDetCheck -> GetXaxis() -> CenterTitle(fCnt);
  hDetCheck -> GetYaxis() -> SetLabelSize(fLbl);
  hDetCheck -> GetYaxis() -> SetLabelFont(fTxt);
  hDetCheck -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDetCheck -> GetYaxis() -> SetTitleFont(fTxt);
  hDetCheck -> GetYaxis() -> SetTitleOffset(fOffY);
  hDetCheck -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set histogram styles." << endl;


  // calculate ratios
  TH1D *hRatioPar   = (TH1D*) hParIn    -> Clone();
  TH1D *hRatioDet   = (TH1D*) hDetIn    -> Clone();
  TH1D *hRatioCheck = (TH1D*) hDetCheck -> Clone();
  hRatioPar   -> SetName("hRatioParInVsOut");
  hRatioDet   -> SetName("hRatioDetInVsOut");
  hRatioCheck -> SetName("hRatioDetOutVsCheck");
  hRatioPar   -> Reset("ICES");
  hRatioDet   -> Reset("ICES");
  hRatioCheck -> Reset("ICES");
  hRatioPar   -> Divide(hParIn, hParOut, 1., 1.);
  hRatioDet   -> Divide(hDetIn, hDetOut, 1., 1.);
  hRatioCheck -> Divide(hDetCheck, hDetOut, 1., 1.);

  // set ratio styles
  const Float_t fTitlR(0.074);
  const Float_t fLabR(0.056);
  const Float_t fOffXR(1.1);
  const Float_t fOffYR(0.7);
  hRatioPar   -> GetXaxis() -> SetTitleSize(fTitlR);
  hRatioDet   -> GetXaxis() -> SetTitleSize(fTitlR);
  hRatioCheck -> GetXaxis() -> SetTitleSize(fTitlR);
  hRatioPar   -> GetXaxis() -> SetTitleOffset(fOffXR);
  hRatioDet   -> GetXaxis() -> SetTitleOffset(fOffXR);
  hRatioCheck -> GetXaxis() -> SetTitleOffset(fOffXR);
  hRatioPar   -> GetXaxis() -> SetLabelSize(fLabR);
  hRatioDet   -> GetXaxis() -> SetLabelSize(fLabR);
  hRatioCheck -> GetXaxis() -> SetLabelSize(fLabR);
  hRatioPar   -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioDet   -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioCheck -> GetYaxis() -> SetTitle(sTitleRY.Data());
  hRatioPar   -> GetYaxis() -> SetTitleSize(fTitlR);
  hRatioDet   -> GetYaxis() -> SetTitleSize(fTitlR);
  hRatioCheck -> GetYaxis() -> SetTitleSize(fTitlR);
  hRatioPar   -> GetYaxis() -> SetTitleOffset(fOffYR);
  hRatioDet   -> GetYaxis() -> SetTitleOffset(fOffYR);
  hRatioCheck -> GetYaxis() -> SetTitleOffset(fOffYR);
  hRatioPar   -> GetYaxis() -> SetLabelSize(fLabR);
  hRatioDet   -> GetYaxis() -> SetLabelSize(fLabR);
  hRatioCheck -> GetYaxis() -> SetLabelSize(fLabR);
  cout << "    Calculated ratios." << endl;


  // make legends and text box
  const UInt_t  fColL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fAlnL(12);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  TLegend *legP = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legP -> SetFillColor(fColL);
  legP -> SetFillStyle(fFilL);
  legP -> SetLineColor(fColL);
  legP -> SetLineStyle(fLinL);
  legP -> SetTextFont(fTxt);
  legP -> SetTextAlign(fAlnL);
  legP -> AddEntry(hParIn, sLegParIn.Data());
  legP -> AddEntry(hParOut, sLegParOut.Data());

  TLegend *legD = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legD -> SetFillColor(fColL);
  legD -> SetFillStyle(fFilL);
  legD -> SetLineColor(fColL);
  legD -> SetLineStyle(fLinL);
  legD -> SetTextFont(fTxt);
  legD -> SetTextAlign(fAlnL);
  legD -> AddEntry(hDetIn, sLegDetIn.Data());
  legD -> AddEntry(hDetOut, sLegDetOut.Data());

  TLegend *legDC = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legDC -> SetFillColor(fColL);
  legDC -> SetFillStyle(fFilL);
  legDC -> SetLineColor(fColL);
  legDC -> SetLineStyle(fLinL);
  legDC -> SetTextFont(fTxt);
  legDC -> SetTextAlign(fAlnL);
  legDC -> AddEntry(hDetCheck, sLegCheck.Data());
  legDC -> AddEntry(hDetOut, sLegDetOut.Data());

  TLegend *legPtP = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legPtP -> SetFillColor(fColL);
  legPtP -> SetFillStyle(fFilL);
  legPtP -> SetLineColor(fColL);
  legPtP -> SetLineStyle(fLinL);
  legPtP -> SetTextFont(fTxt);
  legPtP -> SetTextAlign(fAlnL);
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    legPtP -> AddEntry(hPtPar[iPtPar], sLegPt[iPtPar].Data());
  }

  TLegend *legPtD = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  legPtD -> SetFillColor(fColL);
  legPtD -> SetFillStyle(fFilL);
  legPtD -> SetLineColor(fColL);
  legPtD -> SetLineStyle(fLinL);
  legPtD -> SetTextFont(fTxt);
  legPtD -> SetTextAlign(fAlnL);
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    legPtD -> AddEntry(hPtDet[iPtPar], sLegPt[iPtPar].Data());
  }

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


  // make plots
  const UInt_t  fWidth(750);
  const UInt_t  fHeight(750);
  const UInt_t  fBigHeight(950);
  const UInt_t  fLogY(1);
  const UInt_t  fGrid(0);
  const Float_t fMargin(0.02);
  const Float_t fSmallMargin(0.01);
  const Float_t fBigMargin(0.2);
  const Float_t xyPad[NPad][NVtx] = {{0., 0., 1., 0.35}, {0., 0.35, 1., 1.}};
  const TString sPadNames[NPad]   = {"pRatio", "pDists"};

  // make canvases
  TCanvas *cExp  = new TCanvas("cOnlyExpos", "", fWidth, fHeight);
  cExp    -> SetTopMargin(fMargin);
  cExp    -> SetRightMargin(fMargin);
  cExp    -> SetGrid(fGrid, fGrid);
  cExp    -> SetLogy(fLogY);
  cExp    -> cd();
  hFitExp -> Draw();
  ptTxt   -> Draw();
  fOut    -> cd();
  cExp    -> Write();
  cExp    -> Close();

  TCanvas *cTot = new TCanvas("cTotalFit", "", fWidth, fHeight);
  cTot    -> SetTopMargin(fMargin);
  cTot    -> SetRightMargin(fMargin);
  cTot    -> SetGrid(fGrid, fGrid);
  cTot    -> SetLogy(fLogY);
  cTot    -> cd();
  hFitTot -> Draw();
  ptTxt   -> Draw();
  fOut    -> cd();
  cTot    -> Write();
  cTot    -> Close();
  cout << "    Made canvases." << endl;

  TCanvas *cPar  = new TCanvas("cParticleInVsOut", "", fWidth, fBigHeight);
  TPad    *pParR = new TPad(sPadNames[0].Data(), "", xyPad[0][0], xyPad[0][1], xyPad[0][2], xyPad[0][3]);
  TPad    *pParD = new TPad(sPadNames[1].Data(), "", xyPad[1][0], xyPad[1][1], xyPad[1][2], xyPad[1][3]);
  pParR     -> SetGrid(fGrid, fGrid);
  pParR     -> SetTopMargin(fSmallMargin);
  pParR     -> SetRightMargin(fMargin);
  pParR     -> SetBottomMargin(fBigMargin);
  pParD     -> SetGrid(fGrid, fGrid);
  pParD     -> SetLogy(fLogY);
  pParD     -> SetTopMargin(fMargin);
  pParD     -> SetRightMargin(fMargin);
  pParD     -> SetBottomMargin(fSmallMargin);
  cPar      -> cd();
  pParR     -> Draw();
  pParD     -> Draw();
  pParR     -> cd();
  hRatioPar -> Draw();
  line      -> Draw();
  pParD     -> cd();
  hParIn    -> Draw();
  hParOut   -> Draw("sames");
  legP      -> Draw();
  ptTxt     -> Draw();
  fOut      -> cd();
  cPar      -> Write();
  cPar      -> Close();

  TCanvas *cDet  = new TCanvas("cDetectorInVsOut", "", fWidth, fBigHeight);
  TPad    *pDetR = new TPad(sPadNames[0].Data(), "", xyPad[0][0], xyPad[0][1], xyPad[0][2], xyPad[0][3]);
  TPad    *pDetD = new TPad(sPadNames[1].Data(), "", xyPad[1][0], xyPad[1][1], xyPad[1][2], xyPad[1][3]);
  pDetR     -> SetGrid(fGrid, fGrid);
  pDetR     -> SetTopMargin(fSmallMargin);
  pDetR     -> SetRightMargin(fMargin);
  pDetR     -> SetBottomMargin(fBigMargin);
  pDetD     -> SetGrid(fGrid, fGrid);
  pDetD     -> SetLogy(fLogY);
  pDetD     -> SetTopMargin(fMargin);
  pDetD     -> SetRightMargin(fMargin);
  pDetD     -> SetBottomMargin(fSmallMargin);
  cDet      -> cd();
  pDetR     -> Draw();
  pDetD     -> Draw();
  pDetR     -> cd();
  hRatioDet -> Draw();
  line      -> Draw();
  pDetD     -> cd();
  hDetIn    -> Draw();
  hDetOut   -> Draw("sames");
  legD      -> Draw();
  ptTxt     -> Draw();
  fOut      -> cd();
  cDet      -> Write();
  cDet      -> Close();

  TCanvas *cDetChk  = new TCanvas("cDetectorOutVsCheck", "", fWidth, fBigHeight);
  TPad    *pDetChkR = new TPad(sPadNames[0].Data(), "", xyPad[0][0], xyPad[0][1], xyPad[0][2], xyPad[0][3]);
  TPad    *pDetChkD = new TPad(sPadNames[1].Data(), "", xyPad[1][0], xyPad[1][1], xyPad[1][2], xyPad[1][3]);
  pDetChkR    -> SetGrid(fGrid, fGrid);
  pDetChkR    -> SetTopMargin(fSmallMargin);
  pDetChkR    -> SetRightMargin(fMargin);
  pDetChkR    -> SetBottomMargin(fBigMargin);
  pDetChkD    -> SetGrid(fGrid, fGrid);
  pDetChkD    -> SetLogy(fLogY);
  pDetChkD    -> SetTopMargin(fMargin);
  pDetChkD    -> SetRightMargin(fMargin);
  pDetChkD    -> SetBottomMargin(fSmallMargin);
  cDetChk     -> cd();
  pDetChkR    -> Draw();
  pDetChkD    -> Draw();
  pDetChkR    -> cd();
  hRatioCheck -> Draw();
  line        -> Draw();
  pDetChkD    -> cd();
  hDetCheck   -> Draw();
  hDetOut     -> Draw("sames");
  legDC       -> Draw();
  ptTxt       -> Draw();
  fOut        -> cd();
  cDetChk     -> Write();
  cDetChk     -> Close();

  TCanvas *cPtPar = new TCanvas("cParticleVsPtPar", "", fWidth, fHeight);
  cPtPar    -> SetGrid(fGrid, fGrid);
  cPtPar    -> SetLogy(fLogY);
  cPtPar    -> SetTopMargin(fMargin);
  cPtPar    -> SetRightMargin(fMargin);
  cPtPar    -> cd();
  hPtPar[0] -> Draw();
  for (UInt_t iPtPar = 1; iPtPar < NPtPar; iPtPar++) {
    hPtPar[iPtPar] -> Draw("same");
  }
  legPtP -> Draw();
  ptTxt  -> Draw();
  fOut   -> cd();
  cPtPar -> Write();
  cPtPar -> Close();

  TCanvas *cPtDet = new TCanvas("cDetectorVsPtPar", "", fWidth, fHeight);
  cPtDet    -> SetGrid(fGrid, fGrid);
  cPtDet    -> SetLogy(fLogY);
  cPtDet    -> SetTopMargin(fMargin);
  cPtDet    -> SetRightMargin(fMargin);
  cPtDet    -> cd();
  hPtDet[0] -> Draw();
  for (UInt_t iPtPar = 1; iPtPar < NPtPar; iPtPar++) {
    hPtDet[iPtPar] -> Draw("same");
  }
  legPtD -> Draw();
  ptTxt  -> Draw();
  fOut   -> cd();
  cPtDet -> Write();
  cPtDet -> Close();
  cout << "    Made plots." << endl;


  // save histograms
  fOut        -> cd();
  hResIn      -> Write();
  hResOut     -> Write();
  hParIn      -> Write();
  hParOut     -> Write();
  hDetIn      -> Write();
  hDetOut     -> Write();
  hRatioPar   -> Write();
  hRatioDet   -> Write();
  hRatioCheck -> Write();
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    hPtPar[iPtPar]  -> Write();
    hPtDet[iPtPar]  -> Write();
    fSmooth[iPtPar] -> Write();
  }
  hFitExp -> Write();
  hFitTot -> Write();
  for (UInt_t iExp = 0; iExp < NExp; iExp++) {
    fExp[iExp] -> Write();
  }
  fTot -> Write();
  if (ApplyEff) {
    hEffIn  -> Write();
    fEff    -> Write();
    fEffTot -> Write();
  }
  fParSample -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut   -> cd();
  fOut   -> Close();
  fResIn -> cd();
  fResIn -> Close();
  fParIn -> cd();
  fParIn -> Close();
  fDetIn -> cd();
  fDetIn -> Close();
  cout << "  Finished quick response matrix smoothing test!\n" << endl;

}

// End ------------------------------------------------------------------------
