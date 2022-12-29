// 'FitTsp.C'
// Derek Anderson
// 12.10.2021
//
// Use this to fit a measured TSP distribution
// with a simulated pi0 and photon TSP
// distribution.
//
// NOTE: binning scheme is assumed to be the
//       same between all histograms.

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;

// global constants
static const UInt_t NPlot(2);
static const UInt_t NPar(2);
static const UInt_t NVtx(4);



// main function
void FitTsp() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning TSP fit script..." << endl;

  // io parameters
  const TString sFileO("simVsMeasTsp.withVtxDistribution_morePi0stats.et1520vz55.d25m3y2022.root");
  const TString sFileD("output/2021/December2021/pp200r9data.forTspFit.et1520r02a005chrg.d30m12y2021.root");
  const TString sFileP("input/ParticleGun/particleGun.forTspCheck_withVtxMoreStats_noTsp.et650x650x920vz55tsp0100pi0.d19m3y2022.root");
  const TString sFileG("input/ParticleGun/particleGun.forTspCheck_withVtxNoTsp.et650x650x920vz55tsp0100gam.d22m2y2022.root");
  const TString sHistD("All/hTrgTspA");
  const TString sHistP("match/hTspMatch1520");
  const TString sHistG("match/hTspMatch1520");

  // plot parameters
  const Float_t xPlotRange[NPlot] = {0., 1.};
  const TString sLabelD("measured data");
  const TString sLabelP("simulated #pi^{0} #bf{[particle gun]}");
  const TString sLabelG("simulated #gamma #bf{[particle gun]}");
  const TString sLabelS("weighted sum of #pi^{0} and #gamma");
  const TString sTrg("E_{T}^{trg} #in (11, 15) GeV");

  // fit parameters
  const UInt_t  progSize(500);
  const Float_t stepSize(0.001);
  const Float_t intRatioMin(stepSize);
  const Float_t intRatioMax(3.00);
  const Float_t xLowTsp(0.08);
  const Float_t xHighTsp(0.2);
  const Float_t xPi0Range[NPlot] = {0., 0.08};
  const Float_t xFitRange[NPlot] = {0., 1.};
  const Float_t xFitMax[NPar]    = {-1., -1.};
  const Float_t pFitMin[NPar]    = {stepSize, stepSize};

  // histogram parameters
  const TString sNameD("hData");
  const TString sNameP("hSimPi0");
  const TString sNameG("hSimGam");
  const TString sCompD("hCompDat");
  const TString sCompP("hCompPi0");
  const TString sCompG("hCompGam");
  const TString sCompS("hCompSum");
  const TString sFitP("hFitPi0");
  const TString sFitG("hFitGam");
  const TString sFitS("hFitSum");
  const TString sTitle("");
  const TString sXtitle("TSP");
  const TString sYtitle("a. u.");
  const UInt_t  fColD(923);
  const UInt_t  fColP(859);
  const UInt_t  fColG(899);
  const UInt_t  fColS(879);
  const UInt_t  fMarD(29);
  const UInt_t  fMarP(24);
  const UInt_t  fMarG(25);
  const UInt_t  fMarS(30);

  // binning parameters
  const Bool_t doBinNormD(false);
  const Bool_t doBinNormP(true);
  const Bool_t doBinNormG(true);
  const Bool_t doRebinD(false);
  const Bool_t doRebinP(false);
  const Bool_t doRebinG(false);
  const UInt_t iBinWidth(17);
  const UInt_t fRebin(2);

  // open files
  TFile *fOut  = new TFile(sFileO.Data(), "recreate");
  TFile *fData = new TFile(sFileD.Data(), "read");
  TFile *fSimP = new TFile(sFileP.Data(), "read");
  TFile *fSimG = new TFile(sFileG.Data(), "read");
  if (!fOut || !fData || !fSimP || !fSimG) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut  = " << fOut  << ", fData = " << fData << "\n"
         << "       fSimP = " << fSimP << ", fSimG = " << fSimG << "\n"
         << endl;
    return; 
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hData = (TH1D*) fData -> Get(sHistD.Data());
  TH1D *hSimP = (TH1D*) fSimP -> Get(sHistP.Data());
  TH1D *hSimG = (TH1D*) fSimG -> Get(sHistG.Data());
  if (!hData || !hSimP || !hSimG) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hData = " << hData << ", hSimP = " << hSimP << ", hSimG = " << hSimG << "\n"
         << endl;
    return;
  }
  hData -> SetName(sNameD.Data());
  hSimP -> SetName(sNameP.Data());
  hSimG -> SetName(sNameG.Data());
  cout << "    Grabbed histograms." << endl;

  // normalize for bin width (if need be)
  const Double_t binWidthD = hData -> GetBinWidth(iBinWidth);
  const Double_t binWidthP = hSimP -> GetBinWidth(iBinWidth);
  const Double_t binWidthG = hSimG -> GetBinWidth(iBinWidth);
  if (doBinNormD) hData -> Scale(1. / binWidthD);
  if (doBinNormP) hSimP -> Scale(1. / binWidthP);
  if (doBinNormG) hSimG -> Scale(1. / binWidthG);

  // rebin histograms (if need be)
  if (doRebinD) {
    hData -> Rebin(fRebin);
    hData -> Scale(1. / (Double_t) fRebin);
  }
  if (doRebinP) {
    hSimP -> Rebin(fRebin);
    hSimP -> Scale(1. / (Double_t) fRebin);
  }
  if (doRebinG) {
    hSimG -> Rebin(fRebin);
    hSimG -> Scale(1. / (Double_t) fRebin);
  }
  if (doRebinD || doRebinP || doRebinG) cout << "    Rebinned histograms." << endl;

  // normalize input
  const Double_t intD = hData -> Integral();
  const Double_t intP = hSimP -> Integral();
  const Double_t intG = hSimG -> Integral();
  if (intD > 0.) hData -> Scale(1. / intD);
  if (intP > 0.) hSimP -> Scale(1. / intP);
  if (intG > 0.) hSimG -> Scale(1. / intG);

  TH1D *hSumN = (TH1D*) hSimP -> Clone();
  hSumN -> SetName("hSumNorm");
  hSumN -> Reset("ICES");
  hSumN -> Add(hSimP, hSimG);
  cout << "    Normalized histograms." << endl;

  // determine parameter maxima
  Double_t pFitMax[NPar];
  for (UInt_t iPar = 0; iPar < NPar; iPar++) {

    // get bin content of data at specified point
    Bool_t   isInFitRange = ((xFitMax[iPar] >= xFitRange[0]) && (xFitMax[iPar] <= xFitRange[1]));
    UInt_t   iFitMax(0);
    Double_t yFitMax(0.);
    if (isInFitRange) {
      iFitMax = hData -> FindBin(xFitMax[iPar]);
      yFitMax = hData -> GetBinContent(iFitMax);
    } else {
      yFitMax = hData -> GetMaximum();
    }

    // get maximum of simulation
    Double_t yParMax(1.);
    if (iPar == 0) {
      yParMax = hSimP -> GetMaximum();
    } else {
      yParMax = hSimG -> GetMaximum();
    }

    // set parameter maximum and announce
    pFitMax[iPar] = yFitMax / yParMax;
    if (iPar == 0) {
      cout << "    Determined fit parameter maxima:\n"
           << "      par[" << iPar << "]: max = " << pFitMax[iPar]
           << endl;
    } else {
      cout << "      par[" << iPar << "]: max = " << pFitMax[iPar] << endl;
    }
  }  // end particle loop

  // create histograms for comparison
  TH1D *hCompD = (TH1D*) hData -> Clone();
  TH1D *hCompP = (TH1D*) hSimP -> Clone();
  TH1D *hCompG = (TH1D*) hSimG -> Clone();
  TH1D *hCompS = (TH1D*) hSimP -> Clone();
  hCompS -> Reset("ICES");
  hCompS -> Add(hCompP, hCompG);
  hCompD -> SetName(sCompD.Data());
  hCompP -> SetName(sCompP.Data());
  hCompG -> SetName(sCompG.Data());
  hCompS -> SetName(sCompS.Data());
  cout << "    Created comparison histograms." << endl;

  // fit tsp
  TH1D     *hFitP;
  TH1D     *hFitG;
  TH1D     *hFitS;
  UInt_t   ndf(0);
  UInt_t   iterP(0);
  Bool_t   isNonzero(true);
  Bool_t   isInRange(true);
  Bool_t   isInParRange(true);
  Double_t pPi0(pFitMin[0]);
  Double_t pGam(pFitMin[1]);
  Double_t pPi0Min(1.);
  Double_t pGamMin(1.);
  Double_t chi2(0.);
  Double_t chi2min(99999.);
  Double_t binCenter(-1.);
  Double_t binValueD(0.);
  Double_t binValueP(0.);
  Double_t numer(0.);
  Double_t denom(0.);
  cout << "    Fitting TSP..." << endl;

  // fit parameter loop
  const UInt_t nBinD = hData -> GetNbinsX();
  const UInt_t nBinP = hSimP -> GetNbinsX();
  do {

    // increment pi0 parameter
    if (iterP != 0) {
      pPi0 += stepSize;
    }
    isInParRange = (pPi0 <= pFitMax[0]);
    iterP++;

    // scale histograms
    hFitP = (TH1D*) hCompP -> Clone();
    hFitP -> Scale(pPi0);

    // calculate chi2/ndf
    ndf  = 0;
    chi2 = 0.;
    for (UInt_t iBin = 1; iBin < (nBinD + 1); iBin++) {
      binCenter = hData -> GetBinCenter(iBin);
      binValueD = hData -> GetBinContent(iBin);
      binValueP = hFitP -> GetBinContent(iBin);
      numer     = TMath::Power(binValueD - binValueP, 2.);
      denom     = binValueP;
      isNonzero = ((binValueD > 0.) && (binValueP > 0.));
      isInRange = ((binCenter > xPi0Range[0]) && (binCenter < xPi0Range[1]));
      if (isNonzero && isInRange) {
        chi2 += numer / denom;
        ndf++;
      }
    }
    chi2 = chi2 / (ndf - 1);

    //chi2 = hCompD -> Chi2Test(hFitP, "CHI2/NDF");
    if (chi2 < chi2min) {
      chi2min = chi2;
      pPi0Min = pPi0;
    }
  } while (isInParRange);  // end fit parameter loop

  // determine gamma parameter
  hFitP = (TH1D*) hSimP  -> Clone();
  hFitG = (TH1D*) hSimG  -> Clone();
  hFitS = (TH1D*) hCompS -> Clone();
  hFitP -> Scale(pPi0Min);
  hFitS -> Reset("ICES");
  hFitS -> Add(hFitP, hFitG);

  const Double_t intPi0 = hFitP -> Integral();
  const Double_t intGam = hFitG -> Integral();
  const Double_t intSum = hFitS -> Integral();
  if (intGam > 0.) {
    pGamMin = (1. - intPi0) / intGam;
  }

  // scale histograms
  hFitP = (TH1D*) hSimP -> Clone();
  hFitG = (TH1D*) hSimG -> Clone();
  hSum  = (TH1D*) hSumN -> Clone();
  hFitP -> SetName(sFitP.Data());
  hFitG -> SetName(sFitG.Data());
  hSum  -> SetName(sFitS.Data());
  hFitP -> Scale(pPi0Min);
  hFitG -> Scale(pGamMin);
  hSum  -> Reset("ICES");
  hSum  -> Add(hFitP, hFitG);

  const UInt_t   iIntStart = hData -> FindBin(xFitRange[0]);
  const UInt_t   iIntStop  = hData -> FindBin(xFitRange[1]);
  const Double_t intFitPi0 = hFitP -> Integral(iIntStart, iIntStop);
  const Double_t intFitGam = hFitG -> Integral(iIntStart, iIntStop);
  cout << "    Fit data: chi2/ndf = " << chi2min << "\n"
       << "      pi0 = " << pPi0Min << ", gamma = " << pGamMin << "\n"
       << "      gamma/pi0 = " << intFitGam / intFitPi0
       << endl;

  // set styles
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fFil(0);
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const Float_t fSizT(0.04);
  const Float_t fSizL(0.04);
  const Float_t fOffX(1.5);
  const Float_t fOffY(1.5);
  hData  -> SetMarkerColor(fColD);
  hData  -> SetMarkerStyle(fMarD);
  hData  -> SetLineColor(fColD);
  hData  -> SetLineStyle(fLin);
  hData  -> SetFillColor(fColD);
  hData  -> SetFillStyle(fFil);
  hData  -> SetTitle(sTitle.Data());
  hData  -> SetTitleFont(fTxt);
  hData  -> GetXaxis() -> SetLabelFont(fTxt);
  hData  -> GetXaxis() -> SetLabelSize(fSizL);
  hData  -> GetXaxis() -> SetTitle(sXtitle.Data());
  hData  -> GetXaxis() -> SetTitleFont(fTxt);
  hData  -> GetXaxis() -> SetTitleSize(fSizT);
  hData  -> GetXaxis() -> SetTitleOffset(fOffX);
  hData  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hData  -> GetXaxis() -> CenterTitle(fCnt);
  hData  -> GetYaxis() -> SetLabelFont(fTxt);
  hData  -> GetYaxis() -> SetLabelSize(fSizL);
  hData  -> GetYaxis() -> SetTitle(sYtitle.Data());
  hData  -> GetYaxis() -> SetTitleFont(fTxt);
  hData  -> GetYaxis() -> SetTitleSize(fSizT);
  hData  -> GetYaxis() -> SetTitleOffset(fOffY);
  hData  -> GetYaxis() -> CenterTitle(fCnt);
  hSimP  -> SetMarkerColor(fColP);
  hSimP  -> SetMarkerStyle(fMarP);
  hSimP  -> SetLineColor(fColP);
  hSimP  -> SetLineStyle(fLin);
  hSimP  -> SetFillColor(fColP);
  hSimP  -> SetFillStyle(fFil);
  hSimP  -> SetTitle(sTitle.Data());
  hSimP  -> SetTitleFont(fTxt);
  hSimP  -> GetXaxis() -> SetLabelFont(fTxt);
  hSimP  -> GetXaxis() -> SetLabelSize(fSizL);
  hSimP  -> GetXaxis() -> SetTitle(sXtitle.Data());
  hSimP  -> GetXaxis() -> SetTitleFont(fTxt);
  hSimP  -> GetXaxis() -> SetTitleSize(fSizT);
  hSimP  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSimP  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hSimP  -> GetXaxis() -> CenterTitle(fCnt);
  hSimP  -> GetYaxis() -> SetLabelFont(fTxt);
  hSimP  -> GetYaxis() -> SetLabelSize(fSizL);
  hSimP  -> GetYaxis() -> SetTitle(sYtitle.Data());
  hSimP  -> GetYaxis() -> SetTitleFont(fTxt);
  hSimP  -> GetYaxis() -> SetTitleSize(fSizT);
  hSimP  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSimP  -> GetYaxis() -> CenterTitle(fCnt);
  hSimG  -> SetMarkerColor(fColG);
  hSimG  -> SetMarkerStyle(fMarG);
  hSimG  -> SetLineColor(fColG);
  hSimG  -> SetLineStyle(fLin);
  hSimG  -> SetFillColor(fColG);
  hSimG  -> SetFillStyle(fFil);
  hSimG  -> SetTitle(sTitle.Data());
  hSimG  -> SetTitleFont(fTxt);
  hSimG  -> GetXaxis() -> SetLabelFont(fTxt);
  hSimG  -> GetXaxis() -> SetLabelSize(fSizL);
  hSimG  -> GetXaxis() -> SetTitle(sXtitle.Data());
  hSimG  -> GetXaxis() -> SetTitleFont(fTxt);
  hSimG  -> GetXaxis() -> SetTitleSize(fSizT);
  hSimG  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSimG  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hSimG  -> GetXaxis() -> CenterTitle(fCnt);
  hSimG  -> GetYaxis() -> SetLabelFont(fTxt);
  hSimG  -> GetYaxis() -> SetLabelSize(fSizL);
  hSimG  -> GetYaxis() -> SetTitle(sYtitle.Data());
  hSimG  -> GetYaxis() -> SetTitleFont(fTxt);
  hSimG  -> GetYaxis() -> SetTitleSize(fSizT);
  hSimG  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSimG  -> GetYaxis() -> CenterTitle(fCnt);
  hCompD -> SetMarkerColor(fColD);
  hCompD -> SetMarkerStyle(fMarD);
  hCompD -> SetLineColor(fColD);
  hCompD -> SetLineStyle(fLin);
  hCompD -> SetFillColor(fColD);
  hCompD -> SetFillStyle(fFil);
  hCompD -> SetTitle(sTitle.Data());
  hCompD -> SetTitleFont(fTxt);
  hCompD -> GetXaxis() -> SetLabelFont(fTxt);
  hCompD -> GetXaxis() -> SetLabelSize(fSizL);
  hCompD -> GetXaxis() -> SetTitle(sXtitle.Data());
  hCompD -> GetXaxis() -> SetTitleFont(fTxt);
  hCompD -> GetXaxis() -> SetTitleSize(fSizT);
  hCompD -> GetXaxis() -> SetTitleOffset(fOffX);
  hCompD -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hCompD -> GetXaxis() -> CenterTitle(fCnt);
  hCompD -> GetYaxis() -> SetLabelFont(fTxt);
  hCompD -> GetYaxis() -> SetLabelSize(fSizL);
  hCompD -> GetYaxis() -> SetTitle(sYtitle.Data());
  hCompD -> GetYaxis() -> SetTitleFont(fTxt);
  hCompD -> GetYaxis() -> SetTitleSize(fSizT);
  hCompD -> GetYaxis() -> SetTitleOffset(fOffY);
  hCompD -> GetYaxis() -> CenterTitle(fCnt);
  hCompP -> SetMarkerColor(fColP);
  hCompP -> SetMarkerStyle(fMarP);
  hCompP -> SetLineColor(fColP);
  hCompP -> SetLineStyle(fLin);
  hCompP -> SetFillColor(fColP);
  hCompP -> SetFillStyle(fFil);
  hCompP -> SetTitle(sTitle.Data());
  hCompP -> SetTitleFont(fTxt);
  hCompP -> GetXaxis() -> SetLabelFont(fTxt);
  hCompP -> GetXaxis() -> SetLabelSize(fSizL);
  hCompP -> GetXaxis() -> SetTitle(sXtitle.Data());
  hCompP -> GetXaxis() -> SetTitleFont(fTxt);
  hCompP -> GetXaxis() -> SetTitleSize(fSizT);
  hCompP -> GetXaxis() -> SetTitleOffset(fOffX);
  hCompP -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hCompP -> GetXaxis() -> CenterTitle(fCnt);
  hCompP -> GetYaxis() -> SetLabelFont(fTxt);
  hCompP -> GetYaxis() -> SetLabelSize(fSizL);
  hCompP -> GetYaxis() -> SetTitle(sYtitle.Data());
  hCompP -> GetYaxis() -> SetTitleFont(fTxt);
  hCompP -> GetYaxis() -> SetTitleSize(fSizT);
  hCompP -> GetYaxis() -> SetTitleOffset(fOffY);
  hCompP -> GetYaxis() -> CenterTitle(fCnt);
  hCompG -> SetMarkerColor(fColG);
  hCompG -> SetMarkerStyle(fMarG);
  hCompG -> SetLineColor(fColG);
  hCompG -> SetLineStyle(fLin);
  hCompG -> SetFillColor(fColG);
  hCompG -> SetFillStyle(fFil);
  hCompG -> SetTitle(sTitle.Data());
  hCompG -> SetTitleFont(fTxt);
  hCompG -> GetXaxis() -> SetLabelFont(fTxt);
  hCompG -> GetXaxis() -> SetLabelSize(fSizL);
  hCompG -> GetXaxis() -> SetTitle(sXtitle.Data());
  hCompG -> GetXaxis() -> SetTitleFont(fTxt);
  hCompG -> GetXaxis() -> SetTitleSize(fSizT);
  hCompG -> GetXaxis() -> SetTitleOffset(fOffX);
  hCompG -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hCompG -> GetXaxis() -> CenterTitle(fCnt);
  hCompG -> GetYaxis() -> SetLabelFont(fTxt);
  hCompG -> GetYaxis() -> SetLabelSize(fSizL);
  hCompG -> GetYaxis() -> SetTitle(sYtitle.Data());
  hCompG -> GetYaxis() -> SetTitleFont(fTxt);
  hCompG -> GetYaxis() -> SetTitleSize(fSizT);
  hCompG -> GetYaxis() -> SetTitleOffset(fOffY);
  hCompG -> GetYaxis() -> CenterTitle(fCnt);
  hCompS -> SetMarkerColor(fColS);
  hCompS -> SetMarkerStyle(fMarS);
  hCompS -> SetLineColor(fColS);
  hCompS -> SetLineStyle(fLin);
  hCompS -> SetFillColor(fColS);
  hCompS -> SetFillStyle(fFil);
  hCompS -> SetTitle(sTitle.Data());
  hCompS -> SetTitleFont(fTxt);
  hCompS -> GetXaxis() -> SetLabelFont(fTxt);
  hCompS -> GetXaxis() -> SetLabelSize(fSizL);
  hCompS -> GetXaxis() -> SetTitle(sXtitle.Data());
  hCompS -> GetXaxis() -> SetTitleFont(fTxt);
  hCompS -> GetXaxis() -> SetTitleSize(fSizT);
  hCompS -> GetXaxis() -> SetTitleOffset(fOffX);
  hCompS -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hCompS -> GetXaxis() -> CenterTitle(fCnt);
  hCompS -> GetYaxis() -> SetLabelFont(fTxt);
  hCompS -> GetYaxis() -> SetLabelSize(fSizL);
  hCompS -> GetYaxis() -> SetTitle(sYtitle.Data());
  hCompS -> GetYaxis() -> SetTitleFont(fTxt);
  hCompS -> GetYaxis() -> SetTitleSize(fSizT);
  hCompS -> GetYaxis() -> SetTitleOffset(fOffY);
  hCompS -> GetYaxis() -> CenterTitle(fCnt);
  hFitP  -> SetMarkerColor(fColP);
  hFitP  -> SetMarkerStyle(fMarP);
  hFitP  -> SetLineColor(fColP);
  hFitP  -> SetLineStyle(fLin);
  hFitP  -> SetFillColor(fColP);
  hFitP  -> SetFillStyle(fFil);
  hFitP  -> SetTitle(sTitle.Data());
  hFitP  -> SetTitleFont(fTxt);
  hFitP  -> GetXaxis() -> SetLabelFont(fTxt);
  hFitP  -> GetXaxis() -> SetLabelSize(fSizL);
  hFitP  -> GetXaxis() -> SetTitle(sXtitle.Data());
  hFitP  -> GetXaxis() -> SetTitleFont(fTxt);
  hFitP  -> GetXaxis() -> SetTitleSize(fSizT);
  hFitP  -> GetXaxis() -> SetTitleOffset(fOffX);
  hFitP  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hFitP  -> GetXaxis() -> CenterTitle(fCnt);
  hFitP  -> GetYaxis() -> SetLabelFont(fTxt);
  hFitP  -> GetYaxis() -> SetLabelSize(fSizL);
  hFitP  -> GetYaxis() -> SetTitle(sYtitle.Data());
  hFitP  -> GetYaxis() -> SetTitleFont(fTxt);
  hFitP  -> GetYaxis() -> SetTitleSize(fSizT);
  hFitP  -> GetYaxis() -> SetTitleOffset(fOffY);
  hFitP  -> GetYaxis() -> CenterTitle(fCnt);
  hFitG  -> SetMarkerColor(fColG);
  hFitG  -> SetMarkerStyle(fMarG);
  hFitG  -> SetLineColor(fColG);
  hFitG  -> SetLineStyle(fLin);
  hFitG  -> SetFillColor(fColG);
  hFitG  -> SetFillStyle(fFil);
  hFitG  -> SetTitle(sTitle.Data());
  hFitG  -> SetTitleFont(fTxt);
  hFitG  -> GetXaxis() -> SetLabelFont(fTxt);
  hFitG  -> GetXaxis() -> SetLabelSize(fSizL);
  hFitG  -> GetXaxis() -> SetTitle(sXtitle.Data());
  hFitG  -> GetXaxis() -> SetTitleFont(fTxt);
  hFitG  -> GetXaxis() -> SetTitleSize(fSizT);
  hFitG  -> GetXaxis() -> SetTitleOffset(fOffX);
  hFitG  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hFitG  -> GetXaxis() -> CenterTitle(fCnt);
  hFitG  -> GetYaxis() -> SetLabelFont(fTxt);
  hFitG  -> GetYaxis() -> SetLabelSize(fSizL);
  hFitG  -> GetYaxis() -> SetTitle(sYtitle.Data());
  hFitG  -> GetYaxis() -> SetTitleFont(fTxt);
  hFitG  -> GetYaxis() -> SetTitleSize(fSizT);
  hFitG  -> GetYaxis() -> SetTitleOffset(fOffY);
  hFitG  -> GetYaxis() -> CenterTitle(fCnt);
  hSum   -> SetMarkerColor(fColS);
  hSum   -> SetMarkerStyle(fMarS);
  hSum   -> SetLineColor(fColS);
  hSum   -> SetLineStyle(fLin);
  hSum   -> SetFillColor(fColS);
  hSum   -> SetFillStyle(fFil);
  hSum   -> SetTitle(sTitle.Data());
  hSum   -> SetTitleFont(fTxt);
  hSum   -> GetXaxis() -> SetLabelFont(fTxt);
  hSum   -> GetXaxis() -> SetLabelSize(fSizL);
  hSum   -> GetXaxis() -> SetTitle(sXtitle.Data());
  hSum   -> GetXaxis() -> SetTitleFont(fTxt);
  hSum   -> GetXaxis() -> SetTitleSize(fSizT);
  hSum   -> GetXaxis() -> SetTitleOffset(fOffX);
  hSum   -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hSum   -> GetXaxis() -> CenterTitle(fCnt);
  hSum   -> GetYaxis() -> SetLabelFont(fTxt);
  hSum   -> GetYaxis() -> SetLabelSize(fSizL);
  hSum   -> GetYaxis() -> SetTitle(sYtitle.Data());
  hSum   -> GetYaxis() -> SetTitleFont(fTxt);
  hSum   -> GetYaxis() -> SetTitleSize(fSizT);
  hSum   -> GetYaxis() -> SetTitleOffset(fOffY);
  hSum   -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // create legend and text box
  const UInt_t  fColL(0);
  const UInt_t  fColT(0);
  const UInt_t  fLinL(0);
  const UInt_t  fLinT(0);
  const UInt_t  fFilL(0);
  const UInt_t  fFilT(0);
  const UInt_t  fAlnL(12);
  const UInt_t  fAlnT(32);
  const Float_t hObjL(0.08);
  const Float_t hTotL(hObjL * (NPar + 3));
  const Float_t yTotL(0.1 + hTotL);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, yTotL};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.2};

  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3], sTrg.Data());
  leg -> SetFillColor(fColL);
  leg -> SetFillStyle(fFilL);
  leg -> SetLineColor(fColL);
  leg -> SetLineStyle(fLinL);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnL);
  leg -> AddEntry(hData, sLabelD.Data(), "pf");
  leg -> AddEntry(hFitP, sLabelP.Data(), "pf");
  leg -> AddEntry(hFitG, sLabelG.Data(), "pf");
  leg -> AddEntry(hSum,  sLabelS.Data(), "pf");

  TPaveText *txt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  txt -> SetFillColor(fColT);
  txt -> SetFillStyle(fFilT);
  txt -> SetLineColor(fColT);
  txt -> SetLineStyle(fLinT);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlnT);
  cout << "    Made legend and text box." << endl;

  // create plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  grid(0);
  const UInt_t  tick(1);
  const UInt_t  logX(0);
  const UInt_t  logY(0);
  const UInt_t  mode(0);
  const UInt_t  frame(0);
  const UInt_t  border(2);
  const Float_t marginT(0.02);
  const Float_t marginR(0.02);
  const Float_t marginB(0.15);
  const Float_t marginL(0.15);

  TCanvas *cPlot = new TCanvas("cDataVsSim", "", width, height);
  cPlot -> SetLogx(logX);
  cPlot -> SetLogy(logY);
  cPlot -> SetGrid(grid, grid);
  cPlot -> SetTicks(tick, tick);
  cPlot -> SetBorderMode(mode);
  cPlot -> SetBorderSize(border);
  cPlot -> SetFrameBorderMode(frame);
  cPlot -> SetTopMargin(marginT);
  cPlot -> SetLeftMargin(marginL);
  cPlot -> SetRightMargin(marginR);
  cPlot -> SetBottomMargin(marginB);
  cPlot -> cd();
  hData -> Draw();
  hFitP -> Draw("hist same p");
  hFitG -> Draw("hist same p");
  hSum  -> Draw("hist same p");
  leg   -> Draw();
  fOut  -> cd();
  cPlot -> Write();
  cPlot -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut   -> cd();
  hData  -> Write();
  hSimP  -> Write();
  hSimG  -> Write();
  hSumN  -> Write();
  hCompD -> Write();
  hCompP -> Write();
  hCompG -> Write();
  hCompS -> Write();
  hFitP  -> Write();
  hFitG  -> Write();
  hSum   -> Write();
  fOut   -> Close();

  // close input files
  fData -> cd();
  fData -> Close();
  fSimP -> cd();
  fSimP -> Close();
  fSimG -> cd();
  fSimG -> Close();
  cout << "  TSP fit script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
