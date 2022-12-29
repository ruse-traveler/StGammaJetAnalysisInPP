// 'MakeParVsRecoEtPlot.C'
// Derek Anderson
// 11.02.2021
//
// Use this to plot the eTreco,
// back-smeared eTreco (fit and
// hist), and eTpar from the
// particle gun simulation.
//
//   TrgType = 0: pi0
//   TrgType = 1: gamma

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NTrg(3);
static const UInt_t NTyp(2);
static const UInt_t NVtx(4);
static const UInt_t NPlot(2);
static const UInt_t NHist(4);

// trigger type
static const UInt_t TrgType(0);



void MakeParVsRecoEtPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning eTreco vs. eTpar plot script..." << endl;

  // file parameters
  const TString sOut("particleGunRecoVsPar.et650x920tsp008pi0.d2m11y2021.root");
  const TString sPar("input/particleGun.full650pi0sample_defaultCuts.et650x650x1520vz55tsp008pi0.d24m2y2021.root");
  const TString sReco[NTrg]  = {"input/triggerWeights_et911pi0weights_full650pi0sample.et650x911vz55tsp008pi0.d24m3y2021.root",
                                "input/triggerWeights_et1115pi0weights_full650pi0sample.et650x1115vz55tsp008pi0.d24m3y2021.root",
                                "input/triggerWeights_et1520pi0weights_full650pi0sample.et650x1520vz55tsp008pi0.d24m3y2021.root"};
  const TString sMatch[NTrg] = {"input/triggerWeights_et911pi0weights_full650pi0sample.et650x911vz55tsp008pi0.d24m3y2021.root",
                                "input/triggerWeights_et1115pi0weights_full650pi0sample.et650x1115vz55tsp008pi0.d24m3y2021.root",
                                "input/triggerWeights_et1520pi0weights_full650pi0sample.et650x1520vz55tsp008pi0.d24m3y2021.root"};
  const TString sBackF[NTrg] = {"input/smearedEtMeas_usingFitForSmear.et911vz55tsp008pi0.d2m11y2021.root",
                                "input/smearedEtMeas_usingFitForSmear.et1115vz55tsp008pi0.d2m11y2021.root",
                                "input/smearedEtMeas_usingFitForSmear.et1520vz55tsp008pi0.d2m11y2021.root"};
  const TString sBackH[NTrg] = {"input/smearedEtMeas_usingHistForSmear.et911vz55tsp008pi0.d2m11y2021.root",
                                "input/smearedEtMeas_usingHistForSmear.et1115vz55tsp008pi0.d2m11y2021.root",
                                "input/smearedEtMeas_usingHistForSmear.et1520vz55tsp008pi0.d2m11y2021.root"};

  // input hist parameters
  const TString sHistP("particle/hEtTrgParAllWeight");
  const TString sHistR("hEtMeas");
  const TString sHistM("hEtPar");
  const TString sHistBF("hEtTrgSmear");
  const TString sHistBH("hEtTrgSmear");
  const TString sNameP("hParticleAll");
  const TString sNameR[NTrg]  = {"hReco911", "hReco1115", "hReco1520"};
  const TString sNameM[NTrg]  = {"hMatch911", "hMatch1115", "hMatch1520"};
  const TString sNameBF[NTrg] = {"hBackSmearFit911", "hBackSmearFit1115", "hBackSmearFit1520"};
  const TString sNameBH[NTrg] = {"hBackSmearDist911", "hBackSmearDist1520", "hBackSmearDist1520"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("E_{T}^{trg} [GeV]");
  const TString sTitleY("a. u.");
  const Float_t xPlotRange[NPlot]  = {0., 50.};
  const Float_t xNorm[NTrg][NPlot] = {{9., 11.}, {11., 15.}, {15., 20.}};
  const Float_t fMarSiz[NHist]     = {1.25, 1.25, 1.25, 1.25};
  const UInt_t  fMarR[NTrg]        = {20, 22, 21};
  const UInt_t  fMarM[NTrg]        = {24, 26, 25};
  const UInt_t  fMarBF[NTrg]       = {2, 2, 2};
  const UInt_t  fMarBH[NTrg]       = {5, 5, 5};

  // text parameters
  const TString sEne("E_{T}^{trg}(reco #otimes sim) = (9, 20) #otimes (6, 50) GeV");
  const TString sPrefR("E_{T}^{reco}[ ");
  const TString sPrefM("E_{T}^{match}[ ");
  const TString sPrefBF("#tilde{E}_{TF}^{reco}[ ");
  const TString sPrefBH("#tilde{E}_{TD}^{reco}[ ");
  const TString sBin[NTrg] = {"9-11 GeV", "11-15 GeV", "15-20 GeV"};

  // trigger parameters
  const TString sTypT[NTyp]         = {"#pi^{0}", "#gamma"};
  const TString sTrgT[NTyp]         = {"#bf{#color[859]{simulated #pi^{0}, weighting applied}}", "#bf{#color[899]{simulated #gamma, weighting applied}}"};
  const TString sCutT[NTyp]         = {"|#eta^{trg}| < 0.9, TSP #in (0, 0.08)", "|#eta^{trg}| < 0.9, TSP #in (0.2, 0.6)"};
  const TString sGridT[NTyp]        = {"grid spacing = 0.6 #pm 0.025", "grid spacing = 0.3 #pm 0.025"};
  const TString sHeadP[NHist]       = {"#bf{reconstructed #pi^{0}}", "#bf{matched simulated #pi^{0}}", "#bf{#pi^{0} back-smeared with fit}", "#bf{#pi^{0} back-smeared with distribution}"};
  const TString sHeadG[NHist]       = {"#bf{reconstructed #gamma}", "#bf{matched simulated #gamma}", "#bf{#gamma back-smeared with fit}", "#bf{#gamma back-smeared with distribution}"};
  const UInt_t  fColRT[NTyp][NTrg]  = {{603, 435, 419}, {635, 619, 603}};
  const UInt_t  fColMT[NTyp][NTrg]  = {{597, 429, 413}, {629, 613, 597}};
  const UInt_t  fColBFT[NTyp][NTrg] = {{854, 834, 814}, {894, 874, 854}};
  const UInt_t  fColBHT[NTyp][NTrg] = {{864, 844, 824}, {904, 884, 864}};

  // parse trigger selection
  TString sTyp;
  TString sTrg;
  TString sCut;
  TString sGrid;
  TString sHead[NHist];
  UInt_t  fColR[NTrg];
  UInt_t  fColM[NTrg];
  UInt_t  fColBF[NTrg];
  UInt_t  fColBH[NTrg];
  switch (TrgType) {
    case 0:
      sTyp  = sTypT[0];
      sTrg  = sTrgT[0];
      sCut  = sCutT[0];
      sGrid = sGridT[0];
      for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
        fColR[iTrg]  = fColRT[0][iTrg];
        fColM[iTrg]  = fColMT[0][iTrg];
        fColBF[iTrg] = fColBFT[0][iTrg];
        fColBH[iTrg] = fColBHT[0][iTrg];
      }
      for (UInt_t iHist = 0; iHist < NHist; iHist++) {
        sHead[iHist] = sHeadP[iHist];
      }
      break;
    case 1:
      sTyp  = sTypT[1];
      sTrg  = sTrgT[1];
      sCut  = sCutT[1];
      sGrid = sGridT[1];
      for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
        fColR[iTrg]  = fColRT[1][iTrg];
        fColM[iTrg]  = fColMT[1][iTrg];
        fColBF[iTrg] = fColBFT[1][iTrg];
        fColBH[iTrg] = fColBHT[1][iTrg];
      }
      for (UInt_t iHist = 0; iHist < NHist; iHist++) {
        sHead[iHist] = sHeadG[iHist];
      }
      break;
    default:
      sTyp  = sTypT[0];
      sTrg  = sTrgT[0];
      sCut  = sCutT[0];
      sGrid = sGridT[0];
      for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
        fColR[iTrg]  = fColRT[0][iTrg];
        fColM[iTrg]  = fColMT[0][iTrg];
        fColBF[iTrg] = fColBFT[0][iTrg];
        fColBH[iTrg] = fColBHT[0][iTrg];
      }
      for (UInt_t iHist = 0; iHist < NHist; iHist++) {
        sHead[iHist] = sHeadP[iHist];
      }
      break;
  }
  cout << "    Parsed trigger selection." << endl;

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fPar = new TFile(sPar.Data(), "read");
  if (!fOut || !fPar) {
    cerr << "PANIC: couldn't open output or particle file!\n"
         << "       fOut = " << fOut << ", fPar = " << fPar << "\n"
         << endl;
    return;
  }

  TFile *fReco[NTrg];
  TFile *fMatch[NTrg];
  TFile *fBackF[NTrg];
  TFile *fBackH[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fReco[iTrg]  = new TFile(sReco[iTrg].Data(), "read");
    fMatch[iTrg] = new TFile(sMatch[iTrg].Data(), "read");
    fBackF[iTrg] = new TFile(sBackF[iTrg].Data(), "read");
    fBackH[iTrg] = new TFile(sBackH[iTrg].Data(), "read");
    if (!fReco[iTrg] || !fMatch[iTrg] || !fBackF[iTrg] || !fBackH[iTrg]) {
      cerr << "PANIC: couldn't open a trigger file (#" << iTrg <<")!\n"
           << "       fReco  = " << fReco[iTrg]  << ", fMatch = " << fMatch[iTrg] << "\n"
           << "       fBackF = " << fBackF[iTrg] << ", fBackH = " << fBackH[iTrg] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hPar = (TH1D*) fPar -> Get(sHistP.Data());
  if (!hPar) {
    cerr << "PANIC: couldn't grab particle histogram!\n" << endl;
    return;
  }
  hPar -> SetName(sNameP.Data());

  TH1D *hReco[NTrg];
  TH1D *hMatch[NTrg];
  TH1D *hBackF[NTrg];
  TH1D *hBackH[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hReco[iTrg]  = (TH1D*) fReco[iTrg]  -> Get(sHistR.Data());
    hMatch[iTrg] = (TH1D*) fMatch[iTrg] -> Get(sHistM.Data());
    hBackF[iTrg] = (TH1D*) fBackF[iTrg] -> Get(sHistBF.Data());
    hBackH[iTrg] = (TH1D*) fBackH[iTrg] -> Get(sHistBH.Data());
    if (!hReco[iTrg] || !hMatch[iTrg] || !hBackF[iTrg] || !hBackH[iTrg]) {
      cerr << "PANIC: couldn't grab trigger histogram (#" << iTrg << ")!\n"
           << "       hReco  = " << hReco[iTrg]  << ", hMatch = " << hMatch[iTrg] << "\n"
           << "       hBackF = " << hBackF[iTrg] << ", hBackH = " << hBackH[iTrg] << "\n"
           << endl;
      return;
    }
    hReco[iTrg]  -> SetName(sNameR[iTrg].Data());
    hMatch[iTrg] -> SetName(sNameM[iTrg].Data());
    hBackF[iTrg] -> SetName(sNameBF[iTrg].Data());
    hBackH[iTrg] -> SetName(sNameBH[iTrg].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // normalize histograms
  const Double_t binShift = 0.005;
  const Double_t intPar   = hPar -> Integral();
  const Bool_t   isFinP   = (intPar > 0.);
  if (isFinP) hPar -> Scale(1. / intPar);

  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    const UInt_t   iStart  = hPar         -> FindBin(xNorm[iTrg][0] + binShift);
    const UInt_t   iStop   = hPar         -> FindBin(xNorm[iTrg][1] - binShift);
    const Double_t intNorm = hPar         -> Integral(iStart, iStop);
    const Double_t intRec  = hReco[iTrg]  -> Integral();
    const Double_t intMat  = hMatch[iTrg] -> Integral();
    const Double_t intBaF  = hBackF[iTrg] -> Integral();
    const Double_t intBaH  = hBackH[iTrg] -> Integral();
    const Bool_t   isFinR  = (intRec > 0.);
    const Bool_t   isFinM  = (intMat > 0.);
    const Bool_t   isFinBF = (intBaF > 0.);
    const Bool_t   isFinBH = (intBaH > 0.);
    if (isFinR)  hReco[iTrg]  -> Scale(intNorm / intRec);
    if (isFinM)  hMatch[iTrg] -> Scale(intNorm / intMat);
    if (isFinBF) hBackF[iTrg] -> Scale(intNorm / intBaF);
    if (isFinBH) hBackH[iTrg] -> Scale(intNorm / intBaH);
  }
  cout << "    Normalized histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fSiz(1.25);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hReco[iTrg]  -> SetMarkerColor(fColR[iTrg]);
    hReco[iTrg]  -> SetMarkerStyle(fMarR[iTrg]);
    hReco[iTrg]  -> SetMarkerSize(fMarSiz[0]);
    hReco[iTrg]  -> SetFillColor(fColR[iTrg]);
    hReco[iTrg]  -> SetFillStyle(fFil);
    hReco[iTrg]  -> SetLineColor(fColR[iTrg]);
    hReco[iTrg]  -> SetLineStyle(fLin);
    hReco[iTrg]  -> SetLineWidth(fWid);
    hReco[iTrg]  -> SetTitle(sTitle.Data());
    hReco[iTrg]  -> SetTitleFont(fTxt);
    hReco[iTrg]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hReco[iTrg]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hReco[iTrg]  -> GetXaxis() -> SetTitleFont(fTxt);
    hReco[iTrg]  -> GetXaxis() -> SetTitleSize(fTit);
    hReco[iTrg]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hReco[iTrg]  -> GetXaxis() -> SetLabelFont(fTxt);
    hReco[iTrg]  -> GetXaxis() -> SetLabelSize(fLab);
    hReco[iTrg]  -> GetXaxis() -> CenterTitle(fCnt);
    hReco[iTrg]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    hReco[iTrg]  -> GetYaxis() -> SetTitleFont(fTxt);
    hReco[iTrg]  -> GetYaxis() -> SetTitleSize(fTit);
    hReco[iTrg]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hReco[iTrg]  -> GetYaxis() -> SetLabelFont(fTxt);
    hReco[iTrg]  -> GetYaxis() -> SetLabelSize(fLab);
    hReco[iTrg]  -> GetYaxis() -> CenterTitle(fCnt);
    hMatch[iTrg] -> SetMarkerColor(fColM[iTrg]);
    hMatch[iTrg] -> SetMarkerStyle(fMarM[iTrg]);
    hMatch[iTrg] -> SetMarkerSize(fMarSiz[1]);
    hMatch[iTrg] -> SetFillColor(fColM[iTrg]);
    hMatch[iTrg] -> SetFillStyle(fFil);
    hMatch[iTrg] -> SetLineColor(fColM[iTrg]);
    hMatch[iTrg] -> SetLineStyle(fLin);
    hMatch[iTrg] -> SetLineWidth(fWid);
    hMatch[iTrg] -> SetTitle(sTitle.Data());
    hMatch[iTrg] -> SetTitleFont(fTxt);
    hMatch[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hMatch[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hMatch[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hMatch[iTrg] -> GetXaxis() -> SetTitleSize(fTit);
    hMatch[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hMatch[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hMatch[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hMatch[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hMatch[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hMatch[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hMatch[iTrg] -> GetYaxis() -> SetTitleSize(fTit);
    hMatch[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hMatch[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hMatch[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hMatch[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hBackF[iTrg] -> SetMarkerColor(fColBF[iTrg]);
    hBackF[iTrg] -> SetMarkerStyle(fMarBF[iTrg]);
    hBackF[iTrg] -> SetMarkerSize(fMarSiz[2]);
    hBackF[iTrg] -> SetFillColor(fColBF[iTrg]);
    hBackF[iTrg] -> SetFillStyle(fFil);
    hBackF[iTrg] -> SetLineColor(fColBF[iTrg]);
    hBackF[iTrg] -> SetLineStyle(fLin);
    hBackF[iTrg] -> SetLineWidth(fWid);
    hBackF[iTrg] -> SetTitle(sTitle.Data());
    hBackF[iTrg] -> SetTitleFont(fTxt);
    hBackF[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hBackF[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hBackF[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hBackF[iTrg] -> GetXaxis() -> SetTitleSize(fTit);
    hBackF[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hBackF[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hBackF[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hBackF[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hBackF[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hBackF[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hBackF[iTrg] -> GetYaxis() -> SetTitleSize(fTit);
    hBackF[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hBackF[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hBackF[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hBackF[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hBackH[iTrg] -> SetMarkerColor(fColBH[iTrg]);
    hBackH[iTrg] -> SetMarkerStyle(fMarBH[iTrg]);
    hBackH[iTrg] -> SetMarkerSize(fMarSiz[3]);
    hBackH[iTrg] -> SetFillColor(fColBH[iTrg]);
    hBackH[iTrg] -> SetFillStyle(fFil);
    hBackH[iTrg] -> SetLineColor(fColBH[iTrg]);
    hBackH[iTrg] -> SetLineStyle(fLin);
    hBackH[iTrg] -> SetLineWidth(fWid);
    hBackH[iTrg] -> SetTitle(sTitle.Data());
    hBackH[iTrg] -> SetTitleFont(fTxt);
    hBackH[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hBackH[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hBackH[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hBackH[iTrg] -> GetXaxis() -> SetTitleSize(fTit);
    hBackH[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hBackH[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hBackH[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hBackH[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hBackH[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hBackH[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hBackH[iTrg] -> GetYaxis() -> SetTitleSize(fTit);
    hBackH[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hBackH[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hBackH[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hBackH[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make labels
  TString sLabelR[NTrg];
  TString sLabelM[NTrg];
  TString sLabelBF[NTrg];
  TString sLabelBH[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    sLabelR[iTrg]  = sPrefR;
    sLabelM[iTrg]  = sPrefM;
    sLabelBF[iTrg] = sPrefBF;
    sLabelBH[iTrg] = sPrefBH;
    sLabelR[iTrg].Append(sBin[iTrg].Data());
    sLabelR[iTrg].Append(" | ");
    sLabelR[iTrg].Append(sTyp.Data());
    sLabelR[iTrg].Append(" ]");
    sLabelM[iTrg].Append(sBin[iTrg].Data());
    sLabelM[iTrg].Append(" | ");
    sLabelM[iTrg].Append(sTyp.Data());
    sLabelM[iTrg].Append(" ]");
    sLabelBF[iTrg].Append(sBin[iTrg].Data());
    sLabelBF[iTrg].Append(" | ");
    sLabelBF[iTrg].Append(sTyp.Data());
    sLabelBF[iTrg].Append(" ]");
    sLabelBH[iTrg].Append(sBin[iTrg].Data());
    sLabelBH[iTrg].Append(" | ");
    sLabelBH[iTrg].Append(sTyp.Data());
    sLabelBH[iTrg].Append(" ]");
  }

  // make legend
  const UInt_t  fNumLe(2);
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  nObjLe((NHist / 2) + (2 * NTrg));
  const Float_t hObjLe(0.05);
  const Float_t yTotLe((nObjLe * hObjLe) + 0.1);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.5, yTotLe};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> SetNColumns(fNumLe);
  leg -> AddEntry((TObject*)0, sHead[0].Data(), "");
  leg -> AddEntry((TObject*)0, sHead[1].Data(), "");
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    leg -> AddEntry(hReco[iTrg], sLabelR[iTrg], "pf");
    leg -> AddEntry(hMatch[iTrg], sLabelM[iTrg], "pf");
  }
  leg -> AddEntry((TObject*)0, sHead[2].Data(), "");
  leg -> AddEntry((TObject*)0, sHead[3].Data(), "");
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    leg -> AddEntry(hBackF[iTrg], sLabelBF[iTrg], "pf");
    leg -> AddEntry(hBackH[iTrg], sLabelBH[iTrg], "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const Float_t fTxtXY[NVtx] = {0.5, 0.1, 0.9, 0.3};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sTrg.Data());
  txt -> AddText(sEne.Data());
  txt -> AddText(sCut.Data());
  txt -> AddText(sGrid.Data());
  cout << "    Made text." << endl;

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
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.15);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot    -> SetGrid(fGrid, fGrid);
  cPlot    -> SetTicks(fTick, fTick);
  cPlot    -> SetLogx(fLogX);
  cPlot    -> SetLogy(fLogY);
  cPlot    -> SetBorderMode(fMode);
  cPlot    -> SetBorderSize(fBord);
  cPlot    -> SetFrameBorderMode(fFrame);
  cPlot    -> SetLeftMargin(fMarginL);
  cPlot    -> SetRightMargin(fMarginR);
  cPlot    -> SetTopMargin(fMarginT);
  cPlot    -> SetBottomMargin(fMarginB);
  cPlot    -> cd();
  hReco[0] -> Draw("hist p");
  for (UInt_t iTrg = 1; iTrg < NTrg; iTrg++) {
    hReco[iTrg] -> Draw("hist p same");
  }
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hMatch[iTrg] -> Draw("hist p same");
  }
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hBackF[iTrg] -> Draw("hist p same");
  }
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hBackH[iTrg] -> Draw("hist p same");
  }
  leg   -> Draw();
  txt   -> Draw();
  fOut  -> cd();
  cPlot -> Write();
  cPlot -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut ->  cd();
  hPar -> Write();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hReco[iTrg]  -> Write();
    hMatch[iTrg] -> Write();
    hBackF[iTrg] -> Write();
    hBackH[iTrg] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fPar -> cd();
  fPar -> Close();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fReco[iTrg]  -> cd();
    fReco[iTrg]  -> Close();
    fMatch[iTrg] -> cd();
    fMatch[iTrg] -> Close();
    fBackF[iTrg] -> cd();
    fBackF[iTrg] -> Close();
    fBackH[iTrg] -> cd();
    fBackH[iTrg] -> Close();
  }
  cout << "  Finished eTreco vs. eTpar plot script!\n" << endl;

}

// End ------------------------------------------------------------------------
