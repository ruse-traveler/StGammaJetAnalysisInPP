// 'MakeUnfoldedPlot.C'
// Derek Anderson
// 07.16.2018
//
// Use this to plot an unfolded distribution
// vs. pythia (or otherwise) distribution(s).
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

using namespace std;

// global constants
static const UInt_t NPad(2);
static const UInt_t NPlot(2);
static const UInt_t NData(2);
static const UInt_t NPythia(1);
static const UInt_t NTrgIds(2);
static const UInt_t NTrgBins(3);
static const UInt_t NTrgs(NTrgIds * NTrgBins);

// trigger parameters
static const UInt_t TrigId(0);
static const UInt_t TrigBin(0);



void MakeUnfoldedPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distribution..." << endl;

  // data io parameters
  const TString sOut("pp200r9unfold.noSmoothOrSecondLevy_pTbinFine.et911r05pi0.d27m8y2021.root");
  const TString sInTot("summedErrors_noSmoothOrSecondLevy_pTbinFine.et911r05pi0.d25m8y2021.root");
  const TString sInStat("summedErrors_noSmoothOrSecondLevy_pTbinFine.et911r05pi0.d25m8y2021.root");
  const TString sHistTot("hTotal");
  const TString sHistStat("hStatistics");

  // io pythia parameters
  const TString sInPy[NPythia]   = {"input/pythia/pp200py8par.forComparison_pTbinFine.et911nTrg100Kpt021Kpi0.r05a065rm1chrg.root"};
  const TString sHistPy[NPythia] = {"Pi0/hJetPtCorrP"};

  // data plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T}^{unfold} [GeV/c]");
  const TString  sTitleY("(1/N^{trg}) d^{2}N^{jet}/d(p_{T}^{unfold} #eta^{jet}) [GeV/c]^{-1}");
  const TString  sNameTot("hUnfoldTotal");
  const TString  sNameStat("hUnfoldStat");
  const TString  sLabelTot("unfolded data [stat. #oplus sys.]");
  const TString  sLabelStat("unfolded data [stat.]");
  const Double_t xPlotRange[NPlot] = {-1., 30.};

  // pythia and ratio plot parameters
  const TString sTitleRY("unfolded / par.-lvl.");
  const TString sNamePy[NPythia]     = {"hPy8noWeight"};
  const TString sNameTotR[NPythia]   = {"hRatioTotal"};
  const TString sNameStatR[NPythia]  = {"hRatioStat"};
  const TString sLabelPy[NPythia]    = {"unweighted Pythia8"};
  const TString sLabelTotR[NPythia]  = {"stat. #oplus sys."};
  const TString sLabelStatR[NPythia] = {"stat."};
  const UInt_t  fLinWid[NPythia]     = {2};
  const UInt_t  fLinSty[NPythia]     = {1};
  const UInt_t  fFilStyR[NPythia]    = {3001};

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // define color schemes
  const UInt_t fColDatT[NTrgs]            = {590, 422, 406, 622, 606, 590};
  const UInt_t fColDatS[NTrgs]            = {859, 839, 819, 899, 879, 859};
  const UInt_t fColPyS[NTrgs][NPythia]    = {603, 435, 419, 635, 619, 603};
  const UInt_t fColRatioT[NTrgs][NPythia] = {592, 424, 408, 624, 608, 592};
  const UInt_t fColRatioS[NTrgs][NPythia] = {603, 435, 419, 635, 619, 603};

  // parse trigger selection
  UInt_t fTrg(0);
  if (TrigId == 0) {
    fTrg = TrigBin;
  } else {
    fTrg = TrigBin + 3;
  }

  // select trigger
  UInt_t  fColTot(0);
  UInt_t  fColStat(0);
  UInt_t  fColPy[NPythia];
  UInt_t  fColTotR[NPythia];
  UInt_t  fColStatR[NPythia];
  TString sTrg("");
  switch (fTrg) {
    case 0:
      sTrg         = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      fColTot      = fColDatT[0];
      fColStat     = fColDatS[0];
      fColPy[0]    = fColPyS[0][0];
      fColPy[1]    = fColPyS[0][1];
      fColTotR[0]  = fColRatioT[0][0];
      fColTotR[1]  = fColRatioT[0][1];
      fColStatR[0] = fColRatioS[0][0];
      fColStatR[1] = fColRatioS[0][1];
      break;
    case 1:
      sTrg         = "#pi^{0} trig., E_{T}^{meas} #in (11, 15) GeV";
      fColTot      = fColDatT[1];
      fColStat     = fColDatS[1];
      fColPy[0]    = fColPyS[1][0];
      fColPy[1]    = fColPyS[1][1];
      fColTotR[0]  = fColRatioT[1][0];
      fColTotR[1]  = fColRatioT[1][1];
      fColStatR[0] = fColRatioS[1][0];
      fColStatR[1] = fColRatioS[1][1];
      break;
    case 2:
      sTrg         = "#pi^{0} trig., E_{T}^{meas} #in (15, 20) GeV";
      fColTot      = fColDatT[2];
      fColStat     = fColDatS[2];
      fColPy[0]    = fColPyS[2][0];
      fColPy[1]    = fColPyS[2][1];
      fColTotR[0]  = fColRatioT[2][0];
      fColTotR[1]  = fColRatioT[2][1];
      fColStatR[0] = fColRatioS[2][0];
      fColStatR[1] = fColRatioS[2][1];
      break;
    case 3:
      sTrg         = "#gamma^{dir} trig., E_{T}^{meas} #in (9, 11) GeV";
      fColTot      = fColDatT[3];
      fColStat     = fColDatS[3];
      fColPy[0]    = fColPyS[3][0];
      fColPy[1]    = fColPyS[3][1];
      fColTotR[0]  = fColRatioT[3][0];
      fColTotR[1]  = fColRatioT[3][1];
      fColStatR[0] = fColRatioS[3][0];
      fColStatR[1] = fColRatioS[3][1];
      break;
    case 4:
      sTrg         = "#gamma^{dir} trig., E_{T}^{meas} #in (11, 15) GeV";
      fColTot      = fColDatT[4];
      fColStat     = fColDatS[4];
      fColPy[0]    = fColPyS[4][0];
      fColPy[1]    = fColPyS[4][1];
      fColTotR[0]  = fColRatioT[4][0];
      fColTotR[1]  = fColRatioT[4][1];
      fColStatR[0] = fColRatioS[4][0];
      fColStatR[1] = fColRatioS[4][1];
      break;
    case 5:
      sTrg         = "#gamma^{dir} trig., E_{T}^{meas} #in (15, 20) GeV";
      fColTot      = fColDatT[5];
      fColStat     = fColDatS[5];
      fColPy[0]    = fColPyS[5][0];
      fColPy[1]    = fColPyS[5][1];
      fColTotR[0]  = fColRatioT[5][0];
      fColTotR[1]  = fColRatioT[5][1];
      fColStatR[0] = fColRatioS[5][0];
      fColStatR[1] = fColRatioS[5][1];
      break;
    default:
      sTrg         = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      fColTot      = fColDatT[0];
      fColStat     = fColDatS[0];
      fColPy[0]    = fColPyS[0][0];
      fColPy[1]    = fColPyS[0][1];
      fColTotR[0]  = fColRatioT[0][0];
      fColTotR[1]  = fColRatioT[0][1];
      fColStatR[0] = fColRatioS[0][0];
      fColStatR[1] = fColRatioS[0][1];
      break;
  }  // end switch-case

  // open data files
  TFile *fOut    = new TFile(sOut.Data(), "recreate");
  TFile *fInTot  = new TFile(sInTot.Data(), "read");
  TFile *fInStat = new TFile(sInStat.Data(), "read");
  if (!fOut || !fInTot || !fInStat) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInTot = " << fInTot << ", fInStat" << fInStat
         << endl;
    return;
  }
  cout << "    Opened data files." << endl;

  // open pythia files
  TFile *fInPy[NPythia];
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {

    fInPy[iPythia] = new TFile(sInPy[iPythia], "read");
    if (!fInPy[iPythia]) {
      cerr << "PANIC: couldn't open a pythia file!\n"
           << "       fInPy[" << iPythia << "] = " << fInPy[iPythia]
           << endl;
      return;
    }

  }  // end pythia loop
  cout << "    Opened pythia files." << endl;

  // grab data histograms
  TH1D *hUnfoldTot  = (TH1D*) fInTot  -> Get(sHistTot.Data());
  TH1D *hUnfoldStat = (TH1D*) fInStat -> Get(sHistStat.Data());
  if (!hUnfoldTot || !hUnfoldStat) {
    cerr << "PANIC: couldn't grab a data histogram!\n"
         << "       hUnfoldTot = " << hUnfoldTot << ", hUnfoldStat = " << hUnfoldStat
         << endl;
    return;
  }
  hUnfoldTot  -> SetName(sNameTot.Data());
  hUnfoldStat -> SetName(sNameStat.Data());
  cout << "    Grabbed data histograms." << endl;

  // grab pythia histograms
  TH1D *hPythia[NPythia];
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {

    hPythia[iPythia] = (TH1D*) fInPy[iPythia] -> Get(sHistPy[iPythia].Data());
    if (!hPythia[iPythia]) {
      cerr << "PANIC: couldn't grab a pythia histogram!\n"
           << "        hPythia[" << iPythia << "] = " << hPythia[iPythia]
           << endl;
      return;
    }
    hPythia[iPythia] -> SetName(sNamePy[iPythia].Data());

  }  // end pythia loop
  cout << "    Grabbed pythia histograms." << endl;

  // calculate ratio(s)
  TH1D *hRatioTot[NPythia];
  TH1D *hRatioStat[NPythia];
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    hRatioTot[iPythia]  = (TH1D*) hUnfoldTot  -> Clone();
    hRatioStat[iPythia] = (TH1D*) hUnfoldStat -> Clone();
    hRatioTot[iPythia]  -> Divide(hUnfoldTot, hPythia[iPythia], 1., 1.);
    hRatioStat[iPythia] -> Divide(hUnfoldStat, hPythia[iPythia], 1., 1.);
    hRatioTot[iPythia]  -> SetName(sNameTotR[iPythia].Data());
    hRatioStat[iPythia] -> SetName(sNameStatR[iPythia].Data());
  }
  cout << "    Calculated ratio(s)." << endl;

  // set styles
  const UInt_t  fMarDat(29);
  const UInt_t  fMarPy(1);
  const UInt_t  fFilDat(3001);
  const UInt_t  fFilPy(0);
  const UInt_t  fLinDat(1);
  const UInt_t  fWidDat(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hUnfoldTot  -> SetMarkerColor(fColTot);
  hUnfoldTot  -> SetMarkerStyle(fMarDat);
  hUnfoldTot  -> SetFillColor(fColTot);
  hUnfoldTot  -> SetFillStyle(fFilDat);
  hUnfoldTot  -> SetLineColor(fColTot);
  hUnfoldTot  -> SetLineStyle(fLinDat);
  hUnfoldTot  -> SetLineWidth(fWidDat);
  hUnfoldTot  -> SetTitle(sTitle.Data());
  hUnfoldTot  -> SetTitleFont(fTxt);
  hUnfoldTot  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hUnfoldTot  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hUnfoldTot  -> GetXaxis() -> SetTitleFont(fTxt);
  hUnfoldTot  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hUnfoldTot  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hUnfoldTot  -> GetXaxis() -> SetLabelFont(fTxt);
  hUnfoldTot  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hUnfoldTot  -> GetXaxis() -> CenterTitle(fCnt);
  hUnfoldTot  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hUnfoldTot  -> GetYaxis() -> SetTitleFont(fTxt);
  hUnfoldTot  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hUnfoldTot  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hUnfoldTot  -> GetYaxis() -> SetLabelFont(fTxt);
  hUnfoldTot  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hUnfoldTot  -> GetYaxis() -> CenterTitle(fCnt);
  hUnfoldStat -> SetMarkerColor(fColStat);
  hUnfoldStat -> SetMarkerStyle(fMarDat);
  hUnfoldStat -> SetFillColor(fColStat);
  hUnfoldStat -> SetFillStyle(fFilDat);
  hUnfoldStat -> SetLineColor(fColStat);
  hUnfoldStat -> SetLineStyle(fLinDat);
  hUnfoldStat -> SetLineWidth(fWidDat);
  hUnfoldStat -> SetTitle(sTitle.Data());
  hUnfoldStat -> SetTitleFont(fTxt);
  hUnfoldStat -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hUnfoldStat -> GetXaxis() -> SetTitle(sTitleX.Data());
  hUnfoldStat -> GetXaxis() -> SetTitleFont(fTxt);
  hUnfoldStat -> GetXaxis() -> SetTitleSize(fTit[1]);
  hUnfoldStat -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hUnfoldStat -> GetXaxis() -> SetLabelFont(fTxt);
  hUnfoldStat -> GetXaxis() -> SetLabelSize(fLab[1]);
  hUnfoldStat -> GetXaxis() -> CenterTitle(fCnt);
  hUnfoldStat -> GetYaxis() -> SetTitle(sTitleY.Data());
  hUnfoldStat -> GetYaxis() -> SetTitleFont(fTxt);
  hUnfoldStat -> GetYaxis() -> SetTitleSize(fTit[1]);
  hUnfoldStat -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hUnfoldStat -> GetYaxis() -> SetLabelFont(fTxt);
  hUnfoldStat -> GetYaxis() -> SetLabelSize(fLab[1]);
  hUnfoldStat -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    hPythia[iPythia]    -> SetMarkerColor(fColPy[iPythia]);
    hPythia[iPythia]    -> SetMarkerStyle(fMarPy);
    hPythia[iPythia]    -> SetFillColor(fColPy[iPythia]);
    hPythia[iPythia]    -> SetFillStyle(fFilPy);
    hPythia[iPythia]    -> SetLineColor(fColPy[iPythia]);
    hPythia[iPythia]    -> SetLineStyle(fLinSty[iPythia]);
    hPythia[iPythia]    -> SetLineWidth(fLinWid[iPythia]);
    hPythia[iPythia]    -> SetTitle(sTitle.Data());
    hPythia[iPythia]    -> SetTitleFont(fTxt);
    hPythia[iPythia]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPythia[iPythia]    -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPythia[iPythia]    -> GetXaxis() -> SetTitleFont(fTxt);
    hPythia[iPythia]    -> GetXaxis() -> SetTitleSize(fTit[1]);
    hPythia[iPythia]    -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hPythia[iPythia]    -> GetXaxis() -> SetLabelFont(fTxt);
    hPythia[iPythia]    -> GetXaxis() -> SetLabelSize(fLab[1]);
    hPythia[iPythia]    -> GetXaxis() -> CenterTitle(fCnt);
    hPythia[iPythia]    -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPythia[iPythia]    -> GetYaxis() -> SetTitleFont(fTxt);
    hPythia[iPythia]    -> GetYaxis() -> SetTitleSize(fTit[1]);
    hPythia[iPythia]    -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hPythia[iPythia]    -> GetYaxis() -> SetLabelFont(fTxt);
    hPythia[iPythia]    -> GetYaxis() -> SetLabelSize(fLab[1]);
    hPythia[iPythia]    -> GetYaxis() -> CenterTitle(fCnt);
    hRatioTot[iPythia]  -> SetMarkerColor(fColTotR[iPythia]);
    hRatioTot[iPythia]  -> SetMarkerStyle(fMarDat);
    hRatioTot[iPythia]  -> SetFillColor(fColTotR[iPythia]);
    hRatioTot[iPythia]  -> SetFillStyle(fFilStyR[iPythia]);
    hRatioTot[iPythia]  -> SetLineColor(fColTotR[iPythia]);
    hRatioTot[iPythia]  -> SetLineStyle(fLinDat);
    hRatioTot[iPythia]  -> SetLineWidth(fWidDat);
    hRatioTot[iPythia]  -> SetTitle(sTitle.Data());
    hRatioTot[iPythia]  -> SetTitleFont(fTxt);
    hRatioTot[iPythia]  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatioTot[iPythia]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatioTot[iPythia]  -> GetXaxis() -> SetTitleFont(fTxt);
    hRatioTot[iPythia]  -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatioTot[iPythia]  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatioTot[iPythia]  -> GetXaxis() -> SetLabelFont(fTxt);
    hRatioTot[iPythia]  -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatioTot[iPythia]  -> GetXaxis() -> CenterTitle(fCnt);
    hRatioTot[iPythia]  -> GetYaxis() -> SetTitle(sTitleRY.Data());
    hRatioTot[iPythia]  -> GetYaxis() -> SetTitleFont(fTxt);
    hRatioTot[iPythia]  -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatioTot[iPythia]  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatioTot[iPythia]  -> GetYaxis() -> SetLabelFont(fTxt);
    hRatioTot[iPythia]  -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatioTot[iPythia]  -> GetYaxis() -> CenterTitle(fCnt);
    hRatioStat[iPythia] -> SetMarkerColor(fColStatR[iPythia]);
    hRatioStat[iPythia] -> SetMarkerStyle(fMarDat);
    hRatioStat[iPythia] -> SetFillColor(fColStatR[iPythia]);
    hRatioStat[iPythia] -> SetFillStyle(fFilStyR[iPythia]);
    hRatioStat[iPythia] -> SetLineColor(fColStatR[iPythia]);
    hRatioStat[iPythia] -> SetLineStyle(fLinDat);
    hRatioStat[iPythia] -> SetLineWidth(fWidDat);
    hRatioStat[iPythia] -> SetTitle(sTitle.Data());
    hRatioStat[iPythia] -> SetTitleFont(fTxt);
    hRatioStat[iPythia] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hRatioStat[iPythia] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatioStat[iPythia] -> GetXaxis() -> SetTitleFont(fTxt);
    hRatioStat[iPythia] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hRatioStat[iPythia] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hRatioStat[iPythia] -> GetXaxis() -> SetLabelFont(fTxt);
    hRatioStat[iPythia] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hRatioStat[iPythia] -> GetXaxis() -> CenterTitle(fCnt);
    hRatioStat[iPythia] -> GetYaxis() -> SetTitle(sTitleRY.Data());
    hRatioStat[iPythia] -> GetYaxis() -> SetTitleFont(fTxt);
    hRatioStat[iPythia] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hRatioStat[iPythia] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hRatioStat[iPythia] -> GetYaxis() -> SetLabelFont(fTxt);
    hRatioStat[iPythia] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hRatioStat[iPythia] -> GetYaxis() -> CenterTitle(fCnt);
  }  // end pythia loop
  cout << "    Set styles." << endl;

  // make legends
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot]  = {0.1, 0.1, 0.3, 0.3};
  const Float_t fLegXYR[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.47};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hUnfoldTot, sLabelTot.Data());
  leg -> AddEntry(hUnfoldStat, sLabelStat.Data());
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    leg -> AddEntry(hPythia[iPythia], sLabelPy[iPythia].Data());
  }  // end pythia loop

  TLegend *legR = new TLegend(fLegXYR[0], fLegXYR[1], fLegXYR[2], fLegXYR[3]);
  legR -> SetFillColor(fColLe);
  legR -> SetFillStyle(fFilLe);
  legR -> SetLineColor(fColLe);
  legR -> SetLineStyle(fLinLe);
  legR -> SetTextFont(fTxt);
  legR -> SetTextAlign(fAln);
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    legR -> AddEntry(hRatioTot[iPythia], sLabelTotR[iPythia].Data());
    legR -> AddEntry(hRatioStat[iPythia], sLabelStatR[iPythia].Data());
  }  // end pythia loop
  cout << "    Made legends." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const Float_t fTxtXY[NPlot * NPlot] = {0.3, 0.1, 0.5, 0.3};
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
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(1);
  const UInt_t  fLinLi(2);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NPlot * NPlot] = {xPlotRange[0], 1., xPlotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

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
  hRatioTot[0]  -> Draw("E5");
  hRatioStat[0] -> Draw("E5 same");
  for (UInt_t iPythia = 1; iPythia < NPythia; iPythia++) {
    hRatioTot[iPythia]  -> Draw("E5 same");
    hRatioStat[iPythia] -> Draw("E5 same");
  }  // end pythia loop
  line          -> Draw();
  legR          -> Draw();
  pPad2         -> cd();
  hUnfoldTot    -> Draw("E5");
  hUnfoldStat   -> Draw("E5 SAME");
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    hPythia[iPythia] -> Draw("SAME");
  }  // end pythia loop
  leg           -> Draw();
  txt           -> Draw();
  fOut          -> cd();
  cPlot         -> Write();
  cPlot         -> Close();
  cout << "    Made plot." << endl;

  // close files
  fOut        -> cd();
  hUnfoldTot  -> Write();
  hUnfoldStat -> Write();
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    hPythia[iPythia]    -> Write();
    hRatioTot[iPythia]  -> Write();
    hRatioStat[iPythia] -> Write();
  }  // end pythia loop
  fOut    -> Close();
  fInTot  -> cd();
  fInTot  -> Close();
  fInStat -> cd();
  fInStat -> Close();
  for (UInt_t iPythia = 0; iPythia < NPythia; iPythia++) {
    fInPy[iPythia] -> cd();
    fInPy[iPythia] -> Close();
  }  // end pythia loop
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
