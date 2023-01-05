// 'CompareTriggerProjections.C'
// Derek Anderson
// 02.25.2022
//
// Use this to compare projections from
// the TES/R matrix against the matched
// eT distribution and weighted pythia
// eT distribution

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"

using namespace std;

// global constants
static const UInt_t NVtx(4);
static const UInt_t NType(5);
static const UInt_t NHist(3);
static const UInt_t NRange(2);
static const Bool_t DoMultiply(false);



void CompareTriggerProjections() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning trigger projection comparison script..." << endl;

  // file parameters
  const TString sOut("projAndSmearDataVsMatchWeightPy8.et650x920gam.d20m3y2022.root");
  const TString sInProj[NHist]   = {"output/March2022/triggerMatrix.withWeights.et650x920vz55tsp0206gam.d6m3y2022.root",
                                    "output/March2022/triggerMatrix.withWeights.et650x920vz55tsp0206gam.d6m3y2022.root",
                                    "output/March2022/triggerMatrix.withWeights.et650x920vz55tsp0206gam.d6m3y2022.root"};
  const TString sInMatch[NHist]  = {"output/January2022/particleGun.full650gamSample_withEtMatrix.et650x650x920vz55tsp0206gam.d17m1y2022.root",
                                    "output/January2022/particleGun.full650gamSample_withEtMatrix.et650x650x920vz55tsp0206gam.d17m1y2022.root",
                                    "output/January2022/particleGun.full650gamSample_withEtMatrix.et650x650x920vz55tsp0206gam.d17m1y2022.root"};
  const TString sInPythia[NHist] = {"input/Weights/weightCheck_eTmatchOverEtInput.et920tsp0206gam.d20m3y2022.root",
                                    "input/Weights/weightCheck_eTmatchOverEtInput.et920tsp0206gam.d20m3y2022.root",
                                    "input/Weights/weightCheck_eTmatchOverEtInput.et920tsp0206gam.d20m3y2022.root"};
  const TString sInData[NHist]   = {"input/Weights/pythiaWeights_usingQtHistToSmear.et911vz55tsp0206gam.d23m3y2021.root",
                                    "input/Weights/pythiaWeights_usingQtHistToSmear.et1115vz55tsp0206gam.d23m3y2021.root",
                                    "input/Weights/pythiaWeights_usingQtHistToSmear.et1520vz55tsp0206gam.d23m3y2021.root"};
  const TString sInPyData[NHist] = {"input/Pythia/pp200py8par.forEtCheck_withWeights_pTbinOne.et911r02gam.d8m3y2022.root",
                                    "input/Pythia/pp200py8par.forEtCheck_withWeights_pTbinOne.et1115r02gam.d8m3y2022.root",
                                    "input/Pythia/pp200py8par.forEtCheck_withWeights_pTbinOne.et1520r02gam.d8m3y2022.root"};

  // hist parameters
  const TString sHistProj[NHist]   = {"hEtProj911",
                                      "hEtProj1115",
                                      "hEtProj1520"};
  const TString sHistMatch[NHist]  = {"match/hEtMatch911",
                                      "match/hEtMatch1115",
                                      "match/hEtMatch1520"};
  const TString sHistPythia[NHist] = {"hEtMultiply_911gam",
                                      "hEtMultiply_1115gam",
                                      "hEtMultiply_1520gam"};
  const TString sHistData[NHist]   = {"hEtData",
                                      "hEtData",
                                      "hEtData"};
  const TString sHistPyData[NHist] = {"Gam/hTrgEtWeightG",
                                      "Gam/hTrgEtWeightG",
                                      "Gam/hTrgEtWeightG"};
  const TString sNameProj[NHist]   = {"hEtProj911",   "hEtProj1115",   "hEtProj1520"};
  const TString sNameMatch[NHist]  = {"hEtMatch911",  "hEtMatch1115",  "hEtMatch1520"};
  const TString sNamePythia[NHist] = {"hEtPythia911", "hEtPythia1115", "hEtPythia1520"};
  const TString sNameWeight[NHist] = {"hEtWeight911", "hEtWeight1115", "hEtWeight1520"};
  const TString sNameData[NHist]   = {"hEtData911",   "hEtData1115",   "hEtData1520"};
  const TString sNamePyData[NHist] = {"hEtPyData911", "hEtPyData1115", "hEtPyData1520"};
  const TString sNameHist[NHist]   = {"hHist911",     "hHist1115",     "hHist1520"};
  const TString sNameType[NType]   = {"hProjType", "hMatchType", "hWeightType", "hDataType", "hPythiaDataType"};

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("E_{T}^{true} [GeV/c]");
  const TString sTitleY("a. u.");
  const Float_t xPlot[NRange]     = {3., 30.};
  const UInt_t  fColProj[NHist]   = {799, 899, 879};
  const UInt_t  fColMatch[NHist]  = {803, 893, 883};
  const UInt_t  fColWeight[NHist] = {798, 898, 878};
  const UInt_t  fColData[NHist]   = {801, 891, 881};
  const UInt_t  fColPyData[NHist] = {794, 894, 874};
  //const UInt_t  fColProj[NHist]   = {859, 839, 819};
  //const UInt_t  fColMatch[NHist]  = {863, 843, 813};
  //const UInt_t  fColWeight[NHist] = {858, 838, 818};
  //const UInt_t  fColData[NHist]   = {861, 841, 811};
  //const UInt_t  fColPyData[NHist] = {854, 834, 814};
  const UInt_t  fColPythia[NHist] = {923, 923, 923};
  const UInt_t  fMarProj[NHist]   = {20, 20, 20};
  const UInt_t  fMarMatch[NHist]  = {24, 24, 24};
  const UInt_t  fMarPythia[NHist] = {29, 29, 29};
  const UInt_t  fMarWeight[NHist] = {28, 28, 28};
  const UInt_t  fMarData[NHist]   = {30, 30, 30};
  const UInt_t  fMarPyData[NHist] = {2, 2, 2};

  // label parameters
  const TString sHeader("#color[899]{Simulated #gamma}");
  //const TString sHeader("#color[859]{Simulated #pi^{0}}");
  const TString sLabelPythia("Unweighted Pythia8");
  const TString sLabelHist[NHist] = {"E_{T}^{select} #in (9, 11) GeV",
                                     "E_{T}^{select} #in (11, 15) GeV",
                                     "E_{T}^{select} #in (15, 20) GeV"};
  const TString sLabelType[NType] = {"Matrix Projections",
                                     "Matched Distributions",
                                     "Pythia8 (E_{T}^{match} weighted)",
                                     "Back-smeared data",
                                     "Pythia8 (data weighted)"};
  const UInt_t  fColType[NType] = {923, 923, 923, 923, 923};
  const UInt_t  fColHist[NHist] = {798, 898, 878};
  //const UInt_t  fColHist[NHist] = {859, 839, 819};
  const UInt_t  fMarType[NType] = {20, 24, 28, 30, 2};
  const UInt_t  fMarHist[NHist] = {1, 1, 1};
  const UInt_t  fFilHist[NHist] = {0, 0, 0};

  // misc parameters
  const UInt_t nRebin[NType]   = {1, 2, 1, 1, 1};
  const Bool_t doRebin[NType]  = {false, true, false, false, false};
  const Bool_t drawHist[NType] = {true, true, true, true, true};

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  TFile *fInProj[NHist];
  TFile *fInMatch[NHist];
  TFile *fInPythia[NHist];
  TFile *fInData[NHist];
  TFile *fInPyData[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInProj[iHist]   = new TFile(sInProj[iHist].Data(), "read");
    fInMatch[iHist]  = new TFile(sInMatch[iHist].Data(), "read");
    fInPythia[iHist] = new TFile(sInPythia[iHist].Data(), "read");
    fInData[iHist]   = new TFile(sInData[iHist].Data(), "read");
    fInPyData[iHist] = new TFile(sInPyData[iHist].Data(), "read");
    if (!fInProj[iHist] || !fInMatch[iHist] || !fInPythia[iHist] || !fInData[iHist] || !fInPyData[iHist]) {
      cerr << "PANIC: couldn't open an input file (#" << iHist << ")!\n"
           << "       fInProj = " << fInProj[iHist] << ", fInMatch = " << fInMatch[iHist] << ", fInPythia = "
           << fInPythia[iHist] << ", finData = " << fInData[iHist] << ", fInPyData = " << fInPyData[iHist] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hProj[NHist];
  TH1D *hMatch[NHist];
  TH1D *hPythia[NHist];
  TH1D *hData[NHist];
  TH1D *hPyData[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hProj[iHist]   = (TH1D*) fInProj[iHist]   -> Get(sHistProj[iHist].Data());
    hMatch[iHist]  = (TH1D*) fInMatch[iHist]  -> Get(sHistMatch[iHist].Data());
    hPythia[iHist] = (TH1D*) fInPythia[iHist] -> Get(sHistPythia[iHist].Data());
    hData[iHist]   = (TH1D*) fInData[iHist]   -> Get(sHistData[iHist].Data());
    hPyData[iHist] = (TH1D*) fInPyData[iHist] -> Get(sHistPyData[iHist].Data());
    if (!hProj[iHist] || !hMatch[iHist] || !hPythia[iHist] || !hData[iHist] || !hPyData[iHist]) {
      cerr << "PANIC: couldn't grab an input histogram (#" << iHist << ")!\n"
           << "       hProj = " << hProj[iHist] << ", hMatch = " << hMatch[iHist] << ", hPythia = " << hPythia[iHist]
           << ", hData = " << hData[iHist] << ", hPyData = " << hPyData[iHist] << "\n"
           << endl;
      return;
    }
    hProj[iHist]   -> SetName(sNameProj[iHist].Data());
    hMatch[iHist]  -> SetName(sNameMatch[iHist].Data());
    hPythia[iHist] -> SetName(sNamePythia[iHist].Data());
    hData[iHist]   -> SetName(sNameData[iHist].Data());
    hPyData[iHist] -> SetName(sNamePyData[iHist].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  Bool_t needToRebin(false);
  for (UInt_t iType = 0; iType < NType; iType++) {
    if (doRebin[iType]) {
      needToRebin = doRebin[iType];
      cout << "    Rebinning histograms:" << endl;
      break;
    }
  }

  Double_t binWidth(1.);
  if (needToRebin) {
    for (UInt_t iType = 0; iType < NHist; iType++) {
      if (doRebin[iType]) {
        for (UInt_t iHist = 0; iHist < NHist; iHist++) {
          switch (iType) {
            case 0:
              hProj[iHist]   -> Rebin(nRebin[iType]);
              binWidth = hProj[iHist]   -> GetBinWidth(1);
              hProj[iHist]   -> Scale(1. / binWidth);
              cout << "      Rebinned projection histogram #" << iHist << "..." << endl;
              break;
            case 1:
              hMatch[iHist]  -> Rebin(nRebin[iType]);
              binWidth = hMatch[iHist]  -> GetBinWidth(1);
              hMatch[iHist]  -> Scale(1. / binWidth);
              cout << "      Rebinned matched histogram #"    << iHist << "..." << endl;
              break;
            case 2:
              hPythia[iHist] -> Rebin(nRebin[iType]);
              binWidth = hPythia[iHist] -> GetBinWidth(1);
              hPythia[iHist] -> Scale(1. / binWidth);
              cout << "      Rebinned projection histogram #" << iHist << "..." << endl;
              break;
            case 3:
              hData[iHist]   -> Rebin(nRebin[iType]);
              binWidth = hData[iHist]   -> GetBinWidth(1);
              hData[iHist]   -> Scale(1. / binWidth);
              cout << "      Rebinned projection histogram #" << iHist << "..." << endl;
              break;
            case 4:
              hPyData[iHist] -> Rebin(nRebin[iType]);
              binWidth = hPyData[iHist] -> GetBinWidth(1);
              hPyData[iHist] -> Scale(1. / binWidth);
              cout << "      Rebinned projection histogram #" << iHist << "..." << endl;
              break;
            default:
              hProj[iHist]   -> Rebin(nRebin[iType]);
              binWidth = hProj[iHist]   -> GetBinWidth(1);
              hProj[iHist]   -> Scale(1. / binWidth);
              cout << "      Rebinned projection histogram #" << iHist << "..." << endl;
              break;
          }
        }  // end histogram loop
      }
    }  // end type loop
    cout << "    Finished rebinning histograms." << endl;
  }

  // normalize pythia
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    const Double_t iPythia = hPythia[iHist] -> Integral();
    const Double_t normP   = 1. / iPythia;
    if (normP > 0.) hPythia[iHist] -> Scale(normP);
  }
  cout << "    Normalized pythia." << endl;

  // apply weighting to pythia
  if (DoMultiply) cout << "    Weighting pythia histograms:" << endl;

  TH1D *hWeight[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hWeight[iHist] = (TH1D*) hPythia[iHist] -> Clone();
    hWeight[iHist] -> SetName(sNameWeight[iHist].Data());
    if (DoMultiply) {
      hWeight[iHist] -> Multiply(hProj[iHist]);
      cout << "      Weighted pythia histogram #" << iHist << "..." << endl;
    }
  }

  // QUICK FIX [03.18.2022]
  //for (UInt_t iHist = 0; iHist < NHist; iHist++) {
  //  const Double_t iMatch = hMatch[iHist] -> Integral();
  //  const Double_t normM  = 1. / iMatch;
  //  if (normM > 0.) hMatch[iHist] -> Scale(normM);
  //}
  //cout << "    Normalized matched distributions." << endl;

  // normalize histograms to matched distributions
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    const Double_t iProj   = hProj[iHist]   -> Integral();
    const Double_t iMatch  = hMatch[iHist]  -> Integral();
    const Double_t iWeight = hWeight[iHist] -> Integral();
    const Double_t iData   = hData[iHist]   -> Integral();
    const Double_t iPyData = hPyData[iHist] -> Integral();
    const Double_t normP   = iProj / iProj;
    const Double_t normM   = iProj / iMatch;
    const Double_t normW   = iProj / iWeight;
    const Double_t normD   = iProj / iData;
    const Double_t normPD  = iProj / iPyData;
    if (normM  > 0.) hMatch[iHist]  -> Scale(normM);
    if (normW  > 0.) hWeight[iHist] -> Scale(normW);
    if (normD  > 0.) hData[iHist]   -> Scale(normD);
    if (normPD > 0.) hPyData[iHist] -> Scale(normPD);
  }
  cout << "    Normalized other histograms." << endl;

  // set styles
  TH1D *hType[NType];
  for (UInt_t iType = 0; iType < NType; iType++) {
    switch (iType) {
      case 0:
        hType[iType] = (TH1D*) hProj[0]   -> Clone();
        hType[iType] -> SetName(sNameType[iType].Data());
        break;
      case 1:
        hType[iType] = (TH1D*) hMatch[0]  -> Clone();
        hType[iType] -> SetName(sNameType[iType].Data());
        break;
      case 2:
        hType[iType] = (TH1D*) hWeight[0] -> Clone();
        hType[iType] -> SetName(sNameType[iType].Data());
        break;
      case 3:
        hType[iType] = (TH1D*) hData[0]   -> Clone();
        hType[iType] -> SetName(sNameType[iType].Data());
        break;
      case 4:
        hType[iType] = (TH1D*) hPyData[0] -> Clone();
        hType[iType] -> SetName(sNameType[iType].Data());
        break;
      default:
        hType[iType] = (TH1D*) hProj[0]   -> Clone();
        hType[iType] -> SetName(sNameType[iType].Data());
        break;
    }
  }  // end type loop

  TH1D *hHist[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hHist[iHist] = (TH1D*) hProj[0] -> Clone();
    hHist[iHist] -> SetName(sNameHist[iHist].Data());
  }

  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fLinH(3);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.3);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hProj[iHist]   -> SetMarkerColor(fColProj[iHist]);
    hProj[iHist]   -> SetMarkerStyle(fMarProj[iHist]);
    hProj[iHist]   -> SetFillColor(fColProj[iHist]);
    hProj[iHist]   -> SetFillStyle(fFil);
    hProj[iHist]   -> SetLineColor(fColProj[iHist]);
    hProj[iHist]   -> SetLineStyle(fLin);
    hProj[iHist]   -> SetLineWidth(fWid);
    hProj[iHist]   -> SetTitle(sTitle.Data());
    hProj[iHist]   -> SetTitleFont(fTxt);
    hProj[iHist]   -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hProj[iHist]   -> GetXaxis() -> SetTitle(sTitleX.Data());
    hProj[iHist]   -> GetXaxis() -> SetTitleFont(fTxt);
    hProj[iHist]   -> GetXaxis() -> SetTitleSize(fTit);
    hProj[iHist]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hProj[iHist]   -> GetXaxis() -> SetLabelFont(fTxt);
    hProj[iHist]   -> GetXaxis() -> SetLabelSize(fLab);
    hProj[iHist]   -> GetXaxis() -> CenterTitle(fCnt);
    hProj[iHist]   -> GetYaxis() -> SetTitle(sTitleY.Data());
    hProj[iHist]   -> GetYaxis() -> SetTitleFont(fTxt);
    hProj[iHist]   -> GetYaxis() -> SetTitleSize(fTit);
    hProj[iHist]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hProj[iHist]   -> GetYaxis() -> SetLabelFont(fTxt);
    hProj[iHist]   -> GetYaxis() -> SetLabelSize(fLab);
    hProj[iHist]   -> GetYaxis() -> CenterTitle(fCnt);
    hMatch[iHist]  -> SetMarkerColor(fColMatch[iHist]);
    hMatch[iHist]  -> SetMarkerStyle(fMarMatch[iHist]);
    hMatch[iHist]  -> SetFillColor(fColMatch[iHist]);
    hMatch[iHist]  -> SetFillStyle(fFil);
    hMatch[iHist]  -> SetLineColor(fColMatch[iHist]);
    hMatch[iHist]  -> SetLineStyle(fLin);
    hMatch[iHist]  -> SetLineWidth(fWid);
    hMatch[iHist]  -> SetTitle(sTitle.Data());
    hMatch[iHist]  -> SetTitleFont(fTxt);
    hMatch[iHist]  -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hMatch[iHist]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hMatch[iHist]  -> GetXaxis() -> SetTitleFont(fTxt);
    hMatch[iHist]  -> GetXaxis() -> SetTitleSize(fTit);
    hMatch[iHist]  -> GetXaxis() -> SetTitleOffset(fOffX);
    hMatch[iHist]  -> GetXaxis() -> SetLabelFont(fTxt);
    hMatch[iHist]  -> GetXaxis() -> SetLabelSize(fLab);
    hMatch[iHist]  -> GetXaxis() -> CenterTitle(fCnt);
    hMatch[iHist]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    hMatch[iHist]  -> GetYaxis() -> SetTitleFont(fTxt);
    hMatch[iHist]  -> GetYaxis() -> SetTitleSize(fTit);
    hMatch[iHist]  -> GetYaxis() -> SetTitleOffset(fOffY);
    hMatch[iHist]  -> GetYaxis() -> SetLabelFont(fTxt);
    hMatch[iHist]  -> GetYaxis() -> SetLabelSize(fLab);
    hMatch[iHist]  -> GetYaxis() -> CenterTitle(fCnt);
    hPythia[iHist] -> SetMarkerColor(fColPythia[iHist]);
    hPythia[iHist] -> SetMarkerStyle(fMarPythia[iHist]);
    hPythia[iHist] -> SetFillColor(fColPythia[iHist]);
    hPythia[iHist] -> SetFillStyle(fFil);
    hPythia[iHist] -> SetLineColor(fColPythia[iHist]);
    hPythia[iHist] -> SetLineStyle(fLin);
    hPythia[iHist] -> SetLineWidth(fWid);
    hPythia[iHist] -> SetTitle(sTitle.Data());
    hPythia[iHist] -> SetTitleFont(fTxt);
    hPythia[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hPythia[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPythia[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hPythia[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hPythia[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPythia[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hPythia[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hPythia[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hPythia[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPythia[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hPythia[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hPythia[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPythia[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hPythia[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hPythia[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hWeight[iHist] -> SetMarkerColor(fColWeight[iHist]);
    hWeight[iHist] -> SetMarkerStyle(fMarWeight[iHist]);
    hWeight[iHist] -> SetFillColor(fColWeight[iHist]);
    hWeight[iHist] -> SetFillStyle(fFil);
    hWeight[iHist] -> SetLineColor(fColWeight[iHist]);
    hWeight[iHist] -> SetLineStyle(fLin);
    hWeight[iHist] -> SetLineWidth(fWid);
    hWeight[iHist] -> SetTitle(sTitle.Data());
    hWeight[iHist] -> SetTitleFont(fTxt);
    hWeight[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hWeight[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hWeight[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hWeight[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hWeight[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hWeight[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hWeight[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hWeight[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hWeight[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hWeight[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hWeight[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hWeight[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hWeight[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hWeight[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hWeight[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hData[iHist]   -> SetMarkerColor(fColData[iHist]);
    hData[iHist]   -> SetMarkerStyle(fMarData[iHist]);
    hData[iHist]   -> SetFillColor(fColData[iHist]);
    hData[iHist]   -> SetFillStyle(fFil);
    hData[iHist]   -> SetLineColor(fColData[iHist]);
    hData[iHist]   -> SetLineStyle(fLin);
    hData[iHist]   -> SetLineWidth(fWid);
    hData[iHist]   -> SetTitle(sTitle.Data());
    hData[iHist]   -> SetTitleFont(fTxt);
    hData[iHist]   -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hData[iHist]   -> GetXaxis() -> SetTitle(sTitleX.Data());
    hData[iHist]   -> GetXaxis() -> SetTitleFont(fTxt);
    hData[iHist]   -> GetXaxis() -> SetTitleSize(fTit);
    hData[iHist]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hData[iHist]   -> GetXaxis() -> SetLabelFont(fTxt);
    hData[iHist]   -> GetXaxis() -> SetLabelSize(fLab);
    hData[iHist]   -> GetXaxis() -> CenterTitle(fCnt);
    hData[iHist]   -> GetYaxis() -> SetTitle(sTitleY.Data());
    hData[iHist]   -> GetYaxis() -> SetTitleFont(fTxt);
    hData[iHist]   -> GetYaxis() -> SetTitleSize(fTit);
    hData[iHist]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hData[iHist]   -> GetYaxis() -> SetLabelFont(fTxt);
    hData[iHist]   -> GetYaxis() -> SetLabelSize(fLab);
    hData[iHist]   -> GetYaxis() -> CenterTitle(fCnt);
    hPyData[iHist] -> SetMarkerColor(fColPyData[iHist]);
    hPyData[iHist] -> SetMarkerStyle(fMarPyData[iHist]);
    hPyData[iHist] -> SetFillColor(fColPyData[iHist]);
    hPyData[iHist] -> SetFillStyle(fFil);
    hPyData[iHist] -> SetLineColor(fColPyData[iHist]);
    hPyData[iHist] -> SetLineStyle(fLin);
    hPyData[iHist] -> SetLineWidth(fWid);
    hPyData[iHist] -> SetTitle(sTitle.Data());
    hPyData[iHist] -> SetTitleFont(fTxt);
    hPyData[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hPyData[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPyData[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hPyData[iHist] -> GetXaxis() -> SetTitleSize(fTit);
    hPyData[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPyData[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hPyData[iHist] -> GetXaxis() -> SetLabelSize(fLab);
    hPyData[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hPyData[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPyData[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hPyData[iHist] -> GetYaxis() -> SetTitleSize(fTit);
    hPyData[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPyData[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hPyData[iHist] -> GetYaxis() -> SetLabelSize(fLab);
    hPyData[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hHist[iHist]   -> SetMarkerColor(fColHist[iHist]);
    hHist[iHist]   -> SetMarkerStyle(fMarHist[iHist]);
    hHist[iHist]   -> SetFillColor(fColHist[iHist]);
    hHist[iHist]   -> SetFillStyle(fFilHist[iHist]);
    hHist[iHist]   -> SetLineColor(fColHist[iHist]);
    hHist[iHist]   -> SetLineStyle(fLinH);
    hHist[iHist]   -> SetLineWidth(fWid);
    hHist[iHist]   -> SetTitle(sTitle.Data());
    hHist[iHist]   -> SetTitleFont(fTxt);
    hHist[iHist]   -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hHist[iHist]   -> GetXaxis() -> SetTitle(sTitleX.Data());
    hHist[iHist]   -> GetXaxis() -> SetTitleFont(fTxt);
    hHist[iHist]   -> GetXaxis() -> SetTitleSize(fTit);
    hHist[iHist]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hHist[iHist]   -> GetXaxis() -> SetLabelFont(fTxt);
    hHist[iHist]   -> GetXaxis() -> SetLabelSize(fLab);
    hHist[iHist]   -> GetXaxis() -> CenterTitle(fCnt);
    hHist[iHist]   -> GetYaxis() -> SetTitle(sTitleY.Data());
    hHist[iHist]   -> GetYaxis() -> SetTitleFont(fTxt);
    hHist[iHist]   -> GetYaxis() -> SetTitleSize(fTit);
    hHist[iHist]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hHist[iHist]   -> GetYaxis() -> SetLabelFont(fTxt);
    hHist[iHist]   -> GetYaxis() -> SetLabelSize(fLab);
    hHist[iHist]   -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iType = 0; iType < NType; iType++) {
    hType[iType] -> SetMarkerColor(fColType[iType]);
    hType[iType] -> SetMarkerStyle(fMarType[iType]);
    hType[iType] -> SetFillColor(fColType[iType]);
    hType[iType] -> SetFillStyle(fFil);
    hType[iType] -> SetLineColor(fColType[iType]);
    hType[iType] -> SetLineStyle(fLin);
    hType[iType] -> SetLineWidth(fWid);
    hType[iType] -> SetTitle(sTitle.Data());
    hType[iType] -> SetTitleFont(fTxt);
    hType[iType] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hType[iType] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hType[iType] -> GetXaxis() -> SetTitleFont(fTxt);
    hType[iType] -> GetXaxis() -> SetTitleSize(fTit);
    hType[iType] -> GetXaxis() -> SetTitleOffset(fOffX);
    hType[iType] -> GetXaxis() -> SetLabelFont(fTxt);
    hType[iType] -> GetXaxis() -> SetLabelSize(fLab);
    hType[iType] -> GetXaxis() -> CenterTitle(fCnt);
    hType[iType] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hType[iType] -> GetYaxis() -> SetTitleFont(fTxt);
    hType[iType] -> GetYaxis() -> SetTitleSize(fTit);
    hType[iType] -> GetYaxis() -> SetTitleOffset(fOffY);
    hType[iType] -> GetYaxis() -> SetLabelFont(fTxt);
    hType[iType] -> GetYaxis() -> SetLabelSize(fLab);
    hType[iType] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set Styles." << endl;

  // count number of histograms to draw
  UInt_t nToDraw(0);
  for (UInt_t iType = 0; iType < NType; iType++) {
    if (drawHist[iType]) nToDraw++;
  }
  cout << "    Counted number of histograms to draw:\n"
       << "      number to draw = " << nToDraw
       << endl;

  // make legends
  const UInt_t  fColLe        = 0;
  const UInt_t  fFilLe        = 0;
  const UInt_t  fLinLe        = 0;
  const UInt_t  nColLe        = 2;
  const UInt_t  nObjLe        = 1 + TMath::Max(NType, nToDraw);
  const UInt_t  nObjLeP       = 2 + NHist;
  const Float_t hObj          = 0.05;
  const Float_t hObjLe        = nObjLe * hObj;
  const Float_t hObjLeP       = nObjLeP * hObj;
  const Float_t yObjLe        = 0.1 + hObjLe;
  const Float_t yObjLeP       = 0.1 + hObjLeP;
  const Float_t fLegXY[NVtx]  = {0.1, 0.1, 0.5, yObjLe};
  const Float_t fLegXyP[NVtx] = {0.1, 0.1, 0.3, yObjLeP};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sHeader.Data());
  leg -> SetNColumns(nColLe);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);

  UInt_t iHistLe(0);
  UInt_t iTypeLe(0);
  Bool_t typeDone(false);
  Bool_t histDone(false);
  Bool_t doneFilling(false);
  do {

    // add type histogram
    if (iTypeLe < NType) {
      if (drawHist[iTypeLe]) leg -> AddEntry(hType[iTypeLe], sLabelType[iTypeLe].Data(), "pf");
    } else {
      leg -> AddEntry((TObject*)0, "", "");
    }
    ++iTypeLe;

    // add eT histogram
    if (iHistLe < NHist) {
      leg -> AddEntry(hHist[iHistLe], sLabelHist[iHistLe].Data(), "f");
    } else {
      leg -> AddEntry((TObject*)0, "", "");
    }
    ++iHistLe;

    // check if done
    typeDone    = (iTypeLe > nToDraw);
    histDone    = (iHistLe > NHist);
    doneFilling = (typeDone && histDone);

  }  while (!doneFilling);

  TLegend *legP = new TLegend(fLegXyP[0], fLegXyP[1], fLegXyP[2], fLegXyP[3], sHeader.Data());
  legP -> SetFillColor(fColLe);
  legP -> SetFillStyle(fFilLe);
  legP -> SetLineColor(fColLe);
  legP -> SetLineStyle(fLinLe);
  legP -> SetTextFont(fTxt);
  legP -> SetTextAlign(fAln);
  legP -> AddEntry(hPythia[0], sLabelPythia.Data(), "pf");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    legP -> AddEntry(hWeight[iHist], sLabelHist[iHist].Data(), "pf");
  }
  cout << "    Made legend." << endl;

  // make plot
  const UInt_t  width(950);
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

  TCanvas *cComp = new TCanvas("cComparison", "", width, height);
  cComp      -> SetGrid(fGrid, fGrid);
  cComp      -> SetTicks(fTick, fTick);
  cComp      -> SetBorderMode(fMode);
  cComp      -> SetBorderSize(fBord);
  cComp      -> SetGrid(fGrid, fGrid);
  cComp      -> SetTicks(fTick, fTick);
  cComp      -> SetLogx(fLogX);
  cComp      -> SetLogy(fLogY);
  cComp      -> SetBorderMode(fMode);
  cComp      -> SetBorderSize(fBord);
  cComp      -> SetFrameBorderMode(fFrame);
  cComp      -> SetLeftMargin(fMarginL);
  cComp      -> SetRightMargin(fMarginR);
  cComp      -> SetTopMargin(fMarginT);
  cComp      -> SetBottomMargin(fMarginB);
  cComp      -> cd();
  if (drawHist[0]) hProj[0]   -> Draw("hist p l");
  if (drawHist[1]) hMatch[0]  -> Draw("hist p l same");
  if (drawHist[2]) hWeight[0] -> Draw("hist p l same");
  if (drawHist[3]) hData[0]   -> Draw("hist p l same");
  if (drawHist[4]) hPyData[0] -> Draw("hist p l same");
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    if (drawHist[0]) hProj[iHist]   -> Draw("hist p l same");
    if (drawHist[1]) hMatch[iHist]  -> Draw("hist p l same");
    if (drawHist[2]) hWeight[iHist] -> Draw("hist p l same");
    if (drawHist[3]) hData[iHist]   -> Draw("hist p l same");
    if (drawHist[4]) hPyData[iHist] -> Draw("hist p l same");
  }
  leg   -> Draw();
  fOut  -> cd();
  cComp -> Write();
  cComp -> Close();

  TCanvas *cWeight = new TCanvas("cWeightVsNotPythia", "", width, height);
  cWeight    -> SetGrid(fGrid, fGrid);
  cWeight    -> SetTicks(fTick, fTick);
  cWeight    -> SetBorderMode(fMode);
  cWeight    -> SetBorderSize(fBord);
  cWeight    -> SetGrid(fGrid, fGrid);
  cWeight    -> SetTicks(fTick, fTick);
  cWeight    -> SetLogx(fLogX);
  cWeight    -> SetLogy(fLogY);
  cWeight    -> SetBorderMode(fMode);
  cWeight    -> SetBorderSize(fBord);
  cWeight    -> SetFrameBorderMode(fFrame);
  cWeight    -> SetLeftMargin(fMarginL);
  cWeight    -> SetRightMargin(fMarginR);
  cWeight    -> SetTopMargin(fMarginT);
  cWeight    -> SetBottomMargin(fMarginB);
  cWeight    -> cd();
  hPythia[0] -> Draw();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hWeight[iHist] -> Draw("same");
  }
  legP    -> Draw();
  fOut    -> cd();
  cWeight -> Write();
  cWeight -> Close();
  cout << "    Made plot." << endl;

  // write histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hProj[iHist]   -> Write();
    hMatch[iHist]  -> Write();
    hPythia[iHist] -> Write();
    hWeight[iHist] -> Write();
    hData[iHist]   -> Write();
    hPyData[iHist] -> Write();
    hHist[iHist]   -> Write();
  }
  for (UInt_t iType = 0; iType < NType; iType++) {
    hType[iType] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fInProj[iHist]   -> cd();
    fInProj[iHist]   -> Close();
    fInMatch[iHist]  -> cd();
    fInMatch[iHist]  -> Close();
    fInPythia[iHist] -> cd();
    fInPythia[iHist] -> Close();
    fInData[iHist]   -> cd();
    fInData[iHist]   -> Close();
    fInPyData[iHist] -> cd();
    fInPyData[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
