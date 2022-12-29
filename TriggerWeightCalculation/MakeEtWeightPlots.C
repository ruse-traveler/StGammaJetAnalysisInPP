// 'MakeEtWeightPlots.C'
// Derek Anderson
// 10.13.2021
//
// Makes plots of TES/R weights
//   TrgTyp = 0: pi0
//   TrgTyp = 1: photon

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;

// global constants
const UInt_t NTrg(3);
const UInt_t NPts(2);
const UInt_t NVtx(4);
const UInt_t NHist((NTrg * 3) + 1);
const UInt_t TrgTyp(1);



void MakeEtWeightPlots() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning trigger weight plot script..." << endl;

  // io parameters
  const TString sOut("weightCheck_qTsmearOverTruePy8.et920x650tsp0206gam.d20m3y2022.root");
  const TString sData("hEtData");
  const TString sPythia("hEtPythia");
  const TString sWeight("hEtWeight");
  const TString sInF[NTrg]   = {"weightCheck_qTsmearOverTruePy8_eTbinOne.et911gam.d20m3y2022.root",
                                "weightCheck_qTsmearOverTruePy8_eTbinOne.et1115gam.d20m3y2022.root",
                                "weightCheck_qTsmearOverTruePy8_eTbinOne.et1520gam.d20m3y2022.root"};
  const TString sInH[NTrg]   = {"weightCheck_qTsmearOverTruePy8_eTbinOne.et911gam.d20m3y2022.root",
                                "weightCheck_qTsmearOverTruePy8_eTbinOne.et1115gam.d20m3y2022.root",
                                "weightCheck_qTsmearOverTruePy8_eTbinOne.et1520gam.d20m3y2022.root"};

  // histogram parameters
  const TString sTitle("");
  const TString sTitleX("E_{T}^{trg} [GeV]");
  const TString sTitleYH("arbitrary");
  const TString sTitleYW("weights [q_{T}^{trg}-smeared / E_{T}^{true}]");
  const TString sNameP("hPythia");
  const TString sNameDF[NTrg] = {"hMatchFit911", "hMatchFit1115", "hMatchFit1520"};
  const TString sNameDH[NTrg] = {"hMatchHist911", "hMatchHist1115", "hMatchHist1520"};
  const TString sNameWF[NTrg] = {"hWeightFit911", "hWeightFit1115", "hWeightFit1520"};
  const TString sNameWH[NTrg] = {"hWeightHist911", "hWeightHist1115", "hWeightHist1520"};

  // style parameters
  const UInt_t  fMarP(29);
  const UInt_t  fColP(923);
  const UInt_t  fMarF[NTrg] = {20, 22, 21};
  const UInt_t  fMarH[NTrg] = {24, 26, 25};
  const Float_t xPlot[NPts] = {4., 40.};

  // text parameters
  const TString sCol("#bf{#color[899]{simulated #gamma}}");
  const TString sLabelP("Particle-level Pythia8 [E_{T}^{trg} > 9 GeV]");
  const TString sLabelF("Matched E_{T}");
  const TString sLabelH("Matched E_{T}");
  const TString sLabelE[NTrg] = {"#bf{E_{T}^{select} #in (9, 11) GeV}", "#bf{E_{T}^{select} #in (11, 15) GeV}", "#bf{E_{T}^{select} #in (15, 20) GeV}"};

  // trigger parameters
  const TString sTrgP("#pi^{0} trig., E_{T}^{select} #in (9, 20) GeV");
  const TString sTrgG("#gamma trig., E_{T}^{select} #in (9, 20) GeV");
  const TString sCutP("|#eta^{trg}| < 0.9, TSP #in (0, 0.08)");
  const TString sCutG("|#eta^{trg}| < 0.9, TSP #in (0.2, 0.6)");
  const UInt_t  fColPF[NTrg] = {859, 839, 819};
  const UInt_t  fColPH[NTrg] = {855, 835, 815};
  const UInt_t  fColGF[NTrg] = {799, 899, 879};
  const UInt_t  fColGH[NTrg] = {795, 895, 875};

  // set trigger
  TString sTrg("");
  TString sCut("");
  UInt_t  fMarF[NTrg];
  UInt_t  fMarH[NTrg];
  UInt_t  fColF[NTrg];
  UInt_t  fColH[NTrg];
  switch (TrgTyp) {
    case 0:
      sTrg = sTrgP;
      sCut = sCutP;
      for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
        fColF[iTrg] = fColPF[iTrg];
        fColH[iTrg] = fColPH[iTrg];
      }
      break;
    case 1:
      sTrg = sTrgG;
      sCut = sCutG;
      for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
        fColF[iTrg] = fColGF[iTrg];
        fColH[iTrg] = fColGH[iTrg];
      }
      break;
    default:
      sTrg = sTrgP;
      sCut = sCutP;
      for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
        fColF[iTrg] = fColPF[iTrg];
        fColH[iTrg] = fColPH[iTrg];
      }
      break;
  }

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  TFile *fInF[NTrg];
  TFile *fInH[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fInF[iTrg] = new TFile(sInF[iTrg].Data(), "read");
    fInH[iTrg] = new TFile(sInH[iTrg].Data(), "read");
    if (!fInF[iTrg] || !fInH[iTrg]) {
      cerr << "PANIC: couldn't open an input file!\n"
           << "       fInF[" << iTrg << "] = " << fInF[iTrg] << ", fInH[" << iTrg << "] = " << fInH[iTrg] << "\n"
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hPythia = (TH1D*) fInF[0] -> Get(sPythia.Data());
  if (!hPythia) {
    cerr << "PANIC: couldn't grab pythia histogram!\n" << endl;
    return;
  }
  hPythia -> SetName(sNameP.Data());

  TH1D *hDataF[NTrg];
  TH1D *hDataH[NTrg];
  TH1D *hWeightF[NTrg];
  TH1D *hWeightH[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hDataF[iTrg]   = (TH1D*) fInF[iTrg] -> Get(sData.Data());
    hDataH[iTrg]   = (TH1D*) fInH[iTrg] -> Get(sData.Data());
    hWeightF[iTrg] = (TH1D*) fInF[iTrg] -> Get(sWeight.Data());
    hWeightH[iTrg] = (TH1D*) fInH[iTrg] -> Get(sWeight.Data());
    if (!hDataF[iTrg] || !hDataH[iTrg] || !hWeightF[iTrg] || !hWeightH[iTrg]) {
      cerr << "PANIC: couldn't grab data or weight histogram!\n"
           << "        hDataF[" << iTrg << "]   = " << hDataF[iTrg] << ", hDataH[" << iTrg << "]     = " << hDataH[iTrg] << "\n"
           << "        hWeightF[" << iTrg << "] = " << hWeightF[iTrg] << ", hWeightH[" << iTrg << "] = " << hWeightH[iTrg] << "\n"
           << endl;
      return;
    }
    hDataF[iTrg]   -> SetName(sNameDF[iTrg].Data());
    hDataH[iTrg]   -> SetName(sNameDH[iTrg].Data());
    hWeightF[iTrg] -> SetName(sNameWF[iTrg].Data());
    hWeightH[iTrg] -> SetName(sNameWH[iTrg].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fTit(0.04);
  const Float_t fLab(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.6);
  hPythia -> SetMarkerColor(fColP);
  hPythia -> SetMarkerStyle(fMarP);
  hPythia -> SetFillColor(fColP);
  hPythia -> SetFillStyle(fFil);
  hPythia -> SetLineColor(fColP);
  hPythia -> SetLineStyle(fLin);
  hPythia -> SetLineWidth(fWid);
  hPythia -> SetTitle(sTitle.Data());
  hPythia -> SetTitleFont(fTxt);
  hPythia -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
  hPythia -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPythia -> GetXaxis() -> SetTitleFont(fTxt);
  hPythia -> GetXaxis() -> SetTitleSize(fTit);
  hPythia -> GetXaxis() -> SetTitleOffset(fOffX);
  hPythia -> GetXaxis() -> SetLabelFont(fTxt);
  hPythia -> GetXaxis() -> SetLabelSize(fLab);
  hPythia -> GetXaxis() -> CenterTitle(fCnt);
  hPythia -> GetYaxis() -> SetTitle(sTitleYH.Data());
  hPythia -> GetYaxis() -> SetTitleFont(fTxt);
  hPythia -> GetYaxis() -> SetTitleSize(fTit);
  hPythia -> GetYaxis() -> SetTitleOffset(fOffY);
  hPythia -> GetYaxis() -> SetLabelFont(fTxt);
  hPythia -> GetYaxis() -> SetLabelSize(fLab);
  hPythia -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hDataF[iTrg]   -> SetMarkerColor(fColF[iTrg]);
    hDataF[iTrg]   -> SetMarkerStyle(fMarF[iTrg]);
    hDataF[iTrg]   -> SetFillColor(fColF[iTrg]);
    hDataF[iTrg]   -> SetFillStyle(fFil);
    hDataF[iTrg]   -> SetLineColor(fColF[iTrg]);
    hDataF[iTrg]   -> SetLineStyle(fLin);
    hDataF[iTrg]   -> SetLineWidth(fWid);
    hDataF[iTrg]   -> SetTitle(sTitle.Data());
    hDataF[iTrg]   -> SetTitleFont(fTxt);
    hDataF[iTrg]   -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hDataF[iTrg]   -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDataF[iTrg]   -> GetXaxis() -> SetTitleFont(fTxt);
    hDataF[iTrg]   -> GetXaxis() -> SetTitleSize(fTit);
    hDataF[iTrg]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hDataF[iTrg]   -> GetXaxis() -> SetLabelFont(fTxt);
    hDataF[iTrg]   -> GetXaxis() -> SetLabelSize(fLab);
    hDataF[iTrg]   -> GetXaxis() -> CenterTitle(fCnt);
    hDataF[iTrg]   -> GetYaxis() -> SetTitle(sTitleYH.Data());
    hDataF[iTrg]   -> GetYaxis() -> SetTitleFont(fTxt);
    hDataF[iTrg]   -> GetYaxis() -> SetTitleSize(fTit);
    hDataF[iTrg]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hDataF[iTrg]   -> GetYaxis() -> SetLabelFont(fTxt);
    hDataF[iTrg]   -> GetYaxis() -> SetLabelSize(fLab);
    hDataF[iTrg]   -> GetYaxis() -> CenterTitle(fCnt);
    hDataH[iTrg]   -> SetMarkerColor(fColH[iTrg]);
    hDataH[iTrg]   -> SetMarkerStyle(fMarH[iTrg]);
    hDataH[iTrg]   -> SetFillColor(fColH[iTrg]);
    hDataH[iTrg]   -> SetFillStyle(fFil);
    hDataH[iTrg]   -> SetLineColor(fColH[iTrg]);
    hDataH[iTrg]   -> SetLineStyle(fLin);
    hDataH[iTrg]   -> SetLineWidth(fWid);
    hDataH[iTrg]   -> SetTitle(sTitle.Data());
    hDataH[iTrg]   -> SetTitleFont(fTxt);
    hDataH[iTrg]   -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hDataH[iTrg]   -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDataH[iTrg]   -> GetXaxis() -> SetTitleFont(fTxt);
    hDataH[iTrg]   -> GetXaxis() -> SetTitleSize(fTit);
    hDataH[iTrg]   -> GetXaxis() -> SetTitleOffset(fOffX);
    hDataH[iTrg]   -> GetXaxis() -> SetLabelFont(fTxt);
    hDataH[iTrg]   -> GetXaxis() -> SetLabelSize(fLab);
    hDataH[iTrg]   -> GetXaxis() -> CenterTitle(fCnt);
    hDataH[iTrg]   -> GetYaxis() -> SetTitle(sTitleYH.Data());
    hDataH[iTrg]   -> GetYaxis() -> SetTitleFont(fTxt);
    hDataH[iTrg]   -> GetYaxis() -> SetTitleSize(fTit);
    hDataH[iTrg]   -> GetYaxis() -> SetTitleOffset(fOffY);
    hDataH[iTrg]   -> GetYaxis() -> SetLabelFont(fTxt);
    hDataH[iTrg]   -> GetYaxis() -> SetLabelSize(fLab);
    hDataH[iTrg]   -> GetYaxis() -> CenterTitle(fCnt);
    hWeightF[iTrg] -> SetMarkerColor(fColF[iTrg]);
    hWeightF[iTrg] -> SetMarkerStyle(fMarF[iTrg]);
    hWeightF[iTrg] -> SetFillColor(fColF[iTrg]);
    hWeightF[iTrg] -> SetFillStyle(fFil);
    hWeightF[iTrg] -> SetLineColor(fColF[iTrg]);
    hWeightF[iTrg] -> SetLineStyle(fLin);
    hWeightF[iTrg] -> SetLineWidth(fWid);
    hWeightF[iTrg] -> SetTitle(sTitle.Data());
    hWeightF[iTrg] -> SetTitleFont(fTxt);
    hWeightF[iTrg] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hWeightF[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hWeightF[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hWeightF[iTrg] -> GetXaxis() -> SetTitleSize(fTit);
    hWeightF[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hWeightF[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hWeightF[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hWeightF[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hWeightF[iTrg] -> GetYaxis() -> SetTitle(sTitleYW.Data());
    hWeightF[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hWeightF[iTrg] -> GetYaxis() -> SetTitleSize(fTit);
    hWeightF[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hWeightF[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hWeightF[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hWeightF[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
    hWeightH[iTrg] -> SetMarkerColor(fColH[iTrg]);
    hWeightH[iTrg] -> SetMarkerStyle(fMarH[iTrg]);
    hWeightH[iTrg] -> SetFillColor(fColH[iTrg]);
    hWeightH[iTrg] -> SetFillStyle(fFil);
    hWeightH[iTrg] -> SetLineColor(fColH[iTrg]);
    hWeightH[iTrg] -> SetLineStyle(fLin);
    hWeightH[iTrg] -> SetLineWidth(fWid);
    hWeightH[iTrg] -> SetTitle(sTitle.Data());
    hWeightH[iTrg] -> SetTitleFont(fTxt);
    hWeightH[iTrg] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hWeightH[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hWeightH[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hWeightH[iTrg] -> GetXaxis() -> SetTitleSize(fTit);
    hWeightH[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hWeightH[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hWeightH[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hWeightH[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hWeightH[iTrg] -> GetYaxis() -> SetTitle(sTitleYW.Data());
    hWeightH[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hWeightH[iTrg] -> GetYaxis() -> SetTitleSize(fTit);
    hWeightH[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hWeightH[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hWeightH[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hWeightH[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fSizH(NHist * 0.05);
  const Float_t fSizW((NHist - 1) * 0.05);
  const Float_t fSizW((NHist - 1) * 0.05);
  const Float_t fTopH(0.1 + fSizH);
  const Float_t fTopW(0.1 + fSizW);
  const Float_t fLegXyH[NVtx] = {0.1, 0.1, 0.3, fTopH};
  const Float_t fLegXyW[NVtx] = {0.1, 0.1, 0.3, fTopW};

  // for headers
  TH1D *hEmpty = (TH1D*) hPythia -> Clone();
  hEmpty -> SetName("hEmpty");
  hEmpty -> Reset("ICES");

  TLegend *legH = new TLegend(fLegXyH[0], fLegXyH[1], fLegXyH[2], fLegXyH[3]);
  TLegend *legW = new TLegend(fLegXyW[0], fLegXyW[1], fLegXyW[2], fLegXyW[3]);
  legH -> SetFillColor(fColLe);
  legH -> SetFillStyle(fFilLe);
  legH -> SetLineColor(fColLe);
  legH -> SetLineStyle(fLinLe);
  legH -> SetTextFont(fTxt);
  legH -> SetTextAlign(fAln);
  legW -> SetFillColor(fColLe);
  legW -> SetFillStyle(fFilLe);
  legW -> SetLineColor(fColLe);
  legW -> SetLineStyle(fLinLe);
  legW -> SetTextFont(fTxt);
  legW -> SetTextAlign(fAln);
  legH -> AddEntry(hPythia, sLabelP.Data(), "pf");
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    legH -> AddEntry(hEmpty, sLabelE[iTrg].Data(), "");
    legH -> AddEntry(hDataF[iTrg], sLabelF.Data(), "pf");
    legH -> AddEntry(hDataH[iTrg], sLabelH.Data(), "pf");
    legW -> AddEntry(hEmpty, sLabelE[iTrg].Data(), "");
    legW -> AddEntry(hDataF[iTrg], sLabelF.Data(), "pf");
    legW -> AddEntry(hDataH[iTrg], sLabelH.Data(), "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, 0.25};
  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sCol.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sCut.Data());
  cout << "    Made text." << endl;

  // make line
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(2);
  const Float_t fLinXY[NVtx] = {xPlot[0], 1., xPlot[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
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

  TCanvas *cHist = new TCanvas("cHist", "", width, height);
  cHist   -> SetGrid(fGrid, fGrid);
  cHist   -> SetTicks(fTick, fTick);
  cHist   -> SetLogx(fLogX);
  cHist   -> SetLogy(fLogY);
  cHist   -> SetBorderMode(fMode);
  cHist   -> SetBorderSize(fBord);
  cHist   -> SetFrameBorderMode(fFrame);
  cHist   -> SetLeftMargin(fMarginL);
  cHist   -> SetRightMargin(fMarginR);
  cHist   -> SetTopMargin(fMarginT);
  cHist   -> SetBottomMargin(fMarginB);
  cHist   -> cd();
  hPythia -> Draw();
  for (Int_t iTrgPlot = (NTrg - 1); iTrgPlot > -1; iTrgPlot--) {
    hDataF[iTrgPlot] -> Draw("same");
    //hDataH[iTrgPlot] -> Draw("same");
  }
  legH  -> Draw();
  txt   -> Draw();
  fOut  -> cd();
  cHist -> Write();
  cHist -> Close();

  TCanvas *cWeight = new TCanvas("cWeight", "", width, height);
  cWeight            -> SetGrid(fGrid, fGrid);
  cWeight            -> SetTicks(fTick, fTick);
  cWeight            -> SetLogx(fLogX);
  cWeight            -> SetLogy(fLogY);
  cWeight            -> SetBorderMode(fMode);
  cWeight            -> SetBorderSize(fBord);
  cWeight            -> SetFrameBorderMode(fFrame);
  cWeight            -> SetLeftMargin(fMarginL);
  cWeight            -> SetRightMargin(fMarginR);
  cWeight            -> SetTopMargin(fMarginT);
  cWeight            -> SetBottomMargin(fMarginB);
  cWeight            -> cd();
  hWeightF[NTrg - 1] -> Draw();
  //hWeightH[NTrg - 1] -> Draw("same");
  for (Int_t iTrgPlot = (NTrg - 2); iTrgPlot > -1; iTrgPlot--) {
    hWeightF[iTrgPlot] -> Draw("same");
    //hWeightH[iTrgPlot] -> Draw("same");
  }
  line    -> Draw();
  legW    -> Draw();
  txt     -> Draw();
  fOut    -> cd();
  cWeight -> Write();
  cWeight -> Close();

  // save histograms
  fOut    -> cd();
  hPythia -> Write();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hDataF[iTrg]   -> Write();
    hDataH[iTrg]   -> Write();
    hWeightF[iTrg] -> Write();
    hWeightH[iTrg] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fInF[iTrg] -> cd();
    fInF[iTrg] -> Close();
    fInH[iTrg] -> cd();
    fInH[iTrg] -> Close();
  }
  cout << "  Finished trigger weight plot script!\n" << endl;

}

// End ------------------------------------------------------------------------
