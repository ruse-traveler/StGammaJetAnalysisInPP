// 'ModifyStatistics.C'
// Derek Anderson
// 08.07.2021
//
// Takes a (simulated) jet distribution
// and smears it to match the statistics
// of data.

#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NVtx(4);
static const UInt_t NPad(2);
static const UInt_t NHist(2);
static const UInt_t NRange(2);
static const UInt_t NTrgIds(2);
static const UInt_t NTrgBins(3);
static const UInt_t NTrgPi0[NTrgBins] = {12869, 4918, 694};
static const UInt_t NTrgGam[NTrgBins] = {15232, 7328, 1522};

// calculation parameters
static const UInt_t  NMc(1);
static const UInt_t  TrigId(0);
static const UInt_t  TrigBin(0);
static const Float_t RJet(0.5);



void ModifyStatistics() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Begining statistics modification calculation..." << endl;

  // io parameters
  const TString sOut("pp200r9ff.modifiedStatsWithDetCutoff2_detLvlRFF_pTbinHuge.et911r05pi0.d2m3y2022.root");
  const TString sDat("input/data/pp200r9data.forUnfolding_pTbinHuge.et911pt0230dca1vz55.r05a065rm1chrg.root");
  const TString sSim("input/embed/pp200r9embed.forCutoffCheck_d2pM11_pTbinHuge.et911r05qt05130.root");
  const TString sDatJet("Pi0/hJetPtCorrP");
  const TString sSimJet("hSumDetRFF");
  const UInt_t  nTrgSim(6681);

  // text parameters
  const TString sSys("Py6#oplusGeant(RFF), #sqrt{s} = 200 GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // hist parameters
  const TString sDatName("hDataJets");
  const TString sDatLabel("measured data");
  const TString sTitleX("p_{T}^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{3}N^{jet}/d(p_{T}^{jet} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("sim. / data");
  const TString sSimName[NHist]    = {"hSimJetInput", "hSimJetOutput"};
  const TString sRatioName[NHist]  = {"hRatioInput", "hRatioOutput"};
  const TString sSimLabels[NHist]  = {"input det.-lvl. Py6#oplusGeant(FF)", "det.-lvl. Py6#oplusGeant(FF) w/ modified stats"};
  const Float_t xPlotRange[NRange] = {-1., 30.};

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fDat = new TFile(sDat.Data(), "read");
  TFile *fSim = new TFile(sSim.Data(), "read");
  if (!fOut || !fDat || !fSim) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fDat = " << fDat << ", fSim = " << fSim
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input histograms
  TH1D *hDatJet = (TH1D*) fDat -> Get(sDatJet.Data());
  TH1D *hSimJet = (TH1D*) fSim -> Get(sSimJet.Data());
  if (!hDatJet || !hSimJet) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hDatJet = " << hDatJet << ", hSimJet = " << hSimJet
         << endl;
    return;
  }
  hDatJet -> SetName(sDatName.Data());
  hSimJet -> SetName(sSimName[0].Data());
  cout << "    Grabbed input histograms." << endl;

  // parse trigger selection
  UInt_t fTrg(0);
  if (TrigId == 0) {
    fTrg = TrigBin;
  } else {
    fTrg = TrigBin + (NTrgs + 1);
  }

  // select trigger
  UInt_t  nTrg(0);
  TString sTrgEt("");
  TString sTrgId("");
  switch (fTrg) {
    case 0:
      nTrg   = NTrgPi0[0];
      sTrgEt = "E_{T}^{trg} #in (9, 11) GeV";
      sTrgId = "#pi^{0} trig.";
      break;
    case 1:
      nTrg   = NTrgPi0[1];
      sTrgEt = "E_{T}^{trg} #in (11, 15) GeV";
      sTrgId = "#pi^{0} trig.";
      break;
    case 2:
      nTrg   = NTrgPi0[2];
      sTrgEt = "E_{T}^{trg} #in (15, 20) GeV";
      sTrgId = "#pi^{0} trig.";
      break;
    case 3:
      nTrg   = NTrgGam[0];
      sTrgEt = "E_{T}^{trg} #in (9, 11) GeV";
      sTrgId = "#gamma_{rich} trig.";
      break;
    case 4:
      nTrg   = NTrgGam[1];
      sTrgEt = "E_{T}^{trg} #in (11, 15) GeV";
      sTrgId = "#gamma_{rich} trig.";
      break;
    case 5:
      nTrg   = NTrgGam[2];
      sTrgEt = "E_{T}^{trg} #in (15, 20) GeV";
      sTrgId = "#gamma_{rich} trig.";
      break;
    default:
      nTrg   = NTrgPi0[0];
      sTrgEt = "E_{T}^{trg} #in (9, 11) GeV";
      sTrgId = "#pi^{0} trig.";
      break;
  }
  TString sTrg(sTrgId.Data());
  sTrg.Append(", ");
  sTrg.Append(sTrgEt.Data());
  cout << "    Determined trigger selection:\n"
       << "      trigId = " << sTrgId.Data() << ", " << sTrgEt.Data() << "\n"
       << "      nTrgDat = " << nTrg << ", nTrgSim = " << nTrgSim
       << endl;

  // set seed
  TDatime  *seed = new TDatime();
  TRandom3 *rand = new TRandom3();
  rand -> SetSeed(seed -> Get());
  cout << "    Initialized randomizer." << endl;

  // initialize output histograms
  TH1D *hSimOutMC;
  TH1D *hSimOut = (TH1D*) hSimJet -> Clone();
  hSimOut -> SetName(sSimName[1].Data());
  hSimOut -> Reset("ICES");
  cout << "    Beginning MC loop..." << endl;

  // mc loop
  const UInt_t  nBins  = hSimJet -> GetNbinsX();
  const Float_t etaBin = 2. * (1. - RJet);
  for (UInt_t iMC = 0; iMC < NMc; iMC++) {

    // announce iteration
    cout << "      Processing iteration " << iMC + 1 << "/" << NMc << "\r" << flush;
    if (iMC + 1 == NMc) cout << endl;

    // initialize iteration
    TString sSimMC(sSimName[1].Data());
    sSimMC.Append("_mc");
    sSimMC += iMC;

    hSimOutMC = (TH1D*) hSimJet -> Clone();
    hSimOutMC -> SetName(sSimMC.Data());
    hSimOutMC -> Reset("ICES");

    // loop over bins
    for (UInt_t iBin = 1; iBin < nBins + 1; iBin++) {
      const Double_t simVal  = hSimJet -> GetBinContent(iBin);
      const Double_t simBin  = hSimJet -> GetBinWidth(iBin);
      const Double_t newMu   = nTrg * simVal * simBin * etaBin;
      const Double_t newVal  = rand -> PoissonD(newMu);
      const Double_t newNorm = newVal / (simBin * etaBin);
      const Double_t newErr  = TMath::Sqrt(newMu);
      hSimOutMC -> SetBinContent(iBin, newNorm);
      hSimOutMC -> SetBinError(iBin, newErr);
    }
    hSimOutMC -> Scale(1. / (Double_t) nTrg);
    hSimOut   -> Add(hSimOutMC);
  }  // end MC loop

  // average output
  if (NMc > 1) {
    hSimOut -> Scale(1. / (Double_t) NMc);
  }
  cout << "    Finished MC loop." << endl;

  // calculate ratios
  TH1D *hRatioIn  = (TH1D*) hSimJet -> Clone();
  TH1D *hRatioOut = (TH1D*) hSimOut -> Clone();
  hRatioIn  -> SetName(sRatioName[0].Data());
  hRatioOut -> SetName(sRatioName[1].Data());
  hRatioIn  -> Divide(hRatioIn, hDatJet, 1., 1.);
  hRatioOut -> Divide(hRatioOut, hDatJet, 1., 1.);
  cout << "    Calculated ratios." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const UInt_t  fColDat(923);
  const UInt_t  fMarDat(29);
  const UInt_t  fColSim[NHist] = {859, 899};
  const UInt_t  fMarSim[NHist] = {24, 25};
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  const TString sTitle("");
  hDatJet   -> SetMarkerColor(fColDat);
  hDatJet   -> SetMarkerStyle(fMarDat);
  hDatJet   -> SetFillColor(fColDat);
  hDatJet   -> SetFillStyle(fFil);
  hDatJet   -> SetLineColor(fColDat);
  hDatJet   -> SetLineStyle(fLin);
  hDatJet   -> SetLineWidth(fWid);
  hDatJet   -> SetTitle(sTitle.Data());
  hDatJet   -> SetTitleFont(fTxt);
  hDatJet   -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hDatJet   -> GetXaxis() -> SetTitle(sTitleX.Data());
  hDatJet   -> GetXaxis() -> SetTitleFont(fTxt);
  hDatJet   -> GetXaxis() -> SetTitleSize(fTit[1]);
  hDatJet   -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hDatJet   -> GetXaxis() -> SetLabelFont(fTxt);
  hDatJet   -> GetXaxis() -> SetLabelSize(fLab[1]);
  hDatJet   -> GetXaxis() -> CenterTitle(fCnt);
  hDatJet   -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDatJet   -> GetYaxis() -> SetTitleFont(fTxt);
  hDatJet   -> GetYaxis() -> SetTitleSize(fTit[1]);
  hDatJet   -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hDatJet   -> GetYaxis() -> SetLabelFont(fTxt);
  hDatJet   -> GetYaxis() -> SetLabelSize(fLab[1]);
  hSimJet   -> SetMarkerColor(fColSim[0]);
  hSimJet   -> SetMarkerStyle(fMarSim[0]);
  hSimJet   -> SetFillColor(fColSim[0]);
  hSimJet   -> SetFillStyle(fFil);
  hSimJet   -> SetLineColor(fColSim[0]);
  hSimJet   -> SetLineStyle(fLin);
  hSimJet   -> SetLineWidth(fWid);
  hSimJet   -> SetTitle(sTitle.Data());
  hSimJet   -> SetTitleFont(fTxt);
  hSimJet   -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hSimJet   -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSimJet   -> GetXaxis() -> SetTitleFont(fTxt);
  hSimJet   -> GetXaxis() -> SetTitleSize(fTit[1]);
  hSimJet   -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hSimJet   -> GetXaxis() -> SetLabelFont(fTxt);
  hSimJet   -> GetXaxis() -> SetLabelSize(fLab[1]);
  hSimJet   -> GetXaxis() -> CenterTitle(fCnt);
  hSimJet   -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSimJet   -> GetYaxis() -> SetTitleFont(fTxt);
  hSimJet   -> GetYaxis() -> SetTitleSize(fTit[1]);
  hSimJet   -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hSimJet   -> GetYaxis() -> SetLabelFont(fTxt);
  hSimJet   -> GetYaxis() -> SetLabelSize(fLab[1]);
  hSimJet   -> GetYaxis() -> CenterTitle(fCnt);
  hSimJet   -> GetYaxis() -> CenterTitle(fCnt);
  hSimOut   -> SetMarkerColor(fColSim[1]);
  hSimOut   -> SetMarkerStyle(fMarSim[1]);
  hSimOut   -> SetFillColor(fColSim[1]);
  hSimOut   -> SetFillStyle(fFil);
  hSimOut   -> SetLineColor(fColSim[1]);
  hSimOut   -> SetLineStyle(fLin);
  hSimOut   -> SetLineWidth(fWid);
  hSimOut   -> SetTitle(sTitle.Data());
  hSimOut   -> SetTitleFont(fTxt);
  hSimOut   -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hSimOut   -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSimOut   -> GetXaxis() -> SetTitleFont(fTxt);
  hSimOut   -> GetXaxis() -> SetTitleSize(fTit[1]);
  hSimOut   -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hSimOut   -> GetXaxis() -> SetLabelFont(fTxt);
  hSimOut   -> GetXaxis() -> SetLabelSize(fLab[1]);
  hSimOut   -> GetXaxis() -> CenterTitle(fCnt);
  hSimOut   -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSimOut   -> GetYaxis() -> SetTitleFont(fTxt);
  hSimOut   -> GetYaxis() -> SetTitleSize(fTit[1]);
  hSimOut   -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hSimOut   -> GetYaxis() -> SetLabelFont(fTxt);
  hSimOut   -> GetYaxis() -> SetLabelSize(fLab[1]);
  hSimOut   -> GetYaxis() -> CenterTitle(fCnt);
  hSimOut   -> GetYaxis() -> CenterTitle(fCnt);
  hRatioIn  -> SetMarkerColor(fColSim[0]);
  hRatioIn  -> SetMarkerStyle(fMarSim[0]);
  hRatioIn  -> SetFillColor(fColSim[0]);
  hRatioIn  -> SetFillStyle(fFil);
  hRatioIn  -> SetLineColor(fColSim[0]);
  hRatioIn  -> SetLineStyle(fLin);
  hRatioIn  -> SetLineWidth(fWid);
  hRatioIn  -> SetTitle(sTitle.Data());
  hRatioIn  -> SetTitleFont(fTxt);
  hRatioIn  -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hRatioIn  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioIn  -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioIn  -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatioIn  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatioIn  -> GetXaxis() -> SetLabelFont(fTxt);
  hRatioIn  -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatioIn  -> GetXaxis() -> CenterTitle(fCnt);
  hRatioIn  -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatioIn  -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioIn  -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatioIn  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatioIn  -> GetYaxis() -> SetLabelFont(fTxt);
  hRatioIn  -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatioIn  -> GetYaxis() -> CenterTitle(fCnt);
  hRatioIn  -> GetYaxis() -> CenterTitle(fCnt);
  hRatioOut -> SetMarkerColor(fColSim[1]);
  hRatioOut -> SetMarkerStyle(fMarSim[1]);
  hRatioOut -> SetFillColor(fColSim[1]);
  hRatioOut -> SetFillStyle(fFil);
  hRatioOut -> SetLineColor(fColSim[1]);
  hRatioOut -> SetLineStyle(fLin);
  hRatioOut -> SetLineWidth(fWid);
  hRatioOut -> SetTitle(sTitle.Data());
  hRatioOut -> SetTitleFont(fTxt);
  hRatioOut -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hRatioOut -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioOut -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioOut -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatioOut -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatioOut -> GetXaxis() -> SetLabelFont(fTxt);
  hRatioOut -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatioOut -> GetXaxis() -> CenterTitle(fCnt);
  hRatioOut -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatioOut -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioOut -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatioOut -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatioOut -> GetYaxis() -> SetLabelFont(fTxt);
  hRatioOut -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatioOut -> GetYaxis() -> CenterTitle(fCnt);
  hRatioOut -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // make legend
  TString sTrgInfo("N^{trg}_{data} = ");
  sTrgInfo += nTrg;
  sTrgInfo.Append(", N^{trg}_{sim} = ");
  sTrgInfo += nTrgSim;

  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, 0.3};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sTrgInfo.Data());
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hDatJet, sDatLabel.Data());
  leg -> AddEntry(hSimJet, sSimLabels[0].Data());
  leg -> AddEntry(hSimOut, sSimLabels[1].Data());
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, 0.3};
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
  const UInt_t  fColLi(923);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NVtx] = {xPlotRange[0], 1., xPlotRange[1], 1.};
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
  const UInt_t  fLogY1(0);
  const UInt_t  fLogY2(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT1(0.005);
  const Float_t fMarginT2(0.02);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.005);
  const Float_t fPadXY1[NVtx] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NVtx] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot     -> SetGrid(fGrid, fGrid);
  cPlot     -> SetTicks(fTick, fTick);
  cPlot     -> SetBorderMode(fMode);
  cPlot     -> SetBorderSize(fBord);
  pPad1     -> SetGrid(fGrid, fGrid);
  pPad1     -> SetTicks(fTick, fTick);
  pPad1     -> SetLogx(fLogX);
  pPad1     -> SetLogy(fLogY1);
  pPad1     -> SetBorderMode(fMode);
  pPad1     -> SetBorderSize(fBord);
  pPad1     -> SetFrameBorderMode(fFrame);
  pPad1     -> SetLeftMargin(fMarginL);
  pPad1     -> SetRightMargin(fMarginR);
  pPad1     -> SetTopMargin(fMarginT1);
  pPad1     -> SetBottomMargin(fMarginB1);
  pPad2     -> SetGrid(fGrid, fGrid);
  pPad2     -> SetTicks(fTick, fTick);
  pPad2     -> SetLogx(fLogX);
  pPad2     -> SetLogy(fLogY2);
  pPad2     -> SetBorderMode(fMode);
  pPad2     -> SetBorderSize(fBord);
  pPad2     -> SetFrameBorderMode(fFrame);
  pPad2     -> SetLeftMargin(fMarginL);
  pPad2     -> SetRightMargin(fMarginR);
  pPad2     -> SetTopMargin(fMarginT2);
  pPad2     -> SetBottomMargin(fMarginB2);
  cPlot     -> cd();
  pPad1     -> Draw();
  pPad2     -> Draw();
  pPad1     -> cd();
  hRatioIn  -> Draw("E1");
  hRatioOut -> Draw("E1 same");
  line      -> Draw();
  pPad2     -> cd();
  hDatJet   -> Draw("E1");
  hSimJet   -> Draw("E1 same");
  hSimOut   -> Draw("E1 same");
  leg       -> Draw();
  txt       -> Draw();
  fOut      -> cd();
  cPlot     -> Write();
  cPlot     -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut      -> cd();
  hDatJet   -> Write();
  hSimJet   -> Write();
  hSimOut   -> Write();
  hRatioIn  -> Write();
  hRatioOut -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fDat -> cd();
  fDat -> Close();
  fSim -> cd();
  fSim -> Close();
  cout << "  Finished calculation!\n" << endl;

}

// End ------------------------------------------------------------------------
