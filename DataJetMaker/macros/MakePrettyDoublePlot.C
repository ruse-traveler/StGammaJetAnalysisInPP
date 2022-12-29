// 'MakePrettyDoublePlot.C'
// Derek Anderson
// 01.07.2019
//
// Use this to make a pretty plot comparing
// two jet distributions (and do the direct
// -photon if you want)
//
// NOTE: the gamma-subtraction assumes that
//       the histogram specified by 'sInA'
//       is the gamma-rich distribution.


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
static const UInt_t NHist(5);
static const UInt_t NRebin(2);
static const UInt_t NRebinVar(16);
static const Bool_t DoRebin(false);
static const Bool_t DoVariableRebin(false);
static const Bool_t DoGammaRichPlot(false);
static const Bool_t DoGammaSubtraction(false);



void MakePrettyDoublePlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting unfolded distribution..." << endl;

  // io parameters
  const TString sOut("dataXpythiaXembedComp.noPythiaSmearAndHighEff.et911pt0230vz55.r02a005rm1chrg.d7m1y2019.root");
  const TString sInA("pp200r9.dataXpythiaXembedComp.et911pt0230vz55.r02a005rm1chrg.d7m1y2019.root");
  const TString sInB("pp200r9embed.dataXpythiaXembedComp.et911pt0230vz55had.r02a005rm1chrg.d7m1y2018.root");
  const TString sInC("pp200py8det.dataXpythiaXembedComp.noSmearAndHighEff.et911nTrg20Kpt0230pi0.r02a005rm1chrg.d21m1y2019.root");
  const TString sHistA("Pi0/hJetPtCorrP");
  const TString sHistB("hNorm");
  const TString sHistC("Pi0/hJetPtCorrP");

  // plot parameters
  const TString sTitle("");
  const TString sNameA("hData");
  const TString sNameB("hEmbed");
  const TString sNameC("hPythia");
  const TString sNameRB("hRatioEmbed");
  const TString sNameRC("hRatioPythia");
  const TString sTitleX("p_{T}^{reco} = p_{T}^{jet} - #rhoA^{jet} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) dN^{jet}/d(p_{T}^{reco} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleR("ratio");
  const TString sLabelA("data [#pi^{0} trigger]");
  const TString sLabelB("embedding [h^{#pm} trigger]");
  const TString sLabelC("pythia [#pi^{0} trigger, no smear and high eff.]");

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sTrg("E_{T}^{trg} #in (9, 11) GeV");
  const TString sJet("anti-k_{T}, R = 0.2");
  const TString sTyp("#bf{charged jets}");

  // subtraction parameters
  const UInt_t   fFilGR(3017);
  const UInt_t   fLinGR(1);
  const UInt_t   fWidGR(1);
  const UInt_t   fMarGR(4);
  const UInt_t   fColGR(879);
  const UInt_t   fColP[NHist] = {1, 860, 810, 860, 810};
  const UInt_t   fColG[NHist] = {1, 896, 921, 896, 921};
  const TString  sFileGam("pp200r9.binByBinTest.et911vz55gam.r02a005rm1chrg.p0m3k0n58t4.root");
  const TString  sFilePi0("pp200r9.binByBinTest.et911vz55pi0.r02a005rm1chrg.p0m3k0n58t4.root");
  const TString  sHistGam("hUnfolded");
  const TString  sHistPi0("hUnfolded");
  const TString  sNameGR("hGammaRich");
  const TString  sLabelGR("unfolded [#gamma^{rich}]");
  const Double_t gammaPurity(0.650);

  // misc parameters
  const Double_t plotRange[NPlot]    = {-3., 35.};
  const Double_t varRebin[NRebinVar] = {0., 1., 2., 3., 4., 5., 7., 9., 11., 13., 15., 18., 21., 24., 27., 30.};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInA = new TFile(sInA.Data(), "read");
  TFile *fInB = new TFile(sInB.Data(), "read");
  TFile *fInC = new TFile(sInC.Data(), "read");
  if (!fOut || !fInA || !fInB || !fInC) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInA = " << fInA << ", fInB = " << fInB << ", fInC = " << fInC
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hHistA = (TH1D*) fInA -> Get(sHistA.Data());
  TH1D *hHistB = (TH1D*) fInB -> Get(sHistB.Data());
  TH1D *hHistC = (TH1D*) fInC -> Get(sHistC.Data());
  if (!hHistA || !hHistB) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hHistA = " << hHistA << ", hHistB = " << hHistB << ", hHistC = " << hHistC
         << endl;
    return;
  }
  hHistA -> SetName(sNameA.Data());
  hHistB -> SetName(sNameB.Data());
  hHistC -> SetName(sNameC.Data());
  cout << "    Grabbed histograms." << endl;


  // do gamma-subtraction (if need be)
  if (DoGammaSubtraction) {
    TFile *fPi0 = new TFile(sFilePi0.Data(), "read");
    if (!fPi0) {
      cerr << "PANIC: couldn't open pi0 file!" << endl;
      return;
    }
    TH1D  *hPi0 = (TH1D*) fPi0    -> Get(sHistPi0.Data());
    TH1D  *hGam = (TH1D*) hHistA -> Clone();
    if (!hPi0 || !hGam) {
      cerr << "PANIC: couldn't grab pi0 or gamma histogram!\n"
           << "       hPi0 = " << hPi0 << ", hGam = " << hGam
           << endl;
      return;
    }
    hPi0   -> Scale(gammaPurity);
    hHistA -> Add(hGam, hPi0, 1., -1.);
    hHistA -> Scale(1. / (1. - gammaPurity));
    fPi0   -> Close();
    cout << "    Did gamma subtraction." << endl;
  }


  // rebin (if need be)
  hHistA -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  hHistB -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  hHistC -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  if (DoRebin) {
    hHistA -> Rebin(NRebin);
    hHistB -> Rebin(NRebin);
    hHistC -> Rebin(NRebin);

    const UInt_t nBinA = hHistA -> GetNbinsX();
    const UInt_t nBinB = hHistB -> GetNbinsX();
    const UInt_t nBinC = hHistC -> GetNbinsX();
    for (UInt_t iBinA = 1; iBinA <( nBinA + 1); iBinA++) {
      const Double_t valA = hHistA -> GetBinContent(iBinA);
      const Double_t errA = hHistA -> GetBinError(iBinA);
      const Double_t sizA = hHistA -> GetBinWidth(iBinA);
      hHistA -> SetBinContent(iBinA, valA / sizA);
      hHistA -> SetBinError(iBinA, errA / sizA);
    }
    for (UInt_t iBinB = 1; iBinB < (nBinB + 1); iBinB++) {
      const Double_t valB = hHistB -> GetBinContent(iBinB);
      const Double_t errB = hHistB -> GetBinError(iBinB);
      const Double_t sizB = hHistB -> GetBinWidth(iBinB);
      hHistB -> SetBinContent(iBinB, valB / sizB);
      hHistB -> SetBinError(iBinB, errB / sizB);
    }
    for (UInt_t iBinC = 1; iBinC < (nBinC + 1); iBinC++) {
      const Double_t valC = hHistC -> GetBinContent(iBinC);
      const Double_t errC = hHistC -> GetBinError(iBinC);
      const Double_t sizC = hHistC -> GetBinWidth(iBinC);
      hHistC -> SetBinContent(iBinC, valC / sizC);
      hHistC -> SetBinError(iBinC, errC / sizC);
    }
    cout << "    Rebinned histograms (uniform)." << endl;
  }
  if (DoVariableRebin) {
    hHistA = (TH1D*) hHistA -> Rebin(NRebinVar - 1, sNameA.Data(), varRebin);
    hHistB = (TH1D*) hHistB -> Rebin(NRebinVar - 1, sNameB.Data(), varRebin);
    hHistC = (TH1D*) hHistC -> Rebin(NRebinVar - 1, sNameC.Data(), varRebin);

    const UInt_t nBinA = hHistA -> GetNbinsX();
    const UInt_t nBinB = hHistB -> GetNbinsX();
    const UInt_t nBinC = hHistC -> GetNbinsX();
    for (UInt_t iBinA = 1; iBinA < (nBinA + 1); iBinA++) {
      const Double_t valA = hHistA -> GetBinContent(iBinA);
      const Double_t errA = hHistA -> GetBinError(iBinA);
      const Double_t sizA = hHistA -> GetBinWidth(iBinA);
      hHistA -> SetBinContent(iBinA, valA / sizA);
      hHistA -> SetBinError(iBinA, errA / sizA);
    }
    for (UInt_t iBinB = 1; iBinB < (nBinB + 1); iBinB++) {
      const Double_t valB = hHistB -> GetBinContent(iBinB);
      const Double_t errB = hHistB -> GetBinError(iBinB);
      const Double_t sizB = hHistB -> GetBinWidth(iBinB);
      hHistB -> SetBinContent(iBinB, valB / sizB);
      hHistB -> SetBinError(iBinB, errB / sizB);
    }
    for (UInt_t iBinC = 1; iBinC < (nBinC + 1); iBinC++) {
      const Double_t valC = hHistC -> GetBinContent(iBinC);
      const Double_t errC = hHistC -> GetBinError(iBinC);
      const Double_t sizC = hHistC -> GetBinWidth(iBinC);
      hHistC -> SetBinContent(iBinC, valC / sizC);
      hHistC -> SetBinError(iBinC, errC / sizC);
    }
    cout << "    Rebinned histograms (variable)." << endl;
  }

  // calculate ratio
  TH1D *hRatioB = (TH1D*) hHistA -> Clone();
  TH1D *hRatioC = (TH1D*) hHistA -> Clone();
  hRatioB -> SetName(sNameRB.Data());
  hRatioC -> SetName(sNameRC.Data());

  const UInt_t goodDivB = hRatioB  -> Divide(hHistA, hHistB, 1., 1.);
  const UInt_t goodDivC = hRatioC  -> Divide(hHistA, hHistC, 1., 1.);
  const UInt_t nBinR    = hHistA -> GetNbinsX();
  if (goodDivB== 0) {
    for (UInt_t iBinR = 0; iBinR < nBinR; iBinR++) {
      const UInt_t   iBinB = hHistB -> FindBin(hHistA -> GetBinCenter(iBinR));
      const Double_t num   = hHistA -> GetBinContent(iBinR);
      const Double_t den   = hHistB -> GetBinContent(iBinB);
      const Double_t nErr  = hHistA -> GetBinError(iBinR);
      const Double_t dErr  = hHistA -> GetBinError(iBinB);
      const Double_t nRel  = nErr / num;
      const Double_t dRel  = dErr / den;
      const Double_t ratio = num / den;
      const Double_t error = ratio * TMath::Sqrt((nRel * nRel) + (dRel * dRel));
      if (den > 0.) {
        hRatioB -> SetBinContent(iBinR, ratio);
        hRatioB -> SetBinError(iBinR, error);
      }
      else {
        hRatioB -> SetBinContent(iBinR, 0.);
        hRatioB -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
  }
  if (goodDivC== 0) {
    for (UInt_t iBinR = 0; iBinR < nBinR; iBinR++) {
      const UInt_t   iBinC = hHistC -> FindBin(hHistA -> GetBinCenter(iBinR));
      const Double_t num   = hHistA -> GetBinContent(iBinR);
      const Double_t den   = hHistC -> GetBinContent(iBinC);
      const Double_t nErr  = hHistA -> GetBinError(iBinR);
      const Double_t dErr  = hHistA -> GetBinError(iBinC);
      const Double_t nRel  = nErr / num;
      const Double_t dRel  = dErr / den;
      const Double_t ratio = num / den;
      const Double_t error = ratio * TMath::Sqrt((nRel * nRel) + (dRel * dRel));
      if (den > 0.) {
        hRatioC -> SetBinContent(iBinR, ratio);
        hRatioC -> SetBinError(iBinR, error);
      }
      else {
        hRatioC -> SetBinContent(iBinR, 0.);
        hRatioC -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
  }
  hRatioB -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  hRatioC -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  cout << "    Calculated ratio." << endl;


  // set styles
  UInt_t fCol[NHist];
  if (DoGammaSubtraction) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      fCol[iHist] = fColG[iHist];
    }
  }
  else {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      fCol[iHist] = fColP[iHist];
    }
  }

  const UInt_t  fMar[NHist] = {29, 24, 25, 24, 25};
  const UInt_t  fFil[NHist] = {0, 0, 0, 0, 0};
  const UInt_t  fLin[NHist] = {1, 1, 1, 1, 1};
  const UInt_t  fWid[NHist] = {1, 1, 1, 1, 1};
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hHistA  -> SetMarkerColor(fCol[0]);
  hHistA  -> SetMarkerStyle(fMar[0]);
  hHistA  -> SetFillColor(fCol[0]);
  hHistA  -> SetFillStyle(fFil[0]);
  hHistA  -> SetLineColor(fCol[0]);
  hHistA  -> SetLineStyle(fLin[0]);
  hHistA  -> SetLineWidth(fWid[0]);
  hHistA  -> SetTitle(sTitle.Data());
  hHistA  -> SetTitleFont(fTxt);
  hHistA  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHistA  -> GetXaxis() -> SetTitleFont(fTxt);
  hHistA  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hHistA  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hHistA  -> GetXaxis() -> SetLabelFont(fTxt);
  hHistA  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hHistA  -> GetXaxis() -> CenterTitle(fCnt);
  hHistA  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hHistA  -> GetYaxis() -> SetTitleFont(fTxt);
  hHistA  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hHistA  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hHistA  -> GetYaxis() -> SetLabelFont(fTxt);
  hHistA  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hHistA  -> GetYaxis() -> CenterTitle(fCnt);
  hHistB  -> SetMarkerColor(fCol[1]);
  hHistB  -> SetMarkerStyle(fMar[1]);
  hHistB  -> SetFillColor(fCol[1]);
  hHistB  -> SetFillStyle(fFil[1]);
  hHistB  -> SetLineColor(fCol[1]);
  hHistB  -> SetLineStyle(fLin[1]);
  hHistB  -> SetLineWidth(fWid[1]);
  hHistB  -> SetTitle(sTitle.Data());
  hHistB  -> SetTitleFont(fTxt);
  hHistB  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHistB  -> GetXaxis() -> SetTitleFont(fTxt);
  hHistB  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hHistB  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hHistB  -> GetXaxis() -> SetLabelFont(fTxt);
  hHistB  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hHistB  -> GetXaxis() -> CenterTitle(fCnt);
  hHistB  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hHistB  -> GetYaxis() -> SetTitleFont(fTxt);
  hHistB  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hHistB  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hHistB  -> GetYaxis() -> SetLabelFont(fTxt);
  hHistB  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hHistB  -> GetYaxis() -> CenterTitle(fCnt);
  hHistC  -> SetMarkerColor(fCol[2]);
  hHistC  -> SetMarkerStyle(fMar[2]);
  hHistC  -> SetFillColor(fCol[2]);
  hHistC  -> SetFillStyle(fFil[2]);
  hHistC  -> SetLineColor(fCol[2]);
  hHistC  -> SetLineStyle(fLin[2]);
  hHistC  -> SetLineWidth(fWid[2]);
  hHistC  -> SetTitle(sTitle.Data());
  hHistC  -> SetTitleFont(fTxt);
  hHistC  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hHistC  -> GetXaxis() -> SetTitleFont(fTxt);
  hHistC  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hHistC  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hHistC  -> GetXaxis() -> SetLabelFont(fTxt);
  hHistC  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hHistC  -> GetXaxis() -> CenterTitle(fCnt);
  hHistC  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hHistC  -> GetYaxis() -> SetTitleFont(fTxt);
  hHistC  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hHistC  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hHistC  -> GetYaxis() -> SetLabelFont(fTxt);
  hHistC  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hHistC  -> GetYaxis() -> CenterTitle(fCnt);
  hRatioB -> SetMarkerColor(fCol[3]);
  hRatioB -> SetMarkerStyle(fMar[3]);
  hRatioB -> SetFillColor(fCol[3]);
  hRatioB -> SetFillStyle(fFil[3]);
  hRatioB -> SetLineColor(fCol[3]);
  hRatioB -> SetLineStyle(fLin[3]);
  hRatioB -> SetLineWidth(fWid[3]);
  hRatioB -> SetTitle(sTitle.Data());
  hRatioB -> SetTitleFont(fTxt);
  hRatioB -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioB -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioB -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatioB -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatioB -> GetXaxis() -> SetLabelFont(fTxt);
  hRatioB -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatioB -> GetXaxis() -> CenterTitle(fCnt);
  hRatioB -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatioB -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioB -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatioB -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatioB -> GetYaxis() -> SetLabelFont(fTxt);
  hRatioB -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatioB -> GetYaxis() -> CenterTitle(fCnt);
  hRatioC -> SetMarkerColor(fCol[4]);
  hRatioC -> SetMarkerStyle(fMar[4]);
  hRatioC -> SetFillColor(fCol[4]);
  hRatioC -> SetFillStyle(fFil[4]);
  hRatioC -> SetLineColor(fCol[4]);
  hRatioC -> SetLineStyle(fLin[4]);
  hRatioC -> SetLineWidth(fWid[4]);
  hRatioC -> SetTitle(sTitle.Data());
  hRatioC -> SetTitleFont(fTxt);
  hRatioC -> GetXaxis() -> SetTitle(sTitleX.Data());
  hRatioC -> GetXaxis() -> SetTitleFont(fTxt);
  hRatioC -> GetXaxis() -> SetTitleSize(fTit[0]);
  hRatioC -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hRatioC -> GetXaxis() -> SetLabelFont(fTxt);
  hRatioC -> GetXaxis() -> SetLabelSize(fLab[0]);
  hRatioC -> GetXaxis() -> CenterTitle(fCnt);
  hRatioC -> GetYaxis() -> SetTitle(sTitleR.Data());
  hRatioC -> GetYaxis() -> SetTitleFont(fTxt);
  hRatioC -> GetYaxis() -> SetTitleSize(fTit[0]);
  hRatioC -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hRatioC -> GetYaxis() -> SetLabelFont(fTxt);
  hRatioC -> GetYaxis() -> SetLabelSize(fLab[0]);
  hRatioC -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t fLegXY[NPlot * NPlot] = {0.1, 0.1, 0.3, 0.3};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  leg -> AddEntry(hHistA, sLabelA.Data());
  leg -> AddEntry(hHistB, sLabelB.Data());
  leg -> AddEntry(hHistC, sLabelC.Data());
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
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
  const Float_t fLinXY[NPlot * NPlot] = {plotRange[0], 1., plotRange[1], 1.};
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
  const Float_t fMarginR(0.05);
  const Float_t fMarginT1(0.);
  const Float_t fMarginT2(0.05);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.);
  const Float_t fPadXY1[NPlot * NPlot] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NPlot * NPlot] = {0., 0.35, 1., 1.};

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  TPad    *pPad1 = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2 = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot   -> SetGrid(fGrid, fGrid);
  cPlot   -> SetTicks(fTick, fTick);
  cPlot   -> SetBorderMode(fMode);
  cPlot   -> SetBorderSize(fBord);
  pPad1   -> SetGrid(fGrid, fGrid);
  pPad1   -> SetTicks(fTick, fTick);
  pPad1   -> SetLogx(fLogX);
  pPad1   -> SetLogy(fLogY);
  pPad1   -> SetBorderMode(fMode);
  pPad1   -> SetBorderSize(fBord);
  pPad1   -> SetFrameBorderMode(fFrame);
  pPad1   -> SetLeftMargin(fMarginL);
  pPad1   -> SetRightMargin(fMarginR);
  pPad1   -> SetTopMargin(fMarginT1);
  pPad1   -> SetBottomMargin(fMarginB1);
  pPad2   -> SetGrid(fGrid, fGrid);
  pPad2   -> SetTicks(fTick, fTick);
  pPad2   -> SetLogx(fLogX);
  pPad2   -> SetLogy(fLogY);
  pPad2   -> SetBorderMode(fMode);
  pPad2   -> SetBorderSize(fBord);
  pPad2   -> SetFrameBorderMode(fFrame);
  pPad2   -> SetLeftMargin(fMarginL);
  pPad2   -> SetRightMargin(fMarginR);
  pPad2   -> SetTopMargin(fMarginT2);
  pPad2   -> SetBottomMargin(fMarginB2);
  cPlot   -> cd();
  pPad1   -> Draw();
  pPad2   -> Draw();
  pPad1   -> cd();
  hRatioB -> Draw("E2");
  hRatioC -> Draw("SAME E2");
  line    -> Draw();
  pPad2   -> cd();
  hHistA  -> Draw("E2");
  hHistB  -> Draw("SAME E2");
  hHistC  -> Draw("SAME E2");
  leg     -> Draw();
  txt     -> Draw();
  fOut    -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;


  // close files
  fOut    -> cd();
  hHistA  -> Write();
  hHistB  -> Write();
  hHistC  -> Write();
  hRatioB -> Write();
  hRatioC -> Write();
  fOut    -> Close();
  fInA    -> cd();
  fInA    -> Close();
  fInB    -> cd();
  fInB    -> Close();
  cout << "  Plot made!\n" << endl;

}

// End ------------------------------------------------------------------------
