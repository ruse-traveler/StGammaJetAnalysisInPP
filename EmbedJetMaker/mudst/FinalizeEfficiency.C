// 'FinalizeEfficiency.C'
// Derek Anderson
// 05.15.2019
//
// Use this to sum the relevant particle-
// and detector-level histograms with a
// given set of weights from embedding.
// From the sums, it calculates the
// corresponding efficiency.
//
// NOTE: Configuration = 0, RFF configuration
//       Configuration = 1, FF configuration
//       Configuration = 2, both (1st RFF, then FF)


#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// constants
static const UInt_t   NHistRFF(8);
static const UInt_t   NHistFF(7);
static const UInt_t   NTotal(10);
static const UInt_t   NFunc(5);
static const UInt_t   NLvl(3);
static const UInt_t   NHistAll(NHistRFF + NHistFF);
static const Double_t NTrgsFF[NHistFF]   = {2., 37., 163., 463., 2419., 8762., 4455.};
static const Double_t NTrgsRFF[NHistRFF] = {1., 5., 36., 204., 628., 2807., 7547., 3855.};
static const Double_t NTrgsAll[NHistAll] = {1., 5., 36., 204., 628., 2807., 7547., 3855., 2., 37., 163., 463., 2419., 8762., 4455.};
static const Double_t WeightsFF[NTotal]  = {1.0, 3.361596e-01, 1.401161e-01, 1.337302e-01, 2.895246e-02, 1.042577e-02, 8.294575e-03, 2.064352e-03, 8.088693e-05, 1.417116e-05};
static const Double_t WeightsRFF[NTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};

// options
static const UInt_t NPav(4);
static const UInt_t NHist(NHistAll);
static const UInt_t Configuration(2);
static const Bool_t DoIntNorm(false);
static const Bool_t UseRootNTrgs(true);



void FinalizeEfficiency() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString sOut("pp200r9embed.effFitTest_pTdetOverPtDet_onlyPtMc0ForMatch.et1220pt0230dca3vz55pi0.d11m5y2020.root");
  const TString sIn[NHist]   = {"pp200r9pt4rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt5rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt7rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt9rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt11rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt15rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt25rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt35rff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt5ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt7ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt9ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt11ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt15ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt25ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root", "pp200r9pt35ff.pTmc0.et1220pt0230dca3vz55pi0.d10m5y2020.root"};
  const TString sPar[NHist]  = {"particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par", "particle/hPtAfterQA_par"};
  const TString sDet[NHist]  = {"detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det", "detector/hPtAfterQA_det"};
  const TString sNorm[NHist] = {"particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par", "particle/hNumTrk_par"};

  // histogram parameters
  const TString sNameSum[NLvl]  = {"hParticleSum", "hDetectorSum", "hRatioSum"};
  const TString sNameNorm[NLvl] = {"hParticleNorm", "hDetectorNorm", "hRatioNorm"};
  const TString sTitle[NLvl]    = {"", "", ""};
  const TString sTitleX[NLvl]   = {"p_{T}^{trk} [GeV/c]", "p_{T}^{trk} [GeV/c]", "p_{T}^{trk} [GeV/c]"};
  const TString sTitleYS[NLvl]  = {"dN^{trk}/dp_{T}^{trk} [GeV/c]^{-1}", "dN^{trk}/dp_{T}^{trk} [GeV/c]^{-1}", "detector / particle"};
  const TString sTitleYN[NLvl]  = {"(1/N^{trg}) dN^{trk}/dp_{T}^{trk} [GeV/c]^{-1}", "(1/N^{trg}) dN^{trk}/dp_{T}^{trk} [GeV/c]^{-1}", "detector / particle"};

  // fit parameters
  const TString  sPavF("f(p_{T}^{trk}) = #epsilon_{0} + #epsilon_{1}e^{-#sigma_{1}p_{T}^{trk}} + #epsilon_{2}e^{-#sigma_{2}(p_{T}^{trk})^{2}}");
  const TString  sFunc("[0]+[1]*TMath::Exp(-1.*[2]*x)+[3]*TMath::Exp(-1.*[4]*x*x)");
  const TString  sSumF("fRatioSum");
  const TString  sNormF("fRatioNorm");
  const Double_t fSysErr(0.04);
  const Double_t fFitMin(0.2);
  const Double_t fFitMax(20.); 
  const Double_t parmGuess[NFunc] = {0.90, -5.0, 15., 0.10, 0.25};
  const TString  sNameParm[NFunc] = {"#epsilon_{0}", "#epsilon_{1}", "#sigma_{1}", "#epsilon_{2}", "#sigma_{2}"};

  // plotting parameters
  const UInt_t  fColF(879);
  const UInt_t  fLinF(1);
  const UInt_t  fWidF(2);
  const UInt_t  fCol[NLvl] = {859, 899, 923};
  const UInt_t  fMar[NLvl] = {20, 24, 20};
  const UInt_t  fLin[NLvl] = {1, 1, 1};
  const UInt_t  fFil[NLvl] = {0, 0, 0};
  const TString sLeg[NLvl] = {"particle level (after QA cuts)", "detector level (no matching, after QA cuts)", "detector / particle"};
  const TString sPav[NPav] = {"pp-collisions, #sqrt{s} = 200 GeV", "real #pi^{0} triggers, E_{T}^{trg} #in (9, 20) GeV/c", "using only recoil tracks", ""};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] = new TFile(sIn[iHist].Data(), "read");
    if (!fIn[iHist]) {
      cerr << "PANIC: couldn't open input file no. " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH1D *hPar[NHist];
  TH1D *hDet[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hPar[iHist] = (TH1D*) fIn[iHist] -> Get(sPar[iHist].Data());
    hDet[iHist] = (TH1D*) fIn[iHist] -> Get(sDet[iHist].Data());
    if (!hPar[iHist]) {
      cerr << "PANIC: couldn't grab particle-level histogram no. " << iHist << "!" << endl;
      return;
    }
    if (!hDet[iHist]) {
      cerr << "PANIC: couldn't grab detector-level histogram no. " << iHist << "!" << endl;
      return;
    }
  }
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  Double_t NormsFF[NHistFF];
  Double_t NormsRFF[NHistRFF];
  Double_t NormsAll[NHistAll];
  TH1D     *hNormFF[NHistFF];
  TH1D     *hNormRFF[NHistRFF];
  TH1D     *hNormAll[NHistAll];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        if (UseRootNTrgs) {
          hNormRFF[iHist] = (TH1D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsRFF[iHist] = hNormRFF[iHist] -> GetEntries();
        }
        else
          NormsRFF[iHist] = NTrgsRFF[iHist];
        break;
      case 1:
        if (UseRootNTrgs) {
          hNormFF[iHist] = (TH1D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsFF[iHist] = hNormFF[iHist] -> GetEntries();
        }
        else
          NormsFF[iHist] = NTrgsFF[iHist];
        break;
      case 2:
        if (UseRootNTrgs) {
          hNormAll[iHist] = (TH1D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsAll[iHist] = hNormAll[iHist] -> GetEntries();
        }
        else
          NormsAll[iHist] = NTrgsAll[iHist];
        break;
    }  // end switch
  }  // end histogram loop
  cout << "    Extracted no. of triggers." << endl;


  // scale histograms
  Double_t normer(0.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        hPar[iHist] -> Scale(WeightsRFF[(NTotal - NHist) + iHist]);
        hDet[iHist] -> Scale(WeightsRFF[(NTotal - NHist) + iHist]);
        normer += NormsRFF[iHist] * WeightsRFF[(NTotal - NHist) + iHist];
        break;
      case 1:
        hPar[iHist] -> Scale(WeightsFF[(NTotal - NHist) + iHist]);
        hDet[iHist] -> Scale(WeightsFF[(NTotal - NHist) + iHist]);
        normer += NormsFF[iHist] * WeightsFF[(NTotal - NHist) + iHist];
        break;
      case 2:
        const Bool_t isRFF = (iHist < NHistRFF);
        if (isRFF) {
          hPar[iHist] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHist]);
          hDet[iHist] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHist]);
          normer += NormsAll[iHist] * WeightsRFF[(NTotal - NHistRFF) + iHist];
        }
        else {
          hPar[iHist] -> Scale(WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)]);
          hDet[iHist] -> Scale(WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)]);
          normer += NormsAll[iHist] * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
        }
        break;
    }
  }
  cout << "    Scaled histograms." << endl;

  // sum histograms
  TH1D *hSumP  = (TH1D*) hPar[0] -> Clone();
  TH1D *hSumD  = (TH1D*) hDet[0] -> Clone();
  TH1D *hNormP = (TH1D*) hPar[0] -> Clone();
  TH1D *hNormD = (TH1D*) hDet[0] -> Clone();
  hSumP  -> SetName(sNameSum[0].Data());
  hSumD  -> SetName(sNameSum[1].Data());
  hNormP -> SetName(sNameNorm[0].Data());
  hNormD -> SetName(sNameNorm[1].Data());
  hSumP  -> Reset("ICE");
  hSumD  -> Reset("ICE");
  hNormP -> Reset("ICE");
  hNormD -> Reset("ICE");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hSumP  -> Add(hPar[iHist]);
    hSumD  -> Add(hDet[iHist]);
    hNormP -> Add(hPar[iHist]);
    hNormD -> Add(hDet[iHist]);
  }

  // normalize sum
  if (DoIntNorm) {
    const Double_t integralP = hNormP -> Integral();
    const Double_t integralD = hNormD -> Integral();
    hNormP -> Scale(1. / integralP);
    hNormD -> Scale(1. / integralD);
    cout << "    Normalized histograms.\n"
         << "      normalization(par) = " << integralP << "\n"
         << "      normalization(det) = " << integralD
         << endl;
  }
  else {
    hNormP -> Scale(1. / normer);
    hNormD -> Scale(1. / normer);
    cout << "    Normalized histograms.\n"
         << "      normalization = " << normer
         << endl;
  }


  // calculate ratios
  const Double_t weightP(1.);
  const Double_t weightD(1.);

  TH1D *hSumR  = (TH1D*) hSumP  -> Clone();
  TH1D *hNormR = (TH1D*) hNormP -> Clone();
  hSumR  -> Reset("ICE");
  hNormR -> Reset("ICE");
  hSumR  -> SetName(sNameSum[2].Data());
  hNormR -> SetName(sNameNorm[2].Data());
  hSumR  -> Divide(hSumD, hSumP, weightD, weightP);
  hNormR -> Divide(hNormD, hNormP, weightD, weightP);
  cout << "    Calculated ratios." << endl;


  // fit ratios
  TF1 *fSumR  = new TF1(sSumF.Data(), sFunc.Data(), fFitMin, fFitMax);
  TF1 *fNormR = new TF1(sNormF.Data(), sFunc.Data(), fFitMin, fFitMax);
  for (UInt_t iPar = 0; iPar < NFunc; iPar++) {
    fSumR  -> SetParName(iPar, sNameParm[iPar].Data());
    fNormR -> SetParName(iPar, sNameParm[iPar].Data());
    fSumR  -> SetParameter(iPar, parmGuess[iPar]);
    fNormR -> SetParameter(iPar, parmGuess[iPar]);
  }
  fSumR  -> SetLineColor(fColF);
  fNormR -> SetLineColor(fColF);
  fSumR  -> SetLineStyle(fLinF);
  fNormR -> SetLineStyle(fLinF);
  fSumR  -> SetLineWidth(fWidF);
  fNormR -> SetLineWidth(fWidF);
  hSumR  -> Fit(fSumR, "R");
  hNormR -> Fit(fNormR, "R");
  cout << "    Fit ratios." << endl;


  // set styles
  const UInt_t  fCnt(1);
  const UInt_t  fTxt(42);
  const Float_t fTtl(0.04);
  const Float_t fLbl(0.03);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.2);
  // particle histograms
  hSumP  -> SetMarkerColor(fCol[0]);
  hSumP  -> SetMarkerStyle(fMar[0]);
  hSumP  -> SetLineColor(fCol[0]);
  hSumP  -> SetLineStyle(fLin[0]);
  hSumP  -> SetFillColor(fCol[0]);
  hSumP  -> SetFillStyle(fFil[0]);
  hSumP  -> SetTitle(sTitle[0]);
  hSumP  -> SetTitleFont(fTxt);
  hSumP  -> GetXaxis() -> SetLabelFont(fTxt);
  hSumP  -> GetXaxis() -> SetLabelSize(fLbl);
  hSumP  -> GetXaxis() -> SetTitle(sTitleX[0].Data());
  hSumP  -> GetXaxis() -> SetTitleFont(fTxt);
  hSumP  -> GetXaxis() -> SetTitleSize(fTtl);
  hSumP  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSumP  -> GetXaxis() -> CenterTitle(fCnt);
  hSumP  -> GetYaxis() -> SetLabelFont(fTxt);
  hSumP  -> GetYaxis() -> SetLabelSize(fLbl);
  hSumP  -> GetYaxis() -> SetTitle(sTitleYS[0].Data());
  hSumP  -> GetYaxis() -> SetTitleFont(fTxt);
  hSumP  -> GetYaxis() -> SetTitleSize(fTtl);
  hSumP  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSumP  -> GetYaxis() -> CenterTitle(fCnt);
  hNormP -> SetMarkerColor(fCol[0]);
  hNormP -> SetMarkerStyle(fMar[0]);
  hNormP -> SetLineColor(fCol[0]);
  hNormP -> SetLineStyle(fLin[0]);
  hNormP -> SetFillColor(fCol[0]);
  hNormP -> SetFillStyle(fFil[0]);
  hNormP -> SetTitle(sTitle[0]);
  hNormP -> SetTitleFont(fTxt);
  hNormP -> GetXaxis() -> SetLabelFont(fTxt);
  hNormP -> GetXaxis() -> SetLabelSize(fLbl);
  hNormP -> GetXaxis() -> SetTitle(sTitleX[0].Data());
  hNormP -> GetXaxis() -> SetTitleFont(fTxt);
  hNormP -> GetXaxis() -> SetTitleSize(fTtl);
  hNormP -> GetXaxis() -> SetTitleOffset(fOffX);
  hNormP -> GetXaxis() -> CenterTitle(fCnt);
  hNormP -> GetYaxis() -> SetLabelFont(fTxt);
  hNormP -> GetYaxis() -> SetLabelSize(fLbl);
  hNormP -> GetYaxis() -> SetTitle(sTitleYN[0].Data());
  hNormP -> GetYaxis() -> SetTitleFont(fTxt);
  hNormP -> GetYaxis() -> SetTitleSize(fTtl);
  hNormP -> GetYaxis() -> SetTitleOffset(fOffY);
  hNormP -> GetYaxis() -> CenterTitle(fCnt);
  // detector histograms
  hSumD  -> SetMarkerColor(fCol[1]);
  hSumD  -> SetMarkerStyle(fMar[1]);
  hSumD  -> SetLineColor(fCol[1]);
  hSumD  -> SetLineStyle(fLin[1]);
  hSumD  -> SetFillColor(fCol[1]);
  hSumD  -> SetFillStyle(fFil[1]);
  hSumD  -> SetTitle(sTitle[1]);
  hSumD  -> SetTitleFont(fTxt);
  hSumD  -> GetXaxis() -> SetLabelFont(fTxt);
  hSumD  -> GetXaxis() -> SetLabelSize(fLbl);
  hSumD  -> GetXaxis() -> SetTitle(sTitleX[1].Data());
  hSumD  -> GetXaxis() -> SetTitleFont(fTxt);
  hSumD  -> GetXaxis() -> SetTitleSize(fTtl);
  hSumD  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSumD  -> GetXaxis() -> CenterTitle(fCnt);
  hSumD  -> GetYaxis() -> SetLabelFont(fTxt);
  hSumD  -> GetYaxis() -> SetLabelSize(fLbl);
  hSumD  -> GetYaxis() -> SetTitle(sTitleYS[1].Data());
  hSumD  -> GetYaxis() -> SetTitleFont(fTxt);
  hSumD  -> GetYaxis() -> SetTitleSize(fTtl);
  hSumD  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSumD  -> GetYaxis() -> CenterTitle(fCnt);
  hNormD -> SetMarkerColor(fCol[1]);
  hNormD -> SetMarkerStyle(fMar[1]);
  hNormD -> SetLineColor(fCol[1]);
  hNormD -> SetLineStyle(fLin[1]);
  hNormD -> SetFillColor(fCol[1]);
  hNormD -> SetFillStyle(fFil[1]);
  hNormD -> SetTitle(sTitle[1]);
  hNormD -> SetTitleFont(fTxt);
  hNormD -> GetXaxis() -> SetLabelFont(fTxt);
  hNormD -> GetXaxis() -> SetLabelSize(fLbl);
  hNormD -> GetXaxis() -> SetTitle(sTitleX[1].Data());
  hNormD -> GetXaxis() -> SetTitleFont(fTxt);
  hNormD -> GetXaxis() -> SetTitleSize(fTtl);
  hNormD -> GetXaxis() -> SetTitleOffset(fOffX);
  hNormD -> GetXaxis() -> CenterTitle(fCnt);
  hNormD -> GetYaxis() -> SetLabelFont(fTxt);
  hNormD -> GetYaxis() -> SetLabelSize(fLbl);
  hNormD -> GetYaxis() -> SetTitle(sTitleYN[1].Data());
  hNormD -> GetYaxis() -> SetTitleFont(fTxt);
  hNormD -> GetYaxis() -> SetTitleSize(fTtl);
  hNormD -> GetYaxis() -> SetTitleOffset(fOffY);
  hNormD -> GetYaxis() -> CenterTitle(fCnt);
  // ratio histograms
  hSumR  -> SetMarkerColor(fCol[2]);
  hSumR  -> SetMarkerStyle(fMar[2]);
  hSumR  -> SetLineColor(fCol[2]);
  hSumR  -> SetLineStyle(fLin[2]);
  hSumR  -> SetFillColor(fCol[2]);
  hSumR  -> SetFillStyle(fFil[2]);
  hSumR  -> SetTitle(sTitle[2]);
  hSumR  -> SetTitleFont(fTxt);
  hSumR  -> GetXaxis() -> SetLabelFont(fTxt);
  hSumR  -> GetXaxis() -> SetLabelSize(fLbl);
  hSumR  -> GetXaxis() -> SetTitle(sTitleX[2].Data());
  hSumR  -> GetXaxis() -> SetTitleFont(fTxt);
  hSumR  -> GetXaxis() -> SetTitleSize(fTtl);
  hSumR  -> GetXaxis() -> SetTitleOffset(fOffX);
  hSumR  -> GetXaxis() -> CenterTitle(fCnt);
  hSumR  -> GetYaxis() -> SetLabelFont(fTxt);
  hSumR  -> GetYaxis() -> SetLabelSize(fLbl);
  hSumR  -> GetYaxis() -> SetTitle(sTitleYS[2].Data());
  hSumR  -> GetYaxis() -> SetTitleFont(fTxt);
  hSumR  -> GetYaxis() -> SetTitleSize(fTtl);
  hSumR  -> GetYaxis() -> SetTitleOffset(fOffY);
  hSumR  -> GetYaxis() -> CenterTitle(fCnt);
  hNormR -> SetMarkerColor(fCol[2]);
  hNormR -> SetMarkerStyle(fMar[2]);
  hNormR -> SetLineColor(fCol[2]);
  hNormR -> SetLineStyle(fLin[2]);
  hNormR -> SetFillColor(fCol[2]);
  hNormR -> SetFillStyle(fFil[2]);
  hNormR -> SetTitle(sTitle[2]);
  hNormR -> SetTitleFont(fTxt);
  hNormR -> GetXaxis() -> SetLabelFont(fTxt);
  hNormR -> GetXaxis() -> SetLabelSize(fLbl);
  hNormR -> GetXaxis() -> SetTitle(sTitleX[2].Data());
  hNormR -> GetXaxis() -> SetTitleFont(fTxt);
  hNormR -> GetXaxis() -> SetTitleSize(fTtl);
  hNormR -> GetXaxis() -> SetTitleOffset(fOffX);
  hNormR -> GetXaxis() -> CenterTitle(fCnt);
  hNormR -> GetYaxis() -> SetLabelFont(fTxt);
  hNormR -> GetYaxis() -> SetLabelSize(fLbl);
  hNormR -> GetYaxis() -> SetTitle(sTitleYN[2].Data());
  hNormR -> GetYaxis() -> SetTitleFont(fTxt);
  hNormR -> GetYaxis() -> SetTitleSize(fTtl);
  hNormR -> GetYaxis() -> SetTitleOffset(fOffY);
  hNormR -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;


  // make legends
  const UInt_t  fColL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fAln(12);
  const Float_t legX[2] = {0.7, 0.9};
  const Float_t legY[2] = {0.7, 0.9};

  TLegend *lSum  = new TLegend(legX[0], legY[0], legX[1], legY[1]);
  TLegend *lNorm = new TLegend(legX[0], legY[0], legX[1], legY[1]);
  lSum  -> SetFillColor(fColL);
  lSum  -> SetFillStyle(fFilL);
  lSum  -> SetLineColor(fColL);
  lSum  -> SetLineStyle(fLinL);
  lSum  -> SetTextFont(fTxt);
  lSum  -> SetTextAlign(fAln);
  lSum  -> AddEntry(hSumP, sLeg[0].Data());
  lSum  -> AddEntry(hSumD, sLeg[1].Data());
  lSum  -> AddEntry(hSumR, sLeg[2].Data());
  lNorm -> SetFillColor(fColL);
  lNorm -> SetFillStyle(fFilL);
  lNorm -> SetLineColor(fColL);
  lNorm -> SetLineStyle(fLinL);
  lNorm -> SetTextFont(fTxt);
  lNorm -> SetTextAlign(fAln);
  lNorm -> AddEntry(hNormP, sLeg[0].Data());
  lNorm -> AddEntry(hNormD, sLeg[1].Data());
  lNorm -> AddEntry(hNormR, sLeg[2].Data());
  cout << "    Made legends." << endl;


  // make text boxes
  const UInt_t  fColP(0);
  const UInt_t  fLinP(0);
  const UInt_t  fFilP(0);
  const Float_t pavX1[2] = {0.7, 0.9};
  const Float_t pavX2[2] = {0.5, 0.7};
  const Float_t pavY1[2] = {0.1, 0.3};
  const Float_t pavY2[2] = {0.1, 0.3};

  TPaveText *ptSum  = new TPaveText(pavX1[0], pavY1[0], pavX1[1], pavY1[1], "NDC NB");
  TPaveText *ptNorm = new TPaveText(pavX1[0], pavY1[0], pavX1[1], pavY1[1], "NDC NB");
  ptSum  -> SetFillColor(fColP);
  ptSum  -> SetFillStyle(fFilP);
  ptSum  -> SetLineColor(fColP);
  ptSum  -> SetLineStyle(fLinP);
  ptSum  -> SetTextFont(fTxt);
  ptSum  -> SetTextAlign(fAln);
  ptNorm -> SetFillColor(fColP);
  ptNorm -> SetFillStyle(fFilP);
  ptNorm -> SetLineColor(fColP);
  ptNorm -> SetLineStyle(fLinP);
  ptNorm -> SetTextFont(fTxt);
  ptNorm -> SetTextAlign(fAln);
  for (UInt_t iPav = 0; iPav < NPav; iPav++) {
    ptSum  -> AddText(sPav[iPav].Data());
    ptNorm -> AddText(sPav[iPav].Data());
  }

  TPaveText *ptSumF  = new TPaveText(pavX2[0], pavY2[0], pavX2[1], pavY2[1], "NDC NB");
  TPaveText *ptNormF = new TPaveText(pavX2[0], pavY2[0], pavX2[1], pavY2[1], "NDC NB");
  ptSumF  -> SetFillColor(fColP);
  ptSumF  -> SetFillStyle(fFilP);
  ptSumF  -> SetLineColor(fColP);
  ptSumF  -> SetLineStyle(fLinP);
  ptSumF  -> SetTextFont(fTxt);
  ptSumF  -> SetTextAlign(fAln);
  ptNormF -> SetFillColor(fColP);
  ptNormF -> SetFillStyle(fFilP);
  ptNormF -> SetLineColor(fColP);
  ptNormF -> SetLineStyle(fLinP);
  ptNormF -> SetTextFont(fTxt);
  ptNormF -> SetTextAlign(fAln);

  TString sPavWithColF("#color[");
  sPavWithColF += fColF;
  sPavWithColF += "]{";
  sPavWithColF += sPavF.Data();
  sPavWithColF += "}";

  ptSumF  -> AddText(sPavWithColF.Data());
  ptSumF  -> AddText("TEST TEST TEST");
  ptNormF -> AddText(sPavWithColF.Data());
  ptNormF -> AddText("TEST TEST TEST");
  cout << "    Made text boxes." << endl;



  // make lines
  const UInt_t  fColL(1);
  const UInt_t  fLinL(2);
  const UInt_t  nBinRS    = hSumR  -> GetNbinsX();
  const UInt_t  nBinRN    = hNormR -> GetNbinsX();
  const Float_t xStartRS  = hSumR  -> GetBinLowEdge(1);
  const Float_t xStartRN  = hNormR -> GetBinLowEdge(1);
  const Float_t xStopRS   = hSumR  -> GetBinLowEdge(nBinRS + 1);
  const Float_t xStopRN   = hNormR -> GetBinLowEdge(nBinRN + 1);
  const Float_t lineXS[2] = {xStartRS, xStopRS};
  const Float_t lineXN[2] = {xStartRN, xStopRN};
  const Float_t lineY[2]  = {1., 1.};

  TLine *liSum  = new TLine(lineXS[0], lineY[0], lineXS[1], lineY[1]);
  TLine *liNorm = new TLine(lineXN[0], lineY[0], lineXN[1], lineY[1]);
  liSum  -> SetLineColor(fColL);
  liSum  -> SetLineStyle(fLinL);
  liNorm -> SetLineColor(fColL);
  liNorm -> SetLineStyle(fLinL);
  cout << "    Made lines." << endl;


  // make plots
  const UInt_t  offLogY(0);
  const UInt_t  onLogY(1);
  const UInt_t  grid(0);
  const UInt_t  tick(1);
  const UInt_t  mode(0);
  const UInt_t  border(2);
  const UInt_t  width(1500);
  const UInt_t  height(750);
  const Float_t fMarBig(1.18);
  const Float_t fMarSmall(0.02);
  const Float_t padXR[2]  = {0., 0.5};
  const Float_t padXSN[2] = {0.5, 1.};
  const Float_t padYR[2]  = {0., 1.};
  const Float_t padYSN[2] = {0., 1.};

  TCanvas *cSum  = new TCanvas("cSum", "", width, height);
  TPad    *pSumR = new TPad("pRatioSum", "", padXR[0], padYR[0], padXR[1], padYR[1]);
  TPad    *pSum  = new TPad("pSum", "", padXSN[0], padYSN[0], padXSN[1], padYSN[1]);
  pSumR  -> SetLogy(offLogY);
  pSumR  -> SetGrid(grid, grid);
  pSumR  -> SetTicks(tick, tick);
  pSumR  -> SetBorderMode(mode);
  pSumR  -> SetBorderSize(border);
  pSumR  -> SetTopMargin(fMarSmall);
  pSumR  -> SetRightMargin(fMarSmall);
  pSumR  -> SetBottomMargin(fMarBig);
  pSumR  -> SetLeftMargin(fMarBig);
  pSum   -> SetLogy(onLogY);
  pSum   -> SetGrid(grid, grid);
  pSum   -> SetTicks(tick, tick);
  pSum   -> SetBorderMode(mode);
  pSum   -> SetBorderSize(border);
  pSum   -> SetTopMargin(fMarSmall);
  pSum   -> SetRightMargin(fMarSmall);
  pSum   -> SetBottomMargin(fMarBig);
  pSum   -> SetLeftMargin(fMarBig);
  cSum   -> cd();
  pSumR  -> Draw();
  pSum   -> Draw();
  pSumR  -> cd();
  hSumR  -> Draw("E2");
  liSum  -> Draw();
  ptSum  -> Draw();
  ptSumF -> Draw();
  pSum   -> cd();
  hSumP  -> Draw("E2");
  hSumD  -> Draw("E2 SAME");
  lSum   -> Draw();
  fOut   -> cd();
  cSum   -> Write();
  cSum   -> Close();

  TCanvas *cNorm  = new TCanvas("cNorm", "", width, height);
  TPad    *pNormR = new TPad("pRatioNorm", "", padXR[0], padYR[0], padXR[1], padYR[1]);
  TPad    *pNorm  = new TPad("pNorm", "", padXSN[0], padYSN[0], padXSN[1], padYSN[1]);
  pNormR  -> SetLogy(offLogY);
  pNormR  -> SetGrid(grid, grid);
  pNormR  -> SetTicks(tick, tick);
  pNormR  -> SetBorderMode(mode);
  pNormR  -> SetBorderSize(border);
  pNormR  -> SetTopMargin(fMarSmall);
  pNormR  -> SetRightMargin(fMarSmall);
  pNormR  -> SetBottomMargin(fMarBig);
  pNormR  -> SetLeftMargin(fMarBig);
  pNorm   -> SetLogy(onLogY);
  pNorm   -> SetGrid(grid, grid);
  pNorm   -> SetTicks(tick, tick);
  pNorm   -> SetBorderMode(mode);
  pNorm   -> SetBorderSize(border);
  pNorm   -> SetTopMargin(fMarSmall);
  pNorm   -> SetRightMargin(fMarSmall);
  pNorm   -> SetBottomMargin(fMarBig);
  pNorm   -> SetLeftMargin(fMarBig);
  cNorm   -> cd();
  pNormR  -> Draw();
  pNorm   -> Draw();
  pNormR  -> cd();
  hNormR  -> Draw("E2");
  liNorm  -> Draw();
  ptNorm  -> Draw();
  ptNormF -> Draw();
  pNorm   -> cd();
  hNormP  -> Draw("E2");
  hNormD  -> Draw("E2 SAME");
  lNorm   -> Draw();
  fOut    -> cd();
  cNorm   -> Write();
  cNorm   -> Close();
  cout << "    Made plots." << endl;


  // save and close
  fOut   -> cd();
  hSumP  -> Write();
  hSumD  -> Write();
  hSumR  -> Write();
  fSumR  -> Write();
  hNormP -> Write();
  hNormD -> Write();
  hNormR -> Write();
  fNormR -> Write();
  fOut   -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
