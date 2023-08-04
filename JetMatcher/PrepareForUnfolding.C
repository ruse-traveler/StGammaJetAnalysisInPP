// 'PrepareForUnfolding.C'
// Derek Anderson
// 05.24.2018
//
// Use this to create a file which contains
// the response matrix, reconstruction
// efficiency, prior, etc. for the RFF and
// FF embedding configurations.

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"

using namespace std;

// numerical constants
static const UInt_t   NHistFF(7);
static const UInt_t   NHistRFF(8);
static const UInt_t   NTotal(10);
static const Double_t NTrgsFF[NHistFF]   = {2., 25., 111., 357., 2039., 7688., 4063.};
static const Double_t NTrgsRFF[NHistRFF] = {1., 3., 26., 152., 515., 2350., 6681., 3549.};
static const Double_t WeightsFF[NTotal]  = {1.0, 3.361596e-01, 1.401161e-01, 1.337302e-01, 2.895246e-02, 1.042577e-02, 8.294575e-03, 2.064352e-03, 8.088693e-05, 1.417116e-05};
static const Double_t WeightsRFF[NTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};

// io constants
static const Bool_t  VariableBins(true);
static const Bool_t  UseRootNTrgs(true);
static const Float_t JetResParameter(0.2);
static const TString SParticle("ParticleJets/hJetPtCorrP");
static const TString SDetector("MatchJets/hJetPtCorrM");
static const TString SParEff("EventInfo/hParPtCorr");
static const TString SDetEff("EventInfo/hDetPtCorr");
static const TString SResponse("hResponsePtc");
static const TString SProfile("pResponsePtc");
static const TString SNorm("EventInfo/hRefmultP");



void PrepareForUnfolding() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning preparation for unfolding..." << endl;

  // io parameters
  const TString sOut("pp200r9embed.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root");
  const TString SInFF[NHistFF] = {
    "output/pp200r9pt5ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt7ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt9ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt11ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt15ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt25ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt35ff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root"
  };
  const TString SInRFF[NHistRFF] = {
    "output/pp200r9pt4rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt5rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt7rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt9rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt11rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt15rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt25rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root",
    "output/pp200r9pt35rff.forClosureTest_pTbinGiant.et911r02qt05130.d3m8y2023.root"
  };

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInFF[NHistFF];
  TFile *fInRFF[NHistRFF];
  for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
    fInFF[iHistFF] = new TFile(SInFF[iHistFF].Data(), "read");
    if (!fInFF[iHistFF]) {
      cerr << "PANIC: couldn't open FF input no. " << iHistRFF << "!" << endl;
      return;
    }
  }
  for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
    fInRFF[iHistRFF] = new TFile(SInRFF[iHistRFF].Data(), "read");
    if (!fInRFF[iHistRFF]) {
      cerr << "PANIC: couldn't open RFF input no. " << iHistRFF << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH1D     *hJetsFF[NHistFF][2];
  TH1D     *hJetsRFF[NHistRFF][2];
  TH1D     *hJetEffFF[NHistFF][2];
  TH1D     *hJetEffRFF[NHistRFF][2];
  TH2D     *hResponseFF[NHistFF];
  TH2D     *hResponseRFF[NHistRFF];
  TProfile *pResponseFF[NHistFF];
  TProfile *pResponseRFF[NHistRFF];
  for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
    hJetsFF[iHistFF][0]   = (TH1D*)     fInFF[iHistFF] -> Get(SParticle.Data());
    hJetsFF[iHistFF][1]   = (TH1D*)     fInFF[iHistFF] -> Get(SDetector.Data());
    hJetEffFF[iHistFF][0] = (TH1D*)     fInFF[iHistFF] -> Get(SParEff.Data());
    hJetEffFF[iHistFF][1] = (TH1D*)     fInFF[iHistFF] -> Get(SDetEff.Data());
    hResponseFF[iHistFF]  = (TH2D*)     fInFF[iHistFF] -> Get(SResponse.Data());
    pResponseFF[iHistFF]  = (TProfile*) fInFF[iHistFF] -> Get(SProfile.Data());
    if (!hJetsFF[iHistFF][0] || !hJetsFF[iHistFF][1]) {
      cerr << "PANIC: couldn't grab FF distribution!\n"
           << "       hJetPar[" << iHistFF << "] = " << hJetsFF[iHistFF][0]
           << ", hJetDet[" << iHistFF << "] = " << hJetsFF[iHistFF][1]
           << endl;
      return;
    }
    if (!hJetEffFF[iHistFF][0] || !hJetEffFF[iHistFF][1]) {
      cerr << "PANIC: couldn't grab FF efficiency histogram!\n"
           << "       hParEff[" << iHistFF << "] = " << hJetEffFF[iHistFF][0]
           << ", hDetEff[" << iHistFF << "] = " << hJetEffFF[iHistFF][1]
           << endl;
      return;
    }
    if (!hResponseFF[iHistFF] || !pResponseFF[iHistFF]) {
      cerr << "PANIC: couldn't grab FF response matrix!\n"
           << "       hResponse[" << iHistFF << "] = " << hResponseFF[iHistFF]
           << ", pResponse[" << iHistFF << "] = " << pResponseFF[iHistFF]
           << endl;
      return;
    }
  }
  for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
    hJetsRFF[iHistRFF][0]   = (TH1D*)     fInRFF[iHistRFF] -> Get(SParticle.Data());
    hJetsRFF[iHistRFF][1]   = (TH1D*)     fInRFF[iHistRFF] -> Get(SDetector.Data());
    hJetEffRFF[iHistRFF][0] = (TH1D*)     fInRFF[iHistRFF] -> Get(SParEff.Data());
    hJetEffRFF[iHistRFF][1] = (TH1D*)     fInRFF[iHistRFF] -> Get(SDetEff.Data());
    hResponseRFF[iHistRFF]  = (TH2D*)     fInRFF[iHistRFF] -> Get(SResponse.Data());
    pResponseRFF[iHistRFF]  = (TProfile*) fInRFF[iHistRFF] -> Get(SProfile.Data());
    if (!hJetsRFF[iHistRFF][0] || !hJetsRFF[iHistRFF][1]) {
      cerr << "PANIC: couldn't grab RFF distribution!\n"
           << "       hJetPar[" << iHistRFF << "] = " << hJetsRFF[iHistRFF][0]
           << ", hJetDet[" << iHistRFF << "] = " << hJetsFF[iHistRFF][1]
           << endl;
      return;
    }
    if (!hJetEffRFF[iHistRFF][0] || !hJetEffRFF[iHistRFF][1]) {
      cerr << "PANIC: couldn't grab RFF efficiency histogram!\n"
           << "       hParEff[" << iHistRFF << "] = " << hJetEffRFF[iHistRFF][0]
           << ", hDetEff[" << iHistRFF << "] = " << hJetEffRFF[iHistRFF][1]
           << endl;
      return;
    }
    if (!hResponseRFF[iHistRFF] || !pResponseRFF[iHistRFF]) {
      cerr << "PANIC: couldn't grab RFF response matrix!\n"
           << "       hResponse[" << iHistRFF << "] = " << hResponseRFF[iHistRFF]
           << ", pResponse[" << iHistRFF << "] = " << pResponseRFF[iHistRFF]
           << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // correct for bin size (if need be)
  if (VariableBins) {
    for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
      const UInt_t nBinP = hJetsFF[iHistFF][0] -> GetNbinsX();
      const UInt_t nBinD = hJetsFF[iHistFF][1] -> GetNbinsX();
      for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
        const Double_t pWidth  = hJetsFF[iHistFF][0]   -> GetBinWidth(iBinP);
        const Double_t pJetVal = hJetsFF[iHistFF][0]   -> GetBinContent(iBinP);
        const Double_t pEffVal = hJetEffFF[iHistFF][0] -> GetBinContent(iBinP);
        const Double_t pJetErr = hJetsFF[iHistFF][0]   -> GetBinError(iBinP);
        const Double_t pEffErr = hJetEffFF[iHistFF][0] -> GetBinError(iBinP);
        hJetsFF[iHistFF][0]   -> SetBinContent(iBinP, pJetVal / pWidth);
        hJetEffFF[iHistFF][0] -> SetBinContent(iBinP, pEffVal / pWidth);
        hJetsFF[iHistFF][0]   -> SetBinError(iBinP, pJetErr / pWidth);
        hJetEffFF[iHistFF][0] -> SetBinError(iBinP, pEffErr / pWidth);
      }
      for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
        const Double_t dWidth  = hJetsFF[iHistFF][1]   -> GetBinWidth(iBinD);
        const Double_t dJetVal = hJetsFF[iHistFF][1]   -> GetBinContent(iBinD);
        const Double_t dEffVal = hJetEffFF[iHistFF][1] -> GetBinContent(iBinD);
        const Double_t dJetErr = hJetsFF[iHistFF][1]   -> GetBinError(iBinD);
        const Double_t dEffErr = hJetEffFF[iHistFF][1] -> GetBinError(iBinD);
        hJetsFF[iHistFF][1]   -> SetBinContent(iBinD, dJetVal / dWidth);
        hJetEffFF[iHistFF][1] -> SetBinContent(iBinD, dEffVal / dWidth);
        hJetsFF[iHistFF][1]   -> SetBinError(iBinD, dJetErr / dWidth);
        hJetEffFF[iHistFF][1] -> SetBinError(iBinD, dEffErr / dWidth);
      }

      const UInt_t nBinResX = hResponseFF[iHistFF] -> GetNbinsX();
      const UInt_t nBinResY = hResponseFF[iHistFF] -> GetNbinsY();
      for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
        for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
          const Double_t binVal = hResponseFF[iHistFF] -> GetBinContent(iBinResX, iBinResY);
          const Double_t binErr = hResponseFF[iHistFF] -> GetBinError(iBinResX, iBinResY);
          const Double_t xWidth = hResponseFF[iHistFF] -> GetXaxis() -> GetBinWidth(iBinResX);
          const Double_t yWidth = hResponseFF[iHistFF] -> GetYaxis() -> GetBinWidth(iBinResY);
          const Double_t dArea  = xWidth * yWidth;
          hResponseFF[iHistFF] -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
          hResponseFF[iHistFF] -> SetBinError(iBinResX, iBinResY, binErr / dArea);
        }
      }
      const UInt_t nBinResX = pResponseFF[iHistFF] -> GetNbinsX();
      const UInt_t nBinResY = pResponseFF[iHistFF] -> GetNbinsY();
      for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
        for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
          const Double_t binVal = pResponseFF[iHistFF] -> GetBinContent(iBinResX, iBinResY);
          const Double_t binErr = pResponseFF[iHistFF] -> GetBinError(iBinResX, iBinResY);
          const Double_t xWidth = pResponseFF[iHistFF] -> GetXaxis() -> GetBinWidth(iBinResX);
          const Double_t yWidth = pResponseFF[iHistFF] -> GetYaxis() -> GetBinWidth(iBinResY);
          const Double_t dArea  = xWidth * yWidth;
          pResponseFF[iHistFF] -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
          pResponseFF[iHistFF] -> SetBinError(iBinResX, iBinResY, binErr / dArea);
        }
      }
    }  // end FF loop
    for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
      const UInt_t nBinP = hJetsRFF[iHistRFF][0] -> GetNbinsX();
      const UInt_t nBinD = hJetsRFF[iHistRFF][1] -> GetNbinsX();
      for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
        const Double_t pWidth  = hJetsRFF[iHistRFF][0]   -> GetBinWidth(iBinP);
        const Double_t pJetVal = hJetsRFF[iHistRFF][0]   -> GetBinContent(iBinP);
        const Double_t pEffVal = hJetEffRFF[iHistRFF][0] -> GetBinContent(iBinP);
        const Double_t pJetErr = hJetsRFF[iHistRFF][0]   -> GetBinError(iBinP);
        const Double_t pEffErr = hJetEffRFF[iHistRFF][0] -> GetBinError(iBinP);
        hJetsRFF[iHistRFF][0]   -> SetBinContent(iBinP, pJetVal / pWidth);
        hJetEffRFF[iHistRFF][0] -> SetBinContent(iBinP, pEffVal / pWidth);
        hJetsRFF[iHistRFF][0]   -> SetBinError(iBinP, pJetErr / pWidth);
        hJetEffRFF[iHistRFF][0] -> SetBinError(iBinP, pEffErr / pWidth);
      }
      for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
        const Double_t dWidth  = hJetsRFF[iHistRFF][1] -> GetBinWidth(iBinD);
        const Double_t dJetVal = hJetsRFF[iHistRFF][1] -> GetBinContent(iBinD);
        const Double_t dEffVal = hJetEffRFF[iHistRFF][1] -> GetBinContent(iBinD);
        const Double_t dJetErr = hJetsRFF[iHistRFF][1] -> GetBinError(iBinD);
        const Double_t dEffErr = hJetEffRFF[iHistRFF][1] -> GetBinError(iBinD);
        hJetsRFF[iHistRFF][1]   -> SetBinContent(iBinD, dJetVal / dWidth);
        hJetEffRFF[iHistRFF][1] -> SetBinContent(iBinD, dEffVal / dWidth);
        hJetsRFF[iHistRFF][1]   -> SetBinError(iBinD, dJetErr / dWidth);
        hJetEffRFF[iHistRFF][1] -> SetBinError(iBinD, dEffErr / dWidth);
      }

      const UInt_t nBinResX = hResponseRFF[iHistRFF] -> GetNbinsX();
      const UInt_t nBinResY = hResponseRFF[iHistRFF] -> GetNbinsY();
      for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
        for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
          const Double_t binVal = hResponseRFF[iHistRFF] -> GetBinContent(iBinResX, iBinResY);
          const Double_t binErr = hResponseRFF[iHistRFF] -> GetBinError(iBinResX, iBinResY);
          const Double_t xWidth = hResponseRFF[iHistRFF] -> GetXaxis() -> GetBinWidth(iBinResX);
          const Double_t yWidth = hResponseRFF[iHistRFF] -> GetYaxis() -> GetBinWidth(iBinResY);
          const Double_t dArea  = xWidth * yWidth;
          hResponseRFF[iHistRFF] -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
          hResponseRFF[iHistRFF] -> SetBinError(iBinResX, iBinResY, binErr / dArea);
        }
      }
      const UInt_t nBinResX = pResponseRFF[iHistRFF] -> GetNbinsX();
      const UInt_t nBinResY = pResponseRFF[iHistRFF] -> GetNbinsY();
      for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
        for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
          const Double_t binVal = pResponseRFF[iHistRFF] -> GetBinContent(iBinResX, iBinResY);
          const Double_t binErr = pResponseRFF[iHistRFF] -> GetBinError(iBinResX, iBinResY);
          const Double_t xWidth = pResponseRFF[iHistRFF] -> GetXaxis() -> GetBinWidth(iBinResX);
          const Double_t yWidth = pResponseRFF[iHistRFF] -> GetYaxis() -> GetBinWidth(iBinResY);
          const Double_t dArea  = xWidth * yWidth;
          pResponseRFF[iHistRFF] -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
          pResponseRFF[iHistRFF] -> SetBinError(iBinResX, iBinResY, binErr / dArea);
        }
      }
    }  // end RFF loop
  } else {
    const Double_t parWidthFF  = hJetsFF[0][0]  -> GetBinWidth(17);
    const Double_t parWidthRFF = hJetsRFF[0][0] -> GetBinWidth(17);
    const Double_t detWidthFF  = hJetsFF[0][1]  -> GetBinWidth(17);
    const Double_t detWidthRFF = hJetsRFF[0][1] -> GetBinWidth(17);
    for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
      hJetsFF[iHistFF][0]   -> Scale(1. / parWidthFF);
      hJetsFF[iHistFF][1]   -> Scale(1. / detWidthFF);
      hJetEffFF[iHistFF][0] -> Scale(1. / parWidthFF);
      hJetEffFF[iHistFF][1] -> Scale(1. / detWidthFF);
    }
    for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
      hJetsRFF[iHistRFF][0]   -> Scale(1. / parWidthRFF);
      hJetsRFF[iHistRFF][1]   -> Scale(1. / detWidthRFF);
      hJetEffRFF[iHistRFF][0] -> Scale(1. / parWidthRFF);
      hJetEffRFF[iHistRFF][1] -> Scale(1. / detWidthRFF);
    }  
  }
  cout << "    Scaled by bin width." << endl;

  // determine what to normalize with
  Double_t NormsFF[NHistFF];
  Double_t NormsRFF[NHistRFF];
  TH1D     *hNormFF[NHistFF];
  TH1D     *hNormRFF[NHistRFF];
  if (UseRootNTrgs) {
    for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
      hNormFF[iHistFF] = (TH1D*) fInFF[iHistFF] -> Get(SNorm.Data());
      NormsFF[iHistFF] = hNormFF[iHistFF] -> GetEntries();
    }
    for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
      hNormRFF[iHistRFF] = (TH1D*) fInRFF[iHistRFF] -> Get(SNorm.Data());
      NormsRFF[iHistRFF] = hNormRFF[iHistRFF] -> GetEntries();
    }
  } else {
    for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
      NormsFF[iHistFF] = NTrgsFF[iHistFF];
      cout << "CHECK [" << iHistFF << "]: " << NormsFF[iHistFF] << endl;
    }
    for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
      NormsRFF[iHistRFF] = NTrgsRFF[iHistRFF];
      cout << "CHECK [" << iHistRFF << "]: " << NormsRFF[iHistRFF] << endl;
    }
  }
  cout << "    Determined normalization factors." << endl;

  // scale histograms
  Double_t normFF(0.);
  Double_t normRFF(0.);
  Double_t normAll(0.);
  for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
    hJetsFF[iHistFF][0]   -> Scale(WeightsFF[(NTotal - NHistFF) + iHistFF]);
    hJetsFF[iHistFF][1]   -> Scale(WeightsFF[(NTotal - NHistFF) + iHistFF]);
    hJetEffFF[iHistFF][0] -> Scale(WeightsFF[(NTotal - NHistFF) + iHistFF]);
    hJetEffFF[iHistFF][1] -> Scale(WeightsFF[(NTotal - NHistFF) + iHistFF]);
    hResponseFF[iHistFF]  -> Scale(WeightsFF[(NTotal - NHistFF) + iHistFF]);
    pResponseFF[iHistFF]  -> Scale(WeightsFF[(NTotal - NHistFF) + iHistFF]);
    normFF  += NormsFF[iHistFF] * WeightsFF[(NTotal - NHistFF) + iHistFF];
    normAll += NormsFF[iHistFF] * WeightsFF[(NTotal - NHistFF) + iHistFF];
  }
  for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
    hJetsRFF[iHistRFF][0]   -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHistRFF]);
    hJetsRFF[iHistRFF][1]   -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHistRFF]);
    hJetEffRFF[iHistRFF][0] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHistRFF]);
    hJetEffRFF[iHistRFF][1] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHistRFF]);
    hResponseRFF[iHistRFF]  -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHistRFF]);
    pResponseRFF[iHistRFF]  -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHistRFF]);
    normRFF += NormsRFF[iHistRFF] * WeightsRFF[(NTotal - NHistRFF) + iHistRFF];
    normAll += NormsRFF[iHistRFF] * WeightsRFF[(NTotal - NHistRFF) + iHistRFF];
  }
  cout << "    Scaled histograms:\n"
       << "      norm(FF,RFF,All) = (" << normFF << ", " << normRFF << ", " << normAll << ")"
       << endl;

  // for sums
  TH1D     *hJetSumFF[2];
  TH1D     *hJetSumRFF[2];
  TH1D     *hJetSumAll[2];
  TH1D     *hEffSumFF[2];
  TH1D     *hEffSumRFF[2];
  TH1D     *hEffSumAll[2];
  TH2D     *hResSumFF;
  TH2D     *hResSumRFF;
  TH2D     *hResSumAll;
  TProfile *pResSumFF;
  TProfile *pResSumRFF;
  TProfile *pResSumAll;
  // initialize FF sums
  hJetSumFF[0] = (TH1D*)     hJetsFF[0][0]   -> Clone();
  hJetSumFF[1] = (TH1D*)     hJetsFF[0][1]   -> Clone();
  hEffSumFF[0] = (TH1D*)     hJetEffFF[0][0] -> Clone();
  hEffSumFF[1] = (TH1D*)     hJetEffFF[0][1] -> Clone();
  hResSumFF    = (TH2D*)     hResponseFF[0]  -> Clone();
  pResSumFF    = (TProfile*) pResponseFF[0] -> Clone();
  hJetSumFF[0] -> SetName("hSumParFF");
  hJetSumFF[1] -> SetName("hSumDetFF");
  hEffSumFF[0] -> SetName("hEffParFF");
  hEffSumFF[1] -> SetName("hEffDetFF");
  hResSumFF    -> SetName("hResponseFF");
  pResSumFF    -> SetName("pResponseFF");
  hJetSumFF[0] -> Reset("ICE");
  hJetSumFF[1] -> Reset("ICE");
  hEffSumFF[0] -> Reset("ICE");
  hEffSumFF[1] -> Reset("ICE");
  hResSumFF    -> Reset("ICE");
  pResSumFF    -> Reset("ICE");
  // initialize RFF sums
  hJetSumRFF[0] = (TH1D*)     hJetsRFF[0][0]   -> Clone();
  hJetSumRFF[1] = (TH1D*)     hJetsRFF[0][1]   -> Clone();
  hEffSumRFF[0] = (TH1D*)     hJetEffRFF[0][0] -> Clone();
  hEffSumRFF[1] = (TH1D*)     hJetEffRFF[0][1] -> Clone();
  hResSumRFF    = (TH2D*)     hResponseRFF[0]  -> Clone();
  pResSumRFF    = (TProfile*) pResponseRFF[0]  -> Clone();
  hJetSumRFF[0] -> SetName("hSumParRFF");
  hJetSumRFF[1] -> SetName("hSumDetRFF");
  hEffSumRFF[0] -> SetName("hEffParRFF");
  hEffSumRFF[1] -> SetName("hEffDetRFF");
  hResSumRFF    -> SetName("hResponseRFF");
  pResSumRFF    -> SetName("pResponseRFF");
  hJetSumRFF[0] -> Reset("ICE");
  hJetSumRFF[1] -> Reset("ICE");
  hEffSumRFF[0] -> Reset("ICE");
  hEffSumRFF[1] -> Reset("ICE");
  hResSumRFF    -> Reset("ICE");
  pResSumRFF    -> Reset("ICE");
  // initialize All sums
  hJetSumAll[0] = (TH1D*)     hJetsRFF[0][0]   -> Clone();
  hJetSumAll[1] = (TH1D*)     hJetsRFF[0][1]   -> Clone();
  hEffSumAll[0] = (TH1D*)     hJetEffRFF[0][0] -> Clone();
  hEffSumAll[1] = (TH1D*)     hJetEffRFF[0][1] -> Clone();
  hResSumAll    = (TH2D*)     hResponseRFF[0]  -> Clone();
  pResSumAll    = (TProfile*) pResponseRFF[0]  -> Clone();
  hJetSumAll[0] -> SetName("hSumParAll");
  hJetSumAll[1] -> SetName("hSumDetAll");
  hEffSumAll[0] -> SetName("hEffParAll");
  hEffSumAll[1] -> SetName("hEffDetAll");
  hResSumAll    -> SetName("hResponseAll");
  pResSumAll    -> SetName("pResponseAll");
  hJetSumAll[0] -> Reset("ICE");
  hJetSumAll[1] -> Reset("ICE");
  hEffSumAll[0] -> Reset("ICE");
  hEffSumAll[1] -> Reset("ICE");
  hResSumAll    -> Reset("ICE");
  pResSumAll    -> Reset("ICE");
  // sum histograms
  for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
    hJetSumFF[0]  -> Add(hJetsFF[iHistFF][0]);
    hJetSumFF[1]  -> Add(hJetsFF[iHistFF][1]);
    hJetSumAll[0] -> Add(hJetsFF[iHistFF][0]);
    hJetSumAll[1] -> Add(hJetsFF[iHistFF][1]);
    hEffSumFF[0]  -> Add(hJetEffFF[iHistFF][0]);
    hEffSumFF[1]  -> Add(hJetEffFF[iHistFF][1]);
    hEffSumAll[0] -> Add(hJetEffFF[iHistFF][0]);
    hEffSumAll[1] -> Add(hJetEffFF[iHistFF][1]);
    hResSumFF     -> Add(hResponseFF[iHistFF]);
    pResSumFF     -> Add(pResponseFF[iHistFF]);
    hResSumAll    -> Add(hResponseFF[iHistFF]);
    pResSumAll    -> Add(pResponseFF[iHistFF]);
  }
  for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
    hJetSumRFF[0] -> Add(hJetsRFF[iHistRFF][0]);
    hJetSumRFF[1] -> Add(hJetsRFF[iHistRFF][1]);
    hJetSumAll[0] -> Add(hJetsRFF[iHistRFF][0]);
    hJetSumAll[1] -> Add(hJetsRFF[iHistRFF][1]);
    hEffSumRFF[0] -> Add(hJetEffRFF[iHistRFF][0]);
    hEffSumRFF[1] -> Add(hJetEffRFF[iHistRFF][1]);
    hEffSumAll[0] -> Add(hJetEffRFF[iHistRFF][0]);
    hEffSumAll[1] -> Add(hJetEffRFF[iHistRFF][1]);
    hResSumRFF    -> Add(hResponseRFF[iHistRFF]);
    pResSumRFF    -> Add(pResponseRFF[iHistRFF]);
    hResSumAll    -> Add(hResponseRFF[iHistRFF]);
    pResSumAll    -> Add(pResponseRFF[iHistRFF]);
  }
  cout << "    Summed histograms." << endl;

  // for normalization
  const Double_t etaNorm    = 2. * (1. - JetResParameter);
  const Double_t jetNormFF  = normFF * etaNorm;
  const Double_t jetNormRFF = normRFF * etaNorm;
  const Double_t jetNormAll = normAll * etaNorm;

  // normalize FF sums
  hJetSumFF[0] -> Scale(1. / jetNormFF);
  hJetSumFF[1] -> Scale(1. / jetNormFF);
  hEffSumFF[0] -> Scale(1. / jetNormFF);
  hEffSumFF[1] -> Scale(1. / jetNormFF);
  // normalize RFF sums
  hJetSumRFF[0] -> Scale(1. / jetNormRFF);
  hJetSumRFF[1] -> Scale(1. / jetNormRFF);
  hEffSumRFF[0] -> Scale(1. / jetNormRFF);
  hEffSumRFF[1] -> Scale(1. / jetNormRFF);
  // normalize All sums
  hJetSumAll[0] -> Scale(1. / jetNormAll);
  hJetSumAll[1] -> Scale(1. / jetNormAll);
  hEffSumAll[0] -> Scale(1. / jetNormAll);
  hEffSumAll[1] -> Scale(1. / jetNormAll);

  // normalize response matrices
  const UInt_t nXbinsFF  = hResSumFF  -> GetNbinsX();
  const UInt_t nYbinsFF  = hResSumFF  -> GetNbinsY();
  const UInt_t nXbinsRFF = hResSumRFF -> GetNbinsX();
  const UInt_t nYbinsRFF = hResSumRFF -> GetNbinsY();
  const UInt_t nXbinsAll = hResSumAll -> GetNbinsX();
  const UInt_t nYbinsAll = hResSumAll -> GetNbinsY();
  for (UInt_t iBinY = 0; iBinY < nYbinsFF; iBinY++) {
    const Double_t binNormFF = hResSumFF -> Integral(1, nXbinsFF, iBinY, iBinY);
    if (binNormFF != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbinsFF; iBinX++) {
        const Double_t binVal = hResSumFF -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = hResSumFF -> GetBinError(iBinX, iBinY);
        const Double_t xWidth = hResSumFF -> GetXaxis() -> GetBinWidth(iBinX);
        const Double_t yWidth = hResSumFF -> GetYaxis() -> GetBinWidth(iBinY);
        const Double_t dArea  = xWidth * yWidth;
        const Double_t newVal = (binVal / binNormFF);
        const Double_t newErr = (binErr / binNormFF);
        hResSumFF -> SetBinContent(iBinX, iBinY, newVal);
        hResSumFF -> SetBinError(iBinX, iBinY, newErr);
      }
    }
  }
  for (UInt_t iBinY = 0; iBinY < nYbinsRFF; iBinY++) {
    const Double_t binNormRFF = hResSumRFF -> Integral(1, nXbinsRFF, iBinY, iBinY);
    if (binNormRFF != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbinsRFF; iBinX++) {
        const Double_t binVal = hResSumRFF -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = hResSumRFF -> GetBinError(iBinX, iBinY);
        const Double_t xWidth = hResSumRFF -> GetXaxis() -> GetBinWidth(iBinX);
        const Double_t yWidth = hResSumRFF -> GetYaxis() -> GetBinWidth(iBinY);
        const Double_t dArea  = xWidth * yWidth;
        const Double_t newVal = (binVal / binNormRFF);
        const Double_t newErr = (binErr / binNormRFF);
        hResSumRFF -> SetBinContent(iBinX, iBinY, newVal);
        hResSumRFF -> SetBinError(iBinX, iBinY, newErr);
      }
    }
  }
  for (UInt_t iBinY = 0; iBinY < nYbinsAll; iBinY++) {
    const Double_t binNormAll = hResSumAll -> Integral(1, nXbinsAll, iBinY, iBinY);
    if (binNormAll != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbinsAll; iBinX++) {
        const Double_t binVal = hResSumAll -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = hResSumAll -> GetBinError(iBinX, iBinY);
        const Double_t xWidth = hResSumAll -> GetXaxis() -> GetBinWidth(iBinX);
        const Double_t yWidth = hResSumAll -> GetYaxis() -> GetBinWidth(iBinY);
        const Double_t dArea  = xWidth * yWidth;
        const Double_t newVal = (binVal / binNormAll);
        const Double_t newErr = (binErr / binNormAll);
        hResSumAll -> SetBinContent(iBinX, iBinY, newVal);
        hResSumAll -> SetBinError(iBinX, iBinY, newErr);
      }
    }
  }
  cout << "    Normalized histograms." << endl;

  // calculate efficiency
  TH1D *hEffFF  = (TH1D*) hEffSumFF[0]  -> Clone();
  TH1D *hEffRFF = (TH1D*) hEffSumRFF[0] -> Clone();
  TH1D *hEffAll = (TH1D*) hEffSumAll[0] -> Clone();
  hEffFF  -> SetName("hEfficiencyFF");
  hEffRFF -> SetName("hEfficiencyRFF");
  hEffAll -> SetName("hEfficiencyAll");
  hEffFF  -> Reset();
  hEffRFF -> Reset();
  hEffAll -> Reset();
  hEffFF  -> Divide(hEffSumFF[1], hEffSumFF[0], 1., 1.);
  hEffRFF -> Divide(hEffSumRFF[1], hEffSumRFF[0], 1., 1.);
  hEffAll -> Divide(hEffSumAll[1], hEffSumAll[0], 1., 1.);
  cout << "    Calculated efficiencies." << endl;

  // fit efficiencies
  const TString sFit("[0]*(1.-TMath::Exp(-1.*[1]*x))");
  const Float_t guess0(0.9);
  const Float_t guess1(5.); 
  const Float_t xyFit[2] = {0.2, 8.0};

  TF1 *fFitFF  = new TF1("fEffFF", sFit.Data(), xyFit[0], xyFit[1]);
  TF1 *fFitRFF = new TF1("fEffRFF", sFit.Data(), xyFit[0], xyFit[1]);
  TF1 *fFitAll = new TF1("fEffAll", sFit.Data(), xyFit[0], xyFit[1]);
  fFitFF  -> SetParameter(0, guess0);
  fFitFF  -> SetParameter(1, guess1);
  fFitRFF -> SetParameter(0, guess0);
  fFitRFF -> SetParameter(1, guess1);
  fFitAll -> SetParameter(0, guess0);
  fFitAll -> SetParameter(1, guess1);

  hEffFF  -> Fit(fFitFF, "R0");
  hEffRFF -> Fit(fFitRFF, "R0");
  hEffAll -> Fit(fFitAll, "R0");
  hEffFF  -> GetFunction("fEffFF")  -> ResetBit(1<<9);
  hEffRFF -> GetFunction("fEffRFF") -> ResetBit(1<<9);
  hEffAll -> GetFunction("fEffAll") -> ResetBit(1<<9);
  cout << "    Fit efficiencies." << endl;

  // close files
  fOut          -> cd();
  hJetSumFF[0]  -> Write();
  hJetSumFF[1]  -> Write();
  hJetSumRFF[0] -> Write();
  hJetSumRFF[1] -> Write();
  hJetSumAll[0] -> Write();
  hJetSumAll[1] -> Write();
  hEffSumFF[0]  -> Write();
  hEffSumFF[1]  -> Write();
  hEffSumRFF[0] -> Write();
  hEffSumRFF[1] -> Write();
  hEffSumAll[0] -> Write();
  hEffSumAll[1] -> Write();
  hEffFF        -> Write();
  hEffRFF       -> Write();
  hEffAll       -> Write();
  hResSumFF     -> Write();
  hResSumRFF    -> Write();
  hResSumAll    -> Write();
  pResSumFF     -> Write();
  pResSumRFF    -> Write();
  pResSumAll    -> Write();
  fOut          -> Close();
  for (UInt_t iHistFF = 0; iHistFF < NHistFF; iHistFF++) {
    fInFF[iHistFF] -> cd();
    fInFF[iHistFF] -> Close();
  }
  for (UInt_t iHistRFF = 0; iHistRFF < NHistRFF; iHistRFF++) {
    fInRFF[iHistRFF] -> cd();
    fInRFF[iHistRFF] -> Close();
  }
  cout << "  Preparations finished!\n" << endl;

}

// Ending ---------------------------------------------------------------------
