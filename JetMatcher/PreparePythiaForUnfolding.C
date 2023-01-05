// 'PreparePythiaForUnfolding.C'
// Derek Anderson
// 09.07.2018
//
// Use this to create a file which contains
// the response matrix, reconstruction
// efficiency, prior, etc. for Pythia

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"

using namespace std;

// numerical constant(s)
static const Double_t NTrgs(100000);

// io constants
static const Bool_t  VariableBins(true);
static const Bool_t  UseRootNTrgs(true);
static const Float_t JetResParameter(0.5);
static const TString SParticle("ParticleJets/hJetPtCorrP");
static const TString SDetector("MatchJets/hJetPtCorrM");
static const TString SParEff("EventInfo/hParPtCorr");
static const TString SDetEff("EventInfo/hDetPtCorr");
static const TString SResponse("hResponsePtc");
static const TString SProfile("pResponsePtc");
static const TString SNorm("EventInfo/hRefmultP");



void PreparePythiaForUnfolding() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning preparation for unfolding..." << endl;

  // io parameters
  const TString sOut("pp200py8.forJetEffPlot_pTbinFine.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root");
  const TString sIn("pp200py8merge.forJetEffPlot_pTbinFine.et920pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d14m10y2021.root");

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn  = new TFile(sIn.Data(), "read");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH1D     *hJets[2];
  TH1D     *hJetEff[2];
  TH2D     *hResponse;
  TProfile *pResponse;
  hJets[0]   = (TH1D*)     fIn -> Get(SParticle.Data());
  hJets[1]   = (TH1D*)     fIn -> Get(SDetector.Data());
  hJetEff[0] = (TH1D*)     fIn -> Get(SParEff.Data());
  hJetEff[1] = (TH1D*)     fIn -> Get(SDetEff.Data());
  hResponse  = (TH2D*)     fIn -> Get(SResponse.Data());
  pResponse  = (TProfile*) fIn -> Get(SProfile.Data());
  if (!hJets[0] || !hJets[1]) {
    cerr << "PANIC: couldn't grab jet distribution!\n"
         << "       hJetPar = " << hJets[0]
         << ", hJetDet = " << hJets[1]
         << endl;
    return;
  }
  if (!hJetEff[0] || !hJetEff[1]) {
    cerr << "PANIC: couldn't grab efficiency histogram!\n"
         << "       hParEff = " << hJetEff[0]
         << ", hDetEff = " << hJetEff[1]
         << endl;
    return;
  }
  if (!hResponse || !pResponse) {
    cerr << "PANIC: couldn't grab response matrix!\n"
         << "       hResponse = " << hResponse
         << ", pResponse = " << pResponse
         << endl;
    return;
  }
  cout << "    Grabbed histograms." << endl;

  // correct for bin size (if need be)
  if (VariableBins) {
    const UInt_t nBinP = hJets[0] -> GetNbinsX();
    const UInt_t nBinD = hJets[1] -> GetNbinsX();
    for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
      const Double_t pWidth  = hJets[0]   -> GetBinWidth(iBinP);
      const Double_t pJetVal = hJets[0]   -> GetBinContent(iBinP);
      const Double_t pEffVal = hJetEff[0] -> GetBinContent(iBinP);
      const Double_t pJetErr = hJets[0]   -> GetBinError(iBinP);
      const Double_t pEffErr = hJetEff[0] -> GetBinError(iBinP);
      hJets[0]   -> SetBinContent(iBinP, pJetVal / pWidth);
      hJetEff[0] -> SetBinContent(iBinP, pEffVal / pWidth);
      hJets[0]   -> SetBinError(iBinP, pJetErr / pWidth);
      hJetEff[0] -> SetBinError(iBinP, pEffErr / pWidth);
    }
    for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
      const Double_t dWidth  = hJets[1]   -> GetBinWidth(iBinD);
      const Double_t dJetVal = hJets[1]   -> GetBinContent(iBinD);
      const Double_t dEffVal = hJetEff[1] -> GetBinContent(iBinD);
      const Double_t dJetErr = hJets[1]   -> GetBinError(iBinD);
      const Double_t dEffErr = hJetEff[1] -> GetBinError(iBinD);
      hJets[1]   -> SetBinContent(iBinD, dJetVal / dWidth);
      hJetEff[1] -> SetBinContent(iBinD, dEffVal / dWidth);
      hJets[1]   -> SetBinError(iBinD, dJetErr / dWidth);
      hJetEff[1] -> SetBinError(iBinD, dEffErr / dWidth);
    }

    const UInt_t nBinResX = hResponse -> GetNbinsX();
    const UInt_t nBinResY = hResponse -> GetNbinsY();
    for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
      for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
        const Double_t binVal = hResponse -> GetBinContent(iBinResX, iBinResY);
        const Double_t binErr = hResponse -> GetBinError(iBinResX, iBinResY);
        const Double_t xWidth = hResponse -> GetXaxis() -> GetBinWidth(iBinResX);
        const Double_t yWidth = hResponse -> GetYaxis() -> GetBinWidth(iBinResY);
        const Double_t dArea  = xWidth * yWidth;
        hResponse -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
        hResponse -> SetBinError(iBinResX, iBinResY, binErr / dArea);
      }
    }
    const UInt_t nBinResX = pResponse -> GetNbinsX();
    const UInt_t nBinResY = pResponse -> GetNbinsY();
    for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
      for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
        const Double_t binVal = pResponse -> GetBinContent(iBinResX, iBinResY);
        const Double_t binErr = pResponse -> GetBinError(iBinResX, iBinResY);
        const Double_t xWidth = pResponse -> GetXaxis() -> GetBinWidth(iBinResX);
        const Double_t yWidth = pResponse -> GetYaxis() -> GetBinWidth(iBinResY);
        const Double_t dArea  = xWidth * yWidth;
        pResponse -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
        pResponse -> SetBinError(iBinResX, iBinResY, binErr / dArea);
      }
    }
  } else {
    const Double_t parWidth  = hJets[0] -> GetBinWidth(17);
    const Double_t detWidth  = hJets[1] -> GetBinWidth(17);
    hJets[0]   -> Scale(1. / parWidth);
    hJets[1]   -> Scale(1. / detWidth);
    hJetEff[0] -> Scale(1. / parWidth);
    hJetEff[1] -> Scale(1. / detWidth);
  }
  cout << "    Scaled by bin width." << endl;

  // determine what to normalize with
  Double_t norm(1.);
  TH1D     *hNorm;
  if (UseRootNTrgs) {
    hNorm = (TH1D*) fIn -> Get(SNorm.Data());
    norm  = hNorm -> GetEntries();
  } else {
    norm = NTrgs;
  }
  cout << "    Determined normalization factor." << endl;

  // for normalization
  const Double_t etaNorm = 2. * (1. - JetResParameter);
  const Double_t jetNorm = norm * etaNorm;

  // normalize distributions
  hJets[0]   -> Scale(1. / jetNorm);
  hJets[1]   -> Scale(1. / jetNorm);
  hJetEff[0] -> Scale(1. / jetNorm);
  hJetEff[1] -> Scale(1. / jetNorm);

  // normalize response matrices
  const UInt_t nXbins = hResponse -> GetNbinsX();
  const UInt_t nYbins = hResponse -> GetNbinsY();
  for (UInt_t iBinY = 0; iBinY < nYbins; iBinY++) {
    const Double_t binNorm = hResponse -> Integral(1, nXbins, iBinY, iBinY);
    if (binNorm != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbins; iBinX++) {
        const Double_t binVal = hResponse -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = hResponse -> GetBinError(iBinX, iBinY);
        const Double_t xWidth = hResponse -> GetXaxis() -> GetBinWidth(iBinX);
        const Double_t yWidth = hResponse -> GetYaxis() -> GetBinWidth(iBinY);
        const Double_t dArea  = xWidth * yWidth;
        const Double_t newVal = (binVal / binNorm);
        const Double_t newErr = (binErr / binNorm);
        hResponse -> SetBinContent(iBinX, iBinY, newVal);
        hResponse -> SetBinError(iBinX, iBinY, newErr);
      }
    }
  }
  cout << "    Normalized histograms." << endl;

  // calculate efficiency
  TH1D *hEff = (TH1D*) hJetEff[0] -> Clone();
  hEff -> SetName("hEfficiency");
  hEff -> Reset();
  hEff -> Divide(hJetEff[1], hJetEff[0], 1., 1.);
  cout << "    Calculated efficiencies." << endl;

  // fit efficiencies
  const TString sFit("[0]*(1.-TMath::Exp(-1.*[1]*x))");
  const Float_t guess0(0.9);
  const Float_t guess1(5.); 
  const Float_t xyFit[2] = {0.2, 8.0};

  TF1 *fFit = new TF1("fEff", sFit.Data(), xyFit[0], xyFit[1]);
  fFit -> SetParameter(0, guess0);
  fFit -> SetParameter(1, guess1);
  hEff -> Fit(fFit, "R0");
  hEff -> GetFunction("fEff")  -> ResetBit(1<<9);
  cout << "    Fit efficiencies." << endl;

  // set names
  hJets[0]   -> SetName("hParticle");
  hJets[1]   -> SetName("hDetector");
  hJetEff[0] -> SetName("hEffPar");
  hJetEff[1] -> SetName("hEffDet");
  hResponse  -> SetName("hResponse");
  pResponse  -> SetName("pResponse");

  // close files
  fOut       -> cd();
  hJets[0]   -> Write();
  hJets[1]   -> Write();
  hJetEff[0] -> Write();
  hJetEff[1] -> Write();
  hEff       -> Write();
  hResponse  -> Write();
  pResponse  -> Write();
  fOut       -> Close();
  fIn        -> cd();
  fIn        -> Close();
  cout << "  Preparations finished!\n" << endl;

}

// Ending ---------------------------------------------------------------------
