// 'PrepareAuAuForUnfolding.C'
// Derek Anderson
// 02.20.2023
//
// Use this to create a file which contains
// the response matrix, reconstruction
// efficiency, prior, etc. from the Run14
// Au+Au embedding sample

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile.h"

using namespace std;

/* Numbers from Run14 Au+Au Embedding Sample:
  PtHatEdge[NTotal] = {9.,        11.,       15.,       20.,       25.,       30.,       40.,       50.}
  NEvts[NTotal]     = {1399880.,  864472.,   1143185.,  969676.,   732554.,   501934.,   296936.,   101671.}
  NTrgs[NTotal]     = {789.,      3550.,     12324.,    1.,        31757.,    26979.,    19216.,    6559.}
  XSec[NTotal]      = {5.524e-03, 2.203e-03, 3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09}
  XSum              = 8.812835e-03
*/

// numerical constants
static const UInt_t   NHist(8);
static const UInt_t   NTotal(8);
static const Double_t XSum(8.812835e-03);  // sum of x-sections
static const Double_t NEvts[NTotal]   = {1399880.,  864472.,   1143185.,  969676.,   732554.,   501934.,   296936.,   101671.};
static const Double_t NTrgs[NHist]    = {789.,      3550.,     12324.,    1.,        31757.,    26979.,    19216.,    6559.};  // placefolder for pThat = 20 - 25 GeV/c right now
static const Double_t XSec[NTotal]    = {5.524e-03, 2.203e-03, 3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09}; 

// io constants
static const Bool_t  VariableBins(false);
static const Bool_t  UseRootNTrgs(true);
static const Float_t JetResParameter(0.5);
static const TString SParticle("ParticleJets/hJetPtCorrP");
static const TString SDetector("MatchJets/hJetPtCorrM");
static const TString SParEff("EventInfo/hParPtCorr");
static const TString SDetEff("EventInfo/hDetPtCorr");
static const TString SResponse("hResponsePtc");
static const TString SProfile("pResponsePtc");
static const TString SNorm("EventInfo/hRefmultP");



void PrepareAuAuForUnfolding() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning preparation for unfolding...\n"
       << "    Calculating weights:"
       << endl;

  // calculate weights
  Double_t Weights[NTotal];
  for (UInt_t iPtHat = 0; iPtHat < NTotal; iPtHat++) {
    Weights[iPtHat] = XSec[iPtHat] / XSum;
    cout << "      PtHat Bin [" << iPtHat << "]: weight = " << Weights[iPtHat] << endl;
  }

  // io parameters
  const TString sOut("test.root");
  const TString SIn[NHist] = {"/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt911.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt1115.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt1520.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt2025.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt2530.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt3040.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt4050.testingSetup.d16m2y2023.root",
                              "/star/data01/pwg/nihar/EmbeddingPythiaInAuAu/JetMatcher/auau200r14pt50inf.testingSetup.d16m2y2023.root"};

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] = new TFile(SIn[iHist].Data(), "read");
    if (!fIn[iHist]) {
      cerr << "PANIC: couldn't open input no. " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH1D     *hJets[NHist][2];
  TH1D     *hJetEff[NHist][2];
  TH2D     *hResponse[NHist];
  TProfile *pResponse[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
    hJets[iHist][0]   = (TH1D*)     fIn[iHist] -> Get(SParticle.Data());
    hJets[iHist][1]   = (TH1D*)     fIn[iHist] -> Get(SDetector.Data());
    hJetEff[iHist][0] = (TH1D*)     fIn[iHist] -> Get(SParEff.Data());
    hJetEff[iHist][1] = (TH1D*)     fIn[iHist] -> Get(SDetEff.Data());
    hResponse[iHist]  = (TH2D*)     fIn[iHist] -> Get(SResponse.Data());
    pResponse[iHist]  = (TProfile*) fIn[iHist] -> Get(SProfile.Data());
    if (!hJets[iHist][0] || !hJets[iHist][1]) {
      cerr << "PANIC: couldn't grab input distribution!\n"
           << "       hJetPar[" << iHist << "] = " << hJets[iHist][0]
           << ", hJetDet[" << iHist << "] = " << hJetsFF[iHist][1]
           << endl;
      return;
    }
    if (!hJetEff[iHist][0] || !hJetEff[iHist][1]) {
      cerr << "PANIC: couldn't grab input efficiency histogram!\n"
           << "       hParEff[" << iHist << "] = " << hJetEff[iHist][0]
           << ", hDetEff[" << iHist << "] = " << hJetEff[iHist][1]
           << endl;
      return;
    }
    if (!hResponse[iHist] || !pResponse[iHist]) {
      cerr << "PANIC: couldn't grab input response matrix!\n"
           << "       hResponse[" << iHist << "] = " << hResponse[iHist]
           << ", pResponse[" << iHist << "] = " << pResponse[iHist]
           << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // correct for bin size (if need be)
  if (VariableBins) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
      const UInt_t nBinP = hJets[iHist][0] -> GetNbinsX();
      const UInt_t nBinD = hJets[iHist][1] -> GetNbinsX();
      for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
        const Double_t pWidth  = hJets[iHist][0]   -> GetBinWidth(iBinP);
        const Double_t pJetVal = hJets[iHist][0]   -> GetBinContent(iBinP);
        const Double_t pEffVal = hJetEff[iHist][0] -> GetBinContent(iBinP);
        const Double_t pJetErr = hJets[iHist][0]   -> GetBinError(iBinP);
        const Double_t pEffErr = hJetEff[iHist][0] -> GetBinError(iBinP);
        hJets[iHist][0]   -> SetBinContent(iBinP, pJetVal / pWidth);
        hJetEff[iHist][0] -> SetBinContent(iBinP, pEffVal / pWidth);
        hJets[iHist][0]   -> SetBinError(iBinP, pJetErr / pWidth);
        hJetEff[iHist][0] -> SetBinError(iBinP, pEffErr / pWidth);
      }
      for (UInt_t iBinD = 1; iBinD < (nBinD + 1); iBinD++) {
        const Double_t dWidth  = hJets[iHist][1]   -> GetBinWidth(iBinD);
        const Double_t dJetVal = hJets[iHist][1]   -> GetBinContent(iBinD);
        const Double_t dEffVal = hJetEff[iHist][1] -> GetBinContent(iBinD);
        const Double_t dJetErr = hJets[iHist][1]   -> GetBinError(iBinD);
        const Double_t dEffErr = hJetEff[iHist][1] -> GetBinError(iBinD);
        hJets[iHist][1]   -> SetBinContent(iBinD, dJetVal / dWidth);
        hJetEff[iHist][1] -> SetBinContent(iBinD, dEffVal / dWidth);
        hJets[iHist][1]   -> SetBinError(iBinD, dJetErr / dWidth);
        hJetEff[iHist][1] -> SetBinError(iBinD, dEffErr / dWidth);
      }

      const UInt_t nBinResX = hResponse[iHist] -> GetNbinsX();
      const UInt_t nBinResY = hResponse[iHist] -> GetNbinsY();
      for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
        for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
          const Double_t binVal = hResponse[iHist] -> GetBinContent(iBinResX, iBinResY);
          const Double_t binErr = hResponse[iHist] -> GetBinError(iBinResX, iBinResY);
          const Double_t xWidth = hResponse[iHist] -> GetXaxis() -> GetBinWidth(iBinResX);
          const Double_t yWidth = hResponse[iHist] -> GetYaxis() -> GetBinWidth(iBinResY);
          const Double_t dArea  = xWidth * yWidth;
          hResponse[iHist] -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
          hResponse[iHist] -> SetBinError(iBinResX, iBinResY, binErr / dArea);
        }
      }
      const UInt_t nBinResX = pResponse[iHist] -> GetNbinsX();
      const UInt_t nBinResY = pResponse[iHist] -> GetNbinsY();
      for (UInt_t iBinResX = 1; iBinResX < (nBinResX + 1); iBinResX++) {
        for (UInt_t iBinResY = 1; iBinResY < (nBinResY + 1); iBinResY++) {
          const Double_t binVal = pResponse[iHist] -> GetBinContent(iBinResX, iBinResY);
          const Double_t binErr = pResponse[iHist] -> GetBinError(iBinResX, iBinResY);
          const Double_t xWidth = pResponse[iHist] -> GetXaxis() -> GetBinWidth(iBinResX);
          const Double_t yWidth = pResponse[iHist] -> GetYaxis() -> GetBinWidth(iBinResY);
          const Double_t dArea  = xWidth * yWidth;
          pResponse[iHist] -> SetBinContent(iBinResX, iBinResY, binVal / dArea);
          pResponse[iHist] -> SetBinError(iBinResX, iBinResY, binErr / dArea);
        }
      }
    }  // end  loop
  } else {
    const Double_t parWidth = hJets[0][0] -> GetBinWidth(17);
    const Double_t detWidth = hJets[0][1] -> GetBinWidth(17);
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
      hJets[iHist][0]   -> Scale(1. / parWidth);
      hJets[iHist][1]   -> Scale(1. / detWidth);
      hJetEff[iHist][0] -> Scale(1. / parWidth);
      hJetEff[iHist][1] -> Scale(1. / detWidth);
    }  
  }
  cout << "    Scaled by bin width." << endl;

  // determine what to normalize with
  Double_t Norms[NHist];
  TH1D     *hNorm[NHist];
  if (UseRootNTrgs) {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
      hNorm[iHist] = (TH1D*) fIn[iHist] -> Get(SNorm.Data());
      Norms[iHist] = hNorm[iHist] -> GetEntries();
    }
  } else {
    for (UInt_t iHist = 0; iHist < NHist; iHist++) {
      if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
      Norms[iHist] = NTrgs[iHist];
    }
  }
  cout << "    Determined normalization factors." << endl;

  // scale histograms
  Double_t normAll(0.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
    hJets[iHist][0]   -> Scale(Weights[iHist]);
    hJets[iHist][1]   -> Scale(Weights[iHist]);
    hJetEff[iHist][0] -> Scale(Weights[iHist]);
    hJetEff[iHist][1] -> Scale(Weights[iHist]);
    hResponse[iHist]  -> Scale(Weights[iHist]);
    pResponse[iHist]  -> Scale(Weights[iHist]);
    normAll += Norms[iHist] * Weights[iHist];
  }
  cout << "    Scaled histograms:\n"
       << "      norm = " << normAll
       << endl;

  // for sums
  TH1D     *hJetSumAll[2];
  TH1D     *hEffSumAll[2];
  TH2D     *hResSumAll;
  TProfile *pResSumAll;
  // initialize All sums
  hJetSumAll[0] = (TH1D*)     hJets[0][0]   -> Clone();
  hJetSumAll[1] = (TH1D*)     hJets[0][1]   -> Clone();
  hEffSumAll[0] = (TH1D*)     hJetEff[0][0] -> Clone();
  hEffSumAll[1] = (TH1D*)     hJetEff[0][1] -> Clone();
  hResSumAll    = (TH2D*)     hResponse[0]  -> Clone();
  pResSumAll    = (TProfile*) pResponse[0]  -> Clone();
  hJetSumAll[0] -> SetName("hSumPar");
  hJetSumAll[1] -> SetName("hSumDet");
  hEffSumAll[0] -> SetName("hEffPar");
  hEffSumAll[1] -> SetName("hEffDet");
  hResSumAll    -> SetName("hResponse");
  pResSumAll    -> SetName("pResponse");
  hJetSumAll[0] -> Reset("ICE");
  hJetSumAll[1] -> Reset("ICE");
  hEffSumAll[0] -> Reset("ICE");
  hEffSumAll[1] -> Reset("ICE");
  hResSumAll    -> Reset("ICE");
  pResSumAll    -> Reset("ICE");
  // sum histograms
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
    hJetSumAll[0] -> Add(hJets[iHist][0]);
    hJetSumAll[1] -> Add(hJets[iHist][1]);
    hEffSumAll[0] -> Add(hJetEff[iHist][0]);
    hEffSumAll[1] -> Add(hJetEff[iHist][1]);
    hResSumAll    -> Add(hResponse[iHist]);
    pResSumAll    -> Add(pResponse[iHist]);
  }
  cout << "    Summed histograms." << endl;

  // for normalization
  const Double_t etaNorm    = 2. * (1. - JetResParameter);
  const Double_t jetNormAll = normAll * etaNorm;

  // normalize sums
  hJetSumAll[0] -> Scale(1. / jetNormAll);
  hJetSumAll[1] -> Scale(1. / jetNormAll);
  hEffSumAll[0] -> Scale(1. / jetNormAll);
  hEffSumAll[1] -> Scale(1. / jetNormAll);

  // normalize response matrices
  const UInt_t nXbinsAll = hResSumAll -> GetNbinsX();
  const UInt_t nYbinsAll = hResSumAll -> GetNbinsY();
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
  TH1D *hEffAll = (TH1D*) hEffSumAll[0] -> Clone();
  hEffAll -> SetName("hEfficiencyAll");
  hEffAll -> Reset();
  hEffAll -> Divide(hEffSumAll[1], hEffSumAll[0], 1., 1.);
  cout << "    Calculated efficiencies." << endl;

  // fit efficiencies
  const TString sFit("[0]*(1.-TMath::Exp(-1.*[1]*x))");
  const Float_t guess0(0.9);
  const Float_t guess1(5.); 
  const Float_t xyFit[2] = {0.2, 8.0};

  TF1 *fFitAll = new TF1("fEffAll", sFit.Data(), xyFit[0], xyFit[1]);
  fFitAll -> SetParameter(0, guess0);
  fFitAll -> SetParameter(1, guess1);
  hEffAll -> Fit(fFitAll, "R0");
  hEffAll -> GetFunction("fEffAll") -> ResetBit(1<<9);
  cout << "    Fit efficiencies." << endl;

  // close files
  fOut          -> cd();
  hJetSumAll[0] -> Write();
  hJetSumAll[1] -> Write();
  hEffSumAll[0] -> Write();
  hEffSumAll[1] -> Write();
  hEffAll       -> Write();
  hResSumAll    -> Write();
  pResSumAll    -> Write();
  fOut          -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (iHist == 3) continue;  // DELETE WHEN 20 - 25 GEV/C IS READY [02.20.2023]
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Preparations finished!\n" << endl;

}

// Ending ---------------------------------------------------------------------


