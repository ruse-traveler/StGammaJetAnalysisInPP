// 'DoGammaSubtraction.C'
// Derek Anderson
// 08.17.2019
//
// Use this to subtract the pi0
// background from a gamma-rich
// spectrum.


#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace std;


// global constants
static const UInt_t  nTrgEtBins(3);
static const Float_t purityVal[nTrgEtBins] = {0.569875, 0.520370, 0.466597};
static const Float_t sysError[nTrgEtBins]  = {0.0537321, 0.0356863, 0.0677825};



void DoGammaSubtraction() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning subtraction script..." << endl;

  // parameters
  const UInt_t  iTrgEtBin(1);
  const Bool_t  doSysVariation(true);
  const Bool_t  isUpperVariation(true);
  const TString sFileOut("pp200r9data.forPuritySysPcheck_gamSubBefore_pTbinHuge.et1115r02chrg.d10m10y2021.root");
  const TString sFileGam("pp200r9data.forUnfolding_pTbinHuge.et1115r02chrg.d27m9y2021.root");
  const TString sFilePi0("pp200r9data.forUnfolding_pTbinHuge.et1115r02chrg.d27m9y2021.root");
  const TString sHistGam("Gam/hJetPtCorrG");
  const TString sHistPi0("Pi0/hJetPtCorrP");

  // select gamma purity
  Double_t purityToUse(0.);
  Double_t sysErrToUse(0.);
  switch (iTrgEtBin) {
    // 9 - 11 GeV
    case 0:
      purityToUse = purityVal[0];
      sysErrToUse = sysError[0];
      break;
    // 11 - 15 GeV
    case 1:
      purityToUse = purityVal[1];
      sysErrToUse = sysError[1];
      break;
    // 15 - 20 GeV
    case 2:
      purityToUse = purityVal[2];
      sysErrToUse = sysError[2];
      break;
    // default selection
    default:
      purityTouse = purityVal[0];
      sysErrToUse = sysError[0];
      break;
  }
  if (doSysVariation) {
    if (isUpperVariation)
      purityToUse += sysErrToUse; 
    else
      purityToUse -= sysErrToUse;
  }
  const Double_t gammaPurity = purityToUse;
  cout << "    Selected purity value: " << gammaPurity << endl;


  // open files
  TFile *fOut = new TFile(sFileOut.Data(), "recreate");
  TFile *fGam = new TFile(sFileGam.Data(), "read");
  TFile *fPi0 = new TFile(sFilePi0.Data(), "read");
  if (!fOut || !fGam || !fPi0) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fGam = " << fGam << ", fPi0 = " << fPi0
         << endl;
  }
  cout << "    Opened files." << endl;

  // grab files
  TH1D *hGam = (TH1D*) fGam -> Get(sHistGam.Data());
  TH1D *hPi0 = (TH1D*) fPi0 -> Get(sHistPi0.Data());
  if (!hGam || !hPi0) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hGam = " << hGam << ", hPi0 = " << hPi0
         << endl;
  }
  cout << "    Grabbed histograms." << endl;


  // do subtraction
  TH1D *hOut = (TH1D*) hGam -> Clone();
  TH1D *hSub = (TH1D*) hPi0 -> Clone();
  hOut -> Reset("ICE");
  hSub -> Scale(gammaPurity);
  hOut -> Add(hGam, hSub, 1., -1.);
  hOut -> Scale(1. / (1. - gammaPurity));
  cout << "    Did subtraction.\n"
       << "      purity = " << gammaPurity << ", scale = " << (1. / (1. - gammaPurity))
       << endl;

  // set negative bins to 0
  const UInt_t nBins = (UInt_t) hSub -> GetNbinsX();
  for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
    const Double_t content = hSub -> GetBinContent(iBin);
    if (content < 0.) {
      hSub -> SetBinContent(iBin, 0.);
      hSub -> SetBinError(iBin, 0.);
    }
  }
  cout << "    Set negative bins to 0." << endl;

  // set names
  const TString sOutName("hSub");
  const TString sGamName("hGam");
  const TString sPi0Name("hPi0");
  hOut -> SetName(sOutName.Data());
  hGam -> SetName(sGamName.Data());
  hPi0 -> SetName(sPi0Name.Data());
  cout << "    Set names.\n"
       << "      output histogram = " << sOutName.Data() << "\n"
       << "      gamma histogram  = " << sGamName.Data() << "\n"
       << "      pi0 histogram    = " << sPi0Name.Data()
       << endl;


  // close files
  fOut -> cd();
  hGam -> Write();
  hPi0 -> Write();
  hOut -> Write();
  fOut -> Close();
  fGam -> cd();
  fGam -> Close();
  fPi0 -> cd();
  fPi0 -> Close();
  cout << "  Finished subtraction script!\n" << endl;

}

// End ------------------------------------------------------------------------
