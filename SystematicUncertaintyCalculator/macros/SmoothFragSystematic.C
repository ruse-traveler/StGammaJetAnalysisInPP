// 'SmoothFragSystematic.C'
// Derek Anderson
// 01.07.2021
//
// Use this to smooth the fragmentation
// systematic uncertaintiy

#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"

using namespace std;

// global constants
static const UInt_t NQuad(2);
static const UInt_t NRange(2);



void SmoothFragSystematic() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning fragmentation smoothing script..." << endl;

  // io parameters
  const TString sIn("Ratio_Unfolded_data_11-15_5R.root");
  const TString sOut("smoothedFragSys_withFit.et1115r05gam.root");
  const TString sHistGet("hUnfolded");

  // histogram parameters
  const TString sHistIn("hInput");
  const TString sHistSys("hSmoothSys");
  const TString sHistSmooth("hSmooth");

  // calculation parameters
  const Bool_t   useFit(true);
  const UInt_t   nDegree(2);
  const Double_t dZero(0.01);
  const Double_t xSwitch(6.);
  const Double_t xSysRange[NRange] = {0., 30.};

  // function parameters
  const UInt_t   fColConst(1);
  const UInt_t   fColPol(2);
  const TString  sConst("fConstFit");
  const TString  sLine("fLinearSmooth");
  const TString  sQuad("fQuadraticSmooth");
  const TString  sConstFit("pol0");
  const TString  sLineFit("[0]*(x - [1]) + [2]");
  const TString  sQuadFit("[0]*(x*x - [1]*[1]) + [2]*(x - [1]) + [3]");
  const Double_t constGuess(1.);
  const Double_t lineGuess(1.);
  const Double_t quadGuess[NQuad]  = {-1., xSwitch};
  const Double_t xConstFit[NRange] = {xSwitch, 15.};
  const Double_t xPolFit[NRange]   = {0., xSwitch};

  // open files
  TFile *fIn  = new TFile(sIn.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn || !fOut) {
    cerr << "PANIC: couldn't open the input or output file!\n"
         << "       fIn = " << fIn << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input histogram
  TH1D *hIn = (TH1D*) fIn -> Get(sHistGet.Data());
  if (!hIn) {
    cerr << "PANIC: couldn't grab input histogram!\n" << endl;
    return;
  }
  hIn -> SetName(sHistIn.Data());
  cout << "    Grabbed input histogram." << endl;

  // create output histograms
  TH1D *hSys    = (TH1D*) hIn -> Clone();
  TH1D *hSmooth = (TH1D*) hIn -> Clone();
  hSys    -> SetName(sHistSys.Data());
  hSmooth -> SetName(sHistSmooth.Data());
  hSys    -> Reset("ICES");
  cout << "    Created histograms." << endl;

  // fit high pT range with constant
  TF1 *fConst = new TF1(sConst.Data(), sConstFit.Data(), xConstFit[0], xConstFit[1]);
  fConst  -> SetLineColor(fColConst);
  fConst  -> SetParameter(0, constGuess);
  hSmooth -> Fit(sConst.Data(), "R");

  // if close enough to 1, just use 1
  const Double_t yConst = fConst -> GetParameter(0);
  const Double_t dUnity = TMath::Abs(yConst - 1.);

  Double_t yUse(1.);
  Double_t eUse(0.);
  if (dUnity > dZero) {
    yUse = yConst;
    eUse = dUnity;
  }
  cout << "    Fit constant:\n"
       << "      fit constant  = " << yConst << "\n"
       << "      constant used = " << yUse
       << endl;

  // define smoothing function
  TF1 *fSmooth;
  switch (nDegree) {
    case 1:
      fSmooth = new TF1(sLine.Data(), sLineFit.Data(), xPolFit[0], xPolFit[1]);
      fSmooth -> SetLineColor(fColPol);
      fSmooth -> SetParameter(0, lineGuess);
      fSmooth -> FixParameter(1, xSwitch);
      fSmooth -> FixParameter(2, yUse);
      hSmooth -> Fit(sLine.Data(), "R+");
      break;
    case 2:
      fSmooth = new TF1(sQuad.Data(), sQuadFit.Data(), xPolFit[0], xPolFit[1]);
      fSmooth -> SetLineColor(fColPol);
      fSmooth -> SetParameter(0, quadGuess[0]);
      fSmooth -> SetParameter(2, quadGuess[2]);
      fSmooth -> FixParameter(1, xSwitch);
      fSmooth -> FixParameter(3, yUse);
      hSmooth -> Fit(sQuad.Data(), "R+");
      break;
    default:
      fSmooth = new TF1(sLine.Data(), sLineFit.Data(), xPolFit[0], xPolFit[1]);
      fSmooth -> SetLineColor(fColPol);
      fSmooth -> SetParameter(0, lineGuess);
      fSmooth -> FixParameter(1, xSwitch);
      fSmooth -> FixParameter(2, yUse);
      hSmooth -> Fit(sLine.Data(), "R+");
      break;
  }

  // smooth low pT range
  const UInt_t nBins = hIn -> GetNbinsX();
  if (!useFit) {
    for (UInt_t iBin = nBins; iBin > 0; iBin--) {

      const Double_t thisBinCenter = hSmooth -> GetBinCenter(iBin);
      const Double_t nextBinCenter = hSmooth -> GetBinCenter(iBin - 1);
      const Double_t prevBinCenter = hSmooth -> GetBinCenter(iBin + 1);
      const Double_t thisBinValue  = hSmooth -> GetBinContent(iBin);
      const Double_t nextBinValue  = hSmooth -> GetBinContent(iBin - 1);
      const Double_t prevBinValue  = hSmooth -> GetBinContent(iBin + 1);
      const Double_t thisBinError  = hSmooth -> GetBinError(iBin);
      const Double_t nextBinError  = hSmooth -> GetBinError(iBin - 1);
      const Double_t prevBinError  = hSmooth -> GetBinError(iBin + 1);
      const Double_t thisBinPer    = thisBinError / thisBinValue;
      const Double_t nextBinPer    = nextBinError / nextBinValue;
      const Double_t prevBinPer    = prevBinError / prevBinValue;

      // check if in range
      const Bool_t nextBinInRange = (nextBinCenter > xSysRange[0]);
      const Bool_t prevBinInRange = (prevBinCenter < xSwitch);
      const Bool_t thisBinInRange = ((thisBinCenter > xSysRange[0]) && (thisBinCenter < xSwitch));
      if (!nextBinInRange || !thisBinInRange) continue;

      // if at switch, use fit constant
      Double_t prevComp = prevBinValue;
      if (!prevBinInRange) prevComp = yUse;

      // check bin values
      const Bool_t prevLessThanThis = (prevComp < thisBinValue);
      const Bool_t prevLessThanNext = (prevComp < nextBinValue);
      const Bool_t thisLessThanBoth = ((thisBinValue < prevComp) && (thisBinValue < nextBinValue));

      if (thisLessThanBoth) {
        const Double_t newNextValue = thisBinValue;
        const Double_t newNextError = nextBinPer * newNextValue;
        hSmooth -> SetBinContent(iBin - 1, newNextValue);
        hSmooth -> SetBinError(iBin - 1, newNextError);
      } else if (prevLessThanNext) {
        const Double_t newThisValue = prevComp;
        const Double_t newNextValue = prevComp;
        const Double_t newThisError = thisBinPer * newThisValue;
        const Double_t newNextError = nextBinPer * newNextValue;
        hSmooth -> SetBinContent(iBin, newThisValue);
        hSmooth -> SetBinError(iBin, newThisError);
        hSmooth -> SetBinContent(iBin - 1, newNextValue);
        hSmooth -> SetBinError(iBin - 1, newNextError);
      } else if (prevLessThanThis) {
        const Double_t valueDiff    = prevComp - nextBinValue;
        const Double_t valueAdjust  = valueDiff / 2.;
        const Double_t newThisValue = thisBinValue - valueAdjust;
        const Double_t newThisError = newThisValue * thisBinPer;
        hSmooth -> SetBinContent(iBin, newThisValue);
        hSmooth -> SetBinError(iBin, newThisError);
      }
    }  // end bin loop
  }
  cout << "    Smoothed ratio." << endl;

  // set systematic uncertainties
  for (UInt_t iSys = 1; iSys < (nBins + 1); iSys++) {

    const Double_t binCenter    = hSys -> GetBinCenter(iSys);
    const Bool_t   isInSysRange = ((binCenter > xSysRange[0]) && (binCenter < xSysRange[1]));
    const Bool_t   isInSmoothed = (binCenter < xSwitch);
    if (!isInSysRange) continue;

    // determine systematic to use
    Double_t ySys = hSmooth -> GetBinContent(iSys);
    if (useFit) {
      ySys = fSmooth -> Eval(binCenter);
    }
    const Double_t eSys = TMath::Abs(1. - ySys);

    // set systemaic
    hSys -> SetBinContent(iSys, 1.);
    if (isInSmoothed) {
      hSys -> SetBinError(iSys, eSys);
    } else {
      hSys -> SetBinError(iSys, eUse);
    }
  }
  cout << "    Set systematic uncertainties." << endl;

  // save histograms
  fOut    -> cd();
  hIn     -> Write();
  hSys    -> Write();
  hSmooth -> Write();
  fConst  -> Write();
  fSmooth -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished smoothing fragmentation systematic!\n" << endl;

}

// End ------------------------------------------------------------------------
