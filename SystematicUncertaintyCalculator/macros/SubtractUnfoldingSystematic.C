// 'SubtractUnfoldingSystematic.C'
// Derek Anderson
// 02.18.2022
//
// Use this to subtract an unfolding systematic
// uncertainty from the combined purity-plus-
// unfolding systematic uncertainty.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"

// global constants
static const UInt_t NIn(2);
static const UInt_t NTrg(3);
static const UInt_t NHist(NIn + 1);

// trigger parameters
static const UInt_t TrgBin(2);



void SubtractUnfoldingSystematic() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning unfolding systematic subtraction script..." << endl;

  // io parameters
  const TString sOutput("unfoldErrorSub_unfSubFromTot_withCorrectSmoothAndInterpol.et1520r02gam.d19m5y2022.root");
  const TString sInput[NIn]    = {"et1520r02gam_puritySysCalc/unfoldSys.noBkgdVariation_pTbinHuge.et1520r02gam.d12m2y2022.root",
                                  "et1520r02gam_puritySysCalc/unfoldSys.allBkgdAndUnfold_pTbinHuge.et1520r02gam.d27m4y2022.root"};
  const TString sInHist[NIn]   = {"hTotal", "hTotal"};
  const TString sOutTxt[NHist] = {"unfoldErrorSub_unfSubFromTot_withCorrectSmoothAndInterpol.onlyPur.et1520r02gam.d19m5y2022.txt",
                                  "unfoldErrorSub_unfSubFromTot_withCorrectSmoothAndInterpol.unfAndPur.et1520r02gam.d19m5y2022.txt",
                                  "unfoldErrorSub_unfSubFromTot_withCorrectSmoothAndInterpol.onlyUnf.et1520r02gam.d19m5y2022.txt"};

  // hist parameters
  const TString sName[NHist]    = {"hUnfoldOnly", "hUnfoldAndPurity", "hPurityOnly"};
  const TString sUnity[NHist]   = {"hUnityUnfold", "hUnityTotal", "hUnityPurity"};
  const TString sSmooth[NIn]    = {"hSmoothUnfold", "hSmoothTotal"};
  const Float_t xStartTrg[NTrg] = {0., 0., 0.};
  const Float_t xInterpol[NTrg] = {4., 6., 11.};
  const Bool_t  smoothHist[NIn] = {true, true};

  // parse trigger selection
  Float_t xStart(0.);
  Float_t xIntBegin(0.);
  switch (TrgBin) {
    case 0:
      xStart    = xStartTrg[0];
      xIntBegin = xInterpol[0];
      break;
    case 1:
      xStart    = xStartTrg[1];
      xIntBegin = xInterpol[1];
      break;
    case 2:
      xStart    = xStartTrg[2];
      xIntBegin = xInterpol[2];
      break;
    default:
      xStart    = xStartTrg[0];
      xIntBegin = xInterpol[0];
      break;
  }

  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fUnfold = new TFile(sInput[0].Data(), "read");
  TFile *fTotal  = new TFile(sInput[1].Data(), "read");
  if (!fOutput || !fUnfold || !fTotal) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOutput = " << fOutput << "\n"
         << "       fUnfold = " << fUnfold << "\n"
         << "       fTotal  = " << fTotal  << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hUnfold = (TH1D*) fUnfold -> Get(sInHist[0].Data());
  TH1D *hTotal  = (TH1D*) fTotal  -> Get(sInHist[1].Data());
  if (!hUnfold || !hTotal) {
    cerr << "PANIC: couldn't grab an input histogram!\n"
         << "       hUnfold = " << hUnfold << ", hTotal = " << hTotal << "\n"
         << endl;
    return;
  }
  hUnfold -> SetName(sName[0].Data());
  hTotal  -> SetName(sName[1].Data());
  cout << "    Grabbed histograms." << endl;

  // determine calculation range
  const Int_t  iXstart = hUnfold -> FindBin(xStart);
  const UInt_t iFirst  = hUnfold -> FindFirstBinAbove(0.);
  const UInt_t iLast   = hUnfold -> FindLastBinAbove(0.);
  const UInt_t iStop   = hUnfold -> GetNbinsX();

  // set starting bin
  UInt_t iStart(0);
  if (iFirst < iXstart) {
    iStart = iXstart;
  } else {
    iStart = iFirst;
  }

  // smooth systematics
  TH1D *hSmoothUnfold;
  if (smoothHist[0]) {

    // initialize smoothed histogram
    hSmoothUnfold = (TH1D*) hUnfold -> Clone();
    hSmoothUnfold -> SetName(sSmooth[0]);

    // do smoothing
    const UInt_t nUnfold = hUnfold -> GetNbinsX();
    for (UInt_t iUnfold = 2; iUnfold < nUnfold; iUnfold++) {

      // get bin info
      const Double_t thisBinLoc = hUnfold -> GetBinCenter(iUnfold);
      const Double_t thisBinVal = hUnfold -> GetBinContent(iUnfold);
      const Double_t prevBinVal = hUnfold -> GetBinContent(iUnfold - 1);
      const Double_t nextBinVal = hUnfold -> GetBinContent(iUnfold + 1);
      const Double_t thisBinAbs = hUnfold -> GetBinError(iUnfold);
      const Double_t prevBinAbs = hUnfold -> GetBinError(iUnfold - 1);
      const Double_t nextBinAbs = hUnfold -> GetBinError(iUnfold + 1);
      const Double_t thisBinErr = thisBinAbs / thisBinVal;
      const Double_t prevBinErr = prevBinAbs / prevBinVal;
      const Double_t nextBinErr = nextBinAbs / nextBinVal;

      // check bins
      const Bool_t areBinsNonzero = ((thisBinVal > 0.) && (nextBinVal > 0.));
      const Bool_t isInSmoothZone = (thisBinLoc <= xIntBegin); 
      const Bool_t isAboveStart   = (iUnfold >= iStart);
      const Bool_t isBelowStop    = ((iUnfold <= iStop) && (iUnfold <= iLast));
      const Bool_t isInCalcRange  = (isAboveStart && isBelowStop);
      if (!areBinsNonzero || !isInCalcRange) continue;

      // check how sys compare
      const Bool_t prevMoreThanThis = (prevBinErr > thisBinErr);
      const Bool_t prevMoreThanNext = (prevBinErr > nextBinErr);
      const Bool_t thisMoreThanBoth = ((thisBinErr > prevBinErr) && (thisBinErr > nextBinErr));

      // if (x < xInterpol) smooth; else, interpolate
      if (isInSmoothZone) {
        if (prevMoreThanThis) {
          const Double_t newThisErr = thisBinVal * prevBinErr;
          hSmoothUnfold -> SetBinError(iUnfold, newThisErr);
        }
      } else {
        if (thisMoreThanBoth) {
          const Double_t newNextErr = nextBinVal * thisBinErr;
          hSmoothUnfold -> SetBinError(iUnfold + 1, newNextErr); 
        } else if (prevMoreThanNext) {
          const Double_t newThisErr = thisBinVal * prevBinErr;
          const Double_t newNextErr = nextBinVal * prevBinErr;
          hSmoothUnfold -> SetBinError(iUnfold, newThisErr);
          hSmoothUnfold -> SetBinError(iUnfold + 1, newNextErr);
        } else if (prevMoreThanThis) {
          const Double_t errorDiff  = nextBinErr - prevBinErr;
          const Double_t errorAdj   = errorDiff / 2.;
          const Double_t thisBinAdj = prevBinErr + errorAdj;
          const Double_t newThisErr = thisBinVal * thisBinAdj;
          hSmoothUnfold -> SetBinError(iUnfold, newThisErr);
        }
      }  // end smoothing/interpolating
    }  // end bin loop 
    cout << "    Smoothed unfolding histogram." << endl;
  }

  // smooth systematics
  TH1D *hSmoothTotal;
  if (smoothHist[1]) {

    // initialize smoothed histogram
    hSmoothTotal = (TH1D*) hTotal -> Clone();
    hSmoothTotal -> SetName(sSmooth[1]);

    // do smoothing
    const UInt_t nTotal = hTotal -> GetNbinsX();
    for (UInt_t iTotal = 2; iTotal < nTotal; iTotal++) {

      // get bin info
      const Double_t thisBinLoc = hTotal -> GetBinCenter(iTotal);
      const Double_t thisBinVal = hTotal -> GetBinContent(iTotal);
      const Double_t prevBinVal = hTotal -> GetBinContent(iTotal - 1);
      const Double_t nextBinVal = hTotal -> GetBinContent(iTotal + 1);
      const Double_t thisBinAbs = hTotal -> GetBinError(iTotal);
      const Double_t prevBinAbs = hTotal -> GetBinError(iTotal - 1);
      const Double_t nextBinAbs = hTotal -> GetBinError(iTotal + 1);
      const Double_t thisBinErr = thisBinAbs / thisBinVal;
      const Double_t prevBinErr = prevBinAbs / prevBinVal;
      const Double_t nextBinErr = nextBinAbs / nextBinVal;

      // check bins
      const Bool_t areBinsNonzero = ((thisBinVal > 0.) && (nextBinVal > 0.));
      const Bool_t isInSmoothZone = (thisBinLoc < xIntBegin); 
      const Bool_t isAboveStart   = (iTotal >= iStart);
      const Bool_t isBelowStop    = ((iTotal <= iStop) && (iTotal <= iLast));
      const Bool_t isInCalcRange  = (isAboveStart && isBelowStop);
      if (!areBinsNonzero || !isInCalcRange) continue;

      // check how sys compare
      const Bool_t prevMoreThanThis = (prevBinErr > thisBinErr);
      const Bool_t prevMoreThanNext = (prevBinErr > nextBinErr);
      const Bool_t thisMoreThanBoth = ((thisBinErr > prevBinErr) && (thisBinErr > nextBinErr));

      // if (x < xInterpol) smooth; else, interpolate
      if (isInSmoothZone) {
        if (prevMoreThanThis) {
          const Double_t newThisErr = thisBinVal * prevBinErr;
          hSmoothTotal -> SetBinError(iTotal, newThisErr);
        }
      } else {
        if (thisMoreThanBoth) {
          const Double_t newNextErr = nextBinVal * thisBinErr;
          hSmoothTotal -> SetBinError(iTotal + 1, newNextErr); 
        } else if (prevMoreThanNext) {
          const Double_t newThisErr = thisBinVal * prevBinErr;
          const Double_t newNextErr = nextBinVal * prevBinErr;
          hSmoothTotal -> SetBinError(iTotal, newThisErr);
          hSmoothTotal -> SetBinError(iTotal + 1, newNextErr);
        } else if (prevMoreThanThis) {
          const Double_t errorDiff  = nextBinErr - prevBinErr;
          const Double_t errorAdj   = errorDiff / 2.;
          const Double_t thisBinAdj = prevBinErr + errorAdj;
          const Double_t newThisErr = thisBinVal * thisBinAdj;
          hSmoothTotal -> SetBinError(iTotal, newThisErr);
        }
      }  // end smoothing/interpolating
    }  // end bin loop 
    cout << "    Smoothed unfolding-plus-purity histogram." << endl;
  }

  // do subtraction
  TH1D *hPurity      = (TH1D*) hUnfold -> Clone();
  TH1D *hUnityUnfold = (TH1D*) hUnfold -> Clone();
  TH1D *hUnityTotal  = (TH1D*) hTotal  -> Clone();
  TH1D *hUnityPurity = (TH1D*) hTotal  -> Clone();
  hPurity      -> SetName(sName[2].Data());
  hUnityUnfold -> SetName(sUnity[0].Data());
  hUnityTotal  -> SetName(sUnity[1].Data());
  hUnityPurity -> SetName(sUnity[2].Data());
  hPurity      -> Reset("ICES");
  hUnityUnfold -> Reset("ICES");
  hUnityTotal  -> Reset("ICES");
  hUnityPurity -> Reset("ICES");

  const UInt_t nBins = hUnfold -> GetNbinsX();
  for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {

    // get bin content
    Double_t valUnf(0.);
    Double_t errUnf(0.);
    if (smoothHist[0]) {
      valUnf = hSmoothUnfold -> GetBinContent(iBin);
      errUnf = hSmoothUnfold -> GetBinError(iBin);
    } else {
      valUnf = hUnfold -> GetBinContent(iBin);
      errUnf = hUnfold -> GetBinError(iBin);
    }
    const Double_t perUnf = errUnf / valUnf;

    Double_t valTot(0.);
    Double_t errTot(0.);
    if (smoothHist[1]) {
      valTot = hSmoothTotal -> GetBinContent(iBin);
      errTot = hSmoothTotal -> GetBinError(iBin);
    } else {
      valTot = hTotal -> GetBinContent(iBin);
      errTot = hTotal -> GetBinError(iBin);
    }
    const Double_t perTot = errTot / valTot;

    // check if bins are good
    const Bool_t areBinsNonzero   = ((valUnf > 0.) && (valTot > 0.));
    const Bool_t isUnfLessThanTot = (perUnf < perTot);
    const Bool_t isAboveStart     = (iBin >= iStart);
    const Bool_t isBelowStop      = ((iBin <= iStop) && (iBin <= iLast));
    const Bool_t isInCalcRange    = (isAboveStart && isBelowStop);
    if (!areBinsNonzero || !isUnfLessThanTot || !isInCalcRange) continue;

    // calculate difference and set histograms
    const Double_t perDiff = TMath::Sqrt((perTot * perTot) - (perUnf * perUnf));
    const Double_t valDiff = valUnf;
    const Double_t errDiff = perDiff * valDiff;
    hPurity      -> SetBinContent(iBin, valUnf);
    hUnityUnfold -> SetBinContent(iBin, 1.);
    hUnityTotal  -> SetBinContent(iBin, 1.);
    hUnityPurity -> SetBinContent(iBin, 1.);
    hPurity      -> SetBinError(iBin, errDiff);
    hUnityUnfold -> SetBinError(iBin, perUnf);
    hUnityTotal  -> SetBinError(iBin, perTot);
    hUnityPurity -> SetBinError(iBin, perDiff);
  }
  cout << "    Did subtraction." << endl;

  // stream uncertainties to output
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {

    // create stream
    ofstream oTxt(sOutTxt[iHist].Data());
    oTxt << "binNumber binCenter percentErr" << endl;

    // loop over bins
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {

      // check if bin in range
      const Bool_t isAboveStart  = (iBin >= iStart);
      const Bool_t isBelowStop   = ((iBin <= iStop) && (iBin <= iLast));
      const Bool_t isInCalcRange = (isAboveStart && isBelowStop);
      if (!isInCalcRange) continue;

      // stream out values
      Double_t xBin(-1.);
      Double_t eBin(-1.);
      switch (iHist) {
        case 0:
          xBin = hUnityUnfold -> GetBinCenter(iBin);
          eBin = hUnityUnfold -> GetBinError(iBin);
          break;
        case 1:
          xBin = hUnityTotal  -> GetBinCenter(iBin);
          eBin = hUnityTotal  -> GetBinError(iBin);
          break;
        case 2:
          xBin = hUnityPurity -> GetBinCenter(iBin);
          eBin = hUnityPurity -> GetBinError(iBin);
          break;
        default:
          xBin = hUnityUnfold -> GetBinCenter(iBin);
          eBin = hUnityUnfold -> GetBinError(iBin);
          break;
      }
      oTxt << iBin;
      oTxt << " ";
      oTxt << xBin;
      oTxt << " ";
      oTxt << eBin;
      oTxt << endl;
    }
    oTxt.close();
  }
  cout << "    Streamed uncertainties." << endl;

  // save histograms
  fOutput      -> cd();
  hUnfold      -> Write();
  hTotal       -> Write();
  hPurity      -> Write();
  hUnityUnfold -> Write();
  hUnityTotal  -> Write();
  hUnityPurity -> Write();
  if (smoothHist[0]) {
    hSmoothUnfold -> Write();
  }
  if (smoothHist[1]) {
    hSmoothTotal -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fUnfold -> cd();
  fUnfold -> Close();
  fTotal  -> cd();
  fTotal  -> Close();
  cout << "  Unfolding systematic subtraction script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
