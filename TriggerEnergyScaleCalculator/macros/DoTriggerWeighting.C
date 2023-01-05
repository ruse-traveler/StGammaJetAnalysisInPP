// 'DoTriggerWeighting.C'
// Derek Anderson
// 03.07.2022
//
// Use this apply a TES/R matrix
// to an input trigger specturm.

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"

using namespace std;

// global constants
static const UInt_t NTrg(3);
static const UInt_t NTot(11);
static const UInt_t NBin[NTrg] = {2, 4, 5};


void DoTriggerWeighting() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning TES/R matrix application..." << endl;

  // io parameters
  const TString sOutput("weightedPythia8.et650x920pi0.d7m3y2022.root");
  const TString sInTrg("pp200py8par.forEtProjCheck_pTbinOne.et0100r02pi0.d8m3y2022.root");
  const TString sInMat("output/March2022/triggerMatrix.withWeights.et650x920vz55tsp008pi0.d6m3y2022.root");
  const TString sHistTrg("Pi0/hTrgEtP");
  const TString sHistMat("hWeightMatrix");

  // histogram parameters
  const TString sOutTrg("hInputTrigger");
  const TString sOutMat("hInputMatrix");
  const TString sProjPref("hProjection");
  const TString sApplyPref("hAppliedProj");
  const TString sWeightPref("hWeightedTrigEt");
  const TString sWeightSuff[NTrg] = {"911", "1115", "1520"};

  // calculation parameters
  const Bool_t   doTrigNorm(true);
  const Double_t xLowEnd[NTrg]  = {9., 11., 15.};
  const Double_t xHighEnd[NTrg] = {11., 15., 20.};

  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInTrg  = new TFile(sInTrg.Data(), "read");
  TFile *fInMat  = new TFile(sInMat.Data(), "read");
  if (!fOutput || !fInTrg || !fInMat) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOutput = " << fOutput << ", fInTrg = " << fInTrg << ", fInMat = " << fInMat << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hTrig   = (TH1D*) fInTrg -> Get(sHistTrg.Data());
  TH2D *hMatrix = (TH2D*) fInMat -> Get(sHistMat.Data());
  if (!hTrig || !hMatrix) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hTrig = " << hTrig << ", hMatrix = " << hMatrix << "\n"
         << endl;
    return;
  }
  hTrig   -> SetName(sOutTrg.Data());
  hMatrix -> SetName(sOutMat.Data());
  cout << "    Grabbed histograms." << endl;

  // normalize input
  if (doTrigNorm) {
    const Double_t intTrig = hTrig -> Integral();
    hTrig -> Scale(1. / intTrig);
    cout << "    Normalized input trigger spectrum." << endl;
  }

  // check number of bins
  UInt_t nCheck(0);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    nCheck += NBin[iTrg];
  }

  if (nCheck == NTot) {
    cout << "    Preparing for matrix application. Number of bins looks good..." << endl;
  } else {
    cerr << "PANIC: check NTrg, NTot, and NBin! Entries of NBin should add up to NTot!\n"
         << "       NTrg = " << NTrg << ", NTot = " << NTot << ", nCheck = " << nCheck << "\n"
         << endl;
    return;
  }

  // project and apply matrix
  UInt_t  iTot(0);
  UInt_t  nBinsX(0);
  UInt_t  nBinsY(0);
  UInt_t  nBinsP(0);
  UInt_t  nBin[NTrg];
  TH1D   *hProj[NTot];
  TH1D   *hApply[NTot];
  TH1D   *hWeight[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {

    // initialize bin counter
    nBinsX     = (UInt_t) hMatrix -> GetNbinsX();
    nBinsY     = (UInt_t) hTrig   -> GetNbinsX();
    nBin[iTrg] = 0;

    // create weighted name
    TString sNameWeight(sWeightPref.Data());
    sNameWeight.Append(sWeightSuff[iTrg].Data());

    // create weighted histogram
    hWeight[iTrg] = (TH1D*) hTrig -> Clone();
    hWeight[iTrg] -> SetName(sNameWeight.Data());
    hWeight[iTrg] -> Reset("ICES");

    // determine bin range
    for (UInt_t iBinX = 1; iBinX < (nBinsX + 1); iBinX++) {

      const Double_t xBinCenter  = hMatrix -> GetXaxis() -> GetBinCenter(iBinX);
      const Double_t xBinLowEdge = hMatrix -> GetXaxis() -> GetBinLowEdge(iBinX);
      const Bool_t   isInRange   = ((xBinCenter > xLowEnd[iTrg]) && (xBinCenter < xHighEnd[iTrg]));
      if (isInRange) {

        // create projection names
        TString sNameProj(sProjPref.Data());
        TString sNameApply(sApplyPref.Data());
        sNameProj  += iTot;
        sNameApply += iTot;

        // project matrix
        hProj[iTot]  = (TH1D*) hMatrix -> ProjectionY("", iBinX, iBinX, "") -> Clone();
        hApply[iTot] = (TH1D*) hTrig   -> Clone();
        hProj[iTot]  -> SetName(sNameProj.Data());
        hApply[iTot] -> SetName(sNameApply.Data());

        // apply projection
        nBinsP = hProj[iTot]   -> GetNbinsX();
        for (UInt_t iBinY = 1; iBinY < (nBinsY + 1); iBinY++) {
          const Double_t yBinP    = hApply[iTot] -> GetBinCenter(iBinY);
          const UInt_t   iBinP    = hProj[iTot]  -> FindBin(yBinP);
          const Bool_t   isInProj = ((iBinP > 0) && (iBinP < (nBinsP + 1)));
          if (isInProj) {
            const Double_t trigVal = hApply[iTot] -> GetBinContent(iBinY);
            const Double_t projVal = hProj[iTot]  -> GetBinContent(iBinP);
            const Double_t trigErr = hApply[iTot] -> GetBinError(iBinY);
            const Double_t projErr = hProj[iTot]  -> GetBinError(iBinP);
            const Double_t trigRel = trigErr / trigVal;
            const Double_t projRel = projErr / projVal;
            const Double_t multVal = trigVal * projVal;
            const Double_t multRel = TMath::Sqrt((trigRel * trigRel) + (projRel * projRel));
            const Double_t multErr = multVal * multRel;
            hApply[iTot] -> SetBinContent(iBinY, multVal);
            hApply[iTot] -> SetBinError(iBinY, multErr);
          } else {
            hApply[iTot] -> SetBinContent(iBinY, 0.);
            hApply[iTot] -> SetBinError(iBinY, 0.);
          }
        }  // end y-axis loop

        // add to sum and incrememnt counters
        hWeight[iTrg] -> Add(hApply[iTot]);
        nBin[iTrg]++;
        iTot++;
      }
    }  // end x-axis loop

    // check number of bins used
    if (NBin[iTrg] != nBin[iTrg]) {
      cerr << "PANIC: number of bins is off for trigger bin #" << iTrg << "!\n"
           << "       nBinUsed = " << nBin[iTrg] << ", nBinExpected = " << NBin[iTrg] << "\n"
           << endl;
      return;
    }

    // scale sum
    const Double_t scale = 1. / (Double_t) nBin[iTrg];
    if (scale > 0.) {
      hWeight[iTrg] -> Scale(scale);
    }
  }  // end trigger loop

  // announce details
  cout << "    Finished applying matrix:\n"
       << "      total number of bins used = " << iTot
       << endl;

  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    cout << "      number of bins used for bin " << iTrg << " = " << nBin[iTrg] << endl;
  }

  // write histograms
  fOutput -> cd();
  hTrig   -> Write();
  hMatrix -> Write();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hWeight[iTrg] -> Write();
  }
  for (UInt_t iTot = 0; iTot < NTot; iTot++) {
    hProj[iTot]  -> Write();
    hApply[iTot] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  fInTrg  -> cd();
  fInTrg  -> Close();
  fInMat  -> cd();
  fInMat  -> Close();
  cout << "  Finished TES/R matrix application!\n" << endl;

}

// End ------------------------------------------------------------------------
