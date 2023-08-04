// ApplySystematicsToParticleLevel.C
// Derek Anderson
// 08.04.2023
//
// Applies data systematic uncertainties
// from PRC/PRL to a specified jet
// spectrum.
//
// NOTE: assumes errors on input
// histogram are *just* statistical.
//
// Calculation parameters:
//   fType = 0: pi0
//         = 1: gamma-dir
//   fTrig = 0: 9 - 11
//         = 1: 11 - 15
//         = 2: 15 - 20
//   fReso = 0: R = 0.2
//         = 1: R = 0.5

#include <iostream>
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"

using namespace std;

// global constants
static const Ssiz_t NJets(2);
static const Ssiz_t NTrigs(2);
static const Ssiz_t NRange(2);
static const Ssiz_t NJetBins(6);
static const Ssiz_t NTrigBins(3);



void ApplySystematicsToParticleLevel() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Applying systematics..." << endl;

  // i/o parameters
  const TString sIn("input/embed/pp200r9embed.forClosureTest_pTbinGiant.et911r05qt05130.d3m8y2023.root");
  const TString sOut("summedErrorsFF.forClosureTest_sysFromTable_pTbinGiant.et911r05pi0.d4m8y2023.root");
  const TString sInHist("hSumParFF");

  // calculation parameters
  const UInt_t fType(0);
  const UInt_t fTrig(0);
  const UInt_t fReso(1);
  const Bool_t addUnfoldSys(true);
  const Bool_t addFragSys(false);
  const Bool_t addDetSys(true);

  // jet pt bins
  const Double_t ptJetBins[NJetBins][NRange] = {
    {0.,  5.},
    {5.,  10.},
    {10., 15.},
    {15., 20.},
    {20., 25.},
    {25., 30.}
  };

  // pi0 percent uncertainties (from table 5.3 of pp analysis note)
  const Double_t perSysUnfoldPi0[NTrigBins][NJets][NJetBins] = {
    {{0.01, 0.02, 0.02, 0.04, 0.22, 0.22}, {0.04, 0.04, 0.04, 0.07, 0.20, 0.54}},
    {{0.01, 0.01, 0.02, 0.03, 0.20, 0.57}, {0.05, 0.05, 0.05, 0.09, 0.14, 0.40}},
    {{0.01, 0.02, 0.04, 0.10, 0.48, 0.50}, {0.04, 0.09, 0.09, 0.16, 0.22, 0.53}}
  };
  const Double_t perSysFragPi0[NTrigBins][NJets][NJetBins] = {
    {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {0.30, 0.06, 0.00, 0.00, 0.00, 0.00}},
    {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {0.30, 0.06, 0.00, 0.00, 0.00, 0.00}},
    {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {0.30, 0.06, 0.00, 0.00, 0.00, 0.00}}
  };
  const Double_t perSysDetPi0[NTrigBins][NJets][NJetBins] = {
    {{0.07, 0.07, 0.09, 0.12, 0.10, 0.14}, {0.06, 0.06, 0.10, 0.14, 0.15, 0.11}},
    {{0.06, 0.07, 0.08, 0.10, 0.12, 0.10}, {0.07, 0.05, 0.09, 0.11, 0.14, 0.16}},
    {{0.07, 0.06, 0.07, 0.12, 0.02, 0.10}, {0.07, 0.07, 0.09, 0.08, 0.08, 0.13}}
  };

  // gamma percent uncertainties (from table 5.4 of pp analysis note)
  const Double_t perSysUnfoldGam[NTrigBins][NJets][NJetBins] = {
    {{0.06, 0.80, 0.00, 0.00, 0.00, 0.00}, {0.06, 0.82, 0.00, 0.00, 0.00, 0.00}},
    {{0.01, 0.21, 0.90, 0.00, 0.00, 0.00}, {0.06, 0.08, 0.20, 0.00, 0.00, 0.00}},
    {{0.07, 0.13, 0.16, 0.20, 0.00, 0.00}, {0.11, 0.12, 0.16, 0.16, 0.00, 0.00}}
  };
  const Double_t perSysFragGam[NTrigBins][NJets][NJetBins] = {
    {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {0.30, 0.06, 0.00, 0.00, 0.00, 0.00}},
    {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {0.30, 0.06, 0.00, 0.00, 0.00, 0.00}},
    {{0.00, 0.00, 0.00, 0.00, 0.00, 0.00}, {0.30, 0.06, 0.00, 0.00, 0.00, 0.00}}
  };
  const Double_t perSysDetGam[NTrigBins][NJets][NJetBins] = {
    {{0.07, 0.11, 0.00, 0.00, 0.00, 0.00}, {0.06, 0.14, 0.00, 0.00, 0.00, 0.00}},
    {{0.06, 0.09, 0.17, 0.00, 0.00, 0.00}, {0.06, 0.07, 0.12, 0.00, 0.00, 0.00}},
    {{0.04, 0.03, 0.02, 0.08, 0.00, 0.00}, {0.06, 0.04, 0.03, 0.02, 0.00, 0.00}}
  };

  // parse configuration (guard against bad options)
  Ssiz_t iTrigUse  = 0;
  Ssiz_t iResoUse  = 0;
  Bool_t isTrigPi0 = true;
  switch (fType) {
    case 0:
      isTrigPi0 = true;
      break;
    case 1:
      isTrigPi0 = false;
      break;
    default:
      isTrigPi0 = true;
      break;
  }
  switch (fTrig) {
    case 0:
      iTrigUse = 0;
      break;
    case 1:
      iTrigUse = 1;
      break;
    case 2:
      iTrigUse = 2;
      break;
    default:
      iTrigUse = 0;
      break;
  }
  switch (fReso) {
    case 0:
      iResoUse = 0;
      break;
    case 1:
      iResoUse = 1;
      break;
    default:
      isReso = 0;
      break;
  }
  cout << "    Parsed configuration." << endl;

  // open files
  TFile *fIn  = new TFile(sIn.Data(),  "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fIn || !fOut) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fIn = " << fIn << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input histograms
  TH1D *hInJet = (TH1D*) fIn -> Get(sInHist.Data());
  if (!hInJet) {
    cerr << "PANIC: couldn't grab input histogram!\n" << endl;
    return;
  }
  hInJet -> SetNameTitle("hStatistics", "");
  cout << "    Grabbed input jet spectrum." << endl;

  // create output histograms
  TH1D *hOutDet    = (TH1D*) hInJet -> Clone();
  TH1D *hOutFrag   = (TH1D*) hInJet -> Clone();
  TH1D *hOutUnfold = (TH1D*) hInJet -> Clone();
  TH1D *hOutTotal  = (TH1D*) hInJet -> Clone();
  TH1D *hOutCombo  = (TH1D*) hInJet -> Clone();
  hOutDet    -> SetNameTitle("hSysDet",         "");
  hOutFrag   -> SetNameTitle("hSysFrag",        "");
  hOutUnfold -> SetNameTitle("hSysUnfold",      "");
  hOutTotal  -> SetNameTitle("hSysTotal",       "");
  hOutCombo  -> SetNameTitle("hStatPlusTotSys", "");
  hOutDet    -> Reset("ICES");
  hOutFrag   -> Reset("ICES");
  hOutUnfold -> Reset("ICES");
  hOutTotal  -> Reset("ICES");
  hOutCombo  -> Reset("ICES");
  cout << "    Create output histograms.\n"
       << "    Now applying uncertainties..."
       << endl;

  // sum and apply uncertainties
  const Ssiz_t nBins = hInJet -> GetNbinsX();
  for (Ssiz_t iBin = 1; iBin < (nBins + 1); iBin++) {

    // get bin center, yield, and stat error
    const Double_t binCenter  = hInJet -> GetBinCenter(iBin);
    const Double_t binValue   = hInJet -> GetBinContent(iBin);
    const Double_t binStat    = hInJet -> GetBinError(iBin);
    const Double_t binPerStat = binStat / binValue;

    // determine which jet pt bin to use
    Ssiz_t iJetBinToUse = 0;
    Bool_t isInPtRange  = false;
    for (Ssiz_t iJetBin = 0; iJetBin < NJetBins; iJetBin++) {
      if ((binCenter >= ptJetBins[iJetBin][0]) && (binCenter < ptJetBins[iJetBin][1])) {
        iJetBinToUse = iJetBin;
        isInPtRange  = true;
        break;
      }
    }  // end jet pt bin loop

    // if not in range or yield is zero, skip
    const Bool_t isBinNonzero = (binValue > 0.);
    if (!isInPtRange || !isBinNonzero) continue;

    // grab uncertainties
    Double_t perSysUnfoldUse(0.);
    Double_t perSysFragUse(0.);
    Double_t perSysDetUse(0.);
    if (isTrigPi0) {
      perSysUnfoldUse = perSysUnfoldPi0[iTrigUse][iResoUse][iJetBinToUse];
      perSysFragUse   = perSysFragPi0[iTrigUse][iResoUse][iJetBinToUse];
      perSysDetUse    = perSysDetPi0[iTrigUse][iResoUse][iJetBinToUse];
    } else {
      perSysUnfoldUse = perSysUnfoldGam[iTrigUse][iResoUse][iJetBinToUse];
      perSysFragUse   = perSysFragGam[iTrigUse][iResoUse][iJetBinToUse];
      perSysDetUse    = perSysDetGam[iTrigUse][iResoUse][iJetBinToUse];
    }

    // calculate totals
    Double_t perSysTotSquare(0.);
    if (addUnfoldSys) perSysTotSquare += (perSysUnfoldUse * perSysUnfoldUse);
    if (addFragSys)   perSysTotSquare += (perSysFragUse   * perSysFragUse);
    if (addDetSys)    perSysTotSquare += (perSysDetUse    * perSysDetUse);

    const Double_t perSysTot   = TMath::Sqrt(perSysTotSquare);
    const Double_t statPlusTot = TMath::Sqrt((perSysTot * perSysTot) + (binPerStat * binPerStat));

    // set output histograms
    const Double_t sysUnfoldUse = binValue * perSysUnfoldUse;
    const Double_t sysFragUse   = binValue * perSysFragUse;
    const Double_t sysDetUse    = binValue * perSysDetUse;
    const Double_t sysTotUse    = binValue * perSysTot;
    const Double_t sysComboUse  = binValue * statPlusTot;
    if (addUnfoldSys) {
      hOutUnfold -> SetBinContent(iBin, binValue);
      hOutUnfold -> SetBinError(iBin, sysUnfoldUse);
    }
    if (addFragSys) {
      hOutFrag -> SetBinContent(iBin, binValue);
      hOutFrag -> SetBinError(iBin, sysFragUse);
    }
    if (addDetSys) {
      hOutDet -> SetBinContent(iBin, binValue);
      hOutDet -> SetBinError(iBin, sysDetUse);
    }
    hOutTotal -> SetBinContent(iBin, binValue);
    hOutCombo -> SetBinContent(iBin, binValue);
    hOutTotal -> SetBinError(iBin, sysTotUse);
    hOutCombo -> SetBinError(iBin, sysComboUse);
  }  // end histogram bin loop
  cout << "    Finished applying uncertainties." << endl;

  // save histograms
  fOut       -> cd();
  hInJet     -> Write();
  hOutDet    -> Write();
  hOutFrag   -> Write();
  hOutUnfold -> Write();
  hOutTotal  -> Write();
  hOutCombo  -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  fIn  -> cd();
  fIn  -> Close();
  cout << "  Finished applying systematics!\n" << endl;

}

// End ------------------------------------------------------------------------
