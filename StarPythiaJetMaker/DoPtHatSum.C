// 'DoPtHatSum.C'
// Nihar Sahoo, Saskia Mioduszewski, Derek Anderson
// 08.19.2022
//
//  Use this to sum a set of trigger and jet
//  histograms over pThat.
//
//  NOTE: use FSubSample to control which PYTHIA-6
//  sub-sample is being summed, and use UseSmeared
//  to select which histograms (eT-smeared or not)
//  are being summed.
//    FSubSample = 0:     sum eT > 6 sub-sample
//               = 1:     sum eT > 8 sub-sample
//    UseSmeared = true:  sum eT-smeared histograms
//               = false: sum non-eT-smeared histograms

#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"

using namespace std;

// global constants
static const UInt_t NPtHatBins(4);
static const UInt_t FSubSample(0);
static const Bool_t UseSmeared(false);

// x-sections and event fractions
static const Double_t SigPtHat[NPtHatBins]    = {0.46,      0.023,  0.0065, 0.0013};     // should not change
static const Double_t EvntFrcEt6a[NPtHatBins] = {0.0016,    0.0163, 0.0598, 0.1958};     // [6-25 GeV]
static const Double_t EvntFrcEt6b[NPtHatBins] = {0.0015,    0.0168, 0.0633, 0.2212};     // [6-25 GeV, |eta| < 0.9]
static const Double_t EvntFrcEt6c[NPtHatBins] = {0.00011,   0.00268, 0.00945, 0.04442};  // [6-25 GeV, |eta| < 0.9, only pi0]
static const Double_t EvntFrcEt8a[NPtHatBins] = {0.0001,    0.0017, 0.0092, 0.0771};     // [8-25 GeV]
static const Double_t EvntFrcEt8b[NPtHatBins] = {0.0000001, 0.0015, 0.0135, 0.0836};     // [8-25 GeV, |eta| < 0.9]
static const Double_t EvntFrcEt8c[NPtHatBins] = {1.,        1.,     1.,     1.};         // [8-25 GeV, |eta| < 0.9, only pi0]



void DoPtHatSum() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning pThat-summing script..." << endl;

  // output parameters
  const TString sOutPy6("pp200py6et6sum.test3.et625pi0.d23m10y2022.root");
  const TString sOutPy8("pp200py6et8sum.test3.et825pi0.d23m10y2022.root");

  // input parameters
  // >>> NEED TO CHANGE INPUT FILES <<<
  const TString sInPy6[NPtHatBins]   = {"./pp200py6et6pt47.forJetRerun_noCorruptTrees_wave1t2_pTbinOne.et6100pi0.r02rm1chrg.d1m10y2022.root",
                                        "./pp200py6et6pt79.forJetRerun_noCorruptTrees_wave1t4_pTbinOne.et6100pi0.r02rm1chrg.d5m10y2022.root",
                                        "./pp200py6et6pt912.forJetRerun_noCorruptTrees_run1t5_pTbinOne.et6100pi0.r02rm1chrg.d13m10y2022.root",
                                        "./pp200py6et6pt12.forJetRerun_noCorruptTrees_run1t10_pTbinOne.et6100pi0.r02rm1chrg.d1m10y2022.root"};
  const TString sInPy8[NPtHatBins]   = {"./pp200py6et8.forJetRerun_noCorruptTrees_pTbinOne.et6100r02pi0.d23m10y2022.root",
                                        "./pp200py6et8.forJetRerun_noCorruptTrees_pTbinOne.et6100r02pi0.d23m10y2022.root",
                                        "./pp200py6et8.forJetRerun_noCorruptTrees_pTbinOne.et6100r02pi0.d23m10y2022.root",
                                        "./pp200py6et8.forJetRerun_noCorruptTrees_pTbinOne.et6100r02pi0.d23m10y2022.root"};
  // >>> DO NOT NEED TO CHANGE ANYTHING BELOW <<<
  const TString sTrgPy6s[NPtHatBins] = {"Pi0/hTrgEtWeightP_pt47", "Pi0/hTrgEtWeightP_pt79", "Pi0/hTrgEtWeightP_pt912", "Pi0/hTrgEtWeightP_pt12"};
  const TString sTrgPy6n[NPtHatBins] = {"Pi0/hTrgEtP_pt47",       "Pi0/hTrgEtP_79",         "Pi0/hTrgEtP_pt912",       "Pi0/hTrgEtP_pt12"};
  const TString sTrgPy8s[NPtHatBins] = {"Pi0/hTrgEtWeightP_pt47", "Pi0/hTrgEtWeightP_pt79", "Pi0/hTrgEtWeightP_pt912", "Pi0/hTrgEtWeightP_pt12"};
  const TString sTrgPy8n[NPtHatBins] = {"Pi0/hTrgEtP_pt47",       "Pi0/hTrgEtP_pt79",       "Pi0/hTrgEtP_pt912",       "Pi0/hTrgEtP_pt12"};
  const TString sJetPy6s[NPtHatBins] = {"Pi0/hJetPtCorrP_pt47",   "Pi0/hJetPtCorrP_pt79",   "Pi0/hJetPtCorrP_pt912",   "Pi0/hJetPtCorrP_pt12"};
  const TString sJetPy6n[NPtHatBins] = {"Pi0/hJetPtCorrP_pt47",   "Pi0/hJetPtCorrP_pt79",   "Pi0/hJetPtCorrP_pt912",   "Pi0/hJetPtCorrP_pt12"};
  const TString sJetPy8s[NPtHatBins] = {"Pi0/hJetPtCorrP_pt47",   "Pi0/hJetPtCorrP_pt79",   "Pi0/hJetPtCorrP_pt912",   "Pi0/hJetPtCorrP_pt12"};
  const TString sJetPy8n[NPtHatBins] = {"Pi0/hJetPtCorrP_pt47",   "Pi0/hJetPtCorrP_pt79",   "Pi0/hJetPtCorrP_pt912",   "Pi0/hJetPtCorrP_pt12"};
  const TString sNormPy6[NPtHatBins] = {"Pi0/hChkPtHatP_pt47",    "Pi0/hChkPtHatP_pt79",    "Pi0/hChkPtHatP_pt912",    "Pi0/hChkPtHatP_pt12"};
  const TString sNormPy8[NPtHatBins] = {"Pi0/hChkPtHatP_pt47",    "Pi0/hChkPtHatP_pt79",    "Pi0/hChkPtHatP_pt912",    "Pi0/hChkPtHatP_pt12"};

  // pThat parameters
  const TString sTrgInput("hTrgInputPtHat");
  const TString sJetInput("hJetInputPtHat");
  const TString sTrgBinS("hSmearTrgPtHat");
  const TString sTrgBinN("hTrgPtHat");
  const TString sJetBinS("hSmearJetPtHat");
  const TString sJetBinN("hJetPtHat");
  const TString sNorm("hNormPtHat");
  const TString sProb("_prob");
  const TString sCross("_cross");
  const TString sPtHat[NPtHatBins] = {"47", "79", "912", "12"};

  // select event fractions
  Double_t evntFrc[NPtHatBins];
  switch (FSubSample) {
    case 0:
      cout << "    Using 6 - 25 GeV event fractions:" << endl;
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        evntFrc[iPtHatBin] = EvntFrcEt6c[iPtHatBin];
      }
      break;
    case 1:
      cout << "    Using 8 - 25 GeV event fractions:" << endl;
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        evntFrc[iPtHatBin] = EvntFrcEt8c[iPtHatBin];
      }
      break;
    default:
      cout << "    Using 6 - 25 GeV event fractions:" << endl;
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        evntFrc[iPtHatBin] = EvntFrcEt6c[iPtHatBin];
      }
      break;
  }
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    cout << "      EvntFrc[" << iPtHatBin << "] = " << evntFrc[iPtHatBin] << endl;
  }

  // select output file name
  TString sOut("");
  switch (FSubSample) {
    case 0:
      sOut = sOutPy6;
      break;
    case 1:
      sOut = sOutPy8;
      break;
    default:
      sOut = sOutPy6;
      break;
  }
  cout << "    Selected output name:\n"
       << "      sOut = " << sOut.Data()
       << endl;

  // create output file
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }
  cout << "    Created output file." << endl;

  // select input files
  TString sIn[NPtHatBins];
  switch (FSubSample) {
    case 0:
      cout << "    Selected eT > 6 input files:"  << endl;
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        sIn[iPtHatBin] = sInPy6[iPtHatBin];
      }
      break;
    case 1:
      cout << "    Selected eT > 8 input files:"  << endl;
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        sIn[iPtHatBin] = sInPy8[iPtHatBin];
      }
      break;
    case 2:
      cout << "    Selected eT > 6 input files:"  << endl;
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        sIn[iPtHatBin] = sInPy6[iPtHatBin];
      }
      break;
  }
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    cout << "      sIn[" << iPtHatBin << "] = " << sIn[iPtHatBin].Data() << endl;
  }

  // open input files
  TFile *fIn[NPtHatBins];
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    fIn[iPtHatBin] = new TFile(sIn[iPtHatBin].Data(), "read");
    if (!fIn[iPtHatBin]) {
      cerr << "PANIC: couldn't open input file #" << iPtHatBin "!\n" << endl;
      return;
    }
  }
  cout << "    Opened input files." << endl;

  // select input histogram name
  TString sInTrg[NPtHatBins];
  TString sInJet[NPtHatBins];
  TString sInNorm[NPtHatBins];
  switch (FSubSample) {
    case 0:
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        if (UseSmeared) {
          sInTrg[iPtHatBin] = sTrgPy6s[iPtHatBin];
          sInJet[iPtHatBin] = sJetPy6s[iPtHatBin];
        } else {
          sInTrg[iPtHatBin] = sTrgPy6n[iPtHatBin];
          sInJet[iPtHatBin] = sJetPy6n[iPtHatBin];
        }
        sInNorm[iPtHatBin] = sNormPy6[iPtHatBin];
      }
      break;
    case 1:
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        if (UseSmeared) {
          sInTrg[iPtHatBin] = sTrgPy8s[iPtHatBin];
          sInJet[iPtHatBin] = sJetPy8s[iPtHatBin];
        } else {
          sInTrg[iPtHatBin] = sTrgPy8n[iPtHatBin];
          sInJet[iPtHatBin] = sJetPy8n[iPtHatBin];
        }
        sInNorm[iPtHatBin] = sNormPy8[iPtHatBin];
      }
      break;
    default:
      for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
        if (UseSmeared) {
          sInTrg[iPtHatBin] = sTrgPy6s[iPtHatBin];
          sInJet[iPtHatBin] = sJetPy6s[iPtHatBin];
        } else {
          sInTrg[iPtHatBin] = sTrgPy6n[iPtHatBin];
          sInJet[iPtHatBin] = sJetPy6n[iPtHatBin];
        }
        sInNorm[iPtHatBin] = sNormPy6[iPtHatBin];
      }
      break;
  }
  cout << "    Selected input histograms:" << endl;
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    cout << "      sInTrg[" << iPtHatBin << "]  = " << sInTrg[iPtHatBin].Data() << "\n"
         << "      sInJet[" << iPtHatBin << "]  = " << sInJet[iPtHatBin].Data() << "\n"
         << "      sInNorm[" << iPtHatBin << "] = " << sInNorm[iPtHatBin].Data()
         << endl;
  }

  // grab input histograms
  TH1D *hTrgInput[NPtHatBins];
  TH1D *hJetInput[NPtHatBins];
  TH1D *hNorm[NPtHatBins];
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {

    // grab histograms
    hTrgInput[iPtHatBin] = (TH1D*) (fIn[iPtHatBin] -> Get(sInTrg[iPtHatBin].Data())) -> Clone();
    hJetInput[iPtHatBin] = (TH1D*) (fIn[iPtHatBin] -> Get(sInJet[iPtHatBin].Data())) -> Clone();
    hNorm[iPtHatBin]     = (TH1D*) (fIn[iPtHatBin] -> Get(sInTrg[iPtHatBin].Data())) -> Clone();
    if (!hTrgInput[iPtHatBin] || !hJetInput[iPtHatBin] || !hNorm[iPtHatBin]) {
      cerr << "PANIC: couldn't grab an input histogram!\n"
           << "       hTrgInput[" << iPtHatBin << "] = " << hTrgInput[iPtHatBin] << ", hJetInput[" << iPtHatBin << "] = " << hJetInput[iPtHatBin] << "\n"
           << "       hNorm[" << iPtHatBin << "]     = " << hNorm[iPtHatBin]
           << endl;
      return; 
    }

    // create names
    TString sTrgNameI(sTrgInput.Data());
    TString sJetNameI(sJetInput.Data());
    TString sNormName(sNorm.Data());
    sTrgNameI.Append(sPtHat[iPtHatBin]);
    sJetNameI.Append(sPtHat[iPtHatBin]);
    sNormName.Append(sPtHat[iPtHatBin]);
    hTrgInput[iPtHatBin] -> SetName(sTrgNameI.Data());
    hJetInput[iPtHatBin] -> SetName(sJetNameI.Data());
    hNorm[iPtHatBin]     -> SetName(sNormName.Data());
  }
  fOut -> cd();
  cout << "    Grabbed input histograms." << endl;

  // create histograms to sum
  TH1D *hTrgBinP[NPtHatBins];
  TH1D *hTrgBinX[NPtHatBins];
  TH1D *hJetBinP[NPtHatBins];
  TH1D *hJetBinX[NPtHatBins];
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {

    // clone histograms
    hTrgBinP[iPtHatBin] = (TH1D*) hTrgInput[iPtHatBin] -> Clone();
    hTrgBinX[iPtHatBin] = (TH1D*) hTrgInput[iPtHatBin] -> Clone();
    hJetBinP[iPtHatBin] = (TH1D*) hJetInput[iPtHatBin] -> Clone();
    hJetBinX[iPtHatBin] = (TH1D*) hJetInput[iPtHatBin] -> Clone();

    // create names
    TString sTrgNameP("");
    TString sTrgNameX("");
    TString sJetNameP("");
    TString sJetNameX("");
    if (UseSmeared) {
      sTrgNameP.Append(sTrgBinS);
      sTrgNameX.Append(sTrgBinS);
      sJetNameP.Append(sJetBinS);
      sJetNameX.Append(sJetBinS);
    } else {
      sTrgNameP.Append(sTrgBinN);
      sTrgNameX.Append(sTrgBinN);
      sJetNameP.Append(sJetBinN);
      sJetNameX.Append(sJetBinN);
    }
    sTrgNameP.Append(sPtHat[iPtHatBin]);
    sTrgNameX.Append(sPtHat[iPtHatBin]);
    sJetNameP.Append(sPtHat[iPtHatBin]);
    sJetNameX.Append(sPtHat[iPtHatBin]);
    sTrgNameP.Append(sProb.Data());
    sTrgNameX.Append(sCross.Data());
    sJetNameP.Append(sProb.Data());
    sJetNameX.Append(sCross.Data());
    hTrgBinP[iPtHatBin] -> SetName(sTrgNameP.Data());
    hTrgBinX[iPtHatBin] -> SetName(sTrgNameX.Data());
    hJetBinP[iPtHatBin] -> SetName(sJetNameP.Data());
    hJetBinX[iPtHatBin] -> SetName(sJetNameX.Data());
  }
  cout << "    Created histograms to sum.\n"
       << "    Grabbed raw event numbers:"
       << endl;
  
  UInt_t nRawEvt[NPtHatBins];
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    nRawEvt[iPtHatBin] = hNorm[iPtHatBin] -> GetEntries();
    cout << "      nRawEvt[" << iPtHatBin << "] = " << nRawEvt[iPtHatBin] < <endl;
  }
  cout << "    Calculated cross sections:" << endl;

  // calculate cross sections
  Double_t nTrgTot(0.);
  Double_t nEvtTot(0.);
  Double_t sigTot(0.);
  Double_t sigTrg[NPtHatBins];
  Double_t nCorrEvt[NPtHatBins];
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    nCorrEvt[iPtHatBin] = nRawEvt[iPtHatBin] / evntFrc[iPtHatBin];
    //sigTrg[iPtHatBin]   = SigPtHat[iPtHatBin] / nCorrEvt[iPtHatBin];
    sigTrg[iPtHatBin]   = SigPtHat[iPtHatBin] * evntFrc[iPtHatBin];
    nTrgTot            += nRawEvt[iPtHatBin];
    nEvtTot            += nCorrEvt[iPtHatBin];
    sigTot             += SigPtHat[iPtHatBin];
    cout << "      [" << iPtHatBin << "] corrected events = " << nCorrEvt[iPtHatBin] << ", x-section = " << SigPtHat[iPtHatBin]
         << ", trig. x-section = " << sigTrg[iPtHatBin]
         << endl;
  }
  cout << "    Calculated cross-sections:\n"
       << "      Total number of trig.s = " << nTrgTot << "\n"
       << "      Total number of events = " << nEvtTot << "\n"
       << "      Total cross section    = " << sigTot  << "\n"
       << "    Calculated normalization and pThat weights:"
       << endl;

  // calculate normalization and pThat weights
  Double_t trgNorm[NPtHatBins];
  Double_t ptHatWeight[NPtHatBins];
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    trgNorm[iPtHatBin]     = 1. / nRawEvt[iPtHatBin];
    ptHatWeight[iPtHatBin] = sigTrg[iPtHatBin] / sigTot;  // SM: this one is per event yield of triggers and can be compared directly to PYTHIA-8 and btw. 6-25 and 8-25
    cout << "      [" << iPtHatBin << "] normalization = " << trgNorm[iPtHatBin] << ", pThat weight = " << ptHatWeight[iPtHatBin] << endl;
  }

  // scale histograms
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    hTrgBinP[iPtHatBin] -> Scale(ptHatWeight[iPtHatBin] * trgNorm[iPtHatBin]);
    hJetBinP[iPtHatBin] -> Scale(ptHatWeight[iPtHatBin] * trgNorm[iPtHatBin]);
    hTrgBinX[iPtHatBin] -> Scale(sigTrg[iPtHatBin] * trgNorm[iPtHatBin]);
    hJetBinX[iPtHatBin] -> Scale(sigTrg[iPtHatBin] * trgNorm[iPtHatBin]);
  }
  cout << "    Scaled histograms." << endl;

  // sum histograms
  TH1D *hTrgSumP = (TH1D*) hTrgBinP[0] -> Clone();
  TH1D *hTrgSumX = (TH1D*) hTrgBinX[0] -> Clone();
  TH1D *hJetSumP = (TH1D*) hJetBinP[0] -> Clone();
  TH1D *hJetSumX = (TH1D*) hJetBinX[0] -> Clone();
  hTrgSumP -> SetName("hTrgSumProb");
  hTrgSumX -> SetName("hTrgSumCross");
  hJetSumP -> SetName("hJetSumProb");
  hJetSumX -> SetName("hJetSumCross");

  for (UInt_t iPtHatBin = 1; iPtHatBin < NPtHatBins; iPtHatBin++) {
    hTrgSumP -> Add(hTrgBinP[iPtHatBin]);
    hTrgSumX -> Add(hTrgBinX[iPtHatBin]);
    hJetSumP -> Add(hJetBinP[iPtHatBin]);
    hJetSumX -> Add(hJetBinX[iPtHatBin]);
  }
  cout << "    Summed trigger histograms." << endl;

  // normalize jet histograms 
  const Double_t jetNorm = nEvtTot / nTrgTot;
  for (UInt_t iPtHatBin = 1; iPtHatBin < NPtHatBins; iPtHatBin++) {
    hJetBinP[iPtHatBin] -> Scale(jetNorm);
    hJetBinX[iPtHatBin] -> Scale(jetNorm);
  }
  hJetSumP -> Scale(jetNorm);
  hJetSumX -> Scale(jetNorm);
  cout << "    Normalized jet histograms:\n"
       << "      jet norm = " << jetNorm
       << endl;

  // record raw no. of events, corrected no. of events, and stitch factors [Derek, 08.18.2022]
  TH1D *hEvtFracs     = new TH1D("hEvtFracs",     "Event fractions; #hat{p}_{T} bin; f_{evt}",               NPtHatBins, 0., (Double_t) NPtHatBins);
  TH1D *hTrgSigma     = new TH1D("hTrgSigma",     "Triggered x-section; #hat{p}_{T}; #sigma_{trg}",          NPtHatBins, 0., (Double_t) NPtHatBins);
  TH1D *hRawEvtNums   = new TH1D("hRawEvtNums",   "Raw no. of events; #hat{p}_{T}; N_{trg}",                 NPtHatBins, 0., (Double_t) NPtHatBins);
  TH1D *hCorrEvtNums  = new TH1D("hCorrEvtNums",  "Corrected no. of events; #hat{p}_{T}; N_{trg} / f_{evt}", NPtHatBins, 0., (Double_t) NPtHatBins);
  TH1D *hPtHatWeights = new TH1D("hPtHatWeights", "Applied #hat{p}_{T} weights; #hat{p}_{T}; f_{stitch}",    NPtHatBins, 0., (Double_t) NPtHatBins);
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    hEvtFracs     -> SetBinContent(iPtHatBin + 1, evntFrc[iPtHatBin]);
    hTrgSigma     -> SetBinContent(iPtHatBin + 1, sigTrg[iPtHatBin]);
    hRawEvtNums   -> SetBinContent(iPtHatBin + 1, nRawEvt[iPtHatBin]);
    hCorrEvtNums  -> SetBinContent(iPtHatBin + 1, nCorrEvt[iPtHatBin]);
    hPtHatWeights -> SetBinContent(iPtHatBin + 1, ptHatWeight[iPtHatBin]);
  }
  cout << "    Recorded numbers for calculations." << endl;

  // save histograms 
  fOut           -> cd();
  hTrgSumP       -> Write();
  hTrgSumX       -> Write();
  hJetSumP       -> Write();
  hJetSumX       -> Write();
  hEvtFracs      -> Write();
  hTrgSigma      -> Write();
  hRawEvtNums    -> Write();
  hCorrEvtNums   -> Write();
  hPtHatWeights -> Write();
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    hTrgBinP[iPtHatBin]  -> Write();
    hTrgBinX[iPtHatBin]  -> Write();
    hJetBinP[iPtHatBin]  -> Write();
    hJetBinX[iPtHatBin]  -> Write();
    hNorm[iPtHatBin]     -> Write();
    hTrgInput[iPtHatBin] -> Write();
    hJetInput[iPtHatBin] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iPtHatBin = 0; iPtHatBin < NPtHatBins; iPtHatBin++) {
    fIn[iPtHatBin] -> cd();
    fIn[iPtHatBin] -> Close();
  }
  cout << "  Sum finished!\n" << endl;

}

// End ------------------------------------------------------------------------
