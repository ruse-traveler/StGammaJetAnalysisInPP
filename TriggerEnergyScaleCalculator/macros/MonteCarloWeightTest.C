// 'MonteCarloWeightTest.C'
// Derek Anderson
// 01.28.2021
//
// Toy MC to test eTtrg weighting
// proceedure.


#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom3.h"

using namespace std;



void MonteCarloWeightTest() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning Monte Carlo to test weights..." << endl;

  // io parameters
  const TString sOut("test.root");
  //const TString sOut("testMonteCarlo_reusingEventsAndAdjustingWeights_full620pi0sample.et620x1115vz55tsp008pi0.d2m2y2021.root");
  const TString sInEtReco("testingWeights_full620pi0sample.et620x1115vz55tsp008pi0.d28m1y2021.root");
  const TString sInWeight("testingWeights_full620pi0sample.et620x1115vz55tsp008pi0.d28m1y2021.root");
  const TString sInQtTrg("testingWeights_full620pi0sample.et620x1115vz55tsp008pi0.d28m1y2021.root");
  const TString sHistEtReco("hEtMeas");
  const TString sHistWeight("hEtWeights");
  const TString sFuncQtTrg("fQtFit");

  // hist parameters
  const UInt_t  nBool(4);
  const UInt_t  nQtBins(100);
  const UInt_t  nEtBins(200);
  const Float_t boBin[2]  = {-1., 3.};
  const Float_t qTbin[2]  = {0., 10.};
  const Float_t eTbin[2]  = {0., 100.};
  const Float_t eTcalc[2] = {11., 15.};
  const TString sIsInRange("hEtTrueIsInRange");
  const TString sQtTrgSample("hQtTrgSampled");
  const TString sEtRecoSample("hEtRecoSampled");
  const TString sEtTrueSample("hEtTrueSampled");
  const TString sEtTrueWeight("hEtTrueWeighted");

  // mc parameters
  const UInt_t  nEvts(1000000);
  const Bool_t  reuseEvents(true);
  const Bool_t  restrictQt(false);
  const Bool_t  adjustWeightAvg(true);
  const Float_t adjustAvgTo(1.); 


  // set random seed
  TDatime *dateTime = new TDatime();
  Int_t    date     = dateTime -> GetDate();
  Int_t    time     = dateTime -> GetTime();
  Int_t    seed     = date * time;
  gRandom -> SetSeed(seed);
  cout << "    Set random seed:\n"
       << "      seed = " << seed
       << endl;


  // open files
  TFile *fOut    = new TFile(sOut.Data(), "recreate");
  TFile *fInEtReco = new TFile(sInEtReco.Data(), "read");
  TFile *fInWeight = new TFile(sInWeight.Data(), "read");
  TFile *fInQtTrg  = new TFile(sInQtTrg.Data(), "read");
  if (!fOut || !fInEtReco || !fInWeight || !fInQtTrg) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInEtReco = " << fInEtReco << "\n"
         << "       fInWeight = " << fInWeight << ", fInQtTrg = " << fInQtTrg
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input
  TH1D *hEtReco = (TH1D*) fInEtReco -> Get(sHistEtReco.Data());
  TH1D *hWeight = (TH1D*) fInWeight -> Get(sHistWeight.Data());
  TF1  *fQtTrg  = (TF1*)  fInQtTrg  -> Get(sFuncQtTrg.Data());
  if (!hEtReco || !hWeight || !fQtTrg) {
    cerr << "PANIC: couldn't grab an input!\n"
         << "       hEtReco = " << hEtReco << ", hWeight = " << hWeight << ", fQtTrg = " << fQtTrg
         << endl;
    return;
  }
  hEtReco -> SetName("hEtRecoInput");
  hWeight -> SetName("hEtWeightInput");
  fQtTrg  -> SetName("fQtTrgInput");
  cout << "    Grabbed input." << endl;

  // qT parameters
  const Double_t ampQt = fQtTrg -> GetParameter(0);
  const Double_t muQt  = fQtTrg -> GetParameter(1);
  const Double_t sigQt = fQtTrg -> GetParameter(2);
  cout << "    Grabbed qT parameters:\n"
       << "      ampl. = " << ampQt << "\n"
       << "      mu    = " << muQt << "\n"
       << "      sigma = " << sigQt
       << endl;


  // for if restrciting qT throws
  TF1 *fQtRand;
  if (restrictQt) {
    fQtRand = new TF1("fQtRand", "gaus(0)", qTbin[0], qTbin[1]);
    fQtRand -> FixParameter(0, ampQt);
    fQtRand -> FixParameter(1, muQt);
    fQtRand -> FixParameter(2, sigQt);
    cout << "    Declared restricted qT function." << endl;
  }

  // for adjusting input weights
  TH1D *hUseWeight = hWeight -> Clone();
  hUseWeight -> SetName("hEtWeightUsed");

  if (adjustWeightAvg) {
    const UInt_t  iFirstNonzero = hWeight -> FindFirstBinAbove(0.);
    const UInt_t  iLastNonzero  = hWeight -> FindLastBinAbove(0.);
    const Float_t xFirstAvg     = hWeight -> GetBinLowEdge(iFirstNonzero);
    const Float_t xLastAvg      = hWeight -> GetBinLowEdge(iLastNonzero + 1);

    TF1 *fInitAvg = new TF1("fInitAvg", "pol0(0)", xFirstAvg, xLastAvg);
    TF1 *fAdjAvg  = new TF1("fAdjAvg", "pol0(0)", xFirstAvg, xLastAvg);
    fInitAvg   -> SetParameter(0, adjustAvgTo);
    fAdjAvg    -> SetParameter(0, adjustAvgTo);
    fInitAvg   -> SetLineColor(880);
    fAdjAvg    -> SetLineColor(600);
    hUseWeight -> Fit("fInitAvg", "RWW+");

    const Float_t initAvg  = fInitAvg -> GetParameter(0);
    const Float_t adjust   = adjustAvgTo - initAvg;
    for (UInt_t iAdjBin = iFirstNonzero; iAdjBin <= iLastNonzero; iAdjBin++) {
      const Float_t binVal = hUseWeight -> GetBinContent(iAdjBin);
      const Float_t binErr = hUseWeight -> GetBinError(iAdjBin);
      const Float_t perErr = binErr / binVal;
      const Float_t newVal = binVal + adjust;
      const Float_t newErr = perErr * newVal;
      hUseWeight -> SetBinContent(iAdjBin, newVal);
      hUseWeight -> SetBinError(iAdjBin, newErr);
    }
    hUseWeight -> Fit("fAdjAvg", "RWW+");

    const Float_t adjAvg = fAdjAvg -> GetParameter(0);
    cout << "    Adjusted input weight average:\n"
         << "      range adjust  = (" << xFirstAvg << ", " << xLastAvg << ")\n"
         << "      initial avg.  = " << initAvg << "\n"
         << "      adjusted avg. = " << adjAvg
         << endl;
  }


  // set output histograms
  TH1D *hIsInRange   = new TH1D(sIsInRange.Data(), "Is thrown E_{T}^{true} in range?", nBool, boBin[0], boBin[1]);
  TH1D *hQtTrgSamp   = new TH1D(sQtTrgSample.Data(), "Sampled q_{T}^{trg}", nQtBins, qTbin[0], qTbin[1]);
  TH1D *hEtRecSamp   = new TH1D(sEtRecoSample.Data(), "Sampled E_{T}^{reco}", nEtBins, eTbin[0], eTbin[1]);
  TH1D *hEtTruSamp   = new TH1D(sEtTrueSample.Data(), "Sampled E_{T}^{true}", nEtBins, eTbin[0], eTbin[1]);
  TH1D *hEtTruWeight = new TH1D(sEtTrueWeight.Data(), "Weighted E_{T}^{true}", nEtBins, eTbin[0], eTbin[1]);
  hQtTrgSamp   -> Sumw2();
  hEtRecSamp   -> Sumw2();
  hEtTruSamp   -> Sumw2();
  hEtTruWeight -> Sumw2();
  cout << "    Set output histograms.\n"
       << "    Beginning MC loop:"
       << endl;


  // mc loop
  UInt_t nEtInRange(0);
  UInt_t nEtOutRange(0);
  for (UInt_t iEvt = 0; iEvt < nEvts; iEvt++) {

    // progress bar
    cout << "      Processing event " << (iEvt + 1) << "/" << nEvts << "...\r" << flush;
    if ((iEvt + 1) == nEvts) cout << endl;

    // sample eT
    const Double_t eTrec = hEtReco -> GetRandom();
    const Double_t qTmin = TMath::Min(eTrec / eTcalc[0], eTrec / eTcalc[1]);
    const Double_t qTmax = TMath::Max(eTrec / eTcalc[0], eTrec / eTcalc[1]);

    // sample qT and calculate eTtru
    Double_t qTtrg(0.);
    Double_t eTtru(0.);
    Bool_t   isInRange(false);

    // sample qT until eTtru is in range
    if (reuseEvents) {
      do {
        qTtrg     = gRandom -> Gaus(muQt, sigQt);
        eTtru     = eTrec / qTtrg;
        isInRange = ((eTtru > eTcalc[0]) && (eTtru < eTcalc[1]));
      }  while (!isInRange);
    }
    // restrict qTtrg to a particular range
    if (restrictQt) {
      qTtrg     = fQtRand -> GetRandom(qTmin, qTmax);
      eTtru     = eTrec / qTtrg;
      isInRange = ((eTtru > eTcalc[0]) && (eTtru < eTcalc[1]));
    }
    // default calculation
    if (!reuseEvents && !restrictQt) {
      qTtrg     = gRandom -> Gaus(muQt, sigQt);
      eTtru     = eTrec / qTtrg;
      isInRange = ((eTtru > eTcalc[0]) && (eTtru < eTcalc[1]));
    }

    // fill sampled histograms
    hEtRecSamp -> Fill(eTrec);
    hQtTrgSamp -> Fill(qTtrg);
    hEtTruSamp -> Fill(eTtru);

    // grab weights
    const UInt_t   iEtWeight = hUseWeight -> FindBin(eTtru);
    const Double_t eTweight  = hUseWeight -> GetBinContent(iEtWeight);
    if (isInRange) {
      hEtTruWeight -> Fill(eTtru, eTweight);
      hIsInRange   -> Fill(1.5);
      nEtInRange++;
    }
    else {
      hEtTruWeight -> Fill(eTtru, 0.);
      hIsInRange   -> Fill(0.5);
      nEtOutRange++;
    }

  }  // end mc loop
  cout << "    MC loop finished:\n"
       << "      no. of eTtrue values in range     = " << nEtInRange << "\n"
       << "      no. of eTtrue values out of range = " << nEtOutRange
       << endl;


  // close files
  fOut         -> cd();
  hEtReco      -> Write();
  hWeight      -> Write();
  hUseWeight   -> Write();
  fQtTrg       -> Write();
  hQtTrgSamp   -> Write();
  hEtRecSamp   -> Write();
  hEtTruSamp   -> Write();
  hEtTruWeight -> Write();
  hIsInRange   -> Write();
  fOut         -> Close();
  fInEtReco    -> cd();
  fInEtReco    -> Close();
  fInWeight    -> cd();
  fInWeight    -> Close();
  fInQtTrg     -> cd();
  fInQtTrg     -> Close();
  cout << "  Test Monte Carlo finished!\n" << endl;

}

// End ------------------------------------------------------------------------
