// 'CorrectTriggerSpectrum.C'
// Derek Anderson
// 01.19.2021
//
// Takes a measured (or reconstructed)
// trigger eT spectrum and "smears" it
// according to a provided trigger
// energy scale/resolution.  It then
// computes a set of weights to be
// applied to recoil jet spectra


#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom3.h"

using namespace std;



void CorrectTriggerSpectrum() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning trigger spectrum correction..." << endl;

  // io parameters
  const TString sOut("qTsmearPy8.et650x911vz55tsp008pi0.d20m3y2022.root");
  const TString sInEtMeas("input/Pythia/pp200py8par.forPyQtSmearCheck_pTbinOne.et911r02pi0.d20m3y2022.root");
  const TString sInEtPar("output/January2022/particleGun.full650pi0sample_withEtMatrix.et650x650x920vz55tsp008pi0.d17m1y2022.root");
  const TString sInQtTrg("output/January2022/particleGun.full650pi0sample_withEtMatrix.et650x650x920vz55tsp008pi0.d17m1y2022.root");

  // hist parameters
  const TString sHistEtMeas("Pi0/hTrgEtP");
  const TString sHistEtPar("match/hEtMatch911");
  const TString sHistQtTrg("match/hQtMatch911");

  // calculation parameters
  const Float_t ampGuess(1.);
  const Float_t muGuess(0.98);
  const Float_t sigGuess(0.08);
  const Float_t avgGuess(1.);
  const Float_t eTcalc[2]   = {9., 11.};
  const Float_t fitRange[2] = {0.7, 1.2};


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
  TFile *fEtMeas = new TFile(sInEtMeas.Data(), "read");
  TFile *fEtPar  = new TFile(sInEtPar.Data(), "read");
  TFile *fQtTrg  = new TFile(sInQtTrg.Data(), "read");
  if (!fOut || !fEtMeas || !fEtPar || !fQtTrg) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fEtMeas = " << fEtMeas << "\n"
         << "       fEtPar = " << fEtPar << ", fQtTrg = " << fQtTrg
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hEtMeas = (TH1D*) fEtMeas -> Get(sHistEtMeas.Data());
  TH1D *hEtPar  = (TH1D*) fEtPar  -> Get(sHistEtPar.Data());
  TH1D *hQtTrg  = (TH1D*) fQtTrg  -> Get(sHistQtTrg.Data());
  if (!hEtMeas || !hEtPar || !hQtTrg) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hEtMeas = " << hEtMeas << ", hEtPar = " << hEtPar << ", hQtTrg = " << hQtTrg
         << endl;
    return;
  }
  hEtMeas -> SetName("hEtMeas");
  hEtPar  -> SetName("hEtPar");
  hQtTrg  -> SetName("hQtTrg");
  cout << "    Grabbed histograms." << endl;


  // normalize relevant histograms
  const Float_t iEtParNorm = hEtPar -> Integral();
  const Float_t iQtTrgNorm = hQtTrg -> Integral();
  hEtMeas -> Scale(1. / iEtParNorm);
  hEtPar  -> Scale(1. / iEtParNorm);
  hQtTrg  -> Scale(1. / iQtTrgNorm);


  // fit qT histogram
  TF1 *fQtFit = new TF1("fQtFit", "gaus(0)", fitRange[0], fitRange[1]);
  fQtFit -> SetParName(0, "amp");
  fQtFit -> SetParName(1, "#mu");
  fQtFit -> SetParName(2, "#sigma");
  fQtFit -> SetParameter(0, ampGuess);
  fQtFit -> SetParameter(1, muGuess);
  fQtFit -> SetParameter(2, sigGuess);
  fQtFit -> SetLineColor(2);
  hQtTrg -> Fit("fQtFit", "R");

  const Double_t muFit  = fQtFit -> GetParameter(1);
  const Double_t sigFit = fQtFit -> GetParameter(2);
  cout << "    Fit qT histogram:\n"
       << "      range = (" << fitRange[0] << ", " << fitRange[1] << ")\n"
       << "      mean  = " << muFit << "\n"
       << "      sigma = " << sigFit
       << endl;


  // declare histograms for smearing calculation
  TH1D *hEtSmear  = (TH1D*) hEtPar -> Clone();
  TH1D *hQtThrown = (TH1D*) hQtTrg -> Clone();
  hEtSmear  -> Reset("ICES");
  hQtThrown -> Reset("ICES");
  hEtSmear  -> SetNameTitle("hEtSmear_noNorm", "Smeared E_{T}^{meas} (unnormalized)");
  hQtThrown -> SetNameTitle("hQtThrown", "Thrown values of q_{T}^{trg}");

  // smear measured eT
  const UInt_t nEtMeas  = hEtMeas -> GetEntries();
  for (UInt_t iEtSmear = 0; iEtSmear < nEtMeas; iEtSmear++) {

    // smear eT
    const Double_t eTmeas  = hEtMeas -> GetRandom();
    //const Double_t qTrand  = gRandom -> Gaus(muFit, sigFit);
    const Double_t qTrand  = hQtTrg -> GetRandom();
    const Double_t eTsmear = eTmeas / qTrand;

    // fill histograms
    hEtSmear  -> Fill(eTsmear);
    hQtThrown -> Fill(qTrand);

  }  // end eT sampling loop

  // normalize thrown qT histogram
  const Float_t iQtThrown = hQtThrown -> Integral();
  hQtThrown -> Scale(1. / iQtThrown);

  // compare smeared and measured eT
  const Float_t muMeas   = hEtMeas  -> GetMean();
  const Float_t muSmear  = hEtSmear -> GetMean();
  const Float_t rmsMeas  = hEtMeas  -> GetRMS();
  const Float_t rmsSmear = hEtSmear -> GetRMS();
  cout << "    Smeared measured eT:\n"
       << "      mean(measured, smeared) = (" << muMeas << ", " << muSmear << ")\n"
       << "      rms(measured, smeared)  = (" << rmsMeas << ", " << rmsSmear << ")"
       << endl;

  // determine smeared range
  const UInt_t  iMinBin    = hEtSmear -> FindFirstBinAbove(0.);
  const UInt_t  iMaxBin    = hEtSmear -> FindLastBinAbove(0.);
  const Float_t eTsmearMin = hEtSmear -> GetBinLowEdge(iMinBin);
  const Float_t eTsmearMax = hEtSmear -> GetBinLowEdge(iMaxBin + 1);
  cout << "    Smeared eT range:\n"
       << "      eTsmear = (" << eTsmearMin << ", " << eTsmearMax << ")"
       << endl;


  // normalize smeared distribution
  TH1D *hEtSmearN = (TH1D*) hEtSmear -> Clone();
  hEtSmearN -> SetNameTitle("hEtSmear_normedToPar", "Smeared E_{T}^{meas} (normalized to E_{T}^{sim})");

  // determine normalization range
  const Float_t binWidth    = hEtPar -> GetBinWidth(1);
  const Float_t eTnormStart = eTcalc[0] + (binWidth / 2.);
  const Float_t eTnormStop  = eTcalc[1] - (binWidth / 2.);
  const UInt_t  iNormLo     = hEtPar -> FindBin(eTnormStart);
  const UInt_t  iNormHi     = hEtPar -> FindBin(eTnormStop);
  cout << "    Normalization range:\n"
       << "      eTnorm(Start, Stop)   = (" << eTnormStart << ", " << eTnormStop << ")\n"
       << "      iBinNorm(Start, Stop) = (" << iNormLo << ", " << iNormHi << ")"
       << endl;

  // calculate normalization
  const Float_t iEtSmearN = hEtSmearN -> Integral(iNormLo, iNormHi);
  const Float_t iEtParN   = hEtPar    -> Integral(iNormLo, iNormHi);
  const Float_t smearNorm = iEtParN / iEtSmearN;
  hEtSmearN -> Scale(smearNorm);
  cout << "    Normalized smeared histograms:\n"
       << "      norm wrt.simulated eT  = " << smearNorm
       << endl; 


  // compute weights
  TH1D *hEtWeights = (TH1D*) hEtSmear -> Clone();
  TH1D *hEtParCalc = (TH1D*) hEtPar   -> Clone();
  hEtWeights -> SetNameTitle("hEtWeights", "E_{T}^{trg} weights (normalized to E_{T}^{sim})");
  hEtParCalc -> SetNameTitle("hEtParCalc", "E_{T}^{sim} for weight calculation");
  hEtWeights -> Reset("ICES");

  const UInt_t nWeightBins = hEtWeights -> GetNbinsX();
  for (UInt_t iWeightBin = 1; iWeightBin < (nWeightBins + 1); iWeightBin++) {

    // grab bin values
    const Float_t valPar       = hEtPar    -> GetBinContent(iWeightBin);
    const Float_t errPar       = hEtPar    -> GetBinError(iWeightBin);
    const Float_t valSmear     = hEtSmearN -> GetBinContent(iWeightBin);
    const Float_t errSmear     = hEtSmearN -> GetBinError(iWeightBin);
    const Float_t binCenter    = hEtPar    -> GetBinCenter(iWeightBin);
    const Bool_t  binIsInRange = ((binCenter > eTcalc[0]) && (binCenter < eTcalc[1]));

    // calculate weights
    const Float_t valWeight = valPar / valSmear;
    const Float_t errWeight = TMath::Sqrt(((errPar / valPar) * (errPar / valPar)) + ((errSmear / valSmear) * (errSmear / valSmear))) * valWeight;

    // set histograms
    if (binIsInRange && (valSmear > 0.)) {
      hEtWeights -> SetBinContent(iWeightBin, valWeight);
      hEtParCalc -> SetBinContent(iWeightBin, valPar);
      hEtWeights -> SetBinError(iWeightBin, errWeight);
      hEtParCalc -> SetBinError(iWeightBin, errPar);
    }
    else {
      hEtWeights -> SetBinContent(iWeightBin, 0.);
      hEtParCalc -> SetBinContent(iWeightBin, 0.);
      hEtWeights -> SetBinError(iWeightBin, 0.);
      hEtParCalc -> SetBinError(iWeightBin, 0.);
    }

  }  // end bin loop
  cout << "    Computed weights." << endl;

  // fit weights
  TF1 *fAvgWeights = new TF1("fAvgWeights", "pol0(0)", eTsmearMin, eTsmearMax);
  fAvgWeights -> SetParameter(0, avgGuess);
  fAvgWeights -> SetLineColor(2);
  hEtWeights  -> Fit("fAvgWeights", "RWW");


  // close files and save
  fOut        -> cd();
  hEtMeas     -> Write();
  hEtPar      -> Write();
  hQtTrg      -> Write();
  hQtThrown   -> Write();
  hEtSmear    -> Write();
  hEtSmearN   -> Write();
  hEtWeights  -> Write();
  hEtParCalc  -> Write();
  fQtFit      -> Write();
  fAvgWeights -> Write(); 
  fOut        -> Close();
  fEtMeas     -> cd();
  fEtMeas     -> Close();
  fEtPar      -> cd();
  fEtPar      -> Close();
  fQtTrg      -> cd();
  fQtTrg      -> Close();
  cout << "  Trigger spectrum correction finished!\n" << endl;

}

// End ------------------------------------------------------------------------
