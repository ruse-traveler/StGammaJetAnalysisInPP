// 'SmearTriggerEt.C'
// Derek Anderson
// 03.10.2021
//
// Takes a (measured) trigger
// eT spectrum and smears it
// according to a provided
// qTtrg function.

#include <vector>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom3.h"
#include "TDirectory.h"

using namespace std;

// global constants
static const UInt_t NBins(2);
static const Bool_t SmoothQt(true);



void SmearTriggerEt() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning trigger smearing..." << endl;

  // io parameters
  //const TString sOut("smearedEtMeas_usingFitForSmear.et911vz55tsp008pi0.d2m11y2021.root");
  const TString sOut("test.root");
  const TString sTrigIn("input/Weights/triggerWeights_et911pi0weights_full650pi0sample.et650x911vz55tsp008pi0.d24m3y2021.root");
  const TString sSmearIn("input/ParticleGun/test_withTSP_allowshowers.root");
  const TString sEtHistIn("hEtMeas");
  const TString sQtFuncIn("match/hQtMatch911");

  // calculation parameters
  const Bool_t   reuseEvents(false);
  const Double_t eTtrgCorr[2] = {9., 11.};

  // smoothing parameters
  const UInt_t  nPar();
  const TString sQtSmooth("fQtSmooth");
  const TString sFuncSmooth("gaus(0)+pol3(3)");

  // open files
  TFile *fOut     = new TFile(sOut.Data(), "recreate");
  TFile *fTrigIn  = new TFile(sTrigIn.Data(), "read");
  TFile *fSmearIn = new TFile(sSmearIn.Data(), "read");
  if (!fOut || !fTrigIn || !fSmearIn) {
    cerr <<"PANIC: couldn't open file!\n"
         << "      fOut     = " << fOut << "\n"
         << "      fTrigIn  = " << fTrigIn << "\n"
         << "      fSmearIn = " << fSmearIn << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab input
  TH1D *hEtTrg = (TH1D*) fTrigIn  -> Get(sEtHistIn.Data());
  TH1D *hQtTrg = (TH1D*)  fSmearIn -> Get(sQtFuncIn.Data());
  if (!hEtTrg || !hQtTrg) {
    cerr << "PANIC: couldn't grab weight histogram or qT function!\n"
         << "       hEtTrg = " << hEtTrg << "\n"
         << "       hQtTrg = " << hQtTrg
         << endl;
    return;
  }
  hEtTrg -> SetName("hEtTrgInput");
  hQtTrg -> SetName("hQtTrgInput");
  cout << "    Grabbed input histogram and qT function." << endl;

  // declare histograms
  TH1D *hNumTrgRaw;
  TH1D *hNumTrgSmear;
  TH1D *hNumTrgThrown;
  TH1D *hNumInRange;
  TH1D *hNumOutRange;
  TH1D *hEtTrgRaw;
  TH1D *hEtTrgSmear;
  TH1D *hEtTrgThrown;

  // binning
  const UInt_t nNum(3);
  const UInt_t eTnum(200);
  const Float_t nBins[NBins]  = {0., 3.};
  const Float_t eTbins[NBins] = {0., 100.};
  // nTrg histograms
  hNumTrgRaw    = new TH1D("hNumTrgRaw", "no. of raw triggers sampled from input", nNum, nBins[0], nBins[1]);
  hNumTrgSmear  = new TH1D("hNumTrgSmear", "no. of accepted smeared triggers", nNum, nBins[0], nBins[1]);
  hNumTrgThrown = new TH1D("hNumTrgThrown", "no. of thrown smeared triggers", nNum, nBins[0], nBins[1]);
  hNumInRange   = new TH1D("hNumInRange", "no. of thrown triggers in corr. range", nNum, nBins[0], nBins[1]);
  hNumOutRange  = new TH1D("hNumOutRange", "no. of thrown triggers outside corr. range", nNum, nBins[0], nBins[1]);
  // eTtrg histograms
  hEtTrgRaw     = new TH1D("hEtTrgRaw", "raw trigger E_{T} sampled from input", eTnum, eTbins[0], eTbins[1]);
  hEtTrgSmear   = new TH1D("hEtTrgSmear", "accepted smeared trigger E_{T}", eTnum, eTbins[0], eTbins[1]);
  hEtTrgThrown  = new TH1D("hEtTrgThrown", "thrown smeared trigger E_{T}", eTnum, eTbins[0], eTbins[1]);
  // errors
  hEtTrgRaw    -> Sumw2();
  hEtTrgSmear  -> Sumw2();
  hEtTrgThrown -> Sumw2();
  cout << "    Declared histograms." << endl;


  // set random seed
  TDatime *dateTime = new TDatime();
  Int_t    date     = dateTime -> GetDate();
  Int_t    time     = dateTime -> GetTime();
  Int_t    seed     = date * time;
  gRandom -> SetSeed(seed);
  cout << "    Set random seed:\n"
       << "      seed = " << seed
       << endl;

  // qT parameters
  const Double_t ampQt = hQtTrg -> GetMaximum();
  const Double_t muQt  = hQtTrg -> GetMean();
  const Double_t sigQt = hQtTrg -> GetRMS();
  cout << "    Grabbed qT parameters:\n"
       << "      ampl. = " << ampQt << "\n"
       << "      mu    = " << muQt << "\n"
       << "      sigma = " << sigQt
       << endl;

  // smooth qT histograms
  TF1 *fQtSmooth;
  if (SmoothQt) {
    fQtSmooth = new TF1(sQtSmooth.Data(), sFuncSmooth, eTbins[0], eTbins[1]);
    fQtSmooth -> SetParameter(0, ampQt);
    fQtSmooth -> SetParameter(1, muQt);
    fQtSmooth -> SetParameter(2, sigQt);
  }

  // grab no. of triggers (iterations)
  const UInt_t nTrg = hEtTrg -> GetEntries();
  cout << "    Grabbed number of triggers:\n"
       << "      nTrg = " << nTrg
       << endl;


  // back-smear eTtrg
  Double_t qTtrg(0.);
  Double_t eTraw(0.);
  Double_t eTsmear(0.);
  UInt_t   nTrgRaw(0);
  UInt_t   nTrgSmear(0);
  UInt_t   nTrgThrown(0);
  UInt_t   nInRange(0);
  UInt_t   nOutRange(0);
  Bool_t   isInRange(false);
  for (UInt_t iTrg = 0; iTrg < nTrg; iTrg++) {

    cout << "      Processing event " << iTrg + 1 << "/" << nTrg << "...\r" << flush;
    if (iTrg + 1 == nTrg) cout << endl;

    // sample eT
    eTraw = hEtTrg -> GetRandom();
    hEtTrgRaw -> Fill(eTraw);
    nTrgRaw++; 

    // sample qT
    isInRange = false;
    if (reuseEvents) {
      do {
        //qTtrg     = gRandom -> Gaus(muQt, sigQt);
        qTtrg     = hQtTrg -> GetRandom();
        eTsmear   = eTraw / qTtrg;
        isInRange = ((eTsmear > eTtrgCorr[0]) && (eTsmear < eTtrgCorr[1]));
        hEtTrgThrown -> Fill(eTsmear); 
        nTrgThrown++;
      } while (!isInRange);
    } else {
      qTtrg     = gRandom -> Gaus(muQt, sigQt);
      //qTtrg     = hQtTrg -> GetRandom();
      eTsmear   = eTraw / qTtrg;
      isInRange = ((eTsmear > eTtrgCorr[0]) && (eTsmear < eTtrgCorr[1]));
      hEtTrgThrown -> Fill(eTsmear);
      nTrgThrown++;
    }

    // fill histogram and increment counters
    hEtTrgSmear -> Fill(eTsmear);
    nTrgSmear++;
    if (isInRange) {
      nInRange++;
    } else {
      nOutRange++;
    }

  }  // end trigger loop
  cout << "    Finished trigger loop:\n"
       << "      nTrg(Raw, Thrown, Smear) = (" << nTrgRaw << ", " << nTrgThrown << ", " << nTrgSmear << ")\n"
       << "      nTrg(In, Out)Range       = (" << nInRange << ", " << nOutRange << ")"
       << endl;

  // record trigger numbers
  const UInt_t iBin = hNumTrgRaw -> FindBin(1.5);
  hNumTrgRaw    -> SetBinContent(iBin, nTrgRaw);
  hNumTrgSmear  -> SetBinContent(iBin, nTrgSmear);
  hNumTrgThrown -> SetBinContent(iBin, nTrgThrown);
  hNumInRange   -> SetBinContent(iBin, nInRange);
  hNumOutRange  -> SetBinContent(iBin, nOutRange);
  cout << "    Recorded no. of triggers." << endl;


  // save histograms
  fOut          -> cd();
  hNumTrgRaw    -> Write();
  hNumTrgSmear  -> Write();
  hNumTrgThrown -> Write();
  hNumInRange   -> Write();
  hNumOutRange  -> Write();
  hEtTrgRaw     -> Write();
  hEtTrgSmear   -> Write();
  hEtTrgThrown  -> Write();
  hEtTrg        -> Write();
  hQtTrg        -> Write();
  cout << "    Saved histograms." << endl;

  // save and close files
  fOut     -> cd();
  fOut     -> Close();
  fTrigIn  -> cd();
  fTrigIn  -> Close();
  fSmearIn -> cd();
  fSmearIn -> Close();
  cout << "  Jet reweighting finished!\n" << endl;

}

// End ------------------------------------------------------------------------
