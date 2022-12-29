// CalculatePythiaWeights.C
// Derek Anderson
// 03.23.2021
//
// Takes a Pythia eTtrg spectrum
// and a data eTtrg specturm.
// Then smears the data etTrg and
// calculates the corresponding
// weights to apply to Pythia

#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom3.h"

using namespace std;

// global constants
static const UInt_t NBins(2);
static const Bool_t SmoothWeights(false);



void CalculatePythiaWeights() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Calculating weights for Pythia..." << endl;

  // file parameters
  const TString sOut("weightCheck_qTsmearOverTruePy8_eTbinOne.et1520gam.d20m3y2022.root");
  const TString sDataIn("input/Weights/triggerWeights_et1520pi0weights_full650pi0sample.et650x1520vz55tsp008pi0.d24m3y2021.root");  // 1 for each eT bin
  const TString sSmearIn("input/ParticleGun/test_withTSP_allowshowers.root");
  const TString sPythiaIn("input/Pythia/pp200py8par.forPyQtSmearCheck_pTbinOne.et0100r02gam.d20m3y2022.root");  // stays the same for all eT bins

  // histogram parameters
  const TString sDataEtIn("hEtMeas");
  const TString sQtFuncIn("match/hQtMatch1520");  // 1 for eT bin
  const TString sPythiaEtIn("Gam/hTrgEtG");  // should be summed over pThat

  // calculation parameters
  const Bool_t   reuseEvents(false);
  const Double_t eTtrgCorr[2] = {15., 20.};

  // smoothing parameters
  const TString sQtSmooth("fQtSmooth");
  const TString sFuncSmooth("gaus(0)+pol3(3)");

  // open files
  TFile *fOut      = new TFile(sOut.Data(), "recreate");
  TFile *fDataIn   = new TFile(sDataIn.Data(), "read");
  TFile *fSmearIn  = new TFile(sSmearIn.Data(), "read");
  TFile *fPythiaIn = new TFile(sPythiaIn.Data(), "read");
  if (!fOut || !fDataIn || !fSmearIn || !fPythiaIn) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fDataIn   = " << fDataIn << ", fSmearIn = " << fSmearIn << "\n"
         << "       fPythiaIn = " << fPythiaIn << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hEtData   = (TH1D*) fDataIn   -> Get(sDataEtIn.Data());
  TH1D *hQtTrg    = (TH1D*) fSmearIn  -> Get(sQtFuncIn.Data());
  TH1D *hEtPythia = (TH1D*) fPythiaIn -> Get(sPythiaEtIn.Data());
  if (!hEtData || !hQtTrg || !hEtPythia) {
    cerr << "PANIC: couldn't grab an input histogram!\n"
         << "       hEtData = " << hEtData << ", hQtTrg = " << hQtTrg << ", hEtPythia = " << hEtPythia << "\n"
         << endl;
    return;
  }
  hEtData   -> SetName("hEtDataInput");
  hQtTrg    -> SetName("hQtTrgInput");
  hEtPythia -> SetName("hEtPythiaInput");
  cout << "    Grabbed input histograms." << endl;

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
  const UInt_t  nNum(3);
  const UInt_t  eTnum(200);
  const Float_t nBins[NBins]  = {0., 3.};
  const Float_t eTbins[NBins] = {0., 100.};
  // nTrg histograms
  hNumTrgRaw    = new TH1D("hNumTrgRaw",    "no. of raw triggers sampled from input",     nNum, nBins[0], nBins[1]);
  hNumTrgSmear  = new TH1D("hNumTrgSmear",  "no. of accepted smeared triggers",           nNum, nBins[0], nBins[1]);
  hNumTrgThrown = new TH1D("hNumTrgThrown", "no. of thrown smeared triggers",             nNum, nBins[0], nBins[1]);
  hNumInRange   = new TH1D("hNumInRange",   "no. of thrown triggers in corr. range",      nNum, nBins[0], nBins[1]);
  hNumOutRange  = new TH1D("hNumOutRange",  "no. of thrown triggers outside corr. range", nNum, nBins[0], nBins[1]);
  // eTtrg histograms
  hEtTrgRaw     = new TH1D("hEtTrgRaw",     "raw trigger E_{T} sampled from input",       eTnum, eTbins[0], eTbins[1]);
  hEtTrgSmear   = new TH1D("hEtTrgSmear",   "accepted smeared trigger E_{T}",             eTnum, eTbins[0], eTbins[1]);
  hEtTrgThrown  = new TH1D("hEtTrgThrown",  "thrown smeared trigger E_{T}",               eTnum, eTbins[0], eTbins[1]);
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

  // grab no. of triggers (iterations)
  const UInt_t nTrg = hEtData -> GetEntries();
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
    eTraw = hEtData -> GetRandom();
    hEtTrgRaw -> Fill(eTraw);
    nTrgRaw++; 

    // sample qT
    isInRange = false;
    if (reuseEvents) {
      do {
        qTtrg     = hQtTrg -> GetRandom();
        eTsmear   = eTraw / qTtrg;
        isInRange = ((eTsmear > eTtrgCorr[0]) && (eTsmear < eTtrgCorr[1]));
        hEtTrgThrown -> Fill(eTsmear); 
        nTrgThrown++;
      } while (!isInRange);
    } else {
      qTtrg     = hQtTrg -> GetRandom();
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
       << "      nTrg(Raw, Thrown, Smear) = (" << nTrgRaw  << ", " << nTrgThrown << ", " << nTrgSmear << ")\n"
       << "      nTrg(In, Out)Range       = (" << nInRange << ", " << nOutRange  << ")"
       << endl;

  // record trigger numbers
  const UInt_t iBin = hNumTrgRaw -> FindBin(1.5);
  hNumTrgRaw    -> SetBinContent(iBin, nTrgRaw);
  hNumTrgSmear  -> SetBinContent(iBin, nTrgSmear);
  hNumTrgThrown -> SetBinContent(iBin, nTrgThrown);
  hNumInRange   -> SetBinContent(iBin, nInRange);
  hNumOutRange  -> SetBinContent(iBin, nOutRange);
  cout << "    Recorded no. of triggers." << endl;

  // QUICK FIX [03.17.2022]
  //hEtTrgSmear -> Rebin(2);
  //hEtPythia   -> Rebin(2);

  // normalize histograms for weights
  Double_t iData   = hEtTrgSmear -> Integral();
  Double_t iPythia = hEtPythia   -> Integral();
  hEtTrgSmear -> Scale(1. / iData);
  hEtPythia   -> Scale(1. / iPythia);
  cout << "    Normalized input:\n"
       << "      iData = " << iData << ", iPythia = " << iPythia
       << endl;

  // calculate weights
  TH1D *hEtWeight = (TH1D*) hEtPythia -> Clone();
  hEtWeight -> SetNameTitle("hEtWeight", "weights [q_{T}^{trg} smeared / E_{T}^{true}]");
  hEtWeight -> Reset("ICES");

  // loop over data bins
  const UInt_t nBinsData   = hEtTrgSmear -> GetNbinsX();
  const UInt_t nBinsPythia = hEtPythia   -> GetNbinsX(); 
  for (UInt_t iBinD = 1; iBinD < (nBinsData + 1); iBinD++) {

    // locate pythia bin
    const Float_t binCntD = hEtTrgSmear -> GetBinCenter(iBinD);
    const UInt_t  iBinP   = hEtPythia -> FindBin(binCntD);
    const Bool_t  isInPy  = (iBinP < (nBinsPythia + 1));
    if (!isInPy) break;

    // calculate ratio
    const Double_t binValD  = hEtTrgSmear -> GetBinContent(iBinD);
    const Double_t binValP  = hEtPythia   -> GetBinContent(iBinP);
    const Double_t binErrD  = hEtTrgSmear -> GetBinError(iBinD);
    const Double_t binErrP  = hEtPythia   -> GetBinError(iBianP);
    const Double_t binRelD  = binErrD / binValD;
    const Double_t binRelP  = binErrP / binValP;
    const Double_t ratio    = binValD / binValP;
    const Double_t ratioErr = ratio * TMath::Sqrt((binRelD * binRelD) + (binRelP * binRelP));
    const Bool_t   datIsPos = (binValD > 0.);
    const Bool_t   pyIsPos  = (binValP > 0.);
    if (!datIsPos || !pyIsPos) {
      hEtWeight -> SetBinContent(iBinP, 0.);
      hEtWeight -> SetBinError(iBinP, 0.);
    } else {
      hEtWeight -> SetBinContent(iBinP, ratio);
      hEtWeight -> SetBinError(iBinP, ratioErr);
    }
  }  // end bin loop
  cout << "    Calculated weights." << endl;

  // save histograms
  fOut          -> cd();
  hEtData       -> Write();
  hQtTrg        -> Write();
  hEtPythia     -> Write();
  hNumTrgRaw    -> Write();
  hNumTrgSmear  -> Write();
  hNumTrgThrown -> Write();
  hNumInRange   -> Write();
  hNumOutRange  -> Write();
  hEtTrgRaw     -> Write();
  hEtTrgSmear   -> Write();
  hEtTrgThrown  -> Write();
  hEtWeight     -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut      -> Close();
  fDataIn   -> cd();
  fDataIn   -> Close();
  fSmearIn  -> cd();
  fSmearIn  -> Close();
  fPythiaIn -> cd();
  fPythiaIn -> Close();
  cout << "  Finished Pythia weight calculation!\n" << endl;

}

// End ------------------------------------------------------------------------
