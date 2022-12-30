// 'CheckTreeTriggers.C'
// Derek Anderson
// 03.01.2017
//
// Use this to check the triggers found
// by the 'StThirdMaker' class.

#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "TNtuple.h"

using namespace std;


// global constants
static const Double_t eStrMin  = 0.5;
static const Double_t pStrMin  = 0.5;
static const Double_t eTwrMin  = -0.9;
static const Double_t eTwrMax  = 0.9;
static const Double_t adcMax   = 6004;
static const Double_t pTwrMax  = 3.0;
static const Double_t hTrgMax  = 0.9;
static const Double_t eTtrgMin = 8.;
static const Double_t tspMinP  = 0.;
static const Double_t tspMaxP  = 0.08;
static const Double_t tspMinG  = 0.2;
static const Double_t tspMaxG  = 0.6;
static const Double_t pi = TMath::Pi();
// i/o parameters
static const TString sIn("merge.afterFix.r10145070.d12m4y2017.root");
static const TString sOut("tupleTrigs.r10145070.d12m4y2017.root");
static const TString sTupleEvt1("QAAB");
static const TString sTupleEvt2("AngleTAB");
static const TString sTupleEvt3("CAB");
static const TString sTupleEvt4("SMDE1B");
static const TString sTupleEvt5("PiB");


void CheckTupleTriggers() {

  cout << "\n  Beginning tuple-reading script..." << endl;
  gErrorIgnoreLevel = kError;


  TFile  *fOut = new TFile(sOut.Data(), "recreate");
  TFile  *fIn  = new TFile(sIn.Data(), "read");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    assert(fIn);
  }

  TNtuple *tEvt1 = (TNtuple*) fIn -> Get(sTupleEvt1.Data());
  TNtuple *tEvt2 = (TNtuple*) fIn -> Get(sTupleEvt2.Data());
  TNtuple *tEvt3 = (TNtuple*) fIn -> Get(sTupleEvt3.Data());
  TNtuple *tEvt4 = (TNtuple*) fIn -> Get(sTupleEvt4.Data());
  TNtuple *tEvt5 = (TNtuple*) fIn -> Get(sTupleEvt5.Data());
  if (!tEvt1 || !tEvt2 || !tEvt3 || !tEvt4 || !tEvt5)  {
    cerr << "PANIC: couldn't grab an input tuple!" << endl;
    assert(tEvt1 && tEvt2 && tEvt3 && tEvt4 && tEvt5);
  }

  cout << "    Tuples grabbed..." << endl;

  // QAAB leaves
  Float_t EneT0B;
  Float_t MPosB;
  Float_t MNegB;
  Float_t HB;
  Float_t en4B;
  Float_t enp4B;
  Float_t eneT0B;
  Float_t zB;
  Float_t etav1B;
  Float_t phiv1B;
  Float_t didTB;
  Float_t delta1B;
  Float_t dphi1B;
  Float_t Weta1B;
  // AngleTAB leaves
  Float_t moduleTB;
  Float_t EtaEB;
  Float_t sPB;
  Float_t runIdB;
  Float_t evtB;
  Float_t NOPB;
  Float_t Wphi1B;
  Float_t enp7B;
  Float_t enp1B;
  Float_t enp8B;
  Float_t enp0B;
  Float_t en01B;
  Float_t en02B;
  Float_t PTowerB;
  // CAB leaves
  Float_t didEB;
  Float_t ENET0B;
  Float_t XXB;
  Float_t YYB;
  Float_t adc11B;
  Float_t didPB;
  Float_t adcE1B;
  Float_t adcP1B;
  Float_t eTB;
  Float_t phTB;
  Float_t hSSB;
  Float_t hpedSSB;
  Float_t hRMSSSB;
  Float_t TC7B;
  // SMDE1B leaves
  Float_t en03B;
  Float_t en0B;
  Float_t en1B;
  Float_t en2B;
  Float_t en3B;
  Float_t en5B;
  Float_t en6B;
  Float_t en7B;
  Float_t en8B;
  Float_t enp2B;
  Float_t enp3B;
  Float_t enp5B;
  Float_t enp6B;
  // PiB leaves
  Float_t en9B;
  Float_t en10B;
  Float_t en11B;
  Float_t enp03B;
  Float_t enp02B;
  Float_t enp01B;
  Float_t enp9B;
  Float_t enp10B;
  Float_t enp11B;
  

  // set QAAB addresses
  tEvt1 -> SetBranchAddress("EneT0", &EneT0B);
  tEvt1 -> SetBranchAddress("MPos", &MPosB);
  tEvt1 -> SetBranchAddress("MNeg", &MNegB);
  tEvt1 -> SetBranchAddress("H", &HB);
  tEvt1 -> SetBranchAddress("en4", &en4B);
  tEvt1 -> SetBranchAddress("enp4", &enp4B);
  tEvt1 -> SetBranchAddress("eneT0", &eneT0B);
  tEvt1 -> SetBranchAddress("z", &zB);
  tEvt1 -> SetBranchAddress("etav1", &etav1B);
  tEvt1 -> SetBranchAddress("phiv1", &phiv1B);
  tEvt1 -> SetBranchAddress("didT", &didTB);
  tEvt1 -> SetBranchAddress("delta1", &delta1B);
  tEvt1 -> SetBranchAddress("dphi1", &dphi1B);
  tEvt1 -> SetBranchAddress("Weta1", &Weta1B);
  // set AngleTAB addresses
  tEvt2 -> SetBranchAddress("moduleT", &moduleTB);
  tEvt2 -> SetBranchAddress("EtaE", &EtaEB);
  tEvt2 -> SetBranchAddress("sP", &sPB);
  tEvt2 -> SetBranchAddress("runId", &runIdB);
  tEvt2 -> SetBranchAddress("evt", &evtB);
  tEvt2 -> SetBranchAddress("NOP", &NOPB);
  tEvt2 -> SetBranchAddress("Wphi1", &Wphi1B);
  tEvt2 -> SetBranchAddress("enp7", &enp7B);
  tEvt2 -> SetBranchAddress("enp1", &enp1B);
  tEvt2 -> SetBranchAddress("enp8", &enp8B);
  tEvt2 -> SetBranchAddress("enp0", &enp0B);
  tEvt2 -> SetBranchAddress("en01", &en01B);
  tEvt2 -> SetBranchAddress("en02", &en02B);
  tEvt2 -> SetBranchAddress("PTower", &PTowerB);
  // set CAB addresses
  tEvt3 -> SetBranchAddress("didE", &didEB);
  tEvt3 -> SetBranchAddress("ENET0", &ENET0B);
  tEvt3 -> SetBranchAddress("XX", &XXB);
  tEvt3 -> SetBranchAddress("YY", &YYB);
  tEvt3 -> SetBranchAddress("adc11", &adc11B);
  tEvt3 -> SetBranchAddress("didP", &didPB);
  tEvt3 -> SetBranchAddress("adcE1", &adcE1B);
  tEvt3 -> SetBranchAddress("adcP1", &adcP1B);
  tEvt3 -> SetBranchAddress("eT", &eTB);
  tEvt3 -> SetBranchAddress("phT", &phTB);
  tEvt3 -> SetBranchAddress("hSS", &hSSB);
  tEvt3 -> SetBranchAddress("hpedSS", &hpedSSB);
  tEvt3 -> SetBranchAddress("hRMSSS", &hRMSSSB);
  tEvt3 -> SetBranchAddress("TC7", &TC7B);
  // SMDE1B leaves
  tEvt4 -> SetBranchAddress("en03", &en03B);
  tEvt4 -> SetBranchAddress("en0", &en0B);
  tEvt4 -> SetBranchAddress("en1", &en1B);
  tEvt4 -> SetBranchAddress("en2", &en2B);
  tEvt4 -> SetBranchAddress("en3", &en3B);
  tEvt4 -> SetBranchAddress("en5", &en5B);
  tEvt4 -> SetBranchAddress("en6", &en6B);
  tEvt4 -> SetBranchAddress("en7", &en7B);
  tEvt4 -> SetBranchAddress("en8", &en8B);
  tEvt4 -> SetBranchAddress("enp2", &enp2B);
  tEvt4 -> SetBranchAddress("enp3", &enp3B);
  tEvt4 -> SetBranchAddress("enp5", &enp5B);
  tEvt4 -> SetBranchAddress("enp6", &enp6B);
  // PiB leaves
  tEvt5 -> SetBranchAddress("en9", &en9B);
  tEvt5 -> SetBranchAddress("en10", &en10B);
  tEvt5 -> SetBranchAddress("en11", &en11B);
  tEvt5 -> SetBranchAddress("enp03", &enp03B);
  tEvt5 -> SetBranchAddress("enp02", &enp02B);
  tEvt5 -> SetBranchAddress("enp01", &enp01B);
  tEvt5 -> SetBranchAddress("enp9", &enp9B);
  tEvt5 -> SetBranchAddress("enp10", &enp10B);
  tEvt5 -> SetBranchAddress("enp11", &enp11B);

  cout << "    Branches set..." << endl;


  // initialize histograms
  Int_t    nT = 10;
  Int_t    nR = 15;
  Double_t t1 = 0.;
  Double_t t2 = 10.;
  Double_t r1 = 0.;
  Double_t r2 = 1.5;
  TH1D *hNumTrg   = new TH1D("hNumTrg", "No. of triggers per event", nT, t1, t2);
  TH1D *hNumPass  = new TH1D("hNumPass", "No. of accepted triggers per event", nT, t1, t2);
  TH1D *hNumRatio = new TH1D("hNumRatio", "N_{pass} / N_{trg}", nR, r1, r2);

  cout << "    Histograms initialized..." << endl;


  Int_t nTrgTot  = 0;
  Int_t nTrgTot1 = tEvt1 -> GetEntries();
  Int_t nTrgTot2 = tEvt2 -> GetEntries();
  if (nTrgTot1 != nTrgTot2) {
    cerr << "PANIC: 'tEvt1' and 'tEvt2' do not have\n"
         << "       the same number of entries!"
         << endl;
    assert(nTrgTot1 == nTrgTot2);
  }
  else
    nTrgTot = nTrgTot1;

  cout << "    Beginning event loop..." << endl;

  // event loop
  Int_t    nEvent    = -1;
  Int_t    nTriggers = 1;
  Int_t    nTrigPass = 0;
  Int_t    nTrigDiff = 0;
  Bool_t   isSameEvt = false;
  Double_t ratio     = 0;
  for (Int_t i = 0; i < nTrgTot; i++) {

    tEvt1 -> GetEntry(i);
    tEvt2 -> GetEntry(i);
    tEvt3 -> GetEntry(i);
    tEvt4 -> GetEntry(i);
    tEvt5 -> GetEntry(i);
    cout << "      Processing trigger " << (i + 1) << "/" << nTrgTot << "...\r" << flush;
    if ((i + 1) == nTrgTot) cout << endl;

    // check if same event as previous iteration
    if (nEvent != evtB) {
      hNumTrg   -> Fill(nTriggers);
      if (nTrigPass != 0) {
        ratio = (Double_t) nTrigPass / (Double_t) nTriggers;
        hNumPass  -> Fill(nTrigPass);
        hNumRatio -> Fill(ratio);
      }
      nEvent    = evtB;
      nTriggers = 1;
      nTrigPass = 0;
      isSameEvt = false;
    }
    else {
      ++nTriggers;
      isSameEvt = true;
    }

    // event cuts
    Bool_t isInStrCut  = ((en4B >= eStrMin) && (enp4B >= pStrMin));
    Bool_t isInTwrCut  = ((eTB > eTwrMin) || (eTB < eTwrMax));
    Bool_t isInAdcCut  = (adc11B <= adcMax);
    Bool_t isInPTwrCut = (PTowerB < pTwrMax);
    if (!isInStrCut)
      continue;
    if (!isInTwrCut)
      continue;
    if (!isInAdcCut)
      continue;
    if (!isInPTwrCut)
      continue;

    Double_t hTrg  = etav1B;
    Double_t tTrg  = 2. * atan(exp(-1. * hTrg));
    Double_t eTrg  = EneT0B;
    Double_t eTtrg = eTrg * sin(tTrg);

    // calculate TSP
    Double_t tspN = eneT0B;
    Double_t tspD = (en4B + enp4B) * 0.783 + (en3B + enp3B + en5B + enp5B) * 2.21 + (en2B + enp2B + en6B + enp6B) * 6.26
                  + (en1B + enp1B + en7B + enp7B) * 11.51 + (en0B + enp0B + en8B + enp8B) * 17.73 + (en01B + enp01B + en9B + enp9B) * 24.78
                  + (en02B + enp02B + en10B + enp10B) * 32.57 + (en03B + enp03B + en11B + enp11B) * 41.05;
    Double_t tsp  = tspN / tspD;

 
    // trigger cuts
    Bool_t isInEtaCut  = (abs(hTrg) < hTrgMax);
    Bool_t isInEneCut  = (eTtrg > eTtrgMin);
    Bool_t isInTspPCut = ((tsp > tspMinP) && (tsp < tspMaxP));
    Bool_t isInTspGCut = ((tsp > tspMinG) && (tsp < tspMaxG));
    if (!isInEtaCut)
      continue;
    if (!isInEneCut)
      continue;
    if (!isInTspPCut && !isInTspGCut)
      continue;
    
    if (!isSameEvt)
      nTrigPass = 1;
    else
      ++nTrigPass;

  }  // end event loop

  cout << "    Event loop finished!" << endl;


  // close and save
  fOut      -> cd();
  hNumTrg   -> Write();
  hNumPass  -> Write();
  hNumRatio -> Write();
  fOut      -> Close();
  fIn       -> cd();
  fIn       -> Close();

  cout << "  Tuple reading script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
