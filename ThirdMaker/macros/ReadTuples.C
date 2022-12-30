// 'ReadTuples.C'
// Derek Anderson
// 02.27.2017
//
// Use this to read the output produced by the 
// 'StThirdMaker' class.

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
static const Bool_t   useQACut = false;
static const Double_t pi       = TMath::Pi();
// i/o parameters
static const TString sIn("merge.pp200r9.d26m3y2017.root");
static const TString sOut("tuple.noCut.pp200r9.d26m3y2017.root");
static const TString sTupleEvt0("CAB");
static const TString sTupleEvt1("QAAB");
static const TString sTupleEvt2("AngleTAB");
static const TString sTupleEvt3("PiB");
static const TString sTupleEvt4("SMDE1B");
static const TString sTupleTrk0("Track");
static const TString sTupleTrk1("QAA");
static const TString sTupleTrk2("CA");
static const TString sTupleTrk3("New8A");
static const TString sTupleTrk4("Pi");
static const TString sTupleTrk5("SMDE1");
static const TString sTupleTrk6("AngleTA");


// QA cuts
static const Int_t    nFitMin  = 15;
static const Double_t strMin   = 0.5;
static const Double_t detHmax  = 0.9;
static const Double_t adcMax   = 6004;
static const Double_t twrPmax  = 3.0;
static const Double_t tspPmin  = 0.;
static const Double_t tspPmax  = 0.08;
static const Double_t tspGmin  = 0.2;
static const Double_t tspGmax  = 0.6;
static const Double_t eTmin    = 8.0;
static const Double_t nRatMin  = 0.52;
static const Double_t dcaMax   = 3.0;
static const Double_t pTmin    = 0.2;
static const Double_t etaMax   = 1.0;



void ReadTuples(const Bool_t inBatchMode=false) {

  cout << "\n  Beginning tuple-reading script..." << endl;
  gErrorIgnoreLevel = kError;


  TFile  *fOut = new TFile(sOut.Data(), "recreate");
  TFile  *fIn  = new TFile(sIn.Data(), "read");
  if (!fIn) {
    cerr << "PANIC: couldn't open input file!" << endl;
    assert(fIn);
  }

  TNtuple *tEvt0 = (TNtuple*) fIn -> Get(sTupleEvt0.Data());
  TNtuple *tEvt1 = (TNtuple*) fIn -> Get(sTupleEvt1.Data());
  TNtuple *tEvt2 = (TNtuple*) fIn -> Get(sTupleEvt2.Data());
  TNtuple *tEvt3 = (TNtuple*) fIn -> Get(sTupleEvt3.Data());
  TNtuple *tEvt4 = (TNtuple*) fIn -> Get(sTupleEvt4.Data());
  TNtuple *tTrk0 = (TNtuple*) fIn -> Get(sTupleTrk0.Data());
  TNtuple *tTrk1 = (TNtuple*) fIn -> Get(sTupleTrk1.Data());
  TNtuple *tTrk2 = (TNtuple*) fIn -> Get(sTupleTrk2.Data());
  TNtuple *tTrk3 = (TNtuple*) fIn -> Get(sTupleTrk3.Data());
  TNtuple *tTrk4 = (TNtuple*) fIn -> Get(sTupleTrk4.Data());
  TNtuple *tTrk5 = (TNtuple*) fIn -> Get(sTupleTrk5.Data());
  TNtuple *tTrk6 = (TNtuple*) fIn -> Get(sTupleTrk6.Data());
  if (!tEvt0 || !tEvt1 || !tEvt2 || !tEvt3 || !tEvt4) {
    cerr << "PANIC: couldn't grab an input event tuple!" << endl;
    assert(tEvt0);
    assert(tEvt1);
    assert(tEvt2);
    assert(tEvt3);
    assert(tEvt4);
  }
  if (!tTrk0 || !tTrk1 || !tTrk2 || !tTrk3 || !tTrk4 || !tTrk5 || !tTrk6) {
    cerr << "PANIC: couldn't grab an input track tuple!" << endl;
    assert(tTrk0);
    assert(tTrk1);
    assert(tTrk2);
    assert(tTrk3);
    assert(tTrk4);
    assert(tTrk5);
    assert(tTrk6);
  }

  cout << "    Tuples grabbed..." << endl;

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
  // Track_M leaves
  Float_t pTPCM;
  Float_t ptTPCM;
  Float_t dEdxM;
  Float_t EtaTrM;
  Float_t epsiM;
  Float_t chM;
  Float_t dcaM;
  Float_t FpM;
  Float_t eEM;
  Float_t phEM;
  Float_t ETAPM;
  Float_t PHIPM;
  Float_t eTM;
  Float_t phTM;
  // QAA_M leaves
  Float_t EneT0M;
  Float_t MPosM;
  Float_t MNegM;
  Float_t HM;
  Float_t en4M;
  Float_t enp4M;
  Float_t eneT0M;
  Float_t zM;
  Float_t etav1M;
  Float_t phiv1M;
  Float_t didTM;
  Float_t delta1M;
  Float_t dphi1M;
  Float_t Weta1M;
  // CA_M leaves
  Float_t didEM;
  Float_t ENET0M;
  Float_t XXM;
  Float_t YYM;
  Float_t adc11M;
  Float_t didPM;
  Float_t adcE1M;
  Float_t adcP1M;
  Float_t eTM;
  Float_t phTM;
  Float_t hSSM;
  Float_t hpedSSM;
  Float_t hRMSSSM;
  Float_t TC7M;
  // Pi1_M leaves
  Float_t psif0M;
  Float_t psif1M;
  Float_t psifn0M;
  Float_t psifn1M;
  Float_t PpoM;
  Float_t xgM;
  Float_t ygM;
  Float_t zgM;
  Float_t pxM;
  Float_t pyM;
  Float_t pzM;
  // Pi_M leaves
  Float_t en9M;
  Float_t en10M;
  Float_t en11M;
  Float_t enp03M;
  Float_t enp02M;
  Float_t enp01M;
  Float_t enp9M;
  Float_t enp10M;
  Float_t enp11M;
  // SMDE1_M leaves
  Float_t en03M;
  Float_t en0M;
  Float_t en1M;
  Float_t en2M;
  Float_t en3M;
  Float_t en5M;
  Float_t en6M;
  Float_t en7M;
  Float_t en8M;
  Float_t enp2M;
  Float_t enp3M;
  Float_t enp5M;
  Float_t enp6M;
  // AngleTA_M leaves
  Float_t moduleTM;
  Float_t EtaEM;
  Float_t sPM;
  Float_t runIdM;
  Float_t evtM;
  Float_t NOPM;
  Float_t Wphi1M;
  Float_t enp7M;
  Float_t enp1M;
  Float_t enp8M;
  Float_t enp0M;
  Float_t en01M;
  Float_t en02M;
  Float_t PTowerM;

  // set CAB addresses
  tEvt0 -> SetBranchAddress("didE", &didEB);
  tEvt0 -> SetBranchAddress("ENET0", &ENET0B);
  tEvt0 -> SetBranchAddress("XX", &XXB);
  tEvt0 -> SetBranchAddress("YY", &YYB);
  tEvt0 -> SetBranchAddress("adc11", &adc11B);
  tEvt0 -> SetBranchAddress("didP", &didPB);
  tEvt0 -> SetBranchAddress("adcE1", &adcE1B);
  tEvt0 -> SetBranchAddress("adcP1", &adcP1B);
  tEvt0 -> SetBranchAddress("eT", &eTB);
  tEvt0 -> SetBranchAddress("phT", &phTB);
  tEvt0 -> SetBranchAddress("hSS", &hSSB);
  tEvt0 -> SetBranchAddress("hpedSS", &hpedSSB);
  tEvt0 -> SetBranchAddress("hRMSSS", &hRMSSSB);
  tEvt0 -> SetBranchAddress("TC7", &TC7B);
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
  // set PiB addresses
  tEvt3 -> SetBranchAddress("en9", &en9B);
  tEvt3 -> SetBranchAddress("en10", &en10B);
  tEvt3 -> SetBranchAddress("en11", &en11B);
  tEvt3 -> SetBranchAddress("enp03", &enp03B);
  tEvt3 -> SetBranchAddress("enp02", &enp02B);
  tEvt3 -> SetBranchAddress("enp01", &enp01B);
  tEvt3 -> SetBranchAddress("enp9", &enp9B);
  tEvt3 -> SetBranchAddress("enp10", &enp10B);
  tEvt3 -> SetBranchAddress("enp11", &enp11B);
  // set SMDE1B addresses
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
  // set Track_M addresses
  tTrk0 -> SetBranchAddress("pTPC", &pTPCM);
  tTrk0 -> SetBranchAddress("ptTPC", &ptTPCM);
  tTrk0 -> SetBranchAddress("dEdx", &dEdxM);
  tTrk0 -> SetBranchAddress("EtaTr", &EtaTrM);
  tTrk0 -> SetBranchAddress("epsi", &epsiM);
  tTrk0 -> SetBranchAddress("ch", &chM);
  tTrk0 -> SetBranchAddress("dca", &dcaM);
  tTrk0 -> SetBranchAddress("Fp", &FpM);
  tTrk0 -> SetBranchAddress("eE", &eEM);
  tTrk0 -> SetBranchAddress("phE", &phEM);
  tTrk0 -> SetBranchAddress("ETAP", &ETAPM);
  tTrk0 -> SetBranchAddress("PHIP", &PHIPM);
  tTrk0 -> SetBranchAddress("eT", &eTM);
  tTrk0 -> SetBranchAddress("phT", &phTM);
  // set QAA_M addresses
  tTrk1 -> SetBranchAddress("EneT0", &EneT0M);
  tTrk1 -> SetBranchAddress("MPos", &MPosM);
  tTrk1 -> SetBranchAddress("MNeg", &MNegM);
  tTrk1 -> SetBranchAddress("H", &HM);
  tTrk1 -> SetBranchAddress("en4", &en4M);
  tTrk1 -> SetBranchAddress("enp4", &enp4M);
  tTrk1 -> SetBranchAddress("eneT0", &eneT0M);
  tTrk1 -> SetBranchAddress("z", &zM);
  tTrk1 -> SetBranchAddress("etav1", &etav1M);
  tTrk1 -> SetBranchAddress("phiv1", &phiv1M);
  tTrk1 -> SetBranchAddress("didT", &didTM);
  tTrk1 -> SetBranchAddress("delta1", &delta1M);
  tTrk1 -> SetBranchAddress("dphi1", &dphi1M);
  tTrk1 -> SetBranchAddress("Weta1", &Weta1M);
  // set CA_M addresses
  tTrk2 -> SetBranchAddress("didE", &didEM);
  tTrk2 -> SetBranchAddress("ENET0", &ENET0M);
  tTrk2 -> SetBranchAddress("XX", &XXM);
  tTrk2 -> SetBranchAddress("YY", &YYM);
  tTrk2 -> SetBranchAddress("adc11", &adc11M);
  tTrk2 -> SetBranchAddress("didP", &didPM);
  tTrk2 -> SetBranchAddress("adcE1", &adcE1M);
  tTrk2 -> SetBranchAddress("adcP1", &adcP1M);
  tTrk2 -> SetBranchAddress("eT", &eTM);
  tTrk2 -> SetBranchAddress("phT", &phTM);
  tTrk2 -> SetBranchAddress("hSS", &hSSM);
  tTrk2 -> SetBranchAddress("hpedSS", &hpedSSM);
  tTrk2 -> SetBranchAddress("hRMSSS", &hRMSSSM);
  tTrk2 -> SetBranchAddress("TC7", &TC7M);
  // set New8A_M addresses
  tTrk3 -> SetBranchAddress("psif0", &psif0M);
  tTrk3 -> SetBranchAddress("psif1", &psif1M);
  tTrk3 -> SetBranchAddress("psifn0", &psifn0M);
  tTrk3 -> SetBranchAddress("psifn1", &psifn1M);
  tTrk3 -> SetBranchAddress("Ppo", &PpoM);
  tTrk3 -> SetBranchAddress("xg", &xgM);
  tTrk3 -> SetBranchAddress("yg", &ygM);
  tTrk3 -> SetBranchAddress("zg", &zgM);
  tTrk3 -> SetBranchAddress("px", &pxM);
  tTrk3 -> SetBranchAddress("py", &pyM);
  tTrk3 -> SetBranchAddress("pz", &pzM);
  // set Pi_M addresses
  tTrk4 -> SetBranchAddress("en9", &en9M);
  tTrk4 -> SetBranchAddress("en10", &en10M);
  tTrk4 -> SetBranchAddress("en11", &en11M);
  tTrk4 -> SetBranchAddress("enp03", &enp03M);
  tTrk4 -> SetBranchAddress("enp02", &enp02M);
  tTrk4 -> SetBranchAddress("enp01", &enp01M);
  tTrk4 -> SetBranchAddress("enp9", &enp9M);
  tTrk4 -> SetBranchAddress("enp10", &enp10M);
  tTrk4 -> SetBranchAddress("enp11", &enp11M);
  // set SMDE1_M addresses
  tTrk5 -> SetBranchAddress("en03", &en03M);
  tTrk5 -> SetBranchAddress("en0", &en0M);
  tTrk5 -> SetBranchAddress("en1", &en1M);
  tTrk5 -> SetBranchAddress("en2", &en2M);
  tTrk5 -> SetBranchAddress("en3", &en3M);
  tTrk5 -> SetBranchAddress("en5", &en5M);
  tTrk5 -> SetBranchAddress("en6", &en6M);
  tTrk5 -> SetBranchAddress("en7", &en7M);
  tTrk5 -> SetBranchAddress("en8", &en8M);
  tTrk5 -> SetBranchAddress("enp2", &enp2M);
  tTrk5 -> SetBranchAddress("enp3", &enp3M);
  tTrk5 -> SetBranchAddress("enp5", &enp5M);
  tTrk5 -> SetBranchAddress("enp6", &enp6M);
  // set AngleTA_M addresses
  tTrk6 -> SetBranchAddress("moduleT", &moduleTM);
  tTrk6 -> SetBranchAddress("EtaE", &EtaEM);
  tTrk6 -> SetBranchAddress("sP", &sPM);
  tTrk6 -> SetBranchAddress("runId", &runIdM);
  tTrk6 -> SetBranchAddress("evt", &evtM);
  tTrk6 -> SetBranchAddress("NOP", &NOPM);
  tTrk6 -> SetBranchAddress("Wphi1", &Wphi1M);
  tTrk6 -> SetBranchAddress("enp7", &enp7M);
  tTrk6 -> SetBranchAddress("enp1", &enp1M);
  tTrk6 -> SetBranchAddress("enp8", &enp8M);
  tTrk6 -> SetBranchAddress("enp0", &enp0M);
  tTrk6 -> SetBranchAddress("en01", &en01M);
  tTrk6 -> SetBranchAddress("en02", &en02M);
  tTrk6 -> SetBranchAddress("PTower", &PTowerM);

  cout << "    Branches set..." << endl;


  // initialize histograms
  const Int_t    nVz = 2000;
  const Int_t    nTr = 50;
  const Int_t    nPt = 500;
  const Int_t    nF  = 720;
  const Int_t    nH  = 100;
  const Int_t    nDf = 360;
  const Double_t vZ1 = -100.;
  const Double_t vZ2 = 100.;
  const Double_t tr1 = 0.;
  const Double_t tr2 = 50.;
  const Double_t pT1 = 0.;
  const Double_t pT2 = 50.;
  const Double_t f1  = -2. * pi;
  const Double_t f2  = 2. * pi;
  const Double_t h1  = -5.;
  const Double_t h2  = 5.;
  const Double_t dF1 = -1. * (pi / 2.);
  const Double_t dF2 = 3. * (pi / 2.);
  TH1D *hPrimeVz = new TH1D("hPrimeVz", "Primary vtx., z-coordinate", nVz, vZ1, vZ2);
  TH1D *hEtTrg   = new TH1D("hEtTrg", "Trigger E_{T}", nPt, pT1, pT2);
  TH1D *hFtrg    = new TH1D("hFtrg", "Trigger phi", nF, f1, f2);
  TH1D *hHtrg    = new TH1D("hHtrg", "Trigger eta", nH, h1, h2);
  TH1D *hPtTrk   = new TH1D("hPtTrk", "Primary track pT (trigger cuts applied)", nPt, pT1, pT2);
  TH1D *hDfTrk   = new TH1D("hDfTrk", "Primary track delta-phi (trigger cuts applied)", nDf, dF1, dF2);
  TH1D *hHtrk    = new TH1D("hHtrk", "Primary track eta (trigger cuts applied)", nH, h1, h2);
  hPrimeVz -> Sumw2();
  hEtTrg   -> Sumw2();
  hFtrg    -> Sumw2();
  hHtrg    -> Sumw2();
  hPtTrk   -> Sumw2();
  hDfTrk   -> Sumw2();
  hHtrk    -> Sumw2();

  cout << "    Histograms initialized..." << endl;


  // initialize vectors for graphs
  vector<Int_t>    runs;
  vector<Int_t>    evts;
  vector<Int_t>    index;
  vector<Double_t> vZ;


  cout << "    Beginning event loop..." << endl;

  Int_t nTrgs = 0;
  Int_t nEvts = 0;
  Int_t nEvt0 = tEvt0 -> GetEntries();
  Int_t nEvt1 = tEvt1 -> GetEntries();
  Int_t nEvt2 = tEvt2 -> GetEntries();
  Int_t nEvt3 = tEvt3 -> GetEntries();
  Int_t nEvt4 = tEvt4 -> GetEntries();
  if ((nEvt0 != nEvt1) || (nEvt1 != nEvt2) || (nEvt2 != nEvt3) || (nEvt3 != nEvt4)) {
    cerr << "PANIC: event tuples do not have\n"
         << "       the same number of entries!"
         << endl;
    assert(nEvt0 == nEvt1);
    assert(nEvt1 == nEvt2);
    assert(nEvt2 == nEvt3);
    assert(nEvt3 == nEvt4);
  }
  else
    nEvts = nEvt1;

  // event loop
  for (Int_t i = 0; i < nEvts; i++) {

    tEvt0 -> GetEntry(i);
    tEvt1 -> GetEntry(i);
    tEvt2 -> GetEntry(i);
    tEvt3 -> GetEntry(i);
    tEvt4 -> GetEntry(i);
    if (inBatchMode) {
      cout << "      Processing event " << (i + 1) << "/" << nEvts << "..." << endl;
    }
    else {
      cout << "      Processing event " << (i + 1) << "/" << nEvts << "...\r" << flush;
      if ((i + 1) == nEvts) cout << endl;
    }


    Double_t fTrg  = phiv1B;
    Double_t hTrg  = etav1B;
    Double_t tTrg  = 2. * atan(exp(-1. * hTrg));
    Double_t eTrg  = EneT0B;
    Double_t eTtrg = eTrg * sin(tTrg);

    // calculate TSP
    Double_t TSPn = eneT0B;
    Double_t TSPd = (en4B + enp4B) * 0.783 + (en3B + enp3B + en5B + enp5B) * 2.21 + (en2B + enp2B + en6B + enp6B) * 6.26
                  + (en1B + enp1B + en7B + enp7B) * 11.51 + (en0B + enp0B + en8B + enp8B) * 17.73 + (en01B + enp01B + en9B + enp9B) * 24.78
                  + (en02B + enp02B + en10B + enp10B) * 32.57 + (en03B + enp03B + en11B + enp11B) * 41.05;
    if (TSPd == 0) TSPd = -1.;
    Double_t TSP  = TSPn / TSPd;

    // apply trigger cuts, if necessary
    const Bool_t inStrCut = ((en4B >= strMin) && (enp4B >= strMin));
    const Bool_t inDetCut = (abs(eTB) < detHmax);
    const Bool_t inAdcCut = (adc11B <= adcMax);
    const Bool_t inTwrCut = (PTowerB < twrPmax);
    const Bool_t inPi0Cut = ((TSP > tspPmin) && (TSP < tspPmax));
    const Bool_t inGamCut = ((TSP > tspGmin) && (TSP < tspGmax));
    const Bool_t inEtCut  = (eTtrg > eTmin);
    if (useQACut) {
      if (!inStrCut || !inDetCut)
        continue;
      if (!inAdcCut || !inTwrCut)
        continue;
      if (!inPi0Cut && !inGamCut)
        continue;
      if (!inEtCut)
        continue;
    }

    // store event info
    vZ.push_back(zB);
    runs.push_back(runIdB);
    evts.push_back(evtB);
    index.push_back(i);

    hPrimeVz -> Fill(zB);
    hEtTrg   -> Fill(eTtrg);
    hFtrg    -> Fill(fTrg);
    hHtrg    -> Fill(hTrg);
    ++nTrgs;

  }  // end event loop

  cout << "    Event loop finished!\n"
       << "    Beginning track loop..."
       << endl;


  Int_t nTrks = 0;
  Int_t nTrk0 = tTrk0 -> GetEntries();
  Int_t nTrk1 = tTrk1 -> GetEntries();
  Int_t nTrk2 = tTrk2 -> GetEntries();
  Int_t nTrk3 = tTrk3 -> GetEntries();
  Int_t nTrk4 = tTrk4 -> GetEntries();
  Int_t nTrk5 = tTrk5 -> GetEntries();
  Int_t nTrk6 = tTrk6 -> GetEntries();
  if ((nTrk0 != nTrk1) || (nTrk1 != nTrk2) || (nTrk2 != nTrk3) || (nTrk3 != nTrk4) || (nTrk4 != nTrk5) || (nTrk5 != nTrk6)) {
    cerr << "PANIC: the track tuples do not have\n"
         << "       the same number of entries!"
         << endl;
    assert(nTrk0 == nTrk1);
    assert(nTrk1 == nTrk2);
    assert(nTrk2 == nTrk3);
    assert(nTrk3 == nTrk4);
    assert(nTrk4 == nTrk5);
    assert(nTrk5 == nTrk6);
  }
  else
    nTrks = nTrk1;

  // track loop
  for (Int_t i = 0; i < nTrks; i++) {

    tTrk0 -> GetEntry(i);
    tTrk1 -> GetEntry(i);
    tTrk2 -> GetEntry(i);
    tTrk3 -> GetEntry(i);
    tTrk4 -> GetEntry(i);
    tTrk5 -> GetEntry(i);
    tTrk6 -> GetEntry(i);
    if (inBatchMode) {
      cout << "      Processing track " << (i + 1) << "/" << nTrks << "..." << endl;
    }
    else {
      cout << "      Processing track " << (i + 1) << "/" << nTrks << "...\r" << flush;
      if ((i + 1) == nTrks) cout << endl;
    }


    Double_t hTrg  = etav1M;
    Double_t fTrg  = phiv1M;
    Double_t tTrg  = 2. * atan(exp(-1. * hTrg));
    Double_t eTrg  = EneT0M;
    Double_t eTtrg = eTrg * sin(tTrg);
    Double_t pTtrk = ptTPCM;
    Double_t hTrk  = EtaTrM;
    Double_t fTrk  = epsiM;
    Double_t fTrg  = phiv1M;
    Double_t dFtrk = fTrk - fTrg;
    if (dFtrk < dF1) dFtrk += (2. * pi);
    if (dFtrk > dF2) dFtrk -= (2. * pi);

    // calculate TSP
    Double_t TSPn = eneT0M;
    Double_t TSPd = (en4M + enp4M) * 0.783 + (en3M + enp3M + en5M + enp5M) * 2.21 + (en2M + enp2M + en6M + enp6M) * 6.26
                  + (en1M + enp1M + en7M + enp7M) * 11.51 + (en0M + enp0M + en8M + enp8M) * 17.73 + (en01M + enp01M + en9M + enp9M) * 24.78
                  + (en02M + enp02M + en10M + enp10M) * 32.57 + (en03M + enp03M + en11M + enp11M) * 41.05;
    if (TSPd == 0) TSPd = -1.;
    Double_t TSP  = TSPn / TSPd;

    // apply trigger cuts, if necessary
    const Bool_t inStrCut = ((en4M >= strMin) && (enp4M >= strMin));
    const Bool_t inDetCut = (abs(eTM) < detHmax);
    const Bool_t inAdcCut = (adc11M <= adcMax);
    const Bool_t inTwrCut = (PTowerM < twrPmax);
    const Bool_t inPi0Cut = ((TSP > tspPmin) && (TSP < tspPmax));
    const Bool_t inGamCut = ((TSP > tspGmin) && (TSP < tspGmax));
    const Bool_t inEtCut  = (eTtrg > eTmin);
    if (useQACut) {
      if (!inStrCut || !inDetCut)
        continue;
      if (!inAdcCut || !inTwrCut)
        continue;
      if (!inPi0Cut && !inGamCut)
        continue;
      if (!inEtCut)
        continue;
    }

    // apply track cuts, if necessary
    const Bool_t inFitCut = (FpM > nFitMin);
    const Bool_t inRatCut = ((FpM / PpoM) > nRatMin);
    const Bool_t inDcaCut = (dcaM < dcaMax);
    const Bool_t inPtCut  = (pTtrk > pTmin);
    const Bool_t inEtaCut = (abs(hTrk) < etaMax);
    if (useQACut) {
      if (!inFitCut)
        continue;
      if (!inRatCut)
        continue;
      if (!inDcaCut)
        continue;
      if (!inPtCut)
        continue;
      if (!inEtaCut)
        continue;
    }

    // store track info
    hPtTrk -> Fill(pTtrk);
    hDfTrk -> Fill(dFtrk);
    hHtrk  -> Fill(hTrk);

  }  // end track loop

  cout << "    Track loop finished!" << endl;


  // create graphs
  const Int_t    nPts = (Int_t) index.size();
  const Int_t    run1 = runs[0] - 1;
  const Int_t    run2 = runs[nPts - 1] + 1;
  const Int_t    nRun = run2 - run1;
  const Int_t    evt1 = evts[0];
  const Int_t    evt2 = evts[nPts - 1];
  const Int_t    ind1 = index[0];
  const Int_t    ind2 = index[nPts - 1];
  TH1D *hVzVsRun  = new TH1D("hVzVsRun", "V_{z} vs run no.", nRun, run1, run2);
  TH1D *hVzVsEvt  = new TH1D("hVzVsEvt", "V_{z} vs evt no.", nPts, evt1, evt2);
  TH1D *hVzVsInd  = new TH1D("hVzVsInd", "V_{z} vs index", nPts, ind1, ind2);
  TH1D *hRunVsInd = new TH1D("hRunVsInd", "Run no. vs index", nPts, ind1, ind2);
  TH1D *hEvtVsInd = new TH1D("hEvtVsInd", "Evt no. vs index", nPts, ind1, ind2);
  hVzVsRun  -> Sumw2();
  hVzVsEvt  -> Sumw2();
  hVzVsInd  -> Sumw2();
  hRunVsInd -> Sumw2();
  hEvtVsInd -> Sumw2();
  for (Int_t i = 0; i < nPts; i++) {
    Int_t iRun = hVzVsRun -> FindBin(runs[i]);
    Int_t iEvt = hVzVsEvt -> FindBin(evts[i]);
    Int_t iInd = hVzVsInd -> FindBin(index[i]);
    hVzVsRun  -> SetBinContent(iRun, vZ[i]);
    hVzVsRun  -> SetBinError(iRun, 0);
    hVzVsEvt  -> SetBinContent(iEvt, vZ[i]);
    hVzVsEvt  -> SetBinError(iEvt, 0);
    hVzVsInd  -> SetBinContent(iInd, vZ[i]);
    hVzVsInd  -> SetBinError(iInd, 0);
    hRunVsInd -> SetBinContent(iInd, runs[i]);
    hRunVsInd -> SetBinError(iInd, 0);
    hEvtVsInd -> SetBinContent(iInd, evts[i]);
    hEvtVsInd -> SetBinError(iInd, 0);
  }

  cout << "    Graphs created..." << endl;


  // normalize histograms
  const Double_t pTbin  = (pT2 - pT1) / nPt;
  const Double_t dFbin  = (dF2 - dF1) / nDf;
  const Double_t hBin   = (h2 - h1) / nH;
  const Double_t pTnorm = pTbin * nTrgs;
  const Double_t dFnorm = dFbin * nTrgs;
  const Double_t hNorm  = hBin * nTrgs;
  hPrimeVz -> Scale(1.);
  hEtTrg   -> Scale(1.);
  hFtrg    -> Scale(1.);
  hHtrg    -> Scale(1.);
  if (nTrgs != 0.) {
    hPtTrk   -> Scale(1. / pTnorm);
    hDfTrk   -> Scale(1. / dFnorm);
    hHtrk    -> Scale(1. / hNorm);
  }


  // close and save
  fOut      -> cd();
  hPrimeVz  -> Write();
  hEtTrg    -> Write();
  hFtrg     -> Write();
  hHtrg     -> Write();
  hPtTrk    -> Write();
  hDfTrk    -> Write();
  hHtrk     -> Write();
  hVzVsRun  -> Write();
  hVzVsEvt  -> Write();
  hVzVsInd  -> Write();
  hRunVsInd -> Write();
  hEvtVsInd -> Write();
  fOut      -> Close();
  fIn       -> cd();
  fIn       -> Close();

  cout << "  Tuple reading script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
