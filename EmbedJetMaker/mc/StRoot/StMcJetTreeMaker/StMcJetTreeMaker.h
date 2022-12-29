// 'StMcJetTreeMaker.h'
// Derek Anderson, Nihar Sahoo
// 04.30.2018
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets.


#ifndef StMcJetTreeMaker_h
#define StMcJetTreeMaker_h

#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

// TEST [04.02.2020]
#include "TH2.h"

using namespace std;


// global constants
static const UInt_t   NHistQA    = 4;
static const UInt_t   NTrgTypes  = 3;
static const UInt_t   NHadIds    = 18;
static const UInt_t   IdPi0      = 7;
static const UInt_t   IdGamma    = 1;
static const UInt_t   NBadRun    = 45;
static const UInt_t   NHotTwr    = 41;
static const UInt_t   NTrkMax    = 5000;
static const UInt_t   NTwrMax    = 5000;
static const UInt_t   NMatchMax  = 10;
static const Double_t MassPi0    = 0.140;
static const Double_t RadiusBemc = 225.405;



class StMcJetTreeMaker {

  private:

    Int_t    _fCurrent;
    UInt_t   _jetType;
    Bool_t   _isInBatchMode;
    Bool_t   _adjustTrackEff;
    Bool_t   _turnOffDecays;
    TH1D    *_hEvtQA[NHistQA][NTrgTypes];
    TH1D    *_hTrkQA[NHistQA][NTrgTypes];
    TH1D    *_hTwrQA[NHistQA][NTrgTypes];
    TH1D    *_hJetQA[NHistQA][NTrgTypes];
    TTree   *_tMc;
    TTree   *_tJet;
    TFile   *_fOutput;
    Float_t  _effAdjust;
    TString  _sInput;
    TString  _sOutput;
    TRandom *_random;

    // event parameters
    Double_t _rVtxMax;
    Double_t _zVtxMax;
    Double_t _etaTrgMax;
    Double_t _eTtrgMin;
    Double_t _eTtrgMax;
    Double_t _tspPi0Min;
    Double_t _tspPi0Max; 
    Double_t _tspGamMin;
    Double_t _tspGamMax;
    // track parameters
    Double_t _etaTrkMax;
    Double_t _pTtrkMin;
    Double_t _pTtrkMax;
    // jet parameters
    UInt_t   _nRepeat;
    UInt_t   _nRemove;
    Double_t _rJet;
    Double_t _aGhost;
    Double_t _pTjetMin;
    Double_t _etaGhostMax;
    Double_t _etaJetMax;
    Double_t _etaBkgdMax;

    // output event members
    Int_t    _EventIndex;
    Int_t    _NJets;
    Int_t    _RunId;
    Double_t _Refmult;
    Double_t _PartonicPt;
    Double_t _TSP;
    Double_t _TrgEta;
    Double_t _TrgPhi;
    Double_t _TrgEt;
    Double_t _Rho;
    Double_t _Sigma;
    Double_t _Vz;
    // output jet members
    vector<Int_t>    _JetIndex; 
    vector<Int_t>    _JetNCons;
    vector<Double_t> _JetPt;
    vector<Double_t> _JetPtCorr;
    vector<Double_t> _JetEta;
    vector<Double_t> _JetPhi;
    vector<Double_t> _JetE;
    vector<Double_t> _JetArea;
    // cst. members
    vector<vector<Double_t>> _JetConsPt;
    vector<vector<Double_t>> _JetConsEta;
    vector<vector<Double_t>> _JetConsPhi;
    vector<vector<Double_t>> _JetConsE;

    // input (event) leaf types
    Int_t    _mcEventId;
    Int_t    _mcRunId;
    Int_t    _mcNumTrks;
    Double_t _muVx;
    Double_t _muVy;
    Double_t _muVz;
    Double_t _mcVx;
    Double_t _mcVy;
    Double_t _mcVz;
    // input (track) leaf types
    vector<Int_t>    *_mcIdTrk;
    vector<Int_t>    *_mcIdGeant;
    vector<Int_t>    *_mcIdVx;
    vector<Int_t>    *_mcIdVxEnd;
    vector<Int_t>    *_mcIntrVtx;
    vector<Bool_t>   *_mcIsShower;
    vector<Double_t> *_mcCharge;
    vector<Double_t> *_mcRapidity;
    vector<Double_t> *_mcEta;
    vector<Double_t> *_mcPhi;
    vector<Double_t> *_mcPx;
    vector<Double_t> *_mcPy;
    vector<Double_t> *_mcPz;
    vector<Double_t> *_mcPt;
    vector<Double_t> *_mcPtot;
    vector<Double_t> *_mcEnergy;
    // input branches
    TBranch *_bMcEventId;
    TBranch *_bMcRunId;
    TBranch *_bMcNumTrks;
    TBranch *_bMuVx;
    TBranch *_bMuVy;
    TBranch *_bMuVz;
    TBranch *_bMcVx;
    TBranch *_bMcVy;
    TBranch *_bMcVz;
    TBranch *_bMcIdTrk;
    TBranch *_bMcIdGeant;
    TBranch *_bMcIdVx;
    TBranch *_bMcIdVxEnd;
    TBranch *_bMcIntrVtx;
    TBranch *_bMcIsShower;
    TBranch *_bMcCharge;
    TBranch *_bMcRapidity;
    TBranch *_bMcEta;
    TBranch *_bMcPhi;
    TBranch *_bMcPx;
    TBranch *_bMcPy;
    TBranch *_bMcPz;
    TBranch *_bMcPt;
    TBranch *_bMcPtot;
    TBranch *_bMcEnergy;

    // TEST [04.02.2020]
    TH1D *_hZtCstVsPtJetP[5];
    TH1D *_hZtCstVsPtJetG[5];
    TH1D *_hZtCstVsPtJetH[5];
    TH1D *_hZtCstVsNcstP[4];
    TH1D *_hZtCstVsNcstG[4];
    TH1D *_hZtCstVsNcstH[4];
    TH2D *_hNcstVsPtJetP;
    TH2D *_hNcstVsPtJetG;
    TH2D *_hNcstVsPtJetH;

    // private methods
    void     InitializeInputTree(TTree *tree);
    void     InitializeOutputTree(TTree *tree);
    void     InitializeHistograms();
    void     PrintInfo(const UInt_t code, const UInt_t nEvts=0, const UInt_t iEvt=0);
    Bool_t   IsGoodRunID(const UInt_t runID);
    Bool_t   IsGoodEvent(const Double_t rVtx, const Double_t zVtx);
    Bool_t   IsGoodTrigger(const Double_t etaTrg, const Double_t eTtrg, const UInt_t idTrg);
    Bool_t   IsGoodTrack(const Double_t etaTrk, const Double_t pTtrk);
    Bool_t   IsPi0(const UInt_t idGnt);
    Bool_t   IsGamma(const UInt_t idGnt);
    Bool_t   IsHadron(const UInt_t idGnt);
    Bool_t   IsUndecayed(const Int_t idVtx);
    Bool_t   IsFinalState(const Double_t pTpar);
    Long64_t LoadTree(const Long64_t entry);
    Long64_t GetEntry(const Long64_t entry);


  public:

    // ctor, dtor
    StMcJetTreeMaker(const Bool_t batch=0, TTree *tree=0);
    virtual ~StMcJetTreeMaker();

    // public methods
    void SetInputAndOutputFiles(const TString& sInput, const TString& sOuput, const Double_t pTparton=0.);
    void SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax);
    void SetTriggerParameters(const Double_t etaTrgMax, const Double_t eTtrgMin, const Double_t eTtrgMax);
    void SetTrackParameters(const Double_t etaTrkMax, const Double_t pTtrkMin, const Double_t pTtrkMax);
    void SetJetParameters(const UInt_t type, const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t aGhost, const Double_t pTjetMin, const Double_t etaGhostMax, const Double_t etaJetMax, const Double_t etaBkgdMax);
    void AdjustTrackEfficiency(const Bool_t effAdjust, const Float_t adjustment);
    void TurnOffWeakDecays(const Bool_t turnOffDecays);
    void Init();
    void Make(const UInt_t evtFlag=0, const UInt_t trgFlag=0);
    void Finish();

  ClassDef(StMcJetTreeMaker, 1)

};

#endif
#ifdef StMcJetTreeMaker_cxx



StMcJetTreeMaker::StMcJetTreeMaker(const Bool_t batch, TTree *tree) : _tMc(0) {

  _tMc            = tree;
  _isInBatchMode  = batch;
  _adjustTrackEff = false;
  _turnOffDecays  = false;
  _effAdjust      = 0.;

  // clear output vectors
  _JetIndex.clear();
  _JetNCons.clear();
  _JetPt.clear();
  _JetPtCorr.clear();
  _JetEta.clear();
  _JetPhi.clear();
  _JetE.clear();
  _JetArea.clear();  
  _JetConsPt.clear();
  _JetConsEta.clear();
  _JetConsPhi.clear();
  _JetConsE.clear();

  // initialize input
  _mcEventId   = 0;
  _mcRunId     = 0;
  _mcNumTrks   = 0;
  _muVx        = 0.;
  _muVy        = 0.;
  _muVz        = 0.;
  _mcVx        = 0.;
  _mcVy        = 0.;
  _mcVz        = 0.;
  _mcIdTrk     = 0;
  _mcIdGeant   = 0;
  _mcIdVx      = 0;
  _mcIdVxEnd   = 0;
  _mcIntrVtx   = 0;
  _mcIsShower  = 0;
  _mcCharge    = 0;
  _mcRapidity  = 0;
  _mcEta       = 0;
  _mcPhi       = 0;
  _mcPx        = 0;
  _mcPy        = 0;
  _mcPz        = 0;
  _mcPt        = 0;
  _mcPtot      = 0;
  _mcEnergy    = 0;
  _bMcEventId  = 0;
  _bMcRunId    = 0;
  _bMcNumTrks  = 0;
  _bMuVx       = 0;
  _bMuVy       = 0;
  _bMuVz       = 0;
  _bMcVx       = 0;
  _bMcVy       = 0;
  _bMcVz       = 0;
  _bMcIdTrk    = 0;
  _bMcIdGeant  = 0;
  _bMcIdVx     = 0;
  _bMcIdVxEnd  = 0;
  _bMcIntrVtx  = 0;
  _bMcIsShower = 0;
  _bMcCharge   = 0;
  _bMcRapidity = 0;
  _bMcEta      = 0;
  _bMcPhi      = 0;
  _bMcPx       = 0;
  _bMcPy       = 0;
  _bMcPz       = 0;
  _bMcPt       = 0;
  _bMcPtot     = 0;
  _bMcEnergy   = 0;
  PrintInfo(0);

}  // end ctor



StMcJetTreeMaker::~StMcJetTreeMaker() {

  if (!_tMc) return;
  delete _tMc -> GetCurrentFile();

}  // dtor

#endif

// End ------------------------------------------------------------------------
