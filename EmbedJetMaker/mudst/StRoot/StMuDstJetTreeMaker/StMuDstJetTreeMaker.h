// 'StMuDstJetTreeMaker.h'
// Derek Anderson, Nihar Sahoo
// 06.23.2017
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets.


#ifndef StMuDstJetTreeMaker_h
#define StMuDstJetTreeMaker_h

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

// TEST [04.04.2020]
#include "TH2.h"

using namespace std;


// global constants
static const UInt_t   NHistQA    = 4;
static const UInt_t   NTrgTypes  = 4;
static const UInt_t   NTrkTypes  = 2;
static const UInt_t   NBadRun    = 45;
static const UInt_t   NHotTwr    = 302;
static const UInt_t   NTrkMax    = 5000;
static const UInt_t   NTwrMax    = 5000;
static const UInt_t   NMatchMax  = 10;
static const UInt_t   NSigCut    = 3;
static const Double_t MassPi0    = 0.140;
static const Double_t RadiusBemc = 225.405;



class StMuDstJetTreeMaker {

  private:

    Int_t    _fCurrent;
    UInt_t   _jetType;
    UInt_t   _trgFlag;
    Bool_t   _isInBatchMode;
    Bool_t   _adjustTrackEff;
    TH1D    *_hEvtQA[NHistQA][NTrgTypes];
    TH1D    *_hTrkQA[NHistQA][NTrgTypes];
    TH1D    *_hTwrQA[NHistQA][NTrgTypes];
    TH1D    *_hJetQA[NHistQA][NTrgTypes];
    TH1D    *_hEvtPt[NTrgTypes][NTrkTypes];
    TTree   *_tFemto;
    TTree   *_tJet;
    TFile   *_fOutput;
    Float_t  _effAdjust;
    TString  _sInput;
    TString  _sOutput;
    TRandom *_random;

    // event parameters
    Int_t    _adcMax;
    Double_t _rVtxMax;
    Double_t _zVtxMax;
    Double_t _eEtaMin;
    Double_t _ePhiMin;
    Double_t _pProjMax;
    Double_t _etaTrgMax;
    Double_t _eTtrgMin;
    Double_t _eTtrgMax;
    Double_t _tspPi0Min;
    Double_t _tspPi0Max; 
    Double_t _tspGamMin;
    Double_t _tspGamMax;
    // track parameters
    UInt_t   _nFitMin;
    Double_t _rFitMin;
    Double_t _dcaMax;
    Double_t _etaTrkMax;
    Double_t _pTtrkMin;
    Double_t _pTtrkMax;
    // tower parameters
    Double_t _etaTwrMax;
    Double_t _eTwrMin;
    Double_t _eTwrMax;
    Double_t _eCorrMin;
    Double_t _eCorrMax;
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
    // output cst. members
    vector<vector<Double_t>> _JetConsPt;
    vector<vector<Double_t>> _JetConsEta;
    vector<vector<Double_t>> _JetConsPhi;
    vector<vector<Double_t>> _JetConsE;

    // input leaf types
    UInt_t   _fUniqueID;
    UInt_t   _fBits;
    Long64_t _runNumber;
    Long64_t _eventNumber;
    Int_t    _trigID;
    Int_t    _nGlobalTracks;
    Int_t    _nPrimaryTracks;
    Int_t    _refMult;
    Double_t _vpdVz;
    Double_t _xVertex;
    Double_t _yVertex;
    Double_t _zVertex;
    Double_t _bbcZVertex;
    Double_t _zdcCoincidenceRate;
    Double_t _bbcCoincidenceRate;
    Double_t _backgroundRate;
    Double_t _bbcBlueBackgroundRate;
    Double_t _bbcYellowBackgroundRate;
    Double_t _refMultPos;
    Double_t _refMultNeg;
    Double_t _bTOFTrayMultiplicity;
    Int_t    _nVerticies;
    Double_t _MagF;
    Double_t _VrtxRank;
    Int_t    _FlagEvent_TrgTrkMisMtch;
    Float_t  _Etsp;
    Int_t    _ETwrdidT;
    Int_t    _ETwradc11;
    Float_t  _ETwreneT0;
    Float_t  _ETwreT;
    Float_t  _ETwrENET0;
    Float_t  _ETwrphT;
    Float_t  _ETwrPTower;
    Float_t  _ETwrpidTower;
    Int_t    _ETwrmoduleT;
    Float_t  _EClustEneT0;
    Float_t  _EClustetav1;
    Float_t  _EClustphiv1;
    Float_t  _EEstrpen01;
    Float_t  _EEstrpen02;
    Float_t  _EEstrpen03;
    Float_t  _EEstrpen0;
    Float_t  _EEstrpen1;
    Float_t  _EEstrpen2;
    Float_t  _EEstrpen3;
    Float_t  _EEstrpen4;
    Float_t  _EEstrpen5;
    Float_t  _EEstrpen6;
    Float_t  _EEstrpen7;
    Float_t  _EEstrpen8;
    Float_t  _EEstrpen9;
    Float_t  _EEstrpen10;
    Float_t  _EEstrpen11;
    Float_t  _EEstrpen12;
    Float_t  _EEstrpen13;
    Float_t  _EEstrpen14;
    Float_t  _EEstrpen15;
    Int_t    _ETwrdidE;
    Float_t  _EPstripenp01;
    Float_t  _EPstripenp02;
    Float_t  _EPstripenp03;
    Float_t  _EPstripenp0;
    Float_t  _EPstripenp1;
    Float_t  _EPstripenp2;
    Float_t  _EPstripenp3;
    Float_t  _EPstripenp4;
    Float_t  _EPstripenp5;
    Float_t  _EPstripenp6;
    Float_t  _EPstripenp7;
    Float_t  _EPstripenp8;
    Float_t  _EPstripenp9;
    Float_t  _EPstripenp10;
    Float_t  _EPstripenp11;
    Float_t  _EPstripenp12;
    Float_t  _EPstripenp13;
    Float_t  _EPstripenp14;
    Float_t  _EPstripenp15;
    Float_t  _EclustEnnq1;
    Float_t  _EclustEnnq20;
    Float_t  _EclustEnnq19;
    Float_t  _EclustEnpq1;
    Float_t  _EclustEnpq20;
    Float_t  _EclustEnpq19;
    Float_t  _EclustEnpq21;
    Int_t    _PrimaryTrackArray_;
    UInt_t   _PrimaryTrackArray_fUniqueID[NTrkMax];
    UInt_t   _PrimaryTrackArray_fBits[NTrkMax];
    Int_t    _PrimaryTrackArray_nHitsFit[NTrkMax];
    Int_t    _PrimaryTrackArray_nHitsPoss[NTrkMax];
    Int_t    _PrimaryTrackArray_trackFlag[NTrkMax];
    Double_t _PrimaryTrackArray_pZ[NTrkMax];
    Double_t _PrimaryTrackArray_pX[NTrkMax];
    Double_t _PrimaryTrackArray_pY[NTrkMax];
    Double_t _PrimaryTrackArray_pT[NTrkMax];
    Double_t _PrimaryTrackArray_dEdx[NTrkMax];
    Double_t _PrimaryTrackArray_charge[NTrkMax];
    Double_t _PrimaryTrackArray_tofBeta[NTrkMax];
    Double_t _PrimaryTrackArray_eta[NTrkMax];
    Double_t _PrimaryTrackArray_phi[NTrkMax];
    Double_t _PrimaryTrackArray_nSigElectron[NTrkMax];
    Double_t _PrimaryTrackArray_nSigPion[NTrkMax];
    Double_t _PrimaryTrackArray_nSigKaon[NTrkMax];
    Double_t _PrimaryTrackArray_nSigProton[NTrkMax];
    Double_t _PrimaryTrackArray_dcag[NTrkMax];
    Double_t _PrimaryTrackArray_nHits[NTrkMax];
    Double_t _PrimaryTrackArray_dEdxHits[NTrkMax];
    Double_t _PrimaryTrackArray_firstZPoint[NTrkMax];
    Double_t _PrimaryTrackArray_lastZPoint[NTrkMax];
    Double_t _PrimaryTrackArray_tofSigElectron[NTrkMax];
    Double_t _PrimaryTrackArray_tofSigPion[NTrkMax];
    Double_t _PrimaryTrackArray_tofSigKaon[NTrkMax];
    Double_t _PrimaryTrackArray_tofSigProton[NTrkMax];
    Double_t _PrimaryTrackArray_timeOfflight[NTrkMax];
    Double_t _PrimaryTrackArray_pathLength[NTrkMax];
    Int_t    _PrimaryTrackArray_trkIndex[NTrkMax];
    Int_t    _TowerArray_;
    UInt_t   _TowerArray_fUniqueID[NTwrMax];
    UInt_t   _TowerArray_fBits[NTwrMax];
    Int_t    _TowerArray_TwrId[NTwrMax];
    Float_t  _TowerArray_TwrEng[NTwrMax];
    Float_t  _TowerArray_TwrEta[NTwrMax];
    Float_t  _TowerArray_TwrPhi[NTwrMax];
    Float_t  _TowerArray_TwrADC[NTwrMax];
    Float_t  _TowerArray_TwrPed[NTwrMax];
    Float_t  _TowerArray_TwrRMS[NTwrMax];
    Int_t    _TowerArray_TwrMatchIdnex[NTwrMax];
    Int_t    _TowerArray_NoOfmatchedTrk[NTwrMax];
    Float_t  _TowerArray_TwrMatchP[NTwrMax];
    Float_t  _TowerArray_TwrPx[NTwrMax];
    Float_t  _TowerArray_TwrPy[NTwrMax];
    Float_t  _TowerArray_TwrPz[NTwrMax];
    Int_t    _TowerArray_fNAssocTracks[NTwrMax];
    Int_t    _TowerArray_fMatchedTracksArray_[NTwrMax][NMatchMax];
    Float_t  _TowerArray_fMatchedTracksArray_P[NTwrMax][NMatchMax];
    // input branches
    TBranch *_bEventList_fUniqueID;
    TBranch *_bEventList_fBits;
    TBranch *_bEventList_runNumber;
    TBranch *_bEventList_eventNumber;
    TBranch *_bEventList_trigID;
    TBranch *_bEventList_nGlobalTracks;
    TBranch *_bEventList_nPrimaryTracks;
    TBranch *_bEventList_refMult;
    TBranch *_bEventList_vpdVz;
    TBranch *_bEventList_xVertex;
    TBranch *_bEventList_yVertex;
    TBranch *_bEventList_zVertex;
    TBranch *_bEventList_bbcZVertex;
    TBranch *_bEventList_zdcCoincidenceRate;
    TBranch *_bEventList_bbcCoincidenceRate;
    TBranch *_bEventList_backgroundRate;
    TBranch *_bEventList_bbcBlueBackgroundRate;
    TBranch *_bEventList_bbcYellowBackgroundRate;
    TBranch *_bEventList_refMultPos;
    TBranch *_bEventList_refMultNeg;
    TBranch *_bEventList_bTOFTrayMultiplicity;
    TBranch *_bEventList_nVerticies;
    TBranch *_bEventList_MagF;
    TBranch *_bEventList_VrtxRank;
    TBranch *_bEventList_FlagEvent_TrgTrkMisMtch;
    TBranch *_bEventList_Etsp;
    TBranch *_bEventList_ETwrdidT;
    TBranch *_bEventList_ETwradc11;
    TBranch *_bEventList_ETwreneT0;
    TBranch *_bEventList_ETwreT;
    TBranch *_bEventList_ETwrENET0;
    TBranch *_bEventList_ETwrphT;
    TBranch *_bEventList_ETwrPTower;
    TBranch *_bEventList_ETwrpidTower;
    TBranch *_bEventList_ETwrmoduleT;
    TBranch *_bEventList_EClustEneT0;
    TBranch *_bEventList_EClustetav1;
    TBranch *_bEventList_EClustphiv1;
    TBranch *_bEventList_EEstrpen01;
    TBranch *_bEventList_EEstrpen02;
    TBranch *_bEventList_EEstrpen03;
    TBranch *_bEventList_EEstrpen0;
    TBranch *_bEventList_EEstrpen1;
    TBranch *_bEventList_EEstrpen2;
    TBranch *_bEventList_EEstrpen3;
    TBranch *_bEventList_EEstrpen4;
    TBranch *_bEventList_EEstrpen5;
    TBranch *_bEventList_EEstrpen6;
    TBranch *_bEventList_EEstrpen7;
    TBranch *_bEventList_EEstrpen8;
    TBranch *_bEventList_EEstrpen9;
    TBranch *_bEventList_EEstrpen10;
    TBranch *_bEventList_EEstrpen11;
    TBranch *_bEventList_EEstrpen12;
    TBranch *_bEventList_EEstrpen13;
    TBranch *_bEventList_EEstrpen14;
    TBranch *_bEventList_EEstrpen15;
    TBranch *_bEventList_ETwrdidE;
    TBranch *_bEventList_EPstripenp01;
    TBranch *_bEventList_EPstripenp02;
    TBranch *_bEventList_EPstripenp03;
    TBranch *_bEventList_EPstripenp0;
    TBranch *_bEventList_EPstripenp1;
    TBranch *_bEventList_EPstripenp2;
    TBranch *_bEventList_EPstripenp3;
    TBranch *_bEventList_EPstripenp4;
    TBranch *_bEventList_EPstripenp5;
    TBranch *_bEventList_EPstripenp6;
    TBranch *_bEventList_EPstripenp7;
    TBranch *_bEventList_EPstripenp8;
    TBranch *_bEventList_EPstripenp9;
    TBranch *_bEventList_EPstripenp10;
    TBranch *_bEventList_EPstripenp11;
    TBranch *_bEventList_EPstripenp12;
    TBranch *_bEventList_EPstripenp13;
    TBranch *_bEventList_EPstripenp14;
    TBranch *_bEventList_EPstripenp15;
    TBranch *_bEventList_EclustEnnq1;
    TBranch *_bEventList_EclustEnnq20;
    TBranch *_bEventList_EclustEnnq19;
    TBranch *_bEventList_EclustEnpq1;
    TBranch *_bEventList_EclustEnpq20;
    TBranch *_bEventList_EclustEnpq19;
    TBranch *_bEventList_EclustEnpq21;
    TBranch *_bEventList_PrimaryTrackArray_;
    TBranch *_bPrimaryTrackArray_fUniqueID;
    TBranch *_bPrimaryTrackArray_fBits;
    TBranch *_bPrimaryTrackArray_nHitsFit;
    TBranch *_bPrimaryTrackArray_nHitsPoss;
    TBranch *_bPrimaryTrackArray_trackFlag;
    TBranch *_bPrimaryTrackArray_pZ;
    TBranch *_bPrimaryTrackArray_pX;
    TBranch *_bPrimaryTrackArray_pY;
    TBranch *_bPrimaryTrackArray_pT;
    TBranch *_bPrimaryTrackArray_dEdx;
    TBranch *_bPrimaryTrackArray_charge;
    TBranch *_bPrimaryTrackArray_tofBeta;
    TBranch *_bPrimaryTrackArray_eta;
    TBranch *_bPrimaryTrackArray_phi;
    TBranch *_bPrimaryTrackArray_nSigElectron;
    TBranch *_bPrimaryTrackArray_nSigPion;
    TBranch *_bPrimaryTrackArray_nSigKaon;
    TBranch *_bPrimaryTrackArray_nSigProton;
    TBranch *_bPrimaryTrackArray_dcag;
    TBranch *_bPrimaryTrackArray_nHits;
    TBranch *_bPrimaryTrackArray_dEdxHits;
    TBranch *_bPrimaryTrackArray_firstZPoint;
    TBranch *_bPrimaryTrackArray_lastZPoint;
    TBranch *_bPrimaryTrackArray_tofSigElectron;
    TBranch *_bPrimaryTrackArray_tofSigPion;
    TBranch *_bPrimaryTrackArray_tofSigKaon;
    TBranch *_bPrimaryTrackArray_tofSigProton;
    TBranch *_bPrimaryTrackArray_timeOfflight;
    TBranch *_bPrimaryTrackArray_pathLength;
    TBranch *_bPrimaryTrackArray_trkIndex;
    TBranch *_bEventList_TowerArray_;
    TBranch *_bTowerArray_fUniqueID;
    TBranch *_bTowerArray_fBits;
    TBranch *_bTowerArray_TwrId;
    TBranch *_bTowerArray_TwrEng;
    TBranch *_bTowerArray_TwrEta;
    TBranch *_bTowerArray_TwrPhi;
    TBranch *_bTowerArray_TwrADC;
    TBranch *_bTowerArray_TwrPed;
    TBranch *_bTowerArray_TwrRMS;
    TBranch *_bTowerArray_TwrMatchIdnex;
    TBranch *_bTowerArray_NoOfmatchedTrk;
    TBranch *_bTowerArray_TwrMatchP;
    TBranch *_bTowerArray_TwrPx;
    TBranch *_bTowerArray_TwrPy;
    TBranch *_bTowerArray_TwrPz;
    TBranch *_bTowerArray_fNAssocTracks;
    TBranch *_bTowerArray_fMatchedTracksArray_;
    TBranch *_bTowerArray_fMatchedTracksArray_P;

    // TEST [04.04.2020]
    TH1D *_hZtCstVsPtJet[5];
    TH1D *_hZtCstVsNcst[4];
    TH2D *_hNcstVsPtJet;


    // private methods
    void           InitializeInputTree(TTree *tree);
    void           InitializeOutputTree(TTree *tree);
    void           InitializeHistograms();
    void           PrintInfo(const UInt_t code, const UInt_t nEvts=0, const UInt_t iEvt=0);
    void           UsePionTriggeredEvents();
    void           UseHadronTriggeredEvents();
    void           UseNonTriggeredEvents();
    Bool_t         IsGoodRunID(const UInt_t runID);
    Bool_t         IsGoodTowerID(const UInt_t twrID);
    Bool_t         IsGoodEvent(const Double_t rVtx, const Double_t zVtx);
    Bool_t         IsGoodPionTrigger(const Int_t adc, const Double_t eEta, const Double_t ePhi, const Double_t pProj, const Double_t etaTrg, const Double_t eTtrg, const Double_t tsp);
    Bool_t         IsGoodHadronTrigger(const Double_t etaTrg, const Double_t pTtrg, const Double_t nSigPi, const Double_t nSigK, const Double_t nSigP, const Double_t nSigE);
    Bool_t         IsGoodTrack(const UInt_t nFit, const Double_t rFit, const Double_t dca, const Double_t etaTrk, const Double_t pTtrk);
    Bool_t         IsGoodTower(const Double_t etaTwr, const Double_t eTwr, const Double_t eCorr);
    Bool_t         IsPi0(const Double_t tsp);
    Bool_t         IsGamma(const Double_t tsp);
    Bool_t         IsHadron(const Double_t nSigPi, const Double_t nSigK, const Double_t nSigP, const Double_t nSigE);
    Long64_t       LoadTree(const Long64_t entry);
    Long64_t       GetEntry(const Long64_t entry);
    Double_t       GetHadronicCorrection(const Double_t eTwr, const vector<Double_t> pMatchedTrks);
    TLorentzVector GetTowerMomentumVector(const Double_t rBEMC, const Double_t etaTwr, const Double_t phiTwr, const Double_t eTwr, const TVector3& vtx);


  public:

    // ctor, dtor
    StMuDstJetTreeMaker(const Bool_t batch=0, TTree *tree=0);
    virtual ~StMuDstJetTreeMaker();

    // public methods
    void SetInputAndOutputFiles(const TString& sInput, const TString& sOuput, const Double_t pTparton=0.);
    void SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax);
    void SetTriggerParameters(const Int_t adcMax, const Double_t eEtaMin, const Double_t ePhiMin, const Double_t pProjMax, const Double_t etaTrgMax, const Double_t eTtrgMin, const Double_t eTtrgMax, const Double_t tspPi0Min, const Double_t tspPi0Max, const Double_t tspGamMin, const Double_t tspGamMax);
    void SetTrackParameters(const UInt_t nFitMin, const Double_t rFitMin, const Double_t dcaMax, const Double_t etaTrkMax, const Double_t pTtrkMin, const Double_t pTtrkMax);
    void SetTowerParameters(const Double_t etaTwrMax, const Double_t eTwrMin, const Double_t eTwrMax, const Double_t eCorrMin, const Double_t eCorrMax);
    void SetJetParameters(const UInt_t type, const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t aGhost, const Double_t pTjetMin, const Double_t etaGhostMax, const Double_t etaJetMax, const Double_t etaBkgdMax);
    void AdjustTrackEfficiency(const Bool_t effAdjust, const Float_t adjustment);
    void Init(const UInt_t trgFlag=1);
    void Make();
    void Finish();

  ClassDef(StMuDstJetTreeMaker, 1)

};

#endif
#ifdef StMuDstJetTreeMaker_cxx



StMuDstJetTreeMaker::StMuDstJetTreeMaker(const Bool_t batch, TTree *tree) : _tFemto(0) {

  _tFemto         = tree;
  _isInBatchMode  = batch;
  _adjustTrackEff = false;
  _effAdjust      = 0.;
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
  PrintInfo(0);

}  // end ctor



StMuDstJetTreeMaker::~StMuDstJetTreeMaker() {

  if (!_tFemto) return;
  delete _tFemto -> GetCurrentFile();

}  // end dtor

#endif

// End ------------------------------------------------------------------------
