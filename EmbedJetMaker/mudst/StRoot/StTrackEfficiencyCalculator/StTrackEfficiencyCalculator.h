// 'StTrackEfficiencyCalculator.h'
// Derek Anderson
// 01.30.2020
//
// This class produces trees of jets using
// generated particles and reconstructed
// tracks matched to generated particles.
// Also calculates tracking efficiency.


#ifndef StTrackEfficiencyCalculator_h
#define StTrackEfficiencyCalculator_h

#include <vector>
#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include "TBranch.h"
#include "TProfile.h"
#include "TDirectory.h"

using namespace std;


// global constants
static const UInt_t  NLvls(2);
static const UInt_t  NCuts(2);
static const UInt_t  NHistQA(4);
static const UInt_t  NPtBins(48);
static const UInt_t  NChrgIds(16);
static const UInt_t  NPionIds(1);
static const UInt_t  NBadRuns(45);
static const UInt_t  NTrkMax(5000);
static const UInt_t  NTwrMax(5000);
static const UInt_t  NMatchMax(10);
static const UInt_t  NTrkSpecies(2);
static const Float_t MassPi0(0.140);



class StTrackEfficiencyCalculator {

  private:

    TFile   *_fIn;
    TFile   *_fOut;
    TTree   *_tPar;
    TTree   *_tDet;
    TTree   *_tJetPar;
    TTree   *_tJetDet;
    UInt_t  _fTrig;
    Bool_t  _isInBatchMode;
    Bool_t  _useParTrg;
    Bool_t  _useDetTrg;
    Bool_t  _useOnlyRecoilTracks;
    Bool_t  _useSpecificTrackSpecies;
    TString _sIn;
    TString _sOut;
    TString _sPar;
    TString _sDet;

    // jet parameters
    UInt_t   _nRepeat;
    UInt_t   _nRemove;
    Double_t _rJet;
    Double_t _aGhost;
    Double_t _pTjetMin;
    Double_t _etaJetMax;
    Double_t _etaBkgdMax;
    Double_t _etaGhostMax;

    // output particle event members
    Int_t    _pEventIndex;
    Int_t    _pNJets;
    Int_t    _pRunId;
    Double_t _pRefmult;
    Double_t _pPartonicPt;
    Double_t _pTSP;
    Double_t _pTrgEta;
    Double_t _pTrgPhi;
    Double_t _pTrgEt;
    Double_t _pRho;
    Double_t _pSigma;
    Double_t _pVz;
    // output particle jet members
    vector<Int_t>    _pJetIndex; 
    vector<Int_t>    _pJetNCons;
    vector<Double_t> _pJetPt;
    vector<Double_t> _pJetPtCorr;
    vector<Double_t> _pJetEta;
    vector<Double_t> _pJetPhi;
    vector<Double_t> _pJetE;
    vector<Double_t> _pJetArea;
    // output particle cst. members
    vector<vector<Double_t>> _pJetConsPt;
    vector<vector<Double_t>> _pJetConsEta;
    vector<vector<Double_t>> _pJetConsPhi;
    vector<vector<Double_t>> _pJetConsE;

    // output detector event members
    Int_t    _dEventIndex;
    Int_t    _dNJets;
    Int_t    _dRunId;
    Double_t _dRefmult;
    Double_t _dPartonicPt;
    Double_t _dTSP;
    Double_t _dTrgEta;
    Double_t _dTrgPhi;
    Double_t _dTrgEt;
    Double_t _dRho;
    Double_t _dSigma;
    Double_t _dVz;
    // output detector jet members
    vector<Int_t>    _dJetIndex; 
    vector<Int_t>    _dJetNCons;
    vector<Double_t> _dJetPt;
    vector<Double_t> _dJetPtCorr;
    vector<Double_t> _dJetEta;
    vector<Double_t> _dJetPhi;
    vector<Double_t> _dJetE;
    vector<Double_t> _dJetArea;
    // output detector cst. members
    vector<vector<Double_t>> _dJetConsPt;
    vector<vector<Double_t>> _dJetConsEta;
    vector<vector<Double_t>> _dJetConsPhi;
    vector<vector<Double_t>> _dJetConsE;

    // input particle event leaves
    Int_t    mcEventId;
    Int_t    mcRunId;
    Int_t    mcNumTrks;
    Double_t muVx;
    Double_t muVy;
    Double_t muVz;
    Double_t mcVx;
    Double_t mcVy;
    Double_t mcVz;
    // input particle track leaves
    vector<Int_t>    *mcIdTrk;
    vector<Int_t>    *mcIdGeant;
    vector<Int_t>    *mcIdVx;
    vector<Int_t>    *mcIdVxEnd;
    vector<Int_t>    *mcIntrVtx;
    vector<Bool_t>   *mcIsShower;
    vector<Double_t> *mcCharge;
    vector<Double_t> *mcRapidity;
    vector<Double_t> *mcEta;
    vector<Double_t> *mcPhi;
    vector<Double_t> *mcPx;
    vector<Double_t> *mcPy;
    vector<Double_t> *mcPz;
    vector<Double_t> *mcPt;
    vector<Double_t> *mcPtot;
    vector<Double_t> *mcEnergy;

    // input detector event leaves
    UInt_t   fUniqueID;
    UInt_t   fBits;
    Long64_t runNumber;
    Long64_t eventNumber;
    Int_t    trigID;
    Int_t    nGlobalTracks;
    Int_t    nPrimaryTracks;
    Int_t    refMult;
    Double_t vpdVz;
    Double_t xVertex;
    Double_t yVertex;
    Double_t zVertex;
    Double_t bbcZVertex;
    Double_t zdcCoincidenceRate;
    Double_t bbcCoincidenceRate;
    Double_t backgroundRate;
    Double_t bbcBlueBackgroundRate;
    Double_t bbcYellowBackgroundRate;
    Double_t refMultPos;
    Double_t refMultNeg;
    Double_t bTOFTrayMultiplicity;
    Int_t    nVerticies;
    Double_t MagF;
    Double_t VrtxRank;
    Float_t  Etsp;
    Int_t    ETwrdidT;
    Int_t    ETwradc11;
    Float_t  ETwreneT0;
    Float_t  ETwreT;
    Float_t  ETwrENET0;
    Float_t  ETwrphT;
    Float_t  ETwrPTower;
    Float_t  ETwrpidTower;
    Int_t    ETwrmoduleT;
    Float_t  EClustEneT0;
    Float_t  EClustetav1;
    Float_t  EClustphiv1;
    Float_t  EEstrpen01;
    Float_t  EEstrpen02;
    Float_t  EEstrpen03;
    Float_t  EEstrpen0;
    Float_t  EEstrpen1;
    Float_t  EEstrpen2;
    Float_t  EEstrpen3;
    Float_t  EEstrpen4;
    Float_t  EEstrpen5;
    Float_t  EEstrpen6;
    Float_t  EEstrpen7;
    Float_t  EEstrpen8;
    Float_t  EEstrpen9;
    Float_t  EEstrpen10;
    Float_t  EEstrpen11;
    Float_t  EEstrpen12;
    Float_t  EEstrpen13;
    Float_t  EEstrpen14;
    Float_t  EEstrpen15;
    Int_t    ETwrdidE;
    Float_t  EPstripenp01;
    Float_t  EPstripenp02;
    Float_t  EPstripenp03;
    Float_t  EPstripenp0;
    Float_t  EPstripenp1;
    Float_t  EPstripenp2;
    Float_t  EPstripenp3;
    Float_t  EPstripenp4;
    Float_t  EPstripenp5;
    Float_t  EPstripenp6;
    Float_t  EPstripenp7;
    Float_t  EPstripenp8;
    Float_t  EPstripenp9;
    Float_t  EPstripenp10;
    Float_t  EPstripenp11;
    Float_t  EPstripenp12;
    Float_t  EPstripenp13;
    Float_t  EPstripenp14;
    Float_t  EPstripenp15;
    Float_t  EclustEnnq1;
    Float_t  EclustEnnq20;
    Float_t  EclustEnnq19;
    Float_t  EclustEnpq1;
    Float_t  EclustEnpq20;
    Float_t  EclustEnpq19;
    Float_t  EclustEnpq21;
    // input detector track leaves
    Int_t    PrimaryTrackArray_;
    UInt_t   PrimaryTrackArray_fUniqueID[NTrkMax];
    UInt_t   PrimaryTrackArray_fBits[NTrkMax];
    Int_t    PrimaryTrackArray_nHitsFit[NTrkMax];
    Int_t    PrimaryTrackArray_nHitsPoss[NTrkMax];
    Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
    Int_t    PrimaryTrackArray_pdgId[NTrkMax];
    Int_t    PrimaryTrackArray_geantId[NTrkMax];
    Double_t PrimaryTrackArray_pZ[NTrkMax];
    Double_t PrimaryTrackArray_pX[NTrkMax];
    Double_t PrimaryTrackArray_pY[NTrkMax];
    Double_t PrimaryTrackArray_pT[NTrkMax];
    Double_t PrimaryTrackArray_dEdx[NTrkMax];
    Double_t PrimaryTrackArray_charge[NTrkMax];
    Double_t PrimaryTrackArray_tofBeta[NTrkMax];
    Double_t PrimaryTrackArray_eta[NTrkMax];
    Double_t PrimaryTrackArray_phi[NTrkMax];
    Double_t PrimaryTrackArray_nSigElectron[NTrkMax];
    Double_t PrimaryTrackArray_nSigPion[NTrkMax];
    Double_t PrimaryTrackArray_nSigKaon[NTrkMax];
    Double_t PrimaryTrackArray_nSigProton[NTrkMax];
    Double_t PrimaryTrackArray_dcag[NTrkMax];
    Double_t PrimaryTrackArray_nHits[NTrkMax];
    Double_t PrimaryTrackArray_dEdxHits[NTrkMax];
    Double_t PrimaryTrackArray_firstZPoint[NTrkMax];
    Double_t PrimaryTrackArray_lastZPoint[NTrkMax];
    Double_t PrimaryTrackArray_tofSigElectron[NTrkMax];
    Double_t PrimaryTrackArray_tofSigPion[NTrkMax];
    Double_t PrimaryTrackArray_tofSigKaon[NTrkMax];
    Double_t PrimaryTrackArray_tofSigProton[NTrkMax];
    Double_t PrimaryTrackArray_timeOfflight[NTrkMax];
    Double_t PrimaryTrackArray_pathLength[NTrkMax];
    Int_t    PrimaryTrackArray_trkIndex[NTrkMax];
    // input detector tower leaves
    Int_t    TowerArray_;
    UInt_t   TowerArray_fUniqueID[NTwrMax];
    UInt_t   TowerArray_fBits[NTwrMax];
    Int_t    TowerArray_TwrId[NTwrMax];
    Float_t  TowerArray_TwrEng[NTwrMax];
    Float_t  TowerArray_TwrEta[NTwrMax];
    Float_t  TowerArray_TwrPhi[NTwrMax];
    Float_t  TowerArray_TwrADC[NTwrMax];
    Float_t  TowerArray_TwrPed[NTwrMax];
    Float_t  TowerArray_TwrRMS[NTwrMax];
    Int_t    TowerArray_TwrMatchIdnex[NTwrMax];
    Int_t    TowerArray_NoOfmatchedTrk[NTwrMax];
    Float_t  TowerArray_TwrMatchP[NTwrMax];
    Float_t  TowerArray_TwrPx[NTwrMax];
    Float_t  TowerArray_TwrPy[NTwrMax];
    Float_t  TowerArray_TwrPz[NTwrMax];
    Int_t    TowerArray_fNAssocTracks[NTwrMax];
    Int_t    TowerArray_fMatchedTracksArray_[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_P[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_nSigPi[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_nSigK[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_nSigP[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_nSigE[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_dcag[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_eta[NTwrMax][NMatchMax];
    Float_t  TowerArray_fMatchedTracksArray_pT[NTwrMax][NMatchMax];
    Int_t    TowerArray_fMatchedTracksArray_nFit[NTwrMax][NMatchMax];
    Int_t    TowerArray_fMatchedTracksArray_nPos[NTwrMax][NMatchMax];

    // input particle branches
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

    // input detector branches
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
    TBranch *_bPrimaryTrackArray_pdgId;
    TBranch *_bPrimaryTrackArray_geantId;
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
    TBranch *_bTowerArray_fMatchedTracksArray_nSigPi;
    TBranch *_bTowerArray_fMatchedTracksArray_nSigK;
    TBranch *_bTowerArray_fMatchedTracksArray_nSigP;
    TBranch *_bTowerArray_fMatchedTracksArray_nSigE;
    TBranch *_bTowerArray_fMatchedTracksArray_dcag;
    TBranch *_bTowerArray_fMatchedTracksArray_eta;
    TBranch *_bTowerArray_fMatchedTracksArray_pT;
    TBranch *_bTowerArray_fMatchedTracksArray_nFit;
    TBranch *_bTowerArray_fMatchedTracksArray_nPos;

  public:

    // ctor, dtor
    StTrackEfficiencyCalculator(const Bool_t batch=false);
    virtual ~StTrackEfficiencyCalculator();

    // public methods
    void Init(const TString sInput, const TString sOutput, const TString sParTree, const TString sDetTree);
    void Make(const UInt_t trigger, const Bool_t parTrg, const Bool_t detTrg, const Bool_t recoilTrks, const Bool_t trkSpecies);
    void Finish();
    void SetJetParameters(const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t pTjetMin, const Double_t aGhost, const Double_t hJetMax, const Double_t hBkgdMax, const Double_t hGhostMax);

  ClassDef(StTrackEfficiencyCalculator, 1);

};


#endif
#ifdef StTrackEfficiencyCalculator_cxx

StTrackEfficiencyCalculator::StTrackEfficiencyCalculator(const Bool_t batch) {

  _isInBatchMode = batch;
  _pJetIndex.clear();
  _pJetNCons.clear();
  _pJetPt.clear();
  _pJetPtCorr.clear();
  _pJetEta.clear();
  _pJetPhi.clear();
  _pJetE.clear();
  _pJetArea.clear();  
  _pJetConsPt.clear();
  _pJetConsEta.clear();
  _pJetConsPhi.clear();
  _pJetConsE.clear();
  _dJetIndex.clear();
  _dJetNCons.clear();
  _dJetPt.clear();
  _dJetPtCorr.clear();
  _dJetEta.clear();
  _dJetPhi.clear();
  _dJetE.clear();
  _dJetArea.clear();  
  _dJetConsPt.clear();
  _dJetConsEta.clear();
  _dJetConsPhi.clear();
  _dJetConsE.clear();

}  // end ctor



StTrackEfficiencyCalculator::~StTrackEfficiencyCalculator() {

 if (!_tPar && !_tDet) return;
 delete _tPar -> GetCurrentFile();
 delete _tDet -> GetCurrentFile();

}  // end dtor

#endif

// End ------------------------------------------------------------------------
