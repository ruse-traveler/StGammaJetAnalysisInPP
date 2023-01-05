// 'StTriggerEnergyScaleCalculator.h'
// Derek Anderson
// 06.15.2020
//
// This class produces trees of jets using
// generated particles and reconstructed
// tracks triggered on generated pi0s
// matched to reconstructed pi0s.  Also
// calculates the trigger energy scale
// and resolution.
//
// NOTE: assumes particle-level and
//   detector-level input trees were
//   filled in same order.

#ifndef StTriggerEnergyScaleCalculator_h
#define StTriggerEnergyScaleCalculator_h

#include <vector>
#include <iostream>
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TDirectory.h"

using namespace std;

// global constants
static const UInt_t  NEvtMax(100000);
static const UInt_t  NTrkMax(5000);
static const UInt_t  NTwrMax(5000);
static const UInt_t  NMatchMax(10);
static const UInt_t  NHotTwrs(302);
static const UInt_t  NBadRuns(45);
static const UInt_t  NTrgBin(3);
static const UInt_t  NTrgTsp(2);
static const UInt_t  NTrgIds(1);
static const UInt_t  NTrgs(2);
static const Float_t MassPi0(0.14);



class StTriggerEnergyScaleCalculator {

  private:

    // private methods
    void InitTrees();
    void InitHists();
    void SaveTreesAndHists();
    void DoDetFirstMatching();
    void DoParFirstMatching();

    // calculation parameters
    UInt_t _trgType;
    UInt_t _matchOrder;
    Bool_t _useDetTrigger;
    Bool_t _useParTrigger;
    Bool_t _requireMatch;
    Bool_t _excludeNeighbors;
    Bool_t _avoidDoubleCounting;

    // io and root parameters
    TF1     *_fPi0;
    TF1     *_fGam;
    TF1     *_fWeight;
    TFile   *_fOut;
    TFile   *_fPar;
    TFile   *_fDet;
    TTree   *_tPar;
    TTree   *_tMat;
    TTree   *_tDet;
    TTree   *_tJetPar;
    TTree   *_tJetDet; 
    Bool_t  _isInBatchMode;
    TString _sOut;
    TString _sPar;
    TString _sDet;
    TString _sParTree;
    TString _sMatTree;
    TString _sDetTree;

    // event parameters
    Double_t _rVtxMax;
    Double_t _zVtxMax;

    // detector trigger parameters
    Int_t    _adcMax;
    Double_t _eStrMin;
    Double_t _pProjMax;
    Double_t _hTrgDetMax;
    Double_t _eTtrgDetMin;
    Double_t _eTtrgDetMax;
    Double_t _eTbinMin[NTrgBin];
    Double_t _eTbinMax[NTrgBin];
    Double_t _tspPi0[NTrgTsp];
    Double_t _tspGam[NTrgTsp];
    Double_t _tspUse[NTrgTsp];

    // particle trigger parameters
    Double_t _cTrgPar;
    Double_t _dRtrgMaxPi0;
    Double_t _dRtrgMaxGam;
    Double_t _dRtrgMaxUse;
    Double_t _hTrgParMax;
    Double_t _eTtrgMatMin;
    Double_t _eTtrgMatMax;
    Double_t _eTtrgParMin;
    Double_t _eTtrgParMax;
    Double_t _dMinSeparationPi0;
    Double_t _dMinSeparationGam;
    Double_t _dMinSeparationUse;
    Double_t _aFitParms[NTrgs];
    Double_t _bFitParms[NTrgs];
    UInt_t   _idTrgPi0[NTrgIds];
    UInt_t   _idTrgGam[NTrgIds];
    UInt_t   _idTrgUse[NTrgIds];

    // track parameters
    UInt_t   _nFitTrkMin;
    Double_t _rFitTrkMin;
    Double_t _dcaTrkMax;
    Double_t _hTrkMax;
    Double_t _pTtrkParMin;
    Double_t _pTtrkParMax;
    Double_t _pTtrkDetMin;
    Double_t _pTtrkDetMax;

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
    Double_t _pTrgEtPar;
    Double_t _pTrgEtDet;
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
    vector<vector<Double_t> > _pJetConsPt;
    vector<vector<Double_t> > _pJetConsEta;
    vector<vector<Double_t> > _pJetConsPhi;
    vector<vector<Double_t> > _pJetConsE;

    // output detector event members
    Int_t    _dEventIndex;
    Int_t    _dNJets;
    Int_t    _dRunId;
    Double_t _dRefmult;
    Double_t _dPartonicPt;
    Double_t _dTSP;
    Double_t _dTrgEta;
    Double_t _dTrgPhi;
    Double_t _dTrgEtPar;
    Double_t _dTrgEtDet;
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
    vector<vector<Double_t> > _dJetConsPt;
    vector<vector<Double_t> > _dJetConsEta;
    vector<vector<Double_t> > _dJetConsPhi;
    vector<vector<Double_t> > _dJetConsE;

    // input particle event leaves
    Int_t    _pMcEventId;
    Int_t    _pMcRunId;
    Int_t    _pMcNumTrks;
    Double_t _pMuVx;
    Double_t _pMuVy;
    Double_t _pMuVz;
    Double_t _pMcVx;
    Double_t _pMcVy;
    Double_t _pMcVz;
    // input particle track leaves
    vector<Int_t>    *_pMcIdTrk;
    vector<Int_t>    *_pMcIdGeant;
    vector<Int_t>    *_pMcIdVx;
    vector<Int_t>    *_pMcIdVxEnd;
    vector<Int_t>    *_pMcIntrVtx;
    vector<Bool_t>   *_pMcIsShower;
    vector<Double_t> *_pMcCharge;
    vector<Double_t> *_pMcRapidity;
    vector<Double_t> *_pMcEta;
    vector<Double_t> *_pMcPhi;
    vector<Double_t> *_pMcPx;
    vector<Double_t> *_pMcPy;
    vector<Double_t> *_pMcPz;
    vector<Double_t> *_pMcPt;
    vector<Double_t> *_pMcPtot;
    vector<Double_t> *_pMcEnergy;

    // input matched particle event leaves
    Int_t    _mMcEventId;
    Int_t    _mMcRunId;
    Int_t    _mMcNumTrks;
    Double_t _mMuVx;
    Double_t _mMuVy;
    Double_t _mMuVz;
    Double_t _mMcVx;
    Double_t _mMcVy;
    Double_t _mMcVz;
    // input matched particle track leaves
    vector<Int_t>    *_mMcIdTrk;
    vector<Int_t>    *_mMcIdGeant;
    vector<Int_t>    *_mMcIdVx;
    vector<Int_t>    *_mMcIdVxEnd;
    vector<Int_t>    *_mMcIntrVtx;
    vector<Bool_t>   *_mMcIsShower;
    vector<Double_t> *_mMcCharge;
    vector<Double_t> *_mMcRapidity;
    vector<Double_t> *_mMcEta;
    vector<Double_t> *_mMcPhi;
    vector<Double_t> *_mMcPx;
    vector<Double_t> *_mMcPy;
    vector<Double_t> *_mMcPz;
    vector<Double_t> *_mMcPt;
    vector<Double_t> *_mMcPtot;
    vector<Double_t> *_mMcEnergy;

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
    Int_t    FlagEvent_TrgTrkMisMtch;
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
    Float_t  EtspAlt;
    Float_t  EClustSecondTwrEnergy;
    UInt_t   EClustNumTwrIncluded;
    Int_t    EClustSecondTwrIndex;
    Int_t    EEstrpModuleE;
    Int_t    EPstripModuleP;
    Int_t    EClustSecondTwrModule;
    // input detector track leaves
    Int_t    PrimaryTrackArray_;
    UInt_t   PrimaryTrackArray_fUniqueID[NTrkMax];
    UInt_t   PrimaryTrackArray_fBits[NTrkMax];
    Double_t PrimaryTrackArray_nHitsFit[NTrkMax];
    Double_t PrimaryTrackArray_nHitsPoss[NTrkMax];
    Int_t    PrimaryTrackArray_trackFlag[NTrkMax];
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

    // input particle branches
    TBranch *_bParMcEventId;
    TBranch *_bParMcRunId;
    TBranch *_bParMcNumTrks;
    TBranch *_bParMuVx;
    TBranch *_bParMuVy;
    TBranch *_bParMuVz;
    TBranch *_bParMcVx;
    TBranch *_bParMcVy;
    TBranch *_bParMcVz;
    TBranch *_bParMcIdTrk;
    TBranch *_bParMcIdGeant;
    TBranch *_bParMcIdVx;
    TBranch *_bParMcIdVxEnd;
    TBranch *_bParMcIntrVtx;
    TBranch *_bParMcIsShower;
    TBranch *_bParMcCharge;
    TBranch *_bParMcRapidity;
    TBranch *_bParMcEta;
    TBranch *_bParMcPhi;
    TBranch *_bParMcPx;
    TBranch *_bParMcPy;
    TBranch *_bParMcPz;
    TBranch *_bParMcPt;
    TBranch *_bParMcPtot;
    TBranch *_bParMcEnergy;

    // input matched particle branches
    TBranch *_bMatMcEventId;
    TBranch *_bMatMcRunId;
    TBranch *_bMatMcNumTrks;
    TBranch *_bMatMuVx;
    TBranch *_bMatMuVy;
    TBranch *_bMatMuVz;
    TBranch *_bMatMcVx;
    TBranch *_bMatMcVy;
    TBranch *_bMatMcVz;
    TBranch *_bMatMcIdTrk;
    TBranch *_bMatMcIdGeant;
    TBranch *_bMatMcIdVx;
    TBranch *_bMatMcIdVxEnd;
    TBranch *_bMatMcIntrVtx;
    TBranch *_bMatMcIsShower;
    TBranch *_bMatMcCharge;
    TBranch *_bMatMcRapidity;
    TBranch *_bMatMcEta;
    TBranch *_bMatMcPhi;
    TBranch *_bMatMcPx;
    TBranch *_bMatMcPy;
    TBranch *_bMatMcPz;
    TBranch *_bMatMcPt;
    TBranch *_bMatMcPtot;
    TBranch *_bMatMcEnergy;

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
    TBranch *_bEventList_EtspAlt;
    TBranch *_bEventList_EClustNumTwrIncluded;
    TBranch *_bEventList_EClustSecondTwrIndex;
    TBranch *_bEventList_EClustSecondTwrEnergy;
    TBranch *_bEventList_EEstrpModuleE;
    TBranch *_bEventList_EPstripModuleP;
    TBranch *_bEventList_EClustSecondTwrModule;
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

    // calculation histogram members
    TH1D *_hNumMatch;
    TH1D *_hNumVsCut;
    TH1D *_hRatioVsCut;
    TH1D *_hPhiEff;
    TH1D *_hEtaEff;
    TH1D *_hEtEff;
    TH1D *_hPhiForEff[2];
    TH1D *_hEtaForEff[2];
    TH1D *_hEtForEff[2];
    TH1D *_hEtVsCut[7];
    TH1D *_hTspTrg[4];
    TH1D *_hPhiTrg[10];
    TH1D *_hEtaTrg[10];
    TH1D *_hEtTrg[10];
    TH1D *_hDrTrg[4];
    TH1D *_hQtTrg[4];
    TH1D *_hDqTrg[4];
    TH1D *_hEtDet[NTrgBin];
    TH1D *_hEtDetNW[NTrgBin];
    TH1D *_hTspMatch[NTrgBin];
    TH1D *_hTspMatchNW[NTrgBin];
    TH1D *_hEtMatch[NTrgBin];
    TH1D *_hEtMatchNW[NTrgBin];
    TH1D *_hDrMatch[NTrgBin];
    TH1D *_hDrMatchNW[NTrgBin];
    TH1D *_hQtMatch[NTrgBin];
    TH1D *_hQtMatchNW[NTrgBin];
    TH1D *_hDqMatch[NTrgBin];
    TH1D *_hDqMatchNW[NTrgBin];
    TH2D *_hQtParVsEtPar;
    TH2D *_hQtParVsEtParNW;
    TH2D *_hQtParVsEtDet;
    TH2D *_hQtParVsEtDetNW;
    TH2D *_hDrParVsEtPar;
    TH2D *_hDrParVsEtParNW;
    TH2D *_hDrParVsEtDet;
    TH2D *_hDrParVsEtDetNW;
    TH2D *_hDrParVsQtPar;
    TH2D *_hDrParVsQtParNW;
    TH2D *_hPhiMatchVsDet;
    TH2D *_hPhiMatchVsDetNW;
    TH2D *_hEtaMatchVsDet;
    TH2D *_hEtaMatchVsDetNW;
    TH2D *_hEtMatchVsDet;
    TH2D *_hEtMatchVsDetNW;
    TH2D *_hEtMatchVsDetTrg[NTrgBin];
    TH2D *_hEtMatchVsDetTrgNW[NTrgBin];
    TH2D *_hDrMatchVsEtDet;
    TH2D *_hDrMatchVsEtDetNW;
    TH2D *_hQtMatchVsEtMatch;
    TH2D *_hQtMatchVsEtMatchNW;
    TH2D *_hQtMatchVsEtDet;
    TH2D *_hQtMatchVsEtDetNW;
    TH2D *_hQtInvMatchVsEtDet;
    TH2D *_hQtInvMatchVsEtDetNW;
    TH2D *_hDqMatchVsEtDet;
    TH2D *_hDqMatchVsEtDetNW;
    TH3D *_hQtEtMatchVsDet;

    // check histogram members
    TH1D *_hDidT_weird;
    TH1D *_hDidT_normal;
    TH1D *_hDidE_weird;
    TH1D *_hDidE_normal;
    TH1D *_hTrigID_weird;
    TH1D *_hTrigID_normal;
    TH1D *_hModuleT_weird;
    TH1D *_hModuleT_normal;
    TH1D *_hEeneT0_weird;
    TH1D *_hEeneT0_normal;
    TH1D *_hEENET0_weird;
    TH1D *_hEENET0_normal;
    TH1D *_hEEneT0_weird;
    TH1D *_hEEneT0_normal;
    TH1D *_hEtReco_weird;
    TH1D *_hEtReco_normal;
    TH1D *_hEtSim_weird;
    TH1D *_hEtSim_normal;
    TH1D *_hEenToEEnRatio_weird;
    TH1D *_hEenToEEnRatio_normal;
    TH1D *_hNumTwrInClust_weird;
    TH1D *_hNumTwrInClust_normal;
    TH1D *_hEtaModule;
    TH1D *_hPhiModule;
    TH1D *_hTwrModule;
    TH1D *_hSecModule;
    TH1D *_hDiffPhiMod;
    TH1D *_hDiffTwrMod;
    TH1D *_hDiffSecMod;
    TH1D *_hSecondTwrID_weird;
    TH1D *_hSecondTwrID_normal;
    TH1D *_hSecondTwrEnergy_weird;
    TH1D *_hSecondTwrEnergy_normal;
    TH2D *_hSecondVsLeadTwrEnergy_weird;
    TH2D *_hSecondVsLeadTwrEnergy_normal;

    // profile members
    TProfile *_pPhiMatchVsDet;
    TProfile *_pPhiMatchVsDetNW;
    TProfile *_pEtaMatchVsDet;
    TProfile *_pEtaMatchVsDetNW;
    TProfile *_pEtMatchVsDet;
    TProfile *_pEtMatchVsDetNW;
    TProfile *_pEtMatchVsDetTrg[NTrgBin];
    TProfile *_pEtMatchVsDetTrgNW[NTrgBin];
    TProfile *_pDrMatchVsEtDet;
    TProfile *_pDrMatchVsEtDetNW;
    TProfile *_pQtMatchVsEtMatch;
    TProfile *_pQtMatchVsEtMatchNW;
    TProfile *_pQtMatchVsEtDet;
    TProfile *_pQtMatchVsEtDetNW;
    TProfile *_pQtInvMatchVsEtDet;
    TProfile *_pQtInvMatchVsEtDetNW;
    TProfile *_pDqMatchVsEtDet;
    TProfile *_pDqMatchVsEtDetNW;

    // bad run and hot tower lists
    static const UInt_t _badRunList[NBadRuns];
    static const UInt_t _hotTwrList[NHotTwrs];

  public:

    // ctor, dtor
    StTriggerEnergyScaleCalculator(const Bool_t batch=false);
    virtual ~StTriggerEnergyScaleCalculator();

    // general public methods
    void Init(const TString sOutput, const TString sParticle, const TString sDetector, const TString sParTree, const TString sMatTree, const TString sDetTree);
    void Make(const UInt_t trgType, const UInt_t matchOrder, const Bool_t useDetTrigger, const Bool_t useParTrigger, const Bool_t requireMatch, const Bool_t excludeNeighbors, const Bool_t avoidDoubleCounting);
    void Finish();

    // parameter setters
    void SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax);
    void SetDetectorParameters(const Int_t adcMax, const Double_t eStrMin, const Double_t pProjMax, const Double_t hTrgDetMax, const Double_t eTtrgDetMin, const Double_t eTtrgDetMax, const Double_t *eTbinMin, const Double_t *eTbinMax, const Double_t *tspPi0, const Double_t *tspGam);
    void SetParticleParameters(const Double_t cTrgPar, const Double_t dRtrgMaxPi0, const Double_t dRtrgMaxGam, const Double_t hTrgParMax, const Double_t eTtrgMatMin, const Double_t eTtrgMatMax, const Double_t eTtrgParMin, const Double_t eTtrgParMax, const Double_t dMinSeparationPi0, const Double_t dMinSeparationGam, const Double_t *aFitParms, const Double_t *bFitParms, const UInt_t *idTrgPi0, const UInt_t *idTrgGam);
    void SetTrackParameters(const UInt_t nFitTrkMin, const Double_t rFitTrkMin, const Double_t dcaTrkMax, const Double_t hTrkMax, const Double_t pTtrkParMin, const Double_t pTtrkParMax, const Double_t pTtrkDetMin, const Double_t pTtrkDetMax);
    void SetJetParameters(const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t pTjetMin, const Double_t aGhost, const Double_t hJetMax, const Double_t hBkgdMax, const Double_t hGhostMax);

  ClassDef(StTriggerEnergyScaleCalculator, 1);

};



#endif
#ifdef StTriggerEnergyScaleCalculator_cxx

StTriggerEnergyScaleCalculator::StTriggerEnergyScaleCalculator(const Bool_t batch) {

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



StTriggerEnergyScaleCalculator::~StTriggerEnergyScaleCalculator() {

 if (!_tPar && !_tMat && !_tDet) return;
 delete _tPar -> GetCurrentFile();
 delete _tMat -> GetCurrentFile();
 delete _tDet -> GetCurrentFile();

}  // end dtor



// initialize bad run list
const UInt_t StTriggerEnergyScaleCalculator::_badRunList[NBadRuns] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};



// initialize hot tower list
const UInt_t StTriggerEnergyScaleCalculator::_hotTwrList[NHotTwrs] = {34, 106, 113, 160, 266, 267, 275, 280, 282, 286, 287, 293, 410, 504, 533, 541, 555, 561, 562, 594, 615, 616, 629, 633, 637, 638, 647, 650, 653, 657, 671, 673, 743, 789, 790, 791, 792, 806, 809, 810, 811, 812, 813, 814, 821, 822, 823, 824, 829, 830, 831, 832, 837, 841, 842, 843, 844, 846, 849, 850, 851, 852, 857, 875, 897, 899, 903, 939, 953, 954, 956, 993, 1026, 1046, 1048, 1080, 1081, 1100, 1125, 1130, 1132, 1180, 1197, 1198, 1199, 1200, 1207, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1312, 1348, 1353, 1354, 1388, 1407, 1409, 1434, 1448, 1537, 1567, 1574, 1597, 1612, 1654, 1668, 1713, 1762, 1765, 1766, 1877, 1878, 1984, 2032, 2043, 2054, 2073, 2077, 2092, 2093, 2097, 2107, 2162, 2168, 2214, 2305, 2392, 2409, 2415, 2439, 2459, 2589, 2590, 2633, 2652, 2749, 2834, 2961, 2969, 3005, 3017, 3070, 3071, 3186, 3220, 3289, 3360, 3493, 3494, 3495, 3508, 3588, 3604, 3611, 3668, 3678, 3679, 3690, 3692, 3732, 3738, 3838, 3840, 3927, 3945, 4005, 4006, 4013, 4018, 4019, 4053, 4059, 4124, 4331, 4355, 4357, 4458, 4464, 4500, 4677, 4678, 4684, 4768, 360, 493, 779, 1284, 1306, 1337, 1438, 1709, 2027, 2445, 3407, 3720, 4217, 4288, 95, 96, 296, 316, 443, 479, 555, 562, 637, 671, 709, 740, 743, 796, 857, 897, 899, 915, 953, 1130, 1132, 1294, 1318, 1337, 1348, 1359, 1378, 1427, 1429, 1440, 1537, 1563, 1574, 1709, 1763, 1773, 1819, 1854, 1874, 1936, 1938, 2018, 2043, 2098, 2099, 2256, 2259, 2294, 2514, 2520, 2552, 2589, 2598, 2680, 2706, 2799, 2880, 2897, 2917, 2969, 3020, 3028, 3310, 3319, 3375, 3399, 3504, 3539, 3541, 3679, 3690, 3692, 3718, 3719, 3720, 3738, 3806, 3838, 3840, 3928, 4013, 4017, 4038, 4053, 4057, 4058, 4079, 4097, 4099};

#endif

// End ------------------------------------------------------------------------
