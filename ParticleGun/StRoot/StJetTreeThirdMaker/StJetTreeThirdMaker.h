//________________________________________________________________
//This is the Only clustering algorithim for Gamma-hadron correlation
//analysis. This class uses BEMC, BSMD, and TPC information to
//find HighTower info, pi0/photon discrimination info and hadrons
//Tracks projection info to project on to BEMC.
//
//
// STAR published paper: Phys. Rev. C 82, 034909 (2010)
// based on this clustering algorithm.
// Contact: A. Hamed and Saskia Mioduszewski
//
// Date 07/2014: TPC tracks from StMuTrack -  Nihar r. sahoo
//             : globalDca has been implemented
//
//Date  12/2015: "StJetTreeThirdMaker" inherited from "StThirdMaker"
//              It makes "JetTree" for Full Jet Reconstruction
//              for Gamma+Jet analysis.      - Nihar r. Sahoo
//Date  06/2016: all events, tracks and tower attributes are added and 
//                                           - Nihar

//Date 04/20/2017: ____StJetTreeThirdMaker must be submiited for "1 mudst file per process in .xml"                   !!!!!ATTENTION
//                becasue it creates problem for StEvent and StMuDst mismatch in EventId hence results                !!!!!ATTENTION
//                mismatch in trigger(BEMC) and tracks (TPC) to produce "femtodst"                   -Nihar           !!!!!ATTENTION
//Date 07/12/2017: added E/p calculation, courtesy of Nihar.
//                 -Derek
//Date 07/21/2017: added additional branches to calculate E/p for
//                 electrons and hadrons.
//                 -Derek
//Date 12/23/2020: added additional output branches for cluster information
//                 -Derek
//Date 01/10/2021: added additional output brach to store energy of 2nd tower in cluster
//                 -Derek
//Date 01/12/2021: added additional output branch to store module #s of eta, phi strips and 2nd tower
//                 -Derek
//____________________________________________________________


#ifndef STAR_StJetTreeThirdMaker
#define STAR_StJetTreeThirdMaker

#ifndef StMaker_H
#include "StMaker.h"
#include <string>
#include <vector>


#include "StBichsel/Bichsel.h"
#include "StBichsel/dEdxParameterization.h"


#endif


class StEvent;
class StMcEvent;
class TString;
class Bichsel;
class TNtuple;
class TFile;
class TH1F;
class TH2F;
class StPrimaryVertex;

class StTriggerData;
class StTriggerData2007;

class StEmcCluster;
class StEmcCollection;
class StEmcDetector;
class StEmcModule;
class StEmcModuleHitCollection;
class StEmcRawHit;
class StEmcFilter;
//class StThreeVectorF;
class TRandom;
class St_db_Maker;
class StEmcDecoder;
class emcPed_st;
class StMcVertex;
class StMuDstMaker;
class StMuDst;
class StMuEvent;
class StMuTrack;
class StEmcPosition;
class StBemcTables; //v3.14

//class vector<float>;
//_______Class for Tree

class GammaJetEvent;
class GammaJetTrack;
class GammaJetTower;
class GammaJetTowerUtil;
class TClonesArray;
class  TLorentzVector;

#define TowerHVchangeMax 400

const  Int_t max_pTracks = 10000;
const  Int_t max_gTracks = 10000;



class StJetTreeThirdMaker : public StMaker 
{
 private: //BetheBloch      	mBB;
  TRandom*     mRandom;  
  // Private method declaration if any
  float             mEnergyE[150];
  float             mEnergyP[10][15];
  float             mEnergyE2[150];
  float             mEnergyP2[10][15];
 
  float             mEnergyTT[20][2];
  float             mEnergyTT1[20][2];
  float             mEnergyTT2[20][2];
  float             energy[6];
  
  float             mEnergyTT3[20][2];
  float             mEnergyTT14[20][2];
  float             mEnergyTT25[20][2];
  float             energy2[6];
  float ParentTrackP[4800]; float ParticleId[4800];
  float   ParentTrackPOSe[4800];
  float ParentTrackPOSP[4800];
  float             mEnergyTR[20][2];
  Float_t mStatusT[4800];
  float mHighest[8];
  float ReactionPlane[43];
  
  Float_t  bbcPedstal[2][16];
  Double_t bbcGainFac[2][16]; 
  Int_t HVTowerid[TowerHVchangeMax];
  St_db_Maker* mDbMaker;
  Int_t associateTracksWithTowerHits();Int_t associateTracksWithTowerHitsAll();
  Int_t associateTracksWithSmdeHits();
  Int_t associateTracksWithSmdpHits();
  Int_t associateTracksWithTowerHitsP();
  Int_t associateTracksWithSmdeHitsP();
  Int_t associateTracksWithSmdpHitsP();
  Int_t associateHadronsWithBEMC1(); 
  Int_t associateHadronsWithBEMC2(); 
  Int_t associateHadronsWithBEMC3();
  Int_t TpcTracks();  Int_t CalculateReactionPlane();
  void         readBadRunList();
  // Bichsel* m_Bichsel;
  Bool_t     mHt1;
  Bool_t     mHt2;
  Bool_t     mMinb;
  Bool_t     mCentral;Bool_t mJPsi;
  Bool_t mJPsi1;
  
  Bool_t    mL2gamma1;
  Bool_t    mL2gamma2;
  Bool_t mL2gamma1_2010;
  Bool_t mL2gamma2_2010;
  Bool_t mL2gamma0_2011;
  Bool_t mL2gamma1_2011;
  Bool_t mL2gamma2_2011;
  Bool_t    mbht2mb1;
  Bool_t     mbht2mb2;
  Bool_t     mbht2mb3;
  Bool_t     mbht2mb4;
  Bool_t     mbht2mb5;Bool_t     mbht2mb6;
  Bool_t    mHt0; Bool_t    mht4;
  Bool_t    mht4fast;Bool_t    mupsilondAu;

  Bool_t    L2BGamma1_ppRun9, L2EGamma1_ppRun9, L2EGamma2_ppRun9;


  void     getTrgType();        
  StEmcFilter*    mEmcFilter;
  StEvent           *mStEvent; 
  StMcEvent           *mc_event;
  StEmcDecoder* mEmcDecoder;
  StMuDstMaker* mMuDstMaker;
  StMuDst *mStMuDst;
  StMuDst *mDst; //NRS: new pointer to StMuDst 
  StMuEvent *muEvent;
  
  Float_t mCalibPedValues[4800];	//!
  Float_t bbc_phi[16];
  Float_t mCalibPedRMSValues[4800];       
  void TowerHVcahngeList();
  
  float BBC_GetPhi(int iTile);

 protected:
  // Protected method if any
  
 public:
  void          setPrint(Bool_t);
  /* StJetTreeThirdMaker(const char *name, StMuDstMaker* uDstMaker);*/
  StJetTreeThirdMaker(const char *name="ThirdMaker", char* dataType = "");
  virtual       ~StJetTreeThirdMaker();
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  void     Clear(Option_t *option="");
  /*Int_t associateTracksWithSmdeHits();*/
  void SetFileName( char* name){outFile=name;}; 
  Int_t ResizeHit(); 
  void  setDbMaker(St_db_Maker*);
  Bichsel* m_dEdxParameterization;
  StEmcFilter* getEmcFilter() { return mEmcFilter; };

  Bool_t Check_hot_Tower(Int_t TwrId);
  Bool_t Check_hot_EtaStrip(Int_t EtaStpId);


  const StMuTrack* ptrack;  
  StEmcPosition *mPosition;

  StBemcTables    *mTables;


  
  ///////////////////////

  TH1F* mHistBg;
  TH1F* hya;
  int  mEventCounter;
  TFile* File_output;
  char* outFile;

  // [Derek, 07.12.2017]
  TH1D* hEbp_elec;
  TH1D* hEbp_hadr;
  // [Derek, 07.21.2017]
  TH1D *hPt_elec;
  TH1D *hPt_hadr;
  // [Derek, 09.18.2020]
  TH1D *hCutDetector1_before;
  TH1D *hCutDetector1_after;
  TH1D *hCutModuleHB_before;
  TH1D *hCutModuleHB_after;
  TH1D *hCutStatusAll_before;
  TH1D *hCutStatusAll_after;
  TH1D *hCutEnergyHB_before;
  TH1D *hCutEnergyHB_after;
  // [Derek, 11.16.2020]
  TH1D *hNumPi0;
  TH1D *hNumDecay;

  // for efficiency [Derek, 12.02.2020]
  TH1D   *hNumBeforeVsCut;
  TH1D   *hNumAfterVsCut;
  TH1D   *hNumLeft[13][2];
  UInt_t  nLeftBeforeCut[13];
  UInt_t  nLeftAfterCut[13];

  // [Derek, 12.10.2020]
  TH1D *hPtPi0;
  TH1D *hPtTwr;
  TH1D *hDidT;
  TH1D *hDidE;

  ///________
  
  TBranch *eventBranch;
  
  Float_t *EvValues;
  
  TTree *outTree;
  GammaJetEvent *event;

  // for comparing particle to detector level [Derek, 05.20.2020]
  TTree *_tMcTracks;

  // for matching calculation [Derek, 11.05.2020]
  TTree *_tMcMatch;

  // mc event leaves
  Int_t    _McEvtId;
  Int_t    _McRunId;
  Int_t    _McNtrk;
  Double_t _MuVx;
  Double_t _MuVy;
  Double_t _MuVz;
  Double_t _McVx;
  Double_t _McVy;
  Double_t _McVz;
  // mc track leaves
  vector<Int_t>    _McIdTrk;
  vector<Int_t>    _McIdGnt;
  vector<Int_t>    _McIdVx;
  vector<Int_t>    _McIdVxEnd;
  vector<Int_t>    _McIntrVtx;
  vector<Bool_t>   _McIsShower;
  vector<Double_t> _McChrg;
  vector<Double_t> _McRap;
  vector<Double_t> _McEta;
  vector<Double_t> _McPhi;
  vector<Double_t> _McPx;
  vector<Double_t> _McPy;
  vector<Double_t> _McPz;
  vector<Double_t> _McPt;
  vector<Double_t> _McPtot;
  vector<Double_t> _McEne;

  
  
  //_____________for matching tracks
  Int_t *pTrMatchArr;//!
  Float_t *pTriMatchEtaArr;//!
  Float_t *pTrMatchPhiArr;//!
  Int_t *pTrIndexArray;//!
  Float_t *pTrEtaArray;//!
  Float_t *pTrPhiArray;//!
  
  
  ClassDef(StJetTreeThirdMaker, 1)   //StAF chain virtual base class for Makers
};






#endif
