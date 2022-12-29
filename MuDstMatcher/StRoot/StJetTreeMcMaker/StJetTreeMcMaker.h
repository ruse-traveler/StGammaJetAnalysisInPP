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
//Date  12/2015: "StJetTreeMcMaker" inherited from "StThirdMaker"
//              It makes "JetTree" for Full Jet Reconstruction
//              for Gamma+Jet analysis.      - Nihar r. Sahoo
//Date  06/2016: all events, tracks and tower attributes are added and 
//                                           - Nihar
//Date  07/2017: added histograms to calculate E/p for electrons
//               and hadrons.
//                                           - Derek
//Date  07/2017: added additional branches to calculate E/p for
//               electrons and hadrons.
//                                           - Derek
//Date  02/2018: added tree and variables to hold StMuMcTrack
//               variables
//                                           - Derek
//____________________________________________________________


#ifndef STAR_StJetTreeMcMaker
#define STAR_StJetTreeMcMaker

#ifndef StMaker_H
#include "StMaker.h"
#include <string>
#include <vector>
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/filters/StEmcFilter.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StTriggerUtilities/StTriggerSimuMaker.h"
#include "StMcEvent/StMcEventTypes.hh"
#include "StMcEvent/StMcContainers.hh"
#include "StarClassLibrary/StParticleDefinition.hh"

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
//class vector<float>;
//_______Class for Tree

class GammaJetEvent;
class GammaJetTrack;
class GammaJetTower;
class GammaJetTowerUtil;
class TClonesArray;
class  TLorentzVector;

class StMuDst;
class StMuEvent;
class StMuTrack;
class StMuMcTrack;  // [Derek, 02.19.2018] 
class StMuMcVertex;  // [Derek, 02.20.2018]
class StEmcPosition;
class StBemcTables; //v3.14


#define TowerHVchangeMax 400

const  Int_t max_pTracks = 10000;
const  Int_t max_gTracks = 10000;
/*
//#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
*/
//class PseudoJet;

class StJetTreeMcMaker : public StMaker 
{
 private: //BetheBloch      	mBB;
  TRandom*     mRandom;  

  St_db_Maker* mDbMaker;

  StEmcFilter*    mEmcFilter;
  StEvent           *mStEvent; 
  StMcEvent           *mcEvent;
  StEmcDecoder* mEmcDecoder;
  StMuDstMaker* mMuDstMaker;
  StMuEvent *muEvent;
  StMuDst* muDst;

  //  StMuDstMaker* mMuDstMaker;
  StEmcGeom       *mGeom_;
  StEmcGeom       *mGeom;
  StMuEmcCollection *muEmcCol;
  StEmcCollection *mEmcCol;
  StEmcFilter     *mFilter;
  StBemcTables    *mTables;
  Int_t    mCounter;

 protected:
  // Protected method if any
  
 public:
  void          setPrint(Bool_t);
  /* StJetTreeMcMaker(const char *name, StMuDstMaker* uDstMaker);*/
  StJetTreeMcMaker(const char *name="ThirdMaker", char* dataType = "");
  virtual       ~StJetTreeMcMaker();

  Int_t doMCEvent(Int_t RunIdOfEvent);
  Int_t  doMuEvent(Int_t & runIdEvnt);
  virtual Int_t Finish();
  void     Clear(Option_t *option="");
  virtual Int_t Init();
  virtual Int_t Make();
  /*Int_t associateTracksWithSmdeHits();*/
  void SetFileName( char* name){outFile=name;}; 
  Int_t ResizeHit(); 
  void  setDbMaker(St_db_Maker*);
  Bichsel* m_dEdxParameterization;
  StEmcFilter* getEmcFilter() { return mEmcFilter; };

  Bool_t Check_hot_Tower(Int_t TwrId);
  Bool_t Check_hot_EtaStrip(Int_t EtaStpId);

  //  void RunFastJet(vector<PseudoJet> particles);

  


  const StMuTrack* ptrack;  
  StEmcPosition *mPosition;

  //  StBemcTables    *mTables;


  
  ///////////////////////

  TH1F* mHistBg;

  TH1D* hTotalEbemc;
  TH1D* hEbemcTwr;

  // for E/p calculation [Derek, 07.18.2017;
  TH1D *hEbp_elec;
  TH1D *hEbp_hadr;
  // [Derek, 07.20.2017]
  TH1D *hPt_elec;
  TH1D *hPt_hadr;

  TH1F* hya;
  int  mEventCounter;
  TFile* File_output;
  char* outFile;


  // for comparing particle to detector level [Derek, 02.19.2018]
  TTree *_tMcTracks;

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

 
  ///________
  
  TBranch *eventBranch;
  TBranch *mueventBranch;
  
  Float_t *MuEvValues;
  Float_t *EvValues;
  
  TTree *outTree_gnt;
  TTree *outTree_mu;
  GammaJetEvent *event;
  GammaJetEvent *Muevent;

  Float_t mStatusT[4800];
  Float_t mCalibPedValues[4800];        //!
  Float_t mCalibPedRMSValues[4800];
  
  
  //_____________for matching tracks
  Int_t *pTrMatchArr;//!
  Float_t *pTriMatchEtaArr;//!
  Float_t *pTrMatchPhiArr;//!
  Int_t *pTrIndexArray;//!
  Float_t *pTrEtaArray;//!
  Float_t *pTrPhiArray;//!
  
  
  ClassDef(StJetTreeMcMaker, 1)   //StAF chain virtual base class for Makers
};






#endif
