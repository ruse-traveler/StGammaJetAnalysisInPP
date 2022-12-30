#ifndef STAR_StThirdMaker
#define STAR_StThirdMaker

#ifndef StMaker_H
#include "StMaker.h"
#include <string>

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
#define TowerHVchangeMax 400
class StThirdMaker : public StMaker 
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
 Bool_t mL2gamma_pp2009;
Bool_t    mbht2mb1;
Bool_t     mbht2mb2;
Bool_t     mbht2mb3;
Bool_t     mbht2mb4;
Bool_t     mbht2mb5;Bool_t     mbht2mb6;
Bool_t    mHt0; Bool_t    mht4;
Bool_t    mht4fast;Bool_t    mupsilondAu;
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
 /* StThirdMaker(const char *name, StMuDstMaker* uDstMaker);*/
  StThirdMaker(const char *name="ThirdMaker", char* dataType = "");
  virtual       ~StThirdMaker();
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  /*Int_t associateTracksWithSmdeHits();*/
 void SetFileName( char* name){outFile=name;}; 
 Int_t ResizeHit(); 
 void          setDbMaker(St_db_Maker*);
 Bichsel* m_dEdxParameterization;
 StEmcFilter* getEmcFilter() { return mEmcFilter; };

TNtuple*   mTupleBBC;TNtuple*  mTupleTrackS4;TNtuple*  mTupleTrackS5;
TNtuple*  mTupleBBC1;
 TNtuple* mTupleTri2;
 TNtuple* mTupleTri1;

 TNtuple*    mTupleTrackS;
 TNtuple*    mTupleTrackSS;
 TNtuple*    mTupleTrackSSS;

 //###########
 TNtuple*    mTupleTrackSO;
 TNtuple*    mTupleTrackSSO;
 TNtuple*    mTupleTrackSSSO;
 TNtuple* mTupleTrackS4O;
 TNtuple* mTupleTrackS5O;
 TNtuple* mTupleFTPCTrO;

 //##########
  TNtuple*    mTupleRP1;
  TNtuple*    mTupleRP2;
 TNtuple*    mTupleRP3;
 TNtuple*    mTupleRP4;
 TNtuple*    mTupleRP5;
 TNtuple*    mTupleFour;
 TNtuple*    mTupleFour2;
 TNtuple*    mTupleFive; 
 TNtuple*    mTupleFiveA;
 TNtuple*    mTupleSix;
  TNtuple*   mTupleSeven;
 TNtuple*    mTupleEight;
 TNtuple*    mTupleNine;	
  TNtuple*   mTupleTen;
  TNtuple*   mTupleEleven; 
  TNtuple*   mTupleTweleve; 
  TNtuple*   mTupleThirteen;  
  TNtuple*   mTupleFourteen; 
  TNtuple*   mTupleFiveteen; 
  TNtuple*   mTupleTri;
  TNtuple*   mTupleFiveL;
  TNtuple* mTupleSixteen;
  TNtuple* mTupleSeventeen;
 TNtuple*  mTupleEighteen;
  
  TNtuple* mTupleFiveteenB;
  TNtuple* mTupleSixteenB;
 TNtuple* mTupleSeventeenB;
 TNtuple* mTupleEighteenB;
 
  TNtuple* mTupleFiveteenA;
  TNtuple* mTupleSixteenA;
 TNtuple* mTupleSeventeenA;
 TNtuple* mTupleEighteenA;
 
  TNtuple* mTupleFTPCTr;
  TNtuple* mTupleFTPCG;
  TNtuple* mTupleFTPCQ;
  TNtuple* mTupleFTPCQA;
  
            TNtuple*      mTupleCAB;
         TNtuple*        	mTupleQAAB;
	 TNtuple* 	mTupleAngleAB;
	 TNtuple* 	mTupleSMDE1B;
	 TNtuple* 	mTuplePiB;
	 TNtuple* 	mTuplePhoB;
	 TNtuple* 	mTuplePi1B;
	 TNtuple* 	mTuplePi2B;
	   TNtuple*       mTupleNew1;
	   TNtuple*       mTupleNew2;
	   TNtuple*       mTupleNew3;
           TNtuple*       mTupleNew4;
            TNtuple*      mTupleNew5;
            TNtuple*      mTupleNew6;
           TNtuple*       mTupleNew7;
            TNtuple*      mTupleNew8;
  
             TNtuple*     mTupleTrack;
            TNtuple*      mTupleFiveLL;
            TNtuple*      mTupleCA;
            TNtuple*     	mTupleQAA;
	 TNtuple* 	mTupleAngleA;
	 TNtuple* 	mTupleSMDE1;
	 TNtuple* 	mTuplePi;
	 TNtuple* 	mTuplePho;
	 TNtuple* 	mTuplePi1;
	 TNtuple* 	mTuplePi2;
	    TNtuple*      mTupleNew1A;
	    TNtuple*      mTupleNew2A;
	     TNtuple*     mTupleNew3A;
            TNtuple*      mTupleNew4A;
            TNtuple*      mTupleNew5A;
            TNtuple*      mTupleNew6A;
            TNtuple*      mTupleNew7A;
            TNtuple*      mTupleNew8A;
 TNtuple*      mTupleNew9A;
 
 
       TNtuple*          mTupleTrack3;
       TNtuple*          mTupleFiveLL3;
        TNtuple*         mTupleCA3;
        TNtuple*        	mTupleQAA3;
	 TNtuple*	mTupleAngleA3;
	 TNtuple*	mTupleSMDE13;
	 TNtuple*	mTuplePi3;
	 TNtuple*	mTuplePho3;
	 TNtuple*	mTuplePi13;
	 TNtuple*	mTupleNew8A3;
         TNtuple*        mTupleNew9A3;
         TNtuple*        mTupleFiveteenA3;
         TNtuple*        mTupleSixteenA3;
         TNtuple*        mTupleSeventeenA3;
          TNtuple*       mTupleEighteenA3;
            TNtuple*     mTupleNew7A3;
          TNtuple*       mTupleNew6A3;
           TNtuple*      mTupleFTPCQA3;

	   //#_______________ 
 
	   TNtuple* mTupleTrack_M;
           TNtuple* mTupleFiveLL_M;
           TNtuple* mTupleCA_M;
           TNtuple* mTupleQAA_M;
           TNtuple* mTupleAngleA_M;

           TNtuple* mTupleSMDE1_M;
           TNtuple* mTuplePi_M;
           TNtuple* mTuplePho_M;

           TNtuple* mTuplePi1_M;
           TNtuple* mTupleNew8A_M;
           TNtuple* mTupleNew9A_M;

           TNtuple* mTupleFiveteenA_M;

           TNtuple* mTupleSixteenA_M;
           TNtuple* mTupleSeventeenA_M;

           TNtuple* mTupleEighteenA_M;


           TNtuple* mTupleNew7A_M;
           TNtuple* mTupleNew6A_M;

           TNtuple* mTupleFTPCQA_M;
	   TNtuple*  mTupleQATrkInfo;
 
 
 
 
 
TH1F* mHist;
TH1F* mHistBg;
TH1F* hya; 
int        mEventCounter;
TFile* myGraphFile;
 char* outFile;
  
  ClassDef(StThirdMaker, 1)   //StAF chain virtual base class for Makers
};

#endif
