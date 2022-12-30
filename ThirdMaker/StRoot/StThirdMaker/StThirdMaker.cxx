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
// 
//
//____________________________________________________________


#include "StThirdMaker.h"
#include "StChain.h"
#include "StEvent.h"
#include "StEventTypes.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/filters/StEmcFilter.h"
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include <vector>
#include <math.h>
#include "TNtuple.h"
#include <TString.h>
#include "TRandom.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StEventUtilities/StuFtpcRefMult.hh"
#include "StThreeVector.hh"
#include "StEvent/StEnumerations.h"
#include "StMcEvent.hh"
#include "StMcVertex.hh"
#include "StMcTrack.hh"
#include "StMcEventMaker/StMcEventMaker.h"
#include "StMcEventTypes.hh"
#include "StPrimaryVertex.h"
#include "StBTofHeader.h"

#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>

#include "St_db_Maker/St_db_Maker.h"
#include "StDetectorDbMaker/StDetectorDbTriggerID.h"

#include <emcStatus.h>
#include "tables/St_emcStatus_Table.h"

#include "tables/St_emcPed_Table.h"
//#include "L2Result.h"
#include "StTriggerData.h"
#include "StTriggerData2007.h"
//#include "StDaqLib/TRG/trgStructures.h"
#include "StDaqLib/TRG/trgStructures2007.h"
//#include "StDaqLib/TRG/L2pedResults2006.h"
//#include "StDaqLib/TRG/L2gammaResult2007.h"
#include "StDaqLib/EMC/StEmcDecoder.h"
#include "L2gammaResult.h"
#include "StEvent/StTriggerIdCollection.h"

#include "StEvent/StTriggerId.h"


//StMuDst
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StTriggerData.h"
//------------------------------------------------------------------------------------
void StThirdMaker::setDbMaker(St_db_Maker* dbMaker)
{
  mDbMaker = dbMaker;
}
//-----------------------------------------------------------------------------------
StThirdMaker::StThirdMaker(const char *name,char *dataType):StMaker(name){ 
 mDbMaker = 0;
for(int n=0;n<150;n++) mEnergyE[n]=0.0; for(int Q=0;Q<10;Q++) for(int G=0;G<15;G++)  mEnergyP[Q][G]=0.0;

for(int n2=0;n2<150;n2++) mEnergyE2[n2]=0.0; for(int Q2=0;Q2<10;Q2++)for(int G2=0;G2<15;G2++)  mEnergyP2[Q2][G2]=0.0;

for(int D=0;D<20;D++) for(int S=0;S<2;S++)    mEnergyTT[D][S]=0.0; 	  for(int D1=0;D1<20;D1++)for(int S1=0;S1<2;S1++)  mEnergyTT1[D1][S1]=0.0; 
for(int D2=0;D2<20;D2++)for(int S2=0;S2<2;S2++)  mEnergyTT2[D2][S2]=0.0;for(int C=0;C<6;C++)  energy[C]=0.0;

for(int D3=0;D3<20;D3++) for(int S3=0;S3<2;S3++)    mEnergyTT3[D3][S3]=0.0; for(int D14=0;D14<20;D14++)for(int S14=0;S14<2;S14++)  mEnergyTT14[D14][S14]=0.0; 
for(int D25=0;D25<20;D25++)for(int S25=0;S25<2;S25++)  mEnergyTT25[D25][S25]=0.0;for(int C2=0;C2<6;C2++)  energy2[C2]=0.0;

for(int RW=0; RW<4800;RW++) ParentTrackP[RW]=0.0;
//for(int RW1=0; RW1<18000;RW1++) ParentTrackdEdx[RW1]=0.0;
for(int RW2=0; RW2<4800;RW2++) ParticleId[RW2]=0.0;
for(int RW3=0; RW3<4800;RW3++) ParentTrackPOSe[RW3]=0.0;
for(int RW4=0; RW4<4800;RW4++) 	ParentTrackPOSP[RW4]=0.0;
//for(int RW3=0; RW3<18000;RW3++) ParentTrackEta[RW3]=0;
//for(int RW4=0; RW4<18000;RW4++) ParentTrackCh[RW4]=0;
//for(int DDD=0;DDD<8;DDD++) mHighest[DDD]=0.0;
for(int RW5=0; RW5<43;RW5++) ReactionPlane[RW5]=0.0;
for(int GFac=0;GFac<2;GFac++) for(int GFac1=0;GFac1<16;GFac1++) bbcGainFac[GFac][GFac1]=0.0;

//+++++++++++
 mDst = 0;   // NRS: intialization

  outFile="pbtestLastSSs1.root";
}

//_____________________________________________________________________________
StThirdMaker::~StThirdMaker(){
  //
}

//_____________________________________________________________________________
// Init - is a first method the top level StChain 
// calls to initialize all its makers 
Int_t StThirdMaker::Init(){
  // Create Histograms 
 
   mRandom = new TRandom();
  myGraphFile= new TFile(outFile,"RECREATE");
  
  
   mTupleBBC=new TNtuple("BBC","BBC","psie:psiw:psif:BBCAdcCut:psifn:psif0:psif1:psifn0:psifn1:BBCAdcCut0");
   mTupleBBC1=new TNtuple("BBC1","BBC1","iTile:bbcEast:bbcWest");
  
   mTupleTri1=new TNtuple("Tri1","Tri1","noOfT:noOfGo:noOfGoP:noOfGoN:noOfGoPr:noOfEx:noOfVer:noOfPiVer:Meanpt:Meanpt2:MeanEta:MgF");
 
   mTupleTrackS=new TNtuple("TrackS","TrackS","pTPC:ptTPC:dEdx:ch:Ppo:epsi:flag:RefM:runId:z:evt:px:py:pz");
   mTupleTrackSS=new TNtuple("TrackSS","TrackSS","EtaTr:Fp:T:dca:rank:rank1:nV:muRef:psie:psiw:psifn:psifn0:psifn1");
   mTupleTrackSSS=new TNtuple("TrackSSS","TrackSSS","FpT:YY:XX:angleRe:angleRef:psif:v2chT:v2chF:v2chB:BBCAdcCut:BBCAdcCut0:psif0:psif1");
   mTupleTrackS4=new TNtuple("TrackS4","TrackS4","angleRefS1:angleRefS2:angleReS1:angleReS2:angleRept:angleReS1pt:angleReS2pt:ltrkV:lpnt:pdca:dcaD:dcaZ");
   mTupleTrackS5=new TNtuple("TrackS5","TrackS5","angleRepe:angleReS1pe:angleReS2pe:angleRene:angleReS1ne:angleReS2ne:angleRenc:angleReS1nc:angleReS2nc:angleRepc:angleReS1pc:angleReS2pc");
   mTupleFTPCTr=new TNtuple("FTPCTr","FTPCTr","angleRefW:angleRefE:angleRefW1:angleRefE1:angleRefW2:angleRefE2:vpdz:BBCz:ZDCz");

   //############ OLD INFO to CHECK
   mTupleTrackSO=new TNtuple("TrackSO","TrackSO","pTPC:ptTPC:dEdx:ch:Ppo:epsi:flag:RefM:runId:z:evt:px:py:pz");
   mTupleTrackSSO=new TNtuple("TrackSSO","TrackSSO","EtaTr:Fp:T:dca:rank:rank1:nV:muRef:psie:psiw:psifn:psifn0:psifn1");
   mTupleTrackSSSO=new TNtuple("TrackSSSO","TrackSSSO","FpT:YY:XX:angleRe:angleRef:psif:v2chT:v2chF:v2chB:BBCAdcCut:BBCAdcCut0:psif0:psif1");
   mTupleTrackS4O=new TNtuple("TrackS4O","TrackS4O","angleRefS1:angleRefS2:angleReS1:angleReS2:angleRept:angleReS1pt:angleReS2pt:ltrkV:lpnt:pdca:dcaD:dcaZ");
   mTupleTrackS5O=new TNtuple("TrackS5O","TrackS5O","angleRepe:angleReS1pe:angleReS2pe:angleRene:angleReS1ne:angleReS2ne:angleRenc:angleReS1nc:angleReS2nc:angleRepc:angleReS1pc:angleReS2pc");
   mTupleFTPCTrO=new TNtuple("FTPCTrO","FTPCTrO","angleRefW:angleRefE:angleRefW1:angleRefE1:angleRefW2:angleRefE2:vpdz:BBCz:ZDCz");
   //############ OLD INFO to CHECK

   mTupleFour=new TNtuple("Four","Four","Energy:Sub:Module:Eta:Adc:T:evt:z:runId:MNeg:MPos:Id:phiT:Theta");
   mTupleFour2=new TNtuple("Four2","Four2","hu4:hped4:hRMS4:eT:TC4:PPF:pidF:TPe:TPp");
   mTupleFive=new TNtuple("Five","Five","hiTowerB:SubhB:ModulehB:EtahB:AdchB:T:huB:hpedB:runId:IdhB:eThB:phiThB:ThetahB:hRMSB"); 
   mTupleFiveA=new TNtuple("FiveA","FiveA","hiTower:Subh:Moduleh:Etah:Adch:T:evt:z:runId:Idh:eTh:phiTh:Thetah:PPTPC");
   mTupleFiveteen=new TNtuple("Fiveteen","Fiveteen","angleRe:angleReS1:angleReS2:angleRef:angleRefS1:angleRefS2:NumOfTracks:Q1:Q2:XR:YR:XRf:YRf"); 
   mTupleTri=new TNtuple("Tri","Tri","T:evt:z:runId:NOP:XX:YY:MPos:MNeg:allGlob:allPrim:ctbsum:zdcsum");
   mTupleTri2=new TNtuple("Tri2","Tri2","muRef:nV:rank1:angleReS2nc:angleRepc:angleReS1pc:angleReS2pc:vpdz:BBCz:ZDCz");
   mTupleFiveL=new TNtuple("FiveL","FiveL","hiTower:Subh:Moduleh:Etah:Adch:T:evt:z:runId:Idh:eTh:phiTh:Thetah:PPTPC"); 
   mTupleSixteen=new TNtuple("Sixteen","Sixteen","XRS1:XRS2:YRS1:YRS2:Q:Q3:Q4:Q1S1:Q1S2:Q2S1:Q2S2:Q3S1:Q3S2:Q4S1");
   mTupleSeventeen=new TNtuple("Seventeen","Seventeen","Q4S2:XRfS1:XRfS2:YRfS1:YRfS2:XRfW:YRfW:XRfE:YRfE:XRfSW1:YRfSW1:XRfSW2:YRfSW2:XRfSE1");
   mTupleEighteen=new TNtuple("Eighteen","Eighteen","YRfSE1:XRfSE2:YRfSE2:angleRept:angleReS1pt:angleReS2pt:angleRepe:angleReS1pe:angleReS2pe:angleRene:angleReS1ne:angleReS2ne:angleRenc:angleReS1nc");
   mTupleFTPCG=new TNtuple("FTPCG","FTPCG","angleRefW:angleRefE:angleRefW1:angleRefE1:angleRefW2:angleRefE2");
 
   
   mTupleCAB=new TNtuple("CAB","CAB","didE:ENET0:XX:YY:adc11:didP:adcE1:adcP1:eT:phT:hSS:hpedSS:hRMSSS:TC7");
   mTupleQAAB=new TNtuple("QAAB","QAAB","EneT0:MPos:MNeg:H:en4:enp4:eneT0:z:etav1:phiv1:didT:delta1:dphi1:Weta1");
   mTupleAngleAB=new TNtuple("AngleTAB","AngleTAB","moduleT:EtaE:sP:runId:evt:NOP:Wphi1:enp7:enp1:enp8:enp0:en01:en02:PTower");
   mTupleSMDE1B=new TNtuple("SMDE1B","SMDE1B","en3:en5:en2:en6:en1:en7:en0:en8:enp3:enp5:enp2:enp6:en03:en04");
   mTuplePiB=new TNtuple("PiB","PiB","en9:en10:en11:en12:enp01:enp02:enp03:enp04:enp9:enp10:enp11:enp12:pidTower:Enpq1");
   mTuplePhoB=new TNtuple("PhoB","PhoB","sT:THeta:Ennq1:Enpq20:Ennq20:Enpq19:Ennq19:Enpq21:NumOfTracks:muRef:nV:rank1:rank:noOfPiVer");
   mTuplePi1B=new TNtuple("Pi1B","Pi1B","en05:enp05:en06:enp06:en07:enp07:en08:enp08:en13:enp13:en14:enp14:en15:enp15");
		
   
   mTupleNew8=new TNtuple("New8","New8","psif0:psif1:psifn0:psifn1:angleu0v:angleu1v:angleu2v:angleu3v:angleu4v:angleu5v:angleu6v:angleu7v");
   mTupleFiveteenB=new TNtuple("FiveteenB","FiveteenB","angleRe:angleReS1:angleReS2:angleRef:angleRefS1:angleRefS2:NumOfTracks:Q1:Q2:XR:YR:XRf:YRf");
   mTupleSixteenB=new TNtuple("SixteenB","SixteenB","XRS1:XRS2:YRS1:YRS2:Q:Q3:Q4:Q1S1:Q1S2:Q2S1:Q2S2:Q3S1:Q3S2:Q4S1");
   mTupleSeventeenB=new TNtuple("SeventeenB","SeventeenB","Q4S2:XRfS1:XRfS2:YRfS1:YRfS2:XRfW:YRfW:XRfE:YRfE:XRfSW1:YRfSW1:XRfSW2:YRfSW2:XRfSE1");
   mTupleEighteenB=new TNtuple("EighteenB","EighteenB","YRfSE1:XRfSE2:YRfSE2:psie:psiw:psif:BBCAdcCut:vf2:vf2F:psifn:v2B:Ehighest:BBCAdcCut0");
   
   mTupleNew7=new TNtuple("New7","New7","angleRept:angleReS1pt:angleReS2pt:vpdz:BBCz:ZDCz");
   mTupleNew6=new TNtuple("New6","New6","angleRepe:angleReS1pe:angleReS2pe:angleRene:angleReS1ne:angleReS2ne:angleRenc:angleReS1nc:angleReS2nc:angleRepc:angleReS1pc:angleReS2pc");
   mTupleFTPCQ=new TNtuple("FTPCQ","FTPCQ","angleRefW:angleRefE:angleRefW1:angleRefE1:angleRefW2:angleRefE2");
   
   
   
   
  
   mTupleTrack=new TNtuple("Track","Track","pTPC:ptTPC:dEdx:EtaTr:epsi:ch:dca:Fp:eE:phE:ETAP:PHIP:eT:phT");
   mTupleFiveLL=new TNtuple("FiveLL","FiveLL","hiTower:Subh:Moduleh:Etah:Adch:T:evt:z:runId:Idh:eTh:phiTh:Thetah:PPTPC");
   mTupleCA=new TNtuple("CA","CA","didE:ENET0:XX:YY:adc11:didP:adcE1:adcP1:eT:phT:hSS:hpedSS:hRMSSS:TC7");
   mTupleQAA=new TNtuple("QAA","QAA","EneT0:MPos:MNeg:H:en4:enp4:eneT0:z:etav1:phiv1:didT:delta1:dphi1:Weta1");
   mTupleAngleA=new TNtuple("AngleTA","AngleTA","moduleT:EtaE:sP:runId:evt:NOP:Wphi1:enp7:enp1:enp8:enp0:en01:en02:PTower");
   mTupleSMDE1=new TNtuple("SMDE1","SMDE1","en3:en5:en2:en6:en1:en7:en0:en8:enp3:enp5:enp2:enp6:en03:en04");
   mTuplePi=new TNtuple("Pi","Pi","en9:en10:en11:en12:enp01:enp02:enp03:enp04:enp9:enp10:enp11:enp12:pidTower:Enpq1");
   mTuplePho=new TNtuple("Pho","Pho","sT:eff:energ:dphif:energSt:THeta:Ennq1:Enpq20:Ennq20:Enpq19:Ennq19:Enpq21:dphiin:dphiou");
   mTuplePi1=new TNtuple("Pi1","Pi1","en05:enp05:en06:enp06:en07:enp07:en08:enp08:en13:enp13:en14:enp14:en15:enp15");
                  
   mTupleNew8A=new TNtuple("New8A","New8A","psif0:psif1:psifn0:psifn1:Ppo:xg:yg:zg:px:py:pz:pdca:dcaD:dcaZ");
   mTupleNew9A=new TNtuple("New9A","New9A","angleu0v:angleu1v:angleu2v:angleu3v:angleu4v:angleu5v:angleu6v:angleu7v:NumOfTracks:muRef:nV:rank1:rank:noOfPiVer");
   mTupleFiveteenA=new TNtuple("FiveteenA","FiveteenA","angleRe:angleReS1:angleReS2:angleRef:angleRefS1:angleRefS2:NumOfTracks:Q1:Q2:XR:YR:XRf:YRf");
   mTupleSixteenA=new TNtuple("SixteenA","SixteenA","XRS1:XRS2:YRS1:YRS2:Q:Q3:Q4:Q1S1:Q1S2:Q2S1:Q2S2:Q3S1:Q3S2:Q4S1");
   mTupleSeventeenA=new TNtuple("SeventeenA","SeventeenA","Q4S2:XRfS1:XRfS2:YRfS1:YRfS2:XRfW:YRfW:XRfE:YRfE:XRfSW1:YRfSW1:XRfSW2:YRfSW2:XRfSE1");
   mTupleEighteenA=new TNtuple("EighteenA","EighteenA","YRfSE1:XRfSE2:YRfSE2:psie:psiw:psif:BBCAdcCut:vf2:vf2F:psifn:v2B:Ehighest:BBCAdcCut0");
   
   mTupleNew7A=new TNtuple("New7A","New7A","angleRept:angleReS1pt:angleReS2pt:vpdz:BBCz:ZDCz:TSP:eh");
   mTupleNew6A=new TNtuple("New6A","New6A","angleRepe:angleReS1pe:angleReS2pe:angleRene:angleReS1ne:angleReS2ne:angleRenc:angleReS1nc:angleReS2nc:angleRepc:angleReS1pc:angleReS2pc");
   mTupleFTPCQA=new TNtuple("FTPCQA","FTPCQA","angleRefW:angleRefE:angleRefW1:angleRefE1:angleRefW2:angleRefE2");

   mTupleTrack_M=new TNtuple("Track_M","Track_M","pTPC:ptTPC:dEdx:EtaTr:epsi:ch:dca:Fp:eE:phE:ETAP:PHIP:eT:phT");
   mTupleFiveLL_M=new TNtuple("FiveLL_M","FiveLL_M","hiTower:Subh:Moduleh:Etah:Adch:T:evt:z:runId:Idh:eTh:phiTh:Thetah:PPTPC");
   mTupleCA_M=new TNtuple("CA_M","CA_M","didE:ENET0:XX:YY:adc11:didP:adcE1:adcP1:eT:phT:hSS:hpedSS:hRMSSS:TC7");
   mTupleQAA_M=new TNtuple("QAA_M","QAA_M","EneT0:MPos:MNeg:H:en4:enp4:eneT0:z:etav1:phiv1:didT:delta1:dphi1:Weta1");
   mTupleAngleA_M=new TNtuple("AngleTA_M","AngleTA_M","moduleT:EtaE:sP:runId:evt:NOP:Wphi1:enp7:enp1:enp8:enp0:en01:en02:PTower");
   mTupleSMDE1_M=new TNtuple("SMDE1_M","SMDE1_M","en3:en5:en2:en6:en1:en7:en0:en8:enp3:enp5:enp2:enp6:en03:en04");
   mTuplePi_M=new TNtuple("Pi_M","Pi_M","en9:en10:en11:en12:enp01:enp02:enp03:enp04:enp9:enp10:enp11:enp12:pidTower:Enpq1");
   mTuplePho_M=new TNtuple("Pho_M","Pho_M","sT:eff:energ:dphif:energSt:THeta:Ennq1:Enpq20:Ennq20:Enpq19:Ennq19:Enpq21:dphiin:dphiou");
   mTuplePi1_M=new TNtuple("Pi1_M","Pi1_M","en05:enp05:en06:enp06:en07:enp07:en08:enp08:en13:enp13:en14:enp14:en15:enp15");
   mTupleNew8A_M=new TNtuple("New8A_M","New8A_M","psif0:psif1:psifn0:psifn1:Ppo:xg:yg:zg:px:py:pz");
   mTupleNew9A_M=new TNtuple("New9A_M","New9A_M","angleu0v:angleu1v:angleu2v:angleu3v:angleu4v:angleu5v:angleu6v:angleu7v:NumOfTracks:muRef:nV:rank1:rank:noOfPiVer");
   mTupleFiveteenA_M=new TNtuple("FiveteenA_M","FiveteenA_M","angleRe:angleReS1:angleReS2:angleRef:angleRefS1:angleRefS2:NumOfTracks:Q1:Q2:XR:YR:XRf:YRf");
   mTupleSixteenA_M=new TNtuple("SixteenA_M","SixteenA_M","XRS1:XRS2:YRS1:YRS2:Q:Q3:Q4:Q1S1:Q1S2:Q2S1:Q2S2:Q3S1:Q3S2:Q4S1");
   mTupleSeventeenA_M=new TNtuple("SeventeenA_M","SeventeenA_M","Q4S2:XRfS1:XRfS2:YRfS1:YRfS2:XRfW:YRfW:XRfE:YRfE:XRfSW1:YRfSW1:XRfSW2:YRfSW2:XRfSE1");
   mTupleEighteenA_M=new TNtuple("EighteenA_M","EighteenA_M","YRfSE1:XRfSE2:YRfSE2:psie:psiw:psif:BBCAdcCut:vf2:vf2F:psifn:v2B:Ehighest:BBCAdcCut0");
   mTupleNew7A_M=new TNtuple("New7A_M","New7A_M","angleRept:angleReS1pt:angleReS2pt:vpdz:BBCz:ZDCz:TSP:eh");
   mTupleNew6A_M=new TNtuple("New6A_M","New6A_M","angleRepe:angleReS1pe:angleReS2pe:angleRene:angleReS1ne:angleReS2ne:angleRenc:angleReS1nc:angleReS2nc:angleRepc:angleReS1pc:angleReS2pc");
   mTupleFTPCQA_M=new TNtuple("FTPCQA_M","FTPCQA_M","angleRefW:angleRefE:angleRefW1:angleRefE1:angleRefW2:angleRefE2");

   mTupleQATrkInfo = new TNtuple("QATrkInfo","QATrkInfo","sttrk_prim:stmutrk_prim:z:rank:rank1:nV");
   

// m_dEdxParameterization=new Bichsel();

 mMuDstMaker = (StMuDstMaker*)GetMaker("MuDst");
    return StMaker::Init();
}

//_____________________________________________________________________________
/// Make - this method is called in loop for each mStEvent
Int_t StThirdMaker::Make()
{  
StMuEvent *muEvent= mMuDstMaker->muDst()->event();

Int_t nV=mMuDstMaker->muDst()->numberOfPrimaryVertices();
Float_t rank = -99;
if (nV > 0) {  // Check that we have at least one vertex
  rank =mMuDstMaker->muDst()->primaryVertex(0)->ranking();} //cout<<" rank= " << rank <<endl;

Int_t muRef=muEvent->refMult(); 

 mStEvent = (StEvent*)GetInputDS("StEvent");     
 if(!mStEvent)    {
   //cout <<"Make() Can't get Event pointer" << endl; 
   return kStWarn;}	

  float z(0),YY(0),XX(0);int NOP(0),evt(0), runId(0);
    if(!(mStEvent->numberOfPrimaryVertices()>=1 && fabs(mStEvent->primaryVertex()->position().z())<1000.0)) return kStOK;
  NOP=mStEvent->numberOfPrimaryVertices();   Float_t rank1 = -99;
 
  rank1 =mStEvent->primaryVertex(0)->ranking(); //cout<<" rank1= " << rank1 <<endl;
   

  
   z=mStEvent->primaryVertex()->position().z();     XX=mStEvent->primaryVertex()->position().x();   YY=mStEvent->primaryVertex()->position().y();
    cout<<"vz is"<<z<<endl;    StThreeVectorF vertexPosition=mStEvent->primaryVertex()->position();
   cout<<"vertexPosition is"<<vertexPosition<<endl;      evt=mStEvent->id();      runId=mStEvent->runId();
    Int_t MPos=uncorrectedNumberOfPositivePrimaries(*mStEvent);  Int_t MNeg=uncorrectedNumberOfNegativePrimaries(*mStEvent);

    ///___________________________Check______________
    double evnt_muE(0), runId_muE(0);
    evnt_muE=muEvent->eventId(); runId_muE=muEvent->runId();
    double Vx_muE(0), Vy_muE(0), Vz_muE(0);
    StThreeVectorF Vtx=muEvent->primaryVertexPosition();
    Vz_muE = Vtx.z();
    Vx_muE = Vtx.x();
    Vy_muE = Vtx.y();


    cout << "run and evt = (StEvent) " << runId << " " << evt << " V(x:y:z) "<<XX<<" : "<<YY<<" : "<<z<<endl;
    cout << "run and evt = (StmuEvent) " << runId_muE << " " << evnt_muE << "V(x:y:z) "<<Vx_muE<<" : "<<Vy_muE<<" : "<<Vz_muE<<endl;
    ///_________________________________________


    //    cout<<"RefMultPos is"<<MPos <<endl;    cout<<"RefMultNeg is"<<MNeg <<endl;   // if(!((MPos+MNeg)>=222)) return kStOK; 
    //5% Refmult cut(>=)373, 10% 313, 20% 222, 30% 154, 40% 102, 50% 65, 60% 38, 70% 20, 80% 9.................// 
 
   StEventSummary *summary=mStEvent->summary();
    if ( ! summary ){      cout << "StThirdMaker::Make : No event summary" << endl;   return kStWarn; }	
	
    int    noOfT=summary->numberOfTracks() ;    int    noOfGo=summary->numberOfGoodTracks() ;    int    noOfGoP=summary->numberOfGoodTracks(positive) ;
    int    noOfGoN=summary->numberOfGoodTracks(negative) ;    int    noOfGoPr=summary->numberOfGoodPrimaryTracks() ;
    int    noOfEx=summary->numberOfExoticTracks() ;    int    noOfVer=summary->numberOfVertices() ;    int    noOfPiVer=summary->numberOfPileupVertices() ;
    float  Meanpt=summary->meanPt() ;    float  Meanpt2=summary->meanPt2() ;    float  MeanEta=summary->meanEta() ;   
    double MgF=summary->magneticField() ;

    mTupleTri1->Fill(noOfT,noOfGo,noOfGoP,noOfGoN,noOfGoPr,noOfEx,noOfVer,noOfPiVer,Meanpt,Meanpt2,MeanEta,MgF);	
	
	//-------------------------------------------------------------------------------
   StTriggerIdCollection *trgcol = mStEvent->triggerIdCollection();
    if ( ! trgcol ){      cout << "StThirdMaker::Make : No triggerIdCollection" << endl;    } 

   //--------------------------------------------------------------------------------------------

    int allGlob = 0;       int allPrim = 0;    
    StTrack *track;    StTrack *primarytrack;    StSPtrVecTrackNode& nodes = mStEvent->trackNodes();
    for (unsigned int j=0; j<nodes.size(); j++) {        track = nodes[j]->track(global);
        primarytrack=nodes[j]->track(primary);        if (track) allGlob++;       // if (accept(track)) goodGlob++;
        if (primarytrack) allPrim++;        //if (accept(primarytrack)) goodPrim++;    
	}
	cout<<"::allGlob::"<<allGlob<<"::allPrim::"<<allPrim<<endl;
    
     //-----------------------------------------------------  Get the ZDC and CTB data.-----------------------------------
        StTriggerDetectorCollection *theTriggers =mStEvent->triggerDetectorCollection();
    if (!theTriggers){      // good idea to check if the data is available at all
      gMessMgr->Warning() << "StAnalysisMaker::Make : no triggerDetectorCollection" << endm;      return kStOK;    }
  
    StVpdTriggerDetector &theVpd = theTriggers->vpd(); 
    StCtbTriggerDetector &theCtb = theTriggers->ctb();    
    StZdcTriggerDetector &theZdc = theTriggers->zdc(); 
    StBbcTriggerDetector &theBBC = theTriggers->bbc();
   
  
 
   //these 2 lines for 2010 and add 
   //#include "StMuDSTMaker/COMMON/StMuTofHit.h"
//#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
   StBTofHeader* tofHeader = mMuDstMaker->muDst()->btofHeader();
    double vpdz(0);vpdz =tofHeader->vpdVz(0);
   
   //ZDC
   float ZDCz = theZdc.vertexZ(); //SAVE THOSE
   float BBCz= theBBC.zVertex();
   
    Float_t zdce = -1.;     Float_t zdcw = -1.;
    zdce = theZdc.adcSum(east);     zdcw = theZdc.adcSum(west);
    float zdcsum(0);zdcsum=theZdc.adcSum();
    //       cout << " zdce = " << zdce << " zdcw = " << zdcw << " zdcsum = " << zdcsum <<endl;
       //  Sum all CTB counter
    float ctbsum = 0;    for (unsigned int islat=0; islat<theCtb.numberOfSlats(); islat++)
        for (unsigned int itray=0; itray<theCtb.numberOfTrays(); itray++)            ctbsum += theCtb.mips(itray, islat, 0);
    //cout << " ctbsum = " << ctbsum <<endl;

//for(int GFac=0;GFac<2;GFac++) for(int GFac1=0;GFac1<16;GFac1++)
bbcGainFac[0][0] = 0.962841;bbcGainFac[0][1] = 0.979108;bbcGainFac[0][2] = 1.03717;bbcGainFac[0][3] = 0.957205;
bbcGainFac[0][4] = 1.07229;bbcGainFac[0][5] = 0.991379;bbcGainFac[0][6] = 0.730763;bbcGainFac[0][7] = 1.05997;
bbcGainFac[0][8] = 1.04575;bbcGainFac[0][9] = 1.57317;bbcGainFac[0][10] = 1.02486;bbcGainFac[0][11] = 0.80374;
bbcGainFac[0][12] = 0.919484;bbcGainFac[0][13] = 1.02191;bbcGainFac[0][14] = 1.35782;bbcGainFac[0][15] = 0.928025;

bbcGainFac[1][0] = 0.930545;bbcGainFac[1][1] = 1.01045;bbcGainFac[1][2] = 0.988749;bbcGainFac[1][3] = 1.03118;
bbcGainFac[1][4] = 0.943601;bbcGainFac[1][5] = 1.09548;bbcGainFac[1][6] = 0.771513;bbcGainFac[1][7] = 1.00259;
bbcGainFac[1][8] = 0.972948;bbcGainFac[1][9] = 1.48775;bbcGainFac[1][10] = 0.91112;bbcGainFac[1][11] = 0.772221;
bbcGainFac[1][12] = 0.985141;bbcGainFac[1][13] = 1.12275;bbcGainFac[1][14] = 1.42478;bbcGainFac[1][15] = 1.00545;


//bbcGainFac[2][16] = {{0.962841,0.979108,1.03717,0.957205,1.07229,0.991379,0.730763,1.05997,1.04575,1.57317,1.02486,0.80374,0.919484,
//1.02191,1.35782,0.928025},{0.930545,1.01045,0.988749,1.03118,0.943601,1.09548,0.771513,1.00259,0.972948,1.48775,0.91112,0.772221,0.985141,1.12275,1.42478,1.00545}}; 

//BBC
 Float_t bbcEast=0.0;
 Float_t bbcWest=0.0;
//Gain Corrected Distribution
  for (int iTile= 0; iTile<16; iTile++) {
  bbcEast=(float)theBBC.adc(iTile)/bbcGainFac[0][iTile];
  bbcWest=(float)theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
  //  cout << " iTile = " << iTile << " adc_East = " << bbcEast << "  iTile+24 = " << iTile+24 << " adc_Wset = " << bbcWest << endl;
   mTupleBBC1->Fill(iTile,bbcEast,bbcWest);
  
  //cout << " iTile = " << iTile << " adc "<< endl;
  //pFlowEvent->SetBBCAdc(0, iTile, bbcEast);
  //pFlowEvent->SetBBCAdc(1, iTile, bbcWest);
 }
 
int mBBCAdcCut(0),mBBCAdcCut0(0);
 for (int iTile= 0; iTile<16; iTile++) {
           if(theBBC.adc(iTile) > 245 || theBBC.adc(iTile+24)>245) mBBCAdcCut++ ;
           if(theBBC.adc(iTile) > 120 || theBBC.adc(iTile+24)>120) mBBCAdcCut0++;
       
       // return kFALSE;
      
  }
 //  cout << " mBBCAdcCut = " << mBBCAdcCut << endl;
 double twopi=2*3.14285714;
//const float phi_div=twopi/12.;
//float bbc_phi=phi_div; 

for (int iTile= 0; iTile<16; iTile++) bbc_phi[iTile]=0;
//   Float_t  phi[24] = {1.5708,0.523599,5.75959,4.71239,3.66519,2.61799,1.5708,1.5708,0.523599,0,5.75959,4.71239,4.71239,3.66519,3.14159,2.61799,
//                     1.5708,0,4.71239,3.14159,1.5708,0,4.71239,3.14159};

 //  Float_t phi[16] = {1.5708,0.5236,5.75959,4.71239,3.66519,2.61799,1.5708,1.5708,0.523599,0.0,5.75959,
 //  4.71239,4.71239,3.66519,3.14159,2.61799};


for (int iTile= 0; iTile<16; iTile++){ 
if (iTile==0)bbc_phi[iTile]=3*twopi/12.; //tile 1 correct 90
if (iTile==1)bbc_phi[iTile]=twopi/12.;//tile 2 correct 30
if (iTile==2)bbc_phi[iTile]=-1*twopi/12.;//tile 3 correct -30
if (iTile==3)bbc_phi[iTile]=-3*twopi/12.;//tile 4 correct -90
if (iTile==4)bbc_phi[iTile]=-5*twopi/12.;//tile 5 correct -150
if (iTile==5)bbc_phi[iTile]=5*twopi/12.;//tile 6 correct 150
if (iTile==6)bbc_phi[iTile]= (gRandom->Rndm()>0.5) ? 2*twopi/12.:4*twopi/12.; //tile 7 correct 60 or 120
if (iTile==7)bbc_phi[iTile]=3*twopi/12.;//tile 8 correct 90
if (iTile==8)bbc_phi[iTile]=twopi/12.;//tile 9 correct 30
if (iTile==9)bbc_phi[iTile]=0.; //tile 10 correct 0
if (iTile==10)bbc_phi[iTile]=-1*twopi/12.;//tile 11 correct -30
if (iTile==11)bbc_phi[iTile]=(gRandom->Rndm()>0.5) ? -2*twopi/12.:-4*twopi/12.;//tile 12 correct -60 or -120
if (iTile==12)bbc_phi[iTile]=-3*twopi/12.;//tile 13 correct -90
if (iTile==13)bbc_phi[iTile]=-5*twopi/12.;//tile 14 correct -150
if (iTile==14)bbc_phi[iTile]=twopi/2.0;//tile 15 correct 0
if (iTile==15)bbc_phi[iTile]=5*twopi/12.;//tile 16 correct 150
//if(bbc_phi[iTile]<0.0) bbc_phi[iTile] +=2*3.14;
//if(bbc_phi[iTile]>2*3.14) bbc_phi[iTile] -=2*3.14;
}

//for (int iTile= 0; iTile<16; iTile++) cout << " bbc_phi[iTile] " << bbc_phi[iTile] <<endl;

   Float_t eXsum=0., eYsum=0., eWgt=0., psi_e(0);
   Float_t wXsum=0., wYsum=0., wWgt=0., psi_w(0);
  
   Float_t eXsum0=0., eYsum0=0., eWgt0=0.;
   Float_t wXsum0=0., wYsum0=0., wWgt0=0.;
   Float_t eXsum1=0., eYsum1=0., eWgt1=0.;
   Float_t wXsum1=0., wYsum1=0., wWgt1=0.;
 
 //psi angle from east BBC
 for(int iTile = 0; iTile < 16; iTile++) {
      eXsum += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
      eYsum += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
      eWgt += theBBC.adc(iTile)/bbcGainFac[0][iTile];
    }
   
    psi_e = atan2((eWgt>0.) ? eYsum/eWgt:0.,(eWgt>0.) ? eXsum/eWgt:0.);

  //if (psi_e < 0.)  psi_e += 2.*3.14; 
  
    //  cout << " psi_e = " << psi_e << endl;
  
//----------------------------------------------------------------------


  //psi angle from west BBC

     for(int iTile = 0; iTile <16; iTile++) {
      wXsum += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      wYsum += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      wWgt += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];

         }

    psi_w = atan2((wWgt>0.) ? wYsum/wWgt:0.,(wWgt>0.) ? wXsum/wWgt:0.);

 // if (psi_w < 0.) psi_w += 2.*3.14; 
 
    // cout << " psi_w = " << psi_w << endl;
  
  //psi angle from  BBC

    for(int iTile = 0; iTile < 16; iTile++) {
      eXsum += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
      eYsum += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
      eWgt += theBBC.adc(iTile)/bbcGainFac[0][iTile];
      
      wXsum += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      wYsum += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      wWgt += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      
      Int_t g	= gRandom->Rndm(1) > 0.5 ? 0:1; 
      if(g==0) {eXsum0 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eYsum0 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eWgt0 += theBBC.adc(iTile)/bbcGainFac[0][iTile];
                wXsum0 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wYsum0 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wWgt0 += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];}
      
      
       if(g==1) {eXsum1 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eYsum1 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eWgt1 += theBBC.adc(iTile)/bbcGainFac[0][iTile];
                wXsum1 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wYsum1 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wWgt1 += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];}
        
       }

 Float_t psi_f = atan2((eWgt>0. && wWgt>0.) ?  eYsum/eWgt + wYsum/wWgt:0., (eWgt>0. && wWgt>0.) ? eXsum/eWgt + wXsum/wWgt:0.);
 Float_t psi_f0 = atan2((eWgt0>0. && wWgt0>0.) ?  eYsum0/eWgt0 + wYsum0/wWgt0:0.,(eWgt0>0. && wWgt0>0.) ? eXsum0/eWgt0 + wXsum0/wWgt0:0.);
 Float_t psi_f1 = atan2((eWgt1>0. && wWgt1>0.) ?  eYsum1/eWgt1 + wYsum1/wWgt1:0.,(eWgt1>0. && wWgt1>0.) ? eXsum1/eWgt1 + wXsum1/wWgt1:0.);
 // if (psi_f < 0.) { psi_f += twopi; }
 
 //cout << " psi_f = " << psi_f << endl;
//-----------------------------------------------------------------------

//psi angle from  BBC

    for(int iTile = 0; iTile < 16; iTile++) {
      eXsum += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
      eYsum += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
      eWgt += theBBC.adc(iTile)/bbcGainFac[0][iTile];
      wXsum += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      wYsum += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      wWgt += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
      
      Int_t g	= gRandom->Rndm(1) > 0.5 ? 0:1; 
      if(g==0) {eXsum0 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eYsum0 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eWgt0 += theBBC.adc(iTile)/bbcGainFac[0][iTile];
                wXsum0 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wYsum0 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wWgt0 += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];}
      
      
       if(g==1) {eXsum1 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eYsum1 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile)/bbcGainFac[0][iTile];
                eWgt1 += theBBC.adc(iTile)/bbcGainFac[0][iTile];
                wXsum1 += cos(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wYsum1 += sin(2*bbc_phi[iTile])*theBBC.adc(iTile+24)/bbcGainFac[1][iTile];
                wWgt1 += theBBC.adc(iTile+24)/bbcGainFac[1][iTile];}
      
      

       }

 Float_t psi_fn  = atan2((eWgt>0. && wWgt>0.) ?  eYsum/eWgt - wYsum/wWgt:0., (eWgt>0. && wWgt>0.) ? eXsum/eWgt - wXsum/wWgt:0.);
 Float_t psi_fn0 = atan2((eWgt0>0. && wWgt0>0.) ?  eYsum0/eWgt0 - wYsum0/wWgt0:0.,(eWgt0>0. && wWgt0>0.) ? eXsum0/eWgt0 - wXsum0/wWgt0:0.);
 Float_t psi_fn1 = atan2((eWgt1>0. && wWgt1>0.) ?  eYsum1/eWgt1 - wYsum1/wWgt1:0.,(eWgt1>0. && wWgt1>0.) ? eXsum1/eWgt1 - wXsum1/wWgt1:0.);


 // if (psi_f < 0.) { psi_f += twopi; }
 
 //cout << " psi_fn = " << psi_fn << endl;



mTupleBBC->Fill(psi_e,psi_w,psi_f,mBBCAdcCut,psi_fn,psi_f0,psi_f1,psi_fn0,psi_fn1,mBBCAdcCut0);
/*
Float_t StFlowEvent::BBC_GetPhi(int iTile) {
  //get phi of BBC tile 
//Pmt  PmtNo[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};


//   Float_t  phi[24] = {1.5708,0.523599,5.75959,4.71239,3.66519,2.61799,1.5708,1.5708,0.523599,0,5.75959,4.71239,4.71239,3.66519,3.14159,2.61799,
//                     1.5708,0,4.71239,3.14159,1.5708,0,4.71239,3.14159};

  Float_t phi[16] = {1.5708,0.5236,5.75959,4.71239,3.66519,2.61799,1.5708,1.5708,0.523599,0.0,5.75959,4.71239,4.71239,3.66519,3.14159,2.61799};

  return phi[iTile];
  }
*/
  
  //  cout<<"::ctbsum::"<<ctbsum<<endl;            // CTB
  //  cout<<"::theZdc.adcSum::"<<theZdc.adcSum()<<endl;   // ZDC
 //----------------------------------------------------------------------------------------------------------------------
                                                            getTrgType();
							    
                UInt_t T(0);if(mMinb &&(!mCentral) &&(!mHt1) &&(!mHt2)) T=1;if(mHt1&&(!mCentral) &&(!mMinb) &&(!mHt2)) T=2;			if(mHt2 &&(!mCentral) &&(!mMinb) &&(!mHt1)) T=3;
		if(mCentral &&(!mHt2) &&(!mMinb) &&(!mHt1)) T=4;if(mMinb && mCentral) T=5;if(mMinb && mHt1) T=6;if(mMinb && mHt2) T=7;	
		if(mHt1 && mHt2) T=8;if(mHt1 && mCentral) T=9;if(mHt2 && mCentral) T=10;if(mMinb && mCentral&&mHt1) T=11;
		if(mMinb && mCentral&&mHt2) T=12;if(mCentral&&mHt1&&mHt2) T=13;	if(mMinb && mHt1&&mHt2) T=14;if(mMinb && mCentral&&mHt1&&mHt2) T=15;	
		if(mJPsi&&mJPsi1) T=16;	if(mJPsi&&(!mJPsi1)) T=17;if(mJPsi1&&(!mJPsi)) T=18;		
		//------------------------------AuAu2007-------------------------
		if(mL2gamma1 && (!mL2gamma2) && (!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) T=20;
		if(mL2gamma2 && (!mL2gamma1) && (!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) T=21;
		if(mbht2mb1 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) T=22;
		if(mbht2mb2 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) T=23;
		if(mbht2mb3 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb4) &&(!mbht2mb5)) T=24;
		if(mbht2mb4 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb5)) T=25;
		if(mbht2mb5 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4)) T=26;
		if(mL2gamma1 && mL2gamma2) T=27;if(mL2gamma1 && mbht2mb1) T=28;	if(mL2gamma1 && mbht2mb2) T=29;
		if(mL2gamma1 && mbht2mb3) T=30;	if(mL2gamma1 && mbht2mb4) T=31;	if(mL2gamma1 && mbht2mb5) T=32;
            	if(mL2gamma2 && mbht2mb1) T=33;	if(mL2gamma2 && mbht2mb2) T=34;	if(mL2gamma2 && mbht2mb3) T=35;
		if(mL2gamma2 && mbht2mb4) T=36;	if(mL2gamma2 && mbht2mb5) T=37;	if(mbht2mb1 && mbht2mb2) T=38;		
		if(mbht2mb1 && mbht2mb3) T=39;	if(mbht2mb1 && mbht2mb4) T=40;	if(mbht2mb1 && mbht2mb5) T=41;		
		if(mbht2mb2 && mbht2mb3) T=42;	if(mbht2mb2 && mbht2mb4) T=43;	if(mbht2mb2 && mbht2mb5) T=44;
		if(mbht2mb3 && mbht2mb4) T=45;	if(mbht2mb3 && mbht2mb5) T=46;	if(mbht2mb4 && mbht2mb5) T=47;
		
		if(mHt0) T=48 ; if(mht4) T=49; if(mht4fast) T=50; if(mupsilondAu) T=51;	if(mHt0 && (!mHt1) && (!mHt2) &&(!mht4) &&(!mht4fast) &&(!mupsilondAu)) T=52;
		if((!mHt0) && mHt1 && (!mHt2) &&(!mht4) &&(!mht4fast) &&(!mupsilondAu)) T=53;if((!mHt0) && (!mHt1) && mHt2 &&(!mht4) &&(!mht4fast) &&(!mupsilondAu)) T=54;
		if((!mHt0) && (!mHt1) && (!mHt2) && mht4 &&(!mht4fast) &&(!mupsilondAu)) T=55;	if((!mHt0) && (!mHt1) && (!mHt2) && (!mht4) && mht4fast &&(!mupsilondAu)) T=56;
		if((!mHt0) && (!mHt1) && (!mHt2) && (!mht4) && (!mht4fast) && mupsilondAu) T=57;if(mHt0 && mHt1) T=58; if(mHt0 && mHt2) T=59; if(mHt0 && mht4) T=60;
		if(mHt1 && mHt2) T=61; if(mHt1 && mht4) T=62; 		if(mHt2 && mht4) T=63; if(mbht2mb6) T=64;
		  //-----------------------------L2gamma2010----------------
		if(mL2gamma1_2010) T=100;
		if(mL2gamma2_2010) T=101;
		//-----------------------------L2gamma2011----------------
		if(mL2gamma0_2011) T=200;
		if(mL2gamma1_2011) T=201;
		if(mL2gamma2_2011) T=202;
		//----------------------------------------------------------
		//______L2gamma for pp2009
                if(mL2gamma_pp2009) T=102;

		//if(!(T <16 || T > 47)) return kStOK; //dAu2008 hitower data
		//	                                		cout<<"T is"<<T<<endl; 
//------------------------------------------Get Recation Plane information------------------------------------------------------------------							

//Do BBC and smdZDC
 StSPtrVecTrackNode& trackNodes=mStEvent->trackNodes();
//  StTrack* track;
    int w(0),HT(0),F1(0),F2(0),F3(0),F4(0);      
    float XR(0),YR(0), XRS1(0), XRS2(0),YRS1(0),YRS2(0),angleRe(0),angleReS1(0),angleReS2(0); 
    int QQ(0),Q1(0),Q2(0),Q3(0),Q4(0),Q1S1(0),Q1S2(0),Q2S1(0),Q2S2(0),Q3S1(0),Q3S2(0),Q4S1(0),Q4S2(0);
   float XRf(0),YRf(0), XRfS1(0), XRfS2(0),YRfS1(0),YRfS2(0),angleRef(0),angleRefS1(0),angleRefS2(0),
      Xr(0),Yr(0),Xrf(0),Yrf(0),XrfE(0),YrfE(0),XrfW(0),YrfW(0),
     Xrp(0),Yrp(0),Xrphi(0),Yrphi(0),Xrfp(0),Yrfp(0),Xrfphi(0),Yrfphi(0),
    XrfpW(0),YrfpW(0),XrfphiW(0),YrfphiW(0),XrfpE(0),YrfpE(0),XrfphiE(0),YrfphiE(0),
       XRfW(0),YRfW(0),XRfE(0),YRfE(0),XRfSW1(0),YRfSW1(0),XRfSW2(0),YRfSW2(0),XRfSE1(0),YRfSE1(0),XRfSE2(0),YRfSE2(0),
    XRP(0),YRP(0), XRPS1(0), XRPS2(0),YRPS1(0),YRPS2(0),angleReP(0),angleReS1P(0),angleReS2P(0), 
   XRPHI(0),YRPHI(0), XRPHIS1(0), XRPHIS2(0),YRPHIS1(0),YRPHIS2(0),angleRePHI(0),angleReS1PHI(0),angleReS2PHI(0), 
    XRFP(0),YRFP(0), XRFPS1(0), XRFPS2(0),YRFPS1(0),YRFPS2(0),angleRePf(0),angleReS1FP(0),angleReS2FP(0),
   XRFPHI(0),YRFPHI(0), XRFPHIS1(0), XRFPHIS2(0),YRFPHIS1(0),YRFPHIS2(0),angleRePHIf(0),angleReS1PHIf(0),angleReS2PHIf(0), 
     XRFPW(0),YRFPW(0), XRFPS1W(0), XRFPS2W(0),YRFPS1W(0),YRFPS2W(0),
   XRFPHIW(0),YRFPHIW(0), XRFPHIS1W(0), XRFPHIS2W(0),YRFPHIS1W(0),YRFPHIS2W(0), 
     XRFPE(0),YRFPE(0), XRFPS1E(0), XRFPS2E(0),YRFPS1E(0),YRFPS2E(0), 
   XRFPHIE(0),YRFPHIE(0), XRFPHIS1E(0), XRFPHIS2E(0),YRFPHIS1E(0),YRFPHIS2E(0);
   
   int F1pe(0),Q1pe(0),Q1S1pe(0),Q1S2pe(0);
    float Xrpe(0),Yrpe(0),XRpe(0),YRpe(0),XRS1pe(0),YRS1pe(0),XRS2pe(0),YRS2pe(0),angleReS1pe(0),angleReS2pe(0),angleRepe(0);
    
    int F1ne(0),Q1ne(0),Q1S1ne(0),Q1S2ne(0);
    float Xrne(0),Yrne(0),XRne(0),YRne(0),XRS1ne(0),YRS1ne(0),XRS2ne(0),YRS2ne(0),angleReS1ne(0),angleReS2ne(0),angleRene(0);
    
    int F1nc(0),Q1nc(0),Q1S1nc(0),Q1S2nc(0);
    float Xrnc(0),Yrnc(0),XRnc(0),YRnc(0),XRS1nc(0),YRS1nc(0),XRS2nc(0),YRS2nc(0),angleReS1nc(0),angleReS2nc(0),angleRenc(0);
    
    int F1pc(0),Q1pc(0),Q1S1pc(0),Q1S2pc(0);
    float Xrpc(0),Yrpc(0),XRpc(0),YRpc(0),XRS1pc(0),YRS1pc(0),XRS2pc(0),YRS2pc(0),angleReS1pc(0),angleReS2pc(0),angleRepc(0); 
   
    float Xrpt(0),Yrpt(0),XRpt(0),YRpt(0),XRS1pt(0),YRS1pt(0),XRS2pt(0),YRS2pt(0),angleRept(0),angleReS1pt(0),angleReS2pt(0);
 float PT(0),PT1(0),PT2(0);
 float angleRefW(0),angleRefE(0),angleRefW1(0),angleRefE1(0),angleRefW2(0),angleRefE2(0);
 float XRfW1(0),YRfW1(0),XRfE1(0),YRfE1(0),XRfW2(0),YRfW2(0),XRfE2(0),YRfE2(0);
 
  for(size_t nodeIndex=0;nodeIndex<trackNodes.size();nodeIndex++)
  {    size_t numberOfTracksInNode=trackNodes[nodeIndex]->entries(primary); 
 
    for( size_t trackIndex=0;trackIndex<numberOfTracksInNode;trackIndex++)
    {      track=trackNodes[nodeIndex]->track(primary,trackIndex);
               //new
	       if (!(track && track->flag()>=0)) continue;
              float nFP(0.),nFPW(0.),nFPE(0.),nFPT(0.),Ppo(0.);
	       nFP=track->fitTraits().numberOfFitPoints(kTpcId);//Int_t Ppo=track->numberOfPossiblePoints();
	        nFPT=track->fitTraits().numberOfFitPoints();
	        nFPW=track->fitTraits().numberOfFitPoints(kFtpcWestId);
		 nFPE=track->fitTraits().numberOfFitPoints(kFtpcEastId);
	       Ppo=track->numberOfPossiblePoints();
	       float dca=track->impactParameter(); 
	      
	      //if(track->fitTraits().numberOfFitPoints(kTpcId)<15) continue;      
	       StTrackGeometry* geometry=track->geometry();       // float p=geometry->momentum().mag();
       float Eta=geometry->momentum().pseudoRapidity(); float epsi=geometry->momentum().phi();       
       float pt=geometry->momentum().perp();short ch=geometry->charge();
      if(!(pt>0.1 && pt<=2)) continue;
       w=numberOfTracksInNode;      if(w==0) return kStOK;          
       HT=1; QQ+=HT;    Int_t g	= gRandom->Rndm(1) > 0.5 ? 0:1; 
     
       //     cout << " Eta =" <<  Eta << " nFP= "<<  nFP<<" nFP/Ppo=" <<  nFP/Ppo << " dca= " <<dca <<" pt = "<<pt<< " Ppo =" <<Ppo<<endl;
  
  //Save it in other loops, take care of epsi sign
   //------------------------------- Full TPC
     if(Eta<=1. && Eta>=-1.&&nFP>=16&&nFP/Ppo>0.52&&dca<2.0&&pt>0.15) 
     
       { //      cout << "here is TPC"<< endl; 
       F1=1; Xr=cos(2*epsi);  Yr=sin(2*epsi);  XR+=Xr; YR+=Yr;   Q1+=F1;
       Xrpt=pt*cos(2*epsi);  Yrpt=pt*sin(2*epsi);  XRpt+=Xrpt; YRpt+=Yrpt;PT+=pt;
   //subEvent
   if(g==0)  {XRS1+=Xr; YRS1+=Yr;  XRS1pt+=Xrpt; YRS1pt+=Yrpt; Q1S1+=F1;PT1+=pt;}      
   if(g==1)  {XRS2+=Xr; YRS2+=Yr;  XRS2pt+=Xrpt; YRS2pt+=Yrpt; Q1S2+=F1;PT2+=pt;}    
   
    angleRe=atan2(YR,XR);     if(PT>0)angleRept=atan2(YRpt/PT,XRpt/PT); 
    
    if(g==0) angleReS1=atan2(YRS1,XRS1);     if(g==1) angleReS2=atan2(YRS2,XRS2);
    if(g==0&&PT1>0) angleReS1pt=atan2(YRS1pt/PT1,XRS1pt/PT1);     if(g==1&&PT2>0) angleReS2pt=atan2(YRS2pt/PT2,XRS2pt/PT2);
    }
    //-------------------------------
   
    //---------------------------------------TPC+ve eta
    if(Eta<=1. && Eta>=0.&&nFP>=16&&nFP/Ppo>0.52&&dca<2.0&&pt>0.15) 
     
      { //      cout << "here is TPC"<< endl; 
F1pe=1; Xrpe=cos(2*epsi);  Yrpe=sin(2*epsi);  XRpe+=Xrpe; YRpe+=Yrpe;   Q1pe+=F1pe;
       
   //subEvent
   if(g==0)  {XRS1pe+=Xrpe; YRS1pe+=Yrpe;   Q1S1pe+=F1pe;}      if(g==1)  {XRS2pe+=Xrpe; YRS2pe+=Yrpe;   Q1S2pe+=F1pe;}     
   angleRepe=atan2(YRpe,XRpe);   
   
    if(g==0) angleReS1pe=atan2(YRS1pe,XRS1pe);     if(g==1) angleReS2pe=atan2(YRS2pe,XRS2pe);
    
    }
    
       //-------------------------------
    
    //--------------------------------------TPC-ve eta
     if(Eta>=-1. && Eta<=0.&&nFP>=16&&nFP/Ppo>0.52&&dca<2.0&&pt>0.15) 
     
       { //      cout << "here is TPC"<< endl;
       F1ne=1; Xrne=cos(2*epsi);  Yrne=sin(2*epsi);  XRne+=Xr; YRne+=Yrne;   Q1ne+=F1;
       
   //subEvent
   if(g==0)  {XRS1ne+=Xrne; YRS1ne+=Yrne;   Q1S1ne+=F1ne;}      if(g==1)  {XRS2ne+=Xrne; YRS2ne+=Yrne;   Q1S2ne+=F1ne;}  
      angleRene=atan2(YRne,XRne);   
   
    if(g==0) angleReS1ne=atan2(YRS1ne,XRS1ne);     if(g==1) angleReS2ne=atan2(YRS2ne,XRS2ne);
    
    }
     //-------------------------------
    
    //------------------------------------TPC+ve charge
    if(Eta<=1. && Eta>=-1.&&nFP>=16&&nFP/Ppo>0.52&&dca<2.0&&pt>0.15&&ch>0) 
     
      { //      cout << "here is TPC"<< endl;
       F1pc=1; Xrpc=cos(2*epsi);  Yrpc=sin(2*epsi);  XRpc+=Xrpc; YRpc+=Yrpc;   Q1pc+=F1pc;
       
   //subEvent
   if(g==0)  {XRS1pc+=Xrpc; YRS1pc+=Yrpc;   Q1S1pc+=F1pc;}      if(g==1)  {XRS2pc+=Xrpc; YRS2pc+=Yrpc;   Q1S2pc+=F1pc;}     
   angleRepc=atan2(YRpc,XRpc);   
   
    if(g==0) angleReS1pc=atan2(YRS1pc,XRS1pc);     if(g==1) angleReS2pc=atan2(YRS2pc,XRS2pc);
    
    }
     //-------------------------------
     
    //----------------------------------TPC-ve charge
    if(Eta<=1. && Eta>=-1.&&nFP>=16&&nFP/Ppo>0.52&&dca<2.0&&pt>0.15&&ch<0) 
     
      { //      cout << "here is TPC"<< endl;
       F1nc=1; Xrnc=cos(2*epsi);  Yrnc=sin(2*epsi);  XRnc+=Xrnc; YRnc+=Yrnc;   Q1nc+=F1nc;
       
   //subEvent
   if(g==0)  {XRS1nc+=Xrnc; YRS1nc+=Yrnc;   Q1S1nc+=F1nc;}      if(g==1)  {XRS2nc+=Xrnc; YRS2nc+=Yrnc;   Q1S2nc+=F1nc;}     
   angleRenc=atan2(YRnc,XRnc);   
   
    if(g==0) angleReS1nc=atan2(YRS1nc,XRS1nc);     if(g==1) angleReS2nc=atan2(YRS2nc,XRS2nc);
    
    }
    
    //-------------------------------
       //FTpc
     if(((Eta>=2.5 && Eta<=4.0)||(Eta<=-2.5&&Eta>=-4.0))&& (nFPW>=5||nFPE>=5)&&(nFPW<11||nFPE<11)&&
     ((nFPW/Ppo>0.52&&nFPW/Ppo<1.05)||(nFPE/Ppo>0.52&&nFPE/Ppo<1.05))&&dca<3)
     { 
     
       //cout << "here is FTPC"<< endl;
     F2=1; Xrf=cos(2*epsi); Yrf=sin(2*epsi);       
     
     XRf+=Xrf; YRf+=Yrf;  Q2+=F2;
   
   //subEvent
   if(g==0)  {XRfS1+=Xrf;  YRfS1+=Yrf;   Q2S1+=F2;}      if(g==1)  {XRfS2+=Xrf;  YRfS2+=Yrf;   Q2S2+=F2;}	   angleRef=atan2(YRf,XRf);    
	   
	   if(g==0)  angleRefS1=atan2(YRfS1,XRfS1);   	   if(g==1)  angleRefS2=atan2(YRfS2,XRfS2);
	   
     }
     
   
    //FTPC West
    if((Eta>=2.5 && Eta<=4.0&& nFPW>=5&&nFPW<11&&nFPW/Ppo>0.52&&nFPW/Ppo<1.05&&dca<3))
    { 
      //cout << "here is FTPC west"<< endl; 
    F3=1; XrfW=cos(2*epsi);   YrfW=sin(2*epsi);     
     
    XRfW+=XrfW; YRfW+=YrfW; Q3+=F3;
   angleRefW=atan2(YRfW,XRfW);
   //subEvent
   if(g==0)  {XRfW1+=XrfW; YRfW1+=YrfW;  Q3S1+=F3;}     if(g==1)  {XRfW2+=XrfW; YRfW2+=YrfW;  Q3S2+=F3;} 
   if(g==0)  angleRefW1=atan2(YRfW1,XRfW1);    if(g==1)   angleRefW2=atan2(YRfW2,XRfW2); 
   
   }
     
      //FTPC EAST
     if((Eta<=-2.5&&Eta>=-4.0&& nFPE>=5&&nFPE<11&&nFPE/Ppo>0.52&&nFPE/Ppo<1.05&&dca<3))
     { 
       //cout << "here is FTPC east"<< endl;
     F4=1; XrfE=cos(2*epsi); YrfE=sin(2*epsi);     
           
     XRfE+=XrfE; YRfE+=YrfE; Q4+=F4;
   angleRefE=atan2(YRfE,XRfE);
   //subEvent
   if(g==0)  {XRfE1+=XrfE; YRfE1+=YrfE;   Q4S1+=F4;}      if(g==1)  {XRfE2+=XrfE; YRfE2+=YrfE;   Q4S2+=F4;}
   if(g==0)  angleRefE1=atan2(YRfE1,XRfE1);    if(g==1)   angleRefE2=atan2(YRfE2,XRfE2); 
     
     }
   
   
   }
 }
  

  //StSPtrVecTrackNode& trackNodes=mStEvent->trackNodes();
  // StTrack* track;
    int NumOfTracks(0);
  //--------------------------------------------------Track QA--------------------------------
  Int_t nPrimary2= 0;
  nPrimary2 = mMuDstMaker->muDst()->primaryTracks()->GetEntries();

  cout<<" No. of Primary trk: "<<nPrimary2<<endl;
  
  for ( Int_t itrk = 0; itrk < nPrimary2; itrk++ )
    {
      float pTPC(0), ptTPC(0),dEdx(0),EtaTr(0),epsi(0); short ch(0);
      ch=mMuDstMaker->muDst()->primaryTracks(itrk)->charge();        
      pTPC=mMuDstMaker->muDst()->primaryTracks(itrk)->p().mag();  
      ptTPC=mMuDstMaker->muDst()->primaryTracks(itrk)->pt(); 

      EtaTr=mMuDstMaker->muDst()->primaryTracks(itrk)->eta();
      epsi=mMuDstMaker->muDst()->primaryTracks(itrk)->phi();  
      Int_t Fp=mMuDstMaker->muDst()->primaryTracks(itrk)->nHitsFit();
      
      float xg(0); float yg(0); float zg(0);   
      Int_t Ppo=mMuDstMaker->muDst()->primaryTracks(itrk)->nHitsPoss();

      //      if(!(ptTPC<50.&&ptTPC>1.2)) continue;
      //      if(!((Fp/Ppo) > 0.52) ) continue;
      //      if(!(Fp>15)) continue; 

      Int_t FpT(0); Int_t FpW(0); Int_t FpE(0);


      
      float       px = pTPC*cos( epsi );
      float       py = pTPC*sin( epsi );
      float       pz = sqrt( pTPC*pTPC - ptTPC*ptTPC );
      int flag(0);
      
      float ltrkV=mMuDstMaker->muDst()->primaryTracks(itrk)->length(); //Returns length of track (cm) from primary vertex to last measured point. 
      float lpnt=mMuDstMaker->muDst()->primaryTracks(itrk)->lengthMeasured(); //Returns length of track (cm) from first to last measured point. 
      
      UInt_t RefM=MPos+MNeg;   
      
      
      float dca=mMuDstMaker->muDst()->primaryTracks(itrk)->dcaGlobal().mag(); //Returns 3D distance of closest approach to primary vertex of associated global track. 
      float pdca=mMuDstMaker->muDst()->primaryTracks(itrk)->dca().mag();  //Returns 3D distance of closest approach to primary vertex. 
      float dcaD=mMuDstMaker->muDst()->primaryTracks(itrk)->dcaD(); //Signed radial component of global DCA (projected). 
      float dcaZ=mMuDstMaker->muDst()->primaryTracks(itrk)->dcaZ(); //Z component of global DCA. 

      
      //if(dca>3) continue;
      dEdx=mMuDstMaker->muDst()->primaryTracks(itrk)->dEdx();  
      
      // cout<<"::dEdx::"<<dEdx<<"::ptTPC::"<<ptTPC<<"::NumOfTracks:"<<NumOfTracks<<endl;
      float v2chT(0),v2chF(0), v2chB(0); 
      v2chT=cos(2*(epsi)-angleRe); v2chF=cos(2*(epsi)-angleRef);
      v2chB=cos(2*(epsi)-psi_f);
      
      mTupleTrackS->Fill(pTPC,ptTPC,dEdx,ch,Ppo,epsi,flag,RefM,runId,z,evt,px,py,pz);
      mTupleTrackSS->Fill(EtaTr,Fp,T,dca,rank,rank1,nV,muRef,psi_e,psi_w,psi_fn,psi_fn0,psi_fn1);
      //      mTupleTrackSSS->Fill(FpT,YY,XX,angleRe,angleRef,psi_f,v2chT,v2chF,v2chB,mBBCAdcCut,mBBCAdcCut0,psi_f0,psi_f1);
      mTupleTrackS4->Fill(angleRefS1,angleRefS2,angleReS1,angleReS2,angleRept,angleReS1pt,angleReS2pt,ltrkV,lpnt,pdca,dcaD,dcaZ);
      //      mTupleTrackS5->Fill(angleRepe,angleReS1pe,angleReS2pe,angleRene,angleReS1ne,angleReS2ne,angleRenc,angleReS1nc,angleReS2nc,angleRepc,angleReS1pc,angleReS2pc);
      //      mTupleFTPCTr->Fill(angleRefW,angleRefE,angleRefW1,angleRefE1,angleRefW2,angleRefE2,vpdz,BBCz,ZDCz);
    }

  //________#######################  TRACK INfo from "OLd StThridMaker" (StTrack)##############
 int NumOfTracks1(0);
//--------------------------------------------------Track QA--------------------------------
 
 for(size_t nodeIndex=0;nodeIndex<trackNodes.size();nodeIndex++)
  {float pTPC(0), ptTPC(0),dEdx(0),EtaTr(0),epsi(0); short ch(0);
   size_t numberOfTracksInNode=trackNodes[nodeIndex]->entries(primary);


  cout<<" No. of Primary trk (stTrack): "<<numberOfTracksInNode<<endl;

    for( size_t trackIndex=0;trackIndex<numberOfTracksInNode;trackIndex++)
  {   track=trackNodes[nodeIndex]->track(primary,trackIndex); 
       if (!(track && track->flag()>=0)) continue;
       StTrackGeometry* geometry=track->geometry();        ch=geometry->charge();        
       pTPC=geometry->momentum().mag();         ptTPC=geometry->momentum().perp(); 
       NumOfTracks1++;

       //      if(!(ptTPC<50.&&ptTPC>.1)) continue;
       //      if(!(ptTPC<50.&&ptTPC>1.2)) continue;
      // if(!(ptTPC<50.)) continue; run fine w/ minbias 
       //if(!(ptTPC>=.1&&ptTPC<2)) continue; 
        EtaTr=geometry->momentum().pseudoRapidity();epsi=geometry->momentum().phi();  
       Int_t Fp=track->fitTraits().numberOfFitPoints(kTpcId);    
       //if(Fp<20) continue;//kFtpcWestId//kFtpcEastId
       Int_t FpT=track->fitTraits().numberOfFitPoints();
       Int_t FpW=track->fitTraits().numberOfFitPoints(kFtpcWestId);
       Int_t FpE=track->fitTraits().numberOfFitPoints(kFtpcEastId);
                     

       Int_t Ppo=track->numberOfPossiblePoints(); float xg=geometry->origin().x();
       float yg=geometry->origin().y();float zg=geometry->origin().z(); float px=geometry->momentum().x();
       float py=geometry->momentum().y();float pz=geometry->momentum().z(); float lg=track->length();
      int flag= track->flag();
        UInt_t RefM=MPos+MNeg;   
       StPtrVecTrackPidTraits traits=track->pidTraits(kTpcId); float dca=track->impactParameter(); 
       //if(dca>3) continue;
       StDedxPidTraits* pid; 
          
       for(UInt_t V=0;V<traits.size();V++)        pid=dynamic_cast<StDedxPidTraits*>(traits[V]);
          if(pid&&pid->method()==kTruncatedMeanId)           dEdx=pid->mean();  
// cout<<"::dEdx::"<<dEdx<<"::ptTPC::"<<ptTPC<<"::NumOfTracks:"<<NumOfTracks<<endl;
               float v2chT(0),v2chF(0), v2chB(0); v2chT=cos(2*(epsi)-angleRe); v2chF=cos(2*(epsi)-angleRef);
	       v2chB=cos(2*(epsi)-psi_f);
	      		
	       mTupleTrackSO->Fill(pTPC,ptTPC,dEdx,ch,Ppo,epsi,flag,RefM,runId,z,evt,px,py,pz);
	       mTupleTrackSSO->Fill(EtaTr,Fp,T,dca,rank,rank1,nV,muRef,psi_e,psi_w,psi_fn,psi_fn0,psi_fn1);
	       mTupleTrackSSSO->Fill(FpT,YY,XX,angleRe,angleRef,psi_f,v2chT,v2chF,v2chB,mBBCAdcCut,mBBCAdcCut0,psi_f0,psi_f1);
	       //	       mTupleTrackS4O->Fill(angleRefS1,angleRefS2,angleReS1,angleReS2,angleRept,angleReS1pt,angleReS2pt);
	       //	       mTupleTrackS5O->Fill(angleRepe,angleReS1pe,angleReS2pe,angleRene,angleReS1ne,angleReS2ne,angleRenc,angleReS1nc,angleReS2nc,angleRepc,angleReS1pc,angleReS2pc);
	       //	       mTupleFTPCTrO->Fill(angleRefW,angleRefE,angleRefW1,angleRefE1,angleRefW2,angleRefE2,vpdz,BBCz,ZDCz);
  }
  }


 mTupleQATrkInfo->Fill(NumOfTracks1,nPrimary2,z,rank,rank1,nV);   
 /////###################### Old info keeping end##############

   
  //   cout<<"::NumOfTracks:"<<NumOfTracks<<endl;	
//---------------------------------------------------------------------------------------------
 for (Int_t Iu2 = 1;Iu2 <= 4800;Iu2++) mStatusT[Iu2]=0; for (Int_t iu = 0;iu < 4800;iu++) {mCalibPedValues[iu] = 0;mCalibPedRMSValues[iu] = 0;}
 St_db_Maker *dbMk = (St_db_Maker*)GetMaker("db");if (dbMk == 0) return kStOK;   TDataSet *db = NULL;  St_emcStatus* run2=NULL;
  TString DbName = "Calibrations/emc/y3bemc";db = dbMk->GetDataBase(DbName); TString tableName2 = "bemcStatus";TString tableName = "bemcPed";//bprs
 run2 = (St_emcStatus*)db->Find(tableName2.Data());emcStatus_st *runst2 = run2->GetTable();
	St_emcPed *ped = (St_emcPed*)db->Find(tableName.Data());	emcPed_st *pedst = ped->GetTable();
		int Iu3(0);    for (Iu3 = 1;Iu3 <=4800;Iu3++) { mStatusT[Iu3-1] = runst2[0].Status[Iu3-1];}
		 int iu4(0);  for (iu4 = 0;iu4 <4800;iu4++) { mCalibPedValues[iu4] = ((Float_t)pedst[0].AdcPedestal[iu4])/100.0;}
		 for (iu4 = 0;iu4 <4800;iu4++) { mCalibPedRMSValues[iu4] = ((Float_t)pedst[0].AdcPedestalRMS[iu4])/100.0;}
	
   StEmcCollection* emcCollection = mStEvent->emcCollection();    if(!emcCollection)//return kStOK;
    {
      //cout<<"There is no emcCollection"<<endl;      
 return kStWarn;}
       
   StEmcGeom * mGeom0=StEmcGeom::instance("bemc"); StEmcGeom * mGeom2=StEmcGeom::instance("bsmde");     StEmcGeom * mGeom3=StEmcGeom::instance("bsmdp");
   //   associateTracksWithTowerHits();    associateTracksWithSmdeHits();   associateTracksWithSmdpHits();   associateTracksWithTowerHitsAll();
 
 //----------------------------------------
 
 TowerHVcahngeList();  
 
 //------------------------------------------------------------------------------------------------------------------------------//  
  
  
  
  
  
  
   StEmcDetector* detector1 = emcCollection->detector(kBarrelEmcTowerId);  if (!detector1)   {
     //cout<<"There is no kBarrelEmcTowerId Detector"<<endl; 
return kStWarn;}
 
 //----------------------------------------------------------------Before Track projection------------------------
  float hiTowerB(0);
   float EnergyhB(0);  Int_t SubhB(0);  UInt_t ModulehB(0), EtahB(0),AdchB(0); int TC5(0),IdhB(0),huB(0); float eThB(0),phiThB(0),ThetahB(0),hpedB(0),hRMSB(0);
  for(UInt_t mw1A=1;mw1A<=120;mw1A++)   
  {   StEmcModule* module = detector1->module(mw1A);  StSPtrVecEmcRawHit& hitsEw1A=module->hits(); 
    for(unsigned int rw1A=0;rw1A<hitsEw1A.size();rw1A++) if(hitsEw1A[rw1A])
   { EnergyhB=hitsEw1A[rw1A]->energy(); if(EnergyhB>hiTowerB) hiTowerB=EnergyhB;  

if(EnergyhB==hiTowerB) SubhB=abs(hitsEw1A[rw1A]->sub()); if(EnergyhB==hiTowerB) ModulehB=hitsEw1A[rw1A]->module(); 
if(EnergyhB==hiTowerB) EtahB=hitsEw1A[rw1A]->eta(); if(EnergyhB==hiTowerB) AdchB=hitsEw1A[rw1A]->adc(); 
if(EnergyhB==hiTowerB) mGeom0->getId(ModulehB,EtahB,SubhB,IdhB);if(EnergyhB==hiTowerB) mGeom0->getEtaPhi(IdhB,eThB,phiThB);
if(EnergyhB==hiTowerB) mGeom0->getTheta(ModulehB,EtahB,ThetahB);if(EnergyhB==hiTowerB) huB=mStatusT[IdhB-1];
if(EnergyhB==hiTowerB) hpedB=mCalibPedValues[IdhB-1];if(EnergyhB==hiTowerB) hRMSB=mCalibPedRMSValues[IdhB-1];
if(EnergyhB==hiTowerB) {for (Int_t ku=0; ku<TowerHVchangeMax; ku++){ if (IdhB == HVTowerid[ku]) TC5=1;}}  
   }
 }
mTupleFive->Fill(hiTowerB,SubhB,ModulehB,EtahB,AdchB,T,huB,hpedB,runId,IdhB,eThB,phiThB,ThetahB,hRMSB);
 //---------------------------------------------------------------
//-------------------------------------Project tracks------------------------------------------------------
for(int RW=0; RW<4800;RW++) ParentTrackP[RW]=0.0;
for(int RW2=0; RW2<4800;RW2++) ParticleId[RW2]=0.0;
for(int RW3=0; RW3<4800;RW3++)   ParentTrackPOSe[RW3]=0.0;
for(int RW4=0; RW4<4800;RW4++) 	ParentTrackPOSP[RW4]=0.0;

 associateHadronsWithBEMC1();associateHadronsWithBEMC2();associateHadronsWithBEMC3(); 

 //-------------------------------------------------------------------------
  for(UInt_t mw=1;mw<=120;mw++)   {   StEmcModule* module = detector1->module(mw);  StSPtrVecEmcRawHit& hitsEw=module->hits(); 
    for(unsigned int rw=0;rw<hitsEw.size();rw++) if(hitsEw[rw])
   {  float Energy(0);  Int_t Sub(0);  UInt_t Module(0), Eta(0),Adc(0); int TC4(0),Id(0),hu4(0); float TPe(0),TPp(0),PPF(0),pidF(0),eT(0),phiT(0),Theta(0),hped4(0),hRMS4(0);
 Energy=hitsEw[rw]->energy(); if(Energy<4) continue; Sub=abs(hitsEw[rw]->sub()); Module=hitsEw[rw]->module(); Eta=hitsEw[rw]->eta();   
 Adc=hitsEw[rw]->adc(); mGeom0->getId(Module,Eta,Sub,Id);mGeom0->getEtaPhi(Id,eT,phiT);mGeom0->getTheta(Module,Eta,Theta);
 hu4=mStatusT[Id-1];hped4=mCalibPedValues[Id-1];hRMS4=mCalibPedRMSValues[Id-1]; 
 for (Int_t ku=0; ku<TowerHVchangeMax; ku++) {if (Id == HVTowerid[ku]) TC4=1;} 
 PPF=ParentTrackP[Id-1];pidF=ParticleId[Id-1];TPe=ParentTrackPOSe[Id-1];TPp=ParentTrackPOSP[Id-1];
  if(Energy>=6) mTupleFour->Fill(Energy,Sub,Module,Eta,Adc,T,evt,z,runId,MNeg,MPos,Id,phiT,Theta);
  if(Energy>=6) mTupleFour2->Fill(hu4,hped4,hRMS4,eT,TC4,PPF,pidF,TPe,TPp);
  //if((Adc-hped4-1.5*hRMS4)>24) mTupleFour3->Fill(hu4,hped4,hRMS4,eT,Id,Adc,eT,Theta,phiT,z,T,Energy,MPos,MNeg);
  }
 }
 
//---------------------------------------After Track Projection-----------------------------------
 //high tower
 float hiTower(0); float highestEt(0);
   float Energyh(0),hped(0),hRMS(0);  Int_t hu(0),Subh(0);  UInt_t Moduleh(0), Etah(0),Adch(0); int TC6(0),Idh(0); float TPeh(0),TPph(0),eTh(0),phiTh(0),Thetah(0),PPTPC(0),pidh(0);
    for(UInt_t mw1=1;mw1<=120;mw1++)   {   StEmcModule* module = detector1->module(mw1);  StSPtrVecEmcRawHit& hitsEw1=module->hits(); 
    for(unsigned int rw1=0;rw1<hitsEw1.size();rw1++) if(hitsEw1[rw1])
   { Energyh=hitsEw1[rw1]->energy(); if(Energyh>hiTower) hiTower=Energyh;  
if(Energyh==hiTower) Subh=abs(hitsEw1[rw1]->sub()); if(Energyh==hiTower) Moduleh=hitsEw1[rw1]->module(); 
if(Energyh==hiTower) Etah=hitsEw1[rw1]->eta();  if(Energyh==hiTower) Adch=hitsEw1[rw1]->adc(); 
if(Energyh==hiTower) mGeom0->getId(Moduleh,Etah,Subh,Idh);if(Energyh==hiTower) mGeom0->getEtaPhi(Idh,eTh,phiTh);
if(Energyh==hiTower) mGeom0->getTheta(Moduleh,Etah,Thetah);
if(Energyh==hiTower)PPTPC=ParentTrackP[Idh-1];if(Energyh==hiTower)pidh=ParticleId[Idh-1];
if(Energyh==hiTower)TPeh=ParentTrackPOSe[Idh-1];if(Energyh==hiTower)TPph=ParentTrackPOSP[Idh-1];
 highestEt=hiTower*sin(Thetah);
if(Energyh==hiTower)hu=mStatusT[Idh-1];if(Energyh==hiTower)hped=mCalibPedValues[Idh-1];if(Energyh==hiTower)hRMS=mCalibPedRMSValues[Idh-1];
if(Energyh==hiTower){for (Int_t ku=0; ku<TowerHVchangeMax; ku++) {if (Idh == HVTowerid[ku]) TC6=1;}}
}
}
    //cout<<"highest tower::"<<Idh<<":eta:"<<eTh<<":phi is:"<<phiTh<<":Adch is:"<<Adch<<":hiTower:"<<hiTower<<"::Et::"<<hiTower*sin(Thetah)<<endl;
    //cout<<"::PPTPC::"<<PPTPC<<"::pidh::"<<pidh<<"::status is::"<<hu<<"::pedh is::"<<hped<<"::hRMS is::"<<hRMS<<endl;

mTupleFiveA->Fill(hiTower,Subh,Moduleh,Etah,Adch,T,evt,z,runId,Idh,eTh,phiTh,Thetah,PPTPC);

mTupleFiveteen->Fill(angleRe,angleReS1,angleReS2,angleRef,angleRefS1,angleRefS2,NumOfTracks,Q1,Q2,XR,YR,XRf,YRf);
mTupleSixteen->Fill(XRS1,XRS2,YRS1,YRS2,QQ,Q3,Q4,Q1S1,Q1S2,Q2S1,Q2S2,Q3S1,Q3S2,Q4S1);
mTupleSeventeen->Fill(Q4S2,XRfS1,XRfS2,YRfS1,YRfS2,XRfW,YRfW,XRfE,YRfE,XRfSW1,YRfSW1,XRfSW2,YRfSW2,XRfSE1);
mTupleEighteen->Fill(YRfSE1,XRfSE2,YRfSE2,angleRept,angleReS1pt,angleReS2pt,angleRepe,angleReS1pe,angleReS2pe,angleRene,angleReS1ne,angleReS2ne,angleRenc,angleReS1nc);
mTupleFTPCG->Fill(angleRefW,angleRefE,angleRefW1,angleRefE1,angleRefW2,angleRefE2);

    
      



mTupleTri->Fill(T,evt,z,runId,NOP,XX,YY,MPos,MNeg,allGlob,allPrim,ctbsum,zdcsum);	
mTupleTri2->Fill(muRef,nV,rank1,angleReS2nc,angleRepc,angleReS1pc,angleReS2pc,vpdz,BBCz,ZDCz);
 //add the TPC loop here

//--------------------------------------------Start--------------------------------------------------------
   StEmcDetector* detector = emcCollection->detector(kBarrelSmdEtaStripId);
  if (!detector)   {cout<<"There is no kBarrelEmcTowerId Detector"<<endl; return kStWarn;}
 
  for(UInt_t m=1;m<=120;m++)   {  StEmcModule* moduleE = detector->module(m);   StSPtrVecEmcRawHit& hitsE=moduleE->hits(); 
    for(int n=0;n<150;n++) mEnergyE[n]=0.0;
  for(unsigned int r=0;r<hitsE.size();r++)
  {unsigned int EtaE=hitsE[r]->eta();float eneE=hitsE[r]->energy();int myindex=EtaE-1;mEnergyE[myindex]=eneE;}

  for(unsigned int r=0;r<hitsE.size();r++) if(hitsE[r])
   {  float eE(0),phE(0),x1(0),y1(0),z1(0),eneE(0); Int_t didE(0), sE(0); UInt_t moduleE(0), EtaE(0),adcE1(0);
sE=abs(hitsE[r]->sub()); moduleE=hitsE[r]->module(); EtaE=hitsE[r]->eta();  eneE=hitsE[r]->energy(); adcE1=hitsE[r]->adc(); 
    mGeom2->getId(moduleE,EtaE,sE,didE); mGeom2->getEtaPhi(didE,eE,phE); mGeom2->getXYZ(didE,x1,y1,z1); 
 float en4(0),en5(0),en3(0),en2(0),en6(0),en1(0),en7(0),en0(0),en8(0),en9(0),en10(0),en11(0),en12(0),en01(0),en02(0),en03(0),en04(0),en05(0),en06(0),en07(0),en08(0),en13(0),en14(0),en15(0),en16(0); UInt_t eta4(0);
eta4=hitsE[r]->eta(); en4=mEnergyE[eta4-1];
if(!(en4>=0.5)) continue;
                             
if(eta4<150)     en5=mEnergyE[eta4];  else en5=0;if(eta4-2>=0)    en3=mEnergyE[eta4-2];else en3=0;
if(eta4+1<150)  en6=mEnergyE[eta4+1];else en6=0;if(eta4-3>=0)    en2=mEnergyE[eta4-3];else en2=0;
if(eta4+2<150)  en7=mEnergyE[eta4+2];else en7=0;if(eta4-4>=0)    en1=mEnergyE[eta4-4];else en1=0;
if(eta4+3<150)  en8=mEnergyE[eta4+3];else en8=0;if(eta4-5>=0)    en0=mEnergyE[eta4-5];else en0=0;

if(eta4+4<150)  en9=mEnergyE[eta4+4]; else en9=0;  if(eta4-6>=0)    en01=mEnergyE[eta4-6];else en01=0;
if(eta4+5<150)  en10=mEnergyE[eta4+5];else en10=0;if(eta4-7>=0)    en02=mEnergyE[eta4-7];else en02=0;
if(eta4+6<150)  en11=mEnergyE[eta4+6];else en11=0;if(eta4-8>=0)    en03=mEnergyE[eta4-8];else en03=0;
if(eta4+7<150)  en12=mEnergyE[eta4+7];else en12=0;if(eta4-9>=0)    en04=mEnergyE[eta4-9];else en04=0;

if(eta4+8<150)  en13=mEnergyE[eta4+8]; else en13=0;  if(eta4-10>=0)    en05=mEnergyE[eta4-10];else en05=0;
if(eta4+9<150)  en14=mEnergyE[eta4+9];else en14=0;if(eta4-11>=0)    en06=mEnergyE[eta4-11];else en06=0;
if(eta4+10<150)  en15=mEnergyE[eta4+10];else en15=0;if(eta4-12>=0)    en07=mEnergyE[eta4-12];else en07=0;
if(eta4+11<150)  en16=mEnergyE[eta4+11];else en16=0;if(eta4-13>=0)    en08=mEnergyE[eta4-13];else en08=0;
	                                                       
if(! (en4>=0.5 && en4>en5 && en4>en3 && en4>en6 && en4>en2 && en4>en7 && en4>en1 && en4>en0 && en4>en8 
&&en4>en01&&en4>en9&&en4>en02&&en4>en10&&
en4>en03&&en4>en11/*&&en4>en04&&en4>en12&&en4>en05&&en4>en13&&en4>en06&&en4>en14&&en4>en07&&en4>en15&&en4>en08&&en4>en16*/))continue; 
mGeom2->getId(moduleE,eta4,sE,didE); 
                       
float Eta3(0),Phi3(0),Eta4(0),Phi4(0),Eta5(0),Phi5(0),Eta2(0),Phi2(0),Eta6(0),Phi6(0);//Eta1,Phi1,Eta7,Phi7,Eta0,Phi0,Eta8,Phi8;
Int_t did3(0), did5(0),  did2(0),did6(0), did1(0),did7(0);
 did3=didE-1;  did5=didE+1; did2=didE-2;  did6=didE+2; did1=didE-3; did7=didE+3;//Int_t did0=didE-4; Int_t did8=didE+4;
	                                      
if(!(did3>0 && did5<=18000&&did2>0 && did6<=18000&&did1>0 && did7<=18000)) continue;

mGeom2->getEtaPhi(didE,Eta4,Phi4); mGeom2->getEtaPhi(did5,Eta5,Phi5); mGeom2->getEtaPhi(did3,Eta3,Phi3);
	//cout<<"module no"<<moduleE<<"did4 is"<<didE<<"Eta5 is"<<Eta5<<"::Eta4 is::"<<Eta4<<"::Eta3 is"<<Eta3<<endl;  
						
mGeom2->getEtaPhi(did2,Eta2,Phi2); mGeom2->getEtaPhi(did6,Eta6,Phi6);
 /* mGeom2->getEtaPhi(did1,Eta1,Phi1);  mGeom2->getEtaPhi(did7,Eta7,Phi7);   mGeom2->getEtaPhi(did0,Eta0,Phi0);  mGeom2->getEtaPhi(did8,Eta8,Phi8);*/
	                    
/*
		cout<<"pass SMDE one"<<endl<<"module is"<<moduleE<<endl;      
cout<<"en4 is"<<en4<<"en3"<<en3<<"en5"<<en5<<"en2"<<en2<<"en6"<<en6<<"en1"<<en1<<"en7"<<en7<<"en0"<<en0<<"en8"<<en8<<endl;
cout<<"eta4 is"<<eta4<<"Eta4 is"<<Eta4<<"Eta5 is"<<Eta5<<"Eta3 is"<<Eta3<<endl;cout<<"Phi4 is"<<Phi4<<"Phi5 is"<<Phi5<<"Phi3 is"<<Phi3<<endl;
cout<<"didE is"<<didE<<"did5 is"<<did5<<"did3 is"<<did3<<endl;cout<<"en4 is "<<en4<<endl;  cout<<"eta4 is"<<eta4<<endl;
*/	   
for(Int_t j=3; j<4; j++) 
{  StDetectorId id = static_cast<StDetectorId>(j+kBarrelEmcTowerId);StEmcDetector* detector=emcCollection->detector(id); if(detector)
								 													     
for(UInt_t y=1;y<=120;y++)
{ StEmcModule* module = detector->module(y);if(!(y==m))continue;StSPtrVecEmcRawHit& hits2=module->hits(); 
for(int Q=0;Q<10;Q++) for(int G=0;G<15;G++)	    mEnergyP[Q][G]=0.0;
for(UInt_t x=0;x<hits2.size();x++) if(hits2[x])
{UInt_t eP=hits2[x]->eta();Int_t sP=abs(hits2[x]->sub()); float enerP=hits2[x]->energy();
int Yindex=(sP-1);int myIndex=(eP-1);mEnergyP[myIndex][Yindex]=enerP;}

for(UInt_t x=0;x<hits2.size();x++) if(hits2[x])
{ Int_t didP(0); float PHIP(0),ETAP(0),x2(0),y2(0),z2(0), enerP(0);UInt_t modP(0), eP(0),adcP1(0); Int_t sP(0);
 modP=hits2[x]->module(); eP=hits2[x]->eta(); sP=abs(hits2[x]->sub()); enerP=hits2[x]->energy(); adcP1=hits2[x]->adc();
 mGeom3->getId(modP,eP,sP,didP); mGeom3->getEtaPhi(didP,ETAP,PHIP);  mGeom3->getXYZ(didP,x2,y2,z2); 

float enp4(0),enp5(0),enp3(0),enp2(0),enp6(0),enp1(0),enp7(0),enp0(0),enp8(0),enp9(0),enp10(0),enp11(0),enp12(0),enp01(0),enp02(0),enp03(0),enp04(0),
enp05(0),enp06(0),enp07(0),enp08(0),enp13(0),enp14(0),enp15(0),enp16(0); 
UInt_t sP4(0), eP4(0), sP5(0), sP3(0);
sP4=abs(hits2[x]->sub());eP4=hits2[x]->eta(); sP5=sP4+1; sP3=sP4-1;
 enp4=mEnergyP[eP4-1][sP4-1];
                        
if(!( eE<=(ETAP+0.05) && eE>=(ETAP-0.05))) continue; 
if(!(enp4>=0.5))continue;
									   
if(sP4<15 && eP4>=1 && eP4<=10)   enp5=mEnergyP[eP4-1][sP4];   else enp5=0; if((sP4-2)>=0 && eP4>=1 && eP4<=10)  enp3=mEnergyP[eP4-1][sP4-2]; else enp3=0;
if((sP4+1)<15 && eP4>=1 && eP4<=10) enp6=mEnergyP[eP4-1][sP4+1]; else enp6=0;if((sP4-3)>=0 && eP4>=1 && eP4<=10)  enp2=mEnergyP[eP4-1][sP4-3]; else enp2=0;
if(sP4+2<15 && eP4>=1 && eP4<=10) enp7=mEnergyP[eP4-1][sP4+2]; else enp7=0;if(sP4-4>=0 && eP4>=1 && eP4<=10)  enp1=mEnergyP[eP4-1][sP4-4]; else enp1=0;
if(sP4+3<15 && eP4>=1 && eP4<=10) enp8=mEnergyP[eP4-1][sP4+3]; else enp8=0;if(sP4-5>=0 && eP4>=1 && eP4<=10)  enp0=mEnergyP[eP4-1][sP4-5]; else enp0=0;
if(sP4-6>=0 && eP4>=1 && eP4<=10)  enp01=mEnergyP[eP4-1][sP4-6]; else enp01=0;if(sP4-7>=0 && eP4>=1 && eP4<=10)  enp02=mEnergyP[eP4-1][sP4-7]; else enp02=0;
if(sP4-8>=0 && eP4>=1 && eP4<=10)  enp03=mEnergyP[eP4-1][sP4-8]; else enp03=0;if(sP4-9>=0 && eP4>=1 && eP4<=10)  enp04=mEnergyP[eP4-1][sP4-9]; else enp04=0;
if(sP4+4<15 && eP4>=1 && eP4<=10) enp9=mEnergyP[eP4-1][sP4+4]; else enp9=0;if(sP4+5<15 && eP4>=1 && eP4<=10) enp10=mEnergyP[eP4-1][sP4+5]; else enp10=0;
if(sP4+6<15 && eP4>=1 && eP4<=10) enp11=mEnergyP[eP4-1][sP4+6]; else enp11=0;if(sP4+7<15 && eP4>=1 && eP4<=10) enp12=mEnergyP[eP4-1][sP4+7]; else enp12=0; 
if(sP4-10>=0 && eP4>=1 && eP4<=10)  enp05=mEnergyP[eP4-1][sP4-10]; else enp05=0;if(sP4-11>=0 && eP4>=1 && eP4<=10)  enp06=mEnergyP[eP4-1][sP4-11]; else enp06=0;
if(sP4-12>=0 && eP4>=1 && eP4<=10)  enp07=mEnergyP[eP4-1][sP4-12]; else enp07=0;if(sP4-13>=0 && eP4>=1 && eP4<=10)  enp08=mEnergyP[eP4-1][sP4-13]; else enp08=0;
if(sP4+8<15 && eP4>=1 && eP4<=10) enp13=mEnergyP[eP4-1][sP4+8]; else enp13=0;if(sP4+9<15 && eP4>=1 && eP4<=10) enp14=mEnergyP[eP4-1][sP4+9]; else enp14=0;
if(sP4+10<15 && eP4>=1 && eP4<=10) enp15=mEnergyP[eP4-1][sP4+10]; else enp15=0;if(sP4+11<15 && eP4>=1 && eP4<=10) enp16=mEnergyP[eP4-1][sP4+11]; else enp16=0; 

if (!(enp4>=0.5 && enp4>enp5 && enp4>enp3 && enp4>enp6 && enp4>enp2 && enp4>enp7 && enp4>enp1 && enp4>enp0 && enp4>enp8&&
enp4>enp01&&enp4>enp9&&enp4>enp02&&enp4>enp10&&enp4>enp03&&enp4>enp11/*&&enp4>enp04&&enp4>enp12&&enp4>enp05&&enp4>enp13&&enp4>enp06&&enp4>enp14
&&enp4>enp07&&enp4>enp15&&enp4>enp08&&enp4>enp16*/)) continue; 

 mGeom3->getId(modP,eP,sP4,didP);   
						 
float EtaP3(0),PhiP3(0),EtaP4(0),PhiP4(0),EtaP5(0),PhiP5(0),EtaP2(0),PhiP2(0),EtaP6(0),PhiP6(0);//EtaP1,PhiP1,EtaP7,PhiP7,EtaP0,PhiP0,EtaP8,PhiP8;   
Int_t didP3(0),didP5(0),didP2x(0),didP6(0),didP1(0),didP7(0);				  
didP3=didP-10;  didP5=didP+10;  didP2x=didP-20;  didP6=didP+20;didP1=didP-30; didP7=didP+30; 
//Int_t didP0=didP-40; Int_t didP8=didP+40; 
						  
if(!(didP3>0 && didP5<=18000&&didP2x>0 && didP6<=18000&&didP1>0 && didP7<=18000)) continue;
mGeom3->getEtaPhi(didP,EtaP4,PhiP4); mGeom3->getEtaPhi(didP5,EtaP5,PhiP5);  mGeom3->getEtaPhi(didP3,EtaP3,PhiP3);
mGeom3->getEtaPhi(didP2x,EtaP2,PhiP2); mGeom3->getEtaPhi(didP6,EtaP6,PhiP6); 
/*
		      cout<<"pass smdp one"<<endl<<"module is"<<modP<<endl;
cout<<"enp4 is"<<enp4<<"enp3"<<enp3<<"enp5"<<enp5<<"enp2"<<enp2<<"enp6"<<enp6<<"enp1"<<enp1<<"enp7"<<enp7<<"enp0"<<enp0<<"enp8"<<enp8<<endl;
cout<<"EtaP4 is"<<EtaP4<<"EtaP5 is"<<EtaP5<<"EtaP3 is"<<EtaP3<<endl;cout<<"sP4 is"<<sP4<<"sP5 is"<<sP5<<"sP3 is"<<sP3<<endl;
cout<<"PhiP4 is"<<PhiP4<<"PhiP5 is"<<PhiP5<<"PhiP3 is"<<PhiP3<<endl;cout<<"didP is"<<didP<<"didP5 is"<<didP5<<"didP3 is"<<didP3<<endl;
*/
 for(Int_t j1=0; j1<1; j1++) 
 {  
  StDetectorId id = static_cast<StDetectorId>(j1+kBarrelEmcTowerId);
  StEmcDetector* detector=emcCollection->detector(id); //if(detector)		  
 
 float eT(0),phT(0),XT(0),YT(0),ZT(0); Int_t didT(0); float PTower(0),pidTower(0),TPpT(0),TPeT(0),eneT0(0),hpedSS(0),hRMSSS(0);unsigned int adc11(0);int TC7(0),hSS(0);
 float THeta(0);Int_t sT(0); UInt_t moduleT(0),EtaT(0);
 
 if(detector)	for(unsigned int i=1;i<121;i++)
           {StEmcModule* moduleTg=detector->module(i);if(!(i==m)) continue; 
           StSPtrVecEmcRawHit& hits=moduleTg->hits(); 
	   for(unsigned int k=0;k<hits.size();k++)
       	   {  eneT0=hits[k]->energy(); if(!(eneT0 >=6)) continue;
	   sT=abs(hits[k]->sub());moduleT=hits[k]->module();EtaT=hits[k]->eta(); adc11=hits[k]->adc();
	   
	 mGeom0->getId(moduleT,EtaT,sT,didT); mGeom0->getEtaPhi(didT,eT,phT); mGeom0->getXYZ(didT,XT,YT,ZT); mGeom0->getTheta(moduleT,EtaT,THeta); 
	if (!( moduleE==moduleT && eE-0.0035> (eT-0.025) && eE+0.0035<(eT+0.025))) continue;
        if( !(modP==moduleT && PHIP+0.0035<phT+0.025 &&PHIP-0.0035>phT-0.025)) continue; 
		      	
		hSS=mStatusT[didT-1]; hpedSS=mCalibPedValues[didT-1]; hRMSSS=mCalibPedRMSValues[didT-1];
               for (Int_t ku=0; ku<TowerHVchangeMax; ku++) if (didT == HVTowerid[ku]) TC7=1;  	   
               PTower=ParentTrackP[didT-1];pidTower=ParticleId[didT-1];
	 TPeT=ParentTrackPOSe[didT-1];  TPpT=ParentTrackPOSP[didT-1];


  
float Energyhq(0);  Int_t Subhq(0);  UInt_t Modulehq(0), Etahq(0),Adchq(0); int didTq(0); float eThq(0),phiThq(0),Thetahq(0);
float pedhhq(0),RMShhq(0);Int_t sphhq(0);

float Enpq1(0),Ennq1(0),Enpq20(0),Ennq20(0),Enpq19(0),Ennq19(0),Enpq21(0),Ennq21(0);
int towerpq1(0),towernq1(0),towerpq20(0),towernq20(0),towerpq19(0),towernq19(0),towerpq21(0),towernq21(0);
int Idpq1(0),Idnq1(0),Idpq20(0),Idnq20(0),Idpq19(0),Idnq19(0),Idpq21(0),Idnq21(0);  
unsigned int adcpq1(0),adcnq1(0),adcpq20(0),adcnq20(0), adcpq19(0),adcnq19(0),adcpq21(0),adcnq21(0);
float thetapq1(0),thetanq1(0),thetapq20(0),thetanq20(0),thetapq19(0),thetanq19(0),thetapq21(0),thetanq21(0);
float ppq1(0),pnq1(0),ppq20(0),pnq20(0),ppq19(0),pnq19(0),ppq21(0),pnq21(0);
float pidpq1(0),pidnq1(0),pidpq20(0),pidnq20(0),pidpq19(0),pidnq19(0),pidpq21(0),pidnq21(0);
UInt_t Modpq1(0),Modnq1(0),Modpq20(0),Modnq20(0),Modpq19(0),Modnq19(0),Modpq21(0),Modnq21(0);
Int_t sUbpq1(0),sUbnq1(0),sUbpq20(0),sUbnq20(0),sUbpq19(0),sUbnq19(0),sUbpq21(0),sUbnq21(0);
UInt_t eTapq1(0),eTanq1(0),eTapq20(0),eTanq20(0),eTapq19(0),eTanq19(0),eTapq21(0),eTanq21(0);

int nuOfNiq(0), nextIdq(0), nextIdTq(0),NextIdq(0), NextIdTq(0); 
float nextHiEtq(0),nextHiEq(0);
float nextPedq(0),nextRMSq(0),nextPidq(0),nextPq(0); unsigned int nextAdcq(0); Int_t nextStq(0);
float pedpq1(0),pednq1(0),pedpq19(0),pednq19(0),pedpq21(0),pednq21(0),pedpq20(0),pednq20(0);
float RMSpq1(0), RMSnq1(0),RMSpq19(0),RMSnq19(0),RMSpq21(0),RMSnq21(0),RMSpq20(0),RMSnq20(0); 
Int_t spq1(0),snq1(0),spq19(0),snq19(0),spq21(0),snq21(0),spq20(0),snq20(0);




 for(UInt_t mwq=1;mwq<=120;mwq++)   
  { 
  StEmcModule* module1 = detector->module(mwq);    StSPtrVecEmcRawHit& hitsEwq=module1->hits(); 
for(unsigned int rwq=0;rwq<hitsEwq.size();rwq++) 
{ Modulehq=hitsEwq[rwq]->module(); Subhq=abs(hitsEwq[rwq]->sub()); Etahq=hitsEwq[rwq]->eta();
mGeom0->getId(Modulehq,Etahq,Subhq,didTq);mGeom0->getEtaPhi(didTq,eThq,phiThq);mGeom0->getTheta(Modulehq,Etahq,Thetahq);

if(abs(didT-didTq)!=20&&abs(didT-didTq)!=19&&abs(didT-didTq)!=21&&abs(didT-didTq)!=1&&abs(didT-didTq)!=2380&&abs(didT-didTq)!=2381&&
abs(didT-didTq)!=2379&&eT>0&&eThq>0) continue;
if(abs(didT-didTq)!=20&&abs(didT-didTq)!=19&&abs(didT-didTq)!=21&&abs(didT-didTq)!=1&&abs(didT-didTq)!=2380&&abs(didT-didTq)!=2381&&
abs(didT-didTq)!=2379&&eT<0&&eThq<0) continue;

if(eT>0&&eThq<0&&fabs(fabs(eT)-fabs(eThq))!=0)continue; if(eT<0&&eThq>0&&fabs(fabs(eT)-fabs(eThq))!=0)continue;
if(fabs(fabs(eT)-fabs(eThq))>0.09 ||fabs(fabs(phT)-fabs(phiThq))>0.09) continue;
if(eT>0&&eThq<0&&fabs(fabs(eT)-fabs(eThq))==0&&eT>0.05&&eThq<-0.05)continue;
if(eT<0&&eThq>0&&fabs(fabs(eT)-fabs(eThq))==0&&eThq>0.05&&eT<-0.05)continue;
if(eT>0&&eThq<0&&fabs(fabs(eT)-fabs(eThq))==0&&eT<0.05&&eThq>-0.05&&phT>0&&phiThq<0)continue;//add here to these two conditions tower Id at the edge when phi
if(eT>0&&eThq<0&&fabs(fabs(eT)-fabs(eThq))==0&&eT<0.05&&eThq>-0.05&&phT<0&&phiThq>0)continue;
if(eT<0&&eThq>0&&fabs(fabs(eT)-fabs(eThq))==0&&eT<0.05&&eThq>-0.05&&phT>0&&phiThq<0)continue;
if(eT<0&&eThq>0&&fabs(fabs(eT)-fabs(eThq))==0&&eT<0.05&&eThq>-0.05&&phT<0&&phiThq>0)continue;

Energyhq=hitsEwq[rwq]->energy(); Adchq=hitsEwq[rwq]->adc();sphhq=mStatusT[didTq-1];pedhhq=mCalibPedValues[didTq-1];RMShhq=mCalibPedRMSValues[didTq-1];	 
nuOfNiq++;	 
//cout<<"didTq is::"<<didTq<<"::Energyq is::"<<Energyhq<<"::dphi::"<<fabs(fabs(phT)-fabs(phiThq))<<
//"::deta::"<<fabs(fabs(eT)-fabs(eThq))<<"::phiThq is::"<<phiThq<<"::eThq::"<<eThq<<endl;	 
	 

if(didTq==didT+1&&phT==phiThq)
{Enpq1=Energyhq;Idpq1=didTq;adcpq1=Adchq;thetapq1=Thetahq;towerpq1=didT-Idpq1;
Modpq1=Modulehq;sUbpq1=Subhq;eTapq1=Etahq;pedpq1=pedhhq;RMSpq1=RMShhq;spq1=sphhq;}

if((didTq==didT-1&&phT==phiThq)||(fabs(fabs(phT)-fabs(phiThq))<0.001&&eThq==-1*eT))
{Ennq1=Energyhq;Idnq1=didTq;adcnq1=Adchq;thetanq1=Thetahq;towernq1=didT-Idnq1;
Modnq1=Modulehq;sUbnq1=Subhq;eTanq1=Etahq;pednq1=pedhhq;RMSnq1=RMShhq;snq1=sphhq;}//here add across the sides

if((didTq==didT+20||didTq==didT-2380)&&eT==eThq)
{Enpq20=Energyhq;Idpq20=didTq;adcpq20=Adchq;thetapq20=Thetahq;towerpq20=didT-Idpq20;
Modpq20=Modulehq;sUbpq20=Subhq;eTapq20=Etahq;pedpq20=pedhhq;RMSpq20=RMShhq;spq20=sphhq;}

if((didTq==didT-20||didTq==didT+2380)&&eT==eThq)
{Ennq20=Energyhq;Idnq20=didTq;adcnq20=Adchq;thetanq20=Thetahq;towernq20=didT-Idnq20;
Modnq20=Modulehq;sUbnq20=Subhq;eTanq20=Etahq;pednq20=pedhhq;RMSnq20=RMShhq;snq20=sphhq;}

if((didTq==didT+21||didTq==didT-2381)&&fabs(fabs(eT)-fabs(eThq))<0.09&&fabs(fabs(phT)-fabs(phiThq))<0.09)
{Enpq21=Energyhq;Idpq21=didTq;adcpq21=Adchq;thetapq21=Thetahq;towerpq21=didT-Idpq21;
Modpq21=Modulehq;sUbpq21=Subhq;eTapq21=Etahq;pedpq21=pedhhq;RMSpq21=RMShhq;spq21=sphhq;}

if(((didTq==didT-21||didTq==didT+2381)&&fabs(fabs(eT)-fabs(eThq))<0.09&&fabs(fabs(phT)-fabs(phiThq))<0.09)
||(fabs(fabs(phT)-fabs(phiThq))>0.001&&(fabs(fabs(phT)-fabs(phiThq))<0.09&&fabs(phT)>fabs(phiThq)&&eThq==-1*eT)))
{Ennq21=Energyhq;Idnq21=didTq;adcnq21=Adchq;thetanq21=Thetahq;towernq21=didT-Idnq21;
Modnq21=Modulehq;sUbnq21=Subhq;eTanq21=Etahq;pednq21=pedhhq;RMSnq21=RMShhq;snq21=sphhq;}

if(((didTq==didT+19||didTq==didT-2379)&&fabs(fabs(eT)-fabs(eThq))<0.09&&fabs(fabs(phT)-fabs(phiThq))<0.09)
||(fabs(fabs(phT)-fabs(phiThq))>0.001&&(fabs(fabs(phT)-fabs(phiThq))<0.09&&fabs(phT)<fabs(phiThq)&&eThq==-1*eT)))
{Enpq19=Energyhq;Idpq19=didTq;adcpq19=Adchq;thetapq19=Thetahq;towerpq19=didT-Idpq19;Modpq19=Modulehq;
sUbpq19=Subhq;eTapq19=Etahq;pedpq19=pedhhq;RMSpq19=RMShhq;spq19=sphhq;}

if((didTq==didT-19||didTq==didT+2379)&&fabs(fabs(eT)-fabs(eThq))<0.09&&fabs(fabs(phT)-fabs(phiThq))<0.09)
{Ennq19=Energyhq;Idnq19=didTq;adcnq19=Adchq;thetanq19=Thetahq;towernq19=didT-Idnq19;Modnq19=Modulehq;
sUbnq19=Subhq;eTanq19=Etahq;pednq19=pedhhq;RMSnq19=RMShhq;snq19=sphhq;}




ppq1=ParentTrackP[Idpq1-1];pidpq1=ParticleId[Idpq1-1];pnq1=ParentTrackP[Idnq1-1];pidnq1=ParticleId[Idnq1-1];
ppq20=ParentTrackP[Idpq20-1];pidpq20=ParticleId[Idpq20-1];pnq20=ParentTrackP[Idnq20-1];pidnq20=ParticleId[Idnq20-1];
ppq19=ParentTrackP[Idpq19-1];pidpq19=ParticleId[Idpq19-1];pnq19=ParentTrackP[Idnq19-1];pidnq19=ParticleId[Idnq19-1];
ppq21=ParentTrackP[Idpq21-1];pidpq21=ParticleId[Idpq21-1];pnq21=ParentTrackP[Idnq21-1]; pidnq21=ParticleId[Idnq21-1];

 
if(Enpq1*sin(thetapq1)>nextHiEtq) 
{nextHiEtq=Enpq1*sin(thetapq1);nextIdTq=Idpq1;NextIdTq=1;nextPedq=pedpq1;nextRMSq=RMSpq1;nextPidq=pidpq1;
nextPq=ppq1;nextAdcq=adcpq1;nextStq=spq1;}

if(Ennq1*sin(thetanq1)>nextHiEtq) 
{nextHiEtq=Ennq1*sin(thetanq1);nextIdTq=Idnq1;NextIdTq=-1;nextPedq=pednq1;nextRMSq=RMSnq1;nextPidq=pidnq1;
nextPq=pnq1;nextAdcq=adcnq1;nextStq=snq1;}

if(Enpq21*sin(thetapq21)>nextHiEtq)
{nextHiEtq=Enpq21*sin(thetapq21);nextIdTq=Idpq21;NextIdTq=21;nextPedq=pedpq21;nextRMSq=RMSpq21;nextPidq=pidpq21;
nextPq=ppq21;nextAdcq=adcpq21;nextStq=spq21;}

if(Ennq21*sin(thetanq21)>nextHiEtq)
{nextHiEtq=Ennq21*sin(thetanq21);nextIdTq=Idnq21;NextIdTq=-21;nextPedq=pednq21;nextRMSq=RMSnq21;nextPidq=pidnq21;
nextPq=pnq21;nextAdcq=adcnq21;nextStq=snq21;}

if(Enpq19*sin(thetapq19)>nextHiEtq)
{nextHiEtq=Enpq19*sin(thetapq19);nextIdTq=Idpq19;NextIdTq=19;nextPedq=pedpq19;nextRMSq=RMSpq19;nextPidq=pidpq19;
nextPq=ppq19;nextAdcq=adcpq19;nextStq=spq19;}

if(Ennq19*sin(thetanq19)>nextHiEtq)
{nextHiEtq=Ennq19*sin(thetanq19);nextIdTq=Idnq19;NextIdTq=-19;nextPedq=pednq19;nextRMSq=RMSnq19;nextPidq=pidnq19;
nextPq=pnq19;nextAdcq=adcnq19;nextStq=snq19;}

if(Enpq20*sin(thetapq20)>nextHiEtq)
{nextHiEtq=Enpq20*sin(thetapq20);nextIdTq=Idpq20;NextIdTq=20;nextPedq=pedpq20;nextRMSq=RMSpq20;nextPidq=pidpq20;
nextPq=ppq20;nextAdcq=adcpq20;nextStq=spq20;}

if(Ennq20*sin(thetanq20)>nextHiEtq)
{nextHiEtq=Ennq20*sin(thetanq20);nextIdTq=Idnq20;NextIdTq=-20;nextPedq=pednq20;nextRMSq=RMSnq20;nextPidq=pidnq20;
nextPq=pnq20;nextAdcq=adcnq20;nextStq=snq20;}

if(Enpq1>nextHiEq) {nextHiEq=Enpq1; nextIdq=Idpq1;NextIdq=1;}if(Ennq1>nextHiEq) {nextHiEq=Ennq1; nextIdq=Idnq1;NextIdq=-1;}
if(Enpq21>nextHiEq){ nextHiEq=Enpq21;nextIdq=Idpq21;NextIdq=21;}if(Ennq21>nextHiEq){ nextHiEq=Ennq21;nextIdq=Idnq21;NextIdq=-21;}
if(Enpq19>nextHiEq){ nextHiEq=Enpq19; nextIdq=Idpq19;NextIdq=19;}if(Ennq19>nextHiEq){ nextHiEq=Ennq19; nextIdq=Idnq19;NextIdq=-19;}
if(Enpq20>nextHiEq) {nextHiEq=Enpq20;nextIdq=Idpq20;NextIdq=20;}if(Ennq20>nextHiEq){ nextHiEq=Ennq20;nextIdq=Idnq20;NextIdq=-20;}
}
	 
 }
 

float Ehighest(0),EneT0(0),Enearest(0),EtHighest(0),EtNearest(0);

Ehighest=eneT0+nextHiEq; EtHighest=eneT0*sin(THeta)+nextHiEtq;
	   
// cout<<"eneT0 is"<<eneT0<<endl;
 
//-----------------------------------Find the nearest tower--------------------------------------- 
 StThreeVectorF VECTOR3(x2,y2,z1); const   StThreeVectorF& v3=VECTOR3-vertexPosition;

float u01(0),u02(0),u03(0),u11(0),u12(0),u13(0),u21(0),u22(0),u23(0),u31(0),u32(0),u33(0);
float u41(0),u42(0),u43(0),u51(0),u52(0),u53(0),u61(0),u62(0),u63(0),u71(0),u72(0),u73(0);

if(Idpq1!=0) mGeom0->getXYZ(Idpq1,u01,u02,u03);if(Idnq1!=0) mGeom0->getXYZ(Idnq1,u11,u12,u13);
if(Idpq20!=0) mGeom0->getXYZ(Idpq20,u21,u22,u23);if(Idnq20!=0) mGeom0->getXYZ(Idnq20,u31,u32,u33);
if(Idpq19!=0) mGeom0->getXYZ(Idpq19,u41,u42,u43);if(Idnq19!=0) mGeom0->getXYZ(Idnq19,u51,u52,u53);
if(Idpq21!=0) mGeom0->getXYZ(Idpq21,u61,u62,u63);if(Idnq21!=0) mGeom0->getXYZ(Idnq21,u71,u72,u73);

 StThreeVectorF u0(u01,u02,u03); StThreeVectorF u1(u11,u12,u13); StThreeVectorF u2(u21,u22,u23); StThreeVectorF u3(u31,u32,u33);
 StThreeVectorF u4(u41,u42,u43); StThreeVectorF u5(u51,u52,u53); StThreeVectorF u6(u61,u62,u63); StThreeVectorF u7(u71,u72,u73);

const   StThreeVectorF& U0=u0-vertexPosition;	const   StThreeVectorF& U1=u1-vertexPosition;
const   StThreeVectorF& U2=u2-vertexPosition;   const   StThreeVectorF& U3=u3-vertexPosition;
const   StThreeVectorF& U4=u4-vertexPosition;	const   StThreeVectorF& U5=u5-vertexPosition;
const   StThreeVectorF& U6=u6-vertexPosition;   const   StThreeVectorF& U7=u7-vertexPosition;
//=========================================================================================//
float angleu0v(0),angleu1v(0),angleu2v(0),angleu3v(0),angleu4v(0),angleu5v(0),angleu6v(0),angleu7v(0);   
       
angleu0v=v3.angle(U0); angleu1v=v3.angle(U1);angleu2v=v3.angle(U2);angleu3v=v3.angle(U3);
angleu4v=v3.angle(U4); angleu5v=v3.angle(U5);angleu6v=v3.angle(U6);angleu7v=v3.angle(U7);

//cout<<": angleu0v:"<<angleu0v<<":angleu1v:"<<angleu1v<<":angleu2v:"<<angleu2v<<":angleu3v:"<<angleu3v<<endl;
//cout<<":angleu4v:"<<angleu4v<<":angleu5v:"<<angleu5v<<":angleu6v:"<<angleu6v<<":angleu7v:"<<angleu7v<<endl;

float delta1(0),dphi1(0);		delta1=fabs(eE)-fabs(eT); dphi1=fabs(PHIP)-fabs(phT);
float nearestHiEq(0),nearestP(0),nearestPidq(0);int nearestIdq(0),NearestIdq(0); UInt_t nearestAdc(0);



if(angleu0v<angleu1v&&angleu0v<angleu2v&&angleu0v<angleu3v&&angleu0v<angleu4v&&angleu0v<angleu5v&&angleu0v<angleu6v&&angleu0v<angleu7v) 
{EneT0=eneT0+Enpq1; Enearest=eneT0+Enpq1;EtNearest=eneT0*sin(THeta)+Enpq1*sin(thetapq1);nearestHiEq=Enpq1; nearestIdq=Idpq1;NearestIdq=1;
  nearestAdc=adcpq1;nearestP=ppq1;nearestPidq=pidpq1; //cout<<":the neraest energy is:"<<Enpq1<<endl;
}

if(angleu1v<angleu0v&&angleu1v<angleu2v&&angleu1v<angleu3v&&angleu1v<angleu4v&&angleu1v<angleu5v&&angleu1v<angleu6v&&angleu1v<angleu7v) 
{EneT0=eneT0+Ennq1;  Enearest=eneT0+Ennq1;EtNearest=eneT0*sin(THeta)+Ennq1*sin(thetanq1);nearestHiEq=Ennq1; nearestIdq=Idnq1;NearestIdq=-1;
nearestAdc=adcnq1;nearestP=pnq1;nearestPidq=pidnq1;
//cout<<":the neraest energy is:"<<Ennq1<<endl;
 }

if(angleu2v<angleu1v&&angleu2v<angleu0v&&angleu2v<angleu3v&&angleu2v<angleu4v&&angleu2v<angleu5v&&angleu2v<angleu6v&&angleu2v<angleu7v) 
{EneT0=eneT0+Enpq20; Enearest=eneT0+Enpq20;EtNearest=eneT0*sin(THeta)+Enpq20*sin(thetapq20);nearestHiEq=Enpq20;
nearestIdq=Idpq20;NearestIdq=20;nearestAdc=adcpq20;nearestP=ppq20;nearestPidq=pidpq20;
//cout<<":the neraest energy is:"<<Enpq20<<endl;
}

if(angleu3v<angleu1v&&angleu3v<angleu2v&&angleu3v<angleu0v&&angleu3v<angleu4v&&angleu3v<angleu5v&&angleu3v<angleu6v&&angleu3v<angleu7v ) 
{EneT0=eneT0+Ennq20; Enearest=eneT0+Ennq20;EtNearest=eneT0*sin(THeta)+Ennq20*sin(thetanq20);nearestHiEq=Ennq20;
nearestIdq=Idnq20;NearestIdq=-20;nearestAdc=adcnq20;nearestP=pnq20;nearestPidq=pidnq20; 
//cout<<":the neraest energy is:"<<Ennq20<<endl;
}


if(angleu4v<angleu1v&&angleu4v<angleu2v&&angleu4v<angleu3v&&angleu4v<angleu0v&&angleu4v<angleu5v&&angleu4v<angleu6v&&angleu4v<angleu7v)
{EneT0=eneT0+Enpq19;Enearest=eneT0+Enpq19;EtNearest=eneT0*sin(THeta)+Enpq19*sin(thetapq19);nearestHiEq=Enpq19; nearestIdq=Idpq19;NearestIdq=19;
nearestAdc=adcpq19;nearestP=ppq19;nearestPidq=pidpq19;
//cout<<":the neraest energy is:"<<Enpq19<<endl;
}

if(angleu5v<angleu1v&&angleu5v<angleu2v&&angleu5v<angleu3v&&angleu5v<angleu4v&&angleu5v<angleu0v&&angleu5v<angleu6v&&angleu5v<angleu7v)
{EneT0=eneT0+Ennq19;Enearest=eneT0+Ennq19;EtNearest=eneT0*sin(THeta)+Ennq19*sin(thetanq19);nearestHiEq=Ennq19;
nearestIdq=Idnq19;NearestIdq=-19;nearestAdc=adcnq19;nearestP=pnq19;nearestPidq=pidnq19; 
// cout<<":the neraest energy is:"<<Ennq19<<endl;
}

if(angleu6v<angleu1v&&angleu6v<angleu2v&&angleu6v<angleu3v&&angleu6v<angleu4v&&angleu6v<angleu5v&&angleu6v<angleu0v&&angleu6v<angleu7v)
{EneT0=eneT0+Enpq21;Enearest=eneT0+Enpq21;EtNearest=eneT0*sin(THeta)+Enpq21*sin(thetapq21);nearestHiEq=Enpq21; nearestIdq=Idpq21;NearestIdq=21;
nearestAdc=adcpq21;nearestP=ppq21;nearestPidq=pidpq21; 
// cout<<":the neraest energy is:"<<Enpq21<<endl;
}


if(angleu7v<angleu1v&&angleu7v<angleu2v&&angleu7v<angleu3v&&angleu7v<angleu4v&&angleu7v<angleu5v&&angleu7v<angleu6v&&angleu7v<angleu0v)
{EneT0=eneT0+Ennq21;Enearest=eneT0+Ennq21;EtNearest=eneT0*sin(THeta)+Ennq21*sin(thetanq21);nearestHiEq=Ennq21; nearestIdq=Idnq21;NearestIdq=-21;
nearestAdc=adcnq21;nearestP=pnq21;nearestPidq=pidnq21;
// cout<<":the neraest energy is:"<<Ennq21<<endl;
}

if(fabs(delta1)<=0.018&&fabs(dphi1)<=0.018) EneT0=eneT0;   if(EneT0==0) EneT0=eneT0;
//cout<<"smde Inf::"<<"moduleE"<<m<<"::EtaE::"<<EtaE<<"::sE is::"<<sE<<"::eE is"<<eE<<".."<<Eta5<<""<<Eta3<<""<<Eta2<<""<<Eta6<<endl;
//cout<<"en4 is"<<en4<<"::en3 is::"<<en3<<"::en5::"<<en5<<"::en2 "<<en2<<"::en6"<<en6<<"::en1::"<<en1<<":en7:"<<en7<<":en0"<<en0<<":en8"<<en8<<endl;
//cout<<"smdp Inf::"<<"modP"<<y<<"::ETAP::"<<ETAP<<":sP is"<<sP<<"::PHIP is::"<<PHIP<<".."<<PhiP5<<""<<PhiP3<<""<<PhiP2<<""<<PhiP6<<endl;
//cout<<"enp4 is"<<enp4<<"::enp3 is::"<<enp3<<"::enp5::"<<enp5<<"::enp2 is"<<enp2<<"::enp6 is"<<enp6<<"::enp1::"<<enp1<<":enp7:"<<enp7<<":enp0"<<enp0<<":enp8"<<enp8<<endl;
//cout<<"Towers Inf::"<<"moduleT"<<i<<"::EtaT::"<<EtaT<<"::sT is::"<<sT<<"::phT is:"<<phT<<"didT is"<<didT<<"eT is"<<eT<<endl;			
		float etav1(0),phiv1(0);		etav1 = v3.pseudoRapidity();  		phiv1=v3.phi(); 
	if(en5<0) en5=0;   if(enp5<0) enp5=0;    if(en3<0) en3=0;    if(enp3<0) enp3=0;	if(en6<0) en6=0;   if(enp6<0) enp6=0;    if(en2<0) en2=0;    if(enp2<0) enp2=0;
	if(en7<0) en7=0;   if(enp7<0) enp7=0;    if(en1<0) en1=0;    if(enp1<0) enp1=0;	if(en8<0) en8=0;   if(enp8<0) enp8=0;    if(en0<0) en0=0;    if(enp0<0) enp0=0;
	if(en9<0) en9=0;   if(enp9<0) enp9=0;    if(en01<0) en01=0;  if(enp01<0) enp01=0;if(en10<0) en10=0; if(enp10<0) enp10=0;  if(en02<0) en02=0;if(enp02<0) enp02=0;
	if(en11<0) en11=0; if(enp11<0) enp11=0;  if(en03<0) en03=0;if(enp03<0) enp03=0;	if(en12<0) en12=0; if(enp12<0) enp12=0;  if(en04<0) en04=0;if(enp04<0) enp04=0;
	if(en13<0) en13=0; if(enp13<0) enp13=0;  if(en05<0) en05=0;if(enp05<0) enp05=0;	if(en14<0) en14=0; if(enp14<0) enp14=0;  if(en06<0) en06=0;if(enp06<0) enp06=0;
	if(en15<0) en15=0; if(enp15<0) enp15=0;  if(en07<0) en07=0;if(enp07<0) enp07=0;	if(en16<0) en16=0; if(enp16<0) enp16=0;  if(en08<0) en08=0;if(enp08<0) enp08=0;
		
 float eff(0),eN4(0),eN5(0),eN6(0),eN7(0),eN8(0),eN9(0),eN10(0),eN11(0),eN12(0),eN13(0),eN14(0),eN15(0),eN16(0),eN3(0),eN2(0),eN1(0),eN0(0),
 eN01(0),eN02(0),eN03(0),eN04(0),eN05(0),eN06(0),eN07(0),eN08(0);
 float eNp4(0),eNp5(0),eNp6(0),eNp7(0),eNp8(0),eNp9(0),eNp10(0),eNp11(0),eNp12(0),eNp13(0),eNp14(0),eNp15(0),eNp16(0),eNp3(0),eNp2(0),eNp1(0),eNp0(0),
 eNp01(0),eNp02(0),eNp03(0),eNp04(0),eNp05(0),eNp06(0),eNp07(0),eNp08(0);
float total(0),WEN(0),ENN0(0),ENN1(0),ENN2(0),ENN3(0),ENN4(0),ENN5(0),ENN6(0),ENN7(0),ENN8(0),ENN9(0),ENN10(0),ENN11(0),ENN12(0);


	eN4=en4*en4; eN5=en5*en5; eN6=en6*en6; eN7=en7*en7; eN8=en8*en8; eN9=en9*en9; eN10=en10*en10; eN11=en11*en11; eN12=en12*en12; eN13=en13*en13; eN14=en14*en14;
	eN15=en15*en15;	eN16=en16*en16; 
	eN3=en3*en3; eN2=en2*en2; eN1=en1*en1; eN0=en0*en0; eN01=en01*en01; eN02=en02*en02; eN03=en03*en03; eN04=en04*en04; eN05=en05*en05; eN06=en06*en06; 
	eN07=en07*en07; eN08=en08*en08;
	
	eNp4=enp4*enp4; eNp5=enp5*enp5; eNp6=enp6*enp6; eNp7=enp7*enp7; eNp8=enp8*enp8; eNp9=enp9*enp9; eNp10=enp10*enp10; eNp11=enp11*enp11; eNp12=enp12*enp12; eNp13=enp13*enp13; eNp14=enp14*enp14;
	eNp15=enp15*enp15; eNp16=enp16*enp16; 
	eNp3=enp3*enp3; eNp2=enp2*enp2; eNp1=enp1*enp1; eNp0=enp0*enp0; eNp01=enp01*enp01; eNp02=enp02*enp02; eNp03=enp03*enp03; eNp04=enp04*enp04; eNp05=enp05*enp05; eNp06=enp06*enp06; 
	eNp07=enp07*enp07; eNp08=enp08*enp08;
	
	total=eN4+eN5+eN6+eN7+eN8+eN9+eN10+eN11+eN12+eN13+eN14+eN15+eN16+eN3+eN2+eN1+eN0+eN01+eN02+eN03+eN04+eN05+eN06+eN07+eN08+ 
	eNp4+eNp5+eNp6+eNp7+eNp8+eNp9+eNp10+eNp11+eNp12+eNp13+eNp14+eNp15+eNp16+eNp3+eNp2+eNp1+eNp0+eNp01+eNp02+eNp03+eNp04+eNp05+eNp06+eNp07+eNp08;
	
	ENN0=(en4+enp4)*0.783; 	ENN1=(en3+enp3+en5+enp5)*2.21;	ENN2=(en2+enp2+en6+enp6)*6.26; 	ENN3=(en1+enp1+en7+enp7)*11.51;
	ENN4=(en0+enp0+en8+enp8)*17.73;	ENN5=(en01+enp01+en9+enp9)*24.78;ENN6=(en02+enp02+en10+enp10)*32.57; ENN7=(en03+enp03+en11+enp11)*41.05;
	ENN8=(en04+enp04+en12+enp12)*50.15; ENN9=(en05+enp05+en13+enp13)*59.84;	ENN10=(en06+enp06+en14+enp14)*70.09;
	ENN11=(en07+enp07+en15+enp15)*80.86;ENN12=(en08+enp08+en16+enp16)*92.13;
	
	WEN=ENN0+ENN1+ENN2+ENN3+ENN4+ENN5+ENN6+ENN7+ENN8+ENN9+ENN10+ENN11+ENN12;
	
	
	eff=(sqrt(total))/WEN;
       float Reff(0), efftow(0), Refftow(0); Reff=1/eff; efftow=WEN/eneT0; Refftow=1/efftow;	
		
			float Weta1(0),Wphi1(0), Wboth(0);
		 Weta1=((en0+en1+en2+en3+en5+en6+en7+en8)/(en0+en1+en2+en3+en4+en5+en6+en7+en8));
		 Wphi1=((enp0+enp1+enp2+enp3+enp5+enp6+enp7+enp8)/(enp0+enp1+enp2+enp3+enp4+enp5+enp6+enp7+enp8));
		 Wboth=((en0+en1+en2+en3+en5+en6+en7+en8+enp0+enp1+enp2+enp3+enp5+enp6+enp7+enp8)/(en0+en1+en2+en3+en4+en5+en6+en7+en8+enp0+enp1+enp2+enp3+enp4+enp5+enp6+enp7+enp8));
	
		 //	cout<<"sqrt(total)::"<<sqrt(total)<<"::WEN::"<<WEN<<"::eff::"<<eff<<":Reff:"<<Reff<<"::efftow:"<<efftow<<":Refftow"<<Refftow<<endl;	
		 //	cout<<"Weta1::"<<Weta1<<"::Wphi1::"<<Wphi1<<"::Wboth::"<<Wboth<<endl;
		//----------------------------------------------------Trigger again-------------------------
		UInt_t H(0);
		if(mMinb &&(!mCentral) &&(!mHt1) &&(!mHt2)) H=1;		if(mHt1&&(!mCentral) &&(!mMinb) &&(!mHt2)) H=2;		
		if(mHt2 &&(!mCentral) &&(!mMinb) &&(!mHt1)) H=3;		if(mCentral &&(!mHt2) &&(!mMinb) &&(!mHt1)) H=4;
		if(mMinb && mCentral) H=5;		if(mMinb && mHt1) H=6;		if(mMinb && mHt2) H=7;
		if(mHt1 && mHt2) H=8;		if(mHt1 && mCentral) H=9;		if(mHt2 && mCentral) H=10;
		if(mMinb && mCentral&&mHt1) H=11;		if(mMinb && mCentral&&mHt2) H=12;		if(mCentral&&mHt1&&mHt2) H=13;
		if(mMinb && mHt1&&mHt2) H=14;		if(mMinb && mCentral&&mHt1&&mHt2) H=15;		if(mJPsi&&mJPsi1) H=16;
		if(mJPsi&&(!mJPsi1)) H=17;		if(mJPsi1&&(!mJPsi)) H=18;	
		if(mL2gamma1 && (!mL2gamma2) && (!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) H=20;
		if(mL2gamma2 && (!mL2gamma1) && (!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) H=21;
		if(mbht2mb1 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) H=22;
		if(mbht2mb2 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb3) &&(!mbht2mb4) &&(!mbht2mb5)) H=23;
		if(mbht2mb3 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb4) &&(!mbht2mb5)) H=24;
		if(mbht2mb4 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb5)) H=25;
		if(mbht2mb5 && (!mL2gamma2) && (!mL2gamma1) &&(!mbht2mb1) &&(!mbht2mb2) &&(!mbht2mb3) &&(!mbht2mb4)) H=26;
		if(mL2gamma1 && mL2gamma2) H=27;if(mL2gamma1 && mbht2mb1) H=28;	if(mL2gamma1 && mbht2mb2) H=29;
		if(mL2gamma1 && mbht2mb3) H=30;	if(mL2gamma1 && mbht2mb4) H=31;	if(mL2gamma1 && mbht2mb5) H=32;
     		if(mL2gamma2 && mbht2mb1) H=33;	if(mL2gamma2 && mbht2mb2) H=34;	if(mL2gamma2 && mbht2mb3) H=35;
		if(mL2gamma2 && mbht2mb4) H=36;	if(mL2gamma2 && mbht2mb5) H=37;	if(mbht2mb1 && mbht2mb2) H=38;	
		if(mbht2mb1 && mbht2mb3) H=39;	if(mbht2mb1 && mbht2mb4) H=40;	if(mbht2mb1 && mbht2mb5) H=41;
		if(mbht2mb2 && mbht2mb3) H=42;	if(mbht2mb2 && mbht2mb4) H=43;	if(mbht2mb2 && mbht2mb5) H=44;
		if(mbht2mb3 && mbht2mb4) H=45;	if(mbht2mb3 && mbht2mb5) H=46;	if(mbht2mb4 && mbht2mb5) H=47;
		if(mHt0) H=48 ; if(mht4) H=49; if(mht4fast) H=50; if(mupsilondAu) H=51;
		if(mHt0 && (!mHt1) && (!mHt2) &&(!mht4) &&(!mht4fast) &&(!mupsilondAu)) H=52;
		if((!mHt0) && mHt1 && (!mHt2) &&(!mht4) &&(!mht4fast) &&(!mupsilondAu)) H=53;
		if((!mHt0) && (!mHt1) && mHt2 &&(!mht4) &&(!mht4fast) &&(!mupsilondAu)) H=54;
		if((!mHt0) && (!mHt1) && (!mHt2) && mht4 &&(!mht4fast) &&(!mupsilondAu)) H=55;
		if((!mHt0) && (!mHt1) && (!mHt2) && (!mht4) && mht4fast &&(!mupsilondAu)) H=56;
		if((!mHt0) && (!mHt1) && (!mHt2) && (!mht4) && (!mht4fast) && mupsilondAu) H=57;
		if(mHt0 && mHt1) H=58; if(mHt0 && mHt2) H=59; if(mHt0 && mht4) H=60;
		if(mHt1 && mHt2) H=61; if(mHt1 && mht4) H=62;  
		if(mHt2 && mht4) H=63; if(mbht2mb6) H=64;
		 //-----------------------------L2gamma2010----------------
		if(mL2gamma1_2010) H=100;
		if(mL2gamma2_2010) H=101;
		//-----------------------------L2gamma2011----------------
		if(mL2gamma0_2011) H=200;
		if(mL2gamma1_2011) H=201;
		if(mL2gamma2_2011) H=202;
		//----------------------------------------------------------
		//		                  		cout<<"H is"<<H<<endl;

		//________L2gamma pp2009
                if(mL2gamma_pp2009) H=102;


float ENET0(0),ENET02(0);
ENET0=(Enpq1+Ennq1+Enpq20+Ennq20);
ENET02=(Enpq1+Ennq1+Enpq20+Ennq20+Enpq19+Ennq19+Enpq21+Ennq21);
    
 float vf2(0),vf2F(0),v2B(0);   vf2=cos(2*(phT)-angleRe); vf2F=cos(2*(phT)-angleRef);
        
	        v2B=cos(2*(phT)-psi_f);
		
		float TSPB(0), ehB(0);
	        
	   TSPB=(eneT0/((en4+enp4)*0.783+(en3+enp3+en5+enp5)*2.21+(en2+enp2+en6+enp6)*6.26+(en1+enp1+en7+enp7)*11.51+(en0+enp0+en8+enp8)*17.73+(en01+enp01+en9+enp9)*24.78+(en02+enp02+en10+enp10)*32.57+(en03+enp03+en11+enp11)*41.05));
	        
	   ehB=((Ehighest-2*eneT0)/Ehighest); 


	   cout<<" TSP: "<<TSPB<<"  EneT0: "<<EneT0<<" eT : "<<eT<<endl;
		
		mTupleFiveL->Fill(hiTower,Subh,Moduleh,Etah,Adch,T,evt,z,runId,Idh,eTh,phiTh,Thetah,PPTPC);
                mTupleCAB->Fill(didE,ENET0,XX,YY,adc11,didP,adcE1,adcP1,eT,phT,hSS,hpedSS,hRMSSS,TC7);
               	mTupleQAAB->Fill(EneT0,MPos,MNeg,H,en4,enp4,eneT0,z,etav1,phiv1,didT,delta1,dphi1,Weta1);
		mTupleAngleAB->Fill(moduleT,EtaE,sP,runId,evt,NOP,Wphi1,enp7,enp1,enp8,enp0,en01,en02,PTower);
		mTupleSMDE1B->Fill(en3,en5,en2,en6,en1,en7,en0,en8,enp3,enp5,enp2,enp6,en03,en04);
		mTuplePiB->Fill(en9,en10,en11,en12,enp01,enp02,enp03,enp04,enp9,enp10,enp11,enp12,pidTower,Enpq1);
		mTuplePhoB->Fill(sT,THeta,Ennq1,Enpq20,Ennq20,Enpq19,Ennq19,Enpq21,NumOfTracks,muRef,nV,rank1,rank,noOfPiVer);
		mTuplePi1B->Fill(en05,enp05,en06,enp06,en07,enp07,en08,enp08,en13,enp13,en14,enp14,en15,enp15);
mTupleNew8->Fill(psi_f0,psi_f1,psi_fn0,psi_fn1,angleu0v,angleu1v,angleu2v,angleu3v,angleu4v,angleu5v,angleu6v,angleu7v);

mTupleFiveteenB->Fill(angleRe,angleReS1,angleReS2,angleRef,angleRefS1,angleRefS2,NumOfTracks,Q1,Q2,XR,YR,XRf,YRf);
mTupleSixteenB->Fill(XRS1,XRS2,YRS1,YRS2,QQ,Q3,Q4,Q1S1,Q1S2,Q2S1,Q2S2,Q3S1,Q3S2,Q4S1);
mTupleSeventeenB->Fill(Q4S2,XRfS1,XRfS2,YRfS1,YRfS2,XRfW,YRfW,XRfE,YRfE,XRfSW1,YRfSW1,XRfSW2,YRfSW2,XRfSE1);
mTupleEighteenB->Fill(YRfSE1,XRfSE2,YRfSE2,psi_e,psi_w,psi_f,mBBCAdcCut,vf2,vf2F,psi_fn,v2B,Ehighest,mBBCAdcCut0);

 mTupleNew7->Fill(angleRept,angleReS1pt,angleReS2pt,vpdz,BBCz,ZDCz,TSPB,ehB);
  mTupleNew6->Fill(angleRepe,angleReS1pe,angleReS2pe,angleRene,angleReS1ne,angleReS2ne,angleRenc,angleReS1nc,angleReS2nc,angleRepc,angleReS1pc,angleReS2pc);

mTupleFTPCQ->Fill(angleRefW,angleRefE,angleRefW1,angleRefE1,angleRefW2,angleRefE2);

//REFILL HERE


//-----------------------------------------------TRACKS---------------------------------------------------





 for(size_t nodeIndex=0;nodeIndex<trackNodes.size();nodeIndex++)
   {size_t numberOfTracksInNode=trackNodes[nodeIndex]->entries(primary);

     for( size_t trackIndex=0;trackIndex<numberOfTracksInNode;trackIndex++)
       { //cout<<"the track index is"<<trackIndex<<endl;
	 track=trackNodes[nodeIndex]->track(primary,trackIndex);
	 // if (!(track && track->flag()>0))continue;
	 float pTPC(0), ptTPC(0),dEdx(0),EtaTr(0),epsi(0); short ch(0);
	 float energ(0),energSt(0);


	 StTrackGeometry* geometry=track->geometry();ch=geometry->charge(); 
	 pTPC=geometry->momentum().mag();ptTPC=geometry->momentum().perp();

	 //	 count_trk_sttrak++;
	 //if(!(ptTPC>=0.1&&ptTPC<50.)) continue;
	 EtaTr=geometry->momentum().pseudoRapidity();epsi=geometry->momentum().phi();
	 Int_t Fp=track->fitTraits().numberOfFitPoints(kTpcId);         
	 //if(Fp<20) continue;


	 Int_t Ppo=track->numberOfPossiblePoints(); float xg=geometry->origin().x();
	 float yg=geometry->origin().y();float zg=geometry->origin().z(); float px=geometry->momentum().x();
	 float py=geometry->momentum().y();float pz=geometry->momentum().z(); float lg=track->length(); 
	 int flag= track->flag();

	 StPtrVecTrackPidTraits traits=track->pidTraits(kTpcId); float dca=track->impactParameter(); 
	 //if(dca>3) continue;
	 StDedxPidTraits* pid;


	 cout<<trackIndex<<" th particle info in StTrack: "<<" Fp= "<<Fp<<"  Eta: "<<EtaTr<<" pT: "<<ptTPC<<endl;

	 for(UInt_t V=0;V<traits.size();V++)
	   {pid=dynamic_cast<StDedxPidTraits*>(traits[V]);if(pid&&pid->method()==kTruncatedMeanId)
							    
							    { dEdx=pid->mean();  float dphi(0), dphif(0),dphiin(0),dphiou(0);
							    

							      dphi=(atan2(sin(epsi-phT),cos(epsi-phT))); dphif=dphi;
							      if(dphi<(-3.14/2)) dphif=(2*3.14)+dphi;
							      
							      if(fabs(fabs(phT)-fabs(angleRe))<3.14/4 || fabs(fabs(phT)-fabs(angleRe))>3*3.14/4) dphiin=dphi; 
					      if (dphi<(-3.14/2)) dphiin=(2*3.14)+dphi;//it was&&
					      if(fabs(fabs(phT)-fabs(angleRe))>3.14/4 && fabs(fabs(phT)-fabs(angleRe))<3*3.14/4)dphiou=dphi;
					      if (dphi<(-3.14/2))dphiou=(2*3.14)+dphi;
					      //							      cout<<"::dphif is::"<<dphif<<"::dphiin::"<<dphiin<<"::dphiou:"<<dphiou<<endl;

					      energ=EneT0*sin(2*atan(exp(-1*etav1)));energSt=EneT0*sin(THeta);

					      float TSP, eh;

					      TSP=(eneT0/((en4+enp4)*0.783+(en3+enp3+en5+enp5)*2.21+(en2+enp2+en6+enp6)*6.26+(en1+enp1+en7+enp7)*11.51+(en0+enp0+en8+enp8)*17.73+(en01+enp01+en9+enp9)*24.78+(en02+enp02+en10+enp10)*32.57+(en03+enp03+en11+enp11)*41.05));
					      eh=((Ehighest-2*eneT0)/Ehighest);								  


					      mTupleTrack->Fill(pTPC,ptTPC,dEdx,EtaTr,epsi,ch,dca,Fp,eE,phE,ETAP,PHIP,eT,phT);
					      mTupleCA->Fill(didE,ENET0,XX,YY,adc11,didP,adcE1,adcP1,eT,phT,hSS,hpedSS,hRMSSS,TC7);
					      mTupleQAA->Fill(EneT0,MPos,MNeg,H,en4,enp4,eneT0,z,etav1,phiv1,didT,delta1,dphi1,Weta1);
					      mTupleAngleA->Fill(moduleT,EtaE,sP,runId,evt,NOP,Wphi1,enp7,enp1,enp8,enp0,en01,en02,PTower);
					      mTupleSMDE1->Fill(en3,en5,en2,en6,en1,en7,en0,en8,enp3,enp5,enp2,enp6,en03,en04);
					      mTuplePi->Fill(en9,en10,en11,en12,enp01,enp02,enp03,enp04,enp9,enp10,enp11,enp12,pidTower,Enpq1);
					      mTuplePho->Fill(sT,eff,energ,dphif,energSt,THeta,Ennq1,Enpq20,Ennq20,Enpq19,Ennq19,Enpq21,dphiin,dphiou);
					      mTuplePi1->Fill(en05,enp05,en06,enp06,en07,enp07,en08,enp08,en13,enp13,en14,enp14,en15,enp15);


					      mTupleNew8A->Fill(psi_f0,psi_f1,psi_fn0,psi_fn1,Ppo,xg,yg,zg,px,py,pz);
					      mTupleFiveLL->Fill(hiTower,Subh,Moduleh,Etah,Adch,T,evt,z,runId,Idh,eTh,phiTh,Thetah,PPTPC);
					      mTupleNew9A->Fill(angleu0v,angleu1v,angleu2v,angleu3v,angleu4v,angleu5v,angleu6v,angleu7v,NumOfTracks,muRef,nV,rank1,rank,noOfPiVer);
			  

					      mTupleFiveteenA->Fill(angleRe,angleReS1,angleReS2,angleRef,angleRefS1,angleRefS2,NumOfTracks,Q1,Q2,XR,YR,XRf,YRf);
					      mTupleSixteenA->Fill(XRS1,XRS2,YRS1,YRS2,QQ,Q3,Q4,Q1S1,Q1S2,Q2S1,Q2S2,Q3S1,Q3S2,Q4S1);
					      mTupleSeventeenA->Fill(Q4S2,XRfS1,XRfS2,YRfS1,YRfS2,XRfW,YRfW,XRfE,YRfE,XRfSW1,YRfSW1,XRfSW2,YRfSW2,XRfSE1);
					      mTupleEighteenA->Fill(YRfSE1,XRfSE2,YRfSE2,psi_e,psi_w,psi_f,mBBCAdcCut,vf2,vf2F,psi_fn,v2B,Ehighest,mBBCAdcCut0);

					      mTupleNew7A->Fill(angleRept,angleReS1pt,angleReS2pt,vpdz,BBCz,ZDCz,TSP,eh);
					      mTupleNew6A->Fill(angleRepe,angleReS1pe,angleReS2pe,angleRene,angleReS1ne,angleReS2ne,angleRenc,angleReS1nc,angleReS2nc,angleRepc,angleReS1pc,angleReS2pc);
					      mTupleFTPCQA->Fill(angleRefW,angleRefE,angleRefW1,angleRefE1,angleRefW2,angleRefE2);

							    }
	   }
       }
   }
	 ///////________From StMuTrack_______
Int_t nPrimary1= 0;
//StSPtrVecTrackNode& trackNodes=mStEvent->trackNodes();
 // StTrack* track;
nPrimary1 = mMuDstMaker->muDst()->primaryTracks()->GetEntries();
//cout<<"Hi in primary tracks  "<<nPrimary1<<endl;
 
    for ( Int_t jtrk = 0; jtrk < nPrimary1; jtrk++ )
    {
      float pTPC(0), ptTPC(0),dEdx(0),EtaTr(0),epsi(0); short ch(0);
      float energ(0),energSt(0);
      
	ch=mMuDstMaker->muDst()->primaryTracks(jtrk)->charge();
	
	pTPC=mMuDstMaker->muDst()->primaryTracks(jtrk)->p().mag();     
	ptTPC=mMuDstMaker->muDst()->primaryTracks(jtrk)->pt();; 
    
	EtaTr=mMuDstMaker->muDst()->primaryTracks(jtrk)->eta();
	epsi=mMuDstMaker->muDst()->primaryTracks(jtrk)->phi();  
	float	xg(0); float yg(0); float zg(0);
	Int_t Fp=mMuDstMaker->muDst()->primaryTracks(jtrk)->nHitsFit(); // Here Fit points from TPC  
	//	if(Fp<20) continue;

	Int_t FpT(0);
	Int_t Ppo=mMuDstMaker->muDst()->primaryTracks(jtrk)->nHitsPoss();

	float       px = pTPC*cos( epsi );
	float       py = pTPC*sin( epsi );
	float       pz = sqrt( pTPC*pTPC - ptTPC*ptTPC );  
	int flag(0);
	float dca=mMuDstMaker->muDst()->primaryTracks(jtrk)->dcaGlobal().mag(); //Returns 3D distance of closest approach to primary vertex of associated global track. 
	float pdca=mMuDstMaker->muDst()->primaryTracks(jtrk)->dca().mag();  //Returns 3D distance of closest approach to primary vertex. 
	float dcaD=mMuDstMaker->muDst()->primaryTracks(jtrk)->dcaD(); //Signed radial component of global DCA (projected). 
	float dcaZ=mMuDstMaker->muDst()->primaryTracks(jtrk)->dcaZ(); //Z component of global DCA. 

       	dEdx = mMuDstMaker->muDst()->primaryTracks(jtrk)->dEdx();
    	float dphi=0.;
	float dphif=0.; 
	float dphiin=0.; 
	float dphiou=0.;

	Int_t vrtxIndx(-99999);
	vrtxIndx = mMuDstMaker->muDst()->primaryTracks(jtrk)->vertexIndex();


	 cout<<jtrk<<" th particle info in StMuTrack: "<<" Fp= "<<Fp<<"  Eta: "<<EtaTr<<" pT: "<<ptTPC<<endl;

	//_________________SOME TRACK QUALITY CONDTIONS APPLIED_________
	//       	if(!((Fp/Ppo) > 0.52) ) continue;
	//	if(!(ptTPC<50.&&ptTPC>1.2)) continue; // this is used to minimize the ntuple size
	


	//_______________________________________________________________

    	dphi=(atan2(sin(epsi-phT),cos(epsi-phT))); dphif=dphi; if(dphi<(-3.14/2)) dphif=(2*3.14)+dphi; 
    
	if(fabs(fabs(phT)-fabs(angleRe))<3.14/4 || fabs(fabs(phT)-fabs(angleRe))>3*3.14/4) dphiin=dphi; if (dphi<(-3.14/2)) dphiin=(2*3.14)+dphi;//it was&&
	if(fabs(fabs(phT)-fabs(angleRe))>3.14/4 && fabs(fabs(phT)-fabs(angleRe))<3*3.14/4)dphiou=dphi; if (dphi<(-3.14/2)) dphiou=(2*3.14)+dphi;
	
	//	cout<<"::dphif is::"<<dphif<<"::dphiin::"<<dphiin<<"::dphiou:"<<dphiou<<endl;
	//	cout<<"Hi i am in Primary trck "<<ch<<"  "<<pTPC<<endl;
    
	energ=EneT0*sin(2*atan(exp(-1*etav1)));
	energSt=EneT0*sin(THeta);

	float TSP, eh;
	        
	TSP=(eneT0/((en4+enp4)*0.783+(en3+enp3+en5+enp5)*2.21+(en2+enp2+en6+enp6)*6.26+(en1+enp1+en7+enp7)*11.51+(en0+enp0+en8+enp8)*17.73+(en01+enp01+en9+enp9)*24.78+(en02+enp02+en10+enp10)*32.57+(en03+enp03+en11+enp11)*41.05));
	        
	eh=((Ehighest-2*eneT0)/Ehighest); 
	     


//-------------------------------------------------------------------------------------------------
	mTupleTrack_M->Fill(pTPC,ptTPC,dEdx,EtaTr,epsi,ch,dca,Fp,eE,phE,ETAP,PHIP,eT,phT);
	mTupleCA_M->Fill(didE,ENET0,XX,YY,adc11,didP,adcE1,adcP1,eT,phT,hSS,hpedSS,hRMSSS,TC7);
	mTupleQAA_M->Fill(EneT0,MPos,MNeg,H,en4,enp4,eneT0,z,etav1,phiv1,didT,delta1,dphi1,Weta1);
	mTupleAngleA_M->Fill(moduleT,EtaE,sP,runId,evt,NOP,Wphi1,enp7,enp1,enp8,enp0,en01,en02,PTower);
	mTupleSMDE1_M->Fill(en3,en5,en2,en6,en1,en7,en0,en8,enp3,enp5,enp2,enp6,en03,en04);
	mTuplePi_M->Fill(en9,en10,en11,en12,enp01,enp02,enp03,enp04,enp9,enp10,enp11,enp12,pidTower,Enpq1);
	mTuplePho_M->Fill(sT,eff,energ,dphif,energSt,THeta,Ennq1,Enpq20,Ennq20,Enpq19,Ennq19,Enpq21,dphiin,dphiou);
	mTuplePi1_M->Fill(en05,enp05,en06,enp06,en07,enp07,en08,enp08,en13,enp13,en14,enp14,en15,enp15);
	mTupleNew8A_M->Fill(psi_f0,psi_f1,psi_fn0,psi_fn1,Ppo,xg,yg,zg,px,py,pz,pdca,dcaD,dcaZ);
	mTupleFiveLL_M->Fill(hiTower,Subh,Moduleh,Etah,Adch,T,evt,z,runId,Idh,eTh,phiTh,Thetah,PPTPC);
	mTupleNew9A_M->Fill(angleu0v,angleu1v,angleu2v,angleu3v,angleu4v,angleu5v,angleu6v,angleu7v,NumOfTracks,muRef,nV,rank1,rank,noOfPiVer);
	mTupleFiveteenA_M->Fill(angleRe,angleReS1,angleReS2,angleRef,angleRefS1,angleRefS2,NumOfTracks,Q1,Q2,XR,YR,XRf,YRf);
	mTupleSixteenA_M->Fill(XRS1,XRS2,YRS1,YRS2,QQ,Q3,Q4,Q1S1,Q1S2,Q2S1,Q2S2,Q3S1,Q3S2,Q4S1);
	mTupleSeventeenA_M->Fill(Q4S2,XRfS1,XRfS2,YRfS1,YRfS2,XRfW,YRfW,XRfE,YRfE,XRfSW1,YRfSW1,XRfSW2,YRfSW2,XRfSE1);
	mTupleEighteenA_M->Fill(YRfSE1,XRfSE2,YRfSE2,psi_e,psi_w,psi_f,mBBCAdcCut,vf2,vf2F,psi_fn,v2B,Ehighest,mBBCAdcCut0);

	mTupleNew7A_M->Fill(angleRept,angleReS1pt,angleReS2pt,vpdz,BBCz,ZDCz,TSP,eh,vrtxIndx);
	mTupleNew6A_M->Fill(angleRepe,angleReS1pe,angleReS2pe,angleRene,angleReS1ne,angleReS2ne,angleRenc,angleReS1nc,angleReS2nc,angleRepc,angleReS1pc,angleReS2pc);
	mTupleFTPCQA_M->Fill(angleRefW,angleRefE,angleRefW1,angleRefE1,angleRefW2,angleRefE2);




    }
// track loop end

	   }
	
	   }
 }
 }
		
 }
 
  
 }
	
   }
  }
	    
	 
	    
	  
return kStOk;
//loop: //break;
}
//-----------------------------------------------------------------------------------------------------
Int_t StThirdMaker::Finish()
{
 //write histogram
myGraphFile->Write();
myGraphFile->Close();
 return kStOk;
 }
//-------------------------------------------------------------------------------------------------------
void StThirdMaker::getTrgType() 
{  
  mHt1 = kFALSE;  mHt2 = kFALSE;  mMinb = kFALSE;  mCentral = kFALSE;  mJPsi=kFALSE;   mJPsi1=kFALSE;     mL2gamma1=kFALSE;    mL2gamma2=kFALSE;
   mbht2mb1=kFALSE;    mbht2mb2=kFALSE;    mbht2mb3=kFALSE;    mbht2mb4=kFALSE;    mbht2mb5=kFALSE; mbht2mb6=kFALSE;
    mHt0=kFALSE; mht4=kFALSE; mht4fast=kFALSE; mupsilondAu=kFALSE;mL2gamma1_2010=kFALSE;mL2gamma2_2010=kFALSE;
     mL2gamma0_2011=kFALSE;mL2gamma1_2011=kFALSE;mL2gamma2_2011=kFALSE;
     if (mStEvent->triggerIdCollection()) 
    {      if (mStEvent->triggerIdCollection()->nominal())
      {        if (                // 2003 dAu minimum bias
               mStEvent->triggerIdCollection()->nominal()->isTrigger(2001)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(2003)
            // 2003 pp minimum bias
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(1001) || mStEvent->triggerIdCollection()->nominal()->isTrigger(1000)
            // 2004 AuAu minimum bias
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(15003) || mStEvent->triggerIdCollection()->nominal()->isTrigger(15007)
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(25007)
            // 2004 63 GeV AuAu minimum bias
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(35001)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(35004)
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(35007)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(35009)
	    //2004 pp minumum bias   
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(45010)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(45020)
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(10)
          //cucu2005
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(66007)
	 //pp2006
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(117001)
	  ) mMinb = true;
      
      //dAu 2008   
       if ( mStEvent->triggerIdCollection()->nominal()->isTrigger(210500) //ht0
       || mStEvent->triggerIdCollection()->nominal()->isTrigger(210501)   //ht0 
       ) mHt0 = true;
       
        if (           
             // 2003 dAu HighTower1 trigger
               mStEvent->triggerIdCollection()->nominal()->isTrigger(2201)
            // 2003 pp HighTower1 trigger
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(1201)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(1101)
            // 2004 AuAu HighTower1 trigger
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(15202)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(15203)
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(25203)
	    //2004 pp HighTower1 trigger   
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(45201)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(45512)
	    //cucu2005
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(66201)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(66203) 
	    //pp2005
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(96201)
	    //dAu 2008   
	      || mStEvent->triggerIdCollection()->nominal()->isTrigger(210510) //ht1
	      || mStEvent->triggerIdCollection()->nominal()->isTrigger(210511)  //ht1
           ) mHt1 = true;
         
        if (   
             // 2003 dAu HighTower2 trigger
               mStEvent->triggerIdCollection()->nominal()->isTrigger(2202)
            // 2003 pp HighTower2 trigger
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(1202) || mStEvent->triggerIdCollection()->nominal()->isTrigger(1102)
	    //2004 pp HighTower2 trigger   
	       || mStEvent->triggerIdCollection()->nominal()->isTrigger(45202)   
	    //pp2005
	    || mStEvent->triggerIdCollection()->nominal()->isTrigger(96211) 
	    //pp2006
	    || mStEvent->triggerIdCollection()->nominal()->isTrigger(117211)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(117212)
	    || mStEvent->triggerIdCollection()->nominal()->isTrigger(127212)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(127213)
	    || mStEvent->triggerIdCollection()->nominal()->isTrigger(137213) 
	     //dAu 2008   
	     || mStEvent->triggerIdCollection()->nominal()->isTrigger(210520)//ht2 
	      || mStEvent->triggerIdCollection()->nominal()->isTrigger(210521)//ht2    
	       	       		
           ) mHt2 = true; 
	   //dAu2008
	   if ( mStEvent->triggerIdCollection()->nominal()->isTrigger(210541)//ht4 
	      	   )mht4= true; 
	 if ( mStEvent->triggerIdCollection()->nominal()->isTrigger(210800)//ht4fast 
	      	   )mht4fast= true; 
	if (	mStEvent->triggerIdCollection()->nominal()->isTrigger(210601)  //upsilon 	   
            )mupsilondAu=true;
	    
        if (               // 2004 AuAu central trigger
               mStEvent->triggerIdCollection()->nominal()->isTrigger(15105)
            // 2004 63 GeV AuAu minimum bias
               || mStEvent->triggerIdCollection()->nominal()->isTrigger(35004)|| mStEvent->triggerIdCollection()->nominal()->isTrigger(35009)
           ) mCentral = true;
	   
	 if (                // JPsi trigger
              mStEvent->triggerIdCollection()->nominal()->isTrigger(117602)
            //7090046-7129067
           ) mJPsi = true;   
	   
	  if (                // JPsi trigger
              mStEvent->triggerIdCollection()->nominal()->isTrigger(137603)
            //7133052-7156040
           ) mJPsi1 = true;   
	
	//-------------------------------AuAu2007---------------------
	
	if (                // L2gamma1 trigger first run:8109015	last run:8113068
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200620)
                     ) mL2gamma1 = true; 
	   
	 if (                // L2gamma2 trigger first run:8113102	last run:8177038
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200621)
                     ) mL2gamma2 = true; 
	      
	if (                // bht2mb trigger first run:8094008	last run:8102029
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200211)
                     ) mbht2mb1 = true;    
	   
	if (               // bht2mb trigger first run:8103029	last run:8106004
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200212)
                     ) mbht2mb2 = true;    
	      
	 if (               // bht2mb triggerfirst run: 8106050	last run:8108014
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200220)
                     ) mbht2mb3 = true;    
	  if (              // bht2mb trigger first run:8109015	last run:8113068
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200221)
                     ) mbht2mb4 = true;    
	  if (                // bht2mb trigger first run:8113102	last run:8177038
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200222)
                     ) mbht2mb5 = true;    
  if (                // bht2mb trigger first run:8113102	last run:8177038
              mStEvent->triggerIdCollection()->nominal()->isTrigger(200586)
                     ) mbht2mb6 = true;    

//-------------------------------AuAu2010---------------------
	
	if (                // L2gamma1 trigger first run:8109015	last run:8113068
              mStEvent->triggerIdCollection()->nominal()->isTrigger(260800)
                     ) mL2gamma1_2010 = true; 
	   
	 if (                // L2gamma2 trigger first run:8113102	last run:8177038
              mStEvent->triggerIdCollection()->nominal()->isTrigger(260820)
                     ) mL2gamma2_2010 = true; 
		     
//-------------------------------AuAu2011---------------------
	
	if (                // L2gamma1 trigger first run:12122026	last run:12152016
              mStEvent->triggerIdCollection()->nominal()->isTrigger(20)
                     ) mL2gamma0_2011 = true; 
	   
	 if (                // L2gamma2 trigger first run:12122025	last run:12122025
              mStEvent->triggerIdCollection()->nominal()->isTrigger(260820)
                     ) mL2gamma1_2011 = true; 	
		     
	 if (                // L2gamma2 trigger first run:12153002	last run:12179051
              mStEvent->triggerIdCollection()->nominal()->isTrigger(350069)
                     ) mL2gamma2_2011 = true; 			     	     
	 //__________pp 200 Run09________________

         if(  //L2BGamma trigger fro pp2009
            mStEvent->triggerIdCollection()->nominal()->isTrigger(240620)
              ) mL2gamma_pp2009 = true;


      }   
    } 


}

//----------------------------------------------------------------------------------------
 Int_t StThirdMaker:: associateHadronsWithBEMC3()
{
   
   StEmcCollection* emcCollection = mStEvent->emcCollection(); //Get emcCollection
   if(!emcCollection) //return kStOK;
    {  cout<<"There is no emcCollection"<<endl;       return kStWarn;}
  
   
  Double_t bFld;  StEventSummary* summary=mStEvent->summary();  if(summary){ bFld=summary->magneticField()/10.;}  if (fabs(bFld)<0.01){return kStWarn;}
 
 Int_t mod,eta,sub; StEmcPosition*  pos=new StEmcPosition(); StThreeVectorD position,momentum;
 
 StSPtrVecTrackNode& trackNodes=mStEvent->trackNodes(); StTrack* track;

 StEmcGeom * mGeom0=StEmcGeom::instance("bemc"); StEmcGeom * mGeom2=StEmcGeom::instance("bsmde");      StEmcGeom * mGeom3=StEmcGeom::instance("bsmdp");
  
 for(size_t nodeIndex=0;nodeIndex<trackNodes.size();nodeIndex++)
  {   size_t numberOfTracksInNode=trackNodes[nodeIndex]->entries(global);
 
   for( size_t trackIndex=0;trackIndex<numberOfTracksInNode;trackIndex++)
  { 
  track=trackNodes[nodeIndex]->track(global,trackIndex);   StTrackGeometry* geometry=track->geometry(); 
       short ch=geometry->charge();       float pTPC=geometry->momentum().mag();       
        float ptTPC=geometry->momentum().perp();             Int_t Fp=track->fitTraits().numberOfFitPoints(kTpcId);//if(Fp<20) continue;
     float dca=track->impactParameter(); //if(dca>3) continue;

     //     cout<<" Hello i am projecting ...check dca  "<<dca<<endl;
     
      StPtrVecTrackPidTraits traits=track->pidTraits(kTpcId);       StDedxPidTraits* pid;
       for(UInt_t V=0;V<traits.size();V++)
       { pid=dynamic_cast<StDedxPidTraits*>(traits[V]);
          if(pid&&pid->method()==kTruncatedMeanId)
          { float dEdx=pid->mean();    //if(!(dEdx>0.0000033&&dEdx<0.0000046))  continue; if (!(track && track->flag()>=0)) continue;

                   Bool_t ok=pos->trackOnEmc(&position,&momentum,track,bFld);    if(! ok) continue;
   
   
   
   StEmcDetector* detector = emcCollection->detector(kBarrelSmdEtaStripId);  if (!detector)
   {cout<<"There is no kBarrelSmdEtaStripId Detector"<<endl; return kStWarn;}
    mGeom2->getBin(position.phi(),position.pseudoRapidity(),mod,eta,sub);
    //   cout<<"the track coordinate at bsmde is"<<mod<<endl<<eta<<endl<<sub<<endl<<position.pseudoRapidity()<<endl<<position.phi()<<endl;                                
   // mDist->Fill(mod,eta);
			
for(unsigned int i=1;i<=120;i++)
  {  
   if(fabs(mod)!=i) continue;
      StEmcModule* module=detector->module(i);       StSPtrVecEmcRawHit& hits=module->hits();
	  //  cout<<"the size of hit container is"<<hits.size()<<"in module"<<i<<endl;
      for(unsigned int k=0;k<hits.size();k++) if(hits[k])
	{ 	 unsigned int module=hits[k]->module();	 unsigned int Eta=hits[k]->eta();	//float energyT=hits[k]->energy();
	int s=fabs(hits[k]->sub());int did(0);
	//  cout<<"module no"<<module<<"hit no"<<k<<"Eta is"<<Eta<<"Energy is"<<energyT<<endl; //it was 5 strips 02/11/05
	 if (module==fabs(mod) &&  Eta == fabs(eta))
	{ 	float energyT1=hits[k]->energy();
	mGeom2->getId(module,Eta,s,did);
	   // cout<<"his no"<<k<<"::module no::"<<module<<"::Eta is::"<<Eta<<"::Energy is::"<<energyT1<<endl; 
	
	if((dEdx>0.0000033&&dEdx<0.0000046)&&ch==1)hits[k]->setAdc(6011);
	if((dEdx>0.0000033&&dEdx<0.0000046)&&ch==-1)hits[k]->setAdc(6012);
	
	if((Fp>20&&dca<3&&dEdx>0.0000033&&dEdx<0.0000046)&&ch==1)hits[k]->setAdc(6013);
	if((Fp>20&&dca<3&&dEdx>0.0000033&&dEdx<0.0000046)&&ch==-1)hits[k]->setAdc(6014); 
	
	if(pTPC<1&&ch>=1) hits[k]->setAdc(5998);  if(pTPC<1&&ch<=-1) hits[k]->setAdc(5999);
	if(pTPC>1 &&ch>=1) hits[k]->setAdc(6001);if(pTPC>1 &&ch<=-1) hits[k]->setAdc(6002);
	if(pTPC>2 &&ch>=1) hits[k]->setAdc(6003);if(pTPC>2 &&ch<=-1) hits[k]->setAdc(6004);
	if(pTPC>3 &&ch>=1) hits[k]->setAdc(6005);if(pTPC>3 &&ch<=-1) hits[k]->setAdc(6006);
	if(pTPC>4 &&ch>=1) hits[k]->setAdc(6007);if(pTPC>4 &&ch<=-1) hits[k]->setAdc(6008); 
	if(pTPC>6 &&ch>=1) hits[k]->setAdc(6009);if(pTPC>6 &&ch<=-1) hits[k]->setAdc(6010);
		
	
	
	
	//if ((module==abs(mod) && Eta == abs(eta)))hits[k]->setAdc(6019);
	 	//cout<<"::::didE::"<<did<<":::energy::"<<energyT1<<"::module::"<<module<<"::adc is::"<<hits[k]->adc()<<"::dEdx::"<<dEdx<<"::ch:"<<ch<<"::pTPC::"<<pTPC<<endl;
      }
     }   
    }
 }   
  }
    }
    }               
 
  
delete pos;        
return kStOK;
}    


 Int_t StThirdMaker:: associateHadronsWithBEMC2()
{
   
   StEmcCollection* emcCollection = mStEvent->emcCollection(); //Get emcCollection
   if(!emcCollection) //return kStOK;
    {cout<<"There is no emcCollection"<<endl;       return kStWarn;}
  
   
  Double_t bFld;  StEventSummary* summary=mStEvent->summary();  if(summary){ bFld=summary->magneticField()/10.;}  if (fabs(bFld)<0.01){return kStWarn;}
 
 Int_t mod,eta,sub; StEmcPosition*  pos=new StEmcPosition(); StThreeVectorD position,momentum;
 
 StSPtrVecTrackNode& trackNodes=mStEvent->trackNodes(); StTrack* track;

 StEmcGeom * mGeom0=StEmcGeom::instance("bemc"); StEmcGeom * mGeom2=StEmcGeom::instance("bsmde");      StEmcGeom * mGeom3=StEmcGeom::instance("bsmdp");
  
 for(size_t nodeIndex=0;nodeIndex<trackNodes.size();nodeIndex++)
  {   size_t numberOfTracksInNode=trackNodes[nodeIndex]->entries(global);
 
   for( size_t trackIndex=0;trackIndex<numberOfTracksInNode;trackIndex++)
  { 
  track=trackNodes[nodeIndex]->track(global,trackIndex);   StTrackGeometry* geometry=track->geometry(); 
       short ch=geometry->charge();       float pTPC=geometry->momentum().mag();       
        float ptTPC=geometry->momentum().perp();             Int_t Fp=track->fitTraits().numberOfFitPoints(kTpcId);//if(Fp<20) continue;
     float dca=track->impactParameter(); //if(dca>3) continue;
     
      StPtrVecTrackPidTraits traits=track->pidTraits(kTpcId);       StDedxPidTraits* pid;
       for(UInt_t V=0;V<traits.size();V++)
       { pid=dynamic_cast<StDedxPidTraits*>(traits[V]);
          if(pid&&pid->method()==kTruncatedMeanId)
          { float dEdx=pid->mean();    //if(!(dEdx>0.0000033&&dEdx<0.0000046))  continue; if (!(track && track->flag()>=0)) continue;

                   Bool_t ok=pos->trackOnEmc(&position,&momentum,track,bFld);    if(! ok) continue;
   



  StEmcDetector* detector2 = emcCollection->detector(kBarrelSmdPhiStripId);  if (!detector2)
   {cout<<"There is no kBarrelSmdEtaStripId Detector"<<endl; return kStWarn;}
    mGeom3->getBin(position.phi(),position.pseudoRapidity(),mod,eta,sub);
    			
for(unsigned int i=1;i<=120;i++)
  {  
   if(fabs(mod)!=i) continue;
      StEmcModule* module=detector2->module(i);       StSPtrVecEmcRawHit& hits=module->hits();
	  //  cout<<"the size of hit container is"<<hits.size()<<"in module"<<i<<endl;
      for(unsigned int k=0;k<hits.size();k++) if(hits[k])
	{ 	 unsigned int module=hits[k]->module();	 unsigned int Eta=hits[k]->eta();	//float energyT=hits[k]->energy();
	int s=fabs(hits[k]->sub());int did(0);
	//  cout<<"module no"<<module<<"hit no"<<k<<"Eta is"<<Eta<<"Energy is"<<energyT<<endl; //it was 5 strips 02/11/05
	 if ( module==fabs(mod) && Eta == fabs(eta)&& fabs(s)==fabs(sub))
	{ 	float energyT1=hits[k]->energy();
	mGeom3->getId(module,Eta,s,did);
	
	if((dEdx>0.0000033&&dEdx<0.0000046)&&ch==1)hits[k]->setAdc(6011);
	if((dEdx>0.0000033&&dEdx<0.0000046)&&ch==-1)hits[k]->setAdc(6012);
	
	if((Fp>20&&dca<3&&dEdx>0.0000033&&dEdx<0.0000046)&&ch==1)hits[k]->setAdc(6013);
	if((Fp>20&&dca<3&&dEdx>0.0000033&&dEdx<0.0000046)&&ch==-1)hits[k]->setAdc(6014); 
	
	
	
	if(pTPC<1&&ch>=1) hits[k]->setAdc(5998);  if(pTPC<1&&ch<=-1) hits[k]->setAdc(5999);
	if(pTPC>1 &&ch>=1) hits[k]->setAdc(6001);if(pTPC>1 &&ch<=-1) hits[k]->setAdc(6002);
	if(pTPC>2 &&ch>=1) hits[k]->setAdc(6003);if(pTPC>2 &&ch<=-1) hits[k]->setAdc(6004);
	if(pTPC>3 &&ch>=1) hits[k]->setAdc(6005);if(pTPC>3 &&ch<=-1) hits[k]->setAdc(6006); //here was if(pTPC>2 &&ch<=-1)
	if(pTPC>4 &&ch>=1) hits[k]->setAdc(6007);if(pTPC>4 &&ch<=-1) hits[k]->setAdc(6008); 
	if(pTPC>6 &&ch>=1) hits[k]->setAdc(6009);if(pTPC>6 &&ch<=-1) hits[k]->setAdc(6010);
		
	
		  // cout<<"his no"<<k<<"::module no::"<<module<<"::Eta is::"<<Eta<<"::Energy is::"<<energyT1<<endl; 
	
	
	// if(ch==1) hits[k]->setAdc(6001);  if(ch==-1) hits[k]->setAdc(6002);
	//cout<<"::::didP::"<<did<<":::energy::"<<energyT1<<"::module::"<<module<<"::adc is::"<<hits[k]->adc()<<"::dEdx::"<<dEdx<<"::ch:"<<ch<<"::pTPC::"<<pTPC<<endl;
      }
     }   
    }   
    
  }   
  }
    }
    }               
 
  
delete pos;        
return kStOK;
}   
    
  Int_t StThirdMaker:: associateHadronsWithBEMC1()
{
   
   StEmcCollection* emcCollection = mStEvent->emcCollection(); //Get emcCollection
   if(!emcCollection) //return kStOK;
    {cout<<"There is no emcCollection"<<endl;       return kStWarn;}
  
   for(int RW=0; RW<4800;RW++) ParentTrackP[RW]=0.0;
for(int RW2=0; RW2<4800;RW2++) ParticleId[RW2]=0.0;
for(int RW3=0; RW3<4800;RW3++)ParentTrackPOSe[RW3]=0.0;
for(int RW4=0; RW4<4800;RW4++)ParentTrackPOSP[RW4]=0.0;

  Double_t bFld;  StEventSummary* summary=mStEvent->summary();  if(summary){ bFld=summary->magneticField()/10.;}  if (fabs(bFld)<0.01){return kStWarn;}
 
 Int_t mod,eta,sub; StEmcPosition*  pos=new StEmcPosition(); StThreeVectorD position,momentum;
 
 StSPtrVecTrackNode& trackNodes=mStEvent->trackNodes(); StTrack* track;

 StEmcGeom * mGeom0=StEmcGeom::instance("bemc"); StEmcGeom * mGeom2=StEmcGeom::instance("bsmde");      StEmcGeom * mGeom3=StEmcGeom::instance("bsmdp");
  
 for(size_t nodeIndex=0;nodeIndex<trackNodes.size();nodeIndex++)
  {   size_t numberOfTracksInNode=trackNodes[nodeIndex]->entries(global);
 
   for( size_t trackIndex=0;trackIndex<numberOfTracksInNode;trackIndex++)
  { 
  track=trackNodes[nodeIndex]->track(global,trackIndex);   StTrackGeometry* geometry=track->geometry(); 
       short ch=geometry->charge();       float pTPC=geometry->momentum().mag();       
        float ptTPC=geometry->momentum().perp();             Int_t Fp=track->fitTraits().numberOfFitPoints(kTpcId);//if(Fp<20) continue;
     float dca=track->impactParameter(); //if(dca>3) continue;
     
      StPtrVecTrackPidTraits traits=track->pidTraits(kTpcId);       StDedxPidTraits* pid;
       for(UInt_t V=0;V<traits.size();V++)
       { pid=dynamic_cast<StDedxPidTraits*>(traits[V]);
          if(pid&&pid->method()==kTruncatedMeanId)
          { float dEdx=pid->mean();    //if(!(dEdx>0.0000033&&dEdx<0.0000046))  continue; if (!(track && track->flag()>=0)) continue;

                   Bool_t ok=pos->trackOnEmc(&position,&momentum,track,bFld);    if(! ok) continue;
   
   
    
    

  StEmcDetector* detector3 = emcCollection->detector(kBarrelEmcTowerId);  if (!detector3)
   {cout<<"There is no kBarrelSmdEtaStripId Detector"<<endl; return kStWarn;}
    mGeom0->getBin(position.phi(),position.pseudoRapidity(),mod,eta,sub);
     float POSphi(0), POSeta(0); POSphi=position.phi();POSeta=position.pseudoRapidity();
    
    //   cout<<"the track coordinate at bsmde is"<<mod<<endl<<eta<<endl<<sub<<endl<<position.pseudoRapidity()<<endl<<position.phi()<<endl;                                
   // mDist->Fill(mod,eta);
			
for(unsigned int i=1;i<=120;i++)
  {  
   if(fabs(mod)!=i) continue;
      StEmcModule* module=detector3->module(i);       StSPtrVecEmcRawHit& hits=module->hits();
	  //  cout<<"the size of hit container is"<<hits.size()<<"in module"<<i<<endl;
      for(unsigned int k=0;k<hits.size();k++) if(hits[k])
	{ 	 unsigned int module=hits[k]->module();	 unsigned int Eta=hits[k]->eta();	//float energyT=hits[k]->energy();
	int s=fabs(hits[k]->sub()); int did(0);
	//  cout<<"module no"<<module<<"hit no"<<k<<"Eta is"<<Eta<<"Energy is"<<energyT<<endl; //it was 5 strips 02/11/05
	 if ( module==fabs(mod) && Eta == fabs(eta)&& fabs(s)==fabs(sub))
	{ 	float energyT1=hits[k]->energy(); 
	mGeom0->getId(module,Eta,s,did); 
	   // cout<<"his no"<<k<<"::module no::"<<module<<"::Eta is::"<<Eta<<"::Energy is::"<<energyT1<<endl; 
	   
	   if((dEdx>0.0000033&&dEdx<0.0000046)&&ch==1)hits[k]->setAdc(6011);
	if((dEdx>0.0000033&&dEdx<0.0000046)&&ch==-1)hits[k]->setAdc(6012);
	
	if((Fp>20&&dca<3&&dEdx>0.0000033&&dEdx<0.0000046)&&ch==1)hits[k]->setAdc(6013);
	if((Fp>20&&dca<3&&dEdx>0.0000033&&dEdx<0.0000046)&&ch==-1)hits[k]->setAdc(6014); 
	   
	   
	   
	if(pTPC<1&&ch>=1) hits[k]->setAdc(5998);  if(pTPC<1&&ch<=-1) hits[k]->setAdc(5999);
	if(pTPC>1 &&ch>=1) hits[k]->setAdc(6001);if(pTPC>1 &&ch<=-1) hits[k]->setAdc(6002);
	if(pTPC>2 &&ch>=1) hits[k]->setAdc(6003);if(pTPC>2 &&ch<=-1) hits[k]->setAdc(6004);
	if(pTPC>3 &&ch>=1) hits[k]->setAdc(6005);if(pTPC>3 &&ch<=-1) hits[k]->setAdc(6006);
	if(pTPC>4 &&ch>=1) hits[k]->setAdc(6007);if(pTPC>4 &&ch<=-1) hits[k]->setAdc(6008); 
	if(pTPC>6 &&ch>=1) hits[k]->setAdc(6009);if(pTPC>6 &&ch<=-1) hits[k]->setAdc(6010);
		
		
		 
	// if(ch==1) hits[k]->setAdc(6001);  if(ch==-1) hits[k]->setAdc(6002);//hits[k]->setEnergy(pTPC); 
//if(pTPC>4)cout<<"pTPC is"<<pTPC<<"::Energy is::"<<hits[k]->energy()<<"::did is::"<<did<<"::module::"<<module<<"::adc is::"<<hits[k]->adc()<<"::dEdx::"<<dEdx<<"::ch:"<<ch<<"::pTPC::"<<pTPC<<endl;
	int myiindex=did-1; ParentTrackP[myiindex]=pTPC; ParticleId[myiindex]=dEdx;ParentTrackPOSe[myiindex]=POSeta;
	ParentTrackPOSP[myiindex]=POSphi;
	
	
	
      }
     }   
    }      
    
 
   }   
  }
    }
    }               
 
  
delete pos;        
return kStOK;
}

//-----------------------------------------------------------HV change tower list------------------------------
void StThirdMaker::TowerHVcahngeList()
{
HVTowerid[0] =1;HVTowerid[1] =51;HVTowerid[2] =100;HVTowerid[3] =116;HVTowerid[4] =164;HVTowerid[5] =175;HVTowerid[6] =187;HVTowerid[7] =208;
HVTowerid[8] =221;HVTowerid[9] =234;HVTowerid[10] =315;HVTowerid[11] =341;HVTowerid[12] =353;HVTowerid[13] =371;HVTowerid[14] =424;
HVTowerid[15] =492;HVTowerid[16] =508;HVTowerid[17] =620;HVTowerid[18] =655;HVTowerid[19] =673;HVTowerid[20] =681;HVTowerid[21] =712;
HVTowerid[22] =719;HVTowerid[23] =749;HVTowerid[24] =758;HVTowerid[25] =768;HVTowerid[26] =772;HVTowerid[27] =779;HVTowerid[28] =790;
HVTowerid[29] =801;HVTowerid[30] =808;HVTowerid[31] =839;HVTowerid[32] =855;HVTowerid[33] =858;HVTowerid[34] =859;HVTowerid[35] =880;
HVTowerid[36] =897;HVTowerid[37] =900;HVTowerid[38] =974;HVTowerid[39] =975;HVTowerid[40] =981;HVTowerid[41] =1021;HVTowerid[42] =1023;
HVTowerid[43] =1026;HVTowerid[44] =1033;HVTowerid[45] =1041;HVTowerid[46] =1044;HVTowerid[47] =1045;HVTowerid[48] =1070;HVTowerid[49] =1091;
HVTowerid[50] =1118;HVTowerid[51] =1120;HVTowerid[52] =1128;HVTowerid[53] =1134;HVTowerid[54] =1140;HVTowerid[55] =1151;HVTowerid[56] =1170;
HVTowerid[57] =1174;HVTowerid[58] =1184;HVTowerid[59] =1194;HVTowerid[60] =1200;HVTowerid[61] =1218;HVTowerid[62] =1238;HVTowerid[63] =1241;
HVTowerid[64] =1242;HVTowerid[65] =1244;HVTowerid[66] =1258;HVTowerid[67] =1260;HVTowerid[68] =1263;HVTowerid[69] =1274;HVTowerid[70] =1280;
HVTowerid[71] =1303;HVTowerid[72] =1304;HVTowerid[73] =1312;HVTowerid[74] =1313;HVTowerid[75] =1319;HVTowerid[76] =1325;HVTowerid[77] =1338;
HVTowerid[78] =1340;HVTowerid[79] =1341;HVTowerid[80] =1348;HVTowerid[81] =1354;HVTowerid[82] =1375;HVTowerid[83] =1381;HVTowerid[84] =1385;
HVTowerid[85] =1387;HVTowerid[86] =1388;HVTowerid[87] =1394;HVTowerid[88] =1397;HVTowerid[89] =1407;HVTowerid[90] =1408;HVTowerid[91] =1427;
HVTowerid[92] =1440;HVTowerid[93] =1443;HVTowerid[94] =1487;HVTowerid[95] =1507;HVTowerid[96] =1608;HVTowerid[97] =1728;HVTowerid[98] =1755;
HVTowerid[99] =1772;HVTowerid[100] =1833;HVTowerid[101] =1838;HVTowerid[102] =1856;HVTowerid[103] =1861;HVTowerid[104] =1862;HVTowerid[105] =1864;
HVTowerid[106] =1866;HVTowerid[107] =1867;HVTowerid[108] =1872;HVTowerid[109] =1873;HVTowerid[110] =1880;HVTowerid[111] =1881;HVTowerid[112] =1882;
HVTowerid[113] =1886;HVTowerid[114] =1887;HVTowerid[115] =1891;HVTowerid[116] =1892;HVTowerid[117] =1893;HVTowerid[118] =1895;HVTowerid[119] =1899;
HVTowerid[120] =1900;HVTowerid[121] =1949;HVTowerid[122] =1984;HVTowerid[123] =2033;HVTowerid[124] =2058;HVTowerid[125] =2074;HVTowerid[126] =2075;
HVTowerid[127] =2079;HVTowerid[128] =2085;HVTowerid[129] =2094;HVTowerid[130] =2095;HVTowerid[131] =2105;HVTowerid[132] =2106;HVTowerid[133] =2127;
HVTowerid[134] =2133;HVTowerid[135] =2162;HVTowerid[136] =2168;HVTowerid[137] =2172;HVTowerid[138] =2184;HVTowerid[139] =2196;HVTowerid[140] =2223;
HVTowerid[141] =2243;HVTowerid[142] =2250;HVTowerid[143] =2254;HVTowerid[144] =2262;HVTowerid[145] =2303;HVTowerid[146] =2304;HVTowerid[147] =2306;
HVTowerid[148] =2352;HVTowerid[149] =2374;HVTowerid[150] =2387;HVTowerid[151] =2391;HVTowerid[152] =2392;HVTowerid[153] =2422;HVTowerid[154] =2439;
HVTowerid[155] =2458;HVTowerid[156] =2459;HVTowerid[157] =2470;HVTowerid[158] =2477;HVTowerid[159] =2580;HVTowerid[160] =2583;HVTowerid[161] =2624;
HVTowerid[162] =2658;HVTowerid[163] =2749;HVTowerid[164] =2794;HVTowerid[165] =2842;HVTowerid[166] =2929;HVTowerid[167] =3007;HVTowerid[168] =3025;
HVTowerid[169] =3097;HVTowerid[170] =3220;HVTowerid[171] =3231;HVTowerid[172] =3255;HVTowerid[173] =3259;HVTowerid[174] =3260;HVTowerid[175] =3261;
HVTowerid[176] =3287;HVTowerid[177] =3289;HVTowerid[178] =3290;HVTowerid[179] =3296;HVTowerid[180] =3305;HVTowerid[181] =3333;HVTowerid[182] =3403;
HVTowerid[183] =3452;HVTowerid[184] =3453;HVTowerid[185] =3473;HVTowerid[186] =3495;HVTowerid[187] =3500;HVTowerid[188] =3557;HVTowerid[189] =3602;
HVTowerid[190] =3620;HVTowerid[191] =3626;HVTowerid[192] =3640;HVTowerid[193] =3643;HVTowerid[194] =3653;HVTowerid[195] =3678;HVTowerid[196] =3692;
HVTowerid[197] =3712;HVTowerid[198] =3715;HVTowerid[199] =3723;HVTowerid[200] =3727;HVTowerid[201] =3759;HVTowerid[202] =3761;HVTowerid[203] =3769;
HVTowerid[204] =3803;HVTowerid[205] =3821;HVTowerid[206] =3838;HVTowerid[207] =3984;HVTowerid[208] =3986;HVTowerid[209] =4000;HVTowerid[210] =4013;
HVTowerid[211] =4019;HVTowerid[212] =4020;HVTowerid[213] =4021;HVTowerid[214] =4022;HVTowerid[215] =4023;HVTowerid[216] =4024;HVTowerid[217] =4025;
HVTowerid[218] =4026;HVTowerid[219] =4027;HVTowerid[220] =4028;HVTowerid[221] =4029;HVTowerid[222] =4030;HVTowerid[223] =4031;HVTowerid[224] =4032;
HVTowerid[225] =4033;HVTowerid[226] =4034;HVTowerid[227] =4035;HVTowerid[228] =4036;HVTowerid[229] =4037;HVTowerid[230] =4038;HVTowerid[231] =4039;
HVTowerid[232] =4040;HVTowerid[233] =4041;HVTowerid[234] =4042;HVTowerid[235] =4043;HVTowerid[236] =4044;HVTowerid[237] =4045;HVTowerid[238] =4046;
HVTowerid[239] =4047;HVTowerid[240] =4048;HVTowerid[241] =4049;HVTowerid[242] =4050;HVTowerid[243] =4051;HVTowerid[244] =4052;HVTowerid[245] =4053;
HVTowerid[246] =4054;HVTowerid[247] =4055;HVTowerid[248] =4056;HVTowerid[249] =4057;HVTowerid[250] =4058;HVTowerid[251] =4059;HVTowerid[252] =4060;
HVTowerid[253] =4061;HVTowerid[254] =4062;HVTowerid[255] =4063;HVTowerid[256] =4064;HVTowerid[257] =4065;HVTowerid[258] =4066;HVTowerid[259] =4067;
HVTowerid[260] =4068;HVTowerid[261] =4069;HVTowerid[262] =4070;HVTowerid[263] =4071;HVTowerid[264] =4072;HVTowerid[265] =4073;HVTowerid[266] =4074;
HVTowerid[267] =4075;HVTowerid[268] =4076;HVTowerid[269] =4077;HVTowerid[270] =4078;HVTowerid[271] =4079;HVTowerid[272] =4080;HVTowerid[273] =4081;
HVTowerid[274] =4082;HVTowerid[275] =4083;HVTowerid[276] =4084;HVTowerid[277] =4085;HVTowerid[278] =4086;HVTowerid[279] =4087;HVTowerid[280] =4088;
HVTowerid[281] =4089;HVTowerid[282] =4090;HVTowerid[283] =4091;HVTowerid[284] =4092;HVTowerid[285] =4093;HVTowerid[286] =4094;HVTowerid[287] =4095;
HVTowerid[288] =4096;HVTowerid[289] =4097;HVTowerid[290] =4098;HVTowerid[291] =4099;HVTowerid[292] =4100;HVTowerid[293] =4101;HVTowerid[294] =4102;
HVTowerid[295] =4103;HVTowerid[296] =4104;HVTowerid[297] =4105;HVTowerid[298] =4106;HVTowerid[299] =4107;HVTowerid[300] =4108;HVTowerid[301] =4109;
HVTowerid[302] =4110;HVTowerid[303] =4111;HVTowerid[304] =4112;HVTowerid[305] =4113;HVTowerid[306] =4114;HVTowerid[307] =4115;HVTowerid[308] =4116;
HVTowerid[309] =4117;HVTowerid[310] =4118;HVTowerid[311] =4119;HVTowerid[312] =4120;HVTowerid[313] =4121;HVTowerid[314] =4122;HVTowerid[315] =4123;
HVTowerid[316] =4124;HVTowerid[317] =4125;HVTowerid[318] =4126;HVTowerid[319] =4127;HVTowerid[320] =4128;HVTowerid[312] =4129;HVTowerid[322] =4131;
HVTowerid[323] =4132;HVTowerid[324] =4133;HVTowerid[325] =4134;HVTowerid[326] =4135;HVTowerid[327] =4136;HVTowerid[328] =4137;HVTowerid[329] =4138;
HVTowerid[330] =4139;HVTowerid[331] =4140;HVTowerid[332] =4141;HVTowerid[333] =4142;HVTowerid[334] =4143;HVTowerid[335] =4144;HVTowerid[336] =4145;
HVTowerid[337] =4146;HVTowerid[338] =4147;HVTowerid[339] =4148;HVTowerid[340] =4149;HVTowerid[341] =4150;HVTowerid[342] =4151;HVTowerid[343] =4152;
HVTowerid[344] =4153;HVTowerid[345] =4154;HVTowerid[346] =4155;HVTowerid[347] =4156;HVTowerid[348] =4157;HVTowerid[349] =4158;HVTowerid[350] =4159;
HVTowerid[351] =4160;HVTowerid[352] =4161;HVTowerid[353] =4162;HVTowerid[354] =4163;HVTowerid[355] =4164;HVTowerid[356] =4165;HVTowerid[357] =4166;
HVTowerid[358] =4167;HVTowerid[359] =4168;HVTowerid[360] =4169;HVTowerid[361] =4170;HVTowerid[362] =4172;HVTowerid[363] =4173;HVTowerid[364] =4174;
HVTowerid[365] =4175;HVTowerid[366] =4176;HVTowerid[367] =4177;HVTowerid[368] =4178;HVTowerid[369] =4179;HVTowerid[370] =4180;HVTowerid[371] =4223;
HVTowerid[372] =4240;HVTowerid[373] =4254;HVTowerid[374] =4339;HVTowerid[375] =4357;HVTowerid[376] =4377;HVTowerid[377] =4388;HVTowerid[378] =4426;
HVTowerid[379] =4459;HVTowerid[380] =4464;HVTowerid[381] =4505;HVTowerid[382] =4506;HVTowerid[383] =4507;HVTowerid[384] =4508;HVTowerid[385] =4514;
HVTowerid[386] =4543;HVTowerid[387] =4545;HVTowerid[388] =4546;HVTowerid[389] =4547;HVTowerid[390] =4560;HVTowerid[391] =4565;HVTowerid[392] =4566;
HVTowerid[393] =4567;HVTowerid[394] =4596;HVTowerid[395] =4660;HVTowerid[396] =4684;HVTowerid[397] =4765;HVTowerid[398] =4778;
 return;
}

//--------------------------------add it to .h

Float_t StThirdMaker::BBC_GetPhi(int iTile){
double twopi=2*3.14;
const float phi_div=twopi/12.;
//const float phi_div=Pi/6;
float bbc_phi=phi_div;
switch(iTile) {
case 0: bbc_phi=3*phi_div;
break;
case 1: bbc_phi=phi_div;
break;
case 2: bbc_phi=-1*phi_div;
break;
case 3: bbc_phi=-3*phi_div;
break;
case 4: bbc_phi=-5*phi_div;
break;
case 5: bbc_phi=5*phi_div;
break;
case 6: bbc_phi= (gRandom->Rndm()>0.5) ? 2*phi_div:4*phi_div;//Int_t g	= gRandom->Rndm(1) > 0.5 ? 0:1; 
break;
case 7: bbc_phi=3*phi_div;
break;
case 8: bbc_phi=phi_div;
break;
case 9: bbc_phi=0.;
break;
case 10: bbc_phi=-phi_div;
break;
case 11: bbc_phi=(gRandom->Rndm()>0.5) ? -2*phi_div:-4*phi_div;
break;
case 12: bbc_phi=-3*phi_div;
break;
case 13: bbc_phi=-5*phi_div;
break;
case 14: bbc_phi=twopi/2.0;
break;
case 15: bbc_phi=5*phi_div;
break;
}
if(bbc_phi<0.0) bbc_phi +=2*3.14;
if(bbc_phi>2*3.14) bbc_phi -=2*3.14;
return bbc_phi;
}
/*
/Constant
 Double_t StThirdMaker::bbcGainFac[2][16] = {{0.962841,0.979108,1.03717,0.957205,1.07229,0.991379,0.730763,1.05997,1.04575,1.57317,1.02486,0.80374,0.919484,1.02191,1.35782,0.928025},{0.930545,1.01045,0.988749,1.03118,0.943601,1.09548,0.771513,1.00259,0.972948,1.48775,0.91112,0.772221,0.985141,1.12275,1.42478,1.00545}}; 
 Float_t StThirdMaker::bbcPedstal[2][16]={{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}};
 bbc_nPsiBins      = 64
typedef Double_t BBC_PsiWgt_t[64]; 
//  static Float_t bbcPedstal[2][24]; //all tiles working
//  static Double_t bbcGainFac[2][24]; //all tiles working
  static Float_t  bbcPedstal[2][16];
  static Double_t bbcGainFac[2][16];
  
/StFlowCutEvent

//----------------------------------------------------------------------

Float_t StFlowEvent::BBC_PsiEst() {
  //psi angle from east BBC

  Float_t eXsum=0., eYsum=0., eWgt=0.;

 for(int iTile = 0; iTile < 16; iTile++) {
      eXsum += cos(BBC_GetPhi(iTile))*BBCAdc(0,iTile);
      eYsum += sin(BBC_GetPhi(iTile))*BBCAdc(0,iTile);
      eWgt += BBCAdc(0,iTile);
    }
   
  Float_t psi_e = atan2((eWgt>0.) ? eYsum/eWgt:0.,(eWgt>0.) ? eXsum/eWgt:0.);

  if (psi_e < 0.) { psi_e += twopi; }
    return psi_e;
}

//----------------------------------------------------------------------

Float_t StFlowEvent::BBC_PsiWst() {
  //psi angle from west BBC

 Float_t wXsum=0.,wYsum=0.,wWgt=0.;

// Yadav for BBC

    for(int iTile = 0; iTile <16; iTile++) {
      wXsum += cos(BBC_GetPhi(iTile))*BBCAdc(1,iTile);
      wYsum += sin(BBC_GetPhi(iTile))*BBCAdc(1,iTile);
      wWgt += BBCAdc(1,iTile);

         }

 Float_t psi_w = atan2((wWgt>0.) ? wYsum/wWgt:0.,(wWgt>0.) ? wXsum/wWgt:0.);

  if (psi_w < 0.) { psi_w += twopi; }
  return psi_w;
}
//------------------------------------------------------

//Cut for BBC Adc Saturation ::Added by Yadav 
  StBbcTriggerDetector &BBC=pMuEvent->bbcTriggerDetector();

 for (int iTile= 0; iTile<16; iTile++) {
           if(BBC.adc(iTile) > 245 || BBC.adc(iTile+24)>245) {
         mBBCAdcCut++ ;
        return kFALSE;
      }
  }
   static UInt_t  mBBCAdcCut;             // number of events cut by BBCAdc
StFlowEvent.cxx
TRandom rndm;
const float phi_div=twopi/12.;
Float_t StFlowEvent::mBBCshiftEast_cos[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftEast_sin[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftWest_cos[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftWest_sin[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftFul0_cos[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftFul0_sin[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftFull_cos[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
Float_t StFlowEvent::mBBCshiftFull_sin[8] = {0.,0.,0.,0.,0.,0.,0.,0.};
*/
