#include "GammaJetTrack.h"

//Contains the Member Functions for the TRACK Class                             

ClassImp(GammaJetTrack);

//_____________________________________________________________________         
GammaJetTrack::GammaJetTrack():TObject()
			      , nHitsFit(-999.)
			      , nHitsPoss(-999.)
			      , trackFlag(-999.)
			      , pZ(-999.)
			      , pX(-999.)
			      , pY(-999.)
			      , pT(-999.)
			      , dEdx(-999.)
			      , charge(-999.)
			      , tofBeta(-999.)
			      , eta(-999.)
			      , phi(-999.)
			      , nSigElectron(-999.)//
			      , nSigPion(-999.)
			      , nSigKaon(-999.)
			      , nSigProton(-999.)
			      , dcag(-999.)
			      , nHits(-999.)
			      , dEdxHits(-999.)
			      , firstZPoint(-999.)
  , lastZPoint(-999.)
  , tofSigElectron(-999.)
  , tofSigPion(-999.)
  , tofSigKaon(-999.)
  , tofSigProton(-999.)
  , timeOfflight(-999.)
  ,pathLength(-999.)
  ,trkIndex(-999.)  //track index starts from 1,...... (check  trkIndex_counter in StJetTreeThirdMaker.cxx)
{


}

//_____________________________________________________________________         
GammaJetTrack::~GammaJetTrack(){

  //GammaJetTrack Destructor                                                            

}

//_____________________________________________________________________         
//GammaJetTrack::GammaJetTrack(Float_t *trackarray){
GammaJetTrack::GammaJetTrack(GammaJetTrack  &trackarray)
  :TObject(trackarray)
  , nHitsFit(trackarray.nHitsFit)
  , nHitsPoss(trackarray.nHitsPoss)
  , trackFlag(trackarray.trackFlag)
  , pZ(trackarray.pZ)
  , pX(trackarray.pX)
  , pY(trackarray.pY)
  , pT(trackarray.pT)
  , dEdx(trackarray.dEdx)
  , charge(trackarray.charge)
  , tofBeta(trackarray.tofBeta)
  , eta(trackarray.eta)
  , phi(trackarray.phi)
  , nSigElectron(trackarray.nSigElectron)//
  , nSigPion(trackarray.nSigPion)
  , nSigKaon(trackarray.nSigKaon)
  , nSigProton(trackarray.nSigProton)
  , dcag(trackarray.dcag)
  , nHits(trackarray.nHits)
  , dEdxHits(trackarray.dEdxHits)
  , firstZPoint(trackarray.firstZPoint)
  , lastZPoint(trackarray.lastZPoint)
  , tofSigElectron(trackarray.tofSigElectron)
  , tofSigPion(trackarray.tofSigPion)
  , tofSigKaon(trackarray.tofSigKaon)
  , tofSigProton(trackarray.tofSigProton)
  , timeOfflight(trackarray.timeOfflight)
  ,pathLength(trackarray.pathLength)
  ,trkIndex(trackarray.trkIndex)  
{


}

//_____________________________________________________________________         
void GammaJetTrack::SetTrackAttributes(Float_t *trackarray){
  /*
  //Sets All the Track Attributes - ORDER MATTERS!                              
  nHitsFit       = (Int_t)trackarray[0];
  nHitsPoss      = (Int_t)trackarray[1];
  trackFlag      = (Int_t)trackarray[2];
  pX             = (Double_t)trackarray[3];
  pY             = (Double_t)trackarray[4];
  pZ             = (Double_t)trackarray[5];
  pT             = (Double_t)trackarray[6];
  dEdx           = (Double_t)trackarray[7];
  charge         = (Double_t)trackarray[8];
  tofBeta        = (Double_t)trackarray[9];
  eta            = (Double_t)trackarray[10];
  phi            = (Double_t)trackarray[11];
  nSigElectron   = (Double_t)trackarray[12];
  nSigPion       = (Double_t)trackarray[13];
  nSigKaon       = (Double_t)trackarray[14];
  nSigProton     = (Double_t)trackarray[15];
  dcag           = (Double_t)trackarray[16];
  nHits          = (Double_t)trackarray[17];
  dEdxHits       = (Double_t)trackarray[18];
  firstZPoint    = (Double_t)trackarray[19];
  lastZPoint     = (Double_t)trackarray[20];
  tofSigElectron = (Double_t)trackarray[21];
  tofSigPion     = (Double_t)trackarray[22];
  tofSigKaon     = (Double_t)trackarray[23];
  tofSigProton   = (Double_t)trackarray[24];
  pathLength     = (Double_t)trackarray[25];
  timeOfflight = (Double_t)trackarray[26];
  trkIndex     = (Int_t)trackarray[27];
  */
}

//_____________________________________________________________________         
void GammaJetTrack::PrintTrack(){

}
