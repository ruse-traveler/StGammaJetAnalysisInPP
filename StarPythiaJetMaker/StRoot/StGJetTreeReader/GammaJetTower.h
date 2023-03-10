#ifndef GammaJetTower_H
#define GammaJetTower_H

//Contains the GammaJetTower Class Declaration                                          

#include <stdlib.h>
#include <TROOT.h>
#include "Riostream.h"
#include  <TArrayI.h>
#include  <vector>

const Int_t max_matchTrk =10;

class GammaJetTower:public TObject{

 private:



  Int_t TwrId;
  Float_t TwrEng;
  Float_t TwrEta;
  Float_t TwrPhi;
  Float_t TwrADC;


  //____matching Info
  Int_t TwrMatchIdnex; // info about wheather matched or not ----> Yes: 1  and   No: 0
  Int_t  NoOfmatchedTrk;   // no. of matched track per Twr
  //  Float_t TwrMatchDphi;  // matched tracks Delta Phi = track phi - twr phi  .... (careful abt this cut)
  //  Float_t TwrMatchDEta;  // matched tracks Delta Eta = track Eta - twr Eta   .....(same.....,,////////)
  Float_t TwrMatchP;  // matched tracks momentum
  //  TArrayI TwrMatchP;  // matched tracks momentum
  //  Int_t   TwrMatchIndexOfTrack; // to know which track index from primary track array it is

  //____Tower Momentum Info (Conversion check in GammaJetTowerUtil class)
  Float_t TwrPx;
  Float_t TwrPy;
  Float_t TwrPz;

  Int_t           fNAssocTracks; 
  //TArrayI         fMatchedTracks; 

  Int_t         fMatchedTracksArray_[max_matchTrk]; 
  Float_t         fMatchedTracksArray_P[max_matchTrk]; 
  //  vector<int>   fMatchedTracksArray_; 



 public:
  GammaJetTower();
  GammaJetTower(GammaJetTower &towerArray);
  virtual ~GammaJetTower();
  
  void SetTowerAttributes(Float_t &towerarray);
  void PrintTower();

  Float_t GetTwrId(){return TwrId;}
  void  SetTwrId(Float_t val){TwrId = val;}

  Float_t GetTwrEng(){return TwrEng;}
  void SetTwrEng(Float_t val){ TwrEng= val;}


  Float_t GetTwrEta(){return TwrEta;}
  void SetTwrEta(Float_t val){ TwrEta = val;}

  Float_t GetTwrPhi(){return TwrPhi;}
  void SetTwrPhi(Float_t val){ TwrPhi = val;}

  Float_t GetTwrADC(){return TwrADC;}
  void SetTwrADC(Float_t val){TwrADC = val;}


  Float_t GetTwrPx(){return TwrPx;}
  void SetTwrPx(Float_t val){TwrPx =val;}

  Float_t GetTwrPy(){return TwrPy;}
  void SetTwrPy(Float_t val){TwrPy=val;}

  Float_t GetTwrPz(){return TwrPz;}
  void SetTwrPz(Float_t val){ TwrPz=val;}
  
  Int_t GetTwrMatchIdnex(){return TwrMatchIdnex;}
  void SetTwrMatchIdnex(Int_t val){TwrMatchIdnex=val;}

  Int_t GetNoOfmatchedTrk(){return NoOfmatchedTrk;}
  void SetNoOfmatchedTrk(Int_t val){ NoOfmatchedTrk = val;}

  /*
  Float_t GetTwrMatchDphi(){return TwrMatchDphi;}
  //  void SetTwrMatchDphi(Float_t val){ TwrMatchDphi = val;}
  
  Float_t GetTwrMatchDEta(){return TwrMatchDEta;}


  //TArrayI GetTwrMatchP(){return TwrMatchP;}
  /*Float_t GetTwrMatchP(){return TwrMatchP;}
  Int_t GetTwrMatchIndexOfTrack(){return TwrMatchIndexOfTrack;}
  */

  void AddMatchedTrack(Int_t idx);
  void AddMatchedTrackIndex(Int_t idx) {AddMatchedTrack(idx);}

  Int_t GetNAssocTracks() const {return fNAssocTracks;}
  void  SetNAssocTracks(Int_t val) {fNAssocTracks = val;}


  void SetMatchedTracksArray(Int_t val, Int_t idx){ fMatchedTracksArray_[idx] = val;}
  void SetMatchedTracksArray_P(Float_t val, Int_t idx){fMatchedTracksArray_P[idx] = val;}
  
  //  Int_t   GetMatchedTrackIndex(Int_t idx) {return fMatchedTracks.At(idx);}

  //  const TArrayI *GetMatchedTracks() {return &fMatchedTracks;}
  //  const TArrayI *GetMatchedTrackIndexes() {return &fMatchedTracks;}
    


  ClassDef(GammaJetTower,1);

};

#endif
