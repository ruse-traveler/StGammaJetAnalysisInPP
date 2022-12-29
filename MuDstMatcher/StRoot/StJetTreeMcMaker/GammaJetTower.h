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
  Float_t TwrPed;
  Float_t TwrRMS;


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
  // For E/p calculation [Derek, 07.20.2017]
  Float_t       fMatchedTracksArray_nSigPi[max_matchTrk];
  Float_t       fMatchedTracksArray_nSigK[max_matchTrk];
  Float_t       fMatchedTracksArray_nSigP[max_matchTrk];
  Float_t       fMatchedTracksArray_nSigE[max_matchTrk];
  Float_t       fMatchedTracksArray_dcag[max_matchTrk];
  Float_t       fMatchedTracksArray_eta[max_matchTrk];
  Float_t       fMatchedTracksArray_pT[max_matchTrk];
  Int_t         fMatchedTracksArray_nFit[max_matchTrk];
  Int_t         fMatchedTracksArray_nPos[max_matchTrk];
  //  vector<int>   fMatchedTracksArray_; 



 public:
  GammaJetTower();
  GammaJetTower(GammaJetTower &towerArray);
  virtual ~GammaJetTower();

  void Clear(Option_t *Option = "");
  
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

  Float_t GetTwrPed(){return TwrPed;}
  void SetTwrPed(Float_t val){TwrPed = val;}

  Float_t GetTwrRMS(){return TwrRMS;}
  void SetTwrRMS(Float_t val){TwrRMS = val;}


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
  */

  Float_t GetTwrMatchSumP(){return TwrMatchP;}
  void SetTwrMatchSumP(double val){ TwrMatchP = val;}
  

  void AddMatchedTrack(Int_t idx);
  void AddMatchedTrackIndex(Int_t idx) {AddMatchedTrack(idx);}

  Int_t GetNAssocTracks() const {return fNAssocTracks;}
  void  SetNAssocTracks(Int_t val) {fNAssocTracks = val;}


  void SetMatchedTracksArray(Int_t val, Int_t idx){ fMatchedTracksArray_[idx] = val;}
  void SetMatchedTracksArray_P(Float_t val, Int_t idx){fMatchedTracksArray_P[idx] = val;}

  // for E/p calculation [Derek, 07.20.2017]
  void SetMatchedTracksArray_nSigPi(Float_t val, Int_t idx){fMatchedTracksArray_nSigPi[idx] = val;}
  void SetMatchedTracksArray_nSigK(Float_t val, Int_t idx){fMatchedTracksArray_nSigK[idx] = val;}
  void SetMatchedTracksArray_nSigP(Float_t val, Int_t idx){fMatchedTracksArray_nSigP[idx] = val;}
  void SetMatchedTracksArray_nSigE(Float_t val, Int_t idx){fMatchedTracksArray_nSigE[idx] = val;}
  void SetMatchedTracksArray_dcag(Float_t val, Int_t idx){fMatchedTracksArray_dcag[idx] = val;}
  void SetMatchedTracksArray_eta(Float_t val, Int_t idx){fMatchedTracksArray_eta[idx] = val;}
  void SetMatchedTracksArray_pT(Float_t val, Int_t idx){fMatchedTracksArray_pT[idx] = val;}
  void SetMatchedTracksArray_nFit(Int_t val, Int_t idx){fMatchedTracksArray_nFit[idx] = val;}
  void SetMatchedTracksArray_nPos(Int_t val, Int_t idx){fMatchedTracksArray_nPos[idx] = val;}
  
  //  Int_t   GetMatchedTrackIndex(Int_t idx) {return fMatchedTracks.At(idx);}

  //  const TArrayI *GetMatchedTracks() {return &fMatchedTracks;}
  //  const TArrayI *GetMatchedTrackIndexes() {return &fMatchedTracks;}
    


  ClassDef(GammaJetTower,1);

};

#endif
