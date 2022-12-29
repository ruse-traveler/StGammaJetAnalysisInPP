#ifndef GammaJetTrack_H
#define GammaJetTrack_H

//Contains the GammaJetTrack Class Declaration                                          

#include <stdlib.h>
#include <TROOT.h>
#include "Riostream.h"


class GammaJetTrack:public TObject{

 private:
  Int_t nHitsFit;
  Int_t nHitsPoss;
  Int_t trackFlag;
  Int_t pdgId;  // only for MC StMcTrack
  Int_t geantId;  // only for MC StMcTrack
  Double_t pZ;
  Double_t pX;
  Double_t pY;
  Double_t pT;
  Double_t dEdx;
  Double_t charge;
  Double_t tofBeta;
  Double_t eta;
  Double_t phi;
  Double_t nSigElectron;//
  Double_t nSigPion;
  Double_t nSigKaon;
  Double_t nSigProton;
  Double_t dcag;
  Double_t nHits;
  Double_t dEdxHits;
  Double_t  firstZPoint;
  Double_t  lastZPoint;
  Double_t tofSigElectron;
  Double_t tofSigPion;
  Double_t tofSigKaon;
  Double_t tofSigProton;
  Double_t timeOfflight;
  Double_t pathLength;
  Int_t    trkIndex;  //track index starts from 1,...... (check  trkIndex_counter in StJetTreeThirdMaker.cxx)

 public:
  GammaJetTrack();
  //  GammaJetTrack(Float_t *trackArray);
  GammaJetTrack(GammaJetTrack  &trackarray);
  virtual ~GammaJetTrack();
  void Clear(Option_t *Option = "");
  
  void SetTrackAttributes(Float_t *trackarray);
  void PrintTrack();


  Int_t GetnHitsFit(){return nHitsFit;}
  void  SetnHitsFit(Int_t val){ nHitsFit = val;}
  
  Int_t GetnHitsPoss(){return nHitsPoss;}
  void  SetnHitsPoss(Int_t val){ nHitsPoss =val;}

  Int_t GetTrackFlag(){return trackFlag;}
  void  SetTrackFlag(Int_t val){ trackFlag = val;}

  Int_t GetPdgId(){return pdgId;}
  void  SetPdgId(Int_t val){ pdgId = val;}

  Int_t GetGeantId(){return geantId;}
  void  SetGeantId(Int_t val){ geantId = val;}

  Double_t GetpZ(){return pZ;}
  void     SetpZ(Double_t val){pZ = val;}

  Double_t GetpT(){return pT;}
  void     SetpT(Double_t val){ pT = val;}

  Double_t GetpX(){return pX;}
  void SetpX(Double_t val){ pX = val;}

  Double_t GetpY(){return pY;}
  void     SetpY(Double_t val){ pY = val;}


  Double_t GetdEdx(){return dEdx;}
  void     SetdEdx(Double_t val){ dEdx = val;}

  Double_t GetCharge(){return charge;}
  void      SetCharge(Double_t val){charge = val;}

  Double_t GetTOFBeta(){return tofBeta;}
  void     SetTOFBeta(Double_t val){ tofBeta = val;}


  Double_t GetEta(){return eta;}
  void  SetEta(Double_t val){ eta = val;}


  Double_t GetPhi(){return phi;}
  void SetPhi(Double_t val){ phi = val;}


  Double_t GetnSigElectron(){return nSigElectron;}
  void     SetnSigElectron(Double_t val){ nSigElectron = val;}

  Double_t GetnSigPion(){return nSigPion;}
  void     SetnSigPion(Double_t val){ nSigPion = val;}


  Double_t GetnSigKaon(){return nSigKaon;}
  void SetnSigKaon(Double_t val){ nSigKaon = val;}

  Double_t GetnSigProton(){return nSigProton;}
  void     SetnSigProton(Double_t val){ nSigProton = val;}

  Double_t GetDCAg(){return dcag;}
  void     SetDCAg(Double_t val){dcag = val;}

  Double_t GetnHits(){return nHits;}
  void     SetnHits(Double_t val){ nHits = val;}


  Double_t GetdEdxHits(){return dEdxHits;}
  void     SetdEdxHits(Double_t val){ dEdxHits = val;}

  Double_t GetFirstZPoint(){return firstZPoint;}
  void     SetFirstZPoint(Double_t val){ firstZPoint = val;}

  Double_t GetLastZPoint(){return lastZPoint;}
  void     SetLastZPoint(Double_t val){ lastZPoint = val;}

  Double_t GetTOFSigElectron(){return tofSigElectron;}
  void     SetTOFSigElectron(Double_t val){ tofSigElectron = val;}

  Double_t GetTOFSigPion(){return tofSigPion;}
  void     SetTOFSigPion(Double_t val){tofSigPion = val;}

  Double_t GetTOFSigKaon(){return tofSigKaon;}
  void     SetTOFSigKaon(Double_t val){ tofSigKaon = val;}

  Double_t GetTOFSigProton(){return tofSigProton;}
  void     SetTOFSigProton(Double_t val){tofSigProton = val;}

  Double_t GetPathLength(){return pathLength;}
  void     SetPathLength(Double_t val){ pathLength = val;}

  Double_t GettimeOfflight(){return timeOfflight;}
  void     SettimeOfflight(Double_t val){ timeOfflight = val;}

  Int_t    GettrkIndex(){return trkIndex;}
  void    SettrkIndex(Double_t val){ trkIndex = val;}


  ClassDef(GammaJetTrack,1);

};

#endif
