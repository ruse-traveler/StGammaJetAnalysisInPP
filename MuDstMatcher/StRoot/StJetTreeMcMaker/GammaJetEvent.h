//___________________________
//  Event class to produce  picodst for GammaJet
//  May 30, 2015  :   Nihar R. Sahoo
//_____________________________


#ifndef GammaJetEvent_H
#define GammaJetEvent_H

#include <TROOT.h>
#include "Riostream.h"
#include <TClonesArray.h>
#include "GammaJetTrack.h"
#include "GammaJetTower.h"

//Event Class Definition                                                        
class GammaJetEvent:public TObject{

 private:
  Long64_t runNumber;
  Long64_t eventNumber;
  Int_t trigID;
  Int_t nGlobalTracks;
  Int_t nPrimaryTracks;
  Int_t refMult;
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
  Double_t refMultPos; //
  Double_t refMultNeg; //
  Double_t bTOFTrayMultiplicity;//
  Int_t    nVerticies;//no. of Vertex
  Double_t MagF;  // Magnetic field
  Double_t VrtxRank;  // Vertex ranking

  
  //_______TSP information for event
  Float_t Etsp;
  Int_t ETwrdidT; // Tower Id
  Int_t ETwradc11;  // tower projection adc

  Float_t ETwreneT0; // Tower energy

  Float_t ETwreT; // Tower Eta true
  Float_t ETwrENET0;
  Float_t ETwrphT; // Tower phi true
  Float_t ETwrPTower; // projected track to this twr momentum
  Float_t ETwrpidTower; // projected track to this twr PID (dedx)
  Int_t ETwrmoduleT; // module no of this twr


  Float_t EClustEneT0;   // Cluster Energy
  Float_t EClustetav1; // Cluster Eta
  Float_t EClustphiv1; // Cluster phi



  Float_t EEstrpen01;
  Float_t EEstrpen02;
  Float_t EEstrpen03;

  Float_t EEstrpen0;
  Float_t EEstrpen1;
  Float_t EEstrpen2;
  Float_t EEstrpen3;
  Float_t EEstrpen4;
  Float_t EEstrpen5;
  Float_t EEstrpen6;
  Float_t EEstrpen7;
  Float_t EEstrpen8;
  Float_t EEstrpen9;
  Float_t EEstrpen10;
  Float_t EEstrpen11;
  Float_t EEstrpen12;
  Float_t EEstrpen13;
  Float_t EEstrpen14;
  Float_t EEstrpen15; //eta strip energy from clustering


  Int_t ETwrdidE; // Eta strip Id

  //_______________ Phi Strip Info
  Float_t EPstripenp01;
  Float_t EPstripenp02;
  Float_t EPstripenp03;


  Float_t EPstripenp0;
  Float_t EPstripenp1;
  Float_t EPstripenp2;
  Float_t EPstripenp3;
  Float_t EPstripenp4;
  Float_t EPstripenp5;
  Float_t EPstripenp6;
  Float_t EPstripenp7;
  Float_t EPstripenp8;
  Float_t EPstripenp9;
  Float_t EPstripenp10;
  Float_t EPstripenp11;
  Float_t EPstripenp12;
  Float_t EPstripenp13;
  Float_t EPstripenp14;
  Float_t EPstripenp15; //eta strip energy from clustering


  Float_t EclustEnnq1;
  Float_t EclustEnnq20;
  Float_t EclustEnnq19;

  Float_t EclustEnpq1;
  Float_t EclustEnpq20;
  Float_t EclustEnpq19;
  Float_t EclustEnpq21;


  //_________________________  

  TClonesArray *PrimaryTrackArray;
  TClonesArray *TowerArray;

 public:
  GammaJetEvent();
  GammaJetEvent(Float_t *eventarray);
  virtual ~GammaJetEvent();

  void SetEventAttributes(Float_t *eventarray);
  void PrintEvent(Bool_t printtrack);
  void ResetEvent();
  //void AddTrackArray(TClonesArray *trackarray);
  void AddTrack(Float_t *trackarray, Int_t At);
  void AddTrack(GammaJetTrack *trk, Int_t At);
  //  void AddTower(Float_t *towerarray, Int_t idx);
  void AddTower(GammaJetTower *t, Int_t idx);

  //Gets                                                                        
  Long64_t GetRunNumber(){return runNumber;}
  Long64_t GetEventNumber(){return eventNumber;}
  Int_t    GetTrigID(){return trigID;}
  Int_t    GetNGlobalTracks(){return nGlobalTracks;}
  Int_t    GetNPrimaryTracks(){return nPrimaryTracks;}
  Int_t    GetRefMult(){return refMult;}
  Int_t    GetRefMultPos(){return refMultPos;}
  Int_t    GetRegMultNeg(){return refMultNeg;}
  Int_t    GetnVerticies(){return nVerticies;}
  Double_t GetVPDVz(){return vpdVz;}
  Double_t GetXVertex(){return xVertex;}
  Double_t GetYVertex(){return yVertex;}
  Double_t GetZVertex(){return zVertex;}
  Double_t GetBBCZVertex(){return bbcZVertex;}
  Double_t GetZDCCoincidenceRate(){return zdcCoincidenceRate;}
  Double_t GetBBCCoincidenceRate(){return bbcCoincidenceRate;}
  Double_t GetBackgroundRate(){return backgroundRate;}
  Double_t GetBBCBlueBackgroundRate(){return bbcBlueBackgroundRate;}
  Double_t GetBBCYellowBackgroundRate(){return bbcYellowBackgroundRate;}
  Double_t GetbTOFTrayMultiplicity(){return bTOFTrayMultiplicity;}
  Double_t GetMagneticFiled(){return MagF;}
  Double_t GetVrtxRank(){return VrtxRank;}

  //_________TSP info______
  Float_t GetEtsp(){return Etsp;}
  Int_t GetETwrdidT(){return ETwrdidT;} // Tower Id
  Int_t GetETwradc11(){return ETwradc11;}  // tower projection adc

  Float_t GetETwreneT0(){return ETwreneT0;} // Tower energy

  Float_t GetETwreT(){return ETwreT;} // Tower Eta true
  Float_t GetETwrENET0(){return ETwrENET0;}
  Float_t GetETwrphT(){return ETwrphT;} // Tower phi true
  Float_t GetETwrPTower(){return ETwrPTower;} // projected track to this twr momentum
  Float_t GetETwrpidTower(){return ETwrpidTower;} // projected track to this twr PID (dedx)
  Int_t GetETwrmoduleT(){return ETwrmoduleT;} // module no of this twr


  Float_t GetEClustEneT0(){return EClustEneT0;}   // Cluster Energy
  Float_t GetEClustetav1(){return EClustetav1;} // Cluster Eta
  Float_t GetEClustphiv1(){return EClustphiv1;} // Cluster phi

  Float_t GetEEstrpen01(){return  EEstrpen01;}
  Float_t GetEEstrpen02(){return  EEstrpen02;}
  Float_t GetEEstrpen03(){return  EEstrpen03;}

  Float_t GetEEstrpen0(){return  EEstrpen0;}
  Float_t GetEEstrpen1(){return  EEstrpen1;}
  Float_t GetEEstrpen2(){return  EEstrpen2;}
  Float_t GetEEstrpen3(){return  EEstrpen3;}
  Float_t GetEEstrpen4(){return  EEstrpen4;}
  Float_t GetEEstrpen5(){return  EEstrpen5;}
  Float_t GetEEstrpen6(){return  EEstrpen6;}
  Float_t GetEEstrpen7(){return  EEstrpen7;}
  Float_t GetEEstrpen8(){return  EEstrpen8;}
  Float_t GetEEstrpen9(){return  EEstrpen9;}
  Float_t GetEEstrpen10(){return  EEstrpen10;}
  Float_t GetEEstrpen11(){return  EEstrpen11;}
  Float_t GetEEstrpen12(){return  EEstrpen12;}
  Float_t GetEEstrpen13(){return  EEstrpen13;}
  Float_t GetEEstrpen14(){return  EEstrpen14;}
  Float_t GetEEstrpen15(){return  EEstrpen15;} //eta strip energy from clustering


  Int_t GetETwrdidE(){return  ETwrdidE;} // Eta strip Id

  //_______________ Phi Strip Info
  Float_t GetEPstripenp01(){return  EPstripenp01;}
  Float_t GetEPstripenp02(){return  EPstripenp02;}
  Float_t GetEPstripenp03(){return  EPstripenp03;}


  Float_t GetEPstripenp0(){return  EPstripenp0;}
  Float_t GetEPstripenp1(){return  EPstripenp1;}
  Float_t GetEPstripenp2(){return  EPstripenp2;}
  Float_t GetEPstripenp3(){return  EPstripenp3;}
  Float_t GetEPstripenp4(){return  EPstripenp4;}
  Float_t GetEPstripenp5(){return  EPstripenp5;}
  Float_t GetEPstripenp6(){return  EPstripenp6;}
  Float_t GetEPstripenp7(){return  EPstripenp7;}
  Float_t GetEPstripenp8(){return  EPstripenp8;}
  Float_t GetEPstripenp9(){return  EPstripenp9;}
  Float_t GetEPstripenp10(){return  EPstripenp10;}
  Float_t GetEPstripenp11(){return  EPstripenp11;}
  Float_t GetEPstripenp12(){return  EPstripenp12;}
  Float_t GetEPstripenp13(){return  EPstripenp13;}
  Float_t GetEPstripenp14(){return  EPstripenp14;}
  Float_t GetEPstripenp15(){return  EPstripenp15;} //eta strip energy from clustering


  Float_t GetEclustEnnq1(){return  EclustEnnq1;}
  Float_t GetEclustEnnq20(){return  EclustEnnq20;}
  Float_t GetEclustEnnq19(){return  EclustEnnq19;}

  Float_t GetEclustEnpq1(){return  EclustEnpq1;}
  Float_t GetEclustEnpq20(){return  EclustEnpq20;}
  Float_t GetEclustEnpq19(){return  EclustEnpq19;}
  Float_t GetEclustEnpq21(){return  EclustEnpq21;}




  //TRACK *GetTrack(Int_t i);
  Int_t GetTrackArrayEntries();
  Int_t GetTowerArrayEntries();



  ClassDef(GammaJetEvent,1);

};

#endif
