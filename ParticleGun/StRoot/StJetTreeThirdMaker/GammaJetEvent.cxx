//___________________________
//  Event class to produce  picodst for GammaJet
//  May 30, 2015  :   Nihar R. Sahoo
//_____________________________


#include <TClonesArray.h>
#include "GammaJetEvent.h"
#include "GammaJetTrack.h"
#include "GammaJetTower.h"

ClassImp(GammaJetEvent);

//__________________________________________________________________            
GammaJetEvent::GammaJetEvent()
  :TObject()
  ,PrimaryTrackArray(0)
  ,TowerArray(0)
{

  
  //Initialization to -999

  runNumber       = -999;
  eventNumber     = -999;
  trigID         = -999;
  nGlobalTracks   = -999;
  nPrimaryTracks  = -999;
  refMult         = -999;
  vpdVz           = -999;
  xVertex         = -999;
  yVertex         = -999;
  zVertex         = -999;
  bbcZVertex      = -999;
  zdcCoincidenceRate      = -999;
  bbcCoincidenceRate      = -999;
  backgroundRate          = -999;
  bbcBlueBackgroundRate   = -999;
  bbcYellowBackgroundRate = -999;
  refMultPos              = -999;
  refMultNeg              = -999;
  bTOFTrayMultiplicity    = -999;
  nVerticies              = -999;
  MagF                    = -999;
  VrtxRank                =-999;
  FlagEvent_TrgTrkMisMtch =-999;

  //_______TSP Info
  //_______TSP information for event
   Etsp      =  -999.;
   ETwrdidT      =  -999.; // Tower Id
   ETwradc11      =  -999.;  // tower projection adc

   ETwreneT0      =  -999.; // Tower energy

   ETwreT      =  -999.; // Tower Eta true
   ETwrENET0      =  -999.;
   ETwrphT      =  -999.; // Tower phi true
   ETwrPTower      =  -999.; // projected track to this twr momentum
   ETwrpidTower      =  -999.; // projected track to this twr PID (dedx)
   ETwrmoduleT      =  -999.; // module no of this twr


   EClustEneT0      =  -999.;   // Cluster Energy
   EClustetav1      =  -999.; // Cluster Eta
   EClustphiv1      =  -999.; // Cluster phi

   EEstrpen01   = -999.;
   EEstrpen02   = -999.;
   EEstrpen03   = -999.;

   EEstrpen0      =  -999.;
   EEstrpen1      =  -999.;
   EEstrpen2      =  -999.;
   EEstrpen3      =  -999.;
   EEstrpen4      =  -999.;
   EEstrpen5      =  -999.;
   EEstrpen6      =  -999.;
   EEstrpen7      =  -999.;
   EEstrpen8      =  -999.;
   EEstrpen9      =  -999.;
   EEstrpen10      =  -999.;
   EEstrpen11      =  -999.;
   EEstrpen12      =  -999.;
   EEstrpen13      =  -999.;
   EEstrpen14      =  -999.;
   EEstrpen15      =  -999.; //eta strip energy from clustering


   ETwrdidE      =  -999.; // Eta strip Id

  //_______________ Phi Strip Info
   EPstripenp01      =  -999.;
   EPstripenp02      =  -999.;
   EPstripenp03      =  -999.;

   EPstripenp0      =  -999.;
   EPstripenp1      =  -999.;
   EPstripenp2      =  -999.;
   EPstripenp3      =  -999.;
   EPstripenp4      =  -999.;
   EPstripenp5      =  -999.;
   EPstripenp6      =  -999.;
   EPstripenp7      =  -999.;
   EPstripenp8      =  -999.;
   EPstripenp9      =  -999.;
   EPstripenp10      =  -999.;
   EPstripenp11      =  -999.;
   EPstripenp12      =  -999.;
   EPstripenp13      =  -999.;
   EPstripenp14      =  -999.;
   EPstripenp15      =  -999.; //eta strip energy from clustering


   EclustEnnq1      =  -999.;
   EclustEnnq20      =  -999.;
   EclustEnnq19      =  -999.;

   EclustEnpq1      =  -999.;
   EclustEnpq20      =  -999.;
   EclustEnpq19      =  -999.;
   EclustEnpq21      =  -999.;

   // new output branches [Derek, 12.23.2020/01.10.2021]
   EtspAlt               = -999.;
   EClustSecondTwrEnergy = -999.;
   EClustNumTwrIncluded  = 9999;
   EClustSecondTwrIndex  = -999;

   // new output branches [Derek, 01.12.2021]
   EEstrpModuleE         = -999;
   EPstripModuleP        = -999;
   EClustSecondTwrModule = -999;

  //__________


  PrimaryTrackArray = new TClonesArray("GammaJetTrack",100000);
  TowerArray = new TClonesArray("GammaJetTower",100000);
  
  //  TowerArray(0);

}


//__________________________________________________________________            
GammaJetEvent::GammaJetEvent(Float_t *eventarray){

  //Sets all the Event attributes from the passed in array                      
  runNumber       = (Long64_t)eventarray[0];
  eventNumber     = (Long64_t)eventarray[1];
  trigID         = (Int_t)eventarray[2];
  nGlobalTracks   = (Int_t)eventarray[3];
  nPrimaryTracks  = (Int_t)eventarray[4];
  refMult         = (Int_t)eventarray[5];
  vpdVz           = (Double_t)eventarray[6];
  xVertex         = (Double_t)eventarray[7];
  yVertex         = (Double_t)eventarray[8];
  zVertex         = (Double_t)eventarray[9];
  bbcZVertex      = (Double_t)eventarray[10];
  zdcCoincidenceRate      = (Double_t)eventarray[11];
  bbcCoincidenceRate      = (Double_t)eventarray[12];
  backgroundRate          = (Double_t)eventarray[13];
  bbcBlueBackgroundRate   = (Double_t)eventarray[14];
  bbcYellowBackgroundRate = (Double_t)eventarray[15];
  refMultPos              = (Double_t)eventarray[16];
  refMultNeg              = (Double_t)eventarray[17];
  bTOFTrayMultiplicity    = (Double_t)eventarray[18];
  nVerticies              = (Double_t)eventarray[19];
  MagF                    = (Double_t)eventarray[20];



  //_______TSP information for event
  Etsp      =  (Float_t)eventarray[21];
  ETwrdidT      =  eventarray[22]; // Tower Id
  ETwradc11      =  eventarray[23];  // tower projection adc

  ETwreneT0      =  (Float_t)eventarray[24]; // Tower energy

  ETwreT      =  (Float_t)eventarray[25]; // Tower Eta true
  ETwrENET0      =  (Float_t)eventarray[26];
  ETwrphT      =  (Float_t)eventarray[27]; // Tower phi true
  ETwrPTower      =  (Float_t)eventarray[28]; // projected track to this twr momentum
  ETwrpidTower      =  (Float_t)eventarray[29]; // projected track to this twr PID (dedx)
  ETwrmoduleT      =  (Int_t)eventarray[30]; // module no of this twr


  EClustEneT0      =  (Float_t)eventarray[31];   // Cluster Energy
  EClustetav1      =  (Float_t)eventarray[32]; // Cluster Eta
  EClustphiv1      =  (Float_t)eventarray[33]; // Cluster phi


  EEstrpen0      =  (Float_t)eventarray[34];
  EEstrpen1      =  (Float_t)eventarray[35];
  EEstrpen2      =  (Float_t)eventarray[36];
  EEstrpen3      =  (Float_t)eventarray[37];
  EEstrpen4      =  (Float_t)eventarray[38];
  EEstrpen5      =  (Float_t)eventarray[39];
  EEstrpen6      =  (Float_t)eventarray[40];
  EEstrpen7      =  (Float_t)eventarray[41];
  EEstrpen8      =  (Float_t)eventarray[42];
  EEstrpen9      =  (Float_t)eventarray[43];
  EEstrpen10      =  (Float_t)eventarray[44];
  EEstrpen11      =  (Float_t)eventarray[45];
  EEstrpen12      =  (Float_t)eventarray[46];
  EEstrpen13      =  (Float_t)eventarray[47];
  EEstrpen14      =  (Float_t)eventarray[48];
  EEstrpen15      =  (Float_t)eventarray[49]; //eta strip energy from clustering


  ETwrdidE      =  (Int_t)eventarray[50]; // Eta strip Id

  //_______________ Phi Strip Info
  EPstripenp01      =  (Float_t)eventarray[51];
  EPstripenp02      = (Float_t)eventarray[52];
  EPstripenp03      =  (Float_t)eventarray[53];

  EPstripenp0      =  (Float_t)eventarray[54];
  EPstripenp1      =  (Float_t)eventarray[55];
  EPstripenp2      =  (Float_t)eventarray[56];
  EPstripenp3      =  (Float_t)eventarray[57];
  EPstripenp4      =  (Float_t)eventarray[58];
  EPstripenp5      =  (Float_t)eventarray[59];
  EPstripenp6      =  (Float_t)eventarray[60];
  EPstripenp7      =  (Float_t)eventarray[61];
  EPstripenp8      =  (Float_t)eventarray[62];
  EPstripenp9      =  (Float_t)eventarray[63];
  EPstripenp10      =  (Float_t)eventarray[64];
  EPstripenp11      =  (Float_t)eventarray[65];
  EPstripenp12      =  (Float_t)eventarray[66];
  EPstripenp13      =  (Float_t)eventarray[67];
  EPstripenp14      =  (Float_t)eventarray[68];
  EPstripenp15      =  (Float_t)eventarray[69]; //eta strip energy from clustering


  EclustEnnq1      =  (Float_t)eventarray[70];
  EclustEnnq20      =  (Float_t)eventarray[71];
  EclustEnnq19      =  (Float_t)eventarray[72];

  EclustEnpq1      =  (Float_t)eventarray[73];
  EclustEnpq20      =  (Float_t)eventarray[74];
  EclustEnpq19      =  (Float_t)eventarray[75];
  EclustEnpq21      =  (Float_t)eventarray[76];


  EEstrpen01 =(Float_t)eventarray[77];
  EEstrpen02 =(Float_t)eventarray[78];
  EEstrpen03 =(Float_t)eventarray[79];

  VrtxRank    = (Double_t)eventarray[80];
  FlagEvent_TrgTrkMisMtch    = (Int_t)eventarray[81];

  // new output branches [Derek, 12.23.2020/01.10.2021]
  EtspAlt               = (Float_t) eventarray[82];
  EClustNumTwrIncluded  = (UInt_t)  eventarray[83];
  EClustSecondTwrIndex  = (Int_t)   eventarray[84];
  EClustSecondTwrEnergy = (Float_t) eventarray[85];

  // new output branches [Derek, 01.12.2021]
  EEstrpModuleE         = (Int_t) eventarray[86];
  EPstripModuleP        = (Int_t) eventarray[87];
  EClustSecondTwrModule = (Int_t) eventarray[88];

  //___________


  PrimaryTrackArray = new TClonesArray("GammaJetTrack",10000);
  TowerArray = new TClonesArray("GammaJetTower",10000);

}

//__________________________________________________________________            
GammaJetEvent::~GammaJetEvent(){

  //GammaJetEvent Destructor                                                            
  delete PrimaryTrackArray;
  delete TowerArray;

}

//__________________________________________________________________            
void GammaJetEvent::PrintEvent(Bool_t printtrack){

  //Print Event Information                                                     
  cout <<"--- Event: " <<eventNumber <<" --- " <<"RunNumber: " << runNumber <<" --- " <<"\n"
       <<"  " <<"PrimaryTracks: " <<nPrimaryTracks  <<"\n"
       <<"  " <<"Vx: " <<xVertex <<" Vy: " <<yVertex <<" Vz: " <<zVertex <<"\n"
       <<"  " <<"Entries in Track Array: " <<PrimaryTrackArray->GetEntries() <<"\n";

  //Print Primary Track Information if Requested                                
  if (printtrack){

    for (Int_t i=0; i<PrimaryTrackArray->GetEntries(); i++){
      ((GammaJetTrack*)PrimaryTrackArray->At(i))->PrintTrack();
    }

  }//End Print Tracks                                                           

  cout <<"End Of Event--\n";
}

//__________________________________________________________________            
void GammaJetEvent::ResetEvent(){

  //Clear the Track TClonesArray                                                
  PrimaryTrackArray->Clear("C");
  TowerArray->Clear("C");
}


//__________________________________________________________________            
//void GammaJetEvent::AddTrack(Float_t *trackarray, Int_t idx){
void GammaJetEvent::AddTrack(GammaJetTrack *trackarray, Int_t idx){

  TClonesArray &traks = *PrimaryTrackArray;
  
  GammaJetTrack *tempTrack = (GammaJetTrack*)traks.ConstructedAt(idx);
  //  tempTrack->SetTrackAttributes(trackarray);
  *tempTrack = *trackarray;

}


//___________________________________________________________
//void GammaJetEvent::AddTower(Float_t *towerarray, Int_t id){
void GammaJetEvent::AddTower( GammaJetTower *t, Int_t id){
  
  TClonesArray &towers = *TowerArray;

  GammaJetTower *tower = (GammaJetTower*)towers.ConstructedAt(id);
  *tower = *t;

}

//__________________________________________________________________            
void GammaJetEvent::SetEventAttributes(Float_t *eventarray){

  //Sets all the Event attributes from the passed in array                      
  runNumber       = (Long64_t)eventarray[0];
  eventNumber     = (Long64_t)eventarray[1];
  trigID         = (Int_t)eventarray[2];
  nGlobalTracks   = (Int_t)eventarray[3];
  nPrimaryTracks  = (Int_t)eventarray[4];
  refMult         = (Int_t)eventarray[5];
  vpdVz           = (Double_t)eventarray[6];
  xVertex         = (Double_t)eventarray[7];
  yVertex         = (Double_t)eventarray[8];
  zVertex         = (Double_t)eventarray[9];
  bbcZVertex      = (Double_t)eventarray[10];
  zdcCoincidenceRate      = (Double_t)eventarray[11];
  bbcCoincidenceRate      = (Double_t)eventarray[12];
  backgroundRate          = (Double_t)eventarray[13];
  bbcBlueBackgroundRate   = (Double_t)eventarray[14];
  bbcYellowBackgroundRate = (Double_t)eventarray[15];
  refMultPos              = (Double_t)eventarray[16];
  refMultNeg              = (Double_t)eventarray[17];
  bTOFTrayMultiplicity    = (Double_t)eventarray[18];
  nVerticies              = (Double_t)eventarray[19];
  MagF                    = (Double_t)eventarray[20];

  //_______TSP information for event
  Etsp      =  (Float_t)eventarray[21];
  ETwrdidT      =  eventarray[22]; // Tower Id
  ETwradc11      =  eventarray[23];  // tower projection adc

  ETwreneT0      =  (Float_t)eventarray[24]; // Tower energy

  ETwreT      =  (Float_t)eventarray[25]; // Tower Eta true
  ETwrENET0      =  (Float_t)eventarray[26];
  ETwrphT      =  (Float_t)eventarray[27]; // Tower phi true
  ETwrPTower      =  (Float_t)eventarray[28]; // projected track to this twr momentum
  ETwrpidTower      =  (Float_t)eventarray[29]; // projected track to this twr PID (dedx)
  ETwrmoduleT      =  (Int_t)eventarray[30]; // module no of this twr


  EClustEneT0      =  (Float_t)eventarray[31];   // Cluster Energy
  EClustetav1      =  (Float_t)eventarray[32]; // Cluster Eta
  EClustphiv1      =  (Float_t)eventarray[33]; // Cluster phi


  EEstrpen0      =  (Float_t)eventarray[34];
  EEstrpen1      =  (Float_t)eventarray[35];
  EEstrpen2      =  (Float_t)eventarray[36];
  EEstrpen3      =  (Float_t)eventarray[37];
  EEstrpen4      =  (Float_t)eventarray[38];
  EEstrpen5      =  (Float_t)eventarray[39];
  EEstrpen6      =  (Float_t)eventarray[40];
  EEstrpen7      =  (Float_t)eventarray[41];
  EEstrpen8      =  (Float_t)eventarray[42];
  EEstrpen9      =  (Float_t)eventarray[43];
  EEstrpen10      =  (Float_t)eventarray[44];
  EEstrpen11      =  (Float_t)eventarray[45];
  EEstrpen12      =  (Float_t)eventarray[46];
  EEstrpen13      =  (Float_t)eventarray[47];
  EEstrpen14      =  (Float_t)eventarray[48];
  EEstrpen15      =  (Float_t)eventarray[49]; //eta strip energy from clustering


  ETwrdidE      =  (Int_t)eventarray[50]; // Eta strip Id

  //_______________ Phi Strip Info
  EPstripenp01      =  (Float_t)eventarray[51];
  EPstripenp02      = (Float_t)eventarray[52];
  EPstripenp03      =  (Float_t)eventarray[53];

  EPstripenp0      =  (Float_t)eventarray[54];
  EPstripenp1      =  (Float_t)eventarray[55];
  EPstripenp2      =  (Float_t)eventarray[56];
  EPstripenp3      =  (Float_t)eventarray[57];
  EPstripenp4      =  (Float_t)eventarray[58];
  EPstripenp5      =  (Float_t)eventarray[59];
  EPstripenp6      =  (Float_t)eventarray[60];
  EPstripenp7      =  (Float_t)eventarray[61];
  EPstripenp8      =  (Float_t)eventarray[62];
  EPstripenp9      =  (Float_t)eventarray[63];
  EPstripenp10      =  (Float_t)eventarray[64];
  EPstripenp11      =  (Float_t)eventarray[65];
  EPstripenp12      =  (Float_t)eventarray[66];
  EPstripenp13      =  (Float_t)eventarray[67];
  EPstripenp14      =  (Float_t)eventarray[68];
  EPstripenp15      =  (Float_t)eventarray[69]; //eta strip energy from clustering


  EclustEnnq1      =  (Float_t)eventarray[70];
  EclustEnnq20      =  (Float_t)eventarray[71];
  EclustEnnq19      =  (Float_t)eventarray[72];

  EclustEnpq1      =  (Float_t)eventarray[73];
  EclustEnpq20      =  (Float_t)eventarray[74];
  EclustEnpq19      =  (Float_t)eventarray[75];
  EclustEnpq21      =  (Float_t)eventarray[76];


  EEstrpen01 =(Float_t)eventarray[77];
  EEstrpen02 =(Float_t)eventarray[78];
  EEstrpen03 =(Float_t)eventarray[79];

  VrtxRank    = (Double_t)eventarray[80];
  FlagEvent_TrgTrkMisMtch = (Int_t)eventarray[81];

  // new output branches [Derek, 12.23.2020/01.10.2021]
  EtspAlt               = (Float_t) eventarray[82];
  EClustNumTwrIncluded  = (UInt_t)  eventarray[83];
  EClustSecondTwrIndex  = (Int_t)   eventarray[84];
  EClustSecondTwrEnergy = (Float_t) eventarray[85];

  // new output branches [Derek, 01.12.2021]
  EEstrpModuleE         = (Int_t) eventarray[86];
  EPstripModuleP        = (Int_t) eventarray[87];
  EClustSecondTwrModule = (Int_t) eventarray[88];

  //___________




}


//__________________________________________________________________            

Int_t GammaJetEvent::GetTrackArrayEntries(){
  return PrimaryTrackArray->GetEntries();
}
//__________________________________________________________________            

Int_t GammaJetEvent::GetTowerArrayEntries(){
  return TowerArray->GetEntries();
}
