// Nihar R. Sahoo
// This is used to do some QA on event, BEMC tower and 
// after then produce same "Gfemtodst" 
// 
//  One can use this class to read Gfemtodst and do analysis ..... Jul11 2016
//____________________________________

#define StGJetTreeReader_cxx
#include "StGJetTreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

#include "GammaJetTrack.h"
#include "GammaJetEvent.h"
#include "GammaJetTower.h"
#include "GammaJetTowerUtil.h"
#include <exception>
#include <vector>


using namespace std;


ClassImp( StGJetTreeReader);

void StGJetTreeReader::Make()
{

  Init_Picotree("pp_200_test.root");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   //______________________________
   for (int g = 0; g < 81; g++){ EvValues[g] = -2000; }
   
   
   Int_t evnt_counter=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;



      evnt_counter++;
      if( evnt_counter% 1000 ==0)cout<<"processing event: "<<evnt_counter<<endl;
      
      //      cout<<" RefMult= "<<refMult<<"  TSP= "<< Etsp<<"  Vz= "<<  zVertex
      //	  <<"runNumber= "<<runNumber<<" nPrimaryTracks=  "<<nPrimaryTracks
      //  <<endl;

      //__________Refilling of Event Info
      EvValues[0]   = runNumber;
      EvValues[1]   = eventNumber;
      EvValues[2]   = trigID;  //Trigger Id check getTrgId()
      EvValues[3]   = nGlobalTracks;
      EvValues[4]   = nPrimaryTracks;
      EvValues[5]   = refMult;
      EvValues[6]   = vpdVz;
      EvValues[7]   = xVertex;
      EvValues[8]   = yVertex;
      EvValues[9]   = zVertex;
      EvValues[10]   = bbcZVertex;
      EvValues[11]   = zdcCoincidenceRate;
      EvValues[12]   = bbcCoincidenceRate;
      EvValues[13]   = backgroundRate;
      EvValues[14]   = bbcBlueBackgroundRate;
      EvValues[15]   = bbcYellowBackgroundRate;
      EvValues[16]   = refMultPos;
      EvValues[17]   = refMultNeg;
      EvValues[18]   = bTOFTrayMultiplicity;
      EvValues[19]   = nVerticies;
      EvValues[20]   = MagF;
      EvValues[80]   = VrtxRank;
      

      EvValues[21]= Etsp;    //e_TSP;
      EvValues[22]= ETwrdidT;  //e_didT;
      EvValues[23]= ETwradc11;//e_adc11;
      EvValues[24]= ETwreneT0; //e_eneT0;
      EvValues[25]= ETwreT;    //e_eT;
      EvValues[26]= ETwrENET0;//e_ENET0;
      EvValues[27]= ETwrphT; //e_phT;
      EvValues[28]= ETwrPTower; //e_PTower;
      EvValues[29]= ETwrpidTower; //e_pidTower;
      EvValues[30]= ETwrmoduleT; //e_moduleT;
      EvValues[31]= EClustEneT0; //e_EneT0;
      EvValues[32]= EClustetav1;  //e_etav1;
      EvValues[33]= EClustphiv1; //e_phiv1;
      EvValues[34]= EEstrpen0;//e_en0;
      EvValues[35]= EEstrpen1; //e_en1;
      EvValues[36]= EEstrpen2; //e_en2;
      EvValues[37]= EEstrpen3; //e_en3;
      EvValues[38]= EEstrpen4; //e_en4;
      EvValues[39]= EEstrpen5; //e_en5;
      EvValues[40]= EEstrpen6; //e_en6;
      EvValues[41]= EEstrpen7; //e_en7;
      EvValues[42]= EEstrpen8; //e_en8;
      EvValues[43]= EEstrpen9; //e_en9;
      EvValues[44]= EEstrpen10; //e_en10;
      EvValues[45]= EEstrpen11; //e_en11;
      EvValues[46]= EEstrpen12; //e_en12;
      EvValues[47]= EEstrpen13; //e_en13;
      EvValues[48]= EEstrpen14; //e_en14;
      EvValues[49]= EEstrpen15; //e_en15;
      EvValues[50]= ETwrdidE; //e_didE;


      EvValues[51]= EPstripenp01; //e_enp01;
      EvValues[52]= EPstripenp02;//e_enp02;
      EvValues[53]= EPstripenp03; //e_enp03;
      EvValues[54]= EPstripenp0; //e_enp0;
      EvValues[55]= EPstripenp1; //e_enp1;
      EvValues[56]= EPstripenp2; //e_enp2;
      EvValues[57]= EPstripenp3; //e_enp3;
      EvValues[58]= EPstripenp4; //e_enp4;
      EvValues[59]= EPstripenp5; //e_enp5;


      EvValues[60]= EPstripenp6; //e_enp6;
      EvValues[61]= EPstripenp7; //e_enp7;
      EvValues[62]= EPstripenp8; //e_enp8;
      EvValues[63]= EPstripenp9; //e_enp9;
      EvValues[64]= EPstripenp10; //e_enp10;
      EvValues[65]= EPstripenp11; //e_enp11;
      EvValues[66]= EPstripenp12;//e_enp12;
      EvValues[67]= EPstripenp13; //e_enp13;
      EvValues[68]= EPstripenp14; //e_enp14;
      EvValues[69]= EPstripenp15;//e_enp15;


      EvValues[70]= EclustEnnq1;//e_Ennq1;
      EvValues[71]= EclustEnnq20; //e_Ennq20;
      EvValues[72]= EclustEnnq19;//e_Ennq19;


      EvValues[73]= EclustEnpq1;//e_Enpq1;
      EvValues[74]= EclustEnpq20;//e_Enpq20;
      EvValues[75]= EclustEnpq19; //e_Enpq19;
      EvValues[76]= EclustEnpq21; //e_Enpq21;


      EvValues[77]= EEstrpen01; //e_en01;
      EvValues[78]= EEstrpen02;//e_en02;
      EvValues[79]= EEstrpen03; //e_en03;

      event->SetEventAttributes(EvValues);      

      //_________________Primary Trk Info Refill_____________
      for(int i=0; i<28; i++){TrValues[i]=-9999;}



      GammaJetTrack TrackInfo;

      Int_t track_counter = 0;
      for(Int_t ii=0; ii < PrimaryTrackArray_; ii++)
	{
	  TrackInfo.SetnHitsFit(PrimaryTrackArray_nHitsFit[ii]);
	  TrackInfo.SetnHitsPoss(PrimaryTrackArray_nHitsPoss[ii]);
	  TrackInfo.SetTrackFlag(PrimaryTrackArray_trackFlag[ii]);
	  TrackInfo.SetpZ(PrimaryTrackArray_pZ[ii]);
	  TrackInfo.SetpY(PrimaryTrackArray_pY[ii]);
	  TrackInfo.SetpX(PrimaryTrackArray_pX[ii]);
	  TrackInfo.SetpT(PrimaryTrackArray_pT[ii]);
	  TrackInfo.SetdEdx(PrimaryTrackArray_dEdx[ii]);
	  TrackInfo.SetCharge(PrimaryTrackArray_charge[ii]);
	  TrackInfo.SetTOFBeta(PrimaryTrackArray_tofBeta[ii]);
	  TrackInfo.SetEta(PrimaryTrackArray_eta[ii]);
	  TrackInfo.SetPhi(PrimaryTrackArray_phi[ii]);
	  TrackInfo.SetnSigElectron(PrimaryTrackArray_nSigElectron[ii]);
	  TrackInfo.SetnSigPion(PrimaryTrackArray_nSigPion[ii]);
	  TrackInfo.SetnSigKaon(PrimaryTrackArray_nSigKaon[ii]);
	  TrackInfo.SetnSigProton(PrimaryTrackArray_nSigProton[ii]);
	  TrackInfo.SetDCAg(PrimaryTrackArray_dcag[ii]);
	  TrackInfo.SetnHits(PrimaryTrackArray_nHits[ii]);
	  TrackInfo.SetdEdxHits(PrimaryTrackArray_dEdxHits[ii]);
	  TrackInfo.SetFirstZPoint(PrimaryTrackArray_firstZPoint[ii]);
	  TrackInfo.SetLastZPoint(PrimaryTrackArray_lastZPoint[ii]);
	  TrackInfo.SetTOFSigElectron(PrimaryTrackArray_tofSigElectron[ii]);
	  TrackInfo.SetTOFSigPion(PrimaryTrackArray_tofSigPion[ii]);
	  TrackInfo.SetTOFSigKaon(PrimaryTrackArray_tofSigKaon[ii]);
	  TrackInfo.SetTOFSigProton(PrimaryTrackArray_tofSigProton[ii]);
	  TrackInfo.SetPathLength(PrimaryTrackArray_pathLength[ii]);
	  TrackInfo.SettimeOfflight(PrimaryTrackArray_timeOfflight[ii]);
	  //	  TrackInfo.SettrkIndex(trkIndex_counter);



	  event->AddTrack(&TrackInfo,track_counter);
	  //	  event->AddTrack(TrValues,track_counter);
	  track_counter++;	  

	}// Trk Loop

      
      //______________Tower Info Refill
      GammaJetTower TowerInfo;
      
      Int_t twr_counter = 0;
      vector<Int_t> vec_matchTrkIndex; //to store ptrks index
      vec_matchTrkIndex.clear();

      vector<Float_t> vec_matchTrkP; //to store ptrks p
      vec_matchTrkP.clear();


      for(int it=0; it < TowerArray_; it++)
	{
	  
	  TowerInfo.SetTwrId(TowerArray_TwrId[it]);
	  TowerInfo.SetTwrEng(TowerArray_TwrEng[it]);
	  TowerInfo.SetTwrPhi(TowerArray_TwrPhi[it]);
	  TowerInfo.SetTwrEta(TowerArray_TwrEta[it]);
	  TowerInfo.SetTwrADC(TowerArray_TwrADC[it]);
	  


	  TowerInfo.SetTwrMatchIdnex(TowerArray_TwrMatchIdnex[it]);
	  TowerInfo.SetNoOfmatchedTrk(TowerArray_NoOfmatchedTrk[it]);
	  
	  TowerInfo.SetTwrPx( TowerArray_TwrPx[it]);
	  TowerInfo.SetTwrPy( TowerArray_TwrPy[it]);
	  TowerInfo.SetTwrPz( TowerArray_TwrPz[it]);
	  
	  for(Int_t im=0; im < TowerArray_NoOfmatchedTrk[it]; im++)
	    {
	      TowerInfo.SetMatchedTracksArray(TowerArray_fMatchedTracksArray_[it][im],im);
	      TowerInfo.SetMatchedTracksArray_P(TowerArray_fMatchedTracksArray_P[it][im],im);
	    }




	  event->AddTower(&TowerInfo,twr_counter);
	  twr_counter++;
	  
	}
      
      

      
      


      //_________________
      outTree->Fill();

   } //event loop

   Finish();

} //Make
