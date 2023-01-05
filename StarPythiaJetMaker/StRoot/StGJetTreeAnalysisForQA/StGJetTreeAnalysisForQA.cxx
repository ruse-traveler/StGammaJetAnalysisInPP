
// Nihar R. Sahoo
// This is used to do some QA on event, BEMC tower and 
// after then produce same "Gfemtodst" 
// 
//  One can use this class to read Gfemtodst and do analysis ..... Jul11 2016
//____________________________________

#define StGJetTreeAnalysisForQA_cxx
#include "StGJetTreeAnalysisForQA.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

//#include "GammaJetTrack.h"
//#include "GammaJetEvent.h"
//#include "GammaJetTower.h"
#include "GammaJetTowerUtil.h"
#include "GammaJetQAUtil.h"
#include <exception>
#include <vector>


#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include  "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TLeaf.h>

//#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"


using namespace std;
using namespace fastjet;

//_______Const Defined
static const double Value_Pi = TMath::Pi();
static const double Mass_pi = 0.140;  // 140 MeV
static const double BEMC_radius =  225.405;

static const double PT_TRACK_MAX = 15; 
static const double PT_TRACK_MIN = 0.2; 

static const double ETA_TRACK_MAX = 1.0; 
static const double ETA_TRACK_MIN = -1.0; 

//______Jet Input
static const Double_t Jet_R = 0.3;
static const Int_t Remove_N_hardest = 3;

static const Int_t Jet_type = 0;                     //           "FullJet" = 1  and "ChargedJet" = 0

//TNtuple *ntuple_Trig;

TDatime datime;


//________Histo
TH1D *hTwrEt;
TH1D *hpTTrack;
TH1D *hTwrTotalE;
TProfile *pTwrIdE;
TH2D *h_grefmult_tofmult;
TH2D *h_grefmult_primTrkana;
TH2D *h_Tofmatch_primTrkana;
TH2D *hTofMatch_Refmult;
TH1D *h_TrigTwr_Id;// = new TH1D("h_TrigTwr_Id","TrigTwr_Id",5000,0,5000);
TH1D *h_TrigBSMD_EtaStrpId;// = new TH1D("h_TrigBSMD_EtaStrpId","TrigBSMD_EtaStrpId",20000,0,20000);
TH1D *h_TrigBSMD_PhiStrpId;// = new TH1D("h_TrigBSMD_PhiStrpId","TrigBSMD_PhiStrpId",20000,0,20000);

TH1D *hEbp_elec;
TH1D *hEbp_hadrn;

const   Int_t _N_HISTOGRAMS_ = 3;
TH1D*   HP[_N_HISTOGRAMS_];
TH1D*   HG[_N_HISTOGRAMS_];
TH1D*   HPtrig;
TH1D*   HGtrig;
TH1D*   HP_;  TH1D*   HG_;

TH1D* hPrimTrk_1;
TH1D* hPrimTrk_2;
TH1D* hPrimTrk_3;
TH1D* hPrimTrk_4;
TH1D* hPrimTrk_5;
TH1D* hPrimTrk_6;
TH1D* hPrimTrk_7;
TH1D* hPrimTrk_8;
TH1D* hPrimTrk_9;
TH1D* hPrimTrk_10;
TH1D* hPrimTrk_11;
TH1D* hPrimTrk_12;TH1D* hPrimTrk_0;
TH1D* hPrimTrk_13;
TH1D* hPrimTrk_14;
TH1D* hPrimTrk_15;
TH1D* hPrimTrk_16;
TH1D* hPrimTrk_17;
TH1D* hPrimTrk_18;


TH1D* hTSP;
TH1D* hTSP_pi0;
TH1D* hTSP_g;

TH1D*   hpi0_Et;  TH1D*   hg_Et;

TNtuple *trigNtuple;
TNtuple *trigTrckNtuple;

const double  zT_asso_max  =1.0;        //   Set Upper cut for pT assoc.
const double  zT_asso_min  =0.;  


//__________________

ClassImp( StGJetTreeAnalysisForQA);

void StGJetTreeAnalysisForQA::Make()
{

  TString OutFileName;
  //  OutFileName = "JetRootFiles/pp200_";
  OutFileName = "/star/data01/pwg/nihar/JetReco_pico/Run14_AuAu_WB_prod/output_QA_May24/";
  OutFileName += "AuAu_Run14_upto200_PrimarTrckForMIxing";

  //  OutFileName += "PtDist";
  //  OutFileName += Jet_type;
  //  OutFileName += "_R";
  //  OutFileName += Jet_R*10;
  //  OutFileName += "_Zt";
  //  OutFileName += zT_asso_max*10;
  //  OutFileName += "to";
  //  OutFileName += zT_asso_min*10;
  OutFileName += "_min" ;
  OutFileName += datime.GetMinute();
  OutFileName += "d" ;
  OutFileName +=  datime.GetDay();
  OutFileName += "m" ;
  OutFileName +=  datime.GetMonth();
  OutFileName += "y" ;
  OutFileName +=  datime.GetYear();
  OutFileName +=  ".root";

  cout<<" Output File: "<<OutFileName.Data()<<endl;

  //  Init_Picotree("JetRootFiles/pp200_R03_Charged_2017Jan10.root");
  Init_Picotree(OutFileName.Data());
  //  Init_Picotree("/star/data01/pwg/nihar/JetReco_pico/Run14_AuAu_WB_prod/output_QA_May24/AuAuRun14_jetreco.root");
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Int_t evnt_counter=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
            
      //__________________________________
      Double_t Vr = sqrt(pow(xVertex,2)+pow(yVertex,2));

      Bool_t event_check = Event_Cut(Vr,zVertex);
      if(!(event_check))continue;
      
      if(!( EEstrpen4 >=0.5 && EPstripenp4 >=0.5  && (ETwreT > -0.9 && ETwreT < 0.9 )  ) )continue;
      if(!(ETwradc11 <=6004 && ETwrPTower < 3.0)) continue;
      if(!(EClustetav1 >= -1.0 && EClustetav1 <= 1.0))continue;

      //      if(!( gRefmult > 315))continue;
      //      if(!( gRefmult > 315))continue;

      //      if(!( Etsp < 0.08 || (Etsp > 0.2 && Etsp < 0.6)))continue;

      Double_t Et_Twr = EClustEneT0*sin(2*atan(exp(-1*EClustetav1)));      

      ///___________________
      hTwrEt ->Fill(Et_Twr);


      //if(!(Et_Twr > 8))continue;
      if(!(Et_Twr > 9 && Et_Twr < 30))continue;

      hRefmult->Fill(refMult);      
      //      GammaJetQAUtil *QAInfo = new GammaJetQAUtil("AuAu200_run11");
      h_grefmult_tofmult->Fill(gRefmult,Tofmult);
      //      h_TrigTwr_Id ->Fill(ETwrdidT);
      //      h_TrigBSMD_EtaStrpId ->Fill(ETwrdidE);
      //      h_TrigBSMD_PhiStrpId ->Fill(ETwrdidP);



      hVpdVz_Vz ->Fill(zVertex,vpdVz);
      hVpdVzVzDiff->Fill(vpdVz-zVertex);
      hRefmultPos ->Fill(refMultPos);
      hRefmultNeg ->Fill(refMultNeg);
      hRefmultPos_tofmult ->Fill(refMultPos,Tofmult);
      hRefmultNeg_tofmult ->Fill(refMultNeg,Tofmult);


      
      evnt_counter++;
      if( evnt_counter% 1000 ==0)cout<<"processing event: "<<evnt_counter<<endl;
      if( evnt_counter == 1024000)break;
      //      cout<<"Vz= "<<zVertex<<"  clstEta=    "<<EClustetav1<<"  ETwreneT0=  "<<ETwreneT0<<" EClustEneT0= "<<EClustEneT0<<endl;
                 

      //_________________Primary Trk Info Reading_____________

      Int_t track_counter = 0;
      Int_t track_nonana_counter = 0;
      Int_t track_counter1 = 0;
      Int_t tofmatched =0;

      hprimTrk->Fill(PrimaryTrackArray_);
      
      if( Jet_type == 1 || Jet_type == 0)
	{
	  for(Int_t ii=0; ii < PrimaryTrackArray_; ii++)
	    {
	      double E = 0.;

	      if(PrimaryTrackArray_tofBeta[ii] >0 &&  PrimaryTrackArray_charge[ii] != 0 && fabs(PrimaryTrackArray_eta[ii]) < 0.5 && PrimaryTrackArray_nHitsFit[ii] > 10 && PrimaryTrackArray_dcag[ii] < 3.0) tofmatched++;


	      if(PrimaryTrackArray_nHitsFit[ii] >14){  //
		if(  ((1.0*PrimaryTrackArray_nHitsFit[ii])/(1.0*PrimaryTrackArray_nHitsPoss[ii]) ) > 0.52 ){ ///
		  if(PrimaryTrackArray_dcag[ii]< 1.0){////
    
		    if(!( PrimaryTrackArray_pT[ii] > PT_TRACK_MIN && PrimaryTrackArray_pT[ii] < PT_TRACK_MAX))continue;
		    if(!( PrimaryTrackArray_eta[ii] > ETA_TRACK_MIN && PrimaryTrackArray_eta[ii] < ETA_TRACK_MAX))continue;
		    //if(!( PrimaryTrackArray_eta[ii] > -0.5 && PrimaryTrackArray_eta[ii] < 0.5))continue;


		    double p = sqrt( PrimaryTrackArray_pX[ii]*PrimaryTrackArray_pX[ii] + PrimaryTrackArray_pY[ii]*PrimaryTrackArray_pY[ii] + PrimaryTrackArray_pZ[ii]*PrimaryTrackArray_pZ[ii]);

		    hpTTrack->Fill(PrimaryTrackArray_pT[ii]);
		
		    E = sqrt( (PrimaryTrackArray_pX[ii]*PrimaryTrackArray_pX[ii])+ (PrimaryTrackArray_pY[ii]*PrimaryTrackArray_pY[ii]) + (PrimaryTrackArray_pZ[ii]*PrimaryTrackArray_pZ[ii]) + (Mass_pi*Mass_pi));
		    

		    //___________________________
		      track_counter++;

		  }////
		}///
	      }//

      



	      //	      cout<<"............................ inside track loop....."<<endl;
	    }// Trk Loop
	}// type cond	  
	 //      cout<< "track_counter .... = "<<track_counter<<"  ptrk#  "<<PrimaryTrackArray_<<endl;
      //      cout<<track_counter1<<endl;
      if(!(track_counter > 410))continue;

      h_Tofmatch_primTrkana->Fill(tofmatched,track_counter);
      h_grefmult_primTrkana->Fill(gRefmult,track_counter);

      hTofMatch_Refmult->Fill(tofmatched,refMult);


      hprimTrk_ana->Fill(track_counter);


      if( track_counter >= 410 && track_counter < 460 ){ hPrimTrk_0->Fill(track_counter);}
      if( track_counter >= 460 && track_counter < 510 ){ hPrimTrk_1->Fill(track_counter);}
      if( track_counter >= 510 && track_counter < 560 ){ hPrimTrk_2->Fill(track_counter);}
      if( track_counter >= 560 && track_counter < 610 ){ hPrimTrk_3->Fill(track_counter);}
      if( track_counter >= 610 && track_counter < 660 ){ hPrimTrk_4->Fill(track_counter);}
      if( track_counter >= 660 && track_counter < 710 ){ hPrimTrk_5->Fill(track_counter);}
      if( track_counter >= 710 && track_counter < 760 ){ hPrimTrk_6->Fill(track_counter);}
      if( track_counter >= 760){ hPrimTrk_7->Fill(track_counter);}



#if 0      
      if( !(track_counter != 0))continue;


      //      cout<<tofmatched<<"   "<<refMult<<endl;

	  //______________Tower Info Reading
	  //      Int_t twr_counters=0;
      //	  GammaJetTowerUtil TwrUtil;
      //TVector3 VertexVec(xVertex,yVertex,zVertex);

	  if( Jet_type == 0)
	    {
	      for(int it=0; it < TowerArray_; it++)
		{
		  
		  if(!(TowerArray_TwrEng[it] > 0.2 ))continue;
		  //		  if(!( TowerArray_TwrEta[it] > -1.0 && TowerArray_TwrEta[it] < 1.0))continue;	  

		  if(!(TowerArray_TwrId[it]==ETwrdidT))continue;

		  //		  cout<<"ETwrdidT= "<<ETwrdidT<<"  TowerArray_TwrId= "<<TowerArray_TwrId[it]<<"  Eng= "<<TowerArray_TwrEng[it]<<endl;
		  
		  //____________________________________E/P for each tower and projected track

		  //		  cout<<"Tower: eta= "<<TowerArray_TwrEta[it]<<" phi=  "<<TowerArray_TwrPhi[it]<<endl;
		  for(int j=0; j < 10; j++){

		    if(!(TowerArray_fMatchedTracksArray_[it][j] >=0))continue;//cout<<TowerArray_fMatchedTracksArray_[it][j]<<endl;

		    //		    for(Int_t ii=0; ii < PrimaryTrackArray_; ii++)
			
		      Int_t matck_index = TowerArray_fMatchedTracksArray_[it][j];

		      /*
		      if(!(PrimaryTrackArray_nHitsFit[matck_index] >14))continue;
		      if(!((1.0*PrimaryTrackArray_nHitsFit[matck_index])/(1.0*PrimaryTrackArray_nHitsPoss[matck_index]) ) > 0.52 )continue;
		      if(!(PrimaryTrackArray_dcag[matck_index]< 3.0))continue;
		      */
			
		      double pTrk = sqrt( (PrimaryTrackArray_pX[matck_index]*PrimaryTrackArray_pX[matck_index])+ (PrimaryTrackArray_pY[matck_index]*PrimaryTrackArray_pY[matck_index]) + (PrimaryTrackArray_pZ[matck_index]*PrimaryTrackArray_pZ[matck_index]));

		      if( fabs(PrimaryTrackArray_nSigElectron[matck_index])< 3.0 && fabs(PrimaryTrackArray_nSigProton[matck_index]) > 3.0 && fabs(PrimaryTrackArray_nSigKaon[matck_index]) > 3.0 && fabs(PrimaryTrackArray_nSigPion[matck_index]) > 3.0)
			{
			  //			  cout<<" electron "<<"    E/p= "<<EClustEneT0/pTrk<<endl;
			  hEbp_elec->Fill(EClustEneT0/pTrk);
			}
		      else if( fabs(PrimaryTrackArray_nSigElectron[matck_index])> 3.0 && fabs(PrimaryTrackArray_nSigProton[matck_index]) < 3.0 && fabs(PrimaryTrackArray_nSigKaon[matck_index]) < 3.0 && fabs(PrimaryTrackArray_nSigPion[matck_index]) < 3.0)
			    {
			      //  cout<<" hadron "<<"    E/p= "<<EClustEneT0/pTrk<<endl;
			      hEbp_hadrn->Fill(EClustEneT0/pTrk);

			    }
		      else { 
			//cout<<" Not found hadron/electron "<<endl;
		    }

		    

		}


		  //___________________________________

		}
	    }
#endif
 
	  //	  outTree->Fill();
   } //event loop


   //   File_output = outTree->GetCurrentFile();

   File_output->cd();


   hPrimTrk_1->Write();
   hPrimTrk_2->Write();
   hPrimTrk_3->Write();
   hPrimTrk_4->Write();
   hPrimTrk_5->Write();
   hPrimTrk_6->Write();
   hPrimTrk_7->Write();
   hPrimTrk_8->Write();
   hPrimTrk_9->Write();
   hPrimTrk_10->Write();
   hPrimTrk_11->Write();
   hPrimTrk_12->Write();
   hPrimTrk_0->Write();
   hPrimTrk_13->Write();
   hPrimTrk_14->Write();
   hPrimTrk_15->Write();
   hPrimTrk_16->Write();
   hPrimTrk_17->Write();
   hPrimTrk_18->Write();
   hprimTrk_ana->Write();

#if 0
   h_Tofmatch_primTrkana->Write();
   double scale = 1./(1.0*evnt_counter*0.1*2);
   hpTTrack->Scale(scale);
   hpTTrack->Write();


   hRefmult->Write();

   hRefmultPos ->Write();
   hRefmultNeg ->Write();
   hRefmultPos_tofmult->Write(); 
   hRefmultNeg_tofmult ->Write();
   hVpdVz_Vz ->Write();
   hVpdVzVzDiff->Write();
   h_grefmult_primTrkana->Write();
   HPtrig ->Write();
   HGtrig ->Write();

   hTSP->Write();
   hTSP_pi0->Write();
   hTSP_g->Write();
   
   hpi0_Et->Write();
   hg_Et->Write();
   
   hTwrEt->Write();
   
   hTwrTotalE->Write();
   pTwrIdE->Write();
   h_grefmult_tofmult->Write();
   hTofMatch_Refmult->Write();
   h_TrigTwr_Id->Write();
   h_TrigBSMD_EtaStrpId->Write();
   h_TrigBSMD_PhiStrpId->Write();

   hprimTrk->Write();
   hEbp_elec->Write();
   hEbp_hadrn->Write();
#endif
   //   File_output->Write();
   File_output->Close();

   //   Finish();
   
} //Make



//______________________________________________
Bool_t   StGJetTreeAnalysisForQA::Event_Cut(Double_t Vr, Double_t Vz)
{
  
  if( Vr < 2.0  
      && fabs(Vz) < 70
      && fabs(Vz) >= 0.0001
      
      ){
  return kTRUE;
  } else return kFALSE;
}

//______________________________________________________________
void StGJetTreeAnalysisForQA::Init_Picotree(const char* outFile)
{

  //  ntuple_Trig = new TNtuple("ntuple_Trig","ntuple_Trig","Vz:Ep:refmult:EvntNo:TrgEng:TrgPhi:TrgEta");
  trigNtuple= new TNtuple("trigNtuple","trigNtuple","Vz:tsp:refmult:TrgEng:TrgPhi:TrgEta");
  trigTrckNtuple = new TNtuple("trigTrckNtuple","trigTrckNtuple","Vz:tsp:refmult:TrgEng:TrgPhi:TrgEta:pTtrk:dphif");


  File_output= new TFile(outFile,"RECREATE");
  outTree = new TTree("JetTree","JetTree");
  //outTree -> Branch("Events", &JEvt, 256000, 2);
  outTree->Branch("eventIndex",&EventIndex,"EventIndex/I");
  outTree->Branch("Refmult",&Refmult,"Refmult/D");
  outTree->Branch("NJets",&NJets,"NJets/I");
  outTree->Branch("TSP",&TSP,"TSP/D");
  outTree->Branch("TrgEta",&TrgEta,"TrgEta/D");
  outTree->Branch("TrgPhi",&TrgPhi,"TrgPhi/D");
  outTree->Branch("TrgEt",&TrgEt,"TrgEt/D");
  outTree->Branch("Rho",&Rho,"Rho/D");
  outTree->Branch("Sigma",&Sigma,"Sigma/D");
  outTree->Branch("Vz",&Vz,"Vz/D");

  //______Jet Info
  //  t1->Branch("E", &E);
  outTree->Branch("JetIndex",&JetIndex);
  outTree->Branch("JetNCons",&JetNCons);//etNCons[NJets]/I");
  outTree->Branch("JetPt",&JetPt);//etPt[NJets]/D");
  outTree->Branch("JetPtCorr",&JetPtCorr);//etPtCorr[NJets]/D");
  outTree->Branch("JetEta",&JetEta);//etEta[NJets]/D");
  outTree->Branch("JetPhi",&JetPhi);//etPhi[NJets]/D");
  outTree->Branch("JetE",&JetE);//etE[NJets]/D");
  outTree->Branch("JetArea",&JetArea);//etArea[NJets]/D");


  outTree->Branch("JetConsPt",&JetConsPt);//,"JetConsPt[NJets][400]/D");
  outTree->Branch("JetConsEta",&JetConsEta);//,"JetConsEta[NJets][400]/D");
  outTree->Branch("JetConsPhi",&JetConsPhi);//,"JetConsPhi[NJets][400]/D");
  outTree->Branch("JetConsE",&JetConsE);//,"JetConsE[NJets][JetNCons[400]/D");

  outTree->SetAutoSave(-500000000);  // autosave activated for each 5 MB
  //outTree->SetMaxTreeSize(1000000000);  // new file with name _n.root will be created after 1GB size


  ///________________________________________
  Init_Histograms();



}
//___________________________

void StGJetTreeAnalysisForQA::Init_Histograms()
{

  hVpdVz_Vz =new TH2D("hVpdVz_Vz","",1000,-500,500,1000,-500,500);;
  hVpdVzVzDiff= new TH1D("hVpdVzVzDiff","",1000,-500,500);

  hRefmultPos = new TH1D("hRefmultPos","",1000,0,1000);
  hRefmultNeg= new TH1D("hRefmultNeg","",1000,0,1000);
  hRefmultPos_tofmult = new TH2D("hRefmultPos_tofmult","",1000,0,1000,5000,0,5000);
  hRefmultNeg_tofmult =new TH2D("hRefmultNeg_tofmult","",1000,0,1000,5000,0,5000);

  hprimTrk_ana= new TH1D("hprimTrk_ana","",1000,0,1000);
  hprimTrk= new TH1D("hprimTrk","",2000,0,2000);
  hprimTrk_1= new TH1D("hprimTrk_1","no p cut",2000,0,2000);

  hprimTrk_nonana = new TH1D("hprimTrk_nonana","",2000,0,2000);;
  hprimTrk_nonanapana = new TH1D("hprimTrk_nonanapana","",2000,0,2000);


  hRefmult = new TH1D("hRefmult","",1000,0,1000);
  hTSP = new TH1D("hTSP","",200,0,2);
  hTSP_pi0 = new TH1D("hTSP_pi0","",200,0,2);
  hTSP_g = new TH1D("hTSP_g","",200,0,2);

  hTwrEt = new TH1D("hTwrEt","",1000,0,100);
  hTwrTotalE = new TH1D("hTwrTotalE","",1000,0,1000);
  pTwrIdE= new TProfile("pTwrIdE","",5000,0,5000);

  h_grefmult_primTrkana = new TH2D("h_grefmult_primTrkana","",1000,0,1000,2000,0,2000);
  h_Tofmatch_primTrkana = new TH2D("h_Tofmatch_primTrkana","",1000,0,1000,2000,0,2000);

  hpTTrack = new TH1D("hpTTrack","",500,0,50);hpTTrack->Sumw2();
  hTofMatch_Refmult = new TH2D("hTofMatch_Refmult","",1000,0,1000,1000,0,1000);
  h_grefmult_tofmult = new TH2D("h_grefmult_tofmult","grefmult vs. tofmult",1000,0,1000,5000,0,5000);
  h_TrigTwr_Id = new TH1D("h_TrigTwr_Id","TrigTwr_Id",5000,0,5000);
  h_TrigBSMD_EtaStrpId = new TH1D("h_TrigBSMD_EtaStrpId","TrigBSMD_EtaStrpId",20000,0,20000);
  h_TrigBSMD_PhiStrpId = new TH1D("h_TrigBSMD_PhiStrpId","TrigBSMD_PhiStrpId",20000,0,20000);

  hEbp_elec = new TH1D("hEbp_elec","",1200,0,12);
  hEbp_hadrn= new TH1D("hEbp_hadrn","",1200,0,12);

  hpi0_Et = new TH1D("hpi0_Et","",500,0,100);  hpi0_Et->Sumw2();
  hg_Et = new TH1D("hg_Et","",500,0,100); hg_Et->Sumw2();


  HPtrig = new TH1D("HPtrig","",100,0,100);
  HGtrig = new TH1D("HGtrig","",100,0,100);

  Int_t   XBins  = 70;
  Float_t XMin   = -2.0;
  Float_t XMax   =  5.0;
  HP_ = new TH1D("H0", "Gammajet_pi", XBins, XMin, XMax);
  HG_ = new TH1D("G0", "Gammajet_g", XBins, XMin, XMax);
  
  HP[0] = (TH1D*) HP_->Clone("Pi0");  HP[0]->SetTitle("pi_12_13"); HP[0]->Sumw2();
  HP[1] = (TH1D*) HP_->Clone("Pi1");  HP[1]->SetTitle("pi_13_14");HP[1]->Sumw2();
  HP[2] = (TH1D*) HP_->Clone("Pi2");  HP[2]->SetTitle("pi_14_15");HP[2]->Sumw2();
  
  HG[0] = (TH1D*) HG_->Clone("G0");  HG[0]->SetTitle("g_12_13"); HG[0]->Sumw2();
  HG[1] = (TH1D*) HG_->Clone("G1");  HG[1]->SetTitle("g_13_14"); HG[1]->Sumw2();
  HG[2] = (TH1D*) HG_->Clone("G2");  HG[2]->SetTitle("g_14_15"); HG[2]->Sumw2();


  for(Int_t i=0; i<_N_HISTOGRAMS_; i++)
    { 
      HP[i]->Reset();     HG[i]->Reset(); 
    }

  //TH1D*   hPrimTrk[_N_HISTOGRAMS_P];
  //TH1D*   hPMTrk_;

  Int_t   XPBins  = 1000;
  Float_t XPMin   = 0;
  Float_t XPMax   = 1000;
  hPrimTrk_0 = new TH1D("hPrimTrk_0", "", XPBins, XPMin, XPMax);
  hPrimTrk_1 = new TH1D("hPrimTrk_1", "", XPBins, XPMin, XPMax);
  hPrimTrk_2 = new TH1D("hPrimTrk_2", "", XPBins, XPMin, XPMax);
  hPrimTrk_3 = new TH1D("hPrimTrk_3", "", XPBins, XPMin, XPMax);
  hPrimTrk_4 = new TH1D("hPrimTrk_4", "", XPBins, XPMin, XPMax);
  hPrimTrk_5 = new TH1D("hPrimTrk_5", "", XPBins, XPMin, XPMax);
  hPrimTrk_6 = new TH1D("hPrimTrk_6", "", XPBins, XPMin, XPMax);
  hPrimTrk_7 = new TH1D("hPrimTrk_7", "", XPBins, XPMin, XPMax);
  hPrimTrk_8 = new TH1D("hPrimTrk_8", "", XPBins, XPMin, XPMax);
  hPrimTrk_9 = new TH1D("hPrimTrk_9", "", XPBins, XPMin, XPMax);
  hPrimTrk_10 = new TH1D("hPrimTrk_10", "", XPBins, XPMin, XPMax);
  hPrimTrk_11 = new TH1D("hPrimTrk_11", "", XPBins, XPMin, XPMax);
  hPrimTrk_12 = new TH1D("hPrimTrk_12", "", XPBins, XPMin, XPMax);
  hPrimTrk_13 = new TH1D("hPrimTrk_13", "", XPBins, XPMin, XPMax);
  hPrimTrk_14 = new TH1D("hPrimTrk_14", "", XPBins, XPMin, XPMax);
  hPrimTrk_15 = new TH1D("hPrimTrk_15", "", XPBins, XPMin, XPMax);
  hPrimTrk_16 = new TH1D("hPrimTrk_16", "", XPBins, XPMin, XPMax);
  hPrimTrk_17 = new TH1D("hPrimTrk_17", "", XPBins, XPMin, XPMax);
  hPrimTrk_18 = new TH1D("hPrimTrk_18", "", XPBins, XPMin, XPMax);




}
