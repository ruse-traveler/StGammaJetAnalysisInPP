// Nihar R. Sahoo
// This is used to do some QA on event, BEMC tower and 
// after then produce same "Gfemtodst" 
// 
//  One can use this class to read Gfemtodst and do analysis ..... Jul11 2016
//____________________________________

#define StGJetTreeAnalysis_cxx
#include "StGJetTreeAnalysis.h"
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
static const int    PRIMARY_TRACK_CUT = 400;
static const double INNER_R = 0.05;
static const double OUTER_R = 0.15;
//static const Int_t Jet_type = 0;                     //           "FullJet" = 1  and "ChargedJet" = 0

//static const double Rho_R = 0.2;  // [Derek, 08.28.2022]

TNtuple *ntuple_WB;
TDatime datime;


ClassImp( StGJetTreeAnalysis);

//void StGJetTreeAnalysis::Make()
void StGJetTreeAnalysis::Make(Double_t Jet_R, Int_t Remove_N_hardest, Int_t Jet_type, double PT_TRACK_MAX, double PT_TRACK_MIN, double ETA_TRACK_MAX, double ETA_TRACK_MIN, Int_t REFMULT_CUT, TString lumi, TString sOutPref)  // [Derek, 08.28.2022]
{

  // [Derek, 08.28.2022]
  const double Rho_R = Jet_R;

  TString OutFileName;
  //  OutFileName = "JetRootFiles/pp200_";
  //OutFileName = "/star/data01/pwg/dmawxc/PythiaData/StarPythia/StarPythiaJetMaker/pp200py6";  
  //OutFileName += ".goodFiles6b.et6100pi0.";
  //OutFileName = "/gpfs01/star/pwg_tasks/hp01/output/test/pp200py6et6pt79.forJetRerun_wave1test4";
  OutFileName = sOutPref.Data();  // [Derek,08.28.2022]
  OutFileName += ".et6100pi0.";
  //  OutFileName += "CheckQA_pp_Pythia6StarTune_jettree_new_";
  OutFileName += "r0";
  OutFileName += Jet_R*10;
  OutFileName += "rm";
  OutFileName += Remove_N_hardest;
  OutFileName += "chrg.d" ;
  OutFileName +=  datime.GetDay();
  OutFileName += "m" ;
  OutFileName +=  datime.GetMonth();
  OutFileName += "y" ;
  OutFileName +=  datime.GetYear();
  OutFileName +=  ".root";

  cout<<" Output File: "<<OutFileName.Data()<<endl;

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"//      Job submitted with following setting ........ "<<endl;
  cout<<"//      Jet Radius: " << Jet_R << ", Rho Radius: " << Rho_R << endl;  // [Derek, 08.28.2022]
  cout<<"//      Remove_N_hardest: "<<Remove_N_hardest<<endl;
  cout<<"//      Jet_type: "<<Jet_type<<endl;
  cout<<"//      PT range: ["<<PT_TRACK_MIN<<" to "<<PT_TRACK_MAX<<"] GeV/c"<<endl;
  cout<<"//      ETA range: ["<<ETA_TRACK_MIN<<" to "<<ETA_TRACK_MAX<<"]"<<endl;
  cout<<".............................................................."<<endl;

  //_______________________________________________________________________
  Init_Picotree(OutFileName.Data());

  if (fChain == 0) return;

  // for trying to catch corrupted files [Derek, 08.24.2021]
  /*
  TFile   *fFileCurr;
  TString sFileList(OutFileName.Data());
  sFileList.Remove(5);
  sFileList.Append(".files.txt");
  ofstream out1(sFileList.Data());
  */
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Int_t evnt_counter=0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;

     // for trying to catch corrupted files [Derek, 08.24.2021]
     //fFileCurr = (TFile*) fChain -> GetFile();
     //out1 << fFileCurr -> GetName();

     // to catch corrupted entries [Derek, 08.24.2021]
     nb = fChain -> GetEntry(jentry);
     if (nb < 0) {
       continue;
     } else {
       nbytes += nb;
     }

      //_________________Primary Trk Info Reading_____________
      Int_t track_counter = 0;
      Int_t tofmatched =0;
      double totalPt=0.;
      Double_t annular_pTSum =0;

      evnt_counter++;

      if( evnt_counter% 1000 ==0)cout<<"processing event: "<<evnt_counter<<endl;
      //      if(evnt_counter == 300000)break;


      //___________trigger selection________
      Int_t countNtrig=0;
      double Tet(-99.);
      double Teta(-99.);
      double Tphi(-99.);
      double Tid(-99.);  // [Derek, 08.23.2022]
      Bool_t foundTrig = false;  // [Derek, 08.23.2022]
      Bool_t trig_flag = 0;
      
      vector <double> Trig_Eta;
      vector <double> Trig_Et;
      vector <double> Trig_Phi;
      vector <double> Trig_Id;  // [Derek, 08.23.2022]
      //      cout<<" mParticles_ "<<mParticles_<<endl;
      for(Int_t ii=0; ii < mParticles_; ii++)
	    {
	      // check photon (22) from hard scattering (21)
	      //	      if(mParticles_fPdgCode[ii]==22 && mParticles_fStatusCode[ii] ==21)
	      //	      if(mParticles_fPdgCode[ii]==111)
	      if(mParticles_fPdgCode[ii] == 111)
		{
	      TLorentzVector Tvec;
	      Tvec.SetPxPyPzE(mParticles_fPx[ii],mParticles_fPy[ii],mParticles_fPz[ii],mParticles_fE[ii]);

	      //	      cout<<" Status code: "<<mParticles_fStatusCode[ii]<<endl;
	      double Et1  = Tvec.Et();
	      if (!(Et1 > 6)) {
                continue;
              }

                // [Derek, 08.23.2022]
                const double eTtrig = Tvec.Et();
                const double hTrig  = Tvec.Eta();
                const double fTrig  = Tvec.Phi();
                const double idTrig = (Double_t) mParticles_fPdgCode[ii];
                // [Derek, 09.14.2022]
                if (TMath::Abs(hTrig) < 0.9) {
                  foundTrig = true;
                }
                // uncomment for debugging
                /*
                cout << "CHECK: found trigger in trigger loop!\n"
                     << "       eT = " << eTtrig << ", eta = " << hTrig << ", phi = " << fTrig << "\n"
                     << "       id = " << idTrig << ", foundTrig = " << foundTrig
                     << endl;
                */

	      Tet = Et1;//Tvec.Et();
	      Teta = Tvec.Eta();
	      Tphi = Tvec.Phi();
              Tid  = idTrig;  // [Derek, 08.23.2022]
	      //	      cout<<"pi0: "<<"Et: "<<Et<<" Teta: "<< Teta<<"  Tphi:  "<<Tphi<<"  status code: "<<mParticles_fStatusCode[ii]<<endl;
	      Trig_Eta.push_back(Teta);
	      Trig_Phi.push_back(Tphi);
	      Trig_Phi.push_back(Tet);
              Trig_Id.push_back(Tid);  // [Derek, 08.23.2022]

	      countNtrig++;

                  // the 1st pion in the event is THE trigger [Derek, 08.23.2022/09.14.2022]
                  if (foundTrig) {
                    break;
                  }
		}		 

	    }// Trigger selection loop
      
      //cout<<"Trigger size: "<<Trig_Eta.size()<<endl;

      //if(!(countNtrig ==1))continue;  // This is not an inclusive way of selecting triggers [Derek, 08.23.2022]

      // [Derek, 08.23.2022]
      /* uncomment for debugging
      if (foundTrig) {
        cout << "CHECK: trigger was found!\n"
             << "       eT = " << Tet << ", eta = " << Teta << ", phi = " << Tphi << "\n"
             << "       id = " << Tid << ", foundTrig = " << foundTrig << ", nTrig = " << countNtrig
             << endl;
      }
      */
      if (!foundTrig) continue;


      //      cout<<"mProcessId: "<<mProcessId<<endl;      

	      //cout<<"............................ inside track loop....."<<endl;
	      // Ttrigger cuts
      if(!(Teta >= -0.9 && Teta <= 0.9))continue;
      if(!(Tet > 6))continue;


      if(mProcessId == 11 || mProcessId == 12)hprocessid_1112->Fill(Tet);
      else if(mProcessId == 13)hprocessid_13->Fill(Tet);
      else if(mProcessId == 28)hprocessid_28->Fill(Tet);
      else if(mProcessId == 53)hprocessid_53->Fill(Tet);
      else if(mProcessId == 68)hprocessid_68->Fill(Tet);
      else if(mProcessId == 96)hprocessid_96->Fill(Tet);
      else hprocessid_0->Fill(Tet);
      hprocessid_all->Fill(Tet);

      //	      cout<<mProcessId<<"   "<<mParticles_<<endl;      
      //      cout<<"pi0: "<<"Tet: "<<Tet<<" Teta: "<< Teta<<"  Tphi:  "<<Tphi<<endl;
	      
      hEtTrig->Fill(Tet);      

	      //__________________
	      // for constituent tral loop
	      vector<PseudoJet> particles;
	      particles.clear();
	      
	      Int_t partCounter=0;
	      for(Int_t ij=0; ij < mParticles_; ij++)
		{
		  
		  if(!(mParticles_fStatusCode[ij] != 21))continue;

		  // final state particles
		  if(mParticles_fDaughter[ij][0]== 0 && mParticles_fDaughter[ij][1]==0)
		    {
		      // only p, pbar, K+, K-, pi+, pi-, sigma+, sigma-
		      
		      /* if((mParticles_fPdgCode[ij] != 211 || mParticles_fPdgCode[ij] != -211
                         || mParticles_fPdgCode[ij] != 321 || mParticles_fPdgCode[ij] != -321 ||
			  mParticles_fPdgCode[ij] != 2212 || mParticles_fPdgCode[ij] != -2212 ||
			  mParticles_fPdgCode[ij] != 3222 || mParticles_fPdgCode[ij] != -3112))
			*/

		      /*		      if(mParticles_fPdgCode[ij] != 211 && mParticles_fPdgCode[ij] != -211 && mParticles_fPdgCode[ij] != 2212 && mParticles_fPdgCode[ij] != -2212 && mParticles_fPdgCode[ij] != 321 && mParticles_fPdgCode[ij] != -321){cout<<mParticles_fPdgCode[ij]<<endl;}
		       */
		      if(mParticles_fPdgCode[ij] == 211 || mParticles_fPdgCode[ij] == -211
			 || mParticles_fPdgCode[ij] == 321 || mParticles_fPdgCode[ij] == -321 ||
			 mParticles_fPdgCode[ij] == 2212 || mParticles_fPdgCode[ij] == -2212 ||
			 mParticles_fPdgCode[ij] == 3222 || mParticles_fPdgCode[ij] == -3112)		      
			{
			  TLorentzVector Cvec;
			  Cvec.SetPxPyPzE(mParticles_fPx[ij],mParticles_fPy[ij],mParticles_fPz[ij],mParticles_fE[ij]);
		  
			  double Pt = Cvec.Pt();
			  double eta = Cvec.Eta();
			  double phi = Cvec.Phi();
			  //pT and Eta cut

			  if(!( eta > ETA_TRACK_MIN && eta < ETA_TRACK_MAX))continue;

			  if(!( Pt > PT_TRACK_MIN && Pt < PT_TRACK_MAX))continue;
			  //			 cout<<"PDG: "<<mParticles_fPdgCode[ij]<<" stcode "<<mParticles_fStatusCode[ij]<<" Pt: "<<Pt<<" eta: "<<eta<<" phi: "<<phi<<endl;
			  particles.push_back( PseudoJet(mParticles_fPx[ij],mParticles_fPy[ij],mParticles_fPz[ij],mParticles_fE[ij]));
			  hpTTrack->Fill(Pt);		       
			  
			  partCounter++;

                          // for recoil track check [Derek, 08.23.2022]
                          Double_t dFtrk = phi - Tphi;
                          if (dFtrk < ((-1. * Value_Pi) / 2.)) dFtrk += (2. * Value_Pi);
                          if (dFtrk > ((3. * Value_Pi) / 2.))  dFtrk -= (2. * Value_Pi);

                          const Bool_t isRecoilTrack = ((dFtrk > ((3. * Value_Pi) / 4.)) && (dFtrk < ((5. * Value_Pi) / 4.)));
                          if (isRecoilTrack) {
                            hPtRecoilTrk  -> Fill(Pt);
                            hEtaRecoilTrk -> Fill(eta);
                            hDfRecoilTrk  -> Fill(dFtrk);
                            /* uncomment for debugging
                            cout << "Check: found recoil track!\n"
                                 << "       pT = " << Pt << ", eta = " << eta << ", dF = " << dFtrk
                                 << endl;
                            */
                          }

			} // PID

		    } //  // final state particles		  
		}/// constituent tral loop	  
	      //	      cout<<" partCounter= "<<partCounter<<endl;


  
#if 1
	      //_____________________Reconstruction done here__________
	      	  
	      //__________________Jet Definitions_________
	      JetDefinition jet_def(antikt_algorithm, Jet_R);
	      
      // jet area definition
      //      Double_t ghost_maxrap = 1.0 + Jet_R; // Fiducial cut for background estimation
	      Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
	      //      GhostedAreaSpec area_spec(ghost_maxrap);
	      //AreaDefinition area_def(active_area, area_spec);
	      AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
	      
	      vector<PseudoJet> jets;      
	      // Using  Jet areas in clustering Sequence
	      
	      ClusterSequenceArea cs_hard(particles, jet_def, area_def);
	      double pt_Min = 0.2;
	      vector<PseudoJet> jets_cs = sorted_by_pt(cs_hard.inclusive_jets(pt_Min));
	      
	      Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - Jet_R);        //...... fiducial cut for jets from Alex
	      
	      // For Active area jets after Fiducial cuts
	      jets = Fiducial_cut_selector(jets_cs);
	      
	      /////////////// For Background estimator /////////
	      //      JetDefinition jet_def_bkgd(kt_algorithm, Jet_R);
	      JetDefinition jet_def_bkgd(kt_algorithm, Rho_R);
	      
      
	      // Area Definition for bkgd
	      AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
	      
	      // Selector for bkgd
	      //      Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest));
	      Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest));
	      //      Selector selector = SelectorAbsEtaMax(1.0); // to tes
      
	      // JetMedianBackgroundEstimator definition for bkgd estimation
	      JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
	      
	      //Setting Subtractor
	      Subtractor subtractor(&bkgd_estimator);
      
	      // Default for higher version of fastjet
	      //      subtractor.set_use_rho_m(true);
	      
	      //#if FASTJET_VERSION_NUMBER >= 30100
	      //      subtractor.set_use_rho_m(true);
	      //      subtractor.set_safe_mass(true);
	      //#endif
	      //Setting bkgd_estimator for events particles
	      bkgd_estimator.set_particles(particles);
	      
	      // Estimate jet Rho and sigma
	      
	      Double_t jets_rho = bkgd_estimator.rho();
	      Double_t jets_sigma = bkgd_estimator.sigma();
	      
	      // print out some info
	      //	      cout << "Clustered with " << jet_def.description() << endl;
	      // print the jets
	      
	      Int_t no0f_jets = jets.size();
	      //      cout<<" NO Of JETS# "<<no0f_jets<<endl;
	      //_________Filling Event Info
	      EventIndex      = evnt_counter;
	      //      Refmult         = refMult;
	      Refmult         = partCounter; // Total final state particles only, p,k, pi and anti-part
	      NJets           = no0f_jets;
	      TSP             = Tid;
	      TrgEta          = Teta;
	      TrgPhi          = Tphi;
	      TrgEt           = Tet;
	      Rho             = jets_rho;
	      Sigma           = jets_sigma;
	      Vz              = 0.0;
	      TofMult         =0.;
	      TofMatched      = 0;
	      PrimTrk         = partCounter; // same as Refmult
              PtHat           = mPt;   // [Derek, 09.07.2022]
	      
      

	      //_____clear vector jet info
	      JetIndex.clear();
	      JetPt.clear();
	      JetPtCorr.clear();
	      JetEta.clear();
	      JetPhi.clear();
	      JetE.clear();
	      JetArea.clear();
	      JetNCons.clear();
      

	      //_________clear vector for constituents info
	      JetConsPt.clear();
	      JetConsEta.clear();
	      JetConsPhi.clear();
	      JetConsE.clear();
      
	  
      
      for (Int_t nj = 0; nj < no0f_jets; nj++) {
	
        double jet_pT= jets[nj].perp();
        double jet_pT_sub= jets[nj].perp() - jets[nj].area()*jets_rho;
        double jet_Area = jets[nj].area();
	//        double jet_Rap = jets[nj].rap();
        double jet_PsudoRap = jets[nj].pseudorapidity();
        //      double jet_Phi = jets[nj].phi();  //jet_Phi: 0 to 2*pi
        double jet_Phi = jets[nj].phi_std(); //jet_Phi: -pi to pi
	
	vector<PseudoJet> constituents = jets[nj].constituents();
	Int_t nCon = constituents.size();
	//	cout<< "no. of nCon = "<<nCon<<endl;
	
	//	cout<<" Jet pT: "<<jet_pT<<"  Eta: "<<jet_PsudoRap<<" jet_Phi: "<<jet_Phi<<endl;
	JetIndex.push_back(nj);
	JetPt.push_back(jet_pT);
	JetPtCorr.push_back(jet_pT_sub);
	JetEta.push_back(jet_PsudoRap);
	JetPhi.push_back(jet_Phi);
	JetE.push_back(jets[nj].e());
	JetArea.push_back(jet_Area);
	JetNCons.push_back(nCon);
	

	vector <double> c_pt;
	vector <double> c_eta;
	vector <double> c_phi;
	vector <double> c_e;
	
	c_pt.clear(); 
	c_eta.clear();
	c_phi.clear();
	c_e.clear();

	/*	
#if 0
        for (Int_t jc = 0; jc < nCon; jc++) {

	  //_____________
	  double jet_consti_pt = constituents[jc].perp();
	  double jet_consti_rap = constituents[jc].rap();
	  //      double jet_consti_phi = constituents[j].phi();
	  double jet_consti_phi = constituents[jc].phi_std();
	  //	  double pXcon = constituents[jc].px();
	  //	  double pYcon = constituents[jc].py();
	  //	  double pZcon = constituents[jc].pz();
	  double jet_consti_e  = constituents[jc].e();

          // cout << "    constituent " << j << "’s pt: "<< constituents[j].perp() <<" rap  "<<constituents[j].rap()<<" phi  "<<constituents[j].phi() <<endl;

	  c_pt.push_back(jet_consti_pt);
	  c_eta.push_back(jet_consti_rap);
	  c_phi.push_back(jet_consti_phi);
	  c_e.push_back(jet_consti_e);

        } // for (unsigned j = 0; j < constituents.size(); j++) here

	JetConsPt.push_back(c_pt);
	JetConsEta.push_back(c_eta);
	JetConsPhi.push_back(c_phi);
	JetConsE.push_back(c_e);
#endif
*/

	
      }// Recoil condition
     
      outTree->Fill();


#endif	
	      
   } //event loop


   hpTTrack->Fill(1./evnt_counter);
   
  Finish();
} //Make

//______________________________________________
Bool_t   StGJetTreeAnalysis::Event_Cut(Double_t Vr, Double_t Vz)
{
  
  if( Vr < 2.0  
      && fabs(Vz) < 70
      && fabs(Vz) >= 0.0001
      ){
  return kTRUE;
  } else return kFALSE;
}
//______________________
//Bool_t StGJetTreeAnalysis::IsHighTrkPt_OK(Double_t pt, Double_t dca, Double_t nfit)
Bool_t StGJetTreeAnalysis::IsHighTrkPt_OK(Double_t pt, Double_t nsigPi0, Double_t nfit)
{
  /*
  if( pt >= 14 && nfit > 20){return kTRUE;}
  else if(pt < 14 && nfit > 15){return kTRUE;}
  else return kFALSE;
  */

  if( pt >= 15 && fabs(nsigPi0) < 2.0){return kTRUE;}
  else if(pt < 15 && nfit > 15){return kTRUE;}
  else return kFALSE;
}

//______________________________________________________________
void StGJetTreeAnalysis::Init_Picotree(const char* outFile)
{
  File_output= new TFile(outFile,"RECREATE");

  //  ntuple_WB = new TNtuple("ntuple_WB","ntuple_WB","grefmult_WB:primTrk_WB");

  outTree = new TTree("JetTree","JetTree");
  //outTree -> Branch("Events", &JEvt, 256000, 2);
  outTree->Branch("eventIndex",&EventIndex,"EventIndex/I");

  outTree->Branch("Refmult",&Refmult,"Refmult/D");
  outTree->Branch("NJets",&NJets,"NJets/I");
  outTree->Branch("TofMult",&TofMult,"TofMult/I");
  outTree->Branch("TofMatched",&TofMatched,"TofMatched/I");
  outTree->Branch("PrimTrk",&PrimTrk,"PrimTrk/I");
  outTree->Branch("TSP",&TSP,"TSP/D");
  outTree->Branch("TrgEta",&TrgEta,"TrgEta/D");
  outTree->Branch("TrgPhi",&TrgPhi,"TrgPhi/D");
  outTree->Branch("TrgEt",&TrgEt,"TrgEt/D");
  outTree->Branch("Rho",&Rho,"Rho/D");
  outTree->Branch("Sigma",&Sigma,"Sigma/D");
  outTree->Branch("Vz",&Vz,"Vz/D");
  outTree->Branch("PtHat", &PtHat, "PtHat/D");


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

void StGJetTreeAnalysis::Init_Histograms()
{

  //ProcessId__
  hprocessid_1112 = new TH1D("hprocessid_1112"," qq to qq and qqbar to qqbar",100,0,50);hprocessid_1112->Sumw2();
  hprocessid_13 = new TH1D("hprocessid_13"," qqbar to gg",100,0,50);hprocessid_13->Sumw2();
  hprocessid_28 = new TH1D("hprocessid_28"," qg to qg",100,0,50);hprocessid_28->Sumw2();
  hprocessid_53 = new TH1D("hprocessid_53"," gg to qqbar",100,0,50);hprocessid_53->Sumw2();
  hprocessid_68 = new TH1D("hprocessid_68"," gg to gg",100,0,50);hprocessid_68->Sumw2();
  hprocessid_96 = new TH1D("hprocessid_96"," semihard QCD 2 TO 2",100,0,50);hprocessid_96->Sumw2();
  hprocessid_0 = new TH1D("hprocessid_0"," unrecognized",100,0,50);hprocessid_0->Sumw2();
  hprocessid_all = new TH1D("hprocessid_all"," All",100,0,50);hprocessid_all->Sumw2();
  //___

  hTwrE = new TH1D("hTwrE","",1000,0,1000);
  hTwrTotalE = new TH1D("hTwrTotalE","",1000,0,1000);
  pTwrIdE= new TProfile("pTwrIdE","",5000,0,5000);
  pRunIndex_rho = new TProfile("pRunIndex_rho","",5000,0,5000);
  pDCA_pt = new TProfile("pDCA_pt","",200,0,100);
  hDCA_pt = new TH2D("hDCA_pt","",200,0,100,50,0,5);
  pnfit_pt = new TProfile("pnfit_pt","",200,0,100);
  hnfit_pt = new TH2D("hnfit_pt","",200,0,100,50,0,50);

  ppRunIndex_zdc =new TProfile("ppRunIndex_zdc","",5000,0,5000);
  ppRunIndex_bbc = new TProfile("ppRunIndex_bbc","",5000,0,5000);
  pRunIndex_primtrk = new TProfile("pRunIndex_primtrk","",5000,0,5000);

  hPrimaryTrk = new TH1D("hPrimaryTrk","",2000,0,2000);
  hPrimaryTrk_ana = new TH1D("hPrimaryTrk_ana","",1000,0,1000);
  hPrimary_TotalPt = new TH1D("hPrimary_TotalPt","",4000,0,2000);
  hPrimaryAnaVsTotalPt = new TH2D("hPrimaryAnaVsTotalPt","",1000,0,1000,4000,0,2000);
  hRefmult = new TH1D("hRefmult","",1000,0,1000);

  hgRefmult = new TH1D("hgRefmult","",1000,0,1000);
  hPrimTrk_0 = new TH1D("hPrimTrk_0","",1000,0,1000);
  hPrimTrk_1 = new TH1D("hPrimTrk_1","",1000,0,1000);
  hPrimTrk_2 = new TH1D("hPrimTrk_2","",1000,0,1000);
  hPrimTrk_3 = new TH1D("hPrimTrk_3","",1000,0,1000);
  hPrimTrk_4 = new TH1D("hPrimTrk_4","",1000,0,1000);
  hPrimTrk_5 = new TH1D("hPrimTrk_5","",1000,0,1000);
  hPrimTrk_6 = new TH1D("hPrimTrk_6","",1000,0,1000);
  hPrimTrk_7 = new TH1D("hPrimTrk_7","",1000,0,1000);

  hRadius_pi0                = new TH1D("hRadius_pi0","",1000,0,50);
  hRadius_gamma              = new TH1D("hRadius_gamma","",1000,0,50);
  hAnnularR                  = new TH1D("hAnnularR","",10,0,1);
  hAnnularSumPt        = new TH1D("hAnnularSumPt","",100,0,100);
  hAnnularEtaPhi_pi0   = new TH2D("hAnnularEtaPhi_pi0","",400,-2,2,750,-3.5,3.5);
  hAnnularEtaPhi_gamma = new TH2D("hAnnularEtaPhi_gamma","",400,-2,2,750,-3.5,3.5);
  hAnnularDelEtaDelPhi_pi0 = new TH2D("hAnnularDelEtaDelPhi_pi0","",800,-4,4,750,-7,7);
  hAnnularDelEtaDelPhi_gamma = new TH2D("hAnnularDelEtaDelPhi_gamma","",800,-4,4,750,-7,7);


			  


  hpTTrack = new TH1D("hpTTrack","",250,0,50);hpTTrack->Sumw2();
  hEtTrig = new TH1D("hEtTrig","",250,0,50);hEtTrig->Sumw2();
  hTofMatch_Refmult = new TH2D("hTofMatch_Refmult","",1000,0,1000,1000,0,1000);
  h_grefmult_tofmult = new TH2D("h_grefmult_tofmult","grefmult vs. tofmult",1000,0,1000,5000,0,5000);
  h_TrigTwr_Id = new TH1D("h_TrigTwr_Id","TrigTwr_Id",5000,0,5000);
  h_TrigBSMD_EtaStrpId = new TH1D("h_TrigBSMD_EtaStrpId","TrigBSMD_EtaStrpId",20000,0,20000);
  h_TrigBSMD_PhiStrpId = new TH1D("h_TrigBSMD_PhiStrpId","TrigBSMD_PhiStrpId",20000,0,20000);

  // for recoil track check [Derek, 08.23.2022]
  hPtRecoilTrk  = new TH1D("hPtRecoilTrk", "p_{T} of recoil tracks (|#Delta#varphi^{trk} - #pi| < #pi/4)", 200, 0., 100.);
  hEtaRecoilTrk = new TH1D("hEtaRecoilTrk", "#eta of recoil tracks (|#Delta#varphi^{trk} - #pi| < #pi/4)", 80, -2., 2.);
  hDfRecoilTrk  = new TH1D("hDfRecoilTrk", "#Delta#varphi of recoil tracks (|#Delta#varphi^{trk} - #pi| < #pi/4)", 360, (-2. * Value_Pi), (2. * Value_Pi));
  hPtRecoilTrk  -> Sumw2();
  hEtaRecoilTrk -> Sumw2();
  hDfRecoilTrk  -> Sumw2();

}

//___________________________

//___________________________________________________________________
Bool_t StGJetTreeAnalysis::IsGoodRun(Int_t RunId)
{

  /*
  Int_t GoodRun_200to1800[1580]={
    15084092,15084093,15085001,15085002,15085003,15085004,15085005,15085006,15085007,15085008,15085010,15085012,15085013,15085014,15085015,15085016,15085017,15085018,15085019,15085020,15085021,15085022,15085115,15086001,15086003,15086004,15086005,15086006,15086007,15086009,15086010,15086013,15086014,15086016,15086017,15086018,15086046,15086050,15086051,15086054,15086055,15086058,15086060,15086061,15086062,15086063,15086064,15086065,15086066,15086067,15086068,15086069,15086073,15086074,15086075,15086076,15086077,15086078,15086079,15086082,15087001,15087002,15087004,15087006,15087007,15087008,15087009,15087010,15087011,15087012,15087013,15087016,15087017,15087018,15087019,15087020,15087022,15087036,15087037,15087038,15087039,15087040,15087041,15087042,15087043,15087046,15087047,15087049,15087050,15087055,15087056,15087058,15088003,15088004,15088005,15088006,15089004,15089005,15089006,15089007,15089008,15089009,15089010,15089023,15089025,15089026,15089027,15089028,15089029,15089030,15089031,15089032,15089033,15089034,15089036,15089037,15089039,15089040,15089041,15089042,15089044,15089045,15089051,15089052,15089053,15090002,15090003,15090004,15090005,15090007,15090008,15090010,15090020,15090022,15090038,15090039,15090047,15090048,15090049,15090050,15090053,15090058,15090060,15090062,15090063,15090064,15090065,15090066,15091006,15091007,15091024,15091025,15091026,15091027,15091028,15091029,15091033,15091034,15091036,15091037,15091041,15091042,15091045,15092004,15092005,15092007,15092008,15092011,15092012,15092013,15092077,15092078,15092079,15093001,15093002,15093003,15093004,15093005,15093006,15093007,15093008,15093015,15093016,15093017,15093018,15093031,15093034,15093035,15093036,15093037,15093038,15093039,15093040,15093041,15093043,15093044,15093045,15093046,15093048,15093049,15093050,15093051,15093052,15093053,15093054,15093061,15093062,15093063,15093064,15094001,15094002,15094007,15094008,15094009,15094010,15094011,15094012,15094013,15094014,15094015,15094016,15094017,15094019,15094020,15094056,15094057,15094058,15094059,15094060,15094064,15094065,15094066,15094069,15094071,15094073,15095008,15095009,15095010,15095011,15095012,15095013,15095014,15095016,15095017,15095018,15095019,15095020,15095021,15095022,15095024,15095025,15095026,15095027,15095028,15095040,15095041,15095042,15095045,15095046,15095047,15095048,15095051,15095052,15096012,15096013,15096014,15096015,15096017,15096018,15096019,15096020,15096021,15096022,15096023,15096024,15096025,15096026,15096029,15096030,15096031,15096032,15096050,15096052,15096053,15096054,15096055,15096056,15096057,15096058,15097006,15097007,15097008,15097009,15097010,15097011,15097012,15097013,15097016,15097018,15097019,15097020,15097021,15097022,15097024,15097025,15097026,15097027,15097028,15097029,15097030,15097039,15097040,15097041,15097042,15097043,15097044,15097046,15097050,15097054,15097055,15097056,15097057,15097059,15097061,15097063,15098007,15098008,15098011,15098012,15098013,15098014,15098015,15098016,15098017,15098018,15098019,15098035,15098036,15098037,15098038,15098039,15098040,15098041,15098067,15098068,15098069,15098070,15098071,15098072,15098073,15098074,15100009,15100010,15100011,15100024,15100025,15100027,15100028,15100030,15100031,15100032,15100033,15100035,15100036,15100037,15100038,15100039,15100040,15100100,15100101,15100102,15100103,15100124,15101001,15101002,15101003,15101004,15101005,15101006,15101007,15101008,15101009,15101010,15101011,15101012,15101013,15101014,15101015,15101016,15101017,15101018,15101019,15101020,15101021,15101022,15101023,15101040,15101041,15101042,15101043,15101044,15101045,15101046,15101047,15101048,15101049,15101050,15101051,15102006,15102007,15102008,15102009,15102010,15102011,15102012,15102013,15102014,15102015,15102016,15102017,15102018,15102019,15102020,15102021,15102022,15102023,15102026,15102032,15102033,15102034,15102035,15102036,15102037,15102038,15102039,15102040,15102041,15102042,15102043,15102044,15102046,15102047,15102050,15102056,15102057,15102058,15102059,15102060,15102068,15103014,15103015,15103016,15103017,15103018,15103019,15103021,15103024,15103025,15103026,15103027,15103028,15103029,15103030,15103031,15103042,15103043,15103045,15103046,15103049,15103050,15103051,15103052,15103053,15103054,15103055,15103056,15103057,15103058,15103059,15103060,15103062,15104002,15104003,15104004,15104006,15104007,15104008,15104010,15104011,15104012,15104013,15104017,15104019,15104020,15104021,15104022,15104023,15104024,15104037,15104040,15104042,15104043,15104044,15104052,15104053,15104054,15104055,15104056,15104057,15104058,15104060,15104063,15104065,15104066,15104067,15104068,15104069,15104070,15105001,15105002,15105003,15105004,15105005,15105006,15105008,15105010,15105012,15105013,15105015,15105016,15105017,15105018,15105019,15105020,15105021,15105023,15105024,15105025,15105026,15105028,15105030,15105031,15105032,15105054,15105055,15105056,15105057,15105058,15105061,15105062,15105063,15105064,15105065,15105068,15105070,15105071,15105072,15105073,15106001,15106002,15106005,15106008,15106009,15106010,15106011,15106131,15106132,15106133,15106134,15106135,15106136,15107001,15107002,15107003,15107004,15107005,15107006,15107008,15107009,15107010,15107011,15107013,15107014,15107063,15107064,15107065,15107073,15107074,15107075,15107076,15107078,15107079,15107080,15107081,15107082,15107084,15107085,15107086,15107087,15107089,15107090,15108001,15108003,15108004,15108005,15108007,15108008,15108009,15108010,15108011,15108012,15108013,15108014,15108015,15108016,15108017,15108018,15108019,15108021,15108022,15108023,15108024,15108025,15108027,15108028,15108056,15108057,15108058,15108059,15108063,15108069,15108072,15108073,15108074,15108075,15108076,15108077,15108078,15108079,15108080,15109001,15109002,15109004,15109005,15109006,15109007,15109008,15109009,15109010,15109011,15109012,15109013,15109019,15109020,15109021,15109024,15109025,15109026,15109027,15109028,15109029,15109030,15109031,15109032,15109035,15109036,15109037,15109038,15109039,15109040,15109041,15109042,15109043,15109045,15109046,15109051,15109052,15109054,15109055,15109056,15109057,15109058,15109059,15109060,15109061,15109062,15109063,15109065,15110001,15110002,15110003,15110004,15110005,15110008,15110009,15110010,15110011,15110012,15110013,15110015,15110016,15110017,15110019,15110020,15110021,15110022,15110023,15110024,15110025,15110026,15110027,15110028,15110029,15110030,15110031,15110033,15110034,15110035,15110039,15110040,15110041,15110042,15110043,15110045,15110046,15110047,15110048,15110049,15110050,15110051,15110052,15110053,15110058,15111001,15111002,15111003,15111007,15111008,15111009,15111010,15111011,15111012,15111013,15111014,15111015,15111016,15111050,15111051,15111052,15111053,15111065,15111066,15111067,15111068,15112001,15112005,15112006,15112007,15112009,15112010,15112011,15112012,15112013,15112014,15112019,15112020,15112021,15112022,15112026,15112027,15112028,15112029,15112030,15112031,15112033,15112034,15112035,15112036,15112037,15112038,15112039,15112040,15112041,15112042,15112043,15112044,15112045,15112046,15112047,15112048,15112049,15112050,15112051,15113001,15113004,15113005,15113006,15113009,15113010,15113011,15113012,15113103,15113104,15114001,15114002,15114003,15114010,15114012,15114017,15114018,15114019,15114020,15114021,15114022,15114025,15114026,15114027,15114028,15114034,15114035,15114036,15114037,15114038,15114039,15114040,15114041,15114042,15114043,15114044,15114045,15114048,15114050,15114051,15114052,15114053,15114054,15114055,15114056,15114057,15114058,15115012,15115013,15115014,15115067,15115068,15115069,15115070,15115071,15115072,15115073,15115074,15115075,15115076,15115077,15115078,15115079,15115080,15115082,15115083,15115084,15115085,15115086,15115087,15115088,15116002,15116003,15116004,15116015,15116016,15116017,15116018,15116019,15116020,15116021,15116026,15116029,15116032,15116033,15116034,15116035,15116036,15116037,15116039,15116040,15116041,15116042,15116047,15116048,15116049,15116050,15116051,15116052,15116054,15116055,15116056,15116057,15116058,15116059,15116060,15116061,15116062,15116068,15117003,15117004,15117005,15117050,15117052,15117053,15117054,15117055,15117056,15117057,15117058,15117059,15117060,15117061,15117062,15117063,15117064,15117065,15117066,15117067,15117068,15117069,15117070,15118003,15118004,15118005,15118006,15118007,15118008,15118009,15118011,15118013,15118014,15118015,15118016,15118017,15118018,15118019,15118020,15118021,15118022,15118023,15118024,15118052,15118053,15118054,15118055,15118056,15118057,15118058,15118059,15118060,15118063,15118064,15118065,15118066,15119008,15119009,15119011,15119012,15119013,15119014,15119015,15119018,15119019,15119020,15119022,15119023,15119024,15119029,15119031,15119034,15119035,15119036,15119043,15119045,15119046,15119049,15119051,15119052,15119053,15119054,15119055,15119065,15120001,15120002,15120004,15120005,15120006,15120007,15120008,15120009,15120011,15120105,15120106,15120108,15120109,15120115,15121001,15121002,15121003,15121004,15121005,15121008,15121009,15121012,15121013,15121015,15121016,15121017,15121018,15121061,15121063,15121065,15121066,15121067,15121068,15121070,15121071,15121072,15121076,15121077,15121078,15122003,15122004,15122006,15122008,15122010,15122011,15122017,15122018,15122019,15122020,15122023,15122026,15122030,15122031,15122034,15122039,15122041,15122042,15122043,15122044,15122045,15122049,15122058,15122059,15122060,15122061,15122062,15122063,15122064,15122065,15123001,15123009,15123010,15123011,15123015,15123016,15123017,15123018,15123019,15123020,15123021,15123022,15123027,15123028,15123034,15123035,15123036,15123037,15123044,15123045,15123047,15123048,15123049,15123053,15123054,15124001,15124002,15124003,15124004,15124006,15124008,15124010,15124021,15124022,15124023,15124025,15124026,15124027,15124028,15124031,15124032,15124033,15124035,15124040,15124041,15124042,15124043,15124044,15124049,15124050,15124052,15124053,15124054,15124055,15124056,15124057,15124063,15125001,15125007,15125060,15125061,15125062,15125064,15125065,15125067,15125074,15125075,15126002,15126004,15126006,15126009,15126010,15126011,15126012,15126013,15126015,15126016,15126017,15126021,15126022,15126023,15126036,15126037,15126038,15126039,15126040,15126044,15126045,15126046,15126047,15126048,15126050,15126051,15126052,15126060,15126061,15127001,15127002,15127003,15127004,15127005,15127006,15127007,15127011,15127012,15127013,15128015,15128016,15128017,15128018,15128019,15128020,15128021,15128022,15128024,15128025,15128027,15128028,15128029,15128030,15128031,15128032,15128033,15128048,15128049,15128050,15128051,15128052,15129001,15129002,15129003,15129006,15129007,15129009,15129010,15129013,15129014,15129015,15129016,15129017,15129018,15129022,15129045,15129046,15129047,15129048,15129049,15129050,15129051,15129052,15130001,15130002,15130003,15130004,15130005,15130006,15130007,15130008,15130009,15130011,15130013,15130014,15130015,15130016,15130031,15130032,15130033,15130034,15130035,15130036,15130037,15130038,15130040,15130041,15130042,15130043,15130044,15130045,15130046,15130047,15130048,15130054,15130055,15130056,15131001,15131002,15131003,15131004,15131008,15131009,15131010,15131011,15131012,15131013,15131014,15131030,15131031,15131033,15131034,15131035,15131036,15131039,15131041,15131042,15131043,15131044,15131045,15131046,15131050,15131051,15132002,15132003,15132004,15132005,15132007,15132008,15132009,15132012,15132014,15132015,15132016,15132017,15132018,15132022,15132023,15132024,15132025,15132030,15132031,15132032,15132033,15132034,15132038,15133002,15133004,15133005,15133006,15133010,15133011,15133012,15133017,15133020,15133031,15133032,15133033,15133034,15133035,15133039,15133040,15133041,15133042,15133044,15133047,15133049,15133050,15133051,15133053,15133054,15133055,15134002,15134006,15134007,15134008,15134009,15134010,15134050,15134051,15134053,15134055,15134056,15134057,15134058,15134060,15134062,15134063,15135004,15135005,15135006,15135007,15135008,15135009,15135011,15135012,15135013,15135014,15135015,15135016,15135062,15135063,15135064,15135065,15136001,15136002,15136003,15136004,15136006,15136007,15136008,15136009,15136011,15136012,15136013,15136014,15136015,15136016,15136017,15136033,15136034,15136035,15136036,15136037,15136038,15136039,15136040,15136041,15136042,15137025,15137026,15137027,15137028,15137029,15137030,15137032,15137033,15137034,15137035,15137036,15137037,15137038,15137050,15137051,15137052,15137053,15138035,15138056,15138057,15138058,15138059,15138060,15138061,15138062,15138063,15138064,15138065,15138066,15138067,15138068,15138069,15138070,15138071,15138074,15139003,15139008,15139009,15139010,15139012,15139013,15139016,15139017,15139018,15139020,15139021,15139022,15139026,15139038,15139040,15139041,15139042,15139043,15139044,15139046,15139047,15139050,15139051,15139052,15140002,15140003,15140004,15140005,15140006,15140007,15140011,15140012,15140013,15140014,15140015,15140017,15140018,15140019,15140020,15140021,15140022,15140023,15140025,15140026,15140027,15140028,15140029,15140030,15140031,15140032,15140034,15140035,15140036,15140037,15140045,15140046,15140047,15140048,15140049,15140050,15140051,15140052,15140055,15141001,15141002,15141003,15141004,15141006,15141007,15141009,15141010,15141011,15141012,15141013,15141048,15142003,15142004,15142005,15142008,15142011,15142012,15142015,15142018,15142020,15142021,15142022,15142025,15142055,15142058,15143010,15143011,15143012,15143013,15143026,15143027,15143028,15143029,15143030,15143031,15143032,15143033,15143034,15143035,15143036,15143037,15143038,15143039,15143040,15143049,15143051,15143052,15144001,15144004,15144005,15144006,15144007,15144008,15144009,15144010,15144011,15144012,15144013,15144014,15144015,15144017,15144018,15144028,15144033,15144034,15144035,15144036,15144037,15144039,15144040,15144043,15144044,15144049,15144052,15144054,15144056,15144058,15144060,15144062,15144063,15144064,15144065,15144067,15144068,15144070,15145003,15145004,15145005,15145006,15145007,15145008,15145015,15145016,15145017,15145018,15145019,15145022,15145023,15145024,15145029,15145034,15145035,15145036,15145038,15145056,15145060,15145061,15145062,15145063,15145064};
  */


  Int_t GoodRun_1800abv[523]={15146002,15146004,15146008,15146009,15146010,15146013,15146017,15146020,15146023,15146025,15146026,15148062,15148064,15148065,15148066,15149009,15149010,15149012,15149013,15149015,15149016,15149017,15149062,15149063,15149064,15149067,15149069,15149070,15149071,15149072,15149074,15149075,15149076,15150001,15150004,15150005,15150009,15150010,15150011,15150014,15150015,15150048,15150056,15150057,15150058,15150059,15150060,15150061,15150062,15150063,15150064,15150079,15150080,15150081,15150082,15151005,15151006,15151007,15151008,15151011,15151012,15151013,15151014,15151015,15151017,15151018,15151019,15151020,15151021,15151022,15151023,15151024,15151032,15151033,15151034,15151035,15151036,15151038,15151039,15151040,15151041,15151044,15151045,15151049,15151050,15151051,15151052,15151053,15151054,15151055,15151056,15151058,15151061,15151062,15152001,15152002,15152003,15152004,15152005,15152006,15152007,15152008,15152009,15152010,15152011,15152012,15152013,15152015,15152016,15152017,15152018,15152020,15152031,15152032,15152033,15152034,15152036,15152038,15152039,15152040,15152041,15152042,15152043,15152044,15152046,15152047,15152048,15152049,15152050,15152051,15152052,15152053,15152054,15152055,15152057,15152058,15153001,15153002,15153003,15153004,15153006,15153007,15153008,15153009,15153010,15153011,15153012,15153013,15153014,15153015,15153018,15153019,15153022,15153034,15153035,15153036,15153037,15153038,15153040,15153041,15153042,15153043,15153044,15153045,15153046,15153047,15153048,15153049,15153050,15153052,15153053,15153054,15153055,15153058,15154001,15154002,15154003,15154008,15154009,15154010,15154011,15154012,15154013,15154014,15154015,15154016,15154017,15154018,15154021,15154022,15154023,15154024,15155004,15155005,15155006,15155007,15155008,15155009,15155010,15156001,15156002,15156003,15156004,15156007,15156010,15156015,15156017,15156019,15156020,15156021,15156022,15156023,15156024,15156039,15156040,15156041,15157006,15157007,15157009,15157011,15157012,15157013,15157015,15157016,15157017,15157034,15157039,15157045,15157046,15157047,15157049,15157050,15157053,15157054,15157055,15157056,15157057,15157058,15157059,15157060,15157061,15158001,15158024,15158025,15158026,15158027,15158034,15158035,15158036,15158037,15158038,15158039,15158040,15158057,15158058,15158059,15158060,15158061,15158062,15158063,15158064,15159001,15159002,15159004,15159005,15159006,15159007,15159008,15159009,15159016,15159017,15159018,15159019,15159024,15159025,15159026,15159027,15159028,15159029,15159030,15159031,15159032,15159033,15159034,15159035,15159036,15159039,15159040,15159041,15159042,15159043,15159044,15159045,15159046,15160001,15160002,15160003,15160004,15160005,15160007,15160029,15160030,15160032,15160033,15160034,15160035,15160036,15160037,15160038,15160039,15160040,15160041,15160042,15160043,15160044,15160045,15160046,15160047,15160048,15160050,15161001,15161002,15161006,15161007,15161008,15161009,15161010,15161011,15161012,15161013,15161015,15161016,15161019,15161020,15161021,15161022,15161037,15161051,15161055,15161056,15161057,15161058,15161059,15161060,15161061,15161062,15161063,15161064,15161068,15161069,15161070,15161071,15161072,15162001,15162004,15162005,15162006,15162010,15162011,15162012,15162014,15162018,15162019,15162020,15162021,15162022,15162023,15162024,15162025,15162027,15162028,15162029,15162030,15162031,15162032,15162041,15162043,15162044,15162045,15162046,15162047,15162048,15162049,15162050,15162051,15162053,15162054,15162055,15163001,15163002,15163003,15163005,15163006,15163007,15163008,15163009,15163010,15163019,15163020,15163021,15163022,15163023,15163024,15163025,15163026,15163027,15163028,15163029,15163033,15163035,15163054,15163055,15163056,15163057,15163058,15163059,15163060,15163061,15163062,15164001,15164002,15164023,15164024,15164028,15164030,15164032,15164033,15164034,15164036,15164037,15164039,15164042,15164043,15164044,15164045,15164046,15164047,15164048,15164060,15164061,15164062,15164063,15164064,15164065,15164068,15164069,15164070,15164071,15165001,15165002,15165003,15165008,15165009,15165021,15165022,15165023,15165024,15165025,15165026,15165031,15165035,15165036,15165037,15165038,15165039,15165040,15165041,15165042,15165043,15165044,15165045,15165046,15165047,15165049,15165050,15165051,15165052,15165053,15165054,15165055,15165056,15165057,15165058,15165059,15166004,15166005,15166006,15166007,15166008,15166010,15166013,15166015,15166016,15166017,15166024,15166025,15166026,15166027,15166030,15166031,15166032,15166033,15166034,15166035,15166036,15166040,15166041,15166042,15166043,15166044,15166045,15166046,15166047,15167001,15167002,15167006,15167007,15167008,15167009,15167011,15167012,15167013,15167014};






Bool_t  x  = kFALSE;

//for(Int_t i= 0; i < 1580; i++){if(RunId == GoodRun_200to1800[i])     //badrunlist[i])
for(Int_t i= 0; i < 523; i++){if(RunId == GoodRun_1800abv[i])     //badrunlist[i])
    { x = kTRUE ;break;}}

return x;
}

