#ifndef StGJetTreeAnalysis_h
#define StGJetTreeAnalysis_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TNtuple.h"
//#include "GammaJetEvent.h"
// Header file for the classes stored in the TTree if any.
#include <TObject.h>
#include  <TH1D.h>
#include  <TH2D.h>
#include  <TProfile.h>
#include <vector>
//#include "StJetEvent.h"
#include "GammaJetTowerUtil.h"

class  TLorentzVector;
//class  StJetEvent;
class GammaJetTowerUtil;
class GammaJetQAUtil;

using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxPrimaryTrackArray = 100000;
const Int_t kMaxTowerArray = 10000000;
const Int_t kMaxTowerArray_fMatchedTracksArray = 10000;
const Int_t nJ = 500;
const Int_t nJConst = 800;

class StGJetTreeAnalysis {
 public :
   //TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TChain         *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // an event with some jets
   // StJetEvent JEvt;

   // Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxmParticles = 100000;

   // Declaration of leaf types
 //StPythiaEvent   *PythiaBranch;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           mRunId;
   Int_t           mEventId;
   Int_t           mProcessId;
   Int_t           mTune;
   UInt_t          mVertex_fUniqueID;
   UInt_t          mVertex_fBits;
   Double_t        mVertex_fX;
   Double_t        mVertex_fY;
   Double_t        mVertex_fZ;
   Float_t         mS;
   Float_t         mT;
   Float_t         mU;
   Float_t         mPt;
   Float_t         mCosTheta;
   Float_t         mX1;
   Float_t         mX2;
   Int_t           mMstu72;
   Int_t           mMstu73;
   Int_t           mMstp111;
   Float_t         mPartonALL;
   Float_t         mDF1[34];
   Float_t         mDF2[34];
   Float_t         mF1[2];
   Float_t         mF2[2];
   Int_t           mParticles_;
   UInt_t          mParticles_fUniqueID[kMaxmParticles];   //[mParticles_]
   UInt_t          mParticles_fBits[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineColor[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineStyle[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineWidth[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fPdgCode[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fStatusCode[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fMother[kMaxmParticles][2];   //[mParticles_]
   Int_t           mParticles_fDaughter[kMaxmParticles][2];   //[mParticles_]
   Float_t         mParticles_fWeight[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fCalcMass[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPx[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPy[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPz[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fE[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVx[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVy[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVz[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVt[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPolarTheta[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPolarPhi[kMaxmParticles];   //[mParticles_]


   // List of branches
   TBranch        *b_PythiaBranch_fUniqueID;   //!
   TBranch        *b_PythiaBranch_fBits;   //!
   TBranch        *b_PythiaBranch_mRunId;   //!
   TBranch        *b_PythiaBranch_mEventId;   //!
   TBranch        *b_PythiaBranch_mProcessId;   //!
   TBranch        *b_PythiaBranch_mTune;   //!
   TBranch        *b_PythiaBranch_mVertex_fUniqueID;   //!
   TBranch        *b_PythiaBranch_mVertex_fBits;   //!
   TBranch        *b_PythiaBranch_mVertex_fX;   //!
   TBranch        *b_PythiaBranch_mVertex_fY;   //!
   TBranch        *b_PythiaBranch_mVertex_fZ;   //!
   TBranch        *b_PythiaBranch_mS;   //!
   TBranch        *b_PythiaBranch_mT;   //!
   TBranch        *b_PythiaBranch_mU;   //!
   TBranch        *b_PythiaBranch_mPt;   //!
   TBranch        *b_PythiaBranch_mCosTheta;   //!
   TBranch        *b_PythiaBranch_mX1;   //!
   TBranch        *b_PythiaBranch_mX2;   //!
   TBranch        *b_PythiaBranch_mMstu72;   //!
   TBranch        *b_PythiaBranch_mMstu73;   //!
   TBranch        *b_PythiaBranch_mMstp111;   //!
   TBranch        *b_PythiaBranch_mPartonALL;   //!
   TBranch        *b_PythiaBranch_mDF1;   //!
   TBranch        *b_PythiaBranch_mDF2;   //!
   TBranch        *b_PythiaBranch_mF1;   //!
   TBranch        *b_PythiaBranch_mF2;   //!
   TBranch        *b_PythiaBranch_mParticles_;   //!
   TBranch        *b_mParticles_fUniqueID;   //!
   TBranch        *b_mParticles_fBits;   //!
   TBranch        *b_mParticles_fLineColor;   //!
   TBranch        *b_mParticles_fLineStyle;   //!
   TBranch        *b_mParticles_fLineWidth;   //!
   TBranch        *b_mParticles_fPdgCode;   //!
   TBranch        *b_mParticles_fStatusCode;   //!
   TBranch        *b_mParticles_fMother;   //!
   TBranch        *b_mParticles_fDaughter;   //!
   TBranch        *b_mParticles_fWeight;   //!
   TBranch        *b_mParticles_fCalcMass;   //!
   TBranch        *b_mParticles_fPx;   //!
   TBranch        *b_mParticles_fPy;   //!
   TBranch        *b_mParticles_fPz;   //!
   TBranch        *b_mParticles_fE;   //!
   TBranch        *b_mParticles_fVx;   //!
   TBranch        *b_mParticles_fVy;   //!
   TBranch        *b_mParticles_fVz;   //!
   TBranch        *b_mParticles_fVt;   //!
   TBranch        *b_mParticles_fPolarTheta;   //!
   TBranch        *b_mParticles_fPolarPhi;   //!
   
   //_____________________________
   ///////////////////////
   TFile          *File_output;
   const char     *outFile;
   TBranch        *eventBranch;
   TTree          *outTree;
   const char     *fileBaseName;
   const char     *mFileIndex;
   const char     *mOutDir;


   ///________

   Int_t       EventIndex;
   Double_t    Refmult;
   Int_t       NJets;
   Int_t       TofMult;
   Int_t       TofMatched;
   Double_t    TSP;
   Double_t    TrgEta;
   Double_t    TrgPhi;
   Double_t    TrgEt;
   Double_t    Rho;
   Double_t    Sigma;
   Double_t    Vz;
   Int_t       PrimTrk;
   Double_t    AnuLarSumPt;
   //   Double_t    TSP;
   Double_t    PtHat;  // [Derek, 09.07.2022]

   vector<int> JetIndex; 
   vector<int> JetNCons;
   vector<double> JetPt;
   vector<double> JetPtCorr;
   vector<double> JetEta;
   vector<double> JetPhi;
   vector<double> JetE;// = new vector<double>; 
   vector<double> JetArea;// = new vector<double>; 

   vector < vector<double> > JetConsPt;
   vector <vector<double> > JetConsEta;
   vector <vector <double> > JetConsPhi;
   vector <vector <double> > JetConsE;

   //___________________________
   TProfile *pDCA_pt;
   TH2D *hDCA_pt;
   TProfile *pnfit_pt;
   TH2D *hnfit_pt;
   TNtuple *ntuple_Trig;
   TH1D *hTwrE;
   TH1D *hPrimaryTrk;
   TH1D *hPrimaryTrk_ana;
   TH1D *hpTTrack;
   TH1D *hEtTrig;
   TH1D *hPrimary_TotalPt;
   TH2D *hPrimaryAnaVsTotalPt;
   TH1D *hRefmult;

   // for recoil track check [Derek, 08.23.2022]
   TH1D *hPtRecoilTrk;
   TH1D *hEtaRecoilTrk;
   TH1D *hDfRecoilTrk;

   TH1D *hgRefmult;

   TH1D *hPrimTrk;
   TH1D *hPrimTrk_0;
   TH1D *hPrimTrk_1;
   TH1D *hPrimTrk_2;
   TH1D *hPrimTrk_3;
   TH1D *hPrimTrk_4;
   TH1D *hPrimTrk_5;
   TH1D *hPrimTrk_6;
   TH1D *hPrimTrk_7;

   //Isolation hist
   TH1D *hRadius_pi0;
   TH1D *hRadius_gamma;
   TH1D *hAnnularR;
   TH1D *hAnnularSumPt;
   TH2D *hAnnularEtaPhi_pi0;
   TH2D *hAnnularEtaPhi_gamma;
   TH2D *hAnnularDelEtaDelPhi_pi0;
   TH2D *hAnnularDelEtaDelPhi_gamma;


   TH1D *hTwrTotalE;
   TProfile *pTwrIdE;
   TProfile *pRunIndex_rho;
   TProfile *pRunIndex_primtrk;
   TProfile *ppRunIndex_zdc;
   TProfile *ppRunIndex_bbc;
   TH2D *h_grefmult_tofmult;
   TH2D *hTofMatch_Refmult;
   TH1D *h_TrigTwr_Id;// = new TH1D("h_TrigTwr_Id","TrigTwr_Id",5000,0,5000);
   TH1D *h_TrigBSMD_EtaStrpId;// = new TH1D("h_TrigBSMD_EtaStrpId","TrigBSMD_EtaStrpId",20000,0,20000);
   TH1D *h_TrigBSMD_PhiStrpId;// = new TH1D("h_TrigBSMD_PhiStrpId","TrigBSMD_PhiStrpId",20000,0,20000);
   //____Process ID
   TH1D* hprocessid_1112;
   TH1D* hprocessid_13 ;
   TH1D* hprocessid_28 ;
   TH1D* hprocessid_53 ;
   TH1D* hprocessid_68 ;
   TH1D* hprocessid_96 ;
   TH1D* hprocessid_0 ;
   TH1D* hprocessid_all;

   //____________________________________________
   //StGJetTreeAnalysis(TTree *tree=0);
   StGJetTreeAnalysis(TChain *tree=0, TString sInList="");  // [Derek, 08.28.2022]
   virtual ~StGJetTreeAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Bool_t   Event_Cut(Double_t Vr, Double_t Vz);  /// For event Cut
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   //virtual void     Init(TTree *tree);
   virtual void     Init(TChain *tree);
   virtual void     Init_Picotree(const char *outFile);
   //   virtual void     Make();
   virtual  void     Make(Double_t R, Int_t Remove_N_hardest, Int_t Jet_type, double PT_TRACK_MAX, double PT_TRACK_MIN, double ETA_TRACK_MAX, double ETA_TRACK_MIN, Int_t REFMULT_CUT, TString lumi, TString sOutPref);
   virtual void     Finish();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Bool_t    IsGoodRun(Int_t RunId);
   virtual Double_t GetHadronicCorrection(Double_t TwrE, vector <double> v);
   virtual Double_t Get_TSP( Float_t E,
		       Float_t en4, Float_t enp4,
		       Float_t en3, Float_t enp3, Float_t en5, Float_t enp5,
		       Float_t en2,Float_t enp2, Float_t en6, Float_t enp6,
		       Float_t en1, Float_t enp1, Float_t en7, Float_t enp7,
		       Float_t en0, Float_t enp0, Float_t en8, Float_t enp8,
		       Float_t en01, Float_t enp01, Float_t en9, Float_t enp9,
		       Float_t en02, Float_t enp02, Float_t en10, Float_t enp10,
		       Float_t en03, Float_t enp03, Float_t en11, Float_t enp11);
   virtual Bool_t  IsHighTrkPt_OK(Double_t pt, Double_t dca, Double_t nfit);
   //   virtual void     DoJetReconstruction();
   virtual Double_t GetRadiusInPhaseSpace(Double_t trig_eta, Double_t trig_phi, Double_t trk_eta, Double_t trk_phi);
   virtual void Init_Histograms();

   ClassDef( StGJetTreeAnalysis,1)
};

#endif

#ifdef StGJetTreeAnalysis_cxx
//_____________________________________
Double_t StGJetTreeAnalysis::GetRadiusInPhaseSpace(Double_t trig_eta, Double_t trig_phi, Double_t trk_eta, Double_t trk_phi)
{

  Double_t radius = sqrt( pow(trig_eta - trk_eta,2) + pow(trig_phi - trk_phi,2));
  return radius;
}

//________________________________________________
Double_t StGJetTreeAnalysis::GetHadronicCorrection(Double_t TwrE, vector <double> v)
{
  
  double sum_e = 0.;
  int n = v.size();
  for(int i =0; i < n; i++)    {      sum_e += v.at(i);    }
  
  double Crr = TwrE - sum_e;
  double Corr_E =  Crr < 0 ? 0 : Crr;
  
  return Corr_E;

};
//______________________________________
void StGJetTreeAnalysis::Finish()
{

  //  File_output = ntuple_Trig->GetCurrentFile();
  File_output = outTree->GetCurrentFile();
  File_output->Write();
  File_output->Close();

}
//_____________________________________________________________________
//StGJetTreeAnalysis::StGJetTreeAnalysis(TTree *tree) : fChain(0) 
StGJetTreeAnalysis::StGJetTreeAnalysis(TChain *tree, TString sInList) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/star/data01/pwg/nihar/JetReco_pico/ppFullJet/output/sched_0DB13B7E2A6F3A1A753D077104B1B376_2.root");
     //if (!f || !f->IsOpen()) {
     //  f = new TFile("/star/data01/pwg/nihar/JetReco_pico/ppFullJet/output/sched_0DB13B7E2A6F3A1A753D077104B1B376_2.root");
     //      }
     //      f->GetObject("GTree",tree);
     
     //     TChain 

     TChain *chain = new TChain("PythiaTree"," ");


     //__________________________________________________

     ifstream in1;
     in1.open(sInList.Data());  // [Derek, 08.28.2022]
     //in1.open("./pp200py6et6pt79.forJetRerun_wave1test.list");
     //     in1.open(inputfilelist);

     // uncomment for debugging [Derek, 08.28.2022]
     cout << "[CHECK] input list = " << sInList.Data() << endl;
     //

     TString path;
     Int_t nline=0;
     while (1)
       {

         in1 >> path;
         cout<<path<<endl;

         //_________________
         TFile *f = new TFile(path.Data());
         if(!(f->IsZombie()))
           {
             chain->Add(path.Data());
           }
         else{cout<<path.Data()<<" yes Zombie..."<<endl;}

         if (!in1.good()) break;
         nline++;
       }

     /*
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_0.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_10.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_11.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_12.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_13.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_14.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_15.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_16.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_17.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_18.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_19.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_1.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_20.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_21.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_22.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_23.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_24.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_25.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_26.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_27.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_28.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_29.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_2.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_30.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_31.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_32.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_33.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_34.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_3.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_4.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_5.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_6.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_7.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_8.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_9.root");

     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_35.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_37.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_38.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_39.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_40.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_41.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_42.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_43.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_44.root");

     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_.root");
     chain->Add("/star/data01/pwg/nihar/Pythia6_Htrieg/dijet_rootfiles_t/pythia6_pi0Trig_pTHat4to100_.root");
     */

 


     
     //___________________________________
     tree=chain;


     
   }
   Init(tree);

   //   Init_Picotree("pp_200_test.root");  ");");");
}

StGJetTreeAnalysis::~StGJetTreeAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StGJetTreeAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t StGJetTreeAnalysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

//void StGJetTreeAnalysis::Init(TTree *tree)
void StGJetTreeAnalysis::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_PythiaBranch_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_PythiaBranch_fBits);
   fChain->SetBranchAddress("mRunId", &mRunId, &b_PythiaBranch_mRunId);
   fChain->SetBranchAddress("mEventId", &mEventId, &b_PythiaBranch_mEventId);
   fChain->SetBranchAddress("mProcessId", &mProcessId, &b_PythiaBranch_mProcessId);
   fChain->SetBranchAddress("mTune", &mTune, &b_PythiaBranch_mTune);
   fChain->SetBranchAddress("mVertex.fUniqueID", &mVertex_fUniqueID, &b_PythiaBranch_mVertex_fUniqueID);
   fChain->SetBranchAddress("mVertex.fBits", &mVertex_fBits, &b_PythiaBranch_mVertex_fBits);
   fChain->SetBranchAddress("mVertex.fX", &mVertex_fX, &b_PythiaBranch_mVertex_fX);
   fChain->SetBranchAddress("mVertex.fY", &mVertex_fY, &b_PythiaBranch_mVertex_fY);
   fChain->SetBranchAddress("mVertex.fZ", &mVertex_fZ, &b_PythiaBranch_mVertex_fZ);
   fChain->SetBranchAddress("mS", &mS, &b_PythiaBranch_mS);
   fChain->SetBranchAddress("mT", &mT, &b_PythiaBranch_mT);
   fChain->SetBranchAddress("mU", &mU, &b_PythiaBranch_mU);
   fChain->SetBranchAddress("mPt", &mPt, &b_PythiaBranch_mPt);
   fChain->SetBranchAddress("mCosTheta", &mCosTheta, &b_PythiaBranch_mCosTheta);
   fChain->SetBranchAddress("mX1", &mX1, &b_PythiaBranch_mX1);
   fChain->SetBranchAddress("mX2", &mX2, &b_PythiaBranch_mX2);
   fChain->SetBranchAddress("mMstu72", &mMstu72, &b_PythiaBranch_mMstu72);
   fChain->SetBranchAddress("mMstu73", &mMstu73, &b_PythiaBranch_mMstu73);
   fChain->SetBranchAddress("mMstp111", &mMstp111, &b_PythiaBranch_mMstp111);
   fChain->SetBranchAddress("mPartonALL", &mPartonALL, &b_PythiaBranch_mPartonALL);
   fChain->SetBranchAddress("mDF1[34]", mDF1, &b_PythiaBranch_mDF1);
   fChain->SetBranchAddress("mDF2[34]", mDF2, &b_PythiaBranch_mDF2);
   fChain->SetBranchAddress("mF1[2]", mF1, &b_PythiaBranch_mF1);
   fChain->SetBranchAddress("mF2[2]", mF2, &b_PythiaBranch_mF2);
   fChain->SetBranchAddress("mParticles", &mParticles_, &b_PythiaBranch_mParticles_);
   fChain->SetBranchAddress("mParticles.fUniqueID", mParticles_fUniqueID, &b_mParticles_fUniqueID);
   fChain->SetBranchAddress("mParticles.fBits", mParticles_fBits, &b_mParticles_fBits);
   fChain->SetBranchAddress("mParticles.fLineColor", mParticles_fLineColor, &b_mParticles_fLineColor);
   fChain->SetBranchAddress("mParticles.fLineStyle", mParticles_fLineStyle, &b_mParticles_fLineStyle);
   fChain->SetBranchAddress("mParticles.fLineWidth", mParticles_fLineWidth, &b_mParticles_fLineWidth);
   fChain->SetBranchAddress("mParticles.fPdgCode", mParticles_fPdgCode, &b_mParticles_fPdgCode);
   fChain->SetBranchAddress("mParticles.fStatusCode", mParticles_fStatusCode, &b_mParticles_fStatusCode);
   fChain->SetBranchAddress("mParticles.fMother[2]", mParticles_fMother, &b_mParticles_fMother);
   fChain->SetBranchAddress("mParticles.fDaughter[2]", mParticles_fDaughter, &b_mParticles_fDaughter);
   fChain->SetBranchAddress("mParticles.fWeight", mParticles_fWeight, &b_mParticles_fWeight);
   fChain->SetBranchAddress("mParticles.fCalcMass", mParticles_fCalcMass, &b_mParticles_fCalcMass);
   fChain->SetBranchAddress("mParticles.fPx", mParticles_fPx, &b_mParticles_fPx);
   fChain->SetBranchAddress("mParticles.fPy", mParticles_fPy, &b_mParticles_fPy);
   fChain->SetBranchAddress("mParticles.fPz", mParticles_fPz, &b_mParticles_fPz);
   fChain->SetBranchAddress("mParticles.fE", mParticles_fE, &b_mParticles_fE);
   fChain->SetBranchAddress("mParticles.fVx", mParticles_fVx, &b_mParticles_fVx);
   fChain->SetBranchAddress("mParticles.fVy", mParticles_fVy, &b_mParticles_fVy);
   fChain->SetBranchAddress("mParticles.fVz", mParticles_fVz, &b_mParticles_fVz);
   fChain->SetBranchAddress("mParticles.fVt", mParticles_fVt, &b_mParticles_fVt);
   fChain->SetBranchAddress("mParticles.fPolarTheta", mParticles_fPolarTheta, &b_mParticles_fPolarTheta);
   fChain->SetBranchAddress("mParticles.fPolarPhi", mParticles_fPolarPhi, &b_mParticles_fPolarPhi);
   Notify();
}

Bool_t StGJetTreeAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  cout<<" i am inside ...Notify"<<endl;
   return kTRUE;
}

void StGJetTreeAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t StGJetTreeAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


//____________

//___________________________
Double_t StGJetTreeAnalysis::Get_TSP( Float_t E, Float_t en4, Float_t enp4, Float_t en3, Float_t enp3, Float_t en5, Float_t enp5, Float_t en2,Float_t enp2, Float_t en6, Float_t enp6, Float_t en1, Float_t enp1, Float_t en7, Float_t enp7, Float_t en0, Float_t enp0, Float_t en8, Float_t enp8, Float_t en01, Float_t enp01, Float_t en9, Float_t enp9, Float_t en02, Float_t enp02, Float_t en10, Float_t enp10, Float_t en03, Float_t enp03, Float_t en11, Float_t enp11)
{
  Double_t tsp =(E/((en4+enp4)*0.783+(en3+enp3+en5+enp5)*2.21+(en2+enp2+en6+enp6)*6.26+(en1+enp1+en7+enp7)*11.51+(en0+enp0+en8+enp8)*17.73+(en01+enp01+en9+enp9)*24.78+(en02+enp02+en10+enp10)*32.57+(en03+enp03+en11+enp11)*41.05));


  cout<<E<<"  en4= "<<en4  <<" enp4= "<<enp4
        <<"  en3= "<<en4  <<" enp4= "<<enp4
      <<"  en5= "<<en4  <<" enp4= "<<enp4
      <<"  en2= "<<en4  <<" enp4= "<<enp4
      <<"  en6= "<<en4  <<" enp4= "<<enp4
      <<"  en7= "<<en4  <<" enp4= "<<enp4
      <<"  en8= "<<en4  <<" enp4= "<<enp4
      <<"  en9= "<<en4  <<" enp4= "<<enp4
      <<"  en10= "<<en4  <<" enp4= "<<enp4
      <<"  en11= "<<en4  <<" enp4= "<<enp4
      <<"  en0= "<<en4  <<" enp4= "<<enp4
      <<"  en1= "<<en4  <<" enp4= "<<enp4
      <<"  en03= "<<en03  <<" enp03= "<<enp03
      <<"  en02= "<<en02  <<" enp02= "<<enp02
      <<"  en01= "<<en01  <<" enp01= "<<enp01
      <<endl;
  return tsp;
}




#endif // #ifdef StGJetTreeAnalysis_cxx
