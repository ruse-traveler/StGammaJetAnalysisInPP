//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  5 17:17:19 2016 by ROOT version 5.34/09
// from TTree GTree/GTree
// found on file: /star/data01/pwg/nihar/JetReco_pico/ppFullJet/output/sched_0DB13B7E2A6F3A1A753D077104B1B376_2.root
//////////////////////////////////////////////////////////

#ifndef StGJetTreeReader_h
#define StGJetTreeReader_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "GammaJetEvent.h"
// Header file for the classes stored in the TTree if any.
#include <TObject.h>


class GammaJetEvent;
class GammaJetTrack;
class GammaJetTower;
class GammaJetTowerUtil;
class TClonesArray;
class  TLorentzVector;



using namespace std;

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxPrimaryTrackArray = 100000;
const Int_t kMaxTowerArray = 10000000;
const Int_t kMaxTowerArray_fMatchedTracksArray = 10000;

class StGJetTreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //GammaJetEvent   *EventList;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Long64_t        runNumber;
   Long64_t        eventNumber;
   Int_t           trigID;
   Int_t           nGlobalTracks;
   Int_t           nPrimaryTracks;
   Int_t           refMult;
   Double_t        vpdVz;
   Double_t        xVertex;
   Double_t        yVertex;
   Double_t        zVertex;
   Double_t        bbcZVertex;
   Double_t        zdcCoincidenceRate;
   Double_t        bbcCoincidenceRate;
   Double_t        backgroundRate;
   Double_t        bbcBlueBackgroundRate;
   Double_t        bbcYellowBackgroundRate;
   Double_t        refMultPos;
   Double_t        refMultNeg;
   Double_t        bTOFTrayMultiplicity;
   Int_t           nVerticies;
   Double_t        MagF;
   Double_t        VrtxRank;
   Float_t         Etsp;
   Int_t           ETwrdidT;
   Int_t           ETwradc11;
   Float_t         ETwreneT0;
   Float_t         ETwreT;
   Float_t         ETwrENET0;
   Float_t         ETwrphT;
   Float_t         ETwrPTower;
   Float_t         ETwrpidTower;
   Int_t           ETwrmoduleT;
   Float_t         EClustEneT0;
   Float_t         EClustetav1;
   Float_t         EClustphiv1;
   Float_t         EEstrpen01;
   Float_t         EEstrpen02;
   Float_t         EEstrpen03;
   Float_t         EEstrpen0;
   Float_t         EEstrpen1;
   Float_t         EEstrpen2;
   Float_t         EEstrpen3;
   Float_t         EEstrpen4;
   Float_t         EEstrpen5;
   Float_t         EEstrpen6;
   Float_t         EEstrpen7;
   Float_t         EEstrpen8;
   Float_t         EEstrpen9;
   Float_t         EEstrpen10;
   Float_t         EEstrpen11;
   Float_t         EEstrpen12;
   Float_t         EEstrpen13;
   Float_t         EEstrpen14;
   Float_t         EEstrpen15;
   Int_t           ETwrdidE;
   Float_t         EPstripenp01;
   Float_t         EPstripenp02;
   Float_t         EPstripenp03;
   Float_t         EPstripenp0;
   Float_t         EPstripenp1;
   Float_t         EPstripenp2;
   Float_t         EPstripenp3;
   Float_t         EPstripenp4;
   Float_t         EPstripenp5;
   Float_t         EPstripenp6;
   Float_t         EPstripenp7;
   Float_t         EPstripenp8;
   Float_t         EPstripenp9;
   Float_t         EPstripenp10;
   Float_t         EPstripenp11;
   Float_t         EPstripenp12;
   Float_t         EPstripenp13;
   Float_t         EPstripenp14;
   Float_t         EPstripenp15;
   Float_t         EclustEnnq1;
   Float_t         EclustEnnq20;
   Float_t         EclustEnnq19;
   Float_t         EclustEnpq1;
   Float_t         EclustEnpq20;
   Float_t         EclustEnpq19;
   Float_t         EclustEnpq21;
   Int_t           PrimaryTrackArray_;
   UInt_t          PrimaryTrackArray_fUniqueID[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   UInt_t          PrimaryTrackArray_fBits[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Int_t           PrimaryTrackArray_nHitsFit[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Int_t           PrimaryTrackArray_nHitsPoss[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Int_t           PrimaryTrackArray_trackFlag[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_pZ[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_pX[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_pY[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_pT[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_dEdx[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_charge[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_tofBeta[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_eta[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_phi[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_nSigElectron[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_nSigPion[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_nSigKaon[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_nSigProton[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_dcag[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_nHits[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_dEdxHits[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_firstZPoint[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_lastZPoint[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_tofSigElectron[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_tofSigPion[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_tofSigKaon[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_tofSigProton[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_timeOfflight[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Double_t        PrimaryTrackArray_pathLength[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Int_t           PrimaryTrackArray_trkIndex[kMaxPrimaryTrackArray];   //[PrimaryTrackArray_]
   Int_t           TowerArray_;
   UInt_t          TowerArray_fUniqueID[kMaxTowerArray];   //[TowerArray_]
   UInt_t          TowerArray_fBits[kMaxTowerArray];   //[TowerArray_]
   Int_t           TowerArray_TwrId[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrEng[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrEta[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrPhi[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrADC[kMaxTowerArray];   //[TowerArray_]
   Int_t           TowerArray_TwrMatchIdnex[kMaxTowerArray];   //[TowerArray_]
   Int_t           TowerArray_NoOfmatchedTrk[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrMatchDphi[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrMatchDEta[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrMatchP[kMaxTowerArray];   //[TowerArray_]
   Int_t           TowerArray_TwrMatchIndexOfTrack[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrPx[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrPy[kMaxTowerArray];   //[TowerArray_]
   Float_t         TowerArray_TwrPz[kMaxTowerArray];   //[TowerArray_]
   Int_t           TowerArray_fNAssocTracks[kMaxTowerArray];   //[TowerArray_]
   TArrayI         TowerArray_fMatchedTracks[kMaxTowerArray];
   Int_t           TowerArray_fMatchedTracksArray_[kMaxTowerArray][10];   //[TowerArray_]
   Float_t         TowerArray_fMatchedTracksArray_P[kMaxTowerArray][10];   //[TowerArray_]

   // List of branches
   TBranch        *b_EventList_fUniqueID;   //!
   TBranch        *b_EventList_fBits;   //!
   TBranch        *b_EventList_runNumber;   //!
   TBranch        *b_EventList_eventNumber;   //!
   TBranch        *b_EventList_trigID;   //!
   TBranch        *b_EventList_nGlobalTracks;   //!
   TBranch        *b_EventList_nPrimaryTracks;   //!
   TBranch        *b_EventList_refMult;   //!
   TBranch        *b_EventList_vpdVz;   //!
   TBranch        *b_EventList_xVertex;   //!
   TBranch        *b_EventList_yVertex;   //!
   TBranch        *b_EventList_zVertex;   //!
   TBranch        *b_EventList_bbcZVertex;   //!
   TBranch        *b_EventList_zdcCoincidenceRate;   //!
   TBranch        *b_EventList_bbcCoincidenceRate;   //!
   TBranch        *b_EventList_backgroundRate;   //!
   TBranch        *b_EventList_bbcBlueBackgroundRate;   //!
   TBranch        *b_EventList_bbcYellowBackgroundRate;   //!
   TBranch        *b_EventList_refMultPos;   //!
   TBranch        *b_EventList_refMultNeg;   //!
   TBranch        *b_EventList_bTOFTrayMultiplicity;   //!
   TBranch        *b_EventList_nVerticies;   //!
   TBranch        *b_EventList_MagF;   //!
   TBranch        *b_EventList_VrtxRank;   //!
   TBranch        *b_EventList_Etsp;   //!
   TBranch        *b_EventList_ETwrdidT;   //!
   TBranch        *b_EventList_ETwradc11;   //!
   TBranch        *b_EventList_ETwreneT0;   //!
   TBranch        *b_EventList_ETwreT;   //!
   TBranch        *b_EventList_ETwrENET0;   //!
   TBranch        *b_EventList_ETwrphT;   //!
   TBranch        *b_EventList_ETwrPTower;   //!
   TBranch        *b_EventList_ETwrpidTower;   //!
   TBranch        *b_EventList_ETwrmoduleT;   //!
   TBranch        *b_EventList_EClustEneT0;   //!
   TBranch        *b_EventList_EClustetav1;   //!
   TBranch        *b_EventList_EClustphiv1;   //!
   TBranch        *b_EventList_EEstrpen01;   //!
   TBranch        *b_EventList_EEstrpen02;   //!
   TBranch        *b_EventList_EEstrpen03;   //!
   TBranch        *b_EventList_EEstrpen0;   //!
   TBranch        *b_EventList_EEstrpen1;   //!
   TBranch        *b_EventList_EEstrpen2;   //!
   TBranch        *b_EventList_EEstrpen3;   //!
   TBranch        *b_EventList_EEstrpen4;   //!
   TBranch        *b_EventList_EEstrpen5;   //!
   TBranch        *b_EventList_EEstrpen6;   //!
   TBranch        *b_EventList_EEstrpen7;   //!
   TBranch        *b_EventList_EEstrpen8;   //!
   TBranch        *b_EventList_EEstrpen9;   //!
   TBranch        *b_EventList_EEstrpen10;   //!
   TBranch        *b_EventList_EEstrpen11;   //!
   TBranch        *b_EventList_EEstrpen12;   //!
   TBranch        *b_EventList_EEstrpen13;   //!
   TBranch        *b_EventList_EEstrpen14;   //!
   TBranch        *b_EventList_EEstrpen15;   //!
   TBranch        *b_EventList_ETwrdidE;   //!
   TBranch        *b_EventList_EPstripenp01;   //!
   TBranch        *b_EventList_EPstripenp02;   //!
   TBranch        *b_EventList_EPstripenp03;   //!
   TBranch        *b_EventList_EPstripenp0;   //!
   TBranch        *b_EventList_EPstripenp1;   //!
   TBranch        *b_EventList_EPstripenp2;   //!
   TBranch        *b_EventList_EPstripenp3;   //!
   TBranch        *b_EventList_EPstripenp4;   //!
   TBranch        *b_EventList_EPstripenp5;   //!
   TBranch        *b_EventList_EPstripenp6;   //!
   TBranch        *b_EventList_EPstripenp7;   //!
   TBranch        *b_EventList_EPstripenp8;   //!
   TBranch        *b_EventList_EPstripenp9;   //!
   TBranch        *b_EventList_EPstripenp10;   //!
   TBranch        *b_EventList_EPstripenp11;   //!
   TBranch        *b_EventList_EPstripenp12;   //!
   TBranch        *b_EventList_EPstripenp13;   //!
   TBranch        *b_EventList_EPstripenp14;   //!
   TBranch        *b_EventList_EPstripenp15;   //!
   TBranch        *b_EventList_EclustEnnq1;   //!
   TBranch        *b_EventList_EclustEnnq20;   //!
   TBranch        *b_EventList_EclustEnnq19;   //!
   TBranch        *b_EventList_EclustEnpq1;   //!
   TBranch        *b_EventList_EclustEnpq20;   //!
   TBranch        *b_EventList_EclustEnpq19;   //!
   TBranch        *b_EventList_EclustEnpq21;   //!
   TBranch        *b_EventList_PrimaryTrackArray_;   //!
   TBranch        *b_PrimaryTrackArray_fUniqueID;   //!
   TBranch        *b_PrimaryTrackArray_fBits;   //!
   TBranch        *b_PrimaryTrackArray_nHitsFit;   //!
   TBranch        *b_PrimaryTrackArray_nHitsPoss;   //!
   TBranch        *b_PrimaryTrackArray_trackFlag;   //!
   TBranch        *b_PrimaryTrackArray_pZ;   //!
   TBranch        *b_PrimaryTrackArray_pX;   //!
   TBranch        *b_PrimaryTrackArray_pY;   //!
   TBranch        *b_PrimaryTrackArray_pT;   //!
   TBranch        *b_PrimaryTrackArray_dEdx;   //!
   TBranch        *b_PrimaryTrackArray_charge;   //!
   TBranch        *b_PrimaryTrackArray_tofBeta;   //!
   TBranch        *b_PrimaryTrackArray_eta;   //!
   TBranch        *b_PrimaryTrackArray_phi;   //!
   TBranch        *b_PrimaryTrackArray_nSigElectron;   //!
   TBranch        *b_PrimaryTrackArray_nSigPion;   //!
   TBranch        *b_PrimaryTrackArray_nSigKaon;   //!
   TBranch        *b_PrimaryTrackArray_nSigProton;   //!
   TBranch        *b_PrimaryTrackArray_dcag;   //!
   TBranch        *b_PrimaryTrackArray_nHits;   //!
   TBranch        *b_PrimaryTrackArray_dEdxHits;   //!
   TBranch        *b_PrimaryTrackArray_firstZPoint;   //!
   TBranch        *b_PrimaryTrackArray_lastZPoint;   //!
   TBranch        *b_PrimaryTrackArray_tofSigElectron;   //!
   TBranch        *b_PrimaryTrackArray_tofSigPion;   //!
   TBranch        *b_PrimaryTrackArray_tofSigKaon;   //!
   TBranch        *b_PrimaryTrackArray_tofSigProton;   //!
   TBranch        *b_PrimaryTrackArray_timeOfflight;   //!
   TBranch        *b_PrimaryTrackArray_pathLength;   //!
   TBranch        *b_PrimaryTrackArray_trkIndex;   //!
   TBranch        *b_EventList_TowerArray_;   //!
   TBranch        *b_TowerArray_fUniqueID;   //!
   TBranch        *b_TowerArray_fBits;   //!
   TBranch        *b_TowerArray_TwrId;   //!
   TBranch        *b_TowerArray_TwrEng;   //!
   TBranch        *b_TowerArray_TwrEta;   //!
   TBranch        *b_TowerArray_TwrPhi;   //!
   TBranch        *b_TowerArray_TwrADC;   //!
   TBranch        *b_TowerArray_TwrMatchIdnex;   //!
   TBranch        *b_TowerArray_NoOfmatchedTrk;   //!
   TBranch        *b_TowerArray_TwrMatchDphi;   //!
   TBranch        *b_TowerArray_TwrMatchDEta;   //!
   TBranch        *b_TowerArray_TwrMatchP;   //!
   TBranch        *b_TowerArray_TwrMatchIndexOfTrack;   //!
   TBranch        *b_TowerArray_TwrPx;   //!
   TBranch        *b_TowerArray_TwrPy;   //!
   TBranch        *b_TowerArray_TwrPz;   //!
   TBranch        *b_TowerArray_fNAssocTracks;   //!
   TBranch        *b_TowerArray_fMatchedTracks;   //!
   TBranch        *b_TowerArray_fMatchedTracksArray_;   //!
   TBranch        *b_TowerArray_fMatchedTracksArray_P;   //!



   //_____________________________
   ///////////////////////

   //   TH1F* mHistBg;
   //   TH1F* hya;
   int        mEventCounter;
   TFile* File_output;
   const char* outFile;


   ///________

   TBranch *eventBranch;

   Float_t *EvValues;
   Float_t *TrValues;
   //  Float_t *TwrValues;

   TTree *outTree;
   GammaJetEvent *event;


   const char* fileBaseName;
   const char* mFileIndex;
   const char* mOutDir;
   //___________________________


   StGJetTreeReader(TTree *tree=0);
   virtual ~StGJetTreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Init_Picotree(char* outFile);
   virtual void     Make();
   virtual void     Finish();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   ClassDef( StGJetTreeReader,1)
};

#endif

#ifdef StGJetTreeReader_cxx
//________________________________________________
//void StGJetTreeReader::Init_Picotree(TTree *outTree)
void StGJetTreeReader::Init_Picotree(char* outFile)
{

  File_output= new TFile(outFile,"RECREATE");
  event = new GammaJetEvent();
  outTree = new TTree("GTree","GTree");
  outTree->SetAutoSave(-500000000);  // autosave activated for each 5 MB
  outTree->SetMaxTreeSize(1000000000);  // new file with name _n.root will be created after 1GB size
  eventBranch = outTree->Branch("EventList",&event,1000000);
  EvValues = new float[81];
  TrValues = new float[28];
  


}

//______________________________________
void StGJetTreeReader::Finish()
{
  File_output = outTree->GetCurrentFile();
  File_output->Write();
  File_output->Close();

}
//_____________________________________________________________________
StGJetTreeReader::StGJetTreeReader(TTree *tree) : fChain(0) 
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

     TChain *chain = new TChain("Gfmtodst"," ");



     chain->Add("/star/data01/pwg/nihar/JetReco_pico/ppFullJet/pp_test_1.root");
     chain->Add("/star/data01/pwg/nihar/JetReco_pico/ppFullJet/pp200_test_2.root");
     chain->Add("/star/data01/pwg/nihar/JetReco_pico/ppFullJet/pp200_test_3.root");

     tree=chain;


     
   }
   Init(tree);

   //   Init_Picotree("pp_200_test.root");
}

StGJetTreeReader::~StGJetTreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t StGJetTreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t StGJetTreeReader::LoadTree(Long64_t entry)
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

void StGJetTreeReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_EventList_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_EventList_fBits);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_EventList_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_EventList_eventNumber);
   fChain->SetBranchAddress("trigID", &trigID, &b_EventList_trigID);
   fChain->SetBranchAddress("nGlobalTracks", &nGlobalTracks, &b_EventList_nGlobalTracks);
   fChain->SetBranchAddress("nPrimaryTracks", &nPrimaryTracks, &b_EventList_nPrimaryTracks);
   fChain->SetBranchAddress("refMult", &refMult, &b_EventList_refMult);
   fChain->SetBranchAddress("vpdVz", &vpdVz, &b_EventList_vpdVz);
   fChain->SetBranchAddress("xVertex", &xVertex, &b_EventList_xVertex);
   fChain->SetBranchAddress("yVertex", &yVertex, &b_EventList_yVertex);
   fChain->SetBranchAddress("zVertex", &zVertex, &b_EventList_zVertex);
   fChain->SetBranchAddress("bbcZVertex", &bbcZVertex, &b_EventList_bbcZVertex);
   fChain->SetBranchAddress("zdcCoincidenceRate", &zdcCoincidenceRate, &b_EventList_zdcCoincidenceRate);
   fChain->SetBranchAddress("bbcCoincidenceRate", &bbcCoincidenceRate, &b_EventList_bbcCoincidenceRate);
   fChain->SetBranchAddress("backgroundRate", &backgroundRate, &b_EventList_backgroundRate);
   fChain->SetBranchAddress("bbcBlueBackgroundRate", &bbcBlueBackgroundRate, &b_EventList_bbcBlueBackgroundRate);
   fChain->SetBranchAddress("bbcYellowBackgroundRate", &bbcYellowBackgroundRate, &b_EventList_bbcYellowBackgroundRate);
   fChain->SetBranchAddress("refMultPos", &refMultPos, &b_EventList_refMultPos);
   fChain->SetBranchAddress("refMultNeg", &refMultNeg, &b_EventList_refMultNeg);
   fChain->SetBranchAddress("bTOFTrayMultiplicity", &bTOFTrayMultiplicity, &b_EventList_bTOFTrayMultiplicity);
   fChain->SetBranchAddress("nVerticies", &nVerticies, &b_EventList_nVerticies);
   fChain->SetBranchAddress("MagF", &MagF, &b_EventList_MagF);
   fChain->SetBranchAddress("VrtxRank", &VrtxRank, &b_EventList_VrtxRank);
   fChain->SetBranchAddress("Etsp", &Etsp, &b_EventList_Etsp);
   fChain->SetBranchAddress("ETwrdidT", &ETwrdidT, &b_EventList_ETwrdidT);
   fChain->SetBranchAddress("ETwradc11", &ETwradc11, &b_EventList_ETwradc11);
   fChain->SetBranchAddress("ETwreneT0", &ETwreneT0, &b_EventList_ETwreneT0);
   fChain->SetBranchAddress("ETwreT", &ETwreT, &b_EventList_ETwreT);
   fChain->SetBranchAddress("ETwrENET0", &ETwrENET0, &b_EventList_ETwrENET0);
   fChain->SetBranchAddress("ETwrphT", &ETwrphT, &b_EventList_ETwrphT);
   fChain->SetBranchAddress("ETwrPTower", &ETwrPTower, &b_EventList_ETwrPTower);
   fChain->SetBranchAddress("ETwrpidTower", &ETwrpidTower, &b_EventList_ETwrpidTower);
   fChain->SetBranchAddress("ETwrmoduleT", &ETwrmoduleT, &b_EventList_ETwrmoduleT);
   fChain->SetBranchAddress("EClustEneT0", &EClustEneT0, &b_EventList_EClustEneT0);
   fChain->SetBranchAddress("EClustetav1", &EClustetav1, &b_EventList_EClustetav1);
   fChain->SetBranchAddress("EClustphiv1", &EClustphiv1, &b_EventList_EClustphiv1);
   fChain->SetBranchAddress("EEstrpen01", &EEstrpen01, &b_EventList_EEstrpen01);
   fChain->SetBranchAddress("EEstrpen02", &EEstrpen02, &b_EventList_EEstrpen02);
   fChain->SetBranchAddress("EEstrpen03", &EEstrpen03, &b_EventList_EEstrpen03);
   fChain->SetBranchAddress("EEstrpen0", &EEstrpen0, &b_EventList_EEstrpen0);
   fChain->SetBranchAddress("EEstrpen1", &EEstrpen1, &b_EventList_EEstrpen1);
   fChain->SetBranchAddress("EEstrpen2", &EEstrpen2, &b_EventList_EEstrpen2);
   fChain->SetBranchAddress("EEstrpen3", &EEstrpen3, &b_EventList_EEstrpen3);
   fChain->SetBranchAddress("EEstrpen4", &EEstrpen4, &b_EventList_EEstrpen4);
   fChain->SetBranchAddress("EEstrpen5", &EEstrpen5, &b_EventList_EEstrpen5);
   fChain->SetBranchAddress("EEstrpen6", &EEstrpen6, &b_EventList_EEstrpen6);
   fChain->SetBranchAddress("EEstrpen7", &EEstrpen7, &b_EventList_EEstrpen7);
   fChain->SetBranchAddress("EEstrpen8", &EEstrpen8, &b_EventList_EEstrpen8);
   fChain->SetBranchAddress("EEstrpen9", &EEstrpen9, &b_EventList_EEstrpen9);
   fChain->SetBranchAddress("EEstrpen10", &EEstrpen10, &b_EventList_EEstrpen10);
   fChain->SetBranchAddress("EEstrpen11", &EEstrpen11, &b_EventList_EEstrpen11);
   fChain->SetBranchAddress("EEstrpen12", &EEstrpen12, &b_EventList_EEstrpen12);
   fChain->SetBranchAddress("EEstrpen13", &EEstrpen13, &b_EventList_EEstrpen13);
   fChain->SetBranchAddress("EEstrpen14", &EEstrpen14, &b_EventList_EEstrpen14);
   fChain->SetBranchAddress("EEstrpen15", &EEstrpen15, &b_EventList_EEstrpen15);
   fChain->SetBranchAddress("ETwrdidE", &ETwrdidE, &b_EventList_ETwrdidE);
   fChain->SetBranchAddress("EPstripenp01", &EPstripenp01, &b_EventList_EPstripenp01);
   fChain->SetBranchAddress("EPstripenp02", &EPstripenp02, &b_EventList_EPstripenp02);
   fChain->SetBranchAddress("EPstripenp03", &EPstripenp03, &b_EventList_EPstripenp03);
   fChain->SetBranchAddress("EPstripenp0", &EPstripenp0, &b_EventList_EPstripenp0);
   fChain->SetBranchAddress("EPstripenp1", &EPstripenp1, &b_EventList_EPstripenp1);
   fChain->SetBranchAddress("EPstripenp2", &EPstripenp2, &b_EventList_EPstripenp2);
   fChain->SetBranchAddress("EPstripenp3", &EPstripenp3, &b_EventList_EPstripenp3);
   fChain->SetBranchAddress("EPstripenp4", &EPstripenp4, &b_EventList_EPstripenp4);
   fChain->SetBranchAddress("EPstripenp5", &EPstripenp5, &b_EventList_EPstripenp5);
   fChain->SetBranchAddress("EPstripenp6", &EPstripenp6, &b_EventList_EPstripenp6);
   fChain->SetBranchAddress("EPstripenp7", &EPstripenp7, &b_EventList_EPstripenp7);
   fChain->SetBranchAddress("EPstripenp8", &EPstripenp8, &b_EventList_EPstripenp8);
   fChain->SetBranchAddress("EPstripenp9", &EPstripenp9, &b_EventList_EPstripenp9);
   fChain->SetBranchAddress("EPstripenp10", &EPstripenp10, &b_EventList_EPstripenp10);
   fChain->SetBranchAddress("EPstripenp11", &EPstripenp11, &b_EventList_EPstripenp11);
   fChain->SetBranchAddress("EPstripenp12", &EPstripenp12, &b_EventList_EPstripenp12);
   fChain->SetBranchAddress("EPstripenp13", &EPstripenp13, &b_EventList_EPstripenp13);
   fChain->SetBranchAddress("EPstripenp14", &EPstripenp14, &b_EventList_EPstripenp14);
   fChain->SetBranchAddress("EPstripenp15", &EPstripenp15, &b_EventList_EPstripenp15);
   fChain->SetBranchAddress("EclustEnnq1", &EclustEnnq1, &b_EventList_EclustEnnq1);
   fChain->SetBranchAddress("EclustEnnq20", &EclustEnnq20, &b_EventList_EclustEnnq20);
   fChain->SetBranchAddress("EclustEnnq19", &EclustEnnq19, &b_EventList_EclustEnnq19);
   fChain->SetBranchAddress("EclustEnpq1", &EclustEnpq1, &b_EventList_EclustEnpq1);
   fChain->SetBranchAddress("EclustEnpq20", &EclustEnpq20, &b_EventList_EclustEnpq20);
   fChain->SetBranchAddress("EclustEnpq19", &EclustEnpq19, &b_EventList_EclustEnpq19);
   fChain->SetBranchAddress("EclustEnpq21", &EclustEnpq21, &b_EventList_EclustEnpq21);
   fChain->SetBranchAddress("PrimaryTrackArray", &PrimaryTrackArray_, &b_EventList_PrimaryTrackArray_);
   fChain->SetBranchAddress("PrimaryTrackArray.fUniqueID", PrimaryTrackArray_fUniqueID, &b_PrimaryTrackArray_fUniqueID);
   fChain->SetBranchAddress("PrimaryTrackArray.fBits", PrimaryTrackArray_fBits, &b_PrimaryTrackArray_fBits);
   fChain->SetBranchAddress("PrimaryTrackArray.nHitsFit", PrimaryTrackArray_nHitsFit, &b_PrimaryTrackArray_nHitsFit);
   fChain->SetBranchAddress("PrimaryTrackArray.nHitsPoss", PrimaryTrackArray_nHitsPoss, &b_PrimaryTrackArray_nHitsPoss);
   fChain->SetBranchAddress("PrimaryTrackArray.trackFlag", PrimaryTrackArray_trackFlag, &b_PrimaryTrackArray_trackFlag);
   fChain->SetBranchAddress("PrimaryTrackArray.pZ", PrimaryTrackArray_pZ, &b_PrimaryTrackArray_pZ);
   fChain->SetBranchAddress("PrimaryTrackArray.pX", PrimaryTrackArray_pX, &b_PrimaryTrackArray_pX);
   fChain->SetBranchAddress("PrimaryTrackArray.pY", PrimaryTrackArray_pY, &b_PrimaryTrackArray_pY);
   fChain->SetBranchAddress("PrimaryTrackArray.pT", PrimaryTrackArray_pT, &b_PrimaryTrackArray_pT);
   fChain->SetBranchAddress("PrimaryTrackArray.dEdx", PrimaryTrackArray_dEdx, &b_PrimaryTrackArray_dEdx);
   fChain->SetBranchAddress("PrimaryTrackArray.charge", PrimaryTrackArray_charge, &b_PrimaryTrackArray_charge);
   fChain->SetBranchAddress("PrimaryTrackArray.tofBeta", PrimaryTrackArray_tofBeta, &b_PrimaryTrackArray_tofBeta);
   fChain->SetBranchAddress("PrimaryTrackArray.eta", PrimaryTrackArray_eta, &b_PrimaryTrackArray_eta);
   fChain->SetBranchAddress("PrimaryTrackArray.phi", PrimaryTrackArray_phi, &b_PrimaryTrackArray_phi);
   fChain->SetBranchAddress("PrimaryTrackArray.nSigElectron", PrimaryTrackArray_nSigElectron, &b_PrimaryTrackArray_nSigElectron);
   fChain->SetBranchAddress("PrimaryTrackArray.nSigPion", PrimaryTrackArray_nSigPion, &b_PrimaryTrackArray_nSigPion);
   fChain->SetBranchAddress("PrimaryTrackArray.nSigKaon", PrimaryTrackArray_nSigKaon, &b_PrimaryTrackArray_nSigKaon);
   fChain->SetBranchAddress("PrimaryTrackArray.nSigProton", PrimaryTrackArray_nSigProton, &b_PrimaryTrackArray_nSigProton);
   fChain->SetBranchAddress("PrimaryTrackArray.dcag", PrimaryTrackArray_dcag, &b_PrimaryTrackArray_dcag);
   fChain->SetBranchAddress("PrimaryTrackArray.nHits", PrimaryTrackArray_nHits, &b_PrimaryTrackArray_nHits);
   fChain->SetBranchAddress("PrimaryTrackArray.dEdxHits", PrimaryTrackArray_dEdxHits, &b_PrimaryTrackArray_dEdxHits);
   fChain->SetBranchAddress("PrimaryTrackArray.firstZPoint", PrimaryTrackArray_firstZPoint, &b_PrimaryTrackArray_firstZPoint);
   fChain->SetBranchAddress("PrimaryTrackArray.lastZPoint", PrimaryTrackArray_lastZPoint, &b_PrimaryTrackArray_lastZPoint);
   fChain->SetBranchAddress("PrimaryTrackArray.tofSigElectron", PrimaryTrackArray_tofSigElectron, &b_PrimaryTrackArray_tofSigElectron);
   fChain->SetBranchAddress("PrimaryTrackArray.tofSigPion", PrimaryTrackArray_tofSigPion, &b_PrimaryTrackArray_tofSigPion);
   fChain->SetBranchAddress("PrimaryTrackArray.tofSigKaon", PrimaryTrackArray_tofSigKaon, &b_PrimaryTrackArray_tofSigKaon);
   fChain->SetBranchAddress("PrimaryTrackArray.tofSigProton", PrimaryTrackArray_tofSigProton, &b_PrimaryTrackArray_tofSigProton);
   fChain->SetBranchAddress("PrimaryTrackArray.timeOfflight", PrimaryTrackArray_timeOfflight, &b_PrimaryTrackArray_timeOfflight);
   fChain->SetBranchAddress("PrimaryTrackArray.pathLength", PrimaryTrackArray_pathLength, &b_PrimaryTrackArray_pathLength);
   fChain->SetBranchAddress("PrimaryTrackArray.trkIndex", PrimaryTrackArray_trkIndex, &b_PrimaryTrackArray_trkIndex);
   fChain->SetBranchAddress("TowerArray", &TowerArray_, &b_EventList_TowerArray_);
   fChain->SetBranchAddress("TowerArray.fUniqueID", TowerArray_fUniqueID, &b_TowerArray_fUniqueID);
   fChain->SetBranchAddress("TowerArray.fBits", TowerArray_fBits, &b_TowerArray_fBits);
   fChain->SetBranchAddress("TowerArray.TwrId", TowerArray_TwrId, &b_TowerArray_TwrId);
   fChain->SetBranchAddress("TowerArray.TwrEng", TowerArray_TwrEng, &b_TowerArray_TwrEng);
   fChain->SetBranchAddress("TowerArray.TwrEta", TowerArray_TwrEta, &b_TowerArray_TwrEta);
   fChain->SetBranchAddress("TowerArray.TwrPhi", TowerArray_TwrPhi, &b_TowerArray_TwrPhi);
   fChain->SetBranchAddress("TowerArray.TwrADC", TowerArray_TwrADC, &b_TowerArray_TwrADC);
   fChain->SetBranchAddress("TowerArray.TwrMatchIdnex", TowerArray_TwrMatchIdnex, &b_TowerArray_TwrMatchIdnex);
   fChain->SetBranchAddress("TowerArray.NoOfmatchedTrk", TowerArray_NoOfmatchedTrk, &b_TowerArray_NoOfmatchedTrk);
   fChain->SetBranchAddress("TowerArray.TwrMatchDphi", TowerArray_TwrMatchDphi, &b_TowerArray_TwrMatchDphi);
   fChain->SetBranchAddress("TowerArray.TwrMatchDEta", TowerArray_TwrMatchDEta, &b_TowerArray_TwrMatchDEta);
   fChain->SetBranchAddress("TowerArray.TwrMatchP", TowerArray_TwrMatchP, &b_TowerArray_TwrMatchP);
   fChain->SetBranchAddress("TowerArray.TwrMatchIndexOfTrack", TowerArray_TwrMatchIndexOfTrack, &b_TowerArray_TwrMatchIndexOfTrack);
   fChain->SetBranchAddress("TowerArray.TwrPx", TowerArray_TwrPx, &b_TowerArray_TwrPx);
   fChain->SetBranchAddress("TowerArray.TwrPy", TowerArray_TwrPy, &b_TowerArray_TwrPy);
   fChain->SetBranchAddress("TowerArray.TwrPz", TowerArray_TwrPz, &b_TowerArray_TwrPz);
   fChain->SetBranchAddress("TowerArray.fNAssocTracks", TowerArray_fNAssocTracks, &b_TowerArray_fNAssocTracks);
   fChain->SetBranchAddress("TowerArray.fMatchedTracks", TowerArray_fMatchedTracks, &b_TowerArray_fMatchedTracks);
   fChain->SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", TowerArray_fMatchedTracksArray_, &b_TowerArray_fMatchedTracksArray_);
   fChain->SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", TowerArray_fMatchedTracksArray_P, &b_TowerArray_fMatchedTracksArray_P);
   Notify();
}

Bool_t StGJetTreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  cout<<" i am inside ...Notify"<<endl;
   return kTRUE;
}

void StGJetTreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t StGJetTreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef StGJetTreeReader_cxx
