// 'GammaJetAnalysis.h'
// Derek Anderson
// 01.08.2016
//
// This file contains the defintion of the Analysis class.  It is executed in
// the file 'GammaJetExecute.cc'.  Global variables (notably the various pT
// cuts) are set in 'GammaJetExecute.cc'.  Note that not all the branches in
// the TTree will be relevant for the pythis simulation; they're left in so
// that the TTree can be read in by external macros.
//
// NOTE: 'weird' tracks refer to situations where a TPC track DOESN'T deposit
// all of its energy in the calorimeter.  The energy it does deposit will be
// picked up by the jet-finder as an additional calorimeter hit.  These 'weird'
// tracks are assumed to have the mass of a pion (140 MeV).
//
// Last updated: 08.30.2016

#pragma once

#include "TString.h"
#include "GammaJetEvent.h"
#include "GammaJetTrack.h"
#include "GammaJetVertex.h"
#include "GammaJetParton.h"
#include "GammaJetSystem.h"
#include "GammaJetTrigger.h"



class Analysis {

private:

  int     nEvt;               // no. of events processed
  bool    trig;               // indicates if a trigger has been found
  bool    trig2;              // indicates if secondary trigger has been found
  float   eCut;               // trigger eT (pT) cut
  float   hCut;               // trigger eta cut
  float   trCut;              // "trigger" region endpoints
  float   reCut;              // "recoil" region endpoints
  float   ueCut;              // "underlying evt." region endpoints
  float   a1Cut;              // lower "away" region endpoint
  float   a2Cut;              // upper "away" region endpoint
  float   rIso;               // size of cone for isolation cut
  System  sys;                // misc. info for analysis
  Event   Pevt,  Devt;        // event info
  Vertex  Pvtx,  Dvtx;        // vertex info
  Parton  Pasp,  Pnsp;        // (particle level) hard scatter product info
  Parton  Dasp,  Dnsp;        // (detector level) hard scatter product info
  Trigger Ptrg,  Dtrg;        // trigger info
  vector<Track>  Ptrk, Dtrk;  // track info
  vector<Track>  Wtrk;        // 'weird' tracks
  Pythia8::Event event;       // the event

  // Root members
  TFile    *file;         // output file
  TTree    *Pvent;        // particle level tree
  TTree    *Dvent;        // detector level tree
  TH1D     *pT_trks;
  TH1D     *pT_accTrks;
  TH1D     *zTparton;
  TH1D     *zTns;
  TH1D     *zTas;
  TH2D     *pT_tpcReso;
  TH2D     *E_calReso;
  TH2D     *pT_smear;
  TProfile *pT_tpcResoP;
  TProfile *E_calResoP;
  TProfile *pT_smearP;

  // private methods
  int    sortTrack(const int iTrack, const int iTrig);
  void   fillTrees(const int evtCode);
  void   locateNearSide(const int iTrig);
  void   calculateVtx(const int index, const bool det);
  bool   addNearTrk(const int iTrack, const bool first);
  bool   addAwayTrk(const int iTrack, const bool first);
  bool   addUncorrelatedTrk(const int iTrack);
  bool   sortPartons(const int iTrig);
  double smearPt(const double pt, const int id, const TString& source);

public:

  Analysis() {}

  // public methods
  void init(const float pT, const float eta, const float r, const float dFtr, const float dFre, const float dFue, const float dFa1, const float dFa2, TFile *f);
  void analyze(const int code, Pythia8::Event& evt);
  void finish();

};



// public methods -------------------------------------------------------------


void Analysis::init(const float pT, const float eta, const float r, const float dFtr, const float dFre, const float dFue, const float dFa1, const float dFa2, TFile *f) {

  nEvt  = 0;
  eCut  = pT;
  hCut  = eta;
  trCut = dFtr;
  reCut = dFre;
  ueCut = dFue;
  a1Cut = dFa1;
  a2Cut = dFa2;
  rIso  = r;
  file  = f;

  // make sure vectors are clear
  Ptrk.clear();
  Dtrk.clear();


  // grow some trees
  Pvent = new TTree("ParTree", "A detailed tree of (particle-level) events");
  Dvent = new TTree("DetTree", "A detailed tree of (detector-level) events");

  // (particle) event branches
  Pvent -> Branch("Events_num", &Pevt.Events_num, "Events_num/I");
  Pvent -> Branch("Events_Process", &Pevt.Events_Process, "Events_Process/I");
  Pvent -> Branch("Events_refmult", &Pevt.Events_refmult, "Events_refmult/I");
  Pvent -> Branch("Events_refPos", &Pevt.Events_refPos, "Events_refPos/I");
  Pvent -> Branch("Events_refNeg", &Pevt.Events_refNeg, "Events_refNeg/I");
  Pvent -> Branch("Events_runId", &Pevt.Events_runId, "Events_runId/I");
  Pvent -> Branch("Events_MagF", &Pevt.Events_MagF, "Events_MagF/D");
  Pvent -> Branch("Events_nVertex", &Pevt.Events_nVertex, "Events_nVertex/I");
  Pvent -> Branch("Events_rankV", &Pevt.Events_rankV, "Events_rankV");
  Pvent -> Branch("Events_pvpdz", &Pevt.Events_pvpdz, "Events_pvpdz/D");
  Pvent -> Branch("Events_BBCz", &Pevt.Events_BBCz, "Events_BBCz/D");
  Pvent -> Branch("Events_ZDCz", &Pevt.Events_ZDCz, "Events_ZDCz/D");
  Pvent -> Branch("Events_primVx", &Pevt.Events_primVx, "Events_primVx/D");
  Pvent -> Branch("Events_primVy", &Pevt.Events_primVy, "Events_primVy/D");
  Pvent -> Branch("Events_primVz", &Pevt.Events_primVz, "Events_primVz/D");
  Pvent -> Branch("Events_TrigId", &Pevt.Events_TrigId, "Events_TrigId/I");
  Pvent -> Branch("Events_TrigStat", &Pevt.Events_TrigStat, "Events_TrigStat/I");
  Pvent -> Branch("Events_TrigIndex", &Pevt.Events_TrigIndex, "Events_TrigIndex/I");
  Pvent -> Branch("Events_PartonId_NS", &Pevt.Events_PartonId_NS, "Events_PartonId_NS/I");
  Pvent -> Branch("Events_PartonId_AS", &Pevt.Events_PartonId_AS, "Events_PartonId_AS/I");
  Pvent -> Branch("Events_PartonStat_NS", &Pevt.Events_PartonStat_NS, "Events_PartonStat_NS/I");
  Pvent -> Branch("Events_PartonStat_AS", &Pevt.Events_PartonStat_AS, "Events_PartonStat_AS/I");
  Pvent -> Branch("Events_PartonEta_NS", &Pevt.Events_PartonEta_NS, "Events_PartonEta_NS/F");
  Pvent -> Branch("Events_PartonEta_AS", &Pevt.Events_PartonEta_AS, "Events_PartonEta_AS/F");
  Pvent -> Branch("Events_PartonPhi_NS", &Pevt.Events_PartonPhi_NS, "Events_PartonPhi_NS/F");
  Pvent -> Branch("Events_PartonPhi_AS", &Pevt.Events_PartonPhi_AS, "Events_PartonPhi_AS/F");
  Pvent -> Branch("Events_PartonE_NS", &Pevt.Events_PartonE_NS, "Events_PartonE_NS/F");
  Pvent -> Branch("Events_PartonE_AS", &Pevt.Events_PartonE_AS, "Events_PartonE_AS/F");
  Pvent -> Branch("Events_PartonEt_NS", &Pevt.Events_PartonEt_NS, "Events_PartonEt_NS/F");
  Pvent -> Branch("Events_PartonEt_AS", &Pevt.Events_PartonEt_AS, "Events_PartonEt_AS/F");
  Pvent -> Branch("Events_tsp", &Pevt.Events_tsp, "Events_tsp/F");
  Pvent -> Branch("Events_iso", &Pevt.Events_iso, "Events_iso/F");
  // (particle) tower branches
  Pvent -> Branch("Events_Twr_didT", &Pevt.Events_Twr_didT, "Events_Twr_didT/I");
  Pvent -> Branch("Events_Twr_adc11", &Pevt.Events_Twr_adc11, "Events_Twr_adc11/I");
  Pvent -> Branch("Events_Twr_eneT0", &Pevt.Events_Twr_eneT0, "Events_Twr_eneT0/F");
  Pvent -> Branch("Events_Twr_eT", &Pevt.Events_Twr_eT, "Events_Twr_eT/F");
  Pvent -> Branch("Events_Twr_ENET0", &Pevt.Events_Twr_ENET0, "Events_Twr_ENET0/F");
  Pvent -> Branch("Events_Twr_phT", &Pevt.Events_Twr_phT, "Events_Twr_phT/F");
  Pvent -> Branch("Events_Twr_PTower", &Pevt.Events_Twr_PTower, "Events_Twr_PTower/F");
  Pvent -> Branch("Events_Twr_pidTower", &Pevt.Events_Twr_pidTower, "Events_Twr_pidTower/F");
  Pvent -> Branch("Events_Twr_moduleT", &Pevt.Events_Twr_moduleT, "Events_Twr_moduleT/F");
  // (particle) tower cluster branches
  Pvent -> Branch("Events_Clust_EneT0", &Pevt.Events_Clust_EneT0, "Events_Clust_EneT0/F");
  Pvent -> Branch("Events_Clust_EneR", &Pevt.Events_Clust_EneR, "Events_Clust_EneR/F");
  Pvent -> Branch("Events_Clust_EneH", &Pevt.Events_Clust_EneH, "Events_Clust_EneH/F");
  Pvent -> Branch("Events_Clust_rapv1", &Pevt.Events_Clust_rapv1, "Events_Clust_rapv1/F");
  Pvent -> Branch("Events_Clust_etav1", &Pevt.Events_Clust_etav1, "Events_Clust_etav1/F");
  Pvent -> Branch("Events_Clust_phiv1", &Pevt.Events_Clust_phiv1, "Events_Clust_phiv1/F");
  // (particle) SMD eta strip branches
  Pvent -> Branch("Events_Estrp_en01", &Pevt.Events_Estrp_en01, "Events_Estrp_en01/F");
  Pvent -> Branch("Events_Estrp_en02", &Pevt.Events_Estrp_en02, "Events_Estrp_en02/F");
  Pvent -> Branch("Events_Estrp_en03", &Pevt.Events_Estrp_en03, "Events_Estrp_en03/F");
  Pvent -> Branch("Events_Estrp_en0", &Pevt.Events_Estrp_en0, "Events_Estrp_en0/F");
  Pvent -> Branch("Events_Estrp_en1", &Pevt.Events_Estrp_en1, "Events_Estrp_en1/F");
  Pvent -> Branch("Events_Estrp_en2", &Pevt.Events_Estrp_en2, "Events_Estrp_en2/F");
  Pvent -> Branch("Events_Estrp_en3", &Pevt.Events_Estrp_en3, "Events_Estrp_en3/F");
  Pvent -> Branch("Events_Estrp_en4", &Pevt.Events_Estrp_en4, "Events_Estrp_en4/F");
  Pvent -> Branch("Events_Estrp_en5", &Pevt.Events_Estrp_en5, "Events_Estrp_en5/F");
  Pvent -> Branch("Events_Estrp_en6", &Pevt.Events_Estrp_en6, "Events_Estrp_en6/F");
  Pvent -> Branch("Events_Estrp_en7", &Pevt.Events_Estrp_en7, "Events_Estrp_en7/F");
  Pvent -> Branch("Events_Estrp_en8", &Pevt.Events_Estrp_en8, "Events_Estrp_en8/F");
  Pvent -> Branch("Events_Estrp_en9", &Pevt.Events_Estrp_en9, "Events_Estrp_en9/F");
  Pvent -> Branch("Events_Estrp_en10", &Pevt.Events_Estrp_en10, "Events_Estrp_en10/F");
  Pvent -> Branch("Events_Estrp_en11", &Pevt.Events_Estrp_en11, "Events_Estrp_en11/F");
  Pvent -> Branch("Events_Estrp_en12", &Pevt.Events_Estrp_en12, "Events_Estrp_en12/F");
  Pvent -> Branch("Events_Estrp_en13", &Pevt.Events_Estrp_en13, "Events_Estrp_en13/F");
  Pvent -> Branch("Events_Estrp_en14", &Pevt.Events_Estrp_en14, "Events_Estrp_en14/F");
  Pvent -> Branch("Events_Estrp_en15", &Pevt.Events_Estrp_en15, "Events_Estrp_en15/F");
  // another (particle) tower branch
  Pvent -> Branch("Events_Twr_didE", &Pevt.Events_Twr_didE, "Events_Twr_didE/I");
  // (particle) SMD phi strip branches
  Pvent -> Branch("Events_Pstrip_enp01", &Pevt.Events_Pstrip_enp01, "Events_Pstrip_enp01/F");
  Pvent -> Branch("Events_Pstrip_enp02", &Pevt.Events_Pstrip_enp02, "Events_Pstrip_enp02/F");
  Pvent -> Branch("Events_Pstrip_enp03", &Pevt.Events_Pstrip_enp03, "Events_Pstrip_enp03/F");
  Pvent -> Branch("Events_Pstrip_enp0", &Pevt.Events_Pstrip_enp0, "Events_Pstrip_enp0/F");
  Pvent -> Branch("Events_Pstrip_enp1", &Pevt.Events_Pstrip_enp1, "Events_Pstrip_enp1/F");
  Pvent -> Branch("Events_Pstrip_enp2", &Pevt.Events_Pstrip_enp2, "Events_Pstrip_enp2/F");
  Pvent -> Branch("Events_Pstrip_enp3", &Pevt.Events_Pstrip_enp3, "Events_Pstrip_enp3/F");
  Pvent -> Branch("Events_Pstrip_enp4", &Pevt.Events_Pstrip_enp4, "Events_Pstrip_enp4/F");
  Pvent -> Branch("Events_Pstrip_enp5", &Pevt.Events_Pstrip_enp5, "Events_Pstrip_enp5/F");
  Pvent -> Branch("Events_Pstrip_enp6", &Pevt.Events_Pstrip_enp6, "Events_Pstrip_enp6/F");
  Pvent -> Branch("Events_Pstrip_enp7", &Pevt.Events_Pstrip_enp7, "Events_Pstrip_enp7/F");
  Pvent -> Branch("Events_Pstrip_enp8", &Pevt.Events_Pstrip_enp8, "Events_Pstrip_enp8/F");
  Pvent -> Branch("Events_Pstrip_enp9", &Pevt.Events_Pstrip_enp9, "Events_Pstrip_enp9/F");
  Pvent -> Branch("Events_Pstrip_enp10", &Pevt.Events_Pstrip_enp10, "Events_Pstrip_enp10/F");
  Pvent -> Branch("Events_Pstrip_enp11", &Pevt.Events_Pstrip_enp11, "Events_Pstrip_enp11/F");
  Pvent -> Branch("Events_Pstrip_enp12", &Pevt.Events_Pstrip_enp12, "Events_Pstrip_enp12/F");
  Pvent -> Branch("Events_Pstrip_enp13", &Pevt.Events_Pstrip_enp13, "Events_Pstrip_enp13/F");
  Pvent -> Branch("Events_Pstrip_enp14", &Pevt.Events_Pstrip_enp14, "Events_Pstrip_enp14/F");
  Pvent -> Branch("Events_Pstrip_enp15", &Pevt.Events_Pstrip_enp15, "Events_Pstrip_enp15/F");
  // (particle) SMD cluster branches
  Pvent -> Branch("Events_clust_Ennq1", &Pevt.Events_clust_Ennq1, "Events_clust_Ennq1/F");
  Pvent -> Branch("Events_clust_Ennq20", &Pevt.Events_clust_Ennq20, "Events_clust_Ennq20/F");
  Pvent -> Branch("Events_clust_Ennq19", &Pevt.Events_clust_Ennq19, "Events_clust_Ennq19/F");
  Pvent -> Branch("Events_clust_Enpq1", &Pevt.Events_clust_Enpq1, "Events_clust_Enpq1/F");
  Pvent -> Branch("Events_clust_Enpq20", &Pevt.Events_clust_Enpq20, "Events_clust_Enpq20/F");
  Pvent -> Branch("Events_clust_Enpq19", &Pevt.Events_clust_Enpq19, "Events_clust_Enpq19/F");
  Pvent -> Branch("Events_clust_Enpq21", &Pevt.Events_clust_Enpq21, "Events_clust_Enpq21/F");
  // no. of (particle) primary tracks branch
  Pvent -> Branch("Events_noOfprimaryTrks", &Pevt.Events_noOfprimaryTrks, "Events_noOfprimaryTrks/I");
  // (particle) primary track branches
  Pvent -> Branch("pTracks_pT", &Pevt.pTracks_pT, "pTracks_pT[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_px", &Pevt.pTracks_px, "pTracks_px[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_py", &Pevt.pTracks_py, "pTracks_py[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_pz", &Pevt.pTracks_pz, "pTracks_pz[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_Eta", &Pevt.pTracks_Eta, "pTracks_Eta[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_Phi", &Pevt.pTracks_Phi, "pTracks_Phi[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_dEdx", &Pevt.pTracks_dEdx, "pTracks_dEdx[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_chrg", &Pevt.pTracks_chrg, "pTracks_chrg[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_gdca", &Pevt.pTracks_gdca, "pTracks_gdca[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_Fp", &Pevt.pTracks_Fp, "pTracks_Fp[Events_noOfprimaryTrks]/I");
  Pvent -> Branch("pTracks_Ppo", &Pevt.pTracks_Ppo, "pTracks_Ppo[Events_noOfprimaryTrks]/I");
  Pvent -> Branch("pTracks_nSigPi", &Pevt.pTracks_nSigPi, "pTracks_nSigPi[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_nSigK", &Pevt.pTracks_nSigK, "pTracks_nSigK[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_nSigP", &Pevt.pTracks_nSigP, "pTracks_nSigP[Events_noOfprimaryTrks]/F");
  Pvent -> Branch("pTracks_nSigE", &Pevt.pTracks_nSigE, "pTracks_nSigE[Events_noOfprimaryTrks]/F");

  // (detector) event branches
  Dvent -> Branch("Events_num", &Devt.Events_num, "Events_num/I");
  Dvent -> Branch("Events_Process", &Devt.Events_Process, "Events_Process/I");
  Dvent -> Branch("Events_refmult", &Devt.Events_refmult, "Events_refmult/I");
  Dvent -> Branch("Events_refPos", &Devt.Events_refPos, "Events_refPos/I");
  Dvent -> Branch("Events_refNeg", &Devt.Events_refNeg, "Events_refNeg/I");
  Dvent -> Branch("Events_runId", &Devt.Events_runId, "Events_runId/I");
  Dvent -> Branch("Events_MagF", &Devt.Events_MagF, "Events_MagF/D");
  Dvent -> Branch("Events_nVertex", &Devt.Events_nVertex, "Events_nVertex/I");
  Dvent -> Branch("Events_rankV", &Devt.Events_rankV, "Events_rankV");
  Dvent -> Branch("Events_pvpdz", &Devt.Events_pvpdz, "Events_pvpdz/D");
  Dvent -> Branch("Events_BBCz", &Devt.Events_BBCz, "Events_BBCz/D");
  Dvent -> Branch("Events_ZDCz", &Devt.Events_ZDCz, "Events_ZDCz/D");
  Dvent -> Branch("Events_primVx", &Devt.Events_primVx, "Events_primVx/D");
  Dvent -> Branch("Events_primVy", &Devt.Events_primVy, "Events_primVy/D");
  Dvent -> Branch("Events_primVz", &Devt.Events_primVz, "Events_primVz/D");
  Dvent -> Branch("Events_TrigId", &Devt.Events_TrigId, "Events_TrigId/I");
  Dvent -> Branch("Events_TrigStat", &Devt.Events_TrigStat, "Events_TrigStat/I");
  Dvent -> Branch("Events_TrigIndex", &Devt.Events_TrigIndex, "Events_TrigIndex/I");
  Dvent -> Branch("Events_PartonId_NS", &Devt.Events_PartonId_NS, "Events_PartonId_NS/I");
  Dvent -> Branch("Events_PartonId_AS", &Devt.Events_PartonId_AS, "Events_PartonId_AS/I");
  Dvent -> Branch("Events_PartonStat_NS", &Devt.Events_PartonStat_NS, "Events_PartonStat_NS/I");
  Dvent -> Branch("Events_PartonStat_AS", &Devt.Events_PartonStat_AS, "Events_PartonStat_AS/I");
  Dvent -> Branch("Events_PartonEta_NS", &Devt.Events_PartonEta_NS, "Events_PartonEta_NS/F");
  Dvent -> Branch("Events_PartonEta_AS", &Devt.Events_PartonEta_AS, "Events_PartonEta_AS/F");
  Dvent -> Branch("Events_PartonPhi_NS", &Devt.Events_PartonPhi_NS, "Events_PartonPhi_NS/F");
  Dvent -> Branch("Events_PartonPhi_AS", &Devt.Events_PartonPhi_AS, "Events_PartonPhi_AS/F");
  Dvent -> Branch("Events_PartonE_NS", &Devt.Events_PartonE_NS, "Events_PartonE_NS/F");
  Dvent -> Branch("Events_PartonE_AS", &Devt.Events_PartonE_AS, "Events_PartonE_AS/F");
  Dvent -> Branch("Events_PartonEt_NS", &Devt.Events_PartonEt_NS, "Events_PartonEt_NS/F");
  Dvent -> Branch("Events_PartonEt_AS", &Devt.Events_PartonEt_AS, "Events_PartonEt_AS/F");
  Dvent -> Branch("Events_tsp", &Devt.Events_tsp, "Events_tsp/F");
  Dvent -> Branch("Events_iso", &Devt.Events_iso, "Events_iso/F");
  // (detector) tower branches
  Dvent -> Branch("Events_Twr_didT", &Devt.Events_Twr_didT, "Events_Twr_didT/I");
  Dvent -> Branch("Events_Twr_adc11", &Devt.Events_Twr_adc11, "Events_Twr_adc11/I");
  Dvent -> Branch("Events_Twr_eneT0", &Devt.Events_Twr_eneT0, "Events_Twr_eneT0/F");
  Dvent -> Branch("Events_Twr_eT", &Devt.Events_Twr_eT, "Events_Twr_eT/F");
  Dvent -> Branch("Events_Twr_ENET0", &Devt.Events_Twr_ENET0, "Events_Twr_ENET0/F");
  Dvent -> Branch("Events_Twr_phT", &Devt.Events_Twr_phT, "Events_Twr_phT/F");
  Dvent -> Branch("Events_Twr_PTower", &Devt.Events_Twr_PTower, "Events_Twr_PTower/F");
  Dvent -> Branch("Events_Twr_pidTower", &Devt.Events_Twr_pidTower, "Events_Twr_pidTower/F");
  Dvent -> Branch("Events_Twr_moduleT", &Devt.Events_Twr_moduleT, "Events_Twr_moduleT/F");
  // (detector) tower cluster branches
  Dvent -> Branch("Events_Clust_EneT0", &Devt.Events_Clust_EneT0, "Events_Clust_EneT0/F");
  Dvent -> Branch("Events_Clust_EneR", &Devt.Events_Clust_EneR, "Events_Clust_EneR/F");
  Dvent -> Branch("Events_Clust_EneH", &Devt.Events_Clust_EneH, "Events_Clust_EneH/F");
  Dvent -> Branch("Events_Clust_rapv1", &Devt.Events_Clust_rapv1, "Events_Clust_rapv1/F");
  Dvent -> Branch("Events_Clust_etav1", &Devt.Events_Clust_etav1, "Events_Clust_etav1/F");
  Dvent -> Branch("Events_Clust_phiv1", &Devt.Events_Clust_phiv1, "Events_Clust_phiv1/F");
  // (detector) SMD eta strip branches
  Dvent -> Branch("Events_Estrp_en01", &Devt.Events_Estrp_en01, "Events_Estrp_en01/F");
  Dvent -> Branch("Events_Estrp_en02", &Devt.Events_Estrp_en02, "Events_Estrp_en02/F");
  Dvent -> Branch("Events_Estrp_en03", &Devt.Events_Estrp_en03, "Events_Estrp_en03/F");
  Dvent -> Branch("Events_Estrp_en0", &Devt.Events_Estrp_en0, "Events_Estrp_en0/F");
  Dvent -> Branch("Events_Estrp_en1", &Devt.Events_Estrp_en1, "Events_Estrp_en1/F");
  Dvent -> Branch("Events_Estrp_en2", &Devt.Events_Estrp_en2, "Events_Estrp_en2/F");
  Dvent -> Branch("Events_Estrp_en3", &Devt.Events_Estrp_en3, "Events_Estrp_en3/F");
  Dvent -> Branch("Events_Estrp_en4", &Devt.Events_Estrp_en4, "Events_Estrp_en4/F");
  Dvent -> Branch("Events_Estrp_en5", &Devt.Events_Estrp_en5, "Events_Estrp_en5/F");
  Dvent -> Branch("Events_Estrp_en6", &Devt.Events_Estrp_en6, "Events_Estrp_en6/F");
  Dvent -> Branch("Events_Estrp_en7", &Devt.Events_Estrp_en7, "Events_Estrp_en7/F");
  Dvent -> Branch("Events_Estrp_en8", &Devt.Events_Estrp_en8, "Events_Estrp_en8/F");
  Dvent -> Branch("Events_Estrp_en9", &Devt.Events_Estrp_en9, "Events_Estrp_en9/F");
  Dvent -> Branch("Events_Estrp_en10", &Devt.Events_Estrp_en10, "Events_Estrp_en10/F");
  Dvent -> Branch("Events_Estrp_en11", &Devt.Events_Estrp_en11, "Events_Estrp_en11/F");
  Dvent -> Branch("Events_Estrp_en12", &Devt.Events_Estrp_en12, "Events_Estrp_en12/F");
  Dvent -> Branch("Events_Estrp_en13", &Devt.Events_Estrp_en13, "Events_Estrp_en13/F");
  Dvent -> Branch("Events_Estrp_en14", &Devt.Events_Estrp_en14, "Events_Estrp_en14/F");
  Dvent -> Branch("Events_Estrp_en15", &Devt.Events_Estrp_en15, "Events_Estrp_en15/F");
  // another (detector) tower branch
  Dvent -> Branch("Events_Twr_didE", &Devt.Events_Twr_didE, "Events_Twr_didE/I");
  // (detector) SMD phi strip branches
  Dvent -> Branch("Events_Pstrip_enp01", &Devt.Events_Pstrip_enp01, "Events_Pstrip_enp01/F");
  Dvent -> Branch("Events_Pstrip_enp02", &Devt.Events_Pstrip_enp02, "Events_Pstrip_enp02/F");
  Dvent -> Branch("Events_Pstrip_enp03", &Devt.Events_Pstrip_enp03, "Events_Pstrip_enp03/F");
  Dvent -> Branch("Events_Pstrip_enp0", &Devt.Events_Pstrip_enp0, "Events_Pstrip_enp0/F");
  Dvent -> Branch("Events_Pstrip_enp1", &Devt.Events_Pstrip_enp1, "Events_Pstrip_enp1/F");
  Dvent -> Branch("Events_Pstrip_enp2", &Devt.Events_Pstrip_enp2, "Events_Pstrip_enp2/F");
  Dvent -> Branch("Events_Pstrip_enp3", &Devt.Events_Pstrip_enp3, "Events_Pstrip_enp3/F");
  Dvent -> Branch("Events_Pstrip_enp4", &Devt.Events_Pstrip_enp4, "Events_Pstrip_enp4/F");
  Dvent -> Branch("Events_Pstrip_enp5", &Devt.Events_Pstrip_enp5, "Events_Pstrip_enp5/F");
  Dvent -> Branch("Events_Pstrip_enp6", &Devt.Events_Pstrip_enp6, "Events_Pstrip_enp6/F");
  Dvent -> Branch("Events_Pstrip_enp7", &Devt.Events_Pstrip_enp7, "Events_Pstrip_enp7/F");
  Dvent -> Branch("Events_Pstrip_enp8", &Devt.Events_Pstrip_enp8, "Events_Pstrip_enp8/F");
  Dvent -> Branch("Events_Pstrip_enp9", &Devt.Events_Pstrip_enp9, "Events_Pstrip_enp9/F");
  Dvent -> Branch("Events_Pstrip_enp10", &Devt.Events_Pstrip_enp10, "Events_Pstrip_enp10/F");
  Dvent -> Branch("Events_Pstrip_enp11", &Devt.Events_Pstrip_enp11, "Events_Pstrip_enp11/F");
  Dvent -> Branch("Events_Pstrip_enp12", &Devt.Events_Pstrip_enp12, "Events_Pstrip_enp12/F");
  Dvent -> Branch("Events_Pstrip_enp13", &Devt.Events_Pstrip_enp13, "Events_Pstrip_enp13/F");
  Dvent -> Branch("Events_Pstrip_enp14", &Devt.Events_Pstrip_enp14, "Events_Pstrip_enp14/F");
  Dvent -> Branch("Events_Pstrip_enp15", &Devt.Events_Pstrip_enp15, "Events_Pstrip_enp15/F");
  // (detector) SMD cluster branches
  Dvent -> Branch("Events_clust_Ennq1", &Devt.Events_clust_Ennq1, "Events_clust_Ennq1/F");
  Dvent -> Branch("Events_clust_Ennq20", &Devt.Events_clust_Ennq20, "Events_clust_Ennq20/F");
  Dvent -> Branch("Events_clust_Ennq19", &Devt.Events_clust_Ennq19, "Events_clust_Ennq19/F");
  Dvent -> Branch("Events_clust_Enpq1", &Devt.Events_clust_Enpq1, "Events_clust_Enpq1/F");
  Dvent -> Branch("Events_clust_Enpq20", &Devt.Events_clust_Enpq20, "Events_clust_Enpq20/F");
  Dvent -> Branch("Events_clust_Enpq19", &Devt.Events_clust_Enpq19, "Events_clust_Enpq19/F");
  Dvent -> Branch("Events_clust_Enpq21", &Devt.Events_clust_Enpq21, "Events_clust_Enpq21/F");
  // no. of (detector) primary tracks branch
  Dvent -> Branch("Events_noOfprimaryTrks", &Devt.Events_noOfprimaryTrks, "Events_noOfprimaryTrks/I");
  // (detector) primary track branches
  Dvent -> Branch("pTracks_pT", &Devt.pTracks_pT, "pTracks_pT[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_px", &Devt.pTracks_px, "pTracks_px[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_py", &Devt.pTracks_py, "pTracks_py[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_pz", &Devt.pTracks_pz, "pTracks_pz[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_Eta", &Devt.pTracks_Eta, "pTracks_Eta[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_Phi", &Devt.pTracks_Phi, "pTracks_Phi[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_dEdx", &Devt.pTracks_dEdx, "pTracks_dEdx[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_chrg", &Devt.pTracks_chrg, "pTracks_chrg[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_gdca", &Devt.pTracks_gdca, "pTracks_gdca[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_Fp", &Devt.pTracks_Fp, "pTracks_Fp[Events_noOfprimaryTrks]/I");
  Dvent -> Branch("pTracks_Ppo", &Devt.pTracks_Ppo, "pTracks_Ppo[Events_noOfprimaryTrks]/I");
  Dvent -> Branch("pTracks_nSigPi", &Devt.pTracks_nSigPi, "pTracks_nSigPi[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_nSigK", &Devt.pTracks_nSigK, "pTracks_nSigK[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_nSigP", &Devt.pTracks_nSigP, "pTracks_nSigP[Events_noOfprimaryTrks]/F");
  Dvent -> Branch("pTracks_nSigE", &Devt.pTracks_nSigE, "pTracks_nSigE[Events_noOfprimaryTrks]/F");

  // track histograms (for detector effects) and QA
  pT_trks     = new TH1D("pT_trks", "p_{T} of all (particle-lvl.) tracks", 100, 0., 10.);
  pT_accTrks  = new TH1D("pT_accTrks", "p_{T} of surviving (detector-lvl.) tracks", 100, 0., 10.);
  pT_tpcReso  = new TH2D("pT_tpcReso", "TPC resolution", 100, 0., 10., 100, -0.2, 0.2);
  pT_tpcResoP = new TProfile("pT_tpcResoP", "TPC resolution, profile", 100, 0., 10., -0.5, 0.5, "S");
  E_calReso   = new TH2D("E_calReso", "Cal. resolution", 100, 0., 10., 100, -0.2, 0.2);
  E_calResoP  = new TProfile("E_calResoP", "Cal. resolution, profile", 100, 0., 10., -0.5, 0.5, "S");
  pT_smear    = new TH2D("pT_smear", "Smeared p_{T} vs. true p_{T}", 100, 0., 10., 100, 0., 10.);
  pT_smearP   = new TProfile("pT_smearP", "Smeared p_{T} vs. true p_{T}, profile", 100, 0., 10., 0., 20., "S");
  zTparton    = new TH1D("zTtrig", "p_{T}^{parton}/p_{T}^{trig}", 200, 0., 20.);
  zTparton -> Sumw2();
  zTns        = new TH1D("zTns", "p_{T}^{NS}/p_{T}^{trig}", 200, 0., 20.);
  zTns     -> Sumw2();
  zTas        = new TH1D("zTas", "p_{T}^{AS}/p_{T}^{trig}", 200, 0., 20.);
  zTas     -> Sumw2();

}  // end 'init(float, float, float, TFile)'



void Analysis::analyze(const int code, Pythia8::Event& evt) {

  ++nEvt;
  event      = evt;
  trig       = false;
  trig2      = false;
  sys.reset();
  Ptrk.resize(event.size());
  Dtrk.resize(event.size());

  // seed random number generator
  srand(time(NULL));


  // determine how many photons were produced in hard process
  int nG = 0;
  switch (code) {
    case 201:
      nG = 1;
      break;
    case 202:
      nG = 1;
      break;
    case 203:
      nG = 1;
      break;
    case 204:
      nG = 2;
      break;
    case 205:
      nG = 2;
      break;
    default:
      nG = 0;
      break;
  }


  // find trigger
  int   tID     = 0;
  int   tStat   = 0;
  bool  sort    = false;
  for (int i = 0; i < event.size(); ++i) {

    // locate hard scatter products
    if (event[i].status() == -23) sys.prdx.push_back(i);


    // trigger acceptance
    int id = event[i].idAbs();
    if (!event[i].isFinal()) continue;
    if ((id != 22) && (id != 111)) continue;
    if ((abs(event[i].eta()) < hCut) && (event[i].eT() > eCut)) {

      // make sure we're looking for the correct trigger
      tStat   = event[i].status();
      tID     = event[i].id();
      sys.eID = tID + nG;
      if ((tID == 111)  && (nG != 0.))    continue;
      if ((tID == 22)   && (nG == 0.))    continue;
      if ((tID == 22)   && (tStat != 62)) continue;

      // grab trigger info
      Ptrg.setTrigger(i, event[i]);
      Dtrg.setTrigger(i, event[i]);
      sys.tIndex = i;

      // determine which products are near side and away side
      if (sys.eID == 111) locateNearSide(i);
      sort = sortPartons(i);
      trig = true;
      ++sys.nTrg;
      break;

    }  // end trigger cuts

  }  // end particle loop


  // analyze relevant event
  if (trig) {

    int    id        = 0;
    int    tCode     = 0;
    int    refP      = 0;
    int    refD      = 0;
    int    refW      = 0;
    bool   trkAccept = false;
    float  dR        = 0.;
    double dE        = 0.;
    double resTPC    = 0.;
    double resCal    = 0.;
    for (int i = 0; i < event.size(); ++i) {

      // don't include trigger!
      if (i == sys.tIndex) continue;

      // assign code to track
      if (i > 4) tCode = sortTrack(i, sys.tIndex); 
      if (!event[i].isFinal())   continue;

      // set Cal. and TPC resolution
      resCal = smearPt(event[i].e(), event[i].id(), "cal");
      if (event[i].isCharged())
        resTPC = smearPt(event[i].pT(), event[i].id(), "tpc");


      // particle level tracks will always survive
      Ptrk[refP].aDet = false; 
      dE = Ptrk[refP].setTrack(i, tCode, resTPC, resCal, event[i]);
      if (dE > -1000.)
        trkAccept = true;
      else
        trkAccept = false; 


      if (trkAccept) {

        // store pT of all tracks
        pT_trks -> Fill(Ptrk[refP].aPt);

        Ptrk[refP].aDphi = Ptrk[refP].aPhi - Ptrg.iPhi;
        if (Ptrk[refP].aDphi < -M_PI/2)  Ptrk[refP].aDphi += 2*M_PI;
        if (Ptrk[refP].aDphi > 3*M_PI/2) Ptrk[refP].aDphi -= 2*M_PI;

        // grab track identity
        id = event[i].idAbs();
        switch (id) {
          case 111:
            Ptrk[refP].aSigPi = event[i].id();
            break;
          case 211:
            Ptrk[refP].aSigPi = event[i].id();
            break;
          case 311:
            Ptrk[refP].aSigK  = event[i].id();
            break;
          case 313:
            Ptrk[refP].aSigK  = event[i].id();
            break;
          case 321:
            Ptrk[refP].aSigK  = event[i].id();
            break;
          case 323:
            Ptrk[refP].aSigK  = event[i].id();
            break;
          case 2212:
            Ptrk[refP].aSigP  = event[i].id();
            break;
          default:
            Ptrk[refP].aSigE  = event[i].id();
            break;
        }

        // calculate isolation ratio
        if (tID == 22) {
          dR = Ptrg.deltaR(event[i]);
          if (dR < rIso) Ptrg.addEne(event[i]);
        }

        // calculate no. of vertices / multiplicities
        calculateVtx(i, Ptrk[refP].aDet);
        if (event[i].charge() > 0.) ++sys.posP;
        if (event[i].charge() < 0.) ++sys.negP;
        ++refP;

      }  // end particle track info


      // apply detector smearing
      Dtrk[refD].aDet = true;
      dE = Dtrk[refD].setTrack(i, tCode, resTPC, resCal, event[i]);
      if (dE > -1000.)
        trkAccept = true;
      else
        trkAccept = false;

      E_calReso  -> Fill(event[i].e(), resCal);
      E_calResoP -> Fill(event[i].e(), resCal);
      if (event[i].isCharged()) {
        pT_tpcReso  -> Fill(event[i].pT(), resTPC);
        pT_tpcResoP -> Fill(event[i].pT(), resTPC);
      }
      pT_smear  -> Fill(Ptrk[refP-1].aPt, Dtrk[refD].aPt);
      pT_smearP -> Fill(Ptrk[refP-1].aPt, Dtrk[refD].aPt);


      // if track survives
      if (trkAccept) {

        // store pT of surviving tracks
        pT_accTrks -> Fill(Dtrk[refD].aPt);

        Dtrk[refD].aDphi = Dtrk[refD].aPhi - Dtrg.iPhi;
        if (Dtrk[refD].aDphi < -M_PI/2)  Dtrk[refD].aDphi += 2*M_PI;
        if (Dtrk[refD].aDphi > 3*M_PI/2) Dtrk[refD].aDphi -= 2*M_PI;

        // grab track identity
        id = event[i].idAbs();
        switch (id) {
          case 111:
            Dtrk[refD].aSigPi = event[i].id();
            break;
          case 211:
            Dtrk[refD].aSigPi = event[i].id();
            break;
          case 311:
            Dtrk[refD].aSigK  = event[i].id();
            break;
          case 313:
            Dtrk[refD].aSigK  = event[i].id();
            break;
          case 321:
            Dtrk[refD].aSigK  = event[i].id();
            break;
          case 323:
            Dtrk[refD].aSigK  = event[i].id();
            break;
          case 2212:
            Dtrk[refD].aSigP  = event[i].id();
            break;
          default:
            Dtrk[refD].aSigE  = event[i].id();
            break;
        }

        // calculate isolation ratio
        if (tID == 22) {
          dR = Dtrg.deltaR(event[i]);
          if (dR < rIso) Dtrg.addEne(event[i]);
        }

        // calculate no. of vertices / multiplicities
        calculateVtx(i, Dtrk[refD].aDet);
        if (event[i].charge() > 0.) ++sys.posD;
        if (event[i].charge() < 0.) ++sys.negD;
        ++refD;

      }  // end detector track info


      // for 'weird' tracks (when dE != 0); assume pion mass.
      if (dE > 0.2) {
        Wtrk.push_back(Dtrk[refD-1]);
        Wtrk[refW].setComponents(dE, 0.140);
        ++refW;
      }

    }  // end particle loop


    // add 'weird' tracks to Dtrk
    int nW = Wtrk.size();
    for (int i = 0; i < nW; ++i) {
      Dtrk.push_back(Wtrk[i]);
      ++refD;
    }


    sys.refP = refP;
    sys.refD = refD;
    Ptrg.setIso();
    Dtrg.setIso();
    Pvtx.setPrimeVtx();
    Dvtx.setPrimeVtx();
    fillTrees(code);

  }  // end relevant event


  // save trees and tuples
  if (nEvt % 100000 == 0) {
    Pvent   -> AutoSave();
    Dvent   -> AutoSave();
  }


  // clear particle objects
  Pvtx.clear();
  Ptrk.clear();
  Ptrg.clear();
  Pasp.clear();
  Pnsp.clear();
  // clear detector objects
  Dvtx.clear();
  Wtrk.clear();
  Dtrk.clear();
  Dtrg.clear();
  Dasp.clear();
  Dnsp.clear();

}  // end 'analyze(int, Event&)'



void Analysis::finish() {

  cout << "\n  There were " << sys.nTrg << " triggers!" << endl;

  // write trees / tuples
  Pvent       -> Write();
  Dvent       -> Write();
  // write histograms
  pT_trks     -> Write();
  pT_accTrks  -> Write();
  pT_tpcReso  -> Write();
  pT_tpcResoP -> Write();
  E_calReso   -> Write();
  E_calResoP  -> Write();
  pT_smear    -> Write();
  pT_smearP   -> Write();
  zTparton    -> Write();
  zTns        -> Write();
  zTas        -> Write();
  // close files
  file        -> Close();

}  // end 'finish()'



// private methods ------------------------------------------------------------

int Analysis::sortTrack(const int iTrack, const int iTrig) {

  int   trkCode = -1;
  bool  near    = false;
  bool  away    = false;
  bool  unc     = false;
  bool  TR      = false;
  bool  RE      = false;
  bool  AR1     = false;
  bool  AR2     = false;

  float trkDphi = event[iTrack].phi() - event[iTrig].phi();
  if (trkDphi < 2*M_PI) trkDphi += 2*M_PI;
  if (trkDphi > 2*M_PI) trkDphi -= 2*M_PI;

  // determine what spawned track
  away = addAwayTrk(iTrack, false);
  if (!away) near = addNearTrk(iTrack, false);
  if (!away && !near) unc = addUncorrelatedTrk(iTrack);

  // determine which Delta-Phi region track lies in
  float rPhi  = trkDphi - M_PI;
  float u1Phi = trkDphi - M_PI/2.;
  float u2Phi = trkDphi - 3.*M_PI/2.;
  if ((abs(trkDphi) < trCut) || (abs(trkDphi-2*M_PI) < trCut)) {
    TR = true;
    if (event[iTrack].isFinal()) {
      sys.pTtr += event[iTrack].pT();
      sys.nTR++;
    }
  }
  if (abs(rPhi) < reCut) {
    RE = true;
    if (event[iTrack].isFinal()) {
      sys.pTre += event[iTrack].pT();
      sys.nRE++;
    }
  }
  if ((trkDphi > a1Cut) && (trkDphi < a2Cut)) {
    AR1 = true;
    if (event[iTrack].isFinal()) {
      sys.pTar1 += event[iTrack].pT();
      sys.nAR1++;
      if ((event[iTrack].pT() > 3) && !trig2) {
        trig2       = true;
        sys.t2Index = iTrack;
        sys.t2Pt    = event[iTrack].pT();
        sys.t2Phi   = event[iTrack].phi();
      }
    }
  }
  if ((rPhi > a1Cut) && (rPhi < a2Cut)) {
    AR2 = true;
    if (event[iTrack].isFinal()) {
      sys.pTar2 += event[iTrack].pT();
      sys.nAR2++;
      if ((event[iTrack].pT() > 3) && !trig2) {
        trig2       = true;
        sys.t2Index = iTrack;
        sys.t2Pt    = event[iTrack].pT();
        sys.t2Phi   = event[iTrack].phi();
      }
    }
  }


  // assign code accordingly
  if (near && TR)  trkCode = 0;
  if (near && RE)  trkCode = 1;
  if (near && AR1) trkCode = 2;
  if (near && AR2) trkCode = 3;
  if (away && TR)  trkCode = 10;
  if (away && RE)  trkCode = 11;
  if (away && AR1) trkCode = 12;
  if (away && AR2) trkCode = 13;
  if (unc  && TR)  trkCode = 20;
  if (unc  && RE)  trkCode = 21;
  if (unc  && AR1) trkCode = 22;
  if (unc  && AR2) trkCode = 23;

  float pT = event[iTrack].pT();
  if (event[iTrack].isParton()) zTparton -> Fill(pT/Ptrg.iPt);

  // determine hardest track in AR
  if (event[iTrack].isFinal() && (AR1 || AR2)) {
    if (pT > sys.leadUPt) {
      sys.leadUIndex = iTrack;
      sys.leadUPt    = pT;
      sys.leadUPhi   = event[iTrack].phi();
      if (pT > 3) {
        if ((abs(u1Phi) < ueCut) || (abs(u2Phi) < ueCut))
          sys.uID = 2;
        else
          sys.uID = 1;
      }
    }
  }

  return trkCode;

}  // end 'sortTrack(int, int)'



void Analysis::fillTrees(const int evtCode) {

  // grab (particle) event info
  Pevt.Events_num             = nEvt;
  Pevt.Events_Process         = evtCode;
  Pevt.Events_refmult         = sys.refP;
  Pevt.Events_refPos          = sys.posP;
  Pevt.Events_refNeg          = sys.negP;
  Pevt.Events_runId           = sys.eID;
  Pevt.Events_nVertex         = Pvtx.nVertex;
  Pevt.Events_primVx          = Pvtx.primVx;
  Pevt.Events_primVy          = Pvtx.primVy;
  Pevt.Events_primVz          = Pvtx.primVz;
  Pevt.Events_TrigId          = Ptrg.iID;
  Pevt.Events_TrigStat        = Ptrg.iStat;
  Pevt.Events_TrigIndex       = Ptrg.iIndex;
  Pevt.Events_PartonId_NS     = Pnsp.pID;
  Pevt.Events_PartonId_AS     = Pasp.pID;
  Pevt.Events_PartonStat_NS   = Pnsp.pStat;
  Pevt.Events_PartonStat_AS   = Pasp.pStat;
  Pevt.Events_PartonEta_NS    = Pnsp.pEta;
  Pevt.Events_PartonEta_AS    = Pasp.pEta;
  Pevt.Events_PartonPhi_NS    = Pnsp.pPhi;
  Pevt.Events_PartonPhi_AS    = Pasp.pPhi;
  Pevt.Events_PartonE_NS      = Pnsp.pEne;
  Pevt.Events_PartonE_AS      = Pasp.pEne;
  Pevt.Events_PartonEt_NS     = Pnsp.pEneT;
  Pevt.Events_PartonEt_AS     = Pasp.pEneT;
  Pevt.Events_iso             = Ptrg.iso;
  Pevt.Events_Twr_didT        = sys.leadUIndex;
  Pevt.Events_Twr_adc11       = sys.t2Index;
  Pevt.Events_Twr_eneT0       = sys.leadUPt;
  Pevt.Events_Twr_eT          = sys.t2Pt;
  Pevt.Events_Twr_phT         = sys.leadUPhi;
  Pevt.Events_Twr_PTower      = sys.t2Phi;
  Pevt.Events_Twr_pidTower    = sys.uID;
  Pevt.Events_Clust_EneT0     = Ptrg.iEne;
  Pevt.Events_Clust_EneR      = Ptrg.iEneR;
  Pevt.Events_Clust_EneH      = Ptrg.iEneH;
  Pevt.Events_Clust_rapv1     = Ptrg.iRap;
  Pevt.Events_Clust_etav1     = Ptrg.iEta;
  Pevt.Events_Clust_phiv1     = Ptrg.iPhi;
  Pevt.Events_noOfprimaryTrks = sys.refP;
  for (int i = 0; i < Pevt.Events_noOfprimaryTrks; ++i) {
    Pevt.pTracks_pT[i]     = Ptrk[i].aPt;
    Pevt.pTracks_px[i]     = Ptrk[i].aPx;
    Pevt.pTracks_py[i]     = Ptrk[i].aPy;
    Pevt.pTracks_pz[i]     = Ptrk[i].aPz;
    Pevt.pTracks_Eta[i]    = Ptrk[i].aEta;
    Pevt.pTracks_Phi[i]    = Ptrk[i].aPhi;
    Pevt.pTracks_dEdx[i]   = Ptrk[i].aID;
    Pevt.pTracks_chrg[i]   = Ptrk[i].aChrg;
    Pevt.pTracks_gdca[i]   = Ptrk[i].aMass;
    Pevt.pTracks_Fp[i]     = Ptrk[i].aCode;
    Pevt.pTracks_Ppo[i]    = Ptrk[i].aIndex;
    Pevt.pTracks_nSigPi[i] = Ptrk[i].aSigPi;
    Pevt.pTracks_nSigK[i]  = Ptrk[i].aSigK;
    Pevt.pTracks_nSigP[i]  = Ptrk[i].aSigP;
    Pevt.pTracks_nSigE[i]  = Ptrk[i].aSigE;
  }
  // grab (detector) event info
  Devt.Events_num             = nEvt;
  Devt.Events_Process         = evtCode;
  Devt.Events_refmult         = sys.refD;
  Devt.Events_refPos          = sys.posD;
  Devt.Events_refNeg          = sys.negD;
  Devt.Events_runId           = sys.eID;
  Devt.Events_nVertex         = Dvtx.nVertex;
  Devt.Events_primVx          = Dvtx.primVx;
  Devt.Events_primVy          = Dvtx.primVy;
  Devt.Events_primVz          = Dvtx.primVz;
  Devt.Events_TrigId          = Dtrg.iID;
  Devt.Events_TrigStat        = Dtrg.iStat;
  Devt.Events_TrigIndex       = Dtrg.iIndex;
  Devt.Events_PartonId_NS     = Dnsp.pID;
  Devt.Events_PartonId_AS     = Dasp.pID;
  Devt.Events_PartonStat_NS   = Dnsp.pStat;
  Devt.Events_PartonStat_AS   = Dasp.pStat;
  Devt.Events_PartonEta_NS    = Dnsp.pEta;
  Devt.Events_PartonEta_AS    = Dasp.pEta;
  Devt.Events_PartonPhi_NS    = Dnsp.pPhi;
  Devt.Events_PartonPhi_AS    = Dasp.pPhi;
  Devt.Events_PartonE_NS      = Dnsp.pEne;
  Devt.Events_PartonE_AS      = Dasp.pEne;
  Devt.Events_PartonEt_NS     = Dnsp.pEneT;
  Devt.Events_PartonEt_AS     = Dasp.pEneT;
  Devt.Events_iso             = Dtrg.iso;
  Devt.Events_Twr_didT        = sys.leadUIndex;
  Devt.Events_Twr_adc11       = sys.t2Index;
  Devt.Events_Twr_eneT0       = sys.leadUPt;
  Devt.Events_Twr_eT          = sys.t2Pt;
  Devt.Events_Twr_phT         = sys.leadUPhi;
  Devt.Events_Twr_PTower      = sys.t2Phi;
  Devt.Events_Twr_pidTower    = sys.uID;
  Devt.Events_Clust_EneT0     = Dtrg.iEne;
  Devt.Events_Clust_EneR      = Dtrg.iEneR;
  Devt.Events_Clust_EneH      = Dtrg.iEneH;
  Devt.Events_Clust_rapv1     = Dtrg.iRap;
  Devt.Events_Clust_etav1     = Dtrg.iEta;
  Devt.Events_Clust_phiv1     = Dtrg.iPhi;
  Devt.Events_noOfprimaryTrks = sys.refD;
  for (int i = 0; i < Devt.Events_noOfprimaryTrks; ++i) {
    Devt.pTracks_pT[i]     = Dtrk[i].aPt;
    Devt.pTracks_px[i]     = Dtrk[i].aPx;
    Devt.pTracks_py[i]     = Dtrk[i].aPy;
    Devt.pTracks_pz[i]     = Dtrk[i].aPz;
    Devt.pTracks_Eta[i]    = Dtrk[i].aEta;
    Devt.pTracks_Phi[i]    = Dtrk[i].aPhi;
    Devt.pTracks_dEdx[i]   = Dtrk[i].aID;
    Devt.pTracks_chrg[i]   = Dtrk[i].aChrg;
    Devt.pTracks_gdca[i]   = Dtrk[i].aMass;
    Devt.pTracks_Fp[i]     = Dtrk[i].aCode;
    Devt.pTracks_Ppo[i]    = Dtrk[i].aIndex;
    Devt.pTracks_nSigPi[i] = Dtrk[i].aSigPi;
    Devt.pTracks_nSigK[i]  = Dtrk[i].aSigK;
    Devt.pTracks_nSigP[i]  = Dtrk[i].aSigP;
    Devt.pTracks_nSigE[i]  = Dtrk[i].aSigE;
  }

  // fill trees
  Pvent   -> Fill();
  Dvent   -> Fill();
  zTns    -> Fill(Pnsp.pEneT, Ptrg.iPt);
  zTas    -> Fill(Pasp.pEneT, Ptrg.iPt);

}  // end 'fillTrees(int)'



void Analysis::locateNearSide(const int iTrig) {

  int   nKids, nPrd;
  int   mindex, maxdex;
  int   pIndex, kIndex;
  float kDphi, minDphi, maxDphi;
  vector<int> kids;

  mindex  = 0;
  maxdex  = 0;
  minDphi = M_PI;
  maxDphi = 0.;
  nPrd    = sys.prdx.size();
  for (int i = 0; i < nPrd; ++i) {

    pIndex = sys.prdx[i];
    kids   = event[pIndex].daughterList();
    nKids  = kids.size();
    for (int j = 0; j < nKids; ++j) {

      kIndex = kids[j];
      kDphi  = event[kIndex].phi() - event[iTrig].phi();
      if (kDphi < -M_PI) kDphi += 2*M_PI;
      if (kDphi > M_PI)  kDphi -= 2*M_PI;

      if (abs(kDphi) < minDphi) {
        minDphi = abs(kDphi);
        mindex  = pIndex;
      }
      if (abs(kDphi) > maxDphi) {
        maxDphi = abs(kDphi);
        maxdex  = pIndex;
      }

    }  // end daughter loop

  }  // end product loop

  sys.mindex  = mindex;
  sys.maxdex  = maxdex;
  sys.minDphi = minDphi;
  sys.maxDphi = maxDphi;

}  // end 'locateNearSide(int)'



void Analysis::calculateVtx(const int iTrack, const bool det) {

  bool vert, add;
  Pythia8::Vec4 prod;

  prod = event[iTrack].vProd();
  vert = event[iTrack].hasVertex();
  // detector level
  if (det) {
    if (Dvtx.nVertex == 0) {
      Dvtx.addVtx(prod);
    }
    else if (vert) {
      add = Dvtx.isNewVtx(prod);
      if (add) Dvtx.addVtx(prod);
    }
    Dvtx.asscVtx(vert, prod);
  }
  // particle level
  else {
    if (Pvtx.nVertex == 0) {
      Pvtx.addVtx(prod);
    }
    else if (vert) {
      add = Pvtx.isNewVtx(prod);
      if (add) Pvtx.addVtx(prod);
    }
    Pvtx.asscVtx(vert, prod);
  }

}  // end 'calculateVtx(int, bool)'



bool Analysis::sortPartons(const int iTrig) {

  int  nPrd;
  int  pIndex, nIndex, aIndex;
  int  rando;
  bool nSet, aSet;
  bool near, away;
  bool success;

  aIndex  = 0;
  nIndex  = 0;
  nSet    = false;
  aSet    = false;
  success = false;
  nPrd    = sys.prdx.size();
  for (int i = 0; i < nPrd; ++i) {

    pIndex = sys.prdx[i];
    switch (sys.eID) {

      case 23:
        if (event[pIndex].idAbs() == 22) {
          nSet   = true;
          nIndex = pIndex;
          Pnsp.setParton(pIndex, event[pIndex]);
          Dnsp.setParton(pIndex, event[pIndex]);
          if (event[iTrig].isAncestor(pIndex)) {
            Ptrg.iPrompt = true;
            Dtrg.iPrompt = true;
          }
        }
        else {
          aSet   = true;
          aIndex = pIndex;
          Pasp.setParton(pIndex, event[pIndex]);
          Dasp.setParton(pIndex, event[pIndex]);
        }
        break;

      case 24:
        if (event[iTrig].isAncestor(pIndex)) {
          nSet   = true;
          nIndex = pIndex;
          Pnsp.setParton(pIndex, event[pIndex]);
          Dnsp.setParton(pIndex, event[pIndex]);
          Ptrg.iPrompt = true;
          Dtrg.iPrompt = true;
        }
        else {
          aSet   = true;
          aIndex = pIndex;
          Pasp.setParton(pIndex, event[pIndex]);
          Dasp.setParton(pIndex, event[pIndex]);
        }
        break;

      case 111:
        if (event[pIndex].index() == sys.mindex) {
          nSet   = true;
          nIndex = pIndex;
          Pnsp.setParton(pIndex, event[pIndex]);
          Dnsp.setParton(pIndex, event[pIndex]);
          if (event[iTrig].isAncestor(pIndex)) {
            Ptrg.iPrompt = true;
            Dtrg.iPrompt = true;
          }
        }
        else {
          aSet   = true;
          aIndex = pIndex;
          Pasp.setParton(pIndex, event[pIndex]);
          Dasp.setParton(pIndex, event[pIndex]);
        }
        break;

    }  // end switch(sys.eID)

  }  // end parton loop

  // make sure both partons have been set. If not, randomly assign them.
  if (!nSet || !aSet) {
    rando = rand() % 10;
    if (rando <= 4) {
      nIndex = 5;
      aIndex = 6;
    }
    if (rando > 4) {
      nIndex = 6;
      aIndex = 5;
    }
    Pnsp.setParton(nIndex, event[nIndex]);
    Pasp.setParton(aIndex, event[aIndex]);
    Dnsp.setParton(nIndex, event[nIndex]);
    Dasp.setParton(aIndex, event[aIndex]);
    Ptrg.iPrompt = false;
    Dtrg.iPrompt = false;
    success      = false;
  }
  else if (nSet && aSet) success = true;

  // add NS and AS to nTrks and aTrks respectively
  near       = addNearTrk(nIndex, true);
  away       = addAwayTrk(aIndex, true);
  sys.nIndex = nIndex;
  sys.aIndex = aIndex;

  return success;

}  // end 'sortPartons(int)'



bool Analysis::addNearTrk(const int iTrack, const bool first) {

  int  Nd;
  bool add, nTrk;
  vector<int> d;

  nTrk = false;
  if (first) {
    sys.nTrks.push_back(iTrack);
    nTrk = true;
  }
  else {
    nTrk = sys.hasNmom(event[iTrack]);
    if (!nTrk) return nTrk;
  }

  d  = event[iTrack].daughterList();
  Nd = d.size();
  for (int i = 0; i < Nd; ++i) {
    if (first) {
      sys.nTrks.push_back(d[i]);
    }
    else {
      add = sys.isNewNtrk(d[i]);
      if (add) sys.nTrks.push_back(d[i]);
    }
  }

  return nTrk;

}  // end 'addNearTrk(int, bool)'



bool Analysis::addAwayTrk(const int iTrack, const bool first) {

  int  Nd;
  bool add, aTrk;
  vector<int> d;

  aTrk = false;
  if (first) {
    sys.aTrks.push_back(iTrack);
    aTrk = true;
  }
  else {
    aTrk = sys.hasAmom(event[iTrack]);
    if (!aTrk) return aTrk;
  }

  d  = event[iTrack].daughterList();
  Nd = d.size();
  for (int i = 0; i < Nd; ++i) {
    if (first) {
      sys.aTrks.push_back(d[i]);
    }
    else {
      add = sys.isNewAtrk(d[i]);
      if (add) sys.aTrks.push_back(d[i]);
    }
  }

  return aTrk;

}  // end 'addAwayTrk(int, bool)'



bool Analysis::addUncorrelatedTrk(const int iTrack) {

  bool add = sys.isNewUtrk(iTrack);
  if (add) sys.uTrks.push_back(iTrack);
  return add;

}  // end 'addUncorrelatedTrk(int)'



double Analysis::smearPt(const double pt, const int id, const TString& source) {

  double res = 0.;
  double ene = 0.;
  double exp = 0.;

  // tpc resolution
  if (source == "tpc")
    res = gRandom->Gaus(0.0,(0.01*pt)) + gRandom->Gaus(0.0,0.017);

  // cal. resolution ("pt" here should be total energy)
  if ((source == "cal") && ((abs(id) == 11) || (abs(id) == 22))) {
    ene = pt;
    exp = 1./sqrt(ene);
    res = gRandom->Gaus(0.0,(0.15*exp)) + gRandom->Gaus(0.0,0.015);
  }

  return res;

}  // end 'smearPt(double)'


// End ========================================================================
