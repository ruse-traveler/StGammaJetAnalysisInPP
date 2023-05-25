// 'GammaJetEvent.h'
// Derek Anderson
// 01.08.2016
//
// This file defines the Event class used in the analysis defined by
// 'GammaJetAnalysis.h' and executed by 'GammaJetExecute.cc'.  Note
// that some of the quantities here won't be relevant for the pythia
// simulation; they're left in for completion's sake.  However, the ones
// that are being used have explanations next to them.
//
// Last updated: 06.09.2016

#pragma once

// max. no. of primary tracks
const Int_t NTmax = 100000;


class Event {

public:

  Event();

  // event variables
  Int_t    Events_num;            // event no.
  Int_t    Events_Process;        // type of hard process which occured
  Int_t    Events_BigUETrack;     // if event has a high-pT track in UE
  Int_t    Events_leadUindex;     // evt. record index of hardest track in UE
  Float_t  Events_leadUpT;        // pT of hardest track in UE
  Int_t    Events_refmult;        // no. of charged tracks
  Int_t    Events_refPos;         // no. of positive tracks
  Int_t    Events_refNeg;         // no. of negative tracks
  Int_t    Events_runId;          // id's how NS / AS partons were determined
  Double_t Events_MagF;
  Int_t    Events_nVertex;        // no. of vertices in event 
  Int_t    Events_rankV;
  Double_t Events_pvpdz;
  Double_t Events_BBCz;
  Double_t Events_ZDCz;
  Double_t Events_primVx;         // prime vtx. coordinates
  Double_t Events_primVy;         // should (almost) always be (0,0,0)
  Double_t Events_primVz;
  Int_t    Events_TrigId;         // trigger PID
  Int_t    Events_TrigStat;       // trigger production code
  Int_t    Events_TrigIndex;      // trigger evt. record index
  Int_t    Events_PartonId_NS;    // near side parton / photon pID
  Int_t    Events_PartonId_AS;    // away side parton / photon pID
  Int_t    Events_PartonStat_NS;  // near side parton / photon prod. code
  Int_t    Events_PartonStat_AS;  // away side parton / photon prod. code
  Float_t  Events_PartonEta_NS;   // near side parton / photon eta
  Float_t  Events_PartonEta_AS;   // away side parton / photon eta
  Float_t  Events_PartonPhi_NS;   // near side parton / photon phi
  Float_t  Events_PartonPhi_AS;   // away side parton / photon phi
  Float_t  Events_PartonE_NS;     // near side parton / photon energy (p)
  Float_t  Events_PartonE_AS;     // away side parton / photon energy (p)
  Float_t  Events_PartonEt_NS;    // near side parton / photon eT (pT)
  Float_t  Events_PartonEt_AS;    // away side parton / photon eT (pT)
  Float_t  Events_tsp;
  Float_t  Events_iso;            // ratio of EneH to EneR in cone
  // tower variables
  Int_t    Events_Twr_didT;
  Int_t    Events_Twr_adc11;
  Float_t  Events_Twr_eneT0;
  Float_t  Events_Twr_eT;
  Float_t  Events_Twr_ENET0;
  Float_t  Events_Twr_phT;
  Float_t  Events_Twr_PTower;
  Float_t  Events_Twr_pidTower;
  Int_t    Events_Twr_moduleT;
  // tower cluster variables
  Float_t  Events_Clust_EneT0;  // total energy of trigger
  Float_t  Events_Clust_EneR;   // total energy in a cone about trigger
  Float_t  Events_Clust_EneH;   // total hadronic energy in same cone
  Float_t  Events_Clust_rapv1;  // trigger y
  Float_t  Events_Clust_etav1;  // trigger eta
  Float_t  Events_Clust_phiv1;  // trigger phi
  // SMD eta strip variables
  Int_t    Events_Twr_didE;
  Float_t  Events_Estrp_en03;
  Float_t  Events_Estrp_en02;
  Float_t  Events_Estrp_en01;
  Float_t  Events_Estrp_en0;
  Float_t  Events_Estrp_en1;
  Float_t  Events_Estrp_en2;
  Float_t  Events_Estrp_en3;
  Float_t  Events_Estrp_en4;
  Float_t  Events_Estrp_en5;
  Float_t  Events_Estrp_en6;
  Float_t  Events_Estrp_en7;
  Float_t  Events_Estrp_en8;
  Float_t  Events_Estrp_en9;
  Float_t  Events_Estrp_en10;
  Float_t  Events_Estrp_en11;
  Float_t  Events_Estrp_en12;
  Float_t  Events_Estrp_en13;
  Float_t  Events_Estrp_en14;
  Float_t  Events_Estrp_en15;
  // SMD phi strip variables
  Float_t  Events_Pstrip_enp03;
  Float_t  Events_Pstrip_enp02;
  Float_t  Events_Pstrip_enp01;
  Float_t  Events_Pstrip_enp0;
  Float_t  Events_Pstrip_enp1;
  Float_t  Events_Pstrip_enp2;
  Float_t  Events_Pstrip_enp3;
  Float_t  Events_Pstrip_enp4;
  Float_t  Events_Pstrip_enp5;
  Float_t  Events_Pstrip_enp6;
  Float_t  Events_Pstrip_enp7;
  Float_t  Events_Pstrip_enp8;
  Float_t  Events_Pstrip_enp9;
  Float_t  Events_Pstrip_enp10;
  Float_t  Events_Pstrip_enp11;
  Float_t  Events_Pstrip_enp12;
  Float_t  Events_Pstrip_enp13;
  Float_t  Events_Pstrip_enp14;
  Float_t  Events_Pstrip_enp15;
  // SMD cluster variables
  Float_t  Events_clust_Ennq1;
  Float_t  Events_clust_Ennq20;
  Float_t  Events_clust_Ennq19;
  Float_t  Events_clust_Enpq1;
  Float_t  Events_clust_Enpq20;
  Float_t  Events_clust_Enpq19;
  Float_t  Events_clust_Enpq21;
  // no. of primary tracks
  Int_t    Events_noOfprimaryTrks;  // no. of charged tracks

  // primary track variables
  Float_t  pTracks_pT[NTmax];      // track pT
  Float_t  pTracks_px[NTmax];      // track px comp.
  Float_t  pTracks_py[NTmax];      // track py comp.
  Float_t  pTracks_pz[NTmax];      // track pz comp.
  Float_t  pTracks_Eta[NTmax];     // track eta
  Float_t  pTracks_Phi[NTmax];     // track phi
  Float_t  pTracks_dEdx[NTmax];
  Float_t  pTracks_chrg[NTmax];    // track charge
  Float_t  pTracks_gdca[NTmax];
  Int_t    pTracks_Fp[NTmax];
  Int_t    pTracks_Ppo[NTmax];
  Float_t  pTracks_nSigPi[NTmax];  // indicates if track is a pion
  Float_t  pTracks_nSigK[NTmax];   // indicates if track is a kaon
  Float_t  pTracks_nSigP[NTmax];   // indicates if track is a proton
  Float_t  pTracks_nSigE[NTmax];   // indicates if track is an exotic

};


// member definitions ---------------------------------------------------------

Event::Event() {

  // initialize event variables
  Events_num           = 0;
  Events_Process       = 0;
  Events_BigUETrack    = 0;
  Events_leadUindex    = 0;
  Events_leadUpT       = 0.;
  Events_refmult       = 0;
  Events_refPos        = 0;
  Events_refNeg        = 0;
  Events_runId         = 0;
  Events_MagF          = 0.;
  Events_nVertex       = 0;
  Events_rankV         = 0;
  Events_pvpdz         = 0.;
  Events_BBCz          = 0.;
  Events_ZDCz          = 0.;
  Events_primVx        = 0.;
  Events_primVy        = 0.;
  Events_primVz        = 0.;
  Events_TrigId        = 0;
  Events_TrigStat      = 0;
  Events_TrigIndex     = 0;
  Events_PartonId_NS   = 0;
  Events_PartonId_AS   = 0;
  Events_PartonStat_NS = 0;
  Events_PartonStat_AS = 0;
  Events_PartonE_NS    = 0.;
  Events_PartonE_AS    = 0.;
  Events_PartonEt_NS   = 0.;
  Events_PartonEt_AS   = 0.;
  Events_tsp           = 0.;
  Events_iso           = 0.;
  // initialize tower variables
  Events_Twr_didT     = 0;
  Events_Twr_adc11    = 0;
  Events_Twr_eneT0    = 0.;
  Events_Twr_eT       = 0.;
  Events_Twr_ENET0    = 0.;
  Events_Twr_phT      = 0.;
  Events_Twr_PTower   = 0.;
  Events_Twr_pidTower = 0.;
  Events_Twr_moduleT  = 0;
  // initialize tower cluster variables
  Events_Clust_EneT0 = 0.;
  Events_Clust_EneR  = 0.;
  Events_Clust_EneH  = 0.;
  Events_Clust_etav1 = 0.;
  Events_Clust_phiv1 = 0.;
  // initialize SMD eta strip variables
  Events_Twr_didE   = 0;
  Events_Estrp_en03 = 0.;
  Events_Estrp_en02 = 0.;
  Events_Estrp_en01 = 0.;
  Events_Estrp_en0  = 0.;
  Events_Estrp_en1  = 0.;
  Events_Estrp_en2  = 0.;
  Events_Estrp_en3  = 0.;
  Events_Estrp_en4  = 0.;
  Events_Estrp_en5  = 0.;
  Events_Estrp_en6  = 0.;
  Events_Estrp_en7  = 0.;
  Events_Estrp_en8  = 0.;
  Events_Estrp_en9  = 0.;
  Events_Estrp_en10 = 0.;
  Events_Estrp_en11 = 0.;
  Events_Estrp_en12 = 0.;
  Events_Estrp_en13 = 0.;
  Events_Estrp_en14 = 0.;
  Events_Estrp_en15 = 0.;
  // initialize SMD phi strip variables
  Events_Pstrip_enp03 = 0.;
  Events_Pstrip_enp02 = 0.;
  Events_Pstrip_enp01 = 0.;
  Events_Pstrip_enp0  = 0.;
  Events_Pstrip_enp1  = 0.;
  Events_Pstrip_enp2  = 0.;
  Events_Pstrip_enp3  = 0.;
  Events_Pstrip_enp4  = 0.;
  Events_Pstrip_enp5  = 0.;
  Events_Pstrip_enp6  = 0.;
  Events_Pstrip_enp7  = 0.;
  Events_Pstrip_enp8  = 0.;
  Events_Pstrip_enp9  = 0.;
  Events_Pstrip_enp10 = 0.;
  Events_Pstrip_enp11 = 0.;
  Events_Pstrip_enp12 = 0.;
  Events_Pstrip_enp13 = 0.;
  Events_Pstrip_enp14 = 0.;
  Events_Pstrip_enp15 = 0.;
  // initialize SMD cluster variables
  Events_clust_Ennq1  = 0.;
  Events_clust_Ennq20 = 0.;
  Events_clust_Ennq19 = 0.;
  Events_clust_Enpq1  = 0.;
  Events_clust_Enpq20 = 0.;
  Events_clust_Enpq19 = 0.;
  Events_clust_Enpq21 = 0.;
  // initialize no. of primary tracks
  Events_noOfprimaryTrks = 0;

  // intitialize primary track variables
  for (int i = 0; i < NTmax; ++i) {
    pTracks_pT[i]     = 0.;
    pTracks_px[i]     = 0.;
    pTracks_py[i]     = 0.;
    pTracks_pz[i]     = 0.;
    pTracks_Eta[i]    = 0.;
    pTracks_Phi[i]    = 0.;
    pTracks_dEdx[i]   = 0.;
    pTracks_chrg[i]   = 0.;
    pTracks_gdca[i]   = 0.;
    pTracks_Fp[i]     = 0;
    pTracks_Ppo[i]    = 0;
    pTracks_nSigPi[i] = 0.;
    pTracks_nSigK[i]  = 0.;
    pTracks_nSigP[i]  = 0.;
    pTracks_nSigE[i]  = 0.;
  }

}  // end 'Event()'

// End ------------------------------------------------------------------------
