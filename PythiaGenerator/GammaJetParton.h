// 'GammaJetParton.h'
// Derek Anderson
// 02.01.2016
//
// This file contains the definition of the Parton class (here meaning the
// products of the hard subprocess which generated the event) used in the
// analysis defined by 'GammaJetAnalysis.h' and executed by
// 'GammaJetExecute.h'.  ['p' for Parton.]
//
// Last updated: 05.12.2016

#pragma once


class Parton {

public:

  Parton();

  // methods
  void clear();
  void setParton(const int index, const Pythia8::Particle p);

  Int_t   pID;
  Int_t   pStat;   // production code
  Int_t   pIndex;  // event record index
  Bool_t  pFinal;  // flags if parton is last parton before hadronization
  Float_t pEta;
  Float_t pPhi;
  Float_t pEne;    // total energy
  Float_t pEneT;   // transverse energy

};


// member definitions ---------------------------------------------------------

Parton::Parton() {

  pID    = 0;
  pStat  = 0;
  pIndex = 0;
  pFinal = false;
  pEta   = 0.;
  pPhi   = 0.;
  pEne   = 0.;
  pEneT  = 0.;

}  // end 'Parton()'


void Parton::clear() {

  pID    = 0;
  pStat  = 0;
  pIndex = 0;
  pFinal = false;
  pEta   = 0.;
  pPhi   = 0.;
  pEne   = 0.;
  pEneT  = 0.;

}  // end 'Clear()'


void Parton::setParton(const int index, const Pythia8::Particle p) {

  pIndex = index;
  pID    = p.id();
  pStat  = p.status();
  pEta   = p.eta();
  pPhi   = p.phi();
  pEne   = p.e();
  pEneT  = p.eT();

}  // end 'setParton(Particle)'

// End ------------------------------------------------------------------------
