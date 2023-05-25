// 'GammaJetTrigger.h'
// Derek Anderson
// 01.08.2016
//
// This file contains the defintion of the Trigger class used in the analysis
// defined by 'GammaJetAnalysis.h' and executed by 'GammaJetExecute.h'. ['i'
// for trIgger.]
//
// Last updated: 06.09.2016

#pragma once


class Trigger {

public:

  Trigger();

  // methods
  void  clear();
  void  setIso();
  void  setEneTr();
  void  addEne(const Pythia8::Particle p);
  void  setTrigger(const int index, const Pythia8::Particle i);
  float deltaR(const Pythia8::Particle p);

  Int_t   iID;
  Int_t   iStat;    // production code
  Int_t   iIndex;   // evt. record index
  Bool_t  iPrompt;  // true if trig directly related to hard scatter prod.
  Float_t iRap; 
  Float_t iEta;
  Float_t iPhi;
  Float_t iPt;
  Float_t iEne;     // total energy
  Float_t iEneT;    // transverse energy
  Float_t iEneR;    // total energy within radius rMax
  Float_t iEneH;    // total energy of hadrons within rMax
  Float_t iso;      // ratio of EneH to EneR

};


// member definitions ---------------------------------------------------------

Trigger::Trigger() {

  iID     = 0;
  iStat   = 0;
  iIndex  = 0;
  iPrompt = false;
  iRap    = 0.;
  iEta    = 0.;
  iPhi    = 0.;
  iPt     = 0.;
  iEne    = 0.;
  iEneT   = 0.;
  iEneR   = 0.;
  iEneH   = 0.;
  iso     = 0.;

}  // end 'Trigger()'


void Trigger::clear() {

  iEneR = 0.;
  iEneH = 0.;
  iso   = 0.;

}  // end 'clear()'


void Trigger::setEneTr() {

  float tanArg, sinArg;

  tanArg = exp(-1 * iEta);
  sinArg = 2 * atan(tanArg);
  iEneT  = iEne * sin(sinArg);

}  // end 'setEneTr()'


void Trigger::setIso() {

  if (iEneR == 0.)
    iso = 0.;
  else
    iso = iEneH / iEneR;

}  // end 'setIso()'


void Trigger::addEne(const Pythia8::Particle p) {

  iEneR += p.e();
  if (p.isHadron()) iEneH += p.e();

}  // end 'calcEne(Particle)'


void Trigger::setTrigger(const int index, const Pythia8::Particle i) {

  iID    = i.id();
  iStat  = i.status();
  iIndex = index;
  iRap   = i.y();
  iEta   = i.eta();
  iPhi   = i.phi();
  iPt    = i.pT();
  iEne   = i.e();
  iEneT  = i.eT();

}  // end 'setTrigger(Particle)'


float Trigger::deltaR(const Pythia8::Particle p) {

  float Df, Df2;
  float Dh, Dh2;
  float Dr;

  Dh  = p.eta() - iEta;
  Dh2 = Dh * Dh;
  Df  = p.phi() - iPhi;
  Df2 = Df * Df;
  Dr  = sqrt(Dh2 + Df2);

  return Dr;

}  // end 'deltaR(Particle)'

// End ------------------------------------------------------------------------
