// 'GammaJetTrack.h'
// Derek Anderson
// 01.11.2016
//
// This file contains the definition of the Track class used in the analysis
// defined by 'GammaJetAnalysis.h' and executed by 'GammaJetExecute.cc'.  Note
// that not all quantities here will be relevant for the pythia simulation;
// they're left in for completion's sake.  ['a' for trAck]
//
// NOTE: that we're working in natural units.
//
// Last updated: 08.16.2016

#pragma once

#include "TRandom3.h"


class Track {

public:

  Track();

  // method
  void   setComponents(const double Etot, const double m=-1.);
  double setTrack(const int index, const int code, const double resTPC, const double resCal, const Pythia8::Particle a);

  Int_t   aIndex;   // evt. record index
  Int_t   aCode;    // encodes what produced track and which Delta-Phi region it's in
  Int_t   aID;      // particle PID
  Bool_t  aNear;    // indicates if track descended from NS hard-scatter product
  Bool_t  aAway;    // indicates if track descended from AS hard-scatter product
  Bool_t  aUnc;     // indicates if track didn't descend from either product
  Bool_t  aTR;      // track lies in "Trigger" Delta-Phi region
  Bool_t  aRE;      // track lies in "Recoil" Delta-phi region
  Bool_t  aUE1;     // track lies in upper "Underlying Event" Delta-Phi region
  Bool_t  aUE2;     // track lies in lower "Underlying Event" Delta-Phi region
  Bool_t  aDet;     // indicates if track is particle or detector level
  Float_t aPt;
  Float_t aPx;
  Float_t aPy;
  Float_t aPz;
  Float_t aEne;     // total energy
  Float_t aEneT;    // transverse energy
  Float_t aEta;
  Float_t aPhi;
  Float_t aDphi;
  Float_t aChrg;
  Float_t aMass;
  Float_t aSigPi;   // nonzero if track is pion
  Float_t aSigK;    // nonzero if track is kaon
  Float_t aSigP;    // nonzero if track is proton
  Float_t aSigE;    // nonzero if track is exotic(?)

};


// member defintions ----------------------------------------------------------

Track::Track() {

  aIndex   = 0;
  aCode    = 0;
  aID      = 0;
  aNear    = false;
  aAway    = false;
  aUnc     = false;
  aTR      = false;
  aRE      = false;
  aUE1     = false;
  aUE2     = false;
  aDet     = false;
  aPt      = 0.;
  aPx      = 0.;
  aPy      = 0.;
  aPz      = 0.;
  aEne     = 0.;
  aEneT    = 0.;
  aEta     = 0.;
  aPhi     = 0.;
  aDphi    = 0.;
  aChrg    = 0.;
  aMass    = 0.;
  aSigPi   = 0.;
  aSigK    = 0.;
  aSigP    = 0.;
  aSigE    = 0.;

}  // end 'Track()'



void Track::setComponents(const double Etot, const double m) {

  double Ptot  = 0.;
  double theta = 0.;

  if (m != -1.) {
    Ptot  = sqrt(Etot*Etot - m*m);
    aMass = m;
  }
  else
    Ptot  = sqrt(Etot*Etot - aMass*aMass);
  theta = 2.*exp(-1.*aEta);
  aEne  = Etot;
  aEneT = Etot * sin(theta);
  aPz   = Ptot * cos(theta);
  aPt   = Ptot * sin(theta);
  aPx   = aPt  * cos(aPhi);
  aPy   = aPt  * sin(aPhi);

}  // end 'setComponents(double)'


double Track::setTrack(const int index, const int code, const double resTPC, const double resCal, const Pythia8::Particle a) {

  bool   accept     = true;
  double efficiency = 1.;
  double passEff    = 0.;
  double pTs        = 0.;
  double pXs        = 0.;
  double pYs        = 0.;
  double pZs        = 0.;
  double Ptpc       = 0.;
  double Etpc       = 0.;
  double rando1     = 0.;
  double rando2     = 0.;
  double Edep       = 0.;
  double Pdep       = 0.;
  double dE         = 0.;

  // determine what produced track and its Delta-Phi region
  switch (code) {
    case 0:
      aNear = true;
      aTR   = true;
      break;
    case 1:
      aNear = true;
      aRE   = true;
      break;
    case 2:
      aNear = true;
      aUE1  = true;
      break;
    case 3:
      aNear = true;
      aUE2  = true;
      break;
    case 10:
      aAway = true;
      aTR   = true;
      break;
    case 11:
      aAway = true;
      aRE   = true;
      break;
    case 12:
      aAway = true;
      aUE1  = true;
      break;
    case 13:
      aAway = true;
      aUE2  = true;
      break;
    case 20:
      aUnc  = true;
      aTR   = true;
      break;
    case 21:
      aUnc  = true;
      aRE   = true;
      break;
    case 22:
      aUnc  = true;
      aUE1  = true;
      break;
    case 23:
      aUnc  = true;
      aUE2  = true;
      break;
  }

  // set particle-lvl. properties
  aIndex = index;
  aCode  = code;
  aID    = a.id();
  aPt    = a.pT();
  aPx    = a.px();
  aPy    = a.py();
  aPz    = a.pz();
  aEne   = a.e();
  aEneT  = a.eT();
  aEta   = a.eta();
  aPhi   = a.phi();
  aChrg  = a.charge();
  aMass  = a.m();


  // if detector level, smear relevant quantities
  if (aDet) {

    // TPC response
    if (aChrg != 0.) {
      pTs        = a.pT()*(1.+resTPC);
      pXs        = a.px()*(1.+resTPC);
      pYs        = a.py()*(1.+resTPC);
      pZs        = a.pz()*(1.+resTPC);
      Ptpc       = sqrt(pXs*pXs + pYs*pYs + pZs*pZs);
      Etpc       = sqrt(pXs*pXs + pYs*pYs + pZs*pZs + aMass*aMass);
      efficiency = 0.91 - 0.91 * exp(-4.0*pTs);
      passEff    = gRandom -> Uniform(0., 1.);
      if (passEff > efficiency) accept = false;
    }


    if (accept) {
      // hadron cal. response
      rando1 = gRandom -> Uniform(0.,1.);
      if (a.isHadron() && (rando1 < 0.35)) {
        rando2 = gRandom -> Uniform(0.,1.);
        // MIP
        if (rando2 < 0.6) {
          Edep = gRandom -> Gaus(0.3, 0.08);
          if (Edep < 0.) Edep = 0.;
        }
        // NI
        else {
          Edep = aEne * (gRandom -> Gaus(0.4, 0.2));
          if (Edep < 0.) Edep = 0.;
        }
      }

      // muon cal. response
      if (abs(aID) == 13) {
        Edep = gRandom -> Gaus(0.3, 0.08);
        if (Edep < 0.) Edep = 0.;
      }

      // electron / photon cal. response
      if ((abs(aID) == 11) || (abs(aID) == 22)) {
        Edep = 0.95 * aEne * (gRandom -> Gaus(1.0,resCal));
        if (Edep < 0.) Edep = 0.;
      }


      // store charged detector-lvl. properties
      if (aChrg != 0.) {
        dE = Edep - Ptpc;
        if (dE < 0.2) dE = 0.;
        aEne  = Etpc;
        aEneT = Etpc * sin(2.*atan(exp(-1.*aEta)));
        aPt   = pTs;
        aPx   = pXs;
        aPy   = pYs;
        aPz   = pZs;
      }
      // store neutral detector-lvl. properties
      if (aChrg == 0.) {
        Pdep  = sqrt(Edep*Edep - aMass*aMass);
        aEne  = Edep;
        aEneT = Edep * sin(2.*atan(exp(-1.*aEta)));
        aPt   = Pdep * sin(2.*atan(exp(-1.*aEta)));
        aPx   = aPt  * cos(aPhi);
        aPy   = aPt  * sin(aPhi);
      }

    }  // end 'if(accept)'
    aPt  = pTs;
    aPx  = pXs;
    aPy  = pYs;
    aPz  = pZs;
    aEne = Etpc;

  }  // end smearing


  if (!accept)
    return -1000.;
  else
    return dE;

}  // end 'setTrack(int, double, double, Particle)'

// End ------------------------------------------------------------------------
