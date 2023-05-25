// 'GammaJetSystem.h'
// Derek Anderson
// 05.01.2016
//
// This file contains the definition of the System class.  This collects
// miscellaneous variables and methods which aid the analysis defined in
// 'GammaJetAnalysis.h'.
//
// NOTE: "partons" as used here refer to the products of the hard-scatter
//       which initiated the event unless otherwise noted.  In addition,
//       DeltaPhi is calculated wrt. trigger. 
//
// Last updated: 05.13.2016

#pragma once



class System {

public:

  System();

  // methods
  void reset();
  bool isNewNtrk(const int iTrack);
  bool isNewAtrk(const int iTrack);
  bool isNewUtrk(const int iTrack);
  bool hasNmom(const Pythia8::Particle& p);
  bool hasAmom(const Pythia8::Particle& p);

  // members
  int   eID;          // event ID = trigger ID + no. direct photons (nG)
  int   uID;          // indicates if leading track in UE / AR is > 3 GeV
  int   nTrg;         // no. of triggers found so far
  int   nTR;          // no. of particles inside "trigger" region
  int   nRE;          // no. of particles inside "recoil" region
  int   nAR1;         // no. of particles inside upper "underlying event" region
  int   nAR2;         // no. of particles inside lower "underlying event" region
  int   nIndex;       // index of near-side parton
  int   aIndex;       // index of away-side parton
  int   tIndex;       // index of trigger
  int   t2Index;      // index of secondary trigger
  int   leadUIndex;   // index of leading track in UE / AR
  int   mindex;       // index of parton whose daughter had smallest DeltaPhi
  int   maxdex;       // index of parton whose daughter had largest DeltaPhi
  int   refP;         // particle level multiplicities
  int   posP;
  int   negP;
  int   refD;         // detector level multiplicities
  int   posD;
  int   negD;
  float minDphi;      // smallest daughter DeltaPhi
  float maxDphi;      // largest daughter DeltaPhi
  float t2Phi;        // phi of secondary trigger
  float leadUPhi;     // phi of leading track in UE / AR
  float pTtr;         // total pT in "trigger" DeltaPhi region
  float pTre;         // total pT in "recoil" DeltaPhi region
  float pTar1;        // total pT in upper "underlying event" DeltaPhi region
  float pTar2;        // total pT in lower "underlying event" DeltaPHi region
  float t2Pt;         // pT of secondary trigger
  float leadUPt;      // pT of leading track in UE / AR
  vector<int> prdx;   // indices of partons
  vector<int> nTrks;  // indices of all particles spawned by the near-side parton
  vector<int> aTrks;  // indices of all particles spawned by the away-side parton
  vector<int> uTrks;  // indices of all the rest of the particles

};


// member definitions ---------------------------------------------------------

System::System() {

  eID        = 0;
  uID    = 0;
  nTrg       = 0;
  nTR        = 0;
  nRE        = 0;
  nAR1       = 0;
  nAR2       = 0;
  nIndex     = 0;
  aIndex     = 0;
  tIndex     = 0;
  t2Index    = 0;
  leadUIndex = 0;
  mindex     = 0;
  maxdex     = 0;
  refP       = 0;
  posP       = 0;
  negP       = 0;
  refD       = 0;
  posD       = 0;
  negD       = 0;
  minDphi    = 0.;
  maxDphi    = 0.;
  t2Phi      = 0.;
  leadUPhi   = 0.;
  pTtr       = 0.;
  pTre       = 0.;
  pTar1      = 0.;
  pTar2      = 0.;
  t2Pt       = 0.;
  leadUPt    = 0.;
  prdx.clear();
  nTrks.clear();
  aTrks.clear();
  uTrks.clear();

}  // end 'System()'


void System::reset() {

  eID        = 0;
  uID        = 0;
  nTR        = 0;
  nRE        = 0;
  nTR        = 0;
  nAR1       = 0;
  nAR2       = 0;
  nIndex     = 0;
  aIndex     = 0;
  tIndex     = 0;
  t2Index    = 0;
  leadUIndex = 0;
  mindex     = 0;
  maxdex     = 0;
  refP       = 0;
  posP       = 0;
  negP       = 0;
  refD       = 0;
  posD       = 0;
  negD       = 0;
  minDphi    = 0.;
  maxDphi    = 0.;
  t2Phi      = 0.;
  leadUPhi   = 0.;
  pTtr       = 0.;
  pTre       = 0.;
  pTar1      = 0.;
  pTar2      = 0.;
  t2Pt       = 0.;
  leadUPt    = 0.;
  prdx.clear();
  nTrks.clear();
  aTrks.clear();
  uTrks.clear();

}  // end 'reset()'


bool System::isNewNtrk(const int iTrack) {

  int  Nn;
  bool isNew;

  Nn    = nTrks.size();
  isNew = true;
  for (int i = 0; i < Nn; ++i) {
    if (iTrack == nTrks[i]) {
      isNew = false;
      break;
    }
  }

  return isNew;

}  // end 'isNewNtrk(int)'


bool System::isNewAtrk(const int iTrack) {

  int  Na;
  bool isNew;

  Na    = aTrks.size();
  isNew = true;
  for (int i = 0; i < Na; ++i) {
    if (iTrack == aTrks[i]) {
      isNew = false;
      break;
    }
  }

  return isNew;

}  // end 'isNewAtrk(int)'


bool System::isNewUtrk(const int iTrack) {

  int  Nu;
  bool isNew;

  Nu    = uTrks.size();
  isNew = true;
  for (int i = 0; i < Nu; ++i) {
    if (iTrack == uTrks[i]) {
      isNew = false;
      break;
    }
  }

  return isNew;

}  // end 'isNewUtrk(int)'


bool System::hasNmom(const Pythia8::Particle& p) {

  int  Nn;
  bool Nmom;

  Nn   = nTrks.size();
  Nmom = false;
  for (int i = 0; i < Nn; ++i) {
    Nmom = p.isAncestor(nTrks[i]);
    if (Nmom) break;
  }

  return Nmom;

}  // end 'hasNmom(Particle)'


bool System::hasAmom(const Pythia8::Particle& p) {

  int  Na;
  bool Amom;

  Na    = aTrks.size();
  Amom  = false;
  for (int i = 0; i < Na; ++i) {
    Amom = p.isAncestor(aTrks[i]);
    if (Amom) break;
  }

  return Amom;

}  // end 'hasAmom(Particle)'

// End ------------------------------------------------------------------------
