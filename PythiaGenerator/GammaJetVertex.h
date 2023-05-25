// 'GammaJetVertex.h'
// Derek Anderson
// 01.14.2016
//
// This file contains the defintion of the Vertex class used by the analysis
// defined by 'GammaJetAnalysis.h' and executed by 'GammaJetExecute.cc'.
// This class holds information about all the vertices in an event, such as
// no. of vertices, which is the primary, etc.  Not all the quantities in here
// are relevant for the pythia simulation; they're kept in just for completion's
// sake.
//
// Last updated: 02.01.2016

#pragma once


class Vertex {

public:

  Vertex();

  // methods
  bool isNewVtx(const Pythia8::Vec4 v);
  void addVtx(const Pythia8::Vec4 v);
  void asscVtx(const bool vtx, const Pythia8::Vec4 v);
  void setPrimeVtx();
  void clear();

  // data members
  Int_t   nVertex;
  Float_t primVx;
  Float_t primVy;
  Float_t primVz;
  Float_t primVt;
  vector<Pythia8::Vec4> V;  // list of vertices
  vector<int> N;            // no. of tracks pointing to ith vtx.

};


// member definitions ---------------------------------------------------------

Vertex::Vertex() {

  nVertex = 0;
  primVx  = 0.;
  primVy  = 0.;
  primVz  = 0.;
  primVt  = 0.;
  V.clear();
  N.clear();

}  // end 'Vertex()'


bool Vertex::isNewVtx(const Pythia8::Vec4 v) {

  int  nDiff, nOld;
  bool xSame, ySame, zSame, tSame;
  bool isNew;

  isNew = false;
  nDiff = 0;
  nOld  = V.size();
  for (int i = 0; i < nOld; ++i) {
    xSame = (V[i].px() == v.px());
    ySame = (V[i].py() == v.py());
    zSame = (V[i].pz() == v.pz());
    tSame = (V[i].e()  == v.e()); 
    if (xSame && ySame && zSame && tSame) continue;
    ++nDiff;
  }
  if (nDiff == nOld) isNew = true;

  return isNew;

}  // end 'isNewVtx(Vec4)'


void Vertex::addVtx(const Pythia8::Vec4 v) {

  int n = 0;
  V.push_back(v);
  N.push_back(n);
  ++nVertex;

}  // end 'addVtx(Vec4)'


void Vertex::asscVtx(const bool vtx, const Pythia8::Vec4 v) {

  int  nVtx;
  bool xSame, ySame, zSame, tSame;

  nVtx = V.size();
  for (int i = 0; i < nVtx; ++i) {

    // if track points to origin
    if (!vtx) {
      xSame = (V[i].px() == 0.);
      ySame = (V[i].py() == 0.);
      zSame = (V[i].pz() == 0.);
      tSame = (V[i].e()  == 0.);
      if (xSame && ySame && zSame && tSame) ++N[i];
    }
    // if not origin, determine which vtx
    else {
      xSame = (V[i].px() == v.px());
      ySame = (V[i].py() == v.py());
      zSame = (V[i].pz() == v.pz());
      tSame = (V[i].e()  == v.e());
      if (xSame && ySame && zSame && tSame) ++N[i];
    }

  }

}  // end 'asscVtx(bool, Vec4)'


void Vertex::setPrimeVtx() {

  int nBig, nVtx;

  nBig = 0;
  nVtx = V.size();
  for (int i = 0; i < nVtx; ++i) {

    if (N[i] > nBig) {
      nBig   = N[i];
      primVx = V[i].px();
      primVy = V[i].py();
      primVz = V[i].pz();
      primVt = V[i].e();
    }

  }

}  // end 'setPrimeVtx'


void Vertex::clear() {

  nVertex = 0;
  V.clear();
  N.clear();

}  // end 'clear()'

// End ------------------------------------------------------------------------
