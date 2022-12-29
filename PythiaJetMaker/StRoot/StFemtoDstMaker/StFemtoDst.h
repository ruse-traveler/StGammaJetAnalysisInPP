// 'StFemtoDst.h'
// Derek Anderson, Nihar Sahoo
// 08.04.2016
//
// This is the definition of the 'FemtoDST' object which is
// filled by the 'StFemtoDstMaker' class.  This is a compact
// tree for jet-related analyses.
//
// Last updated: 08.05.2016

#ifndef StFemtoDst_h
#define StFemtoDst_h

#include <vector>
#include "TObject.h"

using namespace std;


class StFemtoDst : public TObject {

public:

  // event variables
  Int_t    EventIndex;
  Double_t Refmult;
  Double_t TSP;
  Double_t TrgEta;
  Double_t TrgPhi;
  Double_t TrgEt;
  Double_t Rho;
  Double_t Sigma;
  Double_t Vz;

  // jet variables
  vector<Double_t> JetPt;
  vector<Double_t> JetNCons;
  vector<Double_t> JetIndex;
  vector<Double_t> JetEta;
  vector<Double_t> JetPhi;
  vector<Double_t> JetE;
  vector<Double_t> JetArea;

  // constituent variables
  vector<vector<Double_t>> JetConsPt;
  vector<vector<Double_t>> JetConsEta;
  vector<vector<Double_t>> JetConsPhi;
  vector<vector<Double_t>> JetConsE;


  // methods
  StFemtoDst()  {}
  ~StFemtoDst() {}
  void Reset() {
    EventIndex = 0;
    Refmult    = 0.;
    TSP        = 0.;
    TrgEta     = 0.;
    TrgPhi     = 0.;
    TrgEt      = 0.;
    Rho        = 0.;
    Sigma      = 0.;
    Vz         = 0.;
    JetPt.clear();
    JetNCons.clear();
    JetIndex.clear();
    JetEta.clear();
    JetPhi.clear();
    JetE.clear();
    JetArea.clear();
    JetConsPt.clear();
    JetConsEta.clear();
    JetConsPhi.clear();
    JetConsE.clear();
  }


  ClassDef(StFemtoDst, 1)

};


#endif

// End ------------------------------------------------------------------------
