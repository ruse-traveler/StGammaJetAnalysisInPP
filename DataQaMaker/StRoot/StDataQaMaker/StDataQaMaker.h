// 'StDataQaMaker.h'
// Derek Anderson
// 09.19.2019
//
// This class produces the data QA
// plots for the neutral-triggered
// pp recoil jet analysis note.

#ifndef StDataQaMaker_h
#define StDataQaMaker_h

#include <vector>
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TBox.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

using namespace std;

// global constants
const UInt_t NTrkMax(5000);
const UInt_t NTwrMax(5000);
const UInt_t NMatchMax(10);
const UInt_t NHotTwr(302);
const UInt_t NBadRuns(45);

// event constants
const UInt_t NTrgCuts(9);
const UInt_t NTrgBins(4);
const UInt_t NTrkBins(4);
const UInt_t NTrgTsp(2);
const UInt_t NVtx(4);

// jet constants
const UInt_t NJetTypes(2);
const UInt_t NJetBins(4);
const UInt_t NResVals(2);



class StDataQaMaker {

  private:

    TString sOutFile;
    TString sInFile;
    TString sInTree;
    Bool_t  isInBatchMode;
    Bool_t  doEventQaPlots;
    Bool_t  doJetQaPlots;

    // private methods
    void MakeEventQaPlots();
    void MakeJetQaPlots();

  public:

    // ctor, dtor
    StDataQaMaker(const Bool_t batch=0);
    virtual ~StDataQaMaker();

    // public methods
    void Init(const TString sOut="", const TString sIn="", const TString sTree="");
    void Make(const Bool_t doEvents=0, const Bool_t doJets=0);
    void Finish();

  ClassDef(StDataQaMaker, 1);

};

#endif
#ifdef StDataQaMaker_cxx



StDataQaMaker::StDataQaMaker(const Bool_t batch) {

  isInBatchMode  = batch;
  cout << "\n  Created 'StDataQaMaker' object." << endl;

}  // end ctor



StDataQaMaker::~StDataQaMaker() {}  // end dtor

#endif

// End ------------------------------------------------------------------------
