// 'StDataQaMaker.cxx'
// Derek Anderson
// 09.19.2019
//
// This class produces the data QA
// plots for the neutral-triggered
// pp recoil jet analysis note.


#define StDataQaMaker_cxx

#include "StDataQaMaker.h"
#include "StDataQaMaker.eventQaPlots.h"
#include "StDataQaMaker.jetQaPlots.h"

using namespace std;

ClassImp(StDataQaMaker);



void StDataQaMaker::Init(const TString sOut, const TString sIn, const TString sTree) {

  sOutFile = sOut;
  sInFile  = sIn;
  sInTree  = sTree;
  cout << "    Initialized 'StDataQaMaker' object:\n"
       << "      Output file = " << sOutFile.Data() << "\n"
       << "      Input file  = " << sInFile.Data() << "\n"
       << "      Input tree  = " << sInTree.Data()
       << endl;

}  // end 'Init(TString, TString)'



void StDataQaMaker::Make(const Bool_t doEvents, const Bool_t doJets) {

  doEventQaPlots = doEvents;
  doJetQaPlots   = doJets;
  if (doEventQaPlots) MakeEventQaPlots();
  if (doJetQaPlots)   MakeJetQaPlots();
  cout << "    Finished!" << endl;

}  // end 'Make(Bool_t, Bool_t)'



void StDataQaMaker::Finish() {

  cout << "  Deleting 'StDataQaMaker' and exiting!\n" << endl;

}  // end 'Finish()'

// End ------------------------------------------------------------------------
