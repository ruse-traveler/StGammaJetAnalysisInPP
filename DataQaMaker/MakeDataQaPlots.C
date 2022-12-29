// 'MakeDataQaPlots.C'
// Derek Anderson
// 09.24.2019
//
// This produces the (data) event
// and jet qa plots.

#include "TString.h"
#include "TSystem.h"

using namespace std;

class StDataQaMaker;

// i/o parameters
static const TString sOutput("jetQaPlots_fixedRanges.et920vz55pt0230.d26m9y2021.root");
static const TString sInput("../FullJetTree/merged/pp200r9.merge.root");
static const TString sTree("Gfmtodst");
static const Bool_t  doEventPlots(false);
static const Bool_t  doJetPlots(true);



void MakeDataQaPlots(const Bool_t isInBatchMode=false) {

  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem -> Load("/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem -> Load("StDataQaMaker");

  StDataQaMaker *dataQaMaker = new StDataQaMaker(isInBatchMode);
  dataQaMaker -> Init(sOutput.Data(), sInput.Data(), sTree.Data());
  dataQaMaker -> Make(doEventPlots, doJetPlots);
  dataQaMaker -> Finish();

}

// End ------------------------------------------------------------------------
