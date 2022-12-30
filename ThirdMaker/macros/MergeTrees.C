// 'MergeTrees.C'
// Derek Anderson
// 02.21.2017
//
// Use this to merge the output from 'StThirdMaker'
// and 'StThirdJetMaker' into TChains.  Run
// 'MakeTreeList.sh' first.

#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TString.h"

using namespace std;


// i/o constants
static const Int_t   nTuples(17);
static const TString sTuple("tuple.r10145070.d14m4y2017.list");
static const TString sTree("test.r10145070.d14m4y2017.list");
static const TString sOut("merge.after.r10145070.d14m4y2017.root");



void MergeTrees() {

  cout << "\n  Beginning merging script..." << endl;
  gErrorIgnoreLevel = kError;


  TFile  *fOut = new TFile(sOut.Data(), "recreate");
  TChain *tTup[nTuples];
  TChain *tTree;

  tTup[0]  = new TChain("QAAB", "QAAB");
  tTup[1]  = new TChain("AngleTAB", "AngleTAB");
  tTup[2]  = new TChain("CAB", "CAB");
  tTup[3]  = new TChain("SMDE1B", "SMDE1B");
  tTup[4]  = new TChain("PiB", "PiB");
  tTup[5]  = new TChain("TrackS", "TrackS");
  tTup[6]  = new TChain("TrackSS", "TrackSS");
  tTup[7]  = new TChain("Track", "Track");
  tTup[8]  = new TChain("QAA", "QAA");
  tTup[9]  = new TChain("Pi", "Pi");
  tTup[10] = new TChain("SMDE1", "SMDE1");
  tTup[11] = new TChain("SMDE1B", "SMDE1B");
  tTup[12] = new TChain("AngleTA", "AngleTA");
  tTup[13] = new TChain("CA", "CA");
  tTup[14] = new TChain("New8A", "New8A");
  tTup[15] = new TChain("FiveLL", "FiveLL");
  tTup[16] = new TChain("Track_M", "Track_M");
  tTree    = new TChain("Gfmtodst", "Gfmtodst");


  cout << "    Streaming filepaths..." << endl;

  ifstream tuples(sTuple.Data());
  ifstream trees(sTree.Data());
  if (!tuples) {
    cerr << "PANIC: couldn't open tuple stream!" << endl;
    assert(tuples);
  }
  if (!trees) {
    cerr << "PANIC: couldn't open tree stream!" << endl;
    assert(trees);
  }

  Int_t  nU(0);
  Int_t  nR(0);
  string uPath;
  string rPath;
  while (getline(tuples, uPath)) {
    TString uFile(uPath);
    for (Int_t i = 0; i < nTuples; i++) {
      tTup[i] -> Add(uFile.Data(), 0);
    }
    nU++;
  }
  while (getline(trees, rPath)) {
    TString rFile(rPath);
    tTree -> Add(rFile.Data(), 0);
    nR++;
  }
  cout << "    Files streamed: " << nU << " tuple files, "
       << nR << " tree files..."
       << endl;

  tuples.close();
  trees.close();


  cout << "    Writing to file..." << endl;

  fOut  -> cd();
  for (Int_t i = 0; i < nTuples; i++) {
    tTup[i] -> Write();
  }
  tTree -> Write();
  fOut  -> Close();


  cout << "  Merging script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
