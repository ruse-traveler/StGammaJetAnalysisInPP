// 'CheckTuples.C'
// Derek Anderson
//
// Creates a chain out of specified tuples and draws some leaves
// to check the output.

#include <cassert>
#include <fstream>
#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

// i/o constants
static const Int_t   nTuples(16);
static const TString sTuple("old.d4m5y2017.list");
static const TString sBad("badder.d4m5y2017.list");
static const TString sOut("old.d4m5y2017.root"); 


void CheckTuples() {

  cout << "\n  Checking ThirdMaker output..." << endl;


  // lower verbosity
  gErrorIgnoreLevel = kFatal;


  // initialize output
  TFile  *fOut = new TFile(sOut.Data(), "recreate");
  TChain *tChain[nTuples];

  tChain[0]  = new TChain("QAAB", "QAAB");
  tChain[1]  = new TChain("AngleTAB", "AngleTAB");
  tChain[2]  = new TChain("CAB", "CAB");
  tChain[3]  = new TChain("SMDE1B", "SMDE1B");
  tChain[4]  = new TChain("PiB", "PiB");
  tChain[5]  = new TChain("TrackS", "TrackS");
  tChain[6]  = new TChain("TrackSS", "TrackSS");
  tChain[7]  = new TChain("Track", "Track");
  tChain[8]  = new TChain("QAA", "QAA");
  tChain[9]  = new TChain("Pi", "Pi");
  tChain[10] = new TChain("SMDE1", "SMDE1");
  tChain[11] = new TChain("SMDE1B", "SMDE1B");
  tChain[12] = new TChain("AngleTA", "AngleTA");
  tChain[13] = new TChain("CA", "CA");
  tChain[14] = new TChain("New8A", "New8A");
  //tChain[15] = new TChain("Track_M", "Track_M");
  tChain[15] = new TChain("FiveLL", "FiveLL");

  // open streams
  ifstream files(sTuple.Data());
  ofstream bad(sBad.Data());
  if (!files) {
    cerr << "PANIC: input file couldn't be opened!" << endl;
    assert(files);
  }


  // merge files
  cout << "    Reading in files..." << endl;

  Int_t  g;
  Int_t  t;
  string f;
  while (files) {
    files >> f;
    TString file(f);

    g = 0;
    t = 0;
    for (Int_t i = 0; i < nTuples; i++) {
      t = tChain[i] -> Add(file.Data(), 0);
      if (t != 0) g++;
    }

    if (g != nTuples) {
      cout << "      Bad file:  '" << f << "'..." << endl;
      bad << f;
      bad << endl;
    }
    else {
      cout << "      Added file '" << f << "'..." << endl;
    }
  }  // end while loop


  // check entries
  cout << "    Drawing plots..." << endl; 

  TCanvas *cRunNumber = new TCanvas("cRunNumber", "Run Number", 700, 500);
  cRunNumber -> cd();
  tChain[1]  -> Draw("runId");
  cRunNumber -> Write();
  cRunNumber -> Close();

  TCanvas *cEtrg = new TCanvas("cEtrg", "Trigger Energy", 700, 500);
  cEtrg     -> cd();
  tChain[0] -> Draw("EneT0");
  cEtrg     -> Write();
  cEtrg     -> Close();


  fOut -> cd();
  for (Int_t i = 0; i < nTuples; i++) {
    tChain[i] -> Write();
  }
  fOut -> Close();


  cout << "  Finished checking!\n" << endl;

}

// End ------------------------------------------------------------------------
