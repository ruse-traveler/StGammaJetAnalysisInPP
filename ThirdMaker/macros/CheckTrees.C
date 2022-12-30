// 'CheckTrees.C'
// Derek Anderson
//
// Creates a chain out of trees and draws some leaves
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

// i/o parameters
static const TString sTree("tree3.d9m6y2017.list");
static const TString sBad("bad3.d9m6y2017.list");
static const TString sOut("tree3.d9m6y2017.root");


void CheckTrees() {

  cout << "\n  Checking ThirdJetMaker output..." << endl;


  // lower verbosity
  gErrorIgnoreLevel = kFatal;


  // initialize output
  TFile  *fOut   = new TFile(sOut.Data(), "recreate");
  TChain *tChain = new TChain("Gfmtodst");

  // open streams
  ifstream files(sTree.Data());
  ofstream bad(sBad.Data());
  if (!files) {
    cerr << "PANIC: input file couldn't be opened!" << endl;
    assert(files);
  }


  // merge files
  cout << "    Reading in files..." << endl;

  Int_t  t;
  string f;
  while (files) {
    files >> f;
    TString file(f);
    t = tChain -> Add(file.Data(), 0);
    if (t == 0) {
      cout << "      Bad file:  '" << f << "'..." << endl;
      bad << f;
      bad << endl;
    }
    else
      cout << "      Added file '" << f << "'..." << endl;
  }


  // check entries
  cout << "    Drawing plots..." << endl; 

  TCanvas *cRunNumber = new TCanvas("cRunNumber", "Run Number", 700, 500);
  cRunNumber -> cd();
  tChain     -> Draw("EventList.runNumber");
  cRunNumber -> Write();
  cRunNumber -> Close();

  TCanvas *cRefmult = new TCanvas("cRefmult", "Refmult", 700, 500);
  cRefmult -> cd();
  tChain   -> Draw("EventList.refMult");
  cRefmult -> Write();
  cRefmult -> Close();


  fOut   -> cd();
  tChain -> Write();
  fOut   -> Close();


  cout << "  Finished checking!\n" << endl;

}

// End ------------------------------------------------------------------------
