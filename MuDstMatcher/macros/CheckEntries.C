// 'CheckEntries.C'
// Derek Anderson
//
// Creates a chain out of matching-output and draws some leaves
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


void CheckEntries() {

  cout << "\nChecking output..." << endl;


  // lower verbosity
  gErrorIgnoreLevel = kFatal;


  // initialize output
  TFile  *fOut   = new TFile("./pt9.checkWithMc.root", "recreate");
  TChain *gChain = new TChain("GfmtoDst_gnt");
  TChain *uChain = new TChain("GfmtoDst_mu");
  TChain *tChain = new TChain("McTracks");

  // open streams
  ifstream files("pt9.checkWithMc.list");
  ofstream bad("pt9.bad.list");
  if (!files) {
    cerr << "PANIC: input file couldn't be opened!" << endl;
    assert(files);
  }


  // merge files
  cout << "  Reading in files..." << endl;

  Int_t  g;
  Int_t  u;
  Int_t  t;
  string f;
  while (files) {
    files >> f;
    TString file(f);
    g = gChain -> Add(file.Data(), 0);
    u = uChain -> Add(file.Data(), 0);
    t = tChain -> Add(file.Data(), 0);
    if ((g == 0) || (u == 0) || (t == 0)) {
      cout << "    Bad file:  '" << f << "'..." << endl;
      bad << f;
      bad << endl;
    }
    else
      cout << "    Added file '" << f << "'..." << endl;
  }


  // check entries
  cout << "  Drawing plots..." << endl; 

  TCanvas *cGntRunNumber = new TCanvas("cGntRunNumber", "Geant Run Number", 700, 500);
  cGntRunNumber -> cd();
  gChain        -> Draw("GntEventList.runNumber");
  cGntRunNumber -> Write();
  cGntRunNumber -> Close();

  TCanvas *cGntRefmult = new TCanvas("cGntRefmult", "Geant Refmult", 700, 500);
  cGntRefmult -> cd();
  gChain      -> Draw("GntEventList.refMult");
  cGntRefmult -> Write();
  cGntRefmult -> Close();

  TCanvas *cMuRunNumber = new TCanvas("cMuRunNumber", "MuDstt Run Number", 700, 500);
  cMuRunNumber -> cd();
  uChain        -> Draw("MuEventList.runNumber");
  cMuRunNumber -> Write();
  cMuRunNumber -> Close();

  TCanvas *cMuRefmult = new TCanvas("cMuRefmult", "MuDst Refmult", 700, 500);
  cMuRefmult -> cd();
  uChain      -> Draw("MuEventList.refMult");
  cMuRefmult -> Write();
  cMuRefmult -> Close();


  fOut   -> cd();
  gChain -> Write();
  uChain -> Write();
  fOut   -> Close();


  cout << "Finished checking!\n" << endl;

}

// End ------------------------------------------------------------------------
