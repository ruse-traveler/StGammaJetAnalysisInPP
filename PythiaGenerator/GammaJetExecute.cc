// 'GammaJetExecute.cc'
// Derek Anderson
// 01.08.2016
//
// This exucutes the analysis defined in 'GammaJetAnalysis.h'.  All input
// settings for Pythia are contained in 'GammaJetExecute.cmnd'.  The
// functions 'timeBegin()', 'timeEnd()', and 'progress(int, int)' are defined
// in 'progress.h'.
//
// NOTE: the four "DeltaPhi Regions" are defined as so:
//
//   "Trigger" region (TR) -- -pi/4 < DeltaPhi < pi/4
//   "Recoil" region (RE)  -- -pi/4 < DeltaPhi - pi < pi/4
//   "Away" region 1 (AR1) -- pi/4 < DeltaPhi < 3pi/4
//   "Away" region 2 (AR2) -- pi/4 < DeltaPhi - pi < 3pi/4
//
// Last updated: 08.30.2016

#include "Pythia8/Pythia.h"
// c libraries
#include <ctime>
#include <cmath>
#include <vector>
// root headers
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TProfile.h"
// user headers
#include "progress.h"
#include "GammaJetAnalysis.h"

using namespace std;
using namespace Pythia8;


// global variables
const float pTmin  = 8.;         // trigger pT cut
const float hMax   = 1.1;        // trigger abs(eta) cut
const float rMax   = 0.15;       // size of isolation cone
const float dPhiTR = M_PI/4.;    // "trigger" region endpoints
const float dPhiRE = M_PI/4.;    // "recoil" region endpoints
const float dPhiUE = M_PI/6.;    // "underlying evt." region endpoints
const float dPhiA1 = M_PI/4.;    // "away" region start
const float dPhiA2 = 3*M_PI/4.;  // "away" region end

TFile *Tfile = new TFile("pp200py8pt5gam.run0.root", "recreate");



int main(int argc, char* argv[]) {

  // check for and confirm correct input file
  if (argc != 2) {
   cerr << "PANIC: Unexpected number of command-line arguments.\n"
        << "       Provide a file name and nothing else!"
        << endl;
    return 1;
  }

  ifstream is(argv[1]);
  if (!is) {
    cerr << "PANIC: Command-line file " << argv[1] << " was not found."
         << endl;
    return 1;
  }

  cout << "Settings will be read from file " << argv[1] << endl;


  // initialize
  Pythia pythia;
  pythia.readFile(argv[1]);
  pythia.init();

  Analysis analysis;
  analysis.init(pTmin, hMax, rMax, dPhiTR, dPhiRE, dPhiUE, dPhiA1, dPhiA2, Tfile);
  timeBegin();

  // number of events, max number of aborts, etc.
  const int nEvent = pythia.mode("Main:numberOfEvents");
  const int nAbort = pythia.mode("Main:timesAllowErrors");
  bool      hasPL  = pythia.flag("PartonLevel:all");


  // event loop
  int code;
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cerr << "PANIC: Event generation aborted prematurely.\n"
           << "       Final event was " << iEvent
           << endl;
      break;
    }

    code = pythia.info.code();
    analysis.analyze(code, (hasPL ? pythia.event : pythia.process));
    //progress(iEvent, nEvent);

  }  // end event loop


  timeEnd();
  pythia.stat();
  analysis.finish();
  return 0;

}


// End ------------------------------------------------------------------------
