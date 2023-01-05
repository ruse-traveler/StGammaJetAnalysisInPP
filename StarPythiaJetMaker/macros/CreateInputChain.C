// 'CreateInputChain.C'
// Derek Anderson
// 06.24.2021
//
// Takes a list of files and
// links them together into
// a TChain.

#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TList.h"
#include "TError.h"
#include "TChain.h"
#include "TString.h"

using namespace std;


void CreateInputChain(const TString sInputList, const TString sOutputFile, const TString sOutputChain, const TString sOutputBad, const TString sOutputGood) {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Creating input TChain..." << endl;

  // create output file and TChain
  TFile  *fOut = new TFile(sOutputFile.Data(), "recreate");
  TChain *cOut = new TChain(sOutputChain.Data(), "");
  if (!fOut) {
    cerr << "PANIC: couldn't create output file!" << endl;
    return;
  }
  cout << "    Created output file." << endl;

  // open input and output streams
  ifstream ifInput(sInputList.Data());
  ofstream ofBad(sOutputBad.Data());
  ofstream ofGood(sOutputGood.Data());
  if (!ifInput) {
    cerr << "PANIC: couldn't open input file list!" << endl;
    return;
  }
  cout << "    Opened input and output streams.\n"
       << "    Linking files..."
       << endl;

  // link files
  Int_t  iRecover(0);
  Int_t  iAdd(0);
  UInt_t iFile(0);
  Bool_t isBad(false);
  string sInput;
  TFile  *fInput;
  while (ifInput) {

    // read in file path
    ifInput >> sInput;
    TString sName(sInput);

    // check file
    fInput   = new TFile(sName.Data(), "read");
    iRecover = fInput -> Recover();
    isBad    = (iRecover == 0);

    // add file
    if (isBad) {
      cout  << "      File '" << sName.Data() << "' no good!" << endl;
      ofBad << sName.Data();
      ofBad << endl;
    } else {
      iAdd = cOut -> Add(sName.Data(), 0);
      if (iAdd == 0) {
        cout  << "      File '" << sName.Data() << "' no good!" << endl;
        ofBad << sName.Data();
        ofBad << endl;
      } else {
        cout << "      File '" << sName.Data() << "' merged." << endl;
        ofGood << sName.Data();
        ofGood << endl;
      }
    }
    iFile++;

  }  // end while loop
  cout << "    Files linked!\n" << endl;

  // save output and close
  fOut -> cd();
  cOut -> Write();
  fOut -> Close();
  cout << "  Created input TChain!\n" << endl;

}

// End ------------------------------------------------------------------------
