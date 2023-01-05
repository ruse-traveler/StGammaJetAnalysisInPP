// 'SumAndStreamErrors.C'
// Derek Anderson
// 04.28.2022
//
// Sum up some histograms and their
// errors why don'cha?  There should
// be NSys + 2 histograms: NSys
// histograms with systematic errors,
// 1 histogram with statistical errors,
// and the last one will be the sum
// of the two.

#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

using namespace std;

// global constants
static const UInt_t  NSys(2);
static const UInt_t  NTot(NSys + 1);
static const UInt_t  NMaxBin(100);
static const Float_t XStart(0.);
static const Float_t XStop(30.);



void SumAndStreamErrors() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Summing errors..." << endl;

  // i/o files
  const TString sOut("summedErrors.forBUR_withPerStat.et911r02pi0.d2m6y2022.root");
  const TString sStat("unfolding/et911r02pi0_rebinUnfold/summedErrors.forLongPaper.et911r02pi0.d9m1y2022.root");
  const TString sStatTxt("summedErrors.forBUR_stat.et911r02pi0.d2m6y2022.txt");
  const TString sSys[NSys] = {"unfolding/et911r02pi0_rebinUnfold/summedErrors.forLongPaper.et911r02pi0.d9m1y2022.root",
                              "unfolding/et911r02pi0_rebinUnfold/summedErrors.forLongPaper.et911r02pi0.d9m1y2022.root"};
  const TString sTxt[NTot] = {"summedErrors.forBUR_unf.et911r02pi0.d2m6y2022.txt",
                              "summedErrors.forBUR_det.et911r02pi0.d2m6y2022.txt",
                              "summedErrors.forBUR_tot.et911r02pi0.d2m6y2022.txt"};

  // histogram names
  const TString sStatIn("hStatistics");
  const TString sStatOut("hStatistics");
  const TString sSysTot("hSystematics");
  const TString sTotOut("hStatPlusSys");
  const TString sPerStat("hPerStat");
  const TString sSysIn[NSys]    = {"hSysUnf", "hSysDet"};
  const TString sSysOut[NSys]   = {"hSysUnf", "hSysDet"};
  const TString sPerOut[NTot]   = {"hPerSysUnf", "hPerSysDet", "hPerSysTot"};
  const TString sStackOut[NSys] = {"hStackUnf", "hStackUnfPlusDet"};

  // open files
  TFile *fOut  = new TFile(sOut.Data(), "recreate");
  TFile *fStat = new TFile(sStat.Data(), "read");
  if (!fStat) {
    cerr << "PANIC: couldn't open statistics file!" << endl;
    return;
  }

  TFile *fSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fSys[iSys] = new TFile(sSys[iSys].Data(), "read");
    if (!fSys[iSys]) {
      cerr << "PANIC: couldn't systematic file " << iSys << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hStat = (TH1D*) fStat -> Get(sStatIn.Data());
  if (!hStat) {
    cerr << "PANIC: couldn't grab statistics histogram!" << endl;
    return;
  }

  TH1D *hSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys] = (TH1D*) fSys[iSys] -> Get(sSysIn[iSys].Data());
    if (!hSys[iSys]) {
      cerr << "PANIC: couldn't grab systematic histogram " << iSys << endl;
      return;
    }
    hSys[iSys] -> SetName(sSysOut[iSys].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // initialize output histograms
  TH1D *hStatOut = (TH1D*) hStat -> Clone();
  TH1D *hPerStat = (TH1D*) hStat -> Clone();
  TH1D *hSysOut  = (TH1D*) hStat -> Clone();
  TH1D *hTotOut  = (TH1D*) hStat -> Clone();
  hStatOut -> SetName(sStatOut.Data());
  hPerStat -> SetName(sPerStat.Data());
  hSysOut  -> SetName(sSysTot.Data());
  hTotOut  -> SetName(sTotOut.Data());
  hStatOut -> Reset("ICES");
  hPerStat -> Reset("ICES");
  hSysOut  -> Reset("ICES");
  hTotOut  -> Reset("ICES");

  TH1D *hPerSys[NTot];
  TH1D *hStackSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hPerSys[iSys]   = (TH1D*) hStat -> Clone();
    hStackSys[iSys] = (TH1D*) hStat -> Clone();
    hPerSys[iSys]   -> SetName(sPerOut[iSys].Data());
    hStackSys[iSys] -> SetName(sStackOut[iSys].Data());
    hPerSys[iSys]   -> Reset("ICES");
    hStackSys[iSys] -> Reset("ICES"); 
  }
  hPerSys[NSys] = (TH1D*) hStat -> Clone();
  hPerSys[NSys] -> SetName(sPerOut[NSys].Data());
  hPerSys[NSys] -> Reset("ICES");

  // for summing errors
  Double_t valStat[NMaxBin];
  Double_t errStat[NMaxBin];
  Double_t perStat[NMaxBin];
  Double_t perSys2[NMaxBin];
  Double_t perTot2[NMaxBin];

  // determine range
  const Int_t  iXstart = hStat -> FindBin(XStart);
  const UInt_t iFirst  = hStat -> FindFirstBinAbove(0.);
  const UInt_t iLast   = hStat -> FindLastBinAbove(0.);
  const UInt_t iStop   = hStat -> GetNbinsX();

  // set starting bin
  UInt_t iStart(0);
  if (iFirst < iXstart) {
    iStart = iXstart;
  } else {
    iStart = iFirst;
  }

  // count up bins
  UInt_t  nBins(0);
  UInt_t  iCount(iStart);
  Bool_t  isInRange(true);
  Float_t binCenter(0.);
  do {
    binCenter = hStat -> GetBinCenter(iCount);
    isInRange = ((binCenter >= XStart) && (binCenter < XStop));
    ++iCount;
    if (isInRange) ++nBins;
  } while (isInRange);
  cout << "    Range and number of bins determined:\n"
       << "      range = (" << XStart << ", " << XStop << "), nBins = " << nBins
       << endl;

  // initialize errors
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    const UInt_t iVal   = iBin + iStart;
    const Bool_t isLast = (iVal > iLast);
    const Bool_t isDone = (iVal > iStop);
    if (isLast || isDone) break;

    const Double_t val  = hStat -> GetBinContent(iVal);
    const Double_t stat = hStat -> GetBinError(iVal);
    const Double_t per  = stat / val;
    valStat[iBin] = val;
    errStat[iBin] = stat;
    perStat[iBin] = per;
    perSys2[iBin] = 0.;
    perTot2[iBin] = per * per;
  }
  cout << "    Initialized errors." << endl;

  // loop over systematics
  ofstream oTot(sTxt[NSys].Data());
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {

    // set systematic start
    UInt_t iFirstSys = hSys[iSys] -> FindFirstBinAbove(0.);
    UInt_t iCheck    = iFirstSys;
    if (iCheck < iXstart) {
      iCheck = iXstart;
    }
    if (iCheck < iStart) {
      iCheck = iStart;
    }

    // check start point
    if (iStart != iCheck) {
      cerr << "WARNING: looks like bin numbering is off in systematic no. " << iSys << "!\n"
           << "         iStart = " << iStart << ", iCheck = " << iCheck
           << endl;
    }

    // loop over bins
    ofstream oTxt(sTxt[iSys].Data());
    for (UInt_t iBin = 0; iBin < nBins; iBin++) {
      const UInt_t iVal   = iBin + iStart;
      const Bool_t isLast = (iVal > iLast);
      const Bool_t isDone = (iVal > iStop);
      if (isLast || isDone) break;

      const Double_t xBin   = hSys[iSys] -> GetBinCenter(iVal);
      const Double_t valAdd = hSys[iSys] -> GetBinContent(iVal);
      const Double_t sysAdd = hSys[iSys] -> GetBinError(iVal);
      const Double_t perAdd = sysAdd / valAdd;
      if (valAdd > 0.) {
        perSys2[iBin] += perAdd * perAdd;
        perTot2[iBin] += perAdd * perAdd;
      }

      // set first line of output stream
      if (iBin == 0) {
        oTxt << "binNumber binCenter percentErr" << endl;
        if ((iSys + 1) == NSys) oTot << "binNumber binCenter percentErr" << endl;
      }

      // stream out values
      oTxt << iBin;
      oTxt << " ";
      oTxt << xBin;
      oTxt << " ";
      oTxt << perAdd;
      oTxt << endl;
      if ((iSys + 1) == NSys) {
        oTot << iBin;
        oTot << " ";
        oTot << xBin;
        oTot << " ";
        oTot << TMath::Sqrt(perSys2[iBin]);
        oTot << endl;
      }

      // set percent histogram values
      hPerSys[iSys]   -> SetBinContent(iVal, 1.);
      hStackSys[iSys] -> SetBinContent(iVal, 1.);
      hPerSys[iSys]   -> SetBinError(iVal, perAdd);
      hStackSys[iSys] -> SetBinError(iVal, TMath::Sqrt(perSys2[iBin]));
      if ((iSys + 1) == NSys) {
        hPerSys[NSys] -> SetBinContent(iVal, 1.);
        hPerSys[NSys] -> SetBinError(iVal, TMath::Sqrt(perSys2[iBin]));
      }
    }  // end bin loop
    oTxt.close();
  }  // end systematic loop
  oTot.close();

  // stream out percent stat errors
  ofstream oStat(sStatTxt.Data());
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    const UInt_t iVal   = iBin + iStart;
    const Bool_t isLast = (iVal > iLast);
    const Bool_t isDone = (iVal > iStop);
    if (isLast || isDone) break;

    const Double_t xBin = hStatOut -> GetBinCenter(iVal);
    if (iBin == 0) {
      oStat << "binNumber binCenter percentErr" << endl;
    }
    oStat << iBin;
    oStat << " ";
    oStat << xBin;
    oStat << " ";
    oStat << perStat[iBin];
    oStat << endl;

    // set percent stat histogram values
    hPerStat -> SetBinContent(iVal, 1.);
    hPerStat -> SetBinError(iVal, perStat[iBin]);
  }  // end bin loop
  oStat.close();
  cout << "    Summed errors." << endl;

  Double_t perSys[NMaxBin];
  Double_t errSys[NMaxBin];
  Double_t perTot[NMaxBin];
  Double_t errTot[NMaxBin];
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    perSys[iBin] = TMath::Sqrt(perSys2[iBin]);
    perTot[iBin] = TMath::Sqrt(perTot2[iBin]);
    errSys[iBin] = perSys[iBin] * valStat[iBin];
    errTot[iBin] = perTot[iBin] * valStat[iBin];
    hStatOut -> SetBinContent(iBin + iStart, valStat[iBin]);
    hStatOut -> SetBinError(iBin + iStart, errStat[iBin]);
    hSysOut  -> SetBinContent(iBin + iStart, valStat[iBin]);
    hSysOut  -> SetBinError(iBin + iStart, errSys[iBin]);
    hTotOut  -> SetBinContent(iBin + iStart, valStat[iBin]);
    hTotOut  -> SetBinError(iBin + iStart, errTot[iBin]);
  }
  cout << "    Calculated errors." << endl;

  // write histograms
  fOut     -> cd();
  hStatOut -> Write();
  hPerStat -> Write();
  hSysOut  -> Write();
  hTotOut  -> Write();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys]      -> Write();
    hPerSys[iSys]   -> Write();
    hStackSys[iSys] -> Write();
  }
  hPerSys[NSys] -> Write();

  // close files
  fOut  -> Close();
  fStat -> cd();
  fStat -> Close();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fSys[iSys] -> cd();
    fSys[iSys] -> Close();
  }
  cout << "  DONE." << endl;

}

// End ------------------------------------------------------------------------
