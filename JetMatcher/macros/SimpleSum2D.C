// 'SimpleSum2D.C'
// Derek Anderson
// 03.31.2020
//
// Use this to sum a set of histograms with
// a given set of weights from embedding.
//
// NOTE: Configuration = 0, RFF configuration
//       Configuration = 1, FF configuration
//       Configuration = 2, both (1st RFF, then FF)


#include <iostream>
#include "TH2.h"
#include "TFile.h"
#include "TString.h"

using namespace std;


// constants
static const UInt_t   NHistRFF(8);
static const UInt_t   NHistFF(7);
static const UInt_t   NTotal(10);
static const UInt_t   NHistAll(NHistRFF + NHistFF);
static const Double_t NTrgsFF[NHistFF]   = {2., 37., 163., 463., 2419., 8762., 4455.};
static const Double_t NTrgsRFF[NHistRFF] = {1., 5., 36., 204., 628., 2807., 7547., 3855.};
static const Double_t NTrgsAll[NHistAll] = {1., 5., 36., 204., 628., 2807., 7547., 3855., 2., 37., 163., 463., 2419., 8762., 4455.};
static const Double_t WeightsFF[NTotal]  = {1.0, 3.361596e-01, 1.401161e-01, 1.337302e-01, 2.895246e-02, 1.042577e-02, 8.294575e-03, 2.064352e-03, 8.088693e-05, 1.417116e-05};
static const Double_t WeightsRFF[NTotal] = {1.0, 3.501425e-01, 1.395103e-01, 1.326444e-01, 2.801546e-02, 1.031377e-02, 8.210314e-03, 1.985107e-03, 8.054588e-05, 1.449037e-05};

// options
static const UInt_t NHist(NHistAll);
static const UInt_t Configuration(2);
static const Bool_t DoBinNorm(true);
static const Bool_t DoIntNorm(true);
static const Bool_t UseRootNTrgs(false);



void SimpleSum2D() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString sOut("jetQtVsPtPar.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.d16m5y2021.root");
  const TString sIn[NHist]   = {"output/pp200r9pt4rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt5rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt7rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt9rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt11rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt15rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt25rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt35rff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt5ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt7ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt9ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt11ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt15ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt25ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root", "output/pp200r9pt35ff.forPtParVsQtCheck.et1520pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.d16m5y2021.root"};
  const TString sHist[NHist] = {"hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt", "hJetPtParVsQt"};
  const TString sNorm[NHist] = {"EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP", "EventInfo/hRefmultP"};


  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fIn[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] = new TFile(sIn[iHist].Data(), "read");
    if (!fIn[iHist]) {
      cerr << "PANIC: couldn't open input file no. " << iHist << "!" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // get histograms
  TH2D *hHist[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hHist[iHist] = (TH2D*) fIn[iHist] -> Get(sHist[iHist].Data());
    if (!hHist[iHist]) {
      cerr << "PANIC: couldn't grab histogram no. " << iHist << "!" << endl;
      return;
    }
  }
  fOut -> cd();
  cout << "    Grabbed histograms." << endl;


  Double_t NormsFF[NHistFF];
  Double_t NormsRFF[NHistRFF];
  Double_t NormsAll[NHistAll];
  TH2D     *hNormFF[NHistFF];
  TH2D     *hNormRFF[NHistRFF];
  TH2D     *hNormAll[NHistAll];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        if (UseRootNTrgs) {
          hNormRFF[iHist] = (TH2D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsRFF[iHist] = hNormRFF[iHist] -> GetEntries();
        }
        else
          NormsRFF[iHist] = NTrgsRFF[iHist];
        break;
      case 1:
        if (UseRootNTrgs) {
          hNormFF[iHist] = (TH2D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsFF[iHist] = hNormFF[iHist] -> GetEntries();
        }
        else
          NormsFF[iHist] = NTrgsFF[iHist];
        break;
      case 2:
        if (UseRootNTrgs) {
          hNormAll[iHist] = (TH2D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsAll[iHist] = hNormAll[iHist] -> GetEntries();
        }
        else
          NormsAll[iHist] = NTrgsAll[iHist];
        break;
    }  // end switch
  }  // end histogram loop
  cout << "    Extracted no. of triggers." << endl;


  // scale histograms
  Double_t normer(0.);
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        hHist[iHist] -> Scale(WeightsRFF[(NTotal - NHist) + iHist]);
        normer += NormsRFF[iHist] * WeightsRFF[(NTotal - NHist) + iHist];
        break;
      case 1:
        hHist[iHist] -> Scale(WeightsFF[(NTotal - NHist) + iHist]);
        normer += NormsFF[iHist] * WeightsFF[(NTotal - NHist) + iHist];
        break;
      case 2:
        const Bool_t isRFF = (iHist < NHistRFF);
        if (isRFF) {
          hHist[iHist] -> Scale(WeightsRFF[(NTotal - NHistRFF) + iHist]);
          normer += NormsAll[iHist] * WeightsRFF[(NTotal - NHistRFF) + iHist];
        }
        else {
          hHist[iHist] -> Scale(WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)]);
          normer += NormsAll[iHist] * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
        }
        break;
    }
  }
  cout << "    Scaled histograms." << endl;

  // sum histograms
  TH2D *hSum  = (TH2D*) hHist[0] -> Clone();
  TH2D *hNorm = (TH2D*) hHist[0] -> Clone();
  hSum  -> SetName("hSum");
  hNorm -> SetName("hNorm");
  hSum  -> Reset("ICE");
  hNorm -> Reset("ICE");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hSum  -> Add(hHist[iHist]);
    hNorm -> Add(hHist[iHist]);
  }

  // normalize by bin width
  if (DoBinNorm) {
    const UInt_t nBinX = hSum -> GetNbinsX();
    const UInt_t nBinY = hSum -> GetNbinsY();
    for (UInt_t iBinX = 1; iBinX < (nBinX + 1); iBinX++) {
      for (UInt_t iBinY = 1; iBinY < (nBinY + 1); iBinY++) {
        const Double_t sumVal  = hSum  -> GetBinContent(iBinX, iBinY);
        const Double_t normVal = hNorm -> GetBinContent(iBinX, iBinY);
        const Double_t sumErr  = hSum  -> GetBinError(iBinX, iBinY);
        const Double_t normErr = hNorm -> GetBinError(iBinX, iBinY);
        const Double_t xWidth  = hSum  -> GetXaxis() -> GetBinWidth(iBinX);
        const Double_t yWidth  = hSum  -> GetYaxis() -> GetBinWidth(iBinY);
        const Double_t dArea   = xWidth * yWidth;
        hSum  -> SetBinContent(iBinX, iBinY, sumVal / dArea);
        hNorm -> SetBinContent(iBinX, iBinY, normVal / dArea);
        hSum  -> SetBinError(iBinX, iBinY, sumErr / dArea);
        hNorm -> SetBinError(iBinX, iBinY, normErr / dArea);
      }
    }
    cout << "    Normalized sum by bin width." << endl;
  }

  // normalize sum
  if (DoIntNorm) {
    const Double_t integral = hNorm -> Integral();
    hNorm -> Scale(1. / integral);
    cout << "    Normalized sum.\n"
         << "      normalization = " << integral
         << endl;
  }
  else {
    hNorm -> Scale(1. / normer);
    cout << "    Normalized sum.\n"
         << "      normalization = " << normer
         << endl;
  }


  // save and close
  fOut  -> cd();
  hSum  -> Write();
  hNorm -> Write();
  fOut  -> Close();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fIn[iHist] -> cd();
    fIn[iHist] -> Close();
  }
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
