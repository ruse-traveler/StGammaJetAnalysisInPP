// 'SimpleSum.C'
// Derek Anderson
// 01.17.2018
//
// Use this to sum a set of histograms with
// a given set of weights from embedding.
//
// NOTE: Configuration = 0, RFF configuration
//       Configuration = 1, FF configuration
//       Configuration = 2, both (1st RFF, then FF)

#include <iostream>
#include "TH1.h"
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
static const UInt_t  NHist(NHistAll);
static const UInt_t  Configuration(2);
static const Bool_t  DoBinNorm(true);
static const Bool_t  DoIntNorm(false);
static const Bool_t  DoJetNorm(true);
static const Bool_t  UseRootNTrgs(true);
static const Float_t JetRes(0.2);



void SimpleSum() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning sum script..." << endl;

  // io parameters
  const TString sOut("fakeJetCheck.jetArea.et911r02qt05130.d6m13y2022.root");
  const TString sIn[NHist]   = {"output/pp200r9pt4rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt5rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt7rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt9rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt11rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt15rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt25rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt35rff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt5ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt7ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt9ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt11ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt15ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt25ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root",
                                "output/pp200r9pt35ff.forFakeJetCheck_pTbinOne.et911r02qt05130.d7m6y2022.root"};
  const TString sHist[NHist] = {"JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ",
                                "JunkJets/hJetAreaJ"};
  const TString sNorm[NHist] = {"EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP",
                                "EventInfo/hRefmultP"};

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
  TH1D *hHist[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hHist[iHist] = (TH1D*) fIn[iHist] -> Get(sHist[iHist].Data());
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
  TH1D     *hNormFF[NHistFF];
  TH1D     *hNormRFF[NHistRFF];
  TH1D     *hNormAll[NHistAll];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    switch (Configuration) {
      case 0:
        if (UseRootNTrgs) {
          hNormRFF[iHist] = (TH1D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsRFF[iHist] = hNormRFF[iHist] -> GetEntries();
        } else {
          NormsRFF[iHist] = NTrgsRFF[iHist];
        }
        break;
      case 1:
        if (UseRootNTrgs) {
          hNormFF[iHist] = (TH1D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsFF[iHist] = hNormFF[iHist] -> GetEntries();
        } else {
          NormsFF[iHist] = NTrgsFF[iHist];
        }
        break;
      case 2:
        if (UseRootNTrgs) {
          hNormAll[iHist] = (TH1D*) fIn[iHist] -> Get(sNorm[iHist].Data());
          NormsAll[iHist] = hNormAll[iHist] -> GetEntries();
        } else {
          NormsAll[iHist] = NTrgsAll[iHist];
        }
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
        } else {
          hHist[iHist] -> Scale(WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)]);
          normer += NormsAll[iHist] * WeightsFF[(NTotal - NHistFF) + (iHist - NHistRFF)];
        }
        break;
    }
  }
  cout << "    Scaled histograms." << endl;

  // sum histograms
  TH1D *hSum  = (TH1D*) hHist[0] -> Clone();
  TH1D *hNorm = (TH1D*) hHist[0] -> Clone();
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
    for (UInt_t iBinX = 1; iBinX < (nBinX + 1); iBinX++) {
      const Double_t sumVal  = hSum  -> GetBinContent(iBinX);
      const Double_t normVal = hNorm -> GetBinContent(iBinX);
      const Double_t sumErr  = hSum  -> GetBinError(iBinX);
      const Double_t normErr = hNorm -> GetBinError(iBinX);
      const Double_t xWidth  = hSum  -> GetXaxis() -> GetBinWidth(iBinX);
      hSum  -> SetBinContent(iBinX, sumVal / xWidth);
      hNorm -> SetBinContent(iBinX, normVal / xWidth);
      hSum  -> SetBinError(iBinX, sumErr / xWidth);
      hNorm -> SetBinError(iBinX, normErr / xWidth);
    }
    cout << "    Normalized sum by bin width." << endl;
  }

  // normalize sum
  if (DoIntNorm) {
    const Double_t integral = hNorm -> Integral();
    hNorm -> Scale(1. / integral);
    cout << "    Normalized histogram.\n"
         << "      normalization = " << integral
         << endl;
  } else {
    hNorm -> Scale(1. / normer);
    cout << "    Normalized histogram.\n"
         << "      normalization = " << normer
         << endl;
  }

  if (DoJetNorm) {
    const Double_t jetNorm = 2. * (1. - JetRes);
    hNorm -> Scale(1. / jetNorm);
    cout << "    Did jet normalization.\n"
         << "      2(1 - Rjet) = " << jetNorm
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
