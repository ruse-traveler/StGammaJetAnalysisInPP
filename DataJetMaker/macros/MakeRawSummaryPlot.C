// 'MakeRawSummaryPlot.C'
// Derek Anderson
// 09.23.2021
//
// Creates a plot summarizing
// the fully corrected recoil
// jet distributions.
//
// NOTE: use the 'JetRes' variable
//       to change how certain spectra
//       are handled.
//         JetRes = 0 -- R = 0.2
//         JetRes = 1 -- R = 0.5

#include <iostream>
#include "TH1.h"
#include "TLine.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

// global constants
static const UInt_t NVtx(4);
static const UInt_t NErr(4);
static const UInt_t NPlot(2);
static const UInt_t NJetRes(2);
static const UInt_t NTrgPi0(2);
static const UInt_t NTrgGam(3);
static const UInt_t NTrgTot(NTrgPi0 + NTrgGam);

// trigger parameters
static const UInt_t JetRes(1);
static const Bool_t IsGammaDir(false);



void MakeRawSummaryPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Plotting raw distributions..." << endl;

  // io parameters
  const TString sOut("gammaVsPionPythiaSummary_betterRange_forQM22.et920r02.d31m3y2022.root");
  const TString sHistPi0("Pi0/hJetPtCorrP");
  const TString sHistGam("Gam/hJetPtCorrG");
  const TString sInPi0[NTrgPi0] = {"../PythiaJetMaker/output/March2022/pp200py8par.forQM22_pTbinOne.et911pt0230r02pi0.d30m3y2022.root",
                                   "../PythiaJetMaker/output/March2022/pp200py8par.forQM22_pTbinOne.et1115pt0230r02pi0.d30m3y2022.root"};
  const TString sInGam[NTrgGam] = {"../PythiaJetMaker/output/March2022/pp200py8par.forQM22_pTbinOne.et911pt0230r02gam.d30m3y2022.root",
                                   "../PythiaJetMaker/output/March2022/pp200py8par.forQM22_pTbinOne.et1115pt0230r02gam.d30m3y2022.root",
                                   "../PythiaJetMaker/output/March2022/pp200py8par.forQM22_pTbinOne.et1520pt0230r02gam.d30m3y2022.root"};

  // text parameters
  const TString sSys("PYTHIA-8, #sqrt{s} = 200 GeV");
  const TString sTyp("#bf{charged jets}");

  // plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T,jet}^{reco,ch} [GeV/c]");
  const TString  sTitleY("(1/N_{trig}) d^{2}N_{jet}/(dp_{T,jet}^{reco,ch} d#eta_{jet}) [GeV/c]^{-1}");
  const Double_t xPlotRange[NPlot] = {-1., 25.};
  const Double_t yPlotRange[NPlot] = {0.000003, 777.};

  // histogram parameters
  const Float_t xMinJet(0.2);
  const Float_t fMarSize(1.5);
  const TString sName[NTrgTot]   = {"h911pi0", "h1115pi0", "h911dir", "h1115dir", "h1520dir"};
  const TString sLabels[NTrgTot] = {"9 - 11 GeV #pi^{0} [x1]", "11 - 15 GeV #pi^{0} [x10]",
                                    "9 - 11 GeV #gamma_{dir} [x1]", "11 - 15 GeV #gamma_{dir} [x10]", "15 - 20 GeV #gamma_{dir} [x100]"};

  // trigger and jet parameters
  const UInt_t  fCol[NJetRes][NTrgTot]      = {{1, 634, 1, 634, 602},
                                               {1, 634, 1, 634, 602}};
  const UInt_t  fMar[NJetRes][NTrgTot]      = {{24, 26, 20, 22, 23},
                                               {24, 26, 20, 22, 23}};
  const Float_t fTrgScale[NJetRes][NTrgTot] = {{1., 10., 1., 10., 100.},
                                               {1., 10., 1., 10., 100.}};
  const Float_t xMaxJet[NJetRes][NTrgTot]   = {{25., 25., 25., 25., 25.},
                                               {25., 25., 25., 25., 25.}};

  // style parameters
  const UInt_t  fLin(1);
  const UInt_t  fWid(1);
  const UInt_t  fFil(0);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.3);

  // select trigger and jets
  UInt_t  fColTrg[NTrgTot];
  UInt_t  fMarTrg[NTrgTot];
  Float_t fScalers[NTrgTot];
  Float_t xCutoff[NTrgTot];
  TString sTrg("");
  TString sJet("");
  switch (JetRes) {
    case 0:
      if (IsGammaDir) {
        sTrg = "#pi^{0}/#gamma_{dir}+jet, anti-k_{T}";
      } else {
        sTrg = "#pi^{0}/#gamma_{rich}+jet, anti-k_{T}";
      }
      sJet = "R = 0.2";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColTrg[iTrg]  = fCol[0][iTrg];
        fMarTrg[iTrg]  = fMar[0][iTrg];
        fScalers[iTrg] = fTrgScale[0][iTrg];
        xCutoff[iTrg]  = xMaxJet[0][iTrg];
      }
      break;
    case 1:
      if (IsGammaDir) {
        sTrg = "#pi^{0}/#gamma_{dir}+jet, anti-k_{T}";
      } else {
        sTrg = "#pi^{0}/#gamma_{rich}+jet, anti-k_{T}";
      }
      sJet = "R = 0.5";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColTrg[iTrg]  = fCol[1][iTrg];
        fMarTrg[iTrg]  = fMar[1][iTrg];
        fScalers[iTrg] = fTrgScale[1][iTrg];
        xCutoff[iTrg]  = xMaxJet[1][iTrg];
      }
      break;
    default:
      if (IsGammaDir) {
        sTrg = "#pi^{0}/#gamma_{dir}+jet, anti-k_{T}";
      } else {
        sTrg = "#pi^{0}/#gamma_{rich}+jet, anti-k_{T}";
      }
      sJet = "R = 0.2";
      for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
        fColTrg[iTrg]  = fCol[0][iTrg];
        fMarTrg[iTrg]  = fMar[0][iTrg];
        fScalers[iTrg] = fTrgScale[0][iTrg];
        xCutoff[iTrg]  = xMaxJet[0][iTrg];
      }
      break;
  }  // end switch-case

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!" << endl;
    return;
  }

  TFile *fInPi0[NTrgPi0];
  TFile *fInGam[NTrgGam];
  for (UInt_t iTrgPi0 = 0; iTrgPi0 < NTrgPi0; iTrgPi0++) {
    fInPi0[iTrgPi0] = new TFile(sInPi0[iTrgPi0].Data(), "read");
    if (!fInPi0[iTrgPi0]) {
      cerr << "PANIC: couldn't open a pi0 file!\n"
           << "       fInPi0[" << iTrgPi0 << "] = " << fInPi0[iTrgPi0]
           << endl;
      return;
    }
  }
  for (UInt_t iTrgGam = 0; iTrgGam < NTrgGam; iTrgGam++) {
    fInGam[iTrgGam] = new TFile(sInGam[iTrgGam].Data(), "read");
    if (!fInGam[iTrgGam]) {
      cerr << "PANIC: couldn't open a gamma file!\n"
           << "       fInGam[" << iTrgGam << "] = " << fInGam[iTrgGam]
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab input graphs
  TH1D *hData[NTrgTot];
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {

    // grab histogram
    if (iTrg < NTrgPi0) {
      hData[iTrg] = (TH1D*) fInPi0[iTrg] -> Get(sHistPi0.Data());
    } else {
      const UInt_t iGamma = iTrg - NTrgPi0;
      hData[iTrg] = (TH1D*) fInGam[iGamma] -> Get(sHistGam.Data());
    }

    // check if good
    if (!hData[iTrg]) {
      cerr << "PANIC: couldn't open an input histogram!\n"
           << "       hData[" << iTrg << "] = " << hData[iTrg]
           << endl;
      return;
    }
    hData[iTrg] -> SetName(sName[iTrg].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // suppress bins < xMinJet or > xMaxJet (for gamma-dir)
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    const UInt_t nBins = hData[iTrg] -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
      const Float_t binCnt     = hData[iTrg] -> GetBinCenter(iBin);
      const Bool_t  isAboveMin = (binCnt > xMinJet);
      const Bool_t  isBelowMax = (binCnt < xCutoff[iTrg]);
      if (!isAboveMin) {
        hData[iTrg] -> SetBinContent(iBin, 0.);
        hData[iTrg] -> SetBinError(iBin, 0.);
      }
      if (IsGammaDir && !isBelowMax) {
        hData[iTrg] -> SetBinContent(iBin, 0.);
        hData[iTrg] -> SetBinError(iBin, 0.);
      }
    }  // end bin loop
  }  // end trigger loop
  cout << "    Suppressed extremal bins." << endl;

  // set styles
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    hData[iTrg] -> SetLineColor(fColTrg[iTrg]);
    hData[iTrg] -> SetLineStyle(fLin);
    hData[iTrg] -> SetLineWidth(fWid);
    hData[iTrg] -> SetFillColor(fColTrg[iTrg]);
    hData[iTrg] -> SetFillStyle(fFil);
    hData[iTrg] -> SetMarkerColor(fColTrg[iTrg]);
    hData[iTrg] -> SetMarkerStyle(fMarTrg[iTrg]);
    hData[iTrg] -> SetMarkerSize(fMarSize);
    hData[iTrg] -> SetTitle(sTitle.Data());
    hData[iTrg] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hData[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hData[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hData[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hData[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hData[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hData[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hData[iTrg] -> GetYaxis() -> SetRangeUser(yPlotRange[0], yPlotRange[1]);
    hData[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hData[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hData[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hData[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hData[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hData[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  fAlnLe(12);
  const Float_t fLegXY[NVtx] = {0.7, 0.7, 0.9, 0.9};
  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAlnLe);
  leg -> AddEntry(hData[0], sLabels[0].Data(), "fp");
  leg -> AddEntry(hData[2], sLabels[3].Data(), "fp");
  leg -> AddEntry(hData[1], sLabels[1].Data(), "fp");
  leg -> AddEntry(hData[3], sLabels[4].Data(), "fp");
  leg -> AddEntry(hData[4], sLabels[2].Data(), "fp");
  cout << "    Made legend." << endl;

  // for QM 2022
  TH1D *hLegPi0 = (TH1D*) hData[0] -> Clone();
  TH1D *hLegGam = (TH1D*) hData[2] -> Clone();
  hLegPi0 -> SetName("hLegPi0");
  hLegGam -> SetName("hLegGam");
  hLegPi0 -> SetFillColor(921);
  hLegGam -> SetFillColor(921);
  hLegPi0 -> SetFillStyle(0);
  hLegGam -> SetFillStyle(0);
  hLegPi0 -> SetLineColor(921);
  hLegGam -> SetLineColor(921);
  hLegPi0 -> SetLineWidth(2);
  hLegGam -> SetLineWidth(2);
  hLegPi0 -> SetMarkerColor(921);
  hLegGam -> SetMarkerColor(921);

  TH1D *hLeg911  = (TH1D*) hData[2] -> Clone();
  TH1D *hLeg1115 = (TH1D*) hData[3] -> Clone();
  TH1D *hLeg1520 = (TH1D*) hData[4] -> Clone();
  hLeg911  -> SetName("hLeg911");
  hLeg1115 -> SetName("hLeg1115");
  hLeg1520 -> SetName("hLeg1520");
  hLeg911  -> SetLineWidth(2);
  hLeg1115 -> SetLineWidth(2);
  hLeg1520 -> SetLineWidth(2);

  TLegend *legQM = new TLegend(0.1, 0.1, 0.5, 0.25);
  legQM -> SetFillColor(0);
  legQM -> SetFillStyle(0);
  legQM -> SetLineColor(0);
  legQM -> SetLineStyle(0);
  legQM -> SetTextFont(42);
  legQM -> SetTextAlign(12);
  legQM -> SetNColumns(2);
  legQM -> AddEntry(hLeg911, "9 < E_{T}^{trig} < 11 GeV [x1]", "f");
  legQM -> AddEntry(hLegPi0, "#pi^{0}+jet", "pf");
  legQM -> AddEntry(hLeg1115, "11 < E_{T}^{trig} < 15 GeV [x10]", "f");
  legQM -> AddEntry(hLegGam, "#gamma_{dir}+jet", "pf");
  legQM -> AddEntry(hLeg1520, "15 < E_{T}^{trig} < 20 GeV [x100]", "f");

  TPaveText *ptQM = new TPaveText(0.5, 0.1, 0.7, 0.2, "NDC NB");
  ptQM -> SetFillColor(0);
  ptQM -> SetFillStyle(0);
  ptQM -> SetLineColor(0);
  ptQM -> SetLineStyle(0);
  ptQM -> SetTextFont(42);
  ptQM -> SetTextAlign(12);
  ptQM -> AddText("PYTHIA-8, p+p 200 GeV");
  ptQM -> AddText("anti-k_{T}, R = 0.2");

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const UInt_t  fAlnTx(32);
  const UInt_t  fAlnLb(22);
  const UInt_t  fFonLb(62);
  const Float_t fTxtTx[NVtx] = {0.5, 0.5, 0.7, 0.7};
  TPaveText *txt = new TPaveText(fTxtTx[0], fTxtTx[1], fTxtTx[2], fTxtTx[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlnTx);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());
  cout << "    Made text box." << endl;

  // scale histograms
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    hData[iTrg] -> Scale(fScalers[iTrg]);
  }
  cout << "    Scaled histograms." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT(0.02);
  const Float_t fMarginB(0.1);
  TCanvas *cSummary = new TCanvas("cRawData", "", width, height);
  cSummary -> SetGrid(fGrid, fGrid);
  cSummary -> SetTicks(fTick, fTick);
  cSummary -> SetLogx(fLogX);
  cSummary -> SetLogy(fLogY);
  cSummary -> SetBorderMode(fMode);
  cSummary -> SetBorderSize(fBord);
  cSummary -> SetFrameBorderMode(fFrame);
  cSummary -> SetLeftMargin(fMarginL);
  cSummary -> SetRightMargin(fMarginR);
  cSummary -> SetTopMargin(fMarginT);
  cSummary -> SetBottomMargin(fMarginB);
  cSummary -> cd();
  for (Int_t iTrgPlot = (NTrgTot - 1); iTrgPlot > -1; iTrgPlot--) {
    if (iTrgPlot == (NTrgTot - 1)) {
      hData[iTrgPlot] -> Draw();
    } else {
      hData[iTrgPlot] -> Draw("same");
    }
  }
  //txt      -> Draw();
  //leg      -> Draw();
  legQM    -> Draw();
  ptQM     -> Draw();
  fOut     -> cd();
  cSummary -> Write();
  cSummary -> Close();
  cout << "    Made summary plot." << endl;

  // save graphs
  fOut -> cd();
  for (UInt_t iTrg = 0; iTrg < NTrgTot; iTrg++) {
    hData[iTrg] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iTrgPi0 = 0; iTrgPi0 < NTrgPi0; iTrgPi0++) {
    fInPi0[iTrgPi0] -> cd();
    fInPi0[iTrgPi0] -> Close();
  }
  for (UInt_t iTrgGam = 0; iTrgGam < NTrgGam; iTrgGam++) {
    fInGam[iTrgGam] -> cd();
    fInGam[iTrgGam] -> Close();
  }
  cout << "  Finished making summary plot!" << endl;

}

// End ------------------------------------------------------------------------
