// 'MakeEtFitPlot.C'
// Derek Anderson
// 11.01.2021
//
// Makes a plot showing various
// eTtrg distributions fit with
// power laws.

#include <iostream>
#include <cstring>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NTrg(2);
static const UInt_t NPar(2);
static const UInt_t NVtx(4);
static const UInt_t NPlot(2);
static const UInt_t NDecMax(5);



void MakeEtFitPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning eTtrg fit plot script..." << endl;

  // io parameters
  const TString sOut("weightFunctionsForParticleGun.et920pi0xGam.d1m11y2021.root");
  const TString sInFile[NTrg] = {"eTtrgFit.dataVsPythia8_et920pi0xDir.d23m10y2020.root", "eTtrgFit.dataVsPythia8_et920pi0xDir.d23m10y2020.root"};
  const TString sInHist[NTrg] = {"hPy8pi0", "hPy8dir"};
  const TString sInName[NTrg] = {"hPythiaPi0", "hPythiaGam"};

  // text parameters
  const TString sSys("#bf{Pythia8, #sqrt{s} = 200 GeV}");
  const TString sFit("f[*](E_{T}^{trg}) = (a #upoint E_{T}^{trg})^{-b}");
  const TString sLabelTrg[NTrg] = {"#pi^{0} trig.", "#gamma trig."};
  const TString sLabelFit[NTrg] = {"f[#pi^{0}](E_{T}^{trg})", "f[#gamma](E_{T}^{trg})"};

  // fit parameters
  const UInt_t  fColFit[NTrg]      = {858, 898};
  const UInt_t  fLinFit[NTrg]      = {1, 1};
  const UInt_t  fWidFit[NTrg]      = {2, 2};
  const Float_t xFit[NPlot]        = {9., 20.};
  const Float_t pGuess[NTrg][NPar] = {{1., 7.}, {1., 5.}};
  const TString sFitName[NTrg]     = {"fFitPi0", "fFitGam"};
  const TString sFitFunc[NTrg]     = {"TMath::Power([0]*x, -1.*[1])", "TMath::Power([0]*x, -1.*[1])"};
  const TString sParName[NPar]     = {"a", "b"};

  // histogram parameters
  const TString sTitle("");
  const TString sTitleX("E_{T}^{trg} [GeV]");
  const TString sTitleY("a. u.");
  const Float_t xPlot[NPlot]  = {8., 21.};
  const Float_t fSizTrg[NTrg] = {1.25, 1.15};
  const UInt_t  fMarTrg[NTrg] = {20, 21};
  const UInt_t  fLinTrg[NTrg] = {1, 1};
  const UInt_t  fWidTrg[NTrg] = {1, 1};
  const UInt_t  fFilTrg[NTrg] = {0, 0};
  const UInt_t  fColTrg[NTrg] = {604, 636};

  // fit parameters

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fOut) {
    cerr << "PANIC: couldn't open output file!\n" << endl;
    return;
  }

  TFile *fIn[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fIn[iTrg] = new TFile(sInFile[iTrg].Data(), "read");
    if (!fIn[iTrg]) {
      cerr << "PANIC: couldn't open input file #" << iTrg << "!\n" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hTrg[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hTrg[iTrg] = (TH1D*) fIn[iTrg] -> Get(sInHist[iTrg].Data());
    if (!hTrg[iTrg]) {
      cerr << "PANIC: couldn't grab input histogram #" << iTrg << "!\n" << endl;
      return;
    }
  }
  cout << "    Grabbed histograms." << endl;

  // normalize histograms
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    const Double_t integral = hTrg[iTrg] -> Integral();
    const Double_t norm     = 1. / integral;
    if (integral > 0.) {
      hTrg[iTrg] -> Scale(norm);
    }
  }
  cout << "    Normalized histograms." << endl;

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hTrg[iTrg] -> SetMarkerColor(fColTrg[iTrg]);
    hTrg[iTrg] -> SetMarkerStyle(fMarTrg[iTrg]);
    hTrg[iTrg] -> SetMarkerSize(fSizTrg[iTrg]);
    hTrg[iTrg] -> SetFillColor(fColTrg[iTrg]);
    hTrg[iTrg] -> SetFillStyle(fFilTrg[iTrg]);
    hTrg[iTrg] -> SetLineColor(fColTrg[iTrg]);
    hTrg[iTrg] -> SetLineStyle(fLinTrg[iTrg]);
    hTrg[iTrg] -> SetLineWidth(fWidTrg[iTrg]);
    hTrg[iTrg] -> SetTitle(sTitle.Data());
    hTrg[iTrg] -> SetTitleFont(fTxt);
    hTrg[iTrg] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hTrg[iTrg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hTrg[iTrg] -> GetXaxis() -> SetTitleFont(fTxt);
    hTrg[iTrg] -> GetXaxis() -> SetTitleSize(fTit);
    hTrg[iTrg] -> GetXaxis() -> SetTitleOffset(fOffX);
    hTrg[iTrg] -> GetXaxis() -> SetLabelFont(fTxt);
    hTrg[iTrg] -> GetXaxis() -> SetLabelSize(fLab);
    hTrg[iTrg] -> GetXaxis() -> CenterTitle(fCnt);
    hTrg[iTrg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hTrg[iTrg] -> GetYaxis() -> SetTitleFont(fTxt);
    hTrg[iTrg] -> GetYaxis() -> SetTitleSize(fTit);
    hTrg[iTrg] -> GetYaxis() -> SetTitleOffset(fOffY);
    hTrg[iTrg] -> GetYaxis() -> SetLabelFont(fTxt);
    hTrg[iTrg] -> GetYaxis() -> SetLabelSize(fLab);
    hTrg[iTrg] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // fit histograms
  TF1     *fFit[NTrg];
  Double_t pFit[NTrg][NPar];
  Double_t pErr[NTrg][NPar];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {

    // initialize functions
    fFit[iTrg] = new TF1(sFitName[iTrg].Data(), sFitFunc[iTrg].Data(), xFit[0], xFit[1]);
    fFit[iTrg] -> SetLineColor(fColFit[iTrg]);
    fFit[iTrg] -> SetLineStyle(fLinFit[iTrg]);
    fFit[iTrg] -> SetLineWidth(fWidFit[iTrg]);
    for (UInt_t iPar = 0; iPar < NPar; iPar++) {
      fFit[iTrg] -> SetParameter(iPar, pGuess[iTrg][iPar]);
      fFit[iTrg] -> SetParName(iPar, sParName[iPar].Data());
    }

    // fit and get parameters
    hTrg[iTrg] -> Fit(fFit[iTrg], "R");
    for (UInt_t iPar = 0; iPar < NPar; iPar++) {
      pFit[iTrg][iPar] = fFit[iTrg] -> GetParameter(iPar);
      pErr[iTrg][iPar] = fFit[iTrg] -> GetParError(iPar);
    }
  }
  cout << "    Fit histograms." << endl;

  // add parameters to labels
  TString sFitWithParLabel[NTrg];
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {

    // prefix and suffix
    TString sPref("#bf{#color[");
    TString sSuff("]}}");
    sPref += fColFit[iTrg];
    sPref += "]{[";

    // add parameters
    TString sPars(sPref.Data());
    for (UInt_t iPar = 0; iPar < NPar; iPar++) {

      // convert floats to strings
      TString sParRaw("");
      TString sErrRaw("");
      sParRaw += pFit[iTrg][iPar];
      sErrRaw += pErr[iTrg][iPar];

      // locate decimals 
      const UInt_t nDecPar = sParRaw.First(".");
      const UInt_t nDecErr = sErrRaw.First(".");
      const UInt_t nTotErr = sErrRaw.Length();

      // determine number of digits
      UInt_t nZeroErr(0);
      for (UInt_t iErr = (nDecErr + 1); iErr < nTotErr; iErr++) {
        const Char_t  digit  = sErrRaw[iErr];
        const TString sDigit = digit;
        const Int_t   isZero = sDigit.CompareTo('0');
        if (isZero == 0) {
          ++nZeroErr;
        } else {
          break;
        }
      }
      const UInt_t nDigPar = nDecPar + (nZeroErr + 2);
      const UInt_t nDigErr = nDecErr + (nZeroErr + 2);

      // create truncated strings
      TString sParTxt("");
      TString sErrTxt("");
      sParTxt.Append(sParRaw, nDigPar);
      sErrTxt.Append(sErrRaw, nDigErr);

      // add to overall string
      sPars += sParName[iPar];
      sPars += " = ";
      sPars += sParTxt;
      sPars += " #pm ";
      sPars += sErrTxt;
      if (iPar != NPar) {
        sPars += ", ";
      }
    }  // end parameter loop
    sPars.Append(sSuff);

    // create label
    sFitWithParLabel[iTrg]  = sLabelFit[iTrg];
    sFitWithParLabel[iTrg] += " ";
    sFitWithParLabel[iTrg] += sPars;

  }  // trigger loop
  cout << "    Added parameters to fit labels." << endl;

  // create legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  nObjLe(2 * NTrg);
  const Float_t hObjLe(0.05);
  const Float_t hTotLe(nObjLe * hObjLe);
  const Float_t yTotLe(0.1 + hTotLe);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yTotLe};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    leg -> AddEntry(hTrg[iTrg], sLabelTrg[iTrg].Data(), "pf");
    leg -> AddEntry(fFit[iTrg], sFitWithParLabel[iTrg].Data(), "l");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NVtx]  = {0.3, 0.1, 0.5, 0.2};
  const Float_t fTxtLin[NVtx] = {0.0, 0.5, 1.0, 0.5};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddLine(fTxtLin[0], fTxtLin[1], fTxtLin[2], fTxtLin[3]);
  txt -> AddText(sFit.Data());
  cout << "    Made text." << endl;

  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
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
  const Float_t fMarginB(0.15);

  TCanvas *cPlot = new TCanvas("cPlot", "", width, height);
  cPlot   -> SetGrid(fGrid, fGrid);
  cPlot   -> SetTicks(fTick, fTick);
  cPlot   -> SetLogx(fLogX);
  cPlot   -> SetLogy(fLogY);
  cPlot   -> SetBorderMode(fMode);
  cPlot   -> SetBorderSize(fBord);
  cPlot   -> SetFrameBorderMode(fFrame);
  cPlot   -> SetLeftMargin(fMarginL);
  cPlot   -> SetRightMargin(fMarginR);
  cPlot   -> SetTopMargin(fMarginT);
  cPlot   -> SetBottomMargin(fMarginB);
  cPlot   -> cd();
  hTrg[0] -> Draw();
  for (UInt_t iTrg = 1; iTrg < NTrg; iTrg++) {
    hTrg[iTrg] -> Draw("same");
  }
  leg   -> Draw();
  txt   -> Draw();
  fOut  -> cd();
  cPlot -> Write();
  cPlot -> Close();
  cout << "    Made plot." << endl;

  // save histograms
  fOut -> cd();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    hTrg[iTrg] -> Write();
    fFit[iTrg] -> Write();
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOut -> cd();
  fOut -> Close();
  for (UInt_t iTrg = 0; iTrg < NTrg; iTrg++) {
    fIn[iTrg] -> cd();
    fIn[iTrg] -> Close();
  }
  cout << "  Finished eT fit plot script!\n" << endl;

}

// End ------------------------------------------------------------------------
