// 'MakeSmoothPriorPlot.C'
// Derek Anderson
// 10.27.2021
//
// Plots priors against the
// functions used to smooth
// them.

#include <iostream>
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t NParMax(13);
static const UInt_t NExpMax(4);
static const UInt_t NTanMax(2);
static const UInt_t NExpPar(2);
static const UInt_t NTanPar(2);
static const UInt_t NPrior(3);
static const UInt_t NPlot(2);
static const UInt_t NVtx(4);



void MakeSmoothPriorPlot() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning smooth prior plot macro..." << endl;

  // io parameters
  const TString sOutput("priorSmoothingFunction.et920r05pi0.d27m10y2021.root");
  const TString sHistPrior("hPrior");
  const TString sFitPrior("fTotalFit");
  const TString sInPrior[NPrior]    = {"pp200r9pi0.default_pTbinHuge.et911r05qt05130.p0m1k3n38t5.root", "pp200r9pi0.default_pTbinHuge.et1115r05qt05130.p0m1k3n38t5.root", "pp200r9pi0.default_pTbinHuge.et1520r05qt05130.p0m1k3n38t5.root"};
  const TString sNamePrior[NPrior]  = {"hPrior911", "hPrior1115", "hPrior1520"};
  const TString sLabelPrior[NPrior] = {"9 - 11 GeV prior", "11 - 15 GeV prior", "15 - 20 GeV prior"};

  // text parameters
  const TString sSys("Py6#oplusGeant, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} trig., E_{T}^{trg} #in (9, 20) GeV");
  const TString sJet("anti-k_{T} algo., R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // plot parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{par} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{3}N^{jet}/d(p_{T}^{par} #eta^{jet}) [GeV/c]^{-1}");
  const Float_t xPlotRange[NPlot] = {-1., 30.};
  const UInt_t  fColPri[NPrior]   = {893, 863, 813};
  const UInt_t  fMarPri[NPrior]   = {20, 21, 22};
  const UInt_t  fLinPri[NPrior]   = {1, 1, 1};
  const UInt_t  fWidPri[NPrior]   = {1, 1, 1};
  const UInt_t  fFilPri[NPrior]   = {0, 0, 0};

  // function parameters
  const Bool_t   useExistingFunctions(true);
  const TString  sNorm("norm");
  const TString  sExp[NExpPar]          = {"c", "b"};
  const TString  sTan[NTanPar]          = {"p_{T}", "a"};
  const TString  sNameFunc[NPrior]      = {"fSmooth911", "fSmooth1115", "fSmooth1520"};
  const TString  sFunction[NPrior]      = {"(expo(0) + expo(2)) * TMath::TanH((x-[4]) / [5]) * [6]",
                                           "(expo(0) + expo(2) + expo(4)) * (TMath::TanH((x-[6]) / [7]) + TMath::TanH((x-[8]) / [9])) * [10]",
                                           "(expo(0) + expo(2) + expo(4) + expo(6)) * (TMath::TanH((x-[8]) / [9]) + TMath::TanH((x-[10]) / [11])) * [12]"};
  const UInt_t   nExp[NPrior]           = {2, 3, 4};
  const UInt_t   nTan[NPrior]           = {1, 2, 2};
  const UInt_t   fColFunc[NPrior]       = {899, 859, 819};
  const UInt_t   fLinFunc[NPrior]       = {1, 1, 1};
  const UInt_t   fWidFunc[NPrior]       = {2, 2, 2};
  const Float_t  xFuncDef[NPlot]        = {0.2, 27.};
  const Double_t fPars[NPrior][NParMax] = {{0.390, -1.17, -1.78, -0.260, -8.12, 10.70, 1., 1., 1., 1., 1., 1., 1.},
                                           {0.390, -1.17, -2.44, -0.140, -1.78, -0.260, 7.14, 13.29, 7.33, -13.34, 1., 1., 1.},
                                           {0.310, -1.18, -1.36, -0.340, -3.81, 0.010, -0.490, -0.220, 9.07, 8.15, 9.28, -8.15, 1.}};

  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  if (!fOutput) {
    cerr << "PANIC: couldn't open file!\n" << endl;
    return;
  }

  TFile *fInPrior[NPrior];
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {
    fInPrior[iPrior] = new TFile(sInPrior[iPrior].Data(), "read");
    if (!fInPrior[iPrior]) {
      cerr << "PANIC: couldn't open prior file #" << iPrior << "!\n" << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab unfolding histograms
  TH1D *hPrior[NPrior];
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {
    hPrior[iPrior] = (TH1D*) fInPrior[iPrior] -> Get(sHistPrior.Data());
    if (!hPrior[iPrior]) {
      cerr << "PANIC: couldn't grab prior histogram #" << iPrior << "!\n" << endl;
      return;
    }
    hPrior[iPrior] -> SetName(sNamePrior[iPrior].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // create functions
  TF1 *fSmooth[NPrior];
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {

    // initialize function
    fSmooth[iPrior] = new TF1(sNameFunc[iPrior].Data(), sFunction[iPrior].Data(), xFuncDef[0], xFuncDef[1]);
    fSmooth[iPrior] -> SetLineColor(fColFunc[iPrior]);
    fSmooth[iPrior] -> SetLineStyle(fLinFunc[iPrior]);
    fSmooth[iPrior] -> SetLineWidth(fWidFunc[iPrior]);

    // set exponential parameters
    const UInt_t nExpTot = NExpPar * nExp[iPrior];
    for (UInt_t iExp = 0; iExp < nExp[iPrior]; iExp++) {

      // get array indices
      const UInt_t iExpConst = NExpPar * iExp;
      const UInt_t iExpSlope = iExpConst + 1;

      // create names
      TString sExpConst(sExp[0].Data());
      TString sExpSlope(sExp[1].Data());
      sExpConst += "_{";
      sExpSlope += "_{";
      sExpConst += iExp;
      sExpSlope += iExp;
      sExpConst += "}";
      sExpSlope += "}";

      // set parameters
      fSmooth[iPrior] -> FixParameter(iExpConst, fPars[iPrior][iExpConst]);
      fSmooth[iPrior] -> FixParameter(iExpSlope, fPars[iPrior][iExpSlope]);
      fSmooth[iPrior] -> SetParName(iExpConst, sExpConst.Data());
      fSmooth[iPrior] -> SetParName(iExpSlope, sExpSlope.Data());
    }

    // set hyper tangent parameters
    const UInt_t nTanTot = NTanPar * nTan[iPrior];
    for (UInt_t iTan = 0; iTan < nTan[iPrior]; iTan++) {

      // get array indices
      const UInt_t iTanZero = (NTanPar * iTan) + nExpTot;
      const UInt_t iTanCurv = iTanZero + 1;

      // create names
      TString sTanZero(sTan[0].Data());
      TString sTanCurv(sTan[1].Data());
      sTanZero += "^{";
      sTanCurv += "_{";
      sTanZero += iTan;
      sTanCurv += iTan;
      sTanZero += "}";
      sTanCurv += "}";

      // set parameters
      fSmooth[iPrior] -> FixParameter(iTanZero, fPars[iPrior][iTanZero]);
      fSmooth[iPrior] -> FixParameter(iTanCurv, fPars[iPrior][iTanCurv]);
      fSmooth[iPrior] -> SetParName(iTanZero, sTanZero.Data());
      fSmooth[iPrior] -> SetParName(iTanCurv, sTanCurv.Data());
    }

    // set normalization
    const UInt_t iNorm = nExpTot + nTanTot;
    fSmooth[iPrior] -> SetParameter(iNorm, fPars[iPrior][iNorm]);
    fSmooth[iPrior] -> SetParName(iNorm, sNorm.Data());
  }
  cout << "    Created functions.\n"
       << "    Normalizing functions:"
       << endl;

  // determine normalization
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {

    // calculate normalization
    const UInt_t   iStart   = hPrior[iPrior]  -> FindBin(xFuncDef[0]);
    const UInt_t   iStop    = hPrior[iPrior]  -> FindBin(xFuncDef[1]);
    const Double_t intHist  = hPrior[iPrior]  -> Integral(iStart, iStop);
    const Double_t intFunc  = fSmooth[iPrior] -> Integral(xFuncDef[0], xFuncDef[1]);
    const Double_t normFunc = intHist / intFunc;

    // set normalization
    const UInt_t iNorm = (NExpPar * nExp[iPrior]) + (NTanPar * nTan[iPrior]);
    fSmooth[iPrior] -> SetParameter(iNorm, normFunc);
  }

  // set styles
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fSiz(1.15);
  const Float_t fLab(0.04);
  const Float_t fTit(0.04);
  const Float_t fOffX(1.1);
  const Float_t fOffY(1.5);
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {
    hPrior[iPrior] -> SetMarkerColor(fColPri[iPrior]);
    hPrior[iPrior] -> SetMarkerStyle(fMarPri[iPrior]);
    hPrior[iPrior] -> SetMarkerSize(fSiz);
    hPrior[iPrior] -> SetFillColor(fColPri[iPrior]);
    hPrior[iPrior] -> SetFillStyle(fFilPri[iPrior]);
    hPrior[iPrior] -> SetLineColor(fColPri[iPrior]);
    hPrior[iPrior] -> SetLineStyle(fLinPri[iPrior]);
    hPrior[iPrior] -> SetLineWidth(fWidPri[iPrior]);
    hPrior[iPrior] -> SetTitle(sTitle.Data());
    hPrior[iPrior] -> SetTitleFont(fTxt);
    hPrior[iPrior] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPrior[iPrior] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPrior[iPrior] -> GetXaxis() -> SetTitleFont(fTxt);
    hPrior[iPrior] -> GetXaxis() -> SetTitleSize(fTit);
    hPrior[iPrior] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPrior[iPrior] -> GetXaxis() -> SetLabelFont(fTxt);
    hPrior[iPrior] -> GetXaxis() -> SetLabelSize(fLab);
    hPrior[iPrior] -> GetXaxis() -> CenterTitle(fCnt);
    hPrior[iPrior] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPrior[iPrior] -> GetYaxis() -> SetTitleFont(fTxt);
    hPrior[iPrior] -> GetYaxis() -> SetTitleSize(fTit);
    hPrior[iPrior] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPrior[iPrior] -> GetYaxis() -> SetLabelFont(fTxt);
    hPrior[iPrior] -> GetYaxis() -> SetLabelSize(fLab);
    hPrior[iPrior] -> GetYaxis() -> CenterTitle(fCnt);
    if (useExistingFunctions) {
      hPrior[iPrior] -> GetFunction(sFitPrior.Data()) -> SetLineColor(fColFunc[iPrior]);
      hPrior[iPrior] -> GetFunction(sFitPrior.Data()) -> SetLineStyle(fLinFunc[iPrior]);
      hPrior[iPrior] -> GetFunction(sFitPrior.Data()) -> SetLineWidth(fWidFunc[iPrior]);
    } else {
      hPrior[iPrior] -> GetFunction(sFitPrior.Data()) -> Delete();
    }
  }
  cout << "    Set styles." << endl;

  // make legend
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const Float_t hObjLe(0.05);
  const Float_t hTotLe(NPrior * hObjLe);
  const Float_t yTopLe(0.1 + hTotLe);
  const Float_t fLegXY[NVtx] = {0.1, 0.1, 0.3, yTopLe};

  TLegend *leg = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3]);
  leg -> SetFillColor(fColLe);
  leg -> SetFillStyle(fFilLe);
  leg -> SetLineColor(fColLe);
  leg -> SetLineStyle(fLinLe);
  leg -> SetTextFont(fTxt);
  leg -> SetTextAlign(fAln);
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {
    leg -> AddEntry(hPrior[iPrior], sLabelPrior[iPrior], "pf");
  }
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NVtx] = {0.3, 0.1, 0.5, 0.3};

  TPaveText *txt = new TPaveText(fTxtXY[0], fTxtXY[1], fTxtXY[2], fTxtXY[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAln);
  txt -> AddText(sSys.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());
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

  TCanvas *cPlot = new TCanvas("cPriorVsSmooth", "", width, height);
  cPlot     -> SetGrid(fGrid, fGrid);
  cPlot     -> SetTicks(fTick, fTick);
  cPlot     -> SetLogx(fLogX);
  cPlot     -> SetLogy(fLogY);
  cPlot     -> SetBorderMode(fMode);
  cPlot     -> SetBorderSize(fBord);
  cPlot     -> SetFrameBorderMode(fFrame);
  cPlot     -> SetLeftMargin(fMarginL);
  cPlot     -> SetRightMargin(fMarginR);
  cPlot     -> SetTopMargin(fMarginT);
  cPlot     -> SetBottomMargin(fMarginB);
  cPlot     -> cd();
  hPrior[0] -> Draw();
  if (!useExistingFunctions) {
    fSmooth[0] -> Draw("same");
  }
  for (UInt_t iPrior = 1; iPrior < NPrior; iPrior++) {
    hPrior[iPrior]  -> Draw("same");
    if (!useExistingFunctions) {
      fSmooth[iPrior] -> Draw("same");
    }
  }
  leg     -> Draw();
  txt     -> Draw();
  fOutput -> cd();
  cPlot   -> Write();
  cPlot   -> Close();
  cout << "    Made plot." << endl;

  // save histograms and functions
  fOutput -> cd();
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {
    hPrior[iPrior]  -> Write();
    if (!useExistingFunctions) {
      fSmooth[iPrior] -> Write();
    }
  }
  cout << "    Saved histograms." << endl;

  // close files
  fOutput -> cd();
  fOutput -> Close();
  for (UInt_t iPrior = 0; iPrior < NPrior; iPrior++) {
    fInPrior[iPrior] -> cd();
    fInPrior[iPrior] -> Close();
  }
  cout << "  Finished smooth prior plot!\n" << endl;

}

// End ------------------------------------------------------------------------
