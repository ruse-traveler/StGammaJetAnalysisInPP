// 'CompareProjections.C'
// Derek Anderson
// 06.05.2021
//
// Creates a series of plots comparing
// the projections of one response
// matrix to another using the output
// PlotResponseMatrixProjections.C


#include <iostream>
#include "TH1.h"
#include "TPad.h"
#include "TFile.h"
#include "TLine.h"
#include "TStyle.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// global constants
//static const UInt_t  NHist(18);
static const UInt_t  NHist(12);
static const UInt_t  NVtx(4);
static const UInt_t  NPts(2);
static const Float_t PtStart(0.2);
//static const Float_t PtStop(34.);
static const Float_t PtStop(38.);



void CompareProjections() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning response matrix projections comparison..." << endl;

  // io parameters
  const TString sInA("pi0TrgResponseProjections_smearWithSmoothRes_pTbinBig.et1520r05pi0.d21m7y2021.root");
  const TString sInB("gamTrgResponseProjections_smearWithSmoothRes_pTbinBig.et1520r05gam.d21m7y2021.root");
  const TString sOut("pi0vsGamTrgResponseProjections_smearWithSmoothRes_pTbinBig.et1520r05.d21m7y2021.root");
  //const TString sHist[NHist]  = {"hProjPt0p2", "hProjPt0p6", "hProjPt1", "hProjPt1p5", "hProjPt2", "hProjPt3", "hProjPt4", "hProjPt5", "hProjPt6p5", "hProjPt8", "hProjPt10", "hProjPt12", "hProjPt14p5", "hProjPt17", "hProjPt20", "hProjPt23", "hProjPt26p5", "hProjPt30"};
  const TString sHist[NHist]  = {"hProjPt0p2", "hProjPt0p6", "hProjPt1", "hProjPt1p5", "hProjPt2", "hProjPt3", "hProjPt5", "hProjPt8", "hProjPt12", "hProjPt17", "hProjPt23", "hProjPt30"};
  //const TString sRatio[NHist] = {"hRatioPt0p2", "hRatioPt0p6", "hRatioPt1", "hRatioPt1p5", "hRatioPt2", "hRatioPt3", "hRatioPt4", "hRatioPt5", "hRatioPt6p5", "hRatioPt8", "hRatioPt10", "hRatioPt12", "hRatioPt14p5", "hRatioPt17", "hRatioPt20", "hRatioPt23", "hRatioPt26p5", "hRatioPt30"};
  const TString sRatio[NHist] = {"hRatioPt0p2", "hRatioPt0p6", "hRatioPt1", "hRatioPt1p5", "hRatioPt2", "hRatioPt3", "hRatioPt5", "hRatioPt8", "hRatioPt12", "hRatioPt17", "hRatioPt23", "hRatioPt30"};
  //const TString sCan[NHist]   = {"cProjPt0p2", "cProjPt0p6", "cProjPt1", "cProjPt1p5", "cProjPt2", "cProjPt3", "cProjPt4", "cProjPt5", "cProjPt6p5", "cProjPt8", "cProjPt10", "cProjPt12", "cProjPt14p5", "cProjPt17", "cProjPt20", "cProjPt23", "cProjPt26p5", "cProjPt30"};
  const TString sCan[NHist]   = {"cProjPt0p2", "cProjPt0p6", "cProjPt1", "cProjPt1p5", "cProjPt2", "cProjPt3", "cProjPt5", "cProjPt8", "cProjPt12", "cProjPt17", "cProjPt23", "cProjPt30"};

  // histogram parameters
  const TString sSuffixA("_pi0");
  const TString sSuffixB("_gam");
  const TString sTitle("");
  const TString sTitleX("p_{T}^{det} [GeV/c]");
  const TString sTitleY("a. u.");
  const TString sTitleR("ratio");
  const Float_t fOffX(1.);
  const Float_t fOffXR(1.1);
  const Float_t fOffY(1.15);
  const Float_t fOffYR(0.62);
  const Float_t fLbl(0.03);
  const Float_t fLblR(0.056);
  const Float_t fTtlR(0.074);
  const Float_t xPlot[NPts]  = {-1., 34.};
  //const UInt_t  fColA[NHist] = {633, 809, 799, 401, 829, 819, 417, 849, 839, 433, 869, 859, 601, 889, 879, 617, 909, 899};
  //const UInt_t  fColB[NHist] = {625, 806, 796, 393, 826, 816, 409, 846, 836, 425, 866, 856, 593, 886, 876, 609, 906, 896};
  const UInt_t  fColA[NHist] = {809, 799, 829, 819, 849, 839, 869, 859, 889, 879, 909, 899};
  const UInt_t  fColB[NHist] = {806, 796, 826, 816, 846, 836, 866, 856, 886, 876, 906, 896};
  //const UInt_t  fMarA[NHist] = {20, 22, 23, 33, 21, 23, 34, 47, 29, 20, 22, 23, 33, 21, 23, 34, 47, 29};
  //const UInt_t  fMarB[NHist] = {24, 26, 32, 27, 25, 42, 28, 46, 30, 24, 26, 32, 27, 25, 42, 28, 46, 30};
  const UInt_t  fMarA[NHist] = {20, 22, 23, 33, 21, 23, 34, 47, 29, 20, 22, 23,};
  const UInt_t  fMarB[NHist] = {24, 26, 32, 27, 25, 42, 28, 46, 30, 24, 26, 32,};
  const UInt_t  fFil(0);
  const UInt_t  fLin(1);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);

  // legend parameters
  const TString sLegA("#pi^{0} trig. R_{ij} (smoothed)");
  const TString sLegB("#gamma_{dir} trig. R_{ij}");
  const TString sLegPrefix("p_{T}^{par} #in ");
  const TString sLegSuffix(" GeV/c");
  //const TString sLegBin[NHist] = {"(0.2, 0.6)", "(0.6, 1)", "(1, 1.5)", "(1.5, 2)", "(2, 3)", "(3, 4)", "(4, 5)", "(5, 6.5)", "(6.5, 8)", "(8, 10)", "(10, 12)", "(12, 14.5)", "(14.5, 17)", "(17, 20)", "(20, 23)", "(23, 26.5)", "(26.5, 30)", "(30, 34)"};
  const TString sLegBin[NHist] = {"(0.2, 0.6)", "(0.6, 1)", "(1, 1.5)", "(1.5, 2)", "(2, 3)", "(3, 5)", "(5, 8)", "(8, 12)", "(12, 17)", "(17, 23)", "(23, 30)", "(30, 38)"};
  //const TString sPadBin[NHist] = {"Pt0p2", "Pt0p6", "Pt1", "Pt1p5", "Pt2", "Pt3", "Pt4", "Pt5", "Pt6p5", "Pt8", "Pt10", "Pt12", "Pt14p5", "Pt17", "Pt20", "Pt23", "Pt26p5", "Pt30"};
  const TString sPadBin[NHist] = {"Pt0p2", "Pt0p6", "Pt1", "Pt1p5", "Pt2", "Pt3", "Pt5", "Pt8", "Pt12", "Pt17", "Pt23", "Pt30"};

  // text parameters
  const TString sSys("Py6#oplusGeant vs. Py8#oplusR^{Py6}_{ij}, #sqrt{s} = 200 GeV");
  const TString sTrg("#pi^{0} vs. #gamma_{dir} trig., E_{T}^{trg} #in (15, 20) GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // open files
  TFile *fInA = new TFile(sInA.Data(), "read");
  TFile *fInB = new TFile(sInB.Data(), "read");
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  if (!fInA || !fInB || !fOut) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fInA = " << fInA << ", fInB = " << fInB << ", fOut = " << fOut << "\n"
         << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hProjA[NHist];
  TH1D *hProjB[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hProjA[iHist] = (TH1D*) fInA -> Get(sHist[iHist].Data());
    hProjB[iHist] = (TH1D*) fInB -> Get(sHist[iHist].Data());
    if (!hProjA[iHist] || !hProjB[iHist]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "    hProjA[" << iHist << "] = " << hProjA[iHist] << ", hProjB[" << iHist << "] = " << hProjB[iHist] << "\n"
           << endl;
      return;
    }

    // set name
    TString sNameA = hProjA[iHist] -> GetName();
    TString sNameB = hProjB[iHist] -> GetName();
    sNameA.Append(sSuffixA);
    sNameB.Append(sSuffixB);
    hProjA[iHist] -> SetName(sNameA.Data());
    hProjB[iHist] -> SetName(sNameB.Data());
  }
  cout << "    Grabbed histograms." << endl;

  // calculate ratios
  TH1D *hRatio[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hRatio[iHist] = (TH1D*) hProjA[iHist] -> Clone();
    hRatio[iHist] -> SetName(sRatio[iHist].Data());
    hRatio[iHist] -> Reset("ICES");
    hRatio[iHist] -> Divide(hProjA[iHist], hProjB[iHist], 1., 1.);
  }
  cout << "    Calculated ratios." << endl;

  // set styles
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hProjA[iHist] -> SetFillColor(fColA[iHist]);
    hProjA[iHist] -> SetFillStyle(fFil);
    hProjA[iHist] -> SetLineColor(fColA[iHist]);
    hProjA[iHist] -> SetLineStyle(fLin);
    hProjA[iHist] -> SetMarkerColor(fColA[iHist]);
    hProjA[iHist] -> SetMarkerStyle(fMarA[iHist]);
    hProjA[iHist] -> SetTitle(sTitle.Data());
    hProjA[iHist] -> SetTitleFont(fTxt);
    hProjA[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hProjA[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
    hProjA[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hProjA[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hProjA[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hProjA[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hProjA[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hProjA[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
    hProjA[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hProjA[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hProjA[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hProjA[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hProjA[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hProjB[iHist] -> SetFillColor(fColB[iHist]);
    hProjB[iHist] -> SetFillStyle(fFil);
    hProjB[iHist] -> SetLineColor(fColB[iHist]);
    hProjB[iHist] -> SetLineStyle(fLin);
    hProjB[iHist] -> SetMarkerColor(fColB[iHist]);
    hProjB[iHist] -> SetMarkerStyle(fMarB[iHist]);
    hProjB[iHist] -> SetTitle(sTitle.Data());
    hProjB[iHist] -> SetTitleFont(fTxt);
    hProjB[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hProjB[iHist] -> GetXaxis() -> SetLabelSize(fLbl);
    hProjB[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hProjB[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hProjB[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hProjB[iHist] -> GetXaxis() -> SetTitleOffset(fOffX);
    hProjB[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hProjB[iHist] -> GetYaxis() -> SetLabelSize(fLbl);
    hProjB[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hProjB[iHist] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hProjB[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hProjB[iHist] -> GetYaxis() -> SetTitleOffset(fOffY);
    hProjB[iHist] -> GetYaxis() -> CenterTitle(fCnt);
    hRatio[iHist] -> SetFillColor(fColA[iHist]);
    hRatio[iHist] -> SetFillStyle(fFil);
    hRatio[iHist] -> SetLineColor(fColA[iHist]);
    hRatio[iHist] -> SetLineStyle(fLin);
    hRatio[iHist] -> SetMarkerColor(fColA[iHist]);
    hRatio[iHist] -> SetMarkerStyle(fMarA[iHist]);
    hRatio[iHist] -> SetTitle(sTitle.Data());
    hRatio[iHist] -> SetTitleFont(fTxt);
    hRatio[iHist] -> GetXaxis() -> SetRangeUser(xPlot[0], xPlot[1]);
    hRatio[iHist] -> GetXaxis() -> SetLabelSize(fLblR);
    hRatio[iHist] -> GetXaxis() -> SetLabelFont(fTxt);
    hRatio[iHist] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hRatio[iHist] -> GetXaxis() -> SetTitleFont(fTxt);
    hRatio[iHist] -> GetXaxis() -> SetTitleSize(fTtlR);
    hRatio[iHist] -> GetXaxis() -> SetTitleOffset(fOffXR);
    hRatio[iHist] -> GetXaxis() -> CenterTitle(fCnt);
    hRatio[iHist] -> GetYaxis() -> SetLabelSize(fLblR);
    hRatio[iHist] -> GetYaxis() -> SetLabelFont(fTxt);
    hRatio[iHist] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hRatio[iHist] -> GetYaxis() -> SetTitleFont(fTxt);
    hRatio[iHist] -> GetYaxis() -> SetTitleSize(fTtlR);
    hRatio[iHist] -> GetYaxis() -> SetTitleOffset(fOffYR);
    hRatio[iHist] -> GetYaxis() -> CenterTitle(fCnt);
  }
  cout << "    Set styles." << endl;

  // make legend and text box
  const UInt_t  fColL(0);
  const UInt_t  fLinL(0);
  const UInt_t  fFilL(0);
  const UInt_t  fAlnL(12);
  const Float_t xyLeg[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyLeA[NVtx] = {0.1, 0.1, 0.3, 0.3};
  const Float_t xyLeB[NVtx] = {0.1, 0.3, 0.3, 0.5};
  const Float_t xyTxt[NVtx] = {0.3, 0.1, 0.5, 0.3};

  TLegend *legA = new TLegend(xyLeA[0], xyLeA[1], xyLeA[2], xyLeA[3], sLegA.Data());
  TLegend *legB = new TLegend(xyLeB[0], xyLeB[1], xyLeB[2], xyLeB[3], sLegB.Data());
  legA -> SetFillColor(fColL);
  legA -> SetFillStyle(fFilL);
  legA -> SetLineColor(fColL);
  legA -> SetLineStyle(fLinL);
  legA -> SetTextFont(fTxt);
  legA -> SetTextAlign(fAlnL);
  legB -> SetFillColor(fColL);
  legB -> SetFillStyle(fFilL);
  legB -> SetLineColor(fColL);
  legB -> SetLineStyle(fLinL);
  legB -> SetTextFont(fTxt);
  legB -> SetTextAlign(fAlnL);

  TLegend *leg[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    // create label
    TString sLeg = sLegPrefix;
    sLeg.Append(sLegBin[iHist]);
    sLeg.Append(sLegSuffix);

    // create legends
    leg[iHist] = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3], sLeg.Data());
    leg[iHist] -> SetFillColor(fColL);
    leg[iHist] -> SetFillStyle(fFilL);
    leg[iHist] -> SetLineColor(fColL);
    leg[iHist] -> SetLineStyle(fLinL);
    leg[iHist] -> SetTextFont(fTxt);
    leg[iHist] -> SetTextAlign(fAlnL);
    leg[iHist] -> AddEntry(hProjA[iHist], sLegA.Data());
    leg[iHist] -> AddEntry(hProjB[iHist], sLegB.Data());
    legA       -> AddEntry(hProjA[iHist], sLeg.Data());
    legB       -> AddEntry(hProjB[iHist], sLeg.Data());
  }

  TPaveText *ptTxt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  ptTxt -> SetFillColor(fColL);
  ptTxt -> SetFillStyle(fFilL);
  ptTxt -> SetLineColor(fColL);
  ptTxt -> SetLineStyle(fLinL);
  ptTxt -> SetTextFont(fTxt);
  ptTxt -> SetTextAlign(fAlnL);
  ptTxt -> AddText(sSys.Data());
  ptTxt -> AddText(sTrg.Data());
  ptTxt -> AddText(sJet.Data());
  ptTxt -> AddText(sTyp.Data());
  cout << "    Made legends and text box." << endl;

  // make line
  const UInt_t  fColL(1);
  const UInt_t  fLinL(9);
  const UInt_t  fWidL(1);
  const Float_t yVal(1.);

  TLine *line = new TLine(xPlot[0], yVal, xPlot[1], yVal);
  line -> SetLineColor(fColL);
  line -> SetLineStyle(fLinL);
  line -> SetLineWidth(fWidL);
  cout << "    Made line." << endl;

  // make plots
  const UInt_t  fWidth(750);
  const UInt_t  fHeight(950);
  const UInt_t  fLogY(1);
  const UInt_t  fLogYR(0);
  const UInt_t  fGrid(0);
  const Float_t fMargin(0.02);
  const Float_t fMarginS(0.005);
  const Float_t fPadXY1[NVtx] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NVtx] = {0., 0.35, 1., 1.};
  const TString sPadProj("pProj");
  const TString sPadRatio("pRatio");

  TCanvas *cProjPt[NHist];
  TPad    *pProjPt[NHist];
  TPad    *pRatioPt[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    // make pad names
    TString sPadP = sPadProj;
    TString sPadR = sPadRatio;
    sPadP.Append(sPadBin[iHist].Data());
    sPadR.Append(sPadBin[iHist].Data());

    // make plot
    cProjPt[iHist]  = new TCanvas(sCan[iHist].Data(), "", fWidth, fHeight);
    pProjPt[iHist]  = new TPad(sPadP.Data(), "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
    pRatioPt[iHist] = new TPad(sPadR.Data(), "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
    pProjPt[iHist]  -> SetGrid(fGrid, fGrid);
    pProjPt[iHist]  -> SetLogy(fLogY);
    pProjPt[iHist]  -> SetTopMargin(fMargin);
    pProjPt[iHist]  -> SetRightMargin(fMargin);
    pProjPt[iHist]  -> SetBottomMargin(fMarginS);
    pRatioPt[iHist] -> SetGrid(fGrid, fGrid);
    pRatioPt[iHist] -> SetLogy(fLogYR);
    pRatioPt[iHist] -> SetTopMargin(fMarginS);
    pRatioPt[iHist] -> SetRightMargin(fMargin);
    cProjPt[iHist]  -> cd();
    pProjPt[iHist]  -> Draw();
    pRatioPt[iHist] -> Draw();
    pProjPt[iHist]  -> cd();
    hProjA[iHist]   -> Draw("hist p l");
    hProjB[iHist]   -> Draw("hist p l sames");
    leg[iHist]      -> Draw();
    ptTxt           -> Draw();
    pRatioPt[iHist] -> cd();
    hRatio[iHist]   -> Draw("hist p");
    line            -> Draw();
    fOut            -> cd();
    cProjPt[iHist]  -> Write();
    cProjPt[iHist]  -> Close();
  }

  TCanvas *cProj  = new TCanvas("cAllProj", "", fWidth, fHeight);
  TPad    *pProj  = new TPad(sPadProj.Data(), "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  TPad    *pRatio = new TPad(sPadRatio.Data(), "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]); 
  pProj     -> SetGrid(fGrid, fGrid);
  pProj     -> SetLogy(fLogY);
  pProj     -> SetTopMargin(fMargin);
  pProj     -> SetRightMargin(fMargin);
  pProj     -> SetBottomMargin(fMarginS);
  pRatio    -> SetGrid(fGrid, fGrid);
  pRatio    -> SetLogy(fLogYR);
  pRatio    -> SetTopMargin(fMarginS);
  pRatio    -> SetRightMargin(fMargin);
  cProj     -> cd();
  pProj     -> Draw();
  pRatio    -> Draw();
  pProj     -> cd();
  hProjA[0] -> Draw("hist p l");
  hProjB[0] -> Draw("hist p l sames");
  for (UInt_t iHist = 1; iHist < NHist; iHist++) { 
    hProjA[iHist] -> Draw("hist p l sames");
    hProjB[iHist] -> Draw("hist p l sames");
  }
  legA      -> Draw();
  legB      -> Draw();
  pRatio    -> cd();
  hRatio[0] -> Draw("hist p");
  for (UInt_t iHist = 1; iHist < NHist; iHist++) {
    hRatio[iHist] -> Draw("hist p same");
  }
  ptTxt -> Draw();
  line  -> Draw();
  fOut  -> cd();
  cProj -> Write();
  cProj -> Close();
  cout << "    Made plots." << endl;

  // save histograms
  fOut -> cd();
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    hProjA[iHist] -> Write();
    hProjB[iHist] -> Write();
    hRatio[iHist] -> Write();
  }

  // close files
  fOut -> cd();
  fOut -> Close();
  fInA -> cd();
  fInA -> Close();
  fInB -> cd();
  fInB -> Close();
  cout << "  Finished response matrix projection comparison.\n" << endl;

}

// End ------------------------------------------------------------------------
