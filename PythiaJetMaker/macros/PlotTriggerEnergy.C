// 'PlotTriggerEnergy.C'
// Derek Anderson
// 02.26.2018
//
// Use this to quickly plot
// the trigger transverse
// energy and compare between
// pi0 and gamma triggers.


#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;


// io parameters
static const TString sOutput("pp200py.eTtrgAvg920.d26m2y2018.root");
static const TString sInPi0("input/py23pi01par.merged.root");
static const TString sInGam("input/py23gam1par.merged.root");
static const TString sChain("ParTree");
static const Bool_t  doEtCut(true);



void PlotTriggerEnergy() {

  gErrorIgnoreLevel = kError;
  cout << "\n  Beginning plotting script..." << endl;


  // open files
  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInPi0  = new TFile(sInPi0.Data(), "read");
  TFile *fInGam  = new TFile(sInGam.Data(), "read");
  if (!fInPi0 || !fInGam) {
    cerr << "PANIC: couldn't open an input file!" << endl;
    return;
  }
  cout << "    Opened files." << endl;

  // grab trees
  TChain *tPi0 = (TChain*) fInPi0 -> Get(sChain.Data());
  TChain *tGam = (TChain*) fInGam -> Get(sChain.Data());
  if (!tPi0 || !tGam) {
    cerr << "PANIC: couldn't grab input tree!" << endl;
    return;
  }
  cout << "    Grabbed trees." << endl;


  // make histograms
  const UInt_t  nBins(16);
  const Float_t bin1(7.);
  const Float_t bin2(23.);
  const TString sPi0("hPi0");
  const TString sGam("hGam");

  TH1D *hPi0 = new TH1D(sPi0.Data(), "", nBins, bin1, bin2);
  TH1D *hGam = new TH1D(sGam.Data(), "", nBins, bin1, bin2);
  hPi0 -> Sumw2();
  hGam -> Sumw2();
  cout << "    Made histograms." << endl;


  // fill histograms
  const TString sLeaf("Events_Clust_EneT0 * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * Events_Clust_etav1)))");
  const TString sCuts("TMath::Abs(Events_Clust_etav1) < 0.9");

  // add eT cut
  TString sCutsPi(sCuts.Data());
  TString sCutsGa(sCuts.Data());
  TString sEtCut(" && (((");
  sEtCut += sLeaf.Data();
  sEtCut += ") > 9.) && ((";
  sEtCut += sLeaf.Data();
  sEtCut += ") < 20.))";
  if (doEtCut) {
    sCutsPi.Append(sEtCut.Data());
    sCutsGa.Append(sEtCut.Data());
  }
  cout << "    Cuts:\n"
       << "      pi0 = " << sCutsPi.Data() << "\n"
       << "      gam = " << sCutsGa.Data()
       << endl; 

  TString sDrawPi(sLeaf.Data());
  TString sDrawGa(sLeaf.Data());
  sDrawPi.Append(">>");
  sDrawGa.Append(">>");
  sDrawPi.Append(sPi0.Data());
  sDrawGa.Append(sGam.Data());

  const Long64_t nPi0  = tPi0 -> Draw(sDrawPi.Data(), sCutsPi.Data());
  const Long64_t nGam  = tGam -> Draw(sDrawGa.Data(), sCutsGa.Data());
  const Double_t piAvg = hPi0 -> GetMean();
  const Double_t gaAvg = hGam -> GetMean();
  cout << "    Filled histograms:\n"
       << "      nPi0 = " << nPi0 << ", nGam = " << nGam << "\n"
       << "      <eTpi0> = " << piAvg << ", <eTgam> = " << gaAvg
       << endl;


  // normalize histograms
  const Double_t dPi0 = hPi0 -> GetBinWidth(17);
  const Double_t dGam = hGam -> GetBinWidth(17);
  hPi0 -> Scale(1. / (Double_t) nPi0);
  hPi0 -> Scale(1. / dPi0);
  hGam -> Scale(1. / (Double_t) nGam);
  hGam -> Scale(1. / dGam);
  cout << "    Normalized histograms." << endl;


  // set styles
  const UInt_t  fColP(858);
  const UInt_t  fColG(898);
  const UInt_t  fMarP(4);
  const UInt_t  fMarG(4);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.3);
  const Float_t fOff(0.9);
  const TString sTitle("Trigger E_{T}");
  const TString sTitleX("E_{T}^{trg}");
  const TString sTitleY("arb. units");
  hPi0 -> SetLineColor(fColP);
  hPi0 -> SetMarkerColor(fColP);
  hPi0 -> SetMarkerStyle(fMarP);
  hPi0 -> SetTitle(sTitle.Data());
  hPi0 -> SetTitleFont(fTxt);
  hPi0 -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPi0 -> GetXaxis() -> SetTitleFont(fTxt);
  hPi0 -> GetXaxis() -> SetTitleOffset(fOff);
  hPi0 -> GetXaxis() -> CenterTitle(fCnt);
  hPi0 -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPi0 -> GetYaxis() -> SetTitleFont(fTxt);
  hPi0 -> GetYaxis() -> SetTitleOffset(fOff);
  hPi0 -> GetYaxis() -> CenterTitle(fCnt);
  hGam -> SetLineColor(fColG);
  hGam -> SetMarkerColor(fColG);
  hGam -> SetMarkerStyle(fMarG);
  hGam -> SetTitle(sTitle.Data());
  hGam -> SetTitleFont(fTxt);
  hGam -> GetXaxis() -> SetTitle(sTitleX.Data());
  hGam -> GetXaxis() -> SetTitleFont(fTxt);
  hGam -> GetXaxis() -> SetTitleOffset(fOff);
  hGam -> GetXaxis() -> CenterTitle(fCnt);
  hGam -> GetYaxis() -> SetTitle(sTitleY.Data());
  hGam -> GetYaxis() -> SetTitleFont(fTxt);
  hGam -> GetYaxis() -> SetTitleOffset(fOff);
  hGam -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;


  // make legend
  const UInt_t  fColL(0);
  const Float_t xyLeg[4] = {0.1, 0.1, 0.3, 0.3};
  const TString sLegP("#pi^{0}");
  const TString sLegG("#gamma^{dir}");
  TLegend *leg = new TLegend(xyLeg[0], xyLeg[1], xyLeg[2], xyLeg[3]);
  leg -> SetFillColor(fColL);
  leg -> SetLineColor(fColL);
  leg -> SetTextFont(fTxt);
  leg -> AddEntry(hPi0, sLegP.Data());
  leg -> AddEntry(hGam, sLegG.Data());
  cout << "    Made legend." << endl;


  // make text
  const UInt_t  nDec(2);
  const UInt_t  fColT(0);
  const UInt_t  fAlign(12);
  const Float_t xyTxt[4] = {0.3, 0.1, 0.5, 0.3};
  const TString sSystem("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sData("Pythia 8");
  const TString sPiAvg("#LTE_{T}^{trg}#GT_{#pi} = ");
  const TString sGaAvg("#LTE_{T}^{trg}#GT_{#gamma} = ");

  TString sPiRaw("");
  TString sGaRaw("");
  TString sPiTxt(sPiAvg.Data());
  TString sGaTxt(sGaAvg.Data());
  sPiRaw += piAvg;
  sGaRaw += gaAvg;

  const UInt_t nRawP = sPiRaw.First(".");
  const UInt_t nRawG = sGaRaw.First(".");
  const UInt_t nTxtP = (nRawP + nDec) + 1;
  const UInt_t nTxtG = (nRawG + nDec) + 1;
  sPiTxt.Append(sPiRaw.Data(), nTxtP);
  sGaTxt.Append(sGaRaw.Data(), nTxtG);

  TString sAvgs("#color[");
  sAvgs += fColP;
  sAvgs += "]{";
  sAvgs += sPiTxt.Data();
  sAvgs += "}, #color[";
  sAvgs += fColG;
  sAvgs += "]{";
  sAvgs += sGaTxt.Data();
  sAvgs += "}";

  TPaveText *txt = new TPaveText(xyTxt[0], xyTxt[1], xyTxt[2], xyTxt[3], "NDC NB");
  txt -> SetFillColor(fColT);
  txt -> SetLineColor(fColT);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlign);
  txt -> AddText(sSystem.Data());
  txt -> AddText(sData.Data());
  txt -> AddText(sAvgs.Data());
  cout << "    Made text." << endl;


  // make plot
  const UInt_t  width(750);
  const UInt_t  height(750);
  const UInt_t  fLog(1);
  const UInt_t  fGrid(0);
  const UInt_t  fTicks(1);
  const TString sCanvas("cEtTrg");
  fOutput -> cd();

  TCanvas *cEtTrg = new TCanvas(sCanvas.Data(), "", width, height);
  cEtTrg -> SetGrid(fGrid, fGrid);
  cEtTrg -> SetTicks(fTicks, fTicks);
  cEtTrg -> SetLogy(fLog);
  cEtTrg -> cd();
  hPi0   -> Draw();
  hGam   -> Draw("same");
  leg    -> Draw();
  txt    -> Draw();
  cEtTrg -> Write();
  cEtTrg -> Close();
  cout << "    Made plot." << endl;


  // close files
  fOutput -> cd();
  hPi0    -> Write();
  hGam    -> Write();
  fOutput -> Close();
  fInPi0  -> cd();
  fInPi0  -> Close();
  fInGam  -> cd();
  fInGam  -> Close();
  cout << "  Finished script!\n" << endl;

}

// End ------------------------------------------------------------------------
