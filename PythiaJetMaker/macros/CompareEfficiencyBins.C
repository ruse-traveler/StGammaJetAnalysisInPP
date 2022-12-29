// 'CompareEfficiencyBins.C'
// Derek Anderson
// 07.04.2017
//
// Use this to compare particle-level Pythia,
// detector-level Pythia (split up into different
// pT bins), and data.


#include <cassert>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// i/o parameters
static const UInt_t  nBins(3);
static const TString sOutput("pythiaXdata.dFjetGam.eff87.r03a02rm1chrg.d4m7y2017.root");
static const TString sInputA("pythiaGP.plots.r03a02rm1chrg.d4m7y2017.root");
static const TString sInputB("pythiaGD.eff87.plots.r03a02rm1chrg.d4m7y2017.root");
static const TString sInputC("../JetData/pp200r9.plots.r03a02rm1chrg.d4m7y2017.root");
static const TString sHistA[nBins] = {"Gam/hBinDeltaPhiG_pT02", "Gam/hBinDeltaPhiG_pT1", "Gam/hBinDeltaPhiG_pT5"};
static const TString sHistB[nBins] = {"Gam/hBinDeltaPhiG_pT02", "Gam/hBinDeltaPhiG_pT1", "Gam/hBinDeltaPhiG_pT5"};
static const TString sHistC[nBins] = {"Gam/hBinDeltaPhiG_pT02", "Gam/hBinDeltaPhiG_pT1", "Gam/hBinDeltaPhiG_pT5"};
// histogram parameters
static const TString sTitle("Jet #Delta#varphi, #gamma^{rich} trigger");
static const TString sTitleX("#Delta#varphi = #varphi_{jet} - #varphi_{trg}");
static const TString sTitleY("(1/N_{trg}) dN_{jet}/d#Delta#varphi");
// legend parameters
static const TString sLegendA("Pythia, particle");
static const TString sLegendB("Pythia, w/ #epsilon_{0} = 87%");
static const TString sLegendC("Run 9 data");
// label parameters
static const TString sLabel1("pp-collisions, #sqrt{s} = 200 GeV");
static const TString sLabel2("E_{T}^{trg} #in (9, 30) GeV, |#eta^{trg}| < 0.9");
static const TString sLabel3("p_{T}^{trk} #in (0.2, 30) GeV/c, |#eta^{trk}| < 1.0");
static const TString sLabel4("Anti-k_{T}, R = 0.3, charged jet");
// canvas parameters
static const TString sCanvas[nBins] = {"cJetDfGam_pT02", "cJetDfGam_pT1", "cJetDfGam_pT5"};



void CompareEfficiencyBins() {

  gErrorIgnoreLevel = kFatal;
  cout << "\n  Beginning comparison script..." << endl;


  TFile *fOutput = new TFile(sOutput.Data(), "recreate");
  TFile *fInputA = new TFile(sInputA.Data(), "read");
  TFile *fInputB = new TFile(sInputB.Data(), "read");
  TFile *fInputC = new TFile(sInputC.Data(), "read");
  if (!fOutput || !fInputA || !fInputB || !fInputC) {
    cerr << "PANIC: couldn't open file!" << endl;
    assert(fOutput);
    assert(fInputA);
    assert(fInputB);
    assert(fInputC);
  }
  cout << "    Files opened." << endl;

  TH1D *hInputA[nBins];
  TH1D *hInputB[nBins];
  TH1D *hInputC[nBins];
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    hInputA[iBin] = (TH1D*) fInputA -> Get(sHistA[iBin].Data());
    hInputB[iBin] = (TH1D*) fInputB -> Get(sHistB[iBin].Data());
    hInputC[iBin] = (TH1D*) fInputC -> Get(sHistC[iBin].Data());
    if (!hInputA[iBin] || !hInputB[iBin] || !hInputC[iBin]) {
      cerr << "PANIC: couldn't grab histogram " << iBin << "!" << endl;
      assert(hInputA[iBin]);
      assert(hInputB[iBin]);
      assert(hInputC[iBin]);
    }
  }
  cout << "    Histograms grabbed." << endl;


  const Int_t    cA  = 890;
  const Int_t    cB  = 850;
  const Int_t    cC  = 810;
  const Int_t    mA  = 20;
  const Int_t    mB  = 24;
  const Int_t    mC  = 7;
  const Int_t    txt = 42;
  const Int_t    cnt = 1;
  const Double_t lab = 0.02;
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    hInputA[iBin] -> SetLineColor(cA);
    hInputA[iBin] -> SetMarkerColor(cA);
    hInputA[iBin] -> SetMarkerStyle(mA);
    hInputA[iBin] -> SetTitleFont(txt);
    hInputA[iBin] -> SetTitle(sTitle.Data());
    hInputA[iBin] -> GetXaxis() -> SetLabelSize(lab);
    hInputA[iBin] -> GetXaxis() -> CenterTitle(cnt);
    hInputA[iBin] -> GetXaxis() -> SetTitleFont(txt);
    hInputA[iBin] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hInputA[iBin] -> GetYaxis() -> SetLabelSize(lab);
    hInputA[iBin] -> GetYaxis() -> CenterTitle(cnt);
    hInputA[iBin] -> GetYaxis() -> SetTitleFont(txt);
    hInputA[iBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hInputB[iBin] -> SetLineColor(cB);
    hInputB[iBin] -> SetMarkerColor(cB);
    hInputB[iBin] -> SetMarkerStyle(mB);
    hInputB[iBin] -> SetTitleFont(txt);
    hInputB[iBin] -> SetTitle(sTitle.Data());
    hInputB[iBin] -> GetXaxis() -> SetLabelSize(lab);
    hInputB[iBin] -> GetXaxis() -> CenterTitle(cnt);
    hInputB[iBin] -> GetXaxis() -> SetTitleFont(txt);
    hInputB[iBin] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hInputB[iBin] -> GetYaxis() -> SetLabelSize(lab);
    hInputB[iBin] -> GetYaxis() -> CenterTitle(cnt);
    hInputB[iBin] -> GetYaxis() -> SetTitleFont(txt);
    hInputB[iBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hInputC[iBin] -> SetLineColor(cC);
    hInputC[iBin] -> SetMarkerColor(cC);
    hInputC[iBin] -> SetMarkerStyle(mC);
    hInputC[iBin] -> SetTitleFont(txt);
    hInputC[iBin] -> SetTitle(sTitle.Data());
    hInputC[iBin] -> GetXaxis() -> SetLabelSize(lab);
    hInputC[iBin] -> GetXaxis() -> CenterTitle(cnt);
    hInputC[iBin] -> GetXaxis() -> SetTitleFont(txt);
    hInputC[iBin] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hInputC[iBin] -> GetYaxis() -> SetLabelSize(lab);
    hInputC[iBin] -> GetYaxis() -> CenterTitle(cnt);
    hInputC[iBin] -> GetYaxis() -> SetTitleFont(txt);
    hInputC[iBin] -> GetYaxis() -> SetTitle(sTitleY.Data());
  }
  cout << "    Styles set." << endl;


  const Int_t    cL  = 0;
  const Int_t    fL  = 0;
  const Int_t    sL  = 0;
  const Double_t x1L = 0.1;
  const Double_t x2L = 0.3;
  const Double_t y1L = 0.1;
  const Double_t y2L = 0.3;
  TLegend *lLegend[nBins];
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    lLegend[iBin] = new TLegend(x1L, y1L, x2L, y2L);
    lLegend[iBin] -> SetFillColor(cL);
    lLegend[iBin] -> SetFillStyle(sL);
    lLegend[iBin] -> SetLineColor(cL);
    lLegend[iBin] -> SetLineStyle(sL);
    lLegend[iBin] -> SetTextFont(txt);
    lLegend[iBin] -> AddEntry(hInputA[iBin], sLegendA.Data());
    lLegend[iBin] -> AddEntry(hInputB[iBin], sLegendB.Data());
    lLegend[iBin] -> AddEntry(hInputC[iBin], sLegendC.Data());
  }
  cout << "    Legends created." << endl;

  const Double_t x1P = 0.3;
  const Double_t x2P = 0.5;
  const Double_t y1P = 0.1;
  const Double_t y2P = 0.3;
  TPaveText *pLabel = new TPaveText(x1P, y1P, x2P, y2P, "NDC NB");
  pLabel -> SetFillColor(cL);
  pLabel -> SetFillStyle(sL);
  pLabel -> SetLineColor(cL);
  pLabel -> SetLineStyle(sL);
  pLabel -> SetTextFont(txt);
  pLabel -> AddText(sLabel1.Data());
  pLabel -> AddText(sLabel2.Data());
  pLabel -> AddText(sLabel3.Data());
  pLabel -> AddText(sLabel4.Data());
  cout << "    Label created." << endl;


  // make plots
  fOutput -> cd();

  const Int_t wC  = 800;
  const Int_t hC  = 800;
  const Int_t grd = 0;
  const Int_t log = 0;
  TCanvas *cPlot[nBins];
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    cPlot[iBin] = new TCanvas(sCanvas[iBin].Data(), "", wC, hC);
    cPlot[iBin]   -> SetGrid(grd, grd);
    hInputA[iBin] -> Draw();
    hInputB[iBin] -> Draw("same");
    hInputC[iBin] -> Draw("same");
    lLegend[iBin] -> Draw();
    pLabel        -> Draw();
    cPlot[iBin]   -> Write();
    cPlot[iBin]   -> Close();
  }
  cout << "    Plots drawn." << endl;


  fOutput -> cd();
  fOutput -> Close();
  fInputA -> cd();
  fInputA -> Close();
  fInputB -> cd();
  fInputB -> Close();
  fInputC -> cd();
  fInputC -> Close();
  cout << "  Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
