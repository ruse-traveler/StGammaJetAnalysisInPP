// 'CompareSpectra.C'
// Derek Anderson
// 12.08.206
//
// This takes various spectra and makes
// some comparison plots.

#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// i/o paths
static const TString sOut("./PythiaVsPpVsAuAu6080.Dec12.root");
static const TString sData("./plots/pp200run9.r03rm1chrg.plots.Dec8.root");
static const TString sPythPG("../JetMacro/Pythia20p.r03a02rm1.gPlots.Dec9.root");
static const TString sPythPP("../JetMacro/Pythia23p.r03a02rm1.pPlots.Dec8.root");
static const TString sPythDG("../JetMacro/Pythia20d.r03a02rm1.gPlots.Dec9.root");
static const TString sPythDP("../JetMacro/Pythia23d.r03a02rm1.pPlots.Dec8.root");
static const TString sAuAu("./input/JetPt.AuAu6080.ChargedHadron9GeV.r03a02rm3.txt");
static const TString sPythAA("./input/JetPt.Pythia.ChargedHadron9GeV.r03a02rm3.txt");
// plot parameters
static const Double_t pTlo = 0.;
static const Double_t yLo  = 0.00001;
static const Double_t pThi = 27.;
static const Double_t yHi  = 5.;



void CompareSpectra() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;

  cout << "\nBeggining comparison script..." << endl;


  TFile *fOut    = new TFile(sOut.Data(), "recreate");
  TFile *fData   = new TFile(sData.Data(), "read");
  TFile *fPythPG = new TFile(sPythPG.Data(), "read");
  TFile *fPythPP = new TFile(sPythPP.Data(), "read");
  TFile *fPythDG = new TFile(sPythDG.Data(), "read");
  TFile *fPythDP = new TFile(sPythDP.Data(), "read");

  // make sure files are open
  Bool_t datIsOpen = (Bool_t) fData;
  Bool_t gamIsOpen = (Bool_t) (fPythPG && fPythDG);
  Bool_t pi0IsOpen = (Bool_t) (fPythPP && fPythDP);
  if (!datIsOpen) {
    cerr << "PANIC: couldn't open data input!" << endl;
    assert(datIsOpen);
  }
  if (!gamIsOpen) {
    cerr << "PANIC: couldn't open Pythia gamma input!" << endl;
    assert(gamIsOpen);
  }
  if (!pi0IsOpen) {
    cerr << "PANIC: couldn't open Pythia pi0 input!" << endl;
    assert(pi0IsOpen);
  }

  cout << "  Files opened..." << endl;

  // grab histograms
  TH1D *hDataG, *hDataP;
  TH1D *hPythPG, *hPythDG;
  TH1D *hPythPP, *hPythDP;
  if (datIsOpen) {
    TH1D *hDg, *hDp;
    fData -> GetObject("Gam/hJetPtRawG", hDg);
    fData -> GetObject("Pi0/hJetPtRawP", hDp);
    hDataG = (TH1D*) hDg -> Clone();
    hDataP = (TH1D*) hDp -> Clone();
    hDataG -> SetNameTitle("hDataG", "Jet p_{T}, #gamma^{rich} (pp)");
    hDataP -> SetNameTitle("hDataP", "Jet p_{T}, #pi^{0} (pp)");
  }
  if (gamIsOpen) {
    TH1D *hPpg, *hPdg;
    fPythPG -> GetObject("Gam/hJetPtRawG", hPpg);
    fPythDG -> GetObject("Gam/hJetPtRawG", hPdg);
    hPythPG = (TH1D*) hPpg -> Clone();
    hPythDG = (TH1D*) hPdg -> Clone();
    hPythPG -> SetNameTitle("hPythPG", "Jet p_{T}, #gamma^{dir} (particle-pythia)");
    hPythDG -> SetNameTitle("hPythDG", "Jet p_{T}, #gamma^{dir} (detector-pythia)");
  }
  if (pi0IsOpen) {
    TH1D *hPpp, *hPdp;
    fPythPP -> GetObject("Pi0/hJetPtRawP", hPpp);
    fPythDP -> GetObject("Pi0/hJetPtRawP", hPdp);
    hPythPP = (TH1D*) hPpp -> Clone();
    hPythDP = (TH1D*) hPdp -> Clone();
    hPythPP -> SetNameTitle("hPythPP", "Jet p_{T}, #pi^{0} (particle-pythia)");
    hPythDP -> SetNameTitle("hPythDP", "Jet p_{T}, #pi^{0} (detector-pythia)");
  }
  fOut -> cd();


  // open AuAu stream
  ifstream auau(sAuAu.Data());
  if (!auau) {
    cerr << "PANIC: input stream couldn't be opened!" << endl;
    assert(auau);
  }

  vector<Double_t> binA;
  vector<Double_t> valA;
  vector<Double_t> eXLA;
  vector<Double_t> eXHA;
  vector<Double_t> eYLA;
  vector<Double_t> eYHA;

  // read file
  Double_t xA  = 0.;
  Double_t yA  = 0.;
  Double_t xLA = 0.;
  Double_t xHA = 0.;
  Double_t yLA = 0.;
  Double_t yHA = 0.;
  while (auau) {
    auau >> xA;
    auau >> yA;
    auau >> xLA;
    auau >> xHA;
    auau >> yLA;
    auau >> yHA;
    binA.push_back(xA);
    valA.push_back(yA);
    eXLA.push_back(xLA);
    eXHA.push_back(xHA);
    eYLA.push_back(yLA);
    eYHA.push_back(yHA);
  }
  auau.close();

  // create histogram
  const Int_t nBinsA = ((Int_t) binA.size()) + 2;
  TH1D *hAuAu = new TH1D("hAuAu", "Jet p_{T}, AuAu 60-80%", nBinsA, pTlo, pThi);
  for (Int_t i = 0; i < nBinsA; i++) {
    Int_t iBin = hAuAu -> FindBin(binA[i]);
    hAuAu -> SetBinContent(iBin, valA[i]);
    hAuAu -> SetBinError(iBin, eYHA[i]);
  }


  // open Pythia (Alex) stream
  ifstream pyth(sPythAA.Data());
  if (!pyth) {
    cerr << "PANIC: input stream couldn't be opened!" << endl;
    assert(pyth);
  }

  vector<Double_t> binP;
  vector<Double_t> valP;

  // read file
  Double_t xP = 0.;
  Double_t yP = 0.;
  while (pyth) {
    pyth >> xP;
    pyth >> yP;
    binP.push_back(xP);
    valP.push_back(yP);
  }
  pyth.close();

  // create histogram
  const Int_t nBinsP = binP.size() + 1; 
  TH1D *hPythAA = new TH1D("hPythAA", "Jet p_{T}, Pythia (Alex)", nBinsP, pTlo, pThi);
  for (Int_t i = 0; i < nBinsP; i++) {
    Int_t iBin = hPythAA -> FindBin(binP[i]);
    hPythAA -> SetBinContent(iBin, valP[i]);
  }
  hPythAA -> SetMarkerStyle(4);

  cout << "  Histograms grabbed..." << endl;


  // set styles
  hAuAu   -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hDataG  -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hDataP  -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hPythPG -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hPythPP -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hPythDG -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hPythDP -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hPythAA -> GetXaxis() -> SetRangeUser(pTlo, pThi);
  hAuAu   -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hDataG  -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hDataP  -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hPythPG -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hPythPP -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hPythDG -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hPythDP -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hPythAA -> GetYaxis() -> SetRangeUser(yLo, yHi);
  hAuAu   -> SetMarkerColor(kBlack);
  hDataG  -> SetMarkerColor(kBlack);
  hDataP  -> SetMarkerColor(kBlack);
  hPythPG -> SetMarkerColor(kBlue);
  hPythPP -> SetMarkerColor(kBlue);
  hPythDG -> SetMarkerColor(kRed);
  hPythDP -> SetMarkerColor(kRed);
  hPythAA -> SetMarkerColor(kViolet);
  hAuAu   -> SetLineColor(kBlack);
  hDataG  -> SetLineColor(kBlack);
  hDataP  -> SetLineColor(kBlack);
  hPythPG -> SetLineColor(kBlue);
  hPythPP -> SetLineColor(kBlue);
  hPythDG -> SetLineColor(kRed);
  hPythDP -> SetLineColor(kRed);
  hPythAA -> SetLineColor(kViolet);


  // create legends
  TLegend *lDgVsPg = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lDgVsPp = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lDpVsPp = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lAuVsPg = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lAuVsPp = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lPaVsPg = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lPaVsPp = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lAuVsG  = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lAuVsP  = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lPpVsG  = new TLegend(0.1, 0.1, 0.3, 0.3);
  TLegend *lPpVsP  = new TLegend(0.1, 0.1, 0.3, 0.3);
  lDgVsPg -> SetLineColor(kWhite);
  lDgVsPp -> SetLineColor(kWhite);
  lDpVsPp -> SetLineColor(kWhite);
  lAuVsPg -> SetLineColor(kWhite);
  lAuVsPp -> SetLineColor(kWhite);
  lPaVsPg -> SetLineColor(kWhite);
  lPaVsPp -> SetLineColor(kWhite);
  lAuVsG  -> SetLineColor(kWhite);
  lPpVsP  -> SetLineColor(kWhite);
  lPpVsG  -> SetLineColor(kWhite);
  lAuVsP  -> SetLineColor(kWhite);
  lDgVsPg -> SetFillColor(kWhite);
  lDgVsPp -> SetFillColor(kWhite);
  lDpVsPp -> SetFillColor(kWhite);
  lAuVsPg -> SetFillColor(kWhite);
  lAuVsPp -> SetFillColor(kWhite);
  lPaVsPg -> SetFillColor(kWhite);
  lPaVsPp -> SetFillColor(kWhite);
  lAuVsG  -> SetFillColor(kWhite);
  lAuVsP  -> SetFillColor(kWhite);
  lPpVsG  -> SetFillColor(kWhite);
  lPpVsP  -> SetFillColor(kWhite);

  lDgVsPg -> AddEntry(hPythPG, "Pythia, particle");
  lDgVsPg -> AddEntry(hPythDG, "Pythia, detector");
  lDgVsPg -> AddEntry(hDataG, "pp data, 200 GeV");
  lDgVsPp -> AddEntry(hPythPP, "Pythia, particle");
  lDgVsPp -> AddEntry(hPythDP, "Pythia, detector");
  lDgVsPp -> AddEntry(hDataG, "pp data, 200 GeV");
  lDpVsPp -> AddEntry(hPythPP, "Pythia, particle");
  lDpVsPp -> AddEntry(hPythDP, "Pythia, detector");
  lDpVsPp -> AddEntry(hDataP, "pp data, 200 GeV");
  lAuVsPg -> AddEntry(hPythPG, "Pythia, particle");
  lAuVsPg -> AddEntry(hPythDG, "Pythia, detector");
  lAuVsPg -> AddEntry(hAuAu, "AuAu data, 60-80%");
  lAuVsPp -> AddEntry(hPythPP, "Pythia, particle");
  lAuVsPp -> AddEntry(hPythDP, "Pythia, detector");
  lAuVsPp -> AddEntry(hAuAu, "AuAu data, 60-80%");
  lPaVsPg -> AddEntry(hPythPG, "Pythia, particle");
  lPaVsPg -> AddEntry(hPythDG, "Pythia, detector");
  lPaVsPg -> AddEntry(hPythAA, "Pythia, h^{#pm}");
  lPaVsPp -> AddEntry(hPythPP, "Pythia, particle");
  lPaVsPp -> AddEntry(hPythDP, "Pythia, detector");
  lPaVsPp -> AddEntry(hPythAA, "Pythia, h^{#pm}");
  lAuVsG  -> AddEntry(hPythPG, "Pythia, particle");
  lAuVsG  -> AddEntry(hPythDG, "Pythia, detector");
  lAuVsG  -> AddEntry(hPythAA, "Pythia, h^{#pm}");
  lAuVsG  -> AddEntry(hAuAu, "AuAu data, 60-80%");
  lAuVsP  -> AddEntry(hPythPP, "Pythia, particle");
  lAuVsP  -> AddEntry(hPythDP, "Pythia, detector");
  lAuVsP  -> AddEntry(hPythAA, "Pythia, h^{#pm}");
  lAuVsP  -> AddEntry(hAuAu, "AuAu data, 60-80%");
  lPpVsG  -> AddEntry(hPythPG, "Pythia, particle");
  lPpVsG  -> AddEntry(hPythDG, "Pythia, detector");
  lPpVsG  -> AddEntry(hPythAA, "Pythia, h^{#pm}");
  lPpVsG  -> AddEntry(hDataP, "pp data, 200 GeV");
  lPpVsP  -> AddEntry(hPythPP, "Pythia, particle");
  lPpVsP  -> AddEntry(hPythDP, "Pythia, detector");
  lPpVsP  -> AddEntry(hPythAA, "Pythia, h^{#pm}");
  lPpVsP  -> AddEntry(hDataP, "pp data, 200 GeV");


  // create labels
  TPaveText *pDatGam = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pGamPi0 = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pDatPi0 = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pHadGam = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pHadPi0 = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pDatHad = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pPythG  = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pPythP  = new TPaveText(0.3, 0.3, 0.5, 0.5, "NDC NB");
  TPaveText *pJets   = new TPaveText(0.5, 0.5, 0.7, 0.7, "NDC NB");
  pDatGam -> SetFillColor(kWhite);
  pDatPi0 -> SetFillColor(kWhite);
  pGamPi0 -> SetFillColor(kWhite);
  pHadGam -> SetFillColor(kWhite);
  pHadPi0 -> SetFillColor(kWhite);
  pDatHad -> SetFillColor(kWhite);
  pPythG  -> SetFillColor(kWhite);
  pPythP  -> SetFillColor(kWhite);
  pJets   -> SetFillColor(kWhite);

  pDatGam -> AddText("Pythia: #gamma^{dir}");
  pDatGam -> AddText("Data: #gamma^{rich}");
  pGamPi0 -> AddText("Pythia: #pi^{0}");
  pGamPi0 -> AddText("Data: #gamma^{rich}");
  pDatPi0 -> AddText("Pythia: #pi^{0}");
  pDatPi0 -> AddText("Data: #pi^{0}");
  pHadGam -> AddText("Pythia: #gamma^{dir}");
  pHadGam -> AddText("Data: charged hadron");
  pHadPi0 -> AddText("Pythia: #pi^{0}");
  pHadPi0 -> AddText("Data: charged hadron");
  pDatHad -> AddText("pp data: #pi^{0}");
  pDatHad -> AddText("AA data: charged hadron");
  pPythG  -> AddText("Pythia (TAMU): #gamma^{dir}");
  pPythG  -> AddText("Pythia (LBNL): h^{#pm}");
  pPythP  -> AddText("Pythia (TAMU): #pi^{0}");
  pPythP  -> AddText("Pythia (LBNL): h^{#pm}");
  pJets   -> AddText("Anti-k_{T}, R = 0.3");
  pJets   -> AddText("A > 0.2, p_{T}^{cst} > 0.2");

  cout << "  Labels created..." << endl;


  // draw plots
  const Double_t width  = 700;
  const Double_t height = 500;
  TCanvas *cDgVsPg = new TCanvas("cDgVsPg", "Gamma-data vs. gamma-sim.", width, height);
  cDgVsPg -> cd();
  cDgVsPg -> SetLogy(1);
  cDgVsPg -> SetGrid(0, 0);
  hPythPG -> SetTitle("Pythia vs. pp-data, #gamma-trigger; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPG -> Draw();
  hPythDG -> Draw("same");
  hDataG  -> Draw("same");
  lDgVsPg -> Draw();
  pDatGam -> Draw();
  pJets   -> Draw();
  cDgVsPg -> Write();
  cDgVsPg -> Close();

  TCanvas *cDgVsPp = new TCanvas("cDgVsPp", "Gamma-data vs. pi0-sim.", width, height);
  cDgVsPp -> cd();
  cDgVsPp -> SetLogy(1);
  cDgVsPp -> SetGrid(0, 0);
  hPythPP -> SetTitle("Pythia vs. pp-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPP -> Draw();
  hPythDP -> Draw("same");
  hDataG  -> Draw("same");
  lDgVsPp -> Draw();
  pGamPi0 -> Draw();
  pJets   -> Draw();
  cDgVsPp -> Write();
  cDgVsPp -> Close();

  TCanvas *cDpVsPp = new TCanvas("cDpVsPp", "Pi0-data vs. pi0-sim.", width, height);
  cDpVsPp -> cd();
  cDpVsPp -> SetLogy(1);
  cDpVsPp -> SetGrid(0, 0);
  hPythPP -> SetTitle("Pythia vs. pp-data, #pi^{0}-trigger; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPP -> Draw();
  hPythDP -> Draw("same");
  hDataP  -> Draw("same");
  lDpVsPp -> Draw();
  pDatPi0 -> Draw();
  pJets   -> Draw();
  cDpVsPp -> Write();
  cDpVsPp -> Close();

  TCanvas *cAuVsPg = new TCanvas("cAuVsPg", "Peripheral-AuAu vs. gamma-sim.", width, height);
  cAuVsPg -> cd();
  cAuVsPg -> SetLogy(1);
  cAuVsPg -> SetGrid(0, 0);
  hPythPG -> SetTitle("Pythia vs. AA-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPG -> Draw();
  hPythDG -> Draw("same");
  hAuAu   -> Draw("same");
  lAuVsPg -> Draw();
  pHadGam -> Draw();
  pJets   -> Draw();
  cAuVsPg -> Write();
  cAuVsPg -> Close();

  TCanvas *cAuVsPp = new TCanvas("cAuVsPp", "Peripheral-AuAu vs. pi0-sim.", width, height);
  cAuVsPp -> cd();
  cAuVsPp -> SetLogy(1);
  cAuVsPp -> SetGrid(0, 0);
  hPythPP -> SetTitle("Pythia vs. AA-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPP -> Draw();
  hPythDP -> Draw("same");
  hAuAu   -> Draw("same");
  lAuVsPp -> Draw();
  pHadPi0 -> Draw();
  pJets   -> Draw();
  cAuVsPp -> Write();
  cAuVsPp -> Close();

  // compare our pythia to alex's pythia
  TCanvas *cPaVsPg = new TCanvas("cPaVsPg", "gamma-sim. vs. pi-sim.", width, height);
  cPaVsPg -> cd();
  cPaVsPg -> SetLogy(1);
  cPaVsPg -> SetGrid(0, 0);
  hPythPG -> SetTitle("Pythia (8, #gamma-trigger) vs. Pythia (6, h^{#pm}-trigger); p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPG -> Draw();
  hPythDG -> Draw("same");
  hPythAA -> Draw("same p");
  lPaVsPg -> Draw();
  pPythG  -> Draw();
  pJets   -> Draw();
  cPaVsPg -> Write();
  cPaVsPg -> Close();

  TCanvas *cPaVsPp = new TCanvas("cPaVsPp", "pi0-sim. vs. pi-sim.", width, height);
  cPaVsPp -> cd();
  cPaVsPp -> SetLogy(1);
  cPaVsPp -> SetGrid(0, 0);
  hPythPP -> SetTitle("Pythia (8, #pi^{0}-trigger) vs. Pythia (6, h^{#pm}); p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPP -> Draw();
  hPythDP -> Draw("same");
  hPythAA -> Draw("same p");
  lPaVsPp -> Draw();
  pPythP  -> Draw();
  pJets   -> Draw();
  cPaVsPp -> Write();
  cPaVsPp -> Close();

  TCanvas *cAuVsG = new TCanvas("cAuVsG", "Peripheral-AuAu vs. gamma- (and hadron-) sim.", width, height);
  cAuVsG  -> cd();
  cAuVsG  -> SetLogy(1);
  cAuVsG  -> SetGrid(0, 0);
  hPythPG -> SetTitle("Pythia vs. AA-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPG -> Draw();
  hPythDG -> Draw("same");
  hPythAA -> Draw("same p");
  hAuAu   -> Draw("same");
  lAuVsG  -> Draw();
  pJets   -> Draw();
  cAuVsG  -> Write();
  cAuVsG  -> Close();

  TCanvas *cAuVsP = new TCanvas("cAuVsP", "Peripheral-AuAu vs. pi0- (and hadron-) sim.", width, height);
  cAuVsP  -> cd();
  cAuVsP  -> SetLogy(1);
  cAuVsP  -> SetGrid(0, 0);
  hPythPP -> SetTitle("Pythia vs. AA-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPP -> Draw();
  hPythDP -> Draw("same");
  hPythAA -> Draw("same p");
  hAuAu   -> Draw("same");
  lAuVsP  -> Draw();
  pJets   -> Draw();
  cAuVsP  -> Write();
  cAuVsP  -> Close();

  TCanvas *cPpVsG = new TCanvas("cPpVsG", "Gamma-data vs. gamma- (and hadron-) sim.", width, height);
  cPpVsG  -> cd();
  cPpVsG  -> SetLogy(1);
  cPpVsG  -> SetGrid(0, 0);
  hPythPG -> SetTitle("Pythia vs. pp-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPG -> Draw();
  hPythDG -> Draw("same");
  hPythAA -> Draw("same p");
  hDataG  -> Draw("same");
  lPpVsG  -> Draw();
  pJets   -> Draw();
  cPpVsG  -> Write();
  cPpVsG  -> Close();

  TCanvas *cPpVsP = new TCanvas("cPpVsP", "Pi0-data vs. pi0- (and hadron-) sim.", width, height);
  cPpVsP  -> cd();
  cPpVsP  -> SetLogy(1);
  cPpVsP  -> SetGrid(0, 0);
  hPythPP -> SetTitle("Pythia vs. pp-data; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hPythPP -> Draw();
  hPythDP -> Draw("same");
  hPythAA -> Draw("same p");
  hDataP  -> Draw("same");
  lPpVsP  -> Draw();
  pJets   -> Draw();
  cPpVsP  -> Write();
  cPpVsP  -> Close();


  // compare pp to AuAu
  hDataP -> SetMarkerStyle(4);
  hDataP -> SetMarkerColor(kMagenta);
  hDataP -> SetLineColor(kMagenta);

  TLegend *lDpVsAu = new TLegend(0.1, 0.1, 0.3, 0.3);
  lDpVsAu -> SetLineColor(kWhite);
  lDpVsAu -> SetFillColor(kWhite);
  lDpVsAu -> AddEntry(hDataP, "pp data, 200 GeV");
  lDpVsAu -> AddEntry(hAuAu, "AuAu data, 60-80%");

  TCanvas *cDpVsAu = new TCanvas("cDpVsAu", "pp-data vs. peripheral-AuAu", width, height);
  cDpVsAu -> cd();
  cDpVsAu -> SetLogy(1);
  cDpVsAu -> SetGrid(0, 0);
  hDataP  -> SetTitle("pp-data vs. peripheral-AuAu; p_{T}^{jet}; (1/N_{trg})dN_{jet}/dp_{T}^{jet}");
  hDataP  -> Draw();
  hAuAu   -> Draw("same");
  lDpVsAu -> Draw();
  pDatHad -> Draw();
  pJets   -> Draw();
  cDpVsAu -> Write();
  cDpVsAu -> Close();

  cout << "  Plots drawn..." << endl;


  // reset a few things
  hDataP  -> SetMarkerStyle(1);
  hDataP  -> SetMarkerColor(kBlack);
  hDataP  -> SetLineColor(kBlack);
  hDataP  -> SetTitle("Jet p_{T}, #pi^{0} (pp);;");
  hPythPG -> SetTitle("Jet p_{T}, #gamma^{dir} (particle-pythia);;");
  hPythPP -> SetTitle("Jet p_{T}, #pi^{0} (particle-pythia);;");


  // save file and close
  fOut    -> cd();
  hAuAu   -> Write();
  hDataG  -> Write();
  hDataP  -> Write();
  hPythPG -> Write();
  hPythPP -> Write();
  hPythDG -> Write();
  hPythDP -> Write();
  hPythAA -> Write();
  fOut    -> Close();

  fData -> cd();
  fData -> Close();
  fPythPG -> cd();
  fPythPG -> Close();
  fPythPP -> cd();
  fPythPP -> Close();
  fPythDG -> cd();
  fPythDG -> Close();
  fPythDP -> cd();
  fPythDP -> Close();


  cout << "Comparison script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
