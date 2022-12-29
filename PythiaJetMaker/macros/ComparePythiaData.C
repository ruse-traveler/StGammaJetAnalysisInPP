// 'ComparePythiaData.C'
// Derek Anderson
//
// Compare the QA plots from the pythia and data
// jet reconstruction macros.

#include <stdio>
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;


// filepaths
static const TString sPythiaP("../JetMacro/detTestPP.r03a02rm1chrg.d29m1y2017.root");
static const TString sPythiaD("../JetMacro/detTestDP.r03a02rm1chrg.d29m1y2017.root");
static const TString sData("../JetData/dataMystery.r03rm1chrg.d1m2y2017.root");
static const TString sOut("detTestPlotsP.r03a02rm1chrg.d2m2y2017.root");


void ComparePythiaData() {

  // lower verbosity
  gErrorIgnoreLevel = kFatal;
  gStyle -> SetOptStat(0);

  cout << "\nBeginning comparison script..." << endl;


  TFile *fOut     = new TFile(sOut.Data(), "recreate");
  TFile *fPythiaP = new TFile(sPythiaP.Data(), "read");
  TFile *fPythiaD = new TFile(sPythiaD.Data(), "read");
  TFile *fData    = new TFile(sData.Data(), "read");
  if (!fPythiaP || !fPythiaD || !fData) {
    cerr << "PANIC: couldn't open one of the input files!" << endl;
    return;
  }
  cout << "  Files opened..." << endl;


  // grab histograms [0 = par. pythia, 1 = det. pythia, 2 = data]
  TH1D *hNtrk[3];
  TH1D *hNtrkG1[3];
  TH1D *hNtrkG2[3];
  TH1D *hTrkPt[3];
  TH1D *hTrkPx[3];
  TH1D *hTrkPy[3];
  TH1D *hTrkPz[3];
  TH1D *hTrkDf[3];
  TH1D *hTrkDfG1[3];
  TH1D *hTrkDfG2[3];
  if (fPythiaP) {
    hNtrk[0]    = (TH1D*) fPythiaP -> Get("QA/hNtrk");
    hNtrkG1[0]  = (TH1D*) fPythiaP -> Get("QA/hNtrkG1");
    hNtrkG2[0]  = (TH1D*) fPythiaP -> Get("QA/hNtrkG2");
    hTrkPt[0]   = (TH1D*) fPythiaP -> Get("QA/hPtTrk");
    hTrkPx[0]   = (TH1D*) fPythiaP -> Get("QA/hPxTrk");
    hTrkPy[0]   = (TH1D*) fPythiaP -> Get("QA/hPyTrk");
    hTrkPz[0]   = (TH1D*) fPythiaP -> Get("QA/hPzTrk");
    hTrkDf[0]   = (TH1D*) fPythiaP -> Get("QA/hDfTrk");
    hTrkDfG1[0] = (TH1D*) fPythiaP -> Get("QA/hDfTrkG1");
    hTrkDfG2[0] = (TH1D*) fPythiaP -> Get("QA/hDfTrkG2");
  }
  if (fPythiaD) {
    hNtrk[1]    = (TH1D*) fPythiaD -> Get("QA/hNtrk");
    hNtrkG1[1]  = (TH1D*) fPythiaD -> Get("QA/hNtrkG1");
    hNtrkG2[1]  = (TH1D*) fPythiaD -> Get("QA/hNtrkG2");
    hTrkPt[1]   = (TH1D*) fPythiaD -> Get("QA/hPtTrk");
    hTrkPx[1]   = (TH1D*) fPythiaP -> Get("QA/hPxTrk");
    hTrkPy[1]   = (TH1D*) fPythiaD -> Get("QA/hPyTrk");
    hTrkPz[1]   = (TH1D*) fPythiaD -> Get("QA/hPzTrk");
    hTrkDf[1]   = (TH1D*) fPythiaD -> Get("QA/hDfTrk");
    hTrkDfG1[1] = (TH1D*) fPythiaD -> Get("QA/hDfTrkG1");
    hTrkDfG2[1] = (TH1D*) fPythiaD -> Get("QA/hDfTrkG2");
  }
  if (fData) {
    hNtrk[2]    = (TH1D*) fData -> Get("hNtrkP");
    hNtrkG1[2]  = (TH1D*) fData -> Get("hNtrkG1P");
    hNtrkG2[2]  = (TH1D*) fData -> Get("hNtrkG2P");
    hTrkPt[2]   = (TH1D*) fData -> Get("hTrkPtP");
    hTrkPx[2]   = (TH1D*) fData -> Get("hTrkPxP");
    hTrkPy[2]   = (TH1D*) fData -> Get("hTrkPyP");
    hTrkPz[2]   = (TH1D*) fData -> Get("hTrkPzP");
    hTrkDf[2]   = (TH1D*) fData -> Get("hTrkDfP");
    hTrkDfG1[2] = (TH1D*) fData -> Get("hTrkDfG1P");
    hTrkDfG2[2] = (TH1D*) fData -> Get("hTrkDfG2P");
  }
  cout << "  Histograms grabbed..." << endl;


  // set styles [blue for par. pythia, red for det. pythia, black for data]
  const Int_t    iColor[3]     = {4, 2, 1};
  const Double_t nRangeX[2]    = {0., 60.};
  const Double_t nRangeY[2]    = {0.00001, 1.};
  const Double_t nG1rangeX[2]  = {0., 30.};
  const Double_t nG1rangeY[2]  = {0.00001, 1.};
  const Double_t nG2rangeX[2]  = {0., 20.};
  const Double_t nG2rangeY[2]  = {0.00001, 1.};
  const Double_t pTrangeX[2]   = {0., 40.};
  const Double_t pTrangeY[2]   = {0.00001, 15.};
  const Double_t pRangeX[2]    = {-40., 40.};
  const Double_t pRangeY[2]    = {0.00001, 15.};
  const Double_t dFrangeY[2]   = {0., 6.};
  const Double_t dFG1rangeY[2] = {0., 4.};
  const Double_t dFG2rangeY[2] = {0., 2.};
  for (Int_t i = 0; i < 3; i++) {
    hNtrk[i]    -> SetLineColor(iColor[i]);
    hNtrk[i]    -> SetMarkerColor(iColor[i]);
    hNtrk[i]    -> SetTitle("Number of tracks");
    hNtrk[i]    -> GetXaxis() -> SetRangeUser(nRangeX[0], nRangeX[1]);
    hNtrk[i]    -> GetXaxis() -> SetLabelSize(0.02);
    hNtrk[i]    -> GetXaxis() -> SetTitle("N_{trk}");
    hNtrk[i]    -> GetXaxis() -> CenterTitle(true);
    hNtrk[i]    -> GetYaxis() -> SetRangeUser(nRangeY[0], nRangeY[1]);
    hNtrk[i]    -> GetYaxis() -> SetLabelSize(0.02);
    hNtrk[i]    -> GetYaxis() -> SetTitle("(1/N_{trg}) counts");
    hNtrk[i]    -> GetYaxis() -> CenterTitle(true);
    hNtrkG1[i]  -> SetLineColor(iColor[i]);
    hNtrkG1[i]  -> SetMarkerColor(iColor[i]);
    hNtrkG1[i]  -> SetTitle("Number of tracks with p_{T}^{trk} > 1 GeV/c");
    hNtrkG1[i]  -> GetXaxis() -> SetRangeUser(nG1rangeX[0], nG1rangeX[1]);
    hNtrkG1[i]  -> GetXaxis() -> SetLabelSize(0.02);
    hNtrkG1[i]  -> GetXaxis() -> SetTitle("N_{trk}(p_{T}^{trk} > 1 GeV/c)");
    hNtrkG1[i]  -> GetXaxis() -> CenterTitle(true);
    hNtrkG1[i]  -> GetYaxis() -> SetRangeUser(nG1rangeY[0], nG1rangeY[1]);
    hNtrkG1[i]  -> GetYaxis() -> SetLabelSize(0.02);
    hNtrkG1[i]  -> GetYaxis() -> SetTitle("(1/N_{trg}) counts");
    hNtrkG1[i]  -> GetYaxis() -> CenterTitle(true);
    hNtrkG2[i]  -> SetLineColor(iColor[i]);
    hNtrkG2[i]  -> SetMarkerColor(iColor[i]);
    hNtrkG2[i]  -> SetTitle("Number of tracks with p_{T}^{trk} > 2 GeV/c");
    hNtrkG2[i]  -> GetXaxis() -> SetRangeUser(nG2rangeX[0], nG2rangeX[1]);
    hNtrkG2[i]  -> GetXaxis() -> SetLabelSize(0.02);
    hNtrkG2[i]  -> GetXaxis() -> SetTitle("N_{trk}(p_{T}^{trk} > 2 GeV/c)");
    hNtrkG2[i]  -> GetXaxis() -> CenterTitle(true);
    hNtrkG2[i]  -> GetYaxis() -> SetRangeUser(nG2rangeY[0], nG2rangeY[1]);
    hNtrkG2[i]  -> GetYaxis() -> SetLabelSize(0.02);
    hNtrkG2[i]  -> GetYaxis() -> SetTitle("(1/N_{trg}) counts");
    hNtrkG2[i]  -> GetYaxis() -> CenterTitle(true);
    hTrkPt[i]   -> SetLineColor(iColor[i]);
    hTrkPt[i]   -> SetMarkerColor(iColor[i]);
    hTrkPt[i]   -> SetTitle("Track p_{T}");
    hTrkPt[i]   -> GetXaxis() -> SetRangeUser(pTrangeX[0], pTrangeX[1]);
    hTrkPt[i]   -> GetXaxis() -> SetLabelSize(0.02);
    hTrkPt[i]   -> GetXaxis() -> SetTitle("p_{T}^{trk}");
    hTrkPt[i]   -> GetXaxis() -> CenterTitle(true);
    hTrkPt[i]   -> GetYaxis() -> SetRangeUser(pTrangeY[0], pTrangeY[1]);
    hTrkPt[i]   -> GetYaxis() -> SetLabelSize(0.02);
    hTrkPt[i]   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/dp_{T}^{trk}");
    hTrkPt[i]   -> GetYaxis() -> CenterTitle(true);
    hTrkPx[i]   -> SetLineColor(iColor[i]);
    hTrkPx[i]   -> SetMarkerColor(iColor[i]);
    hTrkPx[i]   -> SetTitle("Track p_{x}");
    hTrkPx[i]   -> GetXaxis() -> SetRangeUser(pRangeX[0], pRangeX[1]);
    hTrkPx[i]   -> GetXaxis() -> SetLabelSize(0.02);
    hTrkPx[i]   -> GetXaxis() -> SetTitle("p_{x}^{trk}");
    hTrkPx[i]   -> GetXaxis() -> CenterTitle(true);
    hTrkPx[i]   -> GetYaxis() -> SetRangeUser(pRangeY[0], pRangeY[1]);
    hTrkPx[i]   -> GetYaxis() -> SetLabelSize(0.02);
    hTrkPx[i]   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/dp_{x}^{trk}");
    hTrkPx[i]   -> GetYaxis() -> CenterTitle(true);
    hTrkPy[i]   -> SetLineColor(iColor[i]);
    hTrkPy[i]   -> SetMarkerColor(iColor[i]);
    hTrkPy[i]   -> SetTitle("Track p_{y}");
    hTrkPy[i]   -> GetXaxis() -> SetRangeUser(pRangeX[0], pRangeX[1]);
    hTrkPy[i]   -> GetXaxis() -> SetLabelSize(0.02);
    hTrkPy[i]   -> GetXaxis() -> SetTitle("p_{y}^{trk}");
    hTrkPy[i]   -> GetXaxis() -> CenterTitle(true);
    hTrkPy[i]   -> GetYaxis() -> SetRangeUser(pRangeY[0], pRangeY[1]);
    hTrkPy[i]   -> GetYaxis() -> SetLabelSize(0.02);
    hTrkPy[i]   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/dp_{y}^{trk}");
    hTrkPy[i]   -> GetYaxis() -> CenterTitle(true);
    hTrkPz[i]   -> SetLineColor(iColor[i]);
    hTrkPz[i]   -> SetMarkerColor(iColor[i]);
    hTrkPz[i]   -> SetTitle("Track p_{z}");
    hTrkPz[i]   -> GetXaxis() -> SetRangeUser(pRangeX[0], pRangeX[1]);
    hTrkPz[i]   -> GetXaxis() -> SetLabelSize(0.02);
    hTrkPz[i]   -> GetXaxis() -> SetTitle("p_{z}^{trk}");
    hTrkPz[i]   -> GetXaxis() -> CenterTitle(true);
    hTrkPz[i]   -> GetYaxis() -> SetRangeUser(pRangeY[0], pRangeY[1]);
    hTrkPz[i]   -> GetYaxis() -> SetLabelSize(0.02);
    hTrkPz[i]   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/dp_{z}^{trk}");
    hTrkPz[i]   -> GetYaxis() -> CenterTitle(true);
    hTrkDf[i]   -> SetLineColor(iColor[i]);
    hTrkDf[i]   -> SetMarkerColor(iColor[i]);
    hTrkDf[i]   -> SetTitle("Track #Delta#varphi");
    hTrkDf[i]   -> GetXaxis() -> SetTitle("#Delta#varphi_{trk}");
    hTrkDf[i]   -> GetXaxis() -> CenterTitle(true);
    hTrkDf[i]   -> GetXaxis() -> SetLabelSize(0.02);
    hTrkDf[i]   -> GetYaxis() -> SetRangeUser(dFrangeY[0], dFrangeY[1]);
    hTrkDf[i]   -> GetYaxis() -> SetLabelSize(0.02);
    hTrkDf[i]   -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/d#Delta#varphi");
    hTrkDf[i]   -> GetYaxis() -> CenterTitle(true);
    hTrkDfG1[i] -> SetLineColor(iColor[i]);
    hTrkDfG1[i] -> SetMarkerColor(iColor[i]);
    hTrkDfG1[i] -> SetTitle("#Delta#varphi of tracks with p_{T}^{trk} > 1 GeV/c");
    hTrkDfG1[i] -> GetXaxis() -> SetTitle("#Delta#varphi_{trk}(p_{T}^{trk} > 1 GeV/c)");
    hTrkDfG1[i] -> GetXaxis() -> SetLabelSize(0.02);
    hTrkDfG1[i] -> GetXaxis() -> CenterTitle(true);
    hTrkDfG1[i] -> GetYaxis() -> SetRangeUser(dFG1rangeY[0], dFG1rangeY[1]);
    hTrkDfG1[i] -> GetYaxis() -> SetLabelSize(0.02);
    hTrkDfG1[i] -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/d#Delta#varphi");
    hTrkDfG1[i] -> GetYaxis() -> CenterTitle(true);
    hTrkDfG2[i] -> SetLineColor(iColor[i]);
    hTrkDfG2[i] -> SetMarkerColor(iColor[i]);
    hTrkDfG2[i] -> SetTitle("#Delta#varphi of tracks with p_{T}^{trk} > 2 GeV/c");
    hTrkDfG2[i] -> GetXaxis() -> SetTitle("#Delta#varphi_{trk}(p_{T}^{trk} > 2 GeV/c)");
    hTrkDfG2[i] -> GetXaxis() -> SetLabelSize(0.02);
    hTrkDfG2[i] -> GetXaxis() -> CenterTitle(true);
    hTrkDfG2[i] -> GetYaxis() -> SetRangeUser(dFG2rangeY[0], dFG2rangeY[1]);
    hTrkDfG2[i] -> GetYaxis() -> SetLabelSize(0.02);
    hTrkDfG2[i] -> GetYaxis() -> SetTitle("(1/N_{trg}) dN_{trk}/d#Delta#varphi");
    hTrkDfG2[i] -> GetYaxis() -> CenterTitle(true);
  }
  cout << "  Styles set..." << endl;


  // calculate and record integrals
  TString  sInt[10][3];
  Double_t ints[10][3];
  for (Int_t i = 0; i < 3; i++) {
    ints[0][i] = hNtrk[i]    -> Integral();
    ints[1][i] = hNtrkG1[i]  -> Integral();
    ints[2][i] = hNtrkG2[i]  -> Integral();
    ints[3][i] = hTrkPt[i]   -> Integral();
    ints[4][i] = hTrkPx[i]   -> Integral();
    ints[5][i] = hTrkPy[i]   -> Integral();
    ints[6][i] = hTrkPz[i]   -> Integral();
    ints[7][i] = hTrkDf[i]   -> Integral();
    ints[8][i] = hTrkDfG1[i] -> Integral();
    ints[9][i] = hTrkDfG2[i] -> Integral();
  }

  Ssiz_t size = 0;
  Ssiz_t nDig = 0;
  Ssiz_t nDec = 2;
  for (Int_t i = 0; i < 10; i++) {
    for (Int_t j = 0; j < 3; j++) {
      TString string("");
      string += ints[i][j];
      nDig    = string.First(".");
      size    = (nDig + nDec) + 1;
      sInt[i][j].Append(string, size);
    }
  }
  cout << "  Integrals computed..." << endl;


  // create labels
  const TString sStartP("#color[4]{I_{pyth}^{par.} = ");
  const TString sStartD("#color[2]{I_{pyth}^{det.} = ");
  const TString sStartA("#color[1]{I_{data} = ");
  const TString sEnd("}");

  TLegend *lg = new TLegend(0.1, 0.1, 0.3, 0.3);
  lg -> SetLineColor(kWhite);
  lg -> SetFillColor(kWhite);
  lg -> AddEntry(hNtrk[0], "Pythia (particle)");
  lg -> AddEntry(hNtrk[1], "Pythia (Detector)");
  lg -> AddEntry(hNtrk[2], "Data (run 9)");

  TPaveText *pt = new TPaveText(0.3, 0.1, 0.5, 0.3, "NDC NB");
  pt -> SetLineColor(kWhite);
  pt -> SetFillColor(kWhite);
  pt -> AddText("pp collisions, #pi^{0} trigger");
  pt -> AddText("p_{T}^{trg} #in (9, 30) GeV/c, |#eta^{trg}| < 0.9");
  pt -> AddText("p_{T}^{trk} #in (0.2, 30) GeV/c, |#eta^{trk}| < 1");

  TString    sIntN[3];
  TPaveText *pIntN = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntN -> SetLineColor(kWhite);
  pIntN -> SetFillColor(kWhite);
  sIntN[0]  = sStartP;
  sIntN[1]  = sStartD; 
  sIntN[2]  = sStartA;
  sIntN[0] += sInt[0][0];
  sIntN[1] += sInt[0][1];
  sIntN[2] += sInt[0][2];
  sIntN[0] += sEnd;
  sIntN[1] += sEnd;
  sIntN[2] += sEnd;
  pIntN -> AddText(sIntN[0].Data());
  pIntN -> AddText(sIntN[1].Data());
  pIntN -> AddText(sIntN[2].Data());

  TString    sIntNG1[3];
  TPaveText *pIntNG1 = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntNG1 -> SetLineColor(kWhite);
  pIntNG1 -> SetFillColor(kWhite);
  sIntNG1[0]  = sStartP; 
  sIntNG1[1]  = sStartD;
  sIntNG1[2]  = sStartA;
  sIntNG1[0] += sInt[1][0];
  sIntNG1[1] += sInt[1][1];
  sIntNG1[2] += sInt[1][2];
  sIntNG1[0] += sEnd;
  sIntNG1[1] += sEnd;
  sIntNG1[2] += sEnd;
  pIntNG1 -> AddText(sIntNG1[0].Data());
  pIntNG1 -> AddText(sIntNG1[1].Data());
  pIntNG1 -> AddText(sIntNG1[2].Data());

  TString    sIntNG2[3];
  TPaveText *pIntNG2 = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntNG2 -> SetLineColor(kWhite);
  pIntNG2 -> SetFillColor(kWhite);
  sIntNG2[0]  = sStartP; 
  sIntNG2[1]  = sStartD;
  sIntNG2[2]  = sStartA;
  sIntNG2[0] += sInt[2][0];
  sIntNG2[1] += sInt[2][1];
  sIntNG2[2] += sInt[2][2];
  sIntNG2[0] += sEnd;
  sIntNG2[1] += sEnd;
  sIntNG2[2] += sEnd;
  pIntNG2 -> AddText(sIntNG2[0].Data());
  pIntNG2 -> AddText(sIntNG2[1].Data());
  pIntNG2 -> AddText(sIntNG2[2].Data());

  TString    sIntPt[3];
  TPaveText *pIntPt = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntPt -> SetLineColor(kWhite);
  pIntPt -> SetFillColor(kWhite);
  sIntPt[0]  = sStartP; 
  sIntPt[1]  = sStartD;
  sIntPt[2]  = sStartA;
  sIntPt[0] += sInt[3][0];
  sIntPt[1] += sInt[3][1];
  sIntPt[2] += sInt[3][2];
  sIntPt[0] += sEnd;
  sIntPt[1] += sEnd;
  sIntPt[2] += sEnd;
  pIntPt -> AddText(sIntPt[0].Data());
  pIntPt -> AddText(sIntPt[1].Data());
  pIntPt -> AddText(sIntPt[2].Data());

  TString    sIntPx[3];
  TPaveText *pIntPx = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntPx -> SetLineColor(kWhite);
  pIntPx -> SetFillColor(kWhite);
  sIntPx[0]  = sStartP; 
  sIntPx[1]  = sStartD;
  sIntPx[2]  = sStartA;
  sIntPx[0] += sInt[4][0];
  sIntPx[1] += sInt[4][1];
  sIntPx[2] += sInt[4][2];
  sIntPx[0] += sEnd;
  sIntPx[1] += sEnd;
  sIntPx[2] += sEnd;
  pIntPx -> AddText(sIntPx[0].Data());
  pIntPx -> AddText(sIntPx[1].Data());
  pIntPx -> AddText(sIntPx[2].Data());

  TString    sIntPy[3];
  TPaveText *pIntPy = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntPy -> SetLineColor(kWhite);
  pIntPy -> SetFillColor(kWhite);
  sIntPy[0]  = sStartP; 
  sIntPy[1]  = sStartD;
  sIntPy[2]  = sStartA;
  sIntPy[0] += sInt[5][0];
  sIntPy[1] += sInt[5][1];
  sIntPy[2] += sInt[5][2];
  sIntPy[0] += sEnd;
  sIntPy[1] += sEnd;
  sIntPy[2] += sEnd;
  pIntPy -> AddText(sIntPy[0].Data());
  pIntPy -> AddText(sIntPy[1].Data());
  pIntPy -> AddText(sIntPy[2].Data());

  TString    sIntPz[3];
  TPaveText *pIntPz = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntPz -> SetLineColor(kWhite);
  pIntPz -> SetFillColor(kWhite);
  sIntPz[0]  = sStartP; 
  sIntPz[1]  = sStartD;
  sIntPz[2]  = sStartA;
  sIntPz[0] += sInt[6][0];
  sIntPz[1] += sInt[6][1];
  sIntPz[2] += sInt[6][2];
  sIntPz[0] += sEnd;
  sIntPz[1] += sEnd;
  sIntPz[2] += sEnd;
  pIntPz -> AddText(sIntPz[0].Data());
  pIntPz -> AddText(sIntPz[1].Data());
  pIntPz -> AddText(sIntPz[2].Data());

  TString    sIntDf[3];
  TPaveText *pIntDf = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntDf -> SetLineColor(kWhite);
  pIntDf -> SetFillColor(kWhite);
  sIntDf[0]  = sStartP; 
  sIntDf[1]  = sStartD;
  sIntDf[2]  = sStartA;
  sIntDf[0] += sInt[7][0];
  sIntDf[1] += sInt[7][1];
  sIntDf[2] += sInt[7][2];
  sIntDf[0] += sEnd;
  sIntDf[1] += sEnd;
  sIntDf[2] += sEnd;
  pIntDf -> AddText(sIntDf[0].Data());
  pIntDf -> AddText(sIntDf[1].Data());
  pIntDf -> AddText(sIntDf[2].Data());

  TString    sIntDfG1[3];
  TPaveText *pIntDfG1 = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntDfG1 -> SetLineColor(kWhite);
  pIntDfG1 -> SetFillColor(kWhite);
  sIntDfG1[0]  = sStartP; 
  sIntDfG1[1]  = sStartD;
  sIntDfG1[2]  = sStartA;
  sIntDfG1[0] += sInt[8][0];
  sIntDfG1[1] += sInt[8][1];
  sIntDfG1[2] += sInt[8][2];
  sIntDfG1[0] += sEnd;
  sIntDfG1[1] += sEnd;
  sIntDfG1[2] += sEnd;
  pIntDfG1 -> AddText(sIntDfG1[0].Data());
  pIntDfG1 -> AddText(sIntDfG1[1].Data());
  pIntDfG1 -> AddText(sIntDfG1[2].Data());

  TString    sIntDfG2[3];
  TPaveText *pIntDfG2 = new TPaveText(0.5, 0.1, 0.7, 0.3, "NDC NB");
  pIntDfG2 -> SetLineColor(kWhite);
  pIntDfG2 -> SetFillColor(kWhite);
  sIntDfG2[0]  = sStartP; 
  sIntDfG2[1]  = sStartD;
  sIntDfG2[2]  = sStartA;
  sIntDfG2[0] += sInt[9][0];
  sIntDfG2[1] += sInt[9][1];
  sIntDfG2[2] += sInt[9][2];
  sIntDfG2[0] += sEnd;
  sIntDfG2[1] += sEnd;
  sIntDfG2[2] += sEnd;
  pIntDfG2 -> AddText(sIntDfG2[0].Data());
  pIntDfG2 -> AddText(sIntDfG2[1].Data());
  pIntDfG2 -> AddText(sIntDfG2[2].Data());

  cout << "  Labels created..." << endl;


  // change into output file
  fOut -> cd();

  // draw plots
  const Int_t width  = 800;
  const Int_t height = 800;
  TCanvas *cNtrk = new TCanvas("cNtrk", "No. of tracks", width, height);
  cNtrk    -> cd();
  cNtrk    -> SetGrid(0, 0);
  cNtrk    -> SetTickx(1);
  cNtrk    -> SetTicky(1);
  cNtrk    -> SetLogy(1);
  hNtrk[0] -> Draw();
  hNtrk[1] -> Draw("same");
  hNtrk[2] -> Draw("same");
  lg       -> Draw();
  pt       -> Draw();
  pIntN    -> Draw();
  cNtrk    -> Write();
  cNtrk    -> Close();

  TCanvas *cNtrkG1 = new TCanvas("cNtrkG1", "No. of tracks w/ pT>1", width, height);
  cNtrkG1    -> cd();
  cNtrkG1    -> SetGrid(0, 0);
  cNtrkG1    -> SetTickx(1);
  cNtrkG1    -> SetTicky(1);
  cNtrkG1    -> SetLogy(1);
  hNtrkG1[0] -> Draw();
  hNtrkG1[1] -> Draw("same");
  hNtrkG1[2] -> Draw("same");
  lg         -> Draw();
  pt         -> Draw();
  pIntNG1    -> Draw();
  cNtrkG1    -> Write();
  cNtrkG1    -> Close();

  TCanvas *cNtrkG2 = new TCanvas("cNtrkG2", "No. of tracks w/ pT>2", width, height);
  cNtrkG2    -> cd();
  cNtrkG2    -> SetGrid(0, 0);
  cNtrkG2    -> SetTickx(1);
  cNtrkG2    -> SetTicky(1);
  cNtrkG2    -> SetLogy(1);
  hNtrkG2[0] -> Draw();
  hNtrkG2[1] -> Draw("same");
  hNtrkG2[2] -> Draw("same");
  lg         -> Draw();
  pt         -> Draw();
  pIntNG2    -> Draw();
  cNtrkG2    -> Write();
  cNtrkG2    -> Close();

  TCanvas *cTrkPt = new TCanvas("cTrkPt", "Track pT", width, height);
  cTrkPt    -> cd();
  cTrkPt    -> SetGrid(0, 0);
  cTrkPt    -> SetTickx(1);
  cTrkPt    -> SetTicky(1);
  cTrkPt    -> SetLogy(1);
  hTrkPt[0] -> Draw();
  hTrkPt[1] -> Draw("same");
  hTrkPt[2] -> Draw("same");
  lg        -> Draw();
  pt        -> Draw();
  pIntPt    -> Draw();
  cTrkPt    -> Write();
  cTrkPt    -> Close();

  TCanvas *cTrkPx = new TCanvas("cTrkPx", "Track pX", width, height);
  cTrkPx    -> cd();
  cTrkPx    -> SetGrid(0, 0);
  cTrkPx    -> SetTickx(1);
  cTrkPx    -> SetTicky(1);
  cTrkPx    -> SetLogy(1);
  hTrkPx[0] -> Draw();
  hTrkPx[1] -> Draw("same");
  hTrkPx[2] -> Draw("same");
  lg        -> Draw();
  pt        -> Draw();
  pIntPx    -> Draw();
  cTrkPx    -> Write();
  cTrkPx    -> Close();

  TCanvas *cTrkPy = new TCanvas("cTrkPy", "Track pY", width, height);
  cTrkPy    -> cd();
  cTrkPy    -> SetGrid(0, 0);
  cTrkPy    -> SetTickx(1);
  cTrkPy    -> SetTicky(1);
  cTrkPy    -> SetLogy(1);
  hTrkPy[0] -> Draw();
  hTrkPy[1] -> Draw("same");
  hTrkPy[2] -> Draw("same");
  lg        -> Draw();
  pt        -> Draw();
  pIntPy    -> Draw();
  cTrkPy    -> Write();
  cTrkPy    -> Close();

  TCanvas *cTrkPz = new TCanvas("cTrkPz", "Track pZ", width, height);
  cTrkPz    -> cd();
  cTrkPz    -> SetGrid(0, 0);
  cTrkPz    -> SetTickx(1);
  cTrkPz    -> SetTicky(1);
  cTrkPz    -> SetLogy(1);
  hTrkPz[0] -> Draw();
  hTrkPz[1] -> Draw("same");
  hTrkPz[2] -> Draw("same");
  lg        -> Draw();
  pt        -> Draw();
  pIntPz    -> Draw();
  cTrkPz    -> Write();
  cTrkPz    -> Close();

  TCanvas *cTrkDf = new TCanvas("cTrkDf", "Track DeltaPhi", width, height);
  cTrkDf    -> cd();
  cTrkDf    -> SetGrid(0, 0);
  cTrkDf    -> SetTickx(1);
  cTrkDf    -> SetTicky(1);
  hTrkDf[0] -> Draw();
  hTrkDf[1] -> Draw("same");
  hTrkDf[2] -> Draw("same");
  lg        -> Draw();
  pt        -> Draw();
  pIntDf    -> Draw();
  cTrkDf    -> Write();
  cTrkDf    -> Close();

  TCanvas *cTrkDfG1 = new TCanvas("cTrkDfG1", "DeltaPhi of tracks w/ pT>1", width, height);
  cTrkDfG1    -> cd();
  cTrkDfG1    -> SetGrid(0, 0);
  cTrkDfG1    -> SetTickx(1);
  cTrkDfG1    -> SetTicky(1);
  hTrkDfG1[0] -> Draw();
  hTrkDfG1[1] -> Draw("same");
  hTrkDfG1[2] -> Draw("same");
  lg          -> Draw();
  pt          -> Draw();
  pIntDfG1    -> Draw();
  cTrkDfG1    -> Write();
  cTrkDfG1    -> Close();

  TCanvas *cTrkDfG2 = new TCanvas("cTrkDfG2", "DeltaPhi of tracks w/ pT>2", width, height);
  cTrkDfG2    -> cd();
  cTrkDfG2    -> SetGrid(0, 0);
  cTrkDfG2    -> SetTickx(1);
  cTrkDfG2    -> SetTicky(1);
  hTrkDfG2[0] -> Draw();
  hTrkDfG2[1] -> Draw("same");
  hTrkDfG2[2] -> Draw("same");
  lg          -> Draw();
  pt          -> Draw();
  pIntDfG2    -> Draw();
  cTrkDfG2    -> Write();
  cTrkDfG2    -> Close();

  cout << "  Plots drawn..." << endl;


  // set names
  hNtrk[0]    -> SetName("hNtrkP");
  hNtrkG1[0]  -> SetName("hNtrkG1P");
  hNtrkG2[0]  -> SetName("hNtrkG2P");
  hTrkPt[0]   -> SetName("hTrkPtP");
  hTrkPx[0]   -> SetName("hTrkPxP");
  hTrkPy[0]   -> SetName("hTrkPyP");
  hTrkPz[0]   -> SetName("hTrkPzP");
  hTrkDf[0]   -> SetName("hTrkDfP");
  hTrkDfG1[0] -> SetName("hTrkDfG1P");
  hTrkDfG2[0] -> SetName("hTrkDfG2P");
  hNtrk[1]    -> SetName("hNtrkD");
  hNtrkG1[1]  -> SetName("hNtrkG1D");
  hNtrkG2[1]  -> SetName("hNtrkG2D");
  hTrkPt[1]   -> SetName("hTrkPtD");
  hTrkPx[1]   -> SetName("hTrkPxD");
  hTrkPy[1]   -> SetName("hTrkPyD");
  hTrkPz[1]   -> SetName("hTrkPzD");
  hTrkDf[1]   -> SetName("hTrkDfD");
  hTrkDfG1[1] -> SetName("hTrkDfG1D");
  hTrkDfG2[1] -> SetName("hTrkDfG2D");
  hNtrk[2]    -> SetName("hNtrkA");
  hNtrkG1[2]  -> SetName("hNtrkG1A");
  hNtrkG2[2]  -> SetName("hNtrkG2A");
  hTrkPt[2]   -> SetName("hTrkPtA");
  hTrkPx[2]   -> SetName("hTrkPxA");
  hTrkPy[2]   -> SetName("hTrkPyA");
  hTrkPz[2]   -> SetName("hTrkPzA");
  hTrkDf[2]   -> SetName("hTrkDfA");
  hTrkDfG1[2] -> SetName("hTrkDfG1A");
  hTrkDfG2[2] -> SetName("hTrkDfG2A");

  // close files
  fOut -> cd();
  for (Int_t i = 0; i < 3; i++) {
    hNtrk[i]    -> Write();
    hNtrkG1[i]  -> Write();
    hNtrkG2[i]  -> Write();
    hTrkPt[i]   -> Write(); 
    hTrkPx[i]   -> Write(); 
    hTrkPy[i]   -> Write();
    hTrkPz[i]   -> Write();
    hTrkDf[i]   -> Write();
    hTrkDfG1[i] -> Write();
    hTrkDfG2[i] -> Write();
  }
  fOut     -> Close();
  fPythiaP -> cd();
  fPythiaP -> Close();
  fPythiaD -> cd();
  fPythiaD -> Close();
  fData    -> cd();
  fData    -> Close();

  cout << "Script finished!\n" << endl;

}

// End ------------------------------------------------------------------------
