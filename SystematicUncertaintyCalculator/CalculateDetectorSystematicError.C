// 'CalculateDetectorSystematicError.C'
// Derek Anderson
// 07.27.2021
//
// Use this to calculate a systematic
// uncertainty (associated with
// unfolding) wrt. a default or average
// value.
//
// NOTE: use the 'TrigId' and 'TrigBin'
//       variables change the color palette
//       and labels to fit the trigger used.
//       'TrigId' sets the trigger species:
//         TrigId  = 0 -- pi0
//         TrigId  = 1 -- gamma
//       And 'TrgBin' sets the eTtrg bin:
//         TrigBin = 0 -- (9, 11) GeV
//         TrigBin = 1 -- (11, 15) GeV
//         TrigBin = 2 -- (15, 20) GeV

#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TError.h"
#include "TColor.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// histogram constants
static const UInt_t NPad(2);
static const UInt_t NPlot(2);

// systematic constants
static const UInt_t NCat(2);
static const UInt_t NHist(NCat + 1);
static const UInt_t NSysEffP(2);
static const UInt_t NSysResP(2);
static const UInt_t NSysEffG(2);
static const UInt_t NSysResG(2);

// trigger constants
static const UInt_t NPi0(3);
static const UInt_t NGam(3);
static const UInt_t NTrgs(NPi0 + NGam);

// caclulation parameters
static const Bool_t StatsAreCommon(true);
static const Bool_t RemoveExtraBins(true);

// trigger parameters
static const UInt_t NSysEff(NSysEffG);
static const UInt_t NSysRes(NSysResG);
static const UInt_t TrigId(1);
static const UInt_t TrigBin(2);



void CalculateDetectorSystematicError() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Calculating unfolding systematic uncertainty..." << endl;

  // io parameters
  const TString sOut("detectorSys.oneExtraBin_pTbinHuge.et1520r02gam.d20m10y2021.root");
  const TString sHistD("hUnfolded");
  const TString sHistS("hUnfolded");

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");
  const TString sTest("only detector sys. uncertainties");

  // pi0 io and hist parameters
  const TString sInPi0D("et1520r02pi0_rebinUnfold/pp200r9pi0.baseline_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root");
  const TString sInPi0E[NSysEffP]  = {"et1520r02pi0_rebinUnfold/pp200r9pi0.forEffSysM04_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root",
                                      "et1520r02pi0_rebinUnfold/pp200r9pi0.forEffSysP04_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root"};
  const TString sInPi0R[NSysResP]  = {"et1520r02pi0_rebinUnfold/pp200r9pi0.forResSysV1_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root",
                                      "et1520r02pi0_rebinUnfold/pp200r9pi0.forResSysV2_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root"};
  const TString sEffPi0[NSysEffP]  = {"#epsilon_{trk} - 4\%", "#epsilon_{trk} + 4\%"};
  const TString sResPi0[NSysResP]  = {"#Deltap_{T}^{trk}(1)", "#Deltap_{T}^{trk}(2)"};
  const TString sNameVEP[NSysEffP] = {"hVarEffM04", "hVarEffP04"}; 
  const TString sNameVRP[NSysResP] = {"hVarResV1", "hVarResV2"};
  const TString sNameREP[NSysEffP] = {"hRatEffM04", "hRatEffP04"};
  const TString sNameRRP[NSysResP] = {"hRatResV1", "hRatResV2"};

  // gamma io parameters
  const TString sInGamD("et1520r02gam_rebinUnfold/pp200r9gam.baseline_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root");
  const TString sInGamE[NSysEffG]  = {"et1520r02gam_rebinUnfold/pp200r9gam.forEffSysM04_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root",
                                      "et1520r02gam_rebinUnfold/pp200r9gam.forEffSysP04_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root"};
  const TString sInGamR[NSysResG]  = {"et1520r02gam_rebinUnfold/pp200r9gam.forResSysV1_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root",
                                      "et1520r02gam_rebinUnfold/pp200r9gam.forResSysV2_pTbinHuge.et1520r02qt05130.p0m1k4n38t5.root"};
  const TString sEffGam[NSysEffG]  = {"#epsilon_{trk} - 4\%", "#epsilon_{trk} + 4\%"};
  const TString sResGam[NSysResG]  = {"#Deltap_{T}^{trk}(1)", "#Deltap_{T}^{trk}(2)"};
  const TString sNameVEG[NSysEffG] = {"hVarEffM04", "hVarEffP04"};
  const TString sNameVRG[NSysResG] = {"hVarResV1", "hVarResV2"};
  const TString sNameREG[NSysEffG] = {"hRatEffM04", "hRatEffP04"};
  const TString sNameRRG[NSysResG] = {"hRatResV1", "hRatResV2"};

  // pi0 plot parameters
  const TString sTotPi0("total systematic");
  const TString sDefPi0("default [n_{iter}^{*} = 4]");
  const Float_t xCalcPi0[NPi0][NPlot]      = {{0.2, 30.}, {0.2, 30.}, {0.2, 30.}};
  const UInt_t  fMarEffPi0[NSysEffP]       = {21, 23};
  const UInt_t  fMarResPi0[NSysResP]       = {22, 20};
  const UInt_t  fColAvgPi0[NPi0]           = {819, 799, 899};
  const UInt_t  fColEffPi0[NPi0][NSysEffP] = {{859, 879},
                                              {839, 859},
                                              {819, 839}};
  const UInt_t  fColResPi0[NPi0][NSysResP] = {{899, 799},
                                              {879, 899},
                                              {859, 879}};
  const UInt_t  fColSysPi0[NPi0][NHist]    = {{593, 590, 603},
                                              {425, 422, 435},
                                              {409, 406, 419}};

  // gamma plot parameters
  const TString sTotGam("total systematic");
  const TString sDefGam("default [n_{iter}^{*} = 4, central B value]");
  const Float_t xCalcGam[NGam][NPlot]      = {{0.2, 11.}, {0.2, 20.}, {0.2, 25.}};
  const UInt_t  fMarEffGam[NSysEffG]       = {21, 23};
  const UInt_t  fMarResGam[NSysResG]       = {22, 20};
  const UInt_t  fColAvgGam[NGam]           = {859, 869, 819};
  const UInt_t  fColEffGam[NGam][NSysEffG] = {{899, 799},
                                              {879, 809},
                                              {859, 879}};
  const UInt_t  fColResGam[NGam][NSysResG] = {{819, 839},
                                              {829, 849},
                                              {899, 799}};
  const UInt_t  fColSysGam[NGam][NHist]    = {{625, 622, 635},
                                              {609, 606, 619},
                                              {593, 590, 603}};

  // general systematic parameters
  const TString sNameD("hDefault");
  const TString sNameA("hEffAvg");
  const TString sNameAR("hRatEffAvg");
  const TString sNameT("hTotal");
  const TString sNamePT("hPerTotal");
  const TString sLabelA("average of #epsilon_{trk} variations");
  const TString sNameS[NCat]   = {"hSysEff", "hSysRes"};
  const TString sNameP[NCat]   = {"hPerEff", "hPerRes"};
  const TString sNameU[NCat]   = {"hSumEff", "hSumRes"};
  const TString sNameM[NCat]   = {"hSumPerEff", "hSumPerRes"};
  const TString sLabelS[NHist] = {"#epsilon_{trk} systematic", "#Deltap_{T}^{trk} systematic", "total systematic"};
  const TString sLabelU[NCat]  = {"#sigma_{sys}^{det}(#epsilon_{trk})", "#sigma_{sys}^{det}(#epsilon_{trk}) #oplus #sigma_{sys}^{det}(#Deltap_{T}^{trk})"};
  const UInt_t  fFilSys[NCat]  = {3345, 3354};
  const UInt_t  fMarSys[NCat]  = {1, 1};

  // general plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T}^{unfold} [GeV/c]");
  const TString  sTitleY("(1/N^{trg}) d^{3}N^{jet}/d(p_{T}^{unfold} #eta^{jet}) [GeV/c]^{-1}");
  const TString  sTitleR("var. / default");
  const TString  sTitleP("fractional error");
  const TString  sPlot("Plot");
  const Double_t plotRange[NPlot] = {-1., 30.};

  // parse trigger selection
  UInt_t fTrg(0);
  if (TrigId == 0) {
    fTrg = TrigBin;
  } else {
    fTrg = TrigBin + 3;
  }

  // select trigger
  TString sTrg;
  TString sLabelD;
  TString sLabelDS;
  TString sLabelT;
  TString sInD;
  TString sInE[NSysEff];
  TString sInR[NSysRes];
  TString sNameVE[NSysEff];
  TString sNameVR[NSysRes];
  TString sNameRE[NSysEff];
  TString sNameRR[NSysRes];
  TString sLabelE[NSysEff];
  TString sLabelR[NSysRes];
  Float_t calcRange[NPlot];
  UInt_t  fColDef;
  UInt_t  fColAvg;
  UInt_t  fColSys[NCat];
  UInt_t  fColEff[NSysEff];
  UInt_t  fColRes[NSysRes];
  UInt_t  fMarEff[NSysEff];
  UInt_t  fMarRes[NSysRes];
  switch (fTrg) {
    case 0:
      sTrg     = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      sInD     = sInPi0D;
      sLabelD  = sDefPi0;
      sLabelDS = sDefPi0;
      sLabelT  = sTotPi0;
      fColAvg  = fColAvgPi0[0];
      fColDef  = fColSysPi0[0][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffPi0[0][iEff];
        fMarEff[iEff] = fMarEffPi0[iEff];
        sInE[iEff]    = sInPi0E[iEff];
        sNameVE[iEff] = sNameVEP[iEff];
        sNameRE[iEff] = sNameREP[iEff];
        sLabelE[iEff] = sEffPi0[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResPi0[0][iRes];
        fMarRes[iRes] = fMarResPi0[iRes];
        sInR[iRes]    = sInPi0R[iRes];
        sNameVR[iRes] = sNameVRP[iRes];
        sNameRR[iRes] = sNameRRP[iRes];
        sLabelR[iRes] = sResPi0[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysPi0[0][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcPi0[0][iPlot];
      }
      break;
    case 1:
      sTrg     = "#pi^{0} trig., E_{T}^{meas} #in (11, 15) GeV";
      sInD     = sInPi0D;
      sLabelD  = sDefPi0;
      sLabelDS = sDefPi0;
      sLabelT  = sTotPi0;
      fColAvg  = fColAvgPi0[1];
      fColDef  = fColSysPi0[1][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffPi0[1][iEff];
        fMarEff[iEff] = fMarEffPi0[iEff];
        sInE[iEff]    = sInPi0E[iEff];
        sNameVE[iEff] = sNameVEP[iEff];
        sNameRE[iEff] = sNameREP[iEff];
        sLabelE[iEff] = sEffPi0[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResPi0[1][iRes];
        fMarRes[iRes] = fMarResPi0[iRes];
        sInR[iRes]    = sInPi0R[iRes];
        sNameVR[iRes] = sNameVRP[iRes];
        sNameRR[iRes] = sNameRRP[iRes];
        sLabelR[iRes] = sResPi0[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysPi0[1][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcPi0[1][iPlot];
      }
      break;
    case 2:
      sTrg     = "#pi^{0} trig., E_{T}^{meas} #in (15, 20) GeV";
      sInD     = sInPi0D;
      sLabelD  = sDefPi0;
      sLabelDS = sDefPi0;
      sLabelT  = sTotPi0;
      fColAvg  = fColAvgPi0[2];
      fColDef  = fColSysPi0[2][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffPi0[2][iEff];
        fMarEff[iEff] = fMarEffPi0[iEff];
        sInE[iEff]    = sInPi0E[iEff];
        sNameVE[iEff] = sNameVEP[iEff];
        sNameRE[iEff] = sNameREP[iEff];
        sLabelE[iEff] = sEffPi0[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResPi0[2][iRes];
        fMarRes[iRes] = fMarResPi0[iRes];
        sInR[iRes]    = sInPi0R[iRes];
        sNameVR[iRes] = sNameVRP[iRes];
        sNameRR[iRes] = sNameRRP[iRes];
        sLabelR[iRes] = sResPi0[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysPi0[2][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcPi0[2][iPlot];
      }
      break;
    case 3:
      sTrg     = "#gamma_{dir} trig., E_{T}^{meas} #in (9, 11) GeV";
      sInD     = sInGamD;
      sLabelD  = sDefGam;
      sLabelDS = sDefGam;
      sLabelT  = sTotGam;
      fColAvg  = fColAvgGam[0];
      fColDef  = fColSysGam[0][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffGam[0][iEff];
        fMarEff[iEff] = fMarEffGam[iEff];
        sInE[iEff]    = sInGamE[iEff];
        sNameVE[iEff] = sNameVEG[iEff];
        sNameRE[iEff] = sNameREG[iEff];
        sLabelE[iEff] = sEffGam[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResGam[0][iRes];
        fMarRes[iRes] = fMarResGam[iRes];
        sInR[iRes]    = sInGamR[iRes];
        sNameVR[iRes] = sNameVRG[iRes];
        sNameRR[iRes] = sNameRRG[iRes];
        sLabelR[iRes] = sResGam[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysGam[0][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcGam[0][iPlot];
      }
      break;
    case 4:
      sTrg     = "#gamma_{dir} trig., E_{T}^{meas} #in (11, 15) GeV";
      sInD     = sInGamD;
      sLabelD  = sDefGam;
      sLabelDS = sDefGam;
      sLabelT  = sTotGam;
      fColAvg  = fColAvgGam[1];
      fColDef  = fColSysGam[1][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffGam[1][iEff];
        fMarEff[iEff] = fMarEffGam[iEff];
        sInE[iEff]    = sInGamE[iEff];
        sNameVE[iEff] = sNameVEG[iEff];
        sNameRE[iEff] = sNameREG[iEff];
        sLabelE[iEff] = sEffGam[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResGam[1][iRes];
        fMarRes[iRes] = fMarResGam[iRes];
        sInR[iRes]    = sInGamR[iRes];
        sNameVR[iRes] = sNameVRG[iRes];
        sNameRR[iRes] = sNameRRG[iRes];
        sLabelR[iRes] = sResGam[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysGam[1][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcGam[1][iPlot];
      }
      break;
    case 5:
      sTrg     = "#gamma_{dir} trig., E_{T}^{meas} #in (15, 20) GeV";
      sInD     = sInGamD;
      sLabelD  = sDefGam;
      sLabelDS = sDefGam;
      sLabelT  = sTotGam;
      fColAvg  = fColAvgGam[2];
      fColDef  = fColSysGam[2][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffGam[2][iEff];
        fMarEff[iEff] = fMarEffGam[iEff];
        sInE[iEff]    = sInGamR[iEff];
        sNameVE[iEff] = sNameVEG[iEff];
        sNameRE[iEff] = sNameREG[iEff];
        sLabelE[iEff] = sEffGam[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResGam[2][iRes];
        fMarRes[iRes] = fMarResGam[iRes];
        sInR[iRes]    = sInGamR[iRes];
        sNameVR[iRes] = sNameVRG[iRes];
        sNameRR[iRes] = sNameRRG[iRes];
        sLabelR[iRes] = sResGam[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysGam[2][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcGam[2][iPlot];
      }
      break;
    default:
      sTrg     = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      sInD     = sInPi0D;
      sLabelD  = sDefPi0;
      sLabelDS = sDefPi0;
      sLabelT  = sTotPi0;
      fColAvg  = fColAvgPi0[0];
      fColDef  = fColSysPi0[0][NCat];
      for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
        fColEff[iEff] = fColEffPi0[0][iEff];
        fMarEff[iEff] = fMarEffPi0[iEff];
        sInE[iEff]    = sInPi0R[iEff];
        sNameVE[iEff] = sNameVEP[iEff];
        sNameRE[iEff] = sNameREP[iEff];
        sLabelE[iEff] = sEffPi0[iEff];
      }
      for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
        fColRes[iRes] = fColResPi0[0][iRes];
        fMarRes[iRes] = fMarResPi0[iRes];
        sInR[iRes]    = sInPi0R[iRes];
        sNameVR[iRes] = sNameVRP[iRes];
        sNameRR[iRes] = sNameRRP[iRes];
        sLabelR[iRes] = sResPi0[iRes];
      }
      for (UInt_t iCat = 0; iCat < NCat; iCat++) {
        fColSys[iCat] = fColSysPi0[0][iCat];
      }
      for (UInt_t iPlot = 0; iPlot < NPlot; iPlot++) {
        calcRange[iPlot] = xCalcPi0[0][iPlot];
      }
      break;
  }  // end switch-case
  cout << "    Selected trigger." << endl;

  // open files
  TFile *fOut = new TFile(sOut.Data(), "recreate");
  TFile *fInD = new TFile(sInD.Data(), "read");
  if (!fOut || !fInD) {
    cerr << "PANIC: couldn't open a file!\n"
         << "       fOut = " << fOut << ", fInD = " << fInD
         << endl;
    return;
  }

  TFile *fInE[NSysEff];
  TFile *fInR[NSysRes];
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    fInE[iEff] = new TFile(sInE[iEff].Data(), "read");
    if (!fInE[iEff]) {
      cerr << "PANIC: couldn't open a file!\n"
           << "       fInE[" << iEff <<"] = " << fInE[iEff]
           << endl;
      return;
    }
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    fInR[iRes] = new TFile(sInR[iRes].Data(), "read");
    if (!fInR[iRes]) {
      cerr << "PANIC: couldn't open a file!\n"
           << "       fInR[" << iRes <<"] = " << fInR[iRes]
           << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hDefault = (TH1D*) fInD -> Get(sHistD.Data());
  if (!hDefault) {
    cerr << "PANIC: couldn't grab a histogram!\n"
         << "       hDefault = " << hDefault
         << endl;
    return;
  }
  hDefault -> SetName(sNameD.Data());
  hDefault -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);

  TH1D *hVarE[NSysEff];
  TH1D *hVarR[NSysRes];
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    hVarE[iEff] = (TH1D*) fInE[iEff] -> Get(sHistS.Data());
    if (!hVarE[iEff]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hVarE[" << iEff << "] = " << hVarE[iEff]
           << endl;
      return;
    }
    hVarE[iEff] -> SetName(sNameVE[iEff].Data());
    hVarE[iEff] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    hVarR[iRes] = (TH1D*) fInR[iRes] -> Get(sHistS.Data());
    if (!hVarR[iRes]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hVarR[" << iRes << "] = " << hVarR[iRes]
           << endl;
      return;
    }
    hVarR[iRes] -> SetName(sNameVR[iRes].Data());
    hVarR[iRes] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Grabbed histograms." << endl;

  // calculate ratios
  TH1D *hDivE[NSysEff];
  TH1D *hDivR[NSysRes];
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    hDivE[iEff] = (TH1D*) hDefault -> Clone();
    hDivE[iEff] -> SetName(sNameRE[iEff].Data());
    hDivE[iEff] -> Reset("ICE");

    const UInt_t nBinE = hDivE[iEff] -> GetNbinsX();
    for (UInt_t iBinE = 1; iBinE < (nBinE + 1); iBinE++) {
      const Double_t def  = hDefault    -> GetBinContent(iBinE);
      const Double_t var  = hVarE[iEff] -> GetBinContent(iBinE);
      const Double_t rawD = hDefault    -> GetBinError(iBinE);
      const Double_t rawV = hVarE[iEff] -> GetBinError(iBinE);
      const Double_t relD = (rawD / def);
      const Double_t relV = (rawV / var);
      const Double_t divD = (var / def);

      // calculate error differences
      Double_t div = 0.;
      Double_t raw = 0.;
      Double_t rel = 0.;
      if (StatsAreCommon) {
        div = divD;
        raw = rawV;
        rel = raw / def;
      } else {
        div = divD;
        raw = TMath::Sqrt(TMath::Abs((rawV * rawV) - (rawD * rawD)));
        rel = raw / def;
      }  // end if (StatsAreCommon)

      // check for nan's
      const Bool_t isGoodDef = (TMath::Abs(def) > 0.);
      const Bool_t isGoodVar = (TMath::Abs(var) > 0.);

      // assign bin values
      const Double_t valD = div;
      const Double_t errD = valD * rel;
      if (isGoodDef && isGoodVar) {
        hDivE[iEff] -> SetBinContent(iBinE, valD);
        hDivE[iEff] -> SetBinError(iBinE, errD);
      } else {
        hDivE[iEff] -> SetBinContent(iBinE, 0.);
        hDivE[iEff] -> SetBinError(iBinE, 0.);
      }
    }  // end bin loop
    hDivE[iEff] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    hDivR[iRes] = (TH1D*) hDefault -> Clone();
    hDivR[iRes] -> SetName(sNameRR[iRes].Data());
    hDivR[iRes] -> Reset("ICE");

    const UInt_t nBinR = hDivR[iRes] -> GetNbinsX();
    for (UInt_t iBinR = 1; iBinR < (nBinR + 1); iBinR++) {
      const Double_t def  = hDefault    -> GetBinContent(iBinR);
      const Double_t var  = hVarR[iRes] -> GetBinContent(iBinR);
      const Double_t rawD = hDefault    -> GetBinError(iBinR);
      const Double_t rawV = hVarR[iRes] -> GetBinError(iBinR);
      const Double_t relD = (rawD / def);
      const Double_t relV = (rawV / var);
      const Double_t divD = (var / def);

      // calculate error differences
      Double_t div = 0.;
      Double_t raw = 0.;
      Double_t rel = 0.;
      if (StatsAreCommon) {
        div = divD;
        raw = rawV;
        rel = raw / def;
      } else {
        div = divD;
        raw = TMath::Sqrt(TMath::Abs((rawV * rawV) - (rawD * rawD)));
        rel = raw / def;
      }  // end if (StatsAreCommon)

      // check for nan's
      const Bool_t isGoodDef = (TMath::Abs(def) > 0.);
      const Bool_t isGoodVar = (TMath::Abs(var) > 0.);

      // assign bin values
      const Double_t valD = div;
      const Double_t errD = valD * rel;
      if (isGoodDef && isGoodVar) {
        hDivR[iRes] -> SetBinContent(iBinR, valD);
        hDivR[iRes] -> SetBinError(iBinR, errD);
      } else {
        hDivR[iRes] -> SetBinContent(iBinR, 0.);
        hDivR[iRes] -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
    hDivR[iRes] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Calculated ratios." << endl;

  // average efficiencies
  TH1D *hEffAvg = (TH1D*) hDefault -> Clone();
  TH1D *hDivA   = (TH1D*) hDefault -> Clone();
  hEffAvg -> SetName(sNameA.Data());
  hDivA   -> SetName(sNameAR.Data());
  hEffAvg -> Reset("ICES");
  hDivA   -> Reset("ICES");
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    hEffAvg -> Add(hVarE[iEff]);
  }
  hEffAvg -> Scale(1. / (Double_t) NSysEff);
  hDivA   -> Divide(hEffAvg, hDefault, 1., 1.); 
  cout << "    Calculated average of efficiency variations." << endl;

  // calculate systematics
  TH1D *hSys[NHist];
  TH1D *hPer[NHist];
  TH1D *hSum[NCat];
  TH1D *hSumP[NCat];
  for (UInt_t iCat = 0; iCat < NCat; iCat++) {
    hSys[iCat]  = (TH1D*) hDefault -> Clone();
    hPer[iCat]  = (TH1D*) hDefault -> Clone();
    hSum[iCat]  = (TH1D*) hDefault -> Clone();
    hSumP[iCat] = (TH1D*) hDefault -> Clone();
    hSys[iCat]  -> SetName(sNameS[iCat].Data());
    hPer[iCat]  -> SetName(sNameP[iCat].Data());
    hSum[iCat]  -> SetName(sNameU[iCat].Data());
    hSumP[iCat] -> SetName(sNameM[iCat].Data());
  }
  hSys[NCat] = (TH1D*) hDefault -> Clone();
  hPer[NCat] = (TH1D*) hDefault -> Clone();
  hSys[NCat] -> SetName(sNameT.Data());
  hPer[NCat] -> SetName(sNamePT.Data());

  // clear total systematic histograms
  const UInt_t nBinTot = hSys[NCat] -> GetNbinsX();
  for (UInt_t iBinTot = 0; iBinTot < nBinTot; iBinTot++) {
    hPer[NCat] -> SetBinContent(iBinTot, 0.);
    hPer[NCat] -> SetBinError(iBinTot, 0.);
    hSys[NCat] -> SetBinError(iBinTot, 0.);
  }

  // calculate efficiency systematic
  for (UInt_t iBin = 1; iBin < (nBinTot + 1); iBin++) {

    // loop over efficiency variations
    Double_t effSysErr(0.);
    for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {

      // calculate bin differences
      const Double_t avgValE = hEffAvg     -> GetBinContent(iBin);
      const Double_t avgErrE = hEffAvg     -> GetBinError(iBin);
      const Double_t effVal  = hVarE[iEff] -> GetBinContent(iBin);
      const Double_t effErr  = hVarE[iEff] -> GetBinError(iBin);
      const Double_t valDifE = TMath::Abs(avgValE - effVal);
      const Double_t errDifE = TMath::Sqrt(TMath::Abs((avgErrE * avgErrE) - (effErr * effErr)));
      const Bool_t   isNonzero = (effVal > 0.);
      if (!isNonzero) {
        continue;
      }
      if (!StatsAreCommon) {
        hVarE[iEff] -> SetBinError(iBin, errDifE);
      }

      // calculate possible systematic
      Double_t sysValE(0.);
      if (StatsAreCommon) {
        sysValE = valDifE;
      } else {
        if (valDifR > errDifR) {
          sysValE = TMath::Sqrt((valDifE * valDifE) - (errDifE * errDifE));
        } else {
          sysValE = 0.;
        }
      }

      // keep only largest difference
      if (sysValE > effSysErr) {
        effSysErr = sysValE;
      }
    }  // end variation loop

    // get running totals
    const Double_t defVal = hDefault   -> GetBinContent(iBin);
    const Double_t totVal = hSys[NCat] -> GetBinContent(iBin);
    const Double_t totPer = hPer[NCat] -> GetBinError(iBin);
    const Double_t totErr = hSys[NCat] -> GetBinError(iBin);

    // update totals and assign bin values
    Double_t newTotErr = TMath::Sqrt((totErr * totErr) + (effSysErr * effSysErr));
    Double_t newTotPer = totPer;
    Double_t effSysPer = 0.;
    if ((effSysErr > 0.) && (defVal > 0.)) {
      newTotPer = newTotErr / totVal;
      effSysPer = effSysErr / defVal;
    }
    hPer[0]    -> SetBinContent(iBin, 1.);
    hPer[NCat] -> SetBinContent(iBin, 1.);
    hSumP[0]   -> SetBinContent(iBin, 1.);
    hPer[0]    -> SetBinError(iBin, effSysPer);
    hPer[NCat] -> SetBinError(iBin, newTotPer);
    hSys[0]    -> SetBinError(iBin, effSysErr);
    hSys[NCat] -> SetBinError(iBin, newTotErr);
    hSum[0]    -> SetBinError(iBin, newTotErr);
    hSumP[0]   -> SetBinError(iBin, newTotPer);

  } // end bin loop
  cout << "    Calculated efficiency systematic." << endl;

  // calculate resolution systematic
  for (UInt_t iBin = 1; iBin < (nBinTot + 1); iBin++) {

    // loop over resolution variations
    Double_t resSysErr(0.);
    for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {

      // calculate bin differences
      const Double_t defValR   = hDefault    -> GetBinContent(iBin);
      const Double_t defErrR   = hDefault    -> GetBinError(iBin);
      const Double_t resVal    = hVarR[iRes] -> GetBinContent(iBin);
      const Double_t resErr    = hVarR[iRes] -> GetBinError(iBin);
      const Double_t valDifR   = TMath::Abs(defValR - resVal);
      const Double_t errDifR   = TMath::Sqrt(TMath::Abs((defErrR * defErrR) - (resErr * resErr)));
      const Bool_t   isNonzero = (resVal > 0.);
      if (!isNonzero) {
        continue;
      }
      if (!StatsAreCommon) {
        hVarR[iRes] -> SetBinError(iBin, errDifR);
      }

      // calculate possible systematic
      Double_t sysValR(0.);
      if (StatsAreCommon) {
        sysValR = valDifR;
      } else {
        if (valDifR > errDifR) {
          sysValR = TMath::Sqrt((valDifR * valDifR) - (errDifR * errDifR));
        } else {
          sysValR = 0.;
        }
      }

      // keep only largest difference
      if (sysValR > resSysErr) {
        resSysErr = sysValR;
      }
    }  // end variation loop

    // get running totals
    const Double_t defVal = hDefault   -> GetBinContent(iBin);
    const Double_t totVal = hSys[NCat] -> GetBinContent(iBin);
    const Double_t totPer = hPer[NCat] -> GetBinError(iBin);
    const Double_t totErr = hSys[NCat] -> GetBinError(iBin);

    // update totals and assign bin values
    Double_t newTotErr = TMath::Sqrt((totErr * totErr) + (resSysErr * resSysErr));
    Double_t newTotPer = totPer;
    Double_t resSysPer = 0.;
    if ((resSysErr > 0.) && (defVal > 0.)) {
      newTotPer = newTotErr / totVal;
      resSysPer = resSysErr / defVal;
    }
    hPer[1]    -> SetBinContent(iBin, 1.);
    hPer[NCat] -> SetBinContent(iBin, 1.);
    hSumP[1]   -> SetBinContent(iBin, 1.);
    hPer[1]    -> SetBinError(iBin, resSysPer);
    hPer[NCat] -> SetBinError(iBin, newTotPer);
    hSys[1]    -> SetBinError(iBin, resSysErr);
    hSys[NCat] -> SetBinError(iBin, newTotErr);
    hSum[1]    -> SetBinError(iBin, newTotErr);
    hSumP[1]   -> SetBinError(iBin, newTotPer);

  } // end bin loop
  cout << "    Calculated resolution systematic." << endl;

  // zero out extra bins
  if (RemoveExtraBins) {
    const UInt_t nBins = hDefault -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
      const Double_t binCnt    = hDefault -> GetBinCenter(iBin);
      const Bool_t   isInRange = ((binCnt > calcRange[0]) && (binCnt < calcRange[1]));
      if (!isInRange) {
        hDefault -> SetBinContent(iBin, 0.);
        hDefault -> SetBinError(iBin, 0.);
        hEffAvg  -> SetBinContent(iBin, 0.);
        hDivA    -> SetBinContent(iBin, 0.);
        for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
          hVarE[iEff]  -> SetBinContent(iBin, 0.);
          hDivE[iEff]  -> SetBinContent(iBin, 0.);
          hVarE[iEff]  -> SetBinError(iBin, 0.);
          hDivE[iEff]  -> SetBinError(iBin, 0.);
        }
        for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
          hVarR[iRes]  -> SetBinContent(iBin, 0.);
          hDivR[iRes]  -> SetBinContent(iBin, 0.);
          hVarR[iRes]  -> SetBinError(iBin, 0.);
          hDivR[iRes]  -> SetBinError(iBin, 0.);
        }
        for (UInt_t iCat = 0; iCat < NCat; iCat++) {
          hSys[iCat]  -> SetBinContent(iBin, 0.);
          hPer[iCat]  -> SetBinContent(iBin, 0.);
          hSum[iCat]  -> SetBinContent(iBin, 0.);
          hSumP[iCat] -> SetBinContent(iBin, 0.);
          hSys[iCat]  -> SetBinError(iBin, 0.);
          hPer[iCat]  -> SetBinError(iBin, 0.);
          hSum[iCat]  -> SetBinError(iBin, 0.);
          hSumP[iCat] -> SetBinError(iBin, 0.);
        }  // end systematic loop
        hSys[NCat] -> SetBinContent(iBin, 0.);
        hPer[NCat] -> SetBinContent(iBin, 0.);
        hSys[NCat] -> SetBinError(iBin, 0.);
        hPer[NCat] -> SetBinError(iBin, 0.);
      }
    }  // end bin loop
    cout << "    Zeroed out irrelevant bins." << endl;
  }

  // set styles
  const UInt_t  fColT(923);
  const UInt_t  fMarD(20);
  const UInt_t  fMarA(34);
  const UInt_t  fMarU(1);
  const UInt_t  fMarT(1);
  const UInt_t  fFilD(0);
  const UInt_t  fFilV(0);
  const UInt_t  fFilU(1001);
  const UInt_t  fFilT(0);
  const UInt_t  fLinD(1);
  const UInt_t  fLinS(1);
  const UInt_t  fLinT(1);
  const UInt_t  fWidD(1);
  const UInt_t  fWidS(1);
  const UInt_t  fWidT(2);
  const UInt_t  fTxt(42);
  const UInt_t  fAln(12);
  const UInt_t  fCnt(1);
  const Float_t fAlphaD(0.);
  const Float_t fAlphaS(0.5);
  const Float_t fLab[NPad]  = {0.074, 0.04};
  const Float_t fTit[NPad]  = {0.074, 0.04};
  const Float_t fOffX[NPad] = {1.1, 1.};
  const Float_t fOffY[NPad] = {0.7, 1.3};
  hDefault -> SetMarkerColor(fColDef);
  hDefault -> SetMarkerStyle(fMarD);
  hDefault -> SetFillColor(fColDef);
  hDefault -> SetFillStyle(fFilD);
  hDefault -> SetLineColor(fColDef);
  hDefault -> SetLineStyle(fLinD);
  hDefault -> SetLineWidth(fWidD);
  hDefault -> SetTitle(sTitle.Data());
  hDefault -> SetTitleFont(fTxt);
  hDefault -> GetXaxis() -> SetTitle(sTitleX.Data());
  hDefault -> GetXaxis() -> SetTitleFont(fTxt);
  hDefault -> GetXaxis() -> SetTitleSize(fTit[1]);
  hDefault -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hDefault -> GetXaxis() -> SetLabelFont(fTxt);
  hDefault -> GetXaxis() -> SetLabelSize(fLab[1]);
  hDefault -> GetXaxis() -> CenterTitle(fCnt);
  hDefault -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDefault -> GetYaxis() -> SetTitleFont(fTxt);
  hDefault -> GetYaxis() -> SetTitleSize(fTit[1]);
  hDefault -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hDefault -> GetYaxis() -> SetLabelFont(fTxt);
  hDefault -> GetYaxis() -> SetLabelSize(fLab[1]);
  hDefault -> GetYaxis() -> CenterTitle(fCnt);
  hEffAvg  -> SetMarkerColor(fColAvg);
  hEffAvg  -> SetMarkerStyle(fMarA);
  hEffAvg  -> SetFillColor(fColAvg);
  hEffAvg  -> SetFillStyle(fFilD);
  hEffAvg  -> SetLineColor(fColAvg);
  hEffAvg  -> SetLineStyle(fLinD);
  hEffAvg  -> SetLineWidth(fWidD);
  hEffAvg  -> SetTitle(sTitle.Data());
  hEffAvg  -> SetTitleFont(fTxt);
  hEffAvg  -> GetXaxis() -> SetTitle(sTitleX.Data());
  hEffAvg  -> GetXaxis() -> SetTitleFont(fTxt);
  hEffAvg  -> GetXaxis() -> SetTitleSize(fTit[1]);
  hEffAvg  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hEffAvg  -> GetXaxis() -> SetLabelFont(fTxt);
  hEffAvg  -> GetXaxis() -> SetLabelSize(fLab[1]);
  hEffAvg  -> GetXaxis() -> CenterTitle(fCnt);
  hEffAvg  -> GetYaxis() -> SetTitle(sTitleY.Data());
  hEffAvg  -> GetYaxis() -> SetTitleFont(fTxt);
  hEffAvg  -> GetYaxis() -> SetTitleSize(fTit[1]);
  hEffAvg  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hEffAvg  -> GetYaxis() -> SetLabelFont(fTxt);
  hEffAvg  -> GetYaxis() -> SetLabelSize(fLab[1]);
  hEffAvg  -> GetYaxis() -> CenterTitle(fCnt);
  hDivA    -> SetMarkerColor(fColAvg);
  hDivA    -> SetMarkerStyle(fMarA);
  hDivA    -> SetFillColor(fColAvg);
  hDivA    -> SetFillStyle(fFilD);
  hDivA    -> SetLineColor(fColAvg);
  hDivA    -> SetLineStyle(fLinD);
  hDivA    -> SetLineWidth(fWidD);
  hDivA    -> SetTitle(sTitle.Data());
  hDivA    -> SetTitleFont(fTxt);
  hDivA    -> GetXaxis() -> SetTitle(sTitleX.Data());
  hDivA    -> GetXaxis() -> SetTitleFont(fTxt);
  hDivA    -> GetXaxis() -> SetTitleSize(fTit[0]);
  hDivA    -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hDivA    -> GetXaxis() -> SetLabelFont(fTxt);
  hDivA    -> GetXaxis() -> SetLabelSize(fLab[0]);
  hDivA    -> GetXaxis() -> CenterTitle(fCnt);
  hDivA    -> GetYaxis() -> SetTitle(sTitleY.Data());
  hDivA    -> GetYaxis() -> SetTitleFont(fTxt);
  hDivA    -> GetYaxis() -> SetTitleSize(fTit[0]);
  hDivA    -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hDivA    -> GetYaxis() -> SetLabelFont(fTxt);
  hDivA    -> GetYaxis() -> SetLabelSize(fLab[0]);
  hDivA    -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    hVarE[iEff] -> SetMarkerColor(fColEff[iEff]);
    hVarE[iEff] -> SetMarkerStyle(fMarEff[iEff]);
    hVarE[iEff] -> SetFillColor(fColEff[iEff]);
    hVarE[iEff] -> SetFillStyle(fFilV);
    hVarE[iEff] -> SetLineColor(fColEff[iEff]);
    hVarE[iEff] -> SetLineStyle(fLinS);
    hVarE[iEff] -> SetLineWidth(fWidS);
    hVarE[iEff] -> SetTitle(sTitle.Data());
    hVarE[iEff] -> SetTitleFont(fTxt);
    hVarE[iEff] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hVarE[iEff] -> GetXaxis() -> SetTitleFont(fTxt);
    hVarE[iEff] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hVarE[iEff] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hVarE[iEff] -> GetXaxis() -> SetLabelFont(fTxt);
    hVarE[iEff] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hVarE[iEff] -> GetXaxis() -> CenterTitle(fCnt);
    hVarE[iEff] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hVarE[iEff] -> GetYaxis() -> SetTitleFont(fTxt);
    hVarE[iEff] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hVarE[iEff] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hVarE[iEff] -> GetYaxis() -> SetLabelFont(fTxt);
    hVarE[iEff] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hVarE[iEff] -> GetYaxis() -> CenterTitle(fCnt);
    hDivE[iEff] -> SetMarkerColor(fColEff[iEff]);
    hDivE[iEff] -> SetMarkerStyle(fMarEff[iEff]);
    hDivE[iEff] -> SetFillColor(fColEff[iEff]);
    hDivE[iEff] -> SetFillStyle(fFilV);
    hDivE[iEff] -> SetLineColor(fColEff[iEff]);
    hDivE[iEff] -> SetLineStyle(fLinS);
    hDivE[iEff] -> SetLineWidth(fWidS);
    hDivE[iEff] -> SetTitle(sTitle.Data());
    hDivE[iEff] -> SetTitleFont(fTxt);
    hDivE[iEff] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDivE[iEff] -> GetXaxis() -> SetTitleFont(fTxt);
    hDivE[iEff] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hDivE[iEff] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hDivE[iEff] -> GetXaxis() -> SetLabelFont(fTxt);
    hDivE[iEff] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hDivE[iEff] -> GetXaxis() -> CenterTitle(fCnt);
    hDivE[iEff] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hDivE[iEff] -> GetYaxis() -> SetTitleFont(fTxt);
    hDivE[iEff] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hDivE[iEff] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hDivE[iEff] -> GetYaxis() -> SetLabelFont(fTxt);
    hDivE[iEff] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hDivE[iEff] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    hVarR[iRes] -> SetMarkerColor(fColRes[iRes]);
    hVarR[iRes] -> SetMarkerStyle(fMarRes[iRes]);
    hVarR[iRes] -> SetFillColor(fColRes[iRes]);
    hVarR[iRes] -> SetFillStyle(fFilV);
    hVarR[iRes] -> SetLineColor(fColRes[iRes]);
    hVarR[iRes] -> SetLineStyle(fLinS);
    hVarR[iRes] -> SetLineWidth(fWidS);
    hVarR[iRes] -> SetTitle(sTitle.Data());
    hVarR[iRes] -> SetTitleFont(fTxt);
    hVarR[iRes] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hVarR[iRes] -> GetXaxis() -> SetTitleFont(fTxt);
    hVarR[iRes] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hVarR[iRes] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hVarR[iRes] -> GetXaxis() -> SetLabelFont(fTxt);
    hVarR[iRes] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hVarR[iRes] -> GetXaxis() -> CenterTitle(fCnt);
    hVarR[iRes] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hVarR[iRes] -> GetYaxis() -> SetTitleFont(fTxt);
    hVarR[iRes] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hVarR[iRes] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hVarR[iRes] -> GetYaxis() -> SetLabelFont(fTxt);
    hVarR[iRes] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hVarR[iRes] -> GetYaxis() -> CenterTitle(fCnt);
    hDivR[iRes] -> SetMarkerColor(fColRes[iRes]);
    hDivR[iRes] -> SetMarkerStyle(fMarRes[iRes]);
    hDivR[iRes] -> SetFillColor(fColRes[iRes]);
    hDivR[iRes] -> SetFillStyle(fFilV);
    hDivR[iRes] -> SetLineColor(fColRes[iRes]);
    hDivR[iRes] -> SetLineStyle(fLinS);
    hDivR[iRes] -> SetLineWidth(fWidS);
    hDivR[iRes] -> SetTitle(sTitle.Data());
    hDivR[iRes] -> SetTitleFont(fTxt);
    hDivR[iRes] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDivR[iRes] -> GetXaxis() -> SetTitleFont(fTxt);
    hDivR[iRes] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hDivR[iRes] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hDivR[iRes] -> GetXaxis() -> SetLabelFont(fTxt);
    hDivR[iRes] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hDivR[iRes] -> GetXaxis() -> CenterTitle(fCnt);
    hDivR[iRes] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hDivR[iRes] -> GetYaxis() -> SetTitleFont(fTxt);
    hDivR[iRes] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hDivR[iRes] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hDivR[iRes] -> GetYaxis() -> SetLabelFont(fTxt);
    hDivR[iRes] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hDivR[iRes] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iCat = 0; iCat < NCat; iCat++) {
    hSys[iCat]  -> SetMarkerColor(fColSys[iCat]);
    hSys[iCat]  -> SetMarkerStyle(fMarSys[iCat]);
    hSys[iCat]  -> SetFillColor(fColSys[iCat]);
    hSys[iCat]  -> SetFillStyle(fFilSys[iCat]);
    hSys[iCat]  -> SetLineColor(fColSys[iCat]);
    hSys[iCat]  -> SetLineStyle(fLinS);
    hSys[iCat]  -> SetLineWidth(fWidS);
    hSys[iCat]  -> SetTitle(sTitle.Data());
    hSys[iCat]  -> SetTitleFont(fTxt);
    hSys[iCat]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hSys[iCat]  -> GetXaxis() -> SetTitleFont(fTxt);
    hSys[iCat]  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hSys[iCat]  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hSys[iCat]  -> GetXaxis() -> SetLabelFont(fTxt);
    hSys[iCat]  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hSys[iCat]  -> GetXaxis() -> CenterTitle(fCnt);
    hSys[iCat]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    hSys[iCat]  -> GetYaxis() -> SetTitleFont(fTxt);
    hSys[iCat]  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hSys[iCat]  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hSys[iCat]  -> GetYaxis() -> SetLabelFont(fTxt);
    hSys[iCat]  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hSys[iCat]  -> GetYaxis() -> CenterTitle(fCnt);
    hPer[iCat]  -> SetMarkerColor(fColSys[iCat]);
    hPer[iCat]  -> SetMarkerStyle(fMarSys[iCat]);
    hPer[iCat]  -> SetFillColor(fColSys[iCat]);
    hPer[iCat]  -> SetFillStyle(fFilSys[iCat]);
    hPer[iCat]  -> SetLineColor(fColSys[iCat]);
    hPer[iCat]  -> SetLineStyle(fLinS);
    hPer[iCat]  -> SetLineWidth(fWidS);
    hPer[iCat]  -> SetTitle(sTitle.Data());
    hPer[iCat]  -> SetTitleFont(fTxt);
    hPer[iCat]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPer[iCat]  -> GetXaxis() -> SetTitleFont(fTxt);
    hPer[iCat]  -> GetXaxis() -> SetTitleSize(fTit[0]);
    hPer[iCat]  -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hPer[iCat]  -> GetXaxis() -> SetLabelFont(fTxt);
    hPer[iCat]  -> GetXaxis() -> SetLabelSize(fLab[0]);
    hPer[iCat]  -> GetXaxis() -> CenterTitle(fCnt);
    hPer[iCat]  -> GetYaxis() -> SetTitle(sTitleP.Data());
    hPer[iCat]  -> GetYaxis() -> SetTitleFont(fTxt);
    hPer[iCat]  -> GetYaxis() -> SetTitleSize(fTit[0]);
    hPer[iCat]  -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hPer[iCat]  -> GetYaxis() -> SetLabelFont(fTxt);
    hPer[iCat]  -> GetYaxis() -> SetLabelSize(fLab[0]);
    hPer[iCat]  -> GetYaxis() -> CenterTitle(fCnt);
    hSum[iCat]  -> SetMarkerColor(fColSys[iCat]);
    hSum[iCat]  -> SetMarkerStyle(fMarU);
    hSum[iCat]  -> SetFillColor(fColSys[iCat]);
    hSum[iCat]  -> SetFillStyle(fFilU);
    hSum[iCat]  -> SetLineColor(fColSys[iCat]);
    hSum[iCat]  -> SetLineStyle(fLinS);
    hSum[iCat]  -> SetLineWidth(fWidS);
    hSum[iCat]  -> SetTitle(sTitle.Data());
    hSum[iCat]  -> SetTitleFont(fTxt);
    hSum[iCat]  -> GetXaxis() -> SetTitle(sTitleX.Data());
    hSum[iCat]  -> GetXaxis() -> SetTitleFont(fTxt);
    hSum[iCat]  -> GetXaxis() -> SetTitleSize(fTit[1]);
    hSum[iCat]  -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hSum[iCat]  -> GetXaxis() -> SetLabelFont(fTxt);
    hSum[iCat]  -> GetXaxis() -> SetLabelSize(fLab[1]);
    hSum[iCat]  -> GetXaxis() -> CenterTitle(fCnt);
    hSum[iCat]  -> GetYaxis() -> SetTitle(sTitleY.Data());
    hSum[iCat]  -> GetYaxis() -> SetTitleFont(fTxt);
    hSum[iCat]  -> GetYaxis() -> SetTitleSize(fTit[1]);
    hSum[iCat]  -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hSum[iCat]  -> GetYaxis() -> SetLabelFont(fTxt);
    hSum[iCat]  -> GetYaxis() -> SetLabelSize(fLab[1]);
    hSum[iCat]  -> GetYaxis() -> CenterTitle(fCnt);
    hSumP[iCat] -> SetMarkerColor(fColSys[iCat]);
    hSumP[iCat] -> SetMarkerStyle(fMarU);
    hSumP[iCat] -> SetFillColor(fColSys[iCat]);
    hSumP[iCat] -> SetFillStyle(fFilU);
    hSumP[iCat] -> SetLineColor(fColSys[iCat]);
    hSumP[iCat] -> SetLineStyle(fLinS);
    hSumP[iCat] -> SetLineWidth(fWidS);
    hSumP[iCat] -> SetTitle(sTitle.Data());
    hSumP[iCat] -> SetTitleFont(fTxt);
    hSumP[iCat] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hSumP[iCat] -> GetXaxis() -> SetTitleFont(fTxt);
    hSumP[iCat] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hSumP[iCat] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hSumP[iCat] -> GetXaxis() -> SetLabelFont(fTxt);
    hSumP[iCat] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hSumP[iCat] -> GetXaxis() -> CenterTitle(fCnt);
    hSumP[iCat] -> GetYaxis() -> SetTitle(sTitleP.Data());
    hSumP[iCat] -> GetYaxis() -> SetTitleFont(fTxt);
    hSumP[iCat] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hSumP[iCat] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hSumP[iCat] -> GetYaxis() -> SetLabelFont(fTxt);
    hSumP[iCat] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hSumP[iCat] -> GetYaxis() -> CenterTitle(fCnt);
  }
  hSys[NCat] -> SetMarkerColor(fColT);
  hSys[NCat] -> SetMarkerStyle(fMarT);
  hSys[NCat] -> SetFillColor(fColT);
  hSys[NCat] -> SetFillStyle(fFilT);
  hSys[NCat] -> SetLineColor(fColT);
  hSys[NCat] -> SetLineStyle(fLinT);
  hSys[NCat] -> SetLineWidth(fWidT);
  hSys[NCat] -> SetTitle(sTitle.Data());
  hSys[NCat] -> SetTitleFont(fTxt);
  hSys[NCat] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hSys[NCat] -> GetXaxis() -> SetTitleFont(fTxt);
  hSys[NCat] -> GetXaxis() -> SetTitleSize(fTit[1]);
  hSys[NCat] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
  hSys[NCat] -> GetXaxis() -> SetLabelFont(fTxt);
  hSys[NCat] -> GetXaxis() -> SetLabelSize(fLab[1]);
  hSys[NCat] -> GetXaxis() -> CenterTitle(fCnt);
  hSys[NCat] -> GetYaxis() -> SetTitle(sTitleY.Data());
  hSys[NCat] -> GetYaxis() -> SetTitleFont(fTxt);
  hSys[NCat] -> GetYaxis() -> SetTitleSize(fTit[1]);
  hSys[NCat] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
  hSys[NCat] -> GetYaxis() -> SetLabelFont(fTxt);
  hSys[NCat] -> GetYaxis() -> SetLabelSize(fLab[1]);
  hSys[NCat] -> GetYaxis() -> CenterTitle(fCnt);
  hPer[NCat] -> SetMarkerColor(fColT);
  hPer[NCat] -> SetMarkerStyle(fMarT);
  hPer[NCat] -> SetFillColor(fColT);
  hPer[NCat] -> SetFillStyle(fFilT);
  hPer[NCat] -> SetLineColor(fColT);
  hPer[NCat] -> SetLineStyle(fLinT);
  hPer[NCat] -> SetLineWidth(fWidT);
  hPer[NCat] -> SetTitle(sTitle.Data());
  hPer[NCat] -> SetTitleFont(fTxt);
  hPer[NCat] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPer[NCat] -> GetXaxis() -> SetTitleFont(fTxt);
  hPer[NCat] -> GetXaxis() -> SetTitleSize(fTit[0]);
  hPer[NCat] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
  hPer[NCat] -> GetXaxis() -> SetLabelFont(fTxt);
  hPer[NCat] -> GetXaxis() -> SetLabelSize(fLab[0]);
  hPer[NCat] -> GetXaxis() -> CenterTitle(fCnt);
  hPer[NCat] -> GetYaxis() -> SetTitle(sTitleP.Data());
  hPer[NCat] -> GetYaxis() -> SetTitleFont(fTxt);
  hPer[NCat] -> GetYaxis() -> SetTitleSize(fTit[0]);
  hPer[NCat] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
  hPer[NCat] -> GetYaxis() -> SetLabelFont(fTxt);
  hPer[NCat] -> GetYaxis() -> SetLabelSize(fLab[0]);
  hPer[NCat] -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // check if histogram has errors
  Bool_t varEffHasErr[NSysEff];
  Bool_t varResHasErr[NSysRes];
  Bool_t divEffHasErr[NSysEff];
  Bool_t divResHasErr[NSysRes];
  Bool_t sumHasErr[NCat];
  Bool_t sumPHasErr[NCat];
  Bool_t sysHasErr[NHist];
  Bool_t perHasErr[NHist];
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    varEffHasErr[iEff] = false;
    divEffHasErr[iEff] = false;
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    varResHasErr[iRes] = false;
    divResHasErr[iRes] = false;
  }
  for (UInt_t iCat = 0; iCat < NCat; iCat++) {
    sumHasErr[iCat]  = false;
    sumPHasErr[iCat] = false;
    sysHasErr[iCat]  = false;
    perHasErr[iCat]  = false;
  }
  sysHasErr[NCat] = false;
  perHasErr[NCat] = false;

  // lop over bins
  const UInt_t nBin = hSys[NCat] -> GetNbinsX();
  for (UInt_t iBin = 1; iBin < (nBin + 1); iBin++) {

    // efficiency variations
    for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
      const Double_t binCntE = hVarE[iEff] -> GetBinCenter(iBin);
      if ((binCntE < plotRange[0]) || (binCntE > plotRange[1])) continue;

      const Double_t varEffErr = hVarE[iEff] -> GetBinError(iBin);
      const Double_t divEffErr = hDivE[iEff] -> GetBinError(iBin);
      if (varEffErr > 0.) varEffHasErr[iEff] = true;
      if (divEffErr > 0.) divEffHasErr[iEff] = true;
    }

    // resolution variations
    for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
      const Double_t binCntR = hVarR[iRes] -> GetBinCenter(iBin);
      if ((binCntR < plotRange[0]) || (binCntR > plotRange[1])) continue;

      const Double_t varResErr = hVarR[iRes] -> GetBinError(iBin);
      const Double_t divResErr = hDivR[iRes] -> GetBinError(iBin);
      if (varResErr > 0.) varResHasErr[iRes] = true;
      if (divResErr > 0.) divResHasErr[iRes] = true;
    }

    // systematic histograms
    for (UInt_t iCat = 0; iCat < NCat; iCat++) {
      const Double_t binCntS = hSys[iCat] -> GetBinCenter(iBin);
      if ((binCntS < plotRange[0]) || (binCntS > plotRange[1])) continue;

      const Double_t sumErr  = hSum[iCat]  -> GetBinError(iBin);
      const Double_t sumPErr = hSumP[iCat] -> GetBinError(iBin);
      const Double_t sysErr  = hSys[iCat]  -> GetBinError(iBin);
      const Double_t perErr  = hPer[iCat]  -> GetBinError(iBin);
      if (sumErr  > 0.) sumHasErr[iCat]  = true;
      if (sumPErr > 0.) sumPHasErr[iCat] = true;
      if (sysErr  > 0.) sysHasErr[iCat]  = true;
      if (perErr  > 0.) perHasErr[iCat]  = true;
    }

    // total histograms
    const Double_t binCntT = hSys[NCat] -> GetBinCenter(iBin);
    if ((binCntT < plotRange[0]) || (binCntT > plotRange[1])) continue;

    const Double_t sysErr = hSys[NCat] -> GetBinError(iBin);
    const Double_t perErr = hPer[NCat] -> GetBinError(iBin);
    if (sysErr > 0.) sysHasErr[NCat] = true;
    if (perErr > 0.) perHasErr[NCat] = true;
  }  // end bin loop
  cout << "    Checked for errors." << endl;

  // make error bands for ratio plot
  TH1D    *hSumPR[NHist];
  TString sName[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    if (iHist < NCat) {

      // create name
      sName[iHist] = hSumP[iHist] -> GetName();
      sName[iHist].Append(sPlot.Data());

      // create histogram
      hSumPR[iHist] = (TH1D*) hSumP[iHist] -> Clone();
      hSumPR[iHist] -> SetName(sName[iHist].Data());
      hSumPR[iHist] -> GetYaxis() -> SetTitle(sTitleR.Data());
    } else {

      // create name
      sName[iHist] = hPer[iHist] -> GetName();
      sName[iHist].Append(sPlot.Data());

      // create histogram
      hSumPR[iHist] = (TH1D*) hPer[iHist] -> Clone();
      hSumPR[iHist] -> SetName(sName[iHist].Data());
      hSumPR[iHist] -> GetYaxis() -> SetTitle(sTitleR.Data());
    }
  }
  cout << "    Made uncertainty bands for variation plot." << endl;

  // make legends
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  nLinVar(NSysEff + NSysRes + NHist + 1);
  const Float_t yLinVar(0.1 + (nLinVar * 0.05));
  const Float_t fLegXY[NPlot * NPlot]  = {0.1, 0.1, 0.3, 0.5};
  const Float_t fLegVar[NPlot * NPlot] = {0.1, 0.1, 0.3, yLinVar};
  TLegend *lVar = new TLegend(fLegVar[0], fLegVar[1], fLegVar[2], fLegVar[3], sTest.Data());
  TLegend *lSys = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sTest.Data());
  TLegend *lSum = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sTest.Data());
  lVar -> SetFillColor(fColLe);
  lVar -> SetFillStyle(fFilLe);
  lVar -> SetLineColor(fColLe);
  lVar -> SetLineStyle(fLinLe);
  lVar -> SetTextFont(fTxt);
  lVar -> SetTextAlign(fAln);
  lSys -> SetFillColor(fColLe);
  lSys -> SetFillStyle(fFilLe);
  lSys -> SetLineColor(fColLe);
  lSys -> SetLineStyle(fLinLe);
  lSys -> SetTextFont(fTxt);
  lSys -> SetTextAlign(fAln);
  lSum -> SetFillColor(fColLe);
  lSum -> SetFillStyle(fFilLe);
  lSum -> SetLineColor(fColLe);
  lSum -> SetLineStyle(fLinLe);
  lSum -> SetTextFont(fTxt);
  lSum -> SetTextAlign(fAln);
  lVar -> AddEntry(hDefault, sLabelD.Data(), "pf");
  lSys -> AddEntry(hDefault, sLabelDS.Data(), "pf");
  lSum -> AddEntry(hDefault, sLabelDS.Data(), "pf");
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    lVar -> AddEntry(hVarE[iEff], sLabelE[iEff].Data(), "pf");
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    lVar -> AddEntry(hVarR[iRes], sLabelR[iRes].Data(), "pf");
  }
  lVar -> AddEntry(hEffAvg, sLabelA.Data(), "pf");
  for (UInt_t iCat = 0; iCat < NCat; iCat++) {
    lSys -> AddEntry(hSys[iCat], sLabelS[iCat].Data(), "f");
    lSum -> AddEntry(hSum[iCat], sLabelU[iCat].Data(), "f");
    lVar -> AddEntry(hSumPR[iCat], sLabelU[iCat].Data(), "f");
  }
  lSys -> AddEntry(hSys[NCat], sLabelT.Data(), "f");
  lSum -> AddEntry(hSys[NCat], sLabelT.Data(), "f");
  cout << "    Made legend." << endl;

  // make text
  const UInt_t fColTx(0);
  const UInt_t fFilTx(0);
  const UInt_t fLinTx(0);
  const Float_t fTxtXY[NPlot * NPlot] = {0.3, 0.1, 0.5, 0.3};
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

  // make line
  const UInt_t  fColLi(1);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  const Float_t fLinXY[NPlot * NPlot] = {plotRange[0], 1., plotRange[1], 1.};
  TLine *line = new TLine(fLinXY[0], fLinXY[1], fLinXY[2], fLinXY[3]);
  line -> SetLineColor(fColLi);
  line -> SetLineStyle(fLinLi);
  line -> SetLineWidth(fWidLi);
  cout << "    Made line." << endl;

  // make plots
  const UInt_t  width(750);
  const UInt_t  height(950);
  const UInt_t  fMode(0);
  const UInt_t  fBord(2);
  const UInt_t  fGrid(0);
  const UInt_t  fTick(1);
  const UInt_t  fLogX(0);
  const UInt_t  fLogY1(0);
  const UInt_t  fLogY2(1);
  const UInt_t  fFrame(0);
  const Float_t fMarginL(0.15);
  const Float_t fMarginR(0.02);
  const Float_t fMarginT1(0.005);
  const Float_t fMarginT2(0.02);
  const Float_t fMarginB1(0.25);
  const Float_t fMarginB2(0.005);
  const Float_t fPadXY1[NPlot * NPlot] = {0., 0., 1., 0.35};
  const Float_t fPadXY2[NPlot * NPlot] = {0., 0.35, 1., 1.};

  // create variations plot
  TCanvas *cPlot1 = new TCanvas("cVariations", "", width, height);
  TPad    *pPad1  = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2  = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot1       -> SetGrid(fGrid, fGrid);
  cPlot1       -> SetTicks(fTick, fTick);
  cPlot1       -> SetBorderMode(fMode);
  cPlot1       -> SetBorderSize(fBord);
  pPad1        -> SetGrid(fGrid, fGrid);
  pPad1        -> SetTicks(fTick, fTick);
  pPad1        -> SetLogx(fLogX);
  pPad1        -> SetLogy(fLogY1);
  pPad1        -> SetBorderMode(fMode);
  pPad1        -> SetBorderSize(fBord);
  pPad1        -> SetFrameBorderMode(fFrame);
  pPad1        -> SetLeftMargin(fMarginL);
  pPad1        -> SetRightMargin(fMarginR);
  pPad1        -> SetTopMargin(fMarginT1);
  pPad1        -> SetBottomMargin(fMarginB1);
  pPad2        -> SetGrid(fGrid, fGrid);
  pPad2        -> SetTicks(fTick, fTick);
  pPad2        -> SetLogx(fLogX);
  pPad2        -> SetLogy(fLogY2);
  pPad2        -> SetBorderMode(fMode);
  pPad2        -> SetBorderSize(fBord);
  pPad2        -> SetFrameBorderMode(fFrame);
  pPad2        -> SetLeftMargin(fMarginL);
  pPad2        -> SetRightMargin(fMarginR);
  pPad2        -> SetTopMargin(fMarginT2);
  pPad2        -> SetBottomMargin(fMarginB2);
  cPlot1       -> cd();
  pPad1        -> Draw();
  pPad2        -> Draw();
  pPad1        -> cd();
  hSumPR[NCat] -> Draw("E2");
  for (Int_t iPlotCat = (NCat - 1); iPlotCat >= 0; iPlotCat--) {
    hSumPR[iPlotCat] -> Draw("E2 same");
  }
  for (Int_t iPlotRes = (NSysRes - 1); iPlotRes >= 0; iPlotRes--) {
    if (divResHasErr[iPlotRes]) {
      hDivR[iPlotRes] -> Draw("hist pl same");
    }
  }
  for (Int_t iPlotEff = (NSysEff - 1); iPlotEff >= 0; iPlotEff--) {
    if (divEffHasErr[iPlotEff]) {
      hDivE[iPlotEff] -> Draw("hist pl same");
    }
  }
  hDivA    -> Draw("hist pl same");
  line     -> Draw();
  pPad2    -> cd();
  hDefault -> Draw();
  for (Int_t iPlotRes = (NSysRes - 1); iPlotRes >= 0; iPlotRes--) {
    if (varResHasErr[iPlotRes]) {
      hVarR[iPlotRes] -> Draw("hist pl same");
    }
  }
  for (Int_t iPlotEff = (NSysEff - 1); iPlotEff >= 0; iPlotEff--) {
    if (varEffHasErr[iPlotEff]) {
      hVarE[iPlotEff] -> Draw("hist pl same");
    }
  }
  hEffAvg  -> Draw("hist pl same");
  hDefault -> Draw("same");
  lVar     -> Draw();
  txt      -> Draw();
  fOut     -> cd();
  cPlot1   -> Write();
  cPlot1   -> Close();

  // make systematic plot
  TCanvas *cPlot2 = new TCanvas("cSystematics", "", width, height);
  TPad    *pPad3  = new TPad("pPad3", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad4  = new TPad("pPad4", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot2     -> SetGrid(fGrid, fGrid);
  cPlot2     -> SetTicks(fTick, fTick);
  cPlot2     -> SetBorderMode(fMode);
  cPlot2     -> SetBorderSize(fBord);
  pPad3      -> SetGrid(fGrid, fGrid);
  pPad3      -> SetTicks(fTick, fTick);
  pPad3      -> SetLogx(fLogX);
  pPad3      -> SetLogy(fLogY1);
  pPad3      -> SetBorderMode(fMode);
  pPad3      -> SetBorderSize(fBord);
  pPad3      -> SetFrameBorderMode(fFrame);
  pPad3      -> SetLeftMargin(fMarginL);
  pPad3      -> SetRightMargin(fMarginR);
  pPad3      -> SetTopMargin(fMarginT1);
  pPad3      -> SetBottomMargin(fMarginB1);
  pPad4      -> SetGrid(fGrid, fGrid);
  pPad4      -> SetTicks(fTick, fTick);
  pPad4      -> SetLogx(fLogX);
  pPad4      -> SetLogy(fLogY2);
  pPad4      -> SetBorderMode(fMode);
  pPad4      -> SetBorderSize(fBord);
  pPad4      -> SetFrameBorderMode(fFrame);
  pPad4      -> SetLeftMargin(fMarginL);
  pPad4      -> SetRightMargin(fMarginR);
  pPad4      -> SetTopMargin(fMarginT2);
  pPad4      -> SetBottomMargin(fMarginB2);
  cPlot2     -> cd();
  pPad3      -> Draw();
  pPad4      -> Draw();
  pPad3      -> cd();
  hPer[NCat] -> Draw("E2");
  for (Int_t iPlotCat = (NCat - 1); iPlotCat >= 0; iPlotCat--) {
    if (perHasErr[iPlotCat]) {
      hPer[iPlotCat] -> Draw("E2 same");
    }
  }
  hPer[NCat] -> Draw("E2 same");
  line       -> Draw();
  pPad4      -> cd();
  hSys[NCat] -> Draw("E2");
  for (Int_t iPlotCat = (NCat - 1); iPlotCat >= 0; iPlotCat--) {
    if (sysHasErr[iPlotCat]) {
      hSys[iPlotCat] -> Draw("E2 same");
    }
  }
  hSys[NCat] -> Draw("E2 same");
  hDefault   -> Draw("same");
  lSys       -> Draw();
  txt        -> Draw();
  fOut       -> cd();
  cPlot2     -> Write();
  cPlot2     -> Close();

  // make stacked systematic plot
  TCanvas *cPlot3 = new TCanvas("cStacked", "", width, height);
  TPad    *pPad5  = new TPad("pPad5", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad6  = new TPad("pPad6", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cPlot3     -> SetGrid(fGrid, fGrid);
  cPlot3     -> SetTicks(fTick, fTick);
  cPlot3     -> SetBorderMode(fMode);
  cPlot3     -> SetBorderSize(fBord);
  pPad5      -> SetGrid(fGrid, fGrid);
  pPad5      -> SetTicks(fTick, fTick);
  pPad5      -> SetLogx(fLogX);
  pPad5      -> SetLogy(fLogY1);
  pPad5      -> SetBorderMode(fMode);
  pPad5      -> SetBorderSize(fBord);
  pPad5      -> SetFrameBorderMode(fFrame);
  pPad5      -> SetLeftMargin(fMarginL);
  pPad5      -> SetRightMargin(fMarginR);
  pPad5      -> SetTopMargin(fMarginT1);
  pPad5      -> SetBottomMargin(fMarginB1);
  pPad6      -> SetGrid(fGrid, fGrid);
  pPad6      -> SetTicks(fTick, fTick);
  pPad6      -> SetLogx(fLogX);
  pPad6      -> SetLogy(fLogY2);
  pPad6      -> SetBorderMode(fMode);
  pPad6      -> SetBorderSize(fBord);
  pPad6      -> SetFrameBorderMode(fFrame);
  pPad6      -> SetLeftMargin(fMarginL);
  pPad6      -> SetRightMargin(fMarginR);
  pPad6      -> SetTopMargin(fMarginT2);
  pPad6      -> SetBottomMargin(fMarginB2);
  cPlot3     -> cd();
  pPad5      -> Draw();
  pPad6      -> Draw();
  pPad5      -> cd();
  hPer[NCat] -> Draw("E2");
  for (Int_t iPlotCat = (NCat - 1); iPlotCat >= 0; iPlotCat--) {
    if (sumPHasErr[iPlotCat]) {
      hSumP[iPlotCat] -> Draw("E2 same");
    }
  }
  line       -> Draw();
  pPad6      -> cd();
  hSys[NCat] -> Draw("E2");
  for (Int_t iPlotCat = (NCat - 1); iPlotCat >= 0; iPlotCat--) {
    if (sumHasErr[iPlotCat]) {
      hSum[iPlotCat] -> Draw("E2 same");
    }
  }
  hDefault -> Draw("same");
  lSum     -> Draw();
  txt      -> Draw();
  fOut     -> cd();
  cPlot3   -> Write();
  cPlot3   -> Close();
  cout << "    Made plots." << endl;

  // close files
  fOut     -> cd();
  hDefault -> Write();
  hEffAvg  -> Write();
  hDivA    -> Write();
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    hVarE[iEff] -> Write();
    hDivE[iEff] -> Write();
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    hVarR[iRes] -> Write();
    hDivR[iRes] -> Write();
  }
  for (UInt_t iCat = 0; iCat < NCat; iCat++) {
    hSys[iCat]  -> Write();
    hPer[iCat]  -> Write();
    hSum[iCat]  -> Write();
    hSumP[iCat] -> Write();
  }
  hSys[NCat] -> Write();
  hPer[NCat] -> Write();
  fOut       -> Close();
  fInD       -> cd();
  fInD       -> Close();
  for (UInt_t iEff = 0; iEff < NSysEff; iEff++) {
    fInE[iEff] -> cd();
    fInE[iEff] -> Close();
  }
  for (UInt_t iRes = 0; iRes < NSysRes; iRes++) {
    fInR[iRes] -> cd();
    fInR[iRes] -> Close();
  }
  cout << "  Calculating unfolding systematic uncertainty!\n" << endl;

}

// End ------------------------------------------------------------------------
