// 'CalculateUnfoldingSystematicError.C'
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
static const UInt_t NSysRegP(2);
static const UInt_t NSysPriP(3);
static const UInt_t NSysRegG(6);
static const UInt_t NSysPriG(6);

// trigger constants
static const UInt_t NPi0(3);
static const UInt_t NGam(3);
static const UInt_t NTrgs(NPi0 + NGam);

// caclulation parameters
static const Bool_t StatsAreCommon(true);
static const Bool_t RemoveExtraBins(true);

// trigger parameters
static const UInt_t NSysReg(NSysRegG);
static const UInt_t NSysPri(NSysPriG);
static const UInt_t TrigId(1);
static const UInt_t TrigBin(0);



void CalculateUnfoldingSystematicError() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Calculating unfolding systematic uncertainty..." << endl;

  // io parameters
  const TString sOut("unfoldSys.onlyBkgdVariation_pTbinHuge.et911r02gam.d29m4y2022.root");
  const TString sHistD("hUnfolded");
  const TString sHistS("hUnfolded");

  // text parameters
  const TString sSys("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sJet("anti-k_{T}, R = 0.2");
  const TString sTyp("#bf{charged jets}");
  const TString sTest("only unfolding sys. uncertainties");

  // pi0 io and hist parameters
  const TString sInPi0D("et911r02pi0_rebinUnfold/pp200r9pi0.default_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root");
  const TString sInPi0R[NSysRegP]  = {"et911r02pi0_rebinUnfold/pp200r9pi0.forRegSysM1_pTbinHuge.et911r02qt05130.p0m1k2n38t5.root",
                                      "et911r02pi0_rebinUnfold/pp200r9pi0.forRegSysP1_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root"};
  const TString sInPi0P[NSysPriP]  = {"et911r02pi0_rebinUnfold/pp200r9pi0.forDefLevySys_pTbinHuge.et911r02qt05130.p1m1k3n46t4.root",
                                      "et911r02pi0_rebinUnfold/pp200r9pi0.forAltLevySys_pTbinHuge.et911r02qt05130.p1m1k3n61t5.root",
                                      "et911r02pi0_rebinUnfold/pp200r9pi0.forGammaSys_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root"};
  const TString sRegPi0[NSysRegP]  = {"n_{iter}^{*} - 1", "n_{iter}^{*} + 1"};
  const TString sPriPi0[NSysPriP]  = {"levy: n = 4.6, t = 0.4", "levy: n = 6.1, t = 0.5", "#gamma prior"};
  const TString sNameVRP[NSysRegP] = {"hVarRegM1", "hVarRegP1"}; 
  const TString sNameVPP[NSysPriP] = {"hVarDefLevy", "hVarAltLevy", "hVarGamma"};
  const TString sNameRRP[NSysRegP] = {"hRatRegM1", "hRatRegP1"};
  const TString sNameRPP[NSysPriP] = {"hRatDefLevy", "hRatAltLevy", "hRatGamma"};

  // gamma io parameters
  const TString sInGamD("et911r02gam_puritySysCalc/pp200r9gam.default_noBkgdVariation_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root");
  const TString sInGamR[NSysRegG]  = {"et911r02gam_puritySysCalc/pp200r9gam.defaultM_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root",
                                      "et911r02gam_puritySysCalc/pp200r9gam.defaultP_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root",
                                      /*"et911r02gam_puritySysCalc/pp200r9gam.forRegSysM1_noBkgdVariation_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root",*/
                                      "et911r02gam_puritySysCalc/pp200r9gam.forRegSysM1M_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root",
                                      "et911r02gam_puritySysCalc/pp200r9gam.forRegSysM1P_pTbinHuge.et911r02qt05130.p0m1k3n38t5.root",
                                      /*"et911r02gam_puritySysCalc/pp200r9gam.forRegSysP1_noBkgdVariation_pTbinHuge.et911r02qt05130.p0m1k5n38t5.root",*/
                                      "et911r02gam_puritySysCalc/pp200r9gam.forRegSysP1M_pTbinHuge.et911r02qt05130.p0m1k5n38t5.root",
                                      "et911r02gam_puritySysCalc/pp200r9gam.forRegSysP1P_pTbinHuge.et911r02qt05130.p0m1k5n38t5.root"};
  const TString sInGamP[NSysPriG]  = {/*"et911r02gam_puritySysCalc/pp200r9gam.forDefLevySys_noBkgdVariation_pTbinHuge.et911r02qt05130.p1m1k4n118t8.root",*/
                                      "et911r02gam_puritySysCalc/pp200r9gam.forDefLevySysM_pTbinHuge.et911r02qt05130.p1m1k4n118t8.root",
                                      "et911r02gam_puritySysCalc/pp200r9gam.forDefLevySysP_pTbinHuge.et911r02qt05130.p1m1k4n118t8.root",
                                      /*"et911r02gam_puritySysCalc/pp200r9gam.forAltLevySys_noBkgdVariation_pTbinHuge.et911r02qt05130.p1m1k4n258t11.root",*/
                                      "et911r02gam_puritySysCalc/pp200r9gam.forAltLevySysM_pTbinHuge.et911r02qt05130.p1m1k4n258t11.root",
                                      "et911r02gam_puritySysCalc/pp200r9gam.forAltLevySysP_pTbinHuge.et911r02qt05130.p1m1k4n258t11.root",
                                      /*"et911r02gam_puritySysCalc/pp200r9gam.forPi0Sys_noBkgdVariation_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root",*/
                                      "et911r02gam_puritySysCalc/pp200r9gam.forPi0SysM_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root",
                                      "et911r02gam_puritySysCalc/pp200r9gam.forPi0SysP_pTbinHuge.et911r02qt05130.p0m1k4n38t5.root"};
  const TString sRegGam[NSysRegG]  = {"n_{iter}^{*} [B - #deltaB]", "n_{iter}^{*} [B + #deltaB]",
                                      /*"n_{iter}^{*} - 1 [B]",*/ "n_{iter}^{*} - 1 [B - #deltaB]", "n_{iter}^{*} - 1 [B + #deltaB]",
                                      /*"n_{iter}^{*} + 1 [B]",*/ "n_{iter}^{*} + 1 [B - #deltaB]", "n_{iter}^{*} + 1 [B + #deltaB]"};
  const TString sPriGam[NSysPriG]  = {/*"levy: n = 11.8, t = 0.8 [B]",*/ "levy: n = 11.8, t = 0.8 [B - #deltaB]", "levy: n = 11.8, t = 0.8 [B + #deltaB]",
                                      /*"levy: n = 25.8, t = 1.1 [B]",*/ "levy: n = 25.8, t = 1.1 [B - #deltaB]", "levy: n = 25.8, t = 1.1 [B + #deltaB]",
                                      /*"#pi^{0} prior [B]",*/ "#pi^{0} prior [B - #deltaB]", "#pi^{0} prior [B + #deltaB]"};
  const TString sNameVRG[NSysRegG] = {"hVarDefM", "hVarDefP", /*"hVarRegM1",*/ "hVarRegM1M", "hVarRegM1P", /*"hVarRegP1",*/ "hVarRegP1M", "hVarRegP1P"};
  const TString sNameVPG[NSysPriG] = {/*"hVarDefLevy",*/ "hVarDefLevyM", "hVarDefLevyP", /*"hVarAltLevy",*/ "hVarAltLevyM", "hVarAltLevyP", /*"hVarPi0",*/ "hVarPi0M", "hVarPi0P"};
  const TString sNameRRG[NSysRegG] = {"hRatDefM", "hRatDefP", /*"hRatRegM1",*/ "hRatRegM1M", "hRatRegM1P", /*"hRatRegP1",*/ "hRatRegP1M", "hRatRegP1P"};
  const TString sNameRPG[NSysPriG] = {/*"hRatDefLevy",*/ "hRatDefLevyM", "hRatDefLevyP", /*"hRatAltLevy",*/ "hRatAltLevyM", "hRatAltLevyP", /*"hRatPi0",*/ "hRatPi0M", "hRatPi0P"};

  // pi0 plot parameters
  const TString sTotPi0("total systematic");
  const TString sDefPi0("default [n_{iter}^{*} = 4]");
  const Float_t xCalcPi0[NPi0][NPlot]      = {{0.2, 30.}, {0.2, 30.}, {0.2, 30.}};
  const Bool_t  excludeRegP[NSysRegP]      = {false, false};
  const Bool_t  excludePriP[NSysPriP]      = {false, false, false};
  const UInt_t  fMarRegPi0[NSysRegP]       = {21, 23};
  const UInt_t  fMarPriPi0[NSysPriP]       = {22, 20, 34};
  const UInt_t  fColRegPi0[NPi0][NSysRegP] = {{859, 879},
                                              {839, 859},
                                              {819, 839}};
  const UInt_t  fColPriPi0[NPi0][NSysPriP] = {{899, 799, 819},
                                              {879, 899, 799},
                                              {859, 879, 899}};
  const UInt_t  fColSysPi0[NPi0][NHist]    = {{593, 590, 603},
                                              {425, 422, 435},
                                              {409, 406, 419}};

  // gamma plot parameters
  const TString sTotGam("total systematic");
  const TString sDefGam("default [n_{iter}^{*} = 4]");
  const Float_t xCalcGam[NGam][NPlot]      = {{0.2, 11.}, {0.2, 15.}, {0.2, 20.}};
  const Bool_t  excludeRegG[NSysRegG]      = {false, false, false, false, false, false/*, false, false*/};
  const Bool_t  excludePriG[NSysPriG]      = {false, false, false, false, false, false/*, false, false, false*/};
  const UInt_t  fMarRegGam[NSysRegG]       = {21, 25, 23, 22, 33, 34/*, 29, 43*/};
  const UInt_t  fMarPriGam[NSysPriG]       = {20, 24, 32, 26, 27, 28/*, 30, 42, 40*/};
  const UInt_t  fColRegGam[NGam][NSysRegG] = {{899, 633, 809, 799, 401, 829/*, 819, 417*/},
                                              {879, 617, 909, 899, 633, 809/*, 799, 401*/},
                                              {859, 601, 889, 879, 617, 909/*, 899, 633*/}};
  const UInt_t  fColPriGam[NGam][NSysPriG] = {{849, 839, 433, 869, 859, 601/*, 889, 879, 617*/},
                                              {829, 819, 417, 849, 839, 433/*, 869, 859, 601*/},
                                              {809, 799, 401, 829, 799, 401/*, 829, 819, 417*/}};
  const UInt_t  fColSysGam[NGam][NHist]    = {{625, 622, 635},
                                              {609, 606, 619},
                                              {593, 590, 603}};

  // general systematic parameters
  const TString sNameD("hDefault");
  const TString sNameT("hTotal");
  const TString sNamePT("hPerTotal");
  const TString sNameS[NCat]   = {"hSysReg", "hSysPri"};
  const TString sNameP[NCat]   = {"hPerReg", "hPerPri"};
  const TString sNameU[NCat]   = {"hSumReg", "hSumPri"};
  const TString sNameM[NCat]   = {"hSumPerReg", "hSumPerPri"};
  const TString sLabelS[NHist] = {"n_{iter} systematic", "prior systematic", "total systematic"};
  const TString sLabelU[NCat]  = {"#sigma_{sys}^{unf}(n_{iter})", "#sigma_{sys}^{unf}(n_{iter}) #oplus #sigma_{sys}^{unf}(prior)"};
  const UInt_t  fFilSys[NCat]  = {3345, 3354};
  const UInt_t  fMarSys[NCat]  = {1, 1};

  // general plot parameters
  const TString  sTitle("");
  const TString  sTitleX("p_{T}^{unfold} [GeV/c]");
  const TString  sTitleY("(1/N^{trg}) d^{2}N^{jet}/d(p_{T}^{unfold} #eta^{jet}) [GeV/c]^{-1}");
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
  TString sInR[NSysReg];
  TString sInP[NSysPri];
  TString sNameVR[NSysReg];
  TString sNameVP[NSysPri];
  TString sNameRR[NSysReg];
  TString sNameRP[NSysPri];
  TString sLabelR[NSysReg];
  TString sLabelP[NSysPri];
  Float_t calcRange[NPlot];
  Bool_t  excludeReg[NSysReg];
  Bool_t  excludePri[NSysPri];
  UInt_t  fColDef;
  UInt_t  fColSys[NCat];
  UInt_t  fColReg[NSysReg];
  UInt_t  fColPri[NSysPri];
  UInt_t  fMarReg[NSysReg];
  UInt_t  fMarPri[NSysPri];
  switch (fTrg) {
    case 0:
      sTrg     = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      sInD     = sInPi0D;
      sLabelD  = sDefPi0;
      sLabelDS = sDefPi0;
      sLabelT  = sTotPi0;
      fColDef  = fColSysPi0[0][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegP[iReg];
        fColReg[iReg]    = fColRegPi0[0][iReg];
        fMarReg[iReg]    = fMarRegPi0[iReg];
        sInR[iReg]       = sInPi0R[iReg];
        sNameVR[iReg]    = sNameVRP[iReg];
        sNameRR[iReg]    = sNameRRP[iReg];
        sLabelR[iReg]    = sRegPi0[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriP[iPri];
        fColPri[iPri]    = fColPriPi0[0][iPri];
        fMarPri[iPri]    = fMarPriPi0[iPri];
        sInP[iPri]       = sInPi0P[iPri];
        sNameVP[iPri]    = sNameVPP[iPri];
        sNameRP[iPri]    = sNameRPP[iPri];
        sLabelP[iPri]    = sPriPi0[iPri];
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
      fColDef  = fColSysPi0[1][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegP[iReg];
        fColReg[iReg]    = fColRegPi0[1][iReg];
        fMarReg[iReg]    = fMarRegPi0[iReg];
        sInR[iReg]       = sInPi0R[iReg];
        sNameVR[iReg]    = sNameVRP[iReg];
        sNameRR[iReg]    = sNameRRP[iReg];
        sLabelR[iReg]    = sRegPi0[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriP[iPri];
        fColPri[iPri]    = fColPriPi0[1][iPri];
        fMarPri[iPri]    = fMarPriPi0[iPri];
        sInP[iPri]       = sInPi0P[iPri];
        sNameVP[iPri]    = sNameVPP[iPri];
        sNameRP[iPri]    = sNameRPP[iPri];
        sLabelP[iPri]    = sPriPi0[iPri];
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
      fColDef  = fColSysPi0[2][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegP[iReg];
        fColReg[iReg]    = fColRegPi0[2][iReg];
        fMarReg[iReg]    = fMarRegPi0[iReg];
        sInR[iReg]       = sInPi0R[iReg];
        sNameVR[iReg]    = sNameVRP[iReg];
        sNameRR[iReg]    = sNameRRP[iReg];
        sLabelR[iReg]    = sRegPi0[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriP[iPri];
        fColPri[iPri]    = fColPriPi0[2][iPri];
        fMarPri[iPri]    = fMarPriPi0[iPri];
        sInP[iPri]       = sInPi0P[iPri];
        sNameVP[iPri]    = sNameVPP[iPri];
        sNameRP[iPri]    = sNameRPP[iPri];
        sLabelP[iPri]    = sPriPi0[iPri];
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
      fColDef  = fColSysGam[0][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegG[iReg];
        fColReg[iReg]    = fColRegGam[0][iReg];
        fMarReg[iReg]    = fMarRegGam[iReg];
        sInR[iReg]       = sInGamR[iReg];
        sNameVR[iReg]    = sNameVRG[iReg];
        sNameRR[iReg]    = sNameRRG[iReg];
        sLabelR[iReg]    = sRegGam[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriG[iPri];
        fColPri[iPri]    = fColPriGam[0][iPri];
        fMarPri[iPri]    = fMarPriGam[iPri];
        sInP[iPri]       = sInGamP[iPri];
        sNameVP[iPri]    = sNameVPG[iPri];
        sNameRP[iPri]    = sNameRPG[iPri];
        sLabelP[iPri]    = sPriGam[iPri];
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
      fColDef  = fColSysGam[1][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegG[iReg];
        fColReg[iReg]    = fColRegGam[1][iReg];
        fMarReg[iReg]    = fMarRegGam[iReg];
        sInR[iReg]       = sInGamR[iReg];
        sNameVR[iReg]    = sNameVRG[iReg];
        sNameRR[iReg]    = sNameRRG[iReg];
        sLabelR[iReg]    = sRegGam[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriG[iPri];
        fColPri[iPri]    = fColPriGam[1][iPri];
        fMarPri[iPri]    = fMarPriGam[iPri];
        sInP[iPri]       = sInGamP[iPri];
        sNameVP[iPri]    = sNameVPG[iPri];
        sNameRP[iPri]    = sNameRPG[iPri];
        sLabelP[iPri]    = sPriGam[iPri];
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
      fColDef  = fColSysGam[2][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegG[iReg];
        fColReg[iReg]    = fColRegGam[2][iReg];
        fMarReg[iReg]    = fMarRegGam[iReg];
        sInR[iReg]       = sInGamP[iReg];
        sNameVR[iReg]    = sNameVRG[iReg];
        sNameRR[iReg]    = sNameRRG[iReg];
        sLabelR[iReg]    = sRegGam[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriG[iReg];
        fColPri[iPri]    = fColPriGam[2][iPri];
        fMarPri[iPri]    = fMarPriGam[iPri];
        sInP[iPri]       = sInGamP[iPri];
        sNameVP[iPri]    = sNameVPG[iPri];
        sNameRP[iPri]    = sNameRPG[iPri];
        sLabelP[iPri]    = sPriGam[iPri];
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
      fColDef  = fColSysPi0[0][NCat];
      for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
        excludeReg[iReg] = excludeRegP[iReg];
        fColReg[iReg]    = fColRegPi0[0][iReg];
        fMarReg[iReg]    = fMarRegPi0[iReg];
        sInR[iReg]       = sInPi0P[iReg];
        sNameVR[iReg]    = sNameVRP[iReg];
        sNameRR[iReg]    = sNameRRP[iReg];
        sLabelR[iReg]    = sRegPi0[iReg];
      }
      for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
        excludePri[iPri] = excludePriP[iPri];
        fColPri[iPri]    = fColPriPi0[0][iPri];
        fMarPri[iPri]    = fMarPriPi0[iPri];
        sInP[iPri]       = sInPi0P[iPri];
        sNameVP[iPri]    = sNameVPP[iPri];
        sNameRP[iPri]    = sNameRPP[iPri];
        sLabelP[iPri]    = sPriPi0[iPri];
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

  TFile *fInR[NSysReg];
  TFile *fInP[NSysPri];
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    fInR[iReg] = new TFile(sInR[iReg].Data(), "read");
    if (!fInR[iReg]) {
      cerr << "PANIC: couldn't open a file!\n"
           << "       fInR[" << iReg <<"] = " << fInR[iReg]
           << endl;
      return;
    }
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    fInP[iPri] = new TFile(sInP[iPri].Data(), "read");
    if (!fInP[iPri]) {
      cerr << "PANIC: couldn't open a file!\n"
           << "       fInP[" << iPri <<"] = " << fInP[iPri]
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

  TH1D *hVarR[NSysReg];
  TH1D *hVarP[NSysPri];
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    hVarR[iReg] = (TH1D*) fInR[iReg] -> Get(sHistS.Data());
    if (!hVarR[iReg]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hVarR[" << iReg << "] = " << hVarR[iReg]
           << endl;
      return;
    }
    hVarR[iReg] -> SetName(sNameVR[iReg].Data());
    hVarR[iReg] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    hVarP[iPri] = (TH1D*) fInP[iPri] -> Get(sHistS.Data());
    if (!hVarP[iPri]) {
      cerr << "PANIC: couldn't grab a histogram!\n"
           << "       hVarP[" << iPri << "] = " << hVarP[iPri]
           << endl;
      return;
    }
    hVarP[iPri] -> SetName(sNameVP[iPri].Data());
    hVarP[iPri] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Grabbed histograms." << endl;

  // calculate ratios
  TH1D *hDivR[NSysReg];
  TH1D *hDivP[NSysPri];
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    hDivR[iReg] = (TH1D*) hDefault -> Clone();
    hDivR[iReg] -> SetName(sNameRR[iReg].Data());
    hDivR[iReg] -> Reset("ICE");

    const UInt_t nBinR = hDivR[iReg] -> GetNbinsX();
    for (UInt_t iBinR = 1; iBinR < (nBinR + 1); iBinR++) {
      const Double_t def  = hDefault    -> GetBinContent(iBinR);
      const Double_t var  = hVarR[iReg] -> GetBinContent(iBinR);
      const Double_t rawD = hDefault    -> GetBinError(iBinR);
      const Double_t rawV = hVarR[iReg] -> GetBinError(iBinR);
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
        hDivR[iReg] -> SetBinContent(iBinR, valD);
        hDivR[iReg] -> SetBinError(iBinR, errD);
      } else {
        hDivR[iReg] -> SetBinContent(iBinR, 0.);
        hDivR[iReg] -> SetBinError(iBinR, 0.);
      }
    }  // end bin loop
    hDivR[iReg] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    hDivP[iPri] = (TH1D*) hDefault -> Clone();
    hDivP[iPri] -> SetName(sNameRP[iPri].Data());
    hDivP[iPri] -> Reset("ICE");

    const UInt_t nBinP = hDivP[iPri] -> GetNbinsX();
    for (UInt_t iBinP = 1; iBinP < (nBinP + 1); iBinP++) {
      const Double_t def  = hDefault    -> GetBinContent(iBinP);
      const Double_t var  = hVarP[iPri] -> GetBinContent(iBinP);
      const Double_t rawD = hDefault    -> GetBinError(iBinP);
      const Double_t rawV = hVarP[iPri] -> GetBinError(iBinP);
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
        hDivP[iPri] -> SetBinContent(iBinP, valD);
        hDivP[iPri] -> SetBinError(iBinP, errD);
      } else {
        hDivP[iPri] -> SetBinContent(iBinP, 0.);
        hDivP[iPri] -> SetBinError(iBinP, 0.);
      }
    }  // end bin loop
    hDivP[iPri] -> GetXaxis() -> SetRangeUser(plotRange[0], plotRange[1]);
  }
  cout << "    Calculated ratios." << endl;

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

  // calculate regularization systematic
  for (UInt_t iBin = 1; iBin < (nBinTot + 1); iBin++) {

    // loop over regularization variations
    Double_t regSysErr(0.);
    for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {

      // calculate bin differences
      const Double_t defValR = hDefault    -> GetBinContent(iBin);
      const Double_t defErrR = hDefault    -> GetBinError(iBin);
      const Double_t regVal  = hVarR[iReg] -> GetBinContent(iBin);
      const Double_t regErr  = hVarR[iReg] -> GetBinError(iBin);
      const Double_t valDifR = TMath::Abs(defValR - regVal);
      const Double_t errDifR = TMath::Sqrt(TMath::Abs((defErrR * defErrR) - (regErr * regErr)));
      //const Bool_t   isNonzero = (regVal > 0.);
      const Bool_t   isNonzero = true;
      if (!isNonzero || excludeReg[iReg]) {
        continue;
      }
      if (!StatsAreCommon) {
        hVarR[iReg] -> SetBinError(iBin, errDifR);
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
      if (sysValR > regSysErr) {
        regSysErr = sysValR;
      }
    }  // end variation loop

    // get running totals
    const Double_t defVal = hDefault   -> GetBinContent(iBin);
    const Double_t totVal = hSys[NCat] -> GetBinContent(iBin);
    const Double_t totPer = hPer[NCat] -> GetBinError(iBin);
    const Double_t totErr = hSys[NCat] -> GetBinError(iBin);

    // update totals and assign bin values
    Double_t newTotErr = TMath::Sqrt((totErr * totErr) + (regSysErr * regSysErr));
    Double_t newTotPer = totPer;
    Double_t regSysPer = 0.;
    if ((regSysErr > 0.) && (defVal > 0.)) {
      newTotPer = newTotErr / totVal;
      regSysPer = regSysErr / defVal;
    }
    hPer[0]    -> SetBinContent(iBin, 1.);
    hPer[NCat] -> SetBinContent(iBin, 1.);
    hSumP[0]   -> SetBinContent(iBin, 1.);
    hPer[0]    -> SetBinError(iBin, regSysPer);
    hPer[NCat] -> SetBinError(iBin, newTotPer);
    hSys[0]    -> SetBinError(iBin, regSysErr);
    hSys[NCat] -> SetBinError(iBin, newTotErr);
    hSum[0]    -> SetBinError(iBin, newTotErr);
    hSumP[0]   -> SetBinError(iBin, newTotPer);

  } // end bin loop
  cout << "    Calculated regularization systematic." << endl;

  // calculate prior systematic
  for (UInt_t iBin = 1; iBin < (nBinTot + 1); iBin++) {

    // loop over prior variations
    Double_t priSysErr(0.);
    for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {

      // calculate bin differences
      const Double_t defValP   = hDefault    -> GetBinContent(iBin);
      const Double_t defErrP   = hDefault    -> GetBinError(iBin);
      const Double_t priVal    = hVarP[iPri] -> GetBinContent(iBin);
      const Double_t priErr    = hVarP[iPri] -> GetBinError(iBin);
      const Double_t valDifP   = TMath::Abs(defValP - priVal);
      const Double_t errDifP   = TMath::Sqrt(TMath::Abs((defErrP * defErrP) - (priErr * priErr)));
      //const Bool_t   isNonzero = (priVal > 0.);
      const Bool_t   isNonzero = true;
      if (!isNonzero || excludePri[iPri]) {
        continue;
      }
      if (!StatsAreCommon) {
        hVarP[iPri] -> SetBinError(iBin, errDifP);
      }

      // calculate possible systematic
      Double_t sysValP(0.);
      if (StatsAreCommon) {
        sysValP = valDifP;
      } else {
        if (valDifP > errDifP) {
          sysValP = TMath::Sqrt((valDifP * valDifP) - (errDifP * errDifP));
        } else {
          sysValP = 0.;
        }
      }

      // keep only largest difference
      if (sysValP > priSysErr) {
        priSysErr = sysValP;
      }
    }  // end variation loop

    // get running totals
    const Double_t defVal = hDefault   -> GetBinContent(iBin);
    const Double_t totVal = hSys[NCat] -> GetBinContent(iBin);
    const Double_t totPer = hPer[NCat] -> GetBinError(iBin);
    const Double_t totErr = hSys[NCat] -> GetBinError(iBin);

    // update totals and assign bin values
    Double_t newTotErr = TMath::Sqrt((totErr * totErr) + (priSysErr * priSysErr));
    Double_t newTotPer = totPer;
    Double_t priSysPer = 0.;
    if ((priSysErr > 0.) && (defVal > 0.)) {
      newTotPer = newTotErr / totVal;
      priSysPer = priSysErr / defVal;
    }
    hPer[1]    -> SetBinContent(iBin, 1.);
    hPer[NCat] -> SetBinContent(iBin, 1.);
    hSumP[1]   -> SetBinContent(iBin, 1.);
    hPer[1]    -> SetBinError(iBin, priSysPer);
    hPer[NCat] -> SetBinError(iBin, newTotPer);
    hSys[1]    -> SetBinError(iBin, priSysErr);
    hSys[NCat] -> SetBinError(iBin, newTotErr);
    hSum[1]    -> SetBinError(iBin, newTotErr);
    hSumP[1]   -> SetBinError(iBin, newTotPer);

  } // end bin loop
  cout << "    Calculated prior systematic." << endl;

  // zero out extra bins
  if (RemoveExtraBins) {
    const UInt_t nBins = hDefault -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
      const Double_t binCnt    = hDefault -> GetBinCenter(iBin);
      const Bool_t   isInRange = ((binCnt > calcRange[0]) && (binCnt < calcRange[1]));
      if (!isInRange) {
        hDefault -> SetBinContent(iBin, 0.);
        hDefault -> SetBinError(iBin, 0.);
        for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
          hVarR[iReg]  -> SetBinContent(iBin, 0.);
          hDivR[iReg]  -> SetBinContent(iBin, 0.);
          hVarR[iReg]  -> SetBinError(iBin, 0.);
          hDivR[iReg]  -> SetBinError(iBin, 0.);
        }
        for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
          hVarP[iPri]  -> SetBinContent(iBin, 0.);
          hDivP[iPri]  -> SetBinContent(iBin, 0.);
          hVarP[iPri]  -> SetBinError(iBin, 0.);
          hDivP[iPri]  -> SetBinError(iBin, 0.);
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
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    hVarR[iReg] -> SetMarkerColor(fColReg[iReg]);
    hVarR[iReg] -> SetMarkerStyle(fMarReg[iReg]);
    hVarR[iReg] -> SetFillColor(fColReg[iReg]);
    hVarR[iReg] -> SetFillStyle(fFilV);
    hVarR[iReg] -> SetLineColor(fColReg[iReg]);
    hVarR[iReg] -> SetLineStyle(fLinS);
    hVarR[iReg] -> SetLineWidth(fWidS);
    hVarR[iReg] -> SetTitle(sTitle.Data());
    hVarR[iReg] -> SetTitleFont(fTxt);
    hVarR[iReg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hVarR[iReg] -> GetXaxis() -> SetTitleFont(fTxt);
    hVarR[iReg] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hVarR[iReg] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hVarR[iReg] -> GetXaxis() -> SetLabelFont(fTxt);
    hVarR[iReg] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hVarR[iReg] -> GetXaxis() -> CenterTitle(fCnt);
    hVarR[iReg] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hVarR[iReg] -> GetYaxis() -> SetTitleFont(fTxt);
    hVarR[iReg] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hVarR[iReg] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hVarR[iReg] -> GetYaxis() -> SetLabelFont(fTxt);
    hVarR[iReg] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hVarR[iReg] -> GetYaxis() -> CenterTitle(fCnt);
    hDivR[iReg] -> SetMarkerColor(fColReg[iReg]);
    hDivR[iReg] -> SetMarkerStyle(fMarReg[iReg]);
    hDivR[iReg] -> SetFillColor(fColReg[iReg]);
    hDivR[iReg] -> SetFillStyle(fFilV);
    hDivR[iReg] -> SetLineColor(fColReg[iReg]);
    hDivR[iReg] -> SetLineStyle(fLinS);
    hDivR[iReg] -> SetLineWidth(fWidS);
    hDivR[iReg] -> SetTitle(sTitle.Data());
    hDivR[iReg] -> SetTitleFont(fTxt);
    hDivR[iReg] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDivR[iReg] -> GetXaxis() -> SetTitleFont(fTxt);
    hDivR[iReg] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hDivR[iReg] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hDivR[iReg] -> GetXaxis() -> SetLabelFont(fTxt);
    hDivR[iReg] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hDivR[iReg] -> GetXaxis() -> CenterTitle(fCnt);
    hDivR[iReg] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hDivR[iReg] -> GetYaxis() -> SetTitleFont(fTxt);
    hDivR[iReg] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hDivR[iReg] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hDivR[iReg] -> GetYaxis() -> SetLabelFont(fTxt);
    hDivR[iReg] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hDivR[iReg] -> GetYaxis() -> CenterTitle(fCnt);
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    hVarP[iPri] -> SetMarkerColor(fColPri[iPri]);
    hVarP[iPri] -> SetMarkerStyle(fMarPri[iPri]);
    hVarP[iPri] -> SetFillColor(fColPri[iPri]);
    hVarP[iPri] -> SetFillStyle(fFilV);
    hVarP[iPri] -> SetLineColor(fColPri[iPri]);
    hVarP[iPri] -> SetLineStyle(fLinS);
    hVarP[iPri] -> SetLineWidth(fWidS);
    hVarP[iPri] -> SetTitle(sTitle.Data());
    hVarP[iPri] -> SetTitleFont(fTxt);
    hVarP[iPri] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hVarP[iPri] -> GetXaxis() -> SetTitleFont(fTxt);
    hVarP[iPri] -> GetXaxis() -> SetTitleSize(fTit[1]);
    hVarP[iPri] -> GetXaxis() -> SetTitleOffset(fOffX[1]);
    hVarP[iPri] -> GetXaxis() -> SetLabelFont(fTxt);
    hVarP[iPri] -> GetXaxis() -> SetLabelSize(fLab[1]);
    hVarP[iPri] -> GetXaxis() -> CenterTitle(fCnt);
    hVarP[iPri] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hVarP[iPri] -> GetYaxis() -> SetTitleFont(fTxt);
    hVarP[iPri] -> GetYaxis() -> SetTitleSize(fTit[1]);
    hVarP[iPri] -> GetYaxis() -> SetTitleOffset(fOffY[1]);
    hVarP[iPri] -> GetYaxis() -> SetLabelFont(fTxt);
    hVarP[iPri] -> GetYaxis() -> SetLabelSize(fLab[1]);
    hVarP[iPri] -> GetYaxis() -> CenterTitle(fCnt);
    hDivP[iPri] -> SetMarkerColor(fColPri[iPri]);
    hDivP[iPri] -> SetMarkerStyle(fMarPri[iPri]);
    hDivP[iPri] -> SetFillColor(fColPri[iPri]);
    hDivP[iPri] -> SetFillStyle(fFilV);
    hDivP[iPri] -> SetLineColor(fColPri[iPri]);
    hDivP[iPri] -> SetLineStyle(fLinS);
    hDivP[iPri] -> SetLineWidth(fWidS);
    hDivP[iPri] -> SetTitle(sTitle.Data());
    hDivP[iPri] -> SetTitleFont(fTxt);
    hDivP[iPri] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hDivP[iPri] -> GetXaxis() -> SetTitleFont(fTxt);
    hDivP[iPri] -> GetXaxis() -> SetTitleSize(fTit[0]);
    hDivP[iPri] -> GetXaxis() -> SetTitleOffset(fOffX[0]);
    hDivP[iPri] -> GetXaxis() -> SetLabelFont(fTxt);
    hDivP[iPri] -> GetXaxis() -> SetLabelSize(fLab[0]);
    hDivP[iPri] -> GetXaxis() -> CenterTitle(fCnt);
    hDivP[iPri] -> GetYaxis() -> SetTitle(sTitleR.Data());
    hDivP[iPri] -> GetYaxis() -> SetTitleFont(fTxt);
    hDivP[iPri] -> GetYaxis() -> SetTitleSize(fTit[0]);
    hDivP[iPri] -> GetYaxis() -> SetTitleOffset(fOffY[0]);
    hDivP[iPri] -> GetYaxis() -> SetLabelFont(fTxt);
    hDivP[iPri] -> GetYaxis() -> SetLabelSize(fLab[0]);
    hDivP[iPri] -> GetYaxis() -> CenterTitle(fCnt);
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
  Bool_t varRegHasErr[NSysReg];
  Bool_t varPriHasErr[NSysPri];
  Bool_t divRegHasErr[NSysReg];
  Bool_t divPriHasErr[NSysPri];
  Bool_t sumHasErr[NCat];
  Bool_t sumPHasErr[NCat];
  Bool_t sysHasErr[NHist];
  Bool_t perHasErr[NHist];
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    varRegHasErr[iReg] = false;
    divRegHasErr[iReg] = false;
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    varPriHasErr[iPri] = false;
    divPriHasErr[iPri] = false;
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

    // regularization variations
    for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
      const Double_t binCntR = hVarR[iReg] -> GetBinCenter(iBin);
      if ((binCntR < plotRange[0]) || (binCntR > plotRange[1])) continue;

      const Double_t varRegErr = hVarR[iReg] -> GetBinError(iBin);
      const Double_t divRegErr = hDivR[iReg] -> GetBinError(iBin);
      if (varRegErr > 0.) varRegHasErr[iReg] = true;
      if (divRegErr > 0.) divRegHasErr[iReg] = true;
    }

    // prior variations
    for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
      const Double_t binCntP = hVarP[iPri] -> GetBinCenter(iBin);
      if ((binCntP < plotRange[0]) || (binCntP > plotRange[1])) continue;

      const Double_t varPriErr = hVarP[iPri] -> GetBinError(iBin);
      const Double_t divPriErr = hDivP[iPri] -> GetBinError(iBin);
      if (varPriErr > 0.) varPriHasErr[iPri] = true;
      if (divPriErr > 0.) divPriHasErr[iPri] = true;
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
  const UInt_t  nLinVar(NSysReg + NSysPri + NHist);
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
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    lVar -> AddEntry(hVarR[iReg], sLabelR[iReg].Data(), "pf");
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    lVar -> AddEntry(hVarP[iPri], sLabelP[iPri].Data(), "pf");
  }
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
  for (Int_t iPlotPri = (NSysPri - 1); iPlotPri >= 0; iPlotPri--) {
    if (divPriHasErr[iPlotPri]) {
      hDivP[iPlotPri] -> Draw("hist pl same");
    }
  }
  for (Int_t iPlotReg = (NSysReg - 1); iPlotReg >= 0; iPlotReg--) {
    if (divRegHasErr[iPlotReg]) {
      hDivR[iPlotReg] -> Draw("hist pl same");
    }
  }
  line     -> Draw();
  pPad2    -> cd();
  hDefault -> Draw();
  for (Int_t iPlotPri = (NSysPri - 1); iPlotPri >= 0; iPlotPri--) {
    if (varPriHasErr[iPlotPri]) {
      hVarP[iPlotPri] -> Draw("hist pl same");
    }
  }
  for (Int_t iPlotReg = (NSysReg - 1); iPlotReg >= 0; iPlotReg--) {
    if (varRegHasErr[iPlotReg]) {
      hVarR[iPlotReg] -> Draw("hist pl same");
    }
  }
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
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    hVarR[iReg] -> Write();
    hDivR[iReg] -> Write();
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    hVarP[iPri] -> Write();
    hDivP[iPri] -> Write();
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
  for (UInt_t iReg = 0; iReg < NSysReg; iReg++) {
    fInR[iReg] -> cd();
    fInR[iReg] -> Close();
  }
  for (UInt_t iPri = 0; iPri < NSysPri; iPri++) {
    fInP[iPri] -> cd();
    fInP[iPri] -> Close();
  }
  cout << "  Calculating unfolding systematic uncertainty!\n" << endl;

}

// End ------------------------------------------------------------------------
