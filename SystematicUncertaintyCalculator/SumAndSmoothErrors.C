// 'SumAndSmoothErrors.C'
// Derek Anderson
// 04.28.2022
//
// Sum up some histograms and their
// errors why don'cha?  There should be NSys
// + 2 histograms: NSys histograms with
// systematic errors, 1 histogram with
// statistical errors, and the last one will
// be the sum of the two.
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
#include "TPad.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// global constants
static const UInt_t  NSys(2);
static const UInt_t  NPts(2);
static const UInt_t  NVtx(4);
static const UInt_t  NTrgs(6);
static const UInt_t  NHist(3.);
static const UInt_t  NMaxBin(100);
static const Float_t XStart(0.2);
static const Float_t XStop(30.);

// rebinning parameters
static const Bool_t   DoRebin(false);
static const UInt_t   NPtRebin(21);
static const Double_t PtRebin[NPtRebin + 1] = {-12., -8., -5., -3., -2., -1.5, -1., 0., 1., 2., 3., 5., 8., 12., 17., 23., 30., 38., 47., 57., 68., 81.};

// trigger parameters
static const UInt_t TrigId(0);
static const UInt_t TrigBin(0);



void SumAndSmoothErrors() {

  // lower verbosity
  gErrorIgnoreLevel = kError;
  cout << "\n  Summing errors..." << endl;

  // input/output parameters
  const TString sOut("summedErrorsFF.modStats_forClosureTest_pTbinHuge.et911r05pi0.d26m10y2021.root");
  const TString sStat("et911r05ff_rebinClosure/pp200r9ff.modifiedStats_parLvlFF_pTbinHuge.et911r05pi0.d4m10y2021.root");
  const TString sSys[NSys] = {"et911r05pi0_rebinUnfold/unfoldSys.noGamSys_pTbinHuge.et911r05pi0.d7m10y2021.root", "et911r05pi0_rebinUnfold/detectorSys.updatedStyle_pTbinHuge.et911r05pi0.d7m10y2021.root"};

  // input histogram parameters
  const TString sHistStat("hSimJetOutput");
  const TString sHistSys[NSys]   = {"hTotal", "hTotal"};
  const TString sNameSys[NSys]   = {"hSysUnfold", "hSysDet"};
  const TString sPlotSys[NSys]   = {"hPlotSysUnfold", "hPlotSysDet"};
  const TString sPlotTot[NSys]   = {"hPlotTotUnfold", "hPlotTotDet"};
  const TString sLineSys[NSys]   = {"hLineSysUnfold", "hLineSysDet"};
  const TString sLineTot[NSys]   = {"hLineTotUnfold", "hLineTotDet"};
  const Float_t xPlotRange[NPts] = {0., 30.};

  // text parameters
  const TString sCol("pp-collisions, #sqrt{s} = 200 GeV");
  const TString sJet("anti-k_{T}, R = 0.5");
  const TString sTyp("#bf{charged jets}");

  // histogram parameters
  const TString sTitle("");
  const TString sTitleX("p_{T}^{unfold} [GeV/c]");
  const TString sTitleY("(1/N^{trg}) d^{3}N^{jet}/d(p_{T}^{unfold} #eta^{jet}) [GeV/c]^{-1}");
  const TString sTitleL("percent uncertainty");
  const TString sMixSuff("Mix");
  const TString sSmoothSuff("Smooth");
  const TString sLabelStat("#sigma_{stat}");
  const TString sLabelSys[NSys] = {"#sigma_{sys}^{unfold}", "#sigma_{sys}^{unfold} #oplus #sigma_{sys}^{det}"};
  const TString sLabels[NHist]  = {"only #sigma_{stat}", "only #sigma_{sys}", "total uncertainty [#sigma_{stat} vs. #sigma_{sys}]"};
  const Float_t fScales[NHist]  = {1., 10., 100.};
  const Bool_t  useSmooth[NSys] = {true, false};

  // trigger parameters (color schemes and smoothing/interpolation switch)
  const Float_t xInterpol[NTrgs]    = {15., 15., 15., 4., 6., 11.};
  const UInt_t  fColSt[NTrgs]       = {859, 839, 819, 899, 879, 859};
  const UInt_t  fColSy[NTrgs][NSys] = {{593, 590}, {425, 422}, {409, 406}, {625, 622}, {609, 606}, {593, 590}};

  // open files
  TFile *fOut  = new TFile(sOut.Data(), "recreate");
  TFile *fStat = new TFile(sStat.Data(), "read");
  if (!fStat) {
    cerr << "PANIC: couldn't open statistics file!" << endl;
    return;
  }

  TFile *fSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fSys[iSys] = new TFile(sSys[iSys].Data(), "read");
    if (!fSys[iSys]) {
      cerr << "PANIC: couldn't systematic file " << iSys << endl;
      return;
    }
  }
  cout << "    Opened files." << endl;

  // grab histograms
  TH1D *hStat = (TH1D*) fStat -> Get(sHistStat.Data());
  if (!hStat) {
    cerr << "PANIC: couldn't grab statistics histogram!" << endl;
    return;
  }

  TH1D *hSys[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys] = (TH1D*) fSys[iSys] -> Get(sHistSys[iSys].Data());
    if (!hSys[iSys]) {
      cerr << "PANIC: couldn't grab systematic histogram " << iSys << endl;
      return;
    }
    hSys[iSys] -> SetName(sNameSys[iSys].Data());
  }
  cout << "    Grabbed histograms." << endl;

  // rebin histograms
  if (DoRebin) {
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      hSys[iSys] = (TH1D*) hSys[iSys] -> Rebin(NPtRebin, sNameSys[iSys].Data(), PtRebin);
    }
    hStat = (TH1D*) hStat -> Rebin(NPtRebin, sHistStat.Data(), PtRebin);
    cout << "    Rebinned histograms." << endl;
  }

  // for smoothed sys. uncertainties
  TH1D *hSysSmooth[NSys];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {

    // generate names
    TString sNameS = hSys[iSys] -> GetName();
    sNameS.Append(sSmoothSuff.Data());

    // initialize histograms
    hSysSmooth[iSys] = (TH1D*) hSys[iSys] -> Clone();
    hSysSmooth[iSys] -> SetName(sNameS.Data());
  }
  cout << "    Created histograms for smoothing." << endl;

  // parse trigger selection
  UInt_t fTrg(0);
  if (TrigId == 0) {
    fTrg = TrigBin;
  } else {
    fTrg = TrigBin + 3;
  }

  // select trigger
  UInt_t  fColStat(0);
  UInt_t  fColSys[NSys];
  Float_t xIntBegin(0.);
  TString sTrg("");
  switch (fTrg) {
    case 0:
      sTrg      = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      fColStat  = fColSt[0];
      xIntBegin = xInterpol[0];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[0][iSys];
      }
      break;
    case 1:
      sTrg      = "#pi^{0} trig., E_{T}^{meas} #in (11, 15) GeV";
      fColStat  = fColSt[1];
      xIntBegin = xInterpol[1];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[1][iSys];
      }
      break;
    case 2:
      sTrg      = "#pi^{0} trig., E_{T}^{meas} #in (15, 20) GeV";
      fColStat  = fColSt[2];
      xIntBegin = xInterpol[2];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[2][iSys];
      }
      break;
    case 3:
      sTrg      = "#gamma_{dir} trig., E_{T}^{meas} #in (9, 11) GeV";
      fColStat  = fColSt[3];
      xIntBegin = xInterpol[3];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[3][iSys];
      }
      break;
    case 4:
      sTrg      = "#gamma_{dir} trig., E_{T}^{meas} #in (11, 15) GeV";
      fColStat  = fColSt[4];
      xIntBegin = xInterpol[4];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[4][iSys];
      }
      break;
    case 5:
      sTrg      = "#gamma_{dir} trig., E_{T}^{meas} #in (15, 20) GeV";
      fColStat  = fColSt[5];
      xIntBegin = xInterpol[5];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[5][iSys];
      }
      break;
    default:
      sTrg      = "#pi^{0} trig., E_{T}^{meas} #in (9, 11) GeV";
      fColStat  = fColSt[0];
      xIntBegin = xInterpol[0];
      for (UInt_t iSys = 0; iSys < NSys; iSys++) {
        fColSys[iSys] = fColSy[0][iSys];
      }
      break;
  }  // end switch-case
  cout << "    Selected trigger." << endl;

  // create sum histograms
  TH1D *hStatOut = (TH1D*) hStat -> Clone();
  TH1D *hSysOut  = (TH1D*) hStat -> Clone();
  TH1D *hTotOut  = (TH1D*) hStat -> Clone();
  hStatOut -> SetName("hStatistics");
  hSysOut  -> SetName("hSystematics");
  hTotOut  -> SetName("hTotal");
  hStatOut -> Reset("ICE");
  hSysOut  -> Reset("ICE");
  hTotOut  -> Reset("ICE");

  // create plotting histograms
  TH1D *hPlotStat;
  TH1D *hLineStat;
  TH1D *hPlotSys[NSys];
  TH1D *hLineSys[NSys];
  TH1D *hPlotTot[NSys + 1];
  TH1D *hLineTot[NSys + 1];

  hPlotStat      = (TH1D*) hStat -> Clone();
  hLineStat      = (TH1D*) hStat -> Clone();
  hPlotTot[NSys] = (TH1D*) hStat -> Clone();
  hLineTot[NSys] = (TH1D*) hStat -> Clone();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hPlotSys[iSys] = (TH1D*) hSys[iSys] -> Clone();
    hLineSys[iSys] = (TH1D*) hSys[iSys] -> Clone();
    hPlotTot[iSys] = (TH1D*) hSys[iSys] -> Clone();
    hLineTot[iSys] = (TH1D*) hSys[iSys] -> Clone();
    hPlotSys[iSys] -> SetName(sPlotSys[iSys].Data());
    hLineSys[iSys] -> SetName(sLineSys[iSys].Data());
    hPlotTot[iSys] -> SetName(sPlotTot[iSys].Data());
    hLineTot[iSys] -> SetName(sLineTot[iSys].Data());
    hPlotSys[iSys] -> Reset("ICE");
    hLineSys[iSys] -> Reset("ICE");
    hPlotTot[iSys] -> Reset("ICE");
    hLineTot[iSys] -> Reset("ICE");
  }
  hPlotStat      -> SetName("hPlotStat");
  hLineStat      -> SetName("hLineStat");
  hPlotTot[NSys] -> SetName("hPlotTotStat");
  hLineTot[NSys] -> SetName("hLineTotStat");
  hPlotStat      -> Reset("ICE");
  hLineStat      -> Reset("ICE");
  hPlotTot[NSys] -> Reset("ICE");
  hLineTot[NSys] -> Reset("ICE");

  // for smoothed and mix sys. uncertainties
  TH1D *hSysOutMix;
  TH1D *hTotOutMix;
  TH1D *hSysOutSmooth;
  TH1D *hTotOutSmooth;
  TH1D *hPlotSysMix[NSys];
  TH1D *hPlotTotMix[NSys + 1];
  TH1D *hPlotSysSmooth[NSys];
  TH1D *hPlotTotSmooth[NSys + 1];
  TH1D *hLineSysMix[NSys];
  TH1D *hLineTotMix[NSys + 1];
  TH1D *hLineSysSmooth[NSys];
  TH1D *hLineTotSmooth[NSys + 1];
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {

    // generate names
    TString sNamePSM = hPlotSys[iSys] -> GetName();
    TString sNamePSS = hPlotSys[iSys] -> GetName();
    TString sNamePTM = hPlotTot[iSys] -> GetName();
    TString sNamePTS = hPlotTot[iSys] -> GetName();
    TString sNameLSM = hLineSys[iSys] -> GetName();
    TString sNameLSS = hLineSys[iSys] -> GetName();
    TString sNameLTM = hLineTot[iSys] -> GetName();
    TString sNameLTS = hLineTot[iSys] -> GetName();
    sNamePSM.Append(sMixSuff.Data());
    sNamePSS.Append(sSmoothSuff.Data());
    sNamePTM.Append(sMixSuff.Data());
    sNamePTS.Append(sSmoothSuff.Data());
    sNameLSM.Append(sMixSuff.Data());
    sNameLSS.Append(sSmoothSuff.Data());
    sNameLTM.Append(sMixSuff.Data());
    sNameLTS.Append(sSmoothSuff.Data());

    // initialize histograms
    hPlotSysMix[iSys]    = (TH1D*) hPlotSys[iSys] -> Clone();
    hPlotSysSmooth[iSys] = (TH1D*) hPlotSys[iSys] -> Clone();
    hPlotTotMix[iSys]    = (TH1D*) hPlotTot[iSys] -> Clone();
    hPlotTotSmooth[iSys] = (TH1D*) hPlotTot[iSys] -> Clone();
    hLineSysMix[iSys]    = (TH1D*) hLineSys[iSys] -> Clone();
    hLineSysSmooth[iSys] = (TH1D*) hLineSys[iSys] -> Clone();
    hLineTotMix[iSys]    = (TH1D*) hLineTot[iSys] -> Clone();
    hLineTotSmooth[iSys] = (TH1D*) hLineTot[iSys] -> Clone();
    hPlotSysMix[iSys]    -> SetName(sNamePSM.Data());
    hPlotSysSmooth[iSys] -> SetName(sNamePSS.Data());
    hPlotTotMix[iSys]    -> SetName(sNamePTM.Data());
    hPlotTotSmooth[iSys] -> SetName(sNamePTS.Data());
    hLineSysMix[iSys]    -> SetName(sNameLSM.Data());
    hLineSysSmooth[iSys] -> SetName(sNameLSS.Data());
    hLineTotMix[iSys]    -> SetName(sNameLTM.Data());
    hLineTotSmooth[iSys] -> SetName(sNameLTS.Data());
  }

  // generate names
  TString sNameSM   = hSysOut        -> GetName();
  TString sNameSS   = hSysOut        -> GetName();
  TString sNameTM   = hTotOut        -> GetName();
  TString sNameTS   = hTotOut        -> GetName();
  TString sNamePTMn = hPlotTot[NSys] -> GetName();
  TString sNamePTSn = hPlotTot[NSys] -> GetName();
  TString sNameLTMn = hLineTot[NSys] -> GetName();
  TString sNameLTSn = hLineTot[NSys] -> GetName();
  sNameSM.Append(sMixSuff.Data());
  sNameSS.Append(sSmoothSuff.Data());
  sNameTM.Append(sMixSuff.Data());
  sNameTS.Append(sSmoothSuff.Data());
  sNamePTMn.Append(sMixSuff.Data());
  sNamePTSn.Append(sSmoothSuff.Data());
  sNameLTMn.Append(sMixSuff.Data());
  sNameLTSn.Append(sSmoothSuff.Data());

  // initialize histograms
  hSysOutMix           = (TH1D*) hSysOut        -> Clone();
  hSysOutSmooth        = (TH1D*) hSysOut        -> Clone();
  hTotOutMix           = (TH1D*) hTotOut        -> Clone();
  hTotOutSmooth        = (TH1D*) hTotOut        -> Clone();
  hPlotTotMix[NSys]    = (TH1D*) hPlotTot[NSys] -> Clone();
  hPlotTotSmooth[NSys] = (TH1D*) hPlotTot[NSys] -> Clone();
  hLineTotMix[NSys]    = (TH1D*) hLineTot[NSys] -> Clone();
  hLineTotSmooth[NSys] = (TH1D*) hLineTot[NSys] -> Clone();
  hSysOutMix           -> SetName(sNameSM.Data());
  hSysOutSmooth        -> SetName(sNameSS.Data());
  hTotOutMix           -> SetName(sNameTM.Data());
  hTotOutSmooth        -> SetName(sNameTS.Data());
  hPlotTotMix[NSys]    -> SetName(sNamePTMn.Data());
  hPlotTotSmooth[NSys] -> SetName(sNamePTSn.Data());
  hLineTotMix[NSys]    -> SetName(sNameLTMn.Data());
  hLineTotSmooth[NSys] -> SetName(sNameLTSn.Data());
  cout << "    Created histograms." << endl;

  // determine range
  const Int_t  iXstart = hStat -> FindBin(XStart);
  const UInt_t iFirst  = hStat -> FindFirstBinAbove(0.);
  const UInt_t iLast   = hStat -> FindLastBinAbove(0.);
  const UInt_t iStop   = hStat -> GetNbinsX();

  // set starting bin
  UInt_t iStart(0);
  if (iFirst < iXstart) {
    iStart = iXstart;
  } else {
    iStart = iFirst;
  }

  // count up bins
  UInt_t  nBins(0);
  UInt_t  iCount(iStart);
  Bool_t  isInRange(true);
  Float_t binCenter(0.);
  do {
    binCenter = hStat -> GetBinCenter(iCount);
    isInRange = ((binCenter >= XStart) && (binCenter < XStop));
    ++iCount;
    if (isInRange) ++nBins;
  } while (isInRange);
  cout << "    Range and number of bins determined:\n"
       << "      range = (" << XStart << ", " << XStop << "), nBins = " << nBins
       << endl;

  // smooth systematics
  const UInt_t nSmooth = hSysSmooth[0] -> GetNbinsX();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    for (UInt_t iSmooth = 2; iSmooth < nSmooth; iSmooth++) {

      // get bin info
      const Double_t thisBinLoc = hSysSmooth[iSys] -> GetBinCenter(iSmooth);
      const Double_t thisBinVal = hSysSmooth[iSys] -> GetBinContent(iSmooth);
      const Double_t prevBinVal = hSysSmooth[iSys] -> GetBinContent(iSmooth - 1);
      const Double_t nextBinVal = hSysSmooth[iSys] -> GetBinContent(iSmooth + 1);
      const Double_t thisBinAbs = hSysSmooth[iSys] -> GetBinError(iSmooth);
      const Double_t prevBinAbs = hSysSmooth[iSys] -> GetBinError(iSmooth - 1);
      const Double_t nextBinAbs = hSysSmooth[iSys] -> GetBinError(iSmooth + 1);
      const Double_t thisBinErr = thisBinAbs / thisBinVal;
      const Double_t prevBinErr = prevBinAbs / prevBinVal;
      const Double_t nextBinErr = nextBinAbs / nextBinVal;

      // check bins
      const Bool_t areBinsNonzero = ((thisBinVal > 0.) && (nextBinVal > 0.));
      const Bool_t isInSmoothZone = (thisBinLoc < xIntBegin); 
      const Bool_t isAboveStart   = (iSmooth >= iStart);
      const Bool_t isBelowStop    = ((iSmooth <= iStop) && (iSmooth <= iLast));
      const Bool_t isInCalcRange  = (isAboveStart && isBelowStop);
      if (!areBinsNonzero || !isInCalcRange) continue;

      // check how sys compare
      const Bool_t prevMoreThanThis = (prevBinErr > thisBinErr);
      const Bool_t prevMoreThanNext = (prevBinErr > nextBinErr);
      const Bool_t thisMoreThanBoth = ((thisBinErr > prevBinErr) && (thisBinErr > nextBinErr));

      // if (x < xInterpol) smooth; else, interpolate
      if (isInSmoothZone) {
        if (prevMoreThanThis) {
          const Double_t newThisErr = thisBinVal * prevBinErr;
          hSysSmooth[iSys] -> SetBinError(iSmooth, newThisErr);
        }
      } else {
        if (thisMoreThanBoth) {
          const Double_t newNextErr = nextBinVal * thisBinErr;
          hSysSmooth[iSys] -> SetBinError(iSmooth + 1, newNextErr); 
        } else if (prevMoreThanNext) {
          const Double_t newThisErr = thisBinVal * prevBinErr;
          const Double_t newNextErr = nextBinVal * prevBinErr;
          hSysSmooth[iSys] -> SetBinError(iSmooth, newThisErr);
          hSysSmooth[iSys] -> SetBinError(iSmooth + 1, newNextErr);
        } else if (prevMoreThanThis) {
          const Double_t errorDiff  = nextBinErr - prevBinErr;
          const Double_t errorAdj   = errorDiff / 2.;
          const Double_t thisBinAdj = prevBinErr + errorAdj;
          const Double_t newThisErr = thisBinVal * thisBinAdj;
          hSysSmooth[iSys] -> SetBinError(iSmooth, newThisErr);
        }
      }  // end smoothing/interpolating
    }  // end bin loop 
  }  // end sys loop
  cout << "    Smoothed histograms." << endl;

  // initialize errors
  Double_t valStat[NMaxBin];
  Double_t errStat[NMaxBin];
  Double_t perStat[NMaxBin];
  Double_t perSys2R[NMaxBin];
  Double_t perSys2M[NMaxBin];
  Double_t perSys2S[NMaxBin];
  Double_t perTot2R[NMaxBin];
  Double_t perTot2M[NMaxBin];
  Double_t perTot2S[NMaxBin];
  Double_t perSumSys2R[NSys][NMaxBin];
  Double_t perSumSys2M[NSys][NMaxBin];
  Double_t perSumSys2S[NSys][NMaxBin];
  Double_t perSumTot2R[NSys + 1][NMaxBin];
  Double_t perSumTot2M[NSys + 1][NMaxBin];
  Double_t perSumTot2S[NSys + 1][NMaxBin];
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {
    const UInt_t iVal   = iBin + iStart;
    const Bool_t isLast = (iVal > iLast);
    const Bool_t isDone = (iVal > iStop);
    if (isLast || isDone) break;

    const Double_t val  = hStat -> GetBinContent(iVal);
    const Double_t stat = hStat -> GetBinError(iVal);
    const Double_t per  = stat / val;
    valStat[iBin]  = val;
    errStat[iBin]  = stat;
    perStat[iBin]  = per;
    perSys2R[iBin] = 0.;
    perSys2M[iBin] = 0.;
    perSys2S[iBin] = 0.;
    perTot2R[iBin] = 0.;
    perTot2M[iBin] = 0.;
    perTot2S[iBin] = 0.;
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      perSumSys2R[iSys][iBin] = 0.;
      perSumSys2M[iSys][iBin] = 0.;
      perSumSys2S[iSys][iBin] = 0.;
      perSumTot2R[iSys][iBin] = 0.;
      perSumTot2M[iSys][iBin] = 0.;
      perSumTot2S[iSys][iBin] = 0.;
    }
    perSumTot2R[NSys][iBin] = per * per;
    perSumTot2M[NSys][iBin] = per * per;
    perSumTot2S[NSys][iBin] = per * per;
  }
  cout << "    Initialized errors." << endl;

  // loop over systematics
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {

    // set systematic start
    UInt_t iFirstSys = hSys[iSys] -> FindFirstBinAbove(0.);
    UInt_t iCheck    = iFirstSys;
    if (iCheck < iXstart) {
      iCheck = iXstart;
    }
    if (iCheck < iStart) {
      iCheck = iStart;
    }

    // check start point
    if (iStart != iCheck) {
      cerr << "WARNING: looks like bin numbering is off in systematic no. " << iSys << "!\n"
           << "         iStart = " << iStart << ", iCheck = " << iCheck
           << endl;
    }

    // loop over bins
    for (UInt_t iBin = 0; iBin < nBins; iBin++) {
      const UInt_t iVal   = iBin + iStart;
      const Bool_t isLast = (iVal > iLast);
      const Bool_t isDone = (iVal > iStop);
      if (isLast || isDone) break;

      // add uncertainty to total sum
      const Double_t valAddR = hSys[iSys]       -> GetBinContent(iVal);
      const Double_t valAddS = hSysSmooth[iSys] -> GetBinContent(iVal);
      const Double_t sysAddR = hSys[iSys]       -> GetBinError(iVal);
      const Double_t sysAddS = hSysSmooth[iSys] -> GetBinError(iVal);
      const Double_t perAddR = sysAddR / valAddR;
      const Double_t perAddS = sysAddS / valAddS;
      if (valAddR > 0.) {
        perSys2R[iBin] += perAddR * perAddR;
        perTot2R[iBin] += perAddR * perAddR;
      }
      if (valAddS > 0.) {
        perSys2S[iBin] += perAddS * perAddS;
        perTot2S[iBin] += perAddS * perAddS;
      }

      // add to mix uncertainty
      Double_t valAddM(0.);
      Double_t sysAddM(0.);
      if (useSmooth[iSys]) {
        valAddM = hSysSmooth[iSys] -> GetBinContent(iVal);
        sysAddM = hSysSmooth[iSys] -> GetBinError(iVal);
      } else {
        valAddM = hSys[iSys] -> GetBinContent(iVal);
        valAddM = hSys[iSys] -> GetBinError(iVal);
      }

      const Double_t perAddM = sysAddM / valAddM;
      if (valAddM > 0.) {
        perSys2M[iBin] += perAddM * perAddM;
        perTot2M[iBin] += perAddM * perAddM;
      }

      // add uncertainty to breakdown of sums
      if (iSys > 0) {
        perSumSys2R[iSys][iBin] =  perSumSys2R[iSys - 1][iBin];
        perSumSys2M[iSys][iBin] =  perSumSys2M[iSys - 1][iBin];
        perSumSys2S[iSys][iBin] =  perSumSys2S[iSys - 1][iBin];
        perSumTot2R[iSys][iBin] =  perSumTot2R[iSys - 1][iBin];
        perSumTot2M[iSys][iBin] =  perSumTot2M[iSys - 1][iBin];
        perSumTot2S[iSys][iBin] =  perSumTot2S[iSys - 1][iBin];
        perSumSys2R[iSys][iBin] += perAddR * perAddR;
        perSumSys2M[iSys][iBin] += perAddM * perAddM;
        perSumSys2S[iSys][iBin] += perAddS * perAddS;
        perSumTot2R[iSys][iBin] += perAddR * perAddR;
        perSumTot2M[iSys][iBin] += perAddM * perAddM;
        perSumTot2S[iSys][iBin] += perAddS * perAddS;
      } else {
        perSumSys2R[iSys][iBin] += perAddR * perAddR;
        perSumSys2M[iSys][iBin] += perAddM * perAddM;
        perSumSys2S[iSys][iBin] += perAddS * perAddS;
        perSumTot2R[iSys][iBin] += perAddR * perAddR;
        perSumTot2M[iSys][iBin] += perAddM * perAddM;
        perSumTot2S[iSys][iBin] += perAddS * perAddS;
      }
      perSumTot2R[NSys][iBin] = perStat[iBin] * perStat[iBin];
      perSumTot2M[NSys][iBin] = perStat[iBin] * perStat[iBin];
      perSumTot2S[NSys][iBin] = perStat[iBin] * perStat[iBin];
    }  // end bin loop
  }  // end systematic loop
  cout << "    Summed errors." << endl;

  Double_t perSysR[NMaxBin];
  Double_t perSysM[NMaxBin];
  Double_t perSysS[NMaxBin];
  Double_t errSysR[NMaxBin];
  Double_t errSysM[NMaxBin];
  Double_t errSysS[NMaxBin];
  Double_t perTotR[NMaxBin];
  Double_t perTotM[NMaxBin];
  Double_t perTotS[NMaxBin];
  Double_t errTotR[NMaxBin];
  Double_t errTotM[NMaxBin];
  Double_t errTotS[NMaxBin];
  Double_t perSumSysR[NSys][NMaxBin];
  Double_t perSumSysM[NSys][NMaxBin];
  Double_t perSumSysS[NSys][NMaxBin];
  Double_t errSumSysR[NSys][NMaxBin];
  Double_t errSumSysM[NSys][NMaxBin];
  Double_t errSumSysS[NSys][NMaxBin];
  Double_t perSumTotR[NSys + 1][NMaxBin];
  Double_t perSumTotM[NSys + 1][NMaxBin];
  Double_t perSumTotS[NSys + 1][NMaxBin];
  Double_t errSumTotR[NSys + 1][NMaxBin];
  Double_t errSumTotM[NSys + 1][NMaxBin];
  Double_t errSumTotS[NSys + 1][NMaxBin];
  for (UInt_t iBin = 0; iBin < nBins; iBin++) {

    // set total sum values
    perSysR[iBin] = TMath::Sqrt(perSys2R[iBin]);
    perSysM[iBin] = TMath::Sqrt(perSys2M[iBin]);
    perSysS[iBin] = TMath::Sqrt(perSys2S[iBin]);
    perTotR[iBin] = TMath::Sqrt(perTot2R[iBin]);
    perTotM[iBin] = TMath::Sqrt(perTot2M[iBin]);
    perTotS[iBin] = TMath::Sqrt(perTot2S[iBin]);
    errSysR[iBin] = perSysR[iBin] * valStat[iBin];
    errSysM[iBin] = perSysM[iBin] * valStat[iBin];
    errSysS[iBin] = perSysS[iBin] * valStat[iBin];
    errTotR[iBin] = perTotR[iBin] * valStat[iBin];
    errTotM[iBin] = perTotM[iBin] * valStat[iBin];
    errTotS[iBin] = perTotS[iBin] * valStat[iBin];
    hStatOut      -> SetBinContent(iBin + iStart, valStat[iBin]);
    hSysOut       -> SetBinContent(iBin + iStart, valStat[iBin]);
    hSysOutMix    -> SetBinContent(iBin + iStart, valStat[iBin]);
    hSysOutSmooth -> SetBinContent(iBin + iStart, valStat[iBin]);
    hTotOut       -> SetBinContent(iBin + iStart, valStat[iBin]);
    hTotOutMix    -> SetBinContent(iBin + iStart, valStat[iBin]);
    hTotOutSmooth -> SetBinContent(iBin + iStart, valStat[iBin]);
    hStatOut      -> SetBinError(iBin + iStart, errStat[iBin]);
    hSysOut       -> SetBinError(iBin + iStart, errSysR[iBin]);
    hSysOutMix    -> SetBinError(iBin + iStart, errSysM[iBin]);
    hSysOutSmooth -> SetBinError(iBin + iStart, errSysS[iBin]);
    hTotOut       -> SetBinError(iBin + iStart, errTotR[iBin]);
    hTotOutMix    -> SetBinError(iBin + iStart, errTotM[iBin]);
    hTotOutSmooth -> SetBinError(iBin + iStart, errTotS[iBin]);

    // set sum breakdown values
    hPlotStat -> SetBinContent(iBin + iStart, valStat[iBin]);
    hLineStat -> SetBinContent(iBin + iStart, 1.);
    hPlotStat -> SetBinError(iBin + iStart, errStat[iBin]);
    hLineStat -> SetBinError(iBin + iStart, perStat[iBin]);
    for (UInt_t iSys = 0; iSys < NSys; iSys++) {
      perSumSysR[iSys][iBin] = TMath::Sqrt(perSumSys2R[iSys][iBin]);
      perSumSysM[iSys][iBin] = TMath::Sqrt(perSumSys2M[iSys][iBin]);
      perSumSysS[iSys][iBin] = TMath::Sqrt(perSumSys2S[iSys][iBin]);
      perSumTotR[iSys][iBin] = TMath::Sqrt(perSumTot2R[iSys][iBin]);
      perSumTotM[iSys][iBin] = TMath::Sqrt(perSumTot2M[iSys][iBin]);
      perSumTotS[iSys][iBin] = TMath::Sqrt(perSumTot2S[iSys][iBin]);
      errSumSysR[iSys][iBin] = perSumSysR[iSys][iBin] * valStat[iBin];
      errSumSysM[iSys][iBin] = perSumSysM[iSys][iBin] * valStat[iBin];
      errSumSysS[iSys][iBin] = perSumSysS[iSys][iBin] * valStat[iBin];
      errSumTotR[iSys][iBin] = perSumTotR[iSys][iBin] * valStat[iBin];
      errSumTotM[iSys][iBin] = perSumTotM[iSys][iBin] * valStat[iBin];
      errSumTotS[iSys][iBin] = perSumTotS[iSys][iBin] * valStat[iBin];
      hPlotSys[iSys]       -> SetBinContent(iBin + iStart, valStat[iBin]);
      hPlotSysMix[iSys]    -> SetBinContent(iBin + iStart, valStat[iBin]);
      hPlotSysSmooth[iSys] -> SetBinContent(iBin + iStart, valStat[iBin]);
      hPlotTot[iSys]       -> SetBinContent(iBin + iStart, valStat[iBin]);
      hPlotTotMix[iSys]    -> SetBinContent(iBin + iStart, valStat[iBin]);
      hPlotTotSmooth[iSys] -> SetBinContent(iBin + iStart, valStat[iBin]);
      hLineSys[iSys]       -> SetBinContent(iBin + iStart, 1.);
      hLineSysMix[iSys]    -> SetBinContent(iBin + iStart, 1.);
      hLineSysSmooth[iSys] -> SetBinContent(iBin + iStart, 1.);
      hLineTot[iSys]       -> SetBinContent(iBin + iStart, 1.);
      hLineTotMix[iSys]    -> SetBinContent(iBin + iStart, 1.);
      hLineTotSmooth[iSys] -> SetBinContent(iBin + iStart, 1.);
      hPlotSys[iSys]       -> SetBinError(iBin + iStart, errSumSysR[iSys][iBin]);
      hPlotSysMix[iSys]    -> SetBinError(iBin + iStart, errSumSysM[iSys][iBin]);
      hPlotSysSmooth[iSys] -> SetBinError(iBin + iStart, errSumSysS[iSys][iBin]);
      hPlotTot[iSys]       -> SetBinError(iBin + iStart, errSumTotR[iSys][iBin]);
      hPlotTotMix[iSys]    -> SetBinError(iBin + iStart, errSumTotM[iSys][iBin]);
      hPlotTotSmooth[iSys] -> SetBinError(iBin + iStart, errSumTotS[iSys][iBin]);
      hLineSys[iSys]       -> SetBinError(iBin + iStart, perSumSysR[iSys][iBin]);
      hLineSysMix[iSys]    -> SetBinError(iBin + iStart, perSumSysM[iSys][iBin]);
      hLineSysSmooth[iSys] -> SetBinError(iBin + iStart, perSumSysS[iSys][iBin]);
      hLineTot[iSys]       -> SetBinError(iBin + iStart, perSumTotR[iSys][iBin]);
      hLineTotMix[iSys]    -> SetBinError(iBin + iStart, perSumTotM[iSys][iBin]);
      hLineTotSmooth[iSys] -> SetBinError(iBin + iStart, perSumTotS[iSys][iBin]);
    }
    perSumTotR[NSys][iBin] = TMath::Sqrt(perSumTot2R[NSys][iBin]);
    perSumTotS[NSys][iBin] = TMath::Sqrt(perSumTot2S[NSys][iBin]);
    errSumTotR[NSys][iBin] = perSumTotR[NSys][iBin] * valStat[iBin];
    errSumTotS[NSys][iBin] = perSumTotS[NSys][iBin] * valStat[iBin];
    hPlotTot[NSys]       -> SetBinContent(iBin + iStart, valStat[iBin]);
    hPlotTotSmooth[NSys] -> SetBinContent(iBin + iStart, valStat[iBin]);
    hLineTot[NSys]       -> SetBinContent(iBin + iStart, 1.);
    hLineTotSmooth[NSys] -> SetBinContent(iBin + iStart, 1.);
    hPlotTot[NSys]       -> SetBinError(iBin + iStart, errSumTotR[NSys][iBin]);
    hPlotTotSmooth[NSys] -> SetBinError(iBin + iStart, errSumTotS[NSys][iBin]);
    hLineTot[NSys]       -> SetBinError(iBin + iStart, perSumTotR[NSys][iBin]);
    hLineTotSmooth[NSys] -> SetBinError(iBin + iStart, perSumTotS[NSys][iBin]);
  }
  cout << "    Calculated errors." << endl;

  // set styles
  const UInt_t  fMarSt(20);
  const UInt_t  fMarSy(1);
  const UInt_t  fFilSt(0);
  const UInt_t  fFilSy(1001);
  const UInt_t  fLin(1);
  const UInt_t  fWidSt(2);
  const UInt_t  fWidSy(1);
  const UInt_t  fTxt(42);
  const UInt_t  fCnt(1);
  const Float_t fLab(0.03);
  const Float_t fOffX(1.);
  const Float_t fOffY(1.3);
  hPlotStat -> SetMarkerColor(fColStat);
  hPlotStat -> SetMarkerStyle(fMarSt);
  hPlotStat -> SetFillColor(fColStat);
  hPlotStat -> SetFillStyle(fFilSt);
  hPlotStat -> SetLineColor(fColStat);
  hPlotStat -> SetLineStyle(fLin);
  hPlotStat -> SetLineWidth(fWidSt);
  hPlotStat -> SetTitle(sTitle.Data());
  hPlotStat -> SetTitleFont(fTxt);
  hPlotStat -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hPlotStat -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPlotStat -> GetXaxis() -> SetTitleFont(fTxt);
  hPlotStat -> GetXaxis() -> SetTitleOffset(fOffX);
  hPlotStat -> GetXaxis() -> SetLabelFont(fTxt);
  hPlotStat -> GetXaxis() -> SetLabelSize(fLab);
  hPlotStat -> GetXaxis() -> CenterTitle(fCnt);
  hPlotStat -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPlotStat -> GetYaxis() -> SetTitleFont(fTxt);
  hPlotStat -> GetYaxis() -> SetTitleOffset(fOffY);
  hPlotStat -> GetYaxis() -> SetLabelFont(fTxt);
  hPlotStat -> GetYaxis() -> SetLabelSize(fLab);
  hPlotStat -> GetYaxis() -> CenterTitle(fCnt);
  hLineStat -> SetMarkerColor(fColStat);
  hLineStat -> SetMarkerStyle(fMarSt);
  hLineStat -> SetFillColor(fColStat);
  hLineStat -> SetFillStyle(fFilSt);
  hLineStat -> SetLineColor(fColStat);
  hLineStat -> SetLineStyle(fLin);
  hLineStat -> SetLineWidth(fWidSt);
  hLineStat -> SetTitle(sTitle.Data());
  hLineStat -> SetTitleFont(fTxt);
  hLineStat -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hLineStat -> GetXaxis() -> SetTitle(sTitleX.Data());
  hLineStat -> GetXaxis() -> SetTitleFont(fTxt);
  hLineStat -> GetXaxis() -> SetTitleOffset(fOffX);
  hLineStat -> GetXaxis() -> SetLabelFont(fTxt);
  hLineStat -> GetXaxis() -> SetLabelSize(fLab);
  hLineStat -> GetXaxis() -> CenterTitle(fCnt);
  hLineStat -> GetYaxis() -> SetTitle(sTitleL.Data());
  hLineStat -> GetYaxis() -> SetTitleFont(fTxt);
  hLineStat -> GetYaxis() -> SetTitleOffset(fOffY);
  hLineStat -> GetYaxis() -> SetLabelFont(fTxt);
  hLineStat -> GetYaxis() -> SetLabelSize(fLab);
  hLineStat -> GetYaxis() -> CenterTitle(fCnt);
  for (UInt_t iSys = 0; iSys <NSys; iSys++) {
    hPlotSys[iSys]       -> SetMarkerColor(fColSys[iSys]);
    hPlotSys[iSys]       -> SetMarkerStyle(fMarSy);
    hPlotSys[iSys]       -> SetFillColor(fColSys[iSys]);
    hPlotSys[iSys]       -> SetFillStyle(fFilSy);
    hPlotSys[iSys]       -> SetLineColor(fColSys[iSys]);
    hPlotSys[iSys]       -> SetLineStyle(fLin);
    hPlotSys[iSys]       -> SetLineWidth(fWidSy);
    hPlotSys[iSys]       -> SetTitle(sTitle.Data());
    hPlotSys[iSys]       -> SetTitleFont(fTxt);
    hPlotSys[iSys]       -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPlotSys[iSys]       -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPlotSys[iSys]       -> GetXaxis() -> SetTitleFont(fTxt);
    hPlotSys[iSys]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hPlotSys[iSys]       -> GetXaxis() -> SetLabelFont(fTxt);
    hPlotSys[iSys]       -> GetXaxis() -> SetLabelSize(fLab);
    hPlotSys[iSys]       -> GetXaxis() -> CenterTitle(fCnt);
    hPlotSys[iSys]       -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPlotSys[iSys]       -> GetYaxis() -> SetTitleFont(fTxt);
    hPlotSys[iSys]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hPlotSys[iSys]       -> GetYaxis() -> SetLabelFont(fTxt);
    hPlotSys[iSys]       -> GetYaxis() -> SetLabelSize(fLab);
    hPlotSys[iSys]       -> GetYaxis() -> CenterTitle(fCnt);
    hPlotSysMix[iSys]    -> SetMarkerColor(fColSys[iSys]);
    hPlotSysMix[iSys]    -> SetMarkerStyle(fMarSy);
    hPlotSysMix[iSys]    -> SetFillColor(fColSys[iSys]);
    hPlotSysMix[iSys]    -> SetFillStyle(fFilSy);
    hPlotSysMix[iSys]    -> SetLineColor(fColSys[iSys]);
    hPlotSysMix[iSys]    -> SetLineStyle(fLin);
    hPlotSysMix[iSys]    -> SetLineWidth(fWidSy);
    hPlotSysMix[iSys]    -> SetTitle(sTitle.Data());
    hPlotSysMix[iSys]    -> SetTitleFont(fTxt);
    hPlotSysMix[iSys]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPlotSysMix[iSys]    -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPlotSysMix[iSys]    -> GetXaxis() -> SetTitleFont(fTxt);
    hPlotSysMix[iSys]    -> GetXaxis() -> SetTitleOffset(fOffX);
    hPlotSysMix[iSys]    -> GetXaxis() -> SetLabelFont(fTxt);
    hPlotSysMix[iSys]    -> GetXaxis() -> SetLabelSize(fLab);
    hPlotSysMix[iSys]    -> GetXaxis() -> CenterTitle(fCnt);
    hPlotSysMix[iSys]    -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPlotSysMix[iSys]    -> GetYaxis() -> SetTitleFont(fTxt);
    hPlotSysMix[iSys]    -> GetYaxis() -> SetTitleOffset(fOffY);
    hPlotSysMix[iSys]    -> GetYaxis() -> SetLabelFont(fTxt);
    hPlotSysMix[iSys]    -> GetYaxis() -> SetLabelSize(fLab);
    hPlotSysMix[iSys]    -> GetYaxis() -> CenterTitle(fCnt);
    hPlotSysSmooth[iSys] -> SetMarkerColor(fColSys[iSys]);
    hPlotSysSmooth[iSys] -> SetMarkerStyle(fMarSy);
    hPlotSysSmooth[iSys] -> SetFillColor(fColSys[iSys]);
    hPlotSysSmooth[iSys] -> SetFillStyle(fFilSy);
    hPlotSysSmooth[iSys] -> SetLineColor(fColSys[iSys]);
    hPlotSysSmooth[iSys] -> SetLineStyle(fLin);
    hPlotSysSmooth[iSys] -> SetLineWidth(fWidSy);
    hPlotSysSmooth[iSys] -> SetTitle(sTitle.Data());
    hPlotSysSmooth[iSys] -> SetTitleFont(fTxt);
    hPlotSysSmooth[iSys] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPlotSysSmooth[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPlotSysSmooth[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hPlotSysSmooth[iSys] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPlotSysSmooth[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hPlotSysSmooth[iSys] -> GetXaxis() -> SetLabelSize(fLab);
    hPlotSysSmooth[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hPlotSysSmooth[iSys] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPlotSysSmooth[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hPlotSysSmooth[iSys] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPlotSysSmooth[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hPlotSysSmooth[iSys] -> GetYaxis() -> SetLabelSize(fLab);
    hPlotSysSmooth[iSys] -> GetYaxis() -> CenterTitle(fCnt);
    hPlotTot[iSys]       -> SetMarkerColor(fColSys[iSys]);
    hPlotTot[iSys]       -> SetMarkerStyle(fMarSy);
    hPlotTot[iSys]       -> SetFillColor(fColSys[iSys]);
    hPlotTot[iSys]       -> SetFillStyle(fFilSy);
    hPlotTot[iSys]       -> SetLineColor(fColSys[iSys]);
    hPlotTot[iSys]       -> SetLineStyle(fLin);
    hPlotTot[iSys]       -> SetLineWidth(fWidSy);
    hPlotTot[iSys]       -> SetTitle(sTitle.Data());
    hPlotTot[iSys]       -> SetTitleFont(fTxt);
    hPlotTot[iSys]       -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPlotTot[iSys]       -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPlotTot[iSys]       -> GetXaxis() -> SetTitleFont(fTxt);
    hPlotTot[iSys]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hPlotTot[iSys]       -> GetXaxis() -> SetLabelFont(fTxt);
    hPlotTot[iSys]       -> GetXaxis() -> SetLabelSize(fLab);
    hPlotTot[iSys]       -> GetXaxis() -> CenterTitle(fCnt);
    hPlotTot[iSys]       -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPlotTot[iSys]       -> GetYaxis() -> SetTitleFont(fTxt);
    hPlotTot[iSys]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hPlotTot[iSys]       -> GetYaxis() -> SetLabelFont(fTxt);
    hPlotTot[iSys]       -> GetYaxis() -> SetLabelSize(fLab);
    hPlotTot[iSys]       -> GetYaxis() -> CenterTitle(fCnt);
    hPlotTotMix[iSys]    -> SetMarkerColor(fColSys[iSys]);
    hPlotTotMix[iSys]    -> SetMarkerStyle(fMarSy);
    hPlotTotMix[iSys]    -> SetFillColor(fColSys[iSys]);
    hPlotTotMix[iSys]    -> SetFillStyle(fFilSy);
    hPlotTotMix[iSys]    -> SetLineColor(fColSys[iSys]);
    hPlotTotMix[iSys]    -> SetLineStyle(fLin);
    hPlotTotMix[iSys]    -> SetLineWidth(fWidSy);
    hPlotTotMix[iSys]    -> SetTitle(sTitle.Data());
    hPlotTotMix[iSys]    -> SetTitleFont(fTxt);
    hPlotTotMix[iSys]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPlotTotMix[iSys]    -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPlotTotMix[iSys]    -> GetXaxis() -> SetTitleFont(fTxt);
    hPlotTotMix[iSys]    -> GetXaxis() -> SetTitleOffset(fOffX);
    hPlotTotMix[iSys]    -> GetXaxis() -> SetLabelFont(fTxt);
    hPlotTotMix[iSys]    -> GetXaxis() -> SetLabelSize(fLab);
    hPlotTotMix[iSys]    -> GetXaxis() -> CenterTitle(fCnt);
    hPlotTotMix[iSys]    -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPlotTotMix[iSys]    -> GetYaxis() -> SetTitleFont(fTxt);
    hPlotTotMix[iSys]    -> GetYaxis() -> SetTitleOffset(fOffY);
    hPlotTotMix[iSys]    -> GetYaxis() -> SetLabelFont(fTxt);
    hPlotTotMix[iSys]    -> GetYaxis() -> SetLabelSize(fLab);
    hPlotTotMix[iSys]    -> GetYaxis() -> CenterTitle(fCnt);
    hPlotTotSmooth[iSys] -> SetMarkerColor(fColSys[iSys]);
    hPlotTotSmooth[iSys] -> SetMarkerStyle(fMarSy);
    hPlotTotSmooth[iSys] -> SetFillColor(fColSys[iSys]);
    hPlotTotSmooth[iSys] -> SetFillStyle(fFilSy);
    hPlotTotSmooth[iSys] -> SetLineColor(fColSys[iSys]);
    hPlotTotSmooth[iSys] -> SetLineStyle(fLin);
    hPlotTotSmooth[iSys] -> SetLineWidth(fWidSy);
    hPlotTotSmooth[iSys] -> SetTitle(sTitle.Data());
    hPlotTotSmooth[iSys] -> SetTitleFont(fTxt);
    hPlotTotSmooth[iSys] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hPlotTotSmooth[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hPlotTotSmooth[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hPlotTotSmooth[iSys] -> GetXaxis() -> SetTitleOffset(fOffX);
    hPlotTotSmooth[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hPlotTotSmooth[iSys] -> GetXaxis() -> SetLabelSize(fLab);
    hPlotTotSmooth[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hPlotTotSmooth[iSys] -> GetYaxis() -> SetTitle(sTitleY.Data());
    hPlotTotSmooth[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hPlotTotSmooth[iSys] -> GetYaxis() -> SetTitleOffset(fOffY);
    hPlotTotSmooth[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hPlotTotSmooth[iSys] -> GetYaxis() -> SetLabelSize(fLab);
    hPlotTotSmooth[iSys] -> GetYaxis() -> CenterTitle(fCnt);
    hLineSys[iSys]       -> SetMarkerColor(fColSys[iSys]);
    hLineSys[iSys]       -> SetMarkerStyle(fMarSy);
    hLineSys[iSys]       -> SetFillColor(fColSys[iSys]);
    hLineSys[iSys]       -> SetFillStyle(fFilSy);
    hLineSys[iSys]       -> SetLineColor(fColSys[iSys]);
    hLineSys[iSys]       -> SetLineStyle(fLin);
    hLineSys[iSys]       -> SetLineWidth(fWidSy);
    hLineSys[iSys]       -> SetTitle(sTitle.Data());
    hLineSys[iSys]       -> SetTitleFont(fTxt);
    hLineSys[iSys]       -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hLineSys[iSys]       -> GetXaxis() -> SetTitle(sTitleX.Data());
    hLineSys[iSys]       -> GetXaxis() -> SetTitleFont(fTxt);
    hLineSys[iSys]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hLineSys[iSys]       -> GetXaxis() -> SetLabelFont(fTxt);
    hLineSys[iSys]       -> GetXaxis() -> SetLabelSize(fLab);
    hLineSys[iSys]       -> GetXaxis() -> CenterTitle(fCnt);
    hLineSys[iSys]       -> GetYaxis() -> SetTitle(sTitleL.Data());
    hLineSys[iSys]       -> GetYaxis() -> SetTitleFont(fTxt);
    hLineSys[iSys]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hLineSys[iSys]       -> GetYaxis() -> SetLabelFont(fTxt);
    hLineSys[iSys]       -> GetYaxis() -> SetLabelSize(fLab);
    hLineSys[iSys]       -> GetYaxis() -> CenterTitle(fCnt);
    hLineSysMix[iSys]    -> SetMarkerColor(fColSys[iSys]);
    hLineSysMix[iSys]    -> SetMarkerStyle(fMarSy);
    hLineSysMix[iSys]    -> SetFillColor(fColSys[iSys]);
    hLineSysMix[iSys]    -> SetFillStyle(fFilSy);
    hLineSysMix[iSys]    -> SetLineColor(fColSys[iSys]);
    hLineSysMix[iSys]    -> SetLineStyle(fLin);
    hLineSysMix[iSys]    -> SetLineWidth(fWidSy);
    hLineSysMix[iSys]    -> SetTitle(sTitle.Data());
    hLineSysMix[iSys]    -> SetTitleFont(fTxt);
    hLineSysMix[iSys]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hLineSysMix[iSys]    -> GetXaxis() -> SetTitle(sTitleX.Data());
    hLineSysMix[iSys]    -> GetXaxis() -> SetTitleFont(fTxt);
    hLineSysMix[iSys]    -> GetXaxis() -> SetTitleOffset(fOffX);
    hLineSysMix[iSys]    -> GetXaxis() -> SetLabelFont(fTxt);
    hLineSysMix[iSys]    -> GetXaxis() -> SetLabelSize(fLab);
    hLineSysMix[iSys]    -> GetXaxis() -> CenterTitle(fCnt);
    hLineSysMix[iSys]    -> GetYaxis() -> SetTitle(sTitleL.Data());
    hLineSysMix[iSys]    -> GetYaxis() -> SetTitleFont(fTxt);
    hLineSysMix[iSys]    -> GetYaxis() -> SetTitleOffset(fOffY);
    hLineSysMix[iSys]    -> GetYaxis() -> SetLabelFont(fTxt);
    hLineSysMix[iSys]    -> GetYaxis() -> SetLabelSize(fLab);
    hLineSysMix[iSys]    -> GetYaxis() -> CenterTitle(fCnt);
    hLineSysSmooth[iSys] -> SetMarkerColor(fColSys[iSys]);
    hLineSysSmooth[iSys] -> SetMarkerStyle(fMarSy);
    hLineSysSmooth[iSys] -> SetFillColor(fColSys[iSys]);
    hLineSysSmooth[iSys] -> SetFillStyle(fFilSy);
    hLineSysSmooth[iSys] -> SetLineColor(fColSys[iSys]);
    hLineSysSmooth[iSys] -> SetLineStyle(fLin);
    hLineSysSmooth[iSys] -> SetLineWidth(fWidSy);
    hLineSysSmooth[iSys] -> SetTitle(sTitle.Data());
    hLineSysSmooth[iSys] -> SetTitleFont(fTxt);
    hLineSysSmooth[iSys] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hLineSysSmooth[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hLineSysSmooth[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hLineSysSmooth[iSys] -> GetXaxis() -> SetTitleOffset(fOffX);
    hLineSysSmooth[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hLineSysSmooth[iSys] -> GetXaxis() -> SetLabelSize(fLab);
    hLineSysSmooth[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hLineSysSmooth[iSys] -> GetYaxis() -> SetTitle(sTitleL.Data());
    hLineSysSmooth[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hLineSysSmooth[iSys] -> GetYaxis() -> SetTitleOffset(fOffY);
    hLineSysSmooth[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hLineSysSmooth[iSys] -> GetYaxis() -> SetLabelSize(fLab);
    hLineSysSmooth[iSys] -> GetYaxis() -> CenterTitle(fCnt);
    hLineTot[iSys]       -> SetMarkerColor(fColSys[iSys]);
    hLineTot[iSys]       -> SetMarkerStyle(fMarSy);
    hLineTot[iSys]       -> SetFillColor(fColSys[iSys]);
    hLineTot[iSys]       -> SetFillStyle(fFilSy);
    hLineTot[iSys]       -> SetLineColor(fColSys[iSys]);
    hLineTot[iSys]       -> SetLineStyle(fLin);
    hLineTot[iSys]       -> SetLineWidth(fWidSy);
    hLineTot[iSys]       -> SetTitle(sTitle.Data());
    hLineTot[iSys]       -> SetTitleFont(fTxt);
    hLineTot[iSys]       -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hLineTot[iSys]       -> GetXaxis() -> SetTitle(sTitleX.Data());
    hLineTot[iSys]       -> GetXaxis() -> SetTitleFont(fTxt);
    hLineTot[iSys]       -> GetXaxis() -> SetTitleOffset(fOffX);
    hLineTot[iSys]       -> GetXaxis() -> SetLabelFont(fTxt);
    hLineTot[iSys]       -> GetXaxis() -> SetLabelSize(fLab);
    hLineTot[iSys]       -> GetXaxis() -> CenterTitle(fCnt);
    hLineTot[iSys]       -> GetYaxis() -> SetTitle(sTitleL.Data());
    hLineTot[iSys]       -> GetYaxis() -> SetTitleFont(fTxt);
    hLineTot[iSys]       -> GetYaxis() -> SetTitleOffset(fOffY);
    hLineTot[iSys]       -> GetYaxis() -> SetLabelFont(fTxt);
    hLineTot[iSys]       -> GetYaxis() -> SetLabelSize(fLab);
    hLineTot[iSys]       -> GetYaxis() -> CenterTitle(fCnt);
    hLineTotMix[iSys]    -> SetMarkerColor(fColSys[iSys]);
    hLineTotMix[iSys]    -> SetMarkerStyle(fMarSy);
    hLineTotMix[iSys]    -> SetFillColor(fColSys[iSys]);
    hLineTotMix[iSys]    -> SetFillStyle(fFilSy);
    hLineTotMix[iSys]    -> SetLineColor(fColSys[iSys]);
    hLineTotMix[iSys]    -> SetLineStyle(fLin);
    hLineTotMix[iSys]    -> SetLineWidth(fWidSy);
    hLineTotMix[iSys]    -> SetTitle(sTitle.Data());
    hLineTotMix[iSys]    -> SetTitleFont(fTxt);
    hLineTotMix[iSys]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hLineTotMix[iSys]    -> GetXaxis() -> SetTitle(sTitleX.Data());
    hLineTotMix[iSys]    -> GetXaxis() -> SetTitleFont(fTxt);
    hLineTotMix[iSys]    -> GetXaxis() -> SetTitleOffset(fOffX);
    hLineTotMix[iSys]    -> GetXaxis() -> SetLabelFont(fTxt);
    hLineTotMix[iSys]    -> GetXaxis() -> SetLabelSize(fLab);
    hLineTotMix[iSys]    -> GetXaxis() -> CenterTitle(fCnt);
    hLineTotMix[iSys]    -> GetYaxis() -> SetTitle(sTitleL.Data());
    hLineTotMix[iSys]    -> GetYaxis() -> SetTitleFont(fTxt);
    hLineTotMix[iSys]    -> GetYaxis() -> SetTitleOffset(fOffY);
    hLineTotMix[iSys]    -> GetYaxis() -> SetLabelFont(fTxt);
    hLineTotMix[iSys]    -> GetYaxis() -> SetLabelSize(fLab);
    hLineTotMix[iSys]    -> GetYaxis() -> CenterTitle(fCnt);
    hLineTotSmooth[iSys] -> SetMarkerColor(fColSys[iSys]);
    hLineTotSmooth[iSys] -> SetMarkerStyle(fMarSy);
    hLineTotSmooth[iSys] -> SetFillColor(fColSys[iSys]);
    hLineTotSmooth[iSys] -> SetFillStyle(fFilSy);
    hLineTotSmooth[iSys] -> SetLineColor(fColSys[iSys]);
    hLineTotSmooth[iSys] -> SetLineStyle(fLin);
    hLineTotSmooth[iSys] -> SetLineWidth(fWidSy);
    hLineTotSmooth[iSys] -> SetTitle(sTitle.Data());
    hLineTotSmooth[iSys] -> SetTitleFont(fTxt);
    hLineTotSmooth[iSys] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
    hLineTotSmooth[iSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
    hLineTotSmooth[iSys] -> GetXaxis() -> SetTitleFont(fTxt);
    hLineTotSmooth[iSys] -> GetXaxis() -> SetTitleOffset(fOffX);
    hLineTotSmooth[iSys] -> GetXaxis() -> SetLabelFont(fTxt);
    hLineTotSmooth[iSys] -> GetXaxis() -> SetLabelSize(fLab);
    hLineTotSmooth[iSys] -> GetXaxis() -> CenterTitle(fCnt);
    hLineTotSmooth[iSys] -> GetYaxis() -> SetTitle(sTitleL.Data());
    hLineTotSmooth[iSys] -> GetYaxis() -> SetTitleFont(fTxt);
    hLineTotSmooth[iSys] -> GetYaxis() -> SetTitleOffset(fOffY);
    hLineTotSmooth[iSys] -> GetYaxis() -> SetLabelFont(fTxt);
    hLineTotSmooth[iSys] -> GetYaxis() -> SetLabelSize(fLab);
    hLineTotSmooth[iSys] -> GetYaxis() -> CenterTitle(fCnt);
  }
  hPlotTot[NSys]       -> SetMarkerColor(fColStat);
  hPlotTot[NSys]       -> SetMarkerStyle(fMarSt);
  hPlotTot[NSys]       -> SetFillColor(fColStat);
  hPlotTot[NSys]       -> SetFillStyle(fFilSt);
  hPlotTot[NSys]       -> SetLineColor(fColStat);
  hPlotTot[NSys]       -> SetLineStyle(fLin);
  hPlotTot[NSys]       -> SetLineWidth(fWidSt);
  hPlotTot[NSys]       -> SetTitle(sTitle.Data());
  hPlotTot[NSys]       -> SetTitleFont(fTxt);
  hPlotTot[NSys]       -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hPlotTot[NSys]       -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPlotTot[NSys]       -> GetXaxis() -> SetTitleFont(fTxt);
  hPlotTot[NSys]       -> GetXaxis() -> SetTitleOffset(fOffX);
  hPlotTot[NSys]       -> GetXaxis() -> SetLabelFont(fTxt);
  hPlotTot[NSys]       -> GetXaxis() -> SetLabelSize(fLab);
  hPlotTot[NSys]       -> GetXaxis() -> CenterTitle(fCnt);
  hPlotTot[NSys]       -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPlotTot[NSys]       -> GetYaxis() -> SetTitleFont(fTxt);
  hPlotTot[NSys]       -> GetYaxis() -> SetTitleOffset(fOffY);
  hPlotTot[NSys]       -> GetYaxis() -> SetLabelFont(fTxt);
  hPlotTot[NSys]       -> GetYaxis() -> SetLabelSize(fLab);
  hPlotTot[NSys]       -> GetYaxis() -> CenterTitle(fCnt);
  hPlotTotMix[NSys]    -> SetMarkerColor(fColStat);
  hPlotTotMix[NSys]    -> SetMarkerStyle(fMarSt);
  hPlotTotMix[NSys]    -> SetFillColor(fColStat);
  hPlotTotMix[NSys]    -> SetFillStyle(fFilSt);
  hPlotTotMix[NSys]    -> SetLineColor(fColStat);
  hPlotTotMix[NSys]    -> SetLineStyle(fLin);
  hPlotTotMix[NSys]    -> SetLineWidth(fWidSt);
  hPlotTotMix[NSys]    -> SetTitle(sTitle.Data());
  hPlotTotMix[NSys]    -> SetTitleFont(fTxt);
  hPlotTotMix[NSys]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hPlotTotMix[NSys]    -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPlotTotMix[NSys]    -> GetXaxis() -> SetTitleFont(fTxt);
  hPlotTotMix[NSys]    -> GetXaxis() -> SetTitleOffset(fOffX);
  hPlotTotMix[NSys]    -> GetXaxis() -> SetLabelFont(fTxt);
  hPlotTotMix[NSys]    -> GetXaxis() -> SetLabelSize(fLab);
  hPlotTotMix[NSys]    -> GetXaxis() -> CenterTitle(fCnt);
  hPlotTotMix[NSys]    -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPlotTotMix[NSys]    -> GetYaxis() -> SetTitleFont(fTxt);
  hPlotTotMix[NSys]    -> GetYaxis() -> SetTitleOffset(fOffY);
  hPlotTotMix[NSys]    -> GetYaxis() -> SetLabelFont(fTxt);
  hPlotTotMix[NSys]    -> GetYaxis() -> SetLabelSize(fLab);
  hPlotTotMix[NSys]    -> GetYaxis() -> CenterTitle(fCnt);
  hPlotTotSmooth[NSys] -> SetMarkerColor(fColStat);
  hPlotTotSmooth[NSys] -> SetMarkerStyle(fMarSt);
  hPlotTotSmooth[NSys] -> SetFillColor(fColStat);
  hPlotTotSmooth[NSys] -> SetFillStyle(fFilSt);
  hPlotTotSmooth[NSys] -> SetLineColor(fColStat);
  hPlotTotSmooth[NSys] -> SetLineStyle(fLin);
  hPlotTotSmooth[NSys] -> SetLineWidth(fWidSt);
  hPlotTotSmooth[NSys] -> SetTitle(sTitle.Data());
  hPlotTotSmooth[NSys] -> SetTitleFont(fTxt);
  hPlotTotSmooth[NSys] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hPlotTotSmooth[NSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hPlotTotSmooth[NSys] -> GetXaxis() -> SetTitleFont(fTxt);
  hPlotTotSmooth[NSys] -> GetXaxis() -> SetTitleOffset(fOffX);
  hPlotTotSmooth[NSys] -> GetXaxis() -> SetLabelFont(fTxt);
  hPlotTotSmooth[NSys] -> GetXaxis() -> SetLabelSize(fLab);
  hPlotTotSmooth[NSys] -> GetXaxis() -> CenterTitle(fCnt);
  hPlotTotSmooth[NSys] -> GetYaxis() -> SetTitle(sTitleY.Data());
  hPlotTotSmooth[NSys] -> GetYaxis() -> SetTitleFont(fTxt);
  hPlotTotSmooth[NSys] -> GetYaxis() -> SetTitleOffset(fOffY);
  hPlotTotSmooth[NSys] -> GetYaxis() -> SetLabelFont(fTxt);
  hPlotTotSmooth[NSys] -> GetYaxis() -> SetLabelSize(fLab);
  hPlotTotSmooth[NSys] -> GetYaxis() -> CenterTitle(fCnt);
  hLineTot[NSys]       -> SetMarkerColor(fColStat);
  hLineTot[NSys]       -> SetMarkerStyle(fMarSt);
  hLineTot[NSys]       -> SetFillColor(fColStat);
  hLineTot[NSys]       -> SetFillStyle(fFilSt);
  hLineTot[NSys]       -> SetLineColor(fColStat);
  hLineTot[NSys]       -> SetLineStyle(fLin);
  hLineTot[NSys]       -> SetLineWidth(fWidSt);
  hLineTot[NSys]       -> SetTitle(sTitle.Data());
  hLineTot[NSys]       -> SetTitleFont(fTxt);
  hLineTot[NSys]       -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hLineTot[NSys]       -> GetXaxis() -> SetTitle(sTitleX.Data());
  hLineTot[NSys]       -> GetXaxis() -> SetTitleFont(fTxt);
  hLineTot[NSys]       -> GetXaxis() -> SetTitleOffset(fOffX);
  hLineTot[NSys]       -> GetXaxis() -> SetLabelFont(fTxt);
  hLineTot[NSys]       -> GetXaxis() -> SetLabelSize(fLab);
  hLineTot[NSys]       -> GetXaxis() -> CenterTitle(fCnt);
  hLineTot[NSys]       -> GetYaxis() -> SetTitle(sTitleL.Data());
  hLineTot[NSys]       -> GetYaxis() -> SetTitleFont(fTxt);
  hLineTot[NSys]       -> GetYaxis() -> SetTitleOffset(fOffY);
  hLineTot[NSys]       -> GetYaxis() -> SetLabelFont(fTxt);
  hLineTot[NSys]       -> GetYaxis() -> SetLabelSize(fLab);
  hLineTot[NSys]       -> GetYaxis() -> CenterTitle(fCnt);
  hLineTotMix[NSys]    -> SetMarkerColor(fColStat);
  hLineTotMix[NSys]    -> SetMarkerStyle(fMarSt);
  hLineTotMix[NSys]    -> SetFillColor(fColStat);
  hLineTotMix[NSys]    -> SetFillStyle(fFilSt);
  hLineTotMix[NSys]    -> SetLineColor(fColStat);
  hLineTotMix[NSys]    -> SetLineStyle(fLin);
  hLineTotMix[NSys]    -> SetLineWidth(fWidSt);
  hLineTotMix[NSys]    -> SetTitle(sTitle.Data());
  hLineTotMix[NSys]    -> SetTitleFont(fTxt);
  hLineTotMix[NSys]    -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hLineTotMix[NSys]    -> GetXaxis() -> SetTitle(sTitleX.Data());
  hLineTotMix[NSys]    -> GetXaxis() -> SetTitleFont(fTxt);
  hLineTotMix[NSys]    -> GetXaxis() -> SetTitleOffset(fOffX);
  hLineTotMix[NSys]    -> GetXaxis() -> SetLabelFont(fTxt);
  hLineTotMix[NSys]    -> GetXaxis() -> SetLabelSize(fLab);
  hLineTotMix[NSys]    -> GetXaxis() -> CenterTitle(fCnt);
  hLineTotMix[NSys]    -> GetYaxis() -> SetTitle(sTitleL.Data());
  hLineTotMix[NSys]    -> GetYaxis() -> SetTitleFont(fTxt);
  hLineTotMix[NSys]    -> GetYaxis() -> SetTitleOffset(fOffY);
  hLineTotMix[NSys]    -> GetYaxis() -> SetLabelFont(fTxt);
  hLineTotMix[NSys]    -> GetYaxis() -> SetLabelSize(fLab);
  hLineTotMix[NSys]    -> GetYaxis() -> CenterTitle(fCnt);
  hLineTotSmooth[NSys] -> SetMarkerColor(fColStat);
  hLineTotSmooth[NSys] -> SetMarkerStyle(fMarSt);
  hLineTotSmooth[NSys] -> SetFillColor(fColStat);
  hLineTotSmooth[NSys] -> SetFillStyle(fFilSt);
  hLineTotSmooth[NSys] -> SetLineColor(fColStat);
  hLineTotSmooth[NSys] -> SetLineStyle(fLin);
  hLineTotSmooth[NSys] -> SetLineWidth(fWidSt);
  hLineTotSmooth[NSys] -> SetTitle(sTitle.Data());
  hLineTotSmooth[NSys] -> SetTitleFont(fTxt);
  hLineTotSmooth[NSys] -> GetXaxis() -> SetRangeUser(xPlotRange[0], xPlotRange[1]);
  hLineTotSmooth[NSys] -> GetXaxis() -> SetTitle(sTitleX.Data());
  hLineTotSmooth[NSys] -> GetXaxis() -> SetTitleFont(fTxt);
  hLineTotSmooth[NSys] -> GetXaxis() -> SetTitleOffset(fOffX);
  hLineTotSmooth[NSys] -> GetXaxis() -> SetLabelFont(fTxt);
  hLineTotSmooth[NSys] -> GetXaxis() -> SetLabelSize(fLab);
  hLineTotSmooth[NSys] -> GetXaxis() -> CenterTitle(fCnt);
  hLineTotSmooth[NSys] -> GetYaxis() -> SetTitle(sTitleL.Data());
  hLineTotSmooth[NSys] -> GetYaxis() -> SetTitleFont(fTxt);
  hLineTotSmooth[NSys] -> GetYaxis() -> SetTitleOffset(fOffY);
  hLineTotSmooth[NSys] -> GetYaxis() -> SetLabelFont(fTxt);
  hLineTotSmooth[NSys] -> GetYaxis() -> SetLabelSize(fLab);
  hLineTotSmooth[NSys] -> GetYaxis() -> CenterTitle(fCnt);
  cout << "    Set styles." << endl;

  // make legends
  const UInt_t  fColLe(0);
  const UInt_t  fFilLe(0);
  const UInt_t  fLinLe(0);
  const UInt_t  fAlnLe(12);
  const TString sRawErr("Raw Systematic Uncertainties");
  const TString sMixErr("Mixed raw and smoothed sys. uncertainties");
  const TString sSmoothErr("Smoothed Systematic Uncertainties");
  const Float_t fLegXY[NVtx] = {0.7, 0.7, 0.9, 0.9};
  TLegend *legR = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sRawErr.Data());
  legR -> SetFillColor(fColLe);
  legR -> SetFillStyle(fFilLe);
  legR -> SetLineColor(fColLe);
  legR -> SetLineStyle(fLinLe);
  legR -> SetTextFont(fTxt);
  legR -> SetTextAlign(fAlnLe);
  legR -> AddEntry(hPlotStat, sLabelStat.Data(), "pf");
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    legR -> AddEntry(hPlotSys[iSys], sLabelSys[iSys].Data(), "f");
  }

  TLegend *legM = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sMixErr.Data());
  legM -> SetFillColor(fColLe);
  legM -> SetFillStyle(fFilLe);
  legM -> SetLineColor(fColLe);
  legM -> SetLineStyle(fLinLe);
  legM -> SetTextFont(fTxt);
  legM -> SetTextAlign(fAlnLe);
  legM -> AddEntry(hPlotStat, sLabelStat.Data(), "pf");
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    legM -> AddEntry(hPlotSysMix[iSys], sLabelSys[iSys].Data(), "f");
  }

  TLegend *legS = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sSmoothErr.Data());
  legS -> SetFillColor(fColLe);
  legS -> SetFillStyle(fFilLe);
  legS -> SetLineColor(fColLe);
  legS -> SetLineStyle(fLinLe);
  legS -> SetTextFont(fTxt);
  legS -> SetTextAlign(fAlnLe);
  legS -> AddEntry(hPlotStat, sLabelStat.Data(), "pf");
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    legS -> AddEntry(hPlotSysSmooth[iSys], sLabelSys[iSys].Data(), "f");
  }

  TLegend *sysR = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sRawErr.Data());
  sysR -> SetFillColor(fColLe);
  sysR -> SetFillStyle(fFilLe);
  sysR -> SetLineColor(fColLe);
  sysR -> SetLineStyle(fLinLe);
  sysR -> SetTextFont(fTxt);
  sysR -> SetTextAlign(fAlnLe);
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    sysR -> AddEntry(hPlotSys[iSys], sLabelSys[iSys].Data(), "f");
  }

  TLegend *sysM = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sMixErr.Data());
  sysM -> SetFillColor(fColLe);
  sysM -> SetFillStyle(fFilLe);
  sysM -> SetLineColor(fColLe);
  sysM -> SetLineStyle(fLinLe);
  sysM -> SetTextFont(fTxt);
  sysM -> SetTextAlign(fAlnLe);
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    sysM -> AddEntry(hPlotSysMix[iSys], sLabelSys[iSys].Data(), "f");
  }

  TLegend *sysS = new TLegend(fLegXY[0], fLegXY[1], fLegXY[2], fLegXY[3], sSmoothErr.Data());
  sysS -> SetFillColor(fColLe);
  sysS -> SetFillStyle(fFilLe);
  sysS -> SetLineColor(fColLe);
  sysS -> SetLineStyle(fLinLe);
  sysS -> SetTextFont(fTxt);
  sysS -> SetTextAlign(fAlnLe);
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    sysS -> AddEntry(hPlotSysSmooth[iSys], sLabelSys[iSys].Data(), "f");
  }
  cout << "    Made legends." << endl;

  // make text
  const UInt_t  fColTx(0);
  const UInt_t  fFilTx(0);
  const UInt_t  fLinTx(0);
  const UInt_t  fAlnTx(32);
  const UInt_t  fAlnLb(22);
  const UInt_t  fFonLb(62);
  const Float_t fTxtTx[NVtx]        = {0.5, 0.5, 0.7, 0.7};
  const Float_t fTxtLb[NHist][NVtx] = {{0.1, 0.2, 0.3, 0.3}, {0.1, 0.4, 0.3, 0.5}, {0.1, 0.6, 0.3, 0.7}};

  TPaveText *txt = new TPaveText(fTxtTx[0], fTxtTx[1], fTxtTx[2], fTxtTx[3], "NDC NB");
  txt -> SetFillColor(fColTx);
  txt -> SetFillStyle(fFilTx);
  txt -> SetLineColor(fColTx);
  txt -> SetLineStyle(fLinTx);
  txt -> SetTextFont(fTxt);
  txt -> SetTextAlign(fAlnTx);
  txt -> AddText(sCol.Data());
  txt -> AddText(sTrg.Data());
  txt -> AddText(sJet.Data());
  txt -> AddText(sTyp.Data());

  TPaveText *lbl[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    lbl[iHist] = new TPaveText(fTxtLb[iHist][0], fTxtLb[iHist][1], fTxtLb[iHist][2], fTxtLb[iHist][3], "NDC NB");
    lbl[iHist] -> SetFillColor(fColTx);
    lbl[iHist] -> SetFillStyle(fFilTx);
    lbl[iHist] -> SetLineColor(fColTx);
    lbl[iHist] -> SetLineStyle(fLinTx);
    lbl[iHist] -> SetTextFont(fFonLb);
    lbl[iHist] -> SetTextAlign(fAlnLb);
    lbl[iHist] -> AddText(sLabels[iHist].Data());
  }
  cout << "    Made text." << endl;

  // make lines
  const UInt_t  fColLi(1);
  const UInt_t  fLinLi(9);
  const UInt_t  fWidLi(1);
  const UInt_t  fLocLi(1.);
  const Float_t fLinX[NPts] = {xPlotRange[0], xPlotRange[1]};

  TLine   *line[NHist];
  Float_t fLinY[NHist];
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    fLinY[iHist] = fLocLi * fScales[iHist];
    line[iHist]  = new TLine(fLinX[0], fLinY[iHist], fLinX[1], fLinY[iHist]);
    line[iHist]  -> SetLineColor(fColLi);
    line[iHist]  -> SetLineStyle(fLinLi);
    line[iHist]  -> SetLineWidth(fWidLi);
  }
  cout << "    Made lines." << endl;

  // scale histograms
  hPlotStat -> Scale(fScales[0]);
  hLineStat -> Scale(fScales[0]);
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hPlotSys[iSys]       -> Scale(fScales[1]);
    hPlotSysMix[iSys]    -> Scale(fScales[1]);
    hPlotSysSmooth[iSys] -> Scale(fScales[1]);
    hPlotTot[iSys]       -> Scale(fScales[2]);
    hPlotTotMix[iSys]    -> Scale(fScales[2]);
    hPlotTotSmooth[iSys] -> Scale(fScales[2]);
    hLineSys[iSys]       -> Scale(fScales[1]);
    hLineSysMix[iSys]    -> Scale(fScales[1]);
    hLineSysSmooth[iSys] -> Scale(fScales[1]);
    hLineTot[iSys]       -> Scale(fScales[2]);
    hLineTotMix[iSys]    -> Scale(fScales[2]);
    hLineTotSmooth[iSys] -> Scale(fScales[2]);
  }
  hPlotTot[NSys]       -> Scale(fScales[2]);
  hPlotTotMix[NSys]    -> Scale(fScales[2]);
  hPlotTotSmooth[NSys] -> Scale(fScales[2]);
  hLineTot[NSys]       -> Scale(fScales[2]);
  hLineTotMix[NSys]    -> Scale(fScales[2]);
  hLineTotSmooth[NSys] -> Scale(fScales[2]);
  cout << "    Scaled histograms:\n"
       << "      stat. scale = " << fScales[0] << "\n"
       << "      sys. scale  = " << fScales[1] << "\n"
       << "      tot. scale  = " << fScales[2]
       << endl;

  // make plot
  const UInt_t  width(1500);
  const UInt_t  height(950);
  const UInt_t  smallHeight(750);
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
  const Float_t fPadXY1[NVtx] = {0., 0., 0.5, 1.};
  const Float_t fPadXY2[NVtx] = {0.5, 0., 1., 1.};

  TCanvas *cAllErrRaw = new TCanvas("cAllErrRaw", "", width, height);
  TPad    *pPad1      = new TPad("pPad1", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad2      = new TPad("pPad2", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cAllErrRaw -> SetGrid(fGrid, fGrid);
  cAllErrRaw -> SetTicks(fTick, fTick);
  cAllErrRaw -> SetBorderMode(fMode);
  cAllErrRaw -> SetBorderSize(fBord);
  pPad1      -> SetGrid(fGrid, fGrid);
  pPad1      -> SetTicks(fTick, fTick);
  pPad1      -> SetLogx(fLogX);
  pPad1      -> SetLogy(fLogY);
  pPad1      -> SetBorderMode(fMode);
  pPad1      -> SetBorderSize(fBord);
  pPad1      -> SetFrameBorderMode(fFrame);
  pPad1      -> SetLeftMargin(fMarginL);
  pPad1      -> SetRightMargin(fMarginR);
  pPad1      -> SetTopMargin(fMarginT);
  pPad1      -> SetBottomMargin(fMarginB);
  pPad2      -> SetGrid(fGrid, fGrid);
  pPad2      -> SetTicks(fTick, fTick);
  pPad2      -> SetLogx(fLogX);
  pPad2      -> SetLogy(fLogY);
  pPad2      -> SetBorderMode(fMode);
  pPad2      -> SetBorderSize(fBord);
  pPad2      -> SetFrameBorderMode(fFrame);
  pPad2      -> SetLeftMargin(fMarginL);
  pPad2      -> SetRightMargin(fMarginR);
  pPad2      -> SetTopMargin(fMarginT);
  pPad2      -> SetBottomMargin(fMarginB);
  cAllErrRaw -> cd();
  pPad1      -> Draw();
  pPad2      -> Draw();
  pPad1      -> cd();
  hPlotStat  -> Draw("E1 X0");
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    hPlotSys[iPlot] -> Draw("E2 same");
    hPlotTot[iPlot] -> Draw("E2 same");
  }
  hPlotTot[NSys] -> Draw("E1 X0 same");
  txt            -> Draw();
  legR           -> Draw();
  pPad2          -> cd();
  hLineStat      -> Draw("E1 X0");
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    hLineSys[iPlot] -> Draw("E2 same");
    hLineTot[iPlot] -> Draw("E2 same");
  }
  hLineTot[NSys] -> Draw("E1 X0 same");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    lbl[iHist]  -> Draw();
    line[iHist] -> Draw();
  }
  fOut       -> cd();
  cAllErrRaw -> Write();
  cAllErrRaw -> Close();

  TCanvas *cAllErrMix = new TCanvas("cAllErrMix", "", width, height);
  TPad    *pPad3         = new TPad("pPad3", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad4         = new TPad("pPad4", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cAllErrMix -> SetGrid(fGrid, fGrid);
  cAllErrMix -> SetTicks(fTick, fTick);
  cAllErrMix -> SetBorderMode(fMode);
  cAllErrMix -> SetBorderSize(fBord);
  pPad3      -> SetGrid(fGrid, fGrid);
  pPad3      -> SetTicks(fTick, fTick);
  pPad3      -> SetLogx(fLogX);
  pPad3      -> SetLogy(fLogY);
  pPad3      -> SetBorderMode(fMode);
  pPad3      -> SetBorderSize(fBord);
  pPad3      -> SetFrameBorderMode(fFrame);
  pPad3      -> SetLeftMargin(fMarginL);
  pPad3      -> SetRightMargin(fMarginR);
  pPad3      -> SetTopMargin(fMarginT);
  pPad3      -> SetBottomMargin(fMarginB);
  pPad4      -> SetGrid(fGrid, fGrid);
  pPad4      -> SetTicks(fTick, fTick);
  pPad4      -> SetLogx(fLogX);
  pPad4      -> SetLogy(fLogY);
  pPad4      -> SetBorderMode(fMode);
  pPad4      -> SetBorderSize(fBord);
  pPad4      -> SetFrameBorderMode(fFrame);
  pPad4      -> SetLeftMargin(fMarginL);
  pPad4      -> SetRightMargin(fMarginR);
  pPad4      -> SetTopMargin(fMarginT);
  pPad4      -> SetBottomMargin(fMarginB);
  cAllErrMix -> cd();
  pPad3      -> Draw();
  pPad4      -> Draw();
  pPad3      -> cd();
  hPlotStat  -> Draw("E1 X0");
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    hPlotSysMix[iPlot] -> Draw("E2 same");
    hPlotTotMix[iPlot] -> Draw("E2 same");
  }
  hPlotTotMix[NSys] -> Draw("E1 X0 same");
  txt               -> Draw();
  legM              -> Draw();
  pPad4             -> cd();
  hLineStat         -> Draw("E1 X0");
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    hLineSysMix[iPlot] -> Draw("E2 same");
    hLineTotMix[iPlot] -> Draw("E2 same");
  }
  hLineTot[NSys] -> Draw("E1 X0 same");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    lbl[iHist]  -> Draw();
    line[iHist] -> Draw();
  }
  fOut       -> cd();
  cAllErrMix -> Write();
  cAllErrMix -> Close();

  TCanvas *cAllErrSmooth = new TCanvas("cAllErrSmooth", "", width, height);
  TPad    *pPad5         = new TPad("pPad5", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad6         = new TPad("pPad6", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cAllErrSmooth -> SetGrid(fGrid, fGrid);
  cAllErrSmooth -> SetTicks(fTick, fTick);
  cAllErrSmooth -> SetBorderMode(fMode);
  cAllErrSmooth -> SetBorderSize(fBord);
  pPad5         -> SetGrid(fGrid, fGrid);
  pPad5         -> SetTicks(fTick, fTick);
  pPad5         -> SetLogx(fLogX);
  pPad5         -> SetLogy(fLogY);
  pPad5         -> SetBorderMode(fMode);
  pPad5         -> SetBorderSize(fBord);
  pPad5         -> SetFrameBorderMode(fFrame);
  pPad5         -> SetLeftMargin(fMarginL);
  pPad5         -> SetRightMargin(fMarginR);
  pPad5         -> SetTopMargin(fMarginT);
  pPad5         -> SetBottomMargin(fMarginB);
  pPad6         -> SetGrid(fGrid, fGrid);
  pPad6         -> SetTicks(fTick, fTick);
  pPad6         -> SetLogx(fLogX);
  pPad6         -> SetLogy(fLogY);
  pPad6         -> SetBorderMode(fMode);
  pPad6         -> SetBorderSize(fBord);
  pPad6         -> SetFrameBorderMode(fFrame);
  pPad6         -> SetLeftMargin(fMarginL);
  pPad6         -> SetRightMargin(fMarginR);
  pPad6         -> SetTopMargin(fMarginT);
  pPad6         -> SetBottomMargin(fMarginB);
  cAllErrSmooth -> cd();
  pPad5         -> Draw();
  pPad6         -> Draw();
  pPad5         -> cd();
  hPlotStat     -> Draw("E1 X0");
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    hPlotSysSmooth[iPlot] -> Draw("E2 same");
    hPlotTotSmooth[iPlot] -> Draw("E2 same");
  }
  hPlotTotSmooth[NSys] -> Draw("E1 X0 same");
  txt                  -> Draw();
  legS                 -> Draw();
  pPad6                -> cd();
  hLineStat            -> Draw("E1 X0");
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    hLineSysSmooth[iPlot] -> Draw("E2 same");
    hLineTotSmooth[iPlot] -> Draw("E2 same");
  }
  hLineTot[NSys] -> Draw("E1 X0 same");
  for (UInt_t iHist = 0; iHist < NHist; iHist++) {
    lbl[iHist]  -> Draw();
    line[iHist] -> Draw();
  }
  fOut          -> cd();
  cAllErrSmooth -> Write();
  cAllErrSmooth -> Close();
  cout << "    Made compilation plots." << endl;

  // de-scale histograms
  hPlotStat -> Scale(1. / fScales[0]);
  hLineStat -> Scale(1. / fScales[0]);
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hPlotSys[iSys]       -> Scale(1. / fScales[1]);
    hPlotSysMix[iSys]    -> Scale(1. / fScales[1]);
    hPlotSysSmooth[iSys] -> Scale(1. / fScales[1]);
    hPlotTot[iSys]       -> Scale(1. / fScales[2]);
    hPlotTotMix[iSys]    -> Scale(1. / fScales[2]);
    hPlotTotSmooth[iSys] -> Scale(1. / fScales[2]);
    hLineSys[iSys]       -> Scale(1. / fScales[1]);
    hLineSysMix[iSys]    -> Scale(1. / fScales[1]);
    hLineSysSmooth[iSys] -> Scale(1. / fScales[1]);
    hLineTot[iSys]       -> Scale(1. / fScales[2]);
    hLineTotMix[iSys]    -> Scale(1. / fScales[2]);
    hLineTotSmooth[iSys] -> Scale(1. / fScales[2]);
  }
  hPlotTot[NSys]       -> Scale(1. / fScales[2]);
  hPlotTotMix[NSys]    -> Scale(1. / fScales[2]);
  hPlotTotSmooth[NSys] -> Scale(1. / fScales[2]);
  hLineTot[NSys]       -> Scale(1. / fScales[2]);
  hLineTotMix[NSys]    -> Scale(1. / fScales[2]);
  hLineTotSmooth[NSys] -> Scale(1. / fScales[2]);
  cout << "    De-scaled histograms." << endl;

  TCanvas *cOnlySysRaw = new TCanvas("cOnlySysRaw", "", width, smallHeight);
  TPad    *pPad7       = new TPad("pPad7", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad8       = new TPad("pPad8", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cOnlySysRaw -> SetGrid(fGrid, fGrid);
  cOnlySysRaw -> SetTicks(fTick, fTick);
  cOnlySysRaw -> SetBorderMode(fMode);
  cOnlySysRaw -> SetBorderSize(fBord);
  pPad7       -> SetGrid(fGrid, fGrid);
  pPad7       -> SetTicks(fTick, fTick);
  pPad7       -> SetLogx(fLogX);
  pPad7       -> SetLogy(fLogY);
  pPad7       -> SetBorderMode(fMode);
  pPad7       -> SetBorderSize(fBord);
  pPad7       -> SetFrameBorderMode(fFrame);
  pPad7       -> SetLeftMargin(fMarginL);
  pPad7       -> SetRightMargin(fMarginR);
  pPad7       -> SetTopMargin(fMarginT);
  pPad7       -> SetBottomMargin(fMarginB);
  pPad8       -> SetGrid(fGrid, fGrid);
  pPad8       -> SetTicks(fTick, fTick);
  pPad8       -> SetLogx(fLogX);
  pPad8       -> SetLogy(fLogY);
  pPad8       -> SetBorderMode(fMode);
  pPad8       -> SetBorderSize(fBord);
  pPad8       -> SetFrameBorderMode(fFrame);
  pPad8       -> SetLeftMargin(fMarginL);
  pPad8       -> SetRightMargin(fMarginR);
  pPad8       -> SetTopMargin(fMarginT);
  pPad8       -> SetBottomMargin(fMarginB);
  cOnlySysRaw -> cd();
  pPad7       -> Draw();
  pPad8       -> Draw();
  pPad7       -> cd();
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    if (iPlot == (NSys - 1)) {
      hPlotSys[iPlot] -> Draw("E2");
    } else {
      hPlotSys[iPlot] -> Draw("E2 same");
    }
  }
  txt   -> Draw();
  sysR  -> Draw();
  pPad8 -> cd();
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    if (iPlot == (NSys - 1)) {
      hLineSys[iPlot] -> Draw("E2 same");
    } else {
      hLineSys[iPlot] -> Draw("E2 same");
    }
  }
  lbl[1]      -> Draw();
  line[0]     -> Draw();
  fOut        -> cd();
  cOnlySysRaw -> Write();
  cOnlySysRaw -> Close();

  TCanvas *cOnlySysMix = new TCanvas("cOnlySysMix", "", width, smallHeight);
  TPad    *pPad9       = new TPad("pPad9", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad10      = new TPad("pPad10", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cOnlySysMix -> SetGrid(fGrid, fGrid);
  cOnlySysMix -> SetTicks(fTick, fTick);
  cOnlySysMix -> SetBorderMode(fMode);
  cOnlySysMix -> SetBorderSize(fBord);
  pPad9       -> SetGrid(fGrid, fGrid);
  pPad9       -> SetTicks(fTick, fTick);
  pPad9       -> SetLogx(fLogX);
  pPad9       -> SetLogy(fLogY);
  pPad9       -> SetBorderMode(fMode);
  pPad9       -> SetBorderSize(fBord);
  pPad9       -> SetFrameBorderMode(fFrame);
  pPad9       -> SetLeftMargin(fMarginL);
  pPad9       -> SetRightMargin(fMarginR);
  pPad9       -> SetTopMargin(fMarginT);
  pPad9       -> SetBottomMargin(fMarginB);
  pPad10      -> SetGrid(fGrid, fGrid);
  pPad10      -> SetTicks(fTick, fTick);
  pPad10      -> SetLogx(fLogX);
  pPad10      -> SetLogy(fLogY);
  pPad10      -> SetBorderMode(fMode);
  pPad10      -> SetBorderSize(fBord);
  pPad10      -> SetFrameBorderMode(fFrame);
  pPad10      -> SetLeftMargin(fMarginL);
  pPad10      -> SetRightMargin(fMarginR);
  pPad10      -> SetTopMargin(fMarginT);
  pPad10      -> SetBottomMargin(fMarginB);
  cOnlySysMix -> cd();
  pPad9       -> Draw();
  pPad10      -> Draw();
  pPad9       -> cd();
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    if (iPlot == (NSys - 1)) {
      hPlotSysMix[iPlot] -> Draw("E2");
    } else {
      hPlotSysMix[iPlot] -> Draw("E2 same");
    }
  }
  txt    -> Draw();
  sysM   -> Draw();
  pPad10 -> cd();
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    if (iPlot == (NSys - 1)) {
      hLineSysMix[iPlot] -> Draw("E2 same");
    } else {
      hLineSysMix[iPlot] -> Draw("E2 same");
    }
  }
  lbl[1]      -> Draw();
  line[0]     -> Draw();
  fOut        -> cd();
  cOnlySysMix -> Write();
  cOnlySysMix -> Close();

  TCanvas *cOnlySysSmooth = new TCanvas("cOnlySysSmooth", "", width, smallHeight);
  TPad    *pPad11         = new TPad("pPad11", "", fPadXY1[0], fPadXY1[1], fPadXY1[2], fPadXY1[3]);
  TPad    *pPad12         = new TPad("pPad12", "", fPadXY2[0], fPadXY2[1], fPadXY2[2], fPadXY2[3]);
  cOnlySysSmooth -> SetGrid(fGrid, fGrid);
  cOnlySysSmooth -> SetTicks(fTick, fTick);
  cOnlySysSmooth -> SetBorderMode(fMode);
  cOnlySysSmooth -> SetBorderSize(fBord);
  pPad11         -> SetGrid(fGrid, fGrid);
  pPad11         -> SetTicks(fTick, fTick);
  pPad11         -> SetLogx(fLogX);
  pPad11         -> SetLogy(fLogY);
  pPad11         -> SetBorderMode(fMode);
  pPad11         -> SetBorderSize(fBord);
  pPad11         -> SetFrameBorderMode(fFrame);
  pPad11         -> SetLeftMargin(fMarginL);
  pPad11         -> SetRightMargin(fMarginR);
  pPad11         -> SetTopMargin(fMarginT);
  pPad11         -> SetBottomMargin(fMarginB);
  pPad12         -> SetGrid(fGrid, fGrid);
  pPad12         -> SetTicks(fTick, fTick);
  pPad12         -> SetLogx(fLogX);
  pPad12         -> SetLogy(fLogY);
  pPad12         -> SetBorderMode(fMode);
  pPad12         -> SetBorderSize(fBord);
  pPad12         -> SetFrameBorderMode(fFrame);
  pPad12         -> SetLeftMargin(fMarginL);
  pPad12         -> SetRightMargin(fMarginR);
  pPad12         -> SetTopMargin(fMarginT);
  pPad12         -> SetBottomMargin(fMarginB);
  cOnlySysSmooth -> cd();
  pPad11         -> Draw();
  pPad12         -> Draw();
  pPad11         -> cd();
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    if (iPlot == (NSys - 1)) {
      hPlotSysSmooth[iPlot] -> Draw("E2");
    } else {
      hPlotSysSmooth[iPlot] -> Draw("E2 same");
    }
  }
  txt    -> Draw();
  sysS   -> Draw();
  pPad12 -> cd();
  for (Int_t iPlot = NSys - 1; iPlot > -1; iPlot--) {
    if (iPlot == (NSys - 1)) {
      hLineSysSmooth[iPlot] -> Draw("E2 same");
    } else {
      hLineSysSmooth[iPlot] -> Draw("E2 same");
    }
  }
  lbl[1]         -> Draw();
  line[0]        -> Draw();
  fOut           -> cd();
  cOnlySysSmooth -> Write();
  cOnlySysSmooth -> Close();
  cout << "    Made sys.-only plots." << endl;

  // write histograms
  fOut          -> cd();
  hStatOut      -> Write();
  hSysOut       -> Write();
  hSysOutMix    -> Write();
  hSysOutSmooth -> Write();
  hTotOut       -> Write();
  hTotOutMix    -> Write();
  hTotOutSmooth -> Write();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    hSys[iSys]           -> Write();
    hSysSmooth[iSys]     -> Write();
    hPlotSys[iSys]       -> Write();
    hPlotSysMix[iSys]    -> Write();
    hPlotSysSmooth[iSys] -> Write();
    hPlotTot[iSys]       -> Write();
    hPlotTotMix[iSys]    -> Write();
    hPlotTotSmooth[iSys] -> Write();
    hLineSys[iSys]       -> Write();
    hLineSysMix[iSys]    -> Write();
    hLineSysSmooth[iSys] -> Write();
    hLineTot[iSys]       -> Write();
    hLineTotMix[iSys]    -> Write();
    hLineTotSmooth[iSys] -> Write();
  }
  hPlotStat            -> Write();
  hLineStat            -> Write();
  hPlotTot[NSys]       -> Write();
  hPlotTotMix[NSys]    -> Write();
  hPlotTotSmooth[NSys] -> Write();
  hLineTot[NSys]       -> Write();
  hLineTotMix[NSys]    -> Write();
  hLineTotSmooth[NSys] -> Write();
  cout << "    Saved histograms." << endl;

  // close files
  fOut  -> Close();
  fStat -> cd();
  fStat -> Close();
  for (UInt_t iSys = 0; iSys < NSys; iSys++) {
    fSys[iSys] -> cd();
    fSys[iSys] -> Close();
  }
  cout << "  Summing errors complete." << endl;

}

// End ------------------------------------------------------------------------
