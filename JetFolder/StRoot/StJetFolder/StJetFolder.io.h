// 'StJetFolder.io.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates I/O routines.
//
// Last updated: 06.14.2021


#pragma once

using namespace std;



void StJetFolder::SetPrior(const Char_t *pFile, const Char_t *pName, const Bool_t doPriorSmooth, const Bool_t doPriorCutoff) {

  // set prior and smoothing cutoff flags
  _smoothPrior      = doPriorSmooth;
  _applyCutoff      = doPriorCutoff;
  _priorSmoothed    = false;
  _diffPriorCreated = false;

  TFile *fPrior = (TFile*) gROOT -> GetListOfFiles() -> FindObject(pFile);
  if (!fPrior || !fPrior->IsOpen()) {
    fPrior = new TFile(pFile);
  }

  TH1D *hPrior;
  if (fPrior) {
    hPrior = (TH1D*) fPrior -> Get(pName);
  }


  if (hPrior) {
    _hPrior  = (TH1D*) hPrior -> Clone();
    _flag[0] = true;
  }
  else {
    PrintError(0);
    assert(hPrior);
  }

}  // end 'SetPrior(Char_t*, Char_t*, Int_t, Bool_t, Bool_t)'


void StJetFolder::SetSmeared(const Char_t *sFile, const Char_t *sName) {

  TFile *fSmeared = (TFile*) gROOT -> GetListOfFiles() -> FindObject(sFile);
  if (!fSmeared || !fSmeared->IsOpen()) {
    fSmeared = new TFile(sFile);
  }

  TH1D *hSmeared;
  if (fSmeared) {
    hSmeared = (TH1D*) fSmeared -> Get(sName);
  }


  if (hSmeared) {
    _hSmeared = (TH1D*) hSmeared -> Clone();
    _flag[1]  = true;
  }
  else {
    PrintError(1);
    assert(hSmeared);
  }

}  // end 'SetSmeared(Char_t*, Char_t*)'


void StJetFolder::SetMeasured(const Char_t *mFile, const Char_t *mName) {

  TFile *fMeasured = (TFile*) gROOT -> GetListOfFiles() -> FindObject(mFile);
  if (!fMeasured || !fMeasured->IsOpen()) {
    fMeasured = new TFile(mFile);
  }

  TH1D *hMeasured;
  if (fMeasured) {
    hMeasured = (TH1D*) fMeasured -> Get(mName);
  }


  if (hMeasured) {
    _hMeasured = (TH1D*) hMeasured -> Clone();
    _flag[2]   = true;
  }
  else {
    PrintError(2);
    assert(hMeasured);
  }

}  // end 'SetMeasured(Char_t*, Char_t*)'


void StJetFolder::SetResponse(const Char_t *rFile, const Char_t *rName, const Bool_t doResSmooth) {

  _smoothResponse   = doResSmooth;
  _responseSmoothed = false;

  TFile *fResponse = (TFile*) gROOT -> GetListOfFiles() -> FindObject(rFile);
  if (!fResponse || !fResponse->IsOpen()) {
    fResponse = new TFile(rFile);
  }

  TH2D *hResponse;
  if (fResponse) {
    hResponse = (TH2D*) fResponse -> Get(rName);
  }


  if (hResponse) {
    _hResponse = (TH2D*) hResponse -> Clone();
    _flag[3]   = true;
  }
  else {
    PrintError(3);
    assert(hResponse);
  }

}  // end 'SetResponse(Char_t*, Char_t*)'


void StJetFolder::SetEfficiency(const Char_t *eFile, const Char_t *eName, const Bool_t doEffSmooth, const Bool_t doEffFineTune, const Bool_t removeErrors) {

  _smoothEfficiency   = doEffSmooth;
  _fineTuneEfficiency = doEffFineTune;
  _removeErrors       = removeErrors;

  TFile *fEfficiency = (TFile*) gROOT -> GetListOfFiles() -> FindObject(eFile);
  if (!fEfficiency || !fEfficiency->IsOpen()) {
    fEfficiency = new TFile(eFile);
  }

  TH1D *hEfficiency;
  if (fEfficiency) {
    hEfficiency = (TH1D*) fEfficiency -> Get(eName);
  }


  if (hEfficiency) {
    _hEfficiency = (TH1D*) hEfficiency -> Clone();
    _flag[4]     = true;
  } else {
    PrintError(4);
    assert(hEfficiency);
  }

  // fit parameters
  const TString sPolName("fResiduals");
  const TString sEffName("fEfficiency");
  const TString sSmooth("fSmoothEff");
  const TString sFine("fFineTuneEff");
  const TString sEffFunc("[0] * (1 - TMath::Exp(-1.*[1]*x))");
  const TString sEffFine("[0] * (1 - TMath::Exp(-1.*[1]*2)) * pol2(2)");
  const TString sEffParNames[NParEff] = {"#epsilon_{0}", "#sigma"};
  const Float_t xEffRange[NRange]     = {0., 100.};
  const Float_t xFitRange[NRange]     = {0.2, 15.};
  const Float_t pEffGuess[NParEff]    = {0.87, 4.};
  const Float_t pPolGuess[NParPol]    = {1., -1, 0.5};
  const UInt_t  fColEff(2);
  const UInt_t  fColRes(4);

  // initialize efficiency fits
  const Bool_t doSmoothing = (_smoothEfficiency || _fineTuneEfficiency);
  if (doSmoothing) {
    _fSmoothEff = new TF1(sSmooth.Data(), sEffFunc.Data(), xEffRange[0], xEffRange[1]);
    _fSmoothEff -> SetParameter(0, pEffGuess[0]);
    _fSmoothEff -> SetParameter(1, pEffGuess[1]);
    _fSmoothEff -> SetParName(0, sEffParNames[0].Data());
    _fSmoothEff -> SetParName(1, sEffParNames[1].Data());
    _fSmoothEff -> SetLineColor(fColEff);
    if (_fineTuneEfficiency) {
      _fResiduals   = new TF1(sPolName.Data(), "pol2(0)",       xEffRange[0], xEffRange[1]);
      _fFineTuneEff = new TF1(sFine.Data(),    sEffFine.Data(), xEffRange[0], xEffRange[1]);
      _fResiduals   -> SetParameter(0, pPolGuess[0]);
      _fResiduals   -> SetParameter(1, pPolGuess[1]);
      _fResiduals   -> SetParameter(2, pPolGuess[2]);
      _fFineTuneEff -> SetParameter(0, pEffGuess[0]);
      _fFineTuneEff -> SetParameter(1, pEffGuess[1]);
      _fFineTuneEff -> SetParameter(2, pPolGuess[0]);
      _fFineTuneEff -> SetParameter(3, pPolGuess[1]);
      _fFineTuneEff -> SetParameter(4, pPolGuess[2]); 
      _fResiduals   -> SetLineColor(fColRes);
      _fFineTuneEff -> SetLineColor(fColEff);
    }
  }

  // smooth efficiency or remove errors
  if (_smoothEfficiency) {

    // smooth efficiency
    PrintInfo(15);
    _hEfficiency -> Fit(sSmooth.Data(), "", "", xFitRange[0], xFitRange[1]);

    // fine tune smoothing
    if (_fineTuneEfficiency) {

      // divide efficiency by fit
      _hResiduals = (TH1D*) _hEfficiency -> Clone();
      _hResiduals -> SetName("hEffResiduals");
      _hResiduals -> Divide(_fSmoothEff);
      _hResiduals -> Fit(sPolName.Data(), "", "", xFitRange[0], xFitRange[1]);

      // extract parameters and get fit
      const Double_t amp  = _fSmoothEff -> GetParameter(0);
      const Double_t sig  = _fSmoothEff -> GetParameter(1);
      const Double_t res0 = _fResiduals -> GetParameter(0);
      const Double_t res1 = _fResiduals -> GetParameter(1);
      const Double_t res2 = _fResiduals -> GetParameter(2);
      _fFineTuneEff -> SetParameter(0, amp);
      _fFineTuneEff -> SetParameter(1, sig);
      _fFineTuneEff -> SetParameter(2, res0);
      _fFineTuneEff -> SetParameter(3, res1);
      _fFineTuneEff -> SetParameter(4, res2);
      _hEfficiency  -> Fit(sFine.Data(), "", "", xFitRange[0], xFitRange[1]);
    }

    // select efficiency
    if (_fineTuneEfficiency) {
      _fEfficiency = (TF1*) _fFineTuneEff -> Clone();
      _fEfficiency -> SetName(sEffName.Data());
    } else {
      _fEfficiency = (TF1*) _fSmoothEff -> Clone();
      _fEfficiency -> SetName(sEffName.Data());
    }
  } else if (_removeErrors) {
    PrintInfo(16);
    const UInt_t nBins = _hEfficiency -> GetNbinsX();
    for (UInt_t iBin = 1; iBin < (nBins + 1); iBin++) {
      _hEfficiency -> SetBinError(iBin, 0.);
    }  // end bin loop
  }  // end removing bin errors

}  // end 'SetEfficiency(Char_t*, Char_t*, Bool_t, Bool_t)'


void StJetFolder::SetEventInfo(const Int_t beam, const Double_t energy) {

  // no. of digits after decimal
  const Int_t nDecE = 1;


  TString bTxt("");
  TString eStr("");
  TString eTxt("");
  switch (beam) {
    case 0:
      bTxt = "pp collisions, #sqrt{s} = ";
      break;
    case 1:
      bTxt = "AuAu collisions, #sqrt{s_{NN}} = ";
      break;
  }

  eStr += energy;
  ResizeString(eStr, nDecE);
  eTxt.Append(eStr);
  eTxt.Append(" GeV");

  // combine strings
  TString evnt(bTxt);
  evnt.Append(eTxt);
  _sEvnt = new TString(evnt);


  if (_sEvnt) {
    _flag[5] = true;
  }
  else {
    PrintError(5);
    assert(_sEvnt);
  }
  

}  // end 'SetEventInfo(Int_t, Double_t)'


void StJetFolder::SetTriggerInfo(const Int_t trigger, const UInt_t eTtrgIndex, const Double_t hMax) {

  TString tTxt("");
  TString eTxt("");
  switch (trigger) {
    case 0:
      tTxt = "#gamma^{dir} trigger, ";
      break;
    case 1:
      tTxt = "#gamma^{rich} trigger, ";
      break;
    case 2:
      tTxt = "#pi^{0} trigger, ";
      break;
    case 3:
      tTxt = "h^{#pm} trigger, ";
      break;
  }
  _trigger = trigger;


  // determine eTtrg range and create label
  _eTtrgIndex = eTtrgIndex;
  switch (_eTtrgIndex) {
    case 0:
      eTxt.Append("E_{T}^{trg} #in (9, 11) GeV");
      break;
    case 1:
      eTxt.Append("E_{T}^{trg} #in (11, 15) GeV");
      break;
    case 2:
      eTxt.Append("E_{T}^{trg} #in (15, 20) GeV");
      break;
    default:
      eTxt.Append("E_{T}^{trg} #in (9, 11) GeV");
      break;
  }

  // combine strings
  TString trig(tTxt);
  trig.Append(eTxt);
  _sTrig = new TString(trig);


  if (_sTrig) {
    _flag[6] = true;
  }
  else {
    PrintError(6);
    assert(_sTrig);
  }

}  // end 'SetTriggerInfo(Int_t, UInt_t, Double_t)'


void StJetFolder::SetJetInfo(const Int_t type, const Int_t nRM, const Double_t rJet, const Double_t aMin, const Double_t pTmin) {

  const Int_t nDecR = 1;
  const Int_t nDecA = 2;
  const Int_t nDecP = 1;


  TString tTxt("");
  TString nTxt("");
  TString rStr("");
  TString rTxt("");
  TString aStr("");
  TString aTxt("");
  TString pStr("");
  TString pTxt("");
  switch (type) {
    case 0:
      tTxt = "#bf{charged jets}";
      break;
    case 1:
      tTxt = "#bf{full jets};";
      break;
  }

  nTxt += "N_{rm} = ";
  nTxt += nRM;

  rStr += rJet;
  ResizeString(rStr, nDecR);
  rTxt.Append("R = ");
  rTxt.Append(rStr);

  aStr += aMin;
  ResizeString(aStr, nDecA);
  aTxt.Append("A_{jet} > ");
  aTxt.Append(aStr);
  aTxt.Append(", ");

  pStr += pTmin;
  ResizeString(pStr, nDecP);
  pTxt.Append("p_{T}^{cst} > ");
  pTxt.Append(pStr);
  pTxt.Append(", ");

  // combine strings
  _sJet1 = new TString("anti-k_{T}, ");
  _sJet2 = new TString(aTxt);
  _sJet3 = new TString(tTxt);
  _sJet1 -> Append(rTxt);
  _sJet2 -> Append(pTxt);
  _sJet2 -> Append(nTxt);


  Bool_t jetInfoIsSet = (_sJet1 && _sJet2 && _sJet3);
  if (jetInfoIsSet) {
    _flag[7] = true;
  }
  else {
    PrintError(7);
    assert(jetInfoIsSet);
  }

}  // end 'SetJetInfo(Int_t, Int_t, Double_t, Double_t, Double_t)'


void StJetFolder::SetPriorCutoffParameters(const Double_t pTcut, const Double_t aCut) {

  _pTcutoff = pTcut;
  _aCutoff  = aCut;
  PrintInfo(13);

}  // end 'SetPriorCutoffParameters(Double_t, Double_t)'


void StJetFolder::SetPriorParameters(const Int_t prior, const Double_t bPrior, const Double_t mPrior, const Double_t nPrior, const Double_t tPrior, const Double_t pMin, const Double_t pMax) {

  _prior  = prior;
  _bPrior = bPrior;
  _mPrior = mPrior;
  _nPrior = nPrior;
  _tPrior = tPrior;
  _pMin   = pMin;
  _pMax   = pMax;
  if (_prior > 0)
    _differentPrior = true;
  else
    _differentPrior = false;

  // create prior function
  const UInt_t   iStartP = _hPrior -> FindFirstBinAbove(0.);
  const Double_t nBinsP  = _hPrior -> GetNbinsX();
  const Double_t startP  = _hPrior -> GetBinLowEdge(iStartP);
  const Double_t stopP   = _hPrior -> GetBinLowEdge(nBinsP + 1);
  _fLevy        = new TF1("fLevy", StJetFolder::Levy, startP, stopP, 4);
  _fTsallis     = new TF1("fTsallis", StJetFolder::Tsallis, startP, stopP, 3);
  _fExponential = new TF1("fExponential", StJetFolder::Exponential, startP, stopP, 2);
  _fPowerLaw    = new TF1("fPowerLaw", StJetFolder::PowerLaw, startP, stopP, 2);
  _fLandau      = new TF1("fLandau", StJetFolder::Landau, startP, stopP, 2);
  _fLevy        -> SetParameters(_bPrior, _mPrior, _nPrior, _tPrior);
  _fTsallis     -> SetParameters(_bPrior, _nPrior, _tPrior);
  _fExponential -> SetParameters(_bPrior, _tPrior);
  _fPowerLaw    -> SetParameters(_bPrior, _tPrior);
  _fLandau      -> SetParameters(_nPrior, _tPrior);


  _flag[8] = true;

}  // end 'SetPriorParameters(Int_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t)'


void StJetFolder::SetUnfoldParameters(const Int_t method, const Int_t kReg, const Int_t nToy, const Double_t uMax, const Double_t bMax) {

  _method = method;
  _kReg   = kReg;
  _nToy   = nToy;
  _uMax   = uMax;
  _bMax   = bMax;


  _flag[9] = true;

}  // end 'SetUnfoldParameters(Int_t, Int_t, Int_t, Double_t, Double_t)'

// End ------------------------------------------------------------------------

