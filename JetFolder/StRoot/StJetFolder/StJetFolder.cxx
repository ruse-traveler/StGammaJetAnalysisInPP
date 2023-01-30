// 'StJetFolder.cxx'
// Derek Anderson
// 07.15.2018
//
// This class handles the unfolding of a provided spectrum.  This file
// contains the 'Init()', 'Unfold()', 'Backfold()', and 'Finish()'
// routines.   Pearson Coefficient calculation adapted from Rhagav K.
// Elayavalli.
//
// Last updated: 06.18.2021


#define StJetFolder_cxx

// user includes
#include "StJetFolder.h"
#include "StJetFolder.io.h"
#include "StJetFolder.sys.h"
#include "StJetFolder.math.h"
#include "StJetFolder.plot.h"

ClassImp(StJetFolder)

using namespace std;



void StJetFolder::Init() {

  Bool_t inputOK = CheckFlags();
  if (!inputOK) assert(inputOK);

  // rebin histograms according to measured distribution
  const UInt_t nPri      = _hPrior      -> GetNbinsX();
  const UInt_t nSme      = _hSmeared    -> GetNbinsX();
  const UInt_t nMeas     = _hMeasured   -> GetNbinsX();
  const Bool_t rebinPri  = (nPri > nMeas);
  const Bool_t rebinSme  = (nSme > nMeas);
  if (rebinPri)  DoRebinning(_hPrior, _hMeasured);
  if (rebinSme)  DoRebinning(_hSmeared, _hMeasured);

  // create smoothed/alternate response
  if (_smoothPrior) {
    SmoothPrior();
    _response = new RooUnfoldResponse(0, 0, _hResponseDiff);
  }
  else {
    if (_smoothResponse) {
      SmoothPrior();
      SmoothResponse();
      _response = new RooUnfoldResponse(0, 0, _hResponseSmooth);
    }
    if (_differentPrior) {
      InitializePriors();
      _response = new RooUnfoldResponse(0, 0, _hResponseDiff);
    }
  }

  // otherwise create default response
  if (!_smoothPrior && !_smoothResponse && !_differentPrior) {
    _response = new RooUnfoldResponse(0, 0, _hResponse);
  }

  if (_response) {
    PrintInfo(4);
    _flag[10] = true;
  }
  else {
    PrintError(11);
    assert(_response);
  }

}  // end 'Init()'


void StJetFolder::Unfold(Double_t &chi2unfold) {

  PrintInfo(5);

  // do unfolding
  RooUnfoldBayes    *bay;
  RooUnfoldSvd      *svd;
  RooUnfoldBinByBin *bin;
  RooUnfoldTUnfold  *tun;
  RooUnfoldInvert   *inv;
  RooUnfoldErrors   *err;
  TMatrixD          *cov;
  switch (_method) {
    case 0:
      _hUnfolded = (TH1D*) _hMeasured -> Clone("hUnfolded");
      break;
    case 1:
      bay        = new RooUnfoldBayes(_response, _hMeasured, _kReg);
      err        = new RooUnfoldErrors(_nToy, bay);
      cov        = (TMatrixD*) bay -> Ereco().Clone();
      _hUnfolded = (TH1D*)     bay -> Hreco();
      break;
    case 2:
      svd        = new RooUnfoldSvd(_response, _hMeasured, _kReg, _nToy);
      err        = new RooUnfoldErrors(_nToy, svd);
      cov        = (TMatrixD*) svd -> Ereco().Clone();
      _hUnfolded = (TH1D*)     svd -> Hreco();
      break;
    case 3:
      bin        = new RooUnfoldBinByBin(_response, _hMeasured);
      err        = new RooUnfoldErrors(_nToy, bin);
      cov        = (TMatrixD*) bin -> Ereco().Clone();
      _hUnfolded = (TH1D*)     bin -> Hreco();
      break;
    case 4:
      tun        = new RooUnfoldTUnfold(_response, _hMeasured, TUnfold::kRegModeDerivative);
      err        = new RooUnfoldErrors(_nToy, tun);
      cov        = (TMatrixD*) tun -> Ereco().Clone();
      _hUnfolded = (TH1D*)     tun -> Hreco();
      break;
    case 5:
      inv        = new RooUnfoldInvert(_response, _hMeasured);
      err        = new RooUnfoldErrors(_nToy, inv);
      cov        = (TMatrixD*) inv -> Ereco().Clone();
      _hUnfolded = (TH1D*)     inv -> Hreco();
      break;
  }

  // correct for efficiency
  if (_smoothEfficiency)
    _hUnfolded -> Divide(_fEfficiency);
  else
    _hUnfolded -> Divide(_hEfficiency);

  // calculate pearson coef.s
  if (_method != 0)
    _hPearson = GetPearsonCoefficient(cov, _pearsonDebug, "hPearson");
  else
    _hPearson = (TH2D*) _hResponse -> Clone();

  // grab errors and D vectors
  switch (_method) {
    case 0:
      _hUnfoldErrors = (TH1D*) _hUnfolded -> Clone();
      _hSVvector     = (TH1D*) _hUnfolded -> Clone();
      _hDvector      = (TH1D*) _hUnfolded -> Clone();
      break;
    case 1:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) bay -> Hreco();
      _hDvector      = (TH1D*) bay -> Hreco();
      break;
    case 2:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) svd -> Impl() -> GetSV();
      _hDvector      = (TH1D*) svd -> Impl() -> GetD();
      break;
    case 3:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) bin -> Hreco();
      _hDvector      = (TH1D*) bin -> Hreco();
      break;
    case 4:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) tun -> Hreco();
      _hDvector      = (TH1D*) tun -> Hreco();
      break;
    case 5:
      _hUnfoldErrors = (TH1D*) err -> UnfoldingError();
      _hSVvector     = (TH1D*) inv -> Hreco();
      _hDvector      = (TH1D*) inv -> Hreco();
      break;
  }


  // make sure unfolded didn't exceed max bin
  Int_t nU = _hUnfolded -> GetNbinsX();
  for (Int_t i = 1; i < nU + 1; i++) {
    Double_t uBin = _hUnfolded -> GetBinLowEdge(i);
    if (uBin > _uMax) {
      _hUnfolded -> SetBinContent(i, 0.);
      _hUnfolded -> SetBinError(i, 0.);
    }
  }

  // set negative bins to 0 (for gamma-direct triggers)
  if (_trigger == 0) {
    for (Int_t i = 1; i < nU + 1; i++) {
      const Double_t uVal = _hUnfolded -> GetBinContent(i);
      if (uVal < 0.) {
        _hUnfolded -> SetBinContent(i, 0.);
        _hUnfolded -> SetBinError(i, 0.);
      }
    }
  }

  // calculate chi2
  if (_diffPriorCreated) {
    _chi2unfold = CalculateChi2(_hPriorDiff, _hUnfolded);
  }
  else if (_priorSmoothed || _responseSmoothed) {
    _chi2unfold = CalculateChi2(_hPriorSmooth, _hUnfolded);
  }
  else {
    _chi2unfold = CalculateChi2(_hPrior, _hUnfolded);
  }
  chi2unfold  = _chi2unfold;

  PrintInfo(6);

}  // end 'Unfold(Double_t)'


void StJetFolder::Backfold(Double_t &chi2backfold, const Int_t nMcBack) {


  // set no. of mc iterations
  //_nMcBack = nMcBack;
  _nMcBack = _hMeasured -> GetEntries();  // TEST [09.30.2021]
  PrintInfo(7);

  if (_method == 0) {
    // initialize backfolding histograms
    _hBackfolded = (TH1D*) _hUnfolded -> Clone("hBackfolded");
    _hNormalize  = (TH1D*) _hUnfolded -> Clone("hNormalize");

    // apply efficiency
    if (_smoothEfficiency)
      _hBackfolded -> Multiply(_fEfficiency);
    else
      _hBackfolded -> Multiply(_hEfficiency);

    // calculate backfolding chi2
    _chi2backfold = CalculateChi2(_hMeasured, _hBackfolded);
    chi2backfold  = _chi2backfold;
    PrintInfo(8);
    return;
  }

  _hNormalize  = (TH1D*) _hUnfolded -> Clone();
  _hBackfolded = (TH1D*) _hMeasured -> Clone();
  _hNormalize  -> SetNameTitle("hNormalize", "For normalizing backfolded spectrum");
  _hBackfolded -> SetNameTitle("hBackfolded", "Backfolded spectrum");
  _hNormalize  -> Reset("ICE");
  _hBackfolded -> Reset("ICE");

  // create histogram to sample from and apply efficiency
  TH1D *hUnfoldSample = (TH1D*) _hUnfolded -> Clone();
  if (_smoothEfficiency)
    hUnfoldSample -> Multiply(_fEfficiency);
  else
    hUnfoldSample -> Multiply(_hEfficiency);

  // monte-carlo loop
  Double_t u = 0.;
  Double_t b = 0.;
  for (Int_t i = 0; i < _nMcBack; i++) {
    u = hUnfoldSample -> GetRandom();
    b = Smear(u);
    _hNormalize -> Fill(u);
    if (b > -1000.) _hBackfolded -> Fill(b);
  }

  // set negative bins to 0 (for gamma-direct triggers)
  if (_trigger == 0) {
    const UInt_t nB = _hBackfolded -> GetNbinsX();
    const UInt_t nN = _hNormalize  -> GetNbinsX();
    for (UInt_t iB = 1; iB < nB + 1; iB++) {
      const Double_t bVal = _hBackfolded -> GetBinContent(iB);
      if (bVal < 0.) {
        _hBackfolded -> SetBinContent(iB, 0.);
        _hBackfolded -> SetBinError(iB, 0.);
      }
    }
    for (UInt_t iN = 1; iN < nN + 1; iN++) {
      const Double_t nVal = _hNormalize -> GetBinContent(iN);
      if (nVal < 0.) {
        _hNormalize -> SetBinContent(iN, 0.);
        _hNormalize -> SetBinError(iN, 0.);
      }
    }
  }

  // normalize backfolded spectrum
  Double_t iU = hUnfoldSample -> Integral();
  Double_t iN = _hNormalize   -> Integral();
  if (iU > 0.) {
    Double_t scale = iU / iN;
    _hBackfolded -> Scale(scale);
  }

  // calculate chi2
  _chi2backfold = CalculateChi2(_hMeasured, _hBackfolded);
  chi2backfold  = _chi2backfold;

  PrintInfo(8);

}  // end 'Backfold(Double_t, Int_t)'


void StJetFolder::Finish() {

  // calculate ratios
  _hBackVsMeasRatio   = CalculateRatio(_hBackfolded, _hMeasured, "hBackVsMeasRatio");
  _hUnfoldVsMeasRatio = CalculateRatio(_hUnfolded, _hMeasured, "hUnfoldVsMeasRatio");
  if (_diffPriorCreated) {
    _hUnfoldVsPriRatio = CalculateRatio(_hUnfolded, _hPriorDiff, "hUnfoldVsPriDifRatio");
    _hSmearVsMeasRatio = CalculateRatio(_hSmearedDiff, _hMeasured, "hSmearDifVsMeasRatio");
    _hSmearVsPriRatio  = CalculateRatio(_hSmearedDiff, _hPriorDiff, "hSmearDifVsPriDifRatio");
  }
  else if (_priorSmoothed) {
    _hUnfoldVsPriRatio = CalculateRatio(_hUnfolded, _hPriorSmooth, "hUnfoldVsPriSmoothRatio");
    _hSmearVsMeasRatio = CalculateRatio(_hSmearedDiff, _hMeasured, "hSmearDifVsMeasRatio");
    _hSmearVsPriRatio  = CalculateRatio(_hSmearedDiff, _hPriorSmooth, "hSmearDifVsPriSmoothRatio");
  }
  else if (_responseSmoothed) {
    _hUnfoldVsPriRatio = CalculateRatio(_hUnfolded, _hPriorSmooth, "hUnfoldVsPriSmoothRatio");
    _hSmearVsMeasRatio = CalculateRatio(_hSmearedSmooth, _hMeasured, "hSmearSmoothVsMeasRatio");
    _hSmearVsPriRatio  = CalculateRatio(_hSmearedSmooth, _hPriorSmooth, "hSmearSmoothVsPriSmoothRatio");
  }
  else {
    _hUnfoldVsPriRatio = CalculateRatio(_hUnfolded, _hPrior, "hUnfoldVsPriRatio");
    _hSmearVsMeasRatio = CalculateRatio(_hSmeared, _hMeasured, "hSmearVsMeasRatio");
    _hSmearVsPriRatio  = CalculateRatio(_hSmeared, _hPrior, "hSmearVsPriRatio");
  }

  // set names
  _hPrior        -> SetName("hPrior");
  _hSmeared      -> SetName("hSmeared");
  _hMeasured     -> SetName("hMeasured");
  _hUnfolded     -> SetName("hUnfolded");
  _hDvector      -> SetName("hDvector");
  _hSVvector     -> SetName("hSVvector");
  _hUnfoldErrors -> SetName("hUnfoldErrors");
  _hEfficiency   -> SetName("hEfficiency");
  _hPearson      -> SetName("hPearson");
  _hResponse     -> SetName("hResponse");
  if (_diffPriorCreated) {
    _hPriorDiff    -> SetName("hPriorDiff");
    _hSmearedDiff  -> SetName("hSmearedDiff");
    _hResponseDiff -> SetName("hResponseDiff");
  }
  if (_priorSmoothed) {
    _hPriorSmooth  -> SetName("hPriorSmooth");
    _hSmearedDiff  -> SetName("hSmearedDiff");
    _hResponseDiff -> SetName("hResponseDiff");
  }
  if (_responseSmoothed) {
    _hPriorSmooth    -> SetName("hPriorSmooth");
    _hSmearedSmooth  -> SetName("hSmearedSmooth");
    _hResponseSmooth -> SetName("hResponseSmooth");
  }

  CreateLabel();
  CreatePlots();
  PrintInfo(11);

  // save and close file
  _fOut               -> cd();
  _hPrior             -> Write();
  _hSmeared           -> Write();
  _hMeasured          -> Write();
  _hUnfolded          -> Write();
  _hNormalize         -> Write();
  _hBackfolded        -> Write();
  _hBackVsMeasRatio   -> Write();
  _hSmearVsPriRatio   -> Write();
  _hUnfoldVsPriRatio  -> Write();
  _hSmearVsMeasRatio  -> Write();
  _hUnfoldVsMeasRatio -> Write();
  _hPearson           -> Write();
  _hDvector           -> Write();
  _hSVvector          -> Write();
  _hUnfoldErrors      -> Write();
  _hEfficiency        -> Write();
  _hResponse          -> Write();
  if (_diffPriorCreated) {
    _hPriorDiff    -> Write();
    _hSmearedDiff  -> Write();
    _hResponseDiff -> Write();
  }
  if (_priorSmoothed) {
    _hPriorSmooth  -> Write();
    _hSmearedDiff  -> Write();
    _hResponseDiff -> Write();
    _fSmoothPrior  -> Write();
  }
  if (_responseSmoothed) {
    _hPriorSmooth    -> Write();
    _hSmearedSmooth  -> Write();
    _hResponseSmooth -> Write();
    _fSmoothPrior    -> Write();
    _fParToSample    -> Write();
    for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
      _fQtSmooth[iPtPar] -> Write();
    }
  }
  if (_smoothEfficiency || _fineTuneEfficiency) {
    _fEfficiency -> Write();
    _fSmoothEff  -> Write();
    if (_fineTuneEfficiency) {
      _hResiduals   -> Write();
      _fResiduals   -> Write();
      _fFineTuneEff -> Write();
    }
  }
  if (_applyCutoff) {
    _fCutoff -> Write();
  }
  _fOut -> Close();
  PrintInfo(12);

}  // end 'Finish()'

// End ------------------------------------------------------------------------
