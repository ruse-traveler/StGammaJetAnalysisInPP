// 'StJetFolder.sys.h'
// Derek Anderson
// 02.17.2017
//
// This class handles the unfolding of a provided spectrum.  This file
// encapsulates various internal routines (e.g. printing error messages).
//
// Last updated: 06.14.2021


#pragma once

using namespace std;



void StJetFolder::PrintInfo(const Int_t code) {

  switch (code) {
    case 0:
      cout << "\n  Folder created!\n"
           << "    Writing to '" << _fOut -> GetName() << "'"
           << endl;
      break;
    case 1:
      cout << "    Spectra grabbed..." << endl;
      break;
    case 2:
      cout << "    Info set..." << endl;
      break;
    case 3:
      cout << "    Prior parameters set...\n"
           << "      smoothPrior = " << _smoothPrior << ", applyCutoff = " << _applyCutoff << "\n"
           << "      differentPrior = " << _differentPrior << ", nMcSmooth = " << _nMcSmooth << "\n"
           << "      b = " << _bPrior << ", m = " << _mPrior << "\n"
           << "      n = " << _nPrior << ", t = " << _tPrior << "\n"
           << "      pMin = " << _pMin << ", pMax = " << _pMax
           << endl;
      cout << "    Unfolding parameters set...\n"
           << "      method = " << _method << ", k = " << _kReg << ", nToy = " << _nToy << "\n"
           << "      uMax = " << _uMax << ", bMax = " << _bMax
           << endl;
      break;
    case 4:
      cout << "    Folder initialized..." << endl;
      break;
    case 5:
      cout << "    Unfolding..." << endl;
      break;
    case 6:
      cout << "    Unfolding finished!\n"
           << "      Chi2 (unfold) = " << _chi2unfold
           << endl;
      break;
    case 7:
      cout << "    Backfolding...\n"
           << "      nMcBack = " << _nMcBack
           << endl;
      break;
    case 8:
      cout << "    Backfolding finished!\n"
           << "      Chi2 (backfold) = " << _chi2backfold
           << endl;
      break;
    case 9:
      cout << "    Ratios calculated!" << endl;
      break;
    case 10:
      cout << "    Creating plots..." << endl;
      break;
    case 11:
      cout << "    Plots created; saving..." << endl;
      break;
    case 12:
      cout << "  Folding finished!\n" << endl;
      break;
    case 13:
      cout << "    Set prior cutoff parameters:\n"
           << "      pTcut = " << _pTcutoff << ", aCut = " << _aCutoff
           << endl;
      break;
    case 14:
      cout << "    Smoothing prior:" << endl;
      break;
    case 15:
      cout << "    Smoothing efficiency..." << endl;
      break;
    case 16:
      cout << "    Removing errors on efficiency..." << endl;
      break;
    case 17:
      cout << "      Using 2 exponent smoothing..." << endl;
      break;
    case 18:
      cout << "      Using 3 exponent smoothing..." << endl;
      break;
    case 19:
      cout << "      Using 4 exponent smoothing..." << endl;
      break;
    case 20:
      cout << "    Smoothing response." << endl;
      break;
  }

}  // end 'PrintInfo(Int_t)'


void StJetFolder::PrintError(const Int_t code) {

  switch (code) {
    case 0:
      cerr << "PANIC: couldn't grab prior spectrum!" << endl;
      break;
    case 1:
      cerr << "PANIC: couldn't grab smeared spectrum!" << endl;
      break;
    case 2:
      cerr << "PANIC: couldn't grab measured spectrum!" << endl;
      break;
    case 3:
      cerr << "PANIC: couldn't grab response matrix!" << endl;
      break;
    case 4:
      cerr << "PANIC: couldn't grab effiency!" << endl;
      break;
    case 5:
      cerr << "PANIC: event info couldn't be set!" << endl;
      break;
    case 6:
      cerr << "PANIC: trigger info couldn't be set!" << endl;
      break;
    case 7:
      cerr << "PANIC: jet info couldn't be set!" << endl;
      break;
    case 8:
      cerr << "PANIC: at least one specturm not set!" << endl;
      break;
    case 9:
      cerr << "PANIC: some info not set!" << endl;
      break;
    case 10:
      cerr << "PANIC: some parameters not set!" << endl;
      break;
    case 11:
      cerr << "PANIC: folder couldn't be initialized!" << endl;
      break;
    case 12:
      cerr << "PANIC: trying to take ratio of 2 histograms with different dimensions!" << endl;
      break;
    case 13:
      cerr << "PANIC: trying to rebin a histogram with *less* bins than the new binning scheme!" << endl;
      break;
    case 14:
      cerr << "WARNING: 11 - 15 GeV (3 exponent) fit failed! Using default parameters..." << endl;
      break;
  }

}  // end 'PrintInfo(Int_t)'


void StJetFolder::InitializePriors() {

  // for normalization
  TH1D *hPriorNorm = (TH1D*) _hPrior   -> Clone();
  TH1D *hSmearNorm = (TH1D*) _hSmeared -> Clone();
  hSmearNorm -> Reset("ICE");

  // initialize alternate prior and smeared
  _hPriorDiff   = (TH1D*) _hPrior   -> Clone();
  _hSmearedDiff = (TH1D*) _hSmeared -> Clone();
  _hPriorDiff   -> SetTitle("Alternate prior");
  _hSmearedDiff -> SetTitle("Alternate smeared");
  _hPriorDiff   -> Reset("ICE");
  _hSmearedDiff -> Reset("ICE");

  // initialize alternate response
  _hResponseDiff = (TH2D*) _hResponse   -> Clone();
  _hResponseDiff -> SetTitle("Alternate response matrix");
  _hResponseDiff -> Reset("ICE");

  // determine prior bounds
  const UInt_t   iMin  = hPriorNorm -> FindFirstBinAbove(0.);
  const UInt_t   iMax  = hPriorNorm -> FindLastBinAbove(0.);
  const Double_t xMin  = hPriorNorm -> GetBinLowEdge(iMin);
  const Double_t xMax  = hPriorNorm -> GetBinLowEdge(iMax + 1);

  // create prior cutoff function
  if (_applyCutoff) {
    const Double_t xCutMin = TMath::Max(xMin, _pMin);
    const Double_t xCutMax = TMath::Min(xMax, _pMax);
    _fCutoff = new TF1("fCutoff", "tanh((x-[0])/[1])", xCutMin, xCutMax);
    _fCutoff -> SetParameter(0, _pTcutoff);
    _fCutoff -> SetParameter(1, _aCutoff);
  }

  // create particle-level prior
  switch (_prior) {
    case 1:
      _hPriorDiff -> Reset("ICE");
      for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
        const Double_t l = _fLevy -> GetRandom(xMin, xMax);
        _hPriorDiff -> Fill(l);
      }
      break;
    case 2:
      _hPriorDiff -> Reset("ICE");
      for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
        const Double_t t = _fTsallis -> GetRandom(xMin, xMax);
        _hPriorDiff -> Fill(t);
      }
      break;
    case 3:
      _hPriorDiff -> Reset("ICE");
      for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) { 
        const Double_t e = _fExponential -> GetRandom(xMin, xMax);
        _hPriorDiff -> Fill(e);
      }
      break;
    case 4:
      _hPriorDiff -> Reset("ICE");
      for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
        const Double_t p = _fPowerLaw -> GetRandom(xMin, xMax);
        _hPriorDiff -> Fill(p);
      }
      break;
    case 5:
      _hPriorDiff -> Reset("ICE");
     for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
       const Double_t a = _fLandau -> GetRandom(xMin, xMax);
       _hPriorDiff -> Fill(a);
     }
     break;
    default:
      _hPriorDiff -> Reset("ICE");
      for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
        const Double_t d = _fLevy -> GetRandom(xMin, xMax);
        _hPriorDiff -> Fill(d);
      }
      break;
  }

  // apply cutoff if necessary and trim prior
  if (_applyCutoff) {
    _hPriorDiff -> Multiply(_fCutoff);
  }
  ConstrainPriorRange(0);


  // normalize alternate prior
  const Double_t iPar   = hPriorNorm  -> Integral();
  const Double_t iParMC = _hPriorDiff -> Integral();
  const Double_t scaleP = iPar / iParMC;
  if (iParMC > 0.) _hPriorDiff -> Scale(scaleP);

  // normalize by bin width
  const UInt_t nBinsP = _hPriorDiff -> GetNbinsX();
  for (UInt_t iBinP = 1; iBinP < (nBinsP + 1); iBinP++) {
    const Double_t pVal = _hPriorDiff -> GetBinContent(iBinP);
    const Double_t pErr = _hPriorDiff -> GetBinError(iBinP);
    const Double_t pBin = _hPriorDiff -> GetBinWidth(iBinP);
    _hPriorDiff -> SetBinContent(iBinP, pVal / pBin);
    _hPriorDiff -> SetBinError(iBinP, pErr / pBin);
  }


  // create histogram to sample from and apply efficiency
  TH1D *hAltPriorSample = (TH1D*) _hPriorDiff -> Clone();
  if (_smoothEfficiency)
    hAltPriorSample -> Multiply(_fEfficiency);
  else
    hAltPriorSample -> Multiply(_hEfficiency);

  // create detector-level prior and response
  _hSmearedDiff -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t p = hAltPriorSample -> GetRandom();
    const Double_t s = Smear(p);
    // quick fix [02.04.2022]
    //if (s > 2.) {
    if (s > -1000.) {
      const UInt_t   i = hPriorNorm  -> GetBin(p);
      const Double_t o = hPriorNorm  -> GetBinContent(i);
      const Double_t n = _hPriorDiff -> GetBinContent(i);
      const Double_t w = o / n;
      _hSmearedDiff -> Fill(s);
      if ((o == 0.) || (n == 0.))
        _hResponseDiff -> Fill(s, p, 1.);
      else
        _hResponseDiff -> Fill(s, p, w);
    }
    hSmearNorm -> Fill(p);
  }


  // normalize response
  const UInt_t nXbins = _hResponseDiff -> GetNbinsX();
  const UInt_t nYbins = _hResponseDiff -> GetNbinsY();
  for (UInt_t iBinY = 0; iBinY < nYbins; iBinY++) {
    const Double_t binNorm = _hResponseDiff -> Integral(1, nXbins, iBinY, iBinY);
    if (binNorm != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbins; iBinX++) {
        const Double_t binVal = _hResponseDiff -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = _hResponseDiff -> GetBinError(iBinX, iBinY);
        const Double_t newVal = (binVal / binNorm);
        const Double_t newErr = (binErr / binNorm);
        _hResponseDiff -> SetBinContent(iBinX, iBinY, newVal);
        _hResponseDiff -> SetBinError(iBinX, iBinY, newErr);
      }  // end x loop
    }
  }  // end y loop


  const Double_t iNorm  = _hSmeared  -> Integral();
  const Double_t iDetMC = hSmearNorm -> Integral();
  const Double_t scaleS = iNorm / iDetMC;
  if (iNorm > 0.) _hSmearedDiff -> Scale(scaleS);

  // indicate that alternate prior has been created
  _diffPriorCreated = true;

}  // end 'InitializePriors()'


void StJetFolder::SmoothPrior() {

  PrintInfo(14);
  switch (_eTtrgIndex) {
    case 0:
      DoTwoExpoSmoothing();
      break;
    case 1:
      DoThreeExpoSmoothing();
      break;
    case 2:
      DoFourExpoSmoothing();
      break;
    default:
      DoTwoExpoSmoothing();
      break;
  }

}  // end 'SmoothPrior()'


void StJetFolder::DoTwoExpoSmoothing() {

  const UInt_t NExp2(2);
  PrintInfo(17);

  // fit parameters [9 - 11] GeV
  const TString  sFitHistExp("hFitPrior_onlyExpos");
  const TString  sFitHistTot("hFitPrior_totalFit");
  const TString  sTotName("fTotalFit");
  const TString  sTotFunc("(expo(0) + expo(2)) * tanh((x-[4])/[5])");
  const TString  sExpNames[NExp2]          = {"fExpSoft", "fExpHard"};
  const TString  sExpFuncs[NExp2]          = {"expo(0)", "expo(0)"};
  const TString  sParTanhNames[NParTan]    = {"p_{T}^{0}", "a"};
  const Double_t xExpRange[NExp2][NPts]    = {{0.6, 2.}, {2., 40.}};
  const Double_t pExpGuess[NExp2][NParTan] = {{0.25, -1.03}, {-1.76, -0.26}};
  const Double_t pTanhGuess[NParTan]       = {2., 1.};
  const Double_t xTanhRange[NPts]          = {0.1, 57.};
  const UInt_t   fColExp[NExp2]            = {2, 4};
  const UInt_t   fColTot(6);


  // for normalization
  TH1D *hPriorNorm = (TH1D*) _hPrior   -> Clone();
  TH1D *hSmearNorm = (TH1D*) _hSmeared -> Clone();
  hSmearNorm -> Reset("ICE");

  // initialize alternate prior and smeared
  _hPriorSmooth = (TH1D*) _hPrior   -> Clone();
  _hSmearedDiff = (TH1D*) _hSmeared -> Clone();
  _hPriorSmooth -> SetTitle("Smoothed prior");
  _hSmearedDiff -> SetTitle("Alternate smeared (via smoothed prior)");
  _hPriorSmooth -> Reset("ICE");
  _hSmearedDiff -> Reset("ICE");

  // initialize alternate response
  _hResponseDiff = (TH2D*) _hResponse -> Clone();
  _hResponseDiff -> SetTitle("Alternate response matrix (via smoothed prior)");
  _hResponseDiff -> Reset("ICE");

  // determine prior bounds
  const UInt_t   iMin  = hPriorNorm -> FindFirstBinAbove(0.);
  const UInt_t   iMax  = hPriorNorm -> FindLastBinAbove(0.);
  const Double_t xMin  = hPriorNorm -> GetBinLowEdge(iMin);
  const Double_t xMax  = hPriorNorm -> GetBinLowEdge(iMax + 1);

  // create prior cutoff function
  if (_applyCutoff) {
    const Double_t xCutMin = TMath::Max(xMin, _pMin);
    const Double_t xCutMax = TMath::Min(xMax, _pMax);
    _fCutoff = new TF1("fCutoff", "tanh((x-[0])/[1])", xCutMin, xCutMax);
    _fCutoff -> SetParameter(0, _pTcutoff);
    _fCutoff -> SetParameter(1, _aCutoff);
  }


  // declare histograms for fitting
  TH1D *hFitExp = (TH1D*) _hPrior -> Clone();
  hFitExp -> SetName(sFitHistExp.Data());


  // initialize 1st set of exponentials
  TF1 *fExp[NExp2];
  for (UInt_t iExp = 0; iExp < NExp2; iExp++) {
    fExp[iExp] = new TF1(sExpNames[iExp].Data(), sExpFuncs[iExp].Data(), xExpRange[iExp][0], xExpRange[iExp][1]);
    fExp[iExp] -> SetParameter(0, pExpGuess[iExp][0]);
    fExp[iExp] -> SetParameter(1, pExpGuess[iExp][1]);
    fExp[iExp] -> SetLineColor(fColExp[iExp]);
  }

  // determine exponential parameters
  Double_t pExpFit[NExp2][NParTan];
  for (UInt_t iExp = 0; iExp < NExp2; iExp++) {
    if (iExp == 0) {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R");
    }
    else {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R+");
    }
    pExpFit[iExp][0] = fExp[iExp] -> GetParameter(0);
    pExpFit[iExp][1] = fExp[iExp] -> GetParameter(1);
  }

  // save exponent fits to output file
  _fOut   -> cd();
  hFitExp -> Write();

  // initialize total fit function
  _fSmoothPrior = new TF1(sTotName.Data(), sTotFunc.Data(), xTanhRange[0], xTanhRange[1]);
  _fSmoothPrior -> SetLineColor(fColTot);

  const UInt_t iParTanh = 2 * NExp2;
  for (UInt_t iExp = 0; iExp < NExp2; iExp++) {
    const UInt_t  iParExp = 2 * iExp;
    _fSmoothPrior -> FixParameter(iParExp + 0, pExpFit[iExp][0]);
    _fSmoothPrior -> FixParameter(iParExp + 1, pExpFit[iExp][1]);
  }
  _fSmoothPrior -> SetParameter(iParTanh + 0, pTanhGuess[0]);
  _fSmoothPrior -> SetParameter(iParTanh + 1, pTanhGuess[1]);
  _fSmoothPrior -> SetParName(iParTanh + 0, sParTanhNames[0].Data());
  _fSmoothPrior -> SetParName(iParTanh + 1, sParTanhNames[1].Data());

  // fit function
  _hPrior -> Fit(sTotName.Data(), "B", "", xExpRange[0][0], xExpRange[NExp2 - 1][1]);

  // fill prior
  _hPriorSmooth -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t t = _fSmoothPrior -> GetRandom(xMin, xMax);
    _hPriorSmooth -> Fill(t);
  }

  // apply cutoff if necessary and trim prior
  if (_applyCutoff) {
    _hPriorSmooth -> Multiply(_fCutoff);
  }
  ConstrainPriorRange(1);


  // normalize smoothed prior
  const Double_t xParStart = TMath::Max(xMin, xTanhRange[0]);
  const Double_t xParStop  = TMath::Min(xMax, xTanhRange[1]);
  const Double_t intPar    = _fSmoothPrior -> Integral(xParStart, xParStop);
  const Double_t intParMC  = _hPriorSmooth -> Integral();
  const Double_t scaleP    = intPar / intParMC;
  if (intParMC > 0.) _hPriorSmooth -> Scale(scaleP);

  // normalize by bin width
  const UInt_t nBinsP = _hPriorSmooth -> GetNbinsX();
  for (UInt_t iBinP = 1; iBinP < (nBinsP + 1); iBinP++) {
    const Double_t pVal = _hPriorSmooth -> GetBinContent(iBinP);
    const Double_t pErr = _hPriorSmooth -> GetBinError(iBinP);
    const Double_t pBin = _hPriorSmooth -> GetBinWidth(iBinP);
    _hPriorSmooth -> SetBinContent(iBinP, pVal / pBin);
    _hPriorSmooth -> SetBinError(iBinP, pErr / pBin);
  }

  // only need smoothed prior if smoothing response
  if (!_smoothPrior && _smoothResponse) return;


  // create histogram to sample from and apply efficiency
  TH1D *hAltPriorSample = (TH1D*) _hPriorSmooth -> Clone();
  if (_smoothEfficiency)
    hAltPriorSample -> Multiply(_fEfficiency);
  else
    hAltPriorSample -> Multiply(_hEfficiency);

  // create detector-level prior and response
  _hSmearedDiff -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t p = hAltPriorSample -> GetRandom();
    const Double_t s = Smear(p);
    // quick fix [02.04.2022]
    //if (s > 2.) {
    if (s > -1000.) {
      const UInt_t   i = hPriorNorm    -> GetBin(p);
      const Double_t o = hPriorNorm    -> GetBinContent(i);
      const Double_t n = _hPriorSmooth -> GetBinContent(i);
      const Double_t w = o / n;
      _hSmearedDiff -> Fill(s);
      if ((o == 0.) || (n == 0.))
        _hResponseDiff -> Fill(s, p, 1.);
      else
        _hResponseDiff -> Fill(s, p, w);
    }
    hSmearNorm -> Fill(p);
  }


  // normalize response
  const UInt_t nXbins = _hResponseDiff -> GetNbinsX();
  const UInt_t nYbins = _hResponseDiff -> GetNbinsY();
  for (UInt_t iBinY = 0; iBinY < nYbins; iBinY++) {
    const Double_t binNorm = _hResponseDiff -> Integral(1, nXbins, iBinY, iBinY);
    if (binNorm != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbins; iBinX++) {
        const Double_t binVal = _hResponseDiff -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = _hResponseDiff -> GetBinError(iBinX, iBinY);
        const Double_t newVal = (binVal / binNorm);
        const Double_t newErr = (binErr / binNorm);
        _hResponseDiff -> SetBinContent(iBinX, iBinY, newVal);
        _hResponseDiff -> SetBinError(iBinX, iBinY, newErr);
      }  // end x loop
    }
  }  // end y loop


  const UInt_t   iDetStart = _hSmearedDiff -> FindFirstBinAbove(0.);
  const UInt_t   iDetStop  = _hSmearedDiff -> FindLastBinAbove(0.);
  const Double_t intDet    = _hSmeared     -> Integral(iDetStart, iDetStop);
  const Double_t intDetMC  = hSmearNorm    -> Integral();
  const Double_t scaleS    = intDet / intDetMC;
  if (intDetMC > 0.) _hSmearedDiff -> Scale(scaleS);

  // indicate prior has been smoothed
  _priorSmoothed = true;

}  // end 'DoTwoExpoSmoothing()'


void StJetFolder::DoThreeExpoSmoothing() {

  const UInt_t NExp3(3);
  PrintInfo(18);

  // fit parameters [11 - 15 GeV]
  const TString  sFitHistExp("hFitPrior_onlyExpos");
  const TString  sFitHistTot("hFitPrior_totalFit");
  const TString  sTotName("fTotalFit");
  const TString  sTotFunc("(expo(0) + expo(2) + expo(4)) * (tanh((x-[6])/[7]) + tanh((x-[8])/[9]))");
  const TString  sExpNames[NExp3]          = {"fExpSoft", "fExpMed", "fExpHard"};
  const TString  sExpFuncs[NExp3]          = {"expo(0)", "expo(0)", "expo(0)"};
  const TString  sTanhNames[NTan][NParTan] = {{"p_{T}^{0}", "a_{0}"}, {"p_{T}^{1}", "a_{1}"}};
  const Double_t xExpRange[NExp3][NPts]    = {{0.6, 3.}, {3., 12.}, {12., 40.}};
  const Double_t pExpGuess[NExp3][NParTan] = {{0.19, -1.01}, {-2.44, -0.14}, {-1.06, -0.25}};
  const Double_t pTanhGuess[NTan][NParTan] = {{3., 1.}, {12., 1.}};
  const Double_t xTanhRange[NPts]          = {0.1, 57.};
  const UInt_t   fColExp[NExp3]            = {2, 3, 4};
  const UInt_t   fColTot(6);


  // for normalization
  TH1D *hPriorNorm = (TH1D*) _hPrior   -> Clone();
  TH1D *hSmearNorm = (TH1D*) _hSmeared -> Clone();
  hSmearNorm -> Reset("ICE");

  // initialize alternate prior and smeared
  _hPriorSmooth = (TH1D*) _hPrior   -> Clone();
  _hSmearedDiff = (TH1D*) _hSmeared -> Clone();
  _hPriorSmooth -> SetTitle("Smoothed prior");
  _hSmearedDiff -> SetTitle("Alternate smeared (via smoothed prior)");
  _hPriorSmooth -> Reset("ICE");
  _hSmearedDiff -> Reset("ICE");

  // initialize alternate response
  _hResponseDiff = (TH2D*) _hResponse   -> Clone();
  _hResponseDiff -> SetTitle("Alternate response matrix (via smoothed prior)");
  _hResponseDiff -> Reset("ICE");

  // determine prior bounds
  const UInt_t   iMin  = hPriorNorm -> FindFirstBinAbove(0.);
  const UInt_t   iMax  = hPriorNorm -> FindLastBinAbove(0.);
  const Double_t xMin  = hPriorNorm -> GetBinLowEdge(iMin);
  const Double_t xMax  = hPriorNorm -> GetBinLowEdge(iMax + 1);

  // create prior cutoff function
  if (_applyCutoff) {
    const Double_t xCutMin = TMath::Max(xMin, _pMin);
    const Double_t xCutMax = TMath::Min(xMax, _pMax);
    _fCutoff = new TF1("fCutoff", "tanh((x-[0])/[1])", xCutMin, xCutMax);
    _fCutoff -> SetParameter(0, _pTcutoff);
    _fCutoff -> SetParameter(1, _aCutoff);
  }


  // declare histograms for fitting
  TH1D *hFitExp = (TH1D*) _hPrior -> Clone();
  hFitExp -> SetName(sFitHistExp.Data());

  // initialize 1st set of exponentials
  TF1 *fExp[NExp3];
  for (UInt_t iExp = 0; iExp < NExp3; iExp++) {
    fExp[iExp] = new TF1(sExpNames[iExp].Data(), sExpFuncs[iExp].Data(), xExpRange[iExp][0], xExpRange[iExp][1]);
    fExp[iExp] -> SetParameter(0, pExpGuess[iExp][0]);
    fExp[iExp] -> SetParameter(1, pExpGuess[iExp][1]);
    fExp[iExp] -> SetLineColor(fColExp[iExp]);
  }

  // determine exponential parameters
  Double_t pExpFit[NExp3][NParTan];
  for (UInt_t iExp = 0; iExp < NExp3; iExp++) {
    if (iExp == 0) {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R");
    }
    else {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R+");
    }
    pExpFit[iExp][0] = fExp[iExp] -> GetParameter(0);
    pExpFit[iExp][1] = fExp[iExp] -> GetParameter(1);
  }

  // save exponent fits to output file
  _fOut   -> cd();
  hFitExp -> Write();

  // initialize total fit function
  _fSmoothPrior = new TF1(sTotName.Data(), sTotFunc.Data(), xTanhRange[0], xTanhRange[1]);
  _fSmoothPrior -> SetLineColor(fColTot);

  const UInt_t iTanMin = 2 * NExp3;
  for (UInt_t iExp = 0; iExp < NExp3; iExp++) {
    const UInt_t  iParExp = 2 * iExp;
    _fSmoothPrior -> SetParameter(iParExp + 0, pExpFit[iExp][0]);
    _fSmoothPrior -> FixParameter(iParExp + 1, pExpFit[iExp][1]);
  }
  for (UInt_t iTan = 0; iTan < NTan; iTan++) {
    const UInt_t iParTan = iTanMin + (2 * iTan);
    _fSmoothPrior -> SetParameter(iParTan + 0, pTanhGuess[iTan][0]);
    _fSmoothPrior -> SetParameter(iParTan + 1, pTanhGuess[iTan][1]);
    _fSmoothPrior -> SetParName(iParTan + 0, sTanhNames[iTan][0].Data());
    _fSmoothPrior -> SetParName(iParTan + 1, sTanhNames[iTan][1].Data());
  }

  // fit function
  Int_t fitStat(-1);
  fitStat = _hPrior -> Fit(sTotName.Data(), "B", "", xExpRange[0][0], xExpRange[NExp3 - 1][1]);
  fitStat = _hPrior -> Fit(sTotName.Data(), "B", "", xExpRange[0][0], xExpRange[NExp3 - 1][1]);
  if (fitStat != 0) {
    PrintError(14);
    _fSmoothPrior -> SetParameter(0, 3.81);
    _fSmoothPrior -> SetParameter(2, 1.03);
    _fSmoothPrior -> SetParameter(4, -14.53);
    _fSmoothPrior -> FixParameter(6, 7.14);
    _fSmoothPrior -> FixParameter(7, 13.29);
    _fSmoothPrior -> FixParameter(8, 7.33);
    _fSmoothPrior -> FixParameter(9, -13.34);
    _hPrior       -> Fit(sTotName.Data(), "RB");
  }

  // fill prior
  _hPriorSmooth -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t t = _fSmoothPrior -> GetRandom(xMin, xMax);
    _hPriorSmooth -> Fill(t);
  }

  // apply cutoff if necessary and trim prior
  if (_applyCutoff) {
    _hPriorSmooth -> Multiply(_fCutoff);
  }
  ConstrainPriorRange(1);


  // normalize smoothed prior
  const Double_t xParStart = TMath::Max(xMin, xTanhRange[0]);
  const Double_t xParStop  = TMath::Min(xMax, xTanhRange[1]);
  const Double_t intPar    = _fSmoothPrior -> Integral(xParStart, xParStop);
  const Double_t intParMC  = _hPriorSmooth -> Integral();
  const Double_t scaleP    = intPar / intParMC;
  if (intParMC > 0.) _hPriorSmooth -> Scale(scaleP);

  // normalize by bin width
  const UInt_t nBinsP = _hPriorSmooth -> GetNbinsX();
  for (UInt_t iBinP = 1; iBinP < (nBinsP + 1); iBinP++) {
    const Double_t pVal = _hPriorSmooth -> GetBinContent(iBinP);
    const Double_t pErr = _hPriorSmooth -> GetBinError(iBinP);
    const Double_t pBin = _hPriorSmooth -> GetBinWidth(iBinP);
    _hPriorSmooth -> SetBinContent(iBinP, pVal / pBin);
    _hPriorSmooth -> SetBinError(iBinP, pErr / pBin);
  }

  // only need smoothed prior if smoothing response
  if (!_smoothPrior && _smoothResponse) return;


  // create histogram to sample from and apply efficiency
  TH1D *hAltPriorSample = (TH1D*) _hPriorSmooth -> Clone();
  if (_smoothEfficiency)
    hAltPriorSample -> Multiply(_fEfficiency);
  else
    hAltPriorSample -> Multiply(_hEfficiency);

  // create detector-level prior and response
  _hSmearedDiff -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t p = hAltPriorSample -> GetRandom();
    const Double_t s = Smear(p);
    // quick fix [02.04.2022]
    //if (s > 2.) {
    if (s > -1000.) {
      const UInt_t   i = hPriorNorm    -> GetBin(p);
      const Double_t o = hPriorNorm    -> GetBinContent(i);
      const Double_t n = _hPriorSmooth -> GetBinContent(i);
      const Double_t w = o / n;
      _hSmearedDiff -> Fill(s);
      if ((o == 0.) || (n == 0.))
        _hResponseDiff -> Fill(s, p, 1.);
      else
        _hResponseDiff -> Fill(s, p, w);
    }
    hSmearNorm -> Fill(p);
  }


  // normalize response
  const UInt_t nXbins = _hResponseDiff -> GetNbinsX();
  const UInt_t nYbins = _hResponseDiff -> GetNbinsY();
  for (UInt_t iBinY = 0; iBinY < nYbins; iBinY++) {
    const Double_t binNorm = _hResponseDiff -> Integral(1, nXbins, iBinY, iBinY);
    if (binNorm != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbins; iBinX++) {
        const Double_t binVal = _hResponseDiff -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = _hResponseDiff -> GetBinError(iBinX, iBinY);
        const Double_t newVal = (binVal / binNorm);
        const Double_t newErr = (binErr / binNorm);
        _hResponseDiff -> SetBinContent(iBinX, iBinY, newVal);
        _hResponseDiff -> SetBinError(iBinX, iBinY, newErr);
      }  // end x loop
    }
  }  // end y loop


  const UInt_t   iDetStart = _hSmearedDiff -> FindFirstBinAbove(0.);
  const UInt_t   iDetStop  = _hSmearedDiff -> FindLastBinAbove(0.);
  const Double_t intDet    = _hSmeared     -> Integral(iDetStart, iDetStop);
  const Double_t intDetMC  = hSmearNorm    -> Integral();
  const Double_t scaleS    = intDet / intDetMC;
  if (intDetMC > 0.) _hSmearedDiff -> Scale(scaleS);

  // indicate prior has been smoothed
  _priorSmoothed = true;

}  // end 'DoThreeExpoSmoothing()'


void StJetFolder::DoFourExpoSmoothing() {

  const UInt_t NExp4(4);
  PrintInfo(19);

  // fit parameters [15 - 20 GeV]
  const TString  sFitHistExp("hFitPrior_onlyExpos");
  const TString  sFitHistTot("hFitPrior_totalFit");
  const TString  sTotName("fTotalFit");
  const TString  sTotFunc("(expo(0) + expo(2) + expo(4) + expo(6)) * (tanh((x-[8])/[9]) + tanh((x-[10])/[11]))");
  const TString  sExpNames[NExp4]          = {"fExpSoft", "fExpMedSoft", "fExpMedHard", "fExpHard"};
  const TString  sExpFuncs[NExp4]          = {"expo(0)", "expo(0)", "expo(0)", "expo(0)"};
  const TString  sTanhNames[NTan][NParTan] = {{"p_{T}^{0}", "a_{0}"}, {"p_{T}^{1}", "a_{1}"}};
  const Double_t xExpRange[NExp4][NPts]    = {{0.2, 2.}, {2., 5.}, {5., 14.5}, {14.5, 40.}};
  const Double_t pExpGuess[NExp4][NParTan] = {{0.31, -1.18}, {-1.36, -0.39}, {-3.81, -0.01}, {-0.49, -0.22}};
  const Double_t pTanhGuess[NTan][NParTan] = {{7.5, 13.}, {7.5, -13.}};
  const Double_t xTanhRange[NPts]          = {0.1, 57.};
  const UInt_t   fColExp[NExp4]            = {2, 3, 4, 6};
  const UInt_t   fColTot(6);


  // for normalization
  TH1D *hPriorNorm = (TH1D*) _hPrior   -> Clone();
  TH1D *hSmearNorm = (TH1D*) _hSmeared -> Clone();
  hSmearNorm -> Reset("ICE");

  // initialize alternate prior and smeared
  _hPriorSmooth = (TH1D*) _hPrior   -> Clone();
  _hSmearedDiff = (TH1D*) _hSmeared -> Clone();
  _hPriorSmooth -> SetTitle("Smoothed prior");
  _hSmearedDiff -> SetTitle("Alternate smeared (via smoothed prior)");
  _hPriorSmooth -> Reset("ICE");
  _hSmearedDiff -> Reset("ICE");

  // initialize alternate response
  _hResponseDiff = (TH2D*) _hResponse   -> Clone();
  _hResponseDiff -> SetTitle("Alternate response matrix (via smoothed prior)");
  _hResponseDiff -> Reset("ICE");

  // determine prior bounds
  const UInt_t   iMin  = hPriorNorm -> FindFirstBinAbove(0.);
  const UInt_t   iMax  = hPriorNorm -> FindLastBinAbove(0.);
  const Double_t xMin  = hPriorNorm -> GetBinLowEdge(iMin);
  const Double_t xMax  = hPriorNorm -> GetBinLowEdge(iMax + 1);

  // create prior cutoff function
  if (_applyCutoff) {
    const Double_t xCutMin = TMath::Max(xMin, _pMin);
    const Double_t xCutMax = TMath::Min(xMax, _pMax);
    _fCutoff = new TF1("fCutoff", "tanh((x-[0])/[1])", xCutMin, xCutMax);
    _fCutoff -> SetParameter(0, _pTcutoff);
    _fCutoff -> SetParameter(1, _aCutoff);
  }


  // declare histograms for fitting
  TH1D *hFitExp = (TH1D*) _hPrior -> Clone();
  hFitExp -> SetName(sFitHistExp.Data());

  // initialize 1st set of exponentials
  TF1 *fExp[NExp4];
  for (UInt_t iExp = 0; iExp < NExp4; iExp++) {
    fExp[iExp] = new TF1(sExpNames[iExp].Data(), sExpFuncs[iExp].Data(), xExpRange[iExp][0], xExpRange[iExp][1]);
    fExp[iExp] -> SetParameter(0, pExpGuess[iExp][0]);
    fExp[iExp] -> SetParameter(1, pExpGuess[iExp][1]);
    fExp[iExp] -> SetLineColor(fColExp[iExp]);
  }

  // determine exponential parameters
  Double_t pExpFit[NExp4][NParTan];
  for (UInt_t iExp = 0; iExp < NExp4; iExp++) {
    if (iExp == 0) {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R");
    }
    else {
      hFitExp -> Fit(sExpNames[iExp].Data(), "R+");
    }
    pExpFit[iExp][0] = fExp[iExp] -> GetParameter(0);
    pExpFit[iExp][1] = fExp[iExp] -> GetParameter(1);
  }

  // save exponent fits to output file
  _fOut   -> cd();
  hFitExp -> Write();

  // initialize total fit function
  _fSmoothPrior = new TF1(sTotName.Data(), sTotFunc.Data(), xTanhRange[0], xTanhRange[1]);
  _fSmoothPrior -> SetLineColor(fColTot);

  const UInt_t iTanMin = 2 * NExp4;
  for (UInt_t iExp = 0; iExp < NExp4; iExp++) {
    const UInt_t  iParExp = 2 * iExp;
    _fSmoothPrior -> SetParameter(iParExp + 0, pExpFit[iExp][0]);
    _fSmoothPrior -> FixParameter(iParExp + 1, pExpFit[iExp][1]);
  }
  for (UInt_t iTan = 0; iTan < NTan; iTan++) {
    const UInt_t iParTan = iTanMin + (2 * iTan);
    _fSmoothPrior -> SetParameter(iParTan + 0, pTanhGuess[iTan][0]);
    _fSmoothPrior -> SetParameter(iParTan + 1, pTanhGuess[iTan][1]);
    _fSmoothPrior -> SetParName(iParTan + 0, sTanhNames[iTan][0].Data());
    _fSmoothPrior -> SetParName(iParTan + 1, sTanhNames[iTan][1].Data());
  }

  // fit function
  _hPrior -> Fit(sTotName.Data(), "B", "", xExpRange[0][0], xExpRange[NExp4 - 1][1]);

  // fill prior
  _hPriorSmooth -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t t = _fSmoothPrior -> GetRandom(xMin, xMax);
    _hPriorSmooth -> Fill(t);
  }

  // apply cutoff if necessary and trim prior
  if (_applyCutoff) {
    _hPriorSmooth -> Multiply(_fCutoff);
  }
  ConstrainPriorRange(1);


  // normalize smoothed prior
  const Double_t xParStart = TMath::Max(xMin, xTanhRange[0]);
  const Double_t xParStop  = TMath::Min(xMax, xTanhRange[1]);
  const Double_t intPar    = _fSmoothPrior -> Integral(xParStart, xParStop);
  const Double_t intParMC  = _hPriorSmooth -> Integral();
  const Double_t scaleP    = intPar / intParMC;
  if (intParMC > 0.) _hPriorSmooth -> Scale(scaleP);

  // normalize by bin width
  const UInt_t nBinsP = _hPriorSmooth -> GetNbinsX();
  for (UInt_t iBinP = 1; iBinP < (nBinsP + 1); iBinP++) {
    const Double_t pVal = _hPriorSmooth -> GetBinContent(iBinP);
    const Double_t pErr = _hPriorSmooth -> GetBinError(iBinP);
    const Double_t pBin = _hPriorSmooth -> GetBinWidth(iBinP);
    _hPriorSmooth -> SetBinContent(iBinP, pVal / pBin);
    _hPriorSmooth -> SetBinError(iBinP, pErr / pBin);
  }

  // only need smoothed prior if smoothing response
  if (!_smoothPrior && _smoothResponse) return;


  // create histogram to sample from and apply efficiency
  TH1D *hAltPriorSample = (TH1D*) _hPriorSmooth -> Clone();
  if (_smoothEfficiency)
    hAltPriorSample -> Multiply(_fEfficiency);
  else
    hAltPriorSample -> Multiply(_hEfficiency);

  // create detector-level prior and response
  _hSmearedDiff -> Reset("ICE");
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {
    const Double_t p = hAltPriorSample -> GetRandom();
    const Double_t s = Smear(p);
    // quick fix [02.04.2022]
    //if (s > 2.) {
    if (s > -1000.) {
      const UInt_t   i = hPriorNorm    -> GetBin(p);
      const Double_t o = hPriorNorm    -> GetBinContent(i);
      const Double_t n = _hPriorSmooth -> GetBinContent(i);
      const Double_t w = o / n;
      _hSmearedDiff -> Fill(s);
      if ((o == 0.) || (n == 0.))
        _hResponseDiff -> Fill(s, p, 1.);
      else
        _hResponseDiff -> Fill(s, p, w);
    }
    hSmearNorm -> Fill(p);
  }


  // normalize response
  const UInt_t nXbins = _hResponseDiff -> GetNbinsX();
  const UInt_t nYbins = _hResponseDiff -> GetNbinsY();
  for (UInt_t iBinY = 0; iBinY < nYbins; iBinY++) {
    const Double_t binNorm = _hResponseDiff -> Integral(1, nXbins, iBinY, iBinY);
    if (binNorm != 0.) {
      for (UInt_t iBinX = 0; iBinX < nXbins; iBinX++) {
        const Double_t binVal = _hResponseDiff -> GetBinContent(iBinX, iBinY);
        const Double_t binErr = _hResponseDiff -> GetBinError(iBinX, iBinY);
        const Double_t newVal = (binVal / binNorm);
        const Double_t newErr = (binErr / binNorm);
        _hResponseDiff -> SetBinContent(iBinX, iBinY, newVal);
        _hResponseDiff -> SetBinError(iBinX, iBinY, newErr);
      }  // end x loop
    }
  }  // end y loop


  const UInt_t   iDetStart = _hSmearedDiff -> FindFirstBinAbove(0.);
  const UInt_t   iDetStop  = _hSmearedDiff -> FindLastBinAbove(0.);
  const Double_t intDet    = _hSmeared     -> Integral(iDetStart, iDetStop);
  const Double_t intDetMC  = hSmearNorm    -> Integral();
  const Double_t scaleS    = intDet / intDetMC;
  if (intDetMC > 0.) _hSmearedDiff -> Scale(scaleS);

  // indicate prior has been smoothed
  _priorSmoothed = true;

}  // end 'DoFourExpoSmoothing()'


void StJetFolder::ConstrainPriorRange(const UInt_t code) {

  UInt_t   nPbin(0);
  Bool_t   isOutOfRange(false);
  Bool_t   isNegative(false);
  Double_t pVal(0.);
  Double_t pEdgeLo(0.);
  Double_t pEdgeHi(0.);

  // zero-out negative bins and bins outside set range 
  switch (code) {
    case 0:
      nPbin = _hPriorDiff -> GetNbinsX();
      for (UInt_t iPbin = 1; iPbin < (nPbin + 1); iPbin++) {
        pVal         = _hPriorDiff -> GetBinContent(iPbin);
        pEdgeLo      = _hPriorDiff -> GetBinLowEdge(iPbin);
        pEdgeHi      = _hPriorDiff -> GetBinLowEdge(iPbin + 1);
        isOutOfRange = ((pEdgeLo < _pMin) || (pEdgeHi > _pMax));
        isNegative   = (pVal < 0.);
        if (isOutOfRange || isNegative) {
          _hPriorDiff -> SetBinContent(iPbin, 0.);
          _hPriorDiff -> SetBinError(iPbin, 0.);
        }
      }  // end pT bin loop
      break;
    case 1:
      nPbin = _hPriorSmooth -> GetNbinsX();
      for (UInt_t iPbin = 1; iPbin < (nPbin + 1); iPbin++) {
        pVal         = _hPriorSmooth -> GetBinContent(iPbin);
        pEdgeLo      = _hPriorSmooth -> GetBinLowEdge(iPbin);
        pEdgeHi      = _hPriorSmooth -> GetBinLowEdge(iPbin + 1);
        isOutOfRange = ((pEdgeLo < _pMin) || (pEdgeHi > _pMax));
        isNegative   = (pVal < 0.);
        if (isOutOfRange || isNegative) {
          _hPriorSmooth -> SetBinContent(iPbin, 0.);
          _hPriorSmooth -> SetBinError(iPbin, 0.);
        }
      }  // end pT bin loop
      break;
  }  // end switch-case

}  // end 'ConstrainPriorRange()'


void StJetFolder::SmoothResponse() {

  PrintInfo(20);

  // response smoothing parameters
  const TString sFunc("gaus(0) + gaus(3)");
  const TString sQtNames[NPtPar]          = {"fQt_pt0206", "fQt_pt061", "fQt_pt12", "fQt_pt210", "fQt_pt1057"};
  const Float_t xFuncDef[NRange]          = {0., 2.};
  const Float_t xQtRange[NRange]          = {0.5, 1.3};
  const Float_t xPtParLo[NPtPar]          = {0.2, 0.6, 1., 2., 10.};
  const Float_t xPtParHi[NPtPar]          = {0.6, 1., 2., 10., 57.};
  const Float_t fQt911[NPtPar][NParGaus]  = {{0.59, 1.01, 0.06, 0.007, 0.62, 0.06}, {0.68, 1.02, 0.06, 0.04, 0.67, 0.06}, {0.57, 1.02, 0.07, 0.06, 0.68, 0.06}, {0.47, 1.01, 0.07, 0.15, 0.73, 0.08}, {0.48, 0.99, 0.07, 0.12, 0.71, 0.08}};
  const Float_t fQt1115[NPtPar][NParGaus] = {{0.66, 1.00, 0.06, 0.003, 0.60, 0.06}, {0.51, 1.00, 0.06, 0.03, 0.67, 0.06}, {0.53, 1.01, 0.07, 0.12, 0.68, 0.06}, {0.52, 1.00, 0.07, 0.09, 0.71, 0.08}, {0.47, 0.97, 0.07, 0.11, 0.69, 0.08}};
  const Float_t fQt1520[NPtPar][NParGaus] = {{0.53, 1.00, 0.06, 0.003, 0.64, 0.06}, {0.56, 1.00, 0.06, 0.002, 0.60, 0.06}, {0.57, 1.00, 0.07, 0.04, 0.68, 0.08}, {0.48, 1.00, 0.07, 0.14, 0.70, 0.07}, {0.47, 0.97, 0.07, 0.13, 0.68, 0.08}};

  // prior sampling parameters
  const TString sParSampleName("fParForResSmooth");
  const TString sParSample911("(expo(0) + expo(2)) * tanh((x-[4]/[5])) * ([6] * (1 - TMath::Exp(-1.*[7]*x)))");
  const TString sParSample1115("(expo(0) + expo(2) + expo(4)) * (tanh((x-[6])/[7]) + tanh((x-[8])/[9])) * ([10] * (1 - TMath::Exp(-1.*[11]*x)))");
  const TString sParSample1520("(expo(0) + expo(2) + expo(4) + expo(6)) * (tanh((x-[8])/[9]) + tanh((x-[10])/[11])) * ([12] * (1 - TMath::Exp(-1.*[13]*x)))");
  const Double_t xParSampleRange[NPts] = {0.1, 57.};


  // select qT fit parameters
  Float_t fFuncPar[NPtPar][NParGaus];
  switch (_eTtrgIndex) {
    case 0:
      for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
        for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
          fFuncPar[iPtPar][iParGaus] = fQt911[iPtPar][iParGaus];
        }
      }
      break;
    case 1:
      for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
        for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
          fFuncPar[iPtPar][iParGaus] = fQt1115[iPtPar][iParGaus];
        }
      }
      break;
    case 2:
      for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
        for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
          fFuncPar[iPtPar][iParGaus] = fQt1520[iPtPar][iParGaus];
        }
      }
      break;
    default:
      for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
        for (UInt_t iParGaus = 0; iParGaus < NParGaus; iParGaus++) {
          fFuncPar[iPtPar][iParGaus] = fQt911[iPtPar][iParGaus];
        }
      }
      break;
  }  // end switch-case


  // create smoothing functions
  for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
    _fQtSmooth[iPtPar] = new TF1(sQtNames[iPtPar].Data(), sFunc.Data(), xFuncDef[0], xFuncDef[1]);
    for (UInt_t iParam = 0; iParam < NParGaus; iParam++) {
      _fQtSmooth[iPtPar] -> SetParameter(iParam, fFuncPar[iPtPar][iParam]);
    }
  }  // end pTPar bin loop

  // determine fit range
  Double_t xUseRange[NRange] = {0., 0.};
  if (_applyCutoff) {
    xUseRange[0] = TMath::Max(_pTcutoff, xParSampleRange[0]);
  } else {
    xUseRange[0] = xParSampleRange[0];
  }
  xUseRange[1] = xParSampleRange[1];

  // create function to sample from
  switch (_eTtrgIndex) {
    case 0:
      _fParToSample = new TF1(sParSampleName.Data(), sParSample911.Data(), xUseRange[0], xUseRange[1]);
      _fParToSample -> FixParameter(0, _fSmoothPrior -> GetParameter(0));
      _fParToSample -> FixParameter(1, _fSmoothPrior -> GetParameter(1));
      _fParToSample -> FixParameter(2, _fSmoothPrior -> GetParameter(2));
      _fParToSample -> FixParameter(3, _fSmoothPrior -> GetParameter(3));
      _fParToSample -> FixParameter(4, _fSmoothPrior -> GetParameter(4));
      _fParToSample -> FixParameter(5, _fSmoothPrior -> GetParameter(5));
      _fParToSample -> FixParameter(6, _fEfficiency  -> GetParameter(0));
      _fParToSample -> FixParameter(7, _fEfficiency  -> GetParameter(1));
      break;
    case 1:
      _fParToSample = new TF1(sParSampleName.Data(), sParSample1115.Data(), xUseRange[0], xUseRange[1]);
      _fParToSample -> FixParameter(0,  _fSmoothPrior -> GetParameter(0));
      _fParToSample -> FixParameter(1,  _fSmoothPrior -> GetParameter(1));
      _fParToSample -> FixParameter(2,  _fSmoothPrior -> GetParameter(2));
      _fParToSample -> FixParameter(3,  _fSmoothPrior -> GetParameter(3));
      _fParToSample -> FixParameter(4,  _fSmoothPrior -> GetParameter(4));
      _fParToSample -> FixParameter(5,  _fSmoothPrior -> GetParameter(5));
      _fParToSample -> FixParameter(6,  _fSmoothPrior -> GetParameter(6));
      _fParToSample -> FixParameter(7,  _fSmoothPrior -> GetParameter(7));
      _fParToSample -> FixParameter(8,  _fSmoothPrior -> GetParameter(8));
      _fParToSample -> FixParameter(9,  _fSmoothPrior -> GetParameter(9));
      _fParToSample -> FixParameter(10, _fEfficiency  -> GetParameter(0));
      _fParToSample -> FixParameter(11, _fEfficiency  -> GetParameter(1));
      break;
    case 2:
      _fParToSample = new TF1(sParSampleName.Data(), sParSample1520.Data(), xUseRange[0], xUseRange[1]);
      _fParToSample -> FixParameter(0,  _fSmoothPrior -> GetParameter(0));
      _fParToSample -> FixParameter(1,  _fSmoothPrior -> GetParameter(1));
      _fParToSample -> FixParameter(2,  _fSmoothPrior -> GetParameter(2));
      _fParToSample -> FixParameter(3,  _fSmoothPrior -> GetParameter(3));
      _fParToSample -> FixParameter(4,  _fSmoothPrior -> GetParameter(4));
      _fParToSample -> FixParameter(5,  _fSmoothPrior -> GetParameter(5));
      _fParToSample -> FixParameter(6,  _fSmoothPrior -> GetParameter(6));
      _fParToSample -> FixParameter(7,  _fSmoothPrior -> GetParameter(7));
      _fParToSample -> FixParameter(8,  _fSmoothPrior -> GetParameter(8));
      _fParToSample -> FixParameter(9,  _fSmoothPrior -> GetParameter(9));
      _fParToSample -> FixParameter(10, _fSmoothPrior -> GetParameter(10));
      _fParToSample -> FixParameter(11, _fSmoothPrior -> GetParameter(11));
      _fParToSample -> FixParameter(12, _fEfficiency  -> GetParameter(0));
      _fParToSample -> FixParameter(13, _fEfficiency  -> GetParameter(1));
      break;
    default:
      _fParToSample = new TF1(sParSampleName.Data(), sParSample911.Data(), xParSampleRange[0], xParSampleRange[1]);
      _fParToSample -> FixParameter(0, _fSmoothPrior -> GetParameter(0));
      _fParToSample -> FixParameter(1, _fSmoothPrior -> GetParameter(1));
      _fParToSample -> FixParameter(2, _fSmoothPrior -> GetParameter(2));
      _fParToSample -> FixParameter(3, _fSmoothPrior -> GetParameter(3));
      _fParToSample -> FixParameter(4, _fSmoothPrior -> GetParameter(4));
      _fParToSample -> FixParameter(5, _fSmoothPrior -> GetParameter(5));
      _fParToSample -> FixParameter(6, _fEfficiency  -> GetParameter(0));
      _fParToSample -> FixParameter(7, _fEfficiency  -> GetParameter(1));
      break;
  }  // end switch-case


  // initialize smoothed histograms
  _hSmearedSmooth  = (TH1D*) _hSmeared  -> Clone();
  _hResponseSmooth = (TH2D*) _hResponse -> Clone();
  _hSmearedSmooth  -> SetTitle("Smoothed smeared");
  _hResponseSmooth -> SetTitle("Smoothed response matrix");
  _hSmearedSmooth  -> Reset("ICE");
  _hResponseSmooth -> Reset("ICE");


  // mc loop
  for (Int_t iMC = 0; iMC < _nMcSmooth; iMC++) {

    // sample from particle distribution
    Bool_t   isAboveMin(false);
    Double_t pTpar(0.);
    do {
      pTpar      = _fParToSample -> GetRandom();
      isAboveMin = (pTpar > 0.2);
    } while (!isAboveMin);

    // determine pTpar bin and fill histograms
    const UInt_t iPtBig = NPtPar - 1;
    for (UInt_t iPtPar = 0; iPtPar < NPtPar; iPtPar++) {
      if ((pTpar >= xPtParLo[iPtPar]) && (pTpar < xPtParHi[iPtPar])) {
        const Double_t qTjet = _fQtSmooth[iPtPar] -> GetRandom(xQtRange[0], xQtRange[1]);
        const Double_t pTdet = pTpar * qTjet;
        // quick fix [02.04.2022]
        //if (pTdet > 2.) {
        //  _hResponseSmooth -> Fill(pTdet, pTpar);
        //  _hSmearedSmooth  -> Fill(pTdet);
        //}
        _hResponseSmooth -> Fill(pTdet, pTpar);
        _hSmearedSmooth  -> Fill(pTdet);
        break;
      }
      else if (pTpar > xPtParHi[iPtBig]) {
        const Double_t qTjet = _fQtSmooth[iPtBig] -> GetRandom(xQtRange[0], xQtRange[1]);
        const Double_t pTdet = pTpar * qTjet;
        // quick fix [02.04.2022]
        //if (pTdet > 2.) {
        //  _hResponseSmooth -> Fill(pTdet, pTpar);
        //  _hSmearedSmooth  -> Fill(pTdet);
        //}
        _hResponseSmooth -> Fill(pTdet, pTpar);
        _hSmearedSmooth  -> Fill(pTdet);
        break;
      } 
    }  // end pTpar bin loop
  }  // end mc loop


  // normalize response
  const UInt_t nResX = _hResponseSmooth -> GetNbinsX();
  const UInt_t nResY = _hResponseSmooth -> GetNbinsY();
  for (UInt_t iResY = 1; iResY < (nResY + 1); iResY++) {
    const Double_t intResY = _hResponseSmooth -> Integral(1, nResX, iResY, iResY);
    for (UInt_t iResX = 1; iResX < (nResX + 1); iResX++) {
      const Double_t resBin  = _hResponseSmooth -> GetBinContent(iResX, iResY);
      const Double_t resErr  = _hResponseSmooth -> GetBinError(iResX, iResY);
      if (intResY > 0.) {
        _hResponseSmooth -> SetBinContent(iResX, iResY, resBin / intResY);
        _hResponseSmooth -> SetBinError(iResX, iResY, resErr / intResY);
      }
    }  // ned pTdet bin loop
  }  // end pTpar bin loop

  // normalize smoothed smeared hist
  const UInt_t nPtX = _hSmearedSmooth -> GetNbinsX();
  for (UInt_t iPt = 1; iPt < (nPtX + 1); iPt++) {
    const Double_t binDet   = _hSmearedSmooth -> GetBinContent(iPt);
    const Double_t errDet   = _hSmearedSmooth -> GetBinError(iPt);
    const Double_t binWidth = _hSmearedSmooth -> GetBinWidth(iPt);
    _hSmearedSmooth -> SetBinContent(iPt, binDet / binWidth);
    _hSmearedSmooth -> SetBinError(iPt, errDet / binWidth);
  }

  const UInt_t   iDetStart = _hSmearedSmooth -> FindFirstBinAbove(0.);
  const UInt_t   iDetStop  = _hSmearedSmooth -> FindLastBinAbove(0.);
  const Double_t intDet    = _hSmeared       -> Integral(iDetStart, iDetStop);
  const Double_t intDetMC  = _hSmearedSmooth -> Integral();
  const Double_t scaleS    = intDet / intDetMC;
  if (intDetMC > 0.) _hSmearedDiff -> Scale(scaleS);

  // indicate that response is smoothed
  _responseSmoothed = true;

}  // end 'SmoothResponse()'


Bool_t StJetFolder::CheckFlags() {

  // check spectra
  Bool_t spectraOK = true;
  for (Int_t i = 0; i < 5; i++) {
    if (!_flag[i]) {
      spectraOK = false;
      PrintError(8);
      break;
    }
  }
  if (spectraOK)
    PrintInfo(1);

  // check info
  Bool_t infoOK = true;
  for (Int_t i = 5; i < 8; i++) {
    if (!_flag[i]) {
      infoOK = false;
      PrintError(9);
      break;
    }
  }
  if (infoOK)
    PrintInfo(2);

  // check parameters
  Bool_t parametersOK = true;
  for (Int_t i = 8; i < 10; i++) {
    if (!_flag[i]) {
      parametersOK = false;
      PrintError(10);
      break;
    }
  }
  if (parametersOK)
    PrintInfo(3);


  Bool_t inputOK = (spectraOK && infoOK && parametersOK);
  return inputOK;

}  // end 'CheckFlags()'

// End ------------------------------------------------------------------------
