// 'StJetFolder.h'
// Derek Anderson
// 02.14.2017
//
// This class handles the unfolding of a provided spectrum.  Please see below
// for definition of names:
//
//   prior    -- prior spectrum, used to seed unfolding algorithm.
//   smear    -- smeared prior spectrum, the prior spectrum with
//               efficiencies, smearing, etc. applied.
//   measure  -- measured spectrum, to be unfolded.
//   unfold   -- unfolded measured spectrum.
//   backfold -- unfolded spectrum with efficiencies, smearing,
//               etc. applied.
//
// Pearson Coefficient calculation adapted from Rhagav K. Elayavalli.
//
// Last updated: 06.18.2021


#ifndef StJetFolder_h
#define StJetFolder_h

#include <cmath>
#include <cassert>
#include <iostream>
// ROOT includes
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TPaveText.h"
#include "TSVDUnfold.h"
// RooUnfold includes
#include "../RooUnfold/RooUnfoldResponse.h"
#include "../RooUnfold/RooUnfoldBayes.h"
#include "../RooUnfold/RooUnfoldSvd.h"
#include "../RooUnfold/RooUnfoldBinByBin.h"
#include "../RooUnfold/RooUnfoldTUnfold.h"
#include "../RooUnfold/RooUnfoldInvert.h"
#include "../RooUnfold/RooUnfoldErrors.h"

using namespace std;


// global system constants
const Int_t    Nflag = 11;
const Bool_t   Debug = true;
const Double_t Mpion = 0.140;

// global unfolding constants
const Int_t    NBackDef = 10000;
const Double_t UdefMax  = 100.;
const Double_t BdefMax  = 100.;

// global prior constants
const Double_t PdefMin     = -100.;
const Double_t PdefMax     = 100.;
const Double_t PtCutoffDef = 3.;
const Double_t AcutoffDef  = 1.;  

// global smoothing constants
const UInt_t NTan         = 2;
const UInt_t NPts         = 2;
const UInt_t NRange       = 2;
const UInt_t NPtPar       = 5;
const UInt_t NParEff      = 2;
const UInt_t NParTan      = 2;
const UInt_t NParPol      = 3;
const UInt_t NParGaus     = 6;
const UInt_t NSmoothDef   = 5000000;
const Bool_t PriSmoothDef = false;
const Bool_t PriCutoffDef = false;



class StJetFolder {

public:

  StJetFolder(const Char_t *oFile, const Bool_t pearDebug=Debug, const Int_t nMcSmooth=NSmoothDef);
  virtual ~StJetFolder();

  // public methods ('StJetFolder.io.h')
  void SetPrior(const Char_t *pFile, const Char_t *pName, const Bool_t doPriorSmooth=PriSmoothDef, const Bool_t doPriorCutoff=PriCutoffDef);
  void SetSmeared(const Char_t *sFile, const Char_t *sName);
  void SetMeasured(const Char_t *mFile, const Char_t *mName);
  void SetResponse(const Char_t *rFile, const Char_t *rName, const Bool_t doResSmooth);
  void SetEfficiency(const Char_t *eFile, const Char_t *eName, const Bool_t doEffSmooth, const Bool_t doEffFineTune, const Bool_t removeErrors);
  void SetEventInfo(const Int_t beam, const Double_t energy);
  void SetTriggerInfo(const Int_t trigger, const UInt_t eTtrgIndex, const Double_t hMax);
  void SetJetInfo(const Int_t type, const Int_t nRM, const Double_t rJet, const Double_t aMin, const Double_t pTmin);
  void SetPriorCutoffParameters(const Double_t pTcut=PtCutoffDef, const Double_t aCut=AcutoffDef);
  void SetPriorParameters(const Int_t prior, const Double_t bPrior, const Double_t mPrior, const Double_t nPrior, const Double_t tPrior, const Double_t pMin=PdefMin, const Double_t pMax=PdefMax);
  void SetUnfoldParameters(const Int_t method, const Int_t kReg, const Int_t nToy, const Double_t uMax=UdefMax, const Double_t bMax=BdefMax);
  // public methods ('StJetFolder.cxx')
  void Init();
  void Unfold(Double_t &chi2unfold);
  void Backfold(Double_t &chi2backfold, const Int_t nMcBack=NBackDef);
  void Finish();

  // static public methods ('StJetFolder.math.h')
  static Double_t Levy(const Double_t *x, const Double_t *p);
  static Double_t Tsallis(const Double_t *x, const Double_t *p);
  static Double_t Exponential(const Double_t *x, const Double_t *p);
  static Double_t PowerLaw(const Double_t *x, const Double_t *p);
  static Double_t Landau(const Double_t *x, const Double_t *p);
  static Double_t FineTuneEff(const Double_t *x, const Double_t *p);


private:

  // atomic members
  Int_t     _trigger;
  Int_t     _type;
  Int_t     _prior;
  Int_t     _method;
  Int_t     _kReg;
  Int_t     _nMcSmooth;
  Int_t     _nMcBack;
  Int_t     _nToy;
  UInt_t    _eTtrgIndex;
  Bool_t    _smoothResponse;
  Bool_t    _responseSmoothed;
  Bool_t    _smoothEfficiency;
  Bool_t    _fineTuneEfficiency;
  Bool_t    _removeErrors;
  Bool_t    _smoothPrior;
  Bool_t    _priorSmoothed;
  Bool_t    _applyCutoff;
  Bool_t    _differentPrior;
  Bool_t    _diffPriorCreated;
  Bool_t    _pearsonDebug;
  Bool_t    _flag[Nflag];
  Double_t  _pTcutoff;
  Double_t  _aCutoff;
  Double_t  _bPrior;
  Double_t  _mPrior;
  Double_t  _nPrior;
  Double_t  _tPrior;
  Double_t  _chi2unfold;
  Double_t  _chi2backfold;
  Double_t  _pMin;
  Double_t  _pMax;
  Double_t  _uMax;
  Double_t  _bMax;
  // ROOT members
  TF1       *_fResiduals;
  TF1       *_fSmoothEff;
  TF1       *_fFineTuneEff;
  TF1       *_fEfficiency;
  TF1       *_fHighPtEff;
  TF1       *_fQtSmooth[NPtPar];
  TF1       *_fSmoothPrior;
  TF1       *_fParToSample;
  TF1       *_fCutoff;
  TF1       *_fLevy;
  TF1       *_fTsallis;
  TF1       *_fExponential;
  TF1       *_fPowerLaw;
  TF1       *_fLandau;
  TH1D      *_hPrior;
  TH1D      *_hPriorDiff;
  TH1D      *_hPriorSmooth;
  TH1D      *_hSmeared;
  TH1D      *_hSmearedDiff;
  TH1D      *_hSmearedSmooth;
  TH1D      *_hMeasured;
  TH1D      *_hUnfolded;
  TH1D      *_hBackfolded;
  TH1D      *_hNormalize;
  TH1D      *_hBackVsMeasRatio;
  TH1D      *_hUnfoldVsPriRatio;
  TH1D      *_hSmearVsMeasRatio;
  TH1D      *_hUnfoldVsMeasRatio;
  TH1D      *_hSmearVsPriRatio;
  TH1D      *_hDvector;
  TH1D      *_hSVvector;
  TH1D      *_hUnfoldErrors;
  TH1D      *_hEfficiency;
  TH1D      *_hResiduals;
  TH2D      *_hPearson;
  TH2D      *_hResponse;
  TH2D      *_hResponseDiff;
  TH2D      *_hResponseSmooth;
  TFile     *_fOut;
  TString   *_sEvnt;
  TString   *_sTrig;
  TString   *_sJet1;
  TString   *_sJet2;
  TString   *_sJet3;
  TRandom   *_rando;
  TPaveText *_label;
  TPaveText *_pInfo;
  // RooUnfold members
  RooUnfoldResponse *_response;

  // private methods ('StJetFolder.sys.h')
  void     PrintInfo(const Int_t code);
  void     PrintError(const Int_t code);
  void     InitializePriors();
  void     SmoothPrior();
  void     DoTwoExpoSmoothing();
  void     DoThreeExpoSmoothing();
  void     DoFourExpoSmoothing();
  void     ConstrainPriorRange(const UInt_t code);
  void     SmoothResponse();
  Bool_t   CheckFlags();
  // private methods ('StJetFolder.plot.h')
  void     CreateLabel();
  void     CreatePlots();
  void     DoRebinning(TH1 *hToRebin, TH1 *hComparison);
  void     ResizeString(TString &str, const Int_t nDec);
  void     DrawHistogram(TH1 *h, const Char_t *option, const Int_t mColor, const Int_t lColor, const Int_t fColor, const Int_t mStyle, const Int_t lStyle, const Int_t fStyle, const Int_t lWidth, const Double_t mSize, const Double_t fAlpha);
  void     CreateUnfoldInfo();
  // private methods ('StJetFolder.math.h')
  TH1D*    CalculateRatio(const TH1D *hA, const TH1D *hB, const Char_t *rName);
  TH2D*    GetPearsonCoefficient(TMatrixD *mCovMat, Bool_t isInDebugMode=false, TString sHistName="");
  UInt_t   ApplyEff(const Double_t par);
  Double_t Smear(const Double_t yP);
  Double_t CalculateChi2(const TH1D *hA, TH1D *hB);
  


  ClassDef(StJetFolder, 1)

};



#endif
#ifdef StJetFolder_cxx

StJetFolder::StJetFolder(const Char_t *oFile, const Bool_t pearDebug, const Int_t nMcSmooth) {

  _fOut  = new TFile(oFile, "recreate");
  _rando = new TRandom();
  for (Int_t i = 0; i < Nflag; i++) {
    _flag[i] = false;
  }
  _pearsonDebug = pearDebug;
  _pTcutoff     = PtCutoffDef;
  _aCutoff      = AcutoffDef;
  _nMcSmooth    = nMcSmooth;
  PrintInfo(0);

}  // end 'StJetFolder(Char_t*)'


StJetFolder::~StJetFolder() {

}  // end '~StJetFolder()'

#endif

// End ------------------------------------------------------------------------
