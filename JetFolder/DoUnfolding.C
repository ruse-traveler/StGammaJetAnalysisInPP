// 'DoUnfolding.C'
// Derek Anderson
// 12.17.2016
// 
// This macro performs the unfolding / backfolding of a provided jet spectrum
// using the 'StJetFolder' class. Parameters:
// 
//   mIn    -- name of root-file containing measured spectrum
//   rIn    -- name of root-file containing the response matrix
//             and reconstruction efficiency
//   out    -- name of output root-file
//   mName  -- name of measured spectrum's histogram
//   rName  -- name of response matrix's histogram
//   eName  -- name of reconstruction effiency's profile
//   method -- which unfolding algorithm to be used:
//               0 = no unfolding
//               1 = bayesian
//               2 = SVD
//               3 = bin-by-bin
//   kReg   -- regularization parameter (1~5 is usually
//             good, anything more will be too sensitive
//             to statistical fluctuations)
//   nMC    -- no. of iterations of backfolding (~100000 is 
//             usually sufficient)
//   prior  -- which type of prior to be used:
//               0 = pythia
//               1 = Levy
//               2 = Tsallis
//               3 = Exponential
//               4 = Power
//               5 = Landau
//   nPrior -- adjusts how fast Levy fncn. drops off
//             (~5.8 typically produces good results)
//   tPrior -- adjusts slope of Levy fncn. (~0.4
//             typically produces good results)
//
// NOTE: if using a prior option other than 0, then the code will
//       use the provided prior and smeared prior to define the
//       binning of the calculated prior, normalize it, and weight
//       the calculated response accordingly.
//
// Pearson Coefficient calculation adapted from Rhagav K. Elayavalli.

#include <TSystem>
#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TLine.h"
#include "TString.h"
#include "TDatime.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

using namespace std;


class StJetFolder;


// input and output files
static const TString pFile("input/embed/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.root");
static const TString sFile("input/embed/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.root");
static const TString mFile("closure/closureTest_modStats_et911r05rff/modifiedStats.detLvlRFF_withBinAndEtaCorr_pTbinFine.et911r05pi0.d9m8y2021.root");
static const TString eFile("input/embed/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.root");
static const TString rFile("input/embed/pp200r9embed.forUnfolding_pTbinFine.et911pt0230x021Kvz55pi0.r05a065rm1chrg.dr05qt05130.root");
//static const TString oFile("pp200r9rff.forDefLevySys_detCutButNormalRes_pTbinFine.et911r05qt05130");
static const TString oFile("pp200r9rff.testingEffFineTune_smoothWithFineTune.et911r05qt05130");
// input histograms
static const TString pName("hSumParFF");
static const TString sName("hSumDetFF");
static const TString mName("hSimJetOutput");
static const TString eName("hEfficiencyFF");
static const TString rName("hResponseFF");

// unfolding parameters (to loop over)
static const Int_t nM  = 1;
static const Int_t M[] = {1};
static const Int_t nK  = 1;
static const Int_t K[] = {3};

// prior parameters (to loop over)
static const Int_t    nP  = 1;
static const Int_t    nN  = 1;
static const Int_t    nT  = 1;
static const Int_t    P[] = {0};
static const Double_t N[] = {4.6};
static const Double_t T[] = {0.4};

// prior, smoothing, and MC parameters
static const Int_t   nMcBack       = 10000;    // no. of MC iterations for backfolding
static const Int_t   nMcSmooth     = 5000000;  // no. of MC iterations for smoothing/retraining response matrix
static const Bool_t  doEffSmooth   = true;     // smooth efficiency
static const Int_t   doResSmooth   = true;     // smooth response
static const Bool_t  doPriorSmooth = false;    // smooth prior
static const Bool_t  doPriorCutoff = false;    // apply a smooth cutoff to prior
static const Bool_t  doEffFineTune = true;     // fine tune efficiency after smoothing
static const Bool_t  removeEffErr  = false;    // remove errors on efficiency
static const UInt_t  eTtrgIndex    = 0;        // eTrg range (controls how prior is smoothed): 0 = 9 - 11, 1 = 11 - 15, 2 = 15 - 20 
static const UInt_t  trigType      = 2;        // 0 = "gamma-dir", 1 = "gamma-rich", 2 = "pi0", 3 = "h+-" (controls if negative bins are checked for)
static const Float_t pTpriMin      = 0.;
static const Float_t pTpriMax      = 57.;
static const Float_t pTcutoff      = 2.;
static const Float_t aCutoff       = 1.;


// trigger and jet parameter (doesn't impact unfolding)
static const Int_t   beam   = 0;     // 0 = "pp", 1 = "AuAu"
static const Int_t   type   = 0;     // 0 = "charged jets", 1 = "full jets"
static const Float_t energy = 200.;  // sqrt(s)
static const Float_t rJet   = 0.5;

// jet parameters (doesn't impact unfolding)
static const Int_t    nRM     = 1;
static const Double_t aMin    = 0.65;
static const Double_t pTmin   = 0.2;
static const Double_t pTmaxU  = 57.;
static const Double_t pTmaxB  = 30.;
static const Double_t hTrgMax = 0.9;

// misc. parameters (for debugging and such)
static const Int_t    nToy     = 10;       // used to calculate covariances
static const Bool_t   debug    = false;    // debug pearson calculation coefficient
static const Double_t bPrior   = 2.40;     // amplitude for alternate priors
static const Double_t mPrior   = 0.140;    // m-parameter of levy function



void DoUnfolding() {

  gSystem -> Load("../../RooUnfold/libRooUnfold.so");
  gSystem -> Load("StJetFolder");
  
  // lower verbosity
  gErrorIgnoreLevel = kError;

  TDatime start;
  cout << "\nStarting folding: " << start.AsString() << "\n" << endl;


  // create output stream
  const TString sStream(oFile.Data());
  sStream += ".bestFiles.list";

  ofstream bestFiles(sStream.Data());
  if (!bestFiles) {
    cerr << "PANIC: couldn't open output stream!" << endl;
    return;
  }


  // prior loops
  Double_t chi2bestest = 999.;
  TString  bestestFile;
  for (Int_t p = 0; p < nP; p++) {
    for (Int_t n = 0; n < nN; n++) {
      for (Int_t t = 0; t < nT; t++) {

        // don't double count priors...
        const Bool_t isPyth   = (P[p] == 0);
        const Bool_t isExpo   = (P[p] == 3);
        const Bool_t isPowr   = (P[p] == 4);
        const Bool_t isFirstN = (n == 0);
        const Bool_t isFirstT = (t == 0);
        const Bool_t isFirst  = (isFirstN || isFirstT);
        if (isPyth && !isFirst)  continue;
        if (isExpo && !isFirstN) continue;
        if (isPowr && !isFirstN) continue;

        // for file names
        const Double_t nPrior = N[n];
        const Double_t tPrior = T[t];
        const Float_t  Ntxt = N[n] * 10.;
        const Float_t  Ttxt = T[t] * 10.;

        // for recording chi2
        const TString sBayU("hBayUnfold");
        const TString sBayB("hBayBackfold");
        const TString sSvdU("hSvdUnfold");
        const TString sSvdB("hSvdBackfold");
        const TString sChiX("k_{reg}");
        const TString sChiYU("#chi^{2}(unfold, prior)");
        const TString sChiYB("#chi^{2}(backfold, measured)");

        // create performance file
        TString sChi2(oFile.Data());
        sChi2 += ".p";
        sChi2 += P[p];
        sChi2 += "n";
        sChi2 += Ntxt;
        sChi2 += "t";
        sChi2 += Ttxt;
        sChi2 += ".performance.root";

        TFile *fChi2 = new TFile(sChi2.Data(), "recreate");
        TH1D  *hBayUnfold    = new TH1D(sBayU.Data(), "", nK, K[0], K[nK - 1] + 1);
        TH1D  *hBayBackfold  = new TH1D(sBayB.Data(), "", nK, K[0], K[nK - 1] + 1);
        TH1D  *hSvdUnfold    = new TH1D(sSvdU.Data(), "", nK, K[0], K[nK - 1] + 1);
        TH1D  *hSvdBackfold  = new TH1D(sSvdB.Data(), "", nK, K[0], K[nK - 1] + 1);
        hBayUnfold   -> SetTitleFont(42);
        hBayUnfold   -> GetXaxis() -> SetTitle(sChiX.Data());
        hBayUnfold   -> GetXaxis() -> SetTitleOffset(1.);
        hBayUnfold   -> GetXaxis() -> SetTitleFont(42);
        hBayUnfold   -> GetXaxis() -> SetLabelFont(42);
        hBayUnfold   -> GetYaxis() -> SetTitle(sChiYU.Data());
        hBayUnfold   -> GetYaxis() -> SetTitleFont(42);
        hBayUnfold   -> GetYaxis() -> SetLabelFont(42);
        hBayBackfold -> SetTitleFont(42);
        hBayBackfold -> GetXaxis() -> SetTitle(sChiX.Data());
        hBayBackfold -> GetXaxis() -> SetTitleOffset(1.);
        hBayBackfold -> GetXaxis() -> SetTitleFont(42);
        hBayBackfold -> GetXaxis() -> SetLabelFont(42);
        hBayBackfold -> GetYaxis() -> SetTitle(sChiYB.Data());
        hBayBackfold -> GetYaxis() -> SetTitleFont(42);
        hBayBackfold -> GetYaxis() -> SetLabelFont(42);
        hSvdUnfold   -> SetTitleFont(42);
        hSvdUnfold   -> GetXaxis() -> SetTitle(sChiX.Data());
        hSvdUnfold   -> GetXaxis() -> SetTitleOffset(1.);
        hSvdUnfold   -> GetXaxis() -> SetTitleFont(42);
        hSvdUnfold   -> GetXaxis() -> SetLabelFont(42);
        hSvdUnfold   -> GetYaxis() -> SetTitle(sChiYU.Data());
        hSvdUnfold   -> GetYaxis() -> SetTitleFont(42);
        hSvdUnfold   -> GetYaxis() -> SetLabelFont(42);
        hSvdBackfold -> SetTitleFont(42);
        hSvdBackfold -> GetXaxis() -> SetTitle(sChiX.Data());
        hSvdBackfold -> GetXaxis() -> SetTitleOffset(1.);
        hSvdBackfold -> GetXaxis() -> SetTitleFont(42);
        hSvdBackfold -> GetXaxis() -> SetLabelFont(42);
        hSvdBackfold -> GetYaxis() -> SetTitle(sChiYB.Data());
        hSvdBackfold -> GetYaxis() -> SetTitleFont(42);
        hSvdBackfold -> GetYaxis() -> SetLabelFont(42);
        hBayUnfold   -> Sumw2();
        hBayBackfold -> Sumw2();
        hSvdUnfold   -> Sumw2();
        hSvdBackfold -> Sumw2();

        // method and k loop
        Double_t chi2best = 999.;
        TString  bestFile;
        for (Int_t m = 0; m < nM; m++) {
          for (Int_t k = 0; k < nK; k++) {

            const Int_t prior  = P[p];
            const Int_t method = M[m];
            const Int_t kReg   = K[k];

            // don't double count bin-by-bin corrections
            const Bool_t isBinByBin = (method == 3);
            const Bool_t isFirstK   = (k == 0);
            if (isBinByBin && !isFirstK) continue;

            // skip unreasonable kReg
            const Bool_t isBay      = (method == 1);
            const Bool_t isSVD      = (method == 2);
            const Bool_t isGoodBayK = (kReg < 6);
            const Bool_t isGoodSvdK = ((kReg > 5) && (kReg < 12));
            //if (isBay && !isGoodBayK) continue;
            //if (isSVD && !isGoodSvdK) continue;

            // create output name
            TString output(oFile);
            output += ".p";
            output += prior;
            output += "m";
            output += method;
            output += "k";
            output += kReg;
            output += "n";
            output += Ntxt;
            output += "t";
            output += Ttxt;
            output += ".root";

            // create folder
            Double_t    chi2u = 0.;
            Double_t    chi2b = 0.;
            StJetFolder f(output.Data(), debug, nMcSmooth);
            // set spectra
            f.SetPrior(pFile.Data(), pName.Data(), doPriorSmooth, doPriorCutoff);
            f.SetSmeared(sFile.Data(), sName.Data());
            f.SetMeasured(mFile.Data(), mName.Data());
            f.SetResponse(rFile.Data(), rName.Data(), doResSmooth);
            f.SetEfficiency(eFile.Data(), eName.Data(), doEffSmooth, doEffFineTune, removeEffErr);
            // set info and parameters
            f.SetEventInfo(beam, energy);
            f.SetTriggerInfo(trigType, eTtrgIndex, hTrgMax);
            f.SetJetInfo(type, nRM, rJet, aMin, pTmin);
            if (doPriorCutoff) f.SetPriorCutoffParameters(pTcutoff, aCutoff);
            f.SetPriorParameters(prior, bPrior, mPrior, nPrior, tPrior, pTpriMin, pTpriMax);
            f.SetUnfoldParameters(method, kReg, nToy, pTmaxU, pTmaxB);
            // do unfolding
            f.Init();
            f.Unfold(chi2u);
            f.Backfold(chi2b, nMcBack);
            f.Finish();

            const Double_t merit   = TMath::Abs(chi2b - 1);
            const Double_t best    = TMath::Abs(chi2best - 1);
            const Double_t bestest = TMath::Abs(chi2bestest - 1);
            if (merit < best) {
              chi2best = chi2b;
              bestFile = output;
            }
            if (merit < bestest) {
              chi2bestest = chi2b;
              bestestFile = output;
            }

            // record chi2
            UInt_t iReg(0);
            switch (method) {
              case 1:
                iReg = hBayBackfold -> FindBin(kReg);
                hBayUnfold   -> SetBinContent(iReg, chi2u);
                hBayUnfold   -> SetBinError(iReg, 0.);
                hBayBackfold -> SetBinContent(iReg, chi2b);
                hBayBackfold -> SetBinError(iReg, 0.);
                break;
              case 2:
                iReg = hSvdBackfold -> FindBin(kReg);
                hSvdUnfold   -> SetBinContent(iReg, chi2u);
                hSvdUnfold   -> SetBinError(iReg, 0.);
                hSvdBackfold -> SetBinContent(iReg, chi2b);
                hSvdBackfold -> SetBinError(iReg, 0.);
                break;
              default:
                break;
            }

          }  // end k loop
        }  // end method loop


        // make performance plots
        const UInt_t cBay(810);
        const UInt_t cSvd(860);
        hBayUnfold   -> SetLineColor(cBay);
        hBayUnfold   -> SetMarkerColor(cBay);
        hBayBackfold -> SetLineColor(cBay);
        hBayBackfold -> SetMarkerColor(cBay);
        hSvdUnfold   -> SetLineColor(cSvd);
        hSvdUnfold   -> SetMarkerColor(cSvd);
        hSvdBackfold -> SetLineColor(cSvd);
        hSvdBackfold -> SetMarkerColor(cSvd);

        TLegend *lUnfold   = new TLegend(0.1, 0.1, 0.3, 0.3);
        TLegend *lBackfold = new TLegend(0.1, 0.1, 0.3, 0.3);
        lUnfold   -> SetFillColor(0);
        lUnfold   -> SetLineColor(0);
        lUnfold   -> SetTextFont(42);
        lUnfold   -> SetTextAlign(12);
        lUnfold   -> AddEntry(hBayUnfold, "Bayes.");
        lUnfold   -> AddEntry(hSvdUnfold, "SVD");
        lBackfold -> SetFillColor(0);
        lBackfold -> SetLineColor(0);
        lBackfold -> SetTextFont(42);
        lBackfold -> SetTextAlign(12);
        lBackfold -> AddEntry(hBayBackfold, "Bayes.");
        lBackfold -> AddEntry(hSvdBackfold, "SVD");

        TLine *lOne = new TLine(K[0], 1, K[nK - 1] + 1, 1);
        lOne  -> SetLineColor(1);
        lOne  -> SetLineStyle(2);
        fChi2 -> cd();

        TCanvas *cUnfold   = new TCanvas("cUnfold", "", 750, 500);
        TCanvas *cBackfold = new TCanvas("cBackfold", "", 750, 500);
        cUnfold      -> SetGrid(0, 0);
        cBackfold    -> SetGrid(0, 0);
        cUnfold      -> cd();
        hBayUnfold   -> Draw();
        hSvdUnfold   -> Draw("same");
        lUnfold      -> Draw();
        lOne         -> Draw();
        cBackfold    -> cd();
        hBayBackfold -> Draw();
        hSvdBackfold -> Draw("same");
        lBackfold    -> Draw();
        lOne         -> Draw();
        cUnfold      -> Write();
        cUnfold      -> Close();
        cBackfold    -> Write();
        cBackfold    -> Close();

        // save chi2
        fChi2        -> cd();
        hBayUnfold   -> Write();
        hBayBackfold -> Write();
        hSvdUnfold   -> Write();
        hSvdBackfold -> Write();
        fChi2        -> Close();

        // announce winner
        TDatime endPrior;
        cout << "\nFinished folding prior! " << endPrior.AsString() << "\n"
             << "  Best chi2 = " << chi2best << "\n"
             << "  Best file = " << bestFile << "\n"
             << endl;

        // stream winner
        bestFiles << bestFile.Data();
        bestFiles << endl;

      }  // end tPrior loop
    }  // end nPrior loop
  }  // end prior loop


  // announce biggest winner
  TDatime end;
  cout << "\nFinished all folding! " << end.AsString() << "\n"
       << "  Bestest chi2 = " << chi2bestest << "\n"
       << "  Bestest file = " << bestestFile << "\n"
       << endl;
}

// End ------------------------------------------------------------------------
