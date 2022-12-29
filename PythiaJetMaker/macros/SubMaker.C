// 'SubMaker.C'
// Derek Anderson
//
// This grabs a histogram (or histograms) of recoil jet pTcorr and
// "uncorrelated" jet pTcorr and creates a histogram of pTsub
// (pTsub = pTcorr(recoil) - pTcorr(uncorr.)).

#include <TSystem>
#include <iostream>
#include "TH1.h"

using namespace std;


void SubMaker() {

  TFile *oFile = new TFile("pp200aug.r03a02rm1.g.Aug14.root", "recreate");
  TFile *iFile = new TFile("/global/project/projectdirs/star/pwg/starjetc/nrsahoo/Jet_pp_data/treeOfJets/pp200_R03_ChargedJet_Aug.root");

  // retrieve histograms
  TH1D *hR = (TH1D*) iFile -> Get("h_JetPt_gamma_Recoil");
  TH1D *hU = (TH1D*) iFile -> Get("h_JetPt_gamma_UE");

  // get histogram dimensions
  Int_t    nR = hR -> GetNbinsX();
  Double_t r1 = hR -> GetBinLowEdge(1);
  Double_t r2 = hR -> GetBinLowEdge(nR+1); 

  // subtract histograms
  TH1D *hPtCorr_R = (TH1D*) hR -> Clone("hPtRE");
  TH1D *hPtCorr_U = (TH1D*) hU -> Clone("hPtUE");
  TH1D *hPtSub    = new TH1D("hPtSub", "p_{T}^{sub} = p_{T}^{corr}(RE) - p_{T}^{corr}(UE)", nR, r1, r2);
  hPtSub -> Add(hPtCorr_R, hPtCorr_U, 1., -1.);

  // save and close file
  oFile     -> cd();
  hPtCorr_R -> Write();
  hPtCorr_U -> Write();
  hPtSub    -> Write();
  oFile     -> Close();

}

// End ------------------------------------------------------------------------
