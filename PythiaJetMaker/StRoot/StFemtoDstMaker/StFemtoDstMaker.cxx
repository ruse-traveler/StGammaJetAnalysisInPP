// 'StFemtoDstMaker.cxx'
// Derek Anderson
// 08.04.2016
//
// This class takes output from Pythia, finds jets, and stores in them in a
// compact tree.
//
// Last updated: 04.17.2018


#define StFemtoDstMaker_cxx

#include <cmath>
#include <vector>
#include "TDirectory.h"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StFemtoDstMaker.h"

ClassImp(StFemtoDstMaker)

using namespace std;
using namespace fastjet;



void StFemtoDstMaker::SetEfficiency(const TString sEffFile, const TString sEffHist, const Bool_t isHist, const Bool_t isFunc, const Bool_t doHighPtSmooth, const Bool_t doEffAdjust, const Float_t effNudge, const Float_t amplitude, const Float_t sigma) {

  effIsHist = isHist;
  effIsFunc = isFunc;
  if (effIsHist || effIsFunc) {
    // open file with efficiency
    TFile *fEff = (TFile*) gROOT -> GetListOfFiles() -> FindObject(sEffFile.Data());
    if (!fEff || !fEff->IsOpen())
      fEff = new TFile(sEffFile.Data());

    // grab efficiency
    if (effIsHist) {
      fEff -> GetObject(sEffHist.Data(), hPtEff);
      if (!hPtEff)
        cerr << "WARNING: couldn't grab efficiency histogram! Check input!" << endl;
    }
    if (effIsFunc) {
      fEff -> GetObject(sEffHist.Data(), fPtEff);
      if (!fPtEff)
        cerr << "WARNING: couldn't grab efficiency function! Check input!" << endl;
    }
  }
  else {
    effAmplitude = amplitude;
    effSigma     = sigma;

    fPtEff = new TF1("fPtEff", "[0]*(1.-TMath::Exp(-1.*[1]*x))", EffPtMin, EffPtMax);
    fPtEff -> SetParameter(0, effAmplitude);
    fPtEff -> SetParameter(1, effSigma);
  }

  // for smoothing and systematics
  doEffSmooth = doHighPtSmooth;
  doEffSys    = doEffAdjust;
  effAdjust   = effNudge;

}  // end 'SetEfficiency(TString, TString)'



void StFemtoDstMaker::Init(const Int_t nRM, const Int_t tID, const Double_t r, const Double_t a, const Double_t pMin, const Double_t pMax, const Double_t q, const Double_t eMin, const Double_t eMax, const Double_t hTrkMax, const Double_t hTrgMax) {

  Nrm       = nRM;
  trigID    = tID;
  Rjet      = r;
  Amin      = a;
  pTmin     = pMin;
  pTmax     = pMax;
  qMin      = q; 
  eTmin     = eMin;
  eTmax     = eMax;
  EtaTrkMax = hTrkMax;
  EtaTrgMax = hTrgMax;
  oFile     -> cd();

  // initialize output tree
  fDst = new TTree("femtoDST", "A compact tree of jets");
  fDst -> Branch("EventIndex", &femto.EventIndex, "EventIndex/I");
  fDst -> Branch("Refmult", &femto.Refmult, "Refmult/D");
  fDst -> Branch("TSP", &femto.TSP, "TSP/D");
  fDst -> Branch("TrgEta", &femto.TrgEta, "TrgEta/D");
  fDst -> Branch("TrgPhi", &femto.TrgPhi, "TrgPhi/D");
  fDst -> Branch("TrgEt", &femto.TrgEt, "TrgEt/D");
  fDst -> Branch("Rho", &femto.Rho, "Rho/D");
  fDst -> Branch("Sigma", &femto.Vz, "Vz/D");
  fDst -> Branch("JetPt", &femto.JetPt);
  fDst -> Branch("JetNCons", &femto.JetNCons);
  fDst -> Branch("JetIndex", &femto.JetIndex);
  fDst -> Branch("JetEta", &femto.JetEta);
  fDst -> Branch("JetPhi", &femto.JetPhi);
  fDst -> Branch("JetE", &femto.JetE);
  fDst -> Branch("JetArea", &femto.JetArea);
  fDst -> Branch("JetConsPt", &femto.JetConsPt);
  fDst -> Branch("JetConsEta", &femto.JetConsEta);
  fDst -> Branch("JetConsPhi", &femto.JetConsPhi);
  fDst -> Branch("JetConsE", &femto.JetConsE);
  fDst -> SetAutoSave(10000);

  TDirectory *dQA = oFile -> mkdir("QA");

  // initialize QA histograms
  const Int_t    nTSP = 300;
  const Int_t    nPt  = 1200;
  const Int_t    nP   = 240;
  const Int_t    nDf  = 720;
  const Int_t    nH   = 200;
  const Int_t    nR   = 50;
  const Int_t    nA   = 50;
  const Double_t tsp1 = 0.;
  const Double_t tsp2 = 300.;
  const Double_t pT1  = -20.;
  const Double_t pT2  = 100.;
  const Double_t p1   = -60.;
  const Double_t p2   = 60.;
  const Double_t dF1  = -1.*M_PI/2.;
  const Double_t dF2  = 3.*M_PI/2.;
  const Double_t h1   = -10.;
  const Double_t h2   = 10.;
  const Double_t r1   = 0.;
  const Double_t r2   = 5.;
  const Double_t a1   = 0.;
  const Double_t a2   = 5.;
  dQA -> cd();
  // initialize event histograms
  hRefmult    = new TH1D("hRefmult", "Reference mult.", 200, 0., 200.);
  hTSP        = new TH1D("hTSP", "TSP (absolute value of trigger ID)", nTSP, tsp1, tsp2);
  hEtTrg      = new TH1D("hEtTrg", "Trigger e_{T}", nPt, pT1, pT2);
  hFtrg       = new TH1D("hFtrg", "Trigger #varphi", nDf, dF1, dF2);
  hHtrg       = new TH1D("hHtrg", "Trigger #eta", nH, h1, h2);
  hRho        = new TH1D("hRho", "Median background energy density, #rho", nR, r1, r2);
  hNjet       = new TH1D("hNjet", "No. of jets", 100, 0., 100.);
  hNjetG17    = new TH1D("hNjetG17", "No. of jets with p_{T}^{jet}>17 GeV/c", 200, 0., 200.);
  hNjetL17    = new TH1D("hNjetL17", "No. of jets with p_{T}^{jet}<17 GeV/c", 200, 0., 200.);
  hNtrk       = new TH1D("hNtrk", "No. of tracks passing cuts", 200, 0., 200.);
  hNtrkG1     = new TH1D("hNtrkG1", "No. of tracks with p_{T}>1 GeV/c", 200, 0., 200.);
  hNtrkG2     = new TH1D("hNtrkG2", "No. of tracks with p_{T}>2 GeV/c", 200, 0., 200.);
  // initialize track histograms
  hPtTrk      = new TH1D("hPtTrk", "Track p_{T}", nPt, pT1, pT2);
  hPxTrk      = new TH1D("hPxTrk", "Track p_{x}", nP, p1, p2);
  hPyTrk      = new TH1D("hPyTrk", "Track p_{y}", nP, p1, p2);
  hPzTrk      = new TH1D("hPzTrk", "Track p_{z}", nP, p1, p2);
  hDfTrk      = new TH1D("hDfTrk", "Track #Delta#varphi", nDf, dF1, dF2);
  hDfTrkG1    = new TH1D("hDfTrkG1", "#Delta#varphi of tracks with p_{T}>1 GeV/c", nDf, dF1, dF2);
  hDfTrkG2    = new TH1D("hDfTrkG2", "#Delta#varphi of tracks with p_{T}>2 GeV/c", nDf, dF1, dF2);
  hHtrk       = new TH1D("hHtrk", "Track #eta", nH, h1, h2);
  hPtEffPar   = new TH1D("hPtEffPar", "Particle-lvl. track p_{T}", nPt, pT1, pT2);
  hPtEffDet   = new TH1D("hPtEffDet", "Detector-lvl. track p_{T}", nPt, pT1, pT2);
  hPtEffCheck = new TH1D("hPtEffCheck", "Efficiency check", nPt, pT1, pT2);
  hPtRes1D    = new TH1D("hPtRes1D", "", 200, -10., 10.);
  hPtRes2D    = new TH2D("hPtRes2D", "", 300, 0., 30., 200, -10., 10.);
  hPtVsDfTrk  = new TH2D("hPtVsDfTrk", "Track p_{T} vs. #Delta#varphi; #varphi; p_{T}", nDf, dF1, dF2, nPt, pT1, pT2);
  // initialize jet histograms
  hPtReco     = new TH1D("hPtReco", "Jet p_{T}^{raw} (p_{T}^{reco})", nPt, pT1, pT2);
  hPtCorr     = new TH1D("hPtCorr", "Jet p_{T}^{corr}= p_{T}^{raw}-#rhoA", nPt, pT1, pT2);
  hPtRE       = new TH1D("hPtRE", "p_{T}^{corr} for recoil jets", nPt, pT1, pT2);
  hPtUE       = new TH1D("hPtUE", "p_{T}^{corr} for 'uncorrelated' jets", nPt, pT1, pT2);
  hPtSub      = new TH1D("hPtSub", "p_{T}^{corr}(RE) - p_{T}^{corr}(UE)", nPt, pT1, pT2);
  hDfJet      = new TH1D("hDfJet", "Jet #Delta#varphi", nDf, dF1, dF2);
  hHjet       = new TH1D("hHjet", "Jet #eta", nH, h1, h2);
  hAjet       = new TH1D("hAjet", "Jet area", nA, a1, a2);
  // initialize cst. histograms
  hPtCst      = new TH1D("hPtCst", "Constituent p_{T}", nPt, pT1, pT2);
  hDfCst      = new TH1D("hDfCst", "Constituent #Delta#varphi", nDf, dF1, dF2);
  hHcst       = new TH1D("hHcst", "Constituent #eta", nH, h1, h2);
  // set errors
  hRefmult    -> Sumw2();
  hTSP        -> Sumw2();
  hEtTrg      -> Sumw2();
  hFtrg       -> Sumw2();
  hHtrg       -> Sumw2();
  hRho        -> Sumw2();
  hNjet       -> Sumw2();
  hNjetG17    -> Sumw2();
  hNjetL17    -> Sumw2();
  hNtrk       -> Sumw2();
  hNtrkG1     -> Sumw2();
  hNtrkG2     -> Sumw2();
  hPtTrk      -> Sumw2();
  hPxTrk      -> Sumw2();
  hPyTrk      -> Sumw2();
  hPzTrk      -> Sumw2();
  hDfTrk      -> Sumw2();
  hDfTrkG1    -> Sumw2();
  hDfTrkG2    -> Sumw2();
  hHtrk       -> Sumw2();
  hPtEffPar   -> Sumw2();
  hPtEffDet   -> Sumw2();
  hPtEffCheck -> Sumw2();
  hPtRes1D    -> Sumw2();
  hPtRes2D    -> Sumw2();
  hPtVsDfTrk  -> Sumw2();
  hPtReco     -> Sumw2();
  hPtCorr     -> Sumw2();
  hPtRE       -> Sumw2();
  hPtUE       -> Sumw2();
  hPtSub      -> Sumw2();
  hDfJet      -> Sumw2();
  hHjet       -> Sumw2();
  hAjet       -> Sumw2();
  hPtCst      -> Sumw2();
  hDfCst      -> Sumw2();
  hHcst       -> Sumw2();
  oFile       -> cd();

  femto.Reset();

  // announce parameters
  cout << "  Initializing...\n"
       << "    Rjet=" << Rjet << ", Nrm=" << Nrm << ", Amin=" << Amin
       << ", pTmin=" << pTmin << ", pTmax=" << pTmax << "\n" 
       << "    eTmin=" << eTmin << ", eTmax=" << eTmax
       << ", etaTrkMax=" << hTrkMax << ", etaTrgMax=" << hTrgMax
       << endl;

  // TEST [04.29.2020]
  const UInt_t   nDfPtBin   = 60;
  const Double_t dFpTbin[2] = {-1.57, 4.71};
  hDfTrk_pt021 = new TH1D("hDfTrk_pt021", "#Delta#varphi of tracks with p_{T} = 0.2 - 1 GeV/c", nDfPtBin, dFpTbin[0], dFpTbin[1]);
  hDfTrk_pt12  = new TH1D("hDfTrk_pt12", "#Delta#varphi of tracks with p_{T}> = 1 - 2 GeV/c", nDfPtBin, dFpTbin[0], dFpTbin[1]);
  hDfTrk_pt25  = new TH1D("hDfTrk_pt25", "#Delta#varphi of tracks with p_{T} = 2 - 5 GeV/c", nDfPtBin, dFpTbin[0], dFpTbin[1]);
  hDfTrk_pt530 = new TH1D("hDfTrk_pt530", "#Delta#varphi of tracks with p_{T} = 5 - 30 GeV/c", nDfPtBin, dFpTbin[0], dFpTbin[1]);
  hDfTrk_pt021 -> Sumw2();
  hDfTrk_pt12  -> Sumw2();
  hDfTrk_pt25  -> Sumw2();
  hDfTrk_pt530 -> Sumw2();

}  // end 'Init()'



void StFemtoDstMaker::Make(const Int_t nTrgs, const Int_t StartEvt, const Int_t StopEvt) {

  if (fChain == 0) return;


  // for jet finding
  vector<Double_t>  conPt;
  vector<Double_t>  conH;
  vector<Double_t>  conF;
  vector<Double_t>  conE;
  vector<PseudoJet> tracks;
  vector<PseudoJet> jetsCS;
  vector<PseudoJet> jets;

  const Int_t    nRepeat   = 1;
  const Double_t ghostArea = 0.01;
  const Double_t maxGhostH = 1.0 + Rjet;
  const Double_t absEtaMax = 1.0 - Rjet;


  // determine which events to process
  Long64_t nEvt   = fChain -> GetEntriesFast();
  Long64_t iStart = StartEvt;
  Long64_t iStop  = nEvt;
  if ((StopEvt != -1) && (StopEvt < nEvt))
    iStop = StopEvt;
  cout << "  Beginning event loop, processing events " << iStart+1
       << " to " << iStop+1 << "..."
       << endl;

  // how many triggers do you want?
  Int_t nTrgMax = nTrgs;
  if (nTrgMax == -1)
    nTrgMax = nEvt;
  cout << "  Looking for " << nTrgMax << " triggers..." << endl;


  // event loop
  Long64_t nTrg   = 0;
  Long64_t nBytes = 0;
  Long64_t bytes  = 0;
  for (Long64_t i = iStart; i < iStop; i++) {

    Long64_t entry = LoadTree(i);
    if (entry < 0) break;
    bytes   = fChain -> GetEntry(i);
    nBytes += bytes;

    cout << "    Processing event " << i+1 << "/" << nEvt << "..." << endl;
    if (i+1 == nEvt) cout << endl;

    if (nTrg == nTrgMax) {
      cout << "    " << nTrgMax << " triggers found! Ending event loop!" << endl;
      break;
    }

    // TEST [03.16.2022]
    Bool_t isHeavyQuark(false);
    switch (Events_Process) {
      case 121:
        isHeavyQuark = true;
        break;
      case 122:
        isHeavyQuark = true;
        break;
      case 123:
        isHeavyQuark = true;
        break;
      case 124:
        isHeavyQuark = true;
        break;
      default:
        isHeavyQuark = false;
        break;
    }
    //if (isHeavyQuark) continue;


    // trigger cuts
    Double_t TSP   = Events_TrigId;
    Double_t fTrg  = Events_Clust_phiv1;
    Double_t hTrg  = Events_Clust_etav1;
    Double_t tTrg  = 2. * atan(exp(-1. * hTrg));
    Double_t eTrg  = Events_Clust_EneT0;
    Double_t eTtrg = eTrg * sin(tTrg);
    if (abs(TSP) != trigID)
      continue;
    if ((eTtrg < eTmin) || (eTtrg > eTmax))
      continue;
    if (abs(hTrg) > EtaTrgMax)
      continue;

    // store event info
    femto.EventIndex = i;
    femto.Refmult    = Events_refmult;
    femto.TSP        = abs(Events_TrigId);
    femto.TrgEta     = hTrg;
    femto.TrgPhi     = fTrg;
    femto.TrgEt      = eTtrg;
    femto.Vz         = Events_primVz;
    hRefmult -> Fill((Double_t) femto.Refmult);
    hTSP     -> Fill((Double_t) femto.TSP);
    hEtTrg   -> Fill(femto.TrgEt);
    hFtrg    -> Fill(femto.TrgPhi);
    hHtrg    -> Fill(femto.TrgEta);
    ++nTrg;


    // particle loop
    Int_t nTrkG02 = 0;
    Int_t nTrkG1  = 0;
    Int_t nTrkG2  = 0;
    Int_t nTrk    = Events_noOfprimaryTrks;
    for (Int_t j = 0; j < nTrk; j++) {

      Double_t pX   = pTracks_px[j];
      Double_t pY   = pTracks_py[j];
      Double_t pZ   = pTracks_pz[j];
      Double_t pT   = pTracks_pT[j];
      Double_t chrg = pTracks_chrg[j];
      Double_t fTrk = pTracks_Phi[j];
      Double_t hTrk = pTracks_Eta[j];
      Double_t tTrk = 2. * atan(exp(-1. * hTrk));


      // if detector-level, apply response
      if ((level != 0) && (chrg != 0.)) {
        pT = ApplyDetectorResponse(pT);
      }

      if (pT <= -1000.) {
        continue;
      }
      else {
        pX = pT * cos(fTrk);
        pY = pT * sin(fTrk);
        pZ = pT / tan(tTrk);
      }

      // track cuts
      Double_t E = sqrt(pX*pX + pY*pY + pZ*pZ + PionMass*PionMass);
      if ((pT < pTmin) || (pT > pTmax))
        continue;
      if (abs(hTrk) > EtaTrkMax)
        continue;
      if ((type == 0) && (chrg == 0.))
        continue;


      // add track to list of particles
      tracks.push_back(PseudoJet(pX, pY, pZ, E));

      // fill QA histograms
      Double_t dF = fTrk - fTrg;
      if (dF < -1.*M_PI/2.) dF += 2.*M_PI;
      if (dF > 3.*M_PI/2.)  dF -= 2.*M_PI;
      hPtTrk     -> Fill(pT);
      hPxTrk     -> Fill(pX);
      hPyTrk     -> Fill(pY);
      hPzTrk     -> Fill(pZ);
      hDfTrk     -> Fill(dF);
      hHtrk      -> Fill(hTrk);
      hPtVsDfTrk -> Fill(dF, pT);

      ++nTrkG02;
      if (pT > 1.) {
        hDfTrkG1 -> Fill(dF);
        ++nTrkG1;
      }
      if (pT > 2.) {
        hDfTrkG2 -> Fill(dF);
        ++nTrkG2;
      }

      // TEST [04.29.2020]
      if ((pT > 0.2) && (pT < 1.)) {
        hDfTrk_pt021 -> Fill(dF);
      }
      if ((pT > 1.) && (pT < 2.)) {
        hDfTrk_pt12  -> Fill(dF);
      }
      if ((pT > 2.) && (pT < 5.)) {
        hDfTrk_pt25  -> Fill(dF);
      }
      if ((pT > 5.) && (pT < 30.)) {
        hDfTrk_pt530 -> Fill(dF);
      }

    }  // end particle loop

    // fill QA histograms
    hNtrk   -> Fill(nTrkG02);
    hNtrkG1 -> Fill(nTrkG1);
    hNtrkG2 -> Fill(nTrkG2);


    // define jets and jet area
    GhostedAreaSpec areaSpec(maxGhostH, nRepeat, ghostArea);
    AreaDefinition  areaDef(active_area_explicit_ghosts, areaSpec);
    JetDefinition   jetDef(antikt_algorithm, Rjet);

    // cluster jets
    ClusterSequenceArea CS(tracks, jetDef, areaDef);
    jetsCS = sorted_by_pt(CS.inclusive_jets(pTmin));

    // fiducial cut
    Selector fidCut = SelectorAbsEtaMax(absEtaMax);
    jets = fidCut(jetsCS);


    // define bkgd jets and jet area
    GhostedAreaSpec areaSpecBkgd(maxGhostH, nRepeat, ghostArea);
    AreaDefinition  areaDefBkgd(active_area_explicit_ghosts, areaSpecBkgd);
    JetDefinition   jetDefBkgd(kt_algorithm, Rjet);

    // initialize bkgd estimatiors
    Selector bkgdCut = SelectorAbsEtaMax(EtaTrkMax) * (!SelectorNHardest(Nrm));
    JetMedianBackgroundEstimator bkgd(bkgdCut, jetDefBkgd, areaDefBkgd);
    // initialize subtractor
    Subtractor sub(&bkgd);

#if FASTJET_VERSION_NUMBER >= 30100
  sub.set_use_rho_m(true);
  sub.set_safe_mass(true);
#endif

    // estimate bkgd
    bkgd.set_particles(tracks);
    femto.Rho   = bkgd.rho();
    femto.Sigma = bkgd.sigma();


    // jet loop
    Int_t nJetG17 = 0;
    Int_t nJetL17 = 0;
    Int_t nJet    = jets.size();
    for (Int_t j = 0; j < nJet; j++) {

      Double_t pTj = jets[j].perp();
      Double_t hJ  = jets[j].pseudorapidity();
      Double_t fJ  = jets[j].phi();
      Double_t eJ  = jets[j].e();
      Double_t aJ  = jets[j].area();

      Int_t nCon = jets[j].constituents().size();
      femto.JetPt.push_back(pTj);
      femto.JetNCons.push_back(nCon);
      femto.JetIndex.push_back(j);
      femto.JetEta.push_back(hJ);
      femto.JetPhi.push_back(fJ);
      femto.JetArea.push_back(aJ);
      femto.JetE.push_back(eJ);


      // fill QA histograms
      Double_t pTcorr = pTj - femto.Rho*aJ;
      Double_t dFj    = fJ  - fTrg;
      if (dFj < -1.*M_PI/2.) dFj += 2.*M_PI;
      if (dFj > 3.*M_PI/2.)  dFj -= 2.*M_PI;
      hPtReco -> Fill(pTj);
      hPtCorr -> Fill(pTcorr);
      hDfJet  -> Fill(dFj);
      hHjet   -> Fill(hJ);
      hAjet   -> Fill(aJ);

      // count how many jets have pTjet > 17
      if (pTj > 17.)
        ++nJetG17;
      else
        ++nJetL17;

      // grab RE and UE pTcorr spectra
      if (dFj < 0.)      dFj += 2.*M_PI;
      if (dFj > 2.*M_PI) dFj -= 2.*M_PI;
      Double_t fCut = abs(dFj - M_PI);
      if (aJ > Amin) {
        if (fCut < M_PI/4.)
          hPtRE -> Fill(pTcorr);
        else if ((dFj > M_PI/4.) && (dFj < M_PI/2.))
          hPtUE -> Fill(pTcorr);
        else if ((dFj > 3.*M_PI/2.) && (dFj < 7.*M_PI/4.))
          hPtUE -> Fill(pTcorr);
      }

      // constituent loop
      for (Int_t k = 0; k < nCon; k++) {

        Double_t pTc = jets[j].constituents()[k].perp();
        Double_t hC  = jets[j].constituents()[k].pseudorapidity();
        Double_t fC  = jets[j].constituents()[k].phi();
        Double_t eC  = jets[j].constituents()[k].e();

        conPt.push_back(pTc);
        conH.push_back(hC);
        conF.push_back(fC);
        conE.push_back(eC);


        // fill QA histograms
        Double_t dFc = fC - fTrg;
        if (dFc < -1.*M_PI/2.) dFc += 2.*M_PI;
        if (dFc > 3.*M_PI/2.)  dFc -= 2.*M_PI;
        hPtCst -> Fill(pTc);
        hDfCst -> Fill(dFc);
        hHcst  -> Fill(hC);

      }  // end constituent loop


      femto.JetConsPt.push_back(conPt);
      femto.JetConsEta.push_back(conH);
      femto.JetConsPhi.push_back(conF);
      femto.JetConsE.push_back(conE);

      conPt.clear();
      conH.clear();
      conF.clear();
      conE.clear();

    }  // end jet loop


    // fill tree / QA histograms
    fDst     -> Fill();
    hRho     -> Fill(femto.Rho);
    hNjet    -> Fill((Double_t) nJet);
    hNjetG17 -> Fill((Double_t) nJetG17);
    hNjetL17 -> Fill((Double_t) nJetL17);

    // make sure vectors are clear
    tracks.clear();
    jetsCS.clear();
    jets.clear();
    femto.Reset();

  }  // end event loop


  // subtract underlying evt.
  hPtSub -> Add(hPtRE, hPtUE, 1., -1.);

}  // end 'Make()'



void StFemtoDstMaker::Finish() {

  // for checking efficiency
  if (level != 0) hPtEffCheck -> Divide(hPtEffDet, hPtEffPar, 1., 1.);

  // normalize histograms
  const Int_t    nTrg   = hRefmult -> GetEntries();
  const Double_t tspBin = hTSP     -> GetBinWidth(17);
  const Double_t pTbin  = hEtTrg   -> GetBinWidth(17);
  const Double_t pBin   = hPxTrk   -> GetBinWidth(17);
  const Double_t dFbin  = hFtrg    -> GetBinWidth(17);
  const Double_t hBin   = hHtrg    -> GetBinWidth(17);
  const Double_t rBin   = hRho     -> GetBinWidth(17);
  const Double_t aBin   = hAjet    -> GetBinWidth(17);
  hRefmult   -> Scale(1. / nTrg);
  hTSP       -> Scale(1. / nTrg);
  hTSP       -> Scale(1. / tspBin);
  hEtTrg     -> Scale(1. / nTrg);
  hEtTrg     -> Scale(1. / pTbin);
  hFtrg      -> Scale(1. / nTrg);
  hFtrg      -> Scale(1. / dFbin);
  hHtrg      -> Scale(1. / nTrg);
  hHtrg      -> Scale(1. / hBin);
  hRho       -> Scale(1. / nTrg);
  hRho       -> Scale(1. / rBin);
  hNjet      -> Scale(1. / nTrg);
  hNjetG17   -> Scale(1. / nTrg);
  hNjetL17   -> Scale(1. / nTrg);
  hNtrk      -> Scale(1. / nTrg);
  hNtrkG1    -> Scale(1. / nTrg);
  hNtrkG2    -> Scale(1. / nTrg);
  hPtTrk     -> Scale(1. / nTrg);
  hPtTrk     -> Scale(1. / pTbin);
  hPxTrk     -> Scale(1. / nTrg);
  hPxTrk     -> Scale(1. / pBin);
  hPyTrk     -> Scale(1. / nTrg);
  hPyTrk     -> Scale(1. / pBin);
  hPzTrk     -> Scale(1. / nTrg);
  hPzTrk     -> Scale(1. / pBin);
  hDfTrk     -> Scale(1. / nTrg);
  hDfTrk     -> Scale(1. / dFbin);
  hDfTrkG1   -> Scale(1. / nTrg);
  hDfTrkG1   -> Scale(1. / dFbin);
  hDfTrkG2   -> Scale(1. / nTrg);
  hDfTrkG2   -> Scale(1. / dFbin);
  hHtrk      -> Scale(1. / nTrg);
  hHtrk      -> Scale(1. / hBin);
  hPtVsDfTrk -> Scale(1. / nTrg);
  hPtVsDfTrk -> Scale(1. / pTbin);
  hPtVsDfTrk -> Scale(1. / dFbin);
  hPtReco    -> Scale(1. / nTrg);
  hPtReco    -> Scale(1. / pTbin);
  hPtCorr    -> Scale(1. / nTrg);
  hPtCorr    -> Scale(1. / pTbin);
  hPtRE      -> Scale(1. / nTrg);
  hPtRE      -> Scale(1. / pTbin);
  hPtUE      -> Scale(1. / nTrg);
  hPtUE      -> Scale(1. / pTbin);
  hPtSub     -> Scale(1. / nTrg);
  hPtSub     -> Scale(1. / pTbin);
  hDfJet     -> Scale(1. / nTrg);
  hDfJet     -> Scale(1. / dFbin);
  hHjet      -> Scale(1. / nTrg);
  hHjet      -> Scale(1. / hBin);
  hAjet      -> Scale(1. / nTrg);
  hAjet      -> Scale(1. / aBin);
  hPtCst     -> Scale(1. / nTrg);
  hPtCst     -> Scale(1. / pTbin);
  hDfCst     -> Scale(1. / nTrg);
  hDfCst     -> Scale(1. / dFbin);
  hHcst      -> Scale(1. / nTrg);
  hHcst      -> Scale(1. / hBin);
  cout << "  Normalizing QA histograms...\n"
       << "    nTrg=" << nTrg << ", tspBin=" << tspBin << ", pTbin=" << pTbin
       << ", dFbin=" << dFbin << ", hBin=" << hBin << ", rBin="
       << rBin << ", aBin=" << aBin
       << endl;

  // TEST [04.29.2020]
  const Double_t dFpTwidth = hDfTrk_pt021 -> GetBinWidth(17);
  hDfTrk_pt021 -> Scale(1. / nTrg);
  hDfTrk_pt021 -> Scale(1. / dFpTwidth);
  hDfTrk_pt12  -> Scale(1. / nTrg);
  hDfTrk_pt12  -> Scale(1. / dFpTwidth);
  hDfTrk_pt25  -> Scale(1. / nTrg);
  hDfTrk_pt25  -> Scale(1. / dFpTwidth);
  hDfTrk_pt530 -> Scale(1. / nTrg);
  hDfTrk_pt530 -> Scale(1. / dFpTwidth);
  oFile        -> cd();
  hDfTrk_pt021 -> Write();
  hDfTrk_pt12  -> Write();
  hDfTrk_pt25  -> Write();
  hDfTrk_pt530 -> Write();

  // save and close file
  oFile       -> cd("QA");
  hRefmult    -> Write();
  hTSP        -> Write();
  hEtTrg      -> Write();
  hFtrg       -> Write();
  hHtrg       -> Write();
  hRho        -> Write();
  hNjet       -> Write();
  hNjetG17    -> Write();
  hNjetL17    -> Write();
  hNtrk       -> Write();
  hNtrkG1     -> Write();
  hNtrkG2     -> Write();
  hPtTrk      -> Write();
  hPxTrk      -> Write();
  hPyTrk      -> Write();
  hPzTrk      -> Write();
  hDfTrk      -> Write();
  hDfTrkG1    -> Write();
  hDfTrkG2    -> Write();
  hHtrk       -> Write();
  hPtVsDfTrk  -> Write();
  hPtReco     -> Write();
  hPtCorr     -> Write();
  hPtRE       -> Write();
  hPtUE       -> Write();
  hPtSub      -> Write();
  hDfJet      -> Write();
  hHjet       -> Write();
  hAjet       -> Write();
  hPtCst      -> Write();
  hDfCst      -> Write();
  hHcst       -> Write();
  hPtEffPar   -> Write();
  hPtEffDet   -> Write();
  hPtEffCheck -> Write();
  hPtRes1D    -> Write();
  hPtRes2D    -> Write();
  oFile       -> cd();
  fDst        -> Write();
  oFile       -> Close();

  cout << "  Done!\n" << endl;

}  // end 'Finish()'

// End ------------------------------------------------------------------------
