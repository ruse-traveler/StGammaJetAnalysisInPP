// 'StMcJetTreeMaker.utilities.h'
// Derek Anderson, Nihar Sahoo
// 04.30.2018
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets.

#pragma once



void StMcJetTreeMaker::InitializeInputTree(TTree *tree) {

  // set branch addresses and branch pointers
  if (!tree) return;
  _tMc      = tree;
  _fCurrent = -1;
  _tMc -> SetMakeClass(1);
  _tMc -> SetBranchAddress("EventId", &_mcEventId, &_bMcEventId);
  _tMc -> SetBranchAddress("RunId", &_mcRunId, &_bMcRunId);
  _tMc -> SetBranchAddress("NumTrks", &_mcNumTrks, &_bMcNumTrks);
  _tMc -> SetBranchAddress("MuVtxX", &_muVx, &_bMuVx);
  _tMc -> SetBranchAddress("MuVtxY", &_muVy, &_bMuVy);
  _tMc -> SetBranchAddress("MuVtxZ", &_muVz, &_bMuVz);
  _tMc -> SetBranchAddress("McVtxX", &_mcVx, &_bMcVx);
  _tMc -> SetBranchAddress("McVtxY", &_mcVy, &_bMcVy);
  _tMc -> SetBranchAddress("McVtxZ", &_mcVz, &_bMcVz);
  _tMc -> SetBranchAddress("IdTrk", &_mcIdTrk, &_bMcIdTrk);
  _tMc -> SetBranchAddress("IdGeant", &_mcIdGeant, &_bMcIdGeant);
  _tMc -> SetBranchAddress("IdVx", &_mcIdVx, &_bMcIdVx);
  _tMc -> SetBranchAddress("IdVxEnd", &_mcIdVxEnd, &_bMcIdVxEnd);
  _tMc -> SetBranchAddress("IntrVtx", &_mcIntrVtx, &_bMcIntrVtx);
  _tMc -> SetBranchAddress("IsShower", &_mcIsShower, &_bMcIsShower);
  _tMc -> SetBranchAddress("Charge", &_mcCharge, &_bMcCharge);
  _tMc -> SetBranchAddress("Rapidity", &_mcRapidity, &_bMcRapidity);
  _tMc -> SetBranchAddress("Eta", &_mcEta, &_bMcEta);
  _tMc -> SetBranchAddress("Phi", &_mcPhi, &_bMcPhi);
  _tMc -> SetBranchAddress("Px", &_mcPx, &_bMcPx);
  _tMc -> SetBranchAddress("Py", &_mcPy, &_bMcPy);
  _tMc -> SetBranchAddress("Pz", &_mcPz, &_bMcPz);
  _tMc -> SetBranchAddress("Pt", &_mcPt, &_bMcPt);
  _tMc -> SetBranchAddress("Ptot", &_mcPtot, &_bMcPtot);
  _tMc -> SetBranchAddress("Energy", &_mcEnergy, &_bMcEnergy);
  PrintInfo(7);

}  // end 'InitializeInputTree(TTree*)'



void StMcJetTreeMaker::InitializeOutputTree(TTree *tree) {

  _fOutput -> cd();
  _tJet = new TTree("JetTree", "a tree of jets");
  _tJet -> Branch("eventIndex", &_EventIndex, "EventIndex/I");
  _tJet -> Branch("Refmult", &_Refmult, "Refmult/D");
  _tJet -> Branch("RunId", &_RunId, "RunId/I");
  _tJet -> Branch("NJets", &_NJets, "NJets/I");
  _tJet -> Branch("PartonicPt", &_PartonicPt, "PartonicPt/D");
  _tJet -> Branch("TSP", &_TSP, "TSP/D");
  _tJet -> Branch("TrgEta", &_TrgEta, "TrgEta/D");
  _tJet -> Branch("TrgPhi", &_TrgPhi, "TrgPhi/D");
  _tJet -> Branch("TrgEt", &_TrgEt, "TrgEt/D");
  _tJet -> Branch("Rho", &_Rho, "Rho/D");
  _tJet -> Branch("Sigma", &_Sigma, "Sigma/D");
  _tJet -> Branch("Vz", &_Vz, "Vz/D");
  _tJet -> Branch("JetIndex", &_JetIndex);
  _tJet -> Branch("JetNCons", &_JetNCons);
  _tJet -> Branch("JetPt", &_JetPt);
  _tJet -> Branch("JetPtCorr", &_JetPtCorr);
  _tJet -> Branch("JetEta", &_JetEta);
  _tJet -> Branch("JetPhi", &_JetPhi);
  _tJet -> Branch("JetE", &_JetE);
  _tJet -> Branch("JetArea", &_JetArea);
  _tJet -> Branch("JetConsPt", &_JetConsPt);
  _tJet -> Branch("JetConsEta", &_JetConsEta);
  _tJet -> Branch("JetConsPhi", &_JetConsPhi);
  _tJet -> Branch("JetConsE", &_JetConsE);

  // set autosave for every 500 MB
  _tJet -> SetAutoSave(-500000000);
  PrintInfo(8);

}  // end 'InitializeOutputTree(TTree*)'



void StMcJetTreeMaker::InitializeHistograms() {

  const UInt_t   nTrg = 1000000;
  const UInt_t   nTsp = 100;
  const UInt_t   nTrk = 100;
  const UInt_t   nPt  = 240;
  const UInt_t   nDf  = 720;
  const UInt_t   nH   = 200;
  const Double_t trg1 = 0.;
  const Double_t trg2 = 1000000.;
  const Double_t tsp1 = 0.;
  const Double_t tsp2 = 100.;
  const Double_t trk1 = 0.;
  const Double_t trk2 = 100.;
  const Double_t pT1  = -20.;
  const Double_t pT2  = 100.;
  const Double_t dF1  = -1. * TMath::TwoPi();
  const Double_t dF2  = TMath::TwoPi();
  const Double_t h1   = -10.;
  const Double_t h2   = 10.;
  _fOutput -> cd();
  // event QA
  _hEvtQA[0][0] = new TH1D("hNumTrgP", "No. of triggers, #pi^{0}", nTrg, trg1, trg2);
  _hEvtQA[0][1] = new TH1D("hNumTrgG", "No. of triggers, #gamma^{rich}", nTrg, trg1, trg2);
  _hEvtQA[0][2] = new TH1D("hNumTrgH", "No. of triggers, h^{#pm}", nTrg, trg1, trg2);
  _hEvtQA[1][0] = new TH1D("hTrgTspP", "Trigger TSP, #pi^{0}", nTsp, tsp1, tsp2);
  _hEvtQA[1][1] = new TH1D("hTrgTspG", "Trigger TSP, #gamma^{rich}", nTsp, tsp1, tsp2);
  _hEvtQA[1][2] = new TH1D("hTrgTspH", "Trigger TSP, h^{#pm}", nTsp, tsp1, tsp2);
  _hEvtQA[2][0] = new TH1D("hNumTrkP", "No. of primary tracks, #pi^{0}", nTrk, trk1, trk2);
  _hEvtQA[2][1] = new TH1D("hNumTrkG", "No. of primary tracks, #gamma^{rich}", nTrk, trk1, trk2);
  _hEvtQA[2][2] = new TH1D("hNumTrkH", "No. of primary tracks, h^{#pm}", nTrk, trk1, trk2);
  _hEvtQA[3][0] = new TH1D("hTrgEneP", "Trigger energy, #pi^{0}", nPt, pT1, pT2);
  _hEvtQA[3][1] = new TH1D("hTrgEneG", "Trigger energy, #gamma^{rich}", nPt, pT1, pT2);
  _hEvtQA[3][2] = new TH1D("hTrgEneH", "Trigger energy, h^{#pm}", nPt, pT1, pT2);
  _hTrkQA[0][0] = new TH1D("hTrkPtP", "Primary track p_{T}, #pi^{0}", nPt, pT1, pT2);
  _hTrkQA[0][1] = new TH1D("hTrkPtG", "Primary track p_{T}, #gamma^{rich}", nPt, pT1, pT2);
  _hTrkQA[0][2] = new TH1D("hTrkPtH", "Primary track p_{T}, h^{#pm}", nPt, pT1, pT2);
  _hTrkQA[1][0] = new TH1D("hTrkDfP", "Primary track #Delta#varphi, #pi^{0}", nDf, dF1, dF2);
  _hTrkQA[1][1] = new TH1D("hTrkDfG", "Primary track #Delta#varphi, #gamma^{rich}", nDf, dF1, dF2);
  _hTrkQA[1][2] = new TH1D("hTrkDfH", "Primary track #Delta#varphi, h^{#pm}", nDf, dF1, dF2);
  _hTrkQA[2][0] = new TH1D("hTrkEtaP", "Primary track #eta, #pi^{0}", nH, h1, h2);
  _hTrkQA[2][1] = new TH1D("hTrkEtaG", "Primary track #eta, #gamma^{rich}", nH, h1, h2);
  _hTrkQA[2][2] = new TH1D("hTrkEtaH", "Primary track #eta, h^{#pm}", nH, h1, h2);
  _hTrkQA[3][0] = new TH1D("hTrkEneP", "Primary track energy, #pi^{0}", nPt, pT1, pT2);
  _hTrkQA[3][1] = new TH1D("hTrkEneG", "Primary track energy, #gamma^{rich}", nPt, pT1, pT2);
  _hTrkQA[3][2] = new TH1D("hTrkEneH", "Primary track energy, h^{#pm}", nPt, pT1, pT2);
  // tower QA
  _hTwrQA[0][0] = new TH1D("hTwrPtP", "Tower p_{T}, #pi^{0}", nPt, pT1, pT2);
  _hTwrQA[0][1] = new TH1D("hTwrPtG", "Tower p_{T}, #gamma^{rich}", nPt, pT1, pT2);
  _hTwrQA[0][2] = new TH1D("hTwrPtH", "Tower p_{T}, h^{#pm}", nPt, pT1, pT2);
  _hTwrQA[1][0] = new TH1D("hTwrDfP", "Tower #Delta#varphi, #pi^{0}", nDf , dF1, dF2);
  _hTwrQA[1][1] = new TH1D("hTwrDfG", "Tower #Delta#varphi, #gamma^{rich}", nDf , dF1, dF2);
  _hTwrQA[1][2] = new TH1D("hTwrDfH", "Tower #Delta#varphi, h^{#pm}", nDf , dF1, dF2);
  _hTwrQA[2][0] = new TH1D("hTwrEtaP", "Tower #eta, #pi^{0}", nH, h1, h2);
  _hTwrQA[2][1] = new TH1D("hTwrEtaG", "Tower #eta, #gamma^{rich}", nH, h1, h2);
  _hTwrQA[2][2] = new TH1D("hTwrEtaH", "Tower #eta, h^{#pm}", nH, h1, h2);
  _hTwrQA[3][0] = new TH1D("hTwrEneP", "Tower energy, #pi^{0}", nPt, pT1, pT2);
  _hTwrQA[3][1] = new TH1D("hTwrEneG", "Tower energy, #gamma^{rich}", nPt, pT1, pT2);
  _hTwrQA[3][2] = new TH1D("hTwrEneH", "Tower energy, h^{#pm}", nPt, pT1, pT2);
  // jet QA
  _hJetQA[0][0] = new TH1D("hJetPtP", "Jet p_{T}, #pi^{0}", nPt, pT1, pT2);
  _hJetQA[0][1] = new TH1D("hJetPtG", "Jet p_{T}, #gamma^{rich}", nPt, pT1, pT2);
  _hJetQA[0][2] = new TH1D("hJetPtH", "Jet p_{T}, h^{#pm}", nPt, pT1, pT2);
  _hJetQA[1][0] = new TH1D("hJetDfP", "Jet #Delta#varphi, #pi^{0}", nDf, dF1, dF2);
  _hJetQA[1][1] = new TH1D("hJetDfG", "Jet #Delta#varphi, #gamma^{rich}", nDf, dF1, dF2);
  _hJetQA[1][2] = new TH1D("hJetDfH", "Jet #Delta#varphi, h^{#pm}", nDf, dF1, dF2);
  _hJetQA[2][0] = new TH1D("hJetEtaP", "Jet #eta, #pi^{0}", nH, h1, h2);
  _hJetQA[2][1] = new TH1D("hJetEtaG", "Jet #eta, #gamma^{rich}", nH, h1, h2);
  _hJetQA[2][2] = new TH1D("hJetEtaH", "Jet #eta, h^{#pm}", nH, h1, h2);
  _hJetQA[3][0] = new TH1D("hJetEneP", "Jet energy, #pi^{0}", nPt, pT1, pT2);
  _hJetQA[3][1] = new TH1D("hJetEneG", "Jet energy, #gamma^{rich}", nPt, pT1, pT2);
  _hJetQA[3][2] = new TH1D("hJetEneH", "Jet energy, h^{#pm}", nPt, pT1, pT2);
  // errors
  for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
    for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
      _hEvtQA[iHist][iTrg] -> Sumw2();
      _hTrkQA[iHist][iTrg] -> Sumw2();
      _hTwrQA[iHist][iTrg] -> Sumw2();
      _hJetQA[iHist][iTrg] -> Sumw2();
    }
  }
  PrintInfo(9);

  // TEST [04.02.2020]
  const UInt_t   nZtTest(100);
  const UInt_t   nPtTest(100);
  const UInt_t   nNcstTest(100);
  const Double_t zTtest[2]   = {0., 1.};
  const Double_t pTtest[2]   = {0., 50.};
  const Double_t nCstTest[2] = {0., 100.};
  // pi0 test hist.s
  _hZtCstVsPtJetP[0] = new TH1D("hZtCstVsPtJetP_pt021", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (0.2, 1) GeV/c (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetP[1] = new TH1D("hZtCstVsPtJetP_pt12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (1, 2) GeV/c (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetP[2] = new TH1D("hZtCstVsPtJetP_pt25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (2, 5) GeV/c (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetP[3] = new TH1D("hZtCstVsPtJetP_pt510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (5, 10) GeV/c (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetP[4] = new TH1D("hZtCstVsPtJetP_pt10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} > 10 GeV/c (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstP[0] = new TH1D("hZtCstVsNcstP_n12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (1, 2) (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstP[1] = new TH1D("hZtCstVsNcstP_n25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (2, 5) (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstP[2] = new TH1D("hZtCstVsNcstP_n510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (5, 10) (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstP[3] = new TH1D("hZtCstVsNcstP_n10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} > 10 (#pi^{0})", nZtTest, zTtest[0], zTtest[1]);
  _hNcstVsPtJetP    = new TH2D("hNcstVsPtJetP", "N_{cst}^{jet} vs. p_{T}^{jet} (#pi^{0})", nPtTest, pTtest[0], pTtest[1], nNcstTest, nCstTest[0], nCstTest[1]);
  // gam test hist.s
  _hZtCstVsPtJetG[0] = new TH1D("hZtCstVsPtJetG_pt021", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (0.2, 1) GeV/c (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetG[1] = new TH1D("hZtCstVsPtJetG_pt12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (1, 2) GeV/c (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetG[2] = new TH1D("hZtCstVsPtJetG_pt25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (2, 5) GeV/c (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetG[3] = new TH1D("hZtCstVsPtJetG_pt510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (5, 10) GeV/c (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetG[4] = new TH1D("hZtCstVsPtJetG_pt10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} > 10 GeV/c (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstG[0] = new TH1D("hZtCstVsNcstG_n12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (1, 2) (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstG[1] = new TH1D("hZtCstVsNcstG_n25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (2, 5) (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstG[2] = new TH1D("hZtCstVsNcstG_n510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (5, 10) (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstG[3] = new TH1D("hZtCstVsNcstG_n10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} > 10 (#gamma^{rich})", nZtTest, zTtest[0], zTtest[1]);
  _hNcstVsPtJetG    = new TH2D("hNcstVsPtJetG", "N_{cst}^{jet} vs. p_{T}^{jet} (#gamma^{rich})", nPtTest, pTtest[0], pTtest[1], nNcstTest, nCstTest[0], nCstTest[1]);
  // had test hist.s
  _hZtCstVsPtJetH[0] = new TH1D("hZtCstVsPtJetH_pt021", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (0.2, 1) GeV/c (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetH[1] = new TH1D("hZtCstVsPtJetH_pt12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (1, 2) GeV/c (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetH[2] = new TH1D("hZtCstVsPtJetH_pt25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (2, 5) GeV/c (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetH[3] = new TH1D("hZtCstVsPtJetH_pt510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (5, 10) GeV/c (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJetH[4] = new TH1D("hZtCstVsPtJetH_pt10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} > 10 GeV/c (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstH[0] = new TH1D("hZtCstVsNcstH_n12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (1, 2) (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstH[1] = new TH1D("hZtCstVsNcstH_n25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (2, 5) (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstH[2] = new TH1D("hZtCstVsNcstH_n510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (5, 10) (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcstH[3] = new TH1D("hZtCstVsNcstH_n10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} > 10 (h^{#pm})", nZtTest, zTtest[0], zTtest[1]);
  _hNcstVsPtJetH    = new TH2D("hNcstVsPtJetH", "N_{cst}^{jet} vs. p_{T}^{jet} (h^{#pm})", nPtTest, pTtest[0], pTtest[1], nNcstTest, nCstTest[0], nCstTest[1]);
  for (UInt_t iPtTest = 0; iPtTest < 5; iPtTest++) {
    _hZtCstVsPtJetP[iPtTest] -> Sumw2();
    _hZtCstVsPtJetG[iPtTest] -> Sumw2();
    _hZtCstVsPtJetH[iPtTest] -> Sumw2();
  }
  for (UInt_t iNcstTest = 0; iNcstTest < 4; iNcstTest++) {
    _hZtCstVsNcstP[iNcstTest] -> Sumw2();
    _hZtCstVsNcstG[iNcstTest] -> Sumw2();
    _hZtCstVsNcstH[iNcstTest] -> Sumw2();
  }
  _hNcstVsPtJetP -> Sumw2();
  _hNcstVsPtJetG -> Sumw2();
  _hNcstVsPtJetH -> Sumw2();

}  // end 'InitializeHistograms()'



void StMcJetTreeMaker::PrintInfo(const UInt_t code, const UInt_t nEvts, const UInt_t iEvt) {

  switch (code) {
    case 0:
      cout << "\n  Jet-tree maker created!" << endl;
      break;
    case 1:
      cout << "    Input and output set.\n"
           << "      input  -- " << _sInput.Data() << "\n"
           << "      output -- " << _sOutput.Data()
           << endl;
      break;
    case 2:
      cout << "    Event parameters set.\n"
           << "      |rVtx| < " << _rVtxMax << "\n"
           << "      |zVtx| < " << _zVtxMax
           << endl;
      break;
    case 3:
      cout << "    Trigger parameters set.\n"
           << "      |etaTrg| < " << _etaTrgMax << "\n"
           << "      eTtrg    = (" << _eTtrgMin << ", " << _eTtrgMax << ")\n"
           << "      tsp(pi0) = (" << _tspPi0Min << ", " << _tspPi0Max << ")\n"
           << "      tsp(gam) = (" << _tspGamMin << ", " << _tspGamMax << ")"
           << endl;
      break;
    case 4:
      cout << "    Track parameters set.\n"
           << "      |etaTrk| < " << _etaTrkMax << "\n"
           << "      pTtrk    = (" << _pTtrkMin << ", " << _pTtrkMax << ")"
           << endl;
      break;
    case 6:
      cout << "    Jet parameters set.\n"
           << "      jetType    = " << _jetType << "\n"
           << "      nRepeat    = " << _nRepeat << "\n"
           << "      nRemove    = " << _nRemove << "\n"
           << "      rJet       = " << _rJet << "\n"
           << "      aGhost     = " << _aGhost << "\n"
           << "      pTjet      > " << _pTjetMin << "\n"
           << "      |etaGhost| < " << _etaGhostMax << "\n"
           << "      |etaJet|   < " << _etaJetMax << "\n"
           << "      |etaBkgd|  < " << _etaBkgdMax
           << endl;
      break;
    case 7:
      cout << "    Input tree initialized." << endl;
      break;
    case 8:
      cout << "    Output tree initialized." << endl;
      break;
    case 9:
      cout << "    QA histograms initialized." << endl;
      break;
    case 10:
      cout << "    Maker initialized!" << endl;
      break;
    case 11:
      cout << "    Beginning event loop: " << nEvts << " events to process." << endl;
      break;
    case 12:
      if (_isInBatchMode)
        cout << "      processing event " << iEvt + 1 << "/" << nEvts << "..." << endl;
      else {
        cout << "      processing event " << iEvt + 1 << "/" << nEvts << "...\r" << flush;
        if ((iEvt + 1) == nEvts) cout << endl;
      }
      break;
    case 13:
      cout << "    Event loop finished." << endl;
      break;
    case 14:
      cout << "    Histograms normalized." << endl;
      break;
    case 15:
      cout << "    Directories made." << endl;
      break;
    case 16:
      cout << "  Files closed, jet making finished!\n" << endl;
      break;
    case 17:
      cout << "    Adjusting track efficiency by " << _effAdjust << endl;
      break;
    case 18:
      cout << "    Turning off weak decays!" << endl;
      break;
  }

}  // end 'PrintInfo(UInt_t)'



Bool_t StMcJetTreeMaker::IsGoodRunID(const UInt_t runID) {

  // bad run list (pp200 run9)
  const UInt_t badRunList[NBadRun] = {10114082, 10120093, 10159043, 10166054, 10126064, 10128094, 10128102, 10131009, 10131075, 10131087, 10132004, 10135072, 10136036, 10138049, 10140005, 10140011, 10142012, 10142035, 10142093, 10144038, 10144074, 10149008, 10150005, 10151001, 10152010, 10156090, 10157015, 10157053, 10158047, 10160006, 10161006, 10161016, 10161024, 10162007, 10165027, 10165077, 10166024, 10169033, 10170011, 10170029, 10170047, 10171011, 10172054, 10172059, 10172077};

  Bool_t isGoodRun = true;
  for (UInt_t iRun = 0; iRun < NBadRun; iRun++) {
    if (runID == badRunList[iRun]) {
      isGoodRun = false;
      break;
    }
  }
  return isGoodRun;

}  // end 'IsGoodRunID(UInt_t)'



Bool_t StMcJetTreeMaker::IsGoodEvent(const Double_t rVtx, const Double_t zVtx) {

  const Bool_t isInRCut  = (abs(rVtx) < _rVtxMax);
  const Bool_t isInZCut  = (abs(zVtx) < _zVtxMax);
  const Bool_t isGoodEvt = (isInRCut && isInZCut);
  return isGoodEvt;

}  // end 'IsGoodEvent(Double_t, Double_t)'



Bool_t StMcJetTreeMaker::IsGoodTrigger(const Double_t etaTrg, const Double_t eTtrg, const UInt_t idTrg) {

  const Bool_t isInEtaCut = (abs(etaTrg) < _etaTrgMax);
  const Bool_t isInEtCut  = ((eTtrg > _eTtrgMin) && (eTtrg < _eTtrgMax));
  const Bool_t isPi0      = IsPi0(idTrg);
  const Bool_t isGam      = IsGamma(idTrg);;
  const Bool_t isHad      = IsHadron(idTrg);
  const Bool_t isInIdCut  = (isPi0 || isGam || isHad);
  const Bool_t isGoodTrg  = (isInEtaCut && isInEtCut && isInIdCut);
  return isGoodTrg;

}  // end 'IsGoodTrigger(Double_t, Double_t, UInt_t)'



Bool_t StMcJetTreeMaker::IsGoodTrack(const Double_t etaTrk, const Double_t pTtrk) {

  const Bool_t isInEtaCut = (abs(etaTrk) < _etaTrkMax);
  const Bool_t isInPtCut  = ((pTtrk > _pTtrkMin) && (pTtrk < _pTtrkMax));
  const Bool_t isGoodTrk  = (isInEtaCut && isInPtCut);
  return isGoodTrk;

}  // end 'IsGoodTrack(Double_t, Double_t)'



Bool_t StMcJetTreeMaker::IsPi0(const UInt_t idGnt) {

  const Bool_t isPi0 = (idGnt == IdPi0);
  return isPi0;

}  // end 'IsPi0(UInt_t)'



Bool_t StMcJetTreeMaker::IsGamma(const UInt_t idGnt) {

  const Bool_t isGam = (idGnt == IdGamma);
  return isGam;

}  // end 'IsGamma(UInt_t)'



Bool_t StMcJetTreeMaker::IsHadron(const UInt_t idGnt) {

  // consider only charged hadrons
  const UInt_t idHad[NHadIds] = {8, 9, 11, 12, 14, 15, 19, 21, 23, 24, 27, 29, 31, 32, 45, 46};

  // check ID
  Bool_t isHad(false);
  for (UInt_t iHad = 0; iHad < NHadIds; iHad++) {
    if (idGnt == idHad[iHad]) {
      isHad = true;
      break;
    }
  }
  return isHad;

}  // end 'IsHadron(UInt_t)'



Bool_t StMcJetTreeMaker::IsUndecayed(const Int_t idVtx) {

  const Bool_t isIS = (idVtx == 1);
  return isIS;

}  // end 'IsUndecayed(Int_t)'



Bool_t StMcJetTreeMaker::IsFinalState(const Double_t pTpar) {

  const Bool_t isFS = (pTpar > 0);
  return isFS;

}  // end 'IsFinalState(Double_t)'



Long64_t StMcJetTreeMaker::LoadTree(const Long64_t entry) {

  // set the environment to read one entry
  Int_t    fTree;
  Long64_t bytes;
  if (!_tMc)
    bytes = -5;
  else {
    bytes = _tMc -> LoadTree(entry);
    fTree = _tMc -> GetTreeNumber();
    if ((bytes >= 0) && (fTree != _fCurrent))
      _fCurrent = _tMc -> GetTreeNumber();
  }
  return bytes;

}  // end 'LoadTree(Long64_t)'



Long64_t StMcJetTreeMaker::GetEntry(const Long64_t entry) {

  Int_t bytes;
  if (!_tMc)
    bytes = 0;
  else
    bytes = _tMc -> GetEntry(entry);
  return bytes;

}  // end 'GetEntry(Long64_t)'

// End ------------------------------------------------------------------------

