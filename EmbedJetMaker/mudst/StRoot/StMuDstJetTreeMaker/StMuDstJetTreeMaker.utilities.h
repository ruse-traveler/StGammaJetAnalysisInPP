// 'StMuDstJetTreeMaker.utilities.h'
// Derek Anderson, Nihar Sahoo
// 06.23.2017
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets.

#pragma once



void StMuDstJetTreeMaker::InitializeInputTree(TTree *tree) {

  // set branch addresses and branch pointers
  if (!tree) return;
  _tFemto   = tree;
  _fCurrent = -1;
  _tFemto -> SetMakeClass(1);
  _tFemto -> SetBranchAddress("fUniqueID", &_fUniqueID, &_bEventList_fUniqueID);
  _tFemto -> SetBranchAddress("fBits", &_fBits, &_bEventList_fBits);
  _tFemto -> SetBranchAddress("runNumber", &_runNumber, &_bEventList_runNumber);
  _tFemto -> SetBranchAddress("eventNumber", &_eventNumber, &_bEventList_eventNumber);
  _tFemto -> SetBranchAddress("trigID", &_trigID, &_bEventList_trigID);
  _tFemto -> SetBranchAddress("nGlobalTracks", &_nGlobalTracks, &_bEventList_nGlobalTracks);
  _tFemto -> SetBranchAddress("nPrimaryTracks", &_nPrimaryTracks, &_bEventList_nPrimaryTracks);
  _tFemto -> SetBranchAddress("refMult", &_refMult, &_bEventList_refMult);
  _tFemto -> SetBranchAddress("vpdVz", &_vpdVz, &_bEventList_vpdVz);
  _tFemto -> SetBranchAddress("xVertex", &_xVertex, &_bEventList_xVertex);
  _tFemto -> SetBranchAddress("yVertex", &_yVertex, &_bEventList_yVertex);
  _tFemto -> SetBranchAddress("zVertex", &_zVertex, &_bEventList_zVertex);
  _tFemto -> SetBranchAddress("bbcZVertex", &_bbcZVertex, &_bEventList_bbcZVertex);
  _tFemto -> SetBranchAddress("zdcCoincidenceRate", &_zdcCoincidenceRate, &_bEventList_zdcCoincidenceRate);
  _tFemto -> SetBranchAddress("bbcCoincidenceRate", &_bbcCoincidenceRate, &_bEventList_bbcCoincidenceRate);
  _tFemto -> SetBranchAddress("backgroundRate", &_backgroundRate, &_bEventList_backgroundRate);
  _tFemto -> SetBranchAddress("bbcBlueBackgroundRate", &_bbcBlueBackgroundRate, &_bEventList_bbcBlueBackgroundRate);
  _tFemto -> SetBranchAddress("bbcYellowBackgroundRate", &_bbcYellowBackgroundRate, &_bEventList_bbcYellowBackgroundRate);
  _tFemto -> SetBranchAddress("refMultPos", &_refMultPos, &_bEventList_refMultPos);
  _tFemto -> SetBranchAddress("refMultNeg", &_refMultNeg, &_bEventList_refMultNeg);
  _tFemto -> SetBranchAddress("bTOFTrayMultiplicity", &_bTOFTrayMultiplicity, &_bEventList_bTOFTrayMultiplicity);
  _tFemto -> SetBranchAddress("nVerticies", &_nVerticies, &_bEventList_nVerticies);
  _tFemto -> SetBranchAddress("MagF", &_MagF, &_bEventList_MagF);
  _tFemto -> SetBranchAddress("VrtxRank", &_VrtxRank, &_bEventList_VrtxRank);
  _tFemto -> SetBranchAddress("FlagEvent_TrgTrkMisMtch", &_FlagEvent_TrgTrkMisMtch, &_bEventList_FlagEvent_TrgTrkMisMtch);
  _tFemto -> SetBranchAddress("Etsp", &_Etsp, &_bEventList_Etsp);
  _tFemto -> SetBranchAddress("ETwrdidT", &_ETwrdidT, &_bEventList_ETwrdidT);
  _tFemto -> SetBranchAddress("ETwradc11", &_ETwradc11, &_bEventList_ETwradc11);
  _tFemto -> SetBranchAddress("ETwreneT0", &_ETwreneT0, &_bEventList_ETwreneT0);
  _tFemto -> SetBranchAddress("ETwreT", &_ETwreT, &_bEventList_ETwreT);
  _tFemto -> SetBranchAddress("ETwrENET0", &_ETwrENET0, &_bEventList_ETwrENET0);
  _tFemto -> SetBranchAddress("ETwrphT", &_ETwrphT, &_bEventList_ETwrphT);
  _tFemto -> SetBranchAddress("ETwrPTower", &_ETwrPTower, &_bEventList_ETwrPTower);
  _tFemto -> SetBranchAddress("ETwrpidTower", &_ETwrpidTower, &_bEventList_ETwrpidTower);
  _tFemto -> SetBranchAddress("ETwrmoduleT", &_ETwrmoduleT, &_bEventList_ETwrmoduleT);
  _tFemto -> SetBranchAddress("EClustEneT0", &_EClustEneT0, &_bEventList_EClustEneT0);
  _tFemto -> SetBranchAddress("EClustetav1", &_EClustetav1, &_bEventList_EClustetav1);
  _tFemto -> SetBranchAddress("EClustphiv1", &_EClustphiv1, &_bEventList_EClustphiv1);
  _tFemto -> SetBranchAddress("EEstrpen01", &_EEstrpen01, &_bEventList_EEstrpen01);
  _tFemto -> SetBranchAddress("EEstrpen02", &_EEstrpen02, &_bEventList_EEstrpen02);
  _tFemto -> SetBranchAddress("EEstrpen03", &_EEstrpen03, &_bEventList_EEstrpen03);
  _tFemto -> SetBranchAddress("EEstrpen0", &_EEstrpen0, &_bEventList_EEstrpen0);
  _tFemto -> SetBranchAddress("EEstrpen1", &_EEstrpen1, &_bEventList_EEstrpen1);
  _tFemto -> SetBranchAddress("EEstrpen2", &_EEstrpen2, &_bEventList_EEstrpen2);
  _tFemto -> SetBranchAddress("EEstrpen3", &_EEstrpen3, &_bEventList_EEstrpen3);
  _tFemto -> SetBranchAddress("EEstrpen4", &_EEstrpen4, &_bEventList_EEstrpen4);
  _tFemto -> SetBranchAddress("EEstrpen5", &_EEstrpen5, &_bEventList_EEstrpen5);
  _tFemto -> SetBranchAddress("EEstrpen6", &_EEstrpen6, &_bEventList_EEstrpen6);
  _tFemto -> SetBranchAddress("EEstrpen7", &_EEstrpen7, &_bEventList_EEstrpen7);
  _tFemto -> SetBranchAddress("EEstrpen8", &_EEstrpen8, &_bEventList_EEstrpen8);
  _tFemto -> SetBranchAddress("EEstrpen9", &_EEstrpen9, &_bEventList_EEstrpen9);
  _tFemto -> SetBranchAddress("EEstrpen10", &_EEstrpen10, &_bEventList_EEstrpen10);
  _tFemto -> SetBranchAddress("EEstrpen11", &_EEstrpen11, &_bEventList_EEstrpen11);
  _tFemto -> SetBranchAddress("EEstrpen12", &_EEstrpen12, &_bEventList_EEstrpen12);
  _tFemto -> SetBranchAddress("EEstrpen13", &_EEstrpen13, &_bEventList_EEstrpen13);
  _tFemto -> SetBranchAddress("EEstrpen14", &_EEstrpen14, &_bEventList_EEstrpen14);
  _tFemto -> SetBranchAddress("EEstrpen15", &_EEstrpen15, &_bEventList_EEstrpen15);
  _tFemto -> SetBranchAddress("ETwrdidE", &_ETwrdidE, &_bEventList_ETwrdidE);
  _tFemto -> SetBranchAddress("EPstripenp01", &_EPstripenp01, &_bEventList_EPstripenp01);
  _tFemto -> SetBranchAddress("EPstripenp02", &_EPstripenp02, &_bEventList_EPstripenp02);
  _tFemto -> SetBranchAddress("EPstripenp03", &_EPstripenp03, &_bEventList_EPstripenp03);
  _tFemto -> SetBranchAddress("EPstripenp0", &_EPstripenp0, &_bEventList_EPstripenp0);
  _tFemto -> SetBranchAddress("EPstripenp1", &_EPstripenp1, &_bEventList_EPstripenp1);
  _tFemto -> SetBranchAddress("EPstripenp2", &_EPstripenp2, &_bEventList_EPstripenp2);
  _tFemto -> SetBranchAddress("EPstripenp3", &_EPstripenp3, &_bEventList_EPstripenp3);
  _tFemto -> SetBranchAddress("EPstripenp4", &_EPstripenp4, &_bEventList_EPstripenp4);
  _tFemto -> SetBranchAddress("EPstripenp5", &_EPstripenp5, &_bEventList_EPstripenp5);
  _tFemto -> SetBranchAddress("EPstripenp6", &_EPstripenp6, &_bEventList_EPstripenp6);
  _tFemto -> SetBranchAddress("EPstripenp7", &_EPstripenp7, &_bEventList_EPstripenp7);
  _tFemto -> SetBranchAddress("EPstripenp8", &_EPstripenp8, &_bEventList_EPstripenp8);
  _tFemto -> SetBranchAddress("EPstripenp9", &_EPstripenp9, &_bEventList_EPstripenp9);
  _tFemto -> SetBranchAddress("EPstripenp10", &_EPstripenp10, &_bEventList_EPstripenp10);
  _tFemto -> SetBranchAddress("EPstripenp11", &_EPstripenp11, &_bEventList_EPstripenp11);
  _tFemto -> SetBranchAddress("EPstripenp12", &_EPstripenp12, &_bEventList_EPstripenp12);
  _tFemto -> SetBranchAddress("EPstripenp13", &_EPstripenp13, &_bEventList_EPstripenp13);
  _tFemto -> SetBranchAddress("EPstripenp14", &_EPstripenp14, &_bEventList_EPstripenp14);
  _tFemto -> SetBranchAddress("EPstripenp15", &_EPstripenp15, &_bEventList_EPstripenp15);
  _tFemto -> SetBranchAddress("EclustEnnq1", &_EclustEnnq1, &_bEventList_EclustEnnq1);
  _tFemto -> SetBranchAddress("EclustEnnq20", &_EclustEnnq20, &_bEventList_EclustEnnq20);
  _tFemto -> SetBranchAddress("EclustEnnq19", &_EclustEnnq19, &_bEventList_EclustEnnq19);
  _tFemto -> SetBranchAddress("EclustEnpq1", &_EclustEnpq1, &_bEventList_EclustEnpq1);
  _tFemto -> SetBranchAddress("EclustEnpq20", &_EclustEnpq20, &_bEventList_EclustEnpq20);
  _tFemto -> SetBranchAddress("EclustEnpq19", &_EclustEnpq19, &_bEventList_EclustEnpq19);
  _tFemto -> SetBranchAddress("EclustEnpq21", &_EclustEnpq21, &_bEventList_EclustEnpq21);
  _tFemto -> SetBranchAddress("PrimaryTrackArray", &_PrimaryTrackArray_, &_bEventList_PrimaryTrackArray_);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.fUniqueID", _PrimaryTrackArray_fUniqueID, &_bPrimaryTrackArray_fUniqueID);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.fBits", _PrimaryTrackArray_fBits, &_bPrimaryTrackArray_fBits);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nHitsFit", _PrimaryTrackArray_nHitsFit, &_bPrimaryTrackArray_nHitsFit);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nHitsPoss", _PrimaryTrackArray_nHitsPoss, &_bPrimaryTrackArray_nHitsPoss);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.trackFlag", _PrimaryTrackArray_trackFlag, &_bPrimaryTrackArray_trackFlag);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.pZ", _PrimaryTrackArray_pZ, &_bPrimaryTrackArray_pZ);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.pX", _PrimaryTrackArray_pX, &_bPrimaryTrackArray_pX);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.pY", _PrimaryTrackArray_pY, &_bPrimaryTrackArray_pY);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.pT", _PrimaryTrackArray_pT, &_bPrimaryTrackArray_pT);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.dEdx", _PrimaryTrackArray_dEdx, &_bPrimaryTrackArray_dEdx);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.charge", _PrimaryTrackArray_charge, &_bPrimaryTrackArray_charge);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.tofBeta", _PrimaryTrackArray_tofBeta, &_bPrimaryTrackArray_tofBeta);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.eta", _PrimaryTrackArray_eta, &_bPrimaryTrackArray_eta);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.phi", _PrimaryTrackArray_phi, &_bPrimaryTrackArray_phi);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nSigElectron", _PrimaryTrackArray_nSigElectron, &_bPrimaryTrackArray_nSigElectron);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nSigPion", _PrimaryTrackArray_nSigPion, &_bPrimaryTrackArray_nSigPion);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nSigKaon", _PrimaryTrackArray_nSigKaon, &_bPrimaryTrackArray_nSigKaon);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nSigProton", _PrimaryTrackArray_nSigProton, &_bPrimaryTrackArray_nSigProton);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.dcag", _PrimaryTrackArray_dcag, &_bPrimaryTrackArray_dcag);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.nHits", _PrimaryTrackArray_nHits, &_bPrimaryTrackArray_nHits);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.dEdxHits", _PrimaryTrackArray_dEdxHits, &_bPrimaryTrackArray_dEdxHits);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.firstZPoint", _PrimaryTrackArray_firstZPoint, &_bPrimaryTrackArray_firstZPoint);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.lastZPoint", _PrimaryTrackArray_lastZPoint, &_bPrimaryTrackArray_lastZPoint);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.tofSigElectron", _PrimaryTrackArray_tofSigElectron, &_bPrimaryTrackArray_tofSigElectron);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.tofSigPion", _PrimaryTrackArray_tofSigPion, &_bPrimaryTrackArray_tofSigPion);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.tofSigKaon", _PrimaryTrackArray_tofSigKaon, &_bPrimaryTrackArray_tofSigKaon);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.tofSigProton", _PrimaryTrackArray_tofSigProton, &_bPrimaryTrackArray_tofSigProton);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.timeOfflight", _PrimaryTrackArray_timeOfflight, &_bPrimaryTrackArray_timeOfflight);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.pathLength", _PrimaryTrackArray_pathLength, &_bPrimaryTrackArray_pathLength);
  _tFemto -> SetBranchAddress("PrimaryTrackArray.trkIndex", _PrimaryTrackArray_trkIndex, &_bPrimaryTrackArray_trkIndex);
  _tFemto -> SetBranchAddress("TowerArray", &_TowerArray_, &_bEventList_TowerArray_);
  _tFemto -> SetBranchAddress("TowerArray.fUniqueID", _TowerArray_fUniqueID, &_bTowerArray_fUniqueID);
  _tFemto -> SetBranchAddress("TowerArray.fBits", _TowerArray_fBits, &_bTowerArray_fBits);
  _tFemto -> SetBranchAddress("TowerArray.TwrId", _TowerArray_TwrId, &_bTowerArray_TwrId);
  _tFemto -> SetBranchAddress("TowerArray.TwrEng", _TowerArray_TwrEng, &_bTowerArray_TwrEng);
  _tFemto -> SetBranchAddress("TowerArray.TwrEta", _TowerArray_TwrEta, &_bTowerArray_TwrEta);
  _tFemto -> SetBranchAddress("TowerArray.TwrPhi", _TowerArray_TwrPhi, &_bTowerArray_TwrPhi);
  _tFemto -> SetBranchAddress("TowerArray.TwrADC", _TowerArray_TwrADC, &_bTowerArray_TwrADC);
  _tFemto -> SetBranchAddress("TowerArray.TwrPed", _TowerArray_TwrPed, &_bTowerArray_TwrPed);
  _tFemto -> SetBranchAddress("TowerArray.TwrRMS", _TowerArray_TwrRMS, &_bTowerArray_TwrRMS);
  _tFemto -> SetBranchAddress("TowerArray.TwrMatchIdnex", _TowerArray_TwrMatchIdnex, &_bTowerArray_TwrMatchIdnex);
  _tFemto -> SetBranchAddress("TowerArray.NoOfmatchedTrk", _TowerArray_NoOfmatchedTrk, &_bTowerArray_NoOfmatchedTrk);
  _tFemto -> SetBranchAddress("TowerArray.TwrMatchP", _TowerArray_TwrMatchP, &_bTowerArray_TwrMatchP);
  _tFemto -> SetBranchAddress("TowerArray.TwrPx", _TowerArray_TwrPx, &_bTowerArray_TwrPx);
  _tFemto -> SetBranchAddress("TowerArray.TwrPy", _TowerArray_TwrPy, &_bTowerArray_TwrPy);
  _tFemto -> SetBranchAddress("TowerArray.TwrPz", _TowerArray_TwrPz, &_bTowerArray_TwrPz);
  _tFemto -> SetBranchAddress("TowerArray.fNAssocTracks", _TowerArray_fNAssocTracks, &_bTowerArray_fNAssocTracks);
  _tFemto -> SetBranchAddress("TowerArray.fMatchedTracksArray_[10]", _TowerArray_fMatchedTracksArray_, &_bTowerArray_fMatchedTracksArray_);
  _tFemto -> SetBranchAddress("TowerArray.fMatchedTracksArray_P[10]", _TowerArray_fMatchedTracksArray_P, &_bTowerArray_fMatchedTracksArray_P);
  PrintInfo(7);

}  // end 'InitializeInputTree(TTree*)'



void StMuDstJetTreeMaker::InitializeOutputTree(TTree *tree) {

  _fOutput -> cd();
  _tJet = new TTree("JetTree", "a tree of jets");
  _tJet -> Branch("eventIndex", &_EventIndex, "EventIndex/I");
  _tJet -> Branch("Refmult", &_Refmult, "Refmult/D");
  _tJet -> Branch("NJets", &_NJets, "NJets/I");
  _tJet -> Branch("RunId", &_RunId, "RunId/I");
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

  // set autosave for every 5 MB
  _tJet -> SetAutoSave(-500000000);
  PrintInfo(8);

}  // end 'InitializeOutputTree(TTree*)'



void StMuDstJetTreeMaker::InitializeHistograms() {

  const UInt_t   nTrg = 1000000;
  const UInt_t   nTsp = 100;
  const UInt_t   nTrk = 100;
  const UInt_t   nPt  = 240;
  const UInt_t   nDf  = 720;
  const UInt_t   nH   = 200;
  const Double_t trg1 = 0.;
  const Double_t trg2 = 1000000.;
  const Double_t tsp1 = 0.;
  const Double_t tsp2 = 1.;
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
  _hEvtQA[0][3] = new TH1D("hNumTrg0", "No. of triggers, no trg.", nTrg, trg1, trg2);
  _hEvtQA[1][0] = new TH1D("hTrgTspP", "Trigger TSP, #pi^{0}", nTsp, tsp1, tsp2);
  _hEvtQA[1][1] = new TH1D("hTrgTspG", "Trigger TSP, #gamma^{rich}", nTsp, tsp1, tsp2);
  _hEvtQA[1][2] = new TH1D("hTrgTspH", "Trigger TSP, h^{#pm}", nTsp, tsp1, tsp2);
  _hEvtQA[1][3] = new TH1D("hTrgTsp0", "Trigger TSP, no trg.", nTsp, tsp1, tsp2);
  _hEvtQA[2][0] = new TH1D("hNumTrkP", "No. of primary tracks, #pi^{0}", nTrk, trk1, trk2);
  _hEvtQA[2][1] = new TH1D("hNumTrkG", "No. of primary tracks, #gamma^{rich}", nTrk, trk1, trk2);
  _hEvtQA[2][2] = new TH1D("hNumTrkH", "No. of primary tracks, h^{#pm}", nTrk, trk1, trk2);
  _hEvtQA[2][3] = new TH1D("hNumTrk0", "No. of primary tracks, no trg.", nTrk, trk1, trk2);
  _hEvtQA[3][0] = new TH1D("hTrgEneP", "Trigger energy, #pi^{0}", nPt, pT1, pT2);
  _hEvtQA[3][1] = new TH1D("hTrgEneG", "Trigger energy, #gamma^{rich}", nPt, pT1, pT2);
  _hEvtQA[3][2] = new TH1D("hTrgEneH", "Trigger energy, h^{#pm}", nPt, pT1, pT2);
  _hEvtQA[3][3] = new TH1D("hTrgEne0", "Trigger energy, no trg.", nPt, pT1, pT2);
  // track QA
  _hTrkQA[0][0] = new TH1D("hTrkPtP", "Primary track p_{T}, #pi^{0}", nPt, pT1, pT2);
  _hTrkQA[0][1] = new TH1D("hTrkPtG", "Primary track p_{T}, #gamma^{rich}", nPt, pT1, pT2);
  _hTrkQA[0][2] = new TH1D("hTrkPtH", "Primary track p_{T}, h^{#pm}", nPt, pT1, pT2);
  _hTrkQA[0][3] = new TH1D("hTrkPt0", "Primary track p_{T}, no trg.", nPt, pT1, pT2);
  _hTrkQA[1][0] = new TH1D("hTrkDfP", "Primary track #Delta#varphi, #pi^{0}", nDf, dF1, dF2);
  _hTrkQA[1][1] = new TH1D("hTrkDfG", "Primary track #Delta#varphi, #gamma^{rich}", nDf, dF1, dF2);
  _hTrkQA[1][2] = new TH1D("hTrkDfH", "Primary track #Delta#varphi, h^{#pm}", nDf, dF1, dF2);
  _hTrkQA[1][3] = new TH1D("hTrkDf0", "Primary track #Delta#varphi, no trg.", nDf, dF1, dF2);
  _hTrkQA[2][0] = new TH1D("hTrkEtaP", "Primary track #eta, #pi^{0}", nH, h1, h2);
  _hTrkQA[2][1] = new TH1D("hTrkEtaG", "Primary track #eta, #gamma^{rich}", nH, h1, h2);
  _hTrkQA[2][2] = new TH1D("hTrkEtaH", "Primary track #eta, h^{#pm}", nH, h1, h2);
  _hTrkQA[2][3] = new TH1D("hTrkEta0", "Primary track #eta, no trg.", nH, h1, h2);
  _hTrkQA[3][0] = new TH1D("hTrkEneP", "Primary track energy, #pi^{0}", nPt, pT1, pT2);
  _hTrkQA[3][1] = new TH1D("hTrkEneG", "Primary track energy, #gamma^{rich}", nPt, pT1, pT2);
  _hTrkQA[3][2] = new TH1D("hTrkEneH", "Primary track energy, h^{#pm}", nPt, pT1, pT2);
  _hTrkQA[3][3] = new TH1D("hTrkEne0", "Primary track energy, no trg.", nPt, pT1, pT2);
  // tower QA
  _hTwrQA[0][0] = new TH1D("hTwrPtP", "Tower p_{T}, #pi^{0}", nPt, pT1, pT2);
  _hTwrQA[0][1] = new TH1D("hTwrPtG", "Tower p_{T}, #gamma^{rich}", nPt, pT1, pT2);
  _hTwrQA[0][2] = new TH1D("hTwrPtH", "Tower p_{T}, h^{#pm}", nPt, pT1, pT2);
  _hTwrQA[0][3] = new TH1D("hTwrPt0", "Tower p_{T}, no trg.", nPt, pT1, pT2);
  _hTwrQA[1][0] = new TH1D("hTwrDfP", "Tower #Delta#varphi, #pi^{0}", nDf , dF1, dF2);
  _hTwrQA[1][1] = new TH1D("hTwrDfG", "Tower #Delta#varphi, #gamma^{rich}", nDf , dF1, dF2);
  _hTwrQA[1][2] = new TH1D("hTwrDfH", "Tower #Delta#varphi, h^{#pm}", nDf , dF1, dF2);
  _hTwrQA[1][3] = new TH1D("hTwrDf0", "Tower #Delta#varphi, no trg.", nDf , dF1, dF2);
  _hTwrQA[2][0] = new TH1D("hTwrEtaP", "Tower #eta, #pi^{0}", nH, h1, h2);
  _hTwrQA[2][1] = new TH1D("hTwrEtaG", "Tower #eta, #gamma^{rich}", nH, h1, h2);
  _hTwrQA[2][2] = new TH1D("hTwrEtaH", "Tower #eta, h^{#pm}", nH, h1, h2);
  _hTwrQA[2][3] = new TH1D("hTwrEta0", "Tower #eta, no trg.", nH, h1, h2);
  _hTwrQA[3][0] = new TH1D("hTwrEneP", "Tower energy, #pi^{0}", nPt, pT1, pT2);
  _hTwrQA[3][1] = new TH1D("hTwrEneG", "Tower energy, #gamma^{rich}", nPt, pT1, pT2);
  _hTwrQA[3][2] = new TH1D("hTwrEneH", "Tower energy, h^{#pm}", nPt, pT1, pT2);
  _hTwrQA[3][3] = new TH1D("hTwrEne0", "Tower energy, no trg.", nPt, pT1, pT2);
  // jet QA
  _hJetQA[0][0] = new TH1D("hJetPtP", "Jet p_{T}, #pi^{0}", nPt, pT1, pT2);
  _hJetQA[0][1] = new TH1D("hJetPtG", "Jet p_{T}, #gamma^{rich}", nPt, pT1, pT2);
  _hJetQA[0][2] = new TH1D("hJetPtH", "Jet p_{T}, h^{#pm}", nPt, pT1, pT2);
  _hJetQA[0][3] = new TH1D("hJetPt0", "Jet p_{T}, no trg.", nPt, pT1, pT2);
  _hJetQA[1][0] = new TH1D("hJetDfP", "Jet #Delta#varphi, #pi^{0}", nDf, dF1, dF2);
  _hJetQA[1][1] = new TH1D("hJetDfG", "Jet #Delta#varphi, #gamma^{rich}", nDf, dF1, dF2);
  _hJetQA[1][2] = new TH1D("hJetDfH", "Jet #Delta#varphi, h^{#pm}", nDf, dF1, dF2);
  _hJetQA[1][3] = new TH1D("hJetDf0", "Jet #Delta#varphi, no trg.", nDf, dF1, dF2);
  _hJetQA[2][0] = new TH1D("hJetEtaP", "Jet #eta, #pi^{0}", nH, h1, h2);
  _hJetQA[2][1] = new TH1D("hJetEtaG", "Jet #eta, #gamma^{rich}", nH, h1, h2);
  _hJetQA[2][2] = new TH1D("hJetEtaH", "Jet #eta, h^{#pm}", nH, h1, h2);
  _hJetQA[2][3] = new TH1D("hJetEta0", "Jet #eta, no trg.", nH, h1, h2);
  _hJetQA[3][0] = new TH1D("hJetEneP", "Jet energy, #pi^{0}", nPt, pT1, pT2);
  _hJetQA[3][1] = new TH1D("hJetEneG", "Jet energy, #gamma^{rich}", nPt, pT1, pT2);
  _hJetQA[3][2] = new TH1D("hJetEneH", "Jet energy, h^{#pm}", nPt, pT1, pT2);
  _hJetQA[3][3] = new TH1D("hJetEne0", "Jet energy, no trg.", nPt, pT1, pT2);
  // event pT
  _hEvtPt[0][0] = new TH1D("hEvtPtPP", "Event p_{T} (all primary tracks), #pi^{0}", nPt, pT1, pT2);
  _hEvtPt[0][1] = new TH1D("hEvtPtQP", "Event p_{T} (with QA cuts), #pi^{0}", nPt, pT1, pT2);
  _hEvtPt[1][0] = new TH1D("hEvtPtPG", "Event p_{T} (all primary tracks), #gamma^{rich}", nPt, pT1, pT2);
  _hEvtPt[1][1] = new TH1D("hEvtPtQG", "Event p_{T} (with QA cuts), #gamma^{rich}", nPt, pT1, pT2);
  _hEvtPt[2][0] = new TH1D("hEvtPtPH", "Event p_{T} (all primary tracks), h^{#pm}", nPt, pT1, pT2);
  _hEvtPt[2][1] = new TH1D("hEvtPtQH", "Event p_{T} (with QA cuts), h^{#pm}", nPt, pT1, pT2);
  _hEvtPt[3][0] = new TH1D("hEvtPtP0", "Event p_{T} (all primary tracks), no trg.", nPt, pT1, pT2);
  _hEvtPt[3][1] = new TH1D("hEvtPtQ0", "Event p_{T} (with QA cuts), no trg.", nPt, pT1, pT2);
  // errors
  for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
    for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
      _hEvtQA[iHist][iTrg] -> Sumw2();
      _hTrkQA[iHist][iTrg] -> Sumw2();
      _hTwrQA[iHist][iTrg] -> Sumw2();
      _hJetQA[iHist][iTrg] -> Sumw2();
    }
  }
  for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
    for (UInt_t iTrk = 0; iTrk < NTrkTypes; iTrk++) {
      _hEvtPt[iTrg][iTrk] -> Sumw2();
    }
  }
  PrintInfo(9);

  // TEST [04.04.2020]
  const UInt_t   nZtTest(100);
  const UInt_t   nPtTest(100);
  const UInt_t   nNcstTest(100);
  const Double_t zTtest[2]   = {0., 1.};
  const Double_t pTtest[2]   = {0., 50.};
  const Double_t nCstTest[2] = {0., 100.};
  _hZtCstVsPtJet[0] = new TH1D("hZtCstVsPtJet_pt021", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (0.2, 1) GeV/c", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJet[1] = new TH1D("hZtCstVsPtJet_pt12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (1, 2) GeV/c", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJet[2] = new TH1D("hZtCstVsPtJet_pt25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (2, 5) GeV/c", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJet[3] = new TH1D("hZtCstVsPtJet_pt510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} #in (5, 10) GeV/c", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsPtJet[4] = new TH1D("hZtCstVsPtJet_pt10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, p_{T}^{jet} > 10 GeV/c", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcst[0] = new TH1D("hZtCstVsNcst_n12", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (1, 2)", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcst[1] = new TH1D("hZtCstVsNcst_n25", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (2, 5)", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcst[2] = new TH1D("hZtCstVsNcst_n510", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} #in (5, 10)", nZtTest, zTtest[0], zTtest[1]);
  _hZtCstVsNcst[3] = new TH1D("hZtCstVsNcst_n10inf", "z_{T}^{cst} = p_{T}^{cst} / p_{T}^{jet}, N_{cst}^{jet} > 10", nZtTest, zTtest[0], zTtest[1]);
  _hNcstVsPtJet    = new TH2D("hNcstVsPtJet", "N_{cst}^{jet} vs. p_{T}^{jet}", nPtTest, pTtest[0], pTtest[1], nNcstTest, nCstTest[0], nCstTest[1]);
  for (UInt_t iPtTest = 0; iPtTest < 5; iPtTest++) {
    _hZtCstVsPtJet[iPtTest] -> Sumw2();
  }
  for (UInt_t iNcstTest = 0; iNcstTest < 4; iNcstTest++) {
    _hZtCstVsNcst[iNcstTest] -> Sumw2();
  }
  _hNcstVsPtJet -> Sumw2();

}  // end 'InitializeHistograms()'



void StMuDstJetTreeMaker::PrintInfo(const UInt_t code, const UInt_t nEvts, const UInt_t iEvt) {

  switch (code) {
    case 0:
      cout << "\n  Jet-tree maker created!" << endl;
      break;
    case 1:
      cout << "    Input and output set.\n"
           << "      input    -- " << _sInput.Data() << "\n"
           << "      output   -- " << _sOutput.Data() << "\n"
           << "      pTparton -- " << _PartonicPt
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
           << "      adc      <= " << _adcMax << "\n"
           << "      eEta     > " << _eEtaMin << ", ePhi > " << _ePhiMin << "\n"
           << "      pProj    < " << _pProjMax << "\n"
           << "      |etaTrg| < " << _etaTrgMax << "\n"
           << "      eTtrg    = (" << _eTtrgMin << ", " << _eTtrgMax << ")\n"
           << "      tsp(pi0) = (" << _tspPi0Min << ", " << _tspPi0Max << ")\n"
           << "      tsp(gam) = (" << _tspGamMin << ", " << _tspGamMax << ")"
           << endl;
      break;
    case 4:
      cout << "    Track parameters set.\n"
           << "      nFit     > " << _nFitMin << ", rFit > " << _rFitMin << "\n"
           << "      dca      < " << _dcaMax << "\n"
           << "      |etaTrk| < " << _etaTrkMax << "\n"
           << "      pTtrk    = (" << _pTtrkMin << ", " << _pTtrkMax << ")"
           << endl;
      break;
    case 5:
      cout << "    Tower parameters set.\n"
           << "      |etaTwr|     < " << _etaTwrMax << "\n"
           << "      eTwr(uncorr) = (" << _eTwrMin << ", " << _eTwrMax << ")\n"
           << "      eTwr(corr)   = (" << _eCorrMin << ", " << _eCorrMax << ")"
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
  }

}  // end 'PrintInfo(UInt_t)'



Bool_t StMuDstJetTreeMaker::IsGoodRunID(const UInt_t runID) {

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



Bool_t StMuDstJetTreeMaker::IsGoodTowerID(const UInt_t twrID) {

  // hot tower list (pp200 run9)
  const UInt_t hotTwrList[NHotTwr] = {34, 106, 113, 160, 266, 267, 275, 280, 282, 286, 287, 293, 410, 504, 533, 541, 555, 561, 562, 594, 615, 616, 629, 633, 637, 638, 647, 650, 653, 657, 671, 673, 743, 789, 790, 791, 792, 806, 809, 810, 811, 812, 813, 814, 821, 822, 823, 824, 829, 830, 831, 832, 837, 841, 842, 843, 844, 846, 849, 850, 851, 852, 857, 875, 897, 899, 903, 939, 953, 954, 956, 993, 1026, 1046, 1048, 1080, 1081, 1100, 1125, 1130, 1132, 1180, 1197, 1198, 1199, 1200, 1207, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1237, 1238, 1240, 1241, 1242, 1243, 1244, 1257, 1258, 1259, 1260, 1312, 1348, 1353, 1354, 1388, 1407, 1409, 1434, 1448, 1537, 1567, 1574, 1597, 1612, 1654, 1668, 1713, 1762, 1765, 1766, 1877, 1878, 1984, 2032, 2043, 2054, 2073, 2077, 2092, 2093, 2097, 2107, 2162, 2168, 2214, 2305, 2392, 2409, 2415, 2439, 2459, 2589, 2590, 2633, 2652, 2749, 2834, 2961, 2969, 3005, 3017, 3070, 3071, 3186, 3220, 3289, 3360, 3493, 3494, 3495, 3508, 3588, 3604, 3611, 3668, 3678, 3679, 3690, 3692, 3732, 3738, 3838, 3840, 3927, 3945, 4005, 4006, 4013, 4018, 4019, 4053, 4059, 4124, 4331, 4355, 4357, 4458, 4464, 4500, 4677, 4678, 4684, 4768, 360, 493, 779, 1284, 1306, 1337, 1438, 1709, 2027, 2445, 3407, 3720, 4217, 4288, 95, 96, 296, 316, 443, 479, 555, 562, 637, 671, 709, 740, 743, 796, 857, 897, 899, 915, 953, 1130, 1132, 1294, 1318, 1337, 1348, 1359, 1378, 1427, 1429, 1440, 1537, 1563, 1574, 1709, 1763, 1773, 1819, 1854, 1874, 1936, 1938, 2018, 2043, 2098, 2099, 2256, 2259, 2294, 2514, 2520, 2552, 2589, 2598, 2680, 2706, 2799, 2880, 2897, 2917, 2969, 3020, 3028, 3310, 3319, 3375, 3399, 3504, 3539, 3541, 3679, 3690, 3692, 3718, 3719, 3720, 3738, 3806, 3838, 3840, 3928, 4013, 4017, 4038, 4053, 4057, 4058, 4079, 4097, 4099};

  Bool_t isGoodTwr = true;
  for (UInt_t iTwr = 0; iTwr < NHotTwr; iTwr++) {
    if (twrID == hotTwrList[iTwr]) {
      isGoodTwr = false;
      break;
    }
  }
  return isGoodTwr;

}  // end 'IsGoodTowerID(UInt_t)'



Bool_t StMuDstJetTreeMaker::IsGoodEvent(const Double_t rVtx, const Double_t zVtx) {

  const Bool_t isInRCut  = (abs(rVtx) < _rVtxMax);
  const Bool_t isInZCut  = (abs(zVtx) < _zVtxMax);
  const Bool_t isGoodEvt = (isInRCut && isInZCut);
  return isGoodEvt;

}  // end 'IsGoodEvent(Double_t, Double_t)'



Bool_t StMuDstJetTreeMaker::IsGoodPionTrigger(const Int_t adc, const Double_t eEta, const Double_t ePhi, const Double_t pProj, const Double_t etaTrg, const Double_t eTtrg, const Double_t tsp) {

  const Bool_t isInAdcCut  = (adc <= _adcMax);
  const Bool_t isInStrpCut = ((eEta >= _eEtaMin) && (ePhi >= _ePhiMin));
  const Bool_t isInProjCut = (pProj < _pProjMax);
  const Bool_t isInEtaCut  = (abs(etaTrg) < _etaTrgMax);
  const Bool_t isInEtCut   = ((eTtrg > _eTtrgMin) && (eTtrg < _eTtrgMax));
  const Bool_t isInPi0Cut  = ((tsp > _tspPi0Min) && (tsp < _tspPi0Max));
  const Bool_t isInGamCut  = ((tsp > _tspGamMin) && (tsp < _tspGamMax));
  const Bool_t isInTspCut  = (isInPi0Cut || isInGamCut);
  const Bool_t isGoodTrg   = (isInAdcCut && isInStrpCut && isInProjCut && isInEtaCut && isInEtCut && isInTspCut);
  return isGoodTrg;

}  // end 'IsGoodTrigger(Int_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



Bool_t StMuDstJetTreeMaker::IsGoodHadronTrigger(const Double_t etaTrg, const Double_t pTtrg, const Double_t nSigPi, const Double_t nSigK, const Double_t nSigP, const Double_t nSigE) {

  const Bool_t isInEtaCut = (abs(etaTrg) < _etaTrgMax);
  const Bool_t isInPtCut  = ((pTtrg > _eTtrgMin) && (pTtrg < _eTtrgMax));
  const Bool_t isHadron   = IsHadron(nSigPi, nSigK, nSigP, nSigE);
  const Bool_t isGoodTrg  = (isInEtaCut && isInPtCut && isHadron);
  return isGoodTrg;

}  // end 'IsGoodHadronTrigger(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



Bool_t StMuDstJetTreeMaker::IsGoodTrack(const UInt_t nFit, const Double_t rFit, const Double_t dca, const Double_t etaTrk, const Double_t pTtrk) {

  const Bool_t isInFitCut   = (nFit >= _nFitMin);
  const Bool_t isInRatioCut = (rFit >= _rFitMin);
  const Bool_t isInDcaCut   = (dca < _dcaMax);
  const Bool_t isInEtaCut   = (abs(etaTrk) < _etaTrkMax);
  const Bool_t isInPtCut    = ((pTtrk > _pTtrkMin) && (pTtrk < _pTtrkMax));
  const Bool_t isGoodTrk    = (isInFitCut && isInRatioCut && isInDcaCut && isInEtaCut && isInPtCut);
  return isGoodTrk;

}  // end 'IsGoodTrack(UInt_t, Double_t, Double_t, Double_t, Double_t)'



Bool_t StMuDstJetTreeMaker::IsGoodTower(const Double_t etaTwr, const Double_t eTwr, const Double_t eCorr) {

  const Bool_t isInEtaCut  = (abs(etaTwr) < _etaTwrMax);
  const Bool_t isInEneCut  = (eTwr > _eTwrMin);
  const Bool_t isInCorrCut = ((eCorr > _eCorrMin) && (eCorr < _eCorrMax));
  const Bool_t isGoodTwr   = (isInEtaCut && isInEneCut && isInCorrCut);
  return isGoodTwr;

}  // end 'IsGoodTwr(Double_t, Double_t, Double_t)'



Bool_t StMuDstJetTreeMaker::IsPi0(const Double_t tsp) {

  const Bool_t isPi0 = ((tsp > _tspPi0Min) && (tsp < _tspPi0Max));
  return isPi0;

}  // end 'IsPi0(Double_t)'



Bool_t StMuDstJetTreeMaker::IsGamma(const Double_t tsp) {

  const Bool_t isGam = ((tsp > _tspGamMin) && (tsp < _tspGamMax));
  return isGam;

}  // end 'IsGamma(Double_t)'



Bool_t StMuDstJetTreeMaker::IsHadron(const Double_t nSigPi, const Double_t nSigK, const Double_t nSigP, const Double_t nSigE) {

  const Bool_t isPion     = (abs(nSigPi) < NSigCut);
  const Bool_t isKaon     = (abs(nSigK)  < NSigCut);
  const Bool_t isProton   = (abs(nSigP)  < NSigCut);
  const Bool_t isElectron = (abs(nSigE)  < NSigCut);
  const Bool_t isHadron   = ((isPion || isKaon || isProton) && !isElectron);
  return isHadron;

}  // end 'IsHadron(Double_t, Double_t, Double_t, Double_t)'



Long64_t StMuDstJetTreeMaker::LoadTree(const Long64_t entry) {

  // set the environment to read one entry
  Int_t    fTree;
  Long64_t bytes;
  if (!_tFemto)
    bytes = -5;
  else {
    bytes = _tFemto -> LoadTree(entry);
    fTree = _tFemto -> GetTreeNumber();
    if ((bytes >= 0) && (fTree != _fCurrent))
      _fCurrent = _tFemto -> GetTreeNumber();
  }
  return bytes;

}  // end 'LoadTree(Long64_t)'



Long64_t StMuDstJetTreeMaker::GetEntry(const Long64_t entry) {

  Int_t bytes;
  if (!_tFemto)
    bytes = 0;
  else
    bytes = _tFemto -> GetEntry(entry);
  return bytes;

}  // end 'GetEntry(Long64_t)'



Double_t StMuDstJetTreeMaker::GetHadronicCorrection(const Double_t eTwr, const vector<Double_t> pMatchedTrks) {

  Double_t eSum  = 0.;
  Double_t eCorr = 0.;

  // sum momentum of matched tracks
  const UInt_t nTrks = (UInt_t) pMatchedTrks.size();
  for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {
    eSum += pMatchedTrks.at(iTrk);
  }

  // subtract summed momentum
  const Double_t percent = 0.5;
  const Double_t corr    = eTwr - (percent * eSum);
  if (corr < 0.)
    eCorr = 0.;
  else
    eCorr = corr;
  return eCorr;

}  // end 'GetHadronicCorrection(Double_t, vector<Double_t>)'



TLorentzVector StMuDstJetTreeMaker::GetTowerMomentumVector(const Double_t rBEMC, const Double_t etaTwr, const Double_t phiTwr, const Double_t eTwr, const TVector3& vtx) {

  const Double_t mMin = 0.;

  // calculate 3-vector
  TVector3 rTwr;
  TVector3 rVtx;
  TVector3 pTwr;
  rTwr.SetPtEtaPhi(rBEMC, etaTwr, phiTwr);
  rVtx = vtx;
  pTwr = rTwr - rVtx;

  // calculate magnitude
  Double_t pMag;
  if (eTwr > mMin)
    pMag = sqrt((eTwr * eTwr) - (mMin * mMin));
  else
    pMag = eTwr;
  pTwr.SetMag(pMag);

  // return 4-vector
  TLorentzVector qTwr(pTwr.x(), pTwr.y(), pTwr.z(), eTwr);
  return qTwr;

}  // end 'GetTowerMomentumVector(Double_t, Double_t, Double_t, Double_t, TVector3&)'

// End ------------------------------------------------------------------------

