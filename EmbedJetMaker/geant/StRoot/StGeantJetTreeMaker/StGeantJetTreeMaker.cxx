// 'StGeantJetTreeMaker.cxx'
// Derek Anderson, Nihar Sahoo
// 10.03.2017
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets


#define StGeantJetTreeMaker_cxx

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StGeantJetTreeMaker.h"
#include "StGeantJetTreeMaker.utilities.h"

using namespace std;
using namespace fastjet;

ClassImp(StGeantJetTreeMaker);



void StGeantJetTreeMaker::SetInputAndOutputFiles(const TString& sInput, const TString& sOutput, const Double_t pTparton) {

  _sInput     = sInput;
  _sOutput    = sOutput;
  _PartonicPt = pTparton;
  PrintInfo(1);

}  // end 'SetInputAndOutput(TString&, TString&)'



void StGeantJetTreeMaker::SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax) {

  _rVtxMax = rVtxMax;
  _zVtxMax = zVtxMax;
  PrintInfo(2);

}  // end 'StEventParameters(Double_t, Double_t)'



void StGeantJetTreeMaker::SetTriggerParameters(const Double_t etaTrgMax, const Double_t eTtrgMin, const Double_t eTtrgMax) {

  _etaTrgMax = etaTrgMax;
  _eTtrgMin  = eTtrgMin;
  _eTtrgMax  = eTtrgMax;
  PrintInfo(3);

}  // end 'SetTriggerParameters(Double_t, Double_t, Double_t)'



void StGeantJetTreeMaker::SetTrackParameters(const Double_t etaTrkMax, const Double_t pTtrkMin, const Double_t pTtrkMax) {

  _etaTrkMax = etaTrkMax;
  _pTtrkMin  = pTtrkMin;
  _pTtrkMax  = pTtrkMax;
  PrintInfo(4);

}  // end 'SetTrackParameters(UInt_t, Double_t, Double_t)'



void StGeantJetTreeMaker::SetJetParameters(const UInt_t type, const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t aGhost, const Double_t pTjetMin, const Double_t etaGhostMax, const Double_t etaJetMax, const Double_t etaBkgdMax) {

  _jetType     = type;
  _nRepeat     = nRepeat;
  _nRemove     = nRemove;
  _rJet        = rJet;
  _aGhost      = aGhost;
  _pTjetMin    = pTjetMin;
  _etaGhostMax = etaGhostMax;
  _etaJetMax   = etaJetMax;
  _etaBkgdMax  = etaBkgdMax;
  PrintInfo(6);

}  // 'SetJetParameters(UInt_t, UInt_t, UInt_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



void StGeantJetTreeMaker::Init() {

  _fOutput = new TFile(_sOutput.Data(), "recreate");
  if (_tFemto == 0) {
    TFile *fInput = (TFile*) gROOT -> GetListOfFiles() -> FindObject(_sInput.Data());
    if (!fInput || !(fInput -> IsOpen())) {
      fInput = new TFile(_sInput.Data());
    }
    fInput -> GetObject("GfmtoDst_gnt", _tFemto);
  }
  InitializeInputTree(_tFemto);
  InitializeOutputTree(_tJet);
  InitializeHistograms();
  PrintInfo(10);

}  // end 'Init()'



void StGeantJetTreeMaker::Make() {

  if (_tFemto == 0) return;


  const Long64_t nEvts = _tFemto -> GetEntriesFast();
  PrintInfo(11, nEvts);

  // for jet-finding
  vector<Double_t>  matches;
  vector<PseudoJet> particles;
  vector<PseudoJet> jetsCS;
  vector<PseudoJet> jets;
  matches.clear();
  particles.clear();
  jetsCS.clear();
  jets.clear();

  // event loop
  Long64_t nGam   = 0;
  Long64_t nPi0   = 0;
  Long64_t nHad   = 0;
  Long64_t nBytes = 0;
  for (Long64_t iEvt = 0; iEvt < nEvts; iEvt++) {

    const Long64_t entry = LoadTree(iEvt);
    if (entry < 0) break;
    nBytes += _tFemto -> GetEntry(iEvt);
    particles.clear();
    PrintInfo(12, nEvts, iEvt);


    // event info
    const UInt_t   runID = _runNumber;
    const Long64_t nTrks = _nPrimaryTracks;
    const Double_t xVtx  = _xVertex;
    const Double_t yVtx  = _yVertex;
    const Double_t zVtx  = _zVertex;
    const Double_t rVtx  = sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // event cuts
    const Bool_t isGoodRun = IsGoodRunID(runID);
    const Bool_t isGoodEvt = IsGoodEvent(rVtx, zVtx);
    if (!isGoodRun || !isGoodEvt) continue;


    // trigger loop
    UInt_t   idTrg     = 0;
    Bool_t   isGoodTrg = false;
    Double_t hTrg      = 0.;
    Double_t fTrg      = 0.;
    Double_t eTrg      = 0.;
    Double_t eTtrg     = 0.;
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Double_t idTrk = _PrimaryTrackArray_geantId[iTrk];
      const Double_t hTrk  = _PrimaryTrackArray_eta[iTrk];
      const Double_t pXtrk = _PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk = _PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk = _PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk = _PrimaryTrackArray_pT[iTrk];
      const Double_t eTrk  = sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (MassPi0 * MassPi0));
      const Double_t eTtrk = eTrk * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hTrk)));

      // calculate phi
      Double_t fTrk = TMath::ATan(pYtrk / pXtrk);
      if ((pXtrk < 0.) && (pYtrk > 0.)) fTrk += TMath::Pi();
      if ((pXtrk < 0.) && (pYtrk < 0.)) fTrk += TMath::Pi();

      // trigger cuts
      isGoodTrg = IsGoodTrigger(hTrk, eTtrk, idTrk);
      if (isGoodTrg) {
        idTrg = idTrk;
        hTrg  = hTrk;
        fTrg  = fTrk;
        eTrg  = eTrk;
        eTtrg = eTtrk;
        break;
      }

    }  // end trigger loop

    // trigger cuts
    const Bool_t isPi0 = IsPi0(idTrg);
    const Bool_t isGam = IsGamma(idTrg);
    const Bool_t isHad = IsHadron(idTrg);
    if (!isGoodTrg) continue;


    // count triggers
    if (isPi0) {
      _hEvtQA[1][0] -> Fill(idTrg);
      _hEvtQA[2][0] -> Fill(nTrks);
      _hEvtQA[3][0] -> Fill(eTrg);
      nPi0++;
    }
    if (isGam) {
      _hEvtQA[1][1] -> Fill(idTrg);
      _hEvtQA[2][1] -> Fill(nTrks);
      _hEvtQA[3][1] -> Fill(eTrg);
      nGam++;
    }
    if (isHad) {
      _hEvtQA[1][2] -> Fill(idTrg);
      _hEvtQA[2][2] -> Fill(nTrks);
      _hEvtQA[3][2] -> Fill(eTrg);
      nHad++;
    }


    // track loop 
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Double_t qTrk  = _PrimaryTrackArray_charge[iTrk];
      const Double_t hTrk  = _PrimaryTrackArray_eta[iTrk];
      const Double_t pXtrk = _PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk = _PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk = _PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk = _PrimaryTrackArray_pT[iTrk];
      const Double_t eTrk  = sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (MassPi0 * MassPi0));

      // calculate phi
      Double_t fTrk = TMath::ATan(pYtrk / pXtrk);
      if ((pXtrk < 0.) && (pYtrk > 0.)) fTrk += TMath::Pi();
      if ((pXtrk < 0.) && (pYtrk < 0.)) fTrk += TMath::Pi();

      Double_t dFtrk = fTrk - fTrg;
      if (dFtrk < (-1. * TMath::PiOver2()))
        dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))
        dFtrk -= TMath::TwoPi();

      // track cuts
      const Bool_t isCharged = (qTrk != 0.);
      const Bool_t isGoodTrk = IsGoodTrack(hTrk, pTtrk);
      if (!isGoodTrk) continue;


      if (isCharged) {
        if (isPi0) {
          _hTrkQA[0][0] -> Fill(pTtrk);
          _hTrkQA[1][0] -> Fill(dFtrk);
          _hTrkQA[2][0] -> Fill(hTrk);
          _hTrkQA[3][0] -> Fill(eTrk);
        }
        if (isGam) {
          _hTrkQA[0][1] -> Fill(pTtrk);
          _hTrkQA[1][1] -> Fill(dFtrk);
          _hTrkQA[2][1] -> Fill(hTrk);
          _hTrkQA[3][1] -> Fill(eTrk);
        }
        if (isHad) {
          _hTrkQA[0][2] -> Fill(pTtrk);
          _hTrkQA[1][2] -> Fill(dFtrk);
          _hTrkQA[2][2] -> Fill(hTrk);
          _hTrkQA[3][2] -> Fill(eTrk);
        }
        particles.push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));
      }
      else {
        if (isPi0) {
          _hTwrQA[0][0] -> Fill(pTtrk);
          _hTwrQA[1][0] -> Fill(dFtrk);
          _hTwrQA[2][0] -> Fill(hTrk);
          _hTwrQA[3][0] -> Fill(eTrk);
        }
        if (isGam) {
          _hTwrQA[0][1] -> Fill(pTtrk);
          _hTwrQA[1][1] -> Fill(dFtrk);
          _hTwrQA[2][1] -> Fill(hTrk);
          _hTwrQA[3][1] -> Fill(eTrk);
        }
        if (isHad) {
          _hTwrQA[0][2] -> Fill(pTtrk);
          _hTwrQA[1][2] -> Fill(dFtrk);
          _hTwrQA[2][2] -> Fill(hTrk);
          _hTwrQA[3][2] -> Fill(eTrk);
        }
        if (_jetType == 1)
          particles.push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));
      }

    }  // end track loop


    // define jets and jet area
    GhostedAreaSpec areaSpec(_etaGhostMax, _nRepeat, _aGhost);
    AreaDefinition  areaDef(active_area_explicit_ghosts, areaSpec);
    JetDefinition   jetDef(antikt_algorithm, _rJet);

    // cluster jets
    ClusterSequenceArea clusterSeq(particles, jetDef, areaDef);
    jetsCS = sorted_by_pt(clusterSeq.inclusive_jets(_pTjetMin));

    // fiducial cut
    Selector fidCut = SelectorAbsEtaMax(_etaJetMax);
    jets = fidCut(jetsCS);

    // define bkgd jets and jet area
    GhostedAreaSpec areaSpecBkgd(_etaGhostMax, _nRepeat, _aGhost);
    AreaDefinition  areaDefBkgd(active_area_explicit_ghosts, areaSpecBkgd);
    JetDefinition   jetDefBkgd(kt_algorithm, _rJet);

    // initialize bkgd estimators
    Selector bkgdCut = SelectorAbsEtaMax(_etaBkgdMax) * (!SelectorNHardest(_nRemove));
    JetMedianBackgroundEstimator bkgd(bkgdCut, jetDefBkgd, areaDefBkgd);
    Subtractor sub(&bkgd);

#if FASTJET_VERSION_NUMBER >= 30100
  sub.set_use_rho_m(true);
  sub.set_safe_mass(true);
#endif

    // estimate bkgd
    bkgd.set_particles(particles);
    const Double_t rhoJet = bkgd.rho();
    const Double_t sigJet = bkgd.sigma();


    // clear jet info
    _JetIndex.clear();
    _JetNCons.clear();
    _JetPt.clear();
    _JetPtCorr.clear();
    _JetEta.clear();
    _JetPhi.clear();
    _JetE.clear();
    _JetArea.clear();
    _JetConsPt.clear();
    _JetConsEta.clear();
    _JetConsPhi.clear();
    _JetConsE.clear();

    // jet event info
    _EventIndex = _eventNumber;
    _NJets      = (Int_t) jets.size();
    _RunId      = _runNumber;
    _Refmult    = _refMult;
    _TSP        = (Double_t) idTrg;
    _TrgEta     = hTrg;
    _TrgPhi     = fTrg;
    _TrgEt      = eTtrg;
    _Rho        = rhoJet;
    _Sigma      = sigJet;
    _Vz         = _zVertex;

    // jet loop
    const UInt_t nJets = (UInt_t) jets.size();
    for (UInt_t iJet = 0; iJet < nJets; iJet++) {

      // jet info
      const Int_t    nCst   = (Int_t) jets[iJet].constituents().size();
      const Double_t hJet   = jets[iJet].pseudorapidity();
      const Double_t fJet   = jets[iJet].phi_std();
      const Double_t aJet   = jets[iJet].area();
      const Double_t eJet   = jets[iJet].e();
      const Double_t pTjet  = jets[iJet].perp();
      const Double_t pTcorr = pTjet - (aJet * rhoJet);

      Double_t dFjet = fJet - fTrg;
      if (dFjet < (-1. * TMath::PiOver2()))
        dFjet += TMath::TwoPi();
      if (dFjet > (3. * TMath::PiOver2()))
        dFjet -= TMath::TwoPi();


      if (isPi0) {
        _hJetQA[0][0] -> Fill(pTjet);
        _hJetQA[1][0] -> Fill(dFjet);
        _hJetQA[2][0] -> Fill(hJet);
        _hJetQA[3][0] -> Fill(eJet);
      }
      if (isGam) {
        _hJetQA[0][1] -> Fill(pTjet);
        _hJetQA[1][1] -> Fill(dFjet);
        _hJetQA[2][1] -> Fill(hJet);
        _hJetQA[3][1] -> Fill(eJet);
      }
      if (isHad) {
        _hJetQA[0][2] -> Fill(pTjet);
        _hJetQA[1][2] -> Fill(dFjet);
        _hJetQA[2][2] -> Fill(hJet);
        _hJetQA[3][2] -> Fill(eJet);
      }
      _JetIndex.push_back(iJet);
      _JetNCons.push_back(nCst);
      _JetPt.push_back(pTjet);
      _JetPtCorr.push_back(pTcorr);
      _JetEta.push_back(hJet);
      _JetPhi.push_back(fJet);
      _JetE.push_back(eJet);
      _JetArea.push_back(aJet);

    }  // end jet loop

    _tJet -> Fill();

  }  // end event loop

  // record no. of triggers
  _hEvtQA[0][0] -> Fill(nPi0);
  _hEvtQA[0][1] -> Fill(nGam);
  _hEvtQA[0][2] -> Fill(nHad);
  PrintInfo(13);


  // normalize histograms
  for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
    Double_t nTrg;
    switch (iTrg) {
      case 0:
        nTrg = (Double_t) nPi0;
        break;
      case 1:
        nTrg = (Double_t) nGam;
        break;
      case 2:
        nTrg = (Double_t) nHad;
        break;
    }
    for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
      const Double_t trkBin   = _hTrkQA[iHist][iTrg] -> GetBinWidth(17);
      const Double_t twrBin   = _hTwrQA[iHist][iTrg] -> GetBinWidth(17);
      const Double_t jetBin   = _hJetQA[iHist][iTrg] -> GetBinWidth(17);
      const Double_t trkScale = 1. / (nTrg * trkBin);
      const Double_t twrScale = 1. / (nTrg * twrBin);
      const Double_t jetScale = 1. / (nTrg * jetBin);
      _hTrkQA[iHist][iTrg] -> Scale(trkScale);
      _hTwrQA[iHist][iTrg] -> Scale(twrScale);
      _hJetQA[iHist][iTrg] -> Scale(jetScale);
    }
  }
  PrintInfo(14);


}  // end 'Make()'



void StGeantJetTreeMaker::Finish() {

  _fOutput -> cd();
  TDirectory *dHists = _fOutput -> mkdir("QA");
  TDirectory *dPi0QA = dHists   -> mkdir("Pi0");
  TDirectory *dGamQA = dHists   -> mkdir("Gamma");
  TDirectory *dHadQA = dHists   -> mkdir("Hadron");
  PrintInfo(15);

  for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
    for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
      switch (iTrg) {
        case 0:
          dPi0QA -> cd();
          break;
        case 1:
          dGamQA -> cd();
          break;
        case 2:
          dHadQA -> cd();
          break;
      }
      _hEvtQA[iHist][iTrg] -> Write();
      _hTrkQA[iHist][iTrg] -> Write();
      _hTwrQA[iHist][iTrg] -> Write();
      _hJetQA[iHist][iTrg] -> Write();
    }
  }
  _fOutput -> cd();
  _tJet    -> Write();
  _fOutput -> Close();
  _tFemto  -> GetCurrentFile() -> cd();
  _tFemto  -> GetCurrentFile() -> Close();
  PrintInfo(16);

}  // end 'Finish()'

// End ------------------------------------------------------------------------
