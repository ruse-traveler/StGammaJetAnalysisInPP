// 'StMcJetTreeMaker.cxx'
// Derek Anderson, Nihar Sahoo
// 04.30.2018
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets


#define StMcJetTreeMaker_cxx

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StMcJetTreeMaker.h"
#include "StMcJetTreeMaker.utilities.h"

using namespace std;
using namespace fastjet;

ClassImp(StMcJetTreeMaker);



void StMcJetTreeMaker::SetInputAndOutputFiles(const TString& sInput, const TString& sOutput, const Double_t pTparton) {

  _sInput     = sInput;
  _sOutput    = sOutput;
  _PartonicPt = pTparton;
  PrintInfo(1);

}  // end 'SetInputAndOutput(TString&, TString&)'



void StMcJetTreeMaker::SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax) {

  _rVtxMax = rVtxMax;
  _zVtxMax = zVtxMax;
  PrintInfo(2);

}  // end 'StEventParameters(Double_t, Double_t)'



void StMcJetTreeMaker::SetTriggerParameters(const Double_t etaTrgMax, const Double_t eTtrgMin, const Double_t eTtrgMax) {

  _etaTrgMax = etaTrgMax;
  _eTtrgMin  = eTtrgMin;
  _eTtrgMax  = eTtrgMax;
  PrintInfo(3);

}  // end 'SetTriggerParameters(Double_t, Double_t, Double_t)'



void StMcJetTreeMaker::SetTrackParameters(const Double_t etaTrkMax, const Double_t pTtrkMin, const Double_t pTtrkMax) {

  _etaTrkMax = etaTrkMax;
  _pTtrkMin  = pTtrkMin;
  _pTtrkMax  = pTtrkMax;
  PrintInfo(4);

}  // end 'SetTrackParameters(UInt_t, Double_t, Double_t)'



void StMcJetTreeMaker::SetJetParameters(const UInt_t type, const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t aGhost, const Double_t pTjetMin, const Double_t etaGhostMax, const Double_t etaJetMax, const Double_t etaBkgdMax) {

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



void StMcJetTreeMaker::AdjustTrackEfficiency(const Bool_t effAdjust, const Float_t adjustment) {

  _adjustTrackEff = effAdjust;
  _effAdjust      = adjustment;
  _random         = new TRandom();
  PrintInfo(17);

}  // 'AdjustTrackEfficiency(Bool_t, Float_t)'



void StMcJetTreeMaker::TurnOffWeakDecays(const Bool_t turnOffDecays) {

  _turnOffDecays = turnOffDecays;
  PrintInfo(18);

}  // 'TurnOffWeakDecays(Bool_t)'



void StMcJetTreeMaker::Init() {

  _fOutput = new TFile(_sOutput.Data(), "recreate");
  if (_tMc == 0) {
    TFile *fInput = (TFile*) gROOT -> GetListOfFiles() -> FindObject(_sInput.Data());
    if (!fInput || !(fInput -> IsOpen())) {
      fInput = new TFile(_sInput.Data());
    }
    fInput -> GetObject("McTracks", _tMc);
  }
  InitializeInputTree(_tMc);
  InitializeOutputTree(_tJet);
  InitializeHistograms();
  PrintInfo(10);

}  // end 'Init()'



void StMcJetTreeMaker::Make(const UInt_t evtFlag, const UInt_t trgFlag) {

  if (_tMc == 0) return;

  // set event and trigger mode
  Bool_t applyEventCuts(false);
  Bool_t requireTrigger(false);
  switch (evtFlag) {
    case 0:
      applyEventCuts = false;
      break;
    default:
      applyEventCuts = true;
      break;
  }
  switch (trgFlag) {
    case 0:
      requireTrigger = false;
      break;
    default:
      requireTrigger = true;
      break;
  }


  const Long64_t nEvts = _tMc -> GetEntriesFast();
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
    nBytes += _tMc -> GetEntry(iEvt);
    particles.clear();
    PrintInfo(12, nEvts, iEvt);


    // event info
    const UInt_t   runID = _mcRunId;
    const Long64_t nTrks = _mcNumTrks;
    const Double_t xVtx  = _mcVx;
    const Double_t yVtx  = _mcVy;
    const Double_t zVtx  = _mcVz;
    const Double_t rVtx  = sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // event cuts
    const Bool_t isGoodRun  = IsGoodRunID(runID);
    const Bool_t isGoodEvt  = IsGoodEvent(rVtx, zVtx);
    const Bool_t isInEvtCut = (isGoodRun && isGoodEvt);
    if (applyEventCuts && !isInEvtCut) continue;


    // trigger loop
    UInt_t   iTrg     = nTrks + 1;
    UInt_t   idTrg    = 0;
    Bool_t   foundTrg = false;
    Double_t hTrg     = 0.;
    Double_t fTrg     = 0.;
    Double_t eTtrg    = 0.;
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Int_t    idVtx = _mcIdVx     -> at(iTrk);
      const Double_t idTrk = _mcIdGeant  -> at(iTrk);
      const Double_t hTrk  = _mcEta      -> at(iTrk);
      const Double_t pXtrk = _mcPx       -> at(iTrk);
      const Double_t pYtrk = _mcPy       -> at(iTrk);
      const Double_t pZtrk = _mcPz       -> at(iTrk);
      const Double_t pTtrk = _mcPt       -> at(iTrk);
      const Double_t eTrk  = sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (MassPi0 * MassPi0));
      const Double_t eTtrk = eTrk * TMath::Sin(2. * TMath::ATan(TMath::Exp(-1. * hTrk)));

      // calculate phi
      Double_t fTrk = TMath::ATan(pYtrk / pXtrk);
      if ((pXtrk < 0.) && (pYtrk > 0.)) fTrk += TMath::Pi();
      if ((pXtrk < 0.) && (pYtrk < 0.)) fTrk += TMath::Pi();

      // trigger cuts
      const Bool_t isGoodTrg   = IsGoodTrigger(hTrk, eTtrk, idTrk);
      const Bool_t isUndecayed = IsUndecayed(idVtx);
      if (isGoodTrg && isUndecayed) {
        iTrg     = iTrk;
        idTrg    = idTrk;
        hTrg     = hTrk;
        fTrg     = fTrk;
        eTtrg    = eTtrk;
        foundTrg = true;
        break;
      }

    }  // end trigger loop


    // trigger cuts
    const Bool_t isPi0 = IsPi0(idTrg);
    const Bool_t isGam = IsGamma(idTrg);
    const Bool_t isHad = IsHadron(idTrg);
    if (requireTrigger && !foundTrg) continue;

    // count triggers
    if (isPi0) {
      _hEvtQA[1][0] -> Fill(idTrg);
      _hEvtQA[2][0] -> Fill(nTrks);
      _hEvtQA[3][0] -> Fill(eTtrg);
      nPi0++;
    }
    if (isGam) {
      _hEvtQA[1][1] -> Fill(idTrg);
      _hEvtQA[2][1] -> Fill(nTrks);
      _hEvtQA[3][1] -> Fill(eTtrg);
      nGam++;
    }
    if (isHad) {
      _hEvtQA[1][2] -> Fill(idTrg);
      _hEvtQA[2][2] -> Fill(nTrks);
      _hEvtQA[3][2] -> Fill(eTtrg);
      nHad++;
    }


    // track loop 
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const Int_t    idVtx = _mcIdVx    -> at(iTrk);
      const Double_t qTrk  = _mcCharge  -> at(iTrk);
      const Double_t hTrk  = _mcEta     -> at(iTrk);
      const Double_t pXtrk = _mcPx      -> at(iTrk);
      const Double_t pYtrk = _mcPy      -> at(iTrk);
      const Double_t pZtrk = _mcPz      -> at(iTrk);
      const Double_t pTtrk = _mcPt      -> at(iTrk);
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
      const Bool_t isCharged    = (qTrk != 0.);
      const Bool_t isUndecayed  = IsUndecayed(idVtx);
      const Bool_t isFinalState = IsFinalState(pTtrk);
      const Bool_t isGoodTrk    = IsGoodTrack(hTrk, pTtrk);
      const Bool_t isGoodIS     = (isGoodTrk && isUndecayed);
      const Bool_t isGoodFS     = (isGoodTrk && isFinalState);

      // select correct state
      if (_turnOffDecays && !isGoodIS)  continue;
      if (!_turnOffDecays && !isGoodFS) continue;

      // remove trigger
      const Bool_t isTrigger = (iTrk == iTrg);
      if (isTrigger) continue;

      // adjust tracking efficiency (if need be)
      if (_adjustTrackEff && isCharged) {
        const Float_t rando = _random -> Uniform(0., 1.);
        const Bool_t  pass  = (rando > _effAdjust);
        if (!pass) continue;
      }

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
        if (_jetType == 1) {
          particles.push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));
        }
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
    _EventIndex = _mcEventId;
    _NJets      = (Int_t) jets.size();
    _RunId      = _mcRunId;
    _Refmult    = _mcNumTrks;
    _TSP        = (Double_t) idTrg;
    _TrgEta     = hTrg;
    _TrgPhi     = fTrg;
    _TrgEt      = eTtrg;
    _Rho        = rhoJet;
    _Sigma      = sigJet;
    _Vz         = _mcVz;

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

      // TEST [03.30.2020]
      UInt_t nCstReal = 0;
      for (UInt_t iCst = 0; iCst < nCst; iCst++) {
        const Double_t pTcst = jets[iJet].constituents()[iCst].perp();
        if (pTcst > 0.2) nCstReal++;
      }

      // TEST [04.02.2020]
      for (UInt_t iCst = 0; iCst < nCst; iCst++) {
        const Double_t pTcst = jets[iJet].constituents()[iCst].perp();
        const Double_t zTcst = pTcst / pTjet;
        if ((pTcst > 0.2) && (pTjet > 0.2)) {

          // zTcst vs. pTjet
          if ((pTjet >= 0.2) && (pTjet < 1.)) {
            if (isPi0) _hZtCstVsPtJetP[0] -> Fill(zTcst);
            if (isGam) _hZtCstVsPtJetG[0] -> Fill(zTcst);
            if (isHad) _hZtCstVsPtJetH[0] -> Fill(zTcst);
          }
          if ((pTjet >= 1.) && (pTjet < 2.)) {
            if (isPi0) _hZtCstVsPtJetP[1] -> Fill(zTcst);
            if (isGam) _hZtCstVsPtJetG[1] -> Fill(zTcst);
            if (isHad) _hZtCstVsPtJetH[1] -> Fill(zTcst);
          }
          if ((pTjet >= 2.) && (pTjet < 5.)) {
            if (isPi0) _hZtCstVsPtJetP[2] -> Fill(zTcst);
            if (isGam) _hZtCstVsPtJetG[2] -> Fill(zTcst);
            if (isHad) _hZtCstVsPtJetH[2] -> Fill(zTcst);
          }
          if ((pTjet >= 5.) && (pTjet < 10.)) {
            if (isPi0) _hZtCstVsPtJetP[3] -> Fill(zTcst);
            if (isGam) _hZtCstVsPtJetG[3] -> Fill(zTcst);
            if (isHad) _hZtCstVsPtJetH[3] -> Fill(zTcst);
          }
          if (pTjet >= 10.) {
            if (isPi0) _hZtCstVsPtJetP[4] -> Fill(zTcst);
            if (isGam) _hZtCstVsPtJetG[4] -> Fill(zTcst);
            if (isHad) _hZtCstVsPtJetH[4] -> Fill(zTcst);
          }

          // zTcst vs. nCstJet
          if ((nCstReal >= 1) && (nCstReal < 2)) {
            if (isPi0) _hZtCstVsNcstP[0] -> Fill(zTcst);
            if (isGam) _hZtCstVsNcstG[0] -> Fill(zTcst);
            if (isHad) _hZtCstVsNcstH[0] -> Fill(zTcst);
          }
          if ((nCstReal >= 2) && (nCstReal < 5)) {
            if (isPi0) _hZtCstVsNcstP[1] -> Fill(zTcst);
            if (isGam) _hZtCstVsNcstG[1] -> Fill(zTcst);
            if (isHad) _hZtCstVsNcstH[1] -> Fill(zTcst);
          }
          if ((nCstReal >= 5) && (nCstReal < 10)) {
            if (isPi0) _hZtCstVsNcstP[2] -> Fill(zTcst);
            if (isGam) _hZtCstVsNcstG[2] -> Fill(zTcst);
            if (isHad) _hZtCstVsNcstH[2] -> Fill(zTcst);
          }
          if (nCstReal >= 10) {
            if (isPi0) _hZtCstVsNcstP[3] -> Fill(zTcst);
            if (isGam) _hZtCstVsNcstG[3] -> Fill(zTcst);
            if (isHad) _hZtCstVsNcstH[3] -> Fill(zTcst);
          }
        }
      }
      if (isPi0 && (pTjet > 0.2)) _hNcstVsPtJetP -> Fill(pTjet, nCstReal);
      if (isGam && (pTjet > 0.2)) _hNcstVsPtJetG -> Fill(pTjet, nCstReal);
      if (isHad && (pTjet > 0.2)) _hNcstVsPtJetH -> Fill(pTjet, nCstReal);


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
      //_JetNCons.push_back(nCst);
      _JetNCons.push_back(nCstReal);  // TEST [03.30.2020]
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

}  // end 'Make(UInt_t)'



void StMcJetTreeMaker::Finish() {

  _fOutput -> cd();
  TDirectory *dHists = _fOutput -> mkdir("QA");
  TDirectory *dPi0QA = dHists   -> mkdir("Pi0");
  TDirectory *dGamQA = dHists   -> mkdir("Gamma");
  TDirectory *dHadQA = dHists   -> mkdir("Hadron");
  PrintInfo(15);

  // TEST [04.02.2020]
  _fOutput -> cd();
  for (UInt_t iPtTest = 0; iPtTest < 5; iPtTest++) {
    _hZtCstVsPtJetP[iPtTest] -> Write();
    _hZtCstVsPtJetG[iPtTest] -> Write();
    _hZtCstVsPtJetH[iPtTest] -> Write();
  }
  for (UInt_t iNcstTest = 0; iNcstTest < 4; iNcstTest++) {
    _hZtCstVsNcstP[iNcstTest] -> Write();
    _hZtCstVsNcstG[iNcstTest] -> Write();
    _hZtCstVsNcstH[iNcstTest] -> Write();
  }
  _hNcstVsPtJetP -> Write();
  _hNcstVsPtJetG -> Write();
  _hNcstVsPtJetH -> Write();


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
  _tMc     -> GetCurrentFile() -> cd();
  _tMc     -> GetCurrentFile() -> Close();
  PrintInfo(16);

}  // end 'Finish()'

// End ------------------------------------------------------------------------
