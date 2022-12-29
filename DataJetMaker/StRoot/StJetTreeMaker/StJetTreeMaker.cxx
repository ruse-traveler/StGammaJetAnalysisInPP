// 'StJetTreeMaker.cxx'
// Derek Anderson, Nihar Sahoo
// 06.23.2017
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets


#define StJetTreeMaker_cxx

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StJetTreeMaker.h"
#include "StJetTreeMaker.utilities.h"

using namespace std;
using namespace fastjet;

ClassImp(StJetTreeMaker);



void StJetTreeMaker::SetInputAndOutputFiles(const TString& sInput, const TString& sOutput) {

  _sInput  = sInput;
  _sOutput = sOutput;
  PrintInfo(1);

}  // end 'SetInputAndOutput(TString&, TString&)'



void StJetTreeMaker::SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax) {

  _rVtxMax = rVtxMax;
  _zVtxMax = zVtxMax;
  PrintInfo(2);

}  // end 'StEventParameters(Double_t, Double_t)'



void StJetTreeMaker::SetTriggerParameters(const Int_t adcMax, const Double_t eEtaMin, const Double_t ePhiMin, const Double_t pProjMax, const Double_t etaTrgMax, const Double_t eTtrgMin, const Double_t eTtrgMax, const Double_t tspPi0Min, const Double_t tspPi0Max, const Double_t tspGamMin, const Double_t tspGamMax) {

  _adcMax    = adcMax;
  _eEtaMin   = eEtaMin;
  _ePhiMin   = ePhiMin;
  _pProjMax  = pProjMax;
  _etaTrgMax = etaTrgMax;
  _eTtrgMin  = eTtrgMin;
  _eTtrgMax  = eTtrgMax;
  _tspPi0Min = tspPi0Min;
  _tspPi0Max = tspPi0Max;
  _tspGamMin = tspGamMin;
  _tspGamMax = tspGamMax;
  PrintInfo(3);

}  // end 'SetTriggerParameters(Int_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



void StJetTreeMaker::SetTrackParameters(const UInt_t nFitMin, const Double_t rFitMin, const Double_t dcaMax, const Double_t etaTrkMax, const Double_t pTtrkMin, const Double_t pTtrkMax) {

  _nFitMin   = nFitMin;
  _rFitMin   = rFitMin;
  _dcaMax    = dcaMax;
  _etaTrkMax = etaTrkMax;
  _pTtrkMin  = pTtrkMin;
  _pTtrkMax  = pTtrkMax;
  PrintInfo(4);

}  // end 'SetTrackParameters(UInt_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



void StJetTreeMaker::SetTowerParameters(const Double_t etaTwrMax, const Double_t eTwrMin, const Double_t eTwrMax, const Double_t eCorrMin, const Double_t eCorrMax) {

  _etaTwrMax = etaTwrMax;
  _eTwrMin   = eTwrMin;
  _eTwrMax   = eTwrMax;
  _eCorrMin  = eCorrMin;
  _eCorrMax  = eCorrMax;
  PrintInfo(5);

}  // 'SetTowerParameters(Double_t, Double_t, Double_t, Double_t, Double_t)'



void StJetTreeMaker::SetJetParameters(const UInt_t type, const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t aGhost, const Double_t pTjetMin, const Double_t etaGhostMax, const Double_t etaJetMax, const Double_t etaBkgdMax) {

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



void StJetTreeMaker::Init() {

  _fOutput = new TFile(_sOutput.Data(), "recreate");
  if (_tFemto == 0) {
    TFile *fInput = (TFile*) gROOT -> GetListOfFiles() -> FindObject(_sInput.Data());
    if (!fInput || !(fInput -> IsOpen())) {
      fInput = new TFile(_sInput.Data());
    }
    fInput -> GetObject("Gfmtodst", _tFemto);
  }
  InitializeInputTree(_tFemto);
  InitializeOutputTree(_tJet);
  InitializeHistograms();
  PrintInfo(10);

}  // end 'Init()'



void StJetTreeMaker::Make() {

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
    const Long64_t nTwrs = _TowerArray_;
    const Double_t xVtx  = _xVertex;
    const Double_t yVtx  = _yVertex;
    const Double_t zVtx  = _zVertex;
    const Double_t rVtx  = sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // event cuts
    const Bool_t isGoodRun = IsGoodRunID(runID);
    const Bool_t isGoodEvt = IsGoodEvent(rVtx, zVtx);
    if (!isGoodRun || !isGoodEvt) continue;


    // trigger info
    const Int_t    adc   = _ETwradc11;
    const UInt_t   twrID = _ETwrdidT;
    const Double_t tsp   = _Etsp;
    const Double_t eEta4 = _EEstrpen4;
    const Double_t ePhi4 = _EPstripenp4;
    const Double_t pProj = _ETwrPTower;
    const Double_t hTwr  = _ETwreT;
    const Double_t hTrg  = _EClustetav1;
    const Double_t fTrg  = _EClustphiv1;
    const Double_t eTrg  = _EClustEneT0;
    const Double_t eTtrg = eTrg * sin(2. * atan(exp(-1. * hTrg)));

    // trigger cuts
    const Bool_t isPi0     = IsPi0(tsp);
    const Bool_t isGam     = IsGamma(tsp);
    const Bool_t isGoodTwr = IsGoodTowerID(twrID);
    const Bool_t isGoodTrg = IsGoodTrigger(adc, eEta4, ePhi4, pProj, hTwr, hTrg, eTtrg, tsp);
    if (!isGoodTwr || !isGoodTrg) continue;

    // count triggers
    if (isPi0) {
      _hEvtQA[1][0] -> Fill(tsp);
      _hEvtQA[2][0] -> Fill(nTrks);
      _hEvtQA[3][0] -> Fill(eTrg);
      nPi0++;
    }
    if (isGam) {
      _hEvtQA[1][1] -> Fill(tsp);
      _hEvtQA[2][1] -> Fill(nTrks);
      _hEvtQA[3][1] -> Fill(eTrg);
      nGam++;
    }


    // for event pT
    Double_t pTevtP[NTrkTypes];
    Double_t pTevtG[NTrkTypes];
    for (UInt_t iTrkType = 0; iTrkType < NTrkTypes; iTrkType++) {
      pTevtP[iTrkType] = 0.;
      pTevtG[iTrkType] = 0.;
    }

    // track loop
    UInt_t iCstAdd = 0;
    for (UInt_t iTrk = 0; iTrk < nTrks; iTrk++) {

      // track info
      const UInt_t   nFit  = _PrimaryTrackArray_nHitsFit[iTrk];
      const UInt_t   nPoss = _PrimaryTrackArray_nHitsPoss[iTrk];
      const Double_t rFit  = (Double_t) nFit / (Double_t) nPoss;
      const Double_t dca   = _PrimaryTrackArray_dcag[iTrk];
      const Double_t hTrk  = _PrimaryTrackArray_eta[iTrk];
      const Double_t fTrk  = _PrimaryTrackArray_phi[iTrk];
      const Double_t pXtrk = _PrimaryTrackArray_pX[iTrk];
      const Double_t pYtrk = _PrimaryTrackArray_pY[iTrk];
      const Double_t pZtrk = _PrimaryTrackArray_pZ[iTrk];
      const Double_t pTtrk = _PrimaryTrackArray_pT[iTrk];
      const Double_t eTrk  = sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (MassPi0 * MassPi0));

      Double_t dFtrk = fTrk - fTrg;
      if (dFtrk < (-1. * TMath::PiOver2()))
        dFtrk += TMath::TwoPi();
      if (dFtrk > (3. * TMath::PiOver2()))
        dFtrk -= TMath::TwoPi();

      // sum event pT for all primary tracks
      if (isPi0) pTevtP[0] += pTtrk;
      if (isGam) pTevtG[0] += pTtrk;

      // track cuts
      const Bool_t isGoodTrk = IsGoodTrack(nFit, rFit, dca, hTrk, pTtrk);
      if (!isGoodTrk) continue;


      if (isPi0) {
        _hTrkQA[0][0] -> Fill(pTtrk);
        _hTrkQA[1][0] -> Fill(dFtrk);
        _hTrkQA[2][0] -> Fill(hTrk);
        _hTrkQA[3][0] -> Fill(eTrk);
        pTevtP[1] += pTtrk;
      }
      if (isGam) {
        _hTrkQA[0][1] -> Fill(pTtrk);
        _hTrkQA[1][1] -> Fill(dFtrk);
        _hTrkQA[2][1] -> Fill(hTrk);
        _hTrkQA[3][1] -> Fill(eTrk);
        pTevtG[1] += pTtrk;
      }
      particles.push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));
      particles[iCstAdd].set_user_index(FlagChrg);
      iCstAdd++;

    }  // end track loop

    // fill event pT
    for (UInt_t iTrkType = 0; iTrkType < NTrkTypes; iTrkType++) {
      if (isPi0) _hEvtPt[0][iTrkType] -> Fill(pTevtP[iTrkType]);
      if (isGam) _hEvtPt[1][iTrkType] -> Fill(pTevtG[iTrkType]);
    }


    // tower loop
    const TVector3 vVtx(xVtx, yVtx, zVtx);
    for (UInt_t iTwr = 0; iTwr < nTwrs; iTwr++) {

      // tower info
      const Int_t    tID    = _TowerArray_TwrId[iTwr];
      const UInt_t   fMatch = _TowerArray_TwrMatchIdnex[iTwr];
      const UInt_t   nMatch = _TowerArray_NoOfmatchedTrk[iTwr];
      const Double_t hTwr   = _TowerArray_TwrEta[iTwr];
      const Double_t fTwr   = _TowerArray_TwrPhi[iTwr];
      const Double_t eTwr   = _TowerArray_TwrEng[iTwr];

      Double_t eCorr = -999.;
      Double_t dFtwr = fTwr - fTrg;
      if (dFtwr < (-1. * TMath::PiOver2()))
        dFtwr += TMath::TwoPi();
      if (dFtwr > (3. * TMath::PiOver2()))
        dFtwr -= TMath::TwoPi();

      // calculate corrected energy
      const Bool_t isMatched    = (fMatch == 1);
      const Bool_t isNotMatched = (fMatch == 0);
      const Bool_t hasMatches   = (nMatch >= 1);
      const Bool_t hasNoMatches = (nMatch == 0);
      matches.clear();
      if (isMatched && hasMatches) {
        for (UInt_t iMatch = 0; iMatch < nMatch; iMatch++) {
          const Double_t pMatch = _TowerArray_fMatchedTracksArray_P[iTwr][iMatch];
          matches.push_back(pMatch);
        }
        eCorr = GetHadronicCorrection(eTwr, matches);
      }
      else if (isNotMatched && hasNoMatches) {
        eCorr = eTwr;
      }

      // tower cuts
      const Bool_t isTrigger = (tID == twrID);
      const Bool_t isGoodTwr = IsGoodTower(hTwr, eTwr, eCorr);
      if (!isGoodTwr) continue;


      // get tower momentum
      const TLorentzVector vTwr = GetTowerMomentumVector(RadiusBemc, hTwr, fTwr, eCorr, vVtx);

      const Double_t pXtwr = vTwr.Px();
      const Double_t pYtwr = vTwr.Py();
      const Double_t pZtwr = vTwr.Pz();
      const Double_t pTtwr = vTwr.Pt();
      if (isPi0 && !isTrigger) {
        _hTwrQA[0][0] -> Fill(pTtwr);
        _hTwrQA[1][0] -> Fill(dFtwr);
        _hTwrQA[2][0] -> Fill(hTwr);
        _hTwrQA[3][0] -> Fill(eCorr);
      }
      if (isGam && !isTrigger) {
        _hTwrQA[0][1] -> Fill(pTtwr);
        _hTwrQA[1][1] -> Fill(dFtwr);
        _hTwrQA[2][1] -> Fill(hTwr);
        _hTwrQA[3][1] -> Fill(eCorr);
      }
      if (_jetType == 1) {
        particles.push_back(PseudoJet(pXtwr, pYtwr, pZtwr, eCorr));
        particles[iCstAdd].set_user_index(FlagNeut);
        iCstAdd++;
      }

    }  // end tower loop


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
    _JetRt.clear();
    _JetConsPt.clear();
    _JetConsEta.clear();
    _JetConsPhi.clear();
    _JetConsE.clear();

    // jet event info
    _EventIndex = _eventNumber;
    _NJets      = (Int_t) jets.size();
    _Refmult    = _refMult;
    _TSP        = _Etsp;
    _TrgEta     = _EClustetav1;
    _TrgPhi     = _EClustphiv1;
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


      // off axis cone calculation
      Double_t pTsum[2]     = {0., 0.};
      Double_t pTdensity[2] = {0., 0.};
      GetOffAxisTrackPtSum(nTrks, fTrg, dFjet, hJet, pTsum[0], pTsum[1]);
      if (_jetType == 1)
        GetOffAxisTowerPtSum(nTwrs, twrID, fTrg, dFjet, hJet, vVtx, matches, pTsum[0], pTsum[1]);

      const Double_t aCone = (_rJet * _rJet) * TMath::Pi();
      pTdensity[0] = pTsum[0] / aCone;
      pTdensity[1] = pTsum[1] / aCone;

      // rT calculation
      Double_t rTjet  = -1.;
      if (_jetType == 1) {
        Double_t pTchrg = 0.;
        Double_t pTneut = 0.;
        for (UInt_t iCst = 0; iCst < nCst; iCst++) {
          Int_t    flag  = jets[iJet].constituents()[iCst].user_index();
          Double_t pTcst = jets[iJet].constituents()[iCst].perp();
          switch (flag) {
            case 0:
              pTneut += pTcst;
              break;
            case 1:
              pTchrg += pTcst;
              break;
          }
        }  // end constituent loop
        rTjet = pTneut / (pTchrg + pTneut);
      }


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
      _JetIndex.push_back(iJet);
      _JetNCons.push_back(nCst);
      _JetPt.push_back(pTjet);
      _JetPtCorr.push_back(pTcorr);
      _JetEta.push_back(hJet);
      _JetPhi.push_back(fJet);
      _JetE.push_back(eJet);
      _JetArea.push_back(aJet);
      _JetRt.push_back(rTjet);
      _JetPtOffAxisUp.push_back(pTdensity[0]);
      _JetPtOffAxisDown.push_back(pTdensity[1]);

    }  // end jet loop

    _tJet -> Fill();

  }  // end event loop

  // record no. of triggers
  _hEvtQA[0][0] -> Fill(nPi0);
  _hEvtQA[0][1] -> Fill(nGam);
  PrintInfo(13);


  // normalize histograms
  for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
    Double_t nTrg;
    if (iTrg == 0)
      nTrg = (Double_t) nPi0;
    else
      nTrg = (Double_t) nGam;

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
    for (UInt_t iTrk = 0; iTrk < NTrkTypes; iTrk++) {
      const Double_t evtBin   = _hEvtPt[iTrg][iTrk] -> GetBinWidth(17);
      const Double_t evtScale = 1. / (nTrg * evtBin);
      _hEvtPt[iTrg][iTrk] -> Scale(evtScale);
    }
  }
  PrintInfo(14);

}  // end 'Make()'



void StJetTreeMaker::Finish() {

  _fOutput -> cd();
  TDirectory *dHists = _fOutput -> mkdir("QA");
  TDirectory *dPi0QA = dHists   -> mkdir("Pi0");
  TDirectory *dGamQA = dHists   -> mkdir("Gamma");
  PrintInfo(15);

  for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
    if (iTrg == 0)
      dPi0QA -> cd();
    else
      dGamQA -> cd();
    for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
      _hEvtQA[iHist][iTrg] -> Write();
      _hTrkQA[iHist][iTrg] -> Write();
      _hTwrQA[iHist][iTrg] -> Write();
      _hJetQA[iHist][iTrg] -> Write();
    }
    for (UInt_t iTrk = 0; iTrk < NTrkTypes; iTrk++) {
      _hEvtPt[iTrg][iTrk] -> Write();
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
