// 'StMuDstJetTreeMaker.cxx'
// Derek Anderson, Nihar Sahoo
// 06.23.2017
//
// This class reads the 'Gfmtodst' tree
// and produces a tree of jets


#define StMuDstJetTreeMaker_cxx

#include "fastjet/config.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "StMuDstJetTreeMaker.h"
#include "StMuDstJetTreeMaker.utilities.h"

using namespace std;
using namespace fastjet;

ClassImp(StMuDstJetTreeMaker);



void StMuDstJetTreeMaker::SetInputAndOutputFiles(const TString& sInput, const TString& sOutput, const Double_t pTparton) {

  _sInput     = sInput;
  _sOutput    = sOutput;
  _PartonicPt = pTparton;
  PrintInfo(1);

}  // end 'SetInputAndOutput(TString&, TString&)'



void StMuDstJetTreeMaker::SetEventParameters(const Double_t rVtxMax, const Double_t zVtxMax) {

  _rVtxMax = rVtxMax;
  _zVtxMax = zVtxMax;
  PrintInfo(2);

}  // end 'StEventParameters(Double_t, Double_t)'



void StMuDstJetTreeMaker::SetTriggerParameters(const Int_t adcMax, const Double_t eEtaMin, const Double_t ePhiMin, const Double_t pProjMax, const Double_t etaTrgMax, const Double_t eTtrgMin, const Double_t eTtrgMax, const Double_t tspPi0Min, const Double_t tspPi0Max, const Double_t tspGamMin, const Double_t tspGamMax) {

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



void StMuDstJetTreeMaker::SetTrackParameters(const UInt_t nFitMin, const Double_t rFitMin, const Double_t dcaMax, const Double_t etaTrkMax, const Double_t pTtrkMin, const Double_t pTtrkMax) {

  _nFitMin   = nFitMin;
  _rFitMin   = rFitMin;
  _dcaMax    = dcaMax;
  _etaTrkMax = etaTrkMax;
  _pTtrkMin  = pTtrkMin;
  _pTtrkMax  = pTtrkMax;
  PrintInfo(4);

}  // end 'SetTrackParameters(UInt_t, Double_t, Double_t, Double_t, Double_t, Double_t)'



void StMuDstJetTreeMaker::SetTowerParameters(const Double_t etaTwrMax, const Double_t eTwrMin, const Double_t eTwrMax, const Double_t eCorrMin, const Double_t eCorrMax) {

  _etaTwrMax = etaTwrMax;
  _eTwrMin   = eTwrMin;
  _eTwrMax   = eTwrMax;
  _eCorrMin  = eCorrMin;
  _eCorrMax  = eCorrMax;
  PrintInfo(5);

}  // 'SetTowerParameters(Double_t, Double_t, Double_t, Double_t, Double_t)'



void StMuDstJetTreeMaker::SetJetParameters(const UInt_t type, const UInt_t nRepeat, const UInt_t nRemove, const Double_t rJet, const Double_t aGhost, const Double_t pTjetMin, const Double_t etaGhostMax, const Double_t etaJetMax, const Double_t etaBkgdMax) {

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



void StMuDstJetTreeMaker::AdjustTrackEfficiency(const Bool_t effAdjust, const Float_t adjustment) {

  _adjustTrackEff = effAdjust;
  _effAdjust      = adjustment;
  _random         = new TRandom();
  PrintInfo(17);

}  // 'AdjustTrackEfficiency(Bool_t, Float_t)'



void StMuDstJetTreeMaker::Init(const UInt_t trgFlag) {

  _trgFlag = trgFlag;
  _fOutput = new TFile(_sOutput.Data(), "recreate");
  if (_tFemto == 0) {
    TFile *fInput = (TFile*) gROOT -> GetListOfFiles() -> FindObject(_sInput.Data());
    if (!fInput || !(fInput -> IsOpen())) {
      fInput = new TFile(_sInput.Data());
    }
    switch (_trgFlag) {
      case 0:
        fInput -> GetObject("GfmtoDst_mu", _tFemto);
        break;
      case 1:
        fInput -> GetObject("GfmtoDst_mu", _tFemto);
        break;
      case 2:
        fInput -> GetObject("GfmtoDst_mu", _tFemto);
        break;
    }
  }
  InitializeInputTree(_tFemto);
  InitializeOutputTree(_tJet);
  InitializeHistograms();
  PrintInfo(10);

}  // end 'Init()'



void  StMuDstJetTreeMaker::Make() {

  switch (_trgFlag) {
    case 0:
      UsePionTriggeredEvents();
      break;
    case 1:
      UseHadronTriggeredEvents();
      break;
    case 2:
      UseNonTriggeredEvents();
      break;
  }

}  // end 'Make(UInt_t)'



void StMuDstJetTreeMaker::UsePionTriggeredEvents() {

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
    const Bool_t isGoodTrg = IsGoodPionTrigger(adc, eEta4, ePhi4, pProj, hTwr, eTtrg, tsp);
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

      // adjust tracking efficiency (if need be)
      if (_adjustTrackEff) {
        const Float_t rando = _random -> Uniform(0., 1.);
        const Bool_t  pass  = (rando > _effAdjust);
        if (!pass) continue;
      }

      // fill histograms and sum event pT for good tracks
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
      if (_jetType == 1)
        particles.push_back(PseudoJet(pXtwr, pYtwr, pZtwr, eCorr));

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
    _JetConsPt.clear();
    _JetConsEta.clear();
    _JetConsPhi.clear();
    _JetConsE.clear();

    // jet event info
    _EventIndex = _eventNumber;
    _NJets      = (Int_t) jets.size();
    _RunId      = _runNumber;
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

    }  // end jet loop

    _tJet -> Fill();

  }  // end event loop

  // record no. of triggers
  _hEvtQA[0][0] -> Fill(nPi0);
  _hEvtQA[0][1] -> Fill(nGam);
  PrintInfo(13);


  // normalize histograms
  for (UInt_t iTrg = 0; iTrg < (NTrgTypes - 1); iTrg++) {
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
  }
  PrintInfo(14);

}  // end 'UsePionTriggeredEvents()'



void StMuDstJetTreeMaker::UseHadronTriggeredEvents() {

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
    const Long64_t nTwrs = _TowerArray_;
    const Double_t xVtx  = _xVertex;
    const Double_t yVtx  = _yVertex;
    const Double_t zVtx  = _zVertex;
    const Double_t rVtx  = sqrt((xVtx * xVtx) + (yVtx * yVtx));

    // event cuts
    const Bool_t isGoodRun = IsGoodRunID(runID);
    const Bool_t isGoodEvt = IsGoodEvent(rVtx, zVtx);
    if (!isGoodRun || !isGoodEvt) continue;


    // trigger loop
    UInt_t   idTrg(nTrks + 1);
    Bool_t   foundTrg(false);
    Double_t eTrg(0.);
    Double_t fTrg(0.);
    Double_t hTrg(0.);
    Double_t pTtrg(0.);
    Double_t nSigElec(999.);
    for (UInt_t iTrg = 0; iTrg < nTrks; iTrg++) {

      // track info
      const UInt_t   nFit   = _PrimaryTrackArray_nHitsFit[iTrg];
      const UInt_t   nPoss  = _PrimaryTrackArray_nHitsPoss[iTrg];
      const Double_t rFit   = (Double_t) nFit / (Double_t) nPoss;
      const Double_t dca    = _PrimaryTrackArray_dcag[iTrg];
      const Double_t hTrk   = _PrimaryTrackArray_eta[iTrg];
      const Double_t fTrk   = _PrimaryTrackArray_phi[iTrg];
      const Double_t pZtrk  = _PrimaryTrackArray_pZ[iTrg];
      const Double_t pTtrk  = _PrimaryTrackArray_pT[iTrg];
      const Double_t eTrk   = sqrt((pTtrk * pTtrk) + (pZtrk * pZtrk) + (MassPi0 * MassPi0));
      const Double_t nSigPi = _PrimaryTrackArray_nSigPion[iTrg];
      const Double_t nSigK  = _PrimaryTrackArray_nSigKaon[iTrg];
      const Double_t nSigP  = _PrimaryTrackArray_nSigProton[iTrg];
      const Double_t nSigE  = _PrimaryTrackArray_nSigElectron[iTrg];

      // track cuts
      const Bool_t isInFitCut   = (nFit >= _nFitMin);
      const Bool_t isInRatioCut = (rFit >= _rFitMin);
      const Bool_t isInDcaCut   = (dca < _dcaMax);
      const Bool_t isInEtaCut   = (TMath::Abs(hTrk) < _etaTrgMax);
      const Bool_t isInPtCut    = ((pTtrk > _eTtrgMin) && (pTtrk < _eTtrgMax));
      const Bool_t isGoodTrk    = (isInFitCut && isInRatioCut && isInDcaCut && isInEtaCut && isInPtCut);

      // trigger cuts
      const Bool_t isGoodTrg = IsGoodHadronTrigger(hTrk, pTtrk, nSigPi, nSigK, nSigP, nSigE);
      if (!isGoodTrk || !isGoodTrg)
        continue;
      else {
        idTrg    = iTrg;
        eTrg     = eTrk;
        fTrg     = fTrk;
        hTrg     = hTrk;
        pTtrg    = pTtrk;
        nSigElec = nSigE;
        foundTrg = true;
        break;
      }

    }  // end trigger loop

    if (!foundTrg)
      continue;
    else {
      _hEvtQA[1][2] -> Fill(nSigElec);
      _hEvtQA[2][2] -> Fill(nTrks);
      _hEvtQA[3][2] -> Fill(eTrg);
      nHad++;
    }


    // for event pT
    Double_t pTevtH[NTrkTypes];
    for (UInt_t iTrkType = 0; iTrkType < NTrkTypes; iTrkType++) {
      pTevtH[iTrkType] = 0.;
    }

    // track loop 
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
      const Bool_t isTrigger = (iTrk == idTrg);
      if (!isTrigger) pTevtH[0] += pTtrk;

      // track cuts
      const Bool_t isGoodTrk = IsGoodTrack(nFit, rFit, dca, hTrk, pTtrk);
      if (isTrigger || !isGoodTrk) continue;

      // adjust tracking efficiency (if need be)
      if (_adjustTrackEff) {
        const Float_t rando = _random -> Uniform(0., 1.);
        const Bool_t  pass  = (rando > _effAdjust);
        if (!pass) continue;
      } 

      // fill histograms and sum evet pT for good tracks
      pTevtH[1] += pTtrk;
      _hTrkQA[0][2] -> Fill(pTtrk);
      _hTrkQA[1][2] -> Fill(dFtrk);
      _hTrkQA[2][2] -> Fill(hTrk);
      _hTrkQA[3][2] -> Fill(eTrk);
      particles.push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));

    }  // end track loop

    // fill event pT
    for (UInt_t iTrkType = 0; iTrkType < NTrkTypes; iTrkType++) {
      _hEvtPt[2][iTrkType] -> Fill(pTevtH[iTrkType]);
    }


    // tower loop
    const TVector3 vVtx(xVtx, yVtx, zVtx);
    for (UInt_t iTwr = 0; iTwr < nTwrs; iTwr++) {

      // tower info
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
      const Bool_t isGoodTwr = IsGoodTower(hTwr, eTwr, eCorr);
      if (!isGoodTwr) continue;


      // get tower momentum
      const TLorentzVector vTwr = GetTowerMomentumVector(RadiusBemc, hTwr, fTwr, eCorr, vVtx);

      const Double_t pXtwr = vTwr.Px();
      const Double_t pYtwr = vTwr.Py();
      const Double_t pZtwr = vTwr.Pz();
      const Double_t pTtwr = vTwr.Pt();
      _hTwrQA[0][2] -> Fill(pTtwr);
      _hTwrQA[1][2] -> Fill(dFtwr);
      _hTwrQA[2][2] -> Fill(hTwr);
      _hTwrQA[3][2] -> Fill(eCorr);
      if (_jetType == 1)
        particles.push_back(PseudoJet(pXtwr, pYtwr, pZtwr, eCorr));

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
    _JetConsPt.clear();
    _JetConsEta.clear();
    _JetConsPhi.clear();
    _JetConsE.clear();

    // jet event info
    _EventIndex = _eventNumber;
    _NJets      = (Int_t) jets.size();
    _RunId      = _runNumber;
    _Refmult    = _refMult;
    _TSP        = nSigElec;
    _TrgEta     = hTrg;
    _TrgPhi     = fTrg;
    _TrgEt      = pTtrg;
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

      // TEST [03.30.2020]
      UInt_t nCstReal = 0;
      for (UInt_t iCst = 0; iCst < nCst; iCst++) {
        const Double_t pTcst = jets[iJet].constituents()[iCst].perp();
        if (pTcst > 0.2) nCstReal++;
      }

      // TEST [04.04.2020]
      for (UInt_t iCst = 0; iCst < nCst; iCst++) {
        const Double_t pTcst = jets[iJet].constituents()[iCst].perp();
        const Double_t zTcst = pTcst / pTjet;
        if ((pTcst > 0.2) && (pTjet > 0.2)) {

          // zTcst vs. pTjet
          if ((pTjet >= 0.2) && (pTjet < 1.))
            _hZtCstVsPtJet[0] -> Fill(zTcst);
          if ((pTjet >= 1.) && (pTjet < 2.))
            _hZtCstVsPtJet[1] -> Fill(zTcst);
          if ((pTjet >= 2.) && (pTjet < 5.))
            _hZtCstVsPtJet[2] -> Fill(zTcst);
          if ((pTjet >= 5.) && (pTjet < 10.))
            _hZtCstVsPtJet[3] -> Fill(zTcst);
          if (pTjet >= 10.)
            _hZtCstVsPtJet[4] -> Fill(zTcst);

          // zTcst vs. nCstJet
          if ((nCstReal >= 1) && (nCstReal < 2))
            _hZtCstVsNcst[0] -> Fill(zTcst);
          if ((nCstReal >= 2) && (nCstReal < 5))
            _hZtCstVsNcst[1] -> Fill(zTcst);
          if ((nCstReal >= 5) && (nCstReal < 10))
            _hZtCstVsNcst[2] -> Fill(zTcst);
          if (nCstReal >= 10)
            _hZtCstVsNcst[3] -> Fill(zTcst);
        }
      }
      if (pTjet > 0.2) _hNcstVsPtJet -> Fill(pTjet, nCstReal);


      _hJetQA[0][2] -> Fill(pTjet);
      _hJetQA[1][2] -> Fill(dFjet);
      _hJetQA[2][2] -> Fill(hJet);
      _hJetQA[3][2] -> Fill(eJet);
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
  _hEvtQA[0][2] -> Fill(nHad);
  PrintInfo(13);


  // normalize histograms
  for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
    const Double_t trkBin   = _hTrkQA[iHist][2] -> GetBinWidth(17);
    const Double_t twrBin   = _hTwrQA[iHist][2] -> GetBinWidth(17);
    const Double_t jetBin   = _hJetQA[iHist][2] -> GetBinWidth(17);
    const Double_t trkScale = 1. / (nHad * trkBin);
    const Double_t twrScale = 1. / (nHad * twrBin);
    const Double_t jetScale = 1. / (nHad * jetBin);
    _hTrkQA[iHist][2] -> Scale(trkScale);
    _hTwrQA[iHist][2] -> Scale(twrScale);
    _hJetQA[iHist][2] -> Scale(jetScale);
  }
  PrintInfo(14);

}  // end 'UseHadronTriggeredEvents()'



void StMuDstJetTreeMaker::UseNonTriggeredEvents() {

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
  Long64_t nNon   = 0;
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

    // count events
    Double_t eTrg(0.);
    Double_t fTrg(0.);
    Double_t hTrg(0.);
    Double_t pTtrg(0.);
    Double_t nSigElec(999.);
    _hEvtQA[1][3] -> Fill(nSigElec);
    _hEvtQA[2][3] -> Fill(nTrks);
    _hEvtQA[3][3] -> Fill(eTrg);
    nNon++;


    // for event pT
    Double_t pTevt0[NTrkTypes];
    for (UInt_t iTrkType = 0; iTrkType < NTrkTypes; iTrkType++) {
      pTevt0[iTrkType] = 0.;
    }

    // track loop 
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
      pTevt0[0] += pTtrk;

      // track cuts
      const Bool_t isGoodTrk = IsGoodTrack(nFit, rFit, dca, hTrk, pTtrk);
      if (!isGoodTrk) continue;

      // adjust tracking efficiency (if need be)
      if (_adjustTrackEff) {
        const Float_t rando = _random -> Uniform(0., 1.);
        const Bool_t  pass  = (rando > _effAdjust);
        if (!pass) continue;
      } 


      // fill histograms and sum event pT for good tracks
      pTevt0[1] += pTtrk;
      _hTrkQA[0][3] -> Fill(pTtrk);
      _hTrkQA[1][3] -> Fill(dFtrk);
      _hTrkQA[2][3] -> Fill(hTrk);
      _hTrkQA[3][3] -> Fill(eTrk);
      particles.push_back(PseudoJet(pXtrk, pYtrk, pZtrk, eTrk));

    }  // end track loop

    // fill event pT
    for (UInt_t iTrkType = 0; iTrkType < NTrkTypes; iTrkType++) {
      _hEvtPt[3][iTrkType] -> Fill(pTevt0[iTrkType]);
    }


    // tower loop
    const TVector3 vVtx(xVtx, yVtx, zVtx);
    for (UInt_t iTwr = 0; iTwr < nTwrs; iTwr++) {

      // tower info
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
      const Bool_t isGoodTwr = IsGoodTower(hTwr, eTwr, eCorr);
      if (!isGoodTwr) continue;


      // get tower momentum
      const TLorentzVector vTwr = GetTowerMomentumVector(RadiusBemc, hTwr, fTwr, eCorr, vVtx);

      const Double_t pXtwr = vTwr.Px();
      const Double_t pYtwr = vTwr.Py();
      const Double_t pZtwr = vTwr.Pz();
      const Double_t pTtwr = vTwr.Pt();
      _hTwrQA[0][3] -> Fill(pTtwr);
      _hTwrQA[1][3] -> Fill(dFtwr);
      _hTwrQA[2][3] -> Fill(hTwr);
      _hTwrQA[3][3] -> Fill(eCorr);
      if (_jetType == 1)
        particles.push_back(PseudoJet(pXtwr, pYtwr, pZtwr, eCorr));

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
    _JetConsPt.clear();
    _JetConsEta.clear();
    _JetConsPhi.clear();
    _JetConsE.clear();

    // jet event info
    _EventIndex = _eventNumber;
    _NJets      = (Int_t) jets.size();
    _RunId      = _runNumber;
    _Refmult    = _refMult;
    _TSP        = nSigElec;
    _TrgEta     = hTrg;
    _TrgPhi     = fTrg;
    _TrgEt      = pTtrg;
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


      _hJetQA[0][3] -> Fill(pTjet);
      _hJetQA[1][3] -> Fill(dFjet);
      _hJetQA[2][3] -> Fill(hJet);
      _hJetQA[3][3] -> Fill(eJet);
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
  _hEvtQA[0][3] -> Fill(nNon);
  PrintInfo(13);


  // normalize histograms
  for (UInt_t iHist = 0; iHist < NHistQA; iHist++) {
    const Double_t trkBin   = _hTrkQA[iHist][3] -> GetBinWidth(17);
    const Double_t twrBin   = _hTwrQA[iHist][3] -> GetBinWidth(17);
    const Double_t jetBin   = _hJetQA[iHist][3] -> GetBinWidth(17);
    const Double_t trkScale = 1. / (nNon * trkBin);
    const Double_t twrScale = 1. / (nNon * twrBin);
    const Double_t jetScale = 1. / (nNon * jetBin);
    _hTrkQA[iHist][3] -> Scale(trkScale);
    _hTwrQA[iHist][3] -> Scale(twrScale);
    _hJetQA[iHist][3] -> Scale(jetScale);
  }
  PrintInfo(14);

}  // end 'UseNonTriggeredEvents()'



void StMuDstJetTreeMaker::Finish() {

  _fOutput -> cd();
  TDirectory *dHists = _fOutput -> mkdir("QA");
  TDirectory *dPi0QA = dHists   -> mkdir("Pi0");
  TDirectory *dGamQA = dHists   -> mkdir("Gamma");
  TDirectory *dHadQA = dHists   -> mkdir("Hadron");
  TDirectory *dNonQA = dHists   -> mkdir("NoTrigger");
  PrintInfo(15);

  // TEST [04.04.2020]
  _fOutput -> cd();
  for (UInt_t iPtTest = 0; iPtTest < 5; iPtTest++) {
    _hZtCstVsPtJet[iPtTest] -> Write();
  }
  for (UInt_t iNcstTest = 0; iNcstTest < 4; iNcstTest++) {
    _hZtCstVsNcst[iNcstTest] -> Write();
  }
  _hNcstVsPtJet -> Write();

  for (UInt_t iTrg = 0; iTrg < NTrgTypes; iTrg++) {
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
      case 3:
        dNonQA -> cd();
        break;
    }
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
