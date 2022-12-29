//________________________________________________________________
// 
//
//Date 05/10/2016: This is to read geant.root file
//                      -Nihar
//Date 07/18/2017: added E/p calculation in MuDst
//                 files.
//                      -Derek 
//Date  07/20/2017: added additional branches to calculate E/p for
//                  electrons and hadrons.
//                      -Derek
//Date  02/19/2018: added tree and variables to hold StMuMcTrack
//                  variables
//                      -Derek
//____________________________________________________________


#include "StJetTreeMcMaker.h"
#include "StChain.h"
#include "StEvent.h"
#include "StEventTypes.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StEmcUtil/filters/StEmcFilter.h"
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include <vector>
#include "TVector.h"
#include "TVector3.h"
#include <math.h>
#include "TNtuple.h"
#include <TString.h>
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "StEventUtilities/StuRefMult.hh"
#include "StEventUtilities/StuFtpcRefMult.hh"
#include "StThreeVector.hh"
#include "StEvent/StEnumerations.h"
#include "StMcEvent.hh"
#include "StMcVertex.hh"
#include "StMcTrack.hh"
#include "StMcEventMaker/StMcEventMaker.h"
#include "StMcEventTypes.hh"
#include "StPrimaryVertex.h"
#include "StBTofHeader.h"
//#include "StMuEvent.h"
#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>

#include "St_db_Maker/St_db_Maker.h"
#include "StDetectorDbMaker/StDetectorDbTriggerID.h"

#include <emcStatus.h>
#include "tables/St_emcStatus_Table.h"

#include "tables/St_emcPed_Table.h"
//#include "L2Result.h"
#include "StTriggerData.h"
#include "StTriggerData2007.h"
//#include "StDaqLib/TRG/trgStructures.h"
#include "StDaqLib/TRG/trgStructures2007.h"
//#include "StDaqLib/TRG/L2pedResults2006.h"
//#include "StDaqLib/TRG/L2gammaResult2007.h"
#include "StDaqLib/EMC/StEmcDecoder.h"
//#include "L2gammaResult.h"
#include "StEvent/StTriggerIdCollection.h"

#include "StEvent/StTriggerId.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuTofHit.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"  // [Derek, 02.19.2018]
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"  // [Derek, 02.20.2018]
#include "StTriggerData.h"

 
#include "GammaJetTrack.h"
#include "GammaJetEvent.h"
#include "GammaJetTower.h"
#include "GammaJetTowerUtil.h"
#include <exception>

//void muEventInfo(StMuEvent&, const Int_t&);
ClassImp(StJetTreeMcMaker)

// **** FOR DEBUGGING **** //
// * [Derek, 02.20.2017] * //
static const Bool_t debug = true;

//------------------------------------------------------------------------------------
void StJetTreeMcMaker::setDbMaker(St_db_Maker* dbMaker)
{
  mDbMaker = dbMaker;
}
//-----------------------------------------------------------------------------------
StJetTreeMcMaker::StJetTreeMcMaker(const char *name,char *dataType):StMaker(name){ 
  

 pTrMatchArr = 0;
 pTriMatchEtaArr =0;
 pTrMatchPhiArr = 0;
 pTrIndexArray = 0;
 pTrEtaArray =0;
 pTrPhiArray =0;

  

//+++++++++++
 muDst = 0;   // NRS: intialization
 muEvent =0;


  outFile="";
}

//_____________________________________________________________________________
StJetTreeMcMaker::~StJetTreeMcMaker(){
  //
}

//_____________________________________________________________________________
// Init - is a first method the top level StChain 
// calls to initialize all its makers 
Int_t StJetTreeMcMaker::Init(){
  // Create Histograms 

 
  mRandom = new TRandom();
  File_output = new TFile(outFile,"RECREATE");

  event = new GammaJetEvent();
  Muevent = new GammaJetEvent();

  outTree_gnt = new TTree("GfmtoDst_gnt","Gfmtodst_gnt");
  outTree_mu = new TTree("GfmtoDst_mu","Gfmtodst_mu");
  outTree_gnt->SetAutoSave(-500000000);  // autosave activated for each 5 MB
  outTree_mu->SetAutoSave(-500000000);  // autosave activated for each 5 MB
  
  eventBranch = outTree_gnt->Branch("GntEventList",&event,1000000);
  mueventBranch = outTree_mu->Branch("MuEventList",&Muevent,1000000);
  
  EvValues = new float[81]; // for event info
  MuEvValues = new float[81]; // for event info

  mMuDstMaker = (StMuDstMaker*)GetMaker("MuDst");

// for E/p calculation [Derek, 07.18.2017] ************************************

  File_output -> cd();

  const Int_t    nEp = 200;
  const Int_t    nPt = 1000;
  const Double_t eP1 = 0.;
  const Double_t eP2 = 20.;
  const Double_t pT1 = 0.;
  const Double_t pT2 = 100.;
  hEbp_elec = new TH1D("hEbp_elec", "E/p, electrons", nEp, eP1, eP2);
  hEbp_hadr = new TH1D("hEbp_hadr", "E/p, hadrons", nEp, eP1, eP2);
  hPt_elec  = new TH1D("hPt_elec", "pT, electrons", nPt, pT1, pT2);
  hPt_hadr  = new TH1D("hPt_hadr", "pT, hadrons", nPt, pT1, pT2);
  hEbp_elec -> Sumw2();
  hEbp_hadr -> Sumw2();
  hPt_elec  -> Sumw2();
  hPt_hadr  -> Sumw2();

// ***************************************************************************/


  // for comparing particle to detector level [Derek, 02.19.208]
  _tMcTracks = new TTree("McTracks", "a tree for MC tracks");
  _tMcTracks -> Branch("EventId", &_McEvtId, "EventId/I");
  _tMcTracks -> Branch("RunId", &_McRunId, "RunId/I");
  _tMcTracks -> Branch("NumTrks", &_McNtrk, "NumTrks/I");
  _tMcTracks -> Branch("MuVtxX", &_MuVx, "MuVtxX/D");
  _tMcTracks -> Branch("MuVtxY", &_MuVy, "MuVtxY/D");
  _tMcTracks -> Branch("MuVtxZ", &_MuVz, "MuVtxZ/D");
  _tMcTracks -> Branch("McVtxX", &_McVx, "McVtxX/D");
  _tMcTracks -> Branch("McVtxY", &_McVy, "McVtxY/D");
  _tMcTracks -> Branch("McVtxZ", &_McVz, "McVtxZ/D");
  _tMcTracks -> Branch("IdTrk", &_McIdTrk);
  _tMcTracks -> Branch("IdGeant", &_McIdGnt);
  _tMcTracks -> Branch("IdVx", &_McIdVx);
  _tMcTracks -> Branch("IdVxEnd", &_McIdVxEnd);
  _tMcTracks -> Branch("IntrVtx", &_McIntrVtx);
  _tMcTracks -> Branch("IsShower", &_McIsShower);
  _tMcTracks -> Branch("Charge", &_McChrg);
  _tMcTracks -> Branch("Rapidity", &_McRap);
  _tMcTracks -> Branch("Eta", &_McEta);
  _tMcTracks -> Branch("Phi", &_McPhi);
  _tMcTracks -> Branch("Px", &_McPx);
  _tMcTracks -> Branch("Py", &_McPy);
  _tMcTracks -> Branch("Pz", &_McPz);
  _tMcTracks -> Branch("Pt", &_McPt);
  _tMcTracks -> Branch("Ptot", &_McPtot);
  _tMcTracks -> Branch("Energy", &_McEne);
  _tMcTracks -> SetAutoSave(-500000000);

  // initialize event leaves
  _McEvtId = -999;
  _McRunId = -999;
  _McNtrk  = -999;
  _MuVx    = -999.;
  _MuVy    = -999.;
  _MuVz    = -999.;
  _McVx    = -999.;
  _McVy    = -999.;
  _McVz    = -999.;
  // initialize track leaves
  _McIdTrk.clear();
  _McIdGnt.clear();
  _McIdVx.clear();
  _McIdVxEnd.clear();
  _McIntrVtx.clear();
  _McIsShower.clear();
  _McChrg.clear();
  _McRap.clear();
  _McEta.clear();
  _McPhi.clear();
  _McPx.clear();
  _McPy.clear();
  _McPz.clear();
  _McPt.clear();
  _McPtot.clear();
  _McEne.clear();

  return StMaker::Init();
}

//_____________________________________________________________________________
/// Make - this method is called in loop for each mStEvent
Int_t StJetTreeMcMaker::Make()
{


  Int_t runIdOfFile;
  doMuEvent(runIdOfFile);
  
  cout<<"\nRun Id OF this event: "<<runIdOfFile<<endl;
  //__________
  doMCEvent(runIdOfFile);
  

  return kStOk;
}
//_______________________________________
Int_t StJetTreeMcMaker::doMuEvent(Int_t & runIdEvnt)
{

  muDst =  (StMuDst*) GetInputDS("MuDst");
  if(!muDst){
    cout << "No MuDst ....." << endl;
    return kStOK;
  }
  
  muEvent = (StMuEvent*) muDst->event();

  //cout<<" I am Here in MuDSt...."<<endl;
  Int_t refmult = muEvent->refMult();

  Int_t RunId   = muEvent->runId();
  Int_t eventId   = muEvent->eventId();
  //cout<<" RunID = "<<RunId<<"  eventId = " <<eventId<<endl;

  double MgF= 0.1*muEvent->runInfo().magneticField();  
  //Int_t nPrimaryVertex = 0;
  //  nPrimaryVertex = muDst->numberOfPrimaryVertices();

  float Vz= 0., Vx=0., Vy=0., Vxy=0.;

  StThreeVectorF Vtx=muEvent->primaryVertexPosition();

  Vz = Vtx.z();
  Vx = Vtx.x();
  Vy = Vtx.y();

  Vxy = TMath::Sqrt( Vx*Vx + Vy*Vy );

  //  vector<PseudoJet> particles;
  //  particles.clear();

  Int_t nPrimary1= 0;
  nPrimary1 = muDst->primaryTracks()->GetEntries();

  // **** FOR DEBUGGING **** //
  // * [Derek, 02.20.2017] * //
  if (debug) {
    cout << " EVENT INFO:\n"
         << "   RunID = " << RunId << ", EvtID = " << eventId << ", Nprim = " << nPrimary1
         << ", (vX, vY, vZ) = (" << Vx << ", " << Vy << ", " << Vz << ")"
         << endl;
  }
  // *********************** //


  //_________________
  Int_t muevent_counter =0;
  Muevent->ResetEvent();
  muevent_counter++; 

  //Initialize the event variales
  for (int g = 0; g < 81; g++){ MuEvValues[g] = -2000; }


  //  Int_t nV=muDst->numberOfPrimaryVertices();

  //______________________
  pTrMatchArr =new Int_t[nPrimary1];
  pTriMatchEtaArr=new Float_t[nPrimary1];
  pTrMatchPhiArr=new Float_t[nPrimary1];
  pTrEtaArray=new Float_t[nPrimary1];
  pTrPhiArray=new Float_t[nPrimary1];
  pTrIndexArray=new Int_t[nPrimary1];


  for(Int_t jj=0; jj<nPrimary1; jj++){
    pTrMatchArr[jj] = -999.;
    pTriMatchEtaArr[jj] =-999.;
    pTrMatchPhiArr[jj] =-999.;
    pTrEtaArray[jj] =-999.;
    pTrPhiArray[jj] =-999.;
    pTrIndexArray[jj]=-999.;
  }


  //_______________________

  
  mGeom_=StEmcGeom::instance("bemc");
  StEmcGeom * mGeom_=StEmcGeom::instance("bemc");
  if(!mGeom_){
      cout << "No StEmcGeom  ....." << endl;
      return kStOK;
    }
  

  GammaJetTrack MuTrackInfo;

  Int_t  mutrkIndex_counter=0.;
  Int_t  mutrk_counter=0.;
  Int_t matchedTrk=0;
  Int_t unmatchedTrk=0;
  

  // **** FOR DEBUGGING **** //
  // * [Derek, 02.20.2017] * //
  if (debug) {
    cout << " TRACK INFO:" << endl;
  }
  // *********************** //

  for ( Int_t i = 0; i < nPrimary1; i++ ){
    
    mutrkIndex_counter++;


    float    p = muDst->primaryTracks(i)->p().mag();
    float    pT = muDst->primaryTracks(i)->pt();
    float    eta = muDst->primaryTracks(i)->eta();
    float    phi = muDst->primaryTracks(i)->phi();  


    TVector3 trk;
    trk.Clear();

    trk.SetPtEtaPhi(pT,eta,phi);
    float       Px = trk.Px();
    float       Py = trk.Py();
    float       Pz = trk.Pz();





    //if( phi < 0. ) phi += 2.*3.14159;
    int flag = muDst->primaryTracks(i)->flag();    
    float    dEdX = muDst->primaryTracks(i)->dEdx();
    float    nHitdedx = muDst->primaryTracks(i)->nHitsDedx();   
    float    nFitpnts = muDst->primaryTracks(i)->nHitsFit();
    float    dca = muDst->primaryTracks(i)->dcaGlobal().mag();
    float    nFitPoss = muDst->primaryTracks(i)->nHitsPoss();
    float    pdca = muDst->primaryTracks(i)->dca().mag();
    float    Charge =  muDst->primaryTracks(i)->charge();    
    float    pidPi = muDst->primaryTracks(i)->nSigmaPion();
    float    pidK =  muDst->primaryTracks(i)->nSigmaKaon();
    float    pidP = muDst->primaryTracks(i)->nSigmaProton();
    float    pidE = muDst->primaryTracks(i)->nSigmaElectron();

    // for comparing particle to detector level [Derek, 02.19.2018]
    Int_t idTruth = muDst -> primaryTracks(i) -> idTruth();
    Int_t qaTruth = muDst -> primaryTracks(i) -> qaTruth();
    Int_t idParVx = muDst -> primaryTracks(i) -> idParentVx();


    // **** FOR DEBUGGING **** //
    // * [Derek, 02.20.2017] * //
    if (debug) {
      cout << "   TrkID = " << i << ", FP = " << nFitpnts
           << ", eta = " << eta << ", pT = " << pT
           << endl;
    }
    // *********************** //

    MuTrackInfo.Clear();
    MuTrackInfo.SetnHitsFit(nFitpnts);
    MuTrackInfo.SetnHitsPoss(nFitPoss);
    MuTrackInfo.SetTrackFlag(0);
    MuTrackInfo.SetPdgId(0);
    MuTrackInfo.SetGeantId(0);
    MuTrackInfo.SetpZ(Pz);
    MuTrackInfo.SetpY(Py);
    MuTrackInfo.SetpX(Px);
    MuTrackInfo.SetpT(pT);
    MuTrackInfo.SetdEdx(dEdX);
    MuTrackInfo.SetCharge(Charge);
    MuTrackInfo.SetTOFBeta(0);
    MuTrackInfo.SetEta(eta);
    MuTrackInfo.SetPhi(phi);
    MuTrackInfo.SetnSigElectron(pidE);
    MuTrackInfo.SetnSigPion(pidPi);
    MuTrackInfo.SetnSigKaon(pidK);
    MuTrackInfo.SetnSigProton(pidP);
    MuTrackInfo.SetDCAg(dca);
    MuTrackInfo.SetnHits(0);
    MuTrackInfo.SetdEdxHits(nHitdedx);
    MuTrackInfo.SetFirstZPoint(0);
    MuTrackInfo.SetLastZPoint(0);
    MuTrackInfo.SetTOFSigElectron((Double_t) idTruth);
    MuTrackInfo.SetTOFSigPion((Double_t) qaTruth);
    MuTrackInfo.SetTOFSigKaon((Double_t) idParVx);
    MuTrackInfo.SetTOFSigProton(0);
    MuTrackInfo.SetPathLength(0);
    MuTrackInfo.SettimeOfflight(0);

    MuTrackInfo.SettrkIndex(mutrkIndex_counter);

    Muevent->AddTrack(&MuTrackInfo,mutrk_counter);
    mutrk_counter++;


    //_________Projecting on


    ptrack = (const StMuTrack*)mMuDstMaker->muDst()->primaryTracks(i);
    unmatchedTrk++;


    if(!(flag > 0))continue;


    StThreeVectorD momentum,position;

    bool proj_ok= mPosition->projTrack(&position,&momentum,ptrack,(Double_t)MgF);

    Float_t z1,eta1,phi1;
    eta1 =position.pseudoRapidity();
    phi1 =position.phi();
    z1   =position.z();
    Int_t m,e,s;

    Int_t okBin=mGeom_->getBin(phi1,eta1,m,e,s);

    
    if( (s < 0 || okBin ==1 ))continue;

    if(proj_ok){
      pTrMatchArr[i] = 1.;
      pTriMatchEtaArr[i] = eta1;
      pTrMatchPhiArr[i] = phi1;
      pTrEtaArray[i] = eta;
      pTrPhiArray[i] = phi;
      pTrIndexArray[i] =mutrkIndex_counter;

// for E/p calculation [Derek, 07.18.2017] ************************************
#if 0
  // Electron: |NSigmaElectron| < 3 && |NSigmaPion| > 3 && |NSigmaKaon| > 3 && |NSigmaProton| > 3
  if (fabs(nSigE) < 3.0 && fabs(nSigP) > 3.0 && fabs(nSigK) > 3.0 && fabs(nSigPi) > 3.0) {
    pTrElectronIds[i] = 1;
  }
  else if (fabs(nSigE) > 3.0 && fabs(nSigP) < 3.0 && fabs(nSigK) < 3.0 && fabs(nSigPi) < 3.0) {
    pTrElectronIds[i] = 0;
  }
  else {
    pTrElectronIds[i]= 9;
  }
  //pTrElectronIds
#endif
// ***************************************************************************/

      matchedTrk++;
      //cout <<".................matched trak: "<<matchedTrk<<" out of   "<<mutrk_counter<<endl;
  }

  }  // end track loop


  // **** FOR DEBUGGING **** //
  // * [Derek, 02.20.2017] * //
  if (debug) {
    cout << " END TRACK INFO" << endl;
  }
  // *********************** //


  // for comparing particle to detector level [Derek, 02.20.2018]
  _McEvtId = -999;
  _McRunId = -999;
  _McNtrk  = -999;
  _MuVx    = -999.;
  _MuVy    = -999.;
  _MuVz    = -999.;
  _McVx    = -999.;
  _McVy    = -999.;
  _McVz    = -999.;
  _McIdTrk.clear();
  _McIdGnt.clear();
  _McIdVx.clear();
  _McIdVxEnd.clear();
  _McIntrVtx.clear();
  _McIsShower.clear();
  _McChrg.clear();
  _McRap.clear();
  _McEta.clear();
  _McPhi.clear();
  _McPx.clear();
  _McPy.clear();
  _McPz.clear();
  _McPt.clear();
  _McPtot.clear();
  _McEne.clear();

  // mc track loop
  const Int_t nMcTrks = (Int_t) muDst -> mcArray(1) -> GetEntries();
  for (Int_t iMcTrk = 0; iMcTrk < nMcTrks; iMcTrk++) {

    // grab pointer
    const StMuMcTrack *mcTrk = (StMuMcTrack*) muDst -> mcArray(1) -> UncheckedAt(iMcTrk);

    // mc track info
    const Int_t    idTrk    = mcTrk -> Id();
    const Int_t    idGnt    = mcTrk -> GePid();
    const Int_t    idVx     = mcTrk -> IdVx();
    const Int_t    idVxEnd  = mcTrk -> IdVxEnd();
    const Int_t    intrVtx  = mcTrk -> ItrmdVertex();
    const Bool_t   isShower = mcTrk -> IsShower();
    const Double_t charge   = mcTrk -> Charge();
    const Double_t rap      = mcTrk -> Rapidity();
    const Double_t eta      = mcTrk -> Eta();
    const Double_t phi      = mcTrk -> Pxyz().phi();
    const Double_t pX       = mcTrk -> Pxyz().x();
    const Double_t pY       = mcTrk -> Pxyz().y();
    const Double_t pZ       = mcTrk -> Pxyz().z();
    const Double_t pT       = mcTrk -> pT();
    const Double_t pTot     = mcTrk -> Ptot();
    const Double_t energy   = mcTrk -> E();

    // fill vectors
    _McIdTrk.push_back(idTrk);
    _McIdGnt.push_back(idGnt);
    _McIdVx.push_back(idVx);
    _McIdVxEnd.push_back(idVxEnd);
    _McIntrVtx.push_back(intrVtx);
    _McIsShower.push_back(isShower);
    _McChrg.push_back(charge);
    _McRap.push_back(rap);
    _McEta.push_back(eta);
    _McPhi.push_back(phi);
    _McPx.push_back(pX);
    _McPy.push_back(pY);
    _McPz.push_back(pZ);
    _McPt.push_back(pT);
    _McPtot.push_back(pTot);
    _McEne.push_back(energy);

  }  // end mc track loop

  // get mc vertex
  const StMuMcVertex *mcVtx = (StMuMcVertex*) muDst -> mcArray(0) -> UncheckedAt(0);
  const Double_t     mcVx   = mcVtx -> XyzV().x();
  const Double_t     mcVy   = mcVtx -> XyzV().y();
  const Double_t     mcVz   = mcVtx -> XyzV().z();

  // mc event info
  _McEvtId = eventId;
  _McRunId = RunId;
  _McNtrk  = nMcTrks;
  _MuVx    = Vx;
  _MuVy    = Vy;
  _MuVz    = Vz;
  _McVx    = mcVx;
  _McVy    = mcVy;
  _McVz    = mcVz;


  //cout<<" attempting to read BEMC ____________"<<endl;
  //_____________________BEMC___________________
  /*
  for (Int_t Iu2 = 1;Iu2 <= 4800;Iu2++) mStatusT[Iu2]=0;
  for (Int_t iu = 0;iu < 4800;iu++) {mCalibPedValues[iu] = 0;mCalibPedRMSValues[iu] = 0;}
  St_db_Maker *dbMk = (St_db_Maker*)GetMaker("db");
  if (dbMk == 0) return kStOK;
  TDataSet *db = NULL;  St_emcStatus* run2=NULL;
  TString DbName = "Calibrations/emc/y3bemc";db = dbMk->GetDataBase(DbName);
  TString tableName2 = "bemcStatus";TString tableName = "bemcPed";//bprs
  run2 = (St_emcStatus*)db->Find(tableName2.Data());emcStatus_st *runst2 = run2->GetTable();
  St_emcPed *ped = (St_emcPed*)db->Find(tableName.Data());        emcPed_st *pedst = ped->GetTable();
  int Iu3(0);    for (Iu3 = 1;Iu3 <=4800;Iu3++) { mStatusT[Iu3-1] = runst2[0].Status[Iu3-1];}
  int iu4(0);  for (iu4 = 0;iu4 <4800;iu4++) { mCalibPedValues[iu4] = ((Float_t)pedst[0].AdcPedestal[iu4])/100.0;}
  for (iu4 = 0;iu4 <4800;iu4++) { mCalibPedRMSValues[iu4] = ((Float_t)pedst[0].AdcPedestalRMS[iu4])/100.0;}
  */


  //  StBemcTables* tables = mTables;
  //  assert(tables);

  //  Int_t flag=0;
  mEmcCol = (StEmcCollection*)muDst->emcCollection();
  StMuEmcCollection* muEmcCol = muDst->muEmcCollection();


  if(!mEmcCol || !muEmcCol)
    {
      printf("\n***-- no EMC Collection was found --***\n");
      return 1;
    }
  StEmcGeom * mGeom0=StEmcGeom::instance("bemc");

  if(!mGeom0)cout<<" mGeom0 NOT FOUND..."<<endl;
  
  StSPtrVecEmcPoint& container =  mEmcCol->barrelPoints();
  

  double Bemc_radius = mGeom0->Radius();
  
  TVector3 VertexVec(Vx,Vy,Vz);
  GammaJetTower TowerInfo;
  GammaJetTowerUtil TwrUtil;
  
  StEmcDetector* detector1 = mEmcCol->detector(kBarrelEmcTowerId);  if (!detector1)   {
    cout<<"There is no kBarrelEmcTowerId Detector"<<endl;
    return kStWarn;}
  
  //cout<<"Inside attempting to read BEMC ____________"<<endl;    

  
  if(detector1) {

    Int_t statusAll  =-99;
    Int_t IdhB = 0;
    int twr_counter =0;
    int matchtwr_counter =0;

    //cout<<" I am inside BEMC detector1"<<endl;


    for (UInt_t mw1A = 1; mw1A < 121; mw1A++){

      StEmcModule* module = detector1->module(mw1A);
      StSPtrVecEmcRawHit& hitsEw1A=module->hits();

      for(unsigned int rw1A=0;rw1A<hitsEw1A.size();rw1A++)
	{

	  Int_t twrmatchIdex = 0; // default 0 not matched
	  Int_t  SubhB=hitsEw1A[rw1A]->sub();
	  UInt_t   ModulehB=hitsEw1A[rw1A]->module();
	  UInt_t EtahB=hitsEw1A[rw1A]->eta();
	  
	//	mm = (Int_t)emcTowerHits[j]->module();
	//	ee = (Int_t)emcTowerHits[j]->eta();
	//	ss = emcTowerHits[j]->sub();
	  
	if(abs(ModulehB)<=120&&abs(EtahB)<=20&&SubhB<=2)     {

          //cout << "   ---> (3f) Tower Loop" << endl;

          // Looks like 'mGeom' isn't initialized anywhere...
          // Perhaps 'mGeom0' is what was supposed to be here?
          // [Derek, 02.17.2017]

	  //mGeom->getId(ModulehB, EtahB, SubhB, IdhB);
	  mGeom0->getId(ModulehB, EtahB, SubhB, IdhB);

	  //	  tables->getStatus(BTOW, IdhB, statusAll);
	  Int_t statusAll  =-99;
	  //for BTOW  =1
	  Float_t ped, rms;
	  //	  tables->getStatus(1, IdhB, statusAll);
	  //                              tables->getStatus(1, IdhB, statusAll);
	  //                              tables->getPedestal(1,IdhB,0,ped,rms);
	  //                              cout<<" statusAll = "<<statusAll<<" ped= "<<ped<<"  rms=  "<<rms<<endl;
	  // statusAll = mStatusT[IdhB-1];
	  //	  Float_t ped =mCalibPedValues[IdhB-1];
	  //	  Float_t rms =mCalibPedRMSValues[IdhB-1];
	  //                              cout<<" hstatusAll = "<<statusAll <<" hped= "<<ped<<"  hrms=  "<<rms<<endl;
	  
	  //	  if(!(statusAll==1))continue;

	  
	  Float_t    EnergyhB=hitsEw1A[rw1A]->energy();

	  if(!(EnergyhB > 0.2))continue; // constraining lowest twr energy

	  Float_t AdchB=hitsEw1A[rw1A]->adc();
	  Float_t eThB(-99), phiThB(-99), ThetahB(-99.);
	  mGeom0->getEtaPhi(IdhB,eThB,phiThB);
	  mGeom0->getTheta(ModulehB,EtahB,ThetahB);


	  //cout<<" statusAll ="<< statusAll<<"  eta= "<<eThB<<" phi= "<<phiThB<<endl;
	  
	  float dPhi_TwrTrk_[10], dEta_TwrTrk_[10], mTrk_P_[10], m_index[10];
	  for(int i=0; i< 10; i++)
	    {
	      dPhi_TwrTrk_[i]= -99.; dEta_TwrTrk_[i] = -99.; mTrk_P_[i] =-99.; m_index[i]=0;
	      
	    }
	  Float_t dPhi_TwrTrk = -999.;
	  Float_t dEta_TwrTrk = -999.;
	  Float_t matchedTrk_p = -999.;
	  Int_t  matchedTrk_index = -9;
	  Int_t  no_matchTrk_counter =0;
          // for E/p calculation [Derek, 07.20.2017]
          Float_t matched_nSigPi = -999.;
          Float_t matched_nSigK  = -999.;
          Float_t matched_nSigP  = -999.;
          Float_t matched_nSigE  = -999.;
          Float_t matched_dcag   = -999.;
          Float_t matched_eta    = -999.;
          Float_t matched_pT     = -999.;
          Int_t   matched_nFit   = -999;
          Int_t   matched_nPos   = -999;

	  vector<Int_t> vec_matchTrkIndex; //to store ptrks index
	  vec_matchTrkIndex.clear();

	  vector<Float_t> vec_matchTrkP; // to store ptrks p
	  vec_matchTrkP.clear();

          // for E/p calculation [Derek, 07.20.2017]
          vector<Float_t> vec_matchTrkNsigPi;
          vector<Float_t> vec_matchTrkNsigK;
          vector<Float_t> vec_matchTrkNsigP;
          vector<Float_t> vec_matchTrkNsigE;
          vector<Float_t> vec_matchTrkDcag;
          vector<Float_t> vec_matchTrkEta;
          vector<Float_t> vec_matchTrkPt;
          vector<Int_t>   vec_matchTrkNfit;
          vector<Int_t>   vec_matchTrkNpos;
          vec_matchTrkNsigPi.clear();
          vec_matchTrkNsigK.clear();
          vec_matchTrkNsigP.clear();
          vec_matchTrkNsigE.clear();
          vec_matchTrkDcag.clear();
          vec_matchTrkEta.clear();
          vec_matchTrkPt.clear();
          vec_matchTrkNfit.clear();
          vec_matchTrkNpos.clear();

	  Float_t total_matchedTrk_P = 0.;

	 
	  //__________Track pointing to each track
	  //                             cout<<"In side Tower loop #primary  "<<nPrimary1<<endl;
	  for ( Int_t ii = 0; ii < nPrimary1; ii++ )
	    { ////ptrackloop

	      if(pTrMatchArr[ii]==1){ ///


		Float_t mEta, mPhi, tEta, tPhi;
		Int_t mTrkIndex;
		mEta = pTriMatchEtaArr[ii];
		mPhi = pTrMatchPhiArr[ii];
		tEta = pTrEtaArray[ii];
		tPhi = pTrPhiArray[ii];
		//mTrkIndex = pTrIndexArray[ii];

		//____Check dphi and deta____
		dPhi_TwrTrk = mPhi - phiThB;
		dEta_TwrTrk = mEta - eThB;

		if (fabs(dPhi_TwrTrk) > 0.025) continue;
		if (fabs(dEta_TwrTrk) > 0.025) continue;

		// chnage twrmatchIdex index
		twrmatchIdex = 1; // set 1
		
		//__________To check unmatched track diff___
		Float_t trketa_t = mMuDstMaker->muDst()->primaryTracks(ii)->eta();
		Float_t trkphi_t = mMuDstMaker->muDst()->primaryTracks(ii)->phi();
		Float_t dEta_test = trketa_t - tEta;
		Float_t dPhi_test = trkphi_t - tPhi;

		//_____________Momentum of matched trak
		matchedTrk_p = mMuDstMaker->muDst()->primaryTracks(ii)->p().mag();

                // for E/p calculation [Derek, 07.20.2017]
                matched_nSigPi = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaPion();
                matched_nSigK  = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaKaon();
                matched_nSigP  = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaProton();
                matched_nSigE  = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaElectron();
                matched_dcag   = mMuDstMaker -> muDst() -> primaryTracks(ii) -> dcaGlobal().mag();
                matched_eta    = mMuDstMaker -> muDst() -> primaryTracks(ii) -> eta();
                matched_pT     = mMuDstMaker -> muDst() -> primaryTracks(ii) -> pt();
                matched_nFit   = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nHitsFit();
                matched_nPos   = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nHitsPoss();


// for E/p calculation [Derek, 07.18.2017] ************************************
                float  nSigPi_ = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaPion();
                float  nSigK_  = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaKaon();
                float  nSigP_  = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaProton();
                float  nSigE_  = mMuDstMaker -> muDst() -> primaryTracks(ii) -> nSigmaElectron();

//#if 0

  const Float_t matched_rFit = (Float_t) matched_nFit / (Float_t) matched_nPos;
  if (matched_nFit > 15) {
    if (matched_rFit > 0.52) {
      if (matched_dcag < 1.0) {
        if (fabs(matched_eta) < 1.0) {
          if (matched_pT > 0.2) {

            if (fabs(nSigE_) < 3.0 && fabs(nSigP_) > 3.0 && fabs(nSigK_) > 3.0 && fabs(nSigPi_) > 3.0) {
              hEbp_elec -> Fill(EnergyhB / matchedTrk_p);
              hPt_elec  -> Fill(matched_pT);
            }
            else if (fabs(nSigE_) > 3.0 && fabs(nSigP_) < 3.0 && fabs(nSigK_) < 3.0 && fabs(nSigPi_) < 3.0) {
              hEbp_hadr -> Fill(EnergyhB / matchedTrk_p);
              hPt_hadr  -> Fill(matched_pT);
            }
            else {
              //cout << " Not electron/hadron " << endl;
            }  // PID cuts

          }  // pT cut
        }  // eta cut
      }  // dca cut
    }  // nFit/nPoss cut
  }  // nFit cut

//#endif
// ***************************************************************************/

		matchedTrk_index = ii;
		
		if(twrmatchIdex == 1)
		  {
		    total_matchedTrk_P += matchedTrk_p;   //total momentum of all pointed traks
		    no_matchTrk_counter++;
		    //    vec_matchTrkIndex.push_back(mTrkIndex);
		    vec_matchTrkIndex.push_back(matchedTrk_index);
		    vec_matchTrkP.push_back(matchedTrk_p);
                    // for E/p calculation
                    vec_matchTrkNsigPi.push_back(matched_nSigPi);
                    vec_matchTrkNsigK.push_back(matched_nSigK);
                    vec_matchTrkNsigP.push_back(matched_nSigP);
                    vec_matchTrkNsigE.push_back(matched_nSigE);
                    vec_matchTrkDcag.push_back(matched_dcag);
                    vec_matchTrkEta.push_back(matched_eta);
                    vec_matchTrkPt.push_back(matched_pT);
                    vec_matchTrkNfit.push_back(matched_nFit);
                    vec_matchTrkNpos.push_back(matched_nPos);
		  }
		else {
		  vec_matchTrkIndex.push_back(-99);
		  vec_matchTrkP.push_back(-99.);
                  // for E/p calculation
                  vec_matchTrkNsigPi.push_back(-999.);
                  vec_matchTrkNsigK.push_back(-999.);
                  vec_matchTrkNsigP.push_back(-999.);
                  vec_matchTrkNsigE.push_back(-999.);
                  vec_matchTrkDcag.push_back(-99.);
                  vec_matchTrkEta.push_back(-99.);
                  vec_matchTrkPt.push_back(-99.);
                  vec_matchTrkNfit.push_back(-99);
                  vec_matchTrkNpos.push_back(-99);
		}
		
		

	      } //for if(pTrMatchArr[ii]==1)

	    }//for for ( Int_t ii = 0; ii < nPrimary1; ii++ )

	  TowerInfo.Clear();

	  Float_t sum_p =0.;
	  if(twrmatchIdex==1){
	    //                              cout<<" in Tower id= "<<rw1A<<" module# "<<mw1A<<"  matchedIndx= "<<twrmatchIdex<<endl;
	    matchtwr_counter++;

	    for(int i=0; i < vec_matchTrkIndex.size();i++)
	      {
		TowerInfo.SetMatchedTracksArray(vec_matchTrkIndex[i],i);
		TowerInfo.SetMatchedTracksArray_P(vec_matchTrkP[i],i);
		sum_p += vec_matchTrkP[i];
                // for E/p calculation [Derek, 07.20.2017]
                TowerInfo.SetMatchedTracksArray_nSigPi(vec_matchTrkNsigPi[i], i);
                TowerInfo.SetMatchedTracksArray_nSigK(vec_matchTrkNsigK[i], i);
                TowerInfo.SetMatchedTracksArray_nSigP(vec_matchTrkNsigP[i], i);
                TowerInfo.SetMatchedTracksArray_nSigE(vec_matchTrkNsigE[i], i);
                TowerInfo.SetMatchedTracksArray_dcag(vec_matchTrkDcag[i], i);
                TowerInfo.SetMatchedTracksArray_eta(vec_matchTrkEta[i], i);
                TowerInfo.SetMatchedTracksArray_pT(vec_matchTrkPt[i], i);
                TowerInfo.SetMatchedTracksArray_nFit(vec_matchTrkNfit[i], i);
                TowerInfo.SetMatchedTracksArray_nPos(vec_matchTrkNpos[i], i);
	      } //

	  }// if(twrmatchIdex==1


	  TowerInfo.SetTwrId(IdhB);
	  TowerInfo.SetTwrEng(EnergyhB);
	  TowerInfo.SetTwrPhi(phiThB);
	  TowerInfo.SetTwrEta(eThB);
	  TowerInfo.SetTwrADC(AdchB);
	  TowerInfo.SetTwrPed(0);
	  TowerInfo.SetTwrRMS(0);
	  TowerInfo.SetTwrMatchIdnex(twrmatchIdex);
	  TowerInfo.SetNoOfmatchedTrk(no_matchTrk_counter);
	  TowerInfo.SetTwrMatchSumP(sum_p);
	  //GetMomentumVectorForTower(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng, const TVector3& Vrtx)
	  TLorentzVector TwrMomVec = TwrUtil.GetMomentumVectorForTower( Bemc_radius, eThB, phiThB, EnergyhB, VertexVec); ///Find me here

	  float Twr_px = TwrMomVec.Px();
	  float Twr_py = TwrMomVec.Py();
	  float Twr_pz = TwrMomVec.Pz();

	  TowerInfo.SetTwrPx(Twr_px);
	  TowerInfo.SetTwrPy(Twr_py);
	  TowerInfo.SetTwrPz(Twr_pz);

	  Muevent->AddTower(&TowerInfo,twr_counter);
	  twr_counter++;
	 

	}// for if(abs(mm)<=120&&abs(ee)<=20&&ss<=2)

      } //for (UInt_t j = 0; j < emcTowerHits.size(); j++)

      
    } //for (UInt_t i = 1; i < 121; i++)


  } //if(detector1)






  ///____________________________________________________________________________
  //cout<<"MuEvent:  Refmult= "<<refmult<<"  RunId= "<<RunId<<" evntId= "<<eventId<<" V x y z: "<<Vx<<" "<<Vy<<"  "<<Vz<<endl;

  MuEvValues[0]   = RunId;
  MuEvValues[1]   = eventId;
  MuEvValues[2]   = 0;  //Trigger Id check getTrgId()
  MuEvValues[3]   = 0;
  MuEvValues[4]   = mutrkIndex_counter;
  MuEvValues[5]   = refmult;
  MuEvValues[6]   = 0;
  MuEvValues[7]   = Vx;
  MuEvValues[8]   = Vy;
  MuEvValues[9]   = Vz;
  MuEvValues[10]   = 0;
  MuEvValues[11]   = 0;      //zdcConincidenceRate
  MuEvValues[12]   = 0;   //bbcCoincidenceRate
  MuEvValues[13]   = 0;  //backgroundRate
  MuEvValues[14]   = 0;  //bbcBlueBackgroundRate
  MuEvValues[15]   = 0; //bbcYellowBackgroundRate
  MuEvValues[16]   =0;
  MuEvValues[17]   =0;
  MuEvValues[18]   =0 ; //bTOFTrayMultiplicity (need to add)
  MuEvValues[19]   =0;
  MuEvValues[20]   = 0;
  MuEvValues[80]   = 0;

  Muevent->SetEventAttributes(MuEvValues);

  outTree_mu->Fill();

  // for comparing particle to detector tracks [Derek, 02.20.2018]
  _tMcTracks -> Fill();

  
  runIdEvnt =0; // initialize with a non existing runid
  runIdEvnt = RunId;


  return kStOk;
}

///____________________________________________
Int_t StJetTreeMcMaker::doMCEvent(Int_t RunIdOfEvent)
{  

  mcEvent =  (StMcEvent*) GetDataSet("StMcEvent");       

  Int_t event_counter =0;
  
  event->ResetEvent();
  event_counter++;

  GammaJetTrack TrackInfo;  
  //___________________
  //Initialize the event variales
  for (int g = 0; g < 81; g++){ EvValues[g] = -2000; }

  ///______________________________________
  const StPtrVecMcTrack& mcTracks = mcEvent->primaryVertex()->daughters();
  StMcTrackConstIterator mcTrkIter = mcTracks.begin();
  Int_t track_counter = 0;
  Int_t trkIndex_counter=0;

  for ( ; mcTrkIter != mcTracks.end(); ++mcTrkIter) {
    StMcTrack* track = *mcTrkIter;

    trkIndex_counter++;

    Double_t mcTrk_px = track->momentum().x();
    Double_t mcTrk_py = track->momentum().y();
    Double_t mcTrk_pz = track->momentum().z();
    Double_t mcTrk_pt = sqrt(mcTrk_px*mcTrk_px + mcTrk_py*mcTrk_py);
    Double_t mcTrk_eng = track->energy();
    Double_t mcTrk_eta = track->pseudoRapidity();
    Int_t    mcTrk_pdgId = track->pdgId();
    Int_t    mcTrk_geantId = track->geantId();

    Int_t mcTrk_chrg = -9999;
    Int_t mcTrk_dEdX = -9999;
    if( track->particleDefinition()){
      mcTrk_chrg = track->particleDefinition()->charge();
      mcTrk_dEdX = track->particleDefinition()->pdgEncoding();
    }
    else{
      cout << "Particle with no encoding " << endl;
      mcTrk_chrg = 0.;
      mcTrk_dEdX = 0.;

    }
    Int_t mcTrk_NoOfFitHits = track->tpcHits().size();
    
    TrackInfo.Clear();
    TrackInfo.SetnHitsFit(0);
    TrackInfo.SetnHitsPoss(0);
    TrackInfo.SetTrackFlag(0);
    TrackInfo.SetPdgId(mcTrk_pdgId);
    TrackInfo.SetGeantId(mcTrk_geantId);
    TrackInfo.SetpZ(mcTrk_pz);
    TrackInfo.SetpY(mcTrk_py);
    TrackInfo.SetpX(mcTrk_px);
    TrackInfo.SetpT(mcTrk_pt);
    TrackInfo.SetdEdx(mcTrk_dEdX);
    TrackInfo.SetCharge(mcTrk_chrg);
    TrackInfo.SetTOFBeta(0);
    TrackInfo.SetEta(mcTrk_eta);
    TrackInfo.SetPhi(0);
    TrackInfo.SetnSigElectron(0);
    TrackInfo.SetnSigPion(0);
    TrackInfo.SetnSigKaon(0);
    TrackInfo.SetnSigProton(0);
    TrackInfo.SetDCAg(0);
    TrackInfo.SetnHits(0);
    TrackInfo.SetdEdxHits(0);
    TrackInfo.SetFirstZPoint(0);
    TrackInfo.SetLastZPoint(0);
    TrackInfo.SetTOFSigElectron(0);
    TrackInfo.SetTOFSigPion(0);
    TrackInfo.SetTOFSigKaon(0);
    TrackInfo.SetTOFSigProton(0);
    TrackInfo.SetPathLength(0);
    TrackInfo.SettimeOfflight(0);
    TrackInfo.SettrkIndex(trkIndex_counter);

    event->AddTrack(&TrackInfo,track_counter);
    track_counter++;
        cout<<mcTrk_pz<<"  "<<mcTrk_px<<"   "<<mcTrk_py<<"   "<<mcTrk_chrg<<"   "<<mcTrk_dEdX<<"   "<<mcTrk_NoOfFitHits<<endl;
    
    
  }


  Int_t mcEvntNumber = mcEvent->eventNumber();
  Int_t mcRunNumber = mcEvent->runNumber();
  Int_t mcEvntRefmult = mcEvent->eventGeneratorFinalStateTracks();
  //  Int_t mcEvntNtrkCount = nCount;
  Int_t mcEvntNPrimTrk = mcEvent->numberOfPrimaryTracks();
  Double_t mcEvntReactPlane =  mcEvent->phiReactionPlane();
  Double_t mcEvntVrtX = mcEvent->primaryVertex()->position().x();
  Double_t mcEvntVrtY = mcEvent->primaryVertex()->position().y();
  Double_t mcEvntVrtZ = mcEvent->primaryVertex()->position().z();
  
  
  cout<<"Event Info "<<mcEvntNumber<<"  "<<mcRunNumber<<"   "<<mcEvntRefmult<<"  "<<mcEvntNPrimTrk<<"   "<<mcEvntReactPlane<<"  "<<mcEvntVrtX<<"  "<<mcEvntVrtY<<"  "<<mcEvntVrtZ<<endl;


  //  EvValues[0]   = mcRunNumber;
  EvValues[0]   = RunIdOfEvent;
  EvValues[1]   = mcEvntNumber;
  EvValues[2]   = 0;  //Trigger Id check getTrgId()
  EvValues[3]   = 0;
  EvValues[4]   = mcEvntNPrimTrk;
  EvValues[5]   = mcEvntRefmult;
  EvValues[6]   = 0;
  EvValues[7]   = mcEvntVrtX;
  EvValues[8]   = mcEvntVrtY;
  EvValues[9]   = mcEvntVrtZ;
  EvValues[10]   = 0;
  EvValues[11]   = 0;      //zdcConincidenceRate
  EvValues[12]   = 0;   //bbcCoincidenceRate
  EvValues[13]   = 0;  //backgroundRate
  EvValues[14]   = 0;  //bbcBlueBackgroundRate
  EvValues[15]   = 0; //bbcYellowBackgroundRate
  EvValues[16]   =0;
  EvValues[17]   =0;
  EvValues[18]   =0 ; //bTOFTrayMultiplicity (need to add)
  EvValues[19]   =0;
  EvValues[20]   = 0;
  EvValues[80]   = 0;

  event->SetEventAttributes(EvValues);

  outTree_gnt->Fill(); 

  


return kStOk;
 
}
//-----------------------------------------------------------------------------------------------------
Int_t StJetTreeMcMaker::Finish()
{
 //write histogram
  File_output = outTree_gnt->GetCurrentFile();
  File_output = outTree_mu->GetCurrentFile();
  File_output->Write();
  File_output->Close();
  return kStOk;
}

//-----------------------------------------------------------------------------------------------------void
void  StJetTreeMcMaker::Clear(Option_t *opt)
{
  delete [] pTrMatchArr; pTrMatchArr = 0;
  delete []  pTriMatchEtaArr; pTriMatchEtaArr =0;
  delete []  pTrMatchPhiArr; pTrMatchPhiArr = 0;
  delete []  pTrIndexArray; pTrIndexArray = 0;
  delete []  pTrEtaArray; pTrEtaArray =0;
  delete []  pTrPhiArray; pTrPhiArray =0;


  StMaker::Clear();

  // return kStOk;
 }
 
