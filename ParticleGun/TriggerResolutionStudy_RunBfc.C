// 'TriggerResolutionStudy_RunBfc.C'
// Derek Anderson
// 08.30.2020
//
// Runs BFC.  Adapted from code from
// Robert Licenik and Maria Zurek.



void TriggerResolutionStudy_RunBfc(int nEvents = 10, const char* cFzdFile = "kinematics_local.starsim.fzd") {

  // chain options
  const char* cChainOpt = "fzin,gen_T,geomT,sim_T,TpcRS,l0,ry2009d,DbV20120908,pp2009d ITTF BEmcChkStat btof Corr4 OSpaceZ2 OGridLeak3D VFMCE TpxClu -hitfilt TpxClu -VFMinuit VFPPVnoCTB beamLine -hitfilt,TpcMixer,GeantOut,MiniMcMk,McAna,IdTruth,-in,useInTracker,-emcDY2,emcSim sdt20120908";

  // load macro
  gROOT -> LoadMacro("bfc.C");
  if (gClassTable -> GetID("StBFChain") < 0) Load(); 
  bfc(-1, cChainOpt, cFzdFile);

  // run chain
  int iStat = chain - >Init();
  if (iStat) {
    cout << "Chain initialization failed" << endl;
    chain -> Fatal(iStat, "during Init()");
  }

  // detailed listing of makers in the chain
  cout << "Order of makers in BFC:" << endl;
  StMaker::lsMakers(chain);

  // clear chain
  for(Int_t iEvent = 1; iEvent <= nEvents; iEvent++) {
    chain -> Clear();
    Int_t iMake = chain -> Make(iEvent);
    if(iMake % 10 == kStEOF || iMake % 10 == kStFatal) break;
  }
}

// End ------------------------------------------------------------------------
