// 'bemcCalibMacro.C'
// Nihar Sahoo
//
// Use this macro first to create trees from '*.MuDst.root' and
// '*.geant.root'
//
// Usage Note: This has been modified.  The macro assumes that
//   the provided MuDst file is in the same directory as the
//   provided geant file.  Provide the FULL NAME AND PATH for
//   the geant file, but provide ONLY THE NAME of the corresponding
//   MuDst file.  [Derek, 01.03.2018]

// [Derek, 01.03.2018]
#include "TString.h"

/*
void bemcCalibMacro(int nEvents = 1000000,
		    const char* filelist = "pt9_11_10151001_1.MuDst.root",
		    const char* outFile1  = "matchingTest.root",
		    int nFiles = 1,
		    const char* geantfile = "/projecta/projectdirs/starprod/embedding/production2009_200GeV/Jet_pp200_2009.elz14/SL11d_embed/10151001/pt9_11_10151001_1.geant.root")
{
*/
void bemcCalibMacro(int nEvents = 1000000,
		    const char* filelist  = "st_physics_adc_15076102_raw_0000000.MuDst.root",
		    const char* outFile1  = "auau200r14pt911.wave0dir076run15076102_test.root",
		    int nFiles = 1,
		    const char* geantfile = "/star/data105/embedding/AuAu_200_production_2014/Pythia6_pt9_11_100_20192901/P18ih.SL18h/2014/076/15076102/st_physics_adc_15076102_raw_0000000.MuDst.root")
{

  gROOT->Macro("LoadLogger.C");
  gROOT->Macro("loadMuDst.C");
  gSystem->Load("StarMagField.so");
  gSystem->Load("StMagF");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StTriggerUtilities");
  gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StEmcTriggerMaker");

  //gROOT->Macro("loadMuDst.C");
  //  gSystem->Load("/star/u/nihar/Macro_test/Tree_Stuc_v2/.sl64_gcc482/lib/libStJetTreeThirdMaker.so");
  //  gSystem->Load("/star/u/nihar/JetRecon_pp/FullJetTree_pp/.sl64_gcc482/lib/libStJetTreeThirdMaker.so");

  gSystem->Load("/star/data01/pwg/dmawxc/Embedding/Run9pp/NiharMuDstMatcher/.sl73_gcc485/lib/libStJetTreeMcMaker.so");
  //gSystem->Load("/star/u/nihar/EmbeddingResponxMatrx/MuDstMatcher/.sl73_gcc485/lib/libStJetTreeMcMaker.so");

  //gSystem->Load( "StJetTreeMcMaker" );

	
  StChain* chain = new StChain("StChain");

  // I/O maker
  StIOMaker* ioMaker = new StIOMaker; 
  ioMaker->SetFile(geantfile);
  ioMaker->SetIOMode("r");
  ioMaker->SetBranch("*",0,"0");             // Deactivate all branches
  ioMaker->SetBranch("geantBranch",0,"r");   // Activate geant Branch
    
  // StMcEvent maker
  StMcEventMaker* mcEventMaker = new StMcEventMaker;
  mcEventMaker->doPrintEventInfo = false;
  mcEventMaker->doPrintMemoryInfo = false;


  // determine mudst path [Derek, 01.03.2018]
  TString geant(geantfile);
  TString mudst(filelist);
  TString uPath("");
  Ssiz_t  nPath = geant.Last('/');
  Ssiz_t  nName = nPath + 1;
  uPath.Append(geant.Data(), nName);
  uPath.Append(mudst.Data());
  cout << "I/O Info: geant = '" << geant.Data() << "'\n"
       << "          mudst = '" << uPath.Data() << "'"
       << endl;


  // MuDst maker
  StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",uPath.Data(),"",nFiles);

//* I dunno why, but this is causing me errors... [10.26.2016, Derek]
  // star database
  St_db_Maker *dbMaker = new St_db_Maker("StarDb","MySQL:StarDb");
  // Endcap database
  StEEmcDbMaker* eemcb = new StEEmcDbMaker("eemcDb");
  // Barrel ADC to energy maker
  StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
  //  StPreEclMaker *pre_ecl  = new StPreEclMaker();
  //  StEpcMaker *epc         = new StEpcMaker();
//*/

//*
  //get control table so we can turn off BPRS zero-suppression and save hits from "bad" caps
  controlADCtoE_st* control_table = adc->getControlTable();
  control_table->CutOff[1] = -1;
  control_table->CutOffType[1] = 0;
  control_table->DeductPedestal[1] = 2;
  adc->saveAllStEvent(kTRUE);
 
/* might not need... 
  StTriggerSimuMaker* trigsim = new StTriggerSimuMaker();
  trigsim->setMC(2);	// 0=data, 1=simulation, 2=embedding
  trigsim->useBemc();
  trigsim->useEemc();
  //trigsim->useBbc();							//2016.08.04 It turns out bbc was not implemented in simumaker for run12. When I turn this option on, it will just mess thing up and this is why i didn't get right JP2 trigger should fire information! Unfortunately in the data production, I also used this, not sure how would it play out there. but the data on rcf are removed from disk. I cannot go back to it now. It will also affect Au+Au run11 production as well as p+p run12 production
  
  trigsim->useOfflineDB();						// only offlineDB is available on pdsf date 2016.07
  //2016.07.28	let's test difference between off and on line
  //trigsim->bemc->setConfig(StBemcTriggerSimu::kOffline);		// I used Offline for real data. to be consistent  // 2017.11.29 let's test offline
  trigsim->bemc->setConfig(StBemcTriggerSimu::kOnline);		//2016.07.28 let's test Online
*/
  StJetTreeMcMaker* tMcMaker = new StJetTreeMcMaker("tMcMaker");
  tMcMaker->SetFileName(outFile1);


  StMemStat memory;
  memory.PrintMem(NULL);
	
  chain->Init();
  cout<<"chain initialized"<<endl;
	
  TStopwatch total;
  TStopwatch timer;
	
  int i=0;
  while(i<nEvents && chain->Make()==kStOk)
    {
      if(i % 100000 == 0){
	cout<<"done with event "<<i;
	cout<<"\tcpu: "<<timer.CpuTime()<<"\treal: "<<timer.RealTime()<<"\tratio: "<<timer.CpuTime()/timer.RealTime();//<<endl;
	timer.Start();
	memory.PrintMem(NULL);
      }
      i++;
      chain->Clear();
    }
	
  chain->ls(3);
  chain->Finish();
  //  printf("my macro processed %i events in %s",i,nametag);
  cout<<"\tcpu: "<<total.CpuTime()<<"\treal: "<<total.RealTime()<<"\tratio: "<<total.CpuTime()/total.RealTime()<<endl;

  cout << endl;
  cout << "-------------" << endl;
  cout << "(-: Done :-) " << endl;
  cout << "-------------" << endl;
  cout << endl;
}
