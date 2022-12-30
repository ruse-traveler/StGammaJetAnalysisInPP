//_________________
class StJetTreeThirdMaker;

//void doFullJetTree(unsigned long nEvents=2000000000,char* outFile1,char* fName = "*") //100
void doFullJetTree(unsigned long nEvents=2000000000,char* outFile1.root,char* fName = "/star/embed/embedding/JetEmbedding2009/Jet_pp200_2009.elz17/SL11d_embed/10125076/pt9_11_10125076_1.MuDst.root") //100
{ 
  
  char* flag = "";
  bool printOutput = kFALSE;  
  TString inputFile(fName);
  // 1. load shared libraries:
  if (gClassTable->GetID("TTable") < 0)
     {
       gSystem->Load("libStar");
       gSystem->Load("libPhysics");
       
     }
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  //gSystem->Load("StMagF");
  gSystem->Load("StUtilities");  // new addition 22jul99
  gSystem->Load("StTreeMaker");
  gSystem->Load("StIOMaker");
  gSystem->Load("StarClassLibrary");
  // gSystem->Load("StTpcDb");
  // gSystem->Load("StDbUtilities");
  gSystem->Load("StEvent");
  gSystem->Load("StEventUtilities");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StAssociationMaker");
  //gSystem->Load("StMcAnalysisMaker");
  gSystem->Load("StStrangeMuDstMaker");
  gSystem->Load("StEmcUtil");
  gSystem->Load("StMuDSTMaker");
  //gSystem->Load("StBichsel");
  // libraries for the DB
  gSystem->Load("StDaqLib");
  gSystem->Load("StDbLib");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StEpcMaker");

  // 2. load library of your own maker:
  gROOT->Macro("loadMuDst.C");
  //  gSystem->Load("/star/u/nihar/Macro_test/Tree_Stuc_v2/.sl64_gcc482/lib/libStJetTreeThirdMaker.so");
  gSystem->Load("/star/u/dmawxc/JetRecon_pp/FullJetTree/ThirdJetMaker/.sl64_gcc482/lib/libStJetTreeThirdMaker.so");
  // gSystem->Load("/star/u/saskiam/nihar/GammaJet_Test_pp/ThirdJetMaker/.sl64_gcc482/lib/libStJetTreeThirdMaker.so");
  
  gSystem->Load( "StJetTreeThirdMaker" );
   
  // 3. Print some statement
  cout << "Finished loading libraries " << endl;

  // 4. Instantiate an StChain object
  StChain* chain = new StChain("chain");
  chain->SetDebug(1);
  
  // 5. Define an input file
  TString inputDir ="";

  
  if((inputDir!="") && (!inputDir.EndsWith("/")))
    inputDir += "/";
  cout << "inputFile.Data() = " << inputFile.Data() << endl;
  cout << "inputDir.Data() = " << inputDir.Data() << endl;
  
   // StPreEclMaker* preEcl = 0;
   // St_db_Maker* dbMk = 0;
  StDetectorDbMaker* detDbMk=0;
  
   if((inputFile.EndsWith(".list") && strcmp(flag, "")==0) ||
      inputFile.EndsWith(".MuDst.root") || inputFile == "")
     {
       cout << "Processing " << inputFile.Data() << endl;
       StMuDstMaker * muDstMaker =
      	new StMuDstMaker(0, 0, inputDir.Data(), inputFile.Data(), "MuDst.root", 2000, "MuDst");
       StMuDst2StEventMaker* eventMaker = new StMuDst2StEventMaker("eventMaker");
       // 2nd argument is the name of the StMuDstMaker
       St_db_Maker *dbMK=new St_db_Maker("db","MySQL:StarDb");
       //dbMK->SetDateTime(20070101,000001);
       //dbMK->SetFlavor("sim","bprsCalib");
       dbMK->SetFlavor("sim","bsmdeCalib");
       dbMK->SetFlavor("sim","bsmdpCalib");
       //if(!dbMk)   dbMk = new St_db_Maker("StarDb", "MySQL:StarDb","$STAR/StarDb"); 
       // first parameter is the selfname, then go DB directories
       // the order of directories is important:
       // the maker will search tables in this order and the first found will be used
       if(strcmp(flag, "")==0)
	 {
	   //Database interface
	  if(!detDbMk) detDbMk = new StDetectorDbMaker();   
	 }   
     } 
   StIOMaker* ioMaker = 0;
   
    if(inputFile.EndsWith(".list"))
      inputDir.Prepend("@");
    if(inputFile.EndsWith("event.root") || inputDir.BeginsWith("@"))
      {
	StIOMaker* ioMaker = 
	  new StIOMaker("IO", "r", inputDir.Append(inputFile).Data());
	ioMaker->SetIOMode("r");
	ioMaker->SetBranch("*", 0, "0"); //deactivate all branches
	ioMaker->SetBranch("eventBranch", 0, "r");
	St_db_Maker *dbMK=new St_db_Maker("db","MySQL:StarDb");
	// dbMK->SetDateTime(20070101,000001);
	// dbMK->SetFlavor("sim","bsmdeCalib");
	// dbMK->SetFlavor("sim","bsmdpCalib");
	//if(!dbMk)  dbMK=new St_db_Maker("db","MySQL:StarDb");
	if(strcmp(flag, "")==0)
	  {
	    //Database interface
	    if(!detDbMk) detDbMk = new StDetectorDbMaker();   
	  } 
      }
    
    StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
    adc->setPrint(printOutput);
    
    preEcl = new StPreEclMaker();preEcl->setPrint(printOutput);
    epc = new StEpcMaker();
 
    /*controlADCtoE_st* control_table=adc->getControlTable();
     control_table->CutOff[1]=-1;
     control_table->CutOffType[1]=0;*/
    
    // 7. Add your own maker to chain
    StJetTreeThirdMaker* thirdMaker = new StJetTreeThirdMaker("thirdMaker");
    thirdMaker->SetFileName(outFile1);
    
    /*  StPointMaker* pointMaker = new StPointMaker("thirdMaker");
	pointMaker->SetFileName(outFile1);*/
    // 8. Initialize the chain and all makers
    chain->Init(); // This calls Init() for ALL makers
    preEcl->SetClusterConditions("bemc",  4, 4.0, 0.001, 0.1, kFALSE);
    preEcl->SetClusterConditions("bprs",  1, 0.1, 0.001, 0.1, kFALSE);
    preEcl->SetClusterConditions("bsmde", 5, 0.0, 0.001, 0.1, kFALSE);
    preEcl->SetClusterConditions("bsmdp", 5, 0.0, 0.001, 0.1, kFALSE);
    // 9. Event loop
    for (Int_t iev=0; iev<nEvents; iev++) 
      {
	cout << "------------------------------ event "
	     << iev << " -----------------------------" << endl;
	chain->Clear();
	int iret = chain->Make(iev); 
	// This should call the Make() method in ALL makers
	// (including your own maker)
	if (iret||iev==nEvents-1) {
	  // cout << "Bad return code!" << endl;
      break;
	}
      } // Event Loop
    
    chain->Finish(); // This should call the Finish() method in ALL makers
}
// bsub -q star_cas_big -u alexst -L /usr/local/bin/tcsh -o job.log -e err.log root4star -b -q 'pi0.C(200000)'
// for output: bpeek -f

// our triggers
// RUN III --------------------------
// MinBias dAu: 2001, 2003
// HT1 dAu: 2201
// HT2 dAu: 2202
// MinBias pp: 1000
// HT1 pp: 1101
// HT2 pp: 1102
// RUN IV  --------------------------
// MinBias AuAu: 15007 
// HT1 AuAu: 15203 (there's only 1 HighTower threshold for AuAu2004
