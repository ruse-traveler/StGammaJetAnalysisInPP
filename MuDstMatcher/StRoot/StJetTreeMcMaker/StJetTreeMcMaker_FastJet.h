#ifndef STAR_StJetTreeMcMaker
#define STAR_StJetTreeMcMaker

//#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"


void RunFastJet(vector<PseudoJet> particles)
{
  double Jet_R = 0.3;
  
  //__________________Jet Definitions_________
  JetDefinition jet_def(antikt_algorithm, Jet_R);

  // jet area definition
  //      Double_t ghost_maxrap = 1.0 + Jet_R; // Fiducial cut for background estimation
  Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
  //      GhostedAreaSpec area_spec(ghost_maxrap);
  //AreaDefinition area_def(active_area, area_spec);
  AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

  vector<PseudoJet> jets;
  // Using  Jet areas in clustering Sequence

  ClusterSequenceArea cs_hard(particles, jet_def, area_def);
  double pt_Min = 0.2;
  vector<PseudoJet> jets_cs = sorted_by_pt(cs_hard.inclusive_jets(pt_Min));

  Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - Jet_R);        //...... fiducial cut for jets from Alex

  // For Active area jets after Fiducial cuts
  jets = Fiducial_cut_selector(jets_cs);


  /////////////// For Background estimator /////////
  JetDefinition jet_def_bkgd(kt_algorithm, Jet_R);

  // Area Definition for bkgd
  AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

  // Selector for bkgd
  //      Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest));
  Selector selector = SelectorAbsEtaMax(1.0 - Jet_R) * (!SelectorNHardest(Remove_N_hardest));
  //      Selector selector = SelectorAbsEtaMax(1.0); // to tes

  // JetMedianBackgroundEstimator definition for bkgd estimation
  JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);

  //Setting Subtractor
  Subtractor subtractor(&bkgd_estimator);

  // Default for higher version of fastjet
  //      subtractor.set_use_rho_m(true);

  //#if FASTJET_VERSION_NUMBER >= 30100
  //      subtractor.set_use_rho_m(true);
  //      subtractor.set_safe_mass(true);
  //#endif
  //Setting bkgd_estimator for events particles
  bkgd_estimator.set_particles(particles);

  // Estimate jet Rho and sigma

  Double_t jets_rho = bkgd_estimator.rho();
  Double_t jets_sigma = bkgd_estimator.sigma();

  // print out some info
  //      cout << "Clustered with " << jet_def.description() << endl;
  // print the jets

  Int_t no0f_jets = jets.size();


  for (Int_t nj = 0; nj < no0f_jets; nj++) {

    double jet_pT= jets[nj].perp();
    double jet_pT_sub= jets[nj].perp() - jets[nj].area()*jets_rho;
    double jet_Area = jets[nj].area();
    //        double jet_Rap = jets[nj].rap();
    double jet_PsudoRap = jets[nj].pseudorapidity();
    //      double jet_Phi = jets[nj].phi();  //jet_Phi: 0 to 2*pi
    double jet_Phi = jets[nj].phi_std(); //jet_Phi: -pi to pi

    vector<PseudoJet> constituents = jets[nj].constituents();
    Int_t nCon = constituents.size();
    //      cout<< "no. of nCon = "<<nCon<<endl;

  }





}



//__________
#endif 
