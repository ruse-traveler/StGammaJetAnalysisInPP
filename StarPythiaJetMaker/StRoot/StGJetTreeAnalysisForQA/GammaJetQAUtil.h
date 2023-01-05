#ifndef GammaJetQAUtil_H
#define GammaJetQAUtil_H

//Contains the GammaJetTower Class Declaration
//It converts Tower energy to TLorentzVector
// Nihar Sahoo  .... Jun 20, 2016

#include <stdlib.h>
#include <TROOT.h>
#include "Riostream.h"
#include <TObject.h>

#include <TLorentzVector.h>
#include <TVector3.h>


class GammaJetQAUtil:public TObject{

 private:
  double mass_;
  Option_t *opt_dataSet;

 public:
 GammaJetQAUtil(Option_t *Opt = "option"):opt_dataSet("") { opt_dataSet = Opt;}

    TLorentzVector GetMomentumVectorForTower(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng, const TVector3& Vrtx);
    //TLorentzVector GetMomentumVectorForTower(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng);

    Bool_t IsGoodRun(Int_t runId);
    Bool_t IsGoodTowerIdForTrg(Int_t TwrId);
    Int_t GetRunIndex(Int_t runId);



ClassDef(GammaJetQAUtil, 1)
};

#endif
