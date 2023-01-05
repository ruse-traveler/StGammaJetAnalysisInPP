#ifndef GammaJetTowerUtil_H
#define GammaJetTowerUtil_H

//Contains the GammaJetTower Class Declaration
//It converts Tower energy to TLorentzVector
// Nihar Sahoo  .... Jun 20, 2016

#include <stdlib.h>
#include <TROOT.h>
#include "Riostream.h"
#include <TObject.h>

#include <TLorentzVector.h>
#include <TVector3.h>


class GammaJetTowerUtil:public TObject{

 private:
  double mass_;

public:
 GammaJetTowerUtil(double m=0.):mass_(m) { }

    TLorentzVector GetMomentumVectorForTower(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng, const TVector3& Vrtx);
  //  TLorentzVector GetMomentumVectorForTower(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng);

ClassDef(GammaJetTowerUtil, 1)
};

#endif
