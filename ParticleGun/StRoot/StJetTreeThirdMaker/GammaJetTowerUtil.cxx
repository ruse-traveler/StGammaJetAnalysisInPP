#include "GammaJetTowerUtil.h"

//Contains the Member Functions for the TRACK Class

ClassImp(GammaJetTowerUtil);

//_____________________________________________________________________
//TLorentzVector GetMomentumVectorForTower()(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng, TVector3 Vrtx);
TLorentzVector GammaJetTowerUtil::GetMomentumVectorForTower(double BEMC_R, double Twr_eta, double Twr_phi, double Twr_eng, const TVector3& Vrtx)
{
  //  TVector3 Vrtx;
  TVector3 towerLocation;
  towerLocation.SetPtEtaPhi(BEMC_R, Twr_eta, Twr_phi); 

  TVector3 VertexLocation;
  VertexLocation = Vrtx;

  TVector3 momentum = towerLocation - VertexLocation;
  
  double pMag = (Twr_eng > mass_) ? sqrt(pow(Twr_eng,2) - pow(mass_,2)) : Twr_eng;
  momentum.SetMag(pMag);
  TLorentzVector twrV(momentum.x(), momentum.y(), momentum.z(), Twr_eng);

  return twrV;

}
