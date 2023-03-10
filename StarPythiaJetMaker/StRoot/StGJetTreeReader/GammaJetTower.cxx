#include "GammaJetTower.h"

//Contains the Member Functions for the TRACK Class                             

ClassImp(GammaJetTower);

//_____________________________________________________________________         
GammaJetTower::GammaJetTower():TObject()
			      ,TwrId(-999)
			      ,TwrEng ( -999.)
			      ,TwrEta (-999.)
			      ,TwrPhi( -999.)
			      ,TwrADC (-999.)
			      ,TwrMatchIdnex(-999)
			      ,NoOfmatchedTrk(-999)
			      ,TwrPx(-999.)
			      ,TwrPy(-999.)
			      , TwrPz(-999.)
			      , fNAssocTracks(0)
			       //			      , fMatchedTracks(1)
{

  for(int i=0; i < 10; i++){ fMatchedTracksArray_[i]=-9.;   fMatchedTracksArray_P[i]=-9.;}

  
}

//_____________________________________________________________________         
GammaJetTower::~GammaJetTower(){
  
  //GammaJetTower Destructor                                                            
  //  TwrMatchP.clear();
  fNAssocTracks = 0;
  //  fMatchedTracks.Set(1);
  //  fMatchedTracks.Reset(-1);
}

//_____________________________________________________________________         
GammaJetTower::GammaJetTower(GammaJetTower  &towerarray)
  :TObject(towerarray)
  ,TwrId(towerarray.TwrId)
  ,TwrEng(towerarray.TwrEng)
  ,TwrEta(towerarray.TwrEta)
  ,TwrPhi(towerarray.TwrPhi)
  ,TwrADC(towerarray.TwrADC)
  ,TwrMatchIdnex(towerarray.TwrMatchIdnex)
  ,NoOfmatchedTrk(towerarray.NoOfmatchedTrk)
  ,TwrPx(towerarray.TwrPx)
  ,TwrPy(towerarray.TwrPy)
  , TwrPz(towerarray.TwrPz)
  , fNAssocTracks(towerarray.fNAssocTracks)
   //  ,fMatchedTracksArray_(towerarray.fMatchedTracksArray_)
 {
  
   //   fMatchedTracks.Copy(towerarray.fMatchedTracks);
}

//_____________________________________________________________________         
void GammaJetTower::PrintTower(){

  
}

//____________________________
void GammaJetTower::AddMatchedTrack(Int_t idx)
{
  /*
  //
  // Add new index
  // Resize the array if necessary
  //
  if (fNAssocTracks < fMatchedTracks.GetSize())
    {
      fMatchedTracks[fNAssocTracks++] = idx;      
    }
  else
    {
      //Int_t newsize = fMatchedTracks.GetSize() * 2;
      //since the N of matched tracks should drop steeply we can be modest
      //Int_t newsize = fMatchedTracks.GetSize() + 2;
      Int_t newsize = fMatchedTracks.GetSize() + 1;
      //      __DEBUG(3, Form("Resizing indexes array from %d to %d", fMatchedTracks.GetSize(), newsize));
      
      fMatchedTracks.Set(newsize);
      AddMatchedTrack(idx);
    }
  */
}
