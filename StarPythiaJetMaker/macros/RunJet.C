#include <TSystem>
//#include <TCanvas.h>
#include <iostream>
//#endif

using namespace std;

class StGJetTreeAnalysis;


class TString;
class Bichsel;
class TNtuple;
class TFile;
class TH1F;
class TH2F;

static const double PT_TRACK_MAX = 30;//20;// 15;
static const double PT_TRACK_MIN = 0.2;

static const double ETA_TRACK_MAX = 1.0;
static const double ETA_TRACK_MIN = -1.0;

static const Int_t REFMULT_CUT = 315; // unused
//______Jet Input
static const Double_t Jet_R = 0.2;
static const Int_t Remove_N_hardest = 1;

static const Int_t Jet_type = 0;                     //           "FullJet" = 1  and "ChargedJet" = 0

void RunJet(TString lumi)
{

  gSystem->Load( "/opt/star/sl73_gcc485/lib/libfastjet.so");
  gSystem->Load( "/opt/star/sl73_gcc485/lib/libfastjettools.so");
  gSystem->Load("/star/data01/pwg/dmawxc/PythiaData/StarPythia/StarPythiaJetMaker/.sl73_gcc485/lib/libStGJetTreeAnalysis.so");

  //  gSystem->Load( "/opt/star/Xsl64_gcc482/lib/libfastjet.so");
  //  gSystem->Load( "/opt/star/Xsl64_gcc482/lib/libfastjettools.so");
  //  gSystem->Load("StGJetTreeAnalysis");


  cout<<" Macro loaded....."<<endl;
  t= new StGJetTreeAnalysis();
  //  t->Make();
  //  t->Make(R, Int_t Remove_N_hardest,Int_t Jet_type,double PT_TRACK_MAX,double PT_TRACK_MIN, double ETA_TRACK_MAX,double ETA_TRACK_MIN,Int_t REFMULT_CUT);
  t->Make(Jet_R,Remove_N_hardest,Jet_type,PT_TRACK_MAX,PT_TRACK_MIN,ETA_TRACK_MAX,ETA_TRACK_MIN,REFMULT_CUT,lumi);
}
