// 'TriggerResolutionStudy_RunStarsim.C'
// Derek Anderson
// 08.30.2020
//
// Adapted from macro to instantiate the
// Geant3 from within STAR C++ framework
// and get the starsim prompt. To use it
// do:
//
//   root4star starsim.C
//
// Also adapted from code from Robert
// Licenik and Maria Zurek.

// for generating pion grid
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TDatime.h"
#include "TRandom3.h"

using namespace std;

// particle parameters
static const TString sPar   = "pi0";
static const Float_t pTminP = 6.;
static const Float_t pTmaxP = 50.;
static const Float_t pTminC = 0.;
static const Float_t pTmaxC = 1.5;
static const Float_t hMin   = -0.9;
static const Float_t hMax   = 0.9;
static const Float_t fMin   = -1. * TMath::Pi();
static const Float_t fMax   = TMath::Pi();

// event parameters
static const Float_t zRange[2]  = {-55., 55.};
static const Float_t xyRange[2] = {-2., 2.};
static const Float_t muVtx[3]   = {0.34, -0.05, -4.2};
static const Float_t rmsVtx[3]  = {0.06, 0.04, 62.};
static const Bool_t  useVtx     = true;

// number of particles to generate
static const UInt_t nSetPar = 8;
static const UInt_t nSetPiP = 1;
static const UInt_t nSetPiM = 1;
static const UInt_t nMaxPiC = 2;

// parameters for particle grid
static const Float_t aTow    = 0.05;
static const Float_t aGridH  = 40. * aTow;
static const Float_t aGridF  = 11. * aTow;
static const Float_t dMax    = 12. * aTow;
static const Float_t dMaxH   = 18. * aTow;
static const Float_t dMaxF   = 12. * aTow;
static const UInt_t  nHmax   = (UInt_t) TMath::Ceil(TMath::Abs(hMax - hMin)) / aGridH;
static const UInt_t  nFmax   = (UInt_t) TMath::Floor(TMath::Abs(fMax - fMin)) / aGridF;
static const UInt_t  nGrid   = nHmax * nFmax;
static const Bool_t  useGrid = true;
static const Bool_t  debug   = true;



// declare star classes and makers
class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent *event = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary = 0;

class StarKinematics;
StarKinematics *kinematics = 0;

// instantiate functions
TF1 *ptDist  = 0;
TF1 *etaDist = 0;



// define geometry function
void geometry(TString tag, Bool_t agml = true) {

  TString cmd = "DETP GEOM "; cmd += tag;
  if (!geant_maker) geant_maker = (St_geant_Maker*) chain -> GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //if (agml) command("gexec $STAR_LIB/libxgeometry.so");

}  // end 'geometry(TString, Bool_t)'



// define command function
void command(TString cmd) {

  if (!geant_maker) geant_maker = (St_geant_Maker*) chain -> GetMaker("geant");
  geant_maker -> Do(cmd);

}  // end 'command(TString)'



// define trig function
void trig(Int_t n = 1 ) {

  // event loop
  for (Int_t i = 0; i < n; i++ ) {

    // clear the chain from the previous event
    chain -> Clear();

    // generate vertex
    Bool_t   isInX(false);
    Bool_t   isInY(false);
    Bool_t   isInZ(false);
    Bool_t   isInRange(false);
    Double_t vtx[3];
    if (useVtx) {
      do {
        vtx[0]    = gRandom -> Gaus(muVtx[0], rmsVtx[0]);
        vtx[1]    = gRandom -> Gaus(muVtx[1], rmsVtx[1]);
        vtx[2]    = gRandom -> Gaus(muVtx[2], rmsVtx[2]);
        isInX     = ((vtx[0] > xyRange[0]) && (vtx[0] < xyRange[1]));
        isInY     = ((vtx[1] > xyRange[0]) && (vtx[1] < xyRange[1]));
        isInZ     = ((vtx[2] > zRange[0]) && (vtx[2] < zRange[1]));
        isInRange = (isInX && isInY && isInZ);
      } while (!isInRange);
    }

    // create grid of particles
    if (useGrid) {
      if (debug) {
        cout << "\n ==================================== \n"
             << "  Beginning particle grid generation...\n"
             << "    Constants: particle = " << sPar.Data() << "\n"
             << "      eta     = (" << hMin << ", " << hMax << ")\n"
             << "      phi     = (" << fMin << ", " << fMax << ")\n"
             << "      aTow    = " << aTow << ", aGridH = " << aGridH << ", aGridF = " << aGridF << "\n"
             << "      dMaxH   = " << dMaxH << ", dMaxF = " << dMaxF << "\n"
             << "      nHmax   = " << nHmax << ", nFmax = " << nFmax << "\n"
             << "      nGrid   = " << nGrid << "\n"
             << "      zRange  = (" << zRange[0] << ", " << zRange[1] << ")\n"
             << "      xyRange = (" << xyRange[0] << ", " << xyRange[1] << ")"
             << endl;
        if (useVtx) {
          cout << "      muVtx   = (" << muVtx[0] << ", " << muVtx[1] << ", " << muVtx[2] << ")\n"
               << "      rmsVtx  = (" << rmsVtx[0] << ", " << rmsVtx[1] << ", " << rmsVtx[2] << ")\n"
               << "      vertex  = (" << vtx[0] << ", " << vtx[1] << ", " << vtx[2] << ")"
               << endl;
        }
      }  // end debug

      // make sure particle number is valid
      UInt_t nUsePar = nSetPar;
      UInt_t nUsePiP = nSetPiP;
      UInt_t nUsePiM = nSetPiM;
      if ((nSetPar + nSetPiP + nSetPiM) > nGrid) {
        nUsePar = nGrid - nMaxPiC;
      }
      if ((nSetPiP + nSetPiM) > nMaxPiC) {
        nUsePiP = (UInt_t) ((Float_t) nMaxPiC / 2.);
        nUsePiM = (UInt_t) ((Float_t) nMaxPiC / 2.);
      }
      const UInt_t nPar = nUsePar;
      const UInt_t nPiP = nUsePiP;
      const UInt_t nPiM = nUsePiM;
      if (debug) {
        cout << "    Set particle numbers:\n"
             << "      nSetPar = " << nSetPar << ", nPar = " << nPar << "\n"
             << "      nSetPiP = " << nSetPiP << ", nPiP = " << nPiP << "\n"
             << "      nSetPiM = " << nSetPiM << ", nPiM = " << nPiM
             << endl;
      }  // end debug

      // declare arrays [0 = eta, 1 = phi]
      Bool_t isUsed[nGrid];
      Int_t  kGrid[nGrid][2];
      if (debug) cout << "    Declared arrays." << endl;

      // initialize arrays
      UInt_t nSites = 0;
      for (UInt_t iEta = 0; iEta < nHmax; iEta++) {
        for (UInt_t iPhi = 0; iPhi < nFmax; iPhi++) {
          kGrid[nSites][0] = (Int_t) (iEta - ((Float_t) nHmax / 2.));
          kGrid[nSites][1] = (Int_t) (iPhi - ((Float_t) nFmax / 2.));
          isUsed[nSites]   = false;
          nSites++;
        }  // end phi loop
      }  // end eta loop
      if (debug) cout << "    Initialized arrays: " << nSites << " grid sites." << endl;

      // set random seed
      TDatime *dateTime = new TDatime();
      Int_t    date     = dateTime -> GetDate();
      Int_t    time     = dateTime -> GetTime();
      gRandom -> SetSeed(date * time);
      if (debug) cout << "    Set seed: " << date * time << endl;

      // generate particles
      if (debug) cout << "    Generating particles:" << endl;
      UInt_t nFree = nSites;

      // particle loop
      for (UInt_t iPar = 0; iPar < nPar; iPar++) {

        // get free site index
        UInt_t iRand = (UInt_t) gRandom -> Uniform(0., (Float_t) nFree);
        UInt_t iFree = 0;
        UInt_t iUse  = 0;
        for (UInt_t iSite = 0; iSite < nSites; iSite++) {
          if (!isUsed[iSite]) {
            if (iRand == iFree) {
              iUse          = iSite;
              isUsed[iSite] = true;
              break;
            } else {
              iFree++;
            }
          }  // end if (free)
        }  // end site loop

        // calculate eta and phi
        const Int_t   kEta = kGrid[iUse][0];
        const Int_t   kPhi = kGrid[iUse][1];
        const Float_t eta  = kEta * aGridF;
        const Float_t phi  = kPhi * aGridH;

        // calculate ranges
        const Float_t etaMin = TMath::Min(eta - dMaxH, eta + dMaxH);
        const Float_t etaMax = TMath::Max(eta - dMaxH, eta + dMaxH);
        const Float_t phiMin = TMath::Min(phi - dMaxF, phi + dMaxF);
        const Float_t phiMax = TMath::Max(phi - dMaxF, phi + dMaxF);

        // get values
        const Double_t etaUse = gRandom -> Uniform(etaMin, etaMax);
        const Double_t phiUse = gRandom -> Uniform(phiMin, phiMax);
        const Double_t pTuse  = gRandom -> Uniform(pTminP, pTmaxP);
        const Double_t pXuse  = pTuse * TMath::Cos(phiUse);
        const Double_t pYuse  = pTuse * TMath::Sin(phiUse);
        const Double_t pZuse  = pTuse * TMath::SinH(etaUse);

        // declare particle
        if (useVtx) {
          StarGenParticle *par = kinematics -> AddParticle(sPar.Data());
          par -> SetPx(pXuse);
          par -> SetPy(pYuse);
          par -> SetPz(pZuse);
          par -> SetVx(vtx[0]);
          par -> SetVy(vtx[1]);
          par -> SetVz(vtx[2]);
          if (debug) {
            cout << "      Par[" << iPar << "]: free site " << iRand << "/" << nFree << " = site no. " << iUse << "/" << nSites << "\n"
                 << "              k(eta, phi) = (" << kEta << ", " << kPhi << "); (eta, phi) = (" << eta << ", " << phi << ")\n"
                 << "              p(pT, eta, phi) = (" << pTuse << ", " << etaUse << ", " << phiUse << "); p(pX, pY, pZ) = (" << pXuse << ", " << pYuse << ", " << pZuse << ")"
                 << endl;
          }  // end debug
        } else {
          kinematics -> Kine(1, sPar.Data(), pTminP, pTmaxP, etaMin, etaMax, phiMin, phiMax);
          if (debug) {
            cout << "      Par[" << iPar << "]: free site " << iRand << "/" << nFree << " = site no. " << iUse << "/" << nSites << "\n"
                 << "              k(eta, phi) = (" << kEta << ", " << kPhi << "); (eta, phi) = (" << eta << ", " << phi << ")\n"
                 << "              etaRange = (" << etaMin << ", " << etaMax << "); phiRange = (" << phiMin << ", " << phiMax << ")"
                 << endl;
          }  // end debug
        }
        nFree--;

      }  // end particle loop
      if (debug) cout << "    Generating pi+'s:" << endl;

      // pi+ loop
      for (UInt_t iPiP = 0; iPiP < nPiP; iPiP++) {

        // get free site index
        UInt_t iRand = (UInt_t) gRandom -> Uniform(0., (Float_t) nFree);
        UInt_t iFree = 0;
        UInt_t iUse  = 0;
        for (UInt_t iSite = 0; iSite < nSites; iSite++) {
          if (!isUsed[iSite]) {
            if (iRand == iFree) {
              iUse          = iSite;
              isUsed[iSite] = true;
              break;
            } else {
              iFree++;
            }
          }  // end if (free)
        }  // end site loop

        // calculate eta and phi
        const Int_t   kEta = kGrid[iUse][0];
        const Int_t   kPhi = kGrid[iUse][1];
        const Float_t eta  = kEta * aGridH;
        const Float_t phi  = kPhi * aGridF;

        // calculate ranges
        const Float_t etaMin = TMath::Min(eta - dMaxH, eta + dMaxH);
        const Float_t etaMax = TMath::Max(eta - dMaxH, eta + dMaxH);
        const Float_t phiMin = TMath::Min(phi - dMaxF, phi + dMaxF);
        const Float_t phiMax = TMath::Max(phi - dMaxF, phi + dMaxF);

        // get values
        const Double_t etaUse = gRandom -> Uniform(etaMin, etaMax);
        const Double_t phiUse = gRandom -> Uniform(phiMin, phiMax);
        const Double_t pTuse  = gRandom -> Uniform(pTminC, pTmaxC);
        const Double_t pXuse  = pTuse * TMath::Cos(phiUse);
        const Double_t pYuse  = pTuse * TMath::Sin(phiUse);
        const Double_t pZuse  = pTuse * TMath::SinH(etaUse);

        // declare particle
        if (useVtx) {
          StarGenParticle *piP = kinematics -> AddParticle("pi+");
          piP -> SetPx(pXuse);
          piP -> SetPy(pYuse);
          piP -> SetPz(pZuse);
          piP -> SetVx(vtx[0]);
          piP -> SetVy(vtx[1]);
          piP -> SetVz(vtx[2]);
          if (debug) {
            cout << "      Pi+[" << iPiP << "]: free site " << iRand << "/" << nFree << " = site no. " << iUse << "/" << nSites << "\n"
                 << "              k(eta, phi) = (" << kEta << ", " << kPhi << "); (eta, phi) = (" << eta << ", " << phi << ")\n"
                 << "              p(pT, eta, phi) = (" << pTuse << ", " << etaUse << ", " << phiUse << "); p(pX, pY, pZ) = (" << pXuse << ", " << pYuse << ", " << pZuse << ")"
                 << endl;
          }  // end debug
        } else {
          kinematics -> Kine(1, "pi+", pTminC, pTmaxC, etaMin, etaMax, phiMin, phiMax);
          if (debug) {
            cout << "      Pi+[" << iPiP << "]: free site " << iRand << "/" << nFree << " = site no. " << iUse << "/" << nSites << "\n"
                 << "              k(eta, phi) = (" << kEta << ", " << kPhi << "); (eta, phi) = (" << eta << ", " << phi << ")\n"
                 << "              etaRange = (" << etaMin << ", " << etaMax << "); phiRange = (" << phiMin << ", " << phiMax << ")"
                 << endl;
          }  // end debug
        }
        nFree--;

      }  // end pi+ loop
      if (debug) cout << "    Generating pi-'s:" << endl;

      // pi- loop
      for (UInt_t iPiM = 0; iPiM < nPiM; iPiM++) {

        // get free site index
        UInt_t iRand = (UInt_t) gRandom -> Uniform(0., (Float_t) nFree);
        UInt_t iFree = 0;
        UInt_t iUse  = 0;
        for (UInt_t iSite = 0; iSite < nSites; iSite++) {
          if (!isUsed[iSite]) {
            if (iRand == iFree) {
              iUse          = iSite;
              isUsed[iSite] = true;
              break;
            } else {
              iFree++;
            }
          }  // end if (free)
        }  // end site loop

        // calculate eta and phi
        const Int_t   kEta = kGrid[iUse][0];
        const Int_t   kPhi = kGrid[iUse][1];
        const Float_t eta  = kEta * aGridH;
        const Float_t phi  = kPhi * aGridF;

        // calculate ranges
        const Float_t etaMin = TMath::Min(eta - dMaxH, eta + dMaxH);
        const Float_t etaMax = TMath::Max(eta - dMaxH, eta + dMaxH);
        const Float_t phiMin = TMath::Min(phi - dMaxF, phi + dMaxF);
        const Float_t phiMax = TMath::Max(phi - dMaxF, phi + dMaxF);

        // get values
        const Double_t etaUse = gRandom -> Uniform(etaMin, etaMax);
        const Double_t phiUse = gRandom -> Uniform(phiMin, phiMax);
        const Double_t pTuse  = gRandom -> Uniform(pTminC, pTmaxC);
        const Double_t pXuse  = pTuse * TMath::Cos(phiUse);
        const Double_t pYuse  = pTuse * TMath::Sin(phiUse);
        const Double_t pZuse  = pTuse * TMath::SinH(etaUse);

        // declare particle
        if (useVtx) {
          StarGenParticle *piM = kinematics -> AddParticle("pi-");
          piM -> SetPx(pXuse);
          piM -> SetPy(pYuse);
          piM -> SetPz(pZuse);
          piM -> SetVx(vtx[0]);
          piM -> SetVy(vtx[1]);
          piM -> SetVz(vtx[2]);
          if (debug) {
            cout << "      Pi-[" << iPiM << "]: free site " << iRand << "/" << nFree << " = site no. " << iUse << "/" << nSites << "\n"
                 << "              k(eta, phi) = (" << kEta << ", " << kPhi << "); (eta, phi) = (" << eta << ", " << phi << ")\n"
                 << "              p(pT, eta, phi) = (" << pTuse << ", " << etaUse << ", " << phiUse << "); p(pX, pY, pZ) = (" << pXuse << ", " << pYuse << ", " << pZuse << ")"
                 << endl;
          }  // end debug
        } else {
          kinematics -> Kine(1, "pi-", pTminC, pTmaxC, etaMin, etaMax, phiMin, phiMax);
          if (debug) {
            cout << "      Pi-[" << iPiM << "]: free site " << iRand << "/" << nFree << " = site no. " << iUse << "/" << nSites << "\n"
                 << "              k(eta, phi) = (" << kEta << ", " << kPhi << "); (eta, phi) = (" << eta << ", " << phi << ")\n"
                 << "              etaRange = (" << etaMin << ", " << etaMax << "); phiRange = (" << phiMin << ", " << phiMax << ")"
                 << endl;
          }  // end debug
        }
        nFree--;

      }  // end pi- loop
      if (debug) {
        cout << "  Ending particle grid generation.\n"
             << " ================================= \n"
             << endl;
      }  // end debug
    }  // end pion grid generation

    // generate single particle
    if (!useGrid) {
      if (debug) {
        cout << "\n ======================================== \n"
             << "  Beginning particle non-grid generation...\n"
             << "    Constants: particle = " << sPar.Data() << "\n"
             << "      eta   = (" << hMin << ", " << hMax << ")\n"
             << "      phi   = (" << fMin << ", " << fMax << ")\n"
             << "      aTow  = " << aTow << ", aGridF = " << aGridF << ", aGridH = " << aGridH << "\n"
             << "      dMaxF = " << dMaxF << ", dMaxH = " << dMaxH
             << endl;
        if (useVtx) {
          cout << "      muVtx   = (" << muVtx[0] << ", " << muVtx[1] << ", " << muVtx[2] << ")\n"
               << "      rmsVtx  = (" << rmsVtx[0] << ", " << rmsVtx[1] << ", " << rmsVtx[2] << ")\n"
               << "      vertex  = (" << vtx[0] << ", " << vtx[1] << ", " << vtx[2] << ")"
               << endl;
        }
      }  // end debug

      // set random seed
      TDatime *dateTime = new TDatime();
      Int_t    date     = dateTime -> GetDate();
      Int_t    time     = dateTime -> GetTime();
      gRandom -> SetSeed(date * time);
      if (debug) cout << "    Set seed: " << date * time << endl;

      // determine particle site
      const Float_t etaPar = gRandom -> Uniform(hMin, hMax);
      const Float_t phiPar = gRandom -> Uniform(fMin, fMax);

      // calculate particle ranges
      const Float_t etaParLo = TMath::Min(etaPar - dMaxH, etaPar + dMaxH);
      const Float_t etaParHi = TMath::Max(etaPar - dMaxH, etaPar + dMaxH);
      const Float_t phiParLo = TMath::Min(phiPar - dMaxF, phiPar + dMaxF);
      const Float_t phiParHi = TMath::Max(phiPar - dMaxF, phiPar + dMaxF);

      // get values
      const Double_t etaUse = gRandom -> Uniform(etaParLo, etaParHi);
      const Double_t phiUse = gRandom -> Uniform(phiParLo, phiParHi);
      const Double_t pTuse  = gRandom -> Uniform(pTminP, pTmaxP);
      const Double_t pXuse  = pTuse * TMath::Cos(phiUse);
      const Double_t pYuse  = pTuse * TMath::Sin(phiUse);
      const Double_t pZuse  = pTuse * TMath::SinH(etaUse);

      // declare particle
      if (useVtx) {
        StarGenParticle *par = kinematics -> AddParticle(sPar.Data());
        par -> SetPx(pXuse);
        par -> SetPy(pYuse);
        par -> SetPz(pZuse);
        par -> SetVx(vtx[0]);
        par -> SetVy(vtx[1]);
        par -> SetVz(vtx[2]);
        if (debug) {
          cout << "      Par: p(pT, eta, phi) = (" << pTuse << ", " << etaUse << ", " << phiUse << ")\n"
               << "           p(pX, pY, pZ)   = (" << pXuse << ", " << pYuse << ", "<< pZuse << ")"
               << endl;
        }  // end debug
      } else {
        kinematics -> Kine(1, sPar.Data(), pTminP, pTmaxP, etaParLo, etaParHi, phiParLo, phiParHi);
        if (debug) {
          cout << "      Par: pTpar = (" << pTminP << ", " << pTmaxP << "); (eta, phi) = (" << etaPar << ", " << phiPar << ")\n"
               << "           etaRange = (" << etaParLo << ", " << etaParHi << "); phiRange = (" << phiParLo << ", " << phiParHi << ")"
               << endl;
        }  // end debug
      }

      // determine piC site
      Bool_t  foundGoodSite(false);
      Float_t etaPi(0.);
      Float_t phiPi(0.);
      do {
 
        // select site
        etaPi = gRandom -> Uniform(hMin, hMax);
        phiPi = gRandom -> Uniform(fMin, fMax);

        // check if near generated particle
        const Float_t dEta = TMath::Abs(etaPar - etaPi);
        const Float_t dPhi = TMath::Abs(phiPar - phiPi);
        const Float_t dSep = TMath::Sqrt((dEta * dEta) + (dPhi * dPhi));
        if (dSep < dMax) {
          foundGoodSite = false;
        } else {
          foundGoodSite = true;
        }
      }  while (!foundGoodSite);

      // calculate piC ranges
      const Float_t etaPiLo = TMath::Min(etaPi - dMaxH, etaPi + dMaxH);
      const Float_t etaPiHi = TMath::Max(etaPi - dMaxH, etaPi + dMaxH);
      const Float_t phiPiLo = TMath::Min(phiPi - dMaxF, phiPi + dMaxF);
      const Float_t phiPiHi = TMath::Max(phiPi - dMaxF, phiPi + dMaxF);

      // get values
      const Double_t etaPiC = gRandom -> Uniform(etaPiLo, etaPiHi);
      const Double_t phiPiC = gRandom -> Uniform(phiPiLo, phiPiHi);
      const Double_t pTpiC  = gRandom -> Uniform(pTminC, pTmaxC);
      const Double_t pXpiC  = pTpiC * TMath::Cos(phiPiC);
      const Double_t pYpiC  = pTpiC * TMath::Sin(phiPiC);
      const Double_t pZpiC  = pTpiC * TMath::SinH(etaPiC);

      // select charge and declare particles
      const Float_t charge = gRandom -> Uniform(0., 1.);
      if (charge < 0.5) {
        if (useVtx) {
          StarGenParticle *piM = kinematics -> AddParticle("pi-");
          piM -> SetPx(pXpiC);
          piM -> SetPy(pYpiC);
          piM -> SetPz(pZpiC);
          piM -> SetVx(vtx[0]);
          piM -> SetVy(vtx[1]);
          piM -> SetVz(vtx[2]);
          if (debug) {
            cout << "      Pi-: p(pT, eta, phi) = (" << pTpiC << ", " << etaPiC << ", " << phiPiC << ")\n"
                 << "           p(pX, pY, pZ)   = (" << pXpiC << ", " << pYpiC << ", "<< pZpiC << ")"
                 << endl;
          }  // end debug
        } else {
          kinematics -> Kine(1, "pi-", pTminC, pTmaxC, etaPiLo, etaPiHi, phiPiLo, phiPiHi);
          if (debug) {
            cout << "      Pi-: pTpiM = (" << pTminC << ", " << pTmaxC << "); (eta, phi) = (" << etaPi << ", " << phiPi << ")\n"
                 << "           etaRange = (" << etaPiLo << ", " << etaPiHi << "); phiRange = (" << phiPiLo << ", " << phiPiHi << ")"
                 << endl;
          }  // end debug
        }
      } else {
        if (useVtx) {
          StarGenParticle *piP = kinematics -> AddParticle("pi+");
          piP -> SetPx(pXpiC);
          piP -> SetPy(pYpiC);
          piP -> SetPz(pZpiC);
          piP -> SetVx(vtx[0]);
          piP -> SetVy(vtx[1]);
          piP -> SetVz(vtx[2]);
          if (debug) {
            cout << "      Pi+: p(pT, eta, phi) = (" << pTpiC << ", " << etaPiC << ", " << phiPiC << ")\n"
                 << "           p(pX, pY, pZ)   = (" << pXpiC << ", " << pYpiC << ", "<< pZpiC << ")"
                 << endl;
          }  // end debug
        } else {
          kinematics -> Kine(1, "pi+", pTminC, pTmaxC, etaPiLo, etaPiHi, phiPiLo, phiPiHi);
          if (debug) {
            cout << "      Pi+: pTpiP = (" << pTminC << ", " << pTmaxC << "); (eta, phi) = (" << etaPi << ", " << phiPi << ")\n"
                 << "           etaRange = (" << etaPiLo << ", " << etaPiHi << "); phiRange = (" << phiPiLo << ", " << phiPiHi << ")"
                 << endl;
          }  // end debug
        }
      }

      if (debug) {
        cout << "  Ending particle non-grid generation.\n"
             << " ===================================== \n"
             << endl;
      }  // end debug
    }  // end non-grid generation

    // generate the event
    chain -> Make();

    // print the event
    _primary -> event() -> Print();

  }  // end event loop
}  // end 'trig(Int_t)'



// define kinematics function
void Kinematics() {

  // load star kinematics library
  //gSystem->Load("libStarGeneratorPoolPythia6_4_23.so");
  gSystem -> Load("libKinematics.so");
  kinematics = new StarKinematics();

  // add generator to run
  _primary -> AddGenerator(kinematics);

}  // end 'Kinematics()'



// main starsim function
void TriggerResolutionStudy_RunStarsim(Int_t nEvents = 10, Int_t rngSeed = 1234, TString outName = "") {

  // load geometry for geant
  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = "y2009d geant gstar usexgeom agml ";
    bfc(0, simple);
  }

  // load relevant libraries
  gSystem -> Load("libVMC.so");
  gSystem -> Load("StarGeneratorUtil.so");
  gSystem -> Load("StarGeneratorEvent.so");
  gSystem -> Load("StarGeneratorBase.so");
  gSystem -> Load("libMathMore.so");  
  gSystem -> Load("xgeometry.so");

  // set up RNG seed and map all ROOT TRandom here
  StarRandom::seed(rngSeed);
  StarRandom::capture();
 
  // NOTE: Create any primary event generators and
  // insert them before the geant maker

  // create primary event generators
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName(Form("kinematics_%s.starsim.root", outName.Data()));
    chain -> AddBefore("geant", _primary);
  }
  Kinematics();

  // NOTE: initialize primary event generator and
  // all sub makers here

  // initilialize primary event generators
  _primary -> Init();

  // set up geometry and set starsim to use agusread for input
  geometry("y2009d");
  command("gkine -4 0");
  command(Form("gfile o kinematics_%s.starsim.fzd", outName.Data()));
 
  // set up pT distribution
  // [not relevant for now, Derek, 09.17.2020]
  //Double_t pt0 = 3.0;
  //ptDist = new TF1("ptDist", "(x/[0])/(1+(x/[0])^2)^6", 0.0, 15.0);
  //ptDist -> SetParameter(0, pt0);

  // set up eta distribution
  // [not relevent for now, Derek, 09.17.2020]
  //etaDist = new TF1("etaDist","-TMath::Erf(x+2.6)*TMath::Erf(x-2.6)",-1,1);

  // trigger on nevents
  trig(nEvents);

  // make sure that starsim exits properly
  command("call agexit");

}  // end 'TriggerResolutionStudy_RunStarsim(Int_t, Int_t, TString)'

// End -------------------------------------------------------------------------
