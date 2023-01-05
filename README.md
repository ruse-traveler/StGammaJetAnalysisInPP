# StGammaJetAnalysisInPP

This repository compiles the major pieces of code from the STAR $p+p$ $\gamma$+jet analysis. Each module is detailed below, where the order of the list corresponds to the overall flow of the analysis.

  1. **Third Maker:** This reads in STAR $p+p$ data in the form of MuDSTs, and constructs compact ROOT trees (referred to internally as FemtoDSTs) containing $p+p$ collisions tagged by potential energetic $\pi^{0}$ or direct photon based on the [algorithm](https://doi.org/10.1103/PhysRevC.82.034909) developed by Ahmed Hamed and Saskia Mioduszewski. This algorithm was adapted to save output in the FemtoDST format by Nihar Sahoo.
  2. **Data Jet Maker:** This reads in the output of the Third Maker, reconstructs jets from TPC tracks, and saves them to a compact ROOT tree.
  3. **Data QA Maker:** This also reads in the output of the Third Maker, reconstructs jets from TPC tracks, and produces several QA plots of event-, trigger-, track-, and jet-level quantities.
  4. **MuDst Matcher:** This reads in MuDSTs from the Run9 dijet embedding sample, matches them to the corresponding PYTHIA-6 event record, and saves the reconstructed events in FemtoDSTs and the simulated events in small trees.
  5. **Embed Jet Maker:** These three modules read in the output of the MuDST Matcher and reconstructs jets from various stages of the embedding sample.
     - **StGeantJetTreeMaker:** This reconstructs jets from simulated charged particles from the `StMcEvent` collection stored in FemtoDSTs. Ultimately, this was not used for any physics results.
     - **StMcJetTreeMaker:** This reconstructs jets from simulated charged particles from the McArray branch of the MuDSTs stored in a compact tree. 
     - **StMuDstJetTreeMaker:** This reconstructs jets from the reconstructed events of the Run9 embedding sample stored in FemtoDSTs.
     - **StTrackEfficiencyCalculator:** This reads in the same events as the `StMuDstJetTreeMaker` and calculates the tracking efficiency and resolution.
  6. **Pythia Jet Maker:** This reads in $\pi^{0}$- or $\gamma$-triggered PYTHIA-8 events saved in FemtoDSTs and reconstructs jets from charged particles. Also used to apply the calculated tracking efficiency and resolution of `StTrackEfficiencyCalculator` to the PYTHIA-8 events.
  7. **Jet Matcher:** This reads in the output of `StMcJetTreeMaker` and `StMuDstJetTreeMaker` and matches each simulated jet to a reconstructed jet to produce a response matrix and jet-matching efficiency. Also matches the particle- and detector-level output of the Pythia Jet Maker to produce a response matrix and jet-matching efficiency.
  8. **Jet Folder:** This unfolds the measured data (the output of the Data Jet Maker) using the response matrices and jet-matching efficiencies from the Jet Matcher. This module was also used to carry out the closure test of the analysis.
  9. **Systematic Uncertainty Calculator:** This takes the output of the Jet Folder and calculates a systematic uncertainty to be applied to the unfolded data.
  10. **Particle Gun:** This module uses starsim to simulate single $\pi^{0}$ and $\gamma$, pass them through a GEANT-3 simulation of STAR, and reconstruct the output using the Third Maker.
  11. **Trigger Energy Scale Calculator:** This takes the output of the Particle Gun and matches the reconstructed $\pi^{0}$ and $\gamma$ to their simulated counterparts, providing a measure of the energy scale and resolution of our measured triggers.
  12. **Trigger Weight Calculator:** This takes the output of the Trigger Energy Scale Calculator and computes a set of weights which emulate the effect of the trigger energy scale/resolution to be applied to the simulated jet spectra compared against the unfolded data.
  13. **Star Pythia Jet Maker:** This reconstructs jets from charged particles in $\pi^{0}$- and $\gamma$-triggered events simulated by PYTHIA-6 tuned to STAR data. The output of this module was compared against the unfolded data.

The various Jet Maker modules make use of the [FastJet library](http://fastjet.fr) developed by Matteo Cacciari, Gavin Salam, and Gregory Soyez; the Jet Folder modules make use of the [RooUnfold library](https://gitlab.cern.ch/RooUnfold/RooUnfold) developed by Tim Adye; and all modules make use of the [ROOT framework](https://root.cern.ch). Naturally, several modules require the use of [STAR core software](https://github.com/star-bnl).

## **TODO:**
  - [minor] Update driver macro in MuDstMatcher to work with corresponding class.
  - [minor] Add the following packages
    - Pythia8 Generator
  - [major] Add documentation to sub-directories
