<script
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
  type="text/javascript">
</script>

# StGammaJetAnalysisInPP

This repository compiles the major pieces of code from the STAR $pp$ $\gamma$+jet analysis. Each module is detailed below, where the order of the list corresponds to the overall flow of the analysis.

  1. **Third Maker:** Reads in STAR $pp$ data in the form of MuDSTs, and constructs compact ROOT trees (referred to internally as femtoDSTs) containing $pp$ collisions tagged by potential energetic $pi^{0}$ and direct photons.
  2. **Data QA Maker:**
  3. **Data Jet Maker:**
  4. **MuDst Matcher:**
  5. **Embed Jet Maker:**
    - **StGeantJetTreeMaker:**
    - **StMcJetTreeMaker:**
    - **StMuDstJetTreeMaker:**
    - **StTrackEfficiencyCalculator:**
  6. **Jet Matcher:**
  7. **Jet Folder:**
  8. **Pythia Jet Maker:**
  9. **Systematic Uncertainty Calculator:**
  10. **Particle Gun:**
  11. **Trigger Energy Scale Calculator:**
  12. **Trigger Weight Calculator:**
  13. **Star Pythia Jet Maker:**

The various Jet Maker modules make use of the [FastJet library](http://fastjet.fr) developed by Matteo Cacciari, Gavin Salam, and Gregory Soyez; the Jet Folder modules make use of the [RooUnfold library](https://gitlab.cern.ch/RooUnfold/RooUnfold) developed by Tim Adye; and all modules make use of the [ROOT framework](https://root.cern.ch). Naturally, several modules require the use of [STAR core software](https://github.com/star-bnl).

# **TODO:**
  - [minor] Update driver macro in MuDstMatcher to work with corresponding class.
  - [minor] Add the following packages
    - Pythia8 Generator
  - [major] Add documentation to sub-directories
