# Trigger Weight Calculator

This takes the output of the Trigger Energy Scale Calculator and computes a set of weights which emulate the effect of the trigger energy scale/resolution to be applied to the simulated jet spectra compared against the unfolded data.

---

The macro `CalculatePythiaWeights.C` reads in 3 inputs: the measured trigger energy spectrum (from `ReadJetTree.C`), the qT distribution from the Particle Gun (from `StTriggerEnergyScaleCalculator`), and the trigger energy spectrum from PYTHIA (from `ReadJetTreeSim.C`). Uses these, it calculates a set of weights to apply to the corresponding PYTHIA jets. All parameters and input/output files are set in `CalculatePythiaWeights.C`.

The macro `MakeEtWeightPlots.C` produces a plot of the calculated weights for all 3 bins of trigger energy. The other 2 macros were only used in development.

```
root -b -q CalculatePythiaWeights.C
root -b -q MakeEtWeightPlots.C
```
