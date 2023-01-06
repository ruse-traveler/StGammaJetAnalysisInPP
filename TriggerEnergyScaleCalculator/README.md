# Trigger Energy Scale Calculator

This takes the output of the Particle Gun and matches the reconstructed $\pi^{0}$ and $\gamma$ to their simulated counterparts, providing a measure of the energy scale and resolution of our measured triggers.

---

Compile using `cons.sh` (located in`./scripts`) under any STAR version. Input/output files are set in `DoTriggerEnergyScaleCalculation.sh`, and calculation parameters are set in `DoTriggerEnergyScaleCalculation.C`. The class is easily run with `DoTriggerEnergyScaleCalculations.sh`. The macro `PrepareTriggerMatrixAndEfficiency.C` reads in the output of `DoTriggerEnergyScaleCalculation.C` and prepares several key figures summarizing the Trigger Energy Scale and Resolution.

Additional macros and scripts which may be helpful can be found in the `macros` and `scripts` directories.

```
./cons.sh
./DoTriggerEnergyScaleCalculation.sh
root -b -q PrepareTriggerEnergyScaleCalculation.C
```
