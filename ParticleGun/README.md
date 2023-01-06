# Particle Gun

This module uses starsim to simulate single $\pi^{0}$ and $\gamma$, pass them through a GEANT-3 simulation of STAR, and reconstruct the output using the Third Maker.

---

The Particle Gun simulation proceeds in 2 steps: (1) first, grids of single pi0 or photons are simulated using starsim, passed through a GEANT simulation of STAR, and reconstructed using the STAR BFC; (2) and then, triggers are reconstructed using the Third Maker. Step 1 is carried out in the `TriggerResolutionStudy_RunStarsim.C` and and `TriggerResolutionStudy_RunBfc.C` macros, and step 2 is carried out in the `TriggerResolutionStudy_RunThirdmaker.C` macro.

Parameters of the simulation are set in `TriggerResolutionStudy_RunStarsim.C`, the number of events per file (process) is set in `TriggerResolutionStudy_Submit.sh` while the number of files to generate (processes to run) is set in `TriggerResolutionStudy_Submit.xml`. The simulation and BFC are run with `TriggerResolutionStudy_Submit.sh`.

Then the Third Maker is run on the output of step 1 with `TriggerResolutionStudy_RunThirdmaker.sh`, in which the input/output files are set. The Third Maker needs to be compiled under `SL14g` using `cons` before being run.

### Starsim + BFC

```
./TriggerResolutionStudy_Submit.sh
```

### Third Maker

```
starver SL14g
cons
./TriggerResolutionStudy_RunThirdmaker.sh
```
