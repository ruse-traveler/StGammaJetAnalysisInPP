# STAR Pythia Jet Maker

This reconstructs jets from charged particles in $\pi^{0}$- and $\gamma$-triggered events simulated by PYTHIA-6 tuned to STAR data. The output of this module was compared against the unfolded data.

---

Compile using `cons.sh` (located in `./scripts`) under SL14g. Input/output files are set in `RunJetInSerial{02,05}.sh`, and jet finding parameters are set in `RunJet{02,05}.C`. The code is easily run with `RunJetInSerial{02,05}.sh`. The macros `ReadJetTreeSimR{02,05}.C` reads in the output of `RunJet{02,05}.C` respectively and produces several histograms of jet quantities. Additional macros and scripts which may be useful can found in the `macros` and `scripts` directories.

The scripts `RunJetInSerial{02,05}.sh` were set up to run over multiple batches of files simultaneously. Thus one needs to specify 1, 2, 3, or 4 as an argument to the script when running (these correspond to the 4 sets of files which were processed in parallel during the original run of this code).

```
./cons.sh
./RunJetInSerial.sh {1,2,3,4}
```
