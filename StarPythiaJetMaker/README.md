# STAR Pythia Jet Maker

This reconstructs jets from charged particles in $\pi^{0}$- and $\gamma$-triggered events simulated by PYTHIA-6 tuned to STAR data. The output of this module was compared against the unfolded data.

---

Compile using `cons.sh` (located in `./scripts`) under SL14g. Input/output files are set in `RunJetInSerial{02,05}.sh`, and jet finding parameters are set in `RunJet{02,05}.C`. The code is easily run with `RunJetInSerial{02,05}.sh`. Additional macros and scripts which may be useful can found in the `macros` and `scripts` directories.

The scripts `RunJetInSerial{02,05}.sh` were set up to run over multiple batches of files simultaneously. Thus one needs to specify 1, 2, 3, or 4 as an argument to the script when running (these correspond to the 4 sets of files which were processed in parallel during the original run of this code).

In order to enhance statistics, the $\pi^{0}$-triggered STAR-tuned PYTHIA-6 sample was simulated in bins of partonic pT. Thus, these jets need to be combined into a single spectrum after being reconstructed. Histograms are created for each pT bin using `ReadJetTreeSimR{02,05}.C`: trigger and jet acceptance parameters are set in `ReadJetTreeSimR{02,05}.C`, and input/output files are set in `ReadJetTreeSimR{02,05}.sh`, which are used to run the `*.C` files. After running these macros, the jets are summmed using `DoPtHatSum.C`.

```
starver SL14g
./cons.sh
./RunJetInSerial.sh {1,2,3,4}
ReadJetTreeSimR{02,05}.sh
root -b -q DoPtHatSum.C
```

(Note that the last line is only necessary for the $\pi^{0}$-triggered jets.)
