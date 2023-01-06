# MuDst Matcher

This reads in MuDSTs from the Run9 dijet embedding sample, matches them to the corresponding PYTHIA-6 event record, and saves the reconstructed events in FemtoDSTs and the simulated events in small trees.

---

Compile by running `cons` in `SL14g`. The input/output files are set in `SubmitMatching.sh`, and the class is run with `SubmitMatching.sh` (there are no parameters to change). Additional macros and scripts which may be useful can be found in the `macros` and `scripts` directories.

```
starver SL14g
cons
./SubmitMatching.sh
```
