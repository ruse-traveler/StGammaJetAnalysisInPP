# Jet Folder

This unfolds the measured data (the output of the Data Jet Maker) using the response matrices and jet-matching efficiencies from the Jet Matcher. This module was also used to carry out the closure test of the analysis.

---

compile by running `cons` in any star version. Input/output files and unfolding parameters are set in `DoUnfolding.C`. The class can be easily run with `DoUnfolding.sh`. Additional macros and scripts which may be useful can be found in the `macros` and `scripts` directories.

```
cons
./DoUnfolding.sh
```
