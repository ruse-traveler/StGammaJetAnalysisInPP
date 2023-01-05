# Data QA Maker

This reads in the output of the Third Maker, reconstructs jets from TPC tracks, and produces several QA plots of event-, trigger-, track-, and jet-level quantities.

---

Compile by running `cons.sh` in any STAR version. Input/output files and jet finding parameters are set in `MakeDataQaPlots.C`. The class can be easily run with `MakeDataQaPlots.sh`

```
./cons.sh
./MakeDataQaPlots.sh
```
