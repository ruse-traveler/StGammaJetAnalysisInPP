# Pythia Jet Maker

This reads in $\pi^{0}$- or $\gamma$-triggered PYTHIA-8 events saved in FemtoDSTs and reconstructs jets from charged particles. Also used to apply the calculated tracking efficiency and resolution of `StTrackEfficiencyCalculator` to the PYTHIA-8 events. Due to historical reasons, the output jet tree is erroneously referred to as a FemtoDST in the code.

---

Compile using `cons.sh` (located in `./scripts`) in any STAR version. Input/output files are set in `MakeFemtoDst.sh`, and jet finding parameters are set in `MakeFemtoDst.C`. The code is easily run with `MakeFemtoDst.sh`. The macro `ReadJetTreeSim.C` reads in the output of `MakeFemtoDst.C` and produces several histograms of jet quantities. Additional macros and scripts which may be useful can found in the `macros` and `scripts` directories.

`MakeFemtoDst.sh` takes 2 arguments: `particle` vs. `detector`, and `charged` vs. `full` (in that order). The 1st argument tells `StFemtoDstMaker` whether or not to apply the parameterized tracking efficiency and resolution from `StTrackEfficiencyCalculator`, and the 2nd tells `StFemtoDstMaker` to make jets from just charged particles or from both and charged and neutral particles. For example, if one wanted to make detector-level charged jets, one would issue the command `./MakeFemtoDst.sh detector charged`.

```
./cons.sh
./MakeFemtoDst.sh <particle, detector> <charged, full>
``` 
