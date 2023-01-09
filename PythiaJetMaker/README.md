# Pythia Jet Maker

This reads in $\pi^{0}$- or $\gamma$-triggered PYTHIA-8 events saved in FemtoDSTs and reconstructs jets from charged particles. Also used to apply the calculated tracking efficiency and resolution of `StTrackEfficiencyCalculator` to the PYTHIA-8 events. Due to historical reasons, the output jet tree is erroneously referred to as a FemtoDST in the code.

---

Compile using `cons.sh` (located in `./scripts`) in any STAR version. Input/output files are set in `MakeFemtoDst.sh`, and jet finding parameters are set in `MakeFemtoDst.C`. The code is easily run with `MakeFemtoDst.sh`. The macro `ReadJetTreeSim.C` reads in the output of `MakeFemtoDst.C` and produces several histograms of jet quantities. Additional macros and scripts which may be useful can found in the `macros` and `scripts` directories.

`MakeFemtoDst.sh` takes 2 arguments: `particle` vs. `detector`, and `charged` vs. `full` (in that order). The 1st argument tells `StFemtoDstMaker` whether or not to apply the parameterized tracking efficiency and resolution from `StTrackEfficiencyCalculator`, and the 2nd tells `StFemtoDstMaker` to make jets from just charged particles or from both and charged and neutral particles. For example, if one wanted to make detector-level charged jets, one would issue the command `./MakeFemtoDst.sh detector charged`.

```
./cons.sh
./MakeFemtoDst.sh {particle, detector} {charged, full}
```

There are two sets of input for this module, both of which can be found on RCF. The first is a large sample of $\pi^{0}$- and $\gamma$-triggered PYTHIA-8 events with partonic pT greater than 4 GeV/c, and the second is a sample of $\pi^{0}$-triggered PYTHIA-8 events generated in *bins* of partonic pT. The first was used extensively for testing and for comparison against the unfolded data before the STAR-tuned PYTHIA-6 distributions were available, and the second was used to generate the response matrices and jet-matching efficiencies for assessing the tracking systematic uncertainties.

### > 4 GeV/c Input
```
/star/data01/pwg/dmawxc/JetReco_pp/PythiaJetMaker/input/Pythia/*.merged.root
```

### Binned Input
```
/star/data01/pwg/dmawxc/PythiaData/PtHatBins/*.root
```
