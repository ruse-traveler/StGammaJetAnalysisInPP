# Embed Jet Maker (StMuDstJetTreeMaker, StTrackEfficiencyCalculator)

These reconstruct jets from the reconstructed events of the Run9 embedding sample stored in FemtoDSTs and calculates the tracking efficiency and resolution.

---

For both classes, compile with `cons.sh` under any STAR version. For the jet maker, input/output files are set in `MakeMuDstJetTree.sh`, and jet finding parameters are set in `MakeMuDstJetTree.C`. This class is easily run with `MakeMuDstJetTree.sh`. The macro `ReadJetTreeDst.C` reads in the output of `MakeMuDstJetTree.C` and produces several jet-level quantities.

For the efficiency calculator, input/output files are set in `DoEfficiencyCalculation.sh`, and the calculation parameters are set in `DoEfficiencyCalculation.C`. This class is easily run with `DoEfficiencyCalculation.sh`. The macro `FinalizeEfficiency.C` uses the output of `DoEfficiencyCalculation.C` to produce a plot of the tracking efficiency.

The macro `CalculateResolution.C` uses the output of the MuDST Matcher to calculate the tracking resolution. Additional macros and scripts which may be helpful can be found in the `macros` and `scripts` directories.

### StMuDstJetTreeMaker

```
./cons.sh
./MakeMuDstJetTree.sh
```

### StTrackEfficiencyCalculator

```
./cons.sh
./DoEfficiencyCalculation.sh
root -b -q FinalizeEfficiency.C
root -b -q CalculateResolution.C
```
