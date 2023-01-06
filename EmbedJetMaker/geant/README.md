# Embed Jet Maker (StGeantJetTreeMaker)

This reconstructs jets from simulated charged particles from the `StMcEvent` collection stored in FemtoDSTs. Ultimately, this was not used for any physics results.

---

Compile with `cons.sh` under any STAR version. Input/output files and jet finding parameters are set in `MakeGeantJetTree.C`, and the class is easily run with `MakeGeantJetTree.sh`. The macro `ReadJetTreeGnt.C` reads in the output of `MakeGeantJetTree.C` and produces several histograms of jet-level quantities.

```
./cons.sh
./MakeGeantJetTree.sh
```
