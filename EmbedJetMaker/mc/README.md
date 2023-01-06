# Embed Jet Maker (StMcJetTreeMaker)

This reconstructs jets from simulated charged particles from the McArray branch of the MuDSTs stored in a compact tree. 

---

Compile with `cons.sh` under any STAR version. Input/output files are set in `MakeMcJetTree.sh`, and jet finding parameters are set in `MakeMcJetTree.C`. The class is easily run with `MakeMcJetTree.sh`. The macro `ReadJetTreeMc.C` reads in the output of `MakeMcJetTree.C` and produces several jet-level quantities.

```
./cons.sh
./MakeMcJetTree.sh
```
