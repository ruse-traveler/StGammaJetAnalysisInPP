# Data Jet Maker

This reads in the output of the Third Maker, reconstructs jets from TPC tracks, and saves them to a compact ROOT tree.

---

Compile by running `cons.sh` (located in `./scripts`) in any STAR version. Input/output files and jet finding parameters are set in `MakeJetTree.C`. The class can be run easily with `MakeJetTree.sh`. The macro `ReadJetTree.C` reads in the output of `MakeJetTree.C` and produces several histograms of jet quantities. Additional macros and scripts which may be useful can be found in the `macros` and `scripts` directories.

```
./cons.sh
./MakeJetTree.sh
```

The Third Maker output for the entire pp dataset may be on RCF at the following path:

```
/star/data01/pwg/dmawxc/JetReco_pp/FullJetTree/merged/pp200r9.merge.root
```
