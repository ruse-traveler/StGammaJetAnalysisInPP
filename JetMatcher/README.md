# Jet Matcher

This reads in the output of `StMcJetTreeMaker` and `StMuDstJetTreeMaker` and matches each simulated jet to a reconstructed jet to produce a response matrix and jet-matching efficiency. Also matches the particle- and detector-level output of the Pythia Jet Maker to produce a response matrix and jet-matching efficiency.

---

For the output of the Embedding Jet Maker (`StMcJetTreeMaker`/`StMuDstJetTreeMaker`), use `MatchJets.{C,sh}`; and for the particle- and detector-level output of the Pythia Jet Maker (`StFemtoDstMaker`), use `MatchPythiaJets.{C,sh}`. In both cases, the matching parameters are set in the `*.C` file, the input/output files are set in the `*.sh` file, and the matcher is run with the `*.sh` file.

After running the matcher, the normalization of the response matrix is carried out using `PrepareForUnfolding.C` or `PreparePythiaForUnfolding.C`, where the former is to be used on the output of `MatchJets.C` and the latter is to be used on the output of `MatchPythiaJets.C`. Note that since the embedding sample was generated in bins of $\hat{p}_{T}$, `MatchJets.C` must be running over the Embedding Jet Maker output for each $\hat{p}_{T}$ bin (which is automated in `MatchJets.sh`). Then `PrepareForUnfolding.C` sums the output of each bin with the proper weights.

Additional macros which may be useful can be found in the `macros` directory.

### Embedding Jet Maker

```
./MatchJets.sh
root -b -q PrepareForUnfolding.C
```

### Pythia jet Maker
```
./MatchPythiaJets.sh
root -b -q PreparePythiaForUnfolding.C
```
