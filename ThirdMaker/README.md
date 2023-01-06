# Third Maker

This reads in STAR $p+p$ data in the form of MuDSTs, and constructs compact ROOT trees (referred to internally as FemtoDSTs) containing $p+p$ collisions tagged by potential energetic $\pi^{0}$ or direct photon based on the [algorithm](https://doi.org/10.1103/PhysRevC.82.034909) developed by Ahmed Hamed and Saskia Mioduszewski. This algorithm was adapted to save output in the FemtoDST format by Nihar Sahoo.

---

Compile using `cons` in SL14g. Input/output files are set in `doFullJetTree.C` or `GammaJet_tree.xml`, and either can be used to run `StJetTreeThirdMaker`. The former can be used to process a single input file, and the latter can be used to process a large sample using star-submit.

The corresponding source code, driver macro, and XML file for a previous version of the Third Maker are included in this module as well. The previous version saves the output as a set of ROOT TNtuples while the this version saves the output as TTrees (the FemtoDSTs).

```
starver SL14g
cons
root4star -b -q doFullJetTree.C
star-submit GammaJet_tree.xml
```

(Note: do either line 3 *or* line 4.)
