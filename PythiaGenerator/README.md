# Pythia Generator

This generates $\pi^{0}$- or $\gamma$-triggered PYTHIA-8 events and saves them in FemtoDSTs.

 - **GammaJetExecute.cc:** This file runs the analysis code defined by `GammaJetAnalysis.h` (this file contains the event loop); things like the pTtrig cut, output filename, etc. are set in this file.
 - **GammaJetExecute.cmnd:** This file contains all the settings for the simulation like the number of events, pThatMin, etc.
 - **GammaJetAnalysis.h:** This is the analysis performed on each event. This is the main piece of code where the various trees are filled.
 - **GammaJetSystem.h:** This file contains a bunch of miscellaneous variables and methods which aid in the analysis.
 - **GammaJetHistory.h:** This file contains the definition of the 'history' class. It provides a brief overview of the event, as well as a complete list of particles which participated in the event and some of their properties.
 - **GammaJet$\star$.h:** -- The remaining files contain the definitions of the Event, Trigger, Track, Vertex, and Parton classes used in `GammaJetAnalysis.h`. These are just containers and helper functions to make the analysis code easier to follow / modify.

Note that despite all the files being named `GammaJet*`, this module does generate both $\gamma$- *and* $\pi^{0}$-triggered events. The trigger is controlled by uncommenting and commenting the `HardQCD:all = on` and `PromptPhoton:all = on` lines in `GammaJetExecute.cmnd` for $\pi^{0}$ and $\gamma$ triggers respectively. (Make sure *both* aren't commented or uncommented...)

---

If anything *other* than `GammaJetExecute.cmnd` has been modified, recompile the module with

```
'make GammaJetExecute'
```

and the run it with

```
./GammaJetExecute.exe GammaJetExecute.cmnd > <outputfile>
```
