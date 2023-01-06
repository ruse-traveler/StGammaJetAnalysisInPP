# Systematic Uncertainty Calculator

This takes the output of the Jet Folder and calculates a systematic uncertainty to be applied to the unfolded data.

---

`StJetFolder` was run several times to produce a set of systematic variations. The macros in this module read in those variations, calculated the corresponding systematic uncertainties, smooth them where necessary, and sum them to produce a final systematic uncertainty.

The systematic uncertainties are calculated using `CalculateDetectorSystematicError.C` and `CalculateUnfoldingSystematicError.C`, where the former is used to calculate the uncertainty associated with the TPC respnse and the latter used to calculate the uncertainty associated with the unfolding. All parameters and input/output files are set in the macros themselves.

Then the unfolding uncertainties are smoothed, and the detector and systematic uncertainties are summed in the macro `SumAndSmoothErrors.C`. The macro `MakeSummaryPlot.C` produces a plot of all 6 jet spectra (all 3 energy ranges of both $\pi^{0}$ and $\gamma$ triggers) with their systematic and statistical uncertainties for a given $R$.

Additional macros which may be helpful can be found in the `macros` directory.

```
root -b -q Calculate{Detector,Unfolding}SystematicError.C
root -b -q SumAndSmoothErrors.C
root -b -q MakeSummaryPlot.C
```
