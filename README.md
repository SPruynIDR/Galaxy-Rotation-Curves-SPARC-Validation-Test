# Galaxy Rotation Curves: Complete SPARC Validation Test

Can the unexplained motion in galaxy rotation curves be connected to a measurable feature of ordinary matter rather than requiring a separately fitted dark-matter halo for every galaxy?

This repository provides a complete, runnable test of that question using 171 SPARC galaxies and 3,344 rotation-curve measurements.

The test evaluates a single global relation involving baryonic mass and the measured neutral-hydrogen boundary of each galaxy.

The complete notebook runs:

- The fixed global boundary relation
- The baryons-only baseline
- The MOND/RAR benchmark
- Individually fitted NFW dark-matter halos
- RMS and R² comparisons
- AIC and BIC model comparisons
- A 10,000-shuffle control testing whether the real galaxy-to-hydrogen-boundary pairings outperform randomized assignments

The global relation achieves an RMS error of approximately 33.68 km/s.

The tested MOND/RAR benchmark achieves approximately 41.50 km/s.

Individually fitted NFW halos achieve a better raw fit of approximately 25.95 km/s, but use 342 fitted parameters across the 171 galaxies.

The global boundary relation uses three empirical parameters across the entire sample.

AIC favors the individually fitted NFW halos.

BIC favors the three-parameter global relation by approximately 1,008 BIC points.

In the 10,000-shuffle control, none of the randomized galaxy-to-hydrogen-boundary assignments outperform the real pairings.

## ▶️ Run the complete test now

**[CLICK HERE TO OPEN THE COMPLETE TEST IN GOOGLE COLAB](https://colab.research.google.com/github/SPruynIDR/Galaxy-Rotation-Curves-SPARC-Validation-Test/blob/main/SPARC_COMPLETE_ALL_TESTS.ipynb)**

No installation is required.

Open the notebook and select:

**Runtime → Run all**

The notebook runs the complete analysis and displays the results.

## What is being claimed?

The claim is specific and testable:

The unexplained structure in these galaxy rotation curves appears to contain a statistically significant relationship with the measured neutral-hydrogen boundary.

A single low-parameter global relation using that boundary outperforms the tested MOND/RAR benchmark and competes with individually fitted NFW dark-matter halos strongly enough to be favored by BIC.

This repository is provided so the result can be independently checked rather than simply accepted.

Run the test.

Inspect the code.

Check the statistics.

Try to find where it fails.
