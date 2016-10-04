# p-checker: The one-for-all p-value analyzer.

This Shiny app code implements the **p-curve** (Simonsohn, Nelson, & Simmons, 2013; see http://www.p-curve.com), the **R-Index**, the **Test of Insufficient Variance, TIVA** (Schimmack, 2014; see http://www.r-index.org/), and checks whether p-values are reported correctly.

Some code is adapted from Uri Simonsohn (http://www.p-curve.com/Supplement/R/).

The currect stable app is at (http://shinyapps.org/apps/p-checker/).

## Known issues:
- TODO: Clearly separate the inference functions from UI functions
- TODO: Make TIVA computation robust against outliers
- TODO: Make the code prettier, add a lot of annotations