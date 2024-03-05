- make sure default parameter for watershed threshold is 0.01.
- make sure processed$IntensityDF is ordered by m/z value
- make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed
- untargeted metapeak columns should be assigned the name of that m/z value (as.character)
  
- processing output should contain everything needed for shinyTAG app
  - everything in 'processed' (df, coords, correspondence, filtered, untargeted)
  - everything in 'metapeaks' (counts, smooth counts, metapeaks, seeds, limit, propagation_selection)

- estimation criterion for what watershed threshold to use
