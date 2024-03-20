- ~~make sure default parameter for watershed threshold is 0.01.~~

- ~~make sure processed$IntensityDF is ordered by m/z value~~
- **make sure that untargeted metapeaks that are weirdly actually assigned to markers are removed**
- ~~make sure that CorrespondenceMatrix$expected_mz_location is numeric and not character~~

  
- ~~untargeted metapeak columns should be assigned the name of that m/z value (as.character)~~

  
- estimation criterion for what watershed threshold to use


- refactor assignMetapeaks into 3 functions
  - one function `assignMetapeaks` which creates a new dataframe where each row contains the marker and the metapeak that is assigned to it - should be transparent
  - one for computing the stats and updating the dataframe
  - one for generating the images
 
- Output changes:
  - Correspondence matrix:
    - one table containing current information + metapeak information. Missing markers (zero columns) get NA values.
