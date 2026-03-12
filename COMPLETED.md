# Completed Tasks

## Version 0.12.1 — Branch: bioc-fixes

- make sure default parameter for watershed threshold is 0.01
- make sure processed$IntensityDF is ordered by m/z value
- make sure that CorrespondenceMatrix$expected_mz_location is numeric and not character
- untargeted metapeak columns should be assigned the name of that m/z value (as.character)
- Fix EBImage missing from DESCRIPTION and NAMESPACE — already present: EBImage in DESCRIPTION Imports, importFrom(EBImage,Image) and importFrom(EBImage,propagate) in NAMESPACE, and @importFrom tags in asCytoImageList.R and segmentPeakCounts.R
- Fix abind missing from DESCRIPTION and NAMESPACE — already present: abind in DESCRIPTION Imports, importFrom(abind,abind) in NAMESPACE, and @importFrom abind abind in asCytoImageList.R
- Fix DESCRIPTION: Author format, License, LazyData, BiocStyle placement — already correct; removed erroneous @import BiocStyle from preProcess.R and import(BiocStyle) from NAMESPACE
- Fix duplicate Cardinal import in NAMESPACE — no duplicate present; @import Cardinal was already absent from generateMetapeaks.R
