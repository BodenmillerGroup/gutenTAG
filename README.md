# maldi-processing

This R package supports the data handing, processing, quality control and analysis of 
  targeted MALDI-imaging data. The data handling capabilities (import, pre-processing, export) build 
  heavily on the R package Cardinal, designed for general Mass Spectrometry Imaging data. The main 
  functionality of this package includes reliable data processing, extensive quality control and data
  analysis. 
  
  First, data is processing using a novel approach involving the construction of metapeaks
  to account for mass shift and the presence of isotopic peaks. Then image-level and spatial quality
  control statistics, such as computing Geary's C scores can be computed for every marker. Furthermore, 
  spatial analysis can be performed on the high-dimensional images, including constructing variograms 
  and a variety of clustering algorithms. The package also supports dimensionality reduction methods, 
  including PCA, NMF and UMAP. Additionally, the package supports image visualisation from cytomapper
  (and cytoviewer) through the construction of a CytoImageList object. As such, once the data is coerced
  to this format, many of the analysis techniques that are applicable to CytoImageList can be applied.

## Installation

For installing the development version of the package, run the following lines in R:

``` r

install.packages("devtools")
devtools::install_github("BodenmillerGroup/maldi-processing")

```

## Functionality

`maldi-processing` provides an end-to-end workflow for analysing targeted MALDI-imaging data. The main workhorse for
this are the pre-processing and processing steps. Auxilliary functionalities include quality control statistics as well as 
analysis techniques.

#### 1. Pre-processing

The `preprocess` function implements the `Cardinal::normalize`, `Cardinal::smoothSpectra` and `Cardinal::reduceBaseline` functions in series. The output is an `MSImagingExperiment` object, in-line with the Cardinal framework. Alternatively, you can perform pre-processing using the Cardinal functions directly. 
