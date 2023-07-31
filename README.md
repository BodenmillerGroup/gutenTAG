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

To install the development version of the package, install via GitHub:

``` r

install.packages("devtools")
devtools::install_github("BodenmillerGroup/maldi-processing")

```

## Functionality

`maldi-processing` provides an end-to-end workflow for analysing targeted MALDI-imaging data. The main workhorse for
this are the pre-processing and processing steps. Auxilliary functionalities include quality control statistics as well as 
analysis techniques.

#### 1. Read in data

Three input files are required for targeted MALDI-imaging experiments:

- .imzML file (metadata)
- .ibd file (binary data, not required to read in, but should be in the same directory as the metadata file)
- .csv file (Panel)

Currently, there is no function provided in this package to read in the .imzML file. It is instead recommended to simply use the Cardinal `readMSIData` function for easy reading in of the data. 

```r

rawFile <- Cardinal::readMSIData(/path/to/.imzML)

```

The panel can be easily loaded into your R session using the `readPanel` function:

```r

panel <- readPanel(/path/to/panel.csv)

```

#### 2. Pre-processing

The `preprocess` function implements the `Cardinal::normalize`, `Cardinal::smoothSpectra` and `Cardinal::reduceBaseline` functions in series. The output is an `MSImagingExperiment` object, in-line with the Cardinal framework. Alternatively, you can perform pre-processing using the Cardinal functions directly. 

``` r

pre <- preprocess(rawFile, cores = 4)

```

#### 3. Processing

Full processing of the pre-processed MALDI-imaging data can be performed using the `peakDetection`, `metapeakGeneration` and `getIntensityDF` functions in series.

``` r

list_peaks <- peakDetection(pre)
metapeaks <- metapeakGeneration(list_peaks)
processed <- getIntensityDF(metapeaks, refList = panel, pre = pre)

```

The output of processing is a simple object containing 4 elements

- `IntensityDF` An intensity dataframe for all markers that were paired with a metapeak.
- `CorrespondenceMatrix` A targeted correspondence matrix detailing which metapeaks were associated to which mass tags.
- `SpatialCoords` The spatial coordinates for each pixel.
- `Untargeted` A list containing the untargeted intensity dataframe (all metapeaks) and the untargeted correspondence matrix. 


