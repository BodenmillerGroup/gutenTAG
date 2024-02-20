<img src="vignettes/gutenTAG_logo.png" align="right" alt="" width="200" />


# gutenTAG

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
devtools::install_github("BodenmillerGroup/maldi-imaging")

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

Example raw data is stored in the `inst/extdata` directory of this repository. These files can be retrieved using the two chunks of code below.
It is recommended to use the Cardinal `readMSIData` function to easily read in the .imzML data. 

```r
path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
rawFile <- readMSIData(path)
```

The panel can be conveniently loaded into your R session using the `readPanel` function:

```r
panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
panel <- readPanel(path = panel_path)

```

#### 2. Pre-processing

The `preProcess` function implements the `Cardinal::normalize`, `Cardinal::smoothSpectra` and `Cardinal::reduceBaseline` functions in series. The output is an `MSImagingExperiment` object, in-line with the Cardinal framework. Alternatively, you can perform pre-processing using the above Cardinal functions directly. 

``` r

pre <- preProcess(rawFile, cores = 2)

```

#### 3. Processing

Full processing of the pre-processed MALDI-imaging data can be performed using the `peakDetection`, `metapeakGeneration` and `assignMetapeaks` functions in series.

The first step of processing is to identify peaks in the pre-processed spectra. This function is essentially a wrapper for `Cardinal::peakPick`.

``` r

list_peaks <- peakDetection(pre)

```

Once `peakDetection` has been run, `generateMetapeaks` can be run. Metapeaks are single peaks that correspond to clusters of peaks that contain information the same molecular species. The function used to generate these metapeaks, `generateMetapeaks`, is implemented to compensate for the loss of signal through technical shift and isotopic peaks in the spectra. 

``` r

metapeaks <- generateMetapeaks(list_peaks)

```
`assignMetapeaks` is the final function needed to complete the processing workflow. It is responsible for generating the targeted intensity dataframe.

``` r

processed <- assignMetapeaks(metapeaks, refList = panel, pre = pre)

```

The output of processing is a simple object containing 5 elements

- `IntensityDF` An intensity dataframe for all markers that were paired with a metapeak.
- `CorrespondenceMatrix` A targeted correspondence matrix detailing which metapeaks were associated to which mass tags.
- `SpatialCoords` The spatial coordinates for each pixel.
- `FilteredDF` An intensity dataframe for all markers assigned to a metapeak, with background pixels filtered out. Useful for clustering analysis.
- `Untargeted` A list containing the untargeted intensity dataframe (all metapeaks) and the untargeted correspondence matrix. 
