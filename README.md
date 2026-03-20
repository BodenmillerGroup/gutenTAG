<img src="vignettes/imgs/gutenTAG_logo.png" align="right" alt="" width="100" />

# gutenTAG

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

gutenTAG is an R package for the data processing and quality control of targeted MALDI mass spectrometry imaging (MALDI-MSI) data. It builds on [Cardinal](https://bioconductor.org/packages/Cardinal/), a Bioconductor framework for MSI data, and extends it with targeted workflows specific to MALDI-imaging experiments.

The core processing approach uses a novel metapeak construction strategy to account for mass shift and isotopic peaks in the spectra. Downstream, the package provides image-level and spatial quality control statistics, and supports conversion to standard Bioconductor formats such as `SpatialExperiment`, `CytoImageList`, and `AnnData` — enabling dimensionality reduction, clustering, and visualisation using the broader Bioconductor ecosystem. Image visualisation is supported through [cytomapper](https://bioconductor.org/packages/cytomapper/) and [cytoviewer](https://bioconductor.org/packages/cytoviewer/), a Shiny-based interactive viewer for CytoImageList objects that enables interactive exploration of multichannel images.

## Installation

### Bioconductor (recommended)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("gutenTAG")
```

> **Note:** gutenTAG depends on several Bioconductor packages that cannot be installed via `install.packages()`. Always use `BiocManager::install()`.

### Development version (GitHub)

```r
install.packages("devtools")
devtools::install_github("BodenmillerGroup/gutenTAG")
```

## Functionality

`gutenTAG` provides an end-to-end workflow for processing targeted MALDI-imaging data.

### 1. Read in data

Three input files are required for targeted MALDI-imaging experiments:

- `.imzML` file (metadata)
- `.ibd` file (binary spectral data — required, but does not need to be passed explicitly; Cardinal locates it automatically from the `.imzML` path)
- `.csv` file (panel of target masses)

Example raw data is stored in the `inst/extdata` directory of this repository. It is recommended to use the Cardinal `readMSIData` function to read in the `.imzML` data:

```r
path <- system.file("extdata/Example_data.imzML", package = "gutenTAG")
rawFile <- readMSIData(path)
```

The panel can be loaded using `readPanel`:

```r
panel_path <- system.file("extdata/ref_list.csv", package = "gutenTAG")
panel <- readPanel(path = panel_path)
```

### 2. Pre-processing

The `preProcess` function applies `Cardinal::normalize`, `Cardinal::smooth`, and `Cardinal::reduceBaseline` in series. The output is an `MSImagingExperiment` object. Alternatively, you can call the Cardinal functions directly.

```r
pre <- preProcess(rawFile)
```

**Parallelisation:** `preProcess` and `peakDetection` both accept a `BPPARAM` argument for parallel execution via `BiocParallel`. Register a backend once at the start of your session and all subsequent calls will use it automatically:

```r
# Linux / macOS
BiocParallel::register(BiocParallel::MulticoreParam(workers = 4))

# Windows
BiocParallel::register(BiocParallel::SnowParam(workers = 4))
```

See the package vignette for full details.

### 3. Processing

Full processing is performed using `peakDetection`, `generateMetapeaks`, and `assignMetapeaks` in series.

`peakDetection` identifies peaks in the pre-processed spectra (a wrapper for `Cardinal::peakPick`):

```r
list_peaks <- peakDetection(pre)
```

`generateMetapeaks` groups peaks into metapeaks — single composite peaks that correspond to the same molecular species, compensating for technical mass shift and isotopic peaks:

```r
metapeaks <- generateMetapeaks(list_peaks)
```

`assignMetapeaks` matches metapeaks to the target panel and generates the final intensity dataframe:

```r
processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)
```

The output is a list containing 6 elements:

- `IntensityDF` — Intensity dataframe for all markers paired with a metapeak.
- `CorrespondenceMatrix` — Targeted correspondence matrix detailing which metapeaks were associated with which mass tags.
- `SpatialCoords` — Spatial coordinates for each pixel.
- `SummarySpectra` — Summary spectral information across pixels.
- `FilteredDF` — Intensity dataframe with background pixels filtered out. Useful for clustering analysis.
- `AllMetapeaks` — A list with two elements:
  - `AllMetapeaksIntensity` — Intensity dataframe for all detected metapeaks (targeted and untargeted).
  - `AllMetapeaksCorrespondence` — Correspondence matrix for all detected metapeaks.

### 4. Quality control

A detailed quality control report can be found under `inst/qc_report.Rmd`. This runs and explains all exported QC functions from gutenTAG.

### 5. Conversion to Bioconductor formats

The processed output can be converted to standard Bioconductor objects for downstream analysis via the respective ecosystem tools:

```r
# CytoImageList for image visualisation with cytomapper/cytoviewer
cil <- asCytoImageList(processed)

# SpatialExperiment for spatial omics workflows
spe <- asSpatialExperiment(processed)

# AnnData for interoperability with Python/scanpy
ann <- asAnnData(processed)
```

## Citation

If you use gutenTAG in your work, please cite it using our pre-print:

> Zhang M, Abbey J, (2025). *Spatial proteomics with 100+ markers with High Multiplexed MALDI-IHC* bioRxiv (2025). https://www.biorxiv.org/content/10.1101/2025.04.18.649415v1.abstract

```bibtex
@article{,
  title  = {Spatial proteomics with 100+ markers with High Multiplexed MALDI-IHC},
  author = {Mengze Zhang and John Abbey},
  year   = {2025},
  doi = {https://doi.org/10.1101/2025.04.18.649415},
  url = {https://www.biorxiv.org/content/10.1101/2025.04.18.649415v1.abstract},
  journal   = {bioRxiv}
}
```

We will soon have a pre-print out for gutenTAG!

Please also acknowledge the Cardinal framework:

> Bemis et al. (2023). Cardinal v3: a versatile open-source software for mass spectrometry imaging analysis. *Bioinformatics*.


## License

MIT © BodenmillerGroup. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome. Please open an issue at the [issue tracker](https://github.com/BodenmillerGroup/gutenTAG/issues) before submitting a pull request.
