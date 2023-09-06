## ---- echo=FALSE, results="hide"----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE, crop = NULL)
library(BiocStyle)

## ----library, echo=FALSE------------------------------------------------------
library(maldipackage)

## ----read in data-------------------------------------------------------------
imzML_path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
rawFile <- Cardinal::readMSIData(imzML_path)

## ----read in panel------------------------------------------------------------
panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
panel <- maldipackage::readPanel(panel_path)


## ----example data-------------------------------------------------------------
path <- system.file("extdata", package = "maldipackage")
list.files(path, recursive = TRUE)

## ----pre-processing-----------------------------------------------------------
pre <- preProcess(rawFile, cores = 8)

## ----Cardinal pre-processing--------------------------------------------------
pre <- rawFile %>%
    Cardinal::normalize(method = "tic") %>%
    smoothSignal(method = "gaussian", plot = FALSE) %>%
    reduceBaseline(method = "locmin") %>%
    process(BPPARAM = MulticoreParam(workers = 8))


## ----Peak detection-----------------------------------------------------------

detected <- maldipackage::peakDetection(pre, snr = 3, cores = 8)


## ----Metapeak generation------------------------------------------------------

metapeaks <- maldipackage::metapeakGeneration(detected)


## ----Generate final output----------------------------------------------------
processed <- maldipackage::getIntensityDF(metapeaks, refList = panel, pre = pre)


## ----coerce to CytoImageList object-------------------------------------------
library(cytomapper)
cil <- asCytoImageList(processed)


## ----single channel imaging---------------------------------------------------
cytomapper::plotPixels(image = cil, colour_by = "HLA.ABC")


## ----multi-channel imaging----------------------------------------------------
plotPixels(cil, colour_by = c("HLA.ABC", "CD98", "beta.actin"), 
           legend = list(colour_by.title.cex = 0.7), 
           scale_bar = list(length = 4, lwidth = 0.5, margin = c(1,1)),
           image_title = list(margin = c(1, 1)))


## ----cytoviewer---------------------------------------------------------------
library(cytoviewer)
cytoviewer(cil)

## ----pca, echo=TRUE-----------------------------------------------------------

my_pca <- maldipackage::pca(processed$IntensityDF, comp = 5, scree = TRUE)

for (i in 1:ncol(my_pca$x)){
  maldipackage::imageChannel(my_pca$x, coords = processed$SpatialCoords, channel_number = i)
}


## ----nmf----------------------------------------------------------------------
my_nmf <- maldipackage::nmf(processed$IntensityDF, comp = 5)

# visualise the factored dimensions
for (i in 1:ncol(my_pca$x)){
  maldipackage::imageChannel(my_nmf$w, coords = processed$SpatialCoords, channel_number = i)
}

## ----umap---------------------------------------------------------------------
my_umap <-computeUMAP(processed$IntensityDF)


## ----visualise single channel intensities on umap-----------------------------
# create data frame
umap_df <- cbind(my_umap, processed$IntensityDF)

# plot UMAP with channel intensity
library(ggplot2)
ggplot(umap_df) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = HLA.ABC), alpha = 0.6, size = 0.1) +
  scale_color_viridis_c() +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "HLA.DR",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "Intensity")

## ----community detection------------------------------------------------------
clusters <- louvain_cluster(processed$IntensityDF)

## ----visualise single channel intensities on spatial coords-------------------
# create data frame
cluster_df <- cbind(my_umap, clusters)

# plot UMAP with channel intensity
ggplot(cluster_df, aes(x = UMAP1, y = UMAP2, color = clusters)) +
  geom_point(size = 0.1) +
  scale_color_discrete(name = "Cluster Membership") +
  theme_minimal() +
  labs(title = "Cluster Memberships on Spatial Coordinates",
       x = "x",
       y = "y")


## ----visualise single channel intensities on spatial coordinates--------------
# create data frame
cluster_coords <- cbind(processed$SpatialCoords, clusters)

# ggplot spatial distribution of clusters
ggplot(cluster_coords, aes(x = x, y = y, color = clusters)) +
  geom_point(size = 0.1) +
  scale_color_discrete(name="Cluster Membership") +
  #scale_color_brewer(palette = "Set1", name = "Cluster Membership") + # alternative colour set for clustering
  theme_minimal() +
  labs(title = "Cluster Memberships on Spatial Coordinates",
       x = "X Coordinate",
       y = "Y Coordinate")

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

