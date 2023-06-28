# 1046 Normalisation

#### i). Loading packages and functions ####
{cat("Loading of the libraries....  ")
  suppressMessages(library(Cardinal))
  suppressMessages(library(fftw))
  suppressMessages(library(N2R))
  suppressMessages(library(Matrix))
  suppressMessages(library(igraph))
  suppressMessages(library(splus2R))
  suppressMessages(library(EBImage))
  suppressMessages(library(imager))
  suppressMessages(library(FactoMineR))
  suppressMessages(library(hexbin))
  suppressMessages(library(matrixStats))
  suppressMessages(library(mclust))
  suppressMessages(library(dplyr))
  suppressMessages(library(matter))
  cat("... done ! \n")
}

# Load in functions from helper script
source("/mnt/msi_volume/Rscripts/maldi-processing/code/helper_functions.R")

#### ii). Data and annotation loading ####
args <- commandArgs(trailingOnly = TRUE)

cat("Reading the .imzml and .ibd file... ")
Path_to_imzml_file <- args[1]
Path_to_imzml_file <- "/mnt/msi_volume/experiments/JuneforJohn/colon_13_plex/menzha_20230528_53436_colon-50um.imzML"
rawFile <- readMSIData(Path_to_imzml_file)
cat("done ! \n")

# subset for small test image
rawFile <- subsetPixels(rawFile, x<=50, y<= 50)


cat("Subsetting if needed... ")
rawFile <- subsetPixels(rawFile, x<=200, y<= 200)
cat("done ! \n")

cat("Assigning experiment and sample name... ")
Sample_name <- strsplit(Path_to_imzml_file,split = "/",fixed = TRUE)[[1]]
# initialise correct folder as the first string in the sequence
h <- 1
correctFolder <- Sample_name[h]
# while the correct folder is not "experiments", iterate through each string in the sequence until it is correct
if(correctFolder == "experiments"){
  print("Correct on first iteration")
  print(paste("Current folder is ", correctFolder, sep = ""))
  h = h+1
}else{while(correctFolder != "experiments"){
  correctFolder <- Sample_name[h]
  print(paste("Current folder is ", correctFolder, sep = ""))
  h = h+1
  if (correctFolder == "experiments"){
    print(paste("The current folder is ", correctFolder, ".", sep = ""))
  }
}}

Experiment_name <- Sample_name[h]
print(paste("The experiment is ", Experiment_name, ".", sep = ""))
Sample_name <- Sample_name[length(Sample_name)]
Sample_name <- strsplit(Sample_name,split = ".",fixed = TRUE)[[1]]
Sample_name <- Sample_name[1]
print(paste("The sample is ", Sample_name, ".", sep = ""))
cat("done ! \n")

#Path_to_peak_annotation = args[2]
Path_to_peak_annotation <- "/mnt/msi_volume/panels/JuneforJohn/13plex_colon_kidney.csv"
cat("Reading the csv annotation file... ")
peakAnnotation <- read.delim(Path_to_peak_annotation, sep=",", header=TRUE, col.names = c("Name", "FeatureMass"))
# sanity check: columns are correctly named
if (is.numeric(peakAnnotation$Name) == TRUE){
  colnames(peakAnnotation) = c("FeatureMass","Name")
}
peakAnnotation <- peakAnnotation[-3]
peakAnnotation <- dplyr::arrange(peakAnnotation, peakAnnotation$FeatureMass)
cat("done ! \n")

# remove spaces and problematic characters from marker names
for(a in seq_along(peakAnnotation$Name)){
  if (grepl("+", peakAnnotation$Name[a], fixed=TRUE)){
    new_string <- gsub(" ", "", peakAnnotation$Name[a])
    peakAnnotation$Name[a] <- new_string
  }
  else{
    # split strings with spaces into list with individual strings as elements
    new_string <- gsub(" ", "", peakAnnotation$Name[a])
    # replace fullstops with underscores
    new_string <- gsub("-", ".", new_string, fixed = TRUE)
    # if there is a dash at the end of the name, remove it
    if(endsWith(new_string, "-")){
      new_string <- substr(new_string,1, nchar(new_string)-1)
    }
    # replace names in peakAnnotation
    peakAnnotation$Name[a] <- new_string
  }
  # remove slashes
  if (grepl("/", peakAnnotation$Name[a], fixed=TRUE)){
    new_string <- gsub("/", "", peakAnnotation$Name[a])
    peakAnnotation$Name[a] <- new_string
  }


}


#### iii). Pre-Processing ####

cat("Running normalisation and baseline correction... ")
peakPre <- rawFile %>%
  Cardinal::normalize(method = "tic") %>%
  smoothSignal(method = "gaussian", plot=FALSE) %>%
  reduceBaseline(method="locmin") %>%
  process(BPPARAM = MulticoreParam(workers = 4))
cat("...done ! \n")

image_1046 <- iData(peakPre_1046)
index_1046 <- features(peakPre, mz = 1046)

# Find the feature corresponding to 1046
index_1046 <- features(peakPre, mz = 1046)

# Vector of all pixel intensities at mz = 1046
pixel_intensities_1046 <- image[index_1046, ]

# Sanity check
length(image[index_1046, ])

#
range(pixel_intensities_1046)
max(pixel_intensities_1046)


# Test: average all intensity values from 1046 to 1047 and normalise by this

# Find indices of features between 1046 and 1047
index_in_range <- features(peakPre, 1046 < mz & mz < 1047)

# Find features corresponding to indices
pixel_intensity_range <- image[index_in_range, ]

# Average the intensities of pixels at these intensities
avg_intensity_1046 <- colSums(pixel_intensity_range)/nrow(pixel_intensity_range)
# sanity check
length(avg_intensity_1046)

# Find highest intensity value
max_1046 <- max(avg_intensity_1046)

# Normalise all pixels by this value
View(image)
image_1046_norm <- image/max_1046

# Function for 1046 normalisation ####
referenceNorm <- function(x){

  image <- iData(x)
  index_in_range <- features(x, 1046 < mz & mz < 1047)
  pixel_intensity_range <- image[index_in_range, ]
  avg_intensity_1046 <- colSums(pixel_intensity_range) / nrow(pixel_intensity_range)
  max_1046 <- max(avg_intensity_1046)
  image_1046_norm <- image / max_1046

  return(image_1046_norm)

}

peakPre_norm <- referenceNorm(peakPre)

# Assign normalised to an MSImagingExperiment object
peakPre_norm <- MSContinuousImagingSpectraList(peakPre_norm)
fdata <- featureData(peakPre)
pdata <- pixelData(peakPre)

peakPre_1046 <- MSImagingExperiment(imageData = peakPre_norm,
                    featureData = fdata,
                    pixelData = pdata)

Cardinal::image(peakPre_1046, mz = 1000)
Cardinal::image(peakPre, mz = 1000)


