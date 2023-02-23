# MALDI pre-processing, processing, and QC script. Written by Pierre and John

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
cat("... done ! \n")
}

source("mnt/msi_volume/Rscripts/helper_functions.R")

#### ii). Data and annotation loading #### 
args <- commandArgs(trailingOnly = T)

cat("Reading the .imzml and .ibd file... ")
Path_to_imzml_file = args[1]
#Path_to_imzml_file = "mnt/msi_volume/experiments/glycan_experiment_01/menzha_20230210_glycantest3-breastcancer-20um.imzML"
Path_to_imzml_file = "mnt/msi_volume/experiments/QC-report-ImzML/42251-20um-breastcancer.imzML"
rawFile = readMSIData(Path_to_imzml_file)
cat("done ! \n")


cat("Subsetting if needed... ")
rawFile <- subsetPixels(rawFile, x<=200, y<= 250)
cat("done ! \n")

cat("Assigning experiment and sample name... ")
Sample_name = strsplit(Path_to_imzml_file,split = "/",fixed = T)[[1]]
# initialise correct folder as the first string in the sequence
h=1
correctFolder = Sample_name[h]
# while the correct folder is not "experiments", iterate through each string in the sequence until it is correct
if(correctFolder == "experiments"){
  print("Correct on first iteration")
  print(paste("Current folder is ", correctFolder, sep = ""))
  h = h+1
}else{while(correctFolder != "experiments"){
  correctFolder = Sample_name[h]
  print(paste("Current folder is ", correctFolder, sep = ""))
  h = h+1
  if (correctFolder == "experiments"){
    print(paste("The current folder is ", correctFolder, ".", sep = ""))
  }
}}

Experiment_name = Sample_name[h]
print(paste("The experiment is ", Experiment_name, ".", sep = ""))
Sample_name = Sample_name[length(Sample_name)]
Sample_name = strsplit(Sample_name,split = ".",fixed = T)[[1]]
Sample_name = Sample_name[1]
print(paste("The sample is ", Sample_name, ".", sep = ""))
cat("done ! \n")

Path_to_peak_annotation = args[2]
Path_to_peak_annotation = "mnt/msi_volume/panels/SciLS_library/bc_20um.csv"
cat("Reading the csv annotation file... ")
peakAnnotation = read.delim(Path_to_peak_annotation, sep=",", header=T, col.names = c("Name", "FeatureMass"))
# sanity check: columns are correctly named
if (is.numeric(peakAnnotation$Name) == T){
  colnames(peakAnnotation) = c("FeatureMass","Name")
}
peakAnnotation = peakAnnotation[-3]
peakAnnotation = dplyr::arrange(peakAnnotation, peakAnnotation$FeatureMass)
cat("done ! \n")

# remove spaces and problematic characters from marker names
for(a in seq_along(peakAnnotation$Name)){
  if (grepl("+", peakAnnotation$Name[a], fixed=T)){
    new_string = gsub(" ", "", peakAnnotation$Name[a])
    peakAnnotation$Name[a] <- new_string
  }
  else{
    # split strings with spaces into list with individual strings as elements
    new_string = gsub(" ", "-", peakAnnotation$Name[a])
    # replace fullstops with underscores
    new_string = gsub(".", "-", new_string, fixed = T)
    # if there is a dash at the end of the name, remove it
    if(endsWith(new_string, "-")){
      new_string = substr(new_string,1, nchar(new_string)-1)
    }
    # replace names in peakAnnotation
    peakAnnotation$Name[a] <- new_string
  }
  # remove slashes
  if (grepl("/", peakAnnotation$Name[a], fixed=T)){
    new_string = gsub("/", "", peakAnnotation$Name[a])
    peakAnnotation$Name[a] <- new_string
  }
  

}
# shift reference list masses by 0.5mz
#peakAnnotation$FeatureMass = peakAnnotation$FeatureMass +0.5

Path_to_parameter_file = args[3]
Path_to_parameter_file = "mnt/msi_volume/parameter_files/test_param_file.txt"

# load in parameter file
cat("Loading parameter file... ")
Parameters = suppressWarnings(read.table(Path_to_parameter_file,header = F,sep = "="))
Parameters_names = Parameters$V1
Parameters_values = Parameters$V2
names(Parameters_values) = Parameters_names
# define parameters
N_cores = as.numeric(Parameters_values[1])
peak_mz_tolerance = as.numeric(Parameters_values[2])
freqMin = as.numeric(Parameters_values[3])
cat("done ! \n")


#### iii). Pre-Processing #### 

cat("Running normalisation and baseline correction... ")
peakPre <- rawFile %>% 
  Cardinal::normalize(method = "tic") %>%
  smoothSignal(method = "gaussian", plot=FALSE) %>%
  reduceBaseline(method="locmin") %>%
  process()
cat("...done ! \n")

cat("Managing memory... ")
print(mcols(processingData(peakPre))[,-1])
rm(rawFile)
gc()
cat("...done ! \n")


#### iv). Processing (peak detection and peak aggregation) #### 

cat("Extracting pre-processed intensity data... ")
# extract pre-processed intensity dataframe
raw_intensity = peakPre@imageData$data$intensity
cat("...done ! \n")

cat("Defining important variables and constants... ")
# set up some important constants
range_peaks = range(peakPre@featureData@mz)
N_features = length(peakPre@featureData@mz)
N_pixels_total = dim(peakPre)[2]
Threshold_Detection = N_pixels_total*0.01
# all m/z values 
mz_vector = as.data.frame(peakPre@featureData)
Location_pixels = as.data.frame(pData(peakPre))
Location_pixels = Location_pixels[,c(2,3)]
cat("...done ! \n")

cat("Applying peak detection... ")
# apply MAD peak picking function to all pixels
List_peaks = apply(raw_intensity, MARGIN = 2, FUN = function(x) {peakPick.mad(x,SNR = 10,window = 400)})
cat("...done ! \n")

# number of peaks per pixel
N_peak_per_pixel = unlist(lapply(List_peaks,FUN = length))

cat("Getting peak counts ... ")
# re-scale peaks to m/z scale
List_peaks_rescaled = mz_vector[unlist(List_peaks),1]
# only unique peaks
List_peaks_unique = unique(unlist(List_peaks))
# rescale unique peaks to m/z scale
List_peaks_unique_rescaled = mz_vector$mz[List_peaks_unique]
# order rescales unique peaks
List_peaks_unique_rescaled_ordered = List_peaks_unique_rescaled[order(List_peaks_unique_rescaled)]
# how many pixels peaks are present in
Freq_peak_rescaled_ordered = table(factor(List_peaks_rescaled,levels = mz_vector$mz))
Freq_peak_rescaled_ordered = as.numeric(Freq_peak_rescaled_ordered)
cat("...done ! \n")

# We smooth with a Gaussian curve with a sigma of 1mz so that isotopic peaks are merged

cat("Smoothing spectra to merge isotopic peaks... ")
# mz index ratio: range of mz values divided by number of features
mz_index_ratio = (max(mz_vector$mz)-min(mz_vector$mz))/N_features
# inverse of this ratio
sigma_smoothing_isotopic = 1/mz_index_ratio
# number of mz values centered around zero
list_values = seq(-N_features/2, N_features/2, length.out = N_features)
Gaussian_vector = exp(-list_values^2/sigma_smoothing_isotopic^2)
Freq_peak_rescaled_ordered_smooth = stats::convolve(y = Gaussian_vector,x = Freq_peak_rescaled_ordered, conj = TRUE)
# counts of rescaled ordered smooth peaks
Freq_peak_rescaled_ordered_smooth = Freq_peak_rescaled_ordered_smooth[c((N_features/2):N_features,1:(N_features/2-1))]
cat("...done ! \n")

cat("Setting pixels below detection limit to zero... ")
# set values under the pixel detection threshold to zero
Freq_peak_smoothed_thresholded = Freq_peak_rescaled_ordered_smooth
Freq_peak_smoothed_thresholded[Freq_peak_smoothed_thresholded<Threshold_Detection] = 0
cat("...done ! \n")

cat("Finding local maxima... ")
# We want to find local maxima within a range of 3mz 
span_local_maxima = 3/mz_index_ratio
# round up to nearest odd value
Local_maxima_peak_freq = base::which(find_peaks(Freq_peak_smoothed_thresholded,span = span_local_maxima))
# flattened 
Local_maxima_peak_freq_reshaped = matrix(rep(0,length(Freq_peak_rescaled_ordered)),ncol = 1)
Local_maxima_peak_freq_reshaped[Local_maxima_peak_freq]=1:length(Local_maxima_peak_freq)
cat("...done ! \n")

cat("Finding peak maxima, width and center per peak... ")
Freq_peak_smoothed_thresholded = matrix(Freq_peak_smoothed_thresholded,ncol = 1)
Propagation_selection = propagate(Freq_peak_smoothed_thresholded, seeds = Local_maxima_peak_freq_reshaped, mask = Freq_peak_smoothed_thresholded>Threshold_Detection)

# watershed_peak = watershed(as.cimg(Local_maxima_peak_freq_reshaped),as.cimg(matrix(Freq_peak_rescaled_ordered_thresholded,ncol = 1)))
# watershed_peak[Freq_peak_rescaled_ordered_thresholded==0]=0

Peak_location = aggregate.data.frame(data.frame(Location = mz_vector$mz,Freq = Freq_peak_rescaled_ordered),
                                     by=list(as.numeric(Propagation_selection)),FUN = function(x) {sum(x[[1]]*x[[2]])/sum(x[[2]])})
Peak_location_centered = c()
Peak_location_max = c()
Peak_delimitation = c()
for (k in 1:max(Propagation_selection)) {
  # Looking for weigted mean 
  Location_centered_temp = weighted.mean(x=mz_vector$mz[as.numeric(Propagation_selection)==k],
                                         w =Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k] )
  # vector of centered peak locations
  Peak_location_centered = c(Peak_location_centered,Location_centered_temp)
  
  Location_max_temp  = (mz_vector$mz[as.numeric(Propagation_selection)==k])[which.max(Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k])]
  # vevtor of peak maximimum locations
  Peak_location_max = c(Peak_location_max,Location_max_temp)
  
  Cumsum_freq = (cumsum(Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k])/sum(Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k]))
  
  Beginning_peak = max(which(Cumsum_freq < 0.01))
  End_peak = min(which(Cumsum_freq > 0.99))
  if (is.infinite(Beginning_peak)) {
    Beginning_peak = 1
  }
  # beginning and end of peak location
  Peak_delimitation = rbind(Peak_delimitation,c(mz_vector$mz[as.numeric(Propagation_selection)==k][Beginning_peak],mz_vector$mz[as.numeric(Propagation_selection)==k][End_peak]))
}
# width of each peak 
Peak_width = Peak_delimitation[,2]-Peak_delimitation[,1]
# width of each peak on m/z scale 
Peak_width = Peak_width*mz_index_ratio
cat("...done ! \n")

cat("Modifying peakAnnotation df... ")
# create completed peakAnnotation table for adducts
peakAnnotation_extended = peakAnnotation

# create vectors for each adduct ion
FeatureMassHydrogen = peakAnnotation_extended$FeatureMass
FeatureMassSodium = (FeatureMassHydrogen - 1) + 23
FeatureMassAmmonium = (FeatureMassHydrogen - 1) + 18
peakAnnotation_extended = cbind(peakAnnotation_extended, FeatureMassSodium, FeatureMassAmmonium)
# bind the feature masses to the feature names and adduct ions, then bind all features together to get flattened list
peakAnnotation_extended = rbind(cbind(peakAnnotation_extended$Name, peakAnnotation_extended$FeatureMass),
                                cbind(paste(peakAnnotation_extended$Name,"Sodium", sep = "_"), peakAnnotation_extended$FeatureMassSodium),
                                cbind(paste(peakAnnotation_extended$Name,"Ammonium", sep = "_"), peakAnnotation_extended$FeatureMassAmmonium))
# convert to dataframe
peakAnnotation_extended = as.data.frame(peakAnnotation_extended)
colnames(peakAnnotation_extended) = c("FeatureAnnotation","Location")
# convert peak location to numeric
peakAnnotation_extended$Location = as.numeric(peakAnnotation_extended$Location)
cat("...done ! \n")


cat("Computing cross matrix nearest neighbours between detected peaks and metapeaks... ")
# cross compare total reference list with metapeaks (find closest peak across lists)
Mapping_meta = crossKnn(mA = matrix(Peak_location_max, ncol = 1),
                        mB= matrix(peakAnnotation_extended$Location, ncol = 1), k = 10, indexType = "L2")
mz_threshold_association = 1
# remove all mappings below the m/z distance association threshold
Mapping_meta[Mapping_meta>mz_threshold_association]=0
Mapping_meta_cleaned = apply(as.matrix(Mapping_meta),MARGIN = 2,FUN = which_max_modified)
# table comparing the metapeak location with the detected peak location and the associated peak name
cat("...done ! \n")


cat("Assigning marker labels to nearest neighbours... ")
Associated_marker = c()
for (k in 1:length(Peak_location_max)) {
  selected_mz_location = range(mz_vector$mz[as.numeric(Propagation_selection)==k])
  Annotation_assigned = peakAnnotation_extended$FeatureAnnotation[peakAnnotation_extended$Location>selected_mz_location[1] & peakAnnotation_extended$Location<selected_mz_location[2]]
  if (length(Annotation_assigned)==0){
    Annotation_assigned = NA
  }
  # if there are more than one possible annotations, annotate it as such.
  ## what to do with these?? 
  if (length(Annotation_assigned)>1){
    Annotation_assigned = paste(Annotation_assigned,collapse="__OR__")
  }
  Associated_marker = c(Associated_marker,Annotation_assigned)
}
cat("...done ! \n")

cat("Creating Correspondence matrix... ")
Correspondence_matrix = data.frame(mz_location = Peak_location_max,
                                   Expected_mz_location = peakAnnotation_extended$Location[Mapping_meta_cleaned],
                                   Annotation_peaks = peakAnnotation_extended$FeatureAnnotation[Mapping_meta_cleaned])
rownames(peakAnnotation_extended) = peakAnnotation_extended$FeatureAnnotation

#Correspondence_matrix = data.frame(mz_location = Peak_location_max,
 #                                  Expected_mz_location = peakAnnotation_extended$Location[Associated_marker],
  #                                 Annotation_peaks = Associated_marker)
cat("...done ! \n")


cat("Creating Final intensity matrix... ")
Final_intensity_matrix = c()
for (k in 1:nrow(Correspondence_matrix)) {
  # which mz locations are inside the peak (boolean)
  selected_mz_location = mz_vector$mz >= Peak_delimitation[k,1] & mz_vector$mz <= Peak_delimitation[k,2]
  # those mz values
  mz_vector$mz[which(selected_mz_location == T)]
  if (sum(selected_mz_location)>1) {
    intensity_temp = colSums(raw_intensity[selected_mz_location,])
  }
  # if there is only one mz location, the intensity list is only the intensity values of the one mz location
  if (sum(selected_mz_location)==1) {
    intensity_temp = (raw_intensity[selected_mz_location,])
  }
  #intensity_temp = Correspondence_matrix$Annotation_peaks[k]
  Final_intensity_matrix = cbind(Final_intensity_matrix,intensity_temp)
  
}
cat("...done ! \n")

### DEV: GMM analysis on Final intensity matrix

# which of the metapeaks in the final intensity matrix are annotated peaks
Final_intensity_matrix_targeted = Final_intensity_matrix[,!is.na(Correspondence_matrix$Annotation_peaks)]
colnames(Final_intensity_matrix_targeted) = Correspondence_matrix$Annotation_peaks[!is.na(Correspondence_matrix$Annotation_peaks)]
Final_intensity_matrix_targeted = as.data.frame(Final_intensity_matrix_targeted)

cat("... Processing complete ! \n")


#### v). Summary statistics ####

## a). Geary's C

cat("Starting Geary's C score computation... \n ")

cat("Creating spatial weight matrix... ")
# create spatial weighted matrix (adjacency matrix)
spatial_weight_matrix = N2R::Knn(as.matrix(Location_pixels),k = 4,verbose = T,indexType = "L2")
# make sparse matrix
spatial_weight_matrix = as(spatial_weight_matrix,"dgCMatrix")
cat("...done ! \n")

cat("Computing Geary's C score... ")
# Create vector of Geary's C score for each peak
geary_vector = c()
for (i in 1:ncol(Final_intensity_matrix)){
  x = matrix(Final_intensity_matrix[,i],ncol = 1)
  x_squared = x^2
  product_temp_1 = 2*sum(spatial_weight_matrix%*%x_squared)
  product_temp_2 = 2*t(x)%*%spatial_weight_matrix%*%x
  N = length(x)
  W = sum(spatial_weight_matrix)
  Var_X = var(x)*N
  # compure Geary's C score
  Geary_C = as.numeric((N-1)*(product_temp_1-product_temp_2)/(2*Var_X*W))
  # add geary value to geary vector
  geary_vector = c(geary_vector, Geary_C)
}
cat("...done ! \n")

# add geary vector to Correspondance matrix
Correspondence_matrix$Geary = geary_vector
cat("... Finished Geary's C ! \n")

# b). Summary statistics (mean, variance, corrected standard deviation)

cat("Calculating standard statistics for final intensities... ")
### Normalised variance compared to mean
Mean_signal = colMeans(Final_intensity_matrix)
Sd_signal = apply(Final_intensity_matrix,MARGIN = 2,FUN = sd)
Corrected_sd = lm(log(Sd_signal)~log(Mean_signal))
# add peak standard deviation and peak width to correspondance matrix
Correspondence_matrix$Corrected_sd = Corrected_sd$residuals
Correspondence_matrix$PeakWidth = Peak_width
# cumulative signal of each peak in every pixel 
Total_signal = colSums(raw_intensity)
# correlation between the marker intensity and total signal intensity
Correlation_total_signal = apply(Final_intensity_matrix,MARGIN =2 ,FUN = function(x) {cor(log(x+1),log(1+Total_signal))})
Correspondence_matrix$Correlation_total_signal = Correlation_total_signal^2
cat("...done ! \n")

# c). Row statistics for histogram plot
meanColumn <- rowMeans(Final_intensity_matrix)
ticColumn <- rowSums(Final_intensity_matrix)
sdColumn <- rowSds(as.matrix(Final_intensity_matrix[sapply(Final_intensity_matrix, is.numeric)]))
varColumn <- rowSds(as.matrix(Final_intensity_matrix[sapply(Final_intensity_matrix, is.numeric)]))^2
maxColumn <- rowMaxs(as.matrix(Final_intensity_matrix[sapply(Final_intensity_matrix, is.numeric)]))
minColumn <- rowMins(as.matrix(Final_intensity_matrix[sapply(Final_intensity_matrix, is.numeric)]))

# create dataframe
pixelsSummarised = as.data.frame(cbind(ticColumn, meanColumn))

ticSummary = pixelsSummarised$ticColumn
# convert zeros to NA and then remove NA
ticSummary[ticSummary==0] <- NA 
ticSummary = na.omit(ticSummary)

TIC_threshold = Otsu_thresholding(log10(ticSummary))

####  vi). GMM, signal-to-noise and separation score calculations. #### 

# generate empty export table to be filled in
columns = c("MarkerID","SeparationScore", "PositivePixels", "MeanNoiseIntensity", "MeanSignalIntensity") # , "SNR") 
gmmTable = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(gmmTable) = columns
markerID = colnames(Final_intensity_matrix_targeted)

# need to run this outside the pdf so that variables are saved
for (k in 1:ncol(Final_intensity_matrix_targeted)) {
  # remove pixels below the TIC threshold
  x = Final_intensity_matrix_targeted[,k][ticSummary>TIC_threshold]
  
  # keep all pixels
  x = Final_intensity_matrix_targeted[,k]
  # remove pixels with value zero
  x = x[x>0]
  
  # optimal number of Gaussians
  ICL_score = mclust::mclustICL(data = log10(1+ x),G = 1:5,modelNames = c("V"))
  ICL_score = c(ICL_score[[1]],ICL_score[[2]],ICL_score[[3]],ICL_score[[4]],ICL_score[[5]])
  
  # run GMM
  N_gaussian_curves = as.numeric(which.max((ICL_score)))
  #N_gaussian_curves = 2
  GMM_model = Mclust(data = log10(1+x),G = N_gaussian_curves, modelNames = "V")
  
  # get density as data frame and max density for plotting 
  density = densityMclust(x, plot = F) # same as Mclust but has additional parameter density
  #density = densityMclust(GMM_model$data, plot = F) # same as Mclust but has additional parameter density
  density = data.frame(density[17]) # convert to dataframe to get max value
  maxDensity = max(density)
  
  # Define separability score
  
  # compute mean of each Gaussian
  mean1 = GMM_model$parameters$mean[1]
  mean2 = GMM_model$parameters$mean[2]
  
  # compute standard deviation of each Gaussian
  sd1 = sqrt(GMM_model$parameters$variance$sigmasq[1])
  sd2 = sqrt(GMM_model$parameters$variance$sigmasq[2])
  
  # calculate separation score
  meanDiff = abs(mean1 - mean2)
  SNR = mean2/mean1
  
  # proportion positive pixels
  positivePixels = length(which(GMM_model$classification == 2))
  proportionPositive = positivePixels/length(GMM_model$classification)
  pi2 = proportionPositive
  pi1 = 1 - pi2
  
  denominator = sqrt((sd1^2/pi1) + (sd2^2/pi2))
  separationScore = meanDiff/denominator
  
  # Create table for export
  
  # add markerID and separation score to dataframe
  gmmTable[k,1] = as.character(markerID[k])
  gmmTable[k,2] = separationScore
  gmmTable[k,3] = proportionPositive
  gmmTable[k,4] = mean1
  gmmTable[k,5] = mean2
  #exportTable[k,6] = SNR
  
}

# merge separation score with Mean_peak_intensity table
#exportTable = cbind(test_MeanVarIntensity, exportTable)

# relocate marker name to first column
gmmTable = gmmTable %>% relocate(MarkerID)
gmmTable[is.na(gmmTable)] <- 0

# which markers separate
cat(paste(length(gmmTable$MarkerID[which(gmmTable$SeparationScore != 0)]), "out of", length(gmmTable$MarkerID)
          , "markers show signal and noise separation ."))

# which markers didn't work
markers_no_separation = gmmTable$MarkerID[which(gmmTable$SeparationScore == 0)]

cat(paste("Markers with poor signal-noise separation are: ")) # , markers_no_separation[1], ",", markers_no_separation[2], ",", markers_no_separation[3]), ".")
for (i in seq_along(markers_no_separation)){
  cat(paste(markers_no_separation[i]), "\n")
}


####  vii). Generating GMM pdf #### 

cat(paste("Run GMM pdf..."))

# create path for all processed files
Output_path = args[4]
Output_path = "mnt/msi_volume/processed_files"
Output_path = paste(Output_path,"/",sep="")
print(Output_path)

# if directory does not exist, create it
if( !dir.exists(Output_path)) {
  dir.create(Output_path)
}

# create directory for the experiment
experiment_path = paste(Output_path, Experiment_name, sep="")
experiment_path = paste(experiment_path,"/",sep="")
print(experiment_path)

# if directory does not exist, create it
if( !dir.exists(experiment_path)) {
  dir.create(experiment_path)
}

# create folder path for each sample
sample_path = paste(experiment_path,Sample_name,sep="")

# if directory does not exist, create it
if( !dir.exists(sample_path)) {
  dir.create(sample_path)
}

# create GMM pdf path
GMM_pdf_path = paste("GMM_plot_",Sample_name,".pdf",sep ="" )
GMM_pdf_path = paste(sample_path,"/", GMM_pdf_path ,sep = "")

{
  pdf(GMM_pdf_path, height = 30, width = 15, title = "GMM Plots")
  par(las=1, mfrow=c(7,2), mar=c(8,8,12,5)) # first is bottom, second is left, third is top, fourth is right
  for (k in 1:ncol(Final_intensity_matrix_targeted)) {
    #x = Final_intensity_matrix_targeted[,k][log_testTIC>test_thresholdTIC]
    #x = TargetedIntensityDF[,k]
    x = Final_intensity_matrix_targeted[,k]
    x = x[x>0]
    
    # estimation of optimal number of Gaussians
    ICL_score = mclust::mclustICL(data = log10(1+ x),G = 1:5,modelNames = c("V"))
    ICL_score = c(ICL_score[[1]],ICL_score[[2]],ICL_score[[3]],ICL_score[[4]],ICL_score[[5]])
    
    # run GMM
    N_gaussian_curves = as.numeric(which.max((ICL_score)))
    #N_gaussian_curves = 2
    GMM_model = Mclust(data = log10(1+x), G = N_gaussian_curves, modelNames = "V")
    
    # get density as data frame and max density for plotting 
    density = densityMclust(x, plot = F)
    #density = densityMclust(GMM_model$data, plot = F) # same as Mclust but has additional parameter density
    density = data.frame(density[17]) # convert to dataframe to get max value
    maxDensity = max(density)
    
    # Define separability score
    
    # compute mean of each Gaussian
    mean1 = GMM_model$parameters$mean[1]
    mean2 = GMM_model$parameters$mean[2]
    meanDiff = abs(mean1 - mean2)
    
    # compute standard deviation of each Gaussian
    sd1 = sqrt(GMM_model$parameters$variance$sigmasq[1])
    sd2 = sqrt(GMM_model$parameters$variance$sigmasq[2])
    
    # proportion positive pixels
    positivePixels = length(which(GMM_model$classification == 2))
    totalPixels = length(GMM_model$classification)
    proportionPositive = positivePixels/totalPixels
    proportionNegative = 1 - proportionPositive
    pi2 = proportionPositive
    pi1 = proportionNegative
    
    # separability score
    denominator = sqrt((sd1^2/pi1) + (sd2^2/pi2))
    separationScore = meanDiff/denominator
    
    # plot GMM
    hist(log10(1+x),
         breaks = 50,
         freq = F,
         main = paste(colnames(Final_intensity_matrix_targeted[1,])[k]),
         xlab = paste(colnames(Final_intensity_matrix_targeted[1,])[k],"Intensity"),
         cex.main = 1.5, 
         cex.lab = 1.2) # blank xlab for next line
    for (j in 1:N_gaussian_curves) {
      curve(dnorm(x,
                  mean = GMM_model$parameters$mean[j],
                  sd = sqrt(GMM_model$parameter$variance$sigmasq[j]))*GMM_model$parameters$pro[j],
            add = T,
            lwd=2,
            col=rainbow(N_gaussian_curves)[j])
      text(x = max(GMM_model$data), y = (maxDensity - (0.35*maxDensity)), # coordinates of the text
           labels=paste("Separation Score: ",round(separationScore, 3)),
           adj = 1,
           cex=1, 
           font=1)
    }
  }
  dev.off()
}


#### viii). Generating QC pdf #### 

cat(paste("Creating QC report pdf..."))
# create QC pdf path
QC_pdf_path = paste("QC_report_", Sample_name, ".pdf",sep ="" )
QC_pdf_path = paste(sample_path,"/", QC_pdf_path , sep = "")

{
  pdf(QC_pdf_path, height = 30, width = 15, title = "QC Report")
  par(las=1, mfrow=c(4,2), mar=c(8,8,12,5)) # first is bottom, second is left, third is top, fourth is right
  
  
  # i). Figure 1. Mean-variance plot of peaks in correspondence matrix
  plot(Mean_signal, Sd_signal, log="xy", 
       # label annotated peaks in red, unannotated in grey
       pch=21, bg=string.to.colors(is.na(Correspondence_matrix$Annotation_peaks), colors = c("grey","red3")),
       cex=1.5,xlab="Mean signal", ylab="Standard Deviation signal", cex.lab=1.2, 
       main = "Mean-Variance of Annotated Peaks (red) and Untargeted Peaks (grey)",
       cex.main = 1)
  
  
  # ii). Figure 2. Geary's C score box plot
  
  # boxplot comparing values of annotated metapeaks vs. untargetted
  boxplot(geary_vector~is.na(Correspondence_matrix$Annotation_peaks), 
          main = "Targetted vs. Untargetted Geary's C",
          cex.main = 2,
          xlab = "",
          ylab = "Geary's C score", 
          names = c("Annotated Peaks", "Untargetted Peaks"),
          ylim = c(0,1),
          col = c("red3","lightgrey"))
  
  # iii). Figure 3. Corrected standard deviation vs Geary score for each peak
  plot(Correspondence_matrix$Corrected_sd, Correspondence_matrix$Geary,
       pch=21,
       # label annotated peaks in red, unannotated in grey
       bg=string.to.colors(is.na(Correspondence_matrix$Annotation_peaks), colors = c("grey","red3")),
       xlab = "Corrected Standard Deviation",
       ylab = "Geary's C score",
       main = "Geary vs. SD",
       cex.main = 2,
       ylim = c(0,1))
  
  
  # iv). Figure 4. Histogram
  hist(log10(1+ticSummary), 
       100, 
       main = "TIC histogram", 
       cex.main = 2,
       xlab = "Total Ion Count (Log10)",
       xaxs='i',
       yaxs='i',
       #ylim = c(0,max(log10(1+ticSummary)*1.20)),
       cex.lab=1.3)
  abline(v= TIC_threshold, lwd=2, col="red", lty=2)
  
  dev.off()
}
cat(" ...done ! \n")


#### vii). Saving ion images #### 

# save spatial distribution as .png file
cat(paste("Saving ion images for each channel...\n"))

cat(paste("Create path..."))
Spatial_dir = paste(sample_path,"/", "ion_images",sep = "")
# if directory does not exist, create it
if( !dir.exists(Spatial_dir)) {
  dir.create(Spatial_dir)
}


cat(" ...done ! \n")
cat(paste("Saving annotated single channel images...\n"))

# show all ion images 
for (n in 1:length(Correspondence_matrix$Annotation_peaks)){
  if(!is.na(Correspondence_matrix$Annotation_peaks[n])){
    Spatial_path = paste(Spatial_dir,"/", Correspondence_matrix$Annotation_peaks[n], ".png", sep = "")
    print(Spatial_path)
    png(file=Spatial_path, width=1200, height=700)
    Plot_channel(n)
    dev.off()
  }
}

cat(" ...done ! \n")
cat(paste("...Saving ion images complete. \n"))
#dev.new()

#### viii). Exporting processed data #### 

cat("Exporting processed files...")
print("Intensity data")
write.table(raw_intensity, file = paste(sample_path,"/", Sample_name, "_raw_intensity_data.txt", sep = ""), sep="\t", quote = F, row.names = F)
print("Location data")
write.table(Location_pixels, file = paste(sample_path,"/", Sample_name, "_location_data.txt", sep = ""), sep="\t", quote = F, row.names = F)
print("Correspondance Matrix")
Correspondence_matrix = dplyr::arrange(Correspondence_matrix, Correspondence_matrix$Geary)
write.table(Correspondence_matrix, file = paste(sample_path,"/", Sample_name, "_export_table.txt", sep = ""), sep="\t", quote = F, row.names = F)
print("Final Intensity Matrix")
write.table(Final_intensity_matrix, file = paste(sample_path,"/", Sample_name, "_final_intensity_table.txt", sep = ""), sep="\t", quote = F, row.names = F)
print("Processed data")
writeMSIData(peakPre, file = Sample_name, name=paste(Sample_name,"_pre-processed",sep = ""), folder = sample_path,
             mz.type = "32-bit float", intensity.type="32-bit float")
cat(" ... FINISHED ! \n")



