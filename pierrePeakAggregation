install.packages("fftw", repos="http://olafmersmann.github.io/drat")
library(fftw)
library(sparseMatrixStats)
library(N2R)
library(Matrix)
library(igraph)
fftwtools::fftw(data_raw, inverse=0, HermConj=1, n=NULL)


# estimatenoiseMAD function from Cardinal GitHub
{.estimateNoiseMAD <- function(x, blocks=1, fun=mean, tform=diff) {
  mad <- function(y, na.rm) {
    center <- fun(tform(y), na.rm=na.rm)
    fun(abs(y - center), na.rm=na.rm)
  }
  if ( blocks > 1 ) {
    t <- seq_along(x)
    xint <- split_blocks(x, blocks=blocks)
    tint <- split_blocks(t, blocks=blocks)
    noiseval <- sapply(xint, mad, na.rm=TRUE)
    noiseidx <- sapply(tint, mean, na.rm=TRUE)
    noise <- interp1(noiseidx, noiseval, xi=t, method="linear",
                     extrap=fun(noiseval, na.rm=TRUE))
    noise <- supsmu(x=t, y=noise)$y
  } else {
    noise <- mad(x, na.rm=TRUE)
    noise <- rep(noise, length(x))
  }
  noise
}

# cardinal MAD noise estimation function. set SNR to 10
peakPick.mad <- function(x, SNR=10, window=25, blocks=1, fun=mean, tform=diff, ...) {
  noise <- .estimateNoiseMAD(x, blocks=blocks, fun=fun, tform=tform)
  maxs <- (locmax(x, halfWindow =window))
  peaks <- intersect(maxs, which(x / noise >= SNR))
  return(peaks)
}

# define function to convert to mz scale
Convert_to_mz_scale = function(x,range_peaks, N_features) {
  scale_vector = base::seq(range_peaks[1], range_peaks[2], length.out =N_features )
  return(scale_vector[x])
}

# define function for hierarchical clustering with complete linkage
HC_complete_linkage_function = function(x,Threshold_height=2) {
  sub_clustering = 1
  # if there is more than one element in the region, run the distance function
  if (length(x)>1) {
    dist_matrix = dist(x)
    hc_clustering = hclust(dist_matrix,method = "complete")
    sub_clustering = cutree(hc_clustering,h =Threshold_height )
  }
  return(sub_clustering)
}

# define a modified which max function
which_max_modified = function(x) {
  if (sum(x)==0) {
    y = NA
  }
  else {
    x[x==0] = NA
    y = which.min(x)
  }
  return(y)
}

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}

}

### MY FUNCTION ###

pierrePeakAggregation <- function(x, annotation_table, Delta_mz, mz_threshold_association) { # x is a data matrix of normalised and baseline corrected peak intensities before feature space reduction. Delta_mz is a user-defined threshold for the distance between peaks that constitutes a split for h-clustering
  # apply MAD peak picking function to all pixels
  List_peaks = apply(x, MARGIN = 2, FUN = peakPick.mad)
  Concatenated_peaks = unlist(List_peaks)
  
  ## Convert the scale of these peaks to the correct mz range
  range_peaks = range(mz(peakPre))
  N_features = length(mz(peakPre))
  Rescaled_concatenated_peaks = Convert_to_mz_scale(Concatenated_peaks, range_peaks,N_features)
  
  # how many different peaks there are
  Concatenated_peaks_unique = unique(Rescaled_concatenated_peaks)
  # order the peaks
  Concatenated_peaks_ordered = Concatenated_peaks_unique[order(Concatenated_peaks_unique)]
  # calculate difference/distance between each ordered peak
  Diff_peak_vector = diff(Concatenated_peaks_ordered)
  
  # indices of the m/z values where splits should occur
  Location_split = which(Diff_peak_vector>Delta_mz)
  Concatenated_peaks_ordered[Location_split]
  # the m/z locations where splits should be done 
  Location_split = Concatenated_peaks_ordered[Location_split-1] # which indices to grab in the Concatenated_peaks_ordered vector
  
  # which peaks should be in which segments
  Segmentation_peaks = cumsum(Diff_peak_vector>Delta_mz)+1
  Segmentation_peaks = c(1,Segmentation_peaks)
  
  # splits the data (vector) into groups defined by f  
  Splitted_peaks = base::split(Concatenated_peaks_ordered,f = Segmentation_peaks) 
  
  # apply H-clustering on each split group
  Aggregated_HC = lapply(Splitted_peaks,FUN = HC_complete_linkage_function)
  
  ## function for unlisting the aggregated HC list
  Aggregated_HC_unlist = c()
  for (k in 1:length(Aggregated_HC)) {
    # paste together the split index (which split group it is in) and the index of the peak within that split group
    x = paste(rep(k,length(Aggregated_HC[[k]])),Aggregated_HC[[k]])
    # add the new iteration to the list
    Aggregated_HC_unlist =c(Aggregated_HC_unlist,x)
  }
  
  # takes the mean m/z value of each concatenated index (i.e. "1 34")
  Centroided_metapeaks_unweighted = aggregate(Concatenated_peaks_ordered,by=list(Aggregated_HC_unlist),FUN = mean) ## index jumping here!! 
  
  # which m/z values correspond to these indices
  Centroided_metapeaks_unweighted = Centroided_metapeaks_unweighted$x
  # these are the metapeaks
  
  # load in pixel location
  Location_points = as.data.frame(pData(peakPre))
  Location_points = Location_points[,c(2,3)]
  
  ##  metapeak counts
  # nearest neighbours search between the centroid of metapeaks and the detected peaks
  NN_search = N2R::crossKnn(mB= matrix(Centroided_metapeaks_unweighted,ncol = 1), m = matrix(Rescaled_concatenated_peaks, ncol = 1)
                            , k = 10,verbose = T, indexType ="L2")
  
  Metapeak_assignments = apply(as.matrix(NN_search),MARGIN = 2,FUN = which.max)
  # counts for each metapeak and their indices
  Metapeak_assignments_count = table(Metapeak_assignments)
  # only counts, removes indices
  Metapeak_assignments_count = as.numeric(table(Metapeak_assignments))
  
  # Filtering: remove peaks that are in less than 1% of pixels 
  N_pixels_total = dim(peakPre)[2]
  Threshold_Detection = N_pixels_total*0.01
  Selected_peak_location = Centroided_metapeaks_unweighted[Metapeak_assignments_count > Threshold_Detection]
  
  ## Create flattened dataframe for all mass features (including adduct ions)
  
  peakAnnotation_extended = annotation_table
  # create vectors for each adduct ion
  FeatureMassHydrogen = peakAnnotation_extended$FeatureMass
  FeatureMassSodium = (FeatureMassHydrogen - 1) + 23
  FeatureMassPotassium = (FeatureMassHydrogen - 1) + 39
  FeatureMassAmmonium = (FeatureMassHydrogen - 1) + 18
  peakAnnotation_extended = cbind(peakAnnotation_extended, FeatureMassSodium, FeatureMassPotassium, FeatureMassAmmonium)
  # bind the feature masses to the feature names and adduct ions, then bind all features together to get flattened list
  peakAnnotation_extended = rbind(cbind(peakAnnotation_extended$Name, peakAnnotation_extended$FeatureMass),
                                  cbind(paste(peakAnnotation_extended$Name,"Sodium"), peakAnnotation_extended$FeatureMassSodium),
                                  cbind(paste(peakAnnotation_extended$Name,"Ammonium"), peakAnnotation_extended$FeatureMassAmmonium))
  # convert to dataframe
  peakAnnotation_extended = as.data.frame(peakAnnotation_extended)
  colnames(peakAnnotation_extended) = c("FeatureAnnotation","Location")
  # convert peak location to numeric
  peakAnnotation_extended$Location = as.numeric(peakAnnotation_extended$Location)
  
  # cross compare total reference list with metapeaks (find closest peak across lists)
  Mapping_meta = crossKnn(mA = matrix(Selected_peak_location, ncol = 1),
                          mB= matrix(peakAnnotation_extended$Location, ncol = 1), k = 10, indexType = "L2")
  
  
  # here was mz_threshold_definition
  # only keep those above threshold
  Mapping_meta[Mapping_meta>mz_threshold_association]=0
  
  Mapping_meta_cleaned = apply(as.matrix(Mapping_meta),MARGIN = 2,FUN = which_max_modified)
  # table comparing the metapeak location with the detected peak location and the associated peak name
  Correpondance_matrix = data.frame(mz_location = Selected_peak_location,
                                    Expected_mz_location = peakAnnotation_extended$Location[Mapping_meta_cleaned],
                                    Annotation_peaks = peakAnnotation_extended$FeatureAnnotation[Mapping_meta_cleaned])
  
  ## Final Intensity Matrix
  mz_scale_vector = base::seq(range_peaks[1], range_peaks[2], length.out = N_features)
  
  # sum up the intensities of all mz locations surrounding each metapeak
  Final_intensity_matrix = c()
  for (k in 1:nrow(Correpondance_matrix)) {
    # mz location of a given peak
    Metapeak_location_mz_temp = Correpondance_matrix$mz_location[k]
    # which mz locations to sum up (indices)
    mz_location_to_take = which(mz_scale_vector > (Metapeak_location_mz_temp-mz_threshold_association/2) & (mz_scale_vector < Metapeak_location_mz_temp + mz_threshold_association/2))
    # sum up the intensities of these peaks for every pixel
    temporary_signal = colSums(data_raw[mz_location_to_take,])
    # update matrix
    Final_intensity_matrix = cbind(Final_intensity_matrix,temporary_signal)
    #colnames(Final_intensity_matrix) = Metapeak_location_mz_temp
  
  # which of the metapeaks in the final intensity matrix are annotated peaks
  Final_intensity_matrix_targeted = Final_intensity_matrix[,!is.na(Correpondance_matrix$Annotation_peaks)]
  colnames(Final_intensity_matrix_targeted) = Correpondance_matrix$Annotation_peaks[!is.na(Correpondance_matrix$Annotation_peaks)]
  
  # as dataframe
  Final_intensity_matrix_targeted = as.data.frame(Final_intensity_matrix_targeted)
  
  # FINAL STEP: returns the metapeak locations
  test_returnList <- list("Final_Intensity_Matrix" = Final_intensity_matrix, "Targeted_Final_Intensity_Matrix" = Final_intensity_matrix_targeted, "Correspondance_Matrix" =  Correpondance_matrix, 
                     "Metapeaks" = Centroided_metapeaks_unweighted, "Split_Locations" = Location_split,
                     "Coordinates" = Location_points, "Counts" = Metapeak_assignments_count)
  
  # Final_Intensity_Matrix will have the summed intensity of all peaks
  
  return(test_returnList)
}

# test function
test_output <- pierrePeakAggregation(data_raw, annotation_table = peakAnnotation, Delta_mz = 3, mz_threshold_association = 1)
test_output



# import file
Path_to_imzml_file = "mnt/msi_volume/experiments/ImzML-47plex/part-2/colonadenocarcinoma.imzML"
rawFile = readMSIData(Path_to_imzml_file)

# use pre-processed peaks as data
data_raw = peakPre@imageData$data$intensity
plot(data_raw[,1500],type="l")


# apply MAD peak picking function to all pixels
List_peaks = apply(data_raw, MARGIN = 2, FUN = peakPick.mad)
# unlist 
Concatenated_peaks = unlist(List_peaks)

## Convert the scale of these peaks to the correct mz range
range_peaks = range(mz(peakPre))
N_features = length(mz(peakPre))
Rescaled_concatenated_peaks = Convert_to_mz_scale(Concatenated_peaks, range_peaks,N_features)

# how many different peaks there are
Concatenated_peaks_unique = unique(Rescaled_concatenated_peaks)
# order the peaks
Concatenated_peaks_ordered = Concatenated_peaks_unique[order(Concatenated_peaks_unique)]
# calculate the difference/distance between each ordered peak
Diff_peak_vector = diff(Concatenated_peaks_ordered)
# user-defined threshold for the distance between peaks that constitutes a split for h-clustering
Delta_mz = 3
# indices of the m/z values where splits should occur
Location_split = which(Diff_peak_vector>Delta_mz)
Concatenated_peaks_ordered[Location_split]
# the m/z locations where splits should be done 
Location_split = Concatenated_peaks_ordered[Location_split-1] # which indices to grab in the Concatenated_peaks_ordered vector

# plot where splits should be
par(las = 1,mar = c(1, 5, 5, 1))
hist(Rescaled_concatenated_peaks, 2000, xlab="m/z", xlim=range_peaks, xaxs="i")
abline(v=(refList),lwd=0.5,col="red", lty = 2) # hydrogen adduct ions
abline(v=refList+38,lwd=0.5,col="green", lty = 2) # potassium adduct ions
abline(v=refList+18,lwd=0.5,col="blue", lty = 2) # ammonium adduct ions

# which peaks should be in which segments
Segmentation_peaks = cumsum(Diff_peak_vector>Delta_mz)+1
Segmentation_peaks = c(1,Segmentation_peaks)
# plot histogram of peak frequencies accross all spectra

# splits the data (vector) into groups defined by f  
Splited_peaks = base::split(Concatenated_peaks_ordered,f = Segmentation_peaks) 

# Hierarchical clustering of each split group
Aggregated_HC = lapply(Splited_peaks,FUN = HC_complete_linkage_function)

## function for unlisting the aggregated HC list
# initialize empty list
Aggregated_HC_unlist = c()
for (k in 1:length(Aggregated_HC)) {
  # paste together the split index (which split group it is in) and the index of the peak within that split group
  x = paste(rep(k,length(Aggregated_HC[[k]])),Aggregated_HC[[k]])
  # add the new iteration to the list
  Aggregated_HC_unlist =c(Aggregated_HC_unlist,x)
}

# takes the mean m/z value of each concatenated index (i.e. "1 34")
Centroid_metapeaks_unweighted = aggregate(Concatenated_peaks_ordered,by=list(Aggregated_HC_unlist),FUN = mean) ## index jumping here!! 
## ^ Comment: this seems to be ordered first by Aggregated_HC_unlist, then by the m/z values
### ^^ Comment: this is the unweighted mean, i.e, it is doesn't take the amount of peaks at these m/z values into account, as it is built off the unique function

# the list of metaclusters
List_metacluster = Centroid_metapeaks_unweighted$Group.1 # cluster number
Centroid_metapeaks_unweighted = Centroid_metapeaks_unweighted$x # metapeak location

# create a 1 column matrix containing metapeak mz location and cluster number
Assignment_table = matrix(data = Aggregated_HC_unlist,ncol=1)
# assign each detected peak as a name 
rownames(Assignment_table) = Concatenated_peaks_ordered
# assigns each peak to a metapeak cluster. says which metacluster each peak is in
Meta_peaks_assignments = Assignment_table[as.character(Rescaled_concatenated_peaks),]

# counts for each metapeak and their indices
Meta_peaks_assignments_count = as.numeric(table(factor(Meta_peaks_assignments,List_metacluster)))
# plot ordered Metapeak and their counts
plot(Meta_peaks_assignments_count[order(Meta_peaks_assignments_count)])

# plot metapeaks and their counts 
plot(Centroid_metapeaks_unweighted[order(Centroid_metapeaks_unweighted)],Meta_peaks_assignments_count[order(Centroid_metapeaks_unweighted)],type="l",
     xlab = "Metapeak Location",
     ylab = "Counts")
# spectrum of all metapeaks compared to refList and 
abline(v=refList,lwd=2,lty=2,col="red") # hydrogen
abline(v=refList+22,lwd=2,lty=2,col="green") # sodium
abline(v=refList+17,lwd=2,lty=2,col="purple") # ammonium

Refined_peak_location_weighted = aggregate(Rescaled_concatenated_peaks,by=list(Meta_peaks_assignments),FUN = mean)

#### new stuff from here 08.12.2022


# x is the counts of each m/z peak
x = Meta_peaks_assignments_count[order(Centroid_metapeaks_unweighted)]
# character vector of all metapeak masses
names(x) = Centroid_metapeaks_unweighted[order(Centroid_metapeaks_unweighted)] # name each count by it's associated m/z
# generate list of all metapeak values
List_values_to_show = base::seq(range_peaks[1],range_peaks[2],length.out = N_features) # all m/z values detected in order
x_bis = x
x_bis = x[c(names(x),List_values_to_show)]
# name ALL detected mz values 
names(x_bis) = c(names(x),List_values_to_show)
x_bis = x_bis[order(as.numeric(names(x_bis)))]
# set counts that are NA to 0
x_bis[is.na(x_bis)]=0
x_bis
u = as.numeric(names(x_bis))

par(las=1,bty="l")
# peak plot
plot(u,x_bis,type="l",xlab="m/z",ylab="Number of pixels where detected", main = "Tibshirani clustering for peak aggregation")

# zoomed in plot
plot(u,x_bis,type="l",xlim=c(1050,1250),xlab="m/z",ylab="Number of pixels where detected", main = "Tibshirani clustering for peak aggregation")
plot(u,x_bis,type="l",xlim=c(1285,1317),xlab="m/z",ylab="Number of pixels where detected", main = "Tibshirani clustering for peak aggregation")


# filter peaks out that are in less than 1% of pixels 
N_pixels_total = dim(peakPre)[2]
Threshold_Detection = N_pixels_total*0.01
Selected_peak_location = Centroid_metapeaks_unweighted[Meta_peaks_assignments_count>Threshold_Detection]

# include refList lines and tolerance
abline(v=refList,lwd=1,lty=2,col="red") # hydrogen
abline(v=refList+1,lwd=1,lty=2,col="blue") # tolerance of +1Da
abline(v=refList-1,lwd=1,lty=2,col="blue") # tolerance of -1Da

# include threshold line
abline(h=Threshold_Detection,lwd=1,col="grey",lty=2)
abline(v=refList+22,lwd=1,lty=2,col="green") # sodium.5
abline(v=refList+17,lwd=1,lty=2,col="purple") # ammonium
text(x= refList+22+1, y = max(x_bis+400), labels = "Na+", cex = 0.4)
text(x= refList+17+1.5, y = max(x_bis+400), labels = "NH4+", cex = 0.4)
x = Meta_peaks_assignments_count[order(Centroid_metapeaks_unweighted)]

## Create flattened dataframe for all mass features (including adduct ions)
peakAnnotation_extended = peakAnnotation

# create vectors for each adduct ion
FeatureMassHydrogen = peakAnnotation$FeatureMass
FeatureMassSodium = (FeatureMassHydrogen - 1) + 23
FeatureMassPotassium = (FeatureMassHydrogen - 1) + 39
FeatureMassAmmonium = (FeatureMassHydrogen - 1) + 18
# create peakAnnotation df containing adduct ions
peakAnnotation = cbind(peakAnnotation, FeatureMassSodium, FeatureMassPotassium, FeatureMassAmmonium)

# bind the feature masses to the feature names and adduct ions, then bind all features together to get flattened list
peakAnnotation_extended = rbind(cbind(peakAnnotation$Name,peakAnnotation$FeatureMass),
                                cbind(paste(peakAnnotation$Name,"Sodium"),peakAnnotation$FeatureMassSodium),
                                cbind(paste(peakAnnotation$Name,"Ammonium"),peakAnnotation$FeatureMassAmmonium))
# convert to dataframe
peakAnnotation_extended = as.data.frame(peakAnnotation_extended)
# add names
colnames(peakAnnotation_extended) = c("FeatureAnnotation","Location")
# convert peak location to numeric
peakAnnotation_extended$Location = as.numeric(peakAnnotation_extended$Location)

# cross compare total reference list with metapeaks (find closest peak across lists)
Mapping_meta = crossKnn(mA = matrix(Selected_peak_location, ncol = 1),
                        mB= matrix(peakAnnotation_extended$Location,ncol = 1), k = 10, indexType = "L2")
# define mz threshold
mz_threshold_association = 1
# only keep those above threshold
Mapping_meta[Mapping_meta>mz_threshold_association]=0

# create modified which max function
which_max_modified = function(x) {
  if (sum(x)==0) {
    y = NA
  }
  else {
    x[x==0] = NA
    y = which.min(x)
  }
  return(y)
}

Mapping_meta_cleaned = apply(as.matrix(Mapping_meta),MARGIN = 2,FUN = which_max_modified)

# table comparing the metapeak location with the detected peak location and the associated peak name
Correpondance_matrix = data.frame(mz_location = Selected_peak_location,
                                  Expected_mz_location = peakAnnotation_extended$Location[Mapping_meta_cleaned],
                                  Annotation_peaks = peakAnnotation_extended$FeatureAnnotation[Mapping_meta_cleaned])


mz_scale_vector = base::seq(range_peaks[1], range_peaks[2], length.out =N_features )

# sum up the intensities of all mz locations surrounding each metapeak
Final_intensity_matrix = c()
for (k in 1:nrow(Correpondance_matrix)) {
  # mz location of a given peak
  Metapeak_location_mz_temp = Correpondance_matrix$mz_location[k]
  # which mz locations to sum up
  mz_location_to_take = which(mz_scale_vector> (Metapeak_location_mz_temp-mz_threshold_association) & (mz_scale_vector <Metapeak_location_mz_temp+ mz_threshold_association))
  # sum up the intensities of these peaks for every pixel
  signal_temp = colSums(data_raw[mz_location_to_take,])
  # update matrix
  Final_intensity_matrix = cbind(Final_intensity_matrix,signal_temp)
}

#Loading the location part

Location_points = as.data.frame(pData(peakPre))
Location_points = Location_points[,c(2,3)]

# which of the metapeaks in the final intensity matrix are annotated peaks
Final_intensity_matrix_targeted = Final_intensity_matrix[,!is.na(Correpondance_matrix$Annotation_peaks)]
colnames(Final_intensity_matrix_targeted) = Correpondance_matrix$Annotation_peaks[!is.na(Correpondance_matrix$Annotation_peaks)]
Final_intensity_matrix_targeted = as.data.frame(Final_intensity_matrix_targeted)

# correlation heatmap between all markers
pheatmap(cor(log2(Final_intensity_matrix_targeted)),clustering_method = "ward",show_colnames = F)

View(data.frame(apply(log2(Final_intensity_matrix),MARGIN = 2,FUN = var)))

# step by step output
Final_intensity_matrix_targeted
Correpondance_matrix

# function outout
test_output$Targeted_Final_Intensity_Matrix
test_output$Correspondance_Matrix

# spatial plot 
plot(Location_points,pch=21,lwd=0,bg=color_convertion(log2(Final_intensity_matrix[,502])),cex=0.5, main = Correpondance_matrix$Annotation_peaks[400])

# mean-var plot comparing targeted and untargeted
plot(colMeans(log(Final_intensity_matrix)),apply(log2(Final_intensity_matrix),MARGIN = 2,FUN = var),
     pch=21,bg=string.to.colors(is.na(Correpondance_matrix$Annotation_peaks)),
     xlab = "Log Mean",
     ylab = "Log Variance",
     main = "Log Mean-Var plot (blue = targeted, red = untargeted)")
# index
text(colMeans(log(Final_intensity_matrix)), apply(log2(Final_intensity_matrix), MARGIN = 2,FUN = var), labels = 1:ncol(Final_intensity_matrix), cex = 0.7)
text(colMeans(log(Final_intensity_matrix)), apply(log2(Final_intensity_matrix), MARGIN = 2,FUN = var), labels = Correpondance_matrix$Annotation_peaks, cex = 0.7)


hist(log2(Final_intensity_matrix[,27]),100)

Coexpression_score = (Final_intensity_matrix[,137]*Final_intensity_matrix[,24])/sqrt(Final_intensity_matrix[,137]^2+Final_intensity_matrix[,24]^2)
plot(Location_points,pch=21,lwd=0,bg=color_convertion(Coexpression_score),cex=0.5)

# correlation between total signal and peak signal. Low correlation could/should be interesting, high correlation not

hist(corr_total_signal,100)
View(data.frame(corr_total_signal ))

Total_signal = colSums(data_raw)
corr_total_signal = apply(Final_intensity_matrix,MARGIN = 2,FUN = function(x){cor(x,Total_signal)})
Linear_residual_fit=function(x,Total_signal=Total_signal) {
  U = data.frame(Total_signal = log2(Total_signal),Signal = log2(x))
  m = loess(Signal~Total_signal,U,degree = 2)
  return(m$residuals)
}

marker_intensity_vector = Final_intensity_matrix[,16]
plot(marker_intensity_vector,Total_signal,log="xy")

Normalised_data = apply(Final_intensity_matrix,MARGIN = 2,FUN = function(x) {Linear_residual_fit(marker_intensity_vector,Total_signal)})
plot(Normalised_data , Final_intensity_matrix)

# not necessary but maybe useful code
Mean_spectrum = rowMeans(data_raw)
plot(Mean_spectrum, type = "l")
# convert to mz scale
plot(Convert_to_mz_scale(1:length(Mean_spectrum), range_peaks,N_features),Mean_spectrum, type = "l")
# plot reference antibodies on top
abline(v=refList,lwd=1,lty=2,col="red") # hydrogen
# plot zoomed in
plot(Convert_to_mz_scale(1:length(Mean_spectrum), range_peaks,N_features),Mean_spectrum, type = "l",xlim=c(1010,1040))

View(Final_intensity_matrix_targeted)

#y = peakPick.mad(x[,1])

## Convert the scale of these peaks to the correct mz range
range_peaks = range(mz(peakPre))
N_features = length(mz(peakPre))


# define function
Convert_to_mz_scale = function(x,range_peaks, N_features) {
  scale_vector = base::seq(range_peaks[1], range_peaks[2], length.out =N_features )
  return(scale_vector[x])
}

# rescale peaks
Rescaled_concatenated_peaks = Convert_to_mz_scale(Concatenated_peaks, range_peaks,N_features)
# plot 
hist(Rescaled_concatenated_peaks, 2000, xlab="m/z", xlim=range_peaks, xaxs="i")
abline(v=(refList),lwd=0.5,col="red") # hydrogen adduct ions
abline(v=refList+38,lwd=0.5,col="green") # potassium adduct ions
abline(v=refList+18,lwd=0.5,col="blue") # ammonium adduct ions

# add lines for the reference peaks and their adduct ions
hist(log10(Rescaled_concatenated_peaks),2000,xlab="m/z",xaxs="i")
abline(v=log10(refList),lwd=0.5,col="red") # hydrogen 
abline(v=log10(refList+22),lwd=0.5,col="green") # sodium


hc_analysis = hclust(d = dist(Rescaled_concatenated_peaks),method = "complete")


## Spatial Correlation plots between two images

# create image vectors
marker = Final_intensity_matrix_targeted$GLUT1
marker_sodium = Final_intensity_matrix_targeted$`GLUT1 Sodium`
log_marker = log(marker)
log_marker_sodium = log(marker_sodium)

# correlation coefficient
correlation = cor(marker,marker_sodium)
max(log_marker)*0.02
max(log_marker_sodium)*0.9

# plot 
plot(log_marker,log_marker_sodium, pch = 19, col = "lightblue", 
     main = "Spatial Correlation of Glut1 +H and Glut1 +Na", 
     xlab = "Log PR-A/B + H",
     ylab = "Log PR-A/B + Na")

# Regression line
abline(lm(log_marker_sodium ~ log_marker), col = "red", lwd = 2, lty = 2)
text(x = max(log_marker)*0.25, y = max(log_marker_sodium)*0.9,
     labels=paste("Correlation:", round(correlation, digits = 3)),
     adj = 1,
     cex=1, 
     font=1)

