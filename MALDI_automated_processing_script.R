# Title: Fully automated MALDI processing script
# User: John Abbey
# Start: 29.09.2022
# Latest update: 14.10.2022

#0) Basic functions 

# otsu thresholding function
Otsu_thresholding = function(x,number_bins=100) {
  List_bin = quantile(x,base::seq(from=0,to=1,length.out=number_bins))
  Intravariance_vector = c()
  for (k in 1:number_bins) {
    threshold_temp = List_bin[k]
    s = length(x[x<threshold_temp])*var(x[x<threshold_temp]) + length(x[x>threshold_temp])*var(x[x>threshold_temp])
    Intravariance_vector= c(Intravariance_vector,s)
  }
  Selected_values = List_bin[which.min(Intravariance_vector)]
  return(Selected_values)
}

# convert string to a specific colour function
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

# colour conversion function
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

# obtain frequency of each marker as a percentage of total pixels. developed form peakFilter function Cardinal
getFreq <- function(object, freq.min) {
  if ( length(expr) > 0L ) {
    stats <- c("min", "max", "mean", "var", "nnzero")
  } else {
    stats <- c("nnzero")
  }
  summary <- rowStats(spectra(object), stat=stats, drop=FALSE)
  summary[["count"]] <- summary[["nnzero"]]
  summary[["freq"]] <- summary[["count"]] / ncol(object)
  freq <- summary[["freq"]]
  names_freq <- cbind(markerNames, freq)
  return(names_freq)
}



#I)Package loading

cat("Loading of the libraries.... ")
suppressMessages(library(Cardinal))
suppressMessages(library(dplyr))
suppressMessages(library(hexbin))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(matter))
cat("... done ! \n")


#II) Data and annotation loading

args <- commandArgs(trailingOnly = T)

Path_to_imzml_file = args[1]
Path_to_imzml_file = "mnt/msi_volume/imzml-ambergen-5um/tonsil-1.imzML"

cat("Reading the .imzml and .ibd file  ...")
rawFile <- readMSIData(Path_to_imzml_file)
rawFile <- Cardinal::pull(rawFile)
cat("done ! \n")

#Extracting file name

Sample_name = strsplit(Path_to_imzml_file,split = "/",fixed = T)[[1]]
Sample_name = Sample_name[length(Sample_name)]
Sample_name = strsplit(Sample_name,split = ".",fixed = T)[[1]]
Sample_name = Sample_name[1]

Path_to_peak_annotation = args[2]
Path_to_peak_annotation = "mnt/msi_volume/Panels/AmbergenPilotPanel.csv"

cat("Reading the csv annotation file  ...")
peakAnnotation <- read.delim(Path_to_peak_annotation,sep=",")
peakAnnotation <- peakAnnotation[-3]
refList <- peakAnnotation$FeatureMass
markerNames <- peakAnnotation$Name
markerNames[3] <- if(markerNames[3] == "Actin-alpha-SM"){
  markerNames[3] <- "AASM"
}
cat("done ! \n")

# Set parameters (dev only)

Path_to_parameter_file = args[3]
Path_to_parameter_file = "mnt/msi_volume/Parameter_files/test_param_file.txt"

cat("Loading parameter file...")
Parameters = suppressWarnings(read.table(Path_to_parameter_file,header = F,sep = "="))
Parameters_names = Parameters$V1
Parameters_values = Parameters$V2
names(Parameters_values) = Parameters_names

N_cores = as.numeric(Parameters_values[1])
#N_cores = 2
peak_mz_tolerance = as.numeric(Parameters_values[2])
#peak_mz_tolerance = 0.2
freqMin = as.numeric(Parameters_values[3])
#freqMin = 0.05
cat("done ! \n")



##III) Parallelisation

#if (N_cores >= round(0.8*detectCores(),0) ) {
  #stop("Not enough threads available. Check Parameter.txt")
#}


#if (N_cores >1) {
  #cat("Creating parallel environement...")
  #param <- MulticoreParam(workers = N_cores)
  #setCardinalBPPARAM(param)
  
  #getCardinalBPPARAM()
  #setCardinalBPPARAM(SerialParam)
  #cat("... done ! \n")
#}


###IV)Processing

cat("Running normalisation and baseline correction...")
peakPre <- rawFile %>% 
  process(function(x) length(x) * x / sum(x), label="norm", delay=TRUE) %>% 
  #smoothSignal(method = "gaussian", plot=FALSE) %>%
  reduceBaseline(method="locmin") %>%
  process()
cat("...done ! \n")

cat("Running peak detection...")
peakProc <- peakPre %>%
  peakPick(method="mad") %>% 
  peakAlign(ref = refList, tolerance = peak_mz_tolerance, units = "mz", method="diff") %>% 
  peakFilter(freq.min = freqMin) %>%
  process()
cat("...done ! \n")

cat("Running peak binning and intensity computation....")
binRef <- peakProc@featureData@mz

peaksBinned <- peakPre %>%
  peakBin(ref=binRef, tolerance=peak_mz_tolerance, units="mz", type = "area") %>%
  process()
cat("....done ! \n")

# run peak processing again without filtering so that you get full 15 feature dataset for pixel QC 
cat("Running peak detection without filtering...")
peakProc_nofilter <- peakPre %>%
  peakPick(method="mad") %>% 
  peakAlign(ref = refList, tolerance = peak_mz_tolerance, units = "mz", method="diff") %>% 
  process()
cat("done ! \n")

binRef_noFilter <- peakProc_nofilter@featureData@mz

cat("Running peak binning without filtering and intensity computation....")
peaksBinned_noFilter <- peakPre %>%
  peakBin(ref=binRef_noFilter, tolerance=peak_mz_tolerance, units="mz", type = "area") %>%
  process()
cat("....done ! \n")


###V)Reshaping the data


##### Intensity data with filtering

#Extracting the intensity
Intensity_data = iData(peaksBinned,"intensity")
Intensity_data = t(Intensity_data)
#Annotating the peaks
conserved_peak_mz = peaksBinned@featureData@mz
x = peakAnnotation
rownames(x) = x$FeatureMass
name_conserved_channels = x[as.character(conserved_peak_mz),"Name"]
# create table of antibody names
names_table <- as.data.frame(name_conserved_channels)
names(names_table) <- Sample_name
colnames(Intensity_data) = name_conserved_channels
Intensity_data = as.data.frame(Intensity_data)
#Extracting pixel location
Location_data = as.data.frame(pixelData(peaksBinned))
Location_data = Location_data[,c("x","y")]  


######## Intensity data no filtering
#Extracting the intensity
Intensity_data_noFilter = iData(peaksBinned_noFilter,"intensity")
Intensity_data_noFilter = t(Intensity_data_noFilter)
#Annotating the peaks
conserved_peak_mz_noFilter = peaksBinned_noFilter@featureData@mz
x = peakAnnotation
rownames(x) = x$FeatureMass
name_conserved_channels_noFilter = x[as.character(conserved_peak_mz_noFilter),"Name"]
# create table of antibody names
names_table_noFilter <- as.data.frame(name_conserved_channels_noFilter)
names(names_table_noFilter) <- Sample_name
colnames(Intensity_data_noFilter) = name_conserved_channels_noFilter
Intensity_data_noFilter = as.data.frame(Intensity_data_noFilter)
#Extracting pixel location
Location_data_noFilter = as.data.frame(pixelData(peaksBinned_noFilter))
Location_data_noFilter = Location_data_noFilter[,c("x","y")]  



#VI)Quality control tests

#A)Pixel quality control 

pixelsSummarized <- summarizePixels(peaksBinned, FUN =c(tic="sum", mean = "mean", min = "min", max = "max", sd = "sd", var = "var"), as = "DataFrame")
pixelsSummarized = as.data.frame(pixelsSummarized)
ticSummary = pixelsSummarized$tic

TIC_threshold = Otsu_thresholding(log10(ticSummary))

#B)Feature-wise quality control 

Mean_marker_intensity = colMeans(Intensity_data)

#C) Mean-variance plot
Mean_peak_intensity = summarizeFeatures(peaksBinned,FUN = mean)
Mean_peak_intensity = Mean_peak_intensity@featureData
Mean_peak_intensity = as.data.frame(Mean_peak_intensity)
Mean_peak_values <- Mean_peak_intensity$FUN


colour_panel <- c("Red", "Blue", "Green", "Yellow", "Pink", "Brown", "Purple", "Grey", "Black", "Orange", "Turquoise", "Violet", "deeppink", "darkslategray2", "coral4")
# make sure same amount of colours as passed markers
while(length(colour_panel)!= length(name_conserved_channels)){
  colour_panel <- colour_panel[-length(colour_panel)]
}
# add columns to dataframe
Mean_peak_intensity <- cbind(name_conserved_channels, Mean_peak_intensity)
Mean_peak_intensity <- cbind(Mean_peak_intensity, colour_panel)
`colnames<-`(Mean_peak_intensity, c("Marker", "mz", "Mean", "Colours"))

Var_peak_intensity = summarizeFeatures(peaksBinned, FUN = var)
Var_peak_intensity = Var_peak_intensity@featureData
Var_peak_values <- Var_peak_intensity$FUN

#### old stuff from antibody check table
# get markers that didn't pass peak detection
#failed_markers <- setdiff(peakAnnotation$Name, name_conserved_channels)
# make binary vector for passed/failed
#TF_vector <- peakAnnotation$Name %in% name_conserved_channels
#TF_vector <- as.numeric((as.logical(TF_vector)))
#TF_vector <- as.data.frame(TF_vector)
# add column to dataframe
#peakAnnotation <- cbind(peakAnnotation, TF_vector)

#VII) Generating QC pdf

Output_path = args[4]
#Output_path = "mnt/msi_volume/processed_files"
Output_path = paste(Output_path,"/",sep="")

if( !dir.exists(Output_path)) {
  dir.create(Output_path)
}

cat(paste("Creating QC pdf file ..."))

QC_pdf_path = paste("QC_report_",Sample_name,".pdf",sep ="" )
QC_pdf_path = paste(Output_path,QC_pdf_path ,sep = "")

# start PDF 
{
  pdf(QC_pdf_path,height = 20,width = 15, useDingbats = F, title = "QC Report")
  par(las=1,mfrow=c(3,2),mar=c(6,6,6,4))


# i). Figure 1. Histogram
  hist(log10(1+ticSummary), 100, main = "TIC histogram", xlab = "Total Ion Count (Log10)",xaxs='i',yaxs='i',cex.lab=1.3)
  abline(v= TIC_threshold, lwd=2, col="red", lty=2)

# ii) Figure 2. Mean-variance plot
  plot(Mean_peak_values, Var_peak_values, 
       log="xy", 
       main = "Mean-Variance Plot", 
       xlab = "Mean Peak Intensity", ylab = "Variance", 
       col = colour_panel, type = "p", pch = 19)
#par(xpd=FALSE)
  legend("bottomright",
         inset = c(0.01, 0.01), 
         cex = 0.6, 
         bty = "o", 
         legend = name_conserved_channels, 
         text.col = "black",
         col = Mean_peak_intensity$colour_panel, 
         pch = 19, 
         title = "Marker",
         y.intersp = 0.7)

# iii). Figure 3. Spatial distribution of TIC data
  plot(Location_data, pch = 16, bg = color_convertion(ticSummary), col = "white", main = "Spatial distribution of Total Ion Current")

  #plot(Location_data, pch = 21, bg = color_convertion(ticSummary), col = "white", main = "Spatial distribution of Total Ion Current")

# iv). Figure 4. Percentage Signal in Pixels  
  pixelFreq <- getFreq(peakProc_nofilter, freq.min = freqMin)
  pixelFreq <- data.frame(pixelFreq)
  pixelFreq$freq <- as.numeric(pixelFreq$freq)
  pixelFreq <- arrange(pixelFreq, -freq)
  
  #frequencies <- ifelse(frequencies > 0.25, 0.25, frequencies) 
  par(mar=c(4,4,4,4))
  barplot <- barplot(pixelFreq$freq, 
                     main="Percentage Signal in Pixels",
                     names = pixelFreq$markerNames,
                     las = 2,
                     col = ifelse(pixelFreq$freq < freqMin, "red", "grey")
  )
  abline(h = freqMin, lwd=1, col="red", lty=2)
  dev.off()
}
cat("done ! \n")

#dev.new()


#VIII)Exporting processed data

cat("Exporting processed files...")
writeMSIData(peaksBinned, file = Sample_name, name=paste(Sample_name,"_processed",sep = ""),folder = Output_path,
             mz.type = "32-bit float",intensity.type="32-bit float")
write.table(Intensity_data,file = paste(Output_path,"/",Sample_name,"_intensity_data.txt",sep = ""),sep="\t",quote = F,row.names = F)
write.table(Location_data,file = paste(Output_path,"/",Sample_name,"_location_data.txt",sep = ""),sep="\t",quote = F,row.names = F)
cat("... done ! \n")

