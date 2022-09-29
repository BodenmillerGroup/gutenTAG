#0) Basic functions 

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


#I)Package loading

cat("Loading of the libraries.... ")
suppressMessages(library(Cardinal))
suppressMessages(library(dplyr))
suppressMessages(library(hexbin))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
cat("... done ! \n")


#II)Data and annotation loading

args <- commandArgs(trailingOnly = T)

#Path_to_izml_file = args[1]
Path_to_izml_file = "Desktop/MALDI_Imaging_project/test_file_john/Subset_placenta.imzML"

cat("Reading the .izml and .ibd file  ...")
rawFile <- readMSIData(Path_to_izml_file)
rawFile <- Cardinal::pull(rawFile)
cat("done ! \n")

#Extracting file name

Sample_name = strsplit(Path_to_izml_file,split = "/",fixed = T)[[1]]
Sample_name = Sample_name[length(Sample_name)]
Sample_name = strsplit(Sample_name,split = ".",fixed = T)[[1]]
Sample_name = Sample_name[1]

#Path_to_peak_annotation = args[2]
Path_to_peak_annotation = "Desktop/MALDI_Imaging_project/test_file_john/AmbergenPilotPanel.csv"

cat("Reading the csv annotation file  ...")
peakAnnotation <- read.delim(Path_to_peak_annotation,sep=",")
refList <- peakAnnotation$FeatureMass
cat("done ! \n")

#Path_to_parameter_file = args[3]
Path_to_parameter_file = "Desktop/MALDI_Imaging_project/test_file_john/Parameter_processing.txt"

cat("Loading parameter file...")
Parameters = suppressWarnings(read.table(Path_to_parameter_file,header = F,sep = "="))
Parameters_names = Parameters$V1
Parameters_values = Parameters$V2
names(Parameters_values) = Parameters_names

N_cores = as.numeric(Parameters_values["N_cores"])
peak_mz_tolerance = as.numeric(Parameters_values["peak_mz_tolerance"])
cat("done ! \n")

##III)Parallelisation
cat("Creating parallel environement...")

if (N_cores >= round(0.8*detectCores(),0) ) {
  stop("Not enough threads available. Check Parameter.txt")
}

param <- MulticoreParam(workers = N_cores)
setCardinalBPPARAM(param)

getCardinalBPPARAM()
setCardinalBPPARAM(SerialParam)
cat("... done ! \n")


###IV)Processing

cat("Running smoothing and baseline correction...")
peakPre <- rawFile %>% 
  process(function(x) length(x) * x / sum(x), label="norm", delay=TRUE) %>% 
  smoothSignal(method = "gaussian", plot=FALSE) %>%
  reduceBaseline(method="locmin") %>%
  process()
cat("...done ! \n")

cat("Running peak detection...")
peakProc <- peakPre %>%
  peakPick(method="mad") %>% 
  peakAlign(ref = refList, tolerance = peak_mz_tolerance, units = "mz", method="diff") %>% 
  peakFilter(freq.min = 0.01) %>%
  process()
cat("done ! \n")


cat("Running peak binning and intensity computation....")
binRef <- peakProc@featureData@mz

peaksBinned <- peakPre %>%
  peakBin(ref=binRef, tolerance=0.5, units="mz", type = "area") %>%
  process()
cat("....done ! \n")





###V)Reshaping the data

#Extracting the intensity

Intensity_data = iData(peaksBinned,"intensity")
Intensity_data = t(Intensity_data)

#Annotating the peaks
conserved_peak_mz = peaksBinned@featureData@mz
x = peakAnnotation
rownames(x) = x$FeatureMass
name_conserved_channels = x[as.character(conserved_peak_mz),"Name"]

colnames(Intensity_data) = name_conserved_channels
Intensity_data = as.data.frame(Intensity_data)
#Extracting pixel location
Location_data = as.data.frame(pixelData(peaksBinned))
Location_data = Location_data[,c("x","y")]  

#VI)Quality control

#A)Pixel quality control 

pixelsSummarized <- summarizePixels(peaksBinned, FUN =c(tic="sum", mean = "mean", min = "min", max = "max", sd = "sd", var = "var"), as = "DataFrame")
pixelsSummarized = as.data.frame(pixelsSummarized)
ticSummary = pixelsSummarized$tic

TIC_threshold = Otsu_thresholding(log10(ticSummary))




#VII)Generating QC pdf

Output_path = args[4]
Output_path = "Desktop/MALDI_Imaging_project/Output_directory"
Output_path = paste(Output_path,"/",sep="")

if( !dir.exists(Output_path)) {
  dir.create(Output_path)
}

cat(paste("Creating QC pdf file ..."))

QC_pdf_path = paste("QC_report_",Sample_name,".pdf",sep ="" )
QC_pdf_path = paste(Output_path,QC_pdf_path,sep = "/")

pdf(QC_pdf_path,height = 18,width = 12,useDingbats = F) 
par(las=1,mfrow=c(4,3),mar=c(6,6,6,4))

hist(log10(1+ticSummary), 100, main = "TIC histogram", xlab = "Total Ion Count (Log10)",xaxs='i',yaxs='i',cex.lab=1.3)
abline(v= TIC_threshold, lwd=2, col="red", lty=2)

dev.off()
  
#VIII)Exporting processed data


