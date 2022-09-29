# Title: Pre-processing Script for optimal sample prep tonsil data
# User: John Abbey
# Start: 20.09.2022
# Latest update: 28.09.2022

# Step 1. Normalisation.

library(Cardinal)
library(readxl)
library(useful)
library(dplyr)
library(ggplot2)
library(hexbin)
library(RColorBrewer)
library(pheatmap)

setwd("/")


# 1. Data loading

# a). MSI data

rawFile <- readMSIData("mnt/msi_volume/imzml-ambergen-5um/placenta.imzML")
#rawFile <- readMSIData("mnt/msi_volume/optimalsampleprep/20220914_partscan_tonsil_optimalsampleprep_5um-region01.imzML")
rawFile <- Cardinal::pull(rawFile)

peakAnnotation <- read.delim("mnt/msi_volume/Panels/AmbergenPilotPanel.csv",sep=",")
refList <- peakAnnotation$FeatureMass


# 2. Set parellelisation and cores

param <- MulticoreParam(workers = 2)
setCardinalBPPARAM(param)

getCardinalBPPARAM()
setCardinalBPPARAM(SerialParam)

# 3. Normalisation, Signal Smoothing and Baseline Reduction

peakPre <- rawFile %>% 
  process(function(x) length(x) * x / sum(x), label="norm", delay=TRUE) %>% 
  smoothSignal(method = "gaussian", plot=FALSE) %>%
  reduceBaseline(method="locmin") %>%
  process()

#Check if processes are complete
mcols(processingData(peakPre))[,-1]

# sanity check
plot(rawFile, pixel =300, main = "Raw Spectrum", strip = FALSE)
plot(peakPre, pixel =300, main = "Pre-processed", strip = FALSE)


# 3. Apply processing pipeline to pre-processed data (picking, alignment to reference, filtering)

peakProc <- peakPre %>%
  peakPick(method="mad") %>% 
  peakAlign(ref = refList, tolerance = 0.2, units = "mz", method="diff") %>% 
  peakFilter(freq.min = 0.01) %>%
  process()


Cardinal::plot(peakProc, pixel =200, main =  "Processed")
plot(peakProc, pixel = 1000)


# 4. Peak binning. 

binRef <- peakProc@featureData@mz

peaksBinned <- peakPre %>%
  peakBin(ref=binRef, tolerance=0.5, units="mz", type = "area") %>%
  process()

Cardinal::plot(peaksBinned, pixel =4000, main = "Binned", strip = FALSE)



# 5. QC

# a. Image all mz values and save as variables

her2_raw <- Cardinal::image(peaksBinned, mz = 1293.75, plusminus=0.2, main = "CD44", strip = FALSE)
her2_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y,  mz = 1293.75, plusminus=0.2, main = "Log CD44", strip = FALSE)
her2_supp <- Cardinal::image(peaksBinned, mz = 1293.75, plusminus=0.2, contrast.enhance="suppression", main = "Hotspot Suppression HER2", strip = FALSE)

#cd20 <- Cardinal::image(peaksBinned, mz = 997.53, plusminus=0.2, contrast.enhance="suppression")
cd31 <- Cardinal::image(peaksBinned, mz = 1125.63, plusminus=0.2, contrast.enhance="suppression")
aasm <- Cardinal::image(peaksBinned, mz = 1194.66, plusminus=0.2, contrast.enhance="suppression")
cd3 <- Cardinal::image(peaksBinned, mz = 1206.72, plusminus=0.2, contrast.enhance="suppression")
cd44 <- Cardinal::image(peaksBinned, mz = 1216.75, plusminus=0.2, contrast.enhance="suppression", main =  "CD44", strip = FALSE)
b2m <- Cardinal::image(peaksBinned, mz = 1244.93, plusminus=0.2, contrast.enhance="suppression")
panCK <- Cardinal::image(peaksBinned, mz = 1288.72, plusminus=0.2, contrast.enhance="suppression") 
her2 <- Cardinal::image(peaksBinned, mz = 1293.75, plusminus=0.2, contrast.enhance="suppression")
hh3 <- Cardinal::image(peaksBinned, mz = 1320.76, plusminus=0.2, contrast.enhance="suppression")
vim <- Cardinal::image(peaksBinned, mz = 1437.80, plusminus=0.2, contrast.enhance="suppression")
#ck5 <- Cardinal::image(peaksBinned, mz = 1494.82, plusminus=0.2, contrast.enhance="suppression")
ck818 <- Cardinal::image(peaksBinned, mz = 1524.83, plusminus=0.2, contrast.enhance="suppression")
epcad <- Cardinal::image(peaksBinned, mz = 1551.85, plusminus=0.2, contrast.enhance="suppression")
fibronectin <- Cardinal::image(peaksBinned, mz = 1581.85, plusminus=0.2, contrast.enhance="suppression", strip = FALSE)
#gata3 <- Cardinal::image(peaksBinned, mz = 1638.87, plusminus=0.2, contrast.enhance="suppression")

# b. Generate Histograms. Un-list and then convert to numeric, then plot histograms and log histograms.

#cd20Num <- cd20 %>% unlist(cd20) %>% as.numeric(cd20) 
#hist(cd20Num)
#hist(log2(1+cd20Num),100)

cd31Num <- cd31 %>% unlist(cd31) %>% as.numeric(cd31) 
hist(cd31Num)
hist(log2(1+cd31Num),100)

aasmNum <- aasm %>% unlist(aasm) %>% as.numeric(aasm)
hist(aasmNum)
hist(log2(1+aasmNum),100)

cd3Num <- cd3 %>% unlist(cd3) %>% as.numeric(cd3) 
hist(cd3Num)
hist(log2(1+cd3Num),100)

cd44Num <- cd44 %>% unlist(cd44) %>% as.numeric(cd44)
hist(cd44Num)
hist(log2(1+cd44Num),100)

b2mNum <- b2m %>% unlist(b2m) %>% as.numeric(b2m)
hist(b2mNum)
hist(log2(1+b2mNum),100)

panCKNum <- panCK %>% unlist(panCK) %>% as.numeric(panCK)
hist(panCKNum)
hist(log2(1+panCKNum),100)

her2Num <- her2 %>% unlist(her2) %>% as.numeric(her2)
hist(her2Num)
hist(log2(1+her2Num),100)

hh3Num <- hh3 %>% unlist(hh3) %>% as.numeric(hh3)
hist(hh3Num)
hist(log2(1+hh3Num),100)

ck5Num <- ck5 %>% unlist(ck5) %>% as.numeric(ck5)
hist(ck5Num)
loghist <- hist(log2(1+ck5Num),100)

ck818Num <- ck818 %>% unlist(ck818) %>% as.numeric(ck818) 
hist(ck818Num)
hist(log2(1+ck818Num),100)

epcadNum <- epcad %>% unlist(epcad) %>% as.numeric(epcad)
hist(epcadNum)
hist(log2(1+epcadNum),100)

fibronectinNum <- fibronectin %>% unlist(fibronectin) %>% as.numeric(fibronectin)
hist(fibronectinNum)
hist(log2(1+fibronectinNum),100)


# c. Log images

#cd20_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 997.53, plusminus=0.2, contrast.enhance="suppression")
cd31_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1125.63, plusminus=0.2, contrast.enhance="suppression")
aasm_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1194.66, plusminus=0.2, contrast.enhance="suppression")
#aasm_log_smooth <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1194.66, plusminus=0.2, contrast.enhance="suppression", smooth.image = "adaptive")
cd3_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1206.72, plusminus=0.2, contrast.enhance="suppression")
cd44_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1216.75, plusminus=0.2, contrast.enhance="suppression")
b2m_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1244.93, plusminus=0.2, contrast.enhance="suppression")
panCK_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1288.72, plusminus=0.2, contrast.enhance="suppression") 
her2_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1293.75, plusminus=0.2, contrast.enhance="suppression")
hh3_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1320.76, plusminus=0.2, contrast.enhance="suppression")
vim_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1437.80, plusminus=0.2, contrast.enhance="suppression")
#ck5_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1494.82, plusminus=0.2, contrast.enhance="suppression")
ck818_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1524.83, plusminus=0.2, contrast.enhance="suppression")
epcad_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1551.85, plusminus=0.2, contrast.enhance="suppression")
fibronectin_log <- Cardinal::image(peaksBinned, log(1+intensity)~ x*y, mz = 1581.85, plusminus=0.2, contrast.enhance="suppression", strip = FALSE)


# 6. Summary Statistics

# summary statistics as dataframe
pixelsSummarized <- summarizePixels(peaksBinned, FUN =c(tic="sum", mean = "mean", min = "min", max = "max", sd = "sd", var = "var"), as = "DataFrame")

ticSummary <- pixelsSummarized$tic
meanSummary <- pixelsSummarized$mean
maxSummary <- pixelsSummarized$max
minSummary <- pixelsSummarized$min
sdSummary <- pixelsSummarized$sd
varSummary <- pixelsSummarized$var

# summary statistics are imaging experiment
pixelsSummarizedImage <- summarizePixels(peaksBinned, FUN =c(tic="sum", mean = "mean", min = "min", max = "max", sd = "sd", var = "var"), as = "ImagingExperiment")

ticSummaryImage <- pixelsSummarizedImage$tic
meanSummaryImage <- pixelsSummarizedImage$mean
maxSummaryImage <- pixelsSummarizedImage$max
minSummaryImage <- pixelsSummarizedImage$min
sdSummaryImage <- pixelsSummarizedImage$sd
varSummaryImage <- pixelsSummarizedImage$var

Cardinal::image(pixelsSummarizedImage, mz = 3, contrast.enhance="suppression")


# 7. Define threshold for selecting pixels from histogram data

# Otsu threshold function from Pierre
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


# create threshold for tic image
TIC_threshold = Otsu_thresholding(log2(1+ticSummary))

# plot histogram for tic and draw thresholding line 
hist(log2(1+ticSummary), 100, main = "TIC histogram", xlab = "Intensity")
abline(v= TIC_threshold, lwd=2, col="red", lty=2)

# select pixels above threshold line (quality pixels)
Selected_pixels = log2(1+ticSummary) > TIC_threshold

# plot log of tic against log of sd 
plot(log2(ticSummary[Selected_pixels]), log2(sdSummary[Selected_pixels]), xlab = "Log TIC", ylab = "Log SD", main = "TIC vs Standard Deviation")
line = lm(log2(sdSummary[Selected_pixels])~log2(ticSummary[Selected_pixels]))
abline(coef(line), lwd=2, col="red", lty=2)
hist(line$residuals,100, main = "Residuals")

logTIC <- log2(ticSummary[Selected_pixels])
logSD <- log2(sdSummary[Selected_pixels])
data.frame(x = logTIC, y = logSD)

# 6. plot as density plot. import colour palette
my_colours = colorRampPalette(rev(brewer.pal(13,'Spectral')))
hexbinplot(logSD ~ logTIC, style = "colorscale", colramp = my_colours, main = "Density Plot")
#line2 <- lm(logSD~logTIC)
#abline(v=10)

# plot markers against each other
hexbinplot(log2(1+cd3Num[Selected_pixels]) ~ log2(1+cd44Num[Selected_pixels]), style = "colorscale", colramp = my_colours, main = "Density Plot",)
plot(log2(1+cd44Num[Selected_pixels]) ~ log2(1+cd44Num[Selected_pixels]))


# Pierre work
Mean_peak_intensity = summarizeFeatures(peaksBinned,FUN = mean)
Mean_peak_intensity = Mean_peak_intensity@featureData
Mean_peak_intensity = as.data.frame(Mean_peak_intensity)
Mean_peak_values <- Mean_peak_intensity$FUN

# add column of marker names to mean intensity df
passed_marker_names <- c("CD31", "Actin-alpha-SM", "CD3", "CD44", "B2M", "PanCK", "HER2", "HH3", "VIM", "CK8/18", "E/P-Cadherin", "Fibronectin")
colour_panel <- c("Red", "Blue", "Green", "Yellow", "Pink", "Brown", "Purple", "Grey", "Black", "Orange", "Turquoise", "Violet")
Mean_peak_intensity <- cbind(passed_marker_names, Mean_peak_intensity)
Mean_peak_intensity <- cbind(Mean_peak_intensity, colour_panel)
`colnames<-`(Mean_peak_intensity, c("Marker", "mz", "FUN", "Colours"))

Var_peak_intensity = summarizeFeatures(peaksBinned, FUN = var)
Var_peak_intensity = Var_peak_intensity@featureData
Var_peak_values <- Var_peak_intensity$FUN


# mean variance plot

#plot(Mean_peak_values, Var_peak_values, log="xy", main = "Mean-Variance Plot", xlab = "Mean Peak Intensity", ylab = "Variance", col = "black", type = "p", pch = 21, bg = colour_panel)
plot(Mean_peak_values, Var_peak_values, log="xy", main = "Mean-Variance Plot", xlab = "Mean Peak Intensity", ylab = "Variance", col = colour_panel, type = "p", pch = 19)
legend("topright", 
       inset = c(0,0.2), 
       cex = 0.7, 
       bty = "n", 
       legend = Mean_peak_intensity$passed_marker_names, 
       text.col = "black",
       col = Mean_peak_intensity$colour_panel, 
       pch = 19)

#plot(log2(Mean_peak_values), log2(Var_peak_values), main = "Mean-Variance Plot", xlab = "Mean Peak Intensity", ylab = "Variance")
text(Mean_peak_values, (Var_peak_values) , labels = Mean_peak_intensity$passed_marker_names, offset = 0, cex = 0.6)
# add line
mean_var_line <- lm(log2(Mean_peak_values) ~ log2(Var_peak_values))
abline(coef(mean_var_line), lwd=2, col="red", lty=2)

                   
