library(ggplot2)
library(babynames)
library(dplyr)
library(tidyr)
library(reshape2)

# 1. Load inputs ####

path_to_region01 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region01/mean_value.txt"
region01 <-  read.table(path_to_region01)

path_to_region02 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region02/mean_value.txt"
region02 <-  read.table(path_to_region02)

path_to_region03 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region03/mean_value.txt"
region03 <-  read.table(path_to_region03)

path_to_region04 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region04/mean_value.txt"
region04 <-  read.table(path_to_region04)

path_to_region05 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region05/mean_value.txt"
region05 <-  read.table(path_to_region05)

path_to_region06 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region06/mean_value.txt"
region06 <-  read.table(path_to_region06)

path_to_region07 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region07/mean_value.txt"
region07 <-  read.table(path_to_region07)

path_to_region08 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region08/mean_value.txt"
region08 <-  read.table(path_to_region08)

path_to_region09 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region09/mean_value.txt"
region09 <-  read.table(path_to_region09)

path_to_region10 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region10/mean_value.txt"
region10 <-  read.table(path_to_region10)

path_to_region11 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region11/mean_value.txt"
region11 <-  read.table(path_to_region11)

path_to_region12 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region12/mean_value.txt"
region12 <-  read.table(path_to_region12)

path_to_region13 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region13/mean_value.txt"
region13 <-  read.table(path_to_region13)

path_to_region14 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region14/mean_value.txt"
region14 <-  read.table(path_to_region14)

path_to_region15 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region15/mean_value.txt"
region15 <-  read.table(path_to_region15)

path_to_region16 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region16/mean_value.txt"
region16 <-  read.table(path_to_region16)

path_to_region17 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region17/mean_value.txt"
region17 <-  read.table(path_to_region17)

path_to_region18 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region18/mean_value.txt"
region18 <-  read.table(path_to_region18)

path_to_region19 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region19/mean_value.txt"
region19 <-  read.table(path_to_region19)

path_to_region20 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region20/mean_value.txt"
region20 <-  read.table(path_to_region20)

path_to_region21 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region21/mean_value.txt"
region21 <-  read.table(path_to_region21)

path_to_region22 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region22/mean_value.txt"
region22 <-  read.table(path_to_region22)

path_to_region23 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region23/mean_value.txt"
region23 <-  read.table(path_to_region23)

path_to_region24 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region24/mean_value.txt"
region24 <-  read.table(path_to_region24)

path_to_region25 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region25/mean_value.txt"
region25 <-  read.table(path_to_region25)

path_to_region26 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region26/mean_value.txt"
region26 <-  read.table(path_to_region26)

path_to_region27 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region27/mean_value.txt"
region27 <-  read.table(path_to_region27)

path_to_region28 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region28/mean_value.txt"
region28 <-  read.table(path_to_region28)

path_to_region29 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region29/mean_value.txt"
region29 <-  read.table(path_to_region29)

path_to_region30 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/region30/mean_value.txt"
region30 <-  read.table(path_to_region30)



# 2. Create df ####

series1 <- cbind(region01, region02, region03, region04, region05)
series2 <- cbind(region06, region07, region08, region09, region10)
series3 <- cbind(region11, region12, region13, region14, region15)
series4 <- cbind(region16, region17, region18, region19, region20)
series5 <- cbind(region21, region22, region23, region24, region25)
series6 <- cbind(region26, region27, region28, region29, region30)

# make dataframe nice
curate_dataframe <- function(series){

  rownames(series)<- series[,1]

  colnames(series) <- series[1,]

  series <- series[-1,]

  my_names <- c("Peptide")

  series <- series[,!names(series) %in% my_names]
}

# apply to all series
series1 <- curate_dataframe(series1)
series2 <- curate_dataframe(series2)
series3 <- curate_dataframe(series3)
series4 <- curate_dataframe(series4)
series5 <- curate_dataframe(series5)
series6 <- curate_dataframe(series6)

# save rownames as a variable
rownames1 <- rownames(series1)
rownames2 <- rownames(series2)
rownames3 <- rownames(series3)
rownames4 <- rownames(series4)
rownames5 <- rownames(series5)
rownames6 <- rownames(series6)

# function for ammending column names for each dataframe
#ammendColnames <- function(region = region02){

  colnames(region)<- region[1,]
  region <- region[-1,]

}

# create list of all regions
#all_regions <- list(region01, region02, region03, region04, region05, region06, region07, region08, region09, region10,
#                    region11, region12, region13, region14, region15, region16, region17, region18, region19, region20,
 #                   region21, region22, region23, region24, region25, region26, region27, region28, region29, region30)

# ammend colnames for each dataframe
#all_regions <- lapply(all_regions, ammendColnames)


# change columns to numeric
sapply(series1, class)
series1 <- sapply(series1, MARGIN=2, as.numeric)
sapply(series1, class)
rownames(series1) <- rownames1


t(series1)

# add picmoles
data3["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))

data3 <- as_tibble(data3)
#reorder
data3 <- data3[, c("Picomoles", "uncharged_peptide_plus2", "charged_peptide_plus2", "uncharged_peptide_plus1", "charged_peptide_plus1")]

# add id number
id <- as.character(rep(1, times = 5))
id2 <- as.character(rep(2, times = 5))
id3 <- as.character(rep(3, times = 5))
id4 <-as.character(rep(4, times = 5))
data["Replicate"] <- id
data2["Replicate"] <- id2
data3["Replicate"] <- id3
data4["Replicate"] <- id4

uncharged_df <- bind_rows(data, data2)
uncharged_df[,7] <- NULL

charged_df <- bind_rows(data3, data4)

# 3. Plot ####

# plot uncharged peptide by replicate
uncharged_plot <- ggplot(uncharged_df, aes(x = Picomoles)) +
  geom_line(aes(y = uncharged_peptide_plus1, group = Replicate, color = Replicate)) +
  geom_point(aes(y = uncharged_peptide_plus1, group = Replicate, color =  Replicate)) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer() +
  ggtitle("Uncharged peptide") +
  ylab("Intensity")

ggsave(filename = "uncharged_trend.png",
       plot = uncharged_plot,
       device = "png",
       path = "/mnt/msi_volume/processed_files/dilution_experiment/plots/")

# plot charged peptide by replicate
charged_plot <- ggplot(charged_df, aes(x = Picomoles)) +
  geom_line(aes(y = charged_peptide_plus1, group = Replicate, color = Replicate)) +
  geom_point(aes(y = charged_peptide_plus1, group = Replicate, color =  Replicate)) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer() +
  ggtitle("Charged peptide") +
  ylab("Intensity")

ggsave(filename = "charged_trend.png",
       plot = charged_plot,
       device = "png",
       path = "/mnt/msi_volume/processed_files/dilution_experiment/plots/")

