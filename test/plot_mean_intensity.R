library(ggplot2)
library(babynames)
library(dplyr)
library(tidyr)
library(reshape2)

# 1. Load inputs ####
path_to_region03 <- "/mnt/msi_volume/processed_files/dilution_experiment/region03/mean_value.txt"
region03 <-  read.table(path_to_region03)

path_to_region04 <- "/mnt/msi_volume/processed_files/dilution_experiment/region04/mean_value.txt"
region04 <-  read.table(path_to_region04)

path_to_region05 <- "/mnt/msi_volume/processed_files/dilution_experiment/region05/mean_value.txt"
region05 <-  read.table(path_to_region05)

path_to_region06 <- "/mnt/msi_volume/processed_files/dilution_experiment/region06/mean_value.txt"
region06 <-  read.table(path_to_region06)

path_to_region07 <- "/mnt/msi_volume/processed_files/dilution_experiment/region07/mean_value.txt"
region07 <-  read.table(path_to_region07)

path_to_region08 <- "/mnt/msi_volume/processed_files/dilution_experiment/region08/mean_value.txt"
region08 <-  read.table(path_to_region08)

path_to_region09 <- "/mnt/msi_volume/processed_files/dilution_experiment/region09/mean_value.txt"
region09 <-  read.table(path_to_region09)

path_to_region10 <- "/mnt/msi_volume/processed_files/dilution_experiment/region10/mean_value.txt"
region10 <-  read.table(path_to_region10)

path_to_region11 <- "/mnt/msi_volume/processed_files/dilution_experiment/region11/mean_value.txt"
region11 <-  read.table(path_to_region11)

path_to_region12 <- "/mnt/msi_volume/processed_files/dilution_experiment/region12/mean_value.txt"
region12 <-  read.table(path_to_region12)

path_to_region13 <- "/mnt/msi_volume/processed_files/dilution_experiment/region13/mean_value.txt"
region13 <-  read.table(path_to_region13)

path_to_region14 <- "/mnt/msi_volume/processed_files/dilution_experiment/region14/mean_value_region14.txt"
region14 <-  read.table(path_to_region14)

path_to_region15 <- "/mnt/msi_volume/processed_files/dilution_experiment/region15/mean_value_region15.txt"
region15 <-  read.table(path_to_region15)

path_to_region16 <- "/mnt/msi_volume/processed_files/dilution_experiment/region16/mean_value_region16.txt"
region16 <-  read.table(path_to_region16)

path_to_region17 <- "/mnt/msi_volume/processed_files/dilution_experiment/region17/mean_value_region17.txt"
region17 <-  read.table(path_to_region17)


# 2. Create df ####
df2 <- cbind(region09, region10, region11, region12, region13)
df3 <- region08
df4 <- cbind(region14, region15, region16, region17)


data3 <- data.frame(t(df3))

# set row and col names for data3
rownames(data3) <- data3[,1]
data3[,1] <- NULL
colnames(data3) <- data3[1,]
data3 <- data3[-1,]

# change columns to numeric

sapply(data3[,1:4], class)
data3[,1:4] <- data.frame(sapply(data3[,1:4], as.numeric))
sapply(data3[,1:4], class)

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

