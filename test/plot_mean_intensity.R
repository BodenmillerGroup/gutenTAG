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




# 2. Create df ####
df2 <- cbind(region09, region10, region11, region12, region13)
#df["condition"] <- rownames(df)
#df <- as_tibble(df)

data <- data.frame(t(df))
#data2["Moles"] <- c(1e-11, 1e-12, 1e-13, 1e-14, 1e-15)
data["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
data <- as_tibble(data)
#reorder
data <- data[, c("Picomoles", "uncharged_peptide_plus2", "charged_peptide_plus2", "uncharged_peptide_plus1", "charged_peptide_plus1")]

# add id number
id <- as.character(rep(1, times = 5))
id2 <- as.character(rep(2, times = 5))
data["Replicate"] <- id
data2["Replicate"] <- id2


master_df <- bind_rows(data,data2)

# 3. Plot ####

ggplot(master_df, aes(x = Picomoles)) +
  geom_line(aes(y = uncharged_peptide_plus1, group = 1, color = "Uncharged +1")) +
  geom_line(aes(y = charged_peptide_plus1, group = 1, color = "Charged +1")) +
  scale_x_discrete(limits = rev) +
  geom_point(aes(y = uncharged_peptide_plus1)) +
  geom_point(aes(y = charged_peptide_plus1)) +
  scale_fill_brewer() +
  ylab("Intensity")


ggplot(master_df, aes(x = Picomoles)) +
  geom_line(aes(y = uncharged_peptide_plus1, group = Replicate, color = Replicate)) +
  geom_point(aes(y = uncharged_peptide_plus1, group = Replicate, color =  Replicate)) +
  scale_x_discrete(limits = rev) +
  scale_fill_brewer() +
  ggtitle("Uncharged peptide") +
  ylab("Intensity")



