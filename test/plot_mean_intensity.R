{library(ggplot2)
library(babynames)
library(dplyr)
library(tidyr)
library(reshape2)
library(kableExtra)
}

# 1. Load inputs ####

path_to_region01 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region01/mean_value.txt"
region01 <-  read.table(path_to_region01)

path_to_region02 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region02/mean_value.txt"
region02 <-  read.table(path_to_region02)

path_to_region03 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region03/mean_value.txt"
region03 <-  read.table(path_to_region03)

path_to_region04 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region04/mean_value.txt"
region04 <-  read.table(path_to_region04)

path_to_region05 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region05/mean_value.txt"
region05 <-  read.table(path_to_region05)

path_to_region06 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region06/mean_value.txt"
region06 <-  read.table(path_to_region06)

path_to_region07 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region07/mean_value.txt"
region07 <-  read.table(path_to_region07)

path_to_region08 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region08/mean_value.txt"
region08 <-  read.table(path_to_region08)

path_to_region09 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region09/mean_value.txt"
region09 <-  read.table(path_to_region09)

path_to_region10 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region10/mean_value.txt"
region10 <-  read.table(path_to_region10)

path_to_region11 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region11/mean_value.txt"
region11 <-  read.table(path_to_region11)

path_to_region12 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region12/mean_value.txt"
region12 <-  read.table(path_to_region12)

path_to_region13 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region13/mean_value.txt"
region13 <-  read.table(path_to_region13)

path_to_region14 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region14/mean_value.txt"
region14 <-  read.table(path_to_region14)

path_to_region15 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region15/mean_value.txt"
region15 <-  read.table(path_to_region15)

path_to_region16 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region16/mean_value.txt"
region16 <-  read.table(path_to_region16)

path_to_region17 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region17/mean_value.txt"
region17 <-  read.table(path_to_region17)

path_to_region18 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region18/mean_value.txt"
region18 <-  read.table(path_to_region18)

path_to_region19 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region19/mean_value.txt"
region19 <-  read.table(path_to_region19)

path_to_region20 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region20/mean_value.txt"
region20 <-  read.table(path_to_region20)

path_to_region21 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region21/mean_value.txt"
region21 <-  read.table(path_to_region21)

path_to_region22 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region22/mean_value.txt"
region22 <-  read.table(path_to_region22)

path_to_region23 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region23/mean_value.txt"
region23 <-  read.table(path_to_region23)

path_to_region24 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region24/mean_value.txt"
region24 <-  read.table(path_to_region24)

path_to_region25 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region25/mean_value.txt"
region25 <-  read.table(path_to_region25)

path_to_region26 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region26/mean_value.txt"
region26 <-  read.table(path_to_region26)

path_to_region27 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region27/mean_value.txt"
region27 <-  read.table(path_to_region27)

path_to_region28 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region28/mean_value.txt"
region28 <-  read.table(path_to_region28)

path_to_region29 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region29/mean_value.txt"
region29 <-  read.table(path_to_region29)

path_to_region30 <- "/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/region30/mean_value.txt"
region30 <-  read.table(path_to_region30)



# 2. Create df ####

series1 <- cbind(region01, region02, region03, region04, region05)
series2 <- cbind(region06, region07, region08, region09, region10)
series3 <- cbind(region11, region12, region13, region14, region15)
series4 <- cbind(region16, region17, region18, region19, region20)
series5 <- cbind(region21, region22, region23, region24, region25)
series6 <- cbind(region26, region27, region28, region29, region30)

# temp adjustment to series5 because region21 didn't have col labels in loaded df
peptide_vector <- series5[,2]
series5 <- cbind(series5, peptide_vector)
series5 <- series5 %>% relocate(peptide_vector)

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

# change columns to numeric
sapply(series5, class)
series1 <- sapply(series1, MARGIN=2, as.numeric)
series2 <- sapply(series2, MARGIN=2, as.numeric)
series3 <- sapply(series3, MARGIN=2, as.numeric)
series4 <- sapply(series4, MARGIN=2, as.numeric)
series5 <- sapply(series5, MARGIN=2, as.numeric)
series6 <- sapply(series6, MARGIN=2, as.numeric)

# re-apply rownames
rownames(series1) <- rownames1
rownames(series2) <- rownames2
rownames(series3) <- rownames3
rownames(series4) <- rownames4
rownames(series5) <- rownames5
rownames(series6) <- rownames6

sapply(series5, class)


# transpose dataframe so that each peptide is a feature
series1 <- data.frame(t(series1))
series2 <- data.frame(t(series2))
series3 <- data.frame(t(series3))
series4 <- data.frame(t(series4))
series5 <- data.frame(t(series5))
series6 <- data.frame(t(series6))

# add picmoles
series1["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
series2["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
series3["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
series4["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
series5["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
series6["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))

# as tibble
series1 <- as_tibble(series1)
series2 <- as_tibble(series2)
series3 <- as_tibble(series3)
series4 <- as_tibble(series4)
series5 <- as_tibble(series5)
series6 <- as_tibble(series6)

# reorder so that Picomoles is first column
series1 <- series1 %>% relocate(Picomoles)
series2 <- series2 %>% relocate(Picomoles)
series3 <- series3 %>% relocate(Picomoles)
series4 <- series4 %>% relocate(Picomoles)
series5 <- series5 %>% relocate(Picomoles)
series6 <- series6 %>% relocate(Picomoles)

# add id number
#rep1 <- as.character(rep(1, times = 5))
#rep2 <- as.character(rep(2, times = 5))
#rep3 <- as.character(rep(3, times = 5))
#series1["Replicate"] <- rep1
#series2["Replicate"] <- rep1
#series3["Replicate"] <- rep2
#series4["Replicate"] <- rep2
#series5["Replicate"] <- rep3
#series6["Replicate"] <- rep3

# Aggregate series 1, 3, 5 and 2, 4, 6
uncharged_replicates <- bind_rows(series1, series3, series5)
charged_replicates <- bind_rows(series2, series4, series6)

# create dataframe of mean value of replicate
numeric_df_charged <- select_if(charged_replicates, is.numeric)
numeric_df_uncharged <- select_if(uncharged_replicates, is.numeric)

# replicates of each concentrations
pmol_10_charged <- numeric_df_charged[c(1,6,11),]
pmol_1_charged <- numeric_df_charged[c(2,7,12),]
pmol_01_charged <- numeric_df_charged[c(3,8,13),]
pmol_001_charged <- numeric_df_charged[c(4,9,14),]
pmol_0001_charged <- numeric_df_charged[c(5,10,15),]

pmol_10_uncharged <- numeric_df_uncharged[c(1,6,11),]
pmol_1_uncharged <- numeric_df_uncharged[c(2,7,12),]
pmol_01_uncharged <- numeric_df_uncharged[c(3,8,13),]
pmol_001_uncharged <- numeric_df_uncharged[c(4,9,14),]
pmol_0001_uncharged <- numeric_df_uncharged[c(5,10,15),]


# calculate mean for each concentration
mean_10pmol_charged <- colMeans(pmol_10_charged)
mean_1pmol_charged <- colMeans(pmol_1_charged)
mean_01pmol_charged <- colMeans(pmol_01_charged)
mean_001pmol_charged <- colMeans(pmol_001_charged)
mean_0001pmol_charged <- colMeans(pmol_0001_charged)

mean_10pmol_uncharged <- colMeans(pmol_10_uncharged)
mean_1pmol_uncharged <- colMeans(pmol_1_uncharged)
mean_01pmol_uncharged <- colMeans(pmol_01_uncharged)
mean_001pmol_uncharged <- colMeans(pmol_001_uncharged)
mean_0001pmol_uncharged <- colMeans(pmol_0001_uncharged)

# calculate standard error for each replicate
sd_10pmol_charged <- sapply(pmol_10_charged, sd)
sd_1pmol_charged <- sapply(pmol_1_charged, sd)
sd_01pmol_charged <- sapply(pmol_01_charged, sd)
sd_001pmol_charged <- sapply(pmol_001_charged, sd)
sd_0001pmol_charged <- sapply(pmol_0001_charged, sd)

sd_10pmol_uncharged <- sapply(pmol_10_uncharged, sd)
sd_1pmol_uncharged <- sapply(pmol_1_uncharged, sd)
sd_01pmol_uncharged <- sapply(pmol_01_uncharged, sd)
sd_001pmol_uncharged <- sapply(pmol_001_uncharged, sd)
sd_0001pmol_uncharged <- sapply(pmol_0001_uncharged, sd)

# create mean df
mean_charged_df <- bind_rows(mean_10pmol_charged, mean_1pmol_charged, mean_01pmol_charged, mean_001pmol_charged, mean_0001pmol_charged)
mean_charged_df["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
mean_charged_df <- mean_charged_df %>% relocate(Picomoles)

mean_uncharged_df <- bind_rows(mean_10pmol_uncharged, mean_1pmol_uncharged, mean_01pmol_uncharged, mean_001pmol_uncharged, mean_0001pmol_uncharged)
mean_uncharged_df["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
mean_uncharged_df <- mean_uncharged_df %>% relocate(Picomoles)


# create sd df
sd_charged_df <- bind_rows(sd_10pmol_charged, sd_1pmol_charged, sd_01pmol_charged, sd_001pmol_charged, sd_0001pmol_charged)
sd_charged_df["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
sd_charged_df <- sd_charged_df %>% relocate(Picomoles)

sd_uncharged_df <- bind_rows(sd_10pmol_uncharged, sd_1pmol_uncharged, sd_01pmol_uncharged, sd_001pmol_uncharged, sd_0001pmol_uncharged)
sd_uncharged_df["Picomoles"] <- as.character(c(10, 1, 0.1, 0.01, 0.001))
sd_uncharged_df <- sd_uncharged_df %>% relocate(Picomoles)

# standard error (sigma/sqrt(length(n)))


# make mean and sd dataframes long by Picomoles
mean_long_charged <- pivot_longer(mean_charged_df, cols = -Picomoles, values_to = "mean", names_to = "Peptide")
sd_long_charged <- pivot_longer(sd_charged_df, cols = -Picomoles, values_to = "sd", names_to = "Peptide")
mean_long_charged["sd"] <- sd_long_charged[,3]

mean_long_uncharged <- pivot_longer(mean_uncharged_df, cols = -Picomoles, values_to = "mean", names_to = "Peptide")
sd_long_uncharged <- pivot_longer(sd_uncharged_df, cols = -Picomoles, values_to = "sd", names_to = "Peptide")
mean_long_uncharged["sd"] <- sd_long_uncharged[,3]



# 3. Plot ####

experiment_dir = paste("/mnt/msi_volume/processed_files/150523_dilution_experiment/CHCA/", sep = "")
print(experiment_dir)
# if directory does not exist, create it
if( !dir.exists(experiment_dir)) {
  dir.create(experiment_dir)
}

plot_path = paste(experiment_dir, "plot",sep = "")
print(plot_path)
# if directory does not exist, create it
if( !dir.exists(plot_path)) {
  dir.create(plot_path)
}

# ggplot with ribbon
charged_plot <- ggplot(data = mean_long_charged, aes(x = Picomoles, group = Peptide)) +
  geom_line(aes(y = mean, color = Peptide)) +
  geom_point(aes(y = mean, color = Peptide)) +
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = Peptide), alpha = .2) +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.2, ) +
  scale_x_discrete(limits = rev) +
  #theme_bw() +
  ggtitle("Charged peptide") +
  ylab("Mean Intensity")

ggsave(filename = "charged_trend.png",
       plot = charged_plot,
       device = "png",
       path = plot_path)

uncharged_plot <- ggplot(data = mean_long_uncharged, aes(x = Picomoles, group = Peptide)) +
  geom_line(aes(y = mean, color = Peptide)) +
  geom_point(aes(y = mean, color = Peptide)) +
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = Peptide), alpha = .2) +
  scale_x_discrete(limits = rev) +
  #theme_bw() +
  ggtitle("Uncharged peptide") +
  ylab("Mean Intensity")

ggsave(filename = "uncharged_trend.png",
       plot = uncharged_plot,
       device = "png",
       path = plot_path)


# save mean charged intensity table

mean_charged_table_dhb <- mean_charged_df %>%
  kbl() %>%
  kable_styling(latex_options = "condensed", full_width = FALSE) %>%
  save_kable("charged_mean_df.png", zoom = 10)

mean_uncharged_table_dhb <- mean_uncharged_df %>%
  kbl() %>%
  kable_styling(latex_options = "condensed", full_width = FALSE) %>%
  save_kable("uncharged_mean_df.png", zoom = 10)


# save table
table_path = paste(experiment_dir, "table",sep = "")
print(table_path)
# if directory does not exist, create it
if( !dir.exists(table_path)) {
  dir.create(table_path)
}

print(paste(experiment_dir, "/", "mean_charged.txt", sep = ""))
write.table(mean_charged_df, file = paste(experiment_dir, "/", "charged_mean_df.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)
write.table(mean_uncharged_df, file = paste(experiment_dir, "/", "uncharged_mean_df.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)



