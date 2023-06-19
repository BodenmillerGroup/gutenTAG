# Import final intensity matrix
path_to_df <- "//mnt/msi_volume/processed_files/JuneforJohn/menzha_20230526_54242_melanoma-50um/_final_intensity_table.txt"
Final_intensity_matrix <- read.delim(path_to_df, sep="\t", header=TRUE)

# Import Correspondance matrix
path_to_correspondence <- "//mnt/msi_volume/processed_files/JuneforJohn/menzha_20230526_54242_melanoma-50um/_export_table.txt"
Correspondence_matrix <- read.delim(path_to_correspondance, sep="\t", header=TRUE)

# Make targeted matrix
Final_intensity_matrix_targeted <- Final_intensity_matrix[,!is.na(Correspondence_matrix$Annotation_peaks)]
colnames(Final_intensity_matrix_targeted) <- Correspondence_matrix$Annotation_peaks[!is.na(Correspondence_matrix$Annotation_peaks)]
Final_intensity_matrix_targeted <- as.data.frame(Final_intensity_matrix_targeted)

# Import Location Data
path_to_coords <- "//mnt/msi_volume/processed_files/JuneforJohn/menzha_20230526_54242_melanoma-50um/_location_data.txt"
Location_data <- read.delim(path_to_coords, sep="\t", header=TRUE)

