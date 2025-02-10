# Install and load necessary packages
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")

# Call the libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define file path. Make sure to replace backslashes with double-backslashes or front-slashes
file_path <- "D:\\Meller Lab\\Luciferase materials\\Luciferase Analysis\\Analysis 2-7-25\\luciferase 2-7 simplified.xlsx"

# Read the excel file and convert into a dataframe
results <- read_excel(file_path, sheet = "Sheet1", col_names = TRUE)

# Convert all numeric columns to numeric type
results[-1] <- lapply(results[-1], as.numeric)

# Validation of data type 
print(names(results))  # List column names
print(str(results))  # Check data structure

# Clean column names (if necessary)
names(results) <- gsub("\\.\\.\\..*", "", names(results))

# Calculate the average of Fluc blanks (Normally this should work, but if it does not, follow the next following solves)
fluc_blank_avg <- mean(results[results$Assay == "Fluc", c("Blank", "Buffer", "Buffer.1")], na.rm = TRUE)

  # Manually convert all the required columns into numeric data types and run the previous code 
results$Blank <- as.numeric(results$Blank)
results$Buffer <- as.numeric(results$Buffer)
results$Buffer.1 <- as.numeric(results$Buffer.1)

  # If it still doesn't work, check your rows and subset the columns out for a separate calculation 
print(unique(results$Assay))  # Check unique values in the Assay column
fluc_subset <- results[results$Assay == "Fluc", c("Blank", "Buffer", "Buffer.1")]
print(str(fluc_subset)) # Validate the data type. If it shows "num" values, move on to the mean calculation.

fluc_blank_avg <- mean(unlist(fluc_subset), na.rm = TRUE) # worked for me

  # I added the unlist() function to the original code and it worked
fluc_blank_avg <- mean(unlist(results[results$Assay == "Fluc", c("Blank", "Buffer", "Buffer.1")]), na.rm = TRUE)


# Subtract the average Fluc blank from all Fluc values
fluc_values <- setdiff(names(results), c("Assay", "Blank", "Buffer", "Buffer.1")) # Only consider the fluc values
corrected_data <- results # Make a new dataframe for all the corrected Fluc (& Rluc) values 
corrected_data[results$Assay == "Fluc", fluc_values] <- corrected_data[results$Assay == "Fluc", fluc_values] - fluc_blank_avg

# Calculate the average of Rluc blanks
rluc_blank_avg <- mean(unlist(results[results$Assay == "Rluc", c("Blank", "Buffer", "Buffer.1")]), na.rm = TRUE)

# Subtract the average Rluc blank from all Rluc values
rluc_values <- setdiff(names(results), c("Assay", "Blank", "Buffer", "Buffer.1"))
corrected_data[results$Assay == "Rluc", rluc_values] <- corrected_data[results$Assay == "Rluc", rluc_values] - rluc_blank_avg

# This step is only required if you have columns with the same name for both technical replicates (e.g. A1, A1 etc) 
names(corrected_data) <- make.names(names(corrected_data), unique = TRUE)
names(corrected_data) <- gsub("^X", "", names(corrected_data)) # There was an unexpected 'X' before my column names that started with a number, so this code helps cull that out

# Calculate Fluc/Rluc ratios
ratio_values <- setdiff(names(corrected_data), c("Assay", "Blank", "Buffer", "Buffer.1")) # Set aside the blanks and controls
fluc_values <- corrected_data[corrected_data$Assay == "Fluc", ratio_values] 
rluc_values <- corrected_data[corrected_data$Assay == "Rluc", ratio_values]
fluc_rluc_ratios <- fluc_values / rluc_values

# Normalize to control (Control 1a, 1b, 1c - my biological replicates)   

  # Compute mean of all columns starting with "1a", "1b", or "1c". Add the grep() function ONLY IF you have technical replicates that have the same prefixes i.e. in my case - 1a, 1a.1 etc.)
control_avg <- mean(as.numeric(unlist(fluc_rluc_ratios[, grep("^1[abc]", names(fluc_rluc_ratios))])), na.rm = TRUE)

  # Normalize all values by the computed control average
ratio_over_control <- fluc_rluc_ratios / control_avg

# Convert the dataframe into long format for grouping 

  # without labelling the genotypes
ratio_long <- ratio_over_control %>%
  pivot_longer(cols = everything(), names_to = "Replicate", values_to = "Value") %>%
  mutate(
    BioRep = gsub("\\..*", "", Replicate),  # Extract biological replicate name (remove technical suffix)
    Group = case_when(
      grepl("^[1-4][a-c]", BioRep) ~ "Control",
      grepl("^I[1-3]|^J[1-3]|^K[1-3]|^L[1-3]", BioRep) ~ "SMC1 KD",
      grepl("^E[1-3]|^F[1-3]|^G[1-3]|^H[1-3]", BioRep) ~ "CAPH2 KD",
      TRUE ~ "Unknown"
    )
  )

   
  # Labelling the genotypes
ratio_long_2 <- ratio_over_control %>%
  pivot_longer(cols = everything(), names_to = "Replicate", values_to = "Value") %>%
  mutate(
    BioRep = gsub("\\..*", "", Replicate),  # Remove technical replicate suffix
    Genotype = case_when(
      grepl("^1[a-c]|^L[1-3]|^H[1-3]", BioRep) ~ "Fluc",
      grepl("^2[a-c]|^I[1-3]|^E[1-3]", BioRep) ~ "Fluc + 1.688",
      grepl("^3[a-c]|^K[1-3]|^F[1-3]", BioRep) ~ "Fluc + roX1",  
      grepl("^4[a-c]|^J[1-3]|^G[1-3]", BioRep) ~ "Fluc + 1.688 + roX1",  
      TRUE ~ "Unknown"
    ),
    Condition = case_when(
      grepl("^1[a-c]|^2[a-c]|^3[a-c]|^4[a-c]", BioRep) ~ "Control",
      grepl("^I[1-3]|^J[1-3]|^K[1-3]|^L[1-3]", BioRep) ~ "SMC1 KD",
      grepl("^E[1-3]|^F[1-3]|^G[1-3]|^H[1-3]", BioRep) ~ "CAPH2 KD",
      TRUE ~ "Unknown"
    )
  )

# Visualize the technical replicates (if required)

ggplot(ratio_long_2, aes(x = Genotype, y = Value, color = BioRep)) +
  geom_point(size = 3, alpha = 0.7, position = position_jitter(width = 0.3)) +
  facet_wrap(~ Condition) +  # Separate by Control, SMC1 KD, CAPH2 KD
  labs(title = "Dual Luciferase Assay (Technical Replicates by Genotype)",
       x = "Genotype", y = "Luciferase Activity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  # Removes legend
)


# Average technical replicates but keep biological replicates separate
ratio_avg_tech <- ratio_long_2 %>%
  group_by(BioRep, Genotype, Condition) %>%
  summarise(Average_Tech = mean(Value, na.rm = TRUE), .groups = "drop")



# Data Visualization - Scatter Plot
ggplot(ratio_avg_tech, aes(x = Genotype, y = Average_Tech, color = BioRep)) +
  geom_point(size = 4, alpha = 0.8, position = position_jitter(width = 0.15, height = 0)) +
  facet_wrap(~ Condition) +  
  labs(title = "Biological Replicates in Dual Luciferase Assay",
       x = "Genotype", y = "Luciferase Activity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  


