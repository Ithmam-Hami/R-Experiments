# Install the readx1 package and call it

install.packages("readxl")
library(readxl)

# Define the excel file path on computer. Make sure to replace backslashes with double-backslashes or front-slashes

file_path <- "D:\\Meller Lab\\Luciferase materials\\Luciferase Analysis\\Dual Luciferase Reporter Assay System ms 2025.01.15 14_43_14.xlsx"

# Read the excel file and convert into a dataframe

results <- read_excel(file_path, sheet = "Analysis", col_names = TRUE, col_types = rep("text", ncol(read_excel(file_path, sheet = "Analysis"))))

# Calculate the average of Fluc blanks
fluc_blank_avg <- mean(as.numeric(results[results$Assay == "Fluc", c("Blank 1", "Blank 2")]), na.rm = TRUE)

# Subtract the average Fluc blank from all Fluc values

  # Only consider the fluc values
fluc_values <- setdiff(names(results), c("Assay", "Blank 1", "Blank 2"))

# Make a new dataframe for the processed data

corrected_data <- results

# Check the format of the data
str(results)
str(corrected_data)

# If not numeric, convert them
results[fluc_values] <- lapply(results[fluc_values], as.numeric)
corrected_data[fluc_values] <- lapply(corrected_data[fluc_values], as.numeric)
  
  # Subtraction
corrected_data[results$Assay == "Fluc", fluc_values] <- corrected_data[results$Assay == "Fluc", fluc_values] - fluc_blank_avg

# Calculate the average of Rluc blanks
rluc_blank_avg <- mean(as.numeric(results[results$Assay == "Rluc", c("Blank 1", "Blank 2")]), na.rm = TRUE)

# Subtract the average Rluc blank from all Rluc values

  # Only consider the Rluc values
rluc_values <- setdiff(names(results), c("Assay", "Blank 1", "Blank 2"))

  # Subtraction
corrected_data[results$Assay == "Rluc", rluc_values] <- corrected_data[results$Assay == "Rluc", rluc_values] - rluc_blank_avg

# Calculate Fluc/Rluc ratio

 # Consider the needed values only
ratio_values <- setdiff(names(corrected_data), c("Assay", "Blank 1", "Blank 2"))

 # Select the fluc and rluc data only
fluc_values <- corrected_data[corrected_data$Assay == "Fluc", ratio_values]
rluc_values <- corrected_data[corrected_data$Assay == "Rluc", ratio_values]
fluc_rluc_ratios <- fluc_values / rluc_values

# Average the control (no recruiting element or knockdown) 
control_1_avg <- mean(as.numeric(fluc_rluc_ratios[c("Control 1A", "Control 1B")]), na.rm = TRUE)

# Fluc/Rluc ratios over the control
ratio_over_control <- fluc_rluc_ratios / control_1_avg

# Average the Replicates 

install.packages("dplyr")
install.packages("tidyr")
library(dplyr)
library(tidyr)

# Convert the dataframe into long format
ratio_long <- ratio_over_control %>%
  pivot_longer(cols = everything(),                 # Take all the columns into consideration
               names_to = "Replicate",              # Column names are now merged into a new column
               values_to = "Value") %>%             # Values are now merged into a new column
  mutate(Group = gsub("(A|B|1|2)$", "", Replicate)) # Extracts the shared prefixes and eliminates replicate-specific suffixes in the Group column

# Average
averaged_replicates <- ratio_long %>%
  group_by(Group) %>%                               # Group by 'Group' (e.g., Control 1, Control 2, etc.)
  summarize(Average = mean(Value, na.rm = TRUE))    # Calculate the mean for each group

# Data Representation

install.packages("ggplot2")
library(ggplot2)

#Bar Plot

ggplot(averaged_replicates, aes(x= Group, y = Average, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6, fill = "skyblue") +
  labs(
    title = "Dual Luciferase Assay for SMC1 Knockdown",
    x = "Genotypes",
    y = "Luciferase Activity"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Point Plot

ggplot(averaged_replicates, aes(x = Group, y = Average)) +
  geom_point(size = 4, color = "blue") +
  geom_text(aes(label = round(Average, 2)), vjust = -0.5, size = 3) +
  labs(
    title = "Dual Luciferase Assay for SMC1 Knockdown",
    x = "Genotypes",
    y = "Luciferase Activity"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Lollipop Plot

ggplot(averaged_replicates, aes(x = Group, y = Average)) +
  geom_segment(aes(x = Group, xend = Group, y = 0, yend = Average), color = "gray") +
  geom_point(size = 4, color = "red") +
  labs(
    title = "Dual Luciferase Assay for SMC1 Knockdown",
    x = "Genotypes",
    y = "Luciferase Activity"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
