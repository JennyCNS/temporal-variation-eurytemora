#In this script, I explore the dataset ID fluctuating snps and looking at how the changes in number and AF of fluctuating snps varies according to the FDR (.5, .4, .3, .2, .1, .05, .01, .001)

library(dplyr)
library(tidyr)
library(ggplot2)

# Load data
file_path <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim/combined_all_markers_qvalues_with_AF.csv"
af_wide <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

#starting points
#I will not fix everything according to 2009. Better to avoid any misunderstandings.

af_wide <- af_wide %>%
  mutate(
    diff_2009 = `end.2009`- `start.2009`,
    diff_2011 = `end.2011` - `start.2011`,
    diff_2015 = `end.2015` - `start.2015`,
    diff_2022 = `end.2022` - `start.2022`
  )

# Apply renaming to your data frame (adjust 'final_master' to your object name)
af_wide <- af_wide %>%
  rename(
    `2009.start` = start.2009,
    `2009.end`   = end.2009,
    `2011.start` = start.2011,
    `2011.end`   = end.2011,
    `2015.start` = start.2015,
    `2015.end`   = end.2015,
    `2022.start` = start.2022,
    `2022.end`   = end.2022
  )

hist(af_wide$qvalue, 
     breaks = 50, 
     col = "skyblue", 
     main = "Distribution of q-values for all SNPs", 
     xlab = "q-value", 
     ylab = "Frequency")

# Count SNPs per model
snps_per_model <- af_wide %>%
  group_by(model_run) %>%
  summarise(n_snps = n())

# Create the histogram
ggplot(snps_per_model, aes(x = n_snps)) +
  geom_histogram(bins = 30, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

combined_snps_final <- combined_snps %>%
  rowwise() %>%
  mutate(mean_diff = mean(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()
combined_snps_final <- combined_snps_final %>%
  rowwise() %>%
  mutate(median_diff = median(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

#differences 
#mean and median change per model
library(stringr)
library(tidyverse)
af_wide_stats <- af_wide %>%
  rowwise() %>%
  mutate(mean_diff = mean(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

af_wide_stats <- af_wide_stats %>%
  rowwise() %>%
  mutate(median_diff = median(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()
model_results <- af_wide_stats %>%
#Group by that ID
  group_by(model_run) %>%
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(mean_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

head(model_results)

#plot
ggplot(model_results, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")
#plot2
ggplot(model_results, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")


# Identify positive and negative and positive SNPs
# -------------------------
# Identify positive SNPs: end > start for all years
is_negative <- apply(af_wide, 1, function(row){
  (row["2009.end"] < row["2009.start"]) &
  (row["2011.end"] < row["2011.start"]) &
  (row["2015.end"] < row["2015.start"]) &
  (row["2022.end"] < row["2022.start"])
})

is_positive <- apply(af_wide, 1, function(row){
  (row["2009.end"] > row["2009.start"]) &
  (row["2011.end"] > row["2011.start"]) &
  (row["2015.end"] > row["2015.start"]) &
  (row["2022.end"] > row["2022.start"])
})

positive_snps <- af_wide[is_positive, ]
negative_snps <- af_wide[is_negative, ]

nrow(positive_snps)
nrow(negative_snps)

#197853
#199339

#calculate how many significant snps per simulation
combined_snps <- rbind(positive_snps, negative_snps)
library(dplyr)
library(stringr)

combined_snps <- combined_snps %>%
  mutate(marker_model = paste0(MRK, "_", model_run)) %>%
  select(marker_model, everything()) # This moves the new ID to the first column

# 2. Save the final dataset
write.csv(
  combined_snps, 
  file = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim/combined_all_fluctuating_markers_qvalues_with_AF.csv", 
  row.names = FALSE, 
  quote = FALSE
)

#load data
file_path <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim/combined_all_fluctuating_markers_qvalues_with_AF.csv"
combined_snps <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)

#histogram of fluctuating snps

hist(combined_snps$qvalue, 
     breaks = 50, 
     col = "skyblue", 
     main = "Distribution of q-values for Significant SNPs", 
     xlab = "q-value", 
     ylab = "Frequency")

# Add a density line to see the shape
lines(density(combined_snps$qvalue), col = "red", lwd = 2)

# Count SNPs per model
snps_per_model <- combined_snps %>%
  group_by(model_run) %>%
  summarise(n_snps = n())

# Create the histogram
ggplot(snps_per_model, aes(x = n_snps)) +
  geom_histogram(bins = 30, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

####################################################
#test
#plotting 100 random snps from each dataset to check if they doo what I expect

# 1. Sample 100 random markers from each group BEFORE pivoting
set.seed(123) # Use a seed so you get the same 'random' lines every time
pos_subset <- positive_snps %>% sample_n(100)
neg_subset <- negative_snps %>% sample_n(100)

# 1. Transform and ensure the grouping ID exists
norm_long_all <- bind_rows(pos_subset, neg_subset) %>%
  # Create the unique ID if it's missing or got dropped
  mutate(marker_model = paste0(MRK, "_", model_run)) %>% 
  pivot_longer(
    cols = contains("20"), 
    names_to = "timepoint",
    values_to = "AF"
  ) %>%
  mutate(timepoint = factor(timepoint, 
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")))%>%
  filter(!is.na(AF), !is.na(timepoint))
# 2. Plot using marker_model as the group
ggplot(norm_long_all, aes(x = timepoint, y = AF, group = marker_model, color = direction)) +
  geom_line(alpha = 0.3) +
  facet_wrap(~direction) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Trajectories", y = "Allele Frequency")

#okay, they do!

###############skip
#in case I would want to plot the whole thing
norm_long_positive <- positive_snps %>%
  pivot_longer(-marker_model, names_to = "timepoint", values_to = "af") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")),
         trajectory = "Positive")

norm_long_negative <- negative_snps %>%
  pivot_longer(-marker_model, names_to = "timepoint", values_to = "af") %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("2009.start","2009.end",
                                       "2011.start","2011.end",
                                       "2015.start","2015.end",
                                       "2022.start","2022.end")),
         trajectory = "Negative")

norm_long_all <- bind_rows(norm_long_positive, norm_long_negative)

# -------------------------
# Check
# -------------------------
any(is.na(norm_long_all$timepoint))  # should now be FALSE
range(norm_long_all$af)  # all positive values

# -------------------------
# Plot trajectories
# -------------------------
p_traj <- ggplot(norm_long_all, aes(x = timepoint, y = af, group = marker_model, color = trajectory)) +
  geom_line(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Trajectories",
       x = "Timepoint", y = "AF Change")
p_traj


#######################################


# Compute mean/median per fluctuating SNP
# -------------------------

#okay first I will estimate the change of the fluctuating snps per year

combined_snps <- combined_snps %>%
  mutate(
    diff_2009 = abs(diff_2009),
    diff_2011 = abs(diff_2011),
    diff_2015 = abs(diff_2015),
    diff_2022 = abs(diff_2022)
  )

combined_snps_final <- combined_snps %>%
  rowwise() %>%
  mutate(mean_diff = mean(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()
combined_snps_final <- combined_snps_final %>%
  rowwise() %>%
  mutate(median_diff = median(c(diff_2009, diff_2011, diff_2015, diff_2022), na.rm = TRUE)) %>%
  ungroup()

ggplot(combined_snps_final, aes(x = median_diff)) +
  geom_histogram(fill = "steelblue", color = "white", alpha = 0.6, bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median AF Change per SNP", 
       x = "Median AF Change", 
       y = "Count")

ggplot(combined_snps_final, aes(x = mean_diff)) +
  geom_histogram(fill = "steelblue", color = "white", alpha = 0.6, bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean AF Change per SNP", 
       x = "Mean AF Change", 
       y = "Count")

overall_mean <- mean(combined_snps_final$median_diff, na.rm = TRUE)
# 0.05390716
overall_median <- median(combined_snps_final$median_diff, na.rm = TRUE)
#0.04411765

#mean and median change per model
library(stringr)
library(tidyverse)
model_results <- combined_snps_final %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )
head(model_results)

ggplot(model_results, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Median Change",
       y = "Frequency (Number of Models)")

ggplot(model_results, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")


#okay, now lets start adding some significance filtering
#We will try out different FDRs

#first 0.5

model_results_05 <- combined_snps_final %>%
  filter(qvalue < 0.5) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_05)
#489

ggplot(model_results_05, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_05, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_05, aes(x = num_snps)) +
  geom_histogram(bins = 30, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_04$num_snps)
#16133

#now the same for 0.4

model_results_04 <- combined_snps_final %>%
  filter(qvalue < 0.4) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_04)
#485

ggplot(model_results_04, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_04, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_04, aes(x = num_snps)) +
  geom_histogram(bins = 30, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_04$num_snps)
#7010

#now the same for 0.3

model_results_03 <- combined_snps_final %>%
  filter(qvalue < 0.3) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_03)
#470

ggplot(model_results_03, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_03, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_03, aes(x = num_snps)) +
  geom_histogram(bins = 30, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_03$num_snps)
#3430

#we are getting rid of the extreme models!

#now the same for 0.2
model_results_02 <- combined_snps_final %>%
  filter(qvalue < 0.2) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_02)
#443

ggplot(model_results_02, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_02, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_02, aes(x = num_snps)) +
  geom_histogram(bins = 20, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )
sum(model_results_02$num_snps)
#1760

#now the same for 0.1
model_results_01 <- combined_snps_final %>%
  filter(qvalue < 0.1) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_01)
#383

ggplot(model_results_01, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_01, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_01, aes(x = num_snps)) +
  geom_histogram(bins = 10, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_01$num_snps)
#914



#now the same for 0.05
model_results_005 <- combined_snps_final %>%
  filter(qvalue < 0.05) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_005)
#343

ggplot(model_results_005, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_005, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_005, aes(x = num_snps)) +
  geom_histogram(bins = 10, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_005$num_snps)
#683

#now the same for 0.01
model_results_001 <- combined_snps_final %>%
  filter(qvalue < 0.001) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_001)
#235

ggplot(model_results_001, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_001, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_001, aes(x = num_snps)) +
  geom_histogram(bins = 5, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_001$num_snps)
#321


#now the same for 0.001
model_results_0001 <- combined_snps_final %>%
  filter(qvalue < 0.001) %>%
#Extract the model ID (the part after the underscore)
#Group by that ID
  group_by(model_run) %>%
  
#Calculate mean and median of your change column for that group
  summarize(
    model_mean_change = mean(median_diff, na.rm = TRUE),
    model_median_change = median(median_diff, na.rm = TRUE),
    num_snps = n() # Optional: keeps track of how many SNPs are in each model
  )

nrow(model_results_0001)
#277

ggplot(model_results_0001, aes(x = model_mean_change)) +
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Mean Change Across Models",
       x = "Average Mean Change",
       y = "Frequency (Number of Models)")

ggplot(model_results_0001, aes(x = model_median_change)) + # Changed to median
  geom_histogram(fill = "steelblue", color = "white", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Median Change Across Models",
       x = "Model Median Change",
       y = "Frequency (Number of Models)")

# Create the histogram
ggplot(model_results_0001, aes(x = num_snps)) +
  geom_histogram(bins = 10, fill = "forestgreen", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of SNP Counts per Simulation",
    x = "Number of SNPs",
    y = "Number of Models (Simulations)"
  )

sum(model_results_0001$num_snps)
#417
