#Here is the code I used to fix p-values of real data baypass output (qvalues) and test what is resulting when applying different FDR thresholds to the dataset.

library(qvalue)
library(dplyr)
library(tidyr)
library(ggplot2)

df1 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_1_c2-model_summary_contrast.out", h=TRUE)

df2 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_2_c2-model_summary_contrast.out", h=TRUE)

df3 <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-original/original_population_run_3_c2-model_summary_contrast.out", h=TRUE)

pval1 <- df1$log10.1.pval
#range(df1$log10.1.pval)
#plot(df1$log10.1.pval)
pval2 <- df2$log10.1.pval
pval3 <- df3$log10.1.pval


pval1.2 <- 10^(-pval1)
pval2.2 <- 10^(-pval2)
pval3.2 <- 10^(-pval3)

qobj1 <- qvalue(p = pval1.2)
qobj2 <- qvalue(p = pval2.2)
qobj3 <- qvalue(p = pval3.2)

#pdf("q-values.pdf", width=8, height=6)
#par(mfrow = c(3,2))
#plot(qobj1$qvalues)
#plot(qobj2$qvalues)
#plot(qobj3$qvalues)

#now calculate the median and then backtransform to -log10
qvalues1 <- qobj1$qvalues
qvalues2 <- qobj2$qvalues
qvalues3 <- qobj3$qvalues

#combine qvalues
combinedq <- data.frame(q1=qvalues1, q2=qvalues2, q3=qvalues3)

combinedq <- cbind(marker = 1:nrow(combinedq), combinedq)
head(combinedq)

#check for model convergence
q_columns <- grep("^q", names(combinedq), value = TRUE)

# Create an empty matrix to store correlation values
cor_matrix <- matrix(NA, nrow = length(q_columns), ncol = length(q_columns),
                    dimnames = list(q_columns, q_columns))

# Loop through the pairs of columns
for (i in seq_along(q_columns)) {
 for (j in seq_along(q_columns)) {
   cor_matrix[i, j] <- cor(x = combinedq[[q_columns[i]]], 
                           y = combinedq[[q_columns[j]]], 
                           method = "pearson", 
                           use = "complete.obs")
 }
}

# Print the correlation matrix
print(cor_matrix)
#0.68 for all

# median
head(combinedq)
combinedq$median_q <- apply(combinedq[, c("q1", "q2", "q3")], 1, median)
combinedq$log_transformed <- -log10(combinedq$median_q)

#save table
write.csv(
  combinedq,
  file = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-original-output-files/combinedq.csv",
  row.names = FALSE
)


#okay, now I want to check FDR and filter for the ones which are fluctuating
freq_df <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/new-genome-final-vcf/frequency_alt.csv", header=TRUE)

combinedq <- read.table("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-original-output-files/combinedq.csv", header=FALSE)

# Add MRK column
freq_df$MRK <- 1:nrow(freq_df)

# Check
head(freq_df)

#combine two datasets

freq_df_sig <- freq_df %>% 
  left_join(combinedq %>% select (marker, median_q), by= c("MRK"="marker"))

#lets plot all qvalues
#1874896

ggplot(freq_df_sig, aes(x = median_q)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution of Q values",
    x = "Median Q",
    y = "Count"
  ) +
  theme(legend.title = element_blank())


#now lets start exploring the data using different false discovery rates
#fix the dataset
# Keep only the columns you want
freq_df_sig <- freq_df_sig %>%
  select(
    EA_2009_T2.FREQ, EA_2009_T4.FREQ,
    EA_2011_T1.FREQ, EA_2011_T2.FREQ,
    EA_2015_T1.FREQ, EA_2015_T4.FREQ,
    EA_2022_T1.FREQ, EA_2022_T4.FREQ, 
    MRK, median_q
  )

df <- freq_df_sig %>%
  rename(
    `2009.start` = `EA_2009_T2.FREQ`,
    `2009.end` = `EA_2009_T4.FREQ`,
    `2011.start` = `EA_2011_T1.FREQ`,
    `2011.end` = `EA_2011_T2.FREQ`,
    `2015.start` = `EA_2015_T1.FREQ`,
    `2015.end` = `EA_2015_T4.FREQ`,
    `2022.start` = `EA_2022_T1.FREQ`,
    `2022.end` = `EA_2022_T4.FREQ`,
    # Add more as needed
  )

#calculate mean and median change in AF

library(dplyr)

df <- df %>%
  mutate(mean_AF = rowMeans(select(., c("2009.start", "2009.end", "2011.start", "2011.end","2015.start", "2015.end", "2022.start", "2022.end"))),
         na.rm = TRUE)

df <- df %>%
  rowwise() %>%
  mutate(
    median_AF = median(
      c_across(c("2009.start", "2009.end", "2011.start", "2011.end",
                 "2015.start", "2015.end", "2022.start", "2022.end")),
      na.rm = TRUE
    )
  )

#plot
ggplot(df, aes(x = mean_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution of Mean Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#now lets calculate the changes in AF in the whole dataset

df <- df %>%
  mutate(
    diff_2009 = `2009.end`- `2009.start`,
    diff_2011 = `2011.end` - `2011.start`,
    diff_2015 = `2015.end` - `2015.start`,
    diff_2022 = `2022.end` - `2022.start`
  )

df <- df %>%
  mutate(
    diff_2009 = abs(diff_2009),
    diff_2011 = abs(diff_2011),
    diff_2015 = abs(diff_2015),
    diff_2022 = abs(diff_2022)
  )

#mean change in AF
df <- df %>%
  mutate(mean_change_AF = rowMeans(select(., c("diff_2009", "diff_2011", "diff_2015", "diff_2022"))),
         na.rm = TRUE)

#plot
ggplot(df, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())


#okay, lets get the fluctuating snps!

#look at fluctuating ones...d
is_positive <- apply(df, 1, function(row){
  (row["2009.end"] > row["2009.start"]) &
  (row["2011.end"] > row["2011.start"]) &
  (row["2015.end"] > row["2015.start"]) &
  (row["2022.end"] > row["2022.start"])
})

positive_snps <- df[is_positive, ]
#110039

# Subset negative
is_negative <- apply(df, 1, function(row){
  (row["2009.end"] < row["2009.start"]) &
  (row["2011.end"] < row["2011.start"]) &
  (row["2015.end"] < row["2015.start"]) &
  (row["2022.end"] < row["2022.start"])
})

negative_snps <- df[is_negative, ]
#101622

#combine dataset keeping info on positive and negative
positive_snps <- positive_snps %>%
  mutate(type = "positive")

negative_snps <- negative_snps %>%
  mutate(type = "negative")

# Combine the datasets
combined_snps <- bind_rows(positive_snps, negative_snps)
library(matrixStats)

combined_snps <- combined_snps %>%
  mutate(median_change_AF = rowMedians(
    as.matrix(select(., c("diff_2009", "diff_2011", "diff_2015", "diff_2022"))),
    na.rm = TRUE
  ))

#plot mean_change_AF
ggplot(combined_snps, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(combined_snps, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())


#okay, now lets start adding some significance
#############THIS PART HERE HAS A LOT OF REPETITION 
#############I tried for FDR .5, .4, .3, .2, .1 and .05
#select significant markers
sum(combined_snps$median_q < 0.5)
#1343

sig_05 <- subset(combined_snps, median_q <= 0.5)

ggplot(sig_05, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(sig_05, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot snp trajectory

norm_all_sig_05 <- sig_05 %>%
  pivot_longer(
    # Add every column you DON'T want to pivot here:
    cols = -c(MRK, type, median_q, mean_AF, na.rm, starts_with("diff"), contains("change")), 
    names_to = "timepoint",
    values_to = "af"
  )


# -------------------------
# Plot all trajectories
# -------------------------
library(ggplot2)
ggplot(norm_all_sig_05, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Positive SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#calculate mean/median change in AF
mean_val <- mean(sig_05$mean_change_AF, na.rm = TRUE)
mean_val_2 <- mean(sig_05$mean_AF, na.rm = TRUE)
# Calculate Median
median_val <- median(sig_05$mean_change_AF, na.rm = TRUE)
median_val_2 <- median(sig_05$median_change_AF, na.rm = TRUE)

# Print results
print(paste("Mean AF:", mean_val_2))
print(paste("Median AF:", median_val_2))
print(paste("Mean change AF:", mean_val))
print(paste("Median change AF:", median_val))


#now the same for .4 FDR
sig_04 <- subset(combined_snps, median_q <= 0.4)
#436 snps

ggplot(sig_04, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(sig_05, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot snp trajectory

norm_all_sig_04 <- sig_04 %>%
  pivot_longer(
    # Add every column you DON'T want to pivot here:
    cols = -c(MRK, type, median_q, mean_AF, na.rm, starts_with("diff"), contains("change")), 
    names_to = "timepoint",
    values_to = "af"
  )

# -------------------------
# Plot all trajectories
# -------------------------
library(ggplot2)
ggplot(norm_all_sig_04, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#calculate mean/median change in AF
mean_val <- mean(sig_04$mean_change_AF, na.rm = TRUE)
mean_val_2 <- mean(sig_04$mean_AF, na.rm = TRUE)
# Calculate Median
median_val <- median(sig_04$mean_change_AF, na.rm = TRUE)
median_val_2 <- median(sig_04$median_change_AF, na.rm = TRUE)

# Print results
print(paste("Mean AF:", mean_val_2))
print(paste("Median AF:", median_val_2))
print(paste("Mean change AF:", mean_val))
print(paste("Median change AF:", median_val))


#
#now the same for .4 FDR
sig_03 <- subset(combined_snps, median_q <= 0.3)
#183 snps

ggplot(sig_03, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(sig_03, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot snp trajectory

norm_all_sig_03 <- sig_03 %>%
  pivot_longer(
    # Add every column you DON'T want to pivot here:
    cols = -c(MRK, type, median_q, mean_AF, na.rm, starts_with("diff"), contains("change")), 
    names_to = "timepoint",
    values_to = "af"
  )

norm_all_sig_03 <- norm_all_sig_03 %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start", "2009.end",
                                  "2011.start", "2011.end",
                                  "2015.start", "2015.end",
                                  "2022.start", "2022.end"))
  )


# -------------------------
# Plot all trajectories
# -------------------------
library(ggplot2)
ggplot(norm_all_sig_03, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#calculate mean/median change in AF
mean_val <- mean(sig_03$mean_change_AF, na.rm = TRUE)
mean_val_2 <- mean(sig_03$mean_AF, na.rm = TRUE)
# Calculate Median
median_val <- median(sig_03$mean_change_AF, na.rm = TRUE)
median_val_2 <- median(sig_03$median_change_AF, na.rm = TRUE)

# Print results
print(paste("Mean AF:", mean_val_2))
print(paste("Median AF:", median_val_2))
print(paste("Mean change AF:", mean_val))
print(paste("Median change AF:", median_val))

####0.2
sig_02 <- subset(combined_snps, median_q <= 0.2)
#87

ggplot(sig_02, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(sig_02, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot snp trajectory

norm_all_sig_02 <- sig_02 %>%
  pivot_longer(
    # Add every column you DON'T want to pivot here:
    cols = -c(MRK, type, median_q, mean_AF, na.rm, starts_with("diff"), contains("change")), 
    names_to = "timepoint",
    values_to = "af"
  )

norm_all_sig_02 <- norm_all_sig_02 %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start", "2009.end",
                                  "2011.start", "2011.end",
                                  "2015.start", "2015.end",
                                  "2022.start", "2022.end"))
  )


# -------------------------
# Plot all trajectories
# -------------------------
library(ggplot2)
ggplot(norm_all_sig_02, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#calculate mean/median change in AF
mean_val <- mean(sig_02$mean_change_AF, na.rm = TRUE)
mean_val_2 <- mean(sig_02$mean_AF, na.rm = TRUE)
# Calculate Median
median_val <- median(sig_02$mean_change_AF, na.rm = TRUE)
median_val_2 <- median(sig_02$median_change_AF, na.rm = TRUE)

# Print results
print(paste("Mean AF:", mean_val_2))
print(paste("Median AF:", median_val_2))
print(paste("Mean change AF:", mean_val))
print(paste("Median change AF:", median_val))

####0.1
sig_01 <- subset(combined_snps, median_q <= 0.1)
#50

ggplot(sig_01, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(sig_01, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot snp trajectory

norm_all_sig_01 <- sig_01 %>%
  pivot_longer(
    # Add every column you DON'T want to pivot here:
    cols = -c(MRK, type, median_q, mean_AF, na.rm, starts_with("diff"), contains("change")), 
    names_to = "timepoint",
    values_to = "af"
  )

norm_all_sig_01 <- norm_all_sig_01 %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start", "2009.end",
                                  "2011.start", "2011.end",
                                  "2015.start", "2015.end",
                                  "2022.start", "2022.end"))
  )


# -------------------------
# Plot all trajectories
# -------------------------
library(ggplot2)
ggplot(norm_all_sig_01, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#calculate mean/median change in AF
mean_val <- mean(sig_01$mean_change_AF, na.rm = TRUE)
mean_val_2 <- mean(sig_01$mean_AF, na.rm = TRUE)
# Calculate Median
median_val <- median(sig_01$mean_change_AF, na.rm = TRUE)
median_val_2 <- median(sig_01$median_change_AF, na.rm = TRUE)

# Print results
print(paste("Mean AF:", mean_val_2))
print(paste("Median AF:", median_val_2))
print(paste("Mean change AF:", mean_val))
print(paste("Median change AF:", median_val))

####0.05
sig_005 <- subset(combined_snps, median_q <= 0.05)
#32

ggplot(sig_005, aes(x = mean_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Mean AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot median_change_AF
ggplot(sig_005, aes(x = median_change_AF)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution Median AF change Allele Frequencies",
    x = "Mean Allele Frequency",
    y = "Count"
  ) +
  theme(legend.title = element_blank())

#plot snp trajectory

norm_all_sig_005 <- sig_005 %>%
  pivot_longer(
    # Add every column you DON'T want to pivot here:
    cols = -c(MRK, type, median_q, mean_AF, na.rm, starts_with("diff"), contains("change")), 
    names_to = "timepoint",
    values_to = "af"
  )

norm_all_sig_005 <- norm_all_sig_005 %>%
  mutate(
    timepoint = factor(timepoint,
                       levels = c("2009.start", "2009.end",
                                  "2011.start", "2011.end",
                                  "2015.start", "2015.end",
                                  "2022.start", "2022.end"))
  )


# -------------------------
# Plot all trajectories
# -------------------------
library(ggplot2)
ggplot(norm_all_sig_005, aes(x = timepoint, y = af, group = MRK, color = as.factor(MRK))) +
  geom_line(alpha = 0.3) +
#  geom_point(size = 1, position = position_jitter(width = 0.1, height = 0)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "SNP Trajectories (Normalized, Start = 0)",
    x = "Timepoint",
    y = "AF Change"
  ) +
  guides(color = "none")

#calculate mean/median change in AF
mean_val <- mean(sig_005$mean_change_AF, na.rm = TRUE)
mean_val_2 <- mean(sig_005$mean_AF, na.rm = TRUE)
# Calculate Median
median_val <- median(sig_005$mean_change_AF, na.rm = TRUE)
median_val_2 <- median(sig_005$median_change_AF, na.rm = TRUE)

# Print results
print(paste("Mean AF:", mean_val_2))
print(paste("Median AF:", median_val_2))
print(paste("Mean change AF:", mean_val))
print(paste("Median change AF:", median_val))
