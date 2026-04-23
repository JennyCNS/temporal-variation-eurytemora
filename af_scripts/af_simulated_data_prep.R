# In this script we prepare the output data from baypass with information on the 500 neutral simulations, correcting FDR with qvalues, pulling AFs and combining the results

#first step: compile data and get qvalues

input_dir  <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/c2-model-slim"
output_dir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim"

neutral_files <- list.files(path = input_dir, pattern = "^neutral_.*summary_contrast\\.out$", full.names = TRUE)

# 1. Create an empty list to store results from each file
results_list <- list()

for (i in seq_along(neutral_files)) {
  
  full_path <- neutral_files[i]
  current_filename <- basename(full_path)
  
  # Extract the specific model number
  model_run_val <- as.numeric(gsub(".*neutral_([0-9]+)_.*", "\\1", current_filename))
  
  # Read and calculate q-values
  temp_df <- read.table(full_path, header = TRUE)
  pval1 <- 10^(-temp_df$log10.1.pval)
  qobj <- qvalue(p = pval1)
  
  # 2. Store this specific run in the list
  results_list[[i]] <- data.frame(
    MRK = seq_along(qobj$qvalues),
    qvalue = qobj$qvalues,
    model_run = model_run_val
  )
  
  message(paste("Finished model:", model_run_val))
}

# 3. Combine ALL list elements into one big data frame
final_output <- bind_rows(results_list)

# 4. Save the combined file
write.csv(
  final_output,
  file = file.path(output_dir, "combined_all_markers_qvalues.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Check the results
print(table(final_output$model_run))
```
#5528241 lines

Now, I am getting the AF info from all simulations

```r
library(dplyr)

# Define paths
baypass_dir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/slim-simulations-new-genome/simulations/vcfs/baypass_output"
output_dir  <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

baypass_files <- list.files(baypass_dir, pattern = "^neutral_.*\\.baypass$", full.names = TRUE)

message(paste("Found", length(baypass_files), "files to process..."))

# --- 3. Loop through files ---
for (gfile_path in baypass_files) {
  
  fname <- basename(gfile_path)
  
  # Sharp Regex: Extract only the number between 'neutral_' and the next '_'
  model_num_chr <- gsub("^neutral_([0-9]+)_.*", "\\1", fname)
  model_num <- as.numeric(model_num_chr)
  
  # Read the Baypass raw counts
  # Using read.table which handles variable whitespace/tabs automatically
  g_data <- read.table(gfile_path, header = FALSE)
  
  # --- 4. Calculate AF for ALT ---
  # Pairs: Col1=Ref, Col2=Alt, Col3=Ref, Col4=Alt...
  n_pops <- ncol(g_data) / 2
  af_mat <- matrix(nrow = nrow(g_data), ncol = n_pops)
  
  for (p in 1:n_pops) {
    col_alt <- p * 2      # Even columns
    col_ref <- col_alt - 1 # Odd columns
    
    ref_counts <- g_data[[col_ref]]
    alt_counts <- g_data[[col_alt]]
    total_counts <- ref_counts + alt_counts
    
    # Calculate AF (Alt / Total), handling division by zero
    af_mat[, p] <- ifelse(total_counts > 0, alt_counts / total_counts, 0)
  }
  
  # --- 5. Format Dataframe ---
  df_af <- as.data.frame(af_mat)
  colnames(df_af) <- paste0("Pop_", 1:n_pops)
  
  # Add metadata
  df_af$MRK <- 1:nrow(df_af)
  df_af$model_run <- model_num
  
  # --- 6. Save Result ---
  out_name <- paste0("af_table_run_", model_num, ".csv")
  write.csv(df_af, file = file.path(output_dir, out_name), row.names = FALSE, quote = FALSE)
  
  # Optional: Print the first line of the first pop to double check math in console
  if (model_num == 47) {
     actual_af <- g_data[1,2] / (g_data[1,1] + g_data[1,2])
     message(paste("Check Run 47: Row 1, Pop 1 AF is", round(actual_af, 4)))
  }
}

message("Processing complete.")

```

Here, I join the info on AF and qvalues

```R
library(dplyr)

# Paths
output_dir <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim"

# 1. List all AF files
af_files <- list.files(output_dir, pattern = "af_table_run_.*\\.csv$", full.names = TRUE)

# 2. Read and combine all AF files, forcing model_run to character type
message("Combining 500 AF files...")

all_af_data <- lapply(af_files, function(f) {
  df <- read.csv(f)
  df$model_run <- as.character(df$model_run) # Standardize type here
  return(df)
}) %>% bind_rows()

all_af_data <- all_af_data %>%
  rename(
    start.2009 = Pop_1,
    end.2009   = Pop_2,
    start.2011 = Pop_3,
    end.2011   = Pop_4,
    start.2015 = Pop_5,
    end.2015   = Pop_6,
    start.2022 = Pop_7,
    end.2022   = Pop_8
  )

# Load the q-value dataset from your output directory
qval_path <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/baypass-new-genome/output-c2-slim/combined_all_markers_qvalues.csv"
all_qvals <- read.csv(qval_path)


all_qvals$model_run <- as.character(all_qvals$model_run)

# Now the merge will work
final_master <- all_qvals %>%
  inner_join(all_af_data, by = c("MRK", "model_run"))

# Save under the new name
new_filename <- file.path(output_dir, "combined_all_markers_qvalues_with_AF.csv")

write.csv(
  final_master, 
  file = new_filename, 
  row.names = FALSE, 
  quote = FALSE
)
