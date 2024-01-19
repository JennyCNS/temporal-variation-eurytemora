library(devtools)
setwd("~/Jenny-Eurytemora_affinis/baypass")
#install_github("jdstorey/qvalue")
browseVignettes(package = "qvalue")
library(qvalue)
library(dplyr)
library(purrr)


#set working directory
setwd("~/Jenny-Eurytemora_affinis/baypass/AUX-model")
#auxilary model convergence
df1 <- read.table("aux-model_summary_betai.out", header= TRUE)
df2 <- read.table("aux-model2_summary_betai.out", header= TRUE)
df3 <- read.table("aux-model-3_summary_betai.out", header= TRUE)
df4 <- read.table("aux-model-4_summary_betai.out", header= TRUE)
df5 <- read.table("aux-model-5_summary_betai.out", header= TRUE)

head(df1)
head(df2)

par(2,1)
plot(df1$BF.dB)
plot(df2$BF.dB)

#n of significant SNPs (BF > 20)

bf1 <- df1[df1$BF.dB >= 20, ]
nrow(bf1)
#3270
bf2 <- df2[df2$BF.dB >= 20, ]
nrow(bf2)
#3227
bf3 <- df3[df3$BF.dB >= 20, ]
nrow(bf3)
#3354

#list_of_dfs <- list(bf1, bf2, bf3)
#overlap <- reduce(list_of_dfs, inner_join, by = "MRK")
#333
#write.table(overlap, "overlap.txt", sep = "\t", quote = FALSE, row.names = FALSE)

df_list <- list(df1, df2, df3, df4, df5)

# Use lapply to extract the BF.dB column from each data frame
result_df <- data.frame(lapply(df_list, function(df) df$BF.dB))

# Rename the columns if needed
colnames(result_df) <- paste0("bf", 1:5)

bf_combined <- cbind(marker = 1:nrow(result_df), result_df)

#check for convergence
bf_columns <- grep("^bf", names(bf_combined), value = TRUE)

# Create an empty matrix to store correlation values
cor_matrix <- matrix(NA, nrow = length(bf_columns), ncol = length(bf_columns),
                     dimnames = list(bf_columns, bf_columns))

# Loop through the pairs of columns
for (i in seq_along(bf_columns)) {
  for (j in seq_along(bf_columns)) {
    cor_matrix[i, j] <- cor(x = bf_combined[[bf_columns[i]]], 
                            y = bf_combined[[bf_columns[j]]], 
                            method = "pearson", 
                            use = "complete.obs")
  }
}

# Print the correlation matrix
print(cor_matrix)

ggscatter(data = bf_combined, x = "bf1", y = "bf2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "BF-model1", ylab = "BF-model2")

#calculate median for all of them
head(bf_combined)
bf_combined$median_bf <- apply(bf_combined[, c("bf1", "bf2", "bf3", "bf4", "bf5")], 1, median)

plot(bf_combined$median_bf)

#n of significant SNPs (BF > 20)

bf_sig <- bf_combined[bf_combined$median_bf >= 20, ]
nrow(bf_sig)
#951
bf_sig2 <- bf_combined[bf_combined$median_bf >= 15, ]
nrow(bf_sig2)
#2161




#######################################################
#now I'll do the same for c2

setwd("~/Jenny-Eurytemora_affinis/baypass/C2-model")

df1 <- read.table("c2-model_summary_contrast.out", header= TRUE)
df2 <- read.table("c2-model2_summary_contrast.out", header= TRUE)
df3 <- read.table("c2.3-model_summary_contrast.out", header= TRUE)
df4 <- read.table("c2.4-model_summary_contrast.out", header= TRUE)
df5 <- read.table("c2.5-model_summary_contrast.out", header= TRUE)
head(df1)

pval1 <- df1$log10.1.pval
range(df1$log10.1.pval)
plot(data$log10.1.pval)

pval2 <- df2$log10.1.pval
pval3 <- df3$log10.1.pval
pval4 <- df4$log10.1.pval
pval5 <- df5$log10.1.pval
head(pval2)
range(pval2)
pval1.2 <- 10^(-pval1)
pval2.2 <- 10^(-pval2)
pval3.2 <- 10^(-pval3)
pval4.2 <- 10^(-pval4)
pval5.2 <- 10^(-pval5)
range(pval2.2)
head(pval1.2)

pdf("Pvalues-5runs.pdf", width = 8, height = 6)
par(mfrow = c(3,2))
hist(pval1.2)
hist(pval2.2)
hist(pval3.2)
hist(pval4.2)
hist(pval5.2)
dev.off()

qobj <- qvalue(p = pval1.2)
qobj2 <- qvalue(p = pval2.2)
qobj3 <- qvalue(p = pval3.2)
qobj4 <- qvalue(p = pval4.2)
qobj5 <- qvalue(p = pval5.2)
#check if they're alright
plot(qobj2$qvalue)

pdf("q-values.pdf", width=8, height=6)
par(mfrow = c(3,2))
plot(qobj$qvalues)
plot(qobj2$qvalues)
plot(qobj3$qvalues)
plot(qobj4$qvalues)
plot(qobj5$qvalues)
dev.off()

#seem fine

#now calculate the median and then backtransform to -log10
qvalues <- qobj$qvalues
qvalues2 <- qobj2$qvalues
qvalues3 <- qobj3$qvalues
qvalues4 <- qobj4$qvalues
qvalues5 <- qobj5$qvalues

combinedq <- data.frame(q1=qvalues, q2=qvalues2, q3=qvalues3, q4=qvalues4, q5=qvalues5)

head(qvalues)
head(qvalues2)

combinedq <- cbind(marker = 1:nrow(combinedq), combinedq)
head(combinedq)

#check for convergence
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
#0.4 for all

#visualise correlation
#correlation test
install.packages("ggpubr")
library("ggpubr")
str(combined)
ggscatter(data = combinedq, x = "q1", y = "q2", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "c2-model1", ylab = "c2-model2")
ggscatter(data = combinedq, x = "q1", y = "q3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "c2-model1", ylab = "c2-model3")



#what is going ooon? :D
#no convergence at all ! \o/
#so do a meadian of the values? ASK REID IF HE'S OKAY WITH THIS! 

#save output model runs

#The object can be summarized and visualized by:
# Open a text file for writing
sink("output_summary-modelruns.txt")

# Run the summary functions and print the results
sum1 <- summary(qobj)
print(sum1)

sum2 <- summary(qobj2)
print(sum2)

sum3 <- summary(qobj3)
print(sum3)

sum4 <- summary(qobj4)
print(sum4)

sum5 <- summary(qobj5)
print(sum5)

# Close the text file
sink()

pdf("sumplots-model-runs.pdf", width=8, height=6)
par(mfrow= c(3,2))
plot(qobj)
plot(qobj2)
plot(qobj3)
plot(qobj4)
plot(qobj5)
dev.off()

pdf("histo-model-runs.pdf", width=8, height=6)
par(mfrow= c(3,2))
hist(qobj)
hist(qobj2)
hist(qobj3)
hist(qobj4)
hist(qobj5)
dev.off()
##############

head(combinedq)
combinedq$median_q <- apply(combinedq[, c("q1", "q2", "q3", "q4", "q5")], 1, median)
range(combinedq$log_transformed)
combinedq$log_transformed <- -log10(combinedq$median_q)

pdf("qvalues-final.pdf", width=6, height=8)
par(mfrow=c(2,1))
plot(combinedq$median_q, xlab="marker", ylab="qvalue")
abline(h = 0.05, col = "red", lty = 2)
plot(combinedq$log_transformed, xlab="marker", ylab= "q(-log10 scale)")
abline(h = 1.30103, col = "red", lty = 2)
dev.off()

