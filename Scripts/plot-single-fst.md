# re-runing the analyisis to vcf file
```bash
module load R/4.1.1
R
```
```R
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstsingle-snps-newvcf-queue.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))

#means of each colums (are the some higher then others)

mean_fst <- colMeans(fst[, 5:ncol(fst)], na.rm = TRUE)

#scatterplot number of snps per window

mean_fst 
> means
EA_2009_T1.EA_2009_T2 EA_2009_T1.EA_2009_T3 EA_2009_T1.EA_2009_T4
         0.0022117643          0.0021782320          0.0019972159
EA_2009_T1.EA_2011_T1 EA_2009_T1.EA_2011_T2 EA_2009_T1.EA_2015_T1
         0.0021032638          0.0020796221          0.0021146029
EA_2009_T1.EA_2015_T2 EA_2009_T1.EA_2015_T3 EA_2009_T1.EA_2015_T4
         0.0021387906          0.0023764857          0.0020848236
EA_2009_T1.EA_2022_T1 EA_2009_T1.EA_2022_T2 EA_2009_T1.EA_2022_T3
         0.0021692650          0.0022717752          0.0022698521
EA_2009_T1.EA_2022_T4 EA_2009_T2.EA_2009_T3 EA_2009_T2.EA_2009_T4
         0.0021890472          0.0011143184          0.0009553926
EA_2009_T2.EA_2011_T1 EA_2009_T2.EA_2011_T2 EA_2009_T2.EA_2015_T1
         0.0012208570          0.0010285965          0.0010164950
EA_2009_T2.EA_2015_T2 EA_2009_T2.EA_2015_T3 EA_2009_T2.EA_2015_T4
         0.0011170236          0.0016051377          0.0009986369
EA_2009_T2.EA_2022_T1 EA_2009_T2.EA_2022_T2 EA_2009_T2.EA_2022_T3
         0.0011358893          0.0012167867          0.0011948500
EA_2009_T2.EA_2022_T4 EA_2009_T3.EA_2009_T4 EA_2009_T3.EA_2011_T1
         0.0010961251          0.0007104736          0.0010649986
EA_2009_T3.EA_2011_T2 EA_2009_T3.EA_2015_T1 EA_2009_T3.EA_2015_T2
         0.0007819679          0.0007338291          0.0008701361
EA_2009_T3.EA_2015_T3 EA_2009_T3.EA_2015_T4 EA_2009_T3.EA_2022_T1
         0.0014888931          0.0007277233          0.0008983651
EA_2009_T3.EA_2022_T2 EA_2009_T3.EA_2022_T3 EA_2009_T3.EA_2022_T4
         0.0009906515          0.0009526519          0.0008267787
EA_2009_T4.EA_2011_T1 EA_2009_T4.EA_2011_T2 EA_2009_T4.EA_2015_T1
         0.0008297690          0.0005919081          0.0005607921
EA_2009_T4.EA_2015_T2 EA_2009_T4.EA_2015_T3 EA_2009_T4.EA_2015_T4
         0.0006862987          0.0012855092          0.0005493201
EA_2009_T4.EA_2022_T1 EA_2009_T4.EA_2022_T2 EA_2009_T4.EA_2022_T3
         0.0006775511          0.0007576983          0.0007187376
EA_2009_T4.EA_2022_T4 EA_2011_T1.EA_2011_T2 EA_2011_T1.EA_2015_T1
         0.0006200964          0.0008388475          0.0007752428
EA_2011_T1.EA_2015_T2 EA_2011_T1.EA_2015_T3 EA_2011_T1.EA_2015_T4
         0.0008615173          0.0013708276          0.0007455091
EA_2011_T1.EA_2022_T1 EA_2011_T1.EA_2022_T2 EA_2011_T1.EA_2022_T3
         0.0007816151          0.0008773498          0.0008340204
EA_2011_T1.EA_2022_T4 EA_2011_T2.EA_2015_T1 EA_2011_T2.EA_2015_T2
         0.0007739803          0.0005461445          0.0006634203
EA_2011_T2.EA_2015_T3 EA_2011_T2.EA_2015_T4 EA_2011_T2.EA_2022_T1
         0.0013358316          0.0005292419          0.0006461348
EA_2011_T2.EA_2022_T2 EA_2011_T2.EA_2022_T3 EA_2011_T2.EA_2022_T4
         0.0007299168          0.0006968359          0.0005954428
EA_2015_T1.EA_2015_T2 EA_2015_T1.EA_2015_T3 EA_2015_T1.EA_2015_T4
         0.0005338422          0.0012264019          0.0003891341
EA_2015_T1.EA_2022_T1 EA_2015_T1.EA_2022_T2 EA_2015_T1.EA_2022_T3
         0.0005491474          0.0006390626          0.0005818566
EA_2015_T1.EA_2022_T4 EA_2015_T2.EA_2015_T3 EA_2015_T2.EA_2015_T4
         0.0004911598          0.0012947726          0.0005121281
EA_2015_T2.EA_2022_T1 EA_2015_T2.EA_2022_T2 EA_2015_T2.EA_2022_T3
         0.0006568482          0.0007577106          0.0006999701
EA_2015_T2.EA_2022_T4 EA_2015_T3.EA_2015_T4 EA_2015_T3.EA_2022_T1
         0.0006185678          0.0012189289          0.0013398255
EA_2015_T3.EA_2022_T2 EA_2015_T3.EA_2022_T3 EA_2015_T3.EA_2022_T4
         0.0014295481          0.0014076999          0.0013323899
EA_2015_T4.EA_2022_T1 EA_2015_T4.EA_2022_T2 EA_2015_T4.EA_2022_T3
         0.0005112350          0.0006103877          0.0005452260
EA_2015_T4.EA_2022_T4 EA_2022_T1.EA_2022_T2 EA_2022_T1.EA_2022_T3
         0.0004602390          0.0006507826          0.0006013523
EA_2022_T1.EA_2022_T4 EA_2022_T2.EA_2022_T3 EA_2022_T2.EA_2022_T4
         0.0005300033          0.0006770812          0.0006130283
EA_2022_T3.EA_2022_T4
         0.0005440713


unique_values <- as.list(unique(fst$end))
unique_values
#[[1]] 1000
#[[2]]3000
#[[3]]2000
#[[4]] 4000
#[[5]] 5000
#[[6]] 6000
#[[7]] 7000
#[[8]] 8000
#[[9]] 9000
#[[10]] 10000
#[[11]] 11000
#[[12]] 12000
#[[13]] 13000



#boxplot to check variation between each pairwise comparison, this way we can see which ones are higher
fst2 <- fst[,-(1:4)]
fst2 <- as.data.frame(fst2)

# Reshape the data into long format
data_long <- tidyr::gather(fst2, key = "Column", value = "Value")
data_long$Column <- as.factor(data_long$Column)

# Create boxplots using ggplot2
d <- ggplot(data_long, aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplotfst-singlesnps.pdf",d, w=7.5, h=4.7)

 #check how the data looks like with no zeros

dl_filtered <- subset(data_long, rowSums(data_long != 0) > 0)

#nrow(dl_filtered)
#108974754


#boxplot of filtered data

ggplot(dl_filtered , aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#8160 rows removed with missing values


#subset for pairwise comparisons of interest

column_names <- colnames(fst2)
print(column_names)

#somethings funny, lets try individual pops
fst_subset<- fst2[, c("EA_2009_T1.EA_2009_T2")
# Remove zeros from the vector
non_zero_values <- fst_subset[fst_subset != 0]

ggplot(data.frame(Values = non_zero_values), aes(x = "Values", y = Values)) +
  #geom_boxplot(outlier.shape = NA) +  # Hide individual outlier points
  geom_point(aes(color = "red"), position = position_jitter(width = 0.2), show.legend = FALSE) +  # Add jittered individual points as red dots
  scale_color_identity(guide = "none") +  # Hide the color legend
  labs(title = "Fst EA_2009_T1 vs EA_2009_T2 ", x = NULL, y = "EA_2009_T1 vs EA_2009_T2") +
  theme_minimal()

####
#ok this works so I will do this to all the important combinations


fst_subset2<- fst2[, c("EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4", "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4","EA_2011_T1.EA_2011_T2",  "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", 
	"EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4","EA_2022_T1.EA_2022_T2",
	"EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]

#transform to a  long data

data_long <- tidyr::gather(fst_subset2, key = "Column", value = "Value")

#check if data is okay
table(data_long$Column)
#ok all names have the same count

data_long$Column <- as.factor(data_long$Column)

non_zero_dataset <- data_long[data_long$Value != 0, ]

# Print the resulting dataset without zero values
print(non_zero_dataset)

#boxplot

d <- ggplot(non_zero_dataset, aes(x = Column, y = Value)) +
  #geom_point(aes(color = Column), position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_point(aes(color = Column), position = position_jitter(width = 0.4), show.legend = FALSE, alpha = 0.5) +
  labs(title = "Fst single snps", x = NULL, y = "Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/fst-scatter-main-comparisons.pdf",d, w=12, h=6)

#create new variable year

fst_subset <- as.data.frame(non_zero_dataset)
data_clean <- na.omit(fst_subset)
data_clean$Column <- as.factor(data_clean$Column)

names(data_clean)[1] <- "Column"
year<- rep (NA, length(data_clean$Column))
year[grep("2009", data_clean$Column)] <- "2009"
year[grep("2011", data_clean$Column)] <- "2011"
year[grep("2015", data_clean$Column)] <- "2015"
year[grep("2022", data_clean$Column)] <- "2022"

data_clean$year <- year



#scatter plot


d <- ggplot(data_clean, aes(x = Column, y = Value)) +
  #geom_point(aes(color = Column), position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_point(aes(color=Column), position = position_jitter(width = 0.4), show.legend = FALSE, alpha = 0.5) +
  labs(title = "Fst single snps", x = NULL, y = "Fst") +
  theme_bw()+
  scale_x_discrete(labels=c("2009(T1xT2)", "2009(T1xT3)", "2009(T1xT4)", "2009(T2xT3)", "2009(T2xT4)", "2009(T3xT4)",
    "2011(T1xT2)",
    "2015(T1xT2)", "2015(T1xT3)", "2015(T1xT4)", "2015(T2xT3)", "2015(T2xT4)", "2015(T3xT4)",
    "2022(T1xT2)", "2022(T1xT3)", "2022(T1xT4)", "2022(T2xT3)", "2022(T2xT4)", "2022(T3xT4)")) + 
  scale_color_manual(values=c("#D3DDDC","#D3DDDC","#D3DDDC","#D3DDDC","#D3DDDC","#D3DDDC",'#6699CC',"#F2AD00","#F2AD00","#F2AD00","#F2AD00","#F2AD00","#F2AD00","#00A08A","#00A08A","#00A08A","#00A08A","#00A08A","#00A08A"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/scatter-fst-main-comparisons.pdf",d, w=10, h=8)

#violin plot 


d <- ggplot(data_clean, aes(x = Column, y = Value, fill= Column)) +
  geom_violin(show.legend = FALSE) +
  #geom_point(aes(color=Column), position = position_jitter(width = 0.4), show.legend = FALSE, alpha = 0.5) +
  labs(title = "Fst single snps", x = NULL, y = "Fst") +
  theme_bw()+
  scale_x_discrete(labels=c("2009(T1xT2)", "2009(T1xT3)", "2009(T1xT4)", "2009(T2xT3)", "2009(T2xT4)", "2009(T3xT4)",
    "2011(T1xT2)",
    "2015(T1xT2)", "2015(T1xT3)", "2015(T1xT4)", "2015(T2xT3)", "2015(T2xT4)", "2015(T3xT4)",
    "2022(T1xT2)", "2022(T1xT3)", "2022(T1xT4)", "2022(T2xT3)", "2022(T2xT4)", "2022(T3xT4)")) + 
  scale_fill_manual(values=c("#D3DDDC","#D3DDDC","#D3DDDC","#D3DDDC","#D3DDDC","#D3DDDC",'#6699CC',"#F2AD00","#F2AD00","#F2AD00","#F2AD00","#F2AD00","#F2AD00","#00A08A","#00A08A","#00A08A","#00A08A","#00A08A","#00A08A"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/violin-fst-main-comparisons.pdf",d, w=10, h=8)




#trial heatmap


mean_fst_subset <- colMeans(fst_subset, na.rm = TRUE)
mean_fst_subset

EA_2009_T1.EA_2009_T2 EA_2009_T1.EA_2009_T3 EA_2009_T1.EA_2009_T4
         6.522382e-04          9.119229e-04          8.192032e-04
EA_2009_T2.EA_2009_T3 EA_2009_T2.EA_2009_T4 EA_2009_T3.EA_2009_T4
         2.212941e-04          1.727231e-04          9.451288e-05
EA_2011_T1.EA_2011_T2 EA_2015_T1.EA_2015_T2 EA_2015_T1.EA_2015_T3
         3.314538e-04          7.216167e-05          5.353443e-04
EA_2015_T1.EA_2015_T4 EA_2015_T2.EA_2015_T3 EA_2015_T2.EA_2015_T4
         4.000097e-05          5.151115e-04          8.000057e-05
EA_2015_T3.EA_2015_T4 EA_2022_T1.EA_2022_T2 EA_2022_T1.EA_2022_T3
         6.021291e-04          1.041804e-04          1.297219e-04
EA_2022_T1.EA_2022_T4 EA_2022_T2.EA_2022_T3 EA_2022_T2.EA_2022_T4
         1.114017e-04          1.310401e-04          1.054126e-04
EA_2022_T3.EA_2022_T4
         8.650827e-05

heatmap_matrix <- mean_fst_subset %>% rownames_to_column(var = "Pairwise") %>% pivot_longer(cols = -Pairwise, names_to = "Variable", values_to = "Mean")

#fix data cause something went wrong

mean_fst_matrix <- matrix(
  c(
    6.522382e-04, 9.119229e-04, 8.192032e-04,
    2.212941e-04, 1.727231e-04, 9.451288e-05,
    3.314538e-04, 7.216167e-05, 5.353443e-04,
    4.000097e-05, 5.151115e-04, 8.000057e-05,
    6.021291e-04, 1.041804e-04, 1.297219e-04,
    1.114017e-04, 1.310401e-04, 1.054126e-04,
    8.650827e-05),
 nrow = 19,
 ncol=19,   # Number of pairwise comparisons (rows)
 byrow = TRUE,
 dimnames = list(
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4"),
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4")
  	)
  )

# Plot a heatmap using the heatmap() function

# Get the row and column names from dimnames
row_names <- rownames(mean_fst_matrix)
col_names <- colnames(mean_fst_matrix)

# Calculate the number of rows and columns in the matrix
n_rows <- nrow(mean_fst_matrix)
n_cols <- ncol(mean_fst_matrix)

#list diagonal values
 diagonal_values <- diag(mean_fst_matrix)

# Create a data frame with the diagonal values
diagonal_data <- data.frame(
  Pairwise = rownames(mean_fst_matrix),
  Value = diagonal_values
)

# Plot the heatmap of diagonal values
d <- ggplot(diagonal_data, aes(x = Pairwise, y = Pairwise, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

print(heatmap_plot)
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/heatmap-fst-window1000-maincomparisons.pdf",d, w=10, h=8)




#plot full fst
packages <- c("ComplexHeatmap", "grid", "Rcpp", "stats", "base", "plotrix", "wesanderson", "textshape")
install.packages(packages)
#install complexheatmap for older version
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
library("grid")
library("Rcpp")
library("stats")
library("base")
library("plotrix")
library("wesanderson")
library("textshape")
install.packages("pheatmap")
library(pheatmap)

fst_cleaned <- na.omit(fst)
missing_values <- any(is.na(fst))
if (missing_values) {
  cat("Missing values detected in the data.\n")
} else {
  cat("No missing values found in the data.\n")
}


# Remove row names and column names for heatmap
data_for_heatmap <- fst_cleaned[, 5:ncol(fst)]  # Exclude the first column (chrom)

fst_cleaned$countdown <- seq(nrow(fst_cleaned), 1)
fst_cleaned$full_name <- paste0(fst_cleaned$chrom, fst_cleaned$countdown)
rownames(data_for_heatmap) <- fst_cleaned$full_name


# Convert data to a matrix
mat <- as.matrix(data_for_heatmap)

# Create a heatmap using pheatmap
d <- pheatmap(mat,
         scale = "none",  # You can customize scaling if needed
         col = colorRampPalette(c("white", "blue"))(100),
         fontsize_row = 8,
         fontsize_col = 8,
      
         annotation_legend = TRUE, annotation_names_row = FALSE, annotation_names_col = TRUE,
         drop_levels = FALSE, show_rownames = F, show_colnames = T, 
         main = "Heatmap of Data")

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trial-heatmap-fst-window1000.pdf",d, w=10, h=8)

########################snps per bin####################
#####comparison between 1000 and 5000 bin sizes#########
########################################################


fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fst-all-pairwise.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))
fst_cleaned <- na.omit(fst)

sum_smaller_than_10 <- sum(fst_cleaned$snps < 2)

# Print the result
cat("Number of entries smaller than 10:", sum_smaller_than_10, "\n")
#9596
#1349 smaller than 2

#create a new colum for each bin
#fst_cleaned$bins <- paste(fst_cleaned$start, fst_cleaned$end, sep = "-")
fst_cleaned$bins <- paste(fst_cleaned$chrom, fst_cleaned$start, fst_cleaned$end,sep = "-")
fst_cleaned$bins<- factor(fst_cleaned$bins)

#number of bins
num_different_values <- length(unique(fst_cleaned$bins))
num_different_values
#42728 chroms
#2.5% of the bins only have 1 SNPs, perhaps using 5K windows makes more sense 

#plot scatter plot per bin
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trial-snps-per.bin-1000.pdf", height=6,width=10)
plot(fst_cleaned$bins, fst_cleaned$snps, main = "Scatter Plot", xlab = "Bins", ylab = "N SNPs", ylim = c(0, 800), pch = 16, col = "blue",  xaxt = "n")
dev.off()


###########################
#for 5000K window size #
###########################

fst <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstall-parwise-5000.csv", header = TRUE, na.strings = "")
#negative values to zeros

fst <- fst %>% mutate_all(function(x) ifelse(x < 0, 0, x))

#means of each colums (are the some higher then others)

mean_fst <- colMeans(fst[, 5:ncol(fst)], na.rm = TRUE)

#scatterplot number of snps per window

mean_fst 


#boxplot to check variation between each pairwise comparison, this way we can see which ones are higher

fst2 <- as.data.frame(fst2)
#remove missing data

#count NAs
sumnas <- sum(is.na(fst2))
#2411 less missing data than 1K windows

### remove mising data
data_clean<- na.omit(fst2)

# Reshape the data into long format
data_long <- tidyr::gather(data_clean, key = "Column", value = "Value")
data_long$Column <- as.factor(data_long$Column)

# Create boxplots using ggplot2
d <- ggplot(data_long, aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplotfst-5Kwindows.pdf",d, w=7.5, h=4.7)

#subset for pairwise comparisons of interest

column_names <- colnames(fst2)
print(column_names)
 

fst_subset<- data_clean[, c("EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4", "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4","EA_2011_T1.EA_2011_T2",  "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", 
	"EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4","EA_2022_T1.EA_2022_T2",
	"EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]


#create new variable year

fst_subset <- as.data.frame(fst_subset)
#data_clean <- fst_subset
data_long <- tidyr::gather(fst_subset, key = "Column", value = "Value") # gather() function from the tidyr package to convert your dataset to long format.
data_long$Column <- as.factor(data_long$Column)

names(data_long)[1] <- "Column"
year<- rep (NA, length(data_long$Column))
year[grep("2009", data_long$Column)] <- "2009"
year[grep("2011", data_long$Column)] <- "2011"
year[grep("2015", data_long$Column)] <- "2015"
year[grep("2022", data_long$Column)] <- "2022"

data_long$year <- year


#plot

d <- ggplot(data_long, aes(x=Column, y=Value)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplot-fst-main-comparisons-5Kwindows.pdf",d, w=12, h=4.7)

#scater plots
d <- ggplot(data_long, aes(x = Column, y = Value)) +
  geom_point() + #~ Column indicates that each facet corresponds to a different column.
  facet_grid(. ~ year, scales = "free") + # The scales = "free" argument allows each facet to have its own y-axis scale.
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pairwise comparison", y = "Fst")


 ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/boxplot-fst-scatter-main-comparisons-5kwindows.pdf",d, w=12, h=6)


# Create a scatter plot grid using facet_grid


#trial heatmap


mean_fst_subset <- colMeans(fst_subset, na.rm = TRUE)
mean_fst_subset

heatmap_matrix <- mean_fst_subset %>% rownames_to_column(var = "Pairwise") %>% pivot_longer(cols = -Pairwise, names_to = "Variable", values_to = "Mean")

#fix data cause something went wrong

mean_fst_matrix <- matrix(
  c(
     6.251438e-04, 8.726515e-04, 7.699255e-04,
     2.069383e-04, 1.547305e-04, 8.344243e-05,
     3.222143e-04, 6.337224e-05, 5.170126e-04,
     3.289386e-05, 4.931738e-04, 7.707632e-05, 
     5.848110e-04, 9.421378e-05, 1.233745e-04, 
     1.060758e-04, 1.230012e-04, 9.586955e-05, 8.112256e-05),
 nrow = 19,
 ncol=19,   # Number of pairwise comparisons (rows)
 byrow = TRUE,
 dimnames = list(
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4"),
    c(
      "EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4",
      "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4", "EA_2009_T3.EA_2009_T4",
      "EA_2011_T1.EA_2011_T2", "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3",
      "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4",
      "EA_2015_T3.EA_2015_T4", "EA_2022_T1.EA_2022_T2", "EA_2022_T1.EA_2022_T3",
      "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4",
      "EA_2022_T3.EA_2022_T4")
  	)
  )

# Plot a heatmap using the heatmap() function

# Get the row and column names from dimnames
row_names <- rownames(mean_fst_matrix)
col_names <- colnames(mean_fst_matrix)

# Calculate the number of rows and columns in the matrix
n_rows <- nrow(mean_fst_matrix)
n_cols <- ncol(mean_fst_matrix)

#list diagonal values
 diagonal_values <- diag(mean_fst_matrix)

# Create a data frame with the diagonal values
diagonal_data <- data.frame(
  Pairwise = rownames(mean_fst_matrix),
  Value = diagonal_values
)

# Plot the heatmap of diagonal values
d <- ggplot(diagonal_data, aes(x = Pairwise, y = Pairwise, fill = Value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = NULL, y = NULL)

print(heatmap_plot)
ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/heatmap-fst-window5000-maincomparisons.pdf",d, w=10, h=8)


#bin calculations

sum_smaller_than_10 <- sum(data_clean$snps < 10)
sum_smaller_than_10
# Print the result
cat("Number of entries smaller than 10:", sum_smaller_than_10, "\n")
#0
data_clean<- na.omit(fst)

#create a new colum for each bin
#fst_cleaned$bins <- paste(fst_cleaned$start, fst_cleaned$end, sep = "-")
data_clean$bins <- paste(data_clean$chrom, data_clean$start, data_clean$end,sep = "-")
data_clean$bins<- factor(data_clean$bins)

#number of bins
num_different_values <- length(unique(data_clean$bins))
num_different_values
#41161
#0 % of the bins has less than 10 snps per bin

#plot scatter plot per bin
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trial-snps-per.bin-5000.pdf", height=6,width=10)
plot(data_clean$bins, data_clean$snps, main = "Scatter Plot", xlab = "Bins", ylab = "N SNPs", pch = 16, col = "blue",  xaxt = "n")
dev.off()
```

# Manhattan plot of single snp by snp Fst calculation

``` bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
module load R/4.1.1
module load gcc/10.2.0
R
```
```R
library(pcadapt)
library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)
library(tidyr)



#input data must be a matrix with samples as row names and fst values as colums
freq <- read.csv("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstsingle-snps.csv", header = TRUE)

# Count rows with NaN values
sum(apply(freq, 1, function(row) any(is.nan(row))))
freq_clean <- na.omit(freq)
save(freq_clean, file = "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq_clean.RData")

#########################################################
######trial boxplot for relevant pairwise comparisons####
#########################################################

#Processed 41599 chromosomes with 1258282 (non-filtered) positions in 6872734 windows.
#Total filter summary (after applying all sample filters):
#Passed:               1258282

freq_clean<- na.omit(freq_clean)
freq_clean2  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))



fst_subset<- freq_clean2[, c("EA_2009_T1.EA_2009_T2", "EA_2009_T1.EA_2009_T3", "EA_2009_T1.EA_2009_T4", "EA_2009_T2.EA_2009_T3", "EA_2009_T2.EA_2009_T4","EA_2009_T3.EA_2009_T4","EA_2011_T1.EA_2011_T2",  "EA_2015_T1.EA_2015_T2", "EA_2015_T1.EA_2015_T3", 
  "EA_2015_T1.EA_2015_T4", "EA_2015_T2.EA_2015_T3", "EA_2015_T2.EA_2015_T4", "EA_2015_T3.EA_2015_T4","EA_2022_T1.EA_2022_T2",
  "EA_2022_T1.EA_2022_T3", "EA_2022_T1.EA_2022_T4", "EA_2022_T2.EA_2022_T3", "EA_2022_T2.EA_2022_T4", "EA_2022_T3.EA_2022_T4")]
fst2 <- fst_subset[,-(1:4)]
fst2  <- fst2  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
##############################################continuation of analyisis


df_long <- freq_clean %>%
  pivot_longer(cols = 5:95, names_to = "Chrom", values_to = "Value")
df_long$combined <- paste(df_long$chrom, df_long$start)



#ad a numeric column for chrom
data <- df_long %>%
  mutate(number = match(chrom, unique(chrom)))
save("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/formated_table_FSTSINGLE.RData")

#simple plot
x <- seq_along(df_long$chrom)

pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/singlesnpfst.pdf", height=6,width=10)
plot(data$number, data$Value, xlab = "Index", ylab = "Fst", main = "Single SNP Fst", xaxt = "n")
dev.off()

#doesnt work with FST data?

library(qqman)
load("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/longFST.RData")
#plot manhattan plot long data
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/manhattan-singlesnps.pdf", height=6,width=10)
manhattan(data, chr="number", bp="start", snp="combined", p="Value" )
dev.off()

###

load("/gxfs_work1/geomar/smomw573/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/formated_table_FSTSINGLE.RData")
d <- ggplot(data, aes(x=number, y=Value)) +
  geom_point(size=2, shape=23)+
    theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/trialscatterfst.pdf",d, w=7.5, h=4.7)
####
#calculate number of SNPs per chrom



# Define the file path for the new CSV file (same folder)
output_file <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/fstsingle-snps-cleaned.csv"

# Save the cleaned data frame as a new CSV file
write.csv(freq_clean, file = output_file, row.names = FALSE)

cat("Cleaned data frame has been saved as:", output_file, "\n")

#why dont we have an FST for all snps?
nrow(freq_clean)
# FST values 1098611 from 1258282 SNPs in original vcf file, why?


freq_clean <- freq_clean[, !(names(freq_clean) %in% c("start", "end", "snps"))]
freq_clean  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
save(freq_clean, file = "nonanozeros-allsnps.rds")
df <- freq_clean[-1,-1]

count_zero_columns <- sum(rowSums(df) == 0)
print(count_zero_columns)
# 555450

trans<- t(freq_clean)
mtrans <- as.matrix(trans)
mtrans<-matrix(mtrans, ncol= ncol(trans), dimnames =NULL)


#create pcadapt file
output_file <- "/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/pcadapt-input.csv"

# Save the cleaned data frame as a new CSV file
write.table(mtrans, file = output_file, row.names = FALSE, col.names = FALSE, quote=FALSE)

#heatmap FST snp per snp

library("ComplexHeatmap")
library("grid")
library("Rcpp")
library("stats")
library("base")
library("plotrix")
library("wesanderson")
library("textshape")
library(pheatmap)

load ("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq_clean.RData")


freq_clean  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
check <- freq_clean[c(1:100),c(1:95)]

# Remove row names and column names for heatmap
data_for_heatmap <- freq_clean[, 5:ncol(fst)]  # Exclude the first column (chrom)

freq_clean$countdown <- seq(nrow(freq_clean), 1)
freq_clean$full_name <- paste0(freq_clean$chrom, freq_clean$countdown)
rownames(data_for_heatmap) <- freq_clean$full_name

# Convert data to a matrix
mat <- as.matrix(data_for_heatmap)
subset_mat <- mat[,c(1,2,3,64,65,66,86,87,88,89,90,91)]

#freq_clean_sub <- data_for_heatmap[,c(1,2,3,64,65,66,86,87,88,89,90,91)]
#plot mean of single snps across all samples

mean_fst <- rowMeans(data_for_heatmap)
colnames(mean_fst) <- c("chrom", "mean_fst")
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/singlesnpfst_mean.pdf", height=6,width=10)
plot(mean_fst, xlab = "Index", ylab = "Fst", main = "Mean Single SNP Fst")
dev.off()

#plot histogram

#transform data columns into rows
load ("/gxfs_home/geomar/smomw573/work/seasonal_adaptation/analysis/grenedalf-freq/vcf-maf0.01-misfrac0.1-maxcov95/grenedalfst.log/freq_clean.RData")
freq_clean  <- freq_clean  %>% mutate_all(function(x) ifelse(x < 0, 0, x))
slot1 <- freq_clean[, 5:ncol(fst)]  # Exclude the first column (chrom)


#frombefore


new_row_names <- paste0(row.names(slot1), "_", colnames(slot1)[-1])
# Assign the new row names to the dataframe
row.names(data) <- new_row_names

# Reshape the data
reshaped_data <- slot1 %>%
  tibble::rownames_to_column(var = "row_name") %>%
  pivot_longer(-row_name, names_to = "Column", values_to = "Value")
#plot histogram
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/histogram_singlesnp_fst.pdf", height=6,width=10)
hist_data <- hist(reshaped_data$Value, xlab="breaks", ylab= "Frequency", breaks=30)
# Calculate midpoints of each bin
midpoints <- (hist_data$breaks[-1] + hist_data$breaks[-length(hist_data$breaks)]) / 2
# Get counts in each bin
counts <- hist_data$counts
text(hist_data$mids, hist_data$counts, labels = hist_data$counts, pos = 3)
dev.off()

sum(reshaped_data$Value > 0.95)


#freq table
slot1 <- freq_clean[, c(1, (5:ncol(freq_clean)))]
transformed_data <- slot1 %>%
  pivot_longer(cols = -chrom, names_to = "Column", values_to = "Value") %>%
  select(chrom, Column, Value)

freq <- as.data.frame(table(transformed_data$chrom))
sum(transformed_data$Value == 0)
#93805542
sum(transformed_data$Value > 0)
#2872226


#plot SNPs per chrom
pdf(file="/gxfs_home/geomar/smomw573/work/seasonal_adaptation/figures/nsnpsperchrom.pdf", height=6,width=10)
plot(freq$Freq, pch = 16, col = "black", xlab = "Index", ylab = "n SNPs", main = "SNPs per chrom")
dev.off()
mean(freq$Freq)
#2361.277
```

# now I want to estimate the genome wide Fst and plot it
# for this, I am using the poolfstat R package

```bash
srun --pty --x11 --nodes=1 --cpus-per-task=4 --mem=32000 --time=06:00:00 /bin/bash
module load R/4.1.1
module load gcc/10.2.0
module load icu4c/67.1

R
```

enter R environment 


```R
library(poolfstat)
#poolfstat for vcf

library(poolfstat)

#import vcf file
ea.readcount30X <- vcf2pooldata(vcf.file="depth-corrected.recode.vcf.gz",poolsizes=rep(50,14))
# summary of the resulting pooldata object
ea.readcount30X

#ha?
selected.snps.idx <- as.numeric(sub("rs","",rownames(ea.readcount30X@snp.info)))
head(selected.snps.idx)

#estimate genome wide Fst across all the popuatlions
ea.readcount30X.fst<-computeFST(ea.readcount30X)
ea.readcount30X.fst$FST
#genome wide Fst is -0.01960076, so 0

# Block-Jackknife estimation of FST standard-error and confidence intervals:
ea.readcount30X.fst<-computeFST(ea.readcount30X,nsnp.per.bjack.block = 1000, verbose=FALSE)
ea.readcount30X.fst$FST
#same value  -0.01960076

ea.readcount30X.fst$mean.fst #block-jacknife estimate of s.e.
#-0.02009452
ea.readcount30X.fst$se.fst #s.e. of the genome-wide Fst estimate
#2.90547e-05
ea.readcount30X.fst$mean.fst+c(-1.96,1.96)*ea.readcount30X.fst$se.fst
# -0.02015146 -0.02003757


#Computing multi-locus FST to scan the genome over sliding-windows of SNPs

ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=50)
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=100)
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=10)
#we have 42K scaffolds so loosing a lot of data with sliding window size of 50-100 chrom?


#6722 chromosomes scanned (with more than 50 SNPs)
#Average (min-max) Window Sizes 0.3 ( 0.1 - 4.9 ) kb


#931 chromosomes scanned (with more than 100 SNPs)
#Average (min-max) Window Sizes 0.8 ( 0.2 - 7.4 ) kb


#31251 chromosomes scanned (with more than 10 SNPs)
#Average (min-max) Window Sizes 0 ( 0 - 6.9 ) kb

#I will just play around with 100 as I dont really understand
ea.readcount30X.fst<-computeFST(ea.readcount30X,sliding.window.size=100)

plot(ea.readcount30X.fst$sliding.windows.fst$CumulatedPosition/1e6,
     ea.readcount30X.fst$sliding.windows.fst$MultiLocusFst,
     xlab="Cumulated Position (in Mb)",ylab="Muli-locus Fst")
     #col=as.numeric(ea.readcount30X.fst$sliding.windows.fst$Chr),pch=16) Doesnt work as we dont have chromossome numbers
abline(h=ea.readcount30X.fst$FST,lty=2)

head(ea.readcount30X.fst$sliding.windows.fst$CumulatedPosition/1e6)
head(ea.readcount30X.fst$sliding.windows.fst$MultiLocusFst)
head(ea.readcount30X.fst$sliding.windows.fst$Chr)

#Manhattan plot of the multi-locus FST computed over sliding-windows of 50 SNPs on the PoolSeq example data. The dashed line indicates the estimated overall genome-wide FST . The 20 simulated
#chromosomes are represented by alternate colors


#pairwise FST

ea.pairwisefst<-compute.pairwiseFST(ea.readcount30X,verbose=FALSE)
#heatmap
#Heatmap representing the pairwise-population FST matrix of the 14 populations of the 30XPool-Seq example data set
heatmap(ea.pairwisefst)
#it moves pops which are more similar to each other 

#Block-Jackknife estimation of FST standard-error and visualisation of confidence intervals
ea.pairwisefst@PairwiseFSTmatrix
plot(ea.pairwisefst)
```