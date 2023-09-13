

frequency_ind_samples <- read.csv("~/Jenny-Eurytemora_affinis/excel-files/frequency_ind_samples.csv", header = T)

View(frequency_ind_samples)

df <- t(frequency_ind_samples)
str(df)


##work wth subset
df2 <- df[, 1:5]
# Get the second row as the new header
new_header <- unlist(df2[1, ])

# Remove the first row (existing header) from the data frame
df2 <- df2[-1, ]

# Assign the new header to the data frame
colnames(df2) <- new_header

df2 <- df2[-c(1),]

df3<- na.omit(df2)


# Print the first column
print(df2[, 1])


library(ggplot2)
library(tidyr)

###########################
# PCA

##### snps


af <- read.table("~/tonsa_genomics/analysis/filtered_maf_freqs.txt", header=TRUE)

pops <- c("EA_2009_T1.FREQ", "EA_2009_T2.FREQ", "EA_2009_T3.FREQ",
          "EA_2009_T4.FREQ", "EA_2011_T1.FREQ", "EA_2011_T2.FREQ", 
          "EA_2015_T1.FREQ", "EA_2015_T2.FREQ", "EA_2015_T3.FREQ", 
          "EA_2015_T4.FREQ", "EA_2022_T1.FREQ", "EA_2022_T2.FREQ", 
          "EA_2022_T3.FREQ", "EA_2022_T4.FREQ") 

# remove F3 pop

af2 <- af[,grep("F03", colnames(af), invert=T)]
af2 <- af2[,grep("HH_F00", colnames(af2), invert=T)]
pops2 <- pops[grep("F03", pops, invert=T)]
pops2 <- pops2[grep("HH_F00", pops2, invert=T)]
pops <- pops2

varout <- apply(df3[,2:ncol(df3)], 1, var)

freqs <- t(df3[varout != 0,2:ncol(df3)])
colnames(freqs) <-  af$SNP[varout != 0]

nrow(freqs)

##
## plot pca
##

pcaResult <- prcomp(freqs, scale=T)

percentVar <- round(pcaResult$sdev^2/sum(pcaResult$sdev^2)*100, digits=2)

data <- data.frame(id=pops, Line=substr(pops, 1,2),
                   gen=substr(pops, 4,6),
                   PC1 = pcaResult$x[,1],  PC2= pcaResult$x[,2])

data$Line <- c(rep("Founding population", 4),
               rep("Ambient", 4),
               rep("Acidification", 4),
               rep("Warming", 4),
               rep("OWA", 4))

data$Line <- factor(data$Line, levels = c("Founding population","Ambient", "Acidification", "Warming", "OWA"))

data$PC2 <- data$PC2*-1

d <- ggplot(data, aes(PC1, PC2, fill=Line, shape=Line)) +
  geom_point(size=4.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw() +
  # ylim(-30, 23) + xlim(-50, 65)+
  scale_shape_manual(values=c( 21,21,22, 23, 24))+
  scale_color_manual(values=c('black')) +
  #scale_fill_manual(values=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3"),
  scale_fill_manual(values=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A", "#CC3333"),
                    labels = c("Founding population","Ambient", "Acidification",
                               "Warming", "OWA"))+
  #theme(legend.position = c(0.83,0.85),
  #    legend.background = element_blank(),
  #legend.box.background = element_rect(colour = "black"),
  theme(legend.title = element_blank())+
  # theme(legend.text=element_text(size=8))+
  #theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  #        ggtitle("F1")+
  guides(fill=guide_legend(override.aes=list(
    shape=c(21,21,22, 23, 24),
    #fill=c('steelblue1','steelblue','grey45', "darkorchid2", "firebrick3" )),order = 2),
    fill=c("#D3DDDC",'#6699CC',"#F2AD00","#00A08A", "#CC3333")),order = 2),
    shape= FALSE)

ggsave("~/tonsa_genomics/figures/pca_afs_noF3.pdf",d, w=5.5, h=3.7)