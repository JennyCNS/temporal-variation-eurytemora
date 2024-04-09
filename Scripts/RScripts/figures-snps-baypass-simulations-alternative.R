#plotting figures to reid
#poster SAB
#simulation output (500 simulations baypass)
#29.02.2024

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Jenny-Eurytemora_affinis/excel-files")


df<- read.table("final-significance-simulations.txt", header=TRUE)
freq <- read.table(file="\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/freq_new_vcf.txt", header= T)
str(df)
str(freq)
head(df)
sum(df$pvalues_false <= 0.004)
#
af <- freq[, c("EA_2009_T2","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1", "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
af2 <- cbind(df$chrom, af)
af3 <- 1-af
afdec <- cbind(df$chrom, af3)

names(af2)[1] <- c("chrom")
names(afdec)[1] <- c("chrom")
str(af2)
#sigp false
sigpfalse <- df[df$pvalues_false < 0.05, ]

#sigq false
sigqfalse <- df[df$qvalues_false < 0.05, ]

#AF of pvalue significant snps
str(af2)
sigpvals <- merge(sigpfalse, af2, by = "chrom")
sigqvals <- merge(sigqfalse, af2, by = "chrom")

#for decreasing snps
sigpvalsdec <- merge(sigpfalse, afdec, by = "chrom")
sigqvalsdec <- merge(sigqfalse, afdec, by = "chrom")
write.table(sigpvals,"significant_pvalues_simulations.txt", sep="\t", quote=FALSE)
write.table(sigqvals,"significant_qvalues_simulations.txt", sep="\t", quote=FALSE)

#start with pvalues
sigpvals <- sigpvals[, names(af2)]
sigqvals <- sigqvals[, names(af2)]
sigpvalsdec <- sigpvalsdec[, names(afdec)]
sigqvalsdec <- sigqvalsdec[, names(afdec)]

colnames(sigpvals) <- c("chrompos", "2009.start", "2011.start",
                        "2015.start", "2022.start",
                        "2009.end", "2011.end",
                        "2015.end", "2022.end")
colnames(sigqvals) <- c("chrompos", "2009.start", "2011.start",
                        "2015.start", "2022.start",
                        "2009.end", "2011.end",
                        "2015.end", "2022.end")
colnames(sigpvalsdec) <- c("chrompos", "2009.start", "2011.start",
                           "2015.start", "2022.start",
                           "2009.end", "2011.end",
                           "2015.end", "2022.end")
colnames(sigqvalsdec) <- c("chrompos", "2009.start", "2011.start",
                           "2015.start", "2022.start",
                           "2009.end", "2011.end",
                           "2015.end", "2022.end")

#okay, fisrt visualise what each snp is doing normally
#say 2009.start is zero

#add row with numbers so I can ID SNPs


#first check how much they changed
norm <- sigqvals %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

#set start of the year to zero
norm <- norm %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )

#first check how much they changed
normdec <- sigqvalsdec %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

#set start of the year to zero
normdec <- normdec %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )
#select positive snps

result <- norm %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(result)

result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2)
#3278 SNPs
merged <- merge(result2, sigpvals, by = "chrompos", suffixes = c(".result", ".af"))
nrow(merged)
head(merged)
positivesnps <- merged %>%
  select(chrompos, ends_with(".af"))
head(positivesnps)

#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

#positive <- positivesnps %>%
#  mutate(across(all_of(list), ~ . - `2009.start.af`))
#head(positive)
#head(final)
positive <- positivesnps
finalpost <- gather(positive, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
finalpost$year <- factor(finalpost$year, levels = levels)

factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

d <- ggplot(data = finalpost, aes(x = year, y = af)) +
  geom_line(aes(group = chrompos), color = "black", size = 0.8, alpha = 0.1) +
  labs(x = "Time", y = "AF", color = "Year") +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14),   # Adjust x-axis label size
    axis.title.y = element_text(size = 14),   # Adjust y-axis label size
    plot.title = element_text(size = 16)     # Adjust plot title size
  ) +
  scale_x_discrete(labels = factorsnames)
#stat_summary(fun=median, color="dodgerblue3", geom = "line", group=1, size=1.5)
d
ggsave("af-normalizedperyear-pos.pdf",d, w=7, h=5)


#now negative ones

#select positive snps
#inverse AF with af3
resultdec <- normdec %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(resultdec)
result2dec <- resultdec %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2dec)
#3207 SNPs
#flip to positive

mergeddec <- merge(result2dec, sigpvalsdec, by = "chrompos", suffixes = c(".result", ".af"))
nrow(mergeddec)
head(mergeddec)
positivesnpsdec <- mergeddec %>%
  select(chrompos, ends_with(".af"))
head(positivesnpsdec)

#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

positivedec <- positivesnpsdec %>%
  mutate(across(all_of(list), ~ . - `2009.start.af`))
head(positivedec)


finalpostdec <- gather(positivedec, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
finalpostdec$year <- factor(finalpostdec$year, levels = levels)

factorsnames <- c("2009.start", "2009.end",
                  "2011.start", "2011.end",
                  "2015.start", "2015.end",
                  "2022.start", "2022.end"
)

df_combined <- rbind(finalpostdec, finalpost)

d <- ggplot(data = df_combined, aes(x = year, y = af)) +
  geom_line(aes(group = chrompos), color = "black", size = 0.8, alpha = 0.1) +
  labs(x = "Time", y = "AF", color = "Year") +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 14),   # Adjust x-axis label size
    axis.title.y = element_text(size = 14),   # Adjust y-axis label size
    plot.title = element_text(size = 16)     # Adjust plot title size
  ) +
  scale_x_discrete(labels = factorsnames) +
  stat_summary(fun="mean", color="dodgerblue3", geom = "line", group=1, size=1.5)
d
ggsave("af-normalizedperyear-all-pvalues.pdf",d, w=7, h=5)

?stat_summary

d <- ggplot(data = df_combined, aes(x = year, y = af)) +
  geom_line(aes(group = chrompos), color = "gray48", size = 1.2, alpha = 0.1) +
  labs(x = "Time", y = "AF", color = "Year") +
  theme_classic()+
  theme(
    axis.text.x=element_text(size=22),
    axis.text.y=element_text(size=22),
    axis.title.x = element_text(size = 26),   # Adjust x-axis label size
    axis.title.y = element_text(size = 26),   # Adjust y-axis label size
    plot.title = element_text(size = 28)     # Adjust plot title size
  ) +
  scale_x_discrete(labels = factorsnames) +
  stat_summary(fun="mean", color="dodgerblue3", geom = "line", group=1, size=2.5)
d
ggsave("af-normalizedperyear-all-qvalues-nostandard.pdf",d, w=14, h=7)


#plot distribution across genome pvalues

df_chrom <- as.data.frame(sigpvals[,1])
head(df_chrom)
hist(df_chrom)

# Extract chromosome information using regular expressions
chromosomes <- sub("^([^_]*_[^_]*)_.*", "\\1", df_chrom[, 1])
chrom <- as.data.frame(chromosomes)
chrom$count <- 1

chrom_count <- table(chrom$chromosomes)
chrom_count <- sort(chrom_count, decreasing = TRUE)
chrom_count_top10 <- chrom_count[1:10]
# Convert to data frame
chrom_counts_df <- data.frame(chromosomes = names(chrom_count_top10),
                              count = as.numeric(chrom_count_top10))

# Create ggplot
d <- ggplot(chrom_counts_df, aes(x = reorder(chromosomes, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue3", width = 0.7) +
  labs(title = "", x = "Chromosome", y = "Count") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("chromosome_counts_top10.pdf", d, w= 8, h= 6)



##### okay lets plot the manhattan plot
