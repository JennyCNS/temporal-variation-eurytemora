#figures reid

library(dplyr)
library(tidyr)
library(ggplot2)

df <- read.table("\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/positive-sig-seasonal-snps-list.txt", header = TRUE)
glm <- read.table("full-table-glm-181023-2009t2-time-p.txt", header=TRUE)
head(df)
head(glm)

colnames(df) <- c("chrompos", "2009.start", "2011.start",
                  "2015.start", "2022.start",
                  "2009.end", "2011.end",
                  "2015.end", "2022.end")


af<- gather(df, key = "time", value = "af", -chrom)

head(af)
levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
af$time <- factor(af$time, levels = levels)

ggplot(data = af, aes(x = time, y = af, color = time)) +
  geom_line(aes(group = chrom), color = "lightgray", linewidth = 0.5) 
  labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()
  
#wrong plot
  
freq <- read.table(file="\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/freq_new_vcf.txt", header= T)
#chrom <- read.table(file="\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/chroms.txt", header= T)
sig <- read.table(file="\\\\helmholtz.geomar.de/Users$/jnascimento/Daten/Jenny-Eurytemora_affinis/AF-GLM/cov-sig-snps-glm-181023-2009t2-time-p.txt", header= F) 
  
sig$chrompos <- paste(sig$V1, sig$V2, sep="_")
af <- freq[, c("EA_2009_T2","EA_2009_T4" ,"EA_2011_T1", "EA_2011_T2","EA_2015_T1", "EA_2015_T4", "EA_2022_T1", "EA_2022_T4")]
head(af)
nrow(freq)
nrow(af)
head(sig)
chrom <- as.data.frame(glm$chrompos)
head(chrom)

#data with all AFs
af2 <- cbind(chrom, af)
head(af2)

colnames(af2) <- c("chrompos", "2009t2", "2009t4",
                  "2011t1", "2011t2",
                  "2015t1", "2015t4",
                  "2022t1", "2022t4")


#so now I need to pull out the positive ones from this

merged <- merge(df, af2, by = "chrompos")
head(merged)
merged <- merged[, names(af2)]

colnames(merged) <- c("chrompos", "2009.start", "2011.start",
                  "2015.start", "2022.start",
                  "2009.end", "2011.end",
                  "2015.end", "2022.end")
nrow(merged)
#now I want to make them all start at 0 (correct for 2009.start)
columns_to_modify <- names(merged)[-1]

# Subtract 2009t2 from all other columns
normalized_af <- merged %>%
  mutate(across(all_of(columns_to_modify), ~ . - `2009.start`))
head(normalized_af)
head(norm2)

norm2<- gather(normalized_af, key = "time", value = "af", -chrompos)

levels <- c("2009.start", "2009.end", "2011.start", "2011.end", "2015.start", "2015.end", "2022.start", "2022.end")
# Reorder the levels of the factor
norm2$time <- factor(norm2$time, levels = levels)
head(norm2)

ggplot(data = norm2, aes(x = time, y = af)) +
  geom_line(aes(group = chrompos), color = "lightgray", linewidth = 0.5) 
labs(x = "Year / Time Point", y = "AF", color = "Year") +
  ggtitle("") +
  theme_bw()

#############################



##restart




############################
#positives



head(af2)
nrow(af2)
head(af2)

colnames(af2) <- c("chrompos", "2009.start", "2009.end",
                   "2011.start", "2011.end",
                   "2015.start", "2015.end",
                   "2022.start", "2022.end")


norm <- af2 %>%
  mutate(
    '2009.end' = `2009.end` - `2009.start`,
    '2011.end' = `2011.end` - `2011.start`,
    '2015.end' = `2015.end` - `2015.start`,
    '2022.end' = `2022.end` - `2022.start`
  )

norm <- norm %>%
  mutate(
    '2009.start' = 0,
    '2011.start' = 0,
    '2015.start' = 0,
    '2022.start' = 0
  )

head(norm)
head(sig)
normaf2 <- norm %>%
  inner_join(sig, by = c("chrompos" = "chrompos")) %>%
  select(all_of(names(norm))) #do not include columns in sig

head(normaf2)

#positive snps
result <- normaf2 %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . > 0))
head(result)
result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2)

merged <- merge(result2, af2, by = "chrompos", suffixes = c(".result", ".af"))
nrow(merged)
head(merged)

positivesnps <- merged %>%
  select(chrompos, ends_with(".af"))

head(positivesnps)
#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

positive <- positivesnps %>%
  mutate(across(all_of(list), ~ . - `2009.start.af`))
head(positive)
head(final)

finalpost <- gather(positive, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
head
finalpost$year <- factor(finalpost$year, levels = levels)

# Reorder the levels of the factor


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
  scale_x_discrete(labels = factorsnames)+
  stat_summary(fun=median, color="dodgerblue3", geom = "line", group=1, size=1.5)

ggsave("af-normalizedperyear-pos.pdf",d, w=7, h=5)

#335 snps



#negative snps
head(normaf2)
result <- normaf2 %>%
  mutate(across(starts_with("2009.") | starts_with("2011.") | starts_with("2015.") | starts_with("2022."), ~ . < 0))
head(result)
nrow(result)
result2 <- result %>%
  filter_at(vars(ends_with(".end")), all_vars(. == TRUE))
nrow(result2)
head(af2)

merged <- merge(result2, af2, by.x = 1, by.y = 1,  suffixes = c(".result", ".af"))
nrow(merged)
head(merged)

positivesnps <- merged %>%
  select(chrompos, ends_with(".af"))

head(positivesnps)
#now I will correct them all to 2009
list <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")

positive <- positivesnps %>%
  mutate(across(all_of(list), ~ . - `2009.start.af`))
head(positive)

finalpost <- gather(positive, key = "year", value = "af", -chrompos)
levels <- c("2009.start.af", "2009.end.af", "2011.start.af", "2011.end.af", "2015.start.af", "2015.end.af", "2022.start.af", "2022.end.af")
head(finalpost)
finalpost$year <- factor(finalpost$year, levels = levels)

# Reorder the levels of the factor


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
  scale_x_discrete(labels = factorsnames)+
  stat_summary(fun=median, color="dodgerblue3", geom = "line", group=1, size=1.5)

ggsave("af-normalizedperyear-neg.pdf",d, w=7, h=5)

# 393snps
