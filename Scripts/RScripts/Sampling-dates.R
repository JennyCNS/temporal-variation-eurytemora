#plotting sampling dates used for sequencing

install.packages("dplyr")

library(dplyr)
library(readxl)
seqkit <- read_excel("~/Jenny-Eurytemora_affinis/QC-trimmed-rawdata/seqkit.xlsx")
View(seqkit)

df <- slice(seqkit, 1:(n() -1))
tail(df)


#plot sampling dates

install.packages("ggplot2")
library(ggplot2)
str(seqkit)
plot =ggplot(data=df, aes(x=date, y= value)) +
  geom_point() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(x="Date",
       y="Sampling",
       title = "Sampling Eurytemora affinis (2007-2022)") +
  facet_wrap(~year)
plot

ggsave("sampling_dates.png", width = 6, height = 6, units =c("cm"), dpi =3000)
