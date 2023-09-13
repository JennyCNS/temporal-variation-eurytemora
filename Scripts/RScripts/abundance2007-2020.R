#plot abundances
install.packages("openxlsx")
library(openxlsx)
library(ggplot2)
library(dplyr)

af <- read.xlsx("Zooplankton_All.xlsx", sheet="data")

ea <- af[grep("Eurytemora sp.", af[["species.2"]]), ]
levels(ea$station)
str(ea)
ea$station <- as.factor(ea$station)
eanok <- ea[grep("NOK", ea[["station"]]), ]
str(eanok)
ea.final <- eanok[,c(2,4,6,7,8,11)]



ea.final <- ea.final[grep("", ea.final[["year"]]), ]

#plot
ea.final$month <- factor(ea.final$month)
ea.final$year <- factor(ea.final$year)


#transform different copepodites sex stages into 1
ea.final$sex.stage <- gsub("C1-C3|C4-C5", "Cop", ea.final$sex.stage)
ea.final$sex.stage <- gsub("Cop", "Copepodite", ea.final$sex.stage)

custom_colors <- c( "gray", "black", "#6683B5")
custom_colors2 <- c("#6683B5")

plot <- ggplot(ea.final, aes(x = month, y = abundance, fill = sex.stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ year, scales = "free_x") +
  labs(title = "Abundance of Different Life Stages by Month and Year",
       x = "Month",
       y = "Abundance",
       fill = "Sex Stage") +
  scale_fill_manual(values = custom_colors) +  # Assign custom colors
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.5, size = 10),  # Adjust x-axis text size
    axis.text.y = element_text(size = 10),  # Adjust y-axis text size
    axis.title.x = element_text(size = 11),  # Adjust x-axis label size
    axis.title.y = element_text(size = 11, margin = margin(t = 0, r = 10, b = 0, l = 0)),  # Adjust y-axis label size
    strip.background = element_blank(),
    axis.line = element_line(color = "black"),  # Customize x and y axis lines
    panel.background = element_rect(fill = "white"),  # Customize panel background
    plot.background = element_rect(fill = "white"),  # Customize plot background
    legend.background = element_rect(fill = "white", color = NA),  # Remove legend border
    legend.text = element_text(size = 11),  # Set legend text size
    plot.title = element_text(size = 14, hjust = 0.5)  # Adjust title size and center it
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

plot
ggsave("abundance-througout-years.pdf", plot, width = 10, height = 6, device = "pdf")

#not considering lifestages
head(ea.final)


# Assuming your original dataframe is named 'original_df'
# Group by year, month, and sex.stage, then summarize the abundance
summarized_df <- ea.final %>%
  group_by(year, month) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE))

plot.2 <- ggplot(summarized_df, aes(x = month, y = total_abundance)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#6683B5") +
  facet_wrap(~ year, scales = "free_x") +
  labs(title = expression(paste("Abundance of ", italic ("Eurytemora affinis"))),
       x = "Month",
       y = "Abundance") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    strip.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
plot.2

ggsave("abundance-througout-years-merged-stages.pdf", plot.2, width = 10, height = 6, device = "pdf")
