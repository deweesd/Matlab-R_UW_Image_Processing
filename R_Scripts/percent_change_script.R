
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(scales)

diam_df <- read.csv("./Percent_chnge_diam.csv")
colnames(diam_df) <- c("Rat_ID", "Injury", "5x-Rostral", "4x-Rostral", "3x-Rostral", "2x-Rostral", "1x-Rostral", "Center", "1x-Caudal", "2x-Caudal", "3x-Caudal", "4x-Caudal", "5x-Caudal")





## reshape wide format to long format for visualization ## 
updated_df <- diam_df %>% gather(Location, Percent_change, `5x-Rostral`:`5x-Caudal`)

# re-ordering factor levels #
updated_df$Location <- factor(updated_df$Location, levels = c("5x-Rostral", "4x-Rostral", "3x-Rostral", "2x-Rostral", "1x-Rostral", "Center", "1x-Caudal", "2x-Caudal", "3x-Caudal", "4x-Caudal", "5x-Caudal"))

#summary(updated_df)



# Percent Change barplot #
t <- ggplot(updated_df, aes(updated_df$Location, updated_df$Percent_change)) +
  geom_bar(aes(fill = Injury), position = "dodge", stat="identity")


?error.plot


p <- ggbarplot(updated_df, x = "Location", y = "Percent_change", add = "mean_se",
               error.plot = "errorbar",
               color = "Injury", fill = "Injury", palette = c("blue", "red"), 
               position = position_dodge(0.8))+
  theme(
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    legend.title = element_text(colour = "black", size = 14, face = "bold"),
    legend.text = element_text(colour = "black", size = 14, face = "plain"),
    axis.text.x = element_text(colour = "black", size = 14, face = "plain"),
    axis.text.y = element_text(colour = "black", size = 14, face = "plain")) +
  #stat_compare_means(aes(group = Injury ), label = "p.signif", label.y = 2, size = 10, hide.ns = TRUE)+
  labs( x="Location", y = "Percent Change from Baseline")

p

warnings()
