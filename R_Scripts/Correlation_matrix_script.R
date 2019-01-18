getwd()
setwd('/Users/Deweesd/Desktop/ANOVA_MODEL/')
getwd()

library(ggplot2)
library(ggpubr)
library(dplyr)
library(corrplot)
library("Hmisc")
library(gridExtra)
library(grid)
data_frame <- data_frame[-c(13:16),] # remove last 4 rows from df

## Read in files ##
data_frame <- read.csv("./Correlation_DF_Matrix.csv")

colnames(data_frame) <- c("Injury", "BL arrival time delay (sec)", "Acute arrival time delay (sec)", "Chronic arrival time delay (sec)", "BL area deficit size (mm2)", "Acute Area Deficit (mm2)", "Chronic area deficit size (mm2)", "BL BBB", "Week 1 BBB", "Week 8 BBB", "BL ladder walk % error",  "Week 4 ladder walk % error", "Week 8 ladder walk % error", "WK 8 Tortuoisity Distance (mm)" )


data_frame$Injury <- as.numeric(data_frame$Injury)
data_frame$`BL area deficit size (mm2)` <- NULL
data_frame$`BL BBB` <- NULL
data_frame$`Week 4 ladder walk % error` <- NULL


# for specific feature injury paper #
#data_frame$`BL arrival time delay (sec)` <- NULL
#data_frame$`Acute arrival time delay (sec)` <- NULL
#data_frame$`Chronic arrival time delay (sec)` <- NULL
#data_frame$`WK 8 Tortuoisity Distance (mm)` <- NULL
#data_frame$`Chronic area deficit size (mm2)` <- NULL

#Check levels for each index feature
str(data_frame)

#Subset df by injury type
data_frame_mod <- subset(data_frame, Injury =="2")
data_frame_severe <- subset(data_frame, Injury =="3")


par(mfrow=c(2,2))

#cor(data_frame_mod$`Acute Area Deficit (mm2)`, data_frame_mod$`Week 1 BBB`)
plot(data_frame_mod$`Week 1 BBB`, data_frame_mod$`Acute Area Deficit (mm2)`,  main = "Moderate: corr plot showing acute AD vs WK 1 BBB score",
     xlab="Week 1 BBB score", ylab="Acute perfusion deficit (mm2)", pch = 19, col = "blue")
text(10.2, 0.91, label = "r = -0.93", cex = 1.8)
abline(lm(data_frame_mod$`Acute Area Deficit (mm2)`~data_frame_mod$`Week 1 BBB`), col = "black")

#cor(data_frame_mod$`Acute Area Deficit (mm2)`, data_frame_mod$`Week 8 BBB`)
plot(data_frame_mod$`Week 8 BBB`, data_frame_mod$`Acute Area Deficit (mm2)`,  main = "Moderate: corr plot showing acute AD vs WK 8 BBB score",
     xlab="Week 8 BBB score", ylab="Acute perfusion deficit (mm2)", pch = 19, col = "blue")
text(17.8, 0.91, label = "r = -0.52", cex = 1.8)
abline(lm(data_frame_mod$`Acute Area Deficit (mm2)`~data_frame_mod$`Week 8 BBB`), col = "black")

#cor(data_frame_severe$`Acute Area Deficit (mm2)`, data_frame_severe$`Week 1 BBB`)
plot( data_frame_severe$`Week 1 BBB`, data_frame_severe$`Acute Area Deficit (mm2)`, main = "Severe: corr plot showing acute AD vs WK 1 BBB score",
      xlab="Week 1 BBB score", ylab="Acute perfusion deficit(mm2)", pch = 19, col = "red")
text(5, 2, label = "r = 0.35", cex = 1.8)
abline(lm(data_frame_severe$`Acute Area Deficit (mm2)`~data_frame_severe$`Week 1 BBB`), col = "black")

#R_value <- cor(data_frame_severe$`Acute Area Deficit (mm2)`, data_frame_severe$`Week 8 BBB`)
plot( data_frame_severe$`Week 8 BBB`, data_frame_severe$`Acute Area Deficit (mm2)`, main = "Severe: corr plot showing acute AD vs WK 8 BBB score",
      xlab="Week 8 BBB score", ylab="Acute perfusion deficit(mm2)", pch = 19, col = "red")
text(15, 2, label = "r = -0.47", cex = 1.8)
abline(lm(data_frame_severe$`Acute Area Deficit (mm2)`~data_frame_severe$`Week 8 BBB`), col = "black")
  
  
  
#annotate("text", x=15, y=2.0, label = "r == 0.35", parse=T, size=8.5) #adds R^2 value to graph
#rcorr(data_frame_severe$`Acute Area Deficit (mm2)`, data_frame_severe$`Week 1 BBB`, type = "pearson")

















data_frame$Injury <- NULL
data_frame_mod$Injury <- NULL
data_frame_severe$Injury <- NULL


## Setting up raw matrix for each df 
res_complete <- cor(data_frame)
round(res_complete, 2)
res_complete <- rcorr(as.matrix(data_frame))


res_mod <- cor(data_frame_mod)
round(res_mod, 2)
res_mod <- rcorr(as.matrix(data_frame_mod))


res_sev <- cor(data_frame_severe)
round(res_sev, 2)
res_sev <- rcorr(as.matrix(data_frame_severe))

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res_complete<-rcorr(as.matrix(data_frame_mod[,1:9]))
flattenCorrMatrix(res_complete$r, res_complete$P)

res_mod<-rcorr(as.matrix(data_frame_mod[,1:9]))
flattenCorrMatrix(res_mod$r, res_mod$P)


res_sev<-rcorr(as.matrix(data_frame_severe[,1:9]))
flattenCorrMatrix(res_sev$r, res_sev$P)

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
#whiteblack <- c("white", "black")




title <- "Complete Correlation Matrix"
col <- colorRampPalette(c("grey", "white", "green"))
#colnames(M) <- c("BBB", "Time Point", "Area Deficit", "Injury", "Max Time Delay")
mod_corr <- corrplot(res_complete$r, method = "color", 
                     type = "upper", order = "hclust", col = col(200), number.cex = .8,
                     title=title,
                     addCoef.col = "black", # Add coefficient of correlation
                     tl.col = "black", tl.srt = 90, # Text label color and rotation
                     # Combine with significance
                     #p.mat = res_mod$P, sig.level = 0.05, insig = "blank",
                     mar=c(0,0,2,0),
                     # hide correlation coefficient on the principal diagonal
                     diag = TRUE)
#ggsave("Full Correlation matrix.pdf")

#aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
#scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
#guides(color = FALSE, alpha = FALSE)


title <- "Moderate Injury Correlation Matrix"
col <- colorRampPalette(c("grey", "white", "green"))
#colnames(M) <- c("BBB", "Time Point", "Area Deficit", "Injury", "Max Time Delay")
mod_corr <- corrplot(res_mod$r, method = "color", 
         type = "upper", order = "hclus", col = col(200), number.cex = .8,
         title=title,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         #p.mat = res_mod$P, sig.level = 0.05, insig = "blank",
         mar=c(0,0,2,0),
         # hide correlation coefficient on the principal diagonal
         diag = TRUE)
#ggsave("Full Correlation matrix.pdf")

#aes(color = coefficient > 0, alpha = abs(coefficient) > 0.5)) +
  #scale_alpha_manual(values = c("TRUE" = 0.25, "FALSE" = 0)) +
  #guides(color = FALSE, alpha = FALSE)


mod_corr

title <- "Severe Injury Correlation Matrix"
col <- colorRampPalette(c("grey", "white", "green"))
#colnames(M) <- c("BBB", "Time Point", "Area Deficit", "Injury", "Max Time Delay")
sev_corr <- corrplot(res_sev$r, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .8,
         title=title,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         #p.mat = res_sev$P, sig.level = 0.05, insig = "blank",
         mar=c(0,0,2,0),
         # hide correlation coefficient on the principal diagonal
         diag = TRUE)

#grid.arrange(mod_corr, sev_corr)

