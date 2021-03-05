getwd()
setwd('D:/UW/CEUS_Analysis')


library(ggplot2)
library(ggpubr)
library(dplyr)
library(car)
library(plotly)
library(forcats)
library(knitr)
library(formattable)
library(tidyr)
library(gridExtra)
library(grid)
library(lattice)
# Library
library(dygraphs)
library(xts)          # To make the convertion data-frame / xts format
library(tidyverse)
library(lubridate)
#install.packages("readxl")
library(readxl)


## Load and clean Data ##
my_data <- read.csv("Correlation_DF_Matrix_final.csv", header = T)


## Rename cols ##
colnames(my_data) <- c("Injury", "BL Time Delay (sec)", "Acute Time Delay (sec)", "Chronic Time Delay (sec)", "Acute Area Deficit (mm2)", "Chronic Area Deficit (mm2)", "BL BBB", "Acute BBB", "Chronic BBB", "Spared Tissue (mm)", "Chronic Tortuosity", "Caudal BL Angle", "Caudal Acute Angle", "Caudal Diff", "Rostral BL Angle", "Rostral Acute Angle", "Rostral Diff")


## Exploring Data 

update_df <- my_data %>%
  group_by(Injury) %>%
  lapply(function(x) paste("1", x, sep="_")) %>%
  # transmute(Injury, 
  #           recov_thresh = case_when(`Chronic BBB` > 14 ~ "Above",
  #                                    TRUE ~ "Below")) %>%
  mutate_if(is.character, factor)



my_data[1:5]



 

my_data %>%
  group_by(Injury) %>%
  mutate(Injury_num = paste(Injury, row_number(), sep = "_")) %>%
  mutate(num_Injury = paste(row_number(), Injury, sep = "_")) %>%
  ungroup()




library(skimr)
skim(update_df)



## Classification Model ##
library(rsample)
library(tidymodels)


my_data$Injury <- as.factor(my_data$Injury)


glm_fit <- logistic_reg() %>%
  set_engine("glm") %>%
  fit(Injury ~ `Acute Time Delay (sec)`, data = my_data)


## Get the model output - fit method ##
glm_fit


## tidy fit ##
tidy(glm_fit)


## Change injury values to have Injury+1 to each name (Server1, Moderate1, etc).

new_set <- tibble(Injury = unique(my_data$Injury))



lapply(my_data, function(x) paste("1", x, sep="_"))





##### Linear Model #####


acute_fit_AD_BBB <- lm(my_data$`Acute Area Deficit (mm2)` ~ my_data$`Chronic BBB` + my_data$Injury, data = my_data)


## Setting lm variables for predicting pop ##
intercept <- reg_output$coefficients[1]
year_cof <- reg_output$coefficients[2]




## Predict future populations ##
t <- 10:20
predictions <- (t * acute_fit_AD_BBB$coefficients[2]) + acute_fit_AD_BBB$coefficients[1]
predictor_table <- data.frame(t, predictions)
predictor_table



## plot prediction vs actual ##
plot(predict(acute_fit_AD_BBB), my_data$`Acute Area Deficit (mm2)`,
     xlab="predicted", ylab="actual")
abline(a=0, b=1)


## fit ##
fit <- lm(`Chronic BBB` ~ `Acute Area Deficit (mm2)`, data = my_data)
summary(fit)

## plot ##
ggplot(my_data,aes(y=Injury,x=`Acute Area Deficit (mm2)`))+geom_point()+geom_smooth(method="lm")

## multi fit ##
multi_fit <- lm(`Chronic BBB` ~ `Acute Area Deficit (mm2)` + Injury, data = my_data)
summary(multi_fit)


equation1=function(x){coef(multi_fit)[2]*x+coef(multi_fit)[1]}
equation2=function(x){coef(multi_fit)[2]*x+coef(multi_fit)[1]+coef(multi_fit)[3]}

ggplot(my_data,aes(y=`Chronic BBB`,x=`Acute Area Deficit (mm2)`,color=Injury))+geom_point()+
  stat_function(fun=equation1,geom="line",color=scales::hue_pal()(2)[1])+
  stat_function(fun=equation2,geom="line",color=scales::hue_pal()(2)[2])




## multi interactive fit ##
interactive_fit <- lm(`Chronic BBB` ~ `Acute Area Deficit (mm2)`*`Acute Time Delay (sec)`*Injury, data = my_data)
summary(interactive_fit)


## plot ##
ggplot(my_data,aes(y=`Acute Area Deficit (mm2)`,x=`Chronic BBB`,color=`Acute Time Delay (sec)`))+geom_point()+ facet_wrap(~Injury)+stat_smooth(method="lm",se=FALSE)



