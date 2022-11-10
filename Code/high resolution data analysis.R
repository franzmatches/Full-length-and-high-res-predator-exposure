###########################################################
######### high resolution experiment data analysis#########
###########################################################
rm(list=ls())

#load packages


##test test github

#load packages 
library(dplyr)
library(tidyverse)
library(bestNormalize)
library(MuMIn)
library(rstatix)

#load in the treatment x replicate level data
did_mean_data <- read.csv("Data/did_data_means.csv")
hom_mean_data <- read.csv("Data/hom_data_means.csv")
prey_mean_data <- read.csv("Data/prey_data_means.csv")

#add predator treatment column to each then combine data frames
did_mean_data <- did_mean_data %>%
  mutate(predator_treatment = rep("didinium", times = nrow(did_mean_data)))

hom_mean_data <- hom_mean_data %>%
  mutate(predator_treatment = rep("homalozoon", times = nrow(hom_mean_data)))

prey_mean_data <- prey_mean_data %>%
  mutate(predator_treatment = rep("prey", times = nrow(prey_mean_data)))

#combine the data frames
mean_data <- rbind(did_mean_data, hom_mean_data, prey_mean_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon"))) %>% 
  mutate(treatment = factor(treatment, levels = c(15,25)),
         treat_inter = as.factor(interaction(treatment,predator_treatment))) %>% 
  #change the scaling of the time point so that there is an hour gap between time point 0 and the next
  mutate(time_point = case_when(time_point == 0 ~ '0',
                                TRUE ~ paste(time_point + 50))) %>%
  mutate(time_point = as.numeric(time_point))


#create dataset without time_point zero

#visualise the raw data
speed_plot_raw <-
  ggplot(mean_data, aes(x = time_point, y = average_mean_speed, col = as.factor(predator_treatment)))+
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (minutes)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

width_plot_raw <-
  ggplot(mean_data, aes(x = time_point, y = average_mean_width, col = as.factor(predator_treatment)))+
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (minutes)")+
  ylab("Mean width")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

length_plot_raw <-
  ggplot(mean_data, aes(x = time_point, y = max_length, col = as.factor(predator_treatment)))+
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (minutes)")+
  ylab("Max length")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

roundness_plot_raw <- 
  ggplot(mean_data, aes(x = time_point, y = average_roundness, col = as.factor(predator_treatment)))+
  geom_point(alpha = 0.2) +
  geom_smooth(se = FALSE) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (minutes)")+
  ylab("Average roundness")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))


##fitting a gam###


gam1_hr <- gam(average_mean_speed ~ s(time_point) + treatment:predator_treatment + s(time_point,by = treatment) +
                      s(time_point,by = predator_treatment)
                    + s(time_point,by = treat_inter)
                    ,data = mean_data, family = "gaussian")
summary(gam1_hr)

gam.check(gam1_hr)
plot.gam(gam1_hr,pages = 3,all.terms = T)


tt <- cbind(mean_data,
            fit = predict(gam1_hr,type="response",se.fit=T)$fit,
            se = predict(gam1_hr,type="response",se.fit=T)$se.fit)

ggplot(tt, aes(x = time_point, y = average_mean_speed))+
  # geom_smooth(se = FALSE) +
  # geom_line(aes(y=fit)) +
  geom_ribbon(aes(ymin = fit - 1.96*se, ymax =  fit + 1.96*se,fill=as.factor(predator_treatment)))+
  facet_grid(~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

