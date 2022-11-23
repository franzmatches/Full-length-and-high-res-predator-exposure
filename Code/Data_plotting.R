rm(list=ls())
###############################################################
######### predator exposure full length data plotting #########
###############################################################
#load packages
library(dplyr)
library(ggplot2)
library(bestNormalize) #transforming response variables to fit gaussian model
library(MuMIn) #model comparison and simplification
library(rstatix) 
library(glmmTMB) #GLMMs
library(DHARMa) #GLMM diagnostics
library(fitdistrplus) #data to distribution fit
library(plotrix)

#load in the tracked ID data 

did_id_data <- read.csv("Data/did_data_IDs_corrected_max_abundance.csv")
hom_id_data <- read.csv("Data/hom_data_IDs_corrected_max_abundance.csv")
prey_id_data <- read.csv("Data/prey_data_IDs_correct_max_abundance.csv")

#add predator treatment column to each then combine data frames
did_id_data <- did_id_data %>%
  mutate(predator_treatment = rep("didinium", times = nrow(did_id_data)))

hom_id_data <- hom_id_data %>%
  mutate(predator_treatment = rep("homalozoon", times = nrow(hom_id_data)))

prey_id_data <- prey_id_data %>%
  mutate(predator_treatment = rep("prey", times = nrow(prey_id_data)))

  #combine the data frames
id_data <- rbind(did_id_data, hom_id_data, prey_id_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))

#plotting max_abundance
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = mean_speed, group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(max_abundance = max(max_abundance)),
            aes(x = time_point, y = max_abundance, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_max_abundance = mean(max_abundance)),
            aes(x = time_point, y = mean_max_abundance, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "free")+
  xlab("Time (hours)")+
  ylab("Max abundance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()


ggsave("single_replicates_max_abundances_and_treatments_means_max_abundance.pdf", units="in", width=16, height=10) 




#plotting mean speed among replicates and the general treatment means
ggplot()+
    # geom_point(data = id_data, aes(x = time_point, y = mean_speed, group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
    geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(mean_speed_rep = mean(mean_speed)),
            aes(x = time_point, y = mean_speed_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_speed_treat = mean(mean_speed)),
            aes(x = time_point, y = mean_speed_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()
  # theme_classic()
  # theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
  #       panel.background = element_blank(),
  #       axis.line = element_line(colour = "black"),
  #       aspect.ratio = 1,
  #       panel.border = element_rect(fill = NA, colour = "black"))

ggsave("single_replicates_means_and_treatments_means_speed.pdf", units="in", width=16, height=10) 
ll


#mean speed with error bars
data_error_bars<-id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
  summarise(mean_speed_rep = mean(mean_speed),
            se_mean_speed_rep = std.error(mean_speed))

ggplot(data = data_error_bars,
       aes(x = time_point, y = mean_speed_rep, col = as.factor(predator_treatment),
           group = replicate))+
  facet_grid(predator_treatment~treatment, scales = "free")+
  geom_point(alpha = .5, size = 0.5,
             position = position_dodge(.8)) +
  geom_errorbar(aes(ymin = mean_speed_rep - 1.96*se_mean_speed_rep,
                    ymax = mean_speed_rep + 1.96*se_mean_speed_rep), width = .4,
                size = .4,
                position = position_dodge(.8)
                , alpha = .9
                ) +
  geom_line(alpha = .4,
            position = position_dodge(.8)) +
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("single_replicates_means_and_errorbars_speed.pdf", units="in", width=16, height=10) 




#plotting mean speed VARIANCE among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = mean_speed, group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(var_mean_speed_rep = var(mean_speed)),
            aes(x = time_point, y = var_mean_speed_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(var_mean_speed_treat = var(mean_speed)),
            aes(x = time_point, y = var_mean_speed_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "free")+
  xlab("Time (hours)")+
  ylab("Mean speed Variance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()
# theme_classic()
# theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
#       panel.background = element_blank(),
#       axis.line = element_line(colour = "black"),
#       aspect.ratio = 1,
#       panel.border = element_rect(fill = NA, colour = "black"))

ggsave("VARIANCE_single_replicates_means_and_treatments_means_speed.pdf", units="in", width=16, height=10)


#plotting mean_length_um among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = mean_length_um , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(mean_length_rep = mean(mean_length_um)),
            aes(x = time_point, y = mean_length_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_length_treat = mean(mean_length_um)),
            aes(x = time_point, y = mean_length_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean Length")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("single_replicates_means_and_treatments_means_length.pdf.pdf", units="in", width=16, height=10) 

#plotting mean_length_um VARIANCE among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = mean_length_um , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(var_mean_length_rep = var(mean_length_um)),
            aes(x = time_point, y = var_mean_length_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(variance_mean_length_treat = var(mean_length_um)),
            aes(x = time_point, y = variance_mean_length_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, 
             scales = "free"
             )+
  xlab("Time (hours)")+
  ylab("Mean Length Variance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("VARIANCE_single_replicates_means_and_treatments_means_length.pdf.pdf", units="in", width=16, height=10) 


#plotting max_length among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(mean_max_length_rep = mean(max_length)),
            aes(x = time_point, y = mean_max_length_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_max_length_rep_treat = mean(max_length)),
            aes(x = time_point, y = mean_max_length_rep_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Max Length")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("single_replicates_means_and_treatments_means_max_length.pdf", units="in", width=16, height=10) 


#plotting mean_width_um among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(mean_mean_width_um_rep = mean(mean_width_um)),
            aes(x = time_point, y = mean_mean_width_um_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_mean_width_um_rep_treat = mean(mean_width_um)),
            aes(x = time_point, y = mean_mean_width_um_rep_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean Width")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("single_replicates_means_and_treatments_means_width_um.pdf", units="in", width=16, height=10) 


#plotting mean_width_um VARIANCE among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(var_mean_width_um_rep = var(mean_width_um)),
            aes(x = time_point, y = var_mean_width_um_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(var_mean_width_um_rep_treat = var(mean_width_um)),
            aes(x = time_point, y = var_mean_width_um_rep_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "free")+
  xlab("Time (hours)")+
  ylab("Mean Width Variance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("VARIANCE_single_replicates_means_and_treatments_means_width_um.pdf", units="in", width=16, height=10) 


#plotting mean_area_sqrd_um among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(mean_area_sqrd_um_rep = mean(mean_area_sqrd_um)),
            aes(x = time_point, y = mean_area_sqrd_um_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_area_sqrd_um_treat = mean(mean_area_sqrd_um)),
            aes(x = time_point, y = mean_area_sqrd_um_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean area")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("single_replicates_means_and_treatments_means_area.pdf", units="in", width=16, height=10) 


#plotting mean_area_sqrd_um variance among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(var_area_sqrd_um_rep = var(mean_area_sqrd_um)),
            aes(x = time_point, y = var_area_sqrd_um_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(var_area_sqrd_um_treat = var(mean_area_sqrd_um)),
            aes(x = time_point, y = var_area_sqrd_um_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "free")+
  xlab("Time (hours)")+
  ylab("Mean area Variance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("VARIANCE_single_replicates_means_and_treatments_means_area.pdf", units="in", width=16, height=10) 



#plotting mean_roundness VARIANCE among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(var_roundness_rep = var(mean_roundness)),
            aes(x = time_point, y = var_roundness_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(var_roundness_treat = var(mean_roundness)),
            aes(x = time_point, y = var_roundness_treat, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "free")+
  xlab("Time (hours)")+
  ylab("Mean roundness Variance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()

ggsave("VARIANCE_single_replicates_means_and_treatments_means_roundness.pdf", units="in", width=16, height=10) 


