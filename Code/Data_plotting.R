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
library(ggh4x)
library(ggpubr)
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


#remove those data where max_abundance was counted above reasonable numbers
id_data<-id_data %>% filter(max_abundance <35)



#####Plotting parameters through time#########


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




#plotting mean_roundness among replicates and the general treatment means
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = max_length , group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(mean_roundness_rep = mean(mean_roundness)),
            aes(x = time_point, y = mean_roundness_rep, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_roundness_rep = mean(mean_roundness)),
            aes(x = time_point, y = mean_roundness_rep, col = as.factor(predator_treatment)), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment, scales = "free")+
  xlab("Time (hours)")+
  ylab("Mean roundness Variance")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()



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


#-----------------------------------------------------------------------#
####Plotting relationships of traits with the abundance in the patch#####
##let's plot first the time series of abundance
#single treatment mean of abundance
ggplot(data = id_data %>% group_by(time_point,
                                    treatment,
                                    predator_treatment) %>% 
          summarize(mean_max_ab = mean(max_abundance)),
        aes(x = time_point, y = mean_max_ab, col = predator_treatment,fill = predator_treatment))+
  geom_point(alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(alpha = .2)+
  theme_bw()

#see what happens in each replicate  
ggplot(data = id_data %>% group_by(time_point,
                                   treatment,
                                   predator_treatment,
                                   replicate) %>% 
         summarise(max_abundance_rep = max(max_abundance)),
       aes(x = time_point, y = max_abundance_rep,
           col = predator_treatment,fill = predator_treatment))+
  geom_point(alpha = .5)+
  facet_nested_wrap(~treatment*replicate, scales = "free_y")+
  geom_smooth(alpha = .2)+
  # geom_smooth(method = "lm", alpha = .2)+
  theme_bw()
####-----------------using all ID measurements in each replicate---------------------####
ggplot(data = id_data ,
       aes(x = max_abundance, y = mean_length_um, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_jitter(width = 0.5, shape = 1, alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()


#speed
ggplot(data = id_data,
       aes(x = max_abundance, y = mean_speed, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_jitter(width = 0.5, shape = 1, alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()


#area
ggplot(data = id_data,
       aes(x = max_abundance, y = mean_area_sqrd_um, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_jitter(width = 0.5, shape = 1, alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()


#roundness (calculated as roundness = (max(Length_um)/ mean(Width_um)) for each ID)
ggplot(data = id_data,
       aes(x = max_abundance, y = mean_roundness, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_jitter(width = 0.5, shape = 1, alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # geom_smooth(alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()




####-----------------using the across replicates means---------------------####
#length
ggplot(data = id_data %>% group_by(time_point,
                                   treatment,
                                   predator_treatment) %>% 
         summarize(mean_max_ab = mean(max_abundance),
                   mean_mean_length = mean(mean_length_um)),
       aes(x = mean_max_ab, y = mean_mean_length, 
                           col = predator_treatment,
           fill = predator_treatment))+
  geom_point(alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()

glimpse(id_data)


#speed
ggplot(data = id_data %>% group_by(time_point,
                                   treatment,
                                   predator_treatment) %>% 
         summarize(mean_max_ab = mean(max_abundance),
                   mean_mean_speeds = mean(mean_speed)),
       aes(x = mean_max_ab, y = mean_mean_speeds, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_point(alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()


#area
ggplot(data = id_data %>% group_by(time_point,
                                   treatment,
                                   predator_treatment) %>% 
         summarize(mean_max_ab = mean(max_abundance),
                   mean_mean_area = mean(mean_area_sqrd_um)),
       aes(x = mean_max_ab, y = mean_mean_area, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_point(alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()


#roundness (calculated as roundness = (max(Length_um)/ mean(Width_um)) for each ID)
plot_round_avg<-ggplot(data = id_data %>% group_by(time_point,
                                   treatment,
                                   predator_treatment) %>% 
         summarize(mean_max_ab = mean(max_abundance),
                   mean_mean_roundness = mean(mean_roundness)),
       aes(x = mean_max_ab, y = mean_mean_roundness, 
           col = predator_treatment,
           fill = predator_treatment))+
  geom_point(alpha = .5)+
  facet_grid(~treatment)+
  geom_smooth(method = "lm", alpha = .2)+
  # geom_smooth(alpha = .2)+
  # scale_x_continuous(limits = c(0,35))+
  theme_bw()






####Fit a model on the control and scale the treatments based on that###
#15 degrees

stnd_lm_15_fit<-lm(mean_speed_treat~time_point,
            data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_speed_treat = mean(mean_speed)) %>%
              filter(predator_treatment == "prey",
                     treatment == 15))$fitted



stnd_lm_25_fit<-lm(mean_speed_treat~time_point,
               data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
                 summarise(mean_speed_treat = mean(mean_speed)) %>%
                 filter(predator_treatment == "prey",
                        treatment == 25))$fitted

##data wrangling to substract the model fitted on the average mean_speed of the control to each replicate of each treatment
data_reps<-id_data %>% group_by(time_point,
                                treatment,
                                predator_treatment,
                                replicate) %>% 
  summarize(mean_speed_rep = mean(mean_speed)) %>% 
  group_by(predator_treatment,replicate) %>%
  nest() %>%
  mutate(data = purrr::map(data, function(x){
    x <- x %>% 
      arrange(treatment) %>%
      left_join(data.frame(time_point = rep(0:24,2), detrend_ts = c(stnd_lm_15_fit,stnd_lm_25_fit),
                           treatment = c(rep(15,25),rep(25,25))),
                           by = c("time_point","treatment")) %>%
      #group_by(treatment) %>%
      mutate(detrend_speed = mean_speed_rep - detrend_ts)
    return(x)
  })) %>%
  unnest(data)
  

##plot the difference from such model fit over time
ggplot()+
  geom_line(data = data_reps,
            aes(x = time_point, y = detrend_speed, col = as.factor(predator_treatment),
                group = replicate), alpha = .3)+
  facet_grid(predator_treatment~treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw()






