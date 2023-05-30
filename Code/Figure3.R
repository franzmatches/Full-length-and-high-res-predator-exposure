########################################################################################################################
# Figure 3 #
########################################################################################################################
require(brms)
require(tidyverse)
brms_lm_id <- readRDS(file = "Results/brms_lm_id_noar.rds")

did_id_data <- readRDS(file = "Data/didinium_data_IDs_clean.RDS")
hom_id_data <- readRDS(file = "Data/homalozoon_data_IDs_clean.RDS")
prey_id_data <- readRDS(file = "Data/prey_data_IDs_clean.RDS")

#add predator treatment column to each then combine data frames
did_id_data <- did_id_data %>%
  mutate(predator_treatment = "didinium")

hom_id_data <- hom_id_data %>%
  mutate(predator_treatment = "homalozoon")

prey_id_data <- prey_id_data %>%
  mutate(predator_treatment = "prey")

#combine the data frames
id_data_with_predators <- rbind(did_id_data, hom_id_data, prey_id_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))%>%
  mutate(
    # mean_speed_norm = predict(bestNormalize::bestNormalize(mean_speed)), we dont need to normalize
    treatment = factor(treatment, levels = c(15,25)),
    treat_inter = as.factor(interaction(treatment,predator_treatment)))

id_data<-id_data_with_predators %>% filter(Species == "PARcau") %>%
  group_by(ID,treat_inter,replicate,time_point) %>%
  slice_head()

############ 
#Extract conditional effects
############ 

cond_inter_treatment_id<- conditional_effects(brms_lm_id,effects = "mean_length_um:treatment",
                                              conditions = data.frame(
                                                expand.grid(
                                                  time_point = seq(0,24,by=8),
                                                  predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]%>%
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

############ 
#Save figure
############

ggsave("Results/figures/figure3.png",
       ggplot(cond_inter_treatment_id %>%
                mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h")))%>%
                mutate(effect2__ = paste0(effect2__,"\u00B0C")) %>%
                mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","D. nasutum","H. vermiculare"))))+
         geom_point(data = id_data %>%
                      ungroup()%>%
                      dplyr::mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","D. nasutum","H. vermiculare"))) %>%
                      mutate(effect2__ = paste0(treatment,"\u00B0C")) %>%
                      filter(time_point %in% c(0,8,16,24)) %>%
                      mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h"))),
                    aes(x=mean_length_um,y=mean_width_um,col=as.factor(effect2__)), alpha = 0.4)+
         geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
         geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
         facet_grid(cond__~predator_treatment, scales = "fixed") +
         scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
         scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
         xlab(expression(paste("Mean length (", mu,"m)")))+
         ylab(expression(paste("Mean width (", mu,"m)")))+ 
         theme_classic()+
         theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
               strip.text = element_text(face = "bold"),
               strip.text.y.right = element_text(angle = 0),
               strip.text.x = element_text(face = "bold.italic"),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               panel.border = element_rect(fill = NA, colour = "black")),
       width=5,height =5,dpi=300)
