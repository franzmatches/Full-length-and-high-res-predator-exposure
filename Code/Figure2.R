########################################################################################################################
# Figure 2 #
########################################################################################################################
require(brms)
require(tidyverse)
brms_id_weibullk8 <- readRDS(file = "Results/brms_id_noar_weibullk8.rds")

############ 
#Extract conditional effects
############ 

cond_speed_treatment_id <- brms::conditional_effects(brms_id_weibullk8,
                                                     effects = "time_point",
                                                     method = "posterior_epred",
                                                     conditions = (
                                                       data.frame(
                                                         expand.grid(treatment = c(15,25),
                                                                     predator_treatment =  c("prey","didinium","homalozoon"))) %>%
                          dplyr::mutate(treat_inter = interaction(treatment,predator_treatment))))[[1]]

############ 
#Save figure
############

ggsave("Results/figures/figure2.png",
       ggplot(cond_speed_treatment_id %>%
                mutate(treatment = paste0(treatment,"\u00B0C"))%>%
                mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
         geom_point(data = id_data %>%
                      ungroup()%>%
                      mutate(treatment = paste0(treatment,"\u00B0C"))%>%
                      mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
                    aes(x=time_point,y=mean_speed),col="black", alpha = 0.3, shape = 21)+
         geom_line(aes(x = effect1__, y=estimate__),col="black",linewidth=1.5) +
         geom_ribbon(aes(x = time_point,ymin = lower__, ymax =  upper__,fill=as.factor(treatment)),alpha=0.3)+
         facet_grid(treatment~predator_treatment, scales = "fixed")+
         scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
         scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
         theme_classic()+
         xlab("Time (hours)")+
         ylab("Mean speed (mm/s)")+
         theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
               strip.text = element_text(face = "bold"),
               strip.text.y.right = element_text(angle = 0),
               strip.text.x = element_text(face = "bold.italic"),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black"),
               panel.border = element_rect(fill = NA, colour = "black")),
       width=8,height =6,dpi=300)