require(brms)
require(Matrix)
require(tidyverse)
require(camcorder)

brms_lm <- readRDS("Results/brms_lm.rds")

cond_inter_treatment2 <- conditional_effects(brms_lm,effects = "mean_length_um:treatment",
                                             conditions = data.frame(expand.grid(time_point = unique(id_data$time_point),
                                                                                 predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

camcorder::gg_record(
  dir = file.path("/Users/ul20791/Downloads/camcorder_test", "recording100"), # where to save the recording
  device = "png", # device to use to save images
  width = 5,      # width of saved image
  height = 2,     # height of saved image
  units = "in",   # units for width and height
  dpi = 300       # dpi to use when saving image
)


ggplot(cond_inter_treatment2 |>
         mutate(cond__ = paste0("t",time_point)) |>
         filter(cond__ == paste0("t",0)))+
  geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
  geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
  facet_grid(cond__~predator_treatment, scales = "fixed") +
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
 xlab("Mean length")+
  ylab("Mean width")+
  coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",1)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",2)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",3)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",4)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",5)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",6)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",7)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",8)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))


  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",9)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",10)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",11)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",12)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",13)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",14)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",15)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",16)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",17)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",18)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",19)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",20)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",21)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",22)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",23)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black")) 
  
  ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0("t",time_point)) |>
           filter(cond__ == paste0("t",24)))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length")+
    ylab("Mean width")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
  
          panel.border = element_rect(fill = NA, colour = "black")) 
  
camcorder::gg_playback(
  name = file.path("/Users/ul20791/Downloads/camcorder_test","vignette_gif.gif"),
  first_image_duration = 5,
  last_image_duration = 10,
  frame_duration = 0.5,
  last_as_first = FALSE
)

camcorder::gg_stop_recording()

