require(mgcv)

did_id_data <- read.csv("Data/did_data_IDs.csv")
hom_id_data <- read.csv("Data/hom_data_IDs.csv")
prey_id_data <- read.csv("Data/prey_data_IDs.csv")

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
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))%>%
  mutate(mean_speed_norm = predict(bestNormalize::bestNormalize(mean_speed)),
         treatment = factor(treatment, levels = c(15,25)),
         treat_inter = as.factor(interaction(treatment,predator_treatment)))

speed_id_plot <-
  
  ggplot(id_data, aes(x = time_point, y = mean_speed, col = as.factor(predator_treatment)))+
  #geom_point(alpha = 0.2) +
  geom_smooth(method="gam", formula = y ~s(x,bs="tp"),se = FALSE) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

ggpubr::ggarrange(speed_id_plot,tt_plot,legend = "none")
########################################################################################
# Standard Gams
########################################################################################
gam1 <- gam(mean_speed ~ s(time_point) + 
              treatment:predator_treatment + 
              s(time_point,by = treatment) +
              s(time_point,by = predator_treatment) + 
              s(time_point,by = treat_inter) ,
            data = id_data, family = "gaussian")
summary(gam1)
gam.check(gam1)
plot.gam(gam1,pages = 3,all.terms = T,seWithMean = TRUE, shift = coef(gam1)[1],residuals = T)



gam1_rawdata <- gam(mean_speed ~ s(time_point) + treatment:predator_treatment + s(time_point,by = treatment) +
              s(time_point,by = predator_treatment)
              + s(time_point,by = treat_inter)
              ,data = id_data, family = "gaussian")
summary(gam1_rawdata)

gam.check(gam1_rawdata)
plot.gam(gam1_rawdata,pages = 3,all.terms = T)


tt <- cbind(id_data,
            fit = predict(gam1_rawdata,type="response",se.fit=T)$fit,
            se = predict(gam1_rawdata,type="response",se.fit=T)$se.fit)

tt_plot <- ggplot(tt, aes(x = time_point, y = mean_speed))+
  #geom_smooth(se = FALSE) +
  geom_line(aes(y=fit,col=as.factor(predator_treatment))) +
  geom_ribbon(aes(ymin = fit - 1.96*se, ymax =  fit + 1.96*se,fill=as.factor(predator_treatment)),alpha=0.5)+
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  scale_fill_discrete(guide="none")+
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))




p_df <- ggeffects::ggpredict(gam1,terms = c("time_point","treat_inter"),condition = "time_point")

plot(ggeffects::ggpredict(gam1,terms = c("time_point","treat_inter [15.prey, 15.didinium, 15.homalozoon]")),add.data = F)
plot(ggeffects::ggpredict(gam1,terms = c("time_point","treat_inter [25.prey, 25.didinium, 25.homalozoon]")),add.data = F)

plot(ggeffects::ggpredict(gam1,terms = c("time_point","treatment","predator_treatment")),add.data = F)

########################################################################################
# Mixed Effects Gams
########################################################################################
gam2 <- gamm(mean_speed ~ s(time_point) + treatment:predator_treatment + s(time_point,by = treatment) +
              s(time_point,by = predator_treatment) + 
              s(time_point,by = treat_inter) , 
             #random = list(replicate= ~1),
             correlation=corAR1(form=~1|time_point),
            data = id_data, family = "gaussian",method = "REML")

plot.gam(gam2,pages = 3,all.terms = T,seWithMean = TRUE,shift = coef(gam1)[1])

plot.gam(gam2,pages = 3,all.terms = T)
gam.check(gam2)

