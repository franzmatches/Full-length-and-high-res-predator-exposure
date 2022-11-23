#######################################################################
#########predator exposure full length data analysis with gams#########
#######################################################################
#load packages 
library(dplyr)
library(tidyverse)
library(bestNormalize)
library(MuMIn)
library(rstatix)
require(mgcv)

did_id_data <- read.csv("Data/did_data_IDs_corrected_max_abundance.csv")
hom_id_data <- read.csv("Data/hom_data_IDs_corrected_max_abundance.csv")
prey_id_data <- read.csv("Data/prey_data_IDs_correct_max_abundance.csv")

#add predator treatment column to each then combine data frames
did_id_data <- did_id_data %>%
  mutate(predator_treatment = "didinium")

hom_id_data <- hom_id_data %>%
  mutate(predator_treatment = "homalozoon")

prey_id_data <- prey_id_data %>%
  mutate(predator_treatment = "prey")

#combine the data frames
id_data <- rbind(did_id_data, hom_id_data, prey_id_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))%>%
  mutate(mean_speed_norm = predict(bestNormalize::bestNormalize(mean_speed)),
         treatment = factor(treatment, levels = c(15,25)),
         treat_inter = as.factor(interaction(treatment,predator_treatment)))

speed_id_plot <-
  ggplot(id_data, aes(x = time_point, y = mean_speed, col = as.factor(predator_treatment)))+
  #geom_jitter(width = 1) +
  geom_smooth(method="gam", formula = y ~s(x,bs="tp"),se = T) +
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


head(id_data)

means<-id_data %>%
  group_by(treatment, time_point, predator_treatment ) %>%
  summarise(speed=var(mean_roundness))

head(means)

ggplot(means, aes(x=time_point, y=speed, col=predator_treatment))+geom_line()+facet_wrap(~treatment)

ggpubr::ggarrange(speed_id_plot,tt_plot,legend = "none")
########################################################################################
# Standard Gams
########################################################################################
pooled_id_data <- id_data %>% 
  group_by(replicate,treat_inter,treatment, predator_treatment,time_point) %>%
  summarise(avg_speed = mean(mean_speed),
         avg_roundness = mean(mean_roundness),
         avg_length = mean(mean_length_um),
         avg_width = mean(mean_width_um),
         avg_area = mean(mean_area_sqrd_um))  %>%
  mutate(time_point = time_point + 1) %>%
  arrange(time_point,.by_group = T)
  

gam1 <- gam(mean_speed ~ s(time_point) + 
              treatment*predator_treatment + 
              s(time_point,by = treatment) +
              s(time_point,by = predator_treatment) + 
              s(time_point,by = treat_inter) ,
            data = id_data, family = "gaussian")
summary(gam1)
par(mfrow = c(2,2))
gam.check(gam1)
par(mfrow = c(1,1))

##autocorrelation of gam1 residuals
residual_df <- cbind(subset(id_data,select = c("time_point","ID","treatment","replicate","predator_treatment","treat_inter")),
                     "resid" = resid(gam1)) %>%
  #nest_by(treat_inter,replicate) %>%
  nest(data = -c(treat_inter,replicate,ID)) %>%
  mutate(acf_val = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-1),plot = F)$acf),
         acf_lag = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-1),plot = F)$lag),
         ci = 1.96/sqrt(n()-1)) %>%
  unnest(c(acf_val,acf_lag,ci))

ggplot(residual_df,aes(x=acf_lag,y=acf_val)) + 
  geom_segment(aes(yend=0,xend=acf_lag)) +
  geom_hline(aes(yintercept = ci),col="blue",linetype = "dashed")+
  geom_hline(aes(yintercept = -ci),col="blue",linetype = "dashed")+
  facet_grid(replicate~treat_inter)

plot.gam(gam1,pages = 3,all.terms = T,seWithMean = TRUE, shift = coef(gam1)[1],residuals = T)



gam1_rawdata <- gam(mean_speed ~ s(time_point) + 
                      treatment*predator_treatment + 
                      s(time_point,by = treatment) +
              s(time_point,by = predator_treatment)
              + s(time_point,by = treat_inter)
              ,data = id_data, family = "gaussian")
summary(gam1_rawdata)

par(mfrow = c(2,2))
gam.check(gam1_rawdata)
par(mfrow = c(1,1))

plot.gam(gam1_rawdata,pages = 3,all.terms = T)

pred_gam1 <- cbind(id_data,
            fit = predict(gam1,type="response",se.fit=T)$fit,
            se = predict(gam1,type="response",se.fit=T)$se.fit)

pred_gam1_plot <- ggplot(pred_gam1, aes(x = time_point, y = mean_speed))+
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

gam2 <- gam(avg_speed ~ s(time_point) + 
              treatment*predator_treatment + 
              s(time_point,by = treatment) +
              s(time_point,by = predator_treatment) + 
              s(time_point,by = treat_inter) ,
            data = pooled_id_data, family = "gaussian")

par(mfrow = c(2,2))
gam.check(gam2)
par(mfrow = c(1,1))
plot.gam(gam2,pages = 3,all.terms = T,seWithMean = TRUE, shift = coef(gam1)[1],residuals = T)

##autocorrelation of gam1 residuals
residual_df2 <- cbind(subset(pooled_id_data,select = c("time_point","replicate","treat_inter","treatment","predator_treatment")),
                     "resid" = resid(gam2)) %>%
  nest(data = -c(replicate,treat_inter)) %>%
  mutate(acf_val = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-1),plot = F)$acf),
         acf_lag = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-1),plot = F)$lag),
         ci = 1.96/sqrt(n()-1)) %>%
  unnest(c(acf_val,acf_lag,ci))

ggplot(residual_df2,aes(x=acf_lag,y=acf_val)) + 
  geom_segment(aes(yend=0,xend=acf_lag)) +
  #geom_hline(aes(yintercept = ci),col="blue",linetype = "dashed")+
  #geom_hline(aes(yintercept = -ci),col="blue",linetype = "dashed")+
  geom_hline(yintercept = 0.1,col="blue",linetype = "dashed")+
  geom_hline(yintercept = -0.1,col="blue",linetype = "dashed")+
  facet_grid(replicate~treat_inter)


pred_gam2 <- cbind(pooled_id_data,
                   fit = predict(gam2,type="response",se.fit=T)$fit,
                   se = predict(gam2,type="response",se.fit=T)$se.fit)

pred_gam2_plot <- ggplot(pred_gam2, aes(x = time_point, y = avg_speed))+
  #geom_point() +
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

plot(ggeffects::ggpredict(gam1_rawdata,terms = c("time_point","treatment","predator_treatment")),add.data = F)

########################################################################################
# Mixed Effects Gams
########################################################################################
id_data_gamm <- id_data %>%
  mutate(time_point_unique = paste(time_point,treat_inter,ID,sep="_"))

gamm1 <- gamm(mean_speed ~ s(time_point) + treatment:predator_treatment + s(time_point,by = treatment) +
              s(time_point,by = predator_treatment) + 
              s(time_point,by = treat_inter) , 
             random = list(replicate= ~time_point),
            data = id_data_gamm, family = "gaussian",method = "REML")

gamm1 <- gamm4::gamm4(mean_speed ~ s(time_point) + 
                       treatment:predator_treatment + 
                       s(time_point,by = treatment) +
                       s(time_point,by = predator_treatment) + 
                       s(time_point,by = treat_inter), 
             random = ~(time_point|replicate:treat_inter),
             data = id_data, family = "gaussian",REML = T)

summary(gamm1$mer)
par(mfrow = c(2,2))
plot(gamm1$mer)
par(mfrow = c(1,1))

residual_gamm_df <- cbind(subset(id_data,select = c("time_point","ID","treatment","replicate","predator_treatment","treat_inter")),
                     "resid" = resid(gamm1$mer)) %>%
  #nest_by(treat_inter,replicate) %>%
  nest(data = -c(treat_inter,replicate)) %>%
  mutate(acf_val = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-1), plot = F)$acf),
         acf_lag = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-1),plot = F)$lag),
         ci = 1.96/sqrt(n()-1)) %>%
  unnest(c(acf_val,acf_lag,ci))

ggplot(residual_gamm_df,aes(x=acf_lag,y=acf_val)) + 
  geom_segment(aes(yend=0,xend=acf_lag)) +
  geom_hline(aes(yintercept = ci),col="blue",linetype = "dashed")+
  geom_hline(aes(yintercept = -ci),col="blue",linetype = "dashed")+
  facet_grid(replicate~treat_inter)


pred_gamm1_df <- cbind(id_data,
            fit = predict(gamm1$mer,type="response"))

ggplot(pred_gamm1_df, aes(x = time_point, y = mean_speed,group=replicate))+
  #geom_smooth(se = FALSE) +
  geom_line(aes(y=fit,col=as.factor(predator_treatment))) +
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



gamm2 <- gamm(avg_speed ~ s(time_point) , 
              #random = list(replicate= ~1),
              correlation = corAR1(form = ~time_point|replicate:treat_inter),
              data = pooled_id_data, family = "gaussian",method = "REML")

gamm2 <- gamm4::gamm4(avg_speed ~ s(time_point) + 
                        treatment*predator_treatment + 
                        s(time_point,by = treatment) +
                        s(time_point,by = predator_treatment) + 
                        s(time_point,by = treat_inter), 
                      random = ~(time_point|replicate:treat_inter),
                      data = pooled_id_data, family = "gaussian",REML = T)

summary(gamm2$mer)
par(mfrow = c(2,2))
plot(gamm2$mer)
par(mfrow = c(1,1))

residual_gamm2_df <- cbind(subset(pooled_id_data,select = c("time_point","replicate","treatment","predator_treatment","treat_inter")),
                          "resid" = resid(gamm2$mer)) %>%
  nest(data = -c(replicate,treat_inter)) %>%
  mutate(acf_val = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-5), plot = F)$acf),
         acf_lag = purrr::map(data, ~ acf(.x$resid, lag.max = (max(.x$time_point)-5),plot = F)$lag),
         ci = 1.96/sqrt(n()-1)) %>%
  unnest(c(acf_val,acf_lag,ci))

ggplot(residual_gamm2_df,aes(x=acf_lag,y=acf_val)) + 
  geom_segment(aes(yend=0,xend=acf_lag)) +
  geom_hline(aes(yintercept = ci),col="blue",linetype = "dashed")+
  geom_hline(aes(yintercept = -ci),col="blue",linetype = "dashed")+
  facet_grid(replicate~treat_inter)

pred_gamm2_df <- cbind(pooled_id_data,
                       fit = predict(gamm2$mer,type="response"))

ggplot(pred_gamm2_df, aes(x = time_point, y = avg_speed,group=replicate))+
  geom_point(aes(col=as.factor(replicate)),alpha =0.4) +
  #geom_line(aes(y=fit,col=as.factor(predator_treatment))) +
  geom_line(aes(y=fit,col=as.factor(replicate))) +
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

