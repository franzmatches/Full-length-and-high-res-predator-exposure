require(brms)
require(Matrix)
require(tidyverse)
require(bestNormalize)
require(ggh4x)

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

id_data<-id_data_with_predators %>% filter(Species == "PARcau") |>
  group_by(ID,treat_inter,replicate,time_point) |>
  slice_head()

ggplot(id_data|> group_by(treat_inter,replicate,time_point) |> 
         summarise(num_IDs = max(length(unique(ID)))),
       aes(x=time_point,num_IDs)) +
  geom_point()+
  geom_line()+
  facet_wrap(~treat_inter*replicate)+
  theme_bw()


###################################################################################
# Weibull model
###################################################################################
id_speed_prior <- c(prior(normal(0, 1), class = b),
                    prior(exponential(1), class = Intercept, lb = 0),
                    prior(exponential(1), class = sd),
                    prior(exponential(1), class = sds)
                    ,prior(gamma(0.01,0.01),class=shape)
                    #,prior(exponential(1),class = sigma)
)


brms_id_weibull <- brm(bf(mean_speed ~ s(time_point,k=-1) + 
                            treatment*predator_treatment + 
                            s(time_point,by = treatment,k=-1) +
                            s(time_point,by = predator_treatment,k=-1) +
                            s(time_point,by = treat_inter,k=-1) +
                            (1|replicate)),
                       data = id_data,
                       family = weibull(link = "identity"), #exgaussian/shifted_lognormal/lognormal
                       prior = id_speed_prior,
                       chains = 4, 
                       thin =0.0005*10000,
                       cores = 4, 
                       backend = "cmdstanr", 
                       threads = threading(2),
                       iter = 5000, 
                       warmup = 2000, 
                       refresh = 50,
                       control=list(adapt_delta=0.985,max_treedepth = 20),
                       silent = 0)

saveRDS(brms_id_weibull, file = "Results/brms_id_noar_weibull.rds")
brms_id_weibull <- readRDS(file = "Results/brms_id_noar_weibull.rds")
brms_id_weibullk5 <- readRDS(file = "Results/brms_id_noar_weibullk5.rds")
brms_id_weibullk8 <- readRDS(file = "Results/brms_id_noar_weibullk8.rds")

loo::loo(brms_id_weibull) #21682.5
loo::loo(brms_id_weibullk5) #22208.8
loo::loo(brms_id_weibullk8) #21838.6

pp_check(brms_id_weibull,type = "hist")
pp_check(brms_id_weibullk5,type = "hist")
pp_check(brms_id_weibullk8,type = "hist")

#summary information
loo::loo(brms_id_weibull) #analogous to AIC 21682.0
pp_check(brms_id_weibull) #predicted values vs observed (density plot)
pp_check(brms_id_weibull,type = "hist") #predicted values vs observed (histogram)
summary(brms_id_weibull) #full summary table
bayestestR::describe_posterior(brms_id_weibull, ci = 0.95, test="none") #streamlined summary table (for supplementary)

new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = seq(0,24,by=1)) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat_id <- cbind(new_dat,
                       predict(brms_id_weibull,newdata = new_dat, re_formula =NA))


ggplot(global_dat_id |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  #geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(predator_treatment)),alpha=0.5)+
  #geom_line(aes(x = time_point, y=Estimate,col=as.factor(predator_treatment))) +
  geom_point(data = id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
               mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
             aes(x = time_point, y=mean_speed),col="black", alpha = 0.3, shape = 21)+
  geom_line(aes(x = time_point, y=Estimate),col="black",linewidth=1.5) +
  geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(treatment)),alpha=0.3)+
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  #scale_fill_discrete(guide="none")+
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  #scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Time (hours)")+
  ylab("Mean speed (mm/s)")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))



cond_speed_treatment_id <- brms::conditional_effects(brms_id_weibullk8,
                                                     effects = "time_point",
                                                     #method = "posterior_predict",
                                                     method = "posterior_epred",
                                                     conditions = (data.frame(expand.grid(treatment = c(15,25),
                                                                                          #time_point = seq(from = 0,to=24,by=1),
                                                                                          predator_treatment =  c("prey","didinium","homalozoon"))) |>
                                                                     mutate(treat_inter = interaction(treatment,predator_treatment))))[[1]]

ggplot(cond_speed_treatment_id |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_point(data = id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
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
        panel.border = element_rect(fill = NA, colour = "black"))

###################################################################################
# Gamma speed model
###################################################################################
id_speed_prior <- c(prior(normal(0, 1), class = b),
                    prior(exponential(1.5), class = Intercept, lb = 0),
                    prior(exponential(1), class = sds),
                    prior(exponential(1),class = sd)
                    ,prior(gamma(0.01,0.01),class=shape)
                    #,prior(exponential(1),class = sigma)
)


brms_id_gamma <- brm(bf(mean_speed ~ s(time_point,k=8) + 
                          treatment*predator_treatment + 
                          s(time_point,by = treatment,k=8) +
                          s(time_point,by = predator_treatment,k=8) +
                          s(time_point,by = treat_inter,k=8) +
                          (1|replicate)),
                     data = id_data,
                     family = Gamma(), 
                     prior = id_speed_prior,
                     chains = 4,
                     cores = 4, 
                     backend = "cmdstanr", 
                     threads = threading(2),
                     thin =0.0005*10000,
                     iter = 5000, 
                     warmup = 2000, 
                     refresh = 50,
                     silent = 0
                     ,control=list(adapt_delta=0.99,max_treedepth = 20)
)
saveRDS(brms_id_gamma, file = "Results/brms_id_noar_gammak8.rds")
brms_id_gamma <- readRDS( file = "Results/brms_id_noar_gamma.rds")

#summary information
loo::loo(brms_id_gamma) #analogous to AIC 22190.3
pp_check(brms_id_gamma) #predicted values vs observed (density plot)
pp_check(brms_id_gamma,type = "hist") #predicted values vs observed (histogram)
summary(brms_id_gamma) #full summary table
bayestestR::describe_posterior(brms_id_gamma, ci = 0.95, test="none") #streamlined summary table (for supplementary)

new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = seq(0,24,by=1)) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat_id <- cbind(new_dat,
                       predict(brms_id_gamma,newdata = new_dat, re_formula =NA))


ggplot(global_dat_id |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  #geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(predator_treatment)),alpha=0.5)+
  #geom_line(aes(x = time_point, y=Estimate,col=as.factor(predator_treatment))) +
  geom_point(data = id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
               mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
             aes(x = time_point, y=mean_speed),col="black", alpha = 0.3, shape = 21)+
  geom_line(aes(x = time_point, y=Estimate),col="black",linewidth=1.5) +
  geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(treatment)),alpha=0.3)+
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  #scale_fill_discrete(guide="none")+
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  #scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Time (hours)")+
  ylab("Mean speed (mm/s)")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))



cond_speed_treatment_id <- conditional_effects(brms_id_gamma,
                                               effects = "time_point",
                                               #method = "posterior_predict",
                                               method = "posterior_epred",
                                               conditions = (data.frame(expand.grid(treatment = c(15,25),
                                                                                    #time_point = seq(from = 0,to=24,by=1),
                                                                                    predator_treatment =  c("prey","didinium","homalozoon"))) |>
                                                               mutate(treat_inter = interaction(treatment,predator_treatment))))[[1]]

ggplot(cond_speed_treatment_id |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_point(data = id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
               mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
             aes(x=time_point,y=mean_speed),col="black", alpha = 0.3, shape = 21)+
  geom_line(aes(x = effect1__, y=estimate__),col="black",linewidth=1.5) +
  geom_ribbon(aes(x = time_point,ymin = lower__, ymax =  upper__,fill=as.factor(treatment)),alpha=0.3)+
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))


###################################################################################
# Morphology model
###################################################################################

bprior_lm_id <- c(prior(normal(50,10), class = Intercept, lb=0),
                  prior(normal(0, 1), class = b),
                  prior(exponential(1),class = sd),
                  prior(exponential(1),class = sigma))

brms_lm_id<- brm(bf(mean_width_um ~ mean_length_um*time_point*treatment*predator_treatment +
                      (1|replicate)),
                 data = id_data,
                 family = gaussian(), 
                 prior = bprior_lm_id,
                 chains = 4, 
                 thin =0.0005*10000,
                 cores = 4, 
                 iter = 10000, 
                 warmup = 2000, 
                 silent = 0,
                 control=list(adapt_delta=0.975,max_treedepth = 20))


saveRDS(brms_lm_id, file = "Results/brms_lm_id_noar.rds")
brms_lm_id <- readRDS(file = "Results/brms_lm_id_noar.rds")

cond_inter_treatment_id<- conditional_effects(brms_lm_id,effects = "mean_length_um:treatment",
                                              conditions = data.frame(expand.grid(time_point = seq(0,24,by=8),
                                                                                  predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

ggplot(cond_inter_treatment_id |>
         mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h")))|>
         mutate(effect2__ = paste0(effect2__,"\u00B0C")) |>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_point(data = id_data |>
               ungroup()|>
               dplyr::mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))) |>
               mutate(effect2__ = paste0(treatment,"\u00B0C")) |>
               filter(time_point %in% c(0,8,16,24)) |>
               mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h"))),
             aes(x=mean_length_um,y=mean_width_um,col=as.factor(effect2__)), alpha = 0.4)+
  geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
  geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
  facet_grid(cond__~predator_treatment, scales = "fixed") +
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Mean length (\u03BCm)")+
  ylab("Mean width (\u03BCm)")+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))


cond_inter_treatment_id_pred <- conditional_effects(brms_lm_id,effects = "mean_length_um:predator_treatment",
                                                #method = "posterior_predict",
                                              conditions = data.frame(expand.grid(time_point = seq(0,24,by=8),
                                                                                  treatment =  c(15,25))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")


ggplot(cond_inter_treatment_id_pred |>
         mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h")))|>
         mutate(treatment = paste0(treatment,"\u00B0C")) |>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_point(data = id_data |>
               ungroup()|>
               dplyr::mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))) |>
               mutate(treatment = paste0(treatment,"\u00B0C")) |>
               filter(time_point %in% c(0,8,16,24)) |>
               mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h"))),
             aes(x=mean_length_um,y=mean_width_um,col=as.factor(predator_treatment)), alpha = 0.4)+
  geom_line(aes(x = effect1__, y=estimate__,col=as.factor(predator_treatment))) +
  geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(predator_treatment)),alpha=0.5)+
  facet_grid(cond__~treatment, scales = "fixed") +
  scale_fill_manual(name = "Treatment",values = c("#D9CD8D","#D9B8CF","#5B8B8C")) + 
  scale_color_manual(name = "Treatment",values = c("#D9CD8D","#D9B8CF","#5B8B8C")) + 
  xlab("Mean length (\u03BCm)")+
  ylab("Mean width (\u03BCm)")+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))

