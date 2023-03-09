require(brms)
require(Matrix)
require(tidyverse)
require(bestNormalize)
require(ggh4x)


#load data for analysis
did_id_data <- readRDS(file = "Data/didinium_data_IDs_clean.RDS")
hom_id_data <- readRDS(file = "Data/homalozoon_data_IDs_clean.RDS")
prey_id_data <- readRDS(file = "Data/prey_data_IDs_clean.RDS")




##add predator treatment column to each then combine data frames
#didinium
did_id_data <- did_id_data %>%
  mutate(predator_treatment = "didinium")

#plot to see replicates with strange behaviour 
ggplot(data = did_id_data %>% filter(Species == "PARcau") %>% group_by(time_point, replicate, treatment) %>% 
         summarise(max_abundance = unique(max_abundance)),
       aes(x = time_point, y = max_abundance))+
  geom_point()+
  geom_line()+
  facet_wrap(~treatment*replicate)+
  ggtitle("Didinium")+
  theme_bw()

#homalozoon
hom_id_data <- hom_id_data %>%
  mutate(predator_treatment = "homalozoon")

#plot to see replicates with strange behaviour 
ggplot(data = hom_id_data %>% filter(Species == "PARcau") %>% group_by(time_point, replicate, treatment) %>% 
         summarise(max_abundance = unique(max_abundance)),
       aes(x = time_point, y = max_abundance))+
  geom_point()+
  geom_line()+
  facet_wrap(~treatment*replicate)+
  ggtitle("Homalozoon")+
  theme_bw()

#prey
prey_id_data <- prey_id_data %>%
  mutate(predator_treatment = "prey")

#plot to see replicates with strange behaviour 
ggplot(data = prey_id_data %>% filter(Species == "PARcau") %>% group_by(time_point, replicate, treatment) %>% 
         summarise(max_abundance = unique(max_abundance)),
       aes(x = time_point, y = max_abundance))+
  geom_point()+
  geom_line()+
  facet_wrap(~treatment*replicate)+
  ggtitle("Control")+
  theme_bw()


#combine the data frames

#this dataframe still has the abundances and data for the didiniums and homalozoon present
id_data_with_predators <- rbind(did_id_data, hom_id_data, prey_id_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))%>%
  mutate(
    # mean_speed_norm = predict(bestNormalize::bestNormalize(mean_speed)), we dont need to normalize
         treatment = factor(treatment, levels = c(15,25)),
         treat_inter = as.factor(interaction(treatment,predator_treatment)))

# this dataframe will contain just the parameciums, so equivalent of our previous data
id_data<-id_data_with_predators %>% filter(Species == "PARcau")

##let's plot to see where which replicate have strange data

ggplot(data = id_data %>% group_by(time_point, treat_inter, replicate) %>% summarise(max_abundance = unique(max_abundance)),
       aes(x = time_point, y = max_abundance))+
  geom_point()+
  geom_line()+
  facet_wrap(~treat_inter*replicate)+
  theme_bw()

#let's fo at the replicate level and have average values of all the variables 
grouped_id_data <- id_data %>% 
  group_by(treatment,predator_treatment,replicate,treat_inter,time_point) %>%
  summarise(across(c(max_abundance:mean_speed), ~mean(.x))) 

#from here Duncan's data should work with both "id_data" and "grouped_id_data"


write.csv(grouped_id_data, file = "Data/grouped_id_data.csv")

########################################################################################
## Speed Analysis
########################################################################################
tt <- get_prior(bf(mean_speed ~ s(time_point) + 
                     treatment*predator_treatment + 
                     s(time_point,by = treatment) +
                     s(time_point,by = predator_treatment) +
                     s(time_point,by = treat_inter) + 
                     ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
                     (1|replicate/ID)),
                data = id_data,
                family = gaussian())

bprior2 <- c(prior(normal(0, 1), class = b),
            #prior(lkj(1), class = cor),
            prior(normal(0, 2), class = sds),
            prior(exponential(1),class = sigma))

brms3 <- brm(bf(mean_speed ~ s(time_point) + 
                  treatment*predator_treatment + 
                  s(time_point,by = treatment) +
                  s(time_point,by = predator_treatment) +
                  s(time_point,by = treat_inter) + 
                  ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
                  (1|replicate/ID)),
             data = id_data,
             family = gaussian(), 
             prior = bprior2,
             chains = 4, 
             thin =0.0005*10000,
             cores = 10, 
             iter = 10000, 
             warmup = 2000, 
             silent = 0,
             control=list(adapt_delta=0.975,max_treedepth = 15))

saveRDS(brms3, file = "Results/brms3.rds")
brms3<-readRDS("Results/brms3.rds")

summary(brms3)
#bayestestR::describe_posterior(brms3, ci = 0.95, test="none")
pp_check(brms3,type = "dens_overlay")
pp_check(brms3,type = "loo_pit_qq")

brms::mcmc_plot(brms3, 
                #type = "areas",
                type = "intervals",
                prob = 0.95)

new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = seq(0,24,by=1), replicate = NA, ID = NA) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat2 <- cbind(new_dat,
                    predict(brms3,newdata = new_dat, re_formula =NA))

ggplot(global_dat2 |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  #geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(predator_treatment)),alpha=0.5)+
  #geom_line(aes(x = time_point, y=Estimate,col=as.factor(predator_treatment))) +
  geom_line(aes(x = time_point, y=Estimate),col="black") +
  geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(treatment)),alpha=0.3)+
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  #scale_fill_discrete(guide="none")+
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  #scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Time (hours)")+
  ylab("Mean speed  (\u03BCm/s)")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))



conditional_effects(brms3,effects = "time_point")
conditional_effects(brms3,effects = "time_point:treatment")
conditional_effects(brms3,effects = "time_point:predator_treatment")
conditional_effects(brms3,effects = "time_point:treat_inter")

########################################################################################
## Speed Analysis Null
########################################################################################

brms_null <- brm(bf(mean_speed ~ s(time_point) + 
                  treatment*predator_treatment + 
                  s(time_point,by = treatment) +
                  s(time_point,by = predator_treatment) +
                  s(time_point,by = treat_inter)),
             data = id_data,
             family = gaussian(), 
             prior = c(prior(normal(0, 1), class = b)),
             chains = 4, 
             thin =0.0005*4000,
             cores = 8, 
             iter = 4000, 
             warmup = 2000, 
             silent = 0,
             control=list(adapt_delta=0.975,max_treedepth = 20))

saveRDS(brms_null, file = "Results/brms_null.rds")
brms_null <-readRDS("Results/brms_null.rds")
summary(brms_null)
#bayestestR::describe_posterior(brms_null, ci = 0.95, test="none")

brms::mcmc_plot(brms_null, 
                #type = "areas",
                type = "intervals",
                prob = 0.95)
new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = seq(0,24,by=0.1), replicate = NA, ID = NA) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat_null <- cbind(new_dat,
                     predict(brms_null,newdata = new_dat, re_formula =NA))

ggplot(global_dat_null |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_line(aes(x = time_point, y=Estimate),col="black") +
  geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(treatment)),alpha=0.3)+
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  #scale_fill_discrete(guide="none")+
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  #scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Time (hours)")+
  ylab("Mean speed  (\u03BCm/s)")+
  labs(colour = "Predator treatment") +
  theme_bw()+
  # theme_classic()+
  # theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
  #       strip.text.y.right = element_text(angle = 0),
  #       strip.text.x = element_text(face = "bold.italic"),
  #       strip.text = element_text(face = "bold"),
  #       panel.background = element_blank(),
  #       axis.line = element_line(colour = "black"),
  #       panel.border = element_rect(fill = NA, colour = "black"))
  



########################################################################################
## Morphology Analysis
########################################################################################

#potential linear model bayesian for length ~ width
#priors
get_prior(bf(mean_width_um ~ mean_length_um*time_point*treatment*predator_treatment+
               ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
               (1|replicate/ID)),
          data = id_data,
          family = gaussian())

bprior_lm <- c(prior("",class = ar , ub = 1, lb = -1),
               prior(normal(0, 1), class = b),
               prior(exponential(1),class = sd),
               prior(exponential(1),class = sigma))


brms_lm <- brm(bf(mean_width_um ~ mean_length_um*time_point*treatment*predator_treatment+
         ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
         (1|replicate/ID)),
    data = id_data,
    family = gaussian(), 
    prior = bprior_lm,
    chains = 4, 
    thin =0.0005*10000,
    cores = 4, 
    iter = 10000, 
    warmup = 2000, 
    silent = 0,
    control=list(adapt_delta=0.975,max_treedepth = 10))

saveRDS(brms_lm, file = "Results/brms_lm.rds")
brms_lm <- readRDS("Results/brms_lm.rds")

pp_check(brms_lm,type = "dens_overlay")
pp_check(brms_lm,type = "loo_pit_qq")

brms::mcmc_plot(brms_lm, 
                #type = "areas",
                type = "intervals",
                prob = 0.95)

cond_inter_treatment2 <- conditional_effects(brms_lm,effects = "mean_length_um:treatment",
                                             conditions = data.frame(expand.grid(time_point = unique(id_data$time_point),
                                                                                 predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

ggplot(cond_inter_treatment2 |>
         mutate(cond__ = factor(paste0("t",time_point),levels = paste0("t",unique(time_point))))|>
         mutate(effect2__ = paste0(effect2__,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
  geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
  facet_grid(cond__~predator_treatment, scales = "fixed") +
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  scale_color_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Mean length (\u03BCm)")+
  ylab("Mean width (\u03BCm)")+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))


cond_inter_treatment_sub <- conditional_effects(brms_lm,effects = "mean_length_um:treatment",
                                             conditions = data.frame(expand.grid(time_point = seq(0,24,by=8),
                                                                                 predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

ggplot(cond_inter_treatment_sub |>
         mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h")))|>
         mutate(effect2__ = paste0(effect2__,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
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

########################################################################################
## Average across IDs
########################################################################################
group_prior <- c(prior(normal(0, 1), class = b),
                 prior(exponential(1.5), class = Intercept, lb = 0),
             prior(normal(0,0.5), class = ar), 
             prior(exponential(1), class = sds),
             prior(exponential(1),class = sigma))


brms_group <- brm(bf(mean_speed ~ s(time_point) + 
                  treatment*predator_treatment + 
                  s(time_point,by = treatment) +
                  s(time_point,by = predator_treatment) +
                  s(time_point,by = treat_inter) + 
                  ar(time = time_point,gr = replicate:treat_inter,p=1)+
                  (1|replicate)),
             data = grouped_id_data,
             family = gaussian(), 
             prior = group_prior,
             chains = 4, 
             thin =0.0005*10000,
             cores = 4, 
             iter = 10000, 
             warmup = 2000, 
             backend = "cmdstanr", 
             #threads = threading(2),
             silent = 0,
             control=list(adapt_delta=0.975,max_treedepth = 20))

#quick plot of the raw data
ggplot(data = grouped_id_data, aes(x = time_point, y = mean_speed, group = replicate))+
  geom_line()+
  facet_nested_wrap(~treatment*predator_treatment)+
  theme_bw()


#save output of the model
saveRDS(brms_group, file = "Results/brms_group.rds")
coef(brms_group)

#load output of the model
brms_group<-readRDS("Results/brms_group.rds")

pp_check(brms_group,type = "dens_overlay")
pp_check(brms_group,type = "hist")

summary(brms_group)
prior_summary(brms_group)

new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = seq(0,24,by=1), replicate = NA, ID = NA) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat_grouped <- cbind(new_dat,
                     predict(brms_group,newdata = new_dat, re_formula =NA))


ggplot(global_dat_grouped |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  #geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(predator_treatment)),alpha=0.5)+
  #geom_line(aes(x = time_point, y=Estimate,col=as.factor(predator_treatment))) +
  geom_point(data = grouped_id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
               mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
             aes(x = time_point, y=mean_speed),col="black", alpha = 0.3, shape = 21)+
  geom_line(data = grouped_id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
               mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
             aes(x = time_point, y=mean_speed,group = replicate),col="black", alpha = 0.3)+
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



cond_speed_treatment_grouped <- brms::conditional_effects(brms_group,
                                                          effects = "time_point",
                                                          #method = "posterior_predict",
                                                          method = "posterior_epred",
                                                          conditions = (data.frame(expand.grid(treatment = c(15,25),
                                                                                               #time_point = seq(from = 0,to=24,by=1),
                                                                                               predator_treatment =  c("prey","didinium","homalozoon"))) |>
                                                                          mutate(treat_inter = interaction(treatment,predator_treatment))))[[1]]

ggplot(cond_speed_treatment_grouped |>
         mutate(treatment = paste0(treatment,"\u00B0C"))|>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_point(data = grouped_id_data |>
               ungroup()|>
               mutate(treatment = paste0(treatment,"\u00B0C"))|>
               mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
             aes(x=time_point,y=mean_speed),col="black", alpha = 0.3, shape = 21)+
  geom_line(data = grouped_id_data |>
              ungroup()|>
              mutate(treatment = paste0(treatment,"\u00B0C"))|>
              mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))),
            aes(x = time_point, y=mean_speed,group = replicate),col="black", alpha = 0.3)+
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




####Morphology analysis with mean values across IDs##################

bprior_lm_grouped <- c(prior(normal(0,0.5), class = ar), 
                       prior(normal(50,10), class = Intercept, lb=0),
               prior(normal(0, 1), class = b),
               prior(exponential(1),class = sd),
               prior(exponential(1),class = sigma))

bprior_lm_grouped_no_ar <- c(prior(normal(50,10), class = Intercept, lb=0),
                       prior(normal(0, 1), class = b),
                       prior(exponential(1),class = sd),
                       prior(exponential(1),class = sigma))

brms_lm_grouped_prior_no_ar <- brm(bf(mean_width_um ~ mean_length_um*time_point*treatment*predator_treatment+
                                        (1|replicate)),
                       data = grouped_id_data,
                       family = gaussian(), 
                       prior = bprior_lm_grouped_no_ar,
                       chains = 4, 
                       thin =0.0005*10000,
                       cores = 4, 
                       iter = 2000, 
                       warmup = 1000, 
                       # silent = 0,
                       control=list(adapt_delta=0.975,max_treedepth = 20))

pp_check(brms_lm_grouped_prior_chk,type = "bars")

brms_lm_grouped <- brm(bf(mean_width_um ~ mean_length_um*time_point*treatment*predator_treatment+
                    ar(time = time_point,gr = replicate:treat_inter,p=1)+
                    (1|replicate)),
               data = grouped_id_data,
               family = gaussian(), 
               prior = bprior_lm_grouped,
               chains = 4, 
               thin =0.0005*10000,
               cores = 4, 
               iter = 10000, 
               warmup = 2000, 
               backend = "cmdstanr", 
               silent = 0,
               control=list(adapt_delta=0.975,max_treedepth = 20))

saveRDS(brms_lm_grouped, file = "Results/brms_lm_grouped.rds")
brms_lm_grouped <- readRDS("Results/brms_lm_grouped.rds")

pp_check(brms_lm_grouped,type = "dens_overlay")

summary(brms_lm_grouped)
bayestestR::describe_posterior(brms_lm_grouped, ci = 0.95, test="none")

cond_inter_treatment_grouped <- conditional_effects(brms_lm_grouped,effects = "mean_length_um:treatment",
                                                conditions = data.frame(expand.grid(time_point = seq(0,24,by=8),
                                                                                    predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

ggplot(cond_inter_treatment_grouped |>
         mutate(cond__ = factor(paste0(time_point,"h"),levels = c("0h","8h","16h","24h")))|>
         mutate(effect2__ = paste0(effect2__,"\u00B0C")) |>
         mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
  geom_point(data = grouped_id_data |>
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
