require(brms)
require(Matrix)
require(tidyverse)
require(bestNormalize)

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

tt <- get_prior(bf(mean_speed ~ s(time_point) + 
                     treatment*predator_treatment + 
                     s(time_point,by = treatment) +
                     s(time_point,by = predator_treatment) +
                     s(time_point,by = treat_inter) + 
                     ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
                     (1|replicate/ID)),
                data = id_data,
                family = gaussian())

bprior <- c(prior(normal(0, 1), class = b),
            prior(lkj(1), class = cor),
            prior(normal(0, 1), class = sds),
            prior(exponential(1),class = sd),
            prior(exponential(1),class = sigma)
)

brms1 <- brm(bf(mean_speed ~ s(time_point) + 
                  treatment*predator_treatment + 
                  #s(time_point,by = treatment) +
                  #s(time_point,by = predator_treatment) +
                  s(time_point,by = treat_inter) + 
                  #ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
                  (time_point|replicate/ID)),
             data = id_data,
             family = gaussian(), 
             prior = bprior,
             chains = 4, 
             thin =0.0005*4000,
             cores = 4, 
             iter = 4000, 
             warmup = 2000, 
             silent = 0,
             control=list(adapt_delta=0.99,max_treedepth = 11),
             backend = "cmdstanr"
)

brms::mcmc_plot(brms1, 
                type = "areas",
                #type = "intervals",
                prob = 0.95)

pp_check(brms1,type = "loo_pit_qq")
new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = 0:24) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

tt <- predict(brms1,re_formula = ~ (time_point|replicate/ID))
pred_dat <- cbind(id_data,rep_fit = predict(brms1,re_formula = ~ (time_point|replicate/ID))[,1],
                  global_fit = predict(brms1,re_formula =NA)[,1],
                  lwr = predict(brms1,re_formula =NA)[,3],upr = predict(brms1,re_formula =NA)[,4])


pred_dat2 <- cbind(id_data,rep_fit = predict(brms1,re_formula = ~ (time_point|replicate/ID))[,1])

global_dat <- cbind(new_dat,
                    global_fit = predict(brms1,newdata = new_dat, re_formula =NA)[,1],
                    lwr = predict(brms1,newdata = new_dat,re_formula =NA)[,3],
                    upr = predict(brms1,newdata = new_dat,re_formula =NA)[,4])

ggplot(pred_dat)+
  #geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point,y=rep_fit,col=as.factor(predator_treatment),group = interaction(replicate,ID)),alpha = 0.3) +
  geom_line(data = global_dat,aes(x = time_point, y=global_fit,col=as.factor(predator_treatment))) +
  geom_ribbon(data = global_dat,aes(x = time_point,ymin = lwr, ymax =  upr,fill=as.factor(predator_treatment)),alpha=0.5)+
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


brms2 <- brm(bf(mean_speed ~ s(time_point) + 
                  treatment*predator_treatment + 
                  s(time_point,by = treatment) +
                  s(time_point,by = predator_treatment) +
                  s(time_point,by = treat_inter) + 
                  ar(time = time_point,gr = replicate:treat_inter:ID,p=1)+
                  (time_point|replicate/ID)),
             data = id_data,
             family = gaussian(), 
             prior = bprior,
             chains = 4, 
             thin =0.0005*4000,
             cores = 4, 
             iter = 4000, 
             warmup = 2000, 
             silent = 0,
             control=list(adapt_delta=0.99,max_treedepth = 11),
             backend = "cmdstanr"
)



bprior <- c(prior(normal(0, 1), class = b),
            prior(lkj(1), class = cor),
            prior(normal(0, 1), class = sds),
            prior(exponential(1),class = sd),
            prior(exponential(1),class = sigma)
)

bprior2 <- c(prior(normal(0, 1), class = b),
            #prior(lkj(1), class = cor),
            prior(normal(0, 2), class = sds))

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
             thin =0.0005*4000,
             cores = 8, 
             iter = 4000, 
             warmup = 2000, 
             silent = 0,
             control=list(adapt_delta=0.975,max_treedepth = 20))

saveRDS(brms3, file = "Results/brms3.rds")

summary(brms3)
#bayestestR::describe_posterior(brms3, ci = 0.95, test="none")

brms::mcmc_plot(brms3, 
                #type = "areas",
                type = "intervals",
                prob = 0.95)
pred_dat2 <- cbind(id_data,rep_fit = predict(brms3,re_formula = ~ (time_point|replicate/ID))[,1],
                  global_fit = predict(brms3,re_formula =NA)[,1],
                  lwr = predict(brms3,re_formula =NA)[,3],upr = predict(brms3,re_formula =NA)[,4])

new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = 0:24, replicate = NA, ID = NA) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat2 <- cbind(new_dat,
                    predict(brms3,newdata = new_dat, re_formula =NA))

ggplot(global_dat2)+
  geom_line(aes(x = time_point, y=Estimate,col=as.factor(predator_treatment))) +
  geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(predator_treatment)),alpha=0.5)+
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


conditional_effects(brms3,effects = "time_point")
conditional_effects(brms3,effects = "time_point:treatment")
conditional_effects(brms3,effects = "time_point:predator_treatment")
conditional_effects(brms3,effects = "time_point:treat_inter")




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

ggplot(global_dat_null)+
  geom_line(aes(x = time_point, y=Estimate,col=as.factor(predator_treatment))) +
  geom_ribbon(aes(x = time_point,ymin = Q2.5, ymax =  Q97.5,fill=as.factor(predator_treatment)),alpha=0.5)+
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
