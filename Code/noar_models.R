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

id_data<-id_data_with_predators %>% filter(Species == "PARcau")


id_speed_prior <- c(prior(normal(0, 1), class = b),
                 prior(exponential(1.5), class = Intercept, lb = 0),
                 prior(exponential(1), class = sds),
                 prior(exponential(1),class = sigma))


brms_id <- brm(bf(mean_speed ~ s(time_point) + 
                       treatment*predator_treatment + 
                       s(time_point,by = treatment) +
                       s(time_point,by = predator_treatment) +
                       s(time_point,by = treat_inter) +
                    (1|replicate)),
                  data = id_data,
                  family = gaussian(), 
                  prior = id_speed_prior,
                  chains = 4, 
                  thin =0.0005*10000,
                  cores = 4, 
                  iter = 2000, 
                  warmup = 1000, 
                  silent = 0,
                  control=list(adapt_delta=0.975,max_treedepth = 20))

new_dat <- expand.grid(treatment = c(15,25),
                       predator_treatment = c("prey","didinium","homalozoon"),
                       time_point = seq(0,24,by=1)) |>
  mutate(treat_inter = interaction(treatment,predator_treatment))

global_dat_id <- cbind(new_dat,
                            predict(brms_id,newdata = new_dat, re_formula =NA))


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
                       iter = 2000, 
                       warmup = 1000, 
                       silent = 0,
                       control=list(adapt_delta=0.975,max_treedepth = 20))

                 
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
