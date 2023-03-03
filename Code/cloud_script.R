require(brms)
require(Matrix)
require(tidyverse)

did_id_data <- readRDS(file = "/home/opc/didinium_data_IDs_clean.RDS")
hom_id_data <- readRDS(file = "/home/opc/homalozoon_data_IDs_clean.RDS")
prey_id_data <- readRDS(file = "/home/opc/prey_data_IDs_clean.RDS")

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

id_speed_prior <- c(prior(normal(0, 1), class = b),
                    prior(exponential(1), class = Intercept, lb = 0),
                    prior(exponential(1), class = sds),
                    prior(exponential(1),class = sigma))


brms_id <- brm(bf(mean_speed ~ s(time_point,k=5) + 
                    treatment*predator_treatment + 
                    s(time_point,by = treatment,k=5) +
                    s(time_point,by = predator_treatment,k=5) +
                    s(time_point,by = treat_inter,k=5) +
                    (1|replicate)),
               data = id_data,
               family = gaussian(), #exgaussian/shifted_lognormal/lognormal
               #prior = id_speed_prior,
               chains = 4, 
               thin =0.0005*10000,
               cores = 4, 
               iter = 100, 
               #warmup = 100, 
               refresh = 100,
               control=list(adapt_delta=0.975,max_treedepth = 20),
               silent = 0)

saveRDS(brms_id, file = "/home/opc/brms_id_noar.rds")
