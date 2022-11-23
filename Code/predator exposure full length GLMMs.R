###########################################################
######### full length experiment data analysis #########
###########################################################
rm(list=ls())

#load packages
library(tidyverse)
library(bestNormalize) #transforming response variables to fit gaussian model
library(MuMIn) #model comparison and simplification
library(rstatix) 
library(glmmTMB) #GLMMs
library(DHARMa) #GLMM diagnostics
library(fitdistrplus) #data to distribution fit

#load in the tracked ID data 

did_id_data <- read.csv("Data/did_data_IDs_corrected_max_abundance.csv")
hom_id_data <- read.csv("Data/hom_data_IDs_corrected_max_abundance.csv")
prey_id_data <- read.csv("Data/prey_data_IDs_correct_max_abundance.csv")

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
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))

####################################
###### analyse mean speed ##########
####################################

#plot raw data
speed_id_plot <-
  ggplot(id_data, aes(x = time_point, y = mean_speed, col = as.factor(predator_treatment)))+
  geom_smooth(se = FALSE) +
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

#fit speed gaussian GLM with three-way interaction

speed_mod_gaussian <- glm(mean_speed ~ time_point + treatment + predator_treatment +
                               time_point:treatment + time_point:predator_treatment + treatment:predator_treatment +
                               time_point:treatment:predator_treatment,
                             data = id_data, family = "gaussian") 
#model output
summary(speed_mod_gaussian)

#visual model checking
par(mfrow = c(2,2)) #for diagnostic plots
plot(speed_mod_gaussian) #scale-location not straight, qq has heavy tails

#account for non-linearity of time by fitting a spline
#fit splines to account for non-linearity through time
#quadratic spline
speed_mod_gaussian_spline <- glm(mean_speed ~ splines::ns(time_point,2) + treatment + predator_treatment +
                                      splines::ns(time_point,2):treatment + splines::ns(time_point,2):predator_treatment + treatment:predator_treatment +
                                      splines::ns(time_point,2):treatment:predator_treatment,
                                    data = id_data, family = "gaussian")

plot(speed_mod_gaussian_spline) #same problems as above

#fit a cubic spline
speed_mod_gaussian_cbspline <- glm(mean_speed ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                        splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                        splines::ns(time_point,3):treatment:predator_treatment,
                                      data = id_data, family = "gaussian")

plot(speed_mod_gaussian_cbspline)

MuMIn::model.sel(speed_mod_gaussian, speed_mod_gaussian_spline, speed_mod_gaussian_cbspline)
#cubic spline is best fit

#normalise data using bestNormalize and run on cubic spline model
id_data$mean_speed_norm <- predict(bestNormalize(id_data$mean_speed))


speed_norm_mod_gaussian_cbspline <- glm(mean_speed_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                             splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                             splines::ns(time_point,3):treatment:predator_treatment,
                                           data = id_data, family = "gaussian")
plot(speed_norm_mod_gaussian_cbspline) #diagnostics look better, few outliers (4931, 806, 2254) that we could remove

#account for autocorrelation and non-independence of replicates
#by fitting a GLMM
#rename best mod above for simplification
speed_mod_1 <- glm(mean_speed_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                          splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                          splines::ns(time_point,3):treatment:predator_treatment,
                                        data = id_data, family = "gaussian")

#fit GLMM
speed_mod_auto <- glmmTMB::glmmTMB(mean_speed_norm ~ splines::ns(time_point,3) +
                                     treatment + 
                                     predator_treatment + 
                                     splines::ns(time_point,3):predator_treatment +
                                     splines::ns(time_point,3):treatment  +
                                     predator_treatment:treatment +
                                     splines::ns(time_point,3):treatment:predator_treatment +
                                     #bit of code that deals with the autocorrelation
                                     #uncorrelated random intercept and correlated slope within group
                                     ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                   data = id_data, family = "gaussian", REML=F)
summary(speed_mod_auto)
Anova(speed_mod_auto)

#compare data to the distribution fit to it
fit_gaus <- fitdist(id_data$mean_speed_norm,
                    distr = "norm")
plot(fit_gaus) #looks like a good fit

#pseudo r squared
MuMIn::r.squaredGLMM(speed_mod_auto)

#simulate residuals from the model
speed_glmm_sim <- simulateResiduals(speed_mod_auto, n = 1000)
#plot the residuals
plot(speed_glmm_sim)

testOutliers(speed_glmm_sim, plot = TRUE)

#compare AIC vals
MuMIn::model.sel(speed_mod_auto, speed_mod_1) #speed_mod_auto much better fit

#manual model simplification on the new speed mod - start by removing 3-way interaction
speed_mod_auto_2 <- glmmTMB::glmmTMB(mean_speed_norm ~ splines::ns(time_point,3) +
                                       treatment + 
                                       predator_treatment + 
                                       splines::ns(time_point,3):predator_treatment +
                                       splines::ns(time_point,3):treatment  +
                                       predator_treatment:treatment +
                                       #bit of code that deals with the autocorrelation
                                       ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                     data = id_data, family = "gaussian", REML=F)

MuMIn::model.sel(speed_mod_auto, speed_mod_auto_2) #AIC lower with the 3-way interaction so we keep it

#fit model predictions without the random effect
#new dataframe of unique variables required for plotting
plotting_data_s <- unique(id_data[,c('treatment', 'predator_treatment', 'replicate', 'time_point')])

#set the random effects of the model to zero
X_cond <- model.matrix(lme4::nobars(formula(speed_mod_auto)[-2]), plotting_data_s)
beta_cond <- fixef(speed_mod_auto)$cond
pred_cond <- X_cond %*% beta_cond

#transform point estimates of the speed to the response scale and multiply
ilink <- family(speed_mod_auto)$linkinv
pred_speed = ilink(pred_cond)

#load the MASS library
library(MASS)
#set random number generator seed
set.seed(101)

#use posterior predictive simulations to generate upper and lower CIs
#and median predicted speeds
#ignoring variation in the random effects
#conditional
pred_condpar_psim = mvrnorm(1000, mu = beta_cond,Sigma = vcov(speed_mod_auto)$cond)
pred_cond_psim = X_cond %*% t(pred_condpar_psim)

#transform
pred_speed_psim = ilink(pred_cond_psim)

#calculate 95% CIs/
ci_speed = t(apply(pred_speed_psim,1, quantile, c(0.025, 0.975)))
ci_speed = data.frame(ci_speed)

#rename
names(ci_speed) = c("speed_low", "speed_high")

#put into a data frame 
plotting_data_s = data.frame(plotting_data_s, pred_speed, ci_speed)

#plot on the graph 
ggplot(plotting_data_s, aes(x = time_point, y = pred_speed, col = as.factor(predator_treatment)))+
  #      geom_point(data = id_data, aes(x = time_point, y = mean_speed_norm), alpha = 0.1, pch = 1) +
  # geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point, y = pred_speed)) +
  geom_ribbon(aes(x = time_point,
                  ymin = speed_low,
                  ymax = speed_high,
                  group = predator_treatment), alpha = 0.1, linetype = 0) +
  facet_grid(treatment~., scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean speed")+
  labs(colour = "Predator treatment") +
  theme_bw()
  # theme_classic()+
  # theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
  #       panel.background = element_blank(),
  #       axis.line = element_line(colour = "black"),
  #       aspect.ratio = 1,
  #       panel.border = element_rect(fill = NA, colour = "black"))


#######################################
###### analyse mean width ############
#####################################

#fit width gaussian glm with three-way interaction
width_mod_gaussian <- glm(mean_width_um ~ time_point + treatment + predator_treatment +
                               time_point:treatment + time_point:predator_treatment + treatment:predator_treatment +
                               time_point:treatment:predator_treatment,
                             data = id_data, family = "gaussian") 
#model output
summary(width_mod_gaussian)

#visual model checking
par(mfrow = c(2,2)) #for diagnostic plots
plot(width_mod_gaussian) #qq looks a bit skewed

#account for non-linearity of time by fitting a spline
#fit splines to account for non-linearity through time
width_mod_gaussian_spline <- glm(mean_width_um ~ splines::ns(time_point,2) + treatment + predator_treatment +
                                      splines::ns(time_point,2):treatment + splines::ns(time_point,2):predator_treatment + treatment:predator_treatment +
                                      splines::ns(time_point,2):treatment:predator_treatment,
                                    data = id_data, family = "gaussian")

plot(width_mod_gaussian_spline) 

width_mod_gaussian_cbspline <- glm(mean_width_um ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                        splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                        splines::ns(time_point,3):treatment:predator_treatment,
                                      data = id_data, family = "gaussian")

plot(width_mod_gaussian_cbspline)

MuMIn::model.sel(width_mod_gaussian, width_mod_gaussian_spline, width_mod_gaussian_cbspline) #cubic spline best fit

#try to normalise data using best normalize and run on cubic spline model
id_data$mean_width_norm <- predict(bestNormalize(id_data$mean_width))

width_norm_mod_gaussian_cbspline <- glm(mean_width_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                             splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                             splines::ns(time_point,3):treatment:predator_treatment,
                                           data = id_data, family = "gaussian")

plot(width_norm_mod_gaussian_cbspline)

#model simplification by removing 3-way then 2-way interactions
width_mod_2 <- glm(mean_width_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                        splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment,
                      data = id_data, family = "gaussian")
plot(width_mod_2)

MuMIn::model.sel(width_norm_mod_gaussian_cbspline, width_mod_2) #model 2 has higher AIC so we keep model 1

#mixed effects model to account for non-independence and temporal autocorrelation 
#rename best mod from above
width_mod_1 <- glm(mean_width_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                          splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                          splines::ns(time_point,3):treatment:predator_treatment,
                                        data = id_data, family = "gaussian")

width_mod_auto <- glmmTMB::glmmTMB(mean_width_norm ~ splines::ns(time_point,3) +
                                     treatment + 
                                     predator_treatment + 
                                     splines::ns(time_point,3):predator_treatment +
                                     splines::ns(time_point,3):treatment  +
                                     predator_treatment:treatment +
                                     splines::ns(time_point,3):treatment:predator_treatment +
                                     #bit of code that deals with the autocorrelation
                                     #uncorrelated random intercept and correlated slope within group
                                     ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                   data = id_data, family = "gaussian", REML=F)
summary(width_mod_auto)
Anova(width_mod_auto)

#compare data to the distribution fit to it
fit_gaus <- fitdist(id_data$mean_width_norm,
                    distr = "norm")
plot(fit_gaus) #qq still heavy tails

#model diagnostics can't be done the easy way so have to do them manually 
plot(fitted(width_mod_auto), residuals(width_mod_auto))

#psuedo R2 glmm
MuMIn::r.squaredGLMM(width_mod_auto)

#simulate residuals from the model
width_glmm_sim <- simulateResiduals(width_mod_auto, n = 1000)
#plot the residuals
plot(width_glmm_sim)

testOutliers(width_glmm_sim, plot = TRUE)

#manual model simplification on the new speed mod - start by removing 3-way interaction
width_mod_auto_2 <- glmmTMB::glmmTMB(mean_width_norm ~ splines::ns(time_point,3) +
                                       treatment + 
                                       predator_treatment + 
                                       splines::ns(time_point,3):predator_treatment +
                                       splines::ns(time_point,3):treatment  +
                                       predator_treatment:treatment +
                                       #bit of code that deals with the autocorrelation
                                       ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                     data = id_data, family = "gaussian", REML=F)

MuMIn::model.sel(width_mod_auto, width_mod_auto_2) #AIC lower with the 3-way interaction so we keep it

#extract model predictions for plotting
#newdata dataframe of unique variables required for plotting
plotting_data_w <- unique(id_data[,c('treatment', 'predator_treatment', 'replicate', 'time_point')])

#set the random effects of the model to zero
X_cond <- model.matrix(lme4::nobars(formula(width_mod_auto)[-2]), plotting_data_w)
beta_cond <- fixef(width_mod_auto)$cond
#predictions on the link scale generated from the dataframe with random effects set to zero
pred_cond <- X_cond %*% beta_cond

#transform point estimates of the speed to the response scale and multiply
ilink <- family(width_mod_auto)$linkinv
pred_width = ilink(pred_cond)

#set random number generator seed
set.seed(101)

#use posterior predictive simulations to generate upper and lower CIs
#and median predicted speeds
#ignoring variation in the random effects
#conditional
pred_condpar_psim = mvrnorm(1000, mu = beta_cond,Sigma = vcov(width_mod_auto)$cond)
pred_cond_psim = X_cond %*% t(pred_condpar_psim)

#transform
pred_width_psim = ilink(pred_cond_psim)

#calculate 95% CIs/
ci_width = t(apply(pred_width_psim,1, quantile, c(0.025, 0.975)))
ci_width = data.frame(ci_width)

#rename
names(ci_width) = c("width_low", "width_high")

#put into a data frame 
plotting_data_w = data.frame(plotting_data_w, pred_width, ci_width)

#plot on the graph 
ggplot(plotting_data_w, aes(x = time_point, y = pred_width, col = as.factor(predator_treatment)))+
  # geom_point(data = id_data, aes(x = time_point, y = mean_roundness_norm), alpha = 0.1, pch = 1) +
  # geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point, y = pred_width)) +
  geom_ribbon(aes(x = time_point,
                  ymin = width_low,
                  ymax = width_high,
                  group = predator_treatment), alpha = 0.1, linetype = 0) +
  facet_grid(treatment~., scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean width")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

#########################################
######### length ID GLM ################
########################################

#fit length gaussian glm with three-way interaction
length_mod_gaussian <- glm(max_length ~ time_point + treatment + predator_treatment +
                                time_point:treatment + time_point:predator_treatment + treatment:predator_treatment +
                                time_point:treatment:predator_treatment,
                              data = id_data, family = "gaussian") 
#model output
summary(length_mod_gaussian)

#visual model checking
par(mfrow = c(2,2)) #for diagnostic plots
plot(length_mod_gaussian) #scale-location plot not straight line

#account for non-linearity of time by fitting a spline
#fit splines to account for non-linearity through time
length_mod_gaussian_spline <- glm(max_length ~ splines::ns(time_point,2) + treatment + predator_treatment +
                                       splines::ns(time_point,2):treatment + splines::ns(time_point,2):predator_treatment + treatment:predator_treatment +
                                       splines::ns(time_point,2):treatment:predator_treatment,
                                     data = id_data, family = "gaussian")

plot(length_mod_gaussian_spline) #scale_location plot not straight

length_mod_gaussian_cbspline <- glm(max_length ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                         splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                         splines::ns(time_point,3):treatment:predator_treatment,
                                       data = id_data, family = "gaussian")

plot(length_mod_gaussian_cbspline) #scale-location plot worse

MuMIn::model.sel(length_mod_gaussian, length_mod_gaussian_spline, length_mod_gaussian_cbspline) #cubic spline best fit although diagnostics are worse

#normalise data using best normalize and run on cubic spline model
id_data$max_length_norm <- predict(bestNormalize(id_data$max_length))

length_mod_1 <- glm(max_length_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                      splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                      splines::ns(time_point,3):treatment:predator_treatment,
                    data = id_data, family = "gaussian")

plot(length_mod_1) #looks good except 4898

#model simplification by removing 3-way then 2-way interactions
length_mod_2 <- glm(max_length_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                      splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment,
                    data = id_data, family = "gaussian")
plot(length_mod_2)
MuMIn::model.sel(length_mod_1, length_mod_2) #model 2 lower AIC so we continue simplifying

length_mod_3 <- glm(max_length_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                      # splines::ns(time_point,3):treatment + 
                      splines::ns(time_point,3):predator_treatment + 
                      treatment:predator_treatment,
                    data = id_data, family = "gaussian")

length_mod_4 <- glm(max_length_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                      splines::ns(time_point,3):treatment + 
                      #splines::ns(time_point,3):predator_treatment + 
                      treatment:predator_treatment,
                    data = id_data, family = "gaussian")

length_mod_5 <- glm(max_length_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                      splines::ns(time_point,3):treatment + 
                      splines::ns(time_point,3):predator_treatment,
                    #treatment:predator_treatment,
                    data = id_data, family = "gaussian")

AIC(length_mod_1, length_mod_2, length_mod_3, length_mod_4, length_mod_5) #mod 2 still the lowest 

#extract model preds and plot 
id_data$length_preds <- predict(length_mod_2)
id_data$length_se <- 1.96*(predict(length_mod_2, se.fit = T)$se.fit)

ggplot(id_data, aes(x = time_point, y = max_length_norm, col = as.factor(predator_treatment)))+
  #    geom_point(alpha = 0.1) +
  # geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point, y = length_preds)) +
  geom_ribbon(aes(x = time_point,
                  ymin = length_preds - length_se,
                  ymax = length_preds + length_se,
                  group = predator_treatment), alpha = 0.2, linetype = 0) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean length")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

##### length mixed effects model

#model with extra bits, spaced out for clarity
length_mod_auto <- glmmTMB::glmmTMB(max_length_norm ~ splines::ns(time_point,3) +
                                      treatment + 
                                      predator_treatment + 
                                      splines::ns(time_point,3):predator_treatment +
                                      splines::ns(time_point,3):treatment  +
                                      predator_treatment:treatment +
                                      #bit of code that deals with the autocorrelation
                                      #uncorrelated random intercept and correlated slope within group
                                      ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                    data = id_data, family = "gaussian", REML=F)
summary(length_mod_auto)
Anova(length_mod_auto)

#compare fit with other model
MuMIn::model.sel(length_mod_auto, length_mod_1) #auto much lower AIC so continue to check fit 

#compare data to the distribution fit to it
fit_gaus <- fitdist(id_data$max_length_norm,
                    distr = "norm")
plot(fit_gaus) #bottom of qq doesn't look right

#model diagnostics can't be done the easy way so have to do them manually 
plot(fitted(length_mod_auto), residuals(length_mod_auto))

#compare the residuals of both models to check for autocorrelation
acf(resid(length_mod_auto))
acf(resid(length_mod_1))

#pseudo r2
MuMIn::r.squaredGLMM(length_mod_auto)

#simulate residuals from the model
length_glmm_sim <- simulateResiduals(length_mod_auto, n = 1000)
#plot the residuals
plot(length_glmm_sim)

testOutliers(length_glmm_sim, plot = TRUE) #there are outliers at the ends

#try removing some 2-way interactions
length_mod_auto_3 <- glmmTMB::glmmTMB(max_length_norm ~ splines::ns(time_point,3) +
                                        treatment + 
                                        predator_treatment + 
                                        splines::ns(time_point,3):predator_treatment +
                                        splines::ns(time_point,3):treatment  +
                                        #  predator_treatment:treatment +
                                        #bit of code that deals with the autocorrelation
                                        #uncorrelated random intercept and correlated slope within group
                                        ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                      data = id_data, family = "gaussian", REML=F)

length_mod_auto_4 <- glmmTMB::glmmTMB(max_length_norm ~ splines::ns(time_point,3) +
                                        treatment + 
                                        predator_treatment + 
                                        splines::ns(time_point,3):predator_treatment +
                                        # splines::ns(time_point,3):treatment  +
                                        predator_treatment:treatment +
                                        #bit of code that deals with the autocorrelation
                                        #uncorrelated random intercept and correlated slope within group
                                        ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                      data = id_data, family = "gaussian", REML=F)

length_mod_auto_5 <- glmmTMB::glmmTMB(max_length_norm ~ splines::ns(time_point,3) +
                                        treatment + 
                                        predator_treatment + 
                                        # splines::ns(time_point,3):predator_treatment +
                                        splines::ns(time_point,3):treatment  +
                                        predator_treatment:treatment +
                                        #bit of code that deals with the autocorrelation
                                        #uncorrelated random intercept and correlated slope within group
                                        ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                      data = id_data, family = "gaussian", REML=F)

AIC(length_mod_auto, length_mod_auto_3, length_mod_auto_4, length_mod_auto_5)
#all 2-way interactions stay in 

summary(length_mod_auto)
Anova(length_mod_auto)

#extract model predictions for plotting
#newdata dataframe of unique variables required for plotting
plotting_data_l <- unique(id_data[,c('treatment', 'predator_treatment', 'replicate', 'time_point')])

#set the random effects of the model to zero
X_cond <- model.matrix(lme4::nobars(formula(length_mod_auto)[-2]), plotting_data_l)
beta_cond <- fixef(length_mod_auto)$cond
#predictions on the link scale generated from the dataframe with random effects set to zero
pred_cond <- X_cond %*% beta_cond

#transform point estimates of the speed to the response scale and multiply
ilink <- family(length_mod_auto)$linkinv
pred_length = ilink(pred_cond)

#set random number generator seed
set.seed(101)

#use posterior predictive simulations to generate upper and lower CIs
#and median predicted speeds
#ignoring variation in the random effects
#conditional
pred_condpar_psim = mvrnorm(1000, mu = beta_cond,Sigma = vcov(length_mod_auto)$cond)
pred_cond_psim = X_cond %*% t(pred_condpar_psim)

#transform
pred_length_psim = ilink(pred_cond_psim)

#calculate 95% CIs/
ci_length = t(apply(pred_length_psim,1, quantile, c(0.025, 0.975)))
ci_length = data.frame(ci_length)

#rename
names(ci_length) = c("length_low", "length_high")

#put into a data frame 
plotting_data_l = data.frame(plotting_data_l, pred_length, ci_length)

#plot on the graph 
ggplot(plotting_data_l, aes(x = time_point, y = pred_length, col = as.factor(predator_treatment)))+
  #  geom_point(data = id_data, aes(x = time_point, y = mean_roundness_norm), alpha = 0.1, pch = 1) +
  # geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point, y = pred_length)) +
  geom_ribbon(aes(x = time_point,
                  ymin = length_low,
                  ymax = length_high,
                  group = predator_treatment), alpha = 0.1, linetype = 0) +
  facet_grid(treatment~., scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean length")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

#########################################
######### roundness ID GLM ################
########################################

#fit roundness gaussian glm with three-way interaction
roundness_mod_gaussian <- glm(mean_roundness ~ time_point + treatment + predator_treatment +
                                   time_point:treatment + time_point:predator_treatment + treatment:predator_treatment +
                                   time_point:treatment:predator_treatment,
                                 data = id_data, family = "gaussian") 
#model output
summary(roundness_mod_gaussian)

#visual model checking
par(mfrow = c(2,2)) #for diagnostic plots
plot(roundness_mod_gaussian) #weird qq, s-l dipping in middle

#account for non-linearity of time by fitting a spline
#fit splines to account for non-linearity through time
roundness_mod_gaussian_spline <- glm(mean_roundness ~ splines::ns(time_point,2) + treatment + predator_treatment +
                                          splines::ns(time_point,2):treatment + splines::ns(time_point,2):predator_treatment + treatment:predator_treatment +
                                          splines::ns(time_point,2):treatment:predator_treatment,
                                        data = id_data, family = "gaussian")

plot(roundness_mod_gaussian_spline) #same problems as above

roundness_mod_gaussian_cbspline <- glm(mean_roundness ~ splines::ns(time_point,3) + treatment + predator_treatment +
                                            splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                                            splines::ns(time_point,3):treatment:predator_treatment,
                                          data = id_data, family = "gaussian")

plot(roundness_mod_gaussian_cbspline) #bad qq plot

AIC(roundness_mod_gaussian, roundness_mod_gaussian_spline, roundness_mod_gaussian_cbspline) #cubic spline best fit although diagnostics are worse

#normalise data using best normalize and run on cubic spline model
id_data$mean_roundness_norm <- predict(bestNormalize(id_data$mean_roundness))

roundness_mod_1 <- glm(mean_roundness_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                         splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment +
                         splines::ns(time_point,3):treatment:predator_treatment,
                       data = id_data, family = "gaussian")

plot(roundness_mod_1) #looks better but poss still not great

#model simplification by removing 3-way then 2-way interactions
roundness_mod_2 <- glm(mean_roundness_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                         splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment + treatment:predator_treatment,
                       data = id_data, family = "gaussian")
plot(roundness_mod_2)
AIC(roundness_mod_1, roundness_mod_2) #mod 2 is lower so we keep the 3-way interaction out and continue simplifying

roundness_mod_3 <- glm(mean_roundness_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                         splines::ns(time_point,3):treatment + splines::ns(time_point,3):predator_treatment, 
                       #treatment:predator_treatment +
                       #   splines::ns(time_point,3):treatment:predator_treatment,
                       data = id_data, family = "gaussian")

roundness_mod_4 <- glm(mean_roundness_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                         splines::ns(time_point,3):treatment + 
                         #splines::ns(time_point,3):predator_treatment + 
                         treatment:predator_treatment,
                       #   splines::ns(time_point,3):treatment:predator_treatment,
                       data = id_data, family = "gaussian")

roundness_mod_5 <- glm(mean_roundness_norm ~ splines::ns(time_point,3) + treatment + predator_treatment +
                         #splines::ns(time_point,3):treatment + 
                         splines::ns(time_point,3):predator_treatment + 
                         treatment:predator_treatment,
                       #   splines::ns(time_point,3):treatment:predator_treatment,
                       data = id_data, family = "gaussian")

AIC(roundness_mod_1, id_roundness_mod_2, roundness_mod_3, roundness_mod_4, roundness_mod_5) #all 2-way interactions stay

#extract model preds and plot 
id_data$roundness_preds <- predict(roundness_mod_2)
id_data$roundness_se <- 1.96*(predict(roundness_mod_2, se.fit = T)$se.fit)

ggplot(id_data, aes(x = time_point, y = mean_roundness_norm, col = as.factor(predator_treatment)))+
  #    geom_point(alpha = 0.1) +
  # geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point, y = roundness_preds)) +
  geom_ribbon(aes(x = time_point,
                  ymin = roundness_preds - roundness_se,
                  ymax = roundness_preds + roundness_se,
                  group = predator_treatment), alpha = 0.2, linetype = 0) +
  facet_grid(treatment~predator_treatment, scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean roundness")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))

##### roundness mixed effects model

#model with extra bits, spaced out for clarity
roundness_mod_auto <- glmmTMB::glmmTMB(mean_roundness_norm ~ splines::ns(time_point,3) +
                                         treatment + 
                                         predator_treatment + 
                                         splines::ns(time_point,3):predator_treatment +
                                         splines::ns(time_point,3):treatment  +
                                         predator_treatment:treatment +
                                         splines::ns(time_point,3):treatment:predator_treatment +
                                         #bit of code that deals with the autocorrelation
                                         #uncorrelated random intercept and correlated slope within group
                                         ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                       data = id_data, family = "gaussian", REML=F)
summary(roundness_mod_auto)
Anova(roundness_mod_auto)

#compare fit with other model
AIC(roundness_mod_auto, roundness_mod_1) #auto much lower AIC so continue to check fit 

#compare data to the distribution fit to it
fit_gaus <- fitdist(id_data$mean_roundness_norm,
                    distr = "norm")
plot(fit_gaus) #good fit

#residuals vs fitted
plot(fitted(roundness_mod_auto), residuals(roundness_mod_auto))

#compare the residuals of both models to check for autocorrelation
acf(resid(roundness_mod_auto))
acf(resid(roundness_mod_1))

#pseudo r2
MuMIn::r.squaredGLMM(roundness_mod_auto)

#simulate residuals from the model
roundness_glmm_sim <- simulateResiduals(roundness_mod_auto, n = 1000)
#plot the residuals
plot(roundness_glmm_sim)

testOutliers(roundness_glmm_sim, plot = TRUE) #there are outliers at the ends

#checking if model simplification improves fit by removing 3-way interaction
roundness_mod_auto_2 <- glmmTMB::glmmTMB(mean_roundness_norm ~ splines::ns(time_point,3) +
                                           treatment + 
                                           predator_treatment + 
                                           splines::ns(time_point,3):predator_treatment +
                                           splines::ns(time_point,3):treatment  +
                                           predator_treatment:treatment +
                                           #bit of code that deals with the autocorrelation
                                           #uncorrelated random intercept and correlated slope within group
                                           ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                         data = id_data, family = "gaussian", REML=F)
AIC(roundness_mod_auto, roundness_mod_auto_2) #better fit without 3-way interaction so we leave it out

#try removing some 2-way interactions
roundness_mod_auto_3 <- glmmTMB::glmmTMB(mean_roundness_norm ~ splines::ns(time_point,3) +
                                           treatment + 
                                           predator_treatment + 
                                           splines::ns(time_point,3):predator_treatment +
                                           splines::ns(time_point,3):treatment  +
                                           #  predator_treatment:treatment +
                                           #bit of code that deals with the autocorrelation
                                           #uncorrelated random intercept and correlated slope within group
                                           ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                         data = id_data, family = "gaussian", REML=F)

roundness_mod_auto_4 <- glmmTMB::glmmTMB(mean_roundness_norm ~ splines::ns(time_point,3) +
                                           treatment + 
                                           predator_treatment + 
                                           splines::ns(time_point,3):predator_treatment +
                                           # splines::ns(time_point,3):treatment  +
                                           predator_treatment:treatment +
                                           #bit of code that deals with the autocorrelation
                                           #uncorrelated random intercept and correlated slope within group
                                           ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                         data = id_data, family = "gaussian", REML=F)

roundness_mod_auto_5 <- glmmTMB::glmmTMB(mean_roundness_norm ~ splines::ns(time_point,3) +
                                           treatment + 
                                           predator_treatment + 
                                           # splines::ns(time_point,3):predator_treatment +
                                           splines::ns(time_point,3):treatment  +
                                           predator_treatment:treatment +
                                           #bit of code that deals with the autocorrelation
                                           #uncorrelated random intercept and correlated slope within group
                                           ar1(as.factor(time_point) + 0 | replicate) + (1|replicate), 
                                         data = id_data, family = "gaussian", REML=F)

MuMIn::model.sel(roundness_mod_auto, roundness_mod_auto_2, roundness_mod_auto_3, roundness_mod_auto_4, roundness_mod_auto_5)
#all 2-way interactions stay in 

summary(roundness_mod_auto_2)
Anova(roundness_mod_auto_2)

#extract model predictions for plotting
#newdata dataframe of unique variables required for plotting
plotting_data_r <- unique(id_data[,c('treatment', 'predator_treatment', 'replicate', 'time_point')])

#set the random effects of the model to zero
X_cond <- model.matrix(lme4::nobars(formula(roundness_mod_auto_2)[-2]), plotting_data_r)
beta_cond <- fixef(roundness_mod_auto_2)$cond
#predictions on the link scale generated from the dataframe with random effects set to zero
pred_cond <- X_cond %*% beta_cond

#transform point estimates of the speed to the response scale and multiply
ilink <- family(roundness_mod_auto_2)$linkinv
pred_roundness = ilink(pred_cond)

#set random number generator seed
set.seed(101)

#use posterior predictive simulations to generate upper and lower CIs
#and median predicted speeds
#ignoring variation in the random effects
#conditional
pred_condpar_psim = mvrnorm(1000, mu = beta_cond,Sigma = vcov(roundness_mod_auto_2)$cond)
pred_cond_psim = X_cond %*% t(pred_condpar_psim)

#transform
pred_roundness_psim = ilink(pred_cond_psim)

#calculate 95% CIs/
ci_roundness = t(apply(pred_roundness_psim,1, quantile, c(0.025, 0.975)))
ci_roundness = data.frame(ci_roundness)

#rename
names(ci_roundness) = c("roundness_low", "roundness_high")

#put into a data frame 
plotting_data_r = data.frame(plotting_data_r, pred_roundness, ci_roundness)

#plot on the graph 
ggplot(plotting_data_r, aes(x = time_point, y = pred_roundness, col = as.factor(predator_treatment)))+
  #  geom_point(data = id_data, aes(x = time_point, y = mean_roundness_norm), alpha = 0.1, pch = 1) +
  # geom_smooth(se = FALSE) +
  geom_line(aes(x = time_point, y = pred_roundness)) +
  geom_ribbon(aes(x = time_point,
                  ymin = roundness_low,
                  ymax = roundness_high,
                  group = predator_treatment), alpha = 0.1, linetype = 0) +
  facet_grid(treatment~., scales = "fixed")+
  xlab("Time (hours)")+
  ylab("Mean roundness")+
  labs(colour = "Predator treatment") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "black"))





