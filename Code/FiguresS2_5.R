########################################################################################################################
# Supplementary Figures 2-5 #
########################################################################################################################
require(patchwork)
require(brms)
require(ggplot2)
brms_id_weibullk8 <- readRDS(file = "Results/brms_id_noar_weibullk8.rds")
brms_lm_id <- readRDS(file = "Results/brms_lm_id_noar.rds")

############ 
#Extract model trace plots
############ 

png(file = "Results/figures/supplementary_figures/figure_S2a.png",
    width = 6,height = 15,units = "in",res=300)
plot(brms_id_weibullk8,plot = TRUE,theme = theme_bw(),variable =c("^b_","^bs_"),regex = TRUE, N=32,newpage = FALSE)
dev.off()

png(file = "Results/figures/supplementary_figures/figure_S2b.png",
    width = 6,height = 12,units = "in",res=300)
plot(brms_id_weibullk8,plot = TRUE,theme = theme_bw(),variable =c("^sds_","^sd_","shape"),regex = TRUE, N=14,newpage = FALSE)
dev.off()

png(file = "Results/figures/supplementary_figures/figure_S3a.png",
    width = 6,height = 15,units = "in",res=300)
plot(brms_lm_id,plot = TRUE,theme = theme_bw(),variable = paste0("b_",rownames(fixef(brms_lm_id))[1:13]),regex = FALSE, N=31,newpage = FALSE)
dev.off()

png(file = "Results/figures/supplementary_figures/figure_S3b.png",
    width = 6,height = 15,units = "in",res=300)
plot(brms_lm_id,plot = TRUE,theme = theme_bw(),variable =c(paste0("b_",rownames(fixef(brms_lm_id))[14:26]),"sd_replicate__Intercept","sigma"),regex = TRUE, N=13,newpage = FALSE)
dev.off()

############ 
#Plot posterior predictive checks
############ 

ggsave(
  brms::pp_check(brms_id_weibullk8) + 
    theme_bw(),
  filename = "Results/figures/supplementary_figures/figure_S4.png",
  width = 4,height = 4)

ggsave(
  brms::pp_check(brms_lm_id) + 
    theme_bw(),
  filename = "Results/figures/supplementary_figures/figure_S5.png",
  width = 4,height = 4)
