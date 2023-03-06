########################################################################################################################
# Supplementary Figures #
########################################################################################################################
require(patchwork)
require(brms)
require(ggplot2)
brms_id_weibullk8 <- readRDS(file = "Results/brms_id_noar_weibullk8.rds")
brms_lm_id <- readRDS(file = "Results/brms_lm_id_noar.rds")

############ 
#Extract model trace plots
############ 

pdf(file = "Results/figures/supplementary_figures/figure_S1a.pdf",
    width = 6,height = 15)
plot(brms_id_weibullk8,plot = TRUE,theme = theme_bw(),variable =c("^b_","^bs_"),regex = TRUE, N=32,newpage = FALSE)
dev.off()

pdf(file = "Results/figures/supplementary_figures/figure_S1b.pdf",
    width = 6,height = 12)
plot(brms_id_weibullk8,plot = TRUE,theme = theme_bw(),variable =c("^sds_","^sd_","shape"),regex = TRUE, N=14,newpage = FALSE)
dev.off()

pdf(file = "Results/figures/supplementary_figures/figure_S2a.pdf",
    width = 6,height = 15)
plot(brms_lm_id,plot = TRUE,theme = theme_bw(),variable = paste0("b_",rownames(fixef(brms_lm_id))[1:13]),regex = FALSE, N=31,newpage = FALSE)
dev.off()

pdf(file = "Results/figures/supplementary_figures/figure_S2b.pdf",
    width = 6,height = 15)
plot(brms_lm_id,plot = TRUE,theme = theme_bw(),variable =c(paste0("b_",rownames(fixef(brms_lm_id))[14:26]),"sd_replicate__Intercept","sigma"),regex = TRUE, N=13,newpage = FALSE)
dev.off()

############ 
#Plot posterior predictive checks
############ 

ggsave(
  brms::pp_check(brms_id_weibullk8) + 
    theme_bw(),
  filename = "Results/figures/supplementary_figures/figure_S3.pdf",
  width = 4,height = 4)

ggsave(
  brms::pp_check(brms_lm_id) + 
    theme_bw(),
  filename = "Results/figures/supplementary_figures/figure_S4.pdf",
  width = 4,height = 4)
