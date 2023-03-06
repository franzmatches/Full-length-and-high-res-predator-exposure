########################################################################################################################
# Supplementary Tables #
########################################################################################################################
require(tidyverse)
brms_id_weibullk8 <- readRDS(file = "Results/brms_id_noar_weibullk8.rds")
brms_lm_id <- readRDS(file = "Results/brms_lm_id_noar.rds")

############ 
#Extract summary tables
############ 

speed_model_summary_table <- data.frame(rbind(summary(brms_id_weibullk8)$random[[1]],
                                              summary(brms_id_weibullk8)$fixed,
                                              summary(brms_id_weibullk8)$splines)) %>%
  tibble::rownames_to_column("Parameter") %>%
  dplyr::mutate(dplyr::across(Estimate:Tail_ESS,~signif(.,digit=3))) %>%
  dplyr::add_row(Parameter = "random",.before = 1) %>%
  dplyr::add_row(Parameter = "fixed",.before = 3) %>%
  dplyr::add_row(Parameter = "splines",.before = 22) 

shape_model_summary_table <- data.frame(rbind(summary(brms_lm_id)$random[[1]],
                                              summary(brms_lm_id)$fixed,
                                              summary(brms_lm_id)$splines)) %>%
  tibble::rownames_to_column("Parameter") %>%
  dplyr::mutate(dplyr::across(Estimate:Tail_ESS,~signif(.,digit=3))) %>%
  dplyr::add_row(Parameter = "random",.before = 1) %>%
  dplyr::add_row(Parameter = "fixed",.before = 3)

############ 
#Save tables
############ 

readr::write_csv(speed_model_summary_table,
                 file = "Results/figures/supplementary_figures/speed_model_summary_table.csv")

readr::write_csv(shape_model_summary_table,
                 file = "Results/figures/supplementary_figures/shape_model_summary_table.csv")
