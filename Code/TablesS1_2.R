########################################################################################################################
# Supplementary Tables #
########################################################################################################################

require(tidyverse)

#load models
brms_id_weibullk8 <- readRDS(file = "Results/brms_id_noar_weibullk8.rds")
brms_lm_id <- readRDS(file = "Results/brms_lm_id_noar.rds")

############ 
#Extract summary tables
############ 

speed_model_summary_table <- data.frame(rbind(summary(brms_id_weibullk8)$random[[1]],
                                              summary(brms_id_weibullk8)$fixed,
                                              summary(brms_id_weibullk8)$splines)) %>%
  dplyr::mutate(dplyr::across(Estimate:Tail_ESS,~as.numeric(.))) %>% #remove tibble's pillar_number formatting
  tibble::rownames_to_column("Parameter") %>%
  dplyr::mutate(across(Estimate:u.95..CI,~round(.,digits =3))) %>%
  dplyr::mutate(across(Rhat:Tail_ESS,~round(.,digits = 2))) %>%
  dplyr::add_row(Parameter = "random",.before = 1) %>% #format table into sections (random effects, fixed effects & splines)
  dplyr::add_row(Parameter = "fixed",.before = 3) %>%
  dplyr::add_row(Parameter = "splines",.before = 22) 

shape_model_summary_table <- data.frame(rbind(summary(brms_lm_id)$random[[1]],
                                              summary(brms_lm_id)$fixed,
                                              summary(brms_lm_id)$splines)) %>%
  dplyr::mutate(dplyr::across(Estimate:Tail_ESS,~as.numeric(.))) %>% #remove tibble's pillar_number formatting
  tibble::rownames_to_column("Parameter") %>%
  dplyr::mutate(across(Estimate:u.95..CI,~round(.,digits =3))) %>%
  dplyr::mutate(across(Rhat:Tail_ESS,~round(.,digits = 2))) %>%
  dplyr::add_row(Parameter = "random",.before = 1) %>%
  dplyr::add_row(Parameter = "fixed",.before = 3)

############ 
#Save tables
############ 

write.csv(speed_model_summary_table,
                 file = "Results/figures/supplementary_figures/speed_model_summary_table.csv", row.names = FALSE)

write.csv(shape_model_summary_table,
                 file = "Results/figures/supplementary_figures/shape_model_summary_table.csv", row.names = FALSE)
