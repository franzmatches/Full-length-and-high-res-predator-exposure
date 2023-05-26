################################################################
# Figure S1
################################################################
library(tidyverse)

did_id_data <- readRDS(file = "Data/didinium_data_IDs_clean.RDS")
hom_id_data <- readRDS(file = "Data/homalozoon_data_IDs_clean.RDS")
prey_id_data<-readRDS(file = "Data/prey_data_IDs_clean.RDS")
#add predator treatment column to each then combine data frames
did_id_data <- did_id_data %>%
  mutate(predator_treatment = "D. nasutum")

hom_id_data <- hom_id_data %>%
  mutate(predator_treatment = "H. vermiculare")

prey_id_data <- prey_id_data %>%
  mutate(predator_treatment = "Control")

#combine the data frames
id_data <- rbind(did_id_data, hom_id_data, prey_id_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment,
                                     levels = c("Control", "D. nasutum", "H. vermiculare")))

#plotting abundance

ggsave(filename = "Results/figures/supplementary_figures/figure_S1.png",
ggplot()+
  # geom_point(data = id_data, aes(x = time_point, y = mean_speed, group = replicate, col = as.factor(predator_treatment)), alpha = .05)+
  geom_line(data = id_data %>%
              # filter(Species == "PARcau") %>% 
              group_by(time_point, treatment, predator_treatment, replicate) %>% 
              summarise(max_abundance = max(max_abundance)),
            aes(x = time_point, y = max_abundance,
                group = replicate, colour = treatment), alpha = .3)+
  geom_line(data = id_data %>% group_by(time_point, treatment, predator_treatment) %>% 
              summarise(mean_max_abundance = mean(max_abundance)),
            aes(x = time_point, y = mean_max_abundance, color = treatment), size = .7)+
  # geom_smooth(se = FALSE) +
  # geom_jitter(width = 1) +
  facet_grid(predator_treatment~treatment)+
  scale_color_manual(name = "treatment",values = c("#a2d7d8","#de5842")) + 
  xlab("Time (hours)")+
  ylab("Abundance")+
  guides(colour = FALSE) +
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text = element_text(face = "bold"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black")))


