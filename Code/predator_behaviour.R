################################################################
# Predator Behaviour
################################################################

did_id_data <- readRDS(file = "Data/didinium_data_IDs_clean.RDS")
hom_id_data <- readRDS(file = "Data/homalozoon_data_IDs_clean.RDS")

#add predator treatment column to each then combine data frames
did_id_data <- did_id_data %>%
  mutate(predator_treatment = "didinium")

hom_id_data <- hom_id_data %>%
  mutate(predator_treatment = "homalozoon")

#combine the data frames
id_predators <- rbind(did_id_data, hom_id_data) %>%
  #set prey as the first factor for data analysis
  mutate(predator_treatment = factor(predator_treatment, levels = c("prey", "didinium", "homalozoon")))%>%
  mutate(
    treatment = factor(treatment, levels = c(15,25)),
    treat_inter = as.factor(interaction(treatment,predator_treatment))) %>%
  group_by(Species,treatment,predator_treatment,replicate,treat_inter,time_point) %>%
  summarise(across(c(max_abundance:mean_speed), ~mean(.x))) %>%
  mutate(Species = case_when(Species == "DIDnas" ~ "Didinium",
                             Species == "HOMver" ~ "Homalozoon"))

ggplot(data = subset(id_predators,Species != "PARcau"), 
       aes(x=treatment,y=mean_speed,fill=treatment)) +
  geom_boxplot()+
  scale_x_discrete(labels = paste0(unique(id_predators$treatment),"\u00B0C"))+
  facet_grid(~Species)+
  xlab("Temperature treatment")+
  ylab("Mean speed (mm/s)")+
  scale_fill_manual(name = "Treatment",values = c("#a2d7d8","#de5842")) + 
  ggpubr::stat_compare_means( method = "wilcox.test",
                             symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns")))+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        strip.text.y.right = element_text(angle = 0),
        strip.text.x = element_text(face = "bold.italic"),
        strip.text = element_text(face = "bold"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))
