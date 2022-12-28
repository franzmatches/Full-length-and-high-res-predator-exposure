require(brms)
require(Matrix)
require(tidyverse)
require(camcorder)

brms_lm <- readRDS("Results/brms_lm.rds")

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
  mutate(treatment = factor(treatment, levels = c(15,25)),
         treat_inter = as.factor(interaction(treatment,predator_treatment)))

cond_inter_treatment2 <- conditional_effects(brms_lm,effects = "mean_length_um:treatment",
                                             conditions = data.frame(expand.grid(time_point = unique(id_data$time_point),
                                                                                 predator_treatment =  c("prey","didinium","homalozoon"))))[[1]]|>
  mutate(cond__ = paste(time_point,predator_treatment,sep="_"),
         effect = "interaction")

camcorder::gg_record(
  dir = file.path("/Users/ul20791/Downloads/camcorder_test", "recording100"), # where to save the recording
  device = "png", # device to use to save images
  width = 5,      # width of saved image
  height = 2,     # height of saved image
  units = "in",   # units for width and height
  dpi = 300       # dpi to use when saving image
)

for(i in unique(cond_inter_treatment2$time_point)){
  print(
    ggplot(cond_inter_treatment2 |>
           mutate(cond__ = paste0(time_point,"h")) |>
           filter(cond__ == paste0(i,"h")) |>
      mutate(predator_treatment = factor(predator_treatment,labels =  c("Control","Didinium","Homalozoon"))))+
    geom_line(aes(x = effect1__, y=estimate__,col=as.factor(effect2__))) +
    geom_ribbon(aes(x = effect1__,ymin = lower__, ymax =  upper__,fill=as.factor(effect2__)),alpha=0.5)+
    facet_grid(cond__~predator_treatment, scales = "fixed") +
    scale_fill_manual(name = "Temperature \u00B0C",values = c("#a2d7d8","#de5842")) + 
    scale_color_manual(name = "Temperature \u00B0C",values = c("#a2d7d8","#de5842")) + 
    xlab("Mean length (\u03BCm)")+
    ylab("Mean width (\u03BCm)")+
    coord_cartesian(ylim = c(50,140), xlim = c(50,400))+
    theme_classic()+
    theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
          strip.text = element_text(face = "bold"),
          strip.text.y.right = element_text(angle = 0),
          strip.text.x = element_text(face = "bold.italic"),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, colour = "black"))
  )
}
camcorder::gg_playback(
  name = file.path("/Users/ul20791/Downloads/camcorder_test","vignette_gif.gif"),
  first_image_duration = 5,
  last_image_duration = 5,
  frame_duration = 0.25,
  last_as_first = FALSE
)

camcorder::gg_stop_recording()
